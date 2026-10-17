"""Explicit, resumable installation of pinned human references and native Bowtie tools.

Nothing downloads during import, app installation, app updates, or status reads. A completed
release is activated atomically; cached downloads survive interrupted preparation and upgrades.
"""

import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time
import zipfile
from collections.abc import Iterator
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Any, BinaryIO
from urllib.parse import unquote

import httpx

RELEASE = "GRCh38.p14-RefSeq110-v1"
PREFIX = "GCF_000001405.40_GRCh38.p14"
SOURCE = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/{PREFIX}/"
ASSETS = {
    "genomic": ("_genomic.fna.gz", "c30471567037b2b2389d43c908c653e1", 972898531),
    "transcript": ("_rna.fna.gz", "a3306ec93f9b166dd370fe436d671db1", 132922439),
    "annotation": ("_genomic.gff.gz", "d259315d07ca2ecf1d080d7d640cc1e5", 56923273),
    "assembly": ("_assembly_report.txt", "073e01acec1c43398074c61e3628bf45", 79219),
}
TOOLS = {
    "Linux": (
        "1.3.1",
        "https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/"
        "bowtie-1.3.1-linux-x86_64.zip",
        "e538cf41bc1661981ea19dd187d5feea50a48664e6c879bd0fd2bcde80aa91ca",
        52387989,
    ),
    "Windows": (
        "1.2.2",
        "https://github.com/BenLangmead/bowtie/releases/download/v1.2.2_p1/"
        "bowtie-1.2.2-mingw-x86_64.zip",
        "021975c3efc51c5dbfbc001b8a1f13395be7e6ea272592740b476c8056670652",
        31257451,
    ),
}
TSS_PROXY_FEATURES = frozenset({"mRNA", "lnc_RNA", "transcript", "primary_transcript"})

INDEX_BUILD_OPTIONS = {
    "genomic": ("--packed", "--bmaxdivn", "8", "--threads", "2"),
    "transcript": ("--bmaxdivn", "4", "--threads", "2"),
}

INDEX_SUFFIXES = (".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    """Replace a small state file without exposing a partial JSON document."""
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value), encoding="utf-8")
    temporary.replace(path)


def digest(path: Path, algorithm: str) -> str:
    """Hash a file in bounded memory."""
    value = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


@contextmanager
def installation_lock(root: Path) -> Iterator[None]:
    """Hold an OS lock, released even when the app crashes, including on Windows."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / "setup.lock").open("a+b") as stream:
        stream.seek(0, os.SEEK_END)
        if stream.tell() == 0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        try:
            if sys.platform == "win32":
                import msvcrt

                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise ValueError(
                "Reference setup is already running in another app process."
            ) from error
        try:
            yield
        finally:
            stream.seek(0)
            if sys.platform == "win32":
                import msvcrt

                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


class ReferenceManager:
    """Own one background setup task and a durable, versioned reference directory."""

    def __init__(self, root: Path) -> None:
        """Record paths; do not create directories or start any network activity."""
        self.root = root
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._guard = threading.Lock()
        self._state: dict[str, Any] = {}
        self.search_lock = threading.Lock()
        self.shutdown_event = threading.Event()

    def active(self, release: str | None = None) -> tuple[Path, dict[str, Any]]:
        """Return only a fully prepared installation, never a staging directory."""
        try:
            if release is None:
                pointer = json.loads((self.root / "active.json").read_text(encoding="utf-8"))
                if not isinstance(pointer, dict) or not isinstance(pointer.get("release"), str):
                    raise ValueError("Reference activation metadata is invalid. Run setup again.")
                release = pointer["release"]
            if (
                release in (".", "..")
                or release != Path(release).name
                or release.endswith(".preparing")
            ):
                raise ValueError("Invalid reference installation path.")
            folder = self.root / release
            manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
            if (
                not isinstance(manifest, dict)
                or manifest.get("release") != release
                or not isinstance(manifest.get("files"), dict)
                or not all(isinstance(name, str) for name in manifest["files"])
            ):
                raise ValueError("Reference manifest is invalid. Run setup again.")
            expected_platform = platform.system() + "-" + platform.machine()
            if manifest.get("platform", expected_platform).lower() != expected_platform.lower():
                raise ValueError(
                    "The installed search tools belong to another platform. "
                    "Prepare references on this computer."
                )
            for name, size in manifest["files"].items():
                path = folder / name
                if not path.is_file() or path.stat().st_size != size:
                    raise ValueError("Reference files are missing or changed. Run setup again.")
            return folder, manifest
        except (OSError, KeyError, json.JSONDecodeError) as error:
            raise ValueError("Set up the local human references before searching.") from error

    def status(self) -> dict[str, Any]:
        """Return local metadata only, without hashing gigabytes or contacting a provider."""
        installation_issue: str | None = None
        try:
            _, manifest = self.active()
        except ValueError as error:
            manifest = None
            installation_issue = str(error)
        else:
            installation_issue = None
        state = dict(self._state)
        if not state:
            try:
                state = json.loads((self.root / "setup.json").read_text(encoding="utf-8"))
                if state.get("status") == "running":
                    state = {
                        **state,
                        "status": "paused",
                        "message": "Setup was interrupted. Resume to continue.",
                    }
            except (OSError, json.JSONDecodeError):
                state = {"status": "not-installed"}
        if manifest is None and state.get("status") == "ready":
            state = {"status": "error", "message": installation_issue}
        elif manifest is not None and state.get("status") == "not-installed":
            state = {"status": "ready", "message": "Installed references will be reused."}
        supported = platform.system() in TOOLS and platform.machine().lower() in ("amd64", "x86_64")
        if state.get("status") == "running" and state.get("phase_started_at"):
            state["elapsed_seconds"] = int(time.time() - state["phase_started_at"])
        return {
            **state,
            "ready": manifest is not None,
            "installed": manifest,
            "catalog_release": RELEASE,
            "update_available": bool(manifest and manifest["release"] != RELEASE),
            "supported": supported,
            "download_bytes": sum(a[2] for a in ASSETS.values())
            + (TOOLS[platform.system()][3] if supported else 0),
            "required_free_bytes": 20 * 1024**3,
            "directory": str(self.root),
        }

    def start(self) -> None:
        """Start or resume setup explicitly; a ready current release needs no download."""
        with self._guard:
            if self._thread and self._thread.is_alive():
                return
            if not self.status()["supported"]:
                raise ValueError(
                    "Automatic local setup currently supports Windows x64 and Linux x64."
                )
            if self.status()["ready"] and not self.status()["update_available"]:
                return
            self.root.mkdir(parents=True, exist_ok=True)
            self._stop.clear()
            self._state = {"status": "running", "message": "Checking local storage"}
            self._thread = threading.Thread(
                target=self._setup, name="sarna-reference-setup", daemon=True
            )
            self._thread.start()

    def pause(self) -> None:
        """Request a download checkpoint or stop the current index build safely."""
        self._stop.set()

    def close(self) -> None:
        """Stop setup and local searches when the app closes, keeping verified input files."""
        self.pause()
        self.shutdown_event.set()
        if self._thread and self._thread.is_alive():
            self._thread.join(timeout=2)

    def _check(self) -> None:
        if self._stop.is_set():
            raise InterruptedError(
                "Setup paused. Verified downloads are kept; resume when convenient."
            )

    def _report(self, message: str, **detail: str | int) -> None:
        self._state = {
            "status": "running",
            "message": message,
            "phase_started_at": time.time(),
            **detail,
        }
        atomic_json(self.root / "setup.json", self._state)

    def _download(self, url: str, expected: str, size: int, algorithm: str) -> Path:
        cache = self.root / "downloads"
        cache.mkdir(exist_ok=True)
        destination = cache / expected
        if destination.exists():
            self._report("Verifying cached download", file=url.rsplit("/", 1)[-1])
            if destination.stat().st_size == size and digest(destination, algorithm) == expected:
                return destination
            destination.unlink()
        part = destination.with_suffix(".part")
        offset = part.stat().st_size if part.exists() else 0
        if offset > size:
            part.unlink()
            offset = 0
        self._report(
            "Downloading reference data", file=url.rsplit("/", 1)[-1], downloaded=offset, total=size
        )
        if offset < size:
            with (
                httpx.Client(
                    follow_redirects=True, timeout=httpx.Timeout(60, connect=30)
                ) as client,
                client.stream(
                    "GET", url, headers={"Range": f"bytes={offset}-", "Accept-Encoding": "identity"}
                ) as reply,
            ):
                reply.raise_for_status()
                if reply.status_code == 206:
                    expected_range = f"bytes {offset}-{size - 1}/{size}"
                    if reply.headers.get("content-range") != expected_range:
                        raise ValueError("The download server returned an unexpected byte range.")
                else:
                    offset = 0
                updated = time.monotonic()
                if not offset:
                    part.unlink(missing_ok=True)
                with part.open("ab") as output:
                    for chunk in reply.iter_bytes(1024 * 1024):
                        self._check()
                        output.write(chunk)
                        offset += len(chunk)
                        if offset > size:
                            raise ValueError("Download exceeds the pinned reference size.")
                        if time.monotonic() - updated >= 1:
                            self._report(
                                "Downloading reference data",
                                file=url.rsplit("/", 1)[-1],
                                downloaded=offset,
                                total=size,
                            )
                            updated = time.monotonic()
        self._check()
        self._report("Verifying download checksum", file=url.rsplit("/", 1)[-1])
        if part.stat().st_size != size or digest(part, algorithm) != expected:
            part.unlink(missing_ok=True)
            raise ValueError("Download checksum failed. Retry setup to download a verified copy.")
        part.replace(destination)
        return destination

    def command(self, args: list[str], log: Path) -> None:
        """Run a native tool without a shell, with cancellation and a retained diagnostic log."""
        with (
            log.open("wb") as output,
            subprocess.Popen(args, stdout=output, stderr=subprocess.STDOUT) as process,
        ):
            try:
                while process.poll() is None:
                    self._check()
                    self._stop.wait(0.25)
            except BaseException:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                raise
            if process.returncode:
                raise ValueError(f"Local reference tool failed; see {log}.")

    def _setup(self) -> None:
        try:
            with installation_lock(self.root):
                self._prepare()
        except Exception as error:
            self._state = {
                "status": "paused" if isinstance(error, InterruptedError) else "error",
                "message": str(error),
            }
            atomic_json(self.root / "setup.json", self._state)

    def _prepare(self) -> None:
        if shutil.disk_usage(self.root).free < 20 * 1024**3:
            raise ValueError("Reference preparation needs at least 20 GiB of free storage.")
        version, url, checksum, size = TOOLS[platform.system()]
        archive = self._download(url, checksum, size, "sha256")
        stage = self.root / (RELEASE + ".preparing")
        stage.mkdir(exist_ok=True)
        self._report("Preparing verified native search tools")
        with zipfile.ZipFile(archive) as zipped:
            for member in zipped.infolist():
                self._check()
                path = Path(member.filename)
                if (
                    path.is_absolute()
                    or path.drive
                    or ".." in path.parts
                    or "\\" in member.filename
                ):
                    raise ValueError("Unsafe path in tool archive.")
                if not member.is_dir() and not member.filename.startswith("__MACOSX/"):
                    destination = stage / "tools" / path
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with zipped.open(member) as source, destination.open("wb") as output:
                        shutil.copyfileobj(source, output)
                    if os.name != "nt":
                        destination.chmod(0o755)
        extension = ".exe" if os.name == "nt" else ""
        aligner = next((stage / "tools").rglob("bowtie-align-s" + extension))
        builder = next((stage / "tools").rglob("bowtie-build-s" + extension))
        self.command([str(aligner), "--version"], stage / "tool-version.log")
        self._tool_check(aligner, builder, stage)
        assets = {
            key: self._download(SOURCE + PREFIX + suffix, md5, length, "md5")
            for key, (suffix, md5, length) in ASSETS.items()
        }
        # Rebuild derived data after interruption; verified compressed inputs are always reused.
        database = stage / "annotation.sqlite3"
        database.unlink(missing_ok=True)
        with closing(sqlite3.connect(database)) as connection, connection:
            connection.executescript("""
                CREATE TABLE sequences(scope TEXT,accession TEXT,title TEXT,offset INTEGER,
                  bases INTEGER,width INTEGER,length INTEGER,PRIMARY KEY(scope,accession));
                CREATE TABLE transcripts(accession TEXT,contig TEXT,tss INTEGER,
                  strand TEXT,gene TEXT,feature_type TEXT,is_tss_proxy INTEGER);
                CREATE INDEX transcript_accession ON transcripts(accession);
                CREATE INDEX transcript_tss ON transcripts(contig,is_tss_proxy,tss);
                CREATE TABLE contigs(accession TEXT PRIMARY KEY,label TEXT);
            """)
            for scope in ("genomic", "transcript"):
                self._report(f"Preparing {scope} sequences (one-time setup)")
                self._fasta(assets[scope], stage / f"{scope}.fa", scope, connection)
            self._report("Preparing local gene and transcription-start annotations")
            with assets["assembly"].open(encoding="utf-8") as report:
                for line in report:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip().split("\t")
                    if len(fields) >= 10:
                        connection.execute(
                            "INSERT OR REPLACE INTO contigs VALUES(?,?)",
                            (
                                fields[6],
                                fields[2]
                                if fields[1] == "assembled-molecule" and fields[2] != "na"
                                else fields[6],
                            ),
                        )
            self._transcripts(assets["annotation"], connection)
        with self.search_lock:
            for scope in ("genomic", "transcript"):
                self._check()
                self._report(
                    f"Building {scope} search index (one-time setup; this can take several minutes)"
                )
                index = stage / scope
                marker = stage / f"{scope}.index-ready"
                settings = {
                    "source_md5": ASSETS[scope][1],
                    "tool_sha256": checksum,
                    "options": list(INDEX_BUILD_OPTIONS[scope]),
                    "fasta_policy": "upper-60-v1",
                }
                try:
                    completed = json.loads(marker.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    completed = None
                if completed != settings or not all(
                    Path(str(index) + suffix).is_file() for suffix in INDEX_SUFFIXES
                ):
                    for suffix in INDEX_SUFFIXES:
                        Path(str(index) + suffix).unlink(missing_ok=True)
                    marker.unlink(missing_ok=True)
                    self.command(
                        [
                            str(builder),
                            *INDEX_BUILD_OPTIONS[scope],
                            str(stage / f"{scope}.fa"),
                            str(index),
                        ],
                        stage / f"build-{scope}.log",
                    )
                    atomic_json(marker, settings)
                self._check()
            self._report("Checking prepared indexes")
            for scope in ("genomic", "transcript"):
                query = ""
                with (stage / f"{scope}.fa").open(encoding="utf-8") as sequence:
                    for line in sequence:
                        self._check()
                        if not line.startswith(">") and (match := re.search(r"[ACGT]{19}", line)):
                            query = match.group()
                            break
                if not query:
                    raise ValueError("The reference contains no callable 19-base test sequence.")
                self.command(
                    [str(aligner), "-c", "-v", "0", "-k", "1", str(stage / scope), query],
                    stage / f"check-{scope}.log",
                )
                if not any(
                    len(line.split("\t")) >= 8
                    for line in (stage / f"check-{scope}.log")
                    .read_text(encoding="utf-8")
                    .splitlines()
                ):
                    raise ValueError("Prepared index failed its known-reference positive control.")
        self._check()
        manifest = {
            "release": RELEASE,
            "assembly": "GRCh38.p14",
            "accession": "GCF_000001405.40",
            "annotation": "NCBI RefSeq release 110",
            "source": SOURCE,
            "bowtie_version": version,
            "tool_sha256": checksum,
            "tool_url": url,
            "platform": platform.system() + "-" + platform.machine(),
            "index_build": {
                "program": "bowtie-build-s",
                "options_by_scope": {
                    scope: list(options) for scope, options in INDEX_BUILD_OPTIONS.items()
                },
                "format": "small-ebwt",
            },
            "aligner": str(aligner.relative_to(stage)),
            "sources": {
                key: {
                    "url": SOURCE + PREFIX + ASSETS[key][0],
                    "md5": ASSETS[key][1],
                    "sha256": digest(path, "sha256"),
                }
                for key, path in assets.items()
            },
            "files": {
                str(p.relative_to(stage)): p.stat().st_size
                for p in stage.rglob("*")
                if p.is_file() and p.name != "manifest.json"
            },
        }
        manifest["prepared_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        atomic_json(stage / "manifest.json", manifest)
        destination = self.root / RELEASE
        if destination.exists():
            # A damaged installation is preserved for inspection, replaced only after preparation.
            destination.replace(self.root / (RELEASE + f".replaced-{time.time_ns()}"))
        stage.replace(destination)
        atomic_json(self.root / "active.json", {"release": RELEASE})
        self._state = {
            "status": "ready",
            "message": "Local human references are ready. Downloads and indexes will be reused.",
        }
        atomic_json(self.root / "setup.json", self._state)

    def _transcripts(self, source: Path, db: sqlite3.Connection) -> None:
        """Map RNA accessions, separately marking complete transcript models as TSS proxies."""
        with gzip.open(source, "rt", encoding="utf-8") as annotation:
            for line_number, line in enumerate(annotation):
                if line_number % 10000 == 0:
                    self._check()
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 9 or not (
                    fields[2].endswith("RNA") or fields[2] in ("transcript", "primary_transcript")
                ):
                    continue
                attributes = dict(
                    item.split("=", 1) for item in fields[8].split(";") if "=" in item
                )
                accession = unquote(attributes.get("transcript_id", ""))
                if accession:
                    db.execute(
                        "INSERT INTO transcripts VALUES(?,?,?,?,?,?,?)",
                        (
                            accession,
                            fields[0],
                            int(fields[3]) if fields[6] == "+" else int(fields[4]),
                            fields[6],
                            unquote(attributes.get("gene", attributes.get("Name", ""))),
                            fields[2],
                            int(
                                fields[2] in TSS_PROXY_FEATURES
                                and attributes.get("partial") != "true"
                                and fields[6] in ("+", "-")
                            ),
                        ),
                    )

    def _tool_check(self, aligner: Path, builder: Path, stage: Path) -> None:
        """Check native tools against a tiny exhaustive oracle before large downloads."""
        self._report("Checking native tool compatibility and search correctness")
        query = "GTGAAGCCATAGACTCGTA"
        complement = str.maketrans("ACGTN", "TGCAN")
        variants = [query, query.translate(complement)[::-1]]
        for mismatches in (1, 2, 3, 4):
            changed = list(query)
            for position in (2, 7, 12, 17)[:mismatches]:
                changed[position] = next(base for base in "ACGT" if base != changed[position])
            variants.append("".join(changed))
        reference = ("N" * 25).join(variants)
        expected = set()
        for position in range(len(reference) - len(query) + 1):
            segment = reference[position : position + len(query)]
            if set(segment) - set("ACGT"):
                continue
            for strand, oriented in (("+", segment), ("-", segment.translate(complement)[::-1])):
                if sum(a != b for a, b in zip(query, oriented, strict=True)) <= 3:
                    expected.add((strand, position))
        with tempfile.TemporaryDirectory(prefix="tool-check-", dir=stage) as temporary:
            work = Path(temporary)
            fasta = work / "reference.fa"
            reads = work / "query.fa"
            fasta.write_text(f">check\n{reference}\n", encoding="utf-8")
            reads.write_text(f">probe\n{query}\n", encoding="utf-8")
            self.command(
                [str(builder), "--packed", "--threads", "2", str(fasta), str(work / "index")],
                work / "build.log",
            )
            log = work / "align.log"
            self.command(
                [str(aligner), "-f", "-v", "3", "-a", "--best", str(work / "index"), str(reads)],
                log,
            )
            actual = set()
            for line in log.read_text(encoding="utf-8").splitlines():
                fields = line.split("\t")
                if len(fields) < 8:
                    continue
                if fields[0] != "probe" or fields[2] != "check" or fields[1] not in ("+", "-"):
                    raise ValueError(
                        "The native tool compatibility test returned unexpected output."
                    )
                position = int(fields[3])
                if not (set(reference[position : position + len(query)]) - set("ACGT")):
                    actual.add((fields[1], position))
            if actual != expected or len(expected) < 5:
                raise ValueError(
                    "Native search failed its both-strand, 0-3 mismatch correctness test. "
                    "The new reference was not activated."
                )

    def _fasta(self, source: Path, destination: Path, scope: str, db: sqlite3.Connection) -> None:
        """Normalize FASTA into fixed-width records and save byte offsets for exact footprints."""
        with gzip.open(source, "rb") as incoming, destination.open("wb") as output:
            accession, title, length, offset = "", "", 0, 0
            pending = b""
            for index, raw in enumerate(incoming):
                if index % 10000 == 0:
                    self._check()
                if raw.startswith(b">"):
                    if accession:
                        if pending:
                            output.write(pending + b"\n")
                        db.execute(
                            "INSERT INTO sequences VALUES(?,?,?,?,?,?,?)",
                            (scope, accession, title, offset, 60, 61, length),
                        )
                    title = raw[1:].decode("utf-8").strip()
                    accession = title.split()[0]
                    output.write(b">" + accession.encode() + b"\n")
                    offset, length, pending = output.tell(), 0, b""
                else:
                    sequence = raw.strip().upper()
                    length += len(sequence)
                    pending += sequence
                    while len(pending) >= 60:
                        output.write(pending[:60] + b"\n")
                        pending = pending[60:]
            if accession:
                if pending:
                    output.write(pending + b"\n")
                db.execute(
                    "INSERT INTO sequences VALUES(?,?,?,?,?,?,?)",
                    (scope, accession, title, offset, 60, 61, length),
                )


def reference_bases(stream: BinaryIO, offset: int, start: int, length: int) -> str:
    """Read a 0-based interval from a normalized 60-base-line reference FASTA."""
    stream.seek(offset + start // 60 * 61 + start % 60)
    raw = stream.read(length + (start % 60 + length) // 60 + 1)
    sequence = raw.replace(b"\n", b"")[:length].decode("ascii")
    if len(sequence) != length or ">" in sequence:
        raise ValueError("A local alignment extends outside its reference sequence.")
    return sequence
