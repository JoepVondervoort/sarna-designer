"""Check that an installation of saRNA Designer actually works on this machine.

Written for someone with no computational background: it prints one line per check, finishes
with a single PASS or FAIL, and every failure says what to do about it. Run it after installing
and paste the output to whoever gave you the tool.

    sarna-designer-check              # the normal way
    sarna-designer-check --online     # also check that Ensembl can be reached

It deliberately makes **no network calls by default**, so it works on a locked-down machine and
cannot be slow or flaky for the wrong reason. Pass ``--online`` to also check that Ensembl is
reachable, which is needed to look up a gene for the first time.
"""

import argparse
import platform
import sqlite3
import sys
import tempfile
from pathlib import Path

GREEN, RED, YELLOW, RESET = "\033[32m", "\033[31m", "\033[33m", "\033[0m"


class Report:
    """Collects check results and prints them as they happen."""

    def __init__(self, colour: bool) -> None:
        """Start an empty report, with colour unless the terminal cannot show it."""
        self.failures: list[str] = []
        self.warnings: list[str] = []
        self.colour = colour

    def _paint(self, text: str, code: str) -> str:
        """Return ``text`` in colour, or unchanged where colour is unavailable."""
        return f"{code}{text}{RESET}" if self.colour else text

    def ok(self, label: str, detail: str = "") -> None:
        """Record a passing check."""
        print(f"  {self._paint('OK  ', GREEN)} {label}{f' — {detail}' if detail else ''}")

    def fail(self, label: str, fix: str) -> None:
        """Record a failure, with the action that resolves it."""
        print(f"  {self._paint('FAIL', RED)} {label}")
        print(f"       what to do: {fix}")
        self.failures.append(label)

    def warn(self, label: str, detail: str) -> None:
        """Record something worth knowing that does not block use."""
        print(f"  {self._paint('NOTE', YELLOW)} {label} — {detail}")
        self.warnings.append(label)


def check_python(report: Report) -> None:
    """Check the interpreter is the supported version."""
    version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    if sys.version_info[:2] == (3, 11):
        report.ok("Python version", f"{version} on {platform.system()}")
    else:
        report.fail(
            f"Python version is {version}, but this tool needs 3.11",
            "Install with uv, which fetches the right Python for you: see the installation "
            "guide. Do not change your system Python.",
        )


def check_imports(report: Report) -> None:
    """Check the package and its required dependency are importable."""
    try:
        import sarna

        report.ok("saRNA Designer package", f"version {sarna.__version__}")
    except ImportError as error:
        report.fail(
            f"The sarna package is not installed ({error})",
            "Run the install command from the guide again and watch for errors.",
        )
        return
    try:
        import httpx  # noqa: F401

        report.ok("Network library (httpx)")
    except ImportError:
        report.fail(
            "httpx is missing, so gene lookup cannot work",
            "Reinstall with the app extra. The [app] on the end is what installs it: "
            'uv tool install --python 3.11 "sarna-designer[app] @ <where you got it from>"',
        )


def check_engine(report: Report) -> None:
    """Design against a short built-in sequence, with no network and no files."""
    try:
        from sarna.design.order import order_for
        from sarna.design.pipeline import design_sarnas
        from sarna.design.presets import WANG_2015
        from sarna.genome.promoter import Promoter
        from sarna.sequence.nucleic import NucleicSequence

        bases = (
            "GCTTAACGGATTCAAGTTACCGATTACAGGCTAACGTTAGCCATTAAGCGTTACA"
            "GGGTAATTTTGACAGGTCA"
            "ATCGGTTACAGCTTAAGGCTATTACGGCATTAACGGTTACAGCTTAAGGCATTAA"
        )
        result = design_sarnas(
            Promoter(sequence=NucleicSequence.parse(bases)), WANG_2015.resolve()
        )
        generated = len(result.generated)
        if generated == len(bases) - 19 + 1:
            report.ok("Design engine",
                      f"{generated} candidates from a {len(bases)} bp test sequence")
        else:
            report.fail(
                f"The design engine produced {generated} candidates, which is the wrong number",
                "Report this output — it means the installation is damaged, not misconfigured.",
            )
            return
        order = order_for(result.generated[0])
        if order.sense.annotated.endswith("[dT][dT]") and order.sense.length_nt == 21:
            report.ok("Order chemistry", f"21-mers with dTdT, e.g. {order.sense.annotated}")
        else:
            report.fail(
                f"Order chemistry is wrong: {order.sense.annotated}",
                "Report this output.",
            )
    except Exception as error:
        report.fail(
            f"The design engine could not run ({type(error).__name__}: {error})",
            "Report this whole output to whoever gave you the tool.",
        )


def check_storage(report: Report) -> None:
    """Check SQLite can create a database with the settings the run store uses."""
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "check.sqlite3"
            connection = sqlite3.connect(path, timeout=30)
            mode = connection.execute("PRAGMA journal_mode=WAL").fetchone()[0]
            connection.execute("CREATE TABLE t (a TEXT)")
            connection.execute("INSERT INTO t VALUES (?)", ("saved",))
            connection.commit()
            connection.close()
        if str(mode).lower() == "wal":
            report.ok("Local run storage (SQLite)")
        else:
            report.warn(
                "Local run storage",
                f"journal mode is {mode}, not WAL. Saved runs still work; this can happen on a "
                "network drive. Prefer a folder on your own disk.",
            )
    except Exception as error:
        report.fail(
            f"Could not create a local database ({type(error).__name__}: {error})",
            "Choose a data folder on your own disk, not a network or cloud-synced drive.",
        )


def check_text_encoding(report: Report) -> None:
    """Check non-ASCII text survives a write and read, which cp1252 would break."""
    sample = "sodium channel α-subunit — naïve"  # noqa: RUF001 - the point is non-ASCII text
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "encoding.txt"
            path.write_text(sample, encoding="utf-8")
            if path.read_text(encoding="utf-8") == sample:
                report.ok("Text handling", "non-ASCII gene descriptions round-trip correctly")
            else:
                report.fail("Text came back changed after writing", "Report this output.")
    except Exception as error:
        report.fail(f"Could not write a text file ({error})", "Check folder permissions.")


def check_server(report: Report) -> None:
    """Check the local workbench can bind a loopback port."""
    try:
        from sarna.app.service import DesignService
        from sarna.app.web import AppRouter, make_server

        with tempfile.TemporaryDirectory() as tmp:
            router = AppRouter(DesignService(Path(tmp) / "runs"))
            try:
                server = make_server(router, 0)
                port = server.server_port
                server.server_close()
                report.ok("Local web interface", f"can open a browser port (tested on {port})")
            finally:
                router.close()
    except ImportError as error:
        # Reported separately from a port failure, because the remedies share nothing. Sending
        # somebody to look at their firewall when a dependency is missing wastes their
        # afternoon, and this is the likeliest single mistake in the install: leaving the
        # [app] off the end of the command.
        report.fail(
            f"The local web interface is incomplete ({error})",
            "This is a missing dependency, not a network or firewall problem. Reinstall with "
            "the app extra: the [app] on the end of the install command is what adds it.",
        )
    except OSError as error:
        report.fail(
            f"The local web interface could not open a port ({type(error).__name__}: {error})",
            "A firewall or security tool may be blocking local connections. The tool only ever "
            "listens on your own machine (127.0.0.1); allow it if prompted.",
        )
    except Exception as error:
        report.fail(
            f"The local web interface could not start ({type(error).__name__}: {error})",
            "Copy this whole output and send it on; it is written to be useful to whoever "
            "maintains the tool.",
        )


def check_network(report: Report) -> None:
    """Check Ensembl is reachable, which a first-time gene lookup needs."""
    try:
        from sarna.genome.ensembl import EnsemblClient

        client = EnsemblClient()
        gene = client.gene("TP53")
        report.ok("Ensembl connection", f"looked up TP53 ({gene.ensembl_id})")
    except Exception as error:
        report.warn(
            "Ensembl connection",
            f"could not reach Ensembl ({type(error).__name__}). The tool still works on "
            "sequences you paste in, and on genes already cached. Check your internet "
            "connection or institutional proxy.",
        )


def main() -> int:
    """Run every check and return a shell exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", action="store_true",
                        help="also check that Ensembl can be reached")
    parser.add_argument("--no-colour", action="store_true", help="plain text output")
    args = parser.parse_args()

    report = Report(colour=not args.no_colour and sys.stdout.isatty())
    print("\nsaRNA Designer — installation check")
    print(f"  {platform.system()} {platform.release()}, Python {platform.python_version()}\n")
    check_python(report)
    check_imports(report)
    check_engine(report)
    check_storage(report)
    check_text_encoding(report)
    check_server(report)
    if args.online:
        check_network(report)
    else:
        print("  ---- skipped the internet check; add --online to include it")

    print()
    if report.failures:
        print(f"RESULT: FAIL — {len(report.failures)} check(s) did not pass.")
        print("Copy everything above and send it to whoever gave you the tool.")
        return 1
    extra = f" ({len(report.warnings)} note(s) above)" if report.warnings else ""
    print(f"RESULT: PASS — saRNA Designer works on this machine{extra}.")
    print("Start it with:  sarna-designer")
    return 0


if __name__ == "__main__":
    sys.exit(main())
