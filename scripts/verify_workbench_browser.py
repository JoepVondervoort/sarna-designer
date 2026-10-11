"""Check live and archived workbench behaviour in an isolated Chromium/SQLite session.

Requires a Chromium/headless-shell executable, supplied with --chrome. No browser package
or network access is needed by default. --gene SCN1A exercises existing cached Ensembl data;
--gene kcnq2 --online verifies cache-miss recovery and a live design with an empty cache.
"""

import argparse
import base64
import json
import os
import re
import socket
import struct
import subprocess
import time
from collections.abc import Callable
from itertools import pairwise
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread
from typing import Any, cast
from urllib.parse import urlsplit
from urllib.request import urlopen

from sarna import settings
from sarna.app.service import DesignService
from sarna.app.web import AppRouter, make_server
from sarna.genome.coordinates import GenomicInterval
from sarna.offtarget.models import GENOMIC_SCOPE
from sarna.offtarget.ncbi import (
    BlastClient,
    SubmittedSearch,
    parameters_for,
    parse_report,
)
from sarna.offtarget.proximity import BlockedTranscripts

DNA = "GCTCTACTCTCAGTTCCAAACTATGACCTACTATGATACCATAC"


class Browser:
    """Minimal local DevTools client for actual DOM interactions, without UI test dependencies."""

    def __init__(self, address: str) -> None:
        """Connect to the Chromium instance launched by this script."""
        parsed = urlsplit(address)
        self.socket = socket.create_connection((parsed.hostname or "127.0.0.1", parsed.port or 80))
        self.socket.settimeout(60)
        key = base64.b64encode(os.urandom(16)).decode()
        self.socket.sendall((f"GET {parsed.path} HTTP/1.1\r\nHost: {parsed.netloc}\r\n"
                             "Upgrade: websocket\r\nConnection: Upgrade\r\n"
                             f"Sec-WebSocket-Key: {key}\r\nSec-WebSocket-Version: 13\r\n\r\n"
                             ).encode())
        header = b""
        while not header.endswith(b"\r\n\r\n"):
            header += self.read(1)
        if b" 101 " not in header:
            raise RuntimeError(f"DevTools upgrade failed: {header!r}")
        self.ident = 0
        self.exceptions: list[dict[str, Any]] = []

    def read(self, length: int) -> bytes:
        """Read an exact number of bytes, refusing premature disconnection."""
        data = b""
        while len(data) < length:
            chunk = self.socket.recv(length - len(data))
            if not chunk:
                raise EOFError("Chromium disconnected")
            data += chunk
        return data

    def call(self, method: str, **params: object) -> dict[str, Any]:
        """Send one masked WebSocket message and await its matching DevTools result."""
        self.ident += 1
        payload = json.dumps({"id": self.ident, "method": method, "params": params}).encode()
        length = len(payload)
        size = length if length < 126 else 126 if length < 65536 else 127
        header = bytes([129, 128 | size])
        if length >= 126:
            header += struct.pack("!H" if length < 65536 else "!Q", length)
        mask = os.urandom(4)
        self.socket.sendall(header + mask + bytes(c ^ mask[i % 4] for i, c in enumerate(payload)))
        while True:
            first, second = self.read(2)
            if first & 15 != 1:
                raise RuntimeError(f"Unexpected DevTools WebSocket frame: {first}")
            size = second & 127
            if size == 126:
                size = struct.unpack("!H", self.read(2))[0]
            elif size == 127:
                size = struct.unpack("!Q", self.read(8))[0]
            response = json.loads(self.read(size))
            if response.get("method") == "Runtime.exceptionThrown":
                self.exceptions.append(response)
            if response.get("id") == self.ident:
                if "error" in response:
                    raise RuntimeError(response["error"])
                return cast("dict[str, Any]", response["result"])

    def evaluate(self, expression: str) -> object:
        """Evaluate a browser expression, retaining JavaScript exceptions as failures."""
        result = self.call("Runtime.evaluate", expression=expression,
                           awaitPromise=True, returnByValue=True)
        if "exceptionDetails" in result:
            raise RuntimeError(result["exceptionDetails"])
        return result["result"].get("value")

    def wait(self, expression: str) -> None:
        """Wait for a DOM condition with a bounded real-time deadline."""
        message = json.dumps("Timed out: " + expression)
        self.evaluate("new Promise((resolve,reject)=>{const end=Date.now()+45000;"
                      f"function poll(){{if({expression})resolve(true);"
                      f"else if(Date.now()>end)reject(new Error({message}));"
                      "else setTimeout(poll,100)}poll()})")

    def screenshot(self, path: Path) -> None:
        """Save an actual viewport capture when a screenshot directory was requested."""
        data = self.call("Page.captureScreenshot", format="png")["data"]
        path.write_bytes(base64.b64decode(data))


def check_expression(browser: Browser) -> None:
    """Check the expression panel beside the picker, and that it changes no ordering.

    The panel is optional evidence from a second service, so this does not require it to load.
    What it does require is that the picker's order is the annotation's either way, and that a
    median of zero never reads as missing data. Both are claims a reader would act on.
    """
    order_before = browser.evaluate(
        "[...document.getElementById('reference').options].map(o=>o.value).join('|')")
    # The panel is filled by a request the lookup deliberately does not wait for.
    browser.wait("!document.getElementById('expression').hidden")
    panel = browser.evaluate(r"""(()=>{
      const sidebar=document.getElementById('expression');
      const overview=document.getElementById('expression-overview');
      const control=document.getElementById('reference');
      const options=[...control.options];
      const payload=state.expression;
      // A select reserves roughly 28px for its arrow and padding. Measured rather than
      // eyeballed, because the overflow is silent: the browser simply clips the end of the
      // label, which is where the expression value is.
      const probe=document.createElement('span');
      probe.style.cssText='position:absolute;visibility:hidden;white-space:pre;font:'
        + getComputedStyle(control).font;
      document.body.appendChild(probe);
      const room=control.clientWidth-28;
      const tooWide=options.filter(o=>{probe.textContent=o.textContent;
        return probe.offsetWidth>room+2;}).map(o=>o.textContent);
      probe.remove();
      return {
        available: payload?.available !== false,
        text: sidebar.textContent + ' ' + overview.textContent,
        sidebarLength: sidebar.textContent.length,
        order: options.map(o=>o.value).join('|'),
        annotated: options.filter(o=>/brain [0-9<]/.test(o.textContent)).length,
        titled: options.filter(o=>o.title).length,
        tooWide: tooWide,
        transcripts: payload?.transcripts ? Object.keys(payload.transcripts).length : 0,
        keyed: payload?.transcripts && !Array.isArray(payload.transcripts),
      };
    })()""")
    assert isinstance(panel, dict), panel
    assert panel["order"] == order_before, (
        "the expression panel changed the order of the promoter reference picker; annotation "
        "decides that order (ADR 0012) and expression is evidence beside it"
    )
    if not panel["available"]:
        # A GTEx outage or a cold cache offline. It must say so without implying absence.
        assert "says nothing about whether these transcripts are expressed" in str(
            panel["text"]), panel["text"]
        print("    expression: not retrieved, and the panel says so")
        return
    assert panel["keyed"] is True, "the payload is a list, which invites sorting by expression"
    assert panel["transcripts"] > 0, panel
    assert panel["annotated"] > 0, "no option carries its expression summary"
    assert panel["titled"] > 0, "no option carries the untruncated label in its title"
    # The sidebar answers for one transcript and must stay short; the comparison of all of
    # them lives in the gene overview card, which has the width for it.
    assert panel["sidebarLength"] < 420, (
        f"the sidebar panel is {panel['sidebarLength']} characters, which is the crowding it "
        f"was split in two to avoid"
    )
    assert len(panel["tooWide"]) <= 1, (
        f"{len(panel['tooWide'])} promoter reference options are wider than the control, so "
        f"the browser silently clips the end of each: {panel['tooWide'][:4]}"
    )
    text = str(panel["text"])
    # The claim, not its wording: the panel must name AF22, because "is my transcript
    # expressed in AF22" is the question that prompted it, and must say what GTEx can and
    # cannot conclude about it.
    assert "AF22" in text, text
    assert "iPSC-derived neurons" in text, text
    assert "supports or contradicts" in text, text
    # The finding this panel exists for: a transcript measured at zero says measured, not
    # missing. Rendering the two alike turns ignorance into a finding and back again.
    zero = browser.evaluate(
        "Object.values(state.expression.transcripts).find(t=>t.label.includes('0.00'))?.label")
    if zero:
        assert "measured at 0.00" in str(zero), zero
        assert "no data" not in str(zero).lower(), zero
    check_organ_strip(browser)
    print(f"    expression: {panel['transcripts']} transcripts, "
          f"{panel['annotated']} options annotated, {len(panel['tooWide'])} too wide, "
          f"order unchanged")


def check_rejection_legend(browser: Browser) -> None:
    """Check the excluded track's legend measures the marks the track actually draws.

    The legend listed ``rejection_breakdown``, which counts every gate a candidate failed, beside
    a track coloured by the *first* gate that excluded each one, and captioned the pair with the
    first-gate rule. On a 600 bp fixture it printed 67.5 % next to a colour covering 6.0 % of the
    marks, and spent a palette slot on a rule accounting for 7 of 574 while a cause of 27 fell
    into the grey "other". Both numbers were individually true, which is why a green suite and a
    reading of the code did not catch it: the defect was that they answered different questions.

    The guard is therefore that every percentage printed in the legend has to be recomputable
    from the candidates the track draws, and that they add up to the whole.
    """
    measured = cast("dict[str, Any]", browser.evaluate("""(()=>{
      const excluded = state.result.candidates.filter(c => c.rank === null);
      const counts = new Map();
      for (const c of excluded) {
        const first = c.excluded_by[0];
        counts.set(first, (counts.get(first) || 0) + 1);
      }
      const shares = [...counts].sort((a, b) => b[1] - a[1])
        .map(([rule, n]) => Math.round(1000 * n / excluded.length) / 10);
      const legend = document.getElementById('rejection-legend');
      const printed = [...legend.querySelectorAll('span')]
        .map(s => s.textContent.trim())
        .map(t => (t.match(/([\\d.]+)%$/) || [])[1])
        .filter(v => v !== undefined).map(Number);
      return {shares, printed, excluded: excluded.length,
              caption: (legend.querySelector('small') || {}).textContent || ""};
    })()"""))
    excluded = cast("int", measured["excluded"])
    printed = cast("list[float]", measured["printed"])
    shares = cast("list[float]", measured["shares"])
    if not excluded:
        assert not printed, "the legend names exclusions on a run that excluded nothing"
        return
    assert printed, f"{excluded} targets were excluded and the legend printed no shares"
    # The named chips must be the largest actual causes, in order, to the tenth of a percent.
    for index, (shown, actual) in enumerate(zip(printed, shares, strict=False)):
        assert abs(shown - actual) <= 0.1, (
            f"legend entry {index} says {shown}% where the drawn marks are {actual}%"
        )
    total = sum(printed)
    assert abs(total - 100.0) <= 0.5, (
        f"the legend's shares add up to {total}%, so they are not shares of one whole"
    )
    assert "counted once" in measured["caption"], "the legend no longer states its attribution"


def check_identical_duplexes(browser: Browser) -> None:
    """Check the table's twin chips agree with the grouping Python computed.

    The grouping travels in the run so the chip, the detail note and the order CSV cannot come
    to different answers about which rows are one molecule. This asserts the browser is reading
    that field rather than deciding for itself, and that a chip appears exactly where the field
    says there is a twin. A promoter without repeats renders no chips, which is correct and is
    checked as such rather than skipped.
    """
    measured = cast("dict[str, Any]", browser.evaluate("""(()=>{
      const groups = state.result.identical_duplex_groups;
      const rows = [...document.querySelectorAll('#candidate-rows tr[data-candidate]')];
      const wrong = rows.filter(row => {
        const twins = (groups || {})[row.dataset.candidate] || [];
        const chip = row.querySelector('.twin-chip');
        if (!twins.length) return chip !== null;
        return chip === null || !chip.textContent.includes(String(twins.length));
      }).map(row => row.dataset.candidate);
      return {present: groups !== undefined, wrong, rows: rows.length};
    })()"""))
    assert measured["present"], "the run carries no identical_duplex_groups for the table to read"
    assert not measured["wrong"], (
        f"twin chips disagree with the grouping for {measured['wrong']}"
    )


def check_organ_strip(browser: Browser) -> None:
    """Check that an organ symbol claims only what GTEx measured.

    The strip is decoration with a factual job: a symbol says "measured above zero in at least
    one of this organ's sample groups". The ways it could lie are a symbol beside a transcript
    measured at zero, a symbol with no value behind it, and a cultured cell group dressed as a
    body site.
    """
    organs = browser.evaluate(r"""(()=>{
      const payload = state.expression;
      const entries = Object.values(payload.transcripts || {});
      const chips = [...document.querySelectorAll('#expression-overview .organ')];
      const silent = entries.filter(e => e.status === 'measured'
        && (e.all || []).every(m => m.median === 0));
      const everyOrganHasAValue = entries.every(e =>
        (e.organs || []).every(o => o.median > 0 && o.name && o.tissue_name));
      return {
        chips: chips.length,
        untitled: chips.filter(c => !c.title).length,
        withOrgans: entries.filter(e => (e.organs || []).length).length,
        silent: silent.length,
        silentWithOrgans: silent.filter(e => (e.organs || []).length).length,
        notInRelease: entries.filter(e => e.status !== 'measured'
          && (e.organs || []).length).length,
        everyOrganHasAValue,
        culturedFlagged: entries.every(e => (e.organs || []).every(
          o => o.organ_id !== 'cultured' || o.is_body_site === false)),
        emojiDetected: emojiAvailable(),
      };})()""")
    assert isinstance(organs, dict), organs
    assert organs["withOrgans"] > 0, "no transcript shows any organ, so the strip is dead"
    assert organs["chips"] > 0, "the organ strip rendered nothing into the table"
    assert organs["untitled"] == 0, (
        f"{organs['untitled']} organ marks carry no hover text, so a reader cannot find out "
        f"what they stand for or what the value behind them is")
    assert organs["silentWithOrgans"] == 0, (
        f"{organs['silentWithOrgans']} transcripts measured at zero everywhere still show an "
        f"organ, which claims expression that was not measured")
    assert organs["notInRelease"] == 0, (
        "a transcript absent from the annotation release shows an organ, which invents a "
        "measurement for something nobody quantified")
    assert organs["everyOrganHasAValue"] is True, (
        "an organ is shown without the median and sample group behind it")
    assert organs["culturedFlagged"] is True, (
        "GTEx's cultured cell groups are presented as a body site")
    print(f"    organ strip: {organs['chips']} marks over {organs['withOrgans']} transcripts, "
          f"{organs['silent']} silent transcripts show none, "
          f"emoji font {'present' if organs['emojiDetected'] else 'absent, names shown'}")


def check_cell_evidence(browser: Browser) -> None:
    """Check that the cell-line box says what choosing a line does, and what it does not.

    The control looked inert because it nearly is: one of the registry's 24 (line x axis) pairs
    is available. That is worth stating rather than leaving a researcher to infer, and the
    statement has to stay true, so this measures the registry rather than trusting a number
    written into the page.
    """
    totals = browser.evaluate(r"""(()=>{
      const cells = state.catalogue.cells || [];
      const axes = cells.flatMap(c => c.data || []);
      return {lines: cells.length, axes: axes.length,
              available: axes.filter(a => a.available).length};
    })()""")
    assert isinstance(totals, dict), totals
    assert totals["available"] >= 1, (
        "no cell line has any available evidence axis, so the control carries nothing at all"
    )

    def panel(value: str) -> str:
        browser.evaluate(
            "(()=>{const box=document.getElementById('cell-line');"
            f"box.value={json.dumps(value)};"
            "box.dispatchEvent(new Event('input'));})()")
        return str(browser.evaluate("document.getElementById('cell-evidence').innerText"))

    # Whatever is typed, the panel must say the choice does not touch ranking. It does not:
    # nothing downstream of this field reaches a gate, a preference or a rank.
    for typed in ("", "AF22 (lt-NES)", "NotARealLine-7"):
        assert "does not change which candidates pass" in panel(typed), typed

    known = panel("AF22 (lt-NES)")
    assert "CVCL_A5CW" in known, "the accession is missing, so the line cannot be checked"
    assert "AVAILABLE" in known.upper(), known
    assert "NO ROUTE" in known.upper(), (
        "every axis reads as available, which would hide the gaps this registry exists to record"
    )
    assert "GSE108506" in known, "the one available axis does not name its source"

    unknown = panel("NotARealLine-7")
    assert "not in the registry" in unknown, unknown
    assert "No measurement is implied" in unknown, unknown
    assert "AVAILABLE" not in unknown.upper(), (
        "an unknown name shows an availability claim, which would invent evidence for it"
    )
    browser.evaluate(
        "(()=>{const box=document.getElementById('cell-line');box.value='';"
        "box.dispatchEvent(new Event('input'));})()")
    print(f"    cell lines: {totals['available']} of {totals['axes']} evidence axes available "
          f"across {totals['lines']} lines, and the panel names which")


class StubBlast:
    """A BLAST service that answers from a script, so the panel can be checked without NCBI.

    Shaped to produce the cases the panel has to tell apart: one candidate with a full-length
    exact match somewhere else, one that occurs in many places, one with only weak partial
    alignments, and the rest with nothing at all. The genomic script also returns each
    candidate's own locus, which is what the panel must remove rather than count.

    The XML goes through the real parser, so the per-candidate budget and the capped flag are
    exercised here too rather than mocked away.
    """

    def __init__(self) -> None:
        """Start with no submission and no coordinates to place hits against."""
        self.loci: dict[str, tuple[str, int, int]] = {}
        self.submitted: dict[str, str] = {}
        self.scope = ""

    def submit(self, sequences: dict[str, str], scope: str) -> Any:  # noqa: ANN401
        """Record the job exactly as the real client would send it."""
        self.submitted = dict(sequences)
        self.scope = scope
        return SubmittedSearch(request_id="STUB-RID-1", estimated_seconds=0,
                               parameters=parameters_for(scope))

    def fetch(self, request_id: str, parameters: Any) -> Any:  # noqa: ANN401
        """Return a scripted report for whatever was submitted."""
        return parse_report(self._xml(),
                            per_candidate=parameters.max_alignments_per_candidate)

    def close(self) -> None:
        """Nothing to close; the factory hands out one shared stub."""

    def _hit(self, accession: str, title: str, identities: int,
             start: int, frame: int = 1) -> str:
        # Real reports carry the aligned strings, and the panel reads the mismatched positions
        # out of them. Scattering the mismatches rather than parking them at the 3' end keeps
        # the footprint honest: every alignment drawn identically is the failure mode here.
        query = "ACGTACGTACGTACGTACG"
        spread = [3, 11, 7, 16, 1, 14][: 19 - identities]
        midline = "".join(
            " " if index + 1 in spread else "|" for index in range(19))
        subject = "".join(
            "T" if index + 1 in spread else base for index, base in enumerate(query))
        return (f"<Hit><Hit_accession>{accession}</Hit_accession><Hit_def>{title}</Hit_def>"
                f"<Hit_hsps><Hsp><Hsp_identity>{identities}</Hsp_identity>"
                f"<Hsp_align-len>19</Hsp_align-len><Hsp_gaps>0</Hsp_gaps>"
                f"<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>19</Hsp_query-to>"
                f"<Hsp_hit-from>{start}</Hsp_hit-from><Hsp_hit-to>{start + 18}</Hsp_hit-to>"
                f"<Hsp_hit-frame>{frame}</Hsp_hit-frame>"
                f"<Hsp_qseq>{query}</Hsp_qseq><Hsp_hseq>{subject}</Hsp_hseq>"
                f"<Hsp_midline>{midline}</Hsp_midline>"
                f"<Hsp_evalue>{0.004 if identities > 16 else 7.0}</Hsp_evalue>"
                f"<Hsp_bit-score>{identities * 2}</Hsp_bit-score></Hsp></Hit_hsps></Hit>")

    def _xml(self) -> bytes:
        genomic = self.scope == GENOMIC_SCOPE
        iterations = []
        for index, name in enumerate(self.submitted):
            hits = []
            if genomic and name in self.loci:
                contig, start, _end = self.loci[name]
                number = {"X": 23, "Y": 24}.get(
                    contig, int(contig) if contig.isdigit() else 1)
                hits.append(self._hit(f"NC_{int(number):06d}.12",
                                      "Homo sapiens chromosome, own target locus",
                                      19, start + 1))
            # Genomic reports come back on chromosome accessions. Reusing mRNA ones would
            # leave every hit unplaceable, and the distance plot untested.
            def named(slot: int, symbol: str) -> tuple[str, str]:
                return ((f"NC_0000{(slot % 22) + 1:02d}.11",
                         f"Homo sapiens chromosome {(slot % 22) + 1}")
                        if genomic else
                        (f"NM_{100000 + slot}.1",
                         f"Homo sapiens {symbol.lower()} ({symbol}), mRNA"))

            if index == 0:
                accession, title = named(6, "GAPDH")
                hits.append(self._hit(accession, title, 19, 12_455))
                hits += [self._hit(*named(10 + n, f"WEAK{n}"), 13, 40_000 + 100 * n)
                         for n in range(3)]
            elif index == 1:
                hits += [self._hit(*named(n, f"PP{n}"),
                                   17 if n % 2 else 16, 90_000 + 4000 * n) for n in range(7)]
            elif index == 2:
                hits += [self._hit(*named(14, "POX"), 14, 300_900)]
            if not hits:
                continue
            iterations.append(
                f"<Iteration><Iteration_query-def>{name}</Iteration_query-def>"
                f"<Iteration_hits>{''.join(hits)}</Iteration_hits></Iteration>")
        return ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_iterations>"
                + "".join(iterations)
                + "</BlastOutput_iterations></BlastOutput>").encode()


def check_offtarget(browser: Browser, blast: StubBlast, screenshots: Path | None) -> None:
    """Check the off-target panel: batches, the summary, where a hit lands, and the limits.

    The table this replaces printed six counts per candidate and nothing else, so a run of 22
    showed 21 rows of zeros and no indication of where the one real hit was. What matters here
    is that the alignments reach the screen, that a candidate screened clean is distinguishable
    from one never screened, and that the own-locus removal is stated rather than silent.
    """
    candidates = cast(list[dict[str, Any]], browser.evaluate("state.result.candidates") or [])
    blast.loci = {
        c["id"]: (c["genomic"]["contig"], c["genomic"]["start"], c["genomic"]["end"])
        for c in candidates
        if c.get("genomic")
    }
    pool = browser.evaluate(
        "state.result.candidates.filter(c=>c.rank!==null).length")
    assert isinstance(pool, int) and pool > 0

    # A batch of 5, not the whole run. The default used to be a fixed 50, which is what made a
    # genomic search time out after half an hour and come back with nothing.
    browser.evaluate("document.querySelector('[data-batch=\"5\"]').click();"
                     "document.getElementById('screen-transcript').click()")
    browser.wait("document.querySelector('#offtarget-results .ot-table')!==null")
    screened = browser.evaluate("Object.keys(state.screenings.transcript.byCandidate).length")
    assert screened == min(5, pool), f"asked for 5, screened {screened}"

    summary = str(browser.evaluate("document.querySelector('.ot-summary').innerText"))
    assert "Wang 2015's 80%" in summary, summary
    assert "not measured binding" in summary, summary
    if pool >= 2:
        assert "occurring in many places" in summary, summary
    if screened < pool:
        note = str(browser.evaluate("document.getElementById('offtarget-results').innerText"))
        assert "not the same as having been searched clean" in note, (
            "the unscreened candidates are not distinguished from the clean ones"
        )

    # Where it binds, which never reached the screen before: open the first row and read the
    # gene, the matched bases and the footprint out of the DOM rather than out of the record.
    browser.evaluate("document.querySelector('[data-screenrow]').click()")
    browser.wait("document.querySelector('.ot-detail')!==null")
    detail = browser.evaluate(r"""(()=>{
      const row=document.querySelector('.ot-detail');
      return {text:row.innerText,
              footprints:row.querySelectorAll('svg.footprint').length,
              shaded:[...row.querySelectorAll('svg.footprint rect')]
                .filter(r=>r.getAttribute('fill')!=='#e8eeeb').length,
              mismatched:[...row.querySelectorAll('svg.footprint rect')]
                .filter(r=>r.getAttribute('fill')==='#c08a2e').length,
              shapes:new Set([...row.querySelectorAll('svg.footprint')]
                .map(s=>[...s.querySelectorAll('rect')]
                  .map(r=>r.getAttribute('fill')).join(''))).size,
              seed:row.innerText.toLowerCase().includes('seed region')};
    })()""")
    assert isinstance(detail, dict), detail
    assert detail["footprints"] > 0, "no alignment footprint was drawn"
    assert detail["shaded"] > 0, "a footprint with nothing shaded shows no alignment at all"
    # A partial match has to look different from a full one. The first version drew only the
    # aligned span, so 16 of 19 and 19 of 19 were the same picture and the footprint said
    # nothing a reader could not already read in the identity column.
    assert detail["mismatched"] > 0, (
        "no footprint marks a mismatched base, so a partial match is drawn as a perfect one"
    )
    assert detail["shapes"] > 1, (
        f"every footprint in this detail is identical ({detail['shapes']} distinct), so the "
        f"drawing carries no information the identity column did not"
    )
    assert not detail["seed"], (
        "a seed region is marked; Wang 2015 names none and ADR 0023 refuses to import one"
    )
    # Named by gene rather than by accession: an accession alone makes a reader look it up,
    # which is exactly the work the panel exists to save.
    assert re.search(r"\b(PP\d|GAPDH|POX)\b", str(detail["text"])), (
        f"no hit is named by its gene symbol: {str(detail['text'])[:200]}"
    )
    assert "GAPDH" in str(browser.evaluate(
        "document.querySelector('#offtarget-results .ot-table tbody').innerText")), (
        "the strongest match column does not name the gene of the full-length exact match"
    )
    # The promiscuous candidate's list is collapsed rather than drawn in full, and says so.
    collapsed = str(browser.evaluate("document.querySelector('.ot-detail').innerText"))
    if "many places" in collapsed:
        assert "strongest of" in collapsed, (
            "a candidate occurring in many places drew its whole list, which is the wall of "
            "rows this change exists to remove"
        )

    strip = browser.evaluate(
        "document.querySelectorAll('#offtarget-results .strip i').length")
    assert isinstance(strip, int) and strip > 0, "no identity strip was drawn"

    # The next batch is a second dated record merged into the same table, not a replacement.
    if pool > 5:
        browser.evaluate("document.querySelector('[data-batch-mode=\"rest\"]')?.click()")
        browser.wait("Object.keys(state.screenings.transcript.byCandidate).length>5")
        merged = browser.evaluate("state.screenings.transcript.searches.length")
        assert isinstance(merged, int) and merged >= 2, (
            "the second batch replaced the first instead of merging with it"
        )

    # The genomic scope removes the candidate's own locus, or says it could not.
    browser.evaluate("document.getElementById('screen-genomic').click()")
    browser.wait("state.screenings.genomic!==undefined")
    browser.evaluate("document.querySelector('[data-scope-tab=\"genomic\"]')?.click()")
    genomic = str(browser.evaluate("document.getElementById('offtarget-results').innerText"))
    if blast.loci:
        assert "own target loci" in genomic, (
            "the on-target hit was not removed, so every candidate reads as a perfect match"
        )
        perfect = browser.evaluate(
            "Object.values(state.screenings.genomic.byCandidate)"
            ".every(e=>e.row.perfect_hits===0||e.row.on_target_hits===0)")
        assert perfect is True, "an own-locus hit is still being counted as a perfect match"
    else:
        assert "no genomic coordinates" in genomic, (
            "a run without coordinates must say the removal could not be made"
        )
    # Three other places state what has been screened, and all three used to hardcode the
    # absence. A page that has just run a screening and still says none was performed is
    # worse than one that says nothing: it is the duplex summary a reader acts on.
    browser.evaluate("document.querySelector('[data-scope-tab=\"transcript\"]')?.click()")
    for element, where in (("top-result", "the duplex summary"),
                           ("candidate-detail", "the candidate panel")):
        text = str(browser.evaluate(
            f"document.getElementById({json.dumps(element)}).innerText"))
        assert "No off-target screening has been performed" not in text, (
            f"{where} still reports that no screening was performed after one ran"
        )
    check_progress_bar(browser)
    check_human_only(browser)
    if screenshots:
        browser.evaluate("document.querySelector('[data-scope-tab=\"genomic\"]')?.click();"
                         "document.getElementById('offtarget-card')"
                         ".scrollIntoView({block:'start'})")
        time.sleep(0.4)
        browser.screenshot(screenshots / "offtarget-genomic.png")
        browser.evaluate("document.querySelector('.ot-plot')"
                         "?.scrollIntoView({block:'center'})")
        time.sleep(0.4)
        browser.screenshot(screenshots / "offtarget-tss-plot.png")
        browser.evaluate("document.querySelector('[data-scope-tab=\"transcript\"]')?.click();"
                         "document.getElementById('offtarget-card')"
                         ".scrollIntoView({block:'start'})")
        time.sleep(0.4)
        browser.screenshot(screenshots / "offtarget-panel.png")
    print(f"    off-target: {screened} of {pool} screened in the first batch, alignments named "
          f"and placed, own-locus handling stated")


def check_progress_bar(browser: Browser) -> None:
    """The bar must not claim a completion it cannot know.

    NCBI reports a search as started or finished and nothing between; the only forecast is
    RTOE, its own estimate at submission. So the bar is exercised directly with crafted
    states, because what matters is the judgement in it rather than the animation: a real
    fraction where one exists, and an honest refusal where one does not.
    """
    states = browser.evaluate(r"""(()=>{
      const read = job => {
        const html = progressBar(job);
        const host = document.createElement('div');
        host.innerHTML = html;
        const fill = host.querySelector('.bar i');
        return {
          html,
          indeterminate: !!host.querySelector('.bar.indeterminate'),
          width: fill && fill.style.width ? parseFloat(fill.style.width) : null,
          text: host.textContent,
        };
      };
      const running = extra => ({status:'running', progress:extra});
      return {
        early: read(running({phase:'waiting', elapsed_seconds:10, estimate_seconds:100})),
        nearly: read(running({phase:'waiting', elapsed_seconds:99, estimate_seconds:100})),
        over: read(running({phase:'waiting', elapsed_seconds:500, estimate_seconds:100})),
        noEstimate: read(running({phase:'waiting', elapsed_seconds:40, estimate_seconds:0})),
        placing: read(running({phase:'placing', placed:5, to_place:10})),
        idle: read({status:'complete'}),
      };
    })()""")
    assert isinstance(states, dict), states
    early, nearly, over = states["early"], states["nearly"], states["over"]
    assert early["width"] == 10, early
    assert "estimated by NCBI" in early["text"], early["text"]
    # Never past nine tenths on an estimate. A bar that reaches the end and then keeps going
    # reads as a stall, and this one can legitimately run for half an hour.
    assert nearly["width"] is not None and nearly["width"] <= 90, nearly
    assert over["indeterminate"], "an overrun still shows a filling bar"
    assert "longer than" in over["text"], over["text"]
    assert states["noEstimate"]["indeterminate"], states["noEstimate"]
    assert "no percentage to show" in states["noEstimate"]["text"], states["noEstimate"]["text"]
    # Placing is the one phase with a real fraction, so it is the one that may claim it.
    assert states["placing"]["width"] == 50, states["placing"]
    assert states["placing"]["indeterminate"] is False
    assert states["idle"]["html"] == "", "a finished job still draws a progress bar"
    # A search still running at NCBI must offer to be collected, not read as a failure. The
    # id outlives this page, so the only wrong answer is to hide it.
    waiting = browser.evaluate(r"""(()=>{
      const kept = state.screeningJobs;
      state.screeningJobs = {transcript:{status:'waiting', id:'RID-TEST',
        stage:'Still running at NCBI.'}};
      renderScreeningStatus();
      const lane = document.getElementById('screen-bar-transcript');
      const answer = {
        button: !!lane.querySelector('[data-resume="transcript"]'),
        text: lane.textContent,
        tag: document.getElementById('offtarget-status').textContent,
        buttonEnabled: !document.getElementById('screen-transcript').disabled,
      };
      state.screeningJobs = kept;
      renderScreeningStatus();
      return answer;
    })()""")
    assert isinstance(waiting, dict), waiting
    assert waiting["button"], "a search still at NCBI offers no way to collect it"
    assert "Still running at NCBI" in str(waiting["text"]), waiting["text"]
    assert waiting["tag"] == "AT NCBI", (
        f"a search still running at NCBI is labelled {waiting['tag']!r}, which reads as "
        f"finished or failed when it is neither"
    )
    print("    progress bar: estimate capped at 90%, overrun says so, placing is a real "
          "fraction, a search left at NCBI can be collected")


def check_human_only(browser: Browser) -> None:
    """Both scopes are human, in two different ways, and must not be worded alike.

    The genomic assembly holds one species. RefSeq RNA holds every species NCBI has and is
    narrowed by ENTREZ_QUERY, which the URL API does not document. Calling both "human" with
    no further detail would overstate the weaker of the two.
    """
    card = str(browser.evaluate(
        "document.getElementById('offtarget-card').innerText"))
    assert "human sequence only" in card, card[:200]
    for scope, expected in (("transcript", "Entrez query"), ("genomic", "one species")):
        browser.evaluate(
            f"document.querySelector('[data-scope-tab=\"{scope}\"]')?.click()")
        text = str(browser.evaluate(
            "document.getElementById('offtarget-results').innerText"))
        assert "HUMAN ONLY" in text, f"{scope}: no human marker in the results"
        assert expected in text, (
            f"{scope}: the summary does not say how the search was confined to human, so the "
            f"two scopes read as the same guarantee. Got: {text[:300]}"
        )
    database = browser.evaluate(
        "state.screenings.genomic.searches[0].parameters.database")
    assert database == "GPIPE/9606/current/ref_top_level", database
    assert not browser.evaluate(
        "state.screenings.genomic.searches[0].parameters.organism"), (
        "the human assembly was sent with a species filter it does not need"
    )
    # A full-width progress bar dropped into the flex row squeezed the status text to one
    # character per line, which every assertion above still passed. Measure the box.
    squeezed = browser.evaluate(r"""(()=>{
      return ['transcript','genomic'].map(scope => {
        const node = document.getElementById('screen-state-' + scope);
        const box = node.getBoundingClientRect();
        return {scope, text:node.textContent.trim(), width:box.width, height:box.height};
      }).filter(s => s.text && (s.width < 60 || s.height > 40));
    })()""")
    assert squeezed == [], f"the screening status text is squeezed out of shape: {squeezed}"
    print("    human only: stated on the card and in both scopes, each naming how")


def exercise(browser: Browser, url: str, directory: Path, gene: str | None,
             screenshots: Path | None, online: bool = False,
             blast: StubBlast | None = None) -> None:
    """Check draft preservation, rerun, explanations, persistence and an offline archive."""
    browser.call("Runtime.enable")
    browser.call("Emulation.setDeviceMetricsOverride", width=1440, height=1100,
                 deviceScaleFactor=1, mobile=False)
    browser.call("Page.navigate", url=url)
    browser.wait("document.getElementById('preset')?.options.length===4")
    assert browser.evaluate("document.getElementById('gene-symbol').value") == ""
    forced = browser.evaluate("state.catalogue.data_access.offline_forced")
    assert browser.evaluate("document.getElementById('offline').disabled") is forced
    check_cell_evidence(browser)
    if gene:
        browser.evaluate(f"document.getElementById('gene-symbol').value={json.dumps(gene)};"
                         "document.getElementById('offline').checked=true;"
                         "document.getElementById('offline').dispatchEvent(new Event('change'));"
                         "document.getElementById('lookup-gene').click()")
        if online:
            browser.wait("!document.getElementById('error-recovery').hidden")
            assert "Cache-only mode" in str(browser.evaluate(
                "document.getElementById('global-error').textContent"))
            assert "Cache only" in str(browser.evaluate(
                "document.getElementById('data-access-status').textContent"))
            browser.evaluate("document.getElementById('retry-online').click()")
        browser.wait("document.getElementById('reference').value!==''")
        assert browser.evaluate("document.getElementById('gene-symbol').value") == gene.upper()
        assert browser.evaluate("document.getElementById('offline').checked") is not online
        assert browser.evaluate("document.getElementById('error-recovery').hidden") is True
        check_expression(browser)
    else:
        browser.evaluate("document.querySelector('[data-source=sequence]').click();"
                         f"document.getElementById('dna-sequence').value={json.dumps(DNA)}")
    browser.evaluate("document.getElementById('motifs').checked=false;"
                     "document.getElementById('regulatory').checked=false;"
                     "document.getElementById('design-form').requestSubmit()")
    browser.wait("!document.getElementById('result-area').hidden")
    browser.wait("document.querySelector('#candidate-detail .sequence-box')!==null")
    before = browser.evaluate("document.getElementById('metric-passing').textContent")
    original = str(browser.evaluate("state.result.app.run_id"))
    browser.evaluate("document.getElementById('open-settings').click()")
    rule = "tss-window" if gene else "pos19-a"
    browser.evaluate(f"const role=document.querySelector('[data-rule=\"{rule}\"] [data-role]');"
                     "role.value='preference';"
                     "role.dispatchEvent(new Event('change',{bubbles:true}));"
                     "document.querySelector('#settings-dialog .close-dialog').click()")
    assert browser.evaluate("!document.getElementById('draft-notice').hidden") is True
    assert browser.evaluate("document.getElementById('metric-passing').textContent") == before
    browser.evaluate("document.getElementById('design-form').requestSubmit()")
    browser.wait("!document.getElementById('diff-notice').hidden")
    run_id = str(browser.evaluate("state.result.app.run_id"))
    assert run_id != original
    assert browser.evaluate("document.getElementById('global-error').textContent") == ""
    browser.evaluate("document.querySelector('[data-filter=all]').click();"
                     "document.getElementById('select-top').click();"
                     "document.getElementById('zoom-in').click();"
                     "document.getElementById('pan-right').click()")
    browser.evaluate("document.querySelector('#promoter-map [data-tip]').dispatchEvent("
                     "new FocusEvent('focusin',{bubbles:true}))")
    assert browser.evaluate("!document.getElementById('tooltip').hidden") is True
    # Every candidate the map draws must be inside the clip region. A clip ending above the
    # excluded-target track once hid 1,643 of 2,468 rejected candidates completely and left the
    # rest a 2px sliver, while the legend still invited the reader to hover them. The promise
    # that a rejection is never silently dropped is only kept if it is actually on screen.
    clipped = browser.evaluate("""(()=>{
      const clip=document.querySelector('#map-clip rect');
      const top=+clip.getAttribute('y'), bottom=top + +clip.getAttribute('height');
      const marks='#promoter-map rect[data-candidate],#promoter-map rect[data-zoom-from]';
      return [...document.querySelectorAll(marks)]
        .filter(r => +r.getAttribute('y') < top
                  || +r.getAttribute('y') + +r.getAttribute('height') > bottom).length;
    })()""")
    assert clipped == 0, f"{clipped} candidate marks fall outside the promoter map's clip region"
    drawn = browser.evaluate(
        "document.querySelectorAll('#promoter-map rect[data-candidate]').length")
    assert isinstance(drawn, int) and drawn > 0
    check_rejection_legend(browser)
    check_identical_duplexes(browser)
    if blast is not None:
        check_offtarget(browser, blast, screenshots)
    # The ten highest-ranked targets have to be findable among the passing ones. Before this
    # they were drawn in the same colour at the same height as rank 900, so the map could not
    # answer which marks a reader would actually order. Checked on all three channels the
    # drawing uses, because any one alone fails somebody. Checked at full extent, where every
    # candidate is in view: the preceding zoom and pan leave some of the ten off screen, and
    # a mark outside the viewport is correctly not drawn rather than missing.
    browser.evaluate("document.getElementById('zoom-reset').click()")
    top = browser.evaluate(r"""(()=>{
      const fill = r => (r.getAttribute('fill')||'').toLowerCase();
      const y = r => +r.getAttribute('y');
      const marks = [...document.querySelectorAll('#promoter-map rect[data-candidate]')];
      const byId = new Map(state.result.candidates.map(c => [c.id, c]));
      const rankOf = r => byId.get(r.getAttribute('data-candidate'))?.rank ?? null;
      // The passing lane by position, not by colour: every rank now has its own fill, and
      // the selection outline and the excluded row also carry data-candidate.
      const lane = marks.filter(r => fill(r) !== 'none' && y(r) >= 225 && y(r) <= 255);
      const marked = lane.filter(r => +r.getAttribute('height') > 15);
      const labels = [...document.querySelectorAll('#promoter-map text.rank-label')]
        .map(t => +t.textContent);
      // Relative luminance, so "gets paler as rank gets worse" is measured rather than assumed.
      const luminance = hex => {
        const v = i => {
          const c = parseInt(hex.slice(1 + i * 2, 3 + i * 2), 16) / 255;
          return c <= 0.03928 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
        };
        return 0.2126 * v(0) + 0.7152 * v(1) + 0.0722 * v(2);
      };
      const ordered = [...marked].sort((a, b) => rankOf(a) - rankOf(b));
      return {
        expected: state.result.candidates.filter(c => c.rank !== null && c.rank <= TOP_MARKED)
          .length,
        marked: marked.length,
        ranks: marked.map(rankOf),
        distinctFills: new Set(marked.map(fill)).size,
        // Each mark must carry the colour its own rank maps to, not merely some colour.
        mismatched: marked.filter(r => fill(r) !== rankColour(rankOf(r)).fill.toLowerCase())
          .map(rankOf),
        luminances: ordered.map(r => luminance(fill(r))),
        heights: [...new Set(marked.map(r => +r.getAttribute('height')))],
        plainHeights: [...new Set(lane.filter(r => !marked.includes(r))
          .map(r => +r.getAttribute('height')))],
        labels, labelled: new Set(labels).size,
      };})()""")
    assert isinstance(top, dict)
    if top["expected"]:
        assert top["marked"] == top["expected"], (
            f"{top['expected']} targets rank in the top ten but {top['marked']} are drawn as "
            f"such; the rest are indistinguishable from rank 900")
        assert all(isinstance(rank, int) and 1 <= rank <= 10 for rank in top["ranks"]), (
            f"a mark drawn in the top-ten lane does not hold a top-ten rank: {top['ranks']}")
        assert not top["mismatched"], (
            f"these ranks are drawn in another rank's colour: {top['mismatched']}")
        assert top["distinctFills"] == top["marked"], (
            f"{top['marked']} marks share {top['distinctFills']} colours, so the ramp does "
            f"not separate every rank")
        # Rank order must read off the shading. Strictly paler each step, and deliberately
        # one hue: a green-to-red scale would claim rank 1 is good, which the 125 assayed
        # sequences refute (ADR 0006, ADR 0021).
        luminances = top["luminances"]
        assert all(a < b for a, b in pairwise(luminances)), (
            f"the rank shading is not monotone, so a darker mark does not mean a better "
            f"rank: {luminances}")
        # Height carries the same claim as colour, for a reader who cannot use colour.
        assert top["heights"] and max(top["heights"]) > max(top["plainHeights"] or [0]), (
            f"top-ten marks are {top['heights']} high against {top['plainHeights']} for the "
            f"rest; height has to separate them too")
        assert 1 in top["labels"], "rank 1 is not numbered on the map"
        assert top["labelled"] == len(top["labels"]), (
            f"a rank is printed more than once on the map: {top['labels']}")
        assert all(1 <= rank <= 10 for rank in top["labels"]), top["labels"]
        print(f"    rank ramp: {top['marked']} marks, {top['distinctFills']} distinct shades, "
              f"monotone, {len(top['labels'])} numbered")
    # At full-region zoom a target is generated at every offset, so drawing each excluded one
    # produced about three overlapping marks per pixel: a band that claimed to show every
    # rejection and resolved none of them. They are binned instead. The promise that a
    # rejection is never silently dropped now means the bins must account for all of them.
    browser.evaluate("document.getElementById('zoom-reset').click()")
    binned = browser.evaluate(r"""(()=>{
      const bins=[...document.querySelectorAll('#promoter-map rect[data-zoom-from]')];
      const counted=bins.reduce((sum,b)=>{
        const m=(b.getAttribute('data-tip')||'').match(/(\d+) excluded target/);
        return sum + (m?+m[1]:0);}, 0);
      return {bins:bins.length, counted,
              excluded: state.result.candidates.filter(c=>c.rank===null).length,
              span: state.result.promoter.length,
              marks: document.querySelectorAll('#promoter-map rect[data-candidate]').length,
              legend: document.querySelector('#rejection-legend small').textContent,
              binnedMode: state.excludedBinned};})()""")
    assert isinstance(binned, dict)
    # A short pasted sequence is narrower than the separability threshold, so drawing each
    # target individually is the correct behaviour there, not a failure to bin.
    wide = binned["span"] > 280
    if binned["excluded"] and not wide:
        assert binned["binnedMode"] is False, (
            f"a {binned['span']} bp region is narrow enough to draw every target")
        assert binned["marks"] > 0, "no individual targets were drawn in a narrow region"
        assert "individually" in binned["legend"], binned["legend"]
    if binned["excluded"] and wide:
        assert binned["binnedMode"] is True, "the full region should bin excluded targets"
        assert binned["bins"] > 0, "no density columns were drawn"
        assert binned["counted"] == binned["excluded"], (
            f"the bins account for {binned['counted']} excluded targets but the run has "
            f"{binned['excluded']}; some rejection is not represented anywhere")
        # Which mode the track is in is a claim about what the reader is seeing, so it is said.
        assert "same height" in binned["legend"], binned["legend"]
        # Every rejection stays reachable: a column zooms to the targets it counts.
        reached = browser.evaluate("""(()=>{
          const bin=document.querySelector('#promoter-map rect[data-zoom-from]');
          bin.dispatchEvent(new MouseEvent('click',{bubbles:true}));
          return {marks:document.querySelectorAll('#promoter-map rect[data-candidate]').length,
                  binned:state.excludedBinned,
                  legend:document.querySelector('#rejection-legend small').textContent};})()""")
        assert isinstance(reached, dict)
        assert reached["binned"] is False, "clicking a density column did not zoom in far enough"
        assert reached["marks"] > 0, "zooming into a column drew no individual targets"
        assert "individually" in reached["legend"], reached["legend"]
    browser.evaluate("document.getElementById('zoom-reset').click()")

    # The legend's colours must be the map's colours. They were three separate lists and had
    # drifted: the legend showed #e2b772 for the regulatory build against a map drawing
    # #ddb370 and #efcf99, and #a798d2 for motifs against a map drawing #6a5ca8. A legend
    # whose swatches are not the drawing's colours is worse than none, because it is
    # confidently wrong. Both now come from one table, and this checks that they still do.
    legend = browser.evaluate(r"""(()=>{
      const hex = value => {
        const m = value.match(/rgba?\((\d+), ?(\d+), ?(\d+)/);
        return m ? "#" + [1,2,3].map(i=>(+m[i]).toString(16).padStart(2,"0")).join("") : value;
      };
      // The authoritative table, plus the rejection colours, which are generated per run
      // rather than declared. A track with no features in this run still belongs in the
      // legend, so what is checked is the shared source and not what happens to be drawn.
      // The rank ramp is computed from TOP_MARKED rather than listed, so its ten shades come
      // from the same function the map draws with. Listing them here by hand is how a legend
      // and a map drift apart, which is what this whole check exists to prevent.
      const ramp = Array.from({length: TOP_MARKED},
        (_, i) => rankColour(i + 1).fill.toLowerCase());
      const declared = new Set([...Object.values(TRACK_COLOURS).map(c => c.fill.toLowerCase()),
                                ...ramp]);
      const generated = new Set([...REJECTION_COLOURS, OTHER_REJECTION]
        .map(c => c.toLowerCase()));
      const drawn = new Set([...document.querySelectorAll('#promoter-map rect[fill]')]
        .map(r => (r.getAttribute('fill')||'').toLowerCase())
        .filter(f => f && f !== 'none' && f !== 'transparent'));
      const legend = document.getElementById('map-legend');
      legend.open = true;
      const swatches = [...legend.querySelectorAll('.legend-mark i')]
        .map(i => hex(getComputedStyle(i).backgroundColor).toLowerCase());
      return {drawn:[...drawn], swatches, declared:[...declared], generated:[...generated],
              undeclared:[...drawn].filter(f => !declared.has(f) && !generated.has(f)),
              rows:[...legend.querySelectorAll('.legend-table th')].map(t=>t.textContent.trim()),
              collapsible: legend.tagName === 'DETAILS',
              chips:[...legend.querySelectorAll('.legend-chip')].length};})()""")
    assert isinstance(legend, dict), legend
    assert legend["collapsible"] is True, "the legend cannot be collapsed"
    assert legend["chips"] == 7, legend["chips"]
    # Every row of the map is named in the legend, in the order the map draws them.
    assert legend["rows"] == ["Reference / TSS", "Design windows", "CpG islands",
                              "Regulatory build", "TF motif matches", "Passing targets",
                              "Excluded targets"], legend["rows"]
    stray = [s for s in legend["swatches"] if s not in legend["declared"]]
    assert not stray, (
        f"the legend shows {stray}, which is not in the map's colour table "
        f"{sorted(legend['declared'])}; the two have drifted apart again")
    assert not legend["undeclared"], (
        f"the map draws {legend['undeclared']}, which is in neither the shared colour table "
        f"nor the generated rejection palette, so no legend entry can describe it")
    # A promoter must not be the same colour as any other regulatory feature.
    assert "#e8b5a6" in legend["swatches"] and "#efcf99" in legend["swatches"], legend["swatches"]

    # Each factor gets its own colour, and the neutral for factors past the end of the
    # palette must not be one of the palette's own: while it was, SOX2 fell off the end and
    # rendered identically to KLF4 with nothing on screen to tell them apart.
    assert browser.evaluate(
        "MOTIF_COLOURS.includes(TRACK_COLOURS.motif.fill)") is False, (
        "the neutral motif colour is also a per-factor colour, so two factors can collide")
    assert browser.evaluate("new Set(MOTIF_COLOURS).size === MOTIF_COLOURS.length") is True

    # Panning translates the clipped group instead of redrawing 2,600 SVG nodes, because a full
    # redraw measured 14.5 ms against a 16.7 ms frame budget. That shortcut is only legitimate
    # if the translated picture is the one the redraw would have produced, so this compares
    # them: it reads a mark's position, drags, predicts from the transform, releases, and
    # checks the redrawn position against the prediction.
    # Zoomed in, because at full extent there is nowhere to pan and the map says so by not
    # offering the grab cursor at all.
    browser.evaluate("document.getElementById('zoom-in').click()")
    state_span = int(cast("int", browser.evaluate("state.result.promoter.length")))
    drag = browser.evaluate(r"""(()=>{
      const svg=document.getElementById('promoter-map');
      // Track a mark that stays clear of both edges. renderMap clamps a partly-scrolled mark
      // to the viewport and floors its width at 2 units so thin features stay visible, while
      // the drag translates the true rectangle and lets the clip path cut it. Under the clip
      // those agree everywhere except a sub-2-unit sliver at the very edge, so comparing a
      // clamped mark compares the floor rather than the panning arithmetic.
      const inside=r=>{const x=+r.getAttribute('x'), w=+r.getAttribute('width');
                       return x>140 && x+w<958;};
      const candidates=[...svg.querySelectorAll('[data-candidate]')].filter(inside);
      if(!candidates.length) return {error:"no unclamped candidate mark to track"};
      const id=candidates[Math.floor(candidates.length/2)].dataset.candidate;
      const xOf=()=>{const r=svg.querySelector('[data-candidate="'+id+'"]');
                     return r && inside(r) ? +r.getAttribute('x') : null;};
      const box=()=>svg.getBoundingClientRect();
      const send=(target,type,x)=>target.dispatchEvent(new PointerEvent(type,{bubbles:true,
        cancelable:true,clientX:x,clientY:box().top+box().height/2,button:0,buttons:1,
        pointerId:1,isPrimary:true}));
      const b=box(), from=b.left+b.width*0.6, to=b.left+b.width*0.3;
      const before=xOf(), viewportBefore=state.start;
      const idle=getComputedStyle(svg).cursor;
      send(svg,'pointerdown',from);
      const grabbing=getComputedStyle(svg).cursor;
      send(window,'pointermove',to);
      const node=svg.querySelector('#map-pan');
      const tx=node?parseFloat((node.getAttribute('transform')||'translate(0').match(
        /translate\(([-\d.]+)/)[1]):null;
      send(window,'pointerup',to);
      return {before, tx, after:xOf(), grabbing, idle, panned:state.start!==viewportBefore,
              cleared:svg.querySelector('#map-pan').getAttribute('transform')===null,
              inRange: state.start>=0 && state.end<=state.result.promoter.length+1e-6};
    })()""")
    assert isinstance(drag, dict) and "error" not in drag, drag
    assert drag["grabbing"] == "grabbing", "the map gives no grab affordance while dragging"
    assert drag["idle"] == "grab", "a pannable map should offer the grab cursor at rest"
    assert drag["panned"] is True, "dragging the promoter map did not pan it"
    assert drag["cleared"] is True, "the drag transform outlived the drag"
    assert drag["inRange"] is True, "panning left the promoter"
    if drag["before"] is not None and drag["after"] is not None:
        predicted = drag["before"] + drag["tx"]
        assert abs(predicted - drag["after"]) < 0.01, (
            f"the drag drew the mark at {predicted} but the redraw puts it at {drag['after']}; "
            "the fast pan path and the full redraw disagree")
    else:
        # Only acceptable on a region too small to pan a mark clear of the edges.
        assert state_span <= 280, (
            "the tracked mark left the plot area on a region wide enough to keep it inside, "
            "so the pan arithmetic was never actually compared")
    browser.evaluate("document.getElementById('zoom-reset').click()")

    # The measured ranking performance is shown, not left to the footer.
    assert browser.evaluate("!document.getElementById('ranking-performance').hidden") is True
    assert "does not beat chance" in str(browser.evaluate(
        "document.getElementById('ranking-performance').textContent"))
    # A rule row must look openable. `display:flex` on a <summary> drops Chrome's native
    # marker, and for a while these seventeen rows looked like a static checklist while the
    # explanation -- the measured value, the citation, the rationale -- sat one invisible
    # click away. The affordance is now drawn, and this asserts the content behind it.
    opened = str(browser.evaluate("""(()=>{
      // Scoped to the rule list: the guide-prediction disclosure shares the .verdict class.
      const rows=[...document.querySelectorAll('#candidate-detail .rule-verdicts .verdict')];
      if (!rows.length) return "NO RULE ROWS";
      const marker = getComputedStyle(rows[0].querySelector('summary'), '::before').content;
      rows[0].open = true;
      return marker + "|" + rows[0].innerText;})()"""))
    assert not opened.startswith("NO RULE ROWS")
    assert not opened.startswith("none|"), "the rule rows have no disclosure affordance"
    assert "why" in opened.lower()
    # The explanation itself, not just the verdict.
    assert "Why this rule:" in opened and "Changing it:" in opened, opened[:300]
    # The cell line's identity and evidence. All of this reached the browser from the start and
    # none of it was drawn -- the banner printed the name the researcher typed and stopped. A
    # cell line named without its accession cannot be checked against anybody else's record.
    cells = browser.evaluate("""(()=>{
      const banner=document.getElementById('context-banner');
      const unregistered=banner.innerText;
      // Now the registered path, using a panel the server actually sent.
      const known=(state.result.cell_types||[]).find(c=>c.rrid);
      if(!known) return {unregistered, error:"no cell panel carries an accession"};
      renderCellContext(state.result, known);
      banner.querySelectorAll('details').forEach(d=>d.open=true);
      const link=banner.querySelector('.accession');
      return {unregistered, name:known.name, rrid:known.rrid,
              href:link?link.getAttribute('href'):null,
              linkText:link?link.textContent.trim():null,
              axes:[...banner.querySelectorAll('.axis-table tbody th')]
                     .map(t=>t.childNodes[0].textContent.trim()),
              states:[...banner.querySelectorAll('.axis-state')].map(s=>s.textContent.trim()),
              text:banner.innerText};})()""")
    assert isinstance(cells, dict) and "error" not in cells, cells
    # An unrecognised name must say so rather than imply the ranking used a cell context.
    assert "No record is held for this name" in cells["unregistered"], cells["unregistered"][:200]
    # The registered path: accession, resolvable, and every axis named rather than left raw.
    assert cells["rrid"].startswith("CVCL_"), cells["rrid"]
    assert cells["href"] == f"https://www.cellosaurus.org/{cells['rrid']}", cells["href"]
    assert cells["rrid"] in cells["linkText"], cells["linkText"]
    assert cells["axes"] == ["Expression", "Chromatin", "Factor occupancy"], cells["axes"]
    assert set(cells["states"]) <= {"Available", "Not wired in", "No route"}, cells["states"]
    # Saying "Loaded" would credit this run with a measurement it may not have made.
    assert "This describes the cell line, not this run" in cells["text"]

    # The results open on the answer: how many survived and what to order.
    assert browser.evaluate("!document.getElementById('top-result').hidden") is True
    top = str(browser.evaluate("document.getElementById('top-result').textContent"))
    assert "passed every evaluable gate" in top, top[:160]
    assert "[dT][dT]" in top, "the top card must show the order-ready duplex, not the bare core"
    # The duplex is drawn as it anneals, not as three independent 5'->3' strings. The rungs
    # are computed from the bases rather than assumed, so a strand that is not the reverse
    # complement of its partner would draw crosses here rather than quietly look correct.
    ladder = browser.evaluate("""(()=>{
      const rows=[...document.querySelectorAll('#candidate-detail .ladder-row')];
      if (rows.length!==3) return {error:"expected sense, rungs and antisense rows"};
      // Only the core columns: .ov is an overhang and .pad is the gap opposite one, so
      // including them misaligns the two strands by two columns and compares a dT with a blank.
      const core=r=>[...r.querySelectorAll('.bases span.b')].map(s=>s.textContent.trim());
      const rungs=r=>[...r.querySelectorAll('.bases span.rung')].map(s=>s.textContent.trim());
      const roles=rows.filter(r=>r.querySelector('.ladder-label i'))
                      .map(r=>r.querySelector('.ladder-label i').textContent);
      return {sense:core(rows[0]), rungs:rungs(rows[1]), anti:core(rows[2]), roles,
              overhangs:document.querySelectorAll('#candidate-detail .ladder .ov').length};
    })()""")
    assert isinstance(ladder, dict) and "error" not in ladder, ladder
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}
    sense, anti = ladder["sense"], ladder["anti"]
    assert len(sense) == len(anti) == len(ladder["rungs"]), ladder
    # Both directions, because one of them is not enough. Checking only that *drawn* rungs
    # pair is vacuous: a wrong pairing table simply draws fewer rungs, and every rung that
    # survives is still correct. Planting RNA_PAIR = {G:"G", C:"C"} passed that test while the
    # drawing silently marked every G/C column as a mismatch. So each column is required to
    # draw a rung exactly when its two bases genuinely pair.
    blank = ("", "\u00a0")
    for index, (top, bottom) in enumerate(zip(sense, anti, strict=True)):
        drew_rung = ladder["rungs"][index] not in blank and ladder["rungs"][index] != "\u00d7"
        should_pair = pairs.get(top) == bottom
        assert drew_rung == should_pair, (
            f"column {index} pairs {top} with {bottom} "
            f"({'a real pair' if should_pair else 'not a pair'}) but "
            f"{'drew' if drew_rung else 'did not draw'} a rung")
    # And the engine builds the antisense strand as the reverse complement, so on a correct
    # record every column pairs. A cross anywhere means the data itself is wrong.
    assert all(pairs.get(t) == b for t, b in zip(sense, anti, strict=True)), (
        f"the antisense strand is not the reverse complement of the sense strand: "
        f"{''.join(sense)} against {''.join(anti)}")
    # Both 3' dTdT overhangs are drawn: two bases on each of the two strands.
    assert ladder["overhangs"] == 4, ladder["overhangs"]
    assert "guide" in ladder["roles"] and "passenger" in ladder["roles"], ladder["roles"]
    # The subtraction behind the guide call is shown as numbers, not only inside a sentence.
    ends = browser.evaluate(r"""(()=>{
      const boxes=[...document.querySelectorAll('#candidate-detail .end-box b')]
        .map(b=>parseFloat(b.textContent));
      const diff=parseFloat(document.querySelector('#candidate-detail .end-verdict b').textContent);
      const stacks=[...document.querySelectorAll('#candidate-detail .end-box small')]
        .map(s=>s.textContent.split('+').map(t=>parseFloat(t.trim().split(/\s+/)[1])));
      return {boxes, diff, stacks};})()""")
    assert isinstance(ends, dict)
    five, three = ends["boxes"]
    assert abs((three - five) - ends["diff"]) < 0.011, (
        f"the panel shows E(3')={three} and E(5')={five} but a difference of {ends['diff']}")
    for energy, steps in zip(ends["boxes"], ends["stacks"], strict=True):
        assert abs(sum(steps) - energy) < 0.011, (
            f"the stacks {steps} do not sum to the {energy} printed above them")
    # The printable sheet builds from the record, and keeps gates and preferences apart: a
    # passed gate listed as a preference would contradict the rating printed beside it.
    sheet = str(browser.evaluate(
        "(()=>{const r=state.result, best=r.candidates.find(c=>c.rank===1);"
        "document.getElementById('duplex-sheet').innerHTML=duplexSheetHTML(r,best);"
        "return document.getElementById('duplex-sheet').innerText;})()"))
    assert "Order these two oligonucleotides" in sheet
    # Which tube is which, on the rows that name the tubes.
    assert "GUIDE" in sheet and "PASSENGER" in sheet, (
        "the sheet names two strands without saying which is loaded")
    assert "Anneal 1:1" in sheet
    assert sheet.count("[dT][dT]") >= 2, "both strands belong on the sheet"
    assert "Gates passed" in sheet and "Preferences met" in sheet
    # The sheet goes to the bench, so it has to carry the screening that was run rather than
    # the absence it used to hardcode, and it has to carry the date and database with it.
    assert "No off-target screening has been performed" not in sheet, (
        "the bench sheet still reports that no screening was performed after one ran"
    )
    assert "refseq" in sheet, f"the bench sheet omits the database screened against: {sheet[:300]}"
    assert "Wang 2015's 80%" in sheet, "the bench sheet omits the sourced homology count"
    assert "does not beat chance" in sheet, "the sheet must carry the measurement too"
    stated = re.search(r"meeting (\d+) of (\d+) evaluated preferences", sheet)
    listed = re.search(r"Preferences met \((\d+)\)", sheet)
    assert stated and listed, sheet[:400]
    assert stated.group(1) == listed.group(1), (
        f"the sheet says {stated.group(1)} preferences met but lists {listed.group(1)}")
    # The live order export goes to Python over HTTP. Fetch it the way the button does and read
    # what a researcher would actually paste into a vendor form.
    live_order = str(browser.evaluate(
        f"fetch('/api/runs/{run_id}/export.order.csv').then(r=>r.text())"
    ))
    rows = [r for r in live_order.splitlines() if r.strip()]
    assert "sense_order_5to3" in rows[0], rows[0]
    assert len(rows) > 1, "the live order export returned no duplexes"
    assert all("[dT][dT]" in row for row in rows[1:]), rows[1]
    # A screening has now been run against this route, so the order rows must carry its
    # findings with the database and the date rather than still reading "not-performed". The
    # row is what gets pasted into a vendor form, so the caveats have to travel in it.
    screening_cell = rows[1]
    assert "off-target-screening [not-performed]" not in screening_cell, (
        "a screening was run, yet the order row still reports that none was performed"
    )
    assert "refseq" in screening_cell, (
        f"the order row does not name the database it was screened against: {screening_cell}"
    )
    assert "Wang" in screening_cell or "80%" in screening_cell, (
        f"the order row omits the sourced homology count: {screening_cell}"
    )
    # The Excel export, fetched the way the button does and checked as a real package. An
    # .xlsx exists because Excel reads a CSV's separator from the machine's locale, so a
    # colleague on a Dutch Windows opens a comma-separated file with everything in column A.
    book = browser.evaluate(
        f"fetch('/api/runs/{run_id}/export.xlsx').then(async r =>"
        "({status:r.status, type:r.headers.get('content-type'),"
        " size:(await r.arrayBuffer()).byteLength}))"
    )
    assert isinstance(book, dict), book
    assert book["status"] == 200, book
    assert "spreadsheetml.sheet" in str(book["type"]), book
    assert isinstance(book["size"], int) and book["size"] > 2000, book
    assert browser.evaluate(
        "!!document.getElementById('export-xlsx')") is True, "no Excel button"

    if screenshots:
        browser.screenshot(screenshots / "workbench.png")
    browser.evaluate(f"const history=document.getElementById('saved-runs');"
                     f"history.value={json.dumps(original)};"
                     "history.dispatchEvent(new Event('change'))")
    browser.wait(f"state.result.app.run_id==={json.dumps(original)}")
    assert browser.evaluate("document.getElementById('metric-passing').textContent") == before
    if gene:
        assert browser.evaluate("document.getElementById('offline').checked") is not online
        expected = "Online retrieval allowed" if online else "Cache only"
        assert expected in str(browser.evaluate(
            "document.getElementById('data-access-status').textContent"))
    if online:
        browser.evaluate("document.getElementById('offline').checked=true;"
                         "document.getElementById('offline').dispatchEvent(new Event('change'));"
                         "document.getElementById('upstream').value=2001;"
                         "document.getElementById('upstream').dispatchEvent(new Event('input',"
                         "{bubbles:true}));document.getElementById('design-form').requestSubmit()")
        browser.wait("!state.busy && !document.getElementById('error-recovery').hidden")
        assert "sequence/region" in str(browser.evaluate(
            "document.getElementById('global-error').textContent"))
        assert browser.evaluate("state.result.app.run_id") == original
        browser.evaluate("document.getElementById('retry-online').click()")
        browser.wait(f"!state.busy && state.result.app.run_id!=={json.dumps(original)}")
        assert browser.evaluate("state.result.promoter.length") == 2501
        assert browser.evaluate("document.getElementById('global-error').textContent") == ""
    with urlopen(f"{url}api/runs/{run_id}/export.html") as response:
        archive = directory / "archive.html"
        archive.write_bytes(response.read())
    browser.call("Page.navigate", url=archive.as_uri())
    browser.wait("document.body?.classList.contains('frozen') && "
                 "document.querySelector('#candidate-detail .sequence-box')!==null")
    browser.evaluate("document.querySelector('[data-filter=rejected]').click()")
    browser.wait("document.querySelector('#candidate-rows [data-candidate]')!==null")
    browser.evaluate("document.querySelector('#candidate-rows [data-candidate]').click()")
    browser.wait("document.querySelector('#candidate-detail .verdict .outcome.fail')!==null")
    assert browser.evaluate("state.result.candidates.every(c=>Array.isArray(c.rules))") is True
    # A passing candidate in the archive, which is the case the rejected one above cannot
    # cover: the duplex panel is built from result.order, and the archive has no backend to
    # ask for it. If that table were not embedded, the panel would quietly degrade to "this
    # candidate failed a gate" on a candidate that passed every one.
    archived = browser.evaluate(r"""(()=>{
      const best=state.result.candidates.find(c=>c.rank===1);
      if(!best) return {error:"the archive holds no passing candidate"};
      // The table is still filtered to rejected candidates from the check above, so the
      // passing one is not in the DOM to click until the filter is switched back.
      document.querySelector('[data-filter=passing]').click();
      const row=[...document.querySelectorAll('#candidate-rows [data-candidate]')]
        .find(r=>r.dataset.candidate===best.id);
      if(!row) return {error:"rank 1 is not in the passing table"};
      row.click();
      return new Promise(resolve=>setTimeout(()=>{
        const panel=document.querySelector('#candidate-detail .duplex-panel');
        resolve({
          present: !!panel,
          rungs: document.querySelectorAll('#candidate-detail .ladder .rung').length,
          overhangs: document.querySelectorAll('#candidate-detail .ladder .ov').length,
          tubes: [...document.querySelectorAll('#candidate-detail .role-chip')]
                   .map(c=>c.textContent.trim()),
          ends: document.querySelectorAll('#candidate-detail .end-box b').length,
          text: panel?panel.innerText.slice(0,200):"",
        });}, 400));})()""")
    assert isinstance(archived, dict) and "error" not in archived, archived
    assert archived["present"], "the archive lost the duplex panel"
    assert archived["rungs"] == 19, archived["rungs"]
    assert archived["overhangs"] == 4, archived["overhangs"]
    assert "guide" in archived["tubes"] and "passenger" in archived["tubes"], archived["tubes"]
    assert archived["ends"] == 2, "the end energies did not survive the export"
    assert "failed a gate" not in archived["text"], archived["text"]

    browser.evaluate("document.getElementById('run-evidence').click()")
    assert browser.evaluate("document.getElementById('evidence-dialog').open") is True
    browser.evaluate("document.querySelector('#evidence-dialog .close-dialog').click()")
    # The archive has no backend, so its order CSV must come from rows Python stored in the
    # record. Building the text is the whole of what the button does offline.
    order_csv = str(browser.evaluate("orderCsvText(state.result)"))
    header, first = order_csv.splitlines()[0], order_csv.splitlines()[1]
    assert "sense_order_5to3" in header, header
    assert "[dT][dT]" in first, first
    assert "off-target-screening [not-performed]" in first
    assert browser.evaluate(
        "state.result.order.rows.every(r=>r[4].endsWith('[dT][dT]')"
        " && r[5].endsWith('[dT][dT]') && r[9]==='21')"
    ) is True
    browser.call("Emulation.setDeviceMetricsOverride", width=390, height=844,
                 deviceScaleFactor=1, mobile=True)
    assert browser.evaluate("document.documentElement.scrollWidth<=innerWidth") is True
    if screenshots:
        browser.screenshot(screenshots / "archive-mobile.png")
    assert not browser.exceptions, browser.exceptions
    print(f"PASS: {gene or 'pasted DNA'}; draft retained {before} passing; "
          "Python rerun/diff, history, hover, offline rejected verdicts and mobile layout")
    if online:
        print("PASS: empty-cache lookup refused without fetching; visible online retry "
              "retrieved gene, transcripts and promoter; lowercase symbol resolved; "
              "uncached-window retry also completed without losing previous results")


def main() -> None:
    """Launch isolated app/browser processes and always clean up both."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrome", type=Path, required=True)
    parser.add_argument("--gene", help="Optional human gene with cached MANE/promoter responses")
    parser.add_argument("--online", action="store_true",
                        help="Use a new empty cache and test explicit online recovery for --gene")
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    if args.online and (not args.gene or settings.is_offline()):
        parser.error("--online requires --gene and SARNA_OFFLINE unset/off")
    if args.screenshots:
        args.screenshots.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="sarna-browser-") as temporary:
        directory = Path(temporary)
        if args.online:
            os.environ[settings.CACHE_DIR_VARIABLE] = str(directory / "empty-cache")
        blast = StubBlast()
        # Screening is checked against a scripted service, never NCBI: a verification run
        # must not queue a public job, and the panel's job is to render what it is given.
        service = DesignService(
            directory / "runs",
            blast_factory=cast("Callable[[], BlastClient]", lambda: blast),
        )
        service.poll_sleep = lambda _seconds: None
        # Placement is the one part of a genomic screening that would otherwise reach Ensembl.
        # Replaced with two annotated transcripts at known coordinates so the distance plot is
        # drawn from real arithmetic rather than skipped for want of a network.
        def nearby(interval: GenomicInterval) -> list[dict[str, Any]]:
            """Annotated transcripts on a regular grid across whatever is asked for.

            A grid rather than two fixed transcripts because the service snaps lookups to
            megabase blocks and filters the answer back to the hit's own interval: anything
            pinned to the block's start is filtered out again, and every hit comes back
            unplaced. A grid puts a start site within reach of any position, which is what the
            distance plot needs to be exercised at all.
            """
            step, found = 25_000, []
            first = (max(0, interval.start) // step) * step
            for index, start in enumerate(range(first, interval.end + step, step)):
                if start < 1:
                    continue
                found.append({
                    "id": f"ENST{index:011d}", "external_name": f"NEARBY{index % 7}",
                    "start": start, "end": start + 8_000,
                    "strand": -1 if index % 2 else 1,
                })
            return found

        # Wrapped the way the service wraps the real one, so the block cache is exercised too.
        service._transcript_lookup = (  # type: ignore[method-assign]
            lambda: BlockedTranscripts(nearby)
        )
        router = AppRouter(service)
        server = make_server(router, 0)
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        process = subprocess.Popen([
            str(args.chrome), "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-dev-shm-usage", "--disable-background-networking",
            "--remote-debugging-port=0", f"--user-data-dir={directory / 'profile'}", "about:blank",
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        browser = None
        try:
            active_port = directory / "profile" / "DevToolsActivePort"
            deadline = time.monotonic() + 15
            while not active_port.exists():
                if process.poll() is not None or time.monotonic() > deadline:
                    raise RuntimeError("Chromium did not start its local DevTools endpoint")
                time.sleep(0.1)
            port = int(active_port.read_text(encoding="utf-8").splitlines()[0])
            with urlopen(f"http://127.0.0.1:{port}/json/list") as response:
                pages = json.load(response)
            browser = Browser(pages[0]["webSocketDebuggerUrl"])
            exercise(browser, f"http://127.0.0.1:{server.server_port}/", directory,
                     args.gene, args.screenshots, args.online, blast)
        finally:
            if browser:
                browser.socket.close()
            process.terminate()
            process.wait(timeout=10)
            server.shutdown()
            server.server_close()
            router.close()
            thread.join(timeout=5)


if __name__ == "__main__":
    main()
