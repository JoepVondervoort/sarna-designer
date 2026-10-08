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
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread
from typing import Any, cast
from urllib.parse import urlsplit
from urllib.request import urlopen

from sarna import settings
from sarna.app.service import DesignService
from sarna.app.web import AppRouter, make_server

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


def exercise(browser: Browser, url: str, directory: Path, gene: str | None,
             screenshots: Path | None, online: bool = False) -> None:
    """Check draft preservation, rerun, explanations, persistence and an offline archive."""
    browser.call("Runtime.enable")
    browser.call("Emulation.setDeviceMetricsOverride", width=1440, height=1100,
                 deviceScaleFactor=1, mobile=False)
    browser.call("Page.navigate", url=url)
    browser.wait("document.getElementById('preset')?.options.length===4")
    assert browser.evaluate("document.getElementById('gene-symbol').value") == ""
    forced = browser.evaluate("state.catalogue.data_access.offline_forced")
    assert browser.evaluate("document.getElementById('offline').disabled") is forced
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
    assert sheet.count("[dT][dT]") >= 2, "both strands belong on the sheet"
    assert "Gates passed" in sheet and "Preferences met" in sheet
    assert "No off-target screening has been performed" in sheet
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
    assert "off-target-screening [not-performed]" in rows[1]
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
        router = AppRouter(DesignService(directory / "runs"))
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
                     args.gene, args.screenshots, args.online)
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
