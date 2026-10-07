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
      return [...document.querySelectorAll('#promoter-map rect[data-candidate]')]
        .filter(r => +r.getAttribute('y') < top
                  || +r.getAttribute('y') + +r.getAttribute('height') > bottom).length;
    })()""")
    assert clipped == 0, f"{clipped} candidate marks fall outside the promoter map's clip region"
    drawn = browser.evaluate(
        "document.querySelectorAll('#promoter-map rect[data-candidate]').length")
    assert isinstance(drawn, int) and drawn > 0
    # Panning translates the clipped group instead of redrawing 2,600 SVG nodes, because a full
    # redraw measured 14.5 ms against a 16.7 ms frame budget. That shortcut is only legitimate
    # if the translated picture is the one the redraw would have produced, so this compares
    # them: it reads a mark's position, drags, predicts from the transform, releases, and
    # checks the redrawn position against the prediction.
    drag = browser.evaluate(r"""(()=>{
      const svg=document.getElementById('promoter-map');
      const id=(state.result.candidates.find(c=>c.rank!==null)||{}).id;
      if(!id) return {error:"no passing candidate to track"};
      const xOf=()=>{const r=svg.querySelector('[data-candidate="'+id+'"]');
                     return r?+r.getAttribute('x'):null;};
      const box=()=>svg.getBoundingClientRect();
      const send=(target,type,x)=>target.dispatchEvent(new PointerEvent(type,{bubbles:true,
        cancelable:true,clientX:x,clientY:box().top+box().height/2,button:0,buttons:1,
        pointerId:1,isPrimary:true}));
      const b=box(), from=b.left+b.width*0.6, to=b.left+b.width*0.3;
      const before=xOf(), viewportBefore=state.start;
      send(svg,'pointerdown',from);
      const grabbing=getComputedStyle(svg).cursor;
      send(window,'pointermove',to);
      const node=svg.querySelector('#map-pan');
      const tx=node?parseFloat((node.getAttribute('transform')||'translate(0').match(
        /translate\(([-\d.]+)/)[1]):null;
      send(window,'pointerup',to);
      return {before, tx, after:xOf(), grabbing, panned:state.start!==viewportBefore,
              cleared:svg.querySelector('#map-pan').getAttribute('transform')===null,
              inRange: state.start>=0 && state.end<=state.result.promoter.length+1e-6};
    })()""")
    assert isinstance(drag, dict) and "error" not in drag, drag
    assert drag["grabbing"] == "grabbing", "the map gives no grab affordance while dragging"
    assert drag["panned"] is True, "dragging the promoter map did not pan it"
    assert drag["cleared"] is True, "the drag transform outlived the drag"
    assert drag["inRange"] is True, "panning left the promoter"
    if drag["before"] is not None and drag["after"] is not None:
        predicted = drag["before"] + drag["tx"]
        assert abs(predicted - drag["after"]) < 0.01, (
            f"the drag drew the mark at {predicted} but the redraw puts it at {drag['after']}; "
            "the fast pan path and the full redraw disagree")
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
      const text=r=>[...r.querySelectorAll('.bases span')].map(s=>s.textContent.trim());
      const roles=rows.filter(r=>r.querySelector('.ladder-label i'))
                      .map(r=>r.querySelector('.ladder-label i').textContent);
      return {sense:text(rows[0]), rungs:text(rows[1]), anti:text(rows[2]), roles,
              overhangs:document.querySelectorAll('#candidate-detail .ladder .ov').length};
    })()""")
    assert isinstance(ladder, dict) and "error" not in ladder, ladder
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}
    sense, anti = ladder["sense"], ladder["anti"]
    assert len(sense) == len(anti) == len(ladder["rungs"]), ladder
    # Every drawn rung must correspond to a real Watson-Crick pair between the two rows.
    for index, rung in enumerate(ladder["rungs"]):
        if rung not in ("", "\u00a0"):
            top, bottom = sense[index], anti[index]
            assert pairs.get(top) == bottom, (
                f"column {index} draws a rung between {top} and {bottom}, which do not pair")
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
