"""Check live and archived workbench behaviour in an isolated Chromium/SQLite session.

Requires a Chromium/headless-shell executable, supplied with --chrome. No browser package
or network access is needed by default. --gene SCN1A exercises existing cached Ensembl data;
--gene kcnq2 --online verifies cache-miss recovery and a live design with an empty cache.
"""

import argparse
import base64
import json
import os
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
    # The measured ranking performance is shown, not left to the footer.
    assert browser.evaluate("!document.getElementById('ranking-performance').hidden") is True
    assert "does not beat chance" in str(browser.evaluate(
        "document.getElementById('ranking-performance').textContent"))
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
