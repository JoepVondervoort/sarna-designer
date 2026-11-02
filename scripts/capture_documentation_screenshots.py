"""Capture real, unscreened SCN1A views for the repository documentation."""

import argparse
import json
import subprocess
import time
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread
from typing import Any, cast
from urllib.request import urlopen
from zoneinfo import ZoneInfo

from verify_workbench_browser import Browser

from sarna.app.service import DesignService
from sarna.app.web import AppRouter, make_server


def capture(browser: Browser, url: str, output: Path) -> None:
    """Use the current interface and cached public references without submitting screening."""
    browser.call('Runtime.enable')
    browser.call('Emulation.setDeviceMetricsOverride', width=1440, height=1100,
                 deviceScaleFactor=1, mobile=False)
    browser.call('Page.navigate', url=url)
    browser.wait("document.getElementById('preset')?.options.length===4")
    browser.evaluate("document.getElementById('offline').checked=true;"
                     "document.getElementById('offline').dispatchEvent(new Event('change'));"
                     "document.getElementById('gene-symbol').value='SCN1A';"
                     "document.getElementById('lookup-gene').click()")
    browser.wait("state.gene?.symbol==='SCN1A' && "
                 "document.getElementById('reference').value!=='' && "
                 "!document.getElementById('lookup-gene').disabled")
    browser.screenshot(output / 'design-setup.png')
    browser.evaluate("document.getElementById('design-form').requestSubmit()")
    browser.wait("state.result && !state.busy && "
                 "!document.getElementById('design-ready').hidden")
    assert browser.evaluate("document.getElementById('global-error').textContent") == ''
    browser.evaluate("setWorkspaceTab('candidates',{scroll:false});window.scrollTo(0,0)")
    browser.wait("state.investigationRunId===state.result.app.run_id")
    browser.evaluate("new Promise(resolve=>requestAnimationFrame("
                     "()=>requestAnimationFrame(resolve)))")
    browser.screenshot(output / 'promoter-landscape.png')
    facts = cast(dict[str, Any], browser.evaluate("""(()=>{
      const r=state.result, best=r.candidates.find(c=>c.rank===1);
      const {sequence,...promoter}=r.promoter;
      return {engine_version:r.app.engine_version,run_id:r.app.run_id,
        created:r.app.created,gene:r.gene?.symbol,promoter,
        reference:r.app.anchor,request:r.app.request,counts:r.counts,
        configuration:r.configuration.preset,is_custom:r.configuration.is_custom,
        best:best?{id:best.id,tss_span:best.tss_span,target:best.target_site_dna}:null,
        motifs:r.tracks.motifs.length,regulatory_features:r.tracks.regulatory.length};
    })()"""))
    assert facts['is_custom'] is False
    browser.evaluate("""(async()=>{
      await selectCandidate(state.result.candidates.find(c=>c.rank===1).id);
      document.querySelector('[data-inspector-section=duplex]').open=true;
      document.getElementById('candidate-detail').scrollIntoView({block:'start'});
    })()""")
    browser.call('Emulation.setDeviceMetricsOverride', width=1280, height=1250,
                 deviceScaleFactor=1, mobile=False)
    browser.evaluate("new Promise(resolve=>requestAnimationFrame("
                     "()=>requestAnimationFrame(resolve)))")
    browser.evaluate("document.getElementById('candidate-detail').scrollIntoView({block:'start'})")
    browser.screenshot(output / 'duplex-inspection.png')
    browser.call('Emulation.setDeviceMetricsOverride', width=1440, height=1100,
                 deviceScaleFactor=1, mobile=False)
    browser.evaluate("setWorkspaceTab('design',{scroll:false});"
                     "document.getElementById('advanced-design-settings').open=true;"
                     "document.getElementById('open-settings').click();"
                     "document.querySelector('[data-rule-help=\"gc-content\"]').open=true;"
                     "document.querySelector('[data-rule-params=\"gc-content\"]').open=true;"
                     "document.querySelector('[data-rule=\"gc-content\"]')"
                     ".scrollIntoView({block:'center'})")
    browser.screenshot(output / 'design-rules.png')
    browser.evaluate("document.getElementById('settings-dialog').close();"
                     "setWorkspaceTab('candidates',{scroll:false})")
    browser.evaluate("""(()=>{
      const ids=state.result.candidates.filter(c=>c.rank!==null && c.rank<=3).map(c=>c.id);
      for(const id of ids)document.querySelector(`[data-check="${CSS.escape(id)}"]`).click();
      setWorkspaceTab('export',{scroll:false});window.scrollTo(0,0);
    })()""")
    browser.screenshot(output / 'export-options.png')
    browser.evaluate("setWorkspaceTab('candidates',{scroll:false});"
                     "document.querySelector('[data-map-tool=map-search-tool]').click();"
                     "document.getElementById('map-query').value="
                     "state.result.candidates.find(c=>c.rank===1).target_site_dna;"
                     "document.getElementById('map-query-orientation').value='same';"
                     "document.getElementById('map-find').click()")
    browser.wait("!document.getElementById('map-find').disabled && state.mapSearchMatches.length>0")
    browser.evaluate("document.getElementById('map-search-tool').scrollIntoView({block:'center'});"
                     "document.getElementById('map-query').focus({preventScroll:true})")
    browser.screenshot(output / 'sequence-search.png')
    browser.evaluate("setWorkspaceTab('screening',{scroll:false});window.scrollTo(0,0)")
    browser.screenshot(output / 'screening-workspace.png')
    assert not browser.exceptions, browser.exceptions
    facts.update(
        captured_at=datetime.now(ZoneInfo('Europe/Amsterdam')).isoformat(),
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        browser=browser.call('Browser.getVersion')['product'],
        viewport={'width': 1440, 'height': 1100, 'device_scale_factor': 1},
        image_viewports={'duplex-inspection.png': {'width': 1280, 'height': 1250}},
        screening_performed=False,
        files=['design-setup.png', 'promoter-landscape.png', 'duplex-inspection.png',
               'design-rules.png', 'export-options.png', 'sequence-search.png',
               'screening-workspace.png'],
    )
    (output / 'capture.json').write_text(json.dumps(facts, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(facts, indent=2))


def main() -> None:
    """Start an isolated local app and Chromium, capture, and clean up both processes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chrome', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('docs/assets/screenshots'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='sarna-documentation-') as temporary:
        directory = Path(temporary)
        service = DesignService(directory / 'runs')
        router = AppRouter(service)
        server = make_server(router, 0)
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        process = subprocess.Popen([
            str(args.chrome), '--headless', '--no-sandbox', '--disable-gpu',
            '--disable-dev-shm-usage', '--disable-background-networking',
            '--remote-debugging-port=0', f'--user-data-dir={directory / "profile"}', 'about:blank',
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        browser = None
        try:
            active = directory / 'profile' / 'DevToolsActivePort'
            deadline = time.monotonic() + 15
            while not active.exists():
                if process.poll() is not None or time.monotonic() > deadline:
                    raise RuntimeError('Chromium did not start its local DevTools endpoint.')
                time.sleep(0.1)
            port = int(active.read_text(encoding='utf-8').splitlines()[0])
            with urlopen(f'http://127.0.0.1:{port}/json/list') as response:
                pages = json.load(response)
            browser = Browser(pages[0]['webSocketDebuggerUrl'])
            capture(browser, f'http://127.0.0.1:{server.server_port}/', args.output)
            metadata = json.loads(
                (args.output / 'capture.json').read_text(encoding='utf-8'),
            )
            run_id = metadata['run_id']
            assert service.store.screenings_for(run_id) == []
        finally:
            if browser:
                browser.socket.close()
            process.terminate()
            process.wait(timeout=10)
            server.shutdown()
            server.server_close()
            router.close()
            thread.join(timeout=5)


if __name__ == '__main__':
    main()
