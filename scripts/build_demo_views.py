"""Build offline default/custom previews, then serve out/ with Python's HTTP server.

Uses synthetic DNA with no genomic coordinates or TSS. This is an interface demonstration,
not an experiment or an efficacy benchmark. No external lookups or track downloads are needed.

    .venv/bin/python scripts/build_demo_views.py
    .venv/bin/python -m http.server 8000 --bind 127.0.0.1 --directory out
"""

import argparse
import json
from pathlib import Path
from shlex import quote

from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.promoter import Promoter
from sarna.rules.rule import RuleKind
from sarna.sequence.nucleic import NucleicSequence
from sarna.view.html import FALLBACK_TITLE, render
from sarna.view.promoter_view import PromoterView

# A CpG-depleted synthetic sequence also used by the browser regression fixtures.
DEMO_BASES = (
    "GCTACGGTATCATGTCATAAGTTCAATTGACATTGGAATTACCATAGTTAACGTTTAGCATA" * 10
)[:600]


def main(argv: list[str] | None = None) -> int:
    """Write two reproducible demo pages, their JSON records and a preview index."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("out"))
    args = parser.parse_args(argv)
    destination: Path = args.out_dir
    destination.mkdir(parents=True, exist_ok=True)

    promoter = Promoter(sequence=NucleicSequence.parse(DEMO_BASES))
    configurations = {
        "demo-default": WANG_2015.resolve(),
        "demo-custom": WANG_2015.resolve(
            {"gc-content": {"maximum-gc-fraction": 0.65}},
            role_overrides={"pos19-a": RuleKind.PREFERENCE, "tss-window": None},
        ),
    }
    for name, config in configurations.items():
        run = design_sarnas(promoter, config)
        view = PromoterView(result=run)
        page = render(view)
        note = (
            '<p class="notice"><b>Synthetic demonstration.</b> These 600 bases are not a '
            'real gene. No genomic coordinates, TSS, expression or chromatin data are supplied. '
            'Uncheckable criteria remain unknown; candidate ranks do not predict activation. '
            '<a href="index.html">Back to previews</a>.</p>'
        )
        # The note stays outside #app so the page's existing rerender cannot erase it.
        page = page.replace('<div class="wrap">', '<div class="wrap">' + note, 1)
        page = page.replace(FALLBACK_TITLE, f"Synthetic saRNA demo: {name}", 1)
        (destination / f"{name}.html").write_text(page, encoding="utf-8")
        payload = view.to_dict()
        payload["demonstration"] = {"synthetic": True, "purpose": "interface testing"}
        (destination / f"{name}.json").write_text(
            json.dumps(payload, indent=1), encoding="utf-8"
        )
        print(
            f"{name}: {len(run.generated)} generated, {len(run.passing)} passing; "
            f"{config.max_preference_points} configured preferences"
        )

    real_gene = (
        '<li><a href="scn1a.html">SCN1A with Wang defaults</a>, cached real promoter, '
        'MANE Select anchor; histone tracks skipped in this preview.</li>'
        if (destination / "scn1a.html").is_file() else ""
    )
    index = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>saRNA local previews</title>
<style>
body {{font:17px/1.6 system-ui,sans-serif;max-width:820px;margin:48px auto;padding:0 24px;
       background:#f3f4f1;color:#16191c;}}
a {{color:#0f6e76;}} li {{margin:16px 0;}} code {{background:#e6e9e3;padding:2px 5px;}}
</style></head><body>
<h1>saRNA local previews</h1>
<p>These are frozen design records with interactive viewing. Settings are changed in Python
and regenerated; the settings editor and saved-run workbench are still planned.</p>
<ul>
{real_gene}
<li><a href="demo-default.html">Synthetic demo: Wang defaults</a> (all default settings,
with location-dependent checks explicitly unavailable.</li>
<li><a href="demo-custom.html">Synthetic demo: custom settings</a>) upper GC bound 65%,
position-19 as a preference, TSS window deliberately off. See the recorded custom departures.</li>
</ul>
<p>The two synthetic pages use identical DNA. The custom page changes three settings together;
the overall difference cannot be attributed to just one of them.</p>
<h2>Try it</h2>
<ol><li>Click a candidate and read its sequences and rule explanations.</li>
<li>Switch between Passing and All to see rejected candidates and their reasons.</li>
<li>Zoom or drag the promoter track and toggle motif predictions where available.</li>
<li>Read Rules and their sources; compare the baseline, custom roles and disabled rules.</li>
<li>Check that absent cell evidence is described as unavailable rather than negative.</li></ol>
<p>Rerun <code>scripts/build_demo_views.py</code> after changing the example configuration,
then refresh the browser. JSON records are beside each HTML file.</p>
</body></html>
"""
    (destination / "index.html").write_text(index, encoding="utf-8")
    print(f"\nOpen {destination / 'index.html'} in a browser, or serve this directory:")
    print(
        "  .venv/bin/python -m http.server 8000 --bind 127.0.0.1 --directory "
        + quote(str(destination))
    )
    print("  then open http://localhost:8000/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
