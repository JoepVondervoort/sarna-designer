# Legacy scripts — preserved, not maintained

These files are a **historical record and regression oracle**. Do not edit them, do not lint
them, do not refactor them. They are excluded from linting and type checking in
`pyproject.toml`.

| File | Origin | Role |
|---|---|---|
| `sarna_promoter_scanner.py` | `saRNA Promoter Scanner.py` | **Ground truth** for design, filtering and ranking (ADR 0001). Oracle for the `legacy-scanner` preset. |
| `epimarker_checker.py` | `EpiMarker checker .py` | Epigenomic annotation / POS-vs-NON enrichment. Not a design-rule authority. |

`saRNA by Promoter Sequence.py` was supplied as a second file but is **byte-identical** to
`saRNA Promoter Scanner.py`. It is not preserved separately.

## Why these are kept

Per the project brief, legacy code is retired only once (1) functionality is mapped,
(2) regression tests exist, (3) the new engine reproduces expected behaviour, and
(4) differences are explained. None of those four are met yet. Removal requires a new ADR.

## Known defects

Catalogued with severities in [`../docs/scientific/legacy-inventory.md`](../docs/scientific/legacy-inventory.md)
(§7). The scientifically material ones:

- **B1** — minus-strand coordinates are wrong throughout `epimarker_checker.py`
- **B2** — three nearest-neighbour energy parameters are mis-paired in `sarna_promoter_scanner.py`
- **B3** — operator-precedence bug returns a tuple where a score is expected

These are **deliberately not fixed here.** Fixes belong in the new engine, where the behaviour
change is documented and switchable.

## Running them

Neither script runs cleanly in WSL as-is:

- `sarna_promoter_scanner.py` needs `matplotlib` (only for a stray unused `from pylab import *`)
  and reads the promoter from `input()` on stdin.
- `epimarker_checker.py` hardcodes a Windows path (`C:\Users\joep-\...`), calls `os.chdir` on it,
  `pip install`s at runtime, and requires a missing `Designs.xlsx`.

To reproduce legacy Scanner output for regression comparison, use the harness at
`tests/fixtures/` rather than invoking these directly.
