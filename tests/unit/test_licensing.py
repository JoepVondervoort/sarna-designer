"""Licensing invariants, pinned so attribution cannot rot silently.

Apache-2.0 section 4(d) requires the `NOTICE` file to travel with any redistribution, and
`NOTICE` is only useful if it actually lists what is bundled. A dataset added to the wheel
without a line here is a licence violation that no other test would notice, so these assertions
compare the two directly rather than trusting either on its own.

[ADR 0018](../../docs/adr/0018-apache-2-and-bundled-data-licensing.md) records the decisions.
"""

import json
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
NOTICE = (ROOT / "NOTICE").read_text(encoding="utf-8")
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))

# Every data file hatchling ships inside the wheel, by the artifact globs in pyproject.
BUNDLED_DATA = sorted(
    path
    for pattern in ("src/sarna/*/data/*.csv", "src/sarna/*/data/*.toml", "src/sarna/*/data/*.json")
    for path in ROOT.glob(pattern)
)


def test_the_repository_is_apache_two_point_oh() -> None:
    assert PYPROJECT["project"]["license"] == "Apache-2.0"
    licence = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "Apache License" in licence
    assert "Version 2.0, January 2004" in licence
    # The placeholder this replaced granted no rights at all; make its return loud.
    assert "LICENSING UNDECIDED" not in licence
    assert "all rights reserved" not in licence.lower()


def test_the_notice_file_ships_in_the_distribution() -> None:
    # Apache-2.0 4(d): a redistributor must be able to pass NOTICE on, so it has to be in the
    # wheel rather than only in the source tree.
    assert "NOTICE" in PYPROJECT["project"]["license-files"]
    assert "LICENSE" in PYPROJECT["project"]["license-files"]


def test_the_bundled_data_is_actually_findable() -> None:
    # Guards the test below: if the globs stop matching, the attribution check would pass
    # vacuously by comparing two empty sets.
    assert len(BUNDLED_DATA) >= 5, [p.name for p in BUNDLED_DATA]


@pytest.mark.parametrize("path", BUNDLED_DATA, ids=lambda p: p.name)
def test_every_bundled_data_file_is_attributed_in_notice(path: Path) -> None:
    relative = path.relative_to(ROOT).as_posix()
    assert relative in NOTICE, (
        f"{relative} ships in the wheel but NOTICE does not mention it. Add its source and "
        f"licence there before distributing, as Apache-2.0 4(d) requires."
    )


def test_jaspar_is_recorded_as_cc_by_four_not_the_non_commercial_licence() -> None:
    # The repository previously claimed CC BY-NC-SA 4.0. JASPAR's own about page states CC BY
    # 4.0, so the old label was wrong *and more restrictive than reality*, it would have made
    # users refuse permitted uses. This pins the correction in both places it is stated.
    data = json.loads(
        (ROOT / "src/sarna/regulatory/data/jaspar_core_selection.json").read_text(encoding="utf-8")
    )
    assert "CC BY 4.0" in data["licence"]
    assert "NC-SA" not in data["licence"]
    assert data["attribution"]
    assert "CC BY 4.0" in NOTICE
    assert "jaspar.elixir.no" in NOTICE


def test_notice_separates_bundled_data_from_services_queried_at_runtime() -> None:
    # Ensembl and NCBI are queried, not redistributed. Listing them as bundled would claim a
    # redistribution right nobody exercised; omitting them would hide where data leaves the
    # machine. Both sections must exist.
    assert "Bundled third-party data" in NOTICE
    assert "External services used at runtime (not redistributed)" in NOTICE
    assert "rest.ensembl.org" in NOTICE
    assert "blast.ncbi.nlm.nih.gov" in NOTICE


def test_no_web_framework_is_advertised_as_an_extra() -> None:
    # ADR 0015 chose a stdlib HTTP adapter. An `api` extra pulling fastapi and uvicorn told
    # readers there was a web API; there is not, and nothing imported them.
    extras = PYPROJECT["project"]["optional-dependencies"]
    declared = {name for group in extras.values() for name in group}
    assert not any("fastapi" in d or "uvicorn" in d for d in declared), declared
