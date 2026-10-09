"""The version the package reports, the filename the install instructions quote, and the tag.

A colleague installs by copying a literal command containing a wheel filename. If that filename
drifts from the version the build actually produces, the command fails on a file that is not
there, and the person it fails for is the one least able to work out why.

The version also travels further than the install: every saved run records it as
``engine_version`` and the installation check prints it. Two runs produced by different code
reporting the same version is a hole in the record, not a cosmetic issue, which is why the
placeholder this project shipped with is now a real number.
"""

import re
import tomllib
from pathlib import Path

import pytest

import sarna

ROOT = Path(__file__).resolve().parents[2]

DOCUMENTS = ("README.md", "INSTALL.md")

WHEEL = re.compile(r"sarna_designer-([0-9][^-\s\"']*)-py3-none-any\.whl")


@pytest.fixture(scope="module")
def declared_version() -> str:
    """Return the version declared in ``pyproject.toml``."""
    with (ROOT / "pyproject.toml").open("rb") as handle:
        project: dict[str, object] = tomllib.load(handle)["project"]
    version = project["version"]
    assert isinstance(version, str)
    return version


def test_the_installed_package_reports_the_declared_version(declared_version: str) -> None:
    assert sarna.__version__ == declared_version


def test_the_version_is_not_a_placeholder(declared_version: str) -> None:
    # `0.1.0.dev0` never changed across every release handed to a colleague, so no saved run
    # could be traced to the code that produced it.
    assert not declared_version.endswith(".dev0")
    assert declared_version != "0+unknown"


@pytest.mark.parametrize("document", DOCUMENTS)
def test_every_wheel_filename_in_the_instructions_matches_the_build(
    document: str, declared_version: str
) -> None:
    text = (ROOT / document).read_text(encoding="utf-8")
    found = {match.group(1) for match in WHEEL.finditer(text)}
    assert found, f"{document} quotes no wheel filename; the install command needs one"
    assert found == {declared_version}, (
        f"{document} tells a colleague to install "
        f"sarna_designer-{sorted(found)[0]}-py3-none-any.whl, but this build produces "
        f"sarna_designer-{declared_version}-py3-none-any.whl. Update the documents with the "
        f"version, or the copied command fails on a file that is not there."
    )


@pytest.mark.parametrize("document", DOCUMENTS)
def test_the_instructions_say_how_to_update(document: str) -> None:
    # The question a colleague asks second, after "how do I install it".
    text = (ROOT / document).read_text(encoding="utf-8").lower()
    assert "uv tool upgrade sarna-designer" in text, (
        f"{document} does not give the update command"
    )
