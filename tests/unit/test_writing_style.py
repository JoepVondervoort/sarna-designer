"""Check repository editorial conventions and prohibited text encodings.

Documentation, interface text and source comments use consistent punctuation and concrete
language. Historical legacy scripts and generated output are excluded.
"""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

SUFFIXES = frozenset({
    ".md", ".py", ".html", ".css", ".js", ".toml", ".txt", ".cfg", ".ini", ".yml", ".yaml",
    ".json", ".example", ".csv",
})

SKIP_DIRECTORIES = frozenset({
    "legacy", ".venv", ".git", ".cache", "__pycache__", ".mypy_cache", ".pytest_cache",
    "node_modules", "dist", "build", ".ruff_cache", ".sarna-designer",
    # Generated pages, rebuilt from the swept template by scripts/build_demo_views.py and
    # ignored by git. Checking them reports the age of someone's last build, not the source.
    "out",
})

ALLOW = "writing-style: allow"
"""A line marked with this is exempt, for text that is *about* the rule rather than breaking
it. This file and the working notes both have to name the characters they forbid. The marker
is per line and deliberately conspicuous, so an exemption cannot be granted quietly."""

# Built rather than written out, so this file does not itself contain any of the things it
# bans. The test already skips itself, but a plain `git grep` over the repository does not
# know that, and "search for the character and expect nothing" is the check a person actually
# runs. Leaving one literal here made that check report a false positive.
EM_DASH = chr(0x2014)

DISGUISES = (
    ("\\" + "u2014", "a Python or JavaScript escape for an em dash"),
    ("&m" + "dash;", "the HTML entity for an em dash"),
    ("&#" + "8212;", "the numeric HTML entity for an em dash"),
)

PHRASES = {
    "delve": "say what you mean: look at, read, examine",
    "tapestry": "not a word this project has any use for",
    "testament to": "state the evidence instead",
    "seamless": "claims an absence of friction nobody measured",
    "in the realm of": "say in, or for",
    "it is worth noting": "if it is worth noting, note it",
    "it's worth noting": "if it is worth noting, note it",
    "navigating the complexities": "say what is complicated and why",
    "unlock the potential": "say what it does",
    "game-chang": "a claim, not a description",
}
"""Deliberately short. Every entry has no legitimate use in this codebase, which is why it can
be banned outright. Words that are also legitimate technical terms, `harness` for a
test harness and `comprehensive` in a caveat about what a screening is not, are left alone: a
check that cries wolf gets suppressed, and then it checks nothing."""


def text_files() -> list[Path]:
    """Return every human-readable file this rule applies to."""
    found: list[Path] = []
    for path in ROOT.rglob("*"):
        if path.relative_to(ROOT).as_posix() in {
            "CLAUDE.md", "docs/development/continuation-plan.md",
        }:
            continue
        if not path.is_file() or path.suffix.lower() not in SUFFIXES:
            continue
        if SKIP_DIRECTORIES & set(path.relative_to(ROOT).parts):
            continue
        found.append(path)
    return sorted(found)


@pytest.fixture(scope="module")
def files() -> list[Path]:
    return text_files()


def offending_lines(path: Path, needle: str) -> list[str]:
    """Return the lines of ``path`` holding ``needle`` that are not marked exempt."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        shown = path.relative_to(ROOT)
    except ValueError:        # A path outside the tree, which only the marker test passes in.
        shown = path
    return [f"    {shown}:{number}: {line.strip()[:110]}"
            for number, line in enumerate(lines, 1)
            if needle in line and ALLOW not in line]


def scan(files: list[Path], needle: str) -> list[str]:
    """Return every unexempt occurrence of ``needle``, this file excluded.

    This test states the rule, so it necessarily contains what it bans.
    """
    found: list[str] = []
    for path in files:
        if path.resolve() == Path(__file__).resolve():
            continue
        found += offending_lines(path, needle)
    return found


def test_the_sweep_actually_covers_the_repository(files: list[Path]) -> None:
    # A rule that silently matched nothing would pass for ever. These are the files a reader
    # meets first, so if the walk stops finding them the walk is broken, not the repository.
    names = {path.relative_to(ROOT).as_posix() for path in files}
    assert {"README.md", "INSTALL.md", "GLOSSARY.md"} <= names
    assert len(files) > 100, f"only {len(files)} files were scanned, which cannot be right"


def test_no_em_dash_outside_legacy(files: list[Path]) -> None:
    offenders = scan(files, EM_DASH)
    assert not offenders, (
        "Repository style uses commas, colons, full stops and parentheses instead of "
        "em dashes. Replace the character with a comma, a colon, a full stop or "
        "parentheses, whichever the sentence actually wants:\n"
        + "\n".join(offenders[:12])
    )


@pytest.mark.parametrize(("disguise", "what"), DISGUISES)
def test_no_escaped_em_dash(files: list[Path], disguise: str, what: str) -> None:
    # The first sweep missed seven of these in app.js, because grepping for the character
    # does not find the escape that renders as it.
    offenders = scan(files, disguise)
    assert not offenders, (
        f"{disguise} is {what}, and reads to a user exactly like the character:\n"
        + "\n".join(offenders[:12])
    )


def test_no_banned_phrase(files: list[Path]) -> None:
    found: list[str] = []
    for phrase, why in PHRASES.items():
        for path in files:
            if path.resolve() == Path(__file__).resolve():
                continue
            lowered = path.read_text(encoding="utf-8", errors="replace").lower()
            if phrase in lowered:
                found.append(f"  {path.relative_to(ROOT)}: {phrase!r}, {why}")
    assert not found, "Phrasing outside the editorial conventions:\n" + "\n".join(found[:12])


def test_the_exemption_marker_is_honoured_and_narrow(tmp_path: Path) -> None:
    # Planted on purpose: a marker that did not work would make the whole check decorative,
    # and one that worked per file rather than per line would exempt too much.
    sample = tmp_path / "sample.md"
    sample.write_text(f"exempt {EM_DASH} line  {ALLOW}\nplain {EM_DASH} line\n",
                      encoding="utf-8")
    hits = offending_lines(sample, EM_DASH)
    assert len(hits) == 1
    assert "sample.md:2" in hits[0]
