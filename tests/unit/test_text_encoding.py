"""Every text file this project reads or writes must name its encoding.

Python 3.11 picks the *locale* encoding when `open`, `read_text` or `write_text` is called
without one. That is UTF-8 on Linux and macOS and **cp1252 on Windows**, so code that works
here can raise `UnicodeEncodeError` on a colleague's laptop the first time Ensembl returns a
gene description containing a Greek letter, which it does, for many genes.

The failure is invisible in development and arrives as a crash in someone else's hands, so this
is checked by parsing the source rather than by hoping a reviewer spots it. Binary modes are
exempt: they have no encoding to declare.
"""

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEXT_IO = frozenset({"read_text", "write_text", "open"})
SOURCES = sorted(
    path
    for directory in ("src", "scripts")
    for path in (ROOT / directory).rglob("*.py")
)


def _binary_mode(call: ast.Call) -> bool:
    """Return whether this call opens in a binary mode, which takes no encoding."""
    modes: list[str] = [
        a.value for a in call.args if isinstance(a, ast.Constant) and isinstance(a.value, str)
    ]
    keyword = next((k.value for k in call.keywords if k.arg == "mode"), None)
    if isinstance(keyword, ast.Constant) and isinstance(keyword.value, str):
        modes.append(keyword.value)
    return any("b" in mode for mode in modes)


def _offenders(tree: ast.AST) -> list[tuple[int, str]]:
    """Return (line, name) for each text I/O call with no explicit encoding."""
    # ZipFile.open returns a binary ZipExtFile and cannot accept an encoding argument.
    zip_handles = {
        item.optional_vars.id
        for node in ast.walk(tree) if isinstance(node, ast.With)
        for item in node.items
        if isinstance(item.optional_vars, ast.Name)
        and isinstance(item.context_expr, ast.Call)
        and isinstance(item.context_expr.func, ast.Attribute)
        and isinstance(item.context_expr.func.value, ast.Name)
        and item.context_expr.func.value.id == "zipfile"
        and item.context_expr.func.attr == "ZipFile"
    }
    found: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        name = node.func.attr
        if name not in TEXT_IO:
            continue
        # pyBigWig.open and similar third-party openers take no encoding argument.
        if isinstance(node.func.value, ast.Name) and node.func.value.id in {"pyBigWig", "gzip"}:
            continue
        if (name == "open" and isinstance(node.func.value, ast.Name)
                and node.func.value.id in zip_handles):
            continue
        if _binary_mode(node):
            continue
        if any(k.arg == "encoding" for k in node.keywords):
            continue
        found.append((node.lineno, name))
    return found


def test_there_are_sources_to_check() -> None:
    # Without this the parametrised test below could pass by checking nothing at all.
    assert len(SOURCES) > 30, len(SOURCES)


@pytest.mark.parametrize("path", SOURCES, ids=lambda p: p.name)
def test_text_io_names_its_encoding(path: Path) -> None:
    offenders = _offenders(ast.parse(path.read_text(encoding="utf-8")))
    assert not offenders, (
        f"{path.relative_to(ROOT)} calls "
        + ", ".join(f"{name}() at line {line}" for line, name in offenders)
        + ' without encoding="utf-8". On Windows this uses cp1252 and will fail on non-ASCII '
        "text that Ensembl returns for many genes."
    )


def test_cp1252_really_does_fail_on_text_ensembl_returns() -> None:
    # The premise of this whole module, asserted rather than assumed. If a future Python makes
    # UTF-8 the default everywhere (PEP 686), this test is the one to revisit first.
    with pytest.raises(UnicodeEncodeError):
        "sodium channel α-subunit".encode("cp1252")  # noqa: RUF001


def test_the_flow_diagram_and_its_svg_preview_agree() -> None:
    """The two files describe one diagram, and both used to be maintained by hand.

    That is the arrangement that let the promoter map and its legend drift until four of five
    swatches were wrong. The SVG is now generated from the .drawio; this fails if somebody
    edits one and not the other, naming the command that fixes it.
    """
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [sys.executable, str(root / "scripts" / "render_flow_diagram.py"), "--check"],
        capture_output=True, text=True, check=False, cwd=root,
    )
    assert result.returncode == 0, result.stderr or result.stdout
