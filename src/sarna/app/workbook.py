"""Write a minimal `.xlsx` workbook with nothing but the standard library.

A CSV would be less code, and the app already serves two. The reason this exists is that a CSV
is not a file format, it is a convention, and Excel resolves it against the machine's locale:
on a Dutch or German Windows the list separator is ``;``, so a comma-separated export opens
with every column crushed into column A. That is the first thing a colleague would see, and
they would reasonably conclude the tool is broken.

An `.xlsx` has no such ambiguity. It also keeps a sequence a string: left to guess, Excel reads
some values as numbers and will happily turn a target id into a date.

**Deliberately not a spreadsheet library.** No formulas, no styling beyond a bold header, column
widths and a frozen top row, no shared-string table — every cell is an inline string or a
number. That is enough for an export and it keeps `pydantic` the only runtime dependency, which
is what makes :doc:`the install </INSTALL>` two commands.

The parts written here are the minimum an OOXML consumer requires: a content-type map, the
package relationships, a workbook naming its sheets, the workbook's relationships, and one XML
part per sheet.
"""

import re
import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from io import BytesIO
from typing import Any

MAX_SHEET_NAME = 31
"""Excel refuses a longer one, silently in some versions."""

_FORBIDDEN_IN_SHEET_NAME = re.compile(r"[\[\]:*?/\\]")

_CONTROL_CHARACTERS = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
"""Illegal in XML 1.0 even escaped. Stripped rather than escaped, or the file will not open."""

_NAMESPACE = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_RELATIONSHIPS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


@dataclass(frozen=True, slots=True)
class Sheet:
    """One worksheet: a name, a header row, and the rows under it."""

    name: str
    headers: tuple[str, ...]
    rows: tuple[tuple[Any, ...], ...]

    def safe_name(self) -> str:
        """Return the name with the characters Excel refuses removed and the length capped."""
        cleaned = _FORBIDDEN_IN_SHEET_NAME.sub(" ", self.name).strip() or "Sheet"
        return cleaned[:MAX_SHEET_NAME]


def column_letter(index: int) -> str:
    """Return the spreadsheet column letter for a 1-based index: 1 -> A, 27 -> AA."""
    if index < 1:
        raise ValueError(f"column indices are 1-based, got {index}")
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(ord("A") + remainder) + letters
    return letters


def _escape(value: str) -> str:
    """Escape text for XML content, dropping characters XML cannot carry at all."""
    text = _CONTROL_CHARACTERS.sub("", value)
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _cell(reference: str, value: Any, style: int = 0) -> str:  # noqa: ANN401
    """Return one cell, numeric where that is unambiguous and text everywhere else.

    ``bool`` is checked before ``int`` because it is one, and writing ``1`` where a reader
    expects ``TRUE`` loses the distinction. Everything that is not plainly a finite number
    becomes an inline string, which is what keeps a target site a target site rather than
    something Excel reinterprets.
    """
    styled = f' s="{style}"' if style else ""
    if value is None:
        return f'<c r="{reference}"{styled} t="inlineStr"><is><t/></is></c>'
    if isinstance(value, bool):
        word = "TRUE" if value else "FALSE"
        return f'<c r="{reference}"{styled} t="inlineStr"><is><t>{word}</t></is></c>'
    if isinstance(value, int | float):
        number = float(value)
        if number == number and abs(number) != float("inf"):  # NaN and infinities are text
            return f'<c r="{reference}"{styled}><v>{value!r}</v></c>'
    escaped = _escape(str(value))
    # xml:space preserves leading and trailing spaces, which a reader would otherwise trim.
    return (f'<c r="{reference}"{styled} t="inlineStr">'
            f'<is><t xml:space="preserve">{escaped}</t></is></c>')


def _row(number: int, values: tuple[Any, ...], style: int = 0) -> str:
    """Return one row of cells, numbered from 1."""
    cells = "".join(
        _cell(f"{column_letter(column)}{number}", value, style)
        for column, value in enumerate(values, start=1)
    )
    return f'<row r="{number}">{cells}</row>'


def _widths(sheet: Sheet) -> str:
    """Return column widths wide enough to read a 21-mer without dragging every border."""
    if not sheet.headers:
        return ""
    widths = []
    for index, header in enumerate(sheet.headers, start=1):
        longest = max(
            [len(header)] + [len(str(row[index - 1])) for row in sheet.rows if len(row) >= index]
        )
        widths.append(
            f'<col min="{index}" max="{index}" width="{min(max(longest + 2, 10), 60)}" '
            f'customWidth="1"/>'
        )
    return f"<cols>{''.join(widths)}</cols>"


def _sheet_xml(sheet: Sheet) -> str:
    """Return one worksheet part.

    Element order follows the schema: ``dimension``, ``sheetViews``, ``cols``, ``sheetData``,
    ``autoFilter``. Out of order, Excel reports the file as corrupt rather than ignoring it.
    """
    last_column = column_letter(max(1, len(sheet.headers)))
    last_row = len(sheet.rows) + 1
    body = _row(1, sheet.headers, style=1) + "".join(
        _row(number, values) for number, values in enumerate(sheet.rows, start=2)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<worksheet xmlns="{_NAMESPACE}">'
        f'<dimension ref="A1:{last_column}{last_row}"/>'
        '<sheetViews><sheetView workbookViewId="0">'
        # The header stays put while a researcher scrolls a few hundred candidates.
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        "</sheetView></sheetViews>"
        f"{_widths(sheet)}"
        f"<sheetData>{body}</sheetData>"
        f'<autoFilter ref="A1:{last_column}{last_row}"/>'
        "</worksheet>"
    )


_STYLES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
    '<fills count="2"><fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill></fills>'
    '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
    '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
    # A reader that finds no "Normal" style substitutes its own and warns that the workbook has
    # none. openpyxl says so out loud; Excel is quieter about it. Declared rather than relying
    # on either being forgiving.
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
    "</styleSheet>"
)


def write_workbook(sheets: Sequence[Sheet]) -> bytes:
    """Return a complete ``.xlsx`` file as bytes.

    Args:
        sheets: One or more sheets, in the order they should appear.

    Returns:
        The whole package, ready to serve or write to disk.

    Raises:
        ValueError: If no sheet was supplied, or two share a name after cleaning. Excel opens a
            workbook with duplicate sheet names as corrupt, so this refuses rather than
            producing a file that fails later on somebody else's machine.
    """
    if not sheets:
        raise ValueError("a workbook needs at least one sheet")
    names = [sheet.safe_name() for sheet in sheets]
    if len(set(names)) != len(names):
        raise ValueError(f"sheet names must be unique after cleaning, got {names}")

    overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
        f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet'
        f'+xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    sheet_entries = "".join(
        f'<sheet name="{_escape(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, name in enumerate(names, start=1)
    )
    sheet_rels = "".join(
        f'<Relationship Id="rId{index}" Type="{_RELATIONSHIPS}/worksheet" '
        f'Target="worksheets/sheet{index}.xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    styles_id = len(sheets) + 1

    buffer = BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.'
            'relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.'
            'openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.'
            'openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            f"{overrides}</Types>")
        archive.writestr("_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
            'relationships">'
            f'<Relationship Id="rId1" Type="{_RELATIONSHIPS}/officeDocument" '
            'Target="xl/workbook.xml"/></Relationships>')
        archive.writestr("xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<workbook xmlns="{_NAMESPACE}" xmlns:r="{_RELATIONSHIPS}">'
            f"<sheets>{sheet_entries}</sheets></workbook>")
        archive.writestr("xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
            'relationships">'
            f"{sheet_rels}"
            f'<Relationship Id="rId{styles_id}" Type="{_RELATIONSHIPS}/styles" '
            'Target="styles.xml"/></Relationships>')
        archive.writestr("xl/styles.xml", _STYLES)
        for index, sheet in enumerate(sheets, start=1):
            archive.writestr(f"xl/worksheets/sheet{index}.xml", _sheet_xml(sheet))
    return buffer.getvalue()
