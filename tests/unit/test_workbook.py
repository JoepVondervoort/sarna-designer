"""The stdlib ``.xlsx`` writer, read back by a reader that shares none of its code.

The point of these tests is that the grid written is the grid an OOXML consumer gets out. They
therefore do not call anything in :mod:`sarna.app.workbook` to read: :func:`grid_of` parses the
package the way the specification says to, from the content types through the relationships to
the sheet parts, so a writer that agreed only with itself would fail here.
"""

import zipfile
from io import BytesIO
from typing import Any
from xml.etree import ElementTree

import pytest

from sarna.app.workbook import MAX_SHEET_NAME, Sheet, column_letter, write_workbook

MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
PACKAGE_RELS = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def sheets_of(data: bytes) -> dict[str, list[list[Any]]]:
    """Read a workbook the way a consumer does: rels, then workbook, then each sheet part.

    Resolves the sheet name to its part through ``r:id``, rather than assuming ``sheet1.xml``
    is the first sheet, so a mis-wired relationship shows up as a failure here.
    """
    archive = zipfile.ZipFile(BytesIO(data))
    rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    target_of = {
        r.get("Id"): r.get("Target") for r in rels.findall(f"{PACKAGE_RELS}Relationship")
    }
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    relationship = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    read: dict[str, list[list[Any]]] = {}
    for entry in workbook.iter(f"{MAIN}sheet"):
        name = entry.get("name")
        part = target_of[entry.get(relationship)]
        assert name is not None and part is not None
        read[name] = _grid(archive.read(f"xl/{part}"))
    return read


def _grid(part: bytes) -> list[list[Any]]:
    """Return one sheet's cells as a list of rows."""
    sheet = ElementTree.fromstring(part)
    rows = []
    for row in sheet.iter(f"{MAIN}row"):
        values: list[Any] = []
        for cell in row.findall(f"{MAIN}c"):
            if cell.get("t") == "inlineStr":
                text = cell.find(f"{MAIN}is/{MAIN}t")
                values.append("" if text is None or text.text is None else text.text)
            else:
                number = cell.find(f"{MAIN}v")
                values.append(None if number is None else float(number.text or "0"))
        rows.append(values)
    return rows


class TestTheGridSurvivesTheRoundTrip:
    def test_headers_and_rows_come_back_in_order(self) -> None:
        data = write_workbook([Sheet("Candidates", ("rank", "id"), ((1, "saRNA-1"),
                                                                    (2, "saRNA-2")))])
        assert sheets_of(data) == {
            "Candidates": [["rank", "id"], [1.0, "saRNA-1"], [2.0, "saRNA-2"]]
        }

    def test_several_sheets_keep_their_names_and_their_order(self) -> None:
        data = write_workbook([
            Sheet("Candidates", ("a",), (("x",),)),
            Sheet("Run", ("b",), (("y",),)),
            Sheet("Rules", ("c",), (("z",),)),
        ])
        assert list(sheets_of(data)) == ["Candidates", "Run", "Rules"]

    def test_a_sequence_stays_text_rather_than_becoming_a_number(self) -> None:
        # The whole reason for not shipping a CSV. Excel reads an unquoted "1E5" as 100000 and
        # a target id as a date; an inline string has no such latitude.
        data = write_workbook([Sheet("S", ("seq",), (("GUGAAGCCAUAGAUCUGUA",), ("1E5",)))])
        grid = sheets_of(data)["S"]
        assert grid[1] == ["GUGAAGCCAUAGAUCUGUA"]
        assert grid[2] == ["1E5"], "a sequence-like string must not be written as a number"

    def test_numbers_are_numbers_so_a_column_can_be_sorted(self) -> None:
        data = write_workbook([Sheet("S", ("gc",), ((42.1,), (57.9,)))])
        assert sheets_of(data)["S"][1:] == [[42.1], [57.9]]

    def test_a_missing_value_is_empty_rather_than_the_word_none(self) -> None:
        data = write_workbook([Sheet("S", ("rank",), ((None,),))])
        assert sheets_of(data)["S"][1] == [""]

    def test_booleans_are_words_rather_than_ones_and_zeroes(self) -> None:
        # bool is an int in Python, and writing 1 loses the distinction a reader needs.
        data = write_workbook([Sheet("S", ("agree",), ((True,), (False,)))])
        assert sheets_of(data)["S"][1:] == [["TRUE"], ["FALSE"]]


class TestTheFileIsAWellFormedPackage:
    def test_every_part_the_specification_requires_is_present(self) -> None:
        data = write_workbook([Sheet("One", ("a",), (("x",),))])
        names = set(zipfile.ZipFile(BytesIO(data)).namelist())
        assert {"[Content_Types].xml", "_rels/.rels", "xl/workbook.xml",
                "xl/_rels/workbook.xml.rels", "xl/styles.xml",
                "xl/worksheets/sheet1.xml"} <= names

    def test_every_part_is_well_formed_xml(self) -> None:
        data = write_workbook([Sheet("One", ("a&b", "c<d"), (('say "hi"', "x>y"),))])
        archive = zipfile.ZipFile(BytesIO(data))
        for name in archive.namelist():
            ElementTree.fromstring(archive.read(name))  # raises if not well formed

    def test_content_types_declares_every_sheet_part(self) -> None:
        # A sheet part with no declared content type opens as a corrupt workbook.
        data = write_workbook([Sheet(f"S{i}", ("a",), (("x",),)) for i in range(1, 4)])
        declared = ElementTree.fromstring(
            zipfile.ZipFile(BytesIO(data)).read("[Content_Types].xml"))
        parts = {o.get("PartName") for o in declared.iter()
                 if o.tag.endswith("Override")}
        assert {f"/xl/worksheets/sheet{i}.xml" for i in (1, 2, 3)} <= parts

    def test_xml_special_characters_survive_rather_than_breaking_the_file(self) -> None:
        data = write_workbook([Sheet("One", ("note",), (('a & b < c > d "e"',),))])
        assert sheets_of(data)["One"][1] == ['a & b < c > d "e"']

    def test_control_characters_are_dropped_because_xml_cannot_carry_them(self) -> None:
        # Escaping does not help: \x00 is illegal in XML 1.0 in any form, and a file
        # containing one does not open at all.
        data = write_workbook([Sheet("One", ("note",), (("before\x00after",),))])
        assert sheets_of(data)["One"][1] == ["beforeafter"]


    def test_the_header_row_is_styled_on_its_cells_not_on_the_row(self) -> None:
        """A style on ``<row>`` does not reach the cells without ``customFormat``.

        The header rendered unbolded for exactly that reason, and no test here could see it:
        the grid round-tripped perfectly either way. openpyxl reported it, which is why the
        writer is checked against a real OOXML implementation and not only against itself.
        """
        data = write_workbook([Sheet("One", ("a", "b"), (("x", "y"),))])
        sheet = ElementTree.fromstring(
            zipfile.ZipFile(BytesIO(data)).read("xl/worksheets/sheet1.xml"))
        rows = list(sheet.iter(f"{MAIN}row"))
        header_cells = rows[0].findall(f"{MAIN}c")
        assert all(cell.get("s") == "1" for cell in header_cells), "the header is not styled"
        assert rows[1].findall(f"{MAIN}c")[0].get("s") is None, "a body cell must not be bold"

    def test_the_stylesheet_declares_a_default_style(self) -> None:
        # Without a "Normal" cellStyle a reader substitutes its own and warns the workbook has
        # none; Excel is quieter about it than openpyxl, which is how this was found.
        data = write_workbook([Sheet("One", ("a",), (("x",),))])
        styles = zipfile.ZipFile(BytesIO(data)).read("xl/styles.xml").decode()
        assert 'name="Normal"' in styles


class TestTheRefusals:
    def test_a_workbook_with_no_sheets_is_refused(self) -> None:
        with pytest.raises(ValueError, match="at least one sheet"):
            write_workbook([])

    def test_duplicate_sheet_names_are_refused_rather_than_written(self) -> None:
        # Excel opens such a workbook as corrupt, which would surface on a colleague's machine
        # rather than here.
        with pytest.raises(ValueError, match="unique"):
            write_workbook([Sheet("Same", ("a",), ()), Sheet("Same", ("b",), ())])

    def test_an_over_long_name_is_cut_to_what_excel_accepts(self) -> None:
        sheet = Sheet("C" * 50, ("a",), ())
        assert len(sheet.safe_name()) == MAX_SHEET_NAME

    def test_characters_excel_refuses_in_a_name_are_replaced(self) -> None:
        assert "/" not in Sheet("Run 1/2", ("a",), ()).safe_name()
        assert "[" not in Sheet("Run [1]", ("a",), ()).safe_name()


class TestColumnLetters:
    @pytest.mark.parametrize(
        ("index", "letter"), [(1, "A"), (26, "Z"), (27, "AA"), (52, "AZ"), (53, "BA"),
                              (702, "ZZ"), (703, "AAA")]
    )
    def test_the_letters_match_the_spreadsheet_convention(self, index: int, letter: str) -> None:
        # Hand-checked boundaries: the sequence is bijective base-26, not plain base-26, so
        # 26 -> Z and 27 -> AA rather than anything involving a zero digit.
        assert column_letter(index) == letter

    def test_a_zero_or_negative_index_is_refused(self) -> None:
        with pytest.raises(ValueError, match="1-based"):
            column_letter(0)

    def test_a_wide_sheet_still_addresses_its_last_column(self) -> None:
        headers = tuple(f"c{i}" for i in range(1, 60))
        data = write_workbook([Sheet("Wide", headers, (tuple(range(1, 60)),))])
        assert len(sheets_of(data)["Wide"][0]) == 59
