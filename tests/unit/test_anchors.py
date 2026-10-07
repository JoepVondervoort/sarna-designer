"""Tests for promoter anchors and the cell-line registry.

Both exist because the project's own earlier work could not be reproduced without them: chapter 2
designed SCN1A against promoters mapped by Dong 2014 and Long 2008, and the Dong distal region has
no Ensembl transcript within 5.5 kb.

The substance here is that an anchor cannot lie about how well its zero point is established, and
that a cell line cannot present an absence of data as an absence of binding.
"""

import re
from itertools import pairwise
from pathlib import Path

import pytest

from sarna.genome.anchors import (
    AnchorRegistry,
    PromoterAnchor,
    ReferenceKind,
    anchor_paths,
    load_anchors,
)
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.regulatory.evidence import (
    AF22,
    CELL_TYPES,
    HELA,
    CellType,
    DataAxis,
    cell_line_paths,
    cell_type,
    load_cell_types,
)
from sarna.rules.provenance import Provenance

SOURCE = Provenance.citation("Test et al. 2026", "Table 1")
REGISTRY = load_anchors()


def anchor(
    reference_point: int | None = 1500,
    kind: ReferenceKind = ReferenceKind.ANNOTATED_TSS,
    strand: Strand = Strand.PLUS,
    start: int = 1000,
    end: int = 2000,
) -> PromoterAnchor:
    return PromoterAnchor(
        id="test-anchor",
        gene="TEST1",
        interval=GenomicInterval(Assembly.GRCH38, "1", start, end, strand),
        reference_point=reference_point,
        reference_kind=kind,
        provenance=SOURCE,
        note="A test anchor.",
    )


class TestAnAnchorCannotLieAboutItsZeroPoint:
    def test_a_measurable_kind_needs_a_reference_point(self) -> None:
        with pytest.raises(ValueError, match="gives no reference point"):
            anchor(reference_point=None, kind=ReferenceKind.MAPPED_TSS)

    def test_a_reference_point_needs_a_kind_that_is_not_none(self) -> None:
        # Otherwise a caller cannot tell an annotated TSS from a guess.
        with pytest.raises(ValueError, match="state which kind"):
            anchor(reference_point=1500, kind=ReferenceKind.NONE)

    def test_a_reference_point_outside_its_own_region_is_refused(self) -> None:
        with pytest.raises(ValueError, match="outside"):
            anchor(reference_point=5000)

    def test_an_anchor_must_explain_itself(self) -> None:
        with pytest.raises(ValueError, match="indistinguishable from a guess"):
            PromoterAnchor(
                id="bare",
                gene="TEST1",
                interval=GenomicInterval(Assembly.GRCH38, "1", 1000, 2000),
                reference_point=1500,
                reference_kind=ReferenceKind.ANNOTATED_TSS,
                provenance=SOURCE,
                note="   ",
            )

    def test_only_an_annotated_tss_carries_no_caveat(self) -> None:
        assert ReferenceKind.ANNOTATED_TSS.caveat is None
        for kind in (ReferenceKind.MAPPED_TSS, ReferenceKind.INFERRED, ReferenceKind.NONE):
            assert kind.caveat

    def test_an_inferred_zero_point_says_every_offset_inherits_the_uncertainty(self) -> None:
        caveat = ReferenceKind.INFERRED.caveat
        assert caveat is not None
        assert "inferred" in caveat
        assert "uncertainty" in caveat

    def test_a_region_with_no_start_site_is_not_measurable(self) -> None:
        assert not ReferenceKind.NONE.is_measurable
        assert ReferenceKind.ANNOTATED_TSS.is_measurable


@pytest.mark.strand
class TestAnchorWindows:
    def test_a_plus_strand_window_runs_upstream_to_lower_coordinates(self) -> None:
        window = anchor(strand=Strand.PLUS).window(upstream_bp=400, downstream_bp=100)
        assert window.interval.start == 1500 - 400
        assert window.interval.end == 1500 + 100
        assert window.tss_start is not None
        assert window.tss_start.value == -400

    def test_a_minus_strand_window_runs_upstream_to_higher_coordinates(self) -> None:
        # Mirrors Promoter.genomic_interval and EnsemblClient.promoter_for.
        window = anchor(strand=Strand.MINUS).window(upstream_bp=400, downstream_bp=100)
        assert window.interval.end == 1500 + 400 + 1
        assert window.interval.start == 1500 - 100 + 1
        assert window.interval.length == 500

    def test_an_unmeasurable_anchor_yields_the_region_and_no_tss(self) -> None:
        # The honest answer, not a degraded one: TSS rules then report not-applicable rather
        # than measuring from an invented origin.
        unmeasurable = anchor(reference_point=None, kind=ReferenceKind.NONE)
        window = unmeasurable.window(region_padding_bp=250)
        assert window.tss_start is None
        assert window.interval.start == 1000 - 250
        assert window.interval.end == 2000 + 250

    def test_negative_widths_are_refused(self) -> None:
        with pytest.raises(ValueError, match="must not be negative"):
            anchor().window(upstream_bp=-1)

    def test_a_zero_length_measurable_window_is_refused(self) -> None:
        with pytest.raises(ValueError, match="zero length"):
            anchor().window(upstream_bp=0, downstream_bp=0)

    def test_a_window_running_off_the_contig_start_is_refused(self) -> None:
        near_start = anchor(reference_point=50, start=0, end=100)
        with pytest.raises(ValueError, match="before the beginning of contig"):
            near_start.window(upstream_bp=500)

    def test_an_unmeasurable_window_can_be_built_with_no_padding(self) -> None:
        unmeasurable = anchor(reference_point=None, kind=ReferenceKind.NONE)
        window = unmeasurable.window(region_padding_bp=0)
        assert window.interval.length == 1000


class TestTheShippedScn1aAnchors:
    """Chapter 2's three target regions, which are why anchors exist at all."""

    def test_all_three_scn1a_promoters_ship(self) -> None:
        assert {a.id for a in REGISTRY.for_gene("SCN1A")} == {
            "scn1a-long-2008",
            "scn1a-dong-2014-proximal",
            "scn1a-dong-2014-distal",
        }

    def test_the_gene_lookup_is_case_insensitive(self) -> None:
        assert len(REGISTRY.for_gene("scn1a")) == 3

    def test_a_gene_with_no_declared_anchor_returns_empty(self) -> None:
        # Normal: most genes use the annotated transcript path.
        assert REGISTRY.for_gene("GAPDH") == ()

    def test_the_two_annotation_backed_regions_say_so(self) -> None:
        for anchor_id in ("scn1a-long-2008", "scn1a-dong-2014-proximal"):
            assert REGISTRY[anchor_id].reference_kind is ReferenceKind.ANNOTATED_TSS

    def test_the_distal_region_declares_its_zero_point_inferred(self) -> None:
        # The one Ensembl does not annotate, and the one that produced two activators. Its
        # reference point is deduced from where the chapter's designs fell, not reported.
        distal = REGISTRY["scn1a-dong-2014-distal"]
        assert distal.reference_kind is ReferenceKind.INFERRED
        assert "INFERRED" in distal.note
        assert "Dong" in str(distal.provenance)

    def test_the_distal_note_states_the_uncertainty_and_how_to_fix_it(self) -> None:
        note = REGISTRY["scn1a-dong-2014-distal"].note
        assert "410 bp" in note
        assert 'reference_kind = "mapped-tss"' in note

    def test_every_anchor_cites_a_publication(self) -> None:
        for shipped in REGISTRY:
            assert shipped.provenance.source
            assert shipped.provenance.locator

    def test_all_three_are_on_the_minus_strand(self) -> None:
        assert all(a.strand is Strand.MINUS for a in REGISTRY.for_gene("SCN1A"))

    def test_the_regions_do_not_overlap_each_other(self) -> None:
        regions = sorted(
            (a.interval.start, a.interval.end) for a in REGISTRY.for_gene("SCN1A")
        )
        for (_, first_end), (second_start, _) in pairwise(regions):
            assert first_end <= second_start

    def test_the_export_carries_the_caveat(self) -> None:
        exported = REGISTRY["scn1a-dong-2014-distal"].to_dict()
        assert exported["reference_kind"] == "inferred"
        assert exported["caveat"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps([a.to_dict() for a in REGISTRY])


class TestTheRegistry:
    def test_a_duplicate_id_is_refused(self) -> None:
        with pytest.raises(ValueError, match="duplicate anchor id"):
            AnchorRegistry([anchor(), anchor()])

    def test_an_unknown_id_lists_what_exists(self) -> None:
        with pytest.raises(KeyError, match="scn1a-long-2008"):
            REGISTRY["nope"]

    def test_the_registry_reports_its_genes(self) -> None:
        assert REGISTRY.genes == ("SCN1A",)

    def test_membership_and_iteration_work(self) -> None:
        assert "scn1a-long-2008" in REGISTRY
        assert len(list(REGISTRY)) == len(REGISTRY)


class TestUserSuppliedAnchors:
    EXTRA = """
[[anchor]]
id = "mygene-smith-2020"
gene = "MYGENE"
contig = "7"
start = 1000
end = 3000
strand = "+"
reference_point = 2000
reference_kind = "mapped-tss"
source = "Smith et al. 2020"
locator = "J Test 1:1-2"
note = "A user-supplied promoter."
"""

    def test_a_user_file_adds_anchors(self, tmp_path: Path) -> None:
        # The whole point: a colleague adds a promoter without touching the code.
        path = tmp_path / "mine.toml"
        path.write_text(self.EXTRA)
        registry = load_anchors([path])
        assert "mygene-smith-2020" in registry
        assert registry.for_gene("MYGENE")[0].reference_kind is ReferenceKind.MAPPED_TSS

    def test_the_shipped_anchors_are_still_present(self, tmp_path: Path) -> None:
        path = tmp_path / "mine.toml"
        path.write_text(self.EXTRA)
        assert "scn1a-long-2008" in load_anchors([path])

    def test_a_directory_contributes_every_toml(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        (tmp_path / "a.toml").write_text(self.EXTRA)
        monkeypatch.setenv("SARNA_ANCHORS", str(tmp_path))
        assert anchor_paths() == (tmp_path / "a.toml",)

    def test_a_single_file_can_be_named_by_environment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        path = tmp_path / "mine.toml"
        path.write_text(self.EXTRA)
        monkeypatch.setenv("SARNA_ANCHORS", str(path))
        assert "mygene-smith-2020" in load_anchors()

    def test_a_nonexistent_path_is_ignored(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SARNA_ANCHORS", "/no/such/place.toml")
        assert anchor_paths() == ()

    def test_no_environment_variable_means_no_extra_paths(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("SARNA_ANCHORS", raising=False)
        assert anchor_paths() == ()

    def test_malformed_toml_names_the_file(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.toml"
        path.write_text("[[anchor]\nid = ")
        with pytest.raises(ValueError, match="not valid TOML"):
            load_anchors([path])

    def test_a_bad_entry_names_the_file_and_its_position(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.toml"
        path.write_text('[[anchor]]\nid = "x"\ngene = "Y"\n')
        with pytest.raises(ValueError, match="anchor 1"):
            load_anchors([path])

    def test_a_user_anchor_colliding_with_a_shipped_id_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "clash.toml"
        path.write_text(
            self.EXTRA.replace("mygene-smith-2020", "scn1a-long-2008").replace(
                'gene = "MYGENE"', 'gene = "SCN1A"'
            )
        )
        with pytest.raises(ValueError, match="duplicate anchor id"):
            load_anchors([path])


class TestTheCellLineRegistry:
    def test_the_shipped_lines_load_from_data(self) -> None:
        # The two the project has evidence from, plus six declared because the saRNA literature
        # used them or because a missing axis could realistically be answered in them. Pinned by
        # name: adding or removing a line is a decision, not a refactor.
        assert set(CELL_TYPES) == {
            "hela", "af22", "hek293", "hepg2", "k562", "mcf7", "pc3", "sh-sy5y"}
        assert cell_type("hela") is HELA

    def test_only_af22_ships_with_any_data_loaded(self) -> None:
        # Six lines were added for relevance, not for data. If one ever acquires a loaded axis
        # that must be a deliberate change, because it alters what the tool claims to know.
        loaded = {identifier for identifier, cell in CELL_TYPES.items() if cell.available_axes}
        assert loaded == {"af22"}

    def test_proliferation_is_recorded_because_it_bears_on_rnaa(self) -> None:
        assert HELA.is_proliferating
        assert not AF22.is_proliferating

    def test_every_axis_is_declared_for_every_shipped_line(self) -> None:
        # The regression test for a silent key mismatch: tf_occupancy was declared in the file,
        # looked for under a hyphenated name, and never found.
        for cell in CELL_TYPES.values():
            for axis in DataAxis:
                availability = cell.availability(axis)
                assert availability.note, f"{cell.id} declares no note for {axis}"

    def test_af22_has_chromatin_and_nothing_else(self) -> None:
        assert AF22.available_axes == (DataAxis.CHROMATIN,)
        assert AF22.availability(DataAxis.CHROMATIN).source == "GSE108506"

    def test_hela_has_nothing_loaded_but_everything_reachable(self) -> None:
        # The asymmetry worth stating: the tool can learn far more about HeLa than about the
        # disease-relevant line.
        assert HELA.available_axes == ()
        assert all(HELA.availability(axis).is_reachable for axis in DataAxis)

    def test_af22_expression_has_no_route_at_all(self) -> None:
        # An iPSC-derived research line, absent from DepMap and GTEx.
        expression = AF22.availability(DataAxis.EXPRESSION)
        assert not expression.available
        assert not expression.is_reachable

    def test_each_axis_states_the_question_it_answers(self) -> None:
        for axis in DataAxis:
            assert axis.question.endswith("?")

    def test_expression_is_named_as_the_first_question(self) -> None:
        # A saRNA enhances transcription from a locus that already has some.
        assert "transcribed" in DataAxis.EXPRESSION.question

    def test_an_undeclared_axis_defaults_to_nothing_with_no_route(self) -> None:
        bare = CellType(
            id="x",
            name="X",
            description="d",
            is_proliferating=True,
            provenance=SOURCE,
        )
        availability = bare.availability(DataAxis.EXPRESSION)
        assert not availability.available
        assert not availability.is_reachable
        assert "X" in availability.note

    def test_an_unknown_line_explains_how_to_add_one(self) -> None:
        # Deliberately a line this project does not ship: HEK293 used to stand here and now
        # exists, which is why the assertion names something nobody is about to add.
        assert "jurkat" not in CELL_TYPES
        with pytest.raises(KeyError, match="SARNA_CELL_LINES"):
            cell_type("jurkat")

    def test_a_user_file_adds_cell_lines(self, tmp_path: Path) -> None:
        path = tmp_path / "mine.toml"
        path.write_text(
            """
[[cell_line]]
id = "my-ipsc-neurons"
name = "Lab iPSC-derived neurons"
proliferating = false
source = "This laboratory, unpublished"
locator = "internal differentiation protocol, revision 3"
description = "Patient-derived neurons used as the disease-relevant context."

[cell_line.expression]
available = true
source = "DepMap"
note = "RNA-seq available."
"""
        )
        registry = load_cell_types([path])
        assert "my-ipsc-neurons" in registry
        assert registry["my-ipsc-neurons"].available_axes == (DataAxis.EXPRESSION,)
        # A researcher's own line need not carry an accession; the shipped ones must.
        assert registry["my-ipsc-neurons"].rrid == ""

    def test_a_duplicate_cell_line_id_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "clash.toml"
        path.write_text(
            """
[[cell_line]]
id = "hela"
name = "HeLa again"
proliferating = true
source = "x"
locator = "y"
description = "d"
"""
        )
        with pytest.raises(ValueError, match="duplicate cell line id"):
            load_cell_types([path])

    @staticmethod
    def shipped(monkeypatch: pytest.MonkeyPatch) -> dict[str, CellType]:
        """Return only the lines this project ships.

        ``CELL_TYPES`` also carries whatever ``SARNA_CELL_LINES`` points at, and a researcher's
        own line is not held to the obligations below -- they may well be working with cells
        that have no accession yet.
        """
        monkeypatch.delenv("SARNA_CELL_LINES", raising=False)
        return load_cell_types()

    def test_every_shipped_line_is_identified_by_its_cellosaurus_accession(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A cell line name is ambiguous; an RRID is not.

        "HEK293" names several lines that differ, and HeLa has many drifted stocks. A result
        that records only a display name cannot be checked against anybody else's. Each shipped
        line therefore carries its Cellosaurus accession, which resolves to one record with the
        line's origin, its known problems and its cross-references.
        """
        for identifier, cell in self.shipped(monkeypatch).items():
            assert cell.rrid, f"{identifier} ships without a Cellosaurus accession"
            assert re.fullmatch(r"CVCL_[0-9A-Z]+", cell.rrid), (
                f"{identifier} carries {cell.rrid!r}, which is not a Cellosaurus accession")

    def test_every_shipped_line_cites_a_source_and_a_locator(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The same obligation every rule carries: a claim arrives with where it came from.
        for identifier, cell in self.shipped(monkeypatch).items():
            assert cell.provenance.locator, f"{identifier} cites a source with no locator"
            assert len(cell.description) > 80, f"{identifier} has no usable description"

    def test_every_axis_of_every_line_says_where_its_answer_would_come_from(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # "available = false" is a first-class answer, but only when it says what is missing.
        for identifier, cell in self.shipped(monkeypatch).items():
            for axis in DataAxis:
                note = cell.availability(axis).note
                assert note, f"{identifier}/{axis} declares neither data nor a reason"

    def test_the_accession_travels_into_the_export(self) -> None:
        exported = HELA.to_dict()
        assert exported["rrid"] == "CVCL_0030"

    def test_the_export_lists_every_axis(self) -> None:
        exported = AF22.to_dict()
        axes = exported["data"]
        assert isinstance(axes, list)
        assert len(axes) == len(DataAxis)
        assert exported["provenance"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps([cell.to_dict() for cell in CELL_TYPES.values()])


class TestTheTranscriptPathIsAnAnchorToo:
    """``promoter_for`` is a wrapper over ``promoter_for_anchor``, so the arithmetic lives once."""

    def test_a_transcript_becomes_an_annotated_tss_anchor(self) -> None:
        built = PromoterAnchor.from_transcript(
            gene="SCN1A",
            transcript_id="ENST00000674923",
            tss=166_128_019,
            interval=GenomicInterval(
                Assembly.GRCH38, "2", 165_985_813, 166_128_020, Strand.MINUS
            ),
        )
        assert built.reference_kind is ReferenceKind.ANNOTATED_TSS
        assert built.reference_point == 166_128_019
        assert built.reference_kind.caveat is None
        assert "ENST00000674923" in str(built.provenance)

    def test_its_window_matches_what_the_client_produced_before(self) -> None:
        # The behaviour the unification had to preserve: minus strand, 2000 up, 500 down.
        built = PromoterAnchor.from_transcript(
            gene="SCN1A",
            transcript_id="ENST00000674923",
            tss=166_128_019,
            interval=GenomicInterval(
                Assembly.GRCH38, "2", 165_985_813, 166_128_020, Strand.MINUS
            ),
        )
        window = built.window(upstream_bp=2000, downstream_bp=500)
        assert window.interval.start == 166_127_520
        assert window.interval.end == 166_130_020
        assert window.tss_start is not None
        assert window.tss_start.value == -2000


class TestCellLineFileErrors:
    def test_malformed_toml_names_the_file(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.toml"
        path.write_text("[[cell_line]\nid = ")
        with pytest.raises(ValueError, match="not valid TOML"):
            load_cell_types([path])

    def test_a_bad_entry_names_its_position(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.toml"
        path.write_text('[[cell_line]]\nid = "x"\n')
        with pytest.raises(ValueError, match="cell_line 1"):
            load_cell_types([path])

    def test_a_directory_of_files_can_be_named_by_environment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        (tmp_path / "extra.toml").write_text(
            '[[cell_line]]\nid = "x1"\nname = "X1"\nproliferating = true\n'
            'source = "s"\nlocator = "l"\ndescription = "d"\n'
        )
        monkeypatch.setenv("SARNA_CELL_LINES", str(tmp_path))
        assert "x1" in load_cell_types()

    def test_a_single_file_can_be_named_by_environment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        path = tmp_path / "one.toml"
        path.write_text(
            '[[cell_line]]\nid = "x2"\nname = "X2"\nproliferating = false\n'
            'source = "s"\nlocator = "l"\ndescription = "d"\n'
        )
        monkeypatch.setenv("SARNA_CELL_LINES", str(path))
        assert "x2" in load_cell_types()

    def test_a_nonexistent_path_is_ignored(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SARNA_CELL_LINES", "/no/such/file.toml")
        assert cell_line_paths() == ()
