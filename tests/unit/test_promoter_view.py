"""Tests for the promoter view model.

The view's job is to arrange the engine's conclusions in promoter coordinates without deciding
anything, so the tests check three properties: that genomic features land in the right place on
both strands, that nothing is silently dropped, and that every annotation keeps the evidence kind
it arrived with.
"""

from typing import ClassVar

import pytest

from sarna.design.config import DesignConfig
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import Gene, Transcript
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.regulatory.chromatin import (
    GSE108506,
    BackgroundModel,
    ChromatinProfile,
    ChromatinSignal,
    HistoneMark,
)
from sarna.regulatory.evidence import AF22, HELA, CellTypeStatus
from sarna.regulatory.features import (
    RegulatoryAnnotation,
    RegulatoryFeature,
    RegulatoryFeatureKind,
)
from sarna.regulatory.motifs import load_motifs, scan
from sarna.sequence.nucleic import NucleicSequence
from sarna.view.promoter_view import (
    PromoterView,
    alternative_tss_spans,
    cell_panels_for,
)

MOTIFS = {motif.name: motif for motif in load_motifs()}

# A CpG-depleted promoter, so the island detector does not call the whole thing one island.
# See docs/scientific/synthetic-promoter-fixtures.md.
PROMOTER_BASES = ("GCTACGGTATCATGTCATAAGTTCAATTGACATTGGAATTACCATAGTTAACGTTTAGCATA" * 10)[:600]


def promoter_on(strand: Strand) -> Promoter:
    return Promoter(
        sequence=NucleicSequence.parse(PROMOTER_BASES),
        tss_start=TssPosition(-500),
        interval=GenomicInterval(Assembly.GRCH38, "2", 100_000, 100_600, strand),
    )


def transcript_at(ensembl_id: str, tss: int, strand: Strand) -> Transcript:
    low, high = (tss, tss + 1000) if strand is Strand.PLUS else (tss - 1000, tss + 1)
    return Transcript(
        ensembl_id=ensembl_id,
        biotype="protein_coding",
        interval=GenomicInterval(Assembly.GRCH38, "2", low, high, strand),
        is_mane_select=False,
        is_canonical=False,
    )


def view_on(
    strand: Strand, gene: Gene | None = None, transcript: Transcript | None = None
) -> PromoterView:
    promoter = promoter_on(strand)
    result = design_sarnas(promoter, WANG_2015.resolve())
    return PromoterView(result=result, gene=gene, transcript=transcript)


class TestPlacement:
    @pytest.mark.strand
    def test_a_regulatory_feature_lands_at_the_same_place_on_both_strands(self) -> None:
        # The genomic span is the same; its promoter offsets are mirrored. A feature at the low
        # end of the contig is near the promoter's 3' end on the minus strand.
        feature = RegulatoryFeature(
            feature_id="ENSR1",
            kind=RegulatoryFeatureKind.PROMOTER,
            interval=GenomicInterval(Assembly.GRCH38, "2", 100_000, 100_100),
        )
        for strand, expected in ((Strand.PLUS, (0, 100)), (Strand.MINUS, (500, 600))):
            promoter = promoter_on(strand)
            annotation = RegulatoryAnnotation(promoter.interval, (feature,))  # type: ignore[arg-type]
            view = PromoterView(
                result=design_sarnas(promoter, WANG_2015.resolve()), regulatory=annotation
            )
            placed = view.regulatory_spans()
            assert (placed[0].start, placed[0].end) == expected

    def test_a_feature_extending_past_the_edge_is_marked_clipped(self) -> None:
        # So a renderer draws it as reaching the boundary rather than as a complete feature.
        feature = RegulatoryFeature(
            feature_id="ENSR2",
            kind=RegulatoryFeatureKind.ENHANCER,
            interval=GenomicInterval(Assembly.GRCH38, "2", 99_900, 100_050),
        )
        promoter = promoter_on(Strand.PLUS)
        annotation = RegulatoryAnnotation(promoter.interval, (feature,))  # type: ignore[arg-type]
        view = PromoterView(
            result=design_sarnas(promoter, WANG_2015.resolve()), regulatory=annotation
        )
        placed = view.regulatory_spans()[0]
        assert placed.clipped
        assert (placed.start, placed.end) == (0, 50)

    def test_a_feature_entirely_outside_the_promoter_is_omitted(self) -> None:
        feature = RegulatoryFeature(
            feature_id="ENSR3",
            kind=RegulatoryFeatureKind.ENHANCER,
            interval=GenomicInterval(Assembly.GRCH38, "2", 200_000, 200_100),
        )
        promoter = promoter_on(Strand.PLUS)
        annotation = RegulatoryAnnotation(promoter.interval, (feature,))  # type: ignore[arg-type]
        view = PromoterView(
            result=design_sarnas(promoter, WANG_2015.resolve()), regulatory=annotation
        )
        assert view.regulatory_spans() == ()

    def test_no_regulatory_annotation_means_no_spans(self) -> None:
        assert view_on(Strand.PLUS).regulatory_spans() == ()

    def test_the_tss_offset_follows_from_tss_start(self) -> None:
        # tss_start is -500 on a 600 bp promoter, so +1 sits at offset 500.
        for strand in (Strand.PLUS, Strand.MINUS):
            assert view_on(strand).tss_offset == 500

    def test_a_promoter_without_anchors_has_no_tss_offset(self) -> None:
        bare = Promoter(sequence=NucleicSequence.parse(PROMOTER_BASES))
        view = PromoterView(result=design_sarnas(bare, WANG_2015.resolve()))
        assert view.tss_offset is None


class TestAlternativeTranscriptionStartSites:
    """A gene has many TSSs; a diagram showing one hides what the choice cost (ADR 0012)."""

    def test_alternatives_inside_the_window_are_placed(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        others = [
            transcript_at("ENST1", 100_200, Strand.PLUS),
            transcript_at("ENST2", 100_400, Strand.PLUS),
        ]
        spans = alternative_tss_spans(promoter, others, None)
        assert [s.start for s in spans] == [200, 400]

    def test_the_selected_transcript_is_excluded(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        selected = transcript_at("ENST1", 100_200, Strand.PLUS)
        others = [selected, transcript_at("ENST2", 100_400, Strand.PLUS)]
        spans = alternative_tss_spans(promoter, others, selected)
        assert [s.start for s in spans] == [400]

    def test_transcripts_sharing_a_tss_collapse_to_one_span(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        others = [
            transcript_at("ENST1", 100_200, Strand.PLUS),
            transcript_at("ENST2", 100_200, Strand.PLUS),
        ]
        spans = alternative_tss_spans(promoter, others, None)
        assert len(spans) == 1
        assert "ENST1, ENST2" in spans[0].detail

    def test_a_tss_outside_the_window_is_omitted(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        spans = alternative_tss_spans(
            promoter, [transcript_at("ENST9", 500_000, Strand.PLUS)], None
        )
        assert spans == ()

    @pytest.mark.strand
    def test_alternatives_are_mirrored_on_the_minus_strand(self) -> None:
        promoter = promoter_on(Strand.MINUS)
        spans = alternative_tss_spans(
            promoter, [transcript_at("ENST1", 100_200, Strand.MINUS)], None
        )
        # Offset 0 is coordinate 100_599, so 100_200 is 399 bases along.
        assert [s.start for s in spans] == [399]

    def test_the_detail_warns_that_positions_are_not_relative_to_it(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        spans = alternative_tss_spans(
            promoter, [transcript_at("ENST1", 100_200, Strand.PLUS)], None
        )
        assert "not relative to this TSS" in spans[0].detail


class TestMotifTrack:
    SEQUENCE: ClassVar[NucleicSequence] = NucleicSequence.parse(
        ("GGGGCGGGGCCGCCCCTATAAAAGCTACGGTATCATGTCATAAGTT" * 6)[:270]
    )

    def view(self) -> PromoterView:
        promoter = Promoter(sequence=self.SEQUENCE)
        scans = tuple(scan(self.SEQUENCE, m, shuffles=5) for m in MOTIFS.values())
        return PromoterView(
            result=design_sarnas(promoter, WANG_2015.resolve()), motif_scans=scans
        )

    def test_motif_sites_become_placed_spans(self) -> None:
        spans = self.view().motif_spans()
        assert spans
        assert all(span.kind == "motif" for span in spans)

    def test_spans_are_ordered_by_offset(self) -> None:
        starts = [span.start for span in self.view().motif_spans()]
        assert starts == sorted(starts)

    def test_motifs_at_or_below_chance_are_left_off_the_drawing(self) -> None:
        # A GC-rich promoter matches most promoter-factor motifs by composition alone, so drawing
        # every hit would bury the informative ones.
        all_sites = sum(len(s.sites) for s in self.view().motif_scans)
        drawn = len(self.view().motif_spans(minimum_enrichment=1.0))
        assert drawn <= all_sites

    def test_a_high_enrichment_threshold_draws_less(self) -> None:
        view = self.view()
        assert len(view.motif_spans(minimum_enrichment=5.0)) <= len(
            view.motif_spans(minimum_enrichment=1.0)
        )

    def test_every_span_carries_the_prediction_caveat(self) -> None:
        for span in self.view().motif_spans():
            assert "not 'the factor is there'" in span.detail

    def test_the_detail_names_the_matrix_and_its_experiment_type(self) -> None:
        span = self.view().motif_spans()[0]
        assert "MA0" in span.detail

    def test_scans_that_found_nothing_are_still_reported(self) -> None:
        # So a reader can see what was looked for, not only what was found.
        exported = self.view().to_dict()
        assert len(exported["motif_scans"]) == len(MOTIFS)


class _StubChromatinReader:
    """Stands in for the AF22 reader so the view can be tested without 890 MB of bigWigs.

    Declares ``cell_type`` because the real reader does: the view asks the reader whose tracks
    it holds rather than testing identity against AF22.
    """

    marks = (HistoneMark.H3K4ME3, HistoneMark.H3K27ME3)
    cell_type = AF22
    provenance = GSE108506

    def profile(self, interval: GenomicInterval) -> ChromatinProfile:
        """Return a profile with one enriched activating mark."""
        background = BackgroundModel(
            contig=interval.contig,
            mark=HistoneMark.H3K4ME3,
            window_length=interval.length,
            median=1.0,
            percentile_95=3.0,
            windows=250,
            seed=11,
        )
        signal = ChromatinSignal(
            mark=HistoneMark.H3K4ME3,
            cell_type=AF22,
            provenance=GSE108506,
            mean=50.0,
            maximum=90.0,
            background=background,
            source_interval=interval,
        )
        return ChromatinProfile(
            interval=interval, cell_type=AF22, signals=(signal,), provenance=GSE108506
        )

    def signal_track(
        self, interval: GenomicInterval, mark: HistoneMark, bins: int = 200
    ) -> tuple[float, ...] | None:
        """Return a flat track for the marks this stub has."""
        if mark not in self.marks:
            return None
        return tuple(float(i) for i in range(bins))


class TestCellTypePanels:
    def test_only_the_lines_this_project_experimented_in_get_a_panel(self) -> None:
        # The registry also carries lines declared for lookup -- HEK293, HepG2, K562 and the
        # rest -- in which nothing here was measured. A panel for each would be six identical
        # statements of ignorance beside the one line that has chromatin data.
        panels = cell_panels_for()
        assert {panel.cell_type.id for panel in panels} == {"hela", "af22"}
        assert all(panel.cell_type.is_experimental_context for panel in panels)

    def test_a_line_declared_only_for_lookup_can_still_be_asked_for_explicitly(self) -> None:
        from sarna.regulatory.evidence import cell_type

        panels = cell_panels_for([cell_type("sh-sy5y")])
        assert [panel.cell_type.id for panel in panels] == ["sh-sy5y"]
        assert panels[0].status is CellTypeStatus.NOT_ASSAYED

    def test_every_panel_is_honestly_marked_not_assayed(self) -> None:
        # No factor-occupancy data is wired in. Saying so is what stops an empty panel from
        # reading as "nothing binds here".
        for panel in cell_panels_for():
            assert panel.status is CellTypeStatus.NOT_ASSAYED
            assert not panel.status.is_known
            assert panel.absence_note is not None

    def test_the_hela_note_points_at_the_data_that_exists(self) -> None:
        hela = next(p for p in cell_panels_for() if p.cell_type is HELA)
        assert hela.absence_note is not None
        assert "ENCODE" in hela.absence_note

    def test_the_af22_note_explains_why_there_is_no_factor_data(self) -> None:
        af22 = next(p for p in cell_panels_for() if p.cell_type is AF22)
        assert af22.absence_note is not None
        assert "No public ChIP-seq" in af22.absence_note

    def test_a_panel_with_chromatin_still_says_it_is_not_occupancy(self) -> None:
        # The gap that survives loading histone tracks: accessibility is not occupancy, and
        # nothing in this project says which factor is bound in AF22.
        promoter = promoter_on(Strand.PLUS)
        reader = _StubChromatinReader()
        panels = cell_panels_for(
            [AF22], promoter=promoter, chromatin_reader=reader, bins=8  # type: ignore[arg-type]
        )
        panel = panels[0]
        assert panel.has_chromatin
        assert panel.status is not CellTypeStatus.NOT_ASSAYED
        assert panel.absence_note is not None
        assert "factor occupancy is not" in panel.absence_note

    def test_chromatin_signal_tracks_reach_the_panel_for_drawing(self) -> None:
        panels = cell_panels_for(
            [AF22],
            promoter=promoter_on(Strand.PLUS),
            chromatin_reader=_StubChromatinReader(),  # type: ignore[arg-type]
            bins=8,
        )
        tracks = panels[0].signal_tracks
        assert "H3K4me3" in tracks
        assert len(tracks["H3K4me3"]) == 8

    def test_hela_gets_no_chromatin_because_no_hela_tracks_exist(self) -> None:
        panels = cell_panels_for(
            [HELA],
            promoter=promoter_on(Strand.PLUS),
            chromatin_reader=_StubChromatinReader(),  # type: ignore[arg-type]
        )
        assert panels[0].chromatin is None
        assert panels[0].status is CellTypeStatus.NOT_ASSAYED

    def test_a_specific_selection_can_be_requested(self) -> None:
        assert len(cell_panels_for([HELA])) == 1


class TestRejectionBreakdown:
    def test_it_counts_rejections_per_rule_most_first(self) -> None:
        breakdown = view_on(Strand.PLUS).rejection_breakdown()
        counts = [entry["rejected"] for entry in breakdown]
        assert counts == sorted(counts, reverse=True)

    def test_it_reports_a_percentage_of_everything_generated(self) -> None:
        view = view_on(Strand.PLUS)
        for entry in view.rejection_breakdown():
            assert 0.0 <= entry["percent_of_generated"] <= 100.0

    def test_the_position_19_rule_is_among_the_largest(self) -> None:
        # The catalogue claims it cuts roughly three quarters of candidates on its own, written
        # from the paper rather than measured. This is the measurement.
        breakdown = {
            entry["rule"]: entry["percent_of_generated"]
            for entry in view_on(Strand.PLUS).rejection_breakdown()
        }
        assert breakdown["pos19-a"] > 50.0

    def test_a_run_with_no_rejections_has_an_empty_breakdown(self) -> None:
        promoter = Promoter(sequence=NucleicSequence.parse(PROMOTER_BASES))
        # Only preferences selected, so nothing can be excluded.
        config = DesignConfig.of(["target-length-19", "start-1-gc"])
        view = PromoterView(result=design_sarnas(promoter, config))
        assert view.rejection_breakdown() == ()


class TestExport:
    def test_every_generated_candidate_is_present(self) -> None:
        view = view_on(Strand.PLUS)
        exported = view.to_dict()
        assert len(exported["candidates"]) == len(view.result.generated)

    def test_rejected_candidates_can_be_abbreviated_without_being_dropped(self) -> None:
        view = view_on(Strand.PLUS)
        compact = view.to_dict(full_rules_for_rejected=False)
        assert len(compact["candidates"]) == len(view.result.generated)
        abbreviated = [c for c in compact["candidates"] if c.get("abbreviated")]
        assert abbreviated
        assert all(c["excluded_by"] for c in abbreviated)

    def test_an_abbreviated_candidate_keeps_its_reason(self) -> None:
        compact = view_on(Strand.PLUS).to_dict(full_rules_for_rejected=False)
        abbreviated = next(c for c in compact["candidates"] if c.get("abbreviated"))
        assert abbreviated["rejection_summary"]
        assert "rules" not in abbreviated

    def test_passing_candidates_keep_their_full_rule_results(self) -> None:
        compact = view_on(Strand.PLUS).to_dict(full_rules_for_rejected=False)
        passing = [c for c in compact["candidates"] if not c.get("abbreviated")]
        assert passing
        assert all("rules" in c for c in passing)

    def test_abbreviating_makes_the_document_much_smaller(self) -> None:
        import json

        view = view_on(Strand.PLUS)
        full = len(json.dumps(view.to_dict(full_rules_for_rejected=True)))
        compact = len(json.dumps(view.to_dict(full_rules_for_rejected=False)))
        assert compact < full / 2

    def test_the_configuration_that_produced_the_run_travels_with_it(self) -> None:
        exported = view_on(Strand.PLUS).to_dict()
        assert exported["configuration"]["preset"]["id"] == "wang-2015"

    def test_the_promoter_metadata_names_the_assembly_and_strand(self) -> None:
        exported = view_on(Strand.MINUS).to_dict()
        assert exported["promoter"]["assembly"] == "GRCh38"
        assert exported["promoter"]["strand"] == "-"

    def test_unevaluated_rules_are_reported_once_at_the_top(self) -> None:
        bare = Promoter(sequence=NucleicSequence.parse(PROMOTER_BASES))
        view = PromoterView(result=design_sarnas(bare, WANG_2015.resolve()))
        assert "tss-window" in view.to_dict()["unevaluated_rules"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps(view_on(Strand.MINUS).to_dict())

    def test_gene_and_transcript_are_optional(self) -> None:
        exported = view_on(Strand.PLUS).to_dict()
        assert exported["gene"] is None
        assert exported["transcript"] is None

    def test_a_named_gene_and_transcript_are_carried_through(self) -> None:
        gene = Gene(
            ensembl_id="ENSG1",
            symbol="TEST1",
            interval=GenomicInterval(Assembly.GRCH38, "2", 100_000, 101_000),
            description="a test gene",
        )
        transcript = transcript_at("ENST1", 100_500, Strand.PLUS)
        exported = view_on(Strand.PLUS, gene=gene, transcript=transcript).to_dict()
        assert exported["gene"]["symbol"] == "TEST1"
        assert exported["transcript"]["ensembl_id"] == "ENST1"

    def test_a_placed_span_serialises_its_own_fields(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        spans = alternative_tss_spans(
            promoter, [transcript_at("ENST1", 100_200, Strand.PLUS)], None
        )
        as_dict = spans[0].to_dict()
        assert as_dict["start"] == 200
        assert as_dict["kind"] == "alternative-tss"
        assert as_dict["clipped"] is False


class TestCpGIslandTrack:
    def test_islands_are_already_in_promoter_coordinates(self) -> None:
        # Detected on the promoter sequence itself, so no mapping is needed or applied.
        gc_rich = NucleicSequence.parse("GCGCGCGCGCGCATGCGC" * 40)
        promoter = Promoter(sequence=gc_rich)
        view = PromoterView(result=design_sarnas(promoter, WANG_2015.resolve()))
        for span in view.island_spans():
            assert 0 <= span.start < span.end <= promoter.length

    def test_a_run_without_island_detection_has_no_island_track(self) -> None:
        from sarna.design.presets import LEGACY_SCANNER

        promoter = promoter_on(Strand.PLUS)
        view = PromoterView(result=design_sarnas(promoter, LEGACY_SCANNER.resolve()))
        assert view.result.cpg_islands is None
        assert view.island_spans() == ()

    def test_an_island_span_explains_what_it_means_for_a_candidate(self) -> None:
        gc_rich = NucleicSequence.parse("GCGCGCGCGCGCATGCGC" * 40)
        promoter = Promoter(sequence=gc_rich)
        view = PromoterView(result=design_sarnas(promoter, WANG_2015.resolve()))
        spans = view.island_spans()
        if spans:
            assert "excluded by default" in spans[0].detail


class TestCellPanelExport:
    def test_a_panel_serialises_everything_a_renderer_needs(self) -> None:
        panel = cell_panels_for([AF22])[0]
        as_dict = panel.to_dict()
        assert as_dict["id"] == "af22"
        assert as_dict["known"] is False
        assert as_dict["status"] == "not-assayed"
        assert as_dict["proliferating"] is False
        assert as_dict["absence_note"]
        assert as_dict["observations"] == []

    def test_panels_reach_the_view_export(self) -> None:
        promoter = promoter_on(Strand.PLUS)
        view = PromoterView(
            result=design_sarnas(promoter, WANG_2015.resolve()),
            cell_panels=cell_panels_for(),
        )
        exported = view.to_dict()
        assert {panel["id"] for panel in exported["cell_types"]} == {"hela", "af22"}
        assert all(panel["rrid"].startswith("CVCL_") for panel in exported["cell_types"])
        assert all(panel["absence_note"] for panel in exported["cell_types"])


class TestRanksReachTheExport:
    """DesignResult.generated carries no ranks; exporting it directly loses the ordering.

    The failure is silent — a renderer sorting by rank sees ``None`` everywhere and falls back to
    position, so the list looks plausible and is not ranked. Found by checking a generated page.
    """

    def test_generated_candidates_have_no_rank_of_their_own(self) -> None:
        # The premise. Ranking returns copies so it cannot alter candidates a caller already holds.
        result = design_sarnas(promoter_on(Strand.PLUS), WANG_2015.resolve())
        assert all(candidate.rank is None for candidate in result.generated)
        assert any(candidate.rank is not None for candidate in result.ranked)

    def test_the_view_merges_the_ranks_back_in(self) -> None:
        view = view_on(Strand.PLUS)
        ranked = view.ranked_candidates()
        assert len(ranked) == len(view.result.generated)
        with_rank = [candidate for candidate in ranked if candidate.rank is not None]
        assert len(with_rank) == len(view.result.passing)

    def test_position_order_is_preserved(self) -> None:
        offsets = [candidate.offset for candidate in view_on(Strand.PLUS).ranked_candidates()]
        assert offsets == sorted(offsets)

    def test_rejected_candidates_keep_no_rank(self) -> None:
        for candidate in view_on(Strand.PLUS).ranked_candidates():
            if candidate.is_rejected:
                assert candidate.rank is None

    def test_the_export_contains_a_rank_one_candidate(self) -> None:
        exported = view_on(Strand.PLUS).to_dict()
        ranks = [c["rank"] for c in exported["candidates"] if c["rank"] is not None]
        assert 1 in ranks
        assert sorted(ranks) == list(range(1, len(ranks) + 1))
