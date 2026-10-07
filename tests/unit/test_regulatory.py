"""Tests for the regulatory annotation layer.

The substance here is the evidence model, not the arithmetic. A tool that draws a transcription
factor on a promoter is making a claim, and the tests enforce that it cannot make a stronger one
than it holds: a motif match may not name a cell type, an absent annotation must carry its caveat,
and "not assayed" must stay distinguishable from "not bound".
"""

import re
from typing import Any, ClassVar, cast

import pytest

from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import EnsemblClient
from sarna.regulatory.evidence import (
    AF22,
    HELA,
    CellTypeStatus,
    Evidence,
    EvidenceKind,
    aggregated_annotation,
    cell_type,
    motif_prediction,
    observed_in,
)
from sarna.regulatory.features import (
    RegulatoryAnnotation,
    RegulatoryFeatureKind,
    features_by_kind,
    fetch_regulatory_features,
    promoter_support,
)
from sarna.regulatory.motifs import (
    DEFAULT_RELATIVE_THRESHOLD,
    MINIMUM_RAW_SCORE,
    MotifDataType,
    MotifScan,
    _sites_in,
    load_motifs,
    scan,
    scan_all,
)
from sarna.rules.provenance import Provenance
from sarna.sequence.nucleic import NucleicSequence

SOURCE = Provenance.citation("Test et al. 2026", "Table 1")
MOTIFS = load_motifs()
BY_NAME = {motif.name: motif for motif in MOTIFS}


class TestEvidenceCannotOverclaim:
    def test_a_motif_prediction_may_not_name_a_cell_type(self) -> None:
        # The central rule. A sequence match is computed without reference to any cell, so
        # attaching one would imply a measurement nobody made.
        with pytest.raises(ValueError, match="never made"):
            Evidence(
                kind=EvidenceKind.MOTIF_PREDICTION,
                provenance=SOURCE,
                cell_type=HELA,
            )

    def test_an_aggregated_annotation_may_not_name_a_cell_type(self) -> None:
        with pytest.raises(ValueError, match="never made"):
            Evidence(
                kind=EvidenceKind.AGGREGATED_ANNOTATION,
                provenance=SOURCE,
                cell_type=AF22,
            )

    def test_a_cell_type_observation_must_name_its_cell_type(self) -> None:
        with pytest.raises(ValueError, match="must name the cell type"):
            Evidence(kind=EvidenceKind.CELL_TYPE_OBSERVATION, provenance=SOURCE)

    def test_a_cell_type_observation_must_say_what_was_found(self) -> None:
        # "Observed, outcome unknown" is a contradiction.
        with pytest.raises(ValueError, match="contradicts having observed"):
            Evidence(
                kind=EvidenceKind.CELL_TYPE_OBSERVATION,
                provenance=SOURCE,
                cell_type=HELA,
                status=CellTypeStatus.NOT_ASSAYED,
            )

    def test_not_assayed_is_distinguishable_from_not_bound(self) -> None:
        # Rendering these the same way converts ignorance into a negative finding.
        assert CellTypeStatus.NOT_ASSAYED.is_known is False
        assert CellTypeStatus.ASSAYED_NOT_BOUND.is_known is True
        assert CellTypeStatus.OBSERVED_BOUND.is_known is True


class TestEvidenceConstructors:
    def test_a_motif_prediction_is_not_cell_type_specific(self) -> None:
        evidence = motif_prediction(SOURCE)
        assert not evidence.is_cell_type_specific
        assert evidence.cell_type is None
        assert "no cell type implied" in evidence.summary

    def test_an_aggregated_annotation_says_it_is_pooled(self) -> None:
        assert "pooled" in aggregated_annotation(SOURCE).summary

    def test_an_observation_names_the_cell_type_and_the_outcome(self) -> None:
        evidence = observed_in(HELA, SOURCE)
        assert evidence.is_cell_type_specific
        assert "observed bound in HeLa" in evidence.summary

    def test_a_negative_observation_reads_as_one(self) -> None:
        evidence = observed_in(AF22, SOURCE, CellTypeStatus.ASSAYED_NOT_BOUND)
        assert "assayed and not found bound" in evidence.summary

    def test_every_evidence_kind_carries_the_caveat_a_ui_must_show(self) -> None:
        for kind in EvidenceKind:
            assert kind.confidence_note.strip()

    def test_only_a_cell_type_observation_is_cell_type_specific(self) -> None:
        specific = [k for k in EvidenceKind if k.is_cell_type_specific]
        assert specific == [EvidenceKind.CELL_TYPE_OBSERVATION]


class TestCellTypes:
    def test_the_two_cell_types_from_the_screen_are_characterised(self) -> None:
        assert cell_type("hela") is HELA
        assert cell_type("af22") is AF22

    def test_proliferation_is_recorded_because_it_bears_on_rnaa(self) -> None:
        # Chapter 4's largest effect. HeLa divides continuously; AF22 does not.
        assert HELA.is_proliferating
        assert not AF22.is_proliferating

    def test_af22s_description_records_that_nothing_activated_there(self) -> None:
        assert "no saRNA activated any" in AF22.description

    def test_each_cell_type_cites_a_source(self) -> None:
        for cell in (HELA, AF22):
            assert str(cell.provenance)

    def test_an_unknown_cell_type_raises_and_lists_the_known_ones(self) -> None:
        with pytest.raises(KeyError, match="af22, hek293") as raised:
            cell_type("jurkat")
        # The point of the message is that it names what *is* available, so the researcher can
        # see whether their line is spelled differently rather than absent.
        assert "sh-sy5y" in str(raised.value)


class TestMotifMatrices:
    def test_sixteen_curated_matrices_ship(self) -> None:
        assert len(MOTIFS) == 16

    def test_every_matrix_states_why_it_is_in_the_set(self) -> None:
        # A curated subset must justify each member, or it is just an arbitrary subset.
        for motif in MOTIFS:
            assert len(motif.relevance) > 40

    def test_the_founding_saRNA_experiments_factor_is_present(self) -> None:
        # Li et al. 2006 name SP1 sites in the p21 promoter beside the dsP21-322 target.
        assert "SP1" in BY_NAME
        assert "p21" in BY_NAME["SP1"].relevance

    def test_the_neuronal_repressor_is_present_with_its_reason(self) -> None:
        # REST represses SCN1A in non-neuronal cells; the best single candidate explanation for
        # a HeLa-versus-AF22 difference at this locus.
        assert "SCN1A" in BY_NAME["REST"].relevance

    def test_each_matrix_carries_its_accession_and_experiment_type(self) -> None:
        for motif in MOTIFS:
            assert motif.matrix_id.startswith("MA")
            assert isinstance(motif.data_type, MotifDataType)

    def test_provenance_names_the_accession_a_reader_can_look_up(self) -> None:
        provenance = BY_NAME["SP1"].provenance
        assert "MA0079.5" in str(provenance)

    def test_chip_seq_matrices_are_marked_as_in_vivo(self) -> None:
        assert BY_NAME["SP1"].data_type.is_in_vivo
        assert not BY_NAME["E2F1"].data_type.is_in_vivo  # HT-SELEX

    def test_an_unrecognised_experiment_type_becomes_other(self) -> None:
        assert MotifDataType.parse(None) is MotifDataType.OTHER
        assert MotifDataType.parse("Protein binding microarray") is MotifDataType.OTHER

    def test_the_score_range_is_consistent(self) -> None:
        for motif in MOTIFS:
            assert motif.minimum_score < 0 < motif.maximum_score


class TestMotifScoring:
    def test_the_consensus_sequence_scores_at_the_maximum(self) -> None:
        motif = BY_NAME["SP1"]
        consensus = "".join(
            max("ACGT", key=lambda base: motif.counts[base][position])
            for position in range(motif.length)
        )
        assert motif.score(consensus) == pytest.approx(motif.maximum_score)
        assert motif.relative_score(consensus) == pytest.approx(1.0)

    def test_the_worst_sequence_scores_at_the_minimum(self) -> None:
        motif = BY_NAME["SP1"]
        worst = "".join(
            min("ACGT", key=lambda base: motif.counts[base][position])
            for position in range(motif.length)
        )
        assert motif.relative_score(worst) == pytest.approx(0.0)

    def test_a_window_of_the_wrong_length_is_refused(self) -> None:
        with pytest.raises(ValueError, match="but the window is"):
            BY_NAME["SP1"].score("ACGT")

    def test_an_ambiguity_code_is_refused_rather_than_guessed(self) -> None:
        # Scoring N as if it were A would invent evidence.
        motif = BY_NAME["SP1"]
        with pytest.raises(ValueError, match="invent evidence"):
            motif.score("N" * motif.length)


class TestMotifScanning:
    GC_RICH: ClassVar[str] = "GGGGCGGGGCCGCCCCGCCCCGGGGCGGGGTATAAAAGGCGGGGCCGCCCC" * 4

    def sequence(self) -> NucleicSequence:
        return NucleicSequence.parse(self.GC_RICH)

    def test_a_scan_reports_both_strands(self) -> None:
        result = scan(self.sequence(), BY_NAME["SP1"], shuffles=0)
        strands = {site.strand for site in result.sites}
        assert strands <= {Strand.PLUS, Strand.MINUS}
        assert result.sites

    def test_every_site_is_a_prediction_with_no_cell_type(self) -> None:
        result = scan(self.sequence(), BY_NAME["SP1"], shuffles=0)
        for site in result.sites:
            assert site.evidence.kind is EvidenceKind.MOTIF_PREDICTION
            assert site.evidence.cell_type is None

    def test_no_site_has_a_negative_log_odds_score(self) -> None:
        # The raw-score floor. Four of the sixteen matrices admit negative-log-odds matches at a
        # relative score of 0.80, which would mean the match is less likely than background.
        for motif in MOTIFS:
            for site in scan(self.sequence(), motif, shuffles=0).sites:
                assert site.raw_score > MINIMUM_RAW_SCORE

    def test_the_floor_actually_removes_matches_that_pass_the_relative_threshold(self) -> None:
        # NRF1 is the worst case: at 0.80 relative its raw score is about -10.
        nrf1 = BY_NAME["NRF1"]
        threshold_raw = nrf1.minimum_score + DEFAULT_RELATIVE_THRESHOLD * (
            nrf1.maximum_score - nrf1.minimum_score
        )
        assert threshold_raw < 0
        bases = self.GC_RICH
        passing_relative = sum(
            1
            for offset in range(len(bases) - nrf1.length + 1)
            if nrf1.relative_score(bases[offset : offset + nrf1.length])
            >= DEFAULT_RELATIVE_THRESHOLD
        )
        kept = len(scan(self.sequence(), nrf1, shuffles=0).sites)
        assert passing_relative > kept

    def test_a_scan_reports_what_chance_would_produce(self) -> None:
        result = scan(self.sequence(), BY_NAME["SP1"], shuffles=10)
        assert result.expected_hits > 0
        assert result.enrichment is not None

    def test_a_motif_at_or_below_chance_says_so_in_its_summary(self) -> None:
        result = scan(self.sequence(), BY_NAME["CTCF"], shuffles=10)
        if not result.exceeds_chance and result.sites:
            assert "not evidence of anything" in result.summary

    def test_a_planted_tata_box_is_strongly_enriched(self) -> None:
        # The positive control for the whole scanner: TATAAAA appears in the test sequence and
        # should stand far above its shuffled expectation.
        result = scan(self.sequence(), BY_NAME["TBP"], shuffles=10)
        assert result.sites
        assert result.exceeds_chance

    def test_an_empty_result_summarises_as_no_match(self) -> None:
        plain = NucleicSequence.parse("A" * 300)
        result = scan(plain, BY_NAME["CTCF"], shuffles=0)
        assert not result.sites
        assert "No CTCF match" in result.summary

    def test_a_threshold_outside_zero_to_one_is_refused(self) -> None:
        with pytest.raises(ValueError, match=re.escape("relative score in [0, 1]")):
            scan(self.sequence(), BY_NAME["SP1"], threshold=1.5)

    def test_sites_carry_their_own_coordinates(self) -> None:
        site = scan(self.sequence(), BY_NAME["SP1"], shuffles=0).sites[0]
        assert site.end == site.offset + site.motif.length
        assert len(site.matched_bases) == site.motif.length
        assert str(site).startswith("SP1 at ")

    def test_an_ambiguity_code_is_skipped_rather_than_scored(self) -> None:
        # A window containing N is passed over entirely; the rest of the sequence still scans.
        motif = BY_NAME["SP1"]
        bases = "GGGGCGGGG" + "N" * 20 + "GGGGCGGGG"
        assert _sites_in(bases, motif, DEFAULT_RELATIVE_THRESHOLD)
        assert all(
            "N" not in site.matched_bases
            for site in _sites_in(bases, motif, DEFAULT_RELATIVE_THRESHOLD)
        )

    def test_scanning_all_motifs_orders_by_enrichment(self) -> None:
        scans = scan_all(self.sequence(), shuffles=5)
        enrichments = [s.enrichment or 0.0 for s in scans]
        assert enrichments == sorted(enrichments, reverse=True)

    def test_scanning_all_motifs_covers_the_whole_curated_set(self) -> None:
        assert len(scan_all(self.sequence(), shuffles=0)) == len(MOTIFS)

    def test_the_expectation_is_reproducible_from_the_seed(self) -> None:
        first = scan(self.sequence(), BY_NAME["SP1"], shuffles=8, seed=7)
        second = scan(self.sequence(), BY_NAME["SP1"], shuffles=8, seed=7)
        assert first.expected_hits == second.expected_hits

    def test_zero_shuffles_skips_the_estimate(self) -> None:
        result = scan(self.sequence(), BY_NAME["SP1"], shuffles=0)
        assert result.expected_hits == 0.0
        assert result.enrichment is None


class StubEnsembl:
    """An Ensembl client stub serving fixed overlap records."""

    def __init__(self, records: list[dict[str, Any]]) -> None:
        """Record what to serve."""
        self.records = records
        self.assembly = Assembly.GRCH38

    def overlap(self, interval: GenomicInterval, feature: str) -> list[dict[str, Any]]:
        """Return the fixed records."""
        return self.records


SCN1A_MANE_RECORDS: list[dict[str, Any]] = [
    {
        "id": "ENSR2_86TTNX",
        "description": "promoter",
        "seq_region_name": "2",
        "start": 166128005,
        "end": 166128231,
    },
    {
        "id": "ENSR2_86TTNS",
        "description": "enhancer",
        "seq_region_name": "2",
        "start": 166127800,
        "end": 166128004,
    },
    {
        "id": "ENSR2_53G42LK",
        "description": "CTCF_binding_site",
        "seq_region_name": "2",
        "start": 166128090,
        "end": 166128111,
    },
]

REGION = GenomicInterval(Assembly.GRCH38, "2", 166127500, 166129000)


def annotation_for(records: list[dict[str, Any]]) -> RegulatoryAnnotation:
    """Build an annotation from fixed records, without touching the network."""
    client = cast("EnsemblClient", StubEnsembl(records))
    return fetch_regulatory_features(client, REGION)


class TestRegulatoryFeatures:
    def test_the_kind_comes_from_the_description_not_the_feature_type(self) -> None:
        # Ensembl's `feature_type` is always the literal string "regulatory"; the real kind is
        # in `description`. Reading the wrong field makes every feature "other".
        annotation = annotation_for(SCN1A_MANE_RECORDS)
        kinds = {f.kind for f in annotation.features}
        assert kinds == {
            RegulatoryFeatureKind.PROMOTER,
            RegulatoryFeatureKind.ENHANCER,
            RegulatoryFeatureKind.CTCF_BINDING_SITE,
        }

    def test_coordinates_convert_to_0_based_half_open(self) -> None:
        promoter = annotation_for(SCN1A_MANE_RECORDS).of_kind(
            RegulatoryFeatureKind.PROMOTER
        )[0]
        assert promoter.interval.start == 166128004
        assert promoter.interval.end == 166128231

    def test_features_are_returned_in_coordinate_order(self) -> None:
        starts = [f.interval.start for f in annotation_for(SCN1A_MANE_RECORDS).features]
        assert starts == sorted(starts)

    def test_only_promoter_kinds_count_as_promoter_associated(self) -> None:
        assert RegulatoryFeatureKind.PROMOTER.is_promoter_associated
        assert RegulatoryFeatureKind.PROMOTER_FLANKING.is_promoter_associated
        assert not RegulatoryFeatureKind.ENHANCER.is_promoter_associated
        assert not RegulatoryFeatureKind.CTCF_BINDING_SITE.is_promoter_associated

    def test_an_unknown_description_becomes_other_rather_than_failing(self) -> None:
        annotation = annotation_for(
            [{"id": "X", "description": "something_new", "start": 1, "end": 100}]
        )
        assert annotation.features[0].kind is RegulatoryFeatureKind.OTHER

    def test_every_feature_is_an_aggregated_annotation(self) -> None:
        for feature in annotation_for(SCN1A_MANE_RECORDS).features:
            assert feature.evidence.kind is EvidenceKind.AGGREGATED_ANNOTATION
            assert feature.evidence.cell_type is None

    def test_features_group_by_kind_for_one_track_each(self) -> None:
        grouped = features_by_kind(annotation_for(SCN1A_MANE_RECORDS).features)
        assert len(grouped[RegulatoryFeatureKind.PROMOTER]) == 1

    def test_a_feature_renders_readably(self) -> None:
        promoter = annotation_for(SCN1A_MANE_RECORDS).of_kind(
            RegulatoryFeatureKind.PROMOTER
        )[0]
        assert str(promoter).startswith("promoter 2:166128004-166128231")


class TestAbsenceIsNotEvidence:
    def test_an_empty_region_carries_its_caveat(self) -> None:
        # The important case. An interface must not render "nothing annotated" as a clean
        # promoter, so the caveat is a value it has to handle, not a comment in the source.
        annotation = annotation_for([])
        assert annotation.is_empty
        assert annotation.no_features_note is not None
        assert "weak evidence of absence" in annotation.no_features_note

    def test_the_caveat_names_the_cell_types_the_panel_lacks(self) -> None:
        note = annotation_for([]).no_features_note
        assert note is not None
        assert "AF22" in note

    def test_a_populated_region_carries_no_caveat(self) -> None:
        assert annotation_for(SCN1A_MANE_RECORDS).no_features_note is None

    def test_promoter_support_distinguishes_absent_from_elsewhere(self) -> None:
        target = GenomicInterval(Assembly.GRCH38, "2", 166128100, 166128119)
        supported, why = promoter_support(annotation_for(SCN1A_MANE_RECORDS), target)
        assert supported
        assert "not evidence that it is active in yours" in why

    def test_promoter_support_says_nothing_either_way_on_an_empty_region(self) -> None:
        target = GenomicInterval(Assembly.GRCH38, "2", 166128100, 166128119)
        supported, why = promoter_support(annotation_for([]), target)
        assert not supported
        assert "says nothing either way" in why

    def test_promoter_support_flags_a_target_in_a_non_promoter_feature(self) -> None:
        # Overlapping an enhancer rather than a promoter is a different mechanistic setting.
        target = GenomicInterval(Assembly.GRCH38, "2", 166127850, 166127869)
        supported, why = promoter_support(annotation_for(SCN1A_MANE_RECORDS), target)
        assert not supported
        assert "different" in why

    def test_the_export_carries_the_caveat_and_the_provenance(self) -> None:
        exported = annotation_for([]).to_dict()
        assert exported["caveat"]
        assert "Ensembl Regulatory Build" in exported["provenance"]
        assert exported["evidence_kind"] == "aggregated-annotation"

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps(annotation_for(SCN1A_MANE_RECORDS).to_dict())

    def test_the_summary_counts_each_kind(self) -> None:
        summary = annotation_for(SCN1A_MANE_RECORDS).summary
        assert "1 promoter" in summary
        assert "1 enhancer" in summary


class TestSmallAccessors:
    """Accessors a UI needs, pinned so a refactor cannot quietly drop one."""

    def test_a_cell_type_status_renders_as_its_label(self) -> None:
        assert str(CellTypeStatus.NOT_ASSAYED) == "not-assayed"
        assert str(CellTypeStatus.OBSERVED_BOUND) == "observed-bound"

    def test_an_evidence_kind_renders_as_its_label(self) -> None:
        assert str(EvidenceKind.MOTIF_PREDICTION) == "motif-prediction"

    def test_a_cell_type_renders_as_its_display_name(self) -> None:
        assert str(HELA) == "HeLa"

    def test_a_feature_kind_renders_as_its_label(self) -> None:
        assert str(RegulatoryFeatureKind.PROMOTER) == "promoter"

    def test_the_promoters_accessor_returns_promoter_associated_features(self) -> None:
        annotation = annotation_for(SCN1A_MANE_RECORDS)
        assert [f.feature_id for f in annotation.promoters] == ["ENSR2_86TTNX"]

    def test_a_motif_data_type_renders_as_its_label(self) -> None:
        assert str(MotifDataType.CHIP_SEQ) == "ChIP-seq"

    def test_a_scan_below_chance_names_the_expectation(self) -> None:
        # The branch that tells a reader a hit count is explained by composition alone.
        motif = BY_NAME["SP1"]
        result = MotifScan(
            motif=motif,
            sites=(),
            threshold=0.8,
            sequence_length=100,
            expected_hits=5.0,
        )
        assert "No SP1 match" in result.summary
        populated = MotifScan(
            motif=motif,
            sites=scan(NucleicSequence.parse("GGGGCGGGG" * 10), motif, shuffles=0).sites,
            threshold=0.8,
            sequence_length=90,
            expected_hits=1000.0,
        )
        assert "not evidence of anything" in populated.summary
        assert not populated.exceeds_chance
