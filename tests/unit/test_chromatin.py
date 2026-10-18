"""Tests for the AF22 chromatin adapter and the liftover boundary it sits on.

Almost everything here is about refusing to answer. The tracks are GRCh37 and the engine is
GRCh38, and on SCN1A's TSS those differ by 856 kb, so a conversion this adapter cannot guarantee
must produce a stated reason rather than a number, a signal value from the wrong locus is worse
than no value, because it looks like data.

Run against stubs, so the suite needs none of the 890 MB of real tracks. The tests that use the
real files carry the ``tracks`` marker and skip when they are absent.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest

from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.regulatory.chromatin import (
    ENRICHMENT_THRESHOLD,
    GSE108506,
    Af22ChromatinReader,
    BackgroundModel,
    CoordinateConverter,
    HistoneMark,
    LiftoverError,
    LiftoverProblem,
    SignalReader,
    convert_interval,
    find_tracks,
    open_af22_reader,
    track_directory,
)
from sarna.regulatory.evidence import AF22, HELA, CellTypeStatus, EvidenceKind
from sarna.rules.provenance import Provenance

GRCH38_SHIFT = 856_510
"""Measured offset between GRCh38 and hg19 at SCN1A's MANE Select TSS. The size of defect B15."""


def af22_reader(
    readers: Mapping[HistoneMark, SignalReader],
    converter: CoordinateConverter,
    track_assembly: Assembly = Assembly.GRCH37,
) -> Af22ChromatinReader:
    """Build a reader over AF22 stubs, naming the cell line and source the real factory names.

    The constructor requires both, so that a reader can never label another cell line's tracks
    with AF22's accession. These tests are about AF22 behaviour, so they say so once here.
    """
    return Af22ChromatinReader(
        readers,
        converter,
        cell_type=AF22,
        provenance=GSE108506,
        track_assembly=track_assembly,
    )



class StubConverter:
    """A liftover chain that applies a fixed shift, or misbehaves on demand."""

    def __init__(
        self,
        shift: int = GRCH38_SHIFT,
        contig_override: dict[int, str] | None = None,
        strand_override: dict[int, str] | None = None,
        unmapped: frozenset[int] = frozenset(),
        ambiguous: frozenset[int] = frozenset(),
        drift: dict[int, int] | None = None,
    ) -> None:
        """Record how this chain should behave."""
        self.shift = shift
        self.contig_override = contig_override or {}
        self.strand_override = strand_override or {}
        self.unmapped = unmapped
        self.ambiguous = ambiguous
        self.drift = drift or {}

    def convert_coordinate(
        self, chromosome: str, position: int
    ) -> Sequence[tuple[str, int, str, int]] | None:
        """Return the shifted coordinate, or whatever misbehaviour was configured."""
        if position in self.unmapped:
            return []
        target = self.contig_override.get(position, chromosome)
        strand = self.strand_override.get(position, "+")
        moved = position + self.shift + self.drift.get(position, 0)
        if position in self.ambiguous:
            return [(target, moved, strand, 0), (target, moved + 50_000, strand, 0)]
        return [(target, moved, strand, 0)]


class StubTrack:
    """A bigWig stand-in returning a constant inside one span and a constant elsewhere."""

    def __init__(
        self,
        peak_span: tuple[int, int] | None = None,
        peak_value: float = 100.0,
        background_value: float = 1.0,
        contigs: Mapping[str, int] | None = None,
    ) -> None:
        """Record the shape of the signal to serve."""
        self.peak_span = peak_span
        self.peak_value = peak_value
        self.background_value = background_value
        self._contigs = dict(contigs or {"chr2": 243_199_373})

    def chroms(self) -> Mapping[str, int]:
        """Return the contig lengths this track declares."""
        return self._contigs

    def stats(
        self, chromosome: str, start: int, end: int, type: str = "mean", nBins: int = 1
    ) -> Sequence[float | None]:
        """Return the peak value when the span overlaps the peak, else background."""
        if self.peak_span and start < self.peak_span[1] and end > self.peak_span[0]:
            return [self.peak_value]
        return [self.background_value]


PROMOTER_38 = GenomicInterval(Assembly.GRCH38, "2", 166_127_520, 166_130_020, Strand.MINUS)


class TestConversionRefusesWhatItCannotGuarantee:
    def test_a_clean_interval_converts_and_keeps_its_strand(self) -> None:
        converted = convert_interval(PROMOTER_38, StubConverter(), Assembly.GRCH37)
        assert converted.assembly is Assembly.GRCH37
        assert converted.contig == "2"
        assert converted.start == PROMOTER_38.start + GRCH38_SHIFT
        assert converted.length == PROMOTER_38.length
        assert converted.strand is Strand.MINUS

    def test_the_ucsc_chr_prefix_is_added_for_the_query_and_stripped_from_the_result(
        self,
    ) -> None:
        # The tracks are UCSC-named and the domain model is Ensembl-named.
        converted = convert_interval(PROMOTER_38, StubConverter(), Assembly.GRCH37)
        assert not converted.contig.startswith("chr")

    def test_the_mitochondrion_is_refused_outright(self) -> None:
        # hg19 holds NC_001807 and GRCh37 holds rCRS: different sequences, 2 bp apart in length.
        mito = GenomicInterval(Assembly.GRCH38, "MT", 1000, 2000)
        with pytest.raises(LiftoverError) as caught:
            convert_interval(mito, StubConverter(), Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.MITOCHONDRIAL
        assert "different sequences" in str(caught.value)

    def test_an_unmappable_end_is_refused(self) -> None:
        converter = StubConverter(unmapped=frozenset({PROMOTER_38.start}))
        with pytest.raises(LiftoverError) as caught:
            convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.UNMAPPED

    def test_an_ambiguous_position_is_refused_rather_than_resolved(self) -> None:
        # Taking the first candidate is how a pipeline ends up reading a paralogue.
        converter = StubConverter(ambiguous=frozenset({PROMOTER_38.start}))
        with pytest.raises(LiftoverError) as caught:
            convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.UNMAPPED
        assert "maps to 2 places" in str(caught.value)

    def test_a_contig_change_is_refused(self) -> None:
        converter = StubConverter(contig_override={PROMOTER_38.end - 1: "chr7"})
        with pytest.raises(LiftoverError) as caught:
            convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.CONTIG_CHANGED

    def test_a_strand_flip_is_refused(self) -> None:
        converter = StubConverter(strand_override={PROMOTER_38.end - 1: "-"})
        with pytest.raises(LiftoverError) as caught:
            convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.STRAND_FLIPPED

    def test_an_interval_that_changes_length_is_refused(self) -> None:
        # The two spans no longer describe the same locus, so signal averaged over the second
        # answers a different question from the one asked.
        converter = StubConverter(drift={PROMOTER_38.end - 1: 5000})
        with pytest.raises(LiftoverError) as caught:
            convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert caught.value.problem is LiftoverProblem.SPLIT
        assert "do not describe the same locus" in str(caught.value)

    def test_a_few_bases_of_drift_is_tolerated(self) -> None:
        # A small indel still describes the same locus; refusing here would reject most of the
        # genome for no gain.
        converter = StubConverter(drift={PROMOTER_38.end - 1: 10})
        converted = convert_interval(PROMOTER_38, converter, Assembly.GRCH37)
        assert converted.length == PROMOTER_38.length + 10


class TestWhoseDataItIsBelongsToTheReader:
    """A reader cannot label one cell line's tracks with another's accession.

    ``profile()`` used to take ``cell_type`` per call while provenance came from the module-level
    ``GSE108506`` constant, so asking an AF22 reader for a HeLa profile returned HeLa-labelled
    signal citing AF22's accession. The evidence model exists to stop exactly that class of
    mislabelling, so the cell type and its source are now constructor state: a reader reads one
    cell line's tracks, which is a property of the reader and not of the call.
    """

    def tracks(self) -> dict[HistoneMark, StubTrack]:
        peak = (166_127_520 + GRCH38_SHIFT, 166_130_020 + GRCH38_SHIFT)
        return {HistoneMark.H3K4ME3: StubTrack(peak_span=peak, peak_value=100.0)}

    def test_a_signal_cites_the_source_its_reader_declared(self) -> None:
        source = Provenance.citation("ENCODE HeLa-S3 H3K4me3", "ENCSR000AOC")
        reader = Af22ChromatinReader(
            self.tracks(), StubConverter(), cell_type=HELA, provenance=source
        )
        signal = reader.profile(PROMOTER_38).signals[0]
        assert signal.cell_type is HELA
        assert signal.evidence.provenance == source
        assert "GSE108506" not in str(signal.evidence.provenance)

    def test_the_exported_profile_cites_the_same_source(self) -> None:
        source = Provenance.citation("ENCODE HeLa-S3 H3K4me3", "ENCSR000AOC")
        reader = Af22ChromatinReader(
            self.tracks(), StubConverter(), cell_type=HELA, provenance=source
        )
        exported = reader.profile(PROMOTER_38).to_dict()
        assert exported["provenance"] == str(source)
        assert exported["cell_type"] == HELA.id

    def test_a_reader_cannot_be_built_without_saying_whose_data_it_is(self) -> None:
        # Required, not defaulted: a default would let a caller forget, which is how the
        # mislabelling happened in the first place.
        with pytest.raises(TypeError):
            Af22ChromatinReader(self.tracks(), StubConverter())  # type: ignore[call-arg]


class TestTheReaderReportsReasonsNotNumbers:
    def tracks(self, **kwargs: object) -> dict[HistoneMark, StubTrack]:
        peak = (166_127_520 + GRCH38_SHIFT, 166_130_020 + GRCH38_SHIFT)
        return {
            HistoneMark.H3K4ME3: StubTrack(peak_span=peak, peak_value=100.0),
            HistoneMark.H3K27ME3: StubTrack(peak_span=None, background_value=1.0),
        }

    def test_an_unconvertible_interval_yields_a_problem_and_no_signals(self) -> None:
        reader = af22_reader(
            self.tracks(), StubConverter(unmapped=frozenset({PROMOTER_38.start}))
        )
        profile = reader.profile(PROMOTER_38)
        assert not profile.is_readable
        assert profile.problem is LiftoverProblem.UNMAPPED
        assert profile.signals == ()
        assert "could not be read" in profile.summary

    def test_a_readable_interval_yields_one_signal_per_track(self) -> None:
        reader = af22_reader(self.tracks(), StubConverter())
        profile = reader.profile(PROMOTER_38)
        assert profile.is_readable
        assert {s.mark for s in profile.signals} == {
            HistoneMark.H3K4ME3,
            HistoneMark.H3K27ME3,
        }

    def test_an_interval_already_on_the_track_assembly_is_an_error(self) -> None:
        # Being handed one means an assembly was mixed up upstream, which is defect B15 again.
        reader = af22_reader(self.tracks(), StubConverter())
        already = GenomicInterval(Assembly.GRCH37, "2", 1000, 2000)
        with pytest.raises(ValueError, match="mixed up upstream"):
            reader.profile(already)

    def test_a_contig_absent_from_a_track_is_skipped_not_guessed(self) -> None:
        tracks = {HistoneMark.H3K4ME3: StubTrack(contigs={"chr1": 249_250_621})}
        reader = af22_reader(tracks, StubConverter())
        assert reader.profile(PROMOTER_38).signals == ()

    def test_the_signal_records_where_it_was_actually_read(self) -> None:
        # So a reader can check the conversion rather than trust it.
        reader = af22_reader(self.tracks(), StubConverter())
        signal = reader.profile(PROMOTER_38).signals[0]
        assert signal.source_interval.assembly is Assembly.GRCH37
        assert signal.source_interval.start == PROMOTER_38.start + GRCH38_SHIFT

    def test_marks_are_reported_in_a_stable_order(self) -> None:
        reader = af22_reader(self.tracks(), StubConverter())
        order = [s.mark for s in reader.profile(PROMOTER_38).signals]
        assert order == sorted(order, key=lambda m: list(HistoneMark).index(m))


class TestEnrichmentNotRawSignal:
    """Raw values cannot answer "is this mark here": the tracks have no input normalisation."""

    def test_a_peak_is_called_present(self) -> None:
        peak = (166_127_520 + GRCH38_SHIFT, 166_130_020 + GRCH38_SHIFT)
        tracks = {HistoneMark.H3K4ME3: StubTrack(peak_span=peak, peak_value=100.0)}
        reader = af22_reader(tracks, StubConverter())
        signal = reader.profile(PROMOTER_38).signals[0]
        assert signal.fold_over_background == pytest.approx(100.0)
        assert signal.is_present

    def test_background_level_signal_is_not_called_present(self) -> None:
        # The case the raw test got wrong: every mark reads about 1.0 everywhere.
        tracks = {HistoneMark.H3K4ME3: StubTrack(peak_span=None, background_value=1.0)}
        reader = af22_reader(tracks, StubConverter())
        signal = reader.profile(PROMOTER_38).signals[0]
        assert signal.fold_over_background == pytest.approx(1.0)
        assert not signal.is_present

    def test_the_threshold_is_the_documented_one(self) -> None:
        assert ENRICHMENT_THRESHOLD == 5.0

    def test_the_background_is_sampled_once_per_contig_mark_and_width(self) -> None:
        class Counting(StubTrack):
            calls = 0

            def stats(
                self, chromosome: str, start: int, end: int, type: str = "mean", nBins: int = 1
            ) -> Sequence[float | None]:
                Counting.calls += 1
                return super().stats(chromosome, start, end, type, nBins)

        reader = af22_reader({HistoneMark.H3K4ME3: Counting()}, StubConverter())
        reader.profile(PROMOTER_38)
        after_first = Counting.calls
        reader.profile(PROMOTER_38)
        # Two reads for the second profile (mean and max); the 250 background windows are cached.
        assert Counting.calls - after_first == 2

    def test_the_background_model_records_how_it_was_built(self) -> None:
        reader = af22_reader({HistoneMark.H3K4ME3: StubTrack()}, StubConverter())
        background = reader.profile(PROMOTER_38).signals[0].background
        assert background.windows == 250
        assert background.seed == 11
        assert background.window_length == PROMOTER_38.length

    def test_a_zero_background_gives_infinite_rather_than_dividing_by_zero(self) -> None:
        model = BackgroundModel(
            contig="chr2",
            mark=HistoneMark.H3K4ME3,
            window_length=1000,
            median=0.0,
            percentile_95=0.0,
            windows=10,
            seed=1,
        )
        assert model.fold(5.0) == float("inf")


class TestInterpretation:
    def peak_tracks(self, active: bool, repressive: bool) -> dict[HistoneMark, StubTrack]:
        peak = (166_127_520 + GRCH38_SHIFT, 166_130_020 + GRCH38_SHIFT)
        return {
            HistoneMark.H3K4ME3: StubTrack(
                peak_span=peak if active else None, peak_value=100.0
            ),
            HistoneMark.H3K27ME3: StubTrack(
                peak_span=peak if repressive else None, peak_value=100.0
            ),
        }

    def profile(self, active: bool, repressive: bool) -> object:
        reader = af22_reader(self.peak_tracks(active, repressive), StubConverter())
        return reader.profile(PROMOTER_38)

    def test_active_only_reads_as_accessible(self) -> None:
        summary = self.profile(active=True, repressive=False).summary  # type: ignore[attr-defined]
        assert "accessible promoter" in summary

    def test_repressive_only_reads_as_repressed_and_explains_the_consequence(self) -> None:
        summary = self.profile(active=False, repressive=True).summary  # type: ignore[attr-defined]
        assert "repressed promoter" in summary
        assert "failing to activate" in summary

    def test_both_reads_as_bivalent(self) -> None:
        assert "bivalent" in self.profile(active=True, repressive=True).summary  # type: ignore[attr-defined]

    def test_neither_reads_as_closed_and_names_the_mechanism(self) -> None:
        summary = self.profile(active=False, repressive=False).summary  # type: ignore[attr-defined]
        assert "looks closed" in summary
        assert "cannot activate a locus that is not already transcribed" in summary

    def test_every_mark_explains_what_it_indicates(self) -> None:
        for mark in HistoneMark:
            assert len(mark.interpretation) > 40

    def test_marks_are_classified_as_activating_or_repressive_or_neither(self) -> None:
        assert HistoneMark.H3K4ME3.is_activating
        assert HistoneMark.H3K27ME3.is_repressive
        # H3K36me3 marks gene bodies, so it is neither at a promoter.
        assert not HistoneMark.H3K36ME3.is_activating
        assert not HistoneMark.H3K36ME3.is_repressive


class TestEvidenceTier:
    def test_a_present_mark_is_a_cell_type_observation(self) -> None:
        # The only cell-type-specific evidence this project has.
        peak = (166_127_520 + GRCH38_SHIFT, 166_130_020 + GRCH38_SHIFT)
        tracks = {HistoneMark.H3K4ME3: StubTrack(peak_span=peak, peak_value=100.0)}
        reader = af22_reader(tracks, StubConverter())
        evidence = reader.profile(PROMOTER_38).signals[0].evidence
        assert evidence.kind is EvidenceKind.CELL_TYPE_OBSERVATION
        assert evidence.cell_type is AF22
        assert evidence.status is CellTypeStatus.OBSERVED_BOUND

    def test_an_absent_mark_is_assayed_not_bound_rather_than_not_assayed(self) -> None:
        # The track was read and reported nothing, which is a measurement.
        tracks = {HistoneMark.H3K4ME3: StubTrack(peak_span=None)}
        reader = af22_reader(tracks, StubConverter())
        evidence = reader.profile(PROMOTER_38).signals[0].evidence
        assert evidence.status is CellTypeStatus.ASSAYED_NOT_BOUND
        assert evidence.status.is_known

    def test_the_export_carries_the_assembly_caveat(self) -> None:
        reader = af22_reader({HistoneMark.H3K4ME3: StubTrack()}, StubConverter())
        exported = reader.profile(PROMOTER_38).to_dict()
        assert "GRCh37" in exported["assembly_note"]
        assert "B15" in exported["assembly_note"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        reader = af22_reader({HistoneMark.H3K4ME3: StubTrack()}, StubConverter())
        json.dumps(reader.profile(PROMOTER_38).to_dict())

    def test_an_unreadable_profile_exports_its_reason(self) -> None:
        reader = af22_reader(
            {HistoneMark.H3K4ME3: StubTrack()},
            StubConverter(unmapped=frozenset({PROMOTER_38.start})),
        )
        exported = reader.profile(PROMOTER_38).to_dict()
        assert exported["readable"] is False
        assert exported["problem"] == "unmapped"
        assert exported["signals"] == []


class TestTrackDiscovery:
    def test_marks_are_matched_from_the_geo_filenames(self, tmp_path: Path) -> None:
        for name in (
            "GSM2902406_No0028_AF22_H3K4me3_bs50.bw",
            "GSM2902409_No0032_AF22_H3K27ac_bs50.bw",
        ):
            (tmp_path / name).touch()
        found = find_tracks(tmp_path)
        assert set(found) == {HistoneMark.H3K4ME3, HistoneMark.H3K27AC}

    def test_h3k4me1_is_not_mistaken_for_h3k4me3(self, tmp_path: Path) -> None:
        # Substring matching on mark names has exactly one dangerous near-collision.
        (tmp_path / "AF22_H3K4me1_bs50.bw").touch()
        assert set(find_tracks(tmp_path)) == {HistoneMark.H3K4ME1}

    def test_a_missing_directory_yields_nothing_rather_than_raising(self, tmp_path: Path) -> None:
        assert find_tracks(tmp_path / "absent") == {}

    def test_no_tracks_means_no_reader_rather_than_an_error(self, tmp_path: Path) -> None:
        # 890 MB this project does not ship: absence is the normal case and must not stop a run.
        assert open_af22_reader(tmp_path) is None

    def test_the_directory_can_be_overridden_by_environment(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        monkeypatch.setenv("SARNA_AF22_DIR", str(tmp_path))
        assert track_directory() == tmp_path

    def test_the_default_directory_is_used_without_an_override(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("SARNA_AF22_DIR", raising=False)
        assert track_directory().name == "af22_data"


@pytest.fixture(scope="module")
def real_reader() -> Af22ChromatinReader:
    """Open the real AF22 tracks, skipping the whole class when they are absent."""
    if not find_tracks():
        pytest.skip("AF22 tracks not present")
    from pyliftover import LiftOver

    built = open_af22_reader()
    assert built is not None
    assert isinstance(built._converter, LiftOver)
    return built


@pytest.mark.tracks
class TestAgainstTheRealTracks:
    """Skips when the 890 MB of AF22 bigWigs are absent.

    The housekeeping promoters are the positive control. Without them, "SCN1A has no H3K4me3 in
    AF22" is indistinguishable from "the adapter reads the wrong place".
    """

    @staticmethod
    def window(contig: str, tss: int) -> GenomicInterval:
        return GenomicInterval(Assembly.GRCH38, contig, tss - 500, tss + 500, Strand.PLUS)

    def test_gapdhs_promoter_is_strongly_h3k4me3_marked(
        self, real_reader: Af22ChromatinReader
    ) -> None:
        # The control that validates the whole measurement. GAPDH is transcribed in every cell.
        profile = real_reader.profile(self.window("12", 6_534_517))
        signal = profile.signal(HistoneMark.H3K4ME3)
        assert signal is not None
        assert signal.fold_over_background > 20
        assert signal.is_present

    def test_scn1as_promoter_carries_no_enriched_mark_in_af22(
        self, real_reader: Af22ChromatinReader
    ) -> None:
        # The finding. Meaningful only because the control above passes.
        profile = real_reader.profile(self.window("2", 166_128_019))
        assert profile.is_readable
        assert profile.activating_marks == ()
        assert "looks closed" in profile.summary

    def test_the_tracks_are_grch37(self, real_reader: Af22ChromatinReader) -> None:
        assert real_reader.track_assembly is Assembly.GRCH37

    def test_all_six_marks_are_present(self, real_reader: Af22ChromatinReader) -> None:
        assert len(real_reader.marks) == 6


class TestProfileAccessors:
    def test_asking_for_a_mark_with_no_track_returns_none(self) -> None:
        # Distinct from "the mark is absent": no track was available to look in.
        reader = af22_reader({HistoneMark.H3K4ME3: StubTrack()}, StubConverter())
        profile = reader.profile(PROMOTER_38)
        assert profile.signal(HistoneMark.H3K4ME3) is not None
        assert profile.signal(HistoneMark.H3K9ME3) is None

    def test_a_reader_with_no_tracks_says_so_rather_than_reading_nothing(self) -> None:
        reader = af22_reader({}, StubConverter())
        profile = reader.profile(PROMOTER_38)
        assert profile.is_readable
        assert profile.signals == ()
        assert profile.summary == "No AF22 (lt-NES) tracks were available to read."

    def test_a_converted_interval_past_the_contig_end_yields_no_signal(self) -> None:
        # The track's contig is shorter than the converted coordinates, so there is nothing to
        # read. Returning no signal is right; returning zero would be a measurement.
        tracks = {HistoneMark.H3K4ME3: StubTrack(contigs={"chr2": 1000})}
        reader = af22_reader(tracks, StubConverter())
        assert reader.profile(PROMOTER_38).signals == ()


@pytest.mark.strand
class TestSignalTracks:
    """Binned signal for drawing, which must be in the promoter's own 5'-to-3' order."""

    class Ramp(StubTrack):
        """Returns an ascending ramp across the requested bins, so order is visible."""

        def stats(
            self, chromosome: str, start: int, end: int, type: str = "mean", nBins: int = 1
        ) -> Sequence[float | None]:
            return [float(index) for index in range(nBins)]

    def reader(self) -> Af22ChromatinReader:
        return af22_reader({HistoneMark.H3K4ME3: self.Ramp()}, StubConverter())

    def test_a_plus_strand_track_is_returned_in_reference_order(self) -> None:
        plus = GenomicInterval(Assembly.GRCH38, "2", 166_127_520, 166_130_020, Strand.PLUS)
        track = self.reader().signal_track(plus, HistoneMark.H3K4ME3, bins=5)
        assert track == (0.0, 1.0, 2.0, 3.0, 4.0)

    def test_a_minus_strand_track_is_reversed_into_transcript_order(self) -> None:
        # The promoter is drawn 5'->3' in transcript orientation, so a profile read in reference
        # order would be mirrored, the peak would appear at the wrong end of the drawing.
        minus = GenomicInterval(Assembly.GRCH38, "2", 166_127_520, 166_130_020, Strand.MINUS)
        track = self.reader().signal_track(minus, HistoneMark.H3K4ME3, bins=5)
        assert track == (4.0, 3.0, 2.0, 1.0, 0.0)

    def test_a_mark_with_no_track_returns_none_rather_than_zeros(self) -> None:
        # No data is not an absence of signal, and a renderer must be able to tell them apart.
        assert self.reader().signal_track(PROMOTER_38, HistoneMark.H3K9ME3, bins=5) is None

    def test_an_unconvertible_interval_returns_none(self) -> None:
        reader = af22_reader(
            {HistoneMark.H3K4ME3: self.Ramp()},
            StubConverter(unmapped=frozenset({PROMOTER_38.start})),
        )
        assert reader.signal_track(PROMOTER_38, HistoneMark.H3K4ME3, bins=5) is None

    def test_a_contig_absent_from_the_track_returns_none(self) -> None:
        reader = af22_reader(
            {HistoneMark.H3K4ME3: StubTrack(contigs={"chr1": 249_250_621})}, StubConverter()
        )
        assert reader.signal_track(PROMOTER_38, HistoneMark.H3K4ME3, bins=5) is None

    def test_an_interval_past_the_contig_end_returns_none(self) -> None:
        reader = af22_reader(
            {HistoneMark.H3K4ME3: StubTrack(contigs={"chr2": 1000})}, StubConverter()
        )
        assert reader.signal_track(PROMOTER_38, HistoneMark.H3K4ME3, bins=5) is None
