"""Live Ensembl tests, excluded from the default run.

These exist because the offline tests cannot catch the one class of failure that matters most
here: Ensembl changing the shape or content of a response. The MANE annotation bug was exactly
that, `expand=1` alone returns no MANE field, so every gene looked as though it had no MANE
Select transcript, and no amount of recorded-payload testing would have found it.

Run with::

    .venv/bin/python -m pytest -m network

They are slow and depend on an endpoint that intermittently returns HTTP 500, which the client
retries. Failures here mean "check Ensembl", not "the engine is broken".
"""

import pytest

from sarna.design.pipeline import DesignResult, design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.coordinates import Strand
from sarna.genome.ensembl import EnsemblClient, GeneNotFoundError
from sarna.genome.promoter import Promoter

pytestmark = pytest.mark.network

SCN1A_MANE_TSS = 166128019
"""Verified against Ensembl 2026-10-06, GRCh38, transcript ENST00000674923 / NM_001165963.4."""


@pytest.fixture(scope="module")
def client() -> EnsemblClient:
    return EnsemblClient()


class TestLiveLookup:
    def test_scn1a_resolves_to_a_minus_strand_gene_on_chromosome_2(
        self, client: EnsemblClient
    ) -> None:
        gene = client.gene("SCN1A")
        assert gene.ensembl_id == "ENSG00000144285"
        assert gene.contig == "2"
        assert gene.strand is Strand.MINUS

    def test_scn1a_has_a_mane_select_transcript(self, client: EnsemblClient) -> None:
        # The regression test for the mane=1 bug. Without that parameter this raises.
        transcript = client.select_transcript(client.gene("SCN1A"))
        assert transcript.ensembl_id == "ENST00000674923"
        assert transcript.refseq_match == "NM_001165963.4"
        assert transcript.tss_coordinate == SCN1A_MANE_TSS

    def test_scn1a_has_several_distinct_transcription_start_sites(
        self, client: EnsemblClient
    ) -> None:
        # The reason this project refuses to speak of "the promoter" of a gene. Verified at 5
        # distinct protein-coding TSSs on 2026-10-06; asserted loosely so an annotation update
        # does not fail the suite, because the point is "more than one", not the exact count.
        transcripts = client.transcripts(client.gene("SCN1A"))
        coding_tss = {t.tss_coordinate for t in transcripts if t.is_protein_coding}
        assert len(coding_tss) > 1

    def test_an_unknown_symbol_is_reported_as_not_found(self, client: EnsemblClient) -> None:
        with pytest.raises(GeneNotFoundError):
            client.gene("NOTAREALGENESYMBOL")


@pytest.mark.strand
class TestLivePromoter:
    def test_the_promoter_window_reproduces_the_transcript_tss(
        self, client: EnsemblClient
    ) -> None:
        # The cross-check that matters for a minus-strand gene: the promoter computes the TSS
        # from its own interval and tss_start, by different arithmetic from the transcript's.
        transcript = client.select_transcript(client.gene("SCN1A"))
        promoter = client.promoter_for(transcript)
        assert promoter.tss_coordinate == transcript.tss_coordinate

    def test_the_window_has_the_requested_width_and_anchor(
        self, client: EnsemblClient
    ) -> None:
        transcript = client.select_transcript(client.gene("SCN1A"))
        promoter = client.promoter_for(transcript, upstream_bp=1000, downstream_bp=200)
        assert promoter.length == 1200
        assert promoter.tss_start is not None
        assert promoter.tss_start.value == -1000

    def test_upstream_is_toward_higher_coordinates_on_the_minus_strand(
        self, client: EnsemblClient
    ) -> None:
        transcript = client.select_transcript(client.gene("SCN1A"))
        promoter = client.promoter_for(transcript, upstream_bp=1000, downstream_bp=200)
        assert promoter.interval is not None
        assert promoter.interval.end == transcript.tss_coordinate + 1000 + 1

    def test_the_fetched_sequence_is_real_dna_of_the_right_length(
        self, client: EnsemblClient
    ) -> None:
        transcript = client.select_transcript(client.gene("SCN1A"))
        promoter = client.promoter_for(transcript, upstream_bp=500, downstream_bp=100)
        assert len(promoter.sequence) == 600
        assert set(promoter.sequence.bases) <= set("ACGT")


class TestEveryDesignedTargetIsWhereItSaysItIs:
    """The end-to-end coordinate check, against the reference rather than against ourselves.

    Every other test in this project asks the engine to agree with the engine. This one takes a
    finished design, reads the genomic interval the candidate reports, asks Ensembl separately
    what sequence lives there, and requires it to be the target site, reverse-complemented
    when the gene is on the minus strand.

    It is the only test that can catch an off-by-one or an inverted strand in
    ``Promoter.genomic_interval``, ``offset_of`` and ``span_of``, which CLAUDE.md names as
    defect B1's home. A sign error there produces candidates that are internally consistent,
    pass every rule, export cleanly, and point at the wrong piece of DNA, the exact failure
    this project exists to prevent, and one no amount of recorded-payload testing would reveal.

    SCN1A is deliberate: it is on the minus strand, so the flip is exercised rather than
    skipped.
    """

    @staticmethod
    def designs(client: EnsemblClient, symbol: str) -> tuple[Promoter, DesignResult]:
        transcript = client.select_transcript(client.gene(symbol))
        promoter = client.promoter_for(transcript, upstream_bp=2000, downstream_bp=500)
        return promoter, design_sarnas(promoter, WANG_2015.resolve())

    def test_a_minus_strand_gene_places_every_passing_candidate_correctly(
        self, client: EnsemblClient
    ) -> None:
        promoter, result = self.designs(client, "SCN1A")
        assert promoter.interval is not None
        assert promoter.interval.strand is Strand.MINUS, "SCN1A must exercise the flip"
        ranked = result.ranked
        assert ranked, "no candidate passed, so nothing was verified"
        for candidate in ranked:
            interval = candidate.interval
            assert interval is not None, f"{candidate.id} carries no genomic interval"
            reference = client.sequence(interval)
            # The reference genome holds the plus strand; the target site is read in the
            # transcript's direction, so on the minus strand the two are reverse complements.
            expected = (
                reference.reverse_complement() if interval.strand is Strand.MINUS else reference
            )
            assert expected.bases == candidate.target.bases, (
                f"{candidate.id} reports {interval.contig}:{interval.start}-"
                f"{interval.end} on the {interval.strand} strand and the sequence "
                f"{candidate.target.bases}, but GRCh38 holds {expected.bases} there"
            )

    def test_a_plus_strand_gene_places_its_candidates_correctly_too(
        self, client: EnsemblClient
    ) -> None:
        # The flip must not be applied where it does not belong. CHD2 is on the plus strand.
        promoter, result = self.designs(client, "CHD2")
        assert promoter.interval is not None
        assert promoter.interval.strand is Strand.PLUS, "CHD2 must exercise the unflipped path"
        ranked = result.ranked
        assert ranked, "no candidate passed, so nothing was verified"
        for candidate in ranked[:10]:
            assert candidate.interval is not None
            assert client.sequence(candidate.interval).bases == candidate.target.bases

    def test_the_tss_arithmetic_agrees_with_the_coordinates(
        self, client: EnsemblClient
    ) -> None:
        """A candidate's stated distance from the TSS must match its actual coordinates.

        Reported independently of each other: ``tss_span`` comes from offset arithmetic inside
        the promoter, the interval from genomic arithmetic. They can disagree.
        """
        transcript = client.select_transcript(client.gene("SCN1A"))
        _promoter, result = self.designs(client, "SCN1A")
        tss = transcript.tss_coordinate
        for candidate in result.ranked[:10]:
            assert candidate.tss_span is not None
            assert candidate.interval is not None
            # Minus strand: upstream is toward higher coordinates, and the 5' end of the target
            # is the interval's high end. There is no position zero, so -1 abuts +1.
            start = candidate.tss_span[0].value
            offset = -start if start < 0 else start - 1
            assert candidate.interval.end - 1 == tss + offset, (
                f"{candidate.id} says it starts at {start} relative to the TSS at "
                f"{tss}, which is coordinate {tss + offset}, but its interval ends at "
                f"{candidate.interval.end - 1}"
            )
