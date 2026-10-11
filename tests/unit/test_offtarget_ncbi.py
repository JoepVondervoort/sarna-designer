"""The NCBI BLAST client, against recorded replies.

No test here touches the network. The replies below are the shapes NCBI's URL API actually
returns, trimmed; what they pin is that this client reads them correctly and, where it cannot,
says so instead of reporting a clean result.

The recurring risk this file guards is a **false negative**: at 19 nt, a misconfigured search
returns nothing, and "nothing" is indistinguishable from "specific" unless the parameters are
checked. Several tests therefore assert on what is *sent*, not only on what is parsed.
"""

from collections import Counter
from typing import Any

import httpx
import pytest

from sarna.offtarget.models import (
    BY_DATABASE,
    BY_ENTREZ_FILTER,
    GENOMIC_SCOPE,
    TRANSCRIPT_SCOPE,
)
from sarna.offtarget.ncbi import (
    MINIMUM_POLL_SECONDS,
    BlastClient,
    BlastError,
    BlastStillRunning,
    as_fasta,
    candidate_sequences,
    mismatch_positions,
    parameters_for,
    parse_hits,
    parse_report,
)

PUT_REPLY = b"""
<!--QBlastInfoBegin
    RID = ABC123XYZ014
    RTOE = 27
QBlastInfoEnd
-->
"""

WAITING = b"""
QBlastInfoBegin
\tStatus=WAITING
QBlastInfoEnd
"""

READY_WITH_HITS = b"""
QBlastInfoBegin
\tStatus=READY
\tThereAreHits=yes
QBlastInfoEnd
"""

READY_NO_HITS = b"""
QBlastInfoBegin
\tStatus=READY
\tThereAreHits=no
QBlastInfoEnd
"""

UNKNOWN = b"""
QBlastInfoBegin
\tStatus=UNKNOWN
QBlastInfoEnd
"""

REPORT = b"""<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>saRNA-1549</Iteration_query-def>
      <Iteration_hits>
        <Hit>
          <Hit_id>ref|NM_006920.6|</Hit_id>
          <Hit_def>Homo sapiens sodium voltage-gated channel alpha subunit 1</Hit_def>
          <Hit_accession>NM_006920</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>38.2</Hsp_bit-score>
              <Hsp_evalue>0.004</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>19</Hsp_query-to>
              <Hsp_hit-from>455</Hsp_hit-from>
              <Hsp_hit-to>473</Hsp_hit-to>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>19</Hsp_identity>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>19</Hsp_align-len>
            </Hsp>
            <Hsp>
              <Hsp_bit-score>26.3</Hsp_bit-score>
              <Hsp_evalue>12.0</Hsp_evalue>
              <Hsp_query-from>2</Hsp_query-from>
              <Hsp_query-to>15</Hsp_query-to>
              <Hsp_hit-from>980</Hsp_hit-from>
              <Hsp_hit-to>967</Hsp_hit-to>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>13</Hsp_identity>
              <Hsp_gaps>1</Hsp_gaps>
              <Hsp_align-len>14</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
    <Iteration>
      <Iteration_query-def>saRNA-0002 some trailing description</Iteration_query-def>
      <Iteration_hits/>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def client_for(replies: list[bytes]) -> tuple[BlastClient, list[httpx.Request]]:
    """Return a client answering with ``replies`` in order, and the requests it made."""
    seen: list[httpx.Request] = []

    def handler(request: httpx.Request) -> httpx.Response:
        seen.append(request)
        return httpx.Response(200, content=replies[min(len(seen) - 1, len(replies) - 1)])

    transport = httpx.MockTransport(handler)
    return BlastClient(client=httpx.Client(transport=transport)), seen


def form(request: httpx.Request) -> dict[str, str]:
    """Return a posted form as a dictionary."""
    from urllib.parse import parse_qs

    return {k: v[0] for k, v in parse_qs(request.content.decode()).items()}


class TestWhatIsSent:
    """A misconfigured short search returns nothing, which reads as a clean result."""

    def test_short_queries_are_sent_with_parameters_that_can_find_them(self) -> None:
        client, seen = client_for([PUT_REPLY])
        client.submit({"saRNA-001": "GGGTAATTTTGACAGGTCA"}, TRANSCRIPT_SCOPE)
        sent = form(seen[0])
        assert sent["PROGRAM"] == "blastn"
        # Megablast is tuned for long similar sequences and finds nothing at 19 nt.
        assert sent["MEGABLAST"] == "off"
        assert sent["WORD_SIZE"] == "7"
        assert float(sent["EXPECT"]) >= 1000
        # Low-complexity filtering would mask much of a 19-mer.
        assert sent["FILTER"] == "F"
        # Human only. Measured before this was added: a live search of one 19-mer returned 500
        # alignments, 230 of them "near-perfect", in sawfish, plants, grasses and frogs. Real
        # alignments, entirely irrelevant, and enough to swamp every specificity measure.
        assert sent["ENTREZ_QUERY"] == "txid9606[ORGN]"

    def test_the_two_scopes_search_different_databases(self) -> None:
        assert parameters_for(TRANSCRIPT_SCOPE).database == "refseq_rna"
        # The human GRCh38 reference assembly, not the retired `refseq_genomic` (ADR 0024).
        assert parameters_for(GENOMIC_SCOPE).database == "GPIPE/9606/current/ref_top_level"
        with pytest.raises(ValueError, match="unknown search scope"):
            parameters_for("whatever")

    def test_each_scope_records_how_it_was_confined_to_human(self) -> None:
        # Not the same guarantee, and the interface must not word them alike. The genomic
        # assembly holds one species; the transcript database holds every species RefSeq has
        # and is narrowed by a parameter the URL API does not document.
        assert parameters_for(GENOMIC_SCOPE).species_restriction == BY_DATABASE
        assert parameters_for(TRANSCRIPT_SCOPE).species_restriction == BY_ENTREZ_FILTER

    def test_the_genomic_scope_sends_no_entrez_filter_because_it_needs_none(self) -> None:
        client, seen = client_for([PUT_REPLY])
        client.submit({"saRNA-001": "GGGTAATTTTGACAGGTCA"}, GENOMIC_SCOPE)
        sent = form(seen[0])
        assert sent["DATABASE"] == "GPIPE/9606/current/ref_top_level"
        assert "ENTREZ_QUERY" not in sent, (
            "a single-species assembly needs no species filter, and ENTREZ_QUERY is "
            "undocumented, so sending it anyway is a dependency with nothing to gain"
        )

    def test_a_whole_run_goes_in_one_job(self) -> None:
        # One job per candidate would multiply requests by fourteen and risk an IP block for
        # the whole institution.
        client, seen = client_for([PUT_REPLY])
        client.submit({f"saRNA-{i:03d}": "ACGTACGTACGTACGTACG" for i in range(14)},
                      TRANSCRIPT_SCOPE)
        assert len(seen) == 1
        assert form(seen[0])["QUERY"].count(">") == 14

    def test_the_target_dna_is_searched_rather_than_the_rna_core(self) -> None:
        # "Where else does this sequence occur" is a question about DNA.
        records: list[dict[str, Any]] = [
            {"id": "saRNA-2", "rank": 2, "target_site_dna": "AAAACCCCGGGGTTTTACG",
             "sense_strand_rna": "AAAACCCCGGGGUUUUACG"},
            {"id": "saRNA-1", "rank": 1, "target_site_dna": "GGGTAATTTTGACAGGTCA",
             "sense_strand_rna": "GGGUAAUUUUGACAGGUCA"},
            {"id": "saRNA-9", "rank": None, "target_site_dna": "TTTTTTTTTTTTTTTTTTT"},
        ]
        sequences = candidate_sequences(records)
        assert "U" not in "".join(sequences.values())
        # Rejected candidates are not screened: nobody would order them.
        assert list(sequences) == ["saRNA-1", "saRNA-2"]

    def test_an_identifier_with_whitespace_is_refused(self) -> None:
        # BLAST truncates a FASTA header at whitespace, which would silently rename the query
        # and make the hits unattributable.
        with pytest.raises(ValueError, match="whitespace"):
            as_fasta({"saRNA 001": "ACGT"})

    def test_an_empty_search_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no sequences"):
            as_fasta({})


class TestSubmissionAndPolling:
    def test_a_submission_returns_the_request_id_and_the_estimate(self) -> None:
        client, _ = client_for([PUT_REPLY])
        submitted = client.submit({"saRNA-001": "ACGTACGTACGTACGTACG"}, TRANSCRIPT_SCOPE)
        assert submitted.request_id == "ABC123XYZ014"
        assert submitted.estimated_seconds == 27

    def test_a_reply_with_no_request_id_is_an_error_not_an_empty_result(self) -> None:
        client, _ = client_for([b"<html>Your search was rejected</html>"])
        with pytest.raises(BlastError, match="no RID"):
            client.submit({"saRNA-001": "ACGTACGTACGTACGTACG"}, TRANSCRIPT_SCOPE)

    def test_a_waiting_search_is_not_an_error(self) -> None:
        client, _ = client_for([WAITING])
        assert client.status("ABC123") == ("WAITING", False)
        with pytest.raises(BlastStillRunning):
            client.fetch("ABC123", parameters_for(TRANSCRIPT_SCOPE))

    def test_an_expired_request_id_says_so_rather_than_reporting_no_hits(self) -> None:
        client, _ = client_for([UNKNOWN])
        with pytest.raises(BlastError, match="no longer recognises"):
            client.fetch("ABC123", parameters_for(TRANSCRIPT_SCOPE))

    def test_a_finished_search_with_no_hits_returns_nothing_without_fetching_a_report(
        self,
    ) -> None:
        client, seen = client_for([READY_NO_HITS])
        report = client.fetch("ABC123", parameters_for(TRANSCRIPT_SCOPE))
        assert report.hits == ()
        assert report.capped == frozenset()
        assert len(seen) == 1

    def test_the_polling_interval_respects_ncbis_policy(self) -> None:
        assert MINIMUM_POLL_SECONDS >= 15

    def test_an_http_error_explains_that_the_designs_are_unaffected(self) -> None:
        transport = httpx.MockTransport(lambda _: httpx.Response(503))
        client = BlastClient(client=httpx.Client(transport=transport))
        with pytest.raises(BlastError, match="usage policy"):
            client.status("ABC123")

    def test_an_unreachable_service_is_reported_as_such(self) -> None:
        def boom(_: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("no route to host")

        client = BlastClient(client=httpx.Client(transport=httpx.MockTransport(boom)))
        with pytest.raises(BlastError, match="could not reach NCBI"):
            client.status("ABC123")


class TestParsingTheReport:
    def test_every_alignment_becomes_its_own_hit(self) -> None:
        # Two HSPs against one accession are two separate places the duplex could bind.
        hits = parse_hits(REPORT)
        assert len(hits) == 2
        assert {h.accession for h in hits} == {"NM_006920"}

    def test_a_hit_carries_the_numbers_a_researcher_judges_it_by(self) -> None:
        best = parse_hits(REPORT)[0]
        assert best.candidate_id == "saRNA-1549"
        assert best.identities == 19
        assert best.alignment_length == 19
        assert best.gaps == 0
        assert best.evalue == 0.004
        assert best.bit_score == 38.2
        assert best.subject_strand == "plus"
        assert best.is_perfect(19) is True

    def test_a_negative_frame_is_reported_as_the_minus_strand(self) -> None:
        second = parse_hits(REPORT)[1]
        assert second.subject_strand == "minus"
        assert second.gaps == 1
        assert second.is_perfect(19) is False

    def test_hits_are_attributed_to_the_candidate_that_produced_them(self) -> None:
        # A multi-FASTA job returns one Iteration per query; mixing them up would attach one
        # candidate's risk to another.
        assert {h.candidate_id for h in parse_hits(REPORT)} == {"saRNA-1549"}

    def test_a_trailing_description_does_not_become_part_of_the_candidate_id(self) -> None:
        report = REPORT.replace(b"<Iteration_hits/>", b"""<Iteration_hits><Hit>
          <Hit_accession>NM_000001</Hit_accession><Hit_def>x</Hit_def><Hit_hsps><Hsp>
          <Hsp_identity>12</Hsp_identity><Hsp_align-len>12</Hsp_align-len><Hsp_gaps>0</Hsp_gaps>
          <Hsp_evalue>7.0</Hsp_evalue><Hsp_bit-score>24.0</Hsp_bit-score>
          </Hsp></Hit_hsps></Hit></Iteration_hits>""")
        assert "saRNA-0002" in {h.candidate_id for h in parse_hits(report)}

    def test_the_hit_cap_is_honoured(self) -> None:
        assert len(parse_hits(REPORT, max_hits=1)) == 1


def multi_query_report(counts: dict[str, int]) -> bytes:
    """Build a report with a chosen number of alignments for each named query."""
    def alignment(index: int) -> bytes:
        return (
            f"<Hit><Hit_accession>NM_{index:06d}</Hit_accession><Hit_def>x</Hit_def>"
            f"<Hit_hsps><Hsp><Hsp_identity>14</Hsp_identity>"
            f"<Hsp_align-len>19</Hsp_align-len><Hsp_gaps>0</Hsp_gaps>"
            f"<Hsp_evalue>3.0</Hsp_evalue><Hsp_bit-score>24.0</Hsp_bit-score>"
            f"</Hsp></Hit_hsps></Hit>"
        ).encode()

    iterations = b"".join(
        b"<Iteration><Iteration_query-def>" + name.encode() + b"</Iteration_query-def>"
        b"<Iteration_hits>" + b"".join(alignment(i) for i in range(count))
        + b"</Iteration_hits></Iteration>"
        for name, count in counts.items()
    )
    return b"<BlastOutput><BlastOutput_iterations>" + iterations + (
        b"</BlastOutput_iterations></BlastOutput>"
    )


class TestOneQueryCannotStarveTheRest:
    """The defect behind ADR 0023, pinned so it cannot come back.

    A single budget shared across a multi-query report was consumed by whichever query came
    first. In the run that exposed it, one candidate carried 500 alignments and the other 21
    carried zero, and those zeros were printed in the same ink as a real clean result. The
    candidates had not been searched clean. They had not been read.
    """

    def test_a_promiscuous_first_query_does_not_silence_the_later_ones(self) -> None:
        report = multi_query_report({"saRNA-A": 600, "saRNA-B": 4, "saRNA-C": 3})
        parsed = parse_report(report, per_candidate=100)
        found = Counter(h.candidate_id for h in parsed.hits)
        assert found["saRNA-A"] == 100, "the promiscuous query is capped"
        assert found["saRNA-B"] == 4, "and the next query is still read in full"
        assert found["saRNA-C"] == 3

    def test_a_capped_candidate_is_named_so_its_counts_read_as_floors(self) -> None:
        parsed = parse_report(
            multi_query_report({"saRNA-A": 600, "saRNA-B": 4}), per_candidate=100
        )
        # Three outcomes, not two: B found four, A found at least a hundred, and neither of
        # those is "nothing". Only A's counts are floors.
        assert parsed.capped == frozenset({"saRNA-A"})

    def test_a_query_that_exactly_fills_its_budget_is_not_called_capped(self) -> None:
        # Off by one in the honest direction: reporting "there may be more" when there is not
        # would be a false warning, and a check that cries wolf gets ignored.
        parsed = parse_report(multi_query_report({"saRNA-A": 5}), per_candidate=5)
        assert len(parsed.hits) == 5
        assert parsed.capped == frozenset()

    def test_a_budget_of_zero_is_refused_rather_than_reporting_everything_clean(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            parse_report(REPORT, per_candidate=0)

    def test_a_document_declaring_entities_is_refused(self) -> None:
        # xml.etree expands internal entities, so a crafted reply could exhaust memory.
        bomb = b'<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "aaaa">]><BlastOutput/>'
        with pytest.raises(BlastError, match="entities"):
            parse_hits(bomb)

    def test_ncbis_own_doctype_is_accepted(self) -> None:
        # Every genuine reply opens with one. An earlier guard refused all DOCTYPEs and so
        # rejected every real NCBI report -- found only by running a live search.
        real = REPORT.replace(
            b'<?xml version="1.0"?>',
            b'<?xml version="1.0"?>\n<!DOCTYPE BlastOutput PUBLIC '
            b'"-//NCBI//NCBI BlastOutput/EN" '
            b'"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
        )
        assert len(parse_hits(real)) == 2

    def test_unparseable_xml_is_an_error_rather_than_an_empty_result(self) -> None:
        with pytest.raises(BlastError, match="could not parse"):
            parse_hits(b"<BlastOutput><unclosed>")

    def test_a_missing_field_does_not_crash_the_whole_report(self) -> None:
        minimal = b"""<?xml version="1.0"?><BlastOutput><BlastOutput_iterations><Iteration>
          <Iteration_query-def>saRNA-1</Iteration_query-def><Iteration_hits><Hit>
          <Hit_accession>NM_1</Hit_accession><Hit_hsps><Hsp>
          <Hsp_identity>15</Hsp_identity></Hsp></Hit_hsps></Hit></Iteration_hits>
          </Iteration></BlastOutput_iterations></BlastOutput>"""
        only = parse_hits(minimal)[0]
        assert only.identities == 15
        assert only.alignment_length == 0
        assert only.evalue == float("inf")


class TestWhichBasesDidNotMatch:
    """"17 of 19" never says which two, and which two is the question a reader is asking.

    BLAST already sends the aligned strings; these read the mismatched query positions out of
    them so the panel can draw where a partial match falls rather than only how large it is.
    Every expectation below is counted off the strings by hand.
    """

    def test_a_perfect_alignment_has_no_mismatches(self) -> None:
        assert mismatch_positions("ACGTACG", "|||||||", 1) == ()

    def test_the_mismatched_query_positions_are_returned(self) -> None:
        #  position: 1234567
        #     query: ACGTACG
        #   midline: ||| |||   -> column 4 is query position 4
        assert mismatch_positions("ACGTACG", "||| |||", 1) == (4,)

    def test_two_mismatches_in_a_nineteen_base_core(self) -> None:
        query = "ACGTACGTACGTACGTACG"
        midline = "|||||||||||||||||  "
        # Columns 18 and 19 carry no bar, and with no gaps those are query positions 18 and 19.
        assert mismatch_positions(query, midline, 1) == (18, 19)

    def test_positions_are_offset_by_where_the_alignment_starts(self) -> None:
        # An alignment starting at query base 3 must not report its first column as position 1.
        assert mismatch_positions("ACGTA", "|| ||", 3) == (5,)

    def test_a_gap_in_the_query_is_attributed_rather_than_dropped(self) -> None:
        #     query: AC-GT   columns:      1 2 3 4 5
        #   midline: || ||   column 3 is a gap in the query, so it belongs to no position of
        #                    its own and is attributed to position 2, the base before it.
        assert mismatch_positions("AC-GT", "|| ||", 1) == (2,)

    def test_a_report_without_alignment_strings_says_so_rather_than_claiming_a_match(
        self,
    ) -> None:
        # An empty mismatch list drawn as "every base matched" would turn a missing field into
        # a finding, which is the error this whole area keeps making.
        bare = parse_hits(REPORT)[0]
        assert bare.has_alignment_detail is False
        assert bare.mismatch_positions == ()

    def test_a_real_report_carries_the_positions(self) -> None:
        detailed = REPORT.replace(
            b"<Hsp_identity>19</Hsp_identity>",
            b"<Hsp_qseq>ACGTACGTACGTACGTACG</Hsp_qseq>"
            b"<Hsp_hseq>ACGTACGTACGTACGTAAA</Hsp_hseq>"
            b"<Hsp_midline>|||||||||||||||||  </Hsp_midline>"
            b"<Hsp_identity>17</Hsp_identity>",
            1,
        )
        hit = parse_hits(detailed)[0]
        assert hit.has_alignment_detail is True
        assert hit.mismatch_positions == (18, 19)
        assert hit.to_dict()["mismatch_positions"] == [18, 19]


class TestTheCapKeepsTheStrongestNotTheFirst:
    """On a chromosome-level database, document order is nearly the worst order to cap in.

    BLAST groups a report by database sequence, and one chromosome carries thousands of
    alignments. Taking the first N fills the budget with chromosome 1's weak partial matches
    and never reaches the strong alignment on chromosome 7, which is the only one that would
    have changed anybody's mind.
    """

    @staticmethod
    def _report(scores: list[tuple[str, float, int]]) -> bytes:
        hits = "".join(
            f"<Hit><Hit_accession>{accession}</Hit_accession><Hit_def>x</Hit_def>"
            f"<Hit_hsps><Hsp><Hsp_identity>{identities}</Hsp_identity>"
            f"<Hsp_align-len>19</Hsp_align-len><Hsp_gaps>0</Hsp_gaps>"
            f"<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>19</Hsp_query-to>"
            f"<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>19</Hsp_hit-to>"
            f"<Hsp_hit-frame>1</Hsp_hit-frame><Hsp_evalue>1.0</Hsp_evalue>"
            f"<Hsp_bit-score>{score}</Hsp_bit-score></Hsp></Hit_hsps></Hit>"
            for accession, score, identities in scores
        )
        return (
            "<BlastOutput><BlastOutput_iterations><Iteration>"
            "<Iteration_query-def>saRNA-1</Iteration_query-def>"
            f"<Iteration_hits>{hits}</Iteration_hits>"
            "</Iteration></BlastOutput_iterations></BlastOutput>"
        ).encode()

    def test_a_strong_alignment_reported_last_still_survives_the_cap(self) -> None:
        # Three weak alignments first, then the one that matters. A document-order cap of 3
        # drops exactly the alignment a researcher is screening for.
        report = self._report([
            ("NC_000001.11", 26.0, 13),
            ("NC_000001.11", 26.0, 13),
            ("NC_000001.11", 26.0, 13),
            ("NC_000007.14", 38.2, 19),
        ])
        kept = parse_report(report, per_candidate=3)
        assert len(kept.hits) == 3
        assert any(h.identities == 19 for h in kept.hits), (
            "the 19 of 19 match was dropped in favour of three 13 of 19 ones"
        )

    def test_what_is_kept_is_ordered_strongest_first(self) -> None:
        kept = parse_report(self._report([
            ("NC_000001.11", 26.0, 13),
            ("NC_000007.14", 38.2, 19),
            ("NC_000002.12", 30.1, 16),
        ]))
        assert [h.identities for h in kept.hits] == [19, 16, 13]

    def test_dropping_nothing_does_not_mark_the_candidate_capped(self) -> None:
        kept = parse_report(self._report([("NC_000001.11", 26.0, 13)]), per_candidate=5)
        assert kept.capped == frozenset()

    def test_candidates_keep_the_order_blast_reported_them_in(self) -> None:
        report = multi_query_report({"saRNA-B": 2, "saRNA-A": 2})
        assert [h.candidate_id for h in parse_report(report).hits] == [
            "saRNA-B", "saRNA-B", "saRNA-A", "saRNA-A",
        ]
