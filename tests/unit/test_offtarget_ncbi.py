"""The NCBI BLAST client, against recorded replies.

No test here touches the network. The replies below are the shapes NCBI's URL API actually
returns, trimmed; what they pin is that this client reads them correctly and, where it cannot,
says so instead of reporting a clean result.

The recurring risk this file guards is a **false negative**: at 19 nt, a misconfigured search
returns nothing, and "nothing" is indistinguishable from "specific" unless the parameters are
checked. Several tests therefore assert on what is *sent*, not only on what is parsed.
"""

from typing import Any

import httpx
import pytest

from sarna.offtarget.models import GENOMIC_SCOPE, TRANSCRIPT_SCOPE
from sarna.offtarget.ncbi import (
    MINIMUM_POLL_SECONDS,
    BlastClient,
    BlastError,
    BlastStillRunning,
    as_fasta,
    candidate_sequences,
    parameters_for,
    parse_hits,
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

    def test_the_two_scopes_search_different_databases(self) -> None:
        assert parameters_for(TRANSCRIPT_SCOPE).database == "refseq_rna"
        assert parameters_for(GENOMIC_SCOPE).database == "refseq_genomic"
        with pytest.raises(ValueError, match="unknown search scope"):
            parameters_for("whatever")

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
        assert client.fetch("ABC123", parameters_for(TRANSCRIPT_SCOPE)) == ()
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

    def test_a_document_declaring_entities_is_refused(self) -> None:
        # xml.etree expands entities, so a crafted reply could exhaust memory.
        bomb = b'<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "aaaa">]><BlastOutput/>'
        with pytest.raises(BlastError, match="DTD or entities"):
            parse_hits(bomb)

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
