"""The expression route: a second, optional source that must not take the picker down with it.

Two properties matter more than the numbers here. Gene lookup must not depend on GTEx, because
lookup is on the path to designing anything and GTEx is not; and a failed retrieval must read
as a failed retrieval rather than as an absence of expression. Both are tested by making GTEx
fail on purpose.
"""

import json
from pathlib import Path
from typing import Any

import pytest

from sarna.app.expression_view import expression_overview, unavailable
from sarna.app.service import DesignService
from sarna.app.web import AppRouter
from sarna.expression.gtex import GeneNotInGtexError, GtexClient
from sarna.expression.models import ExpressionStatus
from sarna.sources.json_http import SourceCacheMissError, SourceUnavailableError
from tests.unit.test_app import service_for
from tests.unit.test_gtex import DEFAULT_ROUTES, StubTransport


def gtex_for(routes: dict[str, object] | None = None) -> GtexClient:
    """Build a GTEx client backed by the recorded payloads, with retries made instant."""
    return GtexClient(
        client=StubTransport(DEFAULT_ROUTES if routes is None else routes),  # type: ignore[arg-type]
        backoff_seconds=0.0, offline=False,
    )


def service_with_gtex(
    tmp_path: Path, gtex: GtexClient | None = None, failure: Exception | None = None
) -> DesignService:
    """Build a service whose Ensembl and GTEx clients are both recorded, never live."""
    service = service_for(tmp_path)

    def factory(_: bool) -> GtexClient:
        if failure is not None:
            raise failure
        return gtex or gtex_for()

    service.gtex_factory = factory
    return service


class TestTheRouteIsSeparateFromGeneLookup:
    def test_the_expression_route_answers_on_its_own_path(self, tmp_path: Path) -> None:
        router = AppRouter(service_with_gtex(tmp_path))
        status, _, body = router.dispatch("GET", "/api/genes/SCN1A/expression")
        assert status == 200
        assert json.loads(body)["dataset_id"] == "gtex_v10"

    def test_gene_lookup_still_answers_without_touching_gtex(self, tmp_path: Path) -> None:
        # The GTEx factory raises on any call. Lookup must not call it.
        router = AppRouter(service_with_gtex(tmp_path, failure=AssertionError("called")))
        status, _, body = router.dispatch("GET", "/api/genes/SCN1A")
        assert status == 200
        assert json.loads(body)["default_transcript"] == "ENST00000674923"

    def test_an_unknown_subpath_is_a_404_rather_than_a_gene_named_slash(
        self, tmp_path: Path
    ) -> None:
        router = AppRouter(service_with_gtex(tmp_path))
        assert router.dispatch("GET", "/api/genes/SCN1A/nonsense")[0] == 404

    def test_the_offline_flag_reaches_the_expression_route(self, tmp_path: Path) -> None:
        seen: list[bool] = []
        service = service_for(tmp_path)

        def factory(offline: bool) -> GtexClient:
            seen.append(offline)
            return gtex_for()

        service.gtex_factory = factory
        AppRouter(service).dispatch("GET", "/api/genes/SCN1A/expression?offline=true")
        assert seen == [True]


class TestAFailureReadsAsAFailure:
    @pytest.mark.parametrize(
        ("failure", "code"),
        [
            (SourceCacheMissError("not cached"), "offline-cache-miss"),
            (GeneNotInGtexError("no record"), "gtex-no-record"),
            (SourceUnavailableError("portal down"), "gtex-unavailable"),
        ],
    )
    def test_each_failure_is_reported_with_its_own_code_and_no_transcripts(
        self, tmp_path: Path, failure: Exception, code: str
    ) -> None:
        router = AppRouter(service_with_gtex(tmp_path, failure=failure))
        status, _, body = router.dispatch("GET", "/api/genes/SCN1A/expression")
        payload = json.loads(body)
        assert status == 200, "the page is already showing transcripts; this is not fatal"
        assert payload["available"] is False
        assert payload["code"] == code
        assert payload["transcripts"] == {}

    def test_the_note_refuses_to_imply_anything_about_expression(
        self, tmp_path: Path
    ) -> None:
        router = AppRouter(
            service_with_gtex(tmp_path, failure=SourceUnavailableError("down"))
        )
        payload = json.loads(router.dispatch("GET", "/api/genes/SCN1A/expression")[2])
        assert "says nothing about whether these transcripts are expressed" in payload["note"]

    def test_a_gene_ensembl_cannot_resolve_is_still_an_error(self, tmp_path: Path) -> None:
        # The question being wrong is different from the evidence being missing.
        router = AppRouter(service_with_gtex(tmp_path))
        assert router.dispatch("GET", "/api/genes/NOTAGENE/expression")[0] == 404


class TestWhatThePickerIsGiven:
    def payload(self, tmp_path: Path) -> dict[str, Any]:
        router = AppRouter(service_with_gtex(tmp_path))
        decoded = json.loads(router.dispatch("GET", "/api/genes/SCN1A/expression")[2])
        assert isinstance(decoded, dict)
        return decoded

    def test_it_is_keyed_by_transcript_so_the_page_annotates_its_own_list(
        self, tmp_path: Path
    ) -> None:
        # Keyed rather than ordered on purpose: the page already has its list in annotation
        # order, and handing it a sorted list would invite re-sorting by expression.
        transcripts = self.payload(tmp_path)["transcripts"]
        assert isinstance(transcripts, dict)
        assert "ENST00000674923" in transcripts

    def test_a_silent_transcript_says_measured_at_zero_not_no_data(
        self, tmp_path: Path
    ) -> None:
        entry = self.payload(tmp_path)["transcripts"]["ENST00000674923"]
        assert entry["status"] == "measured"
        assert entry["label"] == "measured at 0.00 TPM in all 3 sample groups"
        assert "no data" not in entry["label"].lower()

    def test_an_all_zero_brain_summary_names_no_single_region(self, tmp_path: Path) -> None:
        # Naming whichever region sorted first would read as though it were distinguished.
        entry = self.payload(tmp_path)["transcripts"]["ENST00000674923"]
        assert entry["brain_label"] == "0.00 TPM in every one of the 2"

    def test_an_expressed_transcript_names_its_highest_group(self, tmp_path: Path) -> None:
        entry = self.payload(tmp_path)["transcripts"]["ENST00000303395"]
        assert entry["label"] == "up to 1.63 TPM in Brain - Cortex"
        assert entry["brain_label"] == "1.63 TPM in Brain - Cortex"

    def test_a_transcript_absent_from_the_release_is_labelled_as_such(
        self, tmp_path: Path
    ) -> None:
        entry = self.payload(tmp_path)["transcripts"]["ENST00000409050"]
        assert entry["status"] == "not-in-release"
        assert entry["label"] == "not in this release, so nothing was quantified"
        assert entry["top"] == []

    def test_every_median_carries_its_donor_count(self, tmp_path: Path) -> None:
        entry = self.payload(tmp_path)["transcripts"]["ENST00000303395"]
        cortex = next(m for m in entry["top"] if m["tissue_id"] == "Brain_Cortex")
        assert cortex["sample_count"] == 270

    def test_the_cell_line_note_answers_the_af22_question_directly(
        self, tmp_path: Path
    ) -> None:
        payload = self.payload(tmp_path)
        assert "AF22" in payload["caveats"][0], "the fact must travel with the profile"
        note = payload["cell_line_note"]
        assert "AF22" in note
        assert "statement about GTEx and not about those cells" in note

    def test_the_silent_count_is_reported_beside_the_measured_count(
        self, tmp_path: Path
    ) -> None:
        payload = self.payload(tmp_path)
        assert payload["measured_transcripts"] == 3
        assert payload["silent_transcripts"] == 1

    def test_a_silent_transcript_lists_no_highest_groups(self, tmp_path: Path) -> None:
        # Its first five groups would be an alphabetical slice of zeros, which reads as
        # though those tissues were distinguished. The flag says why the table is empty.
        entry = self.payload(tmp_path)["transcripts"]["ENST00000674923"]
        assert entry["silent"] is True
        assert entry["top"] == []

    def test_an_expressed_transcript_is_not_marked_silent(self, tmp_path: Path) -> None:
        entry = self.payload(tmp_path)["transcripts"]["ENST00000303395"]
        assert entry["silent"] is False
        assert len(entry["top"]) == 3

    def test_the_short_label_fits_a_picker_option(self, tmp_path: Path) -> None:
        # The long label does not: an option already carries an id, a marker and a biotype.
        transcripts = self.payload(tmp_path)["transcripts"]
        assert transcripts["ENST00000303395"]["short_label"] == "brain 1.63 TPM"
        assert transcripts["ENST00000674923"]["short_label"] == "brain 0.00 TPM"
        assert transcripts["ENST00000409050"]["short_label"] == "not in GTEx"
        assert all(len(t["short_label"]) <= 20 for t in transcripts.values())

    def test_the_note_and_the_caveat_complement_rather_than_repeat(
        self, tmp_path: Path
    ) -> None:
        # They first said the same thing about AF22 twice in one list. The caveat states the
        # fact, the note states what may and may not be concluded from it.
        payload = self.payload(tmp_path)
        fact = payload["caveats"][0]
        assert "no cell lines" in fact
        assert "no cell lines" not in payload["cell_line_note"]
        assert "supports or contradicts" in payload["cell_line_note"]

    def test_the_interpretation_is_listed_after_the_fact_it_interprets(
        self, tmp_path: Path
    ) -> None:
        # Rendered from two separate fields, the panel put "not assayed here means exactly
        # that" above the line saying GTEx holds no cell lines.
        listed = self.payload(tmp_path)["caveat_list"]
        fact = next(i for i, c in enumerate(listed) if "no cell lines" in c)
        note = next(i for i, c in enumerate(listed) if "supports or contradicts" in c)
        assert note == fact + 1

    def test_the_caveat_list_holds_every_caveat_plus_the_note(self, tmp_path: Path) -> None:
        payload = self.payload(tmp_path)
        assert len(payload["caveat_list"]) == len(payload["caveats"]) + 1
        assert set(payload["caveats"]) <= set(payload["caveat_list"])

    def test_a_caveat_says_bulk_tissue_averages_over_cell_types(
        self, tmp_path: Path
    ) -> None:
        caveats = " ".join(self.payload(tmp_path)["caveats"])
        assert "average of many cell types" in caveats
        assert "can read as near zero" in caveats


class TestThePresentationNeverRanks:
    def test_the_overview_carries_no_ordering_of_its_own(self, tmp_path: Path) -> None:
        # A list would be an order, and an order next to a picker is a recommendation.
        router = AppRouter(service_with_gtex(tmp_path))
        payload = json.loads(router.dispatch("GET", "/api/genes/SCN1A/expression")[2])
        assert isinstance(payload["transcripts"], dict)

    def test_unavailable_carries_no_transcript_records_at_all(self) -> None:
        # Empty records would render as a column of blanks, which looks like measured absence.
        assert unavailable("down", "gtex-unavailable")["transcripts"] == {}

    def test_the_brain_note_says_the_figure_is_not_an_average(self) -> None:
        profile = gtex_for().expression_for("ENSG00000144285")
        overview = expression_overview(profile)
        assert "not an average" in overview["brain_note"]
        assert "cervical spinal cord" in overview["brain_note"]

    def test_the_overview_keeps_the_status_note_on_every_transcript(self) -> None:
        profile = gtex_for().expression_for("ENSG00000144285", ["ENST00000713692"])
        overview = expression_overview(profile)
        absent = overview["transcripts"]["ENST00000713692"]
        assert absent["status"] == str(ExpressionStatus.NOT_IN_RELEASE)
        assert "not a measurement of zero" in absent["status_note"]
