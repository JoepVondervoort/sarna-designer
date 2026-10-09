"""The ranking's measured performance, as the interface reports it.

The number shown beside a rating is **computed from the bundled chapter 4 dataset**, not a
constant transcribed from a document. A transcribed figure drifts the moment a rule changes and
nobody re-runs the benchmark; a computed one cannot. These tests pin the values anyway, so that
a change to a published threshold breaks a test that names the measurement, which is the point
of measuring it.

Values come from `docs/scientific/empirical-validation.md`, measured over the same 125 assayed
sequences.
"""

import pytest

from sarna.validation.ranking_performance import (
    CHANCE_AUC,
    performance_for_preset,
    performance_summary,
)


class TestTheMeasuredNumbers:
    """Each preset's AUC, pinned. 0.5 is a coin flip."""

    @pytest.mark.parametrize(
        ("preset_id", "expected_auc"),
        [
            ("wang-2015", 0.4633),
            ("strict", 0.4633),
            ("exploratory", 0.5272),
            ("legacy-scanner", 0.4139),
        ],
    )
    def test_auc_matches_the_recorded_measurement(
        self, preset_id: str, expected_auc: float
    ) -> None:
        measured = performance_for_preset(preset_id)
        assert measured is not None
        assert measured.auc == pytest.approx(expected_auc, abs=5e-4)

    def test_no_preset_beats_chance(self) -> None:
        # The finding the interface exists to communicate. If this ever fails, the claim shown
        # to the user is wrong and the UI text must change with it.
        for preset_id in ("wang-2015", "strict", "exploratory", "legacy-scanner"):
            measured = performance_for_preset(preset_id)
            assert measured is not None
            assert not measured.discriminates

    def test_three_of_the_four_fall_below_chance(self) -> None:
        measured = {
            p: performance_for_preset(p)
            for p in ("wang-2015", "strict", "exploratory", "legacy-scanner")
        }
        assert all(m is not None for m in measured.values())
        below = [p for p, m in measured.items() if m is not None and m.auc < CHANCE_AUC]
        assert sorted(below) == ["legacy-scanner", "strict", "wang-2015"]


class TestWhatTheSummaryTellsTheUser:
    def test_the_summary_names_the_dataset_and_its_size(self) -> None:
        summary = performance_summary("wang-2015")
        assert summary is not None
        assert summary["dataset"]["assayed"] == 125
        assert summary["dataset"]["activating"] == 14
        assert "chapter 4" in summary["dataset"]["name"].lower()

    def test_the_summary_states_the_range_restriction(self) -> None:
        # Every sequence was built to pass Wang's gate, so the set can show a gate too strict
        # but never one too permissive. Reporting an AUC without that is misleading.
        summary = performance_summary("wang-2015")
        assert summary is not None
        assert "range-restricted" in summary["caveat"].lower()

    def test_the_summary_says_plainly_that_this_is_not_a_prediction(self) -> None:
        summary = performance_summary("wang-2015")
        assert summary is not None
        assert summary["beats_chance"] is False
        assert "chance" in summary["headline"].lower()

    def test_the_exploratory_preset_declares_its_circularity(self) -> None:
        # exploratory scores highest because two of its rules were fitted to this very data.
        # Showing its AUC without that would make the worst-evidenced preset look the best.
        summary = performance_summary("exploratory")
        assert summary is not None
        assert summary["circular"] is True
        assert "fitted" in summary["circularity_note"].lower()
        assert performance_summary("wang-2015")["circular"] is False  # type: ignore[index]


class TestUnknownAndCustomConfigurations:
    def test_an_unknown_preset_has_no_measurement_rather_than_a_default(self) -> None:
        # Returning 0.5, or silently falling back to Wang's number, would attach a measurement
        # to something that was never measured.
        assert performance_summary("not-a-preset") is None

    def test_the_result_is_cached_so_the_interface_can_ask_freely(self) -> None:
        first = performance_for_preset("wang-2015")
        assert performance_for_preset("wang-2015") is first
