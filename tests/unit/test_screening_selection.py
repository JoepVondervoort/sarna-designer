"""Project-rule boundaries must never be confused with measured biological specificity."""

import pytest

from sarna.offtarget.selection import assess_counts


@pytest.mark.parametrize(
    ("exact", "one_mismatch", "two_mismatches", "expected"),
    [(0, 0, 0, "eligible"), (0, 0, 4, "eligible"), (0, 0, 5, "excluded"),
     (1, 0, 0, "excluded"), (0, 1, 0, "excluded")],
)
def test_selection_rule_boundaries(exact: int, one_mismatch: int, two_mismatches: int,
                                   expected: str):
    counts = {"19": exact, "18": one_mismatch, "17": two_mismatches, "16": 10000}
    assessment = assess_counts(counts, own_locus_known=True)
    assert assessment["status"] == expected
    assert "not an experimentally validated" in assessment["policy"]["evidence"]


def test_unknown_own_locus_is_pending_even_when_exact_matches_exist():
    assert assess_counts({"19": 1, "18": 0, "17": 0, "16": 0},
                         own_locus_known=False)["status"] == "pending"


@pytest.mark.parametrize("counts", [{}, {"19": 0},
    {"19": 0, "18": 0, "17": -1, "16": 0},
    {"19": 0, "18": 0, "17": True, "16": 0}])
def test_missing_or_invalid_counts_are_not_zero_findings(counts: dict[str, int]):
    assert assess_counts(counts, own_locus_known=True)["status"] == "pending"
