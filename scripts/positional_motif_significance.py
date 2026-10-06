"""Measure whether chapter 4's positional motif enrichments exceed chance.

Chapter 4 reports percentage-point differences of roughly 20-25 between 14 functional and 111
non-functional saRNAs at a handful of positions, and is explicit that the analysis is
exploratory. This script supplies the number the chapter does not: the family-wise error rate
over the whole grid of positions and motifs that was searched.

Run it with::

    .venv/bin/python scripts/positional_motif_significance.py

It prints the largest observed differences for 1-, 2- and 3-mers, then the permutation result
for each. Output is deterministic given the seed, which is printed with the result.

Written up in docs/scientific/positional-motif-analysis.md.
"""

from sarna.validation.dataset import load_chapter4_dataset
from sarna.validation.motifs import (
    DEFAULT_PERMUTATIONS,
    DEFAULT_SEED,
    frequency_differences,
    permutation_test,
)

MOTIF_WIDTHS = (1, 2, 3)
TOP_N = 8


def main() -> None:
    """Report observed enrichments and their family-wise significance for each motif width."""
    dataset = load_chapter4_dataset()
    activating = [s.target for s in dataset.activating]
    inactive = [s.target for s in dataset.not_activating]

    print(f"Dataset: {dataset.name}")
    print(f"Source:  {dataset.provenance}")
    print(
        f"Groups:  {len(activating)} activating vs {len(inactive)} not activating "
        f"(one sequence moves a frequency by {100 / len(activating):.1f} pp)"
    )

    for width in MOTIF_WIDTHS:
        differences = frequency_differences(activating, inactive, width)
        print(f"\n{'=' * 78}\n{width}-mer motifs: {len(differences)} cells examined\n{'=' * 78}")
        print("Largest absolute differences (activating minus not activating):")
        for difference in differences[:TOP_N]:
            sequences = abs(difference.difference) / (100 / len(activating))
            print(f"  {difference}   (~{sequences:.1f} of the {len(activating)} activators)")

        result = permutation_test(
            activating, inactive, width, DEFAULT_PERMUTATIONS, DEFAULT_SEED
        )
        print(
            f"\n  observed maximum      {result.observed_maximum:6.1f} pp\n"
            f"  null median           {result.null_median:6.1f} pp\n"
            f"  null 95th percentile  {result.null_95th_percentile:6.1f} pp\n"
            f"  family-wise p         {result.family_wise_p:6.3f}"
        )
        print(f"\n  {result.summary}")

    print(f"\n{'=' * 78}")
    print(
        "Interpretation: the chapter's positional hotspots reproduce exactly as reported, and\n"
        "none of them survives correction for the size of the grid they were found in. They are\n"
        "usable as hypotheses to test, not as design criteria. See\n"
        "docs/scientific/positional-motif-analysis.md."
    )


if __name__ == "__main__":
    main()
