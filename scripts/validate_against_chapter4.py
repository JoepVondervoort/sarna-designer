"""Run every preset against the measured activity data and report what it got right.

The first thing in this project that can show a design rule to be wrong. Run it with::

    .venv/bin/python scripts/validate_against_chapter4.py

Read the output in this order:

1. **Activators the gate rejected.** The one unambiguous finding this dataset supports. A duplex
   that demonstrably worked and that a preset would have discarded is evidence against a
   mandatory rule, whatever the ranking does.
2. **Per-rule pass rates.** Which rules separate the classes, and in which direction. A rule
   passing *fewer* activators than non-activators is pointing the wrong way.
3. **AUC.** Whether the ranking as a whole beats a coin flip. Expect it not to: see the caveats.

Written up in docs/scientific/empirical-validation.md.
"""

from sarna.design.presets import PRESETS
from sarna.validation.benchmark import benchmark
from sarna.validation.dataset import load_chapter4_dataset

TOP_RULES = 6


def main() -> None:
    """Benchmark every shipped preset and print a comparison."""
    dataset = load_chapter4_dataset()
    print(f"Dataset: {dataset.name}")
    print(f"Source:  {dataset.provenance}")
    print(f"Groups:  {len(dataset.activating)} activating, {len(dataset.not_activating)} not\n")
    print("Caveats that bound every number below:")
    for caveat in dataset.caveats:
        print(f"  - {caveat}")

    results = {}
    for preset_id, preset in PRESETS.items():
        result = benchmark(dataset, preset.resolve())
        results[preset_id] = result
        print(f"\n{'=' * 78}\n{preset.name}\n{'=' * 78}")
        print(
            f"  gate rejected {result.rejected_count} of {len(dataset)} sequences, "
            f"including {len(result.rejected_activators)} of "
            f"{len(dataset.activating)} known activators"
        )
        for scored in result.rejected_activators:
            reasons = ", ".join(r.rule_id for r in scored.candidate.excluded_by)
            print(f"      ! {scored.sequence.name} ({scored.sequence.gene}) excluded by {reasons}")
        print(
            f"  mean preference points: {result.mean_points_activating:.2f} among activators "
            f"vs {result.mean_points_not_activating:.2f} among the rest, "
            f"out of {result.max_preference_points}"
        )
        print(f"  AUC {result.auc:.3f}  (p = {result.auc_p_value:.3f})")
        print(f"  {result.summary}")
        print("  rules by separation:")
        for rule in result.rule_discrimination[:TOP_RULES]:
            print(f"      {rule}")

    print(f"\n{'=' * 78}\nSummary\n{'=' * 78}")
    print(f"  {'preset':<16} {'AUC':>6} {'p':>6} {'rejected':>9} {'activators lost':>16}")
    for preset_id, result in results.items():
        print(
            f"  {preset_id:<16} {result.auc:>6.3f} {result.auc_p_value:>6.3f} "
            f"{result.rejected_count:>9} {len(result.rejected_activators):>16}"
        )
    print(
        "\nNo preset's ranking separates activators from non-activators on this set, which is\n"
        "the expected result: every sequence in it was built to pass the gate, so the set holds\n"
        "almost none of the variation the ranking exists to order. What it does establish is\n"
        "which individual gates exclude duplexes that worked. See\n"
        "docs/scientific/empirical-validation.md."
    )


if __name__ == "__main__":
    main()
