"""Measure how often the shipped parameter sets disagree about the guide strand.

ADR 0009 reports every selected table's verdict rather than one. This script quantifies why that
matters, by enumerating every distinguishable case exhaustively instead of sampling a promoter.

Only the terminal 4 nt at each end enter an end-stability sum, so the 11 bases between them are
held constant: all 4^4 x 4^4 = 65 536 end combinations are examined, which covers every input the
calculation can distinguish.

Run from the repository root:

    .venv/bin/python scripts/guide_prediction_divergence.py
"""

from collections import Counter
from itertools import product

from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.ends import DEFAULT_END_LENGTH
from sarna.thermo.guide import assess_asymmetry
from sarna.thermo.tables import DNA_SANTALUCIA_1998, LEGACY_DNA_KJ, RNA_XIA_1998

FILLER = "ACGGATCCAAG"
"""The 11 bases between the two ends. Never read by the calculation; present only for length."""

PAIRS = (
    ("RNA vs legacy script", RNA_XIA_1998, LEGACY_DNA_KJ),
    ("RNA vs DNA parameters", RNA_XIA_1998, DNA_SANTALUCIA_1998),
    ("DNA vs legacy script", DNA_SANTALUCIA_1998, LEGACY_DNA_KJ),
)


def main() -> None:
    """Enumerate every end combination and report each pair's disagreement rate."""
    disagreements: Counter[str] = Counter()
    examples: dict[str, str] = {}
    total = 0

    for head in product("ACGT", repeat=DEFAULT_END_LENGTH):
        for tail in product("ACGT", repeat=DEFAULT_END_LENGTH):
            target = NucleicSequence.parse("".join(head) + FILLER + "".join(tail))
            total += 1
            for label, first, second in PAIRS:
                one = assess_asymmetry(target, first).predicted_guide
                other = assess_asymmetry(target, second).predicted_guide
                if one is not other:
                    disagreements[label] += 1
                    examples.setdefault(label, f"{target.bases}: {one} vs {other}")

    print(f"End combinations examined: {total}\n")
    for label, _, _ in PAIRS:
        count = disagreements[label]
        print(f"{label:24s} {count:6d}  {100 * count / total:5.2f}%")
    print()
    for label, example in examples.items():
        print(f"First disagreement, {label}: {example}")


if __name__ == "__main__":
    main()
