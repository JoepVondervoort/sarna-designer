"""Which candidates would be ordered as the same molecule.

A candidate is a *locus*: a 19 nt target site at one position in the promoter. What gets
synthesised is determined by the target's bases alone, so the same 19-mer occurring at several
positions produces several candidates and exactly one product. On a tandem-repeat fixture every
one of the eight passing candidates carried the same target DNA, and "Order CSV, all 8" would
have bought the same two oligos eight times.

They are deliberately not merged. The position-dependent rules, the 3' flank preferences and the
TSS window, score the same duplex differently at different loci: that same fixture held a 6/7 at
one position and a 5/7 at another. Collapsing them would throw that away and would decide, on the
researcher's behalf, which locus the molecule "really" targets. Naming the twins reports the fact
and leaves the choice where it belongs.
"""

from collections.abc import Iterable


def identical_duplex_groups(
    candidates: Iterable[tuple[str, str]],
) -> dict[str, tuple[str, ...]]:
    """Map each candidate id to the other ids ordering as the same duplex.

    Args:
        candidates: ``(candidate id, target site DNA)`` pairs. The target's bases are the whole
            key, because the ordered strands are derived from them and nothing else.

    Returns:
        Only the candidates that have at least one twin, each mapped to its twins in the order
        the candidates arrived. A candidate with a unique duplex is absent rather than present
        with an empty tuple, so a caller reading ``groups.get(id, ())`` cannot mistake "unique"
        for "not looked at".
    """
    by_target: dict[str, list[str]] = {}
    for candidate_id, target in candidates:
        by_target.setdefault(target, []).append(candidate_id)
    return {
        candidate_id: tuple(other for other in group if other != candidate_id)
        for group in by_target.values()
        if len(group) > 1
        for candidate_id in group
    }
