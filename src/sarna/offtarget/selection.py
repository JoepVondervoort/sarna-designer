"""Conservative project selection policy, separate from experimental evidence and design rank.

Neither the 18/19 exclusion nor the five 17/19 alignment limit is a validated saRNA
specificity threshold. Counts concern reference alignments, not distinct genes or effects.
Both search scopes must pass this policy before the UI offers a screened shortlist.
"""

from collections.abc import Mapping
from typing import Any

POLICY_ID = "conservative-full-length-v1"


def assess_counts(counts: Mapping[str, int], *, own_locus_known: bool) -> dict[str, Any]:
    """Assess complete 19-base enumeration counts; never substitute missing evidence with zero."""
    policy = {
        "id": POLICY_ID,
        "evidence": "Project selection rule; not an experimentally validated safety threshold",
        "exclude_at_identity": 18,
        "exclude_at_17_base_count": 5,
        "count_unit": "reference alignments per scope, including isoforms and alternate records",
    }
    if not own_locus_known:
        return {"status": "pending", "reasons": ["Intended genomic locus is not known."],
                "policy": policy}
    if any(str(n) not in counts for n in range(16, 20)) or any(
        not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in counts.values()
    ):
        return {"status": "pending", "reasons": ["Complete identity counts are unavailable."],
                "policy": policy}
    reasons = []
    close = counts["18"] + counts["19"]
    if close:
        reasons.append(f"{close} full-length match(es) at 18/19 or 19/19 bases.")
    if counts["17"] >= 5:
        reasons.append(f'{counts["17"]} matches at 17/19 bases; project limit is fewer than 5.')
    return {"status": "excluded" if reasons else "eligible", "reasons": reasons,
            "policy": policy}
