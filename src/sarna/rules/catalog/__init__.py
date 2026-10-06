"""The rule catalogue: every design rule this project knows about.

The catalogue is the full set, not a configuration. A preset is a named *selection* from it,
so "which rules exist" and "which rules this run used" stay separately inspectable.

Order here is presentation order: the mandatory gates in Wang et al. 2015 Table 1 order, then
the scored preferences, then the tie-break, then anything exploratory and off by default.
"""

from sarna.rules.catalog.mandatory import MANDATORY_RULES
from sarna.rules.catalog.preferences import (
    EXPLORATORY_RULES,
    PREFERENCE_RULES,
    TIEBREAK_RULES,
)
from sarna.rules.registry import RuleRegistry

CATALOGUE = RuleRegistry(
    [
        *MANDATORY_RULES,
        *PREFERENCE_RULES,
        *TIEBREAK_RULES,
        *EXPLORATORY_RULES,
    ]
)

__all__ = ["CATALOGUE"]
