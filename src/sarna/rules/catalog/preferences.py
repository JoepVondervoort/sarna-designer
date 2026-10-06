"""Scored preferences, the tie-break, and the one rule carried without a source.

Each preference is worth exactly one point and all are weighted equally (ADR 0006), because
Wang et al. publish no weights and inventing some is what made the original score
uninterpretable. Table 1 marks two groups with the word "preferably"; both are scored
per-position, so a candidate meeting three of the four flank positions ranks above one meeting
none instead of being collapsed to a single yes/no.
"""

from sarna.rules.catalog.mandatory import TABLE_1
from sarna.rules.parameters import RuleParameter
from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.rule import Rule, RuleKind

WANG_NARRATIVE = Provenance.citation(
    "Wang et al. 2015", "Results — preferred -500 to -200 bp region (narrative text, not Table 1)"
)

START_GROUP = "Target 5' start G/C"
FLANK_GROUP = "3' flank A/U richness"

_START_RATIONALE = (
    "Wang et al. list a G or C at the first two target positions as preferable. A stable 5' "
    "end on the target is a stable 3' end on the guide, which reinforces the same asymmetry "
    "that selects the guide strand for loading."
)
_START_CONSEQUENCE = (
    "Worth one point. Candidates without it are ranked lower, never excluded, so switching it "
    "off mostly reshuffles the top of the list."
)

_FLANK_RATIONALE = (
    "The four bases immediately 3' of the target are preferably A or T. An A/T-rich flank is "
    "easier to unwind, which helps the duplex engage the target site."
)
_FLANK_CONSEQUENCE = (
    "Worth one point per position, so this group contributes 0 to 4 of the 7 available points. "
    "The original script scored this region too, but weighted it at 1 against 1000 for other "
    "terms, so in practice it almost never changed the ranking."
)


def _preference(
    rule_id: str,
    name: str,
    rationale: str,
    consequence: str,
    group: str,
    evidence_class: EvidenceClass = EvidenceClass.EMPIRICAL,
    provenance: Provenance = TABLE_1,
    parameters: tuple[RuleParameter, ...] = (),
) -> Rule:
    return Rule(
        id=rule_id,
        name=name,
        kind=RuleKind.PREFERENCE,
        evidence_class=evidence_class,
        provenance=provenance,
        rationale=rationale,
        consequence=consequence,
        group=group,
        parameters=parameters,
    )


START_PREFERENCES: tuple[Rule, ...] = tuple(
    _preference(
        f"start-{position}-gc",
        f"Position {position} is G or C",
        _START_RATIONALE,
        _START_CONSEQUENCE,
        START_GROUP,
    )
    for position in (1, 2)
)

FLANK_PREFERENCES: tuple[Rule, ...] = tuple(
    _preference(
        f"flank-{position}-au",
        f"Position {position} is A or T",
        _FLANK_RATIONALE,
        _FLANK_CONSEQUENCE,
        FLANK_GROUP,
    )
    for position in (20, 21, 22, 23)
)

TSS_PREFERRED_WINDOW = _preference(
    "tss-preferred-window",
    "Target lies 200-500 bp upstream of the TSS",
    rationale=(
        "Within the permitted -1000 to -100 bp window, Wang et al. report activity "
        "concentrated between roughly -500 and -200 bp. Targets there are more likely to sit "
        "over the regulatory elements the duplex is meant to act on."
    ),
    consequence=(
        "Worth one point, and it sits inside the mandatory window, so it narrows the ranking "
        "rather than admission. Note this comes from the paper's text rather than Table 1, so "
        "it is on weaker footing than the criteria beside it."
    ),
    group="TSS proximity",
    provenance=WANG_NARRATIVE,
    parameters=(
        RuleParameter(
            name="preferred-window-start-bp",
            default=-500,
            minimum=-10_000,
            maximum=0,
            unit="bp",
            description="Furthest upstream edge of the preferred region, relative to the TSS.",
        ),
        RuleParameter(
            name="preferred-window-end-bp",
            default=-200,
            minimum=-10_000,
            maximum=0,
            unit="bp",
            description="Nearest edge of the preferred region, relative to the TSS.",
        ),
    ),
)

POS7_U = Rule(
    id="pos7-u",
    name="Position 7 is U",
    kind=RuleKind.PREFERENCE,
    evidence_class=EvidenceClass.EXPLORATORY,
    provenance=Provenance.unsourced(
        "Position 7 appears nowhere in Wang et al. 2015, in Table 1 or elsewhere. It may be a "
        "half-remembered Reynolds et al. 2004 siRNA criterion, which concerns position 10 of "
        "the guide strand, but no source has been produced."
    ),
    rationale=(
        "No mechanism is established for this rule. It is carried because the original script "
        "weighted it at 1000, making it one of the two dominant ranking terms, and deleting it "
        "outright would leave historical rankings unexplainable."
    ),
    consequence=(
        "Off by default. Enabling it adds an eighth preference worth one point, like every "
        "other preference — not the dominant term it was in the original script. Treat any "
        "ranking change it causes as untested."
    ),
    enabled_by_default=False,
    group="Exploratory",
)

ASYMMETRY_MARGIN = Rule(
    id="asymmetry-margin",
    name="Thermodynamic asymmetry margin",
    kind=RuleKind.TIEBREAK,
    evidence_class=EvidenceClass.ESTABLISHED,
    provenance=Provenance.citation(
        "Schwarz et al. 2003", "Cell 115(2):199-208, asymmetric RISC assembly"
    ),
    rationale=(
        "Where candidates tie on preference count, the one with the larger stability "
        "difference between the two duplex ends has the more confident strand selection. It "
        "is continuous, so it discriminates where a checklist cannot."
    ),
    consequence=(
        "Affects ordering only. It never changes the score, never excludes a candidate, and "
        "is reported as a tie-break rather than folded into the number."
    ),
)

PREFERENCE_RULES: tuple[Rule, ...] = (
    *START_PREFERENCES,
    *FLANK_PREFERENCES,
    TSS_PREFERRED_WINDOW,
)

TIEBREAK_RULES: tuple[Rule, ...] = (ASYMMETRY_MARGIN,)

EXPLORATORY_RULES: tuple[Rule, ...] = (POS7_U,)
