"""Optional criteria from the siRNA design literature and from thesis chapter 4.

Every rule here is **off by default**. The default preset stays exactly conformant with Wang et
al. 2015, and nothing in this module changes which candidates it admits unless a user selects it.

Two bodies of work motivate them, and they carry different weight.

**Established siRNA design criteria.** Reynolds et al. 2004 and Ui-Tei et al. 2004 derived
sequence criteria from large siRNA panels (180 and 62 duplexes respectively) and these are among
the most replicated rules in the field. The caveat that governs all of them: they were derived for
**knockdown**, where the loaded guide finds a cytoplasmic mRNA and directs its cleavage. Promoter-
targeted activation shares the first half of that path (duplex recognition, strand selection, AGO2
loading) and none of the second (nuclear import, promoter engagement, no cleavage). So a criterion
whose mechanism is *strand selection* has a reasonable claim to transfer, and one whose mechanism
is *target-site accessibility or cleavage geometry* does not. Each rule below says which it is.

**Thesis chapter 4's positional analysis.** Observational, n = 14 activators, and the enrichments
do not survive correction for the grid they were found in, a random relabelling of the same 125
sequences produces more positional structure than the activity labels do
(`docs/scientific/positional-motif-analysis.md`). They are offered so a hypothesis can be stated
and tested deliberately, never because the project believes them (ADR 0011).

One convergence is worth more than either source alone. Reynolds criterion 2 asks for at least
three A/U among target positions 15-19; Ui-Tei criterion iii asks for at least five A/U in the 5'
terminal third of the guide, which is the same stretch of the same molecule; and chapter 4's
stability analysis independently found weak bases enriched at target positions 15-17. Three
analyses, two of them well powered, pointing at one region. That is the strongest optional rule
here, and it is :data:`AU_RICH_GUIDE_5_PRIME`.
"""

from sarna.rules.options import RuleOption
from sarna.rules.parameters import RuleParameter
from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.rule import Rule, RuleKind

REYNOLDS_2004 = Provenance.citation(
    "Reynolds et al. 2004",
    "Nat Biotechnol 22(3):326-330, rational siRNA design, criteria 1-8 (sense strand)",
)

UI_TEI_2004 = Provenance.citation(
    "Ui-Tei et al. 2004",
    "Nucleic Acids Res 32(3):936-948, four sequence conditions for effective RNAi",
)

THESIS_CHAPTER_4 = Provenance.citation(
    "van de Vondervoort 2026",
    "chapter 4, positional and stability motif analysis (Figures 4-5); significance "
    "reassessed in docs/scientific/positional-motif-analysis.md",
)

SIRNA_GROUP = "siRNA design criteria (knockdown-derived)"
THESIS_GROUP = "Chapter 4 positional hypotheses"

_TRANSFER_CAVEAT = (
    "Derived from siRNA knockdown panels, not from promoter-targeted activation. "
)

_THESIS_CAVEAT = (
    "Observational, from 14 activators, and it does not survive correction for the 76-cell "
    "grid it was found in: a random relabelling of the same sequences produces a larger "
    "positional difference than the real labels do. Selecting it states a hypothesis. "
)


def _optional(
    rule_id: str,
    name: str,
    provenance: Provenance,
    evidence_class: EvidenceClass,
    rationale: str,
    consequence: str,
    group: str,
    parameters: tuple[RuleParameter, ...] = (),
    options: tuple[RuleOption, ...] = (),
) -> Rule:
    """Build an off-by-default preference."""
    return Rule(
        id=rule_id,
        name=name,
        kind=RuleKind.PREFERENCE,
        evidence_class=evidence_class,
        provenance=provenance,
        rationale=rationale,
        consequence=consequence,
        enabled_by_default=False,
        group=group,
        parameters=parameters,
        options=options,
    )


AU_RICH_GUIDE_5_PRIME = _optional(
    "au-rich-guide-5-prime",
    "A/T-rich 3' end of target (weak 5' end of guide)",
    UI_TEI_2004,
    EvidenceClass.ESTABLISHED,
    rationale=(
        "The strand whose 5' end is thermodynamically weaker is the strand AGO2 retains as the "
        "guide (Schwarz et al. 2003; Khvorova et al. 2003). The guide's 5' end corresponds to "
        "the target's 3' end, so A/T there biases loading toward the intended strand. Ui-Tei et "
        "al. require at least five A/U in the 5' terminal third of the antisense strand; "
        "Reynolds et al. independently require at least three A/U at sense positions 15-19, the "
        "same stretch. Chapter 4's stability analysis found weak bases enriched at target "
        "positions 15-17 without looking for either criterion."
    ),
    consequence=(
        "Worth one point. The best-supported rule in this module: its mechanism is strand "
        "selection, which activation and knockdown share, and three independent analyses agree "
        "on the region. It overlaps the mandatory asymmetry rule, which measures the same "
        "physics continuously, so expect it to reorder the list modestly rather than "
        "dramatically. Wang's position 18 and 19 rules already constrain two of these seven "
        "positions."
    ),
    group=SIRNA_GROUP,
    parameters=(
        RuleParameter(
            name="window-length-nt",
            default=7,
            minimum=3,
            maximum=10,
            unit="nt",
            description=(
                "How many 3'-terminal target positions to examine. Seven is the 5' terminal "
                "third of a 19 nt guide, which is what Ui-Tei et al. specify."
            ),
        ),
        RuleParameter(
            name="minimum-au-count",
            default=5,
            minimum=1,
            maximum=10,
            unit="nt",
            description=(
                "How many of those positions must be A or T. Five is Ui-Tei et al.'s "
                "threshold; three over positions 15-19 is Reynolds et al.'s, which is the "
                "more permissive of the two."
            ),
        ),
    ),
    options=(
        RuleOption(
            id="ui-tei-2004",
            label="5 of the last 7 (Ui-Tei et al. 2004)",
            provenance=UI_TEI_2004,
            values={"window-length-nt": 7, "minimum-au-count": 5},
            note=(
                "The published criterion: at least five A/U in the 5' terminal third of the "
                "guide. Derived from siRNA knockdown, and on the one saRNA dataset available it "
                "is not met by 7 of 14 measured activators, so at this threshold the rule costs "
                "half of them a point. Kept as the default because it is what the source says."
            ),
            is_default=True,
        ),
        RuleOption(
            id="reynolds-2004",
            label="3 of the last 5 (Reynolds et al. 2004)",
            provenance=REYNOLDS_2004,
            values={"window-length-nt": 5, "minimum-au-count": 3},
            note=(
                "Reynolds et al.'s form of the same criterion, over positions 15-19. Markedly "
                "better calibrated for activation: only 1 of the 14 measured activators misses "
                "it, against 7 under Ui-Tei's threshold. Prefer this one if you intend the rule "
                "to reward activators rather than to restate a knockdown criterion."
            ),
        ),
    ),
)

NO_GC_STRETCH = _optional(
    "no-gc-stretch",
    "No GC stretch longer than nine bases",
    UI_TEI_2004,
    EvidenceClass.ESTABLISHED,
    rationale=(
        "A long uninterrupted G/C tract is thermodynamically stable enough to resist unwinding "
        "and to favour intramolecular structure, either of which keeps the duplex from being "
        "loaded. Ui-Tei et al. list its absence as one of four conditions an effective duplex "
        "must satisfy simultaneously."
    ),
    consequence=(
        "Worth one point, and it almost never fires inside the default 40-60 % GC band: a "
        "ten-base G/C tract in a 19-mer forces GC content above 52 %. Mostly useful when the "
        "GC band has been widened, or with the legacy preset's 65 % ceiling. Distinct from "
        "`no-homopolymer-run`, which counts repeats of one base rather than G/C of either kind."
    ),
    group=SIRNA_GROUP,
    parameters=(
        RuleParameter(
            name="maximum-gc-stretch-nt",
            default=9,
            minimum=4,
            maximum=19,
            unit="nt",
            description=(
                "Longest permitted run of G and C in any combination. Nine is Ui-Tei et al.'s "
                "figure; the rule fails a candidate exceeding it."
            ),
        ),
    ),
)

POS3_A = _optional(
    "pos3-a",
    "Position 3 is A",
    REYNOLDS_2004,
    EvidenceClass.EXPLORATORY,
    rationale=(
        "Reynolds et al. found A at sense position 3 associated with efficient knockdown "
        "across 180 duplexes. No mechanism was established for it, and the position is in the "
        "guide's 3' supplementary region rather than anywhere strand selection acts."
    ),
    consequence=(
        _TRANSFER_CAVEAT + "With no mechanism to argue from, there is no reason to expect it "
        "to carry over to activation, and chapter 4's data gives it no support: position 3 "
        "shows C depleted, not A enriched. Worth one point if you want Reynolds' full set."
    ),
    group=SIRNA_GROUP,
)

POS10_U = _optional(
    "pos10-u",
    "Position 10 is T (U on the duplex)",
    REYNOLDS_2004,
    EvidenceClass.EXPLORATORY,
    rationale=(
        "Position 10 of the guide is where AGO2 cleaves its target, and Reynolds et al. found "
        "U favoured there. The mechanism is cleavage geometry."
    ),
    consequence=(
        _TRANSFER_CAVEAT + "This is the clearest case against transfer: activation does not "
        "proceed by cleavage, so the reason the criterion exists in knockdown does not apply. "
        "Chapter 4's data points the other way, finding G at position 10 its single largest "
        "enrichment (+24.8 pp), see `pos10-g`, which contradicts this rule. Selecting both "
        "guarantees every candidate loses one of the two points."
    ),
    group=SIRNA_GROUP,
)

NO_G_AT_13 = _optional(
    "no-g-at-13",
    "Position 13 is not G",
    REYNOLDS_2004,
    EvidenceClass.EXPLORATORY,
    rationale=(
        "Reynolds et al. applied a penalty for G at sense position 13. As with position 3, the "
        "association was empirical and no mechanism was offered."
    ),
    consequence=(
        _TRANSFER_CAVEAT + "Chapter 4's stability analysis points the opposite way here too, "
        "finding strong (G/C) bases enriched around positions 12-13. Worth one point; include "
        "it for completeness of Reynolds' set rather than for a reason specific to activation."
    ),
    group=SIRNA_GROUP,
)

POS10_G = _optional(
    "pos10-g",
    "Position 10 is G",
    THESIS_CHAPTER_4,
    EvidenceClass.EXPLORATORY,
    rationale=(
        "The largest single positional difference in chapter 4's screen: G at target position "
        "10 appeared in 24.8 percentage points more of the activating sequences than of the "
        "non-activating ones. Equivalent to about three of the fourteen activators."
    ),
    consequence=(
        _THESIS_CAVEAT + "It also directly contradicts `pos10-u` from Reynolds et al., so the "
        "two should not be selected together. If you want to test this properly, design "
        "duplexes that differ at position 10 and nowhere else: the grid then has one cell "
        "instead of 76, and the same data becomes decisive."
    ),
    group=THESIS_GROUP,
)

CENTRAL_STABILITY = _optional(
    "central-stability",
    "Strong (G/C) bases in the target's central region",
    THESIS_CHAPTER_4,
    EvidenceClass.EXPLORATORY,
    rationale=(
        "Chapter 4's stability analysis, which collapses bases to strong (G/C) and weak (A/U), "
        "found strong bases enriched at target positions 10 and 12 and SS pairs around 12-13, "
        "alongside weak bases toward position 15. Read together that is a stability gradient: "
        "a stable core with a weak 3' end, which is what asymmetric loading predicts."
    ),
    consequence=(
        _THESIS_CAVEAT + "The gradient's weak-3'-end half is already covered, better, by "
        "`au-rich-guide-5-prime` and by the mandatory asymmetry rule, so this rule contributes "
        "only the central-stability half and partly double-counts them. Raising the required "
        "count narrows the ranking sharply, because a GC-rich core pushes against the 60 % GC "
        "ceiling from the other direction."
    ),
    group=THESIS_GROUP,
    parameters=(
        RuleParameter(
            name="central-region-start",
            default=10,
            minimum=1,
            maximum=19,
            unit="position",
            description="First target position of the central region, 1-based inclusive.",
        ),
        RuleParameter(
            name="central-region-end",
            default=13,
            minimum=1,
            maximum=19,
            unit="position",
            description="Last target position of the central region, 1-based inclusive.",
        ),
        RuleParameter(
            name="minimum-strong-count",
            default=3,
            minimum=1,
            maximum=19,
            unit="nt",
            description=(
                "How many positions in the region must be G or C. Three of four matches the "
                "pattern chapter 4 describes without demanding all of it."
            ),
        ),
    ),
)


OPTIONAL_RULES: tuple[Rule, ...] = (
    AU_RICH_GUIDE_5_PRIME,
    NO_GC_STRETCH,
    POS3_A,
    POS10_U,
    NO_G_AT_13,
    POS10_G,
    CENTRAL_STABILITY,
)

__all__ = [
    "AU_RICH_GUIDE_5_PRIME",
    "CENTRAL_STABILITY",
    "NO_GC_STRETCH",
    "NO_G_AT_13",
    "OPTIONAL_RULES",
    "POS3_A",
    "POS10_G",
    "POS10_U",
    "REYNOLDS_2004",
    "THESIS_CHAPTER_4",
    "UI_TEI_2004",
]
