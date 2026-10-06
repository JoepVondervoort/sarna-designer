"""The mandatory criteria from Wang et al. 2015, Table 1.

Failing any of these excludes a candidate, and the UI names the rule that excluded it. They are
not scored: an admission requirement is not a preference you can trade off (ADR 0006).

The ``consequence`` text is written for a researcher deciding whether to change a setting, not
for a developer reading the code.
"""

from sarna.rules.parameters import RuleParameter
from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.rule import Rule, RuleKind

TABLE_1 = Provenance.citation("Wang et al. 2015", "Table 1")


def _gate(
    rule_id: str,
    name: str,
    evidence_class: EvidenceClass,
    rationale: str,
    consequence: str,
    parameters: tuple[RuleParameter, ...] = (),
    provenance: Provenance = TABLE_1,
) -> Rule:
    return Rule(
        id=rule_id,
        name=name,
        kind=RuleKind.MANDATORY,
        evidence_class=evidence_class,
        provenance=provenance,
        rationale=rationale,
        consequence=consequence,
        parameters=parameters,
    )


TARGET_LENGTH = _gate(
    "target-length-19",
    "Target length is 19 nucleotides",
    EvidenceClass.ESTABLISHED,
    rationale=(
        "A 19 nt core carrying a 3' UU overhang on both strands is the standard siRNA/saRNA "
        "duplex geometry that RISC loads efficiently. Shorter cores lose specificity; longer "
        "ones are processed inconsistently."
    ),
    consequence=(
        "This defines the candidate window itself, so changing it changes every candidate. "
        "It also moves the goalposts for every position-based rule: 'position 18' and the "
        "20-23 flank only mean what they mean because the target is 19 long."
    ),
    parameters=(
        RuleParameter(
            name="target-length",
            default=19,
            minimum=18,
            maximum=21,
            unit="nt",
            description="Length of the targeted promoter subsequence.",
        ),
    ),
)

GC_CONTENT = _gate(
    "gc-content",
    "GC content between 40 % and 60 %",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "GC content sets duplex stability. Below the band the duplex is too weak and the "
        "region is often repetitive; above it the duplex resists unwinding and RISC loading, "
        "and the target increasingly sits in the CpG-rich promoter core."
    ),
    consequence=(
        "Widening the band admits more candidates, but the ones gained at the top end tend to "
        "fall inside CpG islands. The original script allowed up to 65 %, and that single "
        "difference accounts for much of the gap between the two presets."
    ),
    parameters=(
        RuleParameter(
            name="minimum-gc-fraction",
            default=0.40,
            minimum=0.0,
            maximum=1.0,
            unit="fraction",
            description="Lowest GC fraction a target may have.",
        ),
        RuleParameter(
            name="maximum-gc-fraction",
            default=0.60,
            minimum=0.0,
            maximum=1.0,
            unit="fraction",
            description="Highest GC fraction a target may have.",
        ),
    ),
)

NO_HOMOPOLYMER_RUN = _gate(
    "no-homopolymer-run",
    "No run of five or more identical bases",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "Long single-base runs are synthesised inaccurately, encourage off-target "
        "hybridisation, and can form secondary structures that block strand loading."
    ),
    consequence=(
        "Allowing longer runs admits candidates that frequently fail at synthesis or behave "
        "unpredictably. Tightening it below four rejects a great many candidates for little "
        "gain."
    ),
    parameters=(
        RuleParameter(
            name="maximum-run-length",
            default=4,
            minimum=3,
            maximum=8,
            unit="nt",
            description=(
                "Longest permitted run of one repeated base. The default of 4 rejects any run "
                "of 5 or more."
            ),
        ),
    ),
)

CPG_COUNT = _gate(
    "cpg-count",
    "At most one CpG within the target",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "CpG dinucleotides are the substrate for DNA methylation. A methylated target region "
        "can block the duplex from engaging, making activity depend on the cell's methylation "
        "state rather than on the design."
    ),
    consequence=(
        "Allowing more CpGs opens up CpG-rich promoters at the cost of methylation-dependent "
        "behaviour. The original script tolerated none at all, which is stricter than the "
        "paper and discards otherwise sound candidates."
    ),
    parameters=(
        RuleParameter(
            name="maximum-cpg-count",
            default=1,
            minimum=0,
            maximum=5,
            unit="count",
            description="How many CpG dinucleotides the target may contain.",
        ),
    ),
)

NOT_IN_CPG_ISLAND = _gate(
    "not-in-cpg-island",
    "Target does not lie inside a CpG island",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "CpG islands are densely methylatable and typically the most tightly regulated part "
        "of a promoter. Targets inside one are harder to access and vary more between cell "
        "types. Island boundaries follow the Takai & Jones definition."
    ),
    consequence=(
        "Switching this off admits candidates in the CpG-rich promoter core, which is exactly "
        "where many strong-looking candidates sit. In the original script this rule was "
        "written but commented out, so it never ran."
    ),
    parameters=(
        RuleParameter(
            name="island-minimum-length",
            default=200,
            minimum=100,
            maximum=2000,
            unit="bp",
            description="Shortest window that can qualify as a CpG island.",
        ),
        RuleParameter(
            name="island-minimum-gc-fraction",
            default=0.50,
            minimum=0.0,
            maximum=1.0,
            unit="fraction",
            description="Lowest GC fraction for a window to qualify as a CpG island.",
        ),
        RuleParameter(
            name="island-minimum-observed-expected-cpg",
            default=0.60,
            minimum=0.0,
            maximum=2.0,
            unit="ratio",
            description=(
                "Lowest observed-to-expected CpG ratio for a window to qualify as an island."
            ),
        ),
    ),
)

NO_SIMPLE_REPEATS = _gate(
    "no-simple-repeats",
    "No di- or tri-nucleotide simple repeats",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "Simple repeats recur in many copies genome-wide, so a duplex targeting one is likely "
        "to engage unintended promoters. They also form unstable secondary structures."
    ),
    consequence=(
        "This rule is absent from the original script, so the legacy preset returns "
        "repeat-containing candidates that this one rejects."
    ),
    parameters=(
        RuleParameter(
            name="maximum-tandem-repeat-units",
            default=3,
            minimum=2,
            maximum=8,
            unit="count",
            description=(
                "How many tandem copies of a di- or tri-nucleotide are tolerated before the "
                "target is treated as a simple repeat."
            ),
        ),
    ),
)

TSS_WINDOW = _gate(
    "tss-window",
    "Target lies 100-1000 bp upstream of the TSS",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "Activating duplexes in Wang et al. cluster in the promoter-proximal region upstream "
        "of the transcription start site. Further upstream the sequence is less likely to be "
        "the functional promoter; downstream the duplex overlaps the transcript itself and "
        "tends to silence rather than activate."
    ),
    consequence=(
        "The original script has no concept of a TSS at all and scans whatever sequence is "
        "pasted in, so this is the largest single difference between the presets on a "
        "genome-based run. Widening the window admits candidates with weaker evidence of "
        "promoter activity."
    ),
    parameters=(
        RuleParameter(
            name="window-start-bp",
            default=-1000,
            minimum=-10_000,
            maximum=0,
            unit="bp",
            description="Furthest upstream edge of the permitted window, relative to the TSS.",
        ),
        RuleParameter(
            name="window-end-bp",
            default=-100,
            minimum=-10_000,
            maximum=0,
            unit="bp",
            description="Nearest edge of the permitted window, relative to the TSS.",
        ),
    ),
)

POS18_AT = _gate(
    "pos18-at",
    "Position 18 is A or T",
    EvidenceClass.EMPIRICAL,
    rationale=(
        "An A or T at target position 18 weakens the 3' end of the target, which is the 5' "
        "end of the guide strand. A less stable guide 5' end is preferentially loaded into "
        "RISC, so this biases loading towards the strand that should do the work."
    ),
    consequence=(
        "Relaxing this admits candidates whose passenger strand may be loaded instead, which "
        "dilutes activity and can cause off-target effects from the wrong strand."
    ),
)

POS19_A = _gate(
    "pos19-a",
    "Position 19 is A",
    EvidenceClass.ESTABLISHED,
    rationale=(
        "An A at target position 19 places a U at the 5' end of the guide strand. AGO2 "
        "strongly prefers U there (Ui-Tei et al. 2004), so this is one of the most reliable "
        "levers on which strand gets loaded."
    ),
    consequence=(
        "This is the most restrictive rule in the set, cutting roughly three quarters of "
        "candidates on its own. Turning it off greatly increases the candidate count at a "
        "real cost in strand-selection confidence."
    ),
)

THERMODYNAMIC_ASYMMETRY = _gate(
    "thermodynamic-asymmetry",
    "Guide 5' end is less stable than the passenger 5' end",
    EvidenceClass.ESTABLISHED,
    rationale=(
        "RISC loads whichever strand has the more weakly paired 5' end (Schwarz et al. 2003; "
        "Khvorova et al. 2003). The difference is computed with a nearest-neighbour "
        "thermodynamic model over the terminal base pairs of the duplex."
    ),
    consequence=(
        "Raising the threshold above zero demands a stronger asymmetry, giving fewer "
        "candidates whose guide strand is more confidently the one loaded."
    ),
    parameters=(
        RuleParameter(
            name="minimum-asymmetry-kcal",
            default=0.0,
            minimum=-5.0,
            maximum=5.0,
            unit="kcal/mol",
            description=(
                "Smallest acceptable stability difference between the two duplex ends. Zero "
                "requires only that the guide 5' end is not the more stable of the two."
            ),
        ),
    ),
)

MANDATORY_RULES: tuple[Rule, ...] = (
    TARGET_LENGTH,
    GC_CONTENT,
    NO_HOMOPOLYMER_RUN,
    CPG_COUNT,
    NOT_IN_CPG_ISLAND,
    NO_SIMPLE_REPEATS,
    TSS_WINDOW,
    POS18_AT,
    POS19_A,
    THERMODYNAMIC_ASYMMETRY,
)
