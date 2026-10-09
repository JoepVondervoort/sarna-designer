"""The mandatory criteria from Wang et al. 2015, Table 1.

Failing any of these excludes a candidate, and the UI names the rule that excluded it. They are
not scored: an admission requirement is not a preference you can trade off (ADR 0006).

The ``consequence`` text is written for a researcher deciding whether to change a setting, not
for a developer reading the code.
"""

from sarna.rules.options import RuleOption
from sarna.rules.parameters import RuleParameter
from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.rule import Rule, RuleKind

TABLE_1 = Provenance.citation("Wang et al. 2015", "Table 1")


LEGACY_SCRIPT = Provenance.legacy_script("sarna_promoter_scanner.py")


def _gate(
    rule_id: str,
    name: str,
    evidence_class: EvidenceClass,
    rationale: str,
    consequence: str,
    parameters: tuple[RuleParameter, ...] = (),
    provenance: Provenance = TABLE_1,
    options: tuple[RuleOption, ...] = (),
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
        options=options,
    )


TARGET_LENGTH = _gate(
    "target-length-19",
    "Target length is 19 nucleotides",
    EvidenceClass.ESTABLISHED,
    rationale=(
        "A 19 nt core carrying a 2 nt 3' overhang on both strands is the standard siRNA/saRNA "
        "duplex geometry that RISC loads efficiently. Shorter cores lose specificity; longer "
        "ones are processed inconsistently. This project's own duplexes use dTdT overhangs "
        "rather than ribo-UU, which is a synthesis choice and does not change the 19 nt core."
    ),
    consequence=(
        "This defines the candidate window itself, so changing it changes every candidate. "
        "It also moves the goalposts for every position-based rule: 'position 18' and the "
        "20-23 flank only mean what they mean because the target is 19 long. At 18 nt the "
        "position-19 gate becomes not-applicable, so a run at that length silently loses the "
        "strongest sequence filter the engine has."
    ),
    parameters=(
        RuleParameter(
            name="target-length",
            default=19,
            minimum=18,
            # Capped at 19, not 21. The flank preferences are registered at literal positions
            # 20-23, so a 20 nt target collides with flank-20-au and the run raises. Declaring
            # a range that crashes is worse than declaring a narrow one, and 19 nt is the
            # settled chemistry for this project's duplexes rather than a convenience.
            maximum=19,
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
    options=(
        RuleOption(
            id="wang-2015",
            label="40-60 % (Wang et al. 2015)",
            provenance=TABLE_1,
            values={"minimum-gc-fraction": 0.40, "maximum-gc-fraction": 0.60},
            note=(
                "The band stated in the paper. Recommended unless you are reproducing an "
                "older design."
            ),
            is_default=True,
        ),
        RuleOption(
            id="legacy-scanner",
            label="40-65 % (original script)",
            provenance=LEGACY_SCRIPT,
            values={"minimum-gc-fraction": 0.40, "maximum-gc-fraction": 0.65},
            note=(
                "What the original script used. More permissive at the top end, so it admits "
                "GC-rich candidates the paper would exclude. Choose it to reproduce an "
                "existing design, not to make a new one."
            ),
        ),
        RuleOption(
            id="reynolds-2004",
            label="30-52 % (Reynolds et al. 2004)",
            provenance=Provenance.citation(
                "Reynolds et al. 2004",
                "Nat Biotechnol 22(3):326-330, criterion 1 of the rational design set",
            ),
            values={"minimum-gc-fraction": 0.30, "maximum-gc-fraction": 0.52},
            note=(
                "The band derived from 180 siRNA duplexes, and shifted distinctly lower than "
                "Wang's. It admits AT-rich candidates Wang excludes and rejects the whole "
                "52-60 % range Wang allows, so it is not a tightening of the default but a "
                "different band. Derived for knockdown rather than activation; worth trying if "
                "a promoter is AT-rich enough that the 40 % floor leaves you few candidates, "
                "but it is the least Wang-conformant choice here."
            ),
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
    options=(
        RuleOption(
            id="wang-2015",
            label="At most one CpG",
            provenance=TABLE_1,
            values={"maximum-cpg-count": 1},
            note="The tolerance stated in the paper.",
            is_default=True,
        ),
        RuleOption(
            id="legacy-scanner",
            label="No CpG at all (original script)",
            provenance=LEGACY_SCRIPT,
            values={"maximum-cpg-count": 0},
            note=(
                "Stricter than the paper. The original script rejected any target containing "
                "a CpG, which discards sound candidates in CpG-rich promoters."
            ),
        ),
        RuleOption(
            id="permissive",
            label="Up to two CpGs",
            provenance=Provenance.unsourced(
                "No source states a tolerance above one; offered only to widen a search that "
                "returns too few candidates."
            ),
            values={"maximum-cpg-count": 2},
            note=(
                "Beyond anything the literature supports. Use only to see what a very "
                "restricted promoter would otherwise yield, and treat the extra candidates "
                "as unvalidated."
            ),
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
        "types. Wang et al. say to avoid them but name no detection algorithm, so the island "
        "definition is this project's choice and is stated explicitly below."
    ),
    consequence=(
        "Switching this off admits candidates in the CpG-rich promoter core, which is exactly "
        "where many strong-looking candidates sit. In the original script this rule was "
        "written but commented out, so it never ran, and the legacy preset leaves it off for "
        "that reason. Because the rule excludes candidates, a loosely-defined island costs you "
        "usable candidates: the looser 1987 criteria call more than half their islands on Alu "
        "repeats, which are not islands at all."
    ),
    parameters=(
        RuleParameter(
            name="island-minimum-length",
            default=500,
            minimum=100,
            maximum=2000,
            unit="bp",
            description="Shortest window that can qualify as a CpG island.",
        ),
        RuleParameter(
            name="island-minimum-gc-fraction",
            default=0.55,
            minimum=0.0,
            maximum=1.0,
            unit="fraction",
            description="Lowest GC fraction for a window to qualify as a CpG island.",
        ),
        RuleParameter(
            name="island-minimum-observed-expected-cpg",
            default=0.65,
            minimum=0.0,
            maximum=2.0,
            unit="ratio",
            description=(
                "Lowest observed-to-expected CpG ratio for a window to qualify as an island."
            ),
        ),
    ),
    options=(
        RuleOption(
            id="takai-jones-2002",
            label="At least 500 bp, GC >= 55 %, obs/exp >= 0.65 (Takai & Jones 2002)",
            provenance=Provenance.citation(
                "Takai & Jones 2002", "PNAS 99(6):3740-3745"
            ),
            values={
                "island-minimum-length": 500,
                "island-minimum-gc-fraction": 0.55,
                "island-minimum-observed-expected-cpg": 0.65,
            },
            note=(
                "Recommended. Takai & Jones tightened the original thresholds precisely because "
                "those over-call: on human chromosomes 21 and 22 the 1987 criteria gave 14 062 "
                "islands, over half of them Alu repeats, where these gave 1 101. Since this "
                "rule excludes candidates, the stricter definition means fewer candidates lost "
                "to islands that are not really there. The trade-off is that a genuine island "
                "shorter than 500 bp is missed."
            ),
            is_default=True,
        ),
        RuleOption(
            id="gardiner-garden-frommer-1987",
            label="At least 200 bp, GC > 50 %, obs/exp > 0.60 (Gardiner-Garden & Frommer 1987)",
            provenance=Provenance.citation(
                "Gardiner-Garden & Frommer 1987", "J Mol Biol 196(2):261-282"
            ),
            values={
                "island-minimum-length": 200,
                "island-minimum-gc-fraction": 0.50,
                "island-minimum-observed-expected-cpg": 0.60,
            },
            note=(
                "The original and most widely cited definition, and the thresholds the legacy "
                "script's own detector used. Catches short islands the stricter criteria miss, "
                "at the cost of calling many Alu repeats as islands and so excluding candidates "
                "for sitting in something that is not an island."
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
            maximum=-1,
            unit="bp",
            description="Furthest upstream edge of the permitted window, relative to the TSS.",
        ),
        RuleParameter(
            name="window-end-bp",
            default=-100,
            minimum=-10_000,
            maximum=-1,
            unit="bp",
            description="Nearest edge of the permitted window, relative to the TSS.",
        ),
    ),
    options=(
        RuleOption(
            id="wang-2015",
            label="-1000 to -100 bp (Wang et al. 2015)",
            provenance=TABLE_1,
            values={"window-start-bp": -1000, "window-end-bp": -100},
            is_default=True,
            note=(
                "The window stated in the paper, and the one the four published saRNAs in the "
                "validation set all fall inside: VEGF at -708, TP53 at -394, dsP21-322 at -310 "
                "and VHL-821 at -818 from their MANE Select TSS. Recommended, with the caveat in "
                "the next option."
            ),
        ),
        RuleOption(
            id="vondervoort-2026",
            label="-2000 to -1 bp (upstream part of thesis chapter 4 window)",
            provenance=Provenance.citation(
                "van de Vondervoort 2026",
                "chapter 4, Material and methods, saRNA design, promoter defined as -2000 to "
                "+500 bp of the MANE Select TSS",
            ),
            values={"window-start-bp": -2000, "window-end-bp": -1},
            note=(
                "The thesis defined its promoter as -2000 to +500 bp. This option implements "
                "only its upstream portion, -2000 to -1; it does not reproduce the full thesis "
                "window. Activation was measured outside Wang's window: KLF4 worked at -45, "
                "closer than Wang's -100 edge, and CEBPA (in human clinical trials) at -1088, "
                "just beyond the -1000 edge. Choose this to avoid excluding candidates like "
                "those; expect a longer candidate list, without evidence that the added "
                "candidates are more effective. Downstream targeting is not expressed by this "
                "upstream-window rule. The last upstream base is -1: there is no TSS position 0."
            ),
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
        RuleParameter(
            name="nearest-neighbour-table",
            default="rna-xia-1998",
            choices=("rna-xia-1998", "dna-santalucia-1998", "legacy-dna-kj"),
            unit="table",
            description="Which set of stacking energies to compute the duplex ends with.",
        ),
        RuleParameter(
            name="end-length-nt",
            default=4,
            minimum=2,
            maximum=10,
            unit="nt",
            description=(
                "How many terminal nucleotides count as an 'end' when comparing the two ends' "
                "stability. Four, because the unwinding helicase reads about four bases at a "
                "time (Schwarz et al. 2003). Raising it averages over more of the duplex and "
                "blurs the terminal asymmetry the rule is testing."
            ),
        ),
    ),
    options=(
        RuleOption(
            id="rna-xia-1998",
            label="RNA parameters (Xia et al. 1998)",
            provenance=Provenance.citation(
                "Xia et al. 1998", "Biochemistry 37(42):14719-14735"
            ),
            values={"nearest-neighbour-table": "rna-xia-1998"},
            note=(
                "Recommended. An saRNA duplex is RNA:RNA, so these are the parameters that "
                "describe it."
            ),
            is_default=True,
        ),
        RuleOption(
            id="dna-santalucia-1998",
            label="DNA parameters (SantaLucia 1998)",
            provenance=Provenance.citation("SantaLucia 1998", "PNAS 95(4):1460-1465"),
            values={"nearest-neighbour-table": "dna-santalucia-1998"},
            note=(
                "Computes asymmetry on the DNA-equivalent sequence. Much of the siRNA design "
                "literature follows this convention, so it is offered for comparison, but the "
                "duplex that is actually loaded is RNA."
            ),
        ),
        RuleOption(
            id="legacy-dna-kj",
            label="Original script table (defective)",
            provenance=LEGACY_SCRIPT,
            values={"nearest-neighbour-table": "legacy-dna-kj"},
            note=(
                "DNA values in kJ/mol with three stacks assigned to the wrong parameter "
                "(defect B2), and no initiation or terminal penalty. Preserved only so older "
                "designs can be reproduced exactly. Do not use it for new work."
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
