# ADR 0011: Optional rules are hypotheses, and the default preset never carries one

**Status:** accepted
**Date:** 2026-10-06

## Context

Two bodies of evidence arrived that are not Wang et al. 2015, and both bear on sequence design.

**Thesis chapter 4** screened 10 promoter-targeting duplexes against each of six epilepsy genes
and reported positional nucleotide preferences separating its 14 activators from 111
non-activators: G at target position 10 (+24.8 pp), U at 15, C depleted at 9, and a stability
gradient with strong bases centrally and weak bases toward the 3' end. The chapter states these
are exploratory.

They are weaker than exploratory. Reanalysis
(`docs/scientific/positional-motif-analysis.md`) reproduces every reported value exactly and
then measures what the chapter does not: a max-statistic permutation test over the grid that was
searched. At every motif width the observed maximum sits **below the median of the null**, a
random relabelling of the same 125 sequences produces more positional structure than the real
activity labels do. Family-wise p = 0.94, 0.67 and 0.73 for 1-, 2- and 3-mers.

**The siRNA design literature** offers criteria derived from far larger panels: Reynolds et al.
2004 (180 duplexes) and Ui-Tei et al. 2004 (62 duplexes, four conditions an effective duplex must
satisfy simultaneously). These are well replicated, but they were derived for *knockdown*.

So: a lot of material a user might reasonably want to apply, of sharply unequal quality, none of
it from the paper this project treats as canonical.

The forces in tension:

- Researchers need to investigate alternative criteria. Recording those choices in the tool
  preserves provenance that would otherwise be lost in external analyses.
- Adding a rule to the default ranking changes every design the project recommends. The evidence
  above does not come close to justifying that.
- "Optional" is not self-explaining. A rule offered without its cost is a trap (the existing
  `RuleOption.note` contract already says so).

## Decision

**Every rule from outside Wang et al. 2015 is catalogued, off by default, as a `PREFERENCE`, and
states in its `consequence` why a user would decline it.**

Four parts, each load-bearing:

**Catalogued, not hidden.** They live in `src/sarna/rules/catalog/optional.py` with full
provenance, so selecting one is a recorded, exportable decision rather than an edit to a number.

**Off by default.** `wang-2015` selects none of them. Its denominator stays 7. A test asserts
this directly, because the failure mode is silent: adding a default-on rule to the catalogue
would change every ranking the project produces without any diff touching the preset.

**Preference, never mandatory.** An optional rule can reorder candidates; it can never exclude
one. This bounds the damage a bad hypothesis can do: a user who selects every optional rule still
sees the same candidate set as the default, differently sorted. A test asserts no optional rule
is a gate.

**The `consequence` field carries the argument against.** Not a disclaimer, the specific reason.
For `pos10-u`: activation does not proceed by cleavage, so the geometric reason the criterion
exists in knockdown does not apply. For `pos10-g`: it fails correction, and it contradicts
`pos10-u`. For `central-stability`: it partly double-counts the mandatory asymmetry rule.

### Evidence class is assigned by mechanism, not by citation count

A rule is `ESTABLISHED` only where the mechanism it invokes is one activation and knockdown
*share*. Both pathways begin with duplex recognition, strand selection and AGO2 loading; they
diverge after that, and only knockdown cleaves anything.

| Rule | Class | Why |
|---|---|---|
| `au-rich-guide-5-prime` | ESTABLISHED | Strand selection (shared mechanism, three agreeing sources |
| `no-gc-stretch` | ESTABLISHED | Duplex stability and secondary structure) shared |
| `pos3-a` | EXPLORATORY | No mechanism offered by its source |
| `pos10-u` | EXPLORATORY | Mechanism is cleavage geometry, which activation does not use |
| `no-g-at-13` | EXPLORATORY | No mechanism offered by its source |
| `pos10-g` | EXPLORATORY | Fails correction for the grid it was found in |
| `central-stability` | EXPLORATORY | Same, and overlaps the asymmetry gate |

`au-rich-guide-5-prime` is the one rule here with a genuine claim on a user's attention, and it is
worth saying why: Reynolds criterion 2 (≥3 A/U at sense positions 15-19), Ui-Tei criterion iii
(≥5 A/U in the guide's 5' terminal third) and chapter 4's stability analysis all identify the same
stretch of the same molecule, and the first two were derived from panels large enough to mean
something. Chapter 4 found it without looking for it.

### Two new presets

`strict` promotes Wang's own softer statements to requirements, the −500 to −200 bp region the
paper reports activity in becomes the admission window, no CpG is tolerated, asymmetry must exceed
0.5 kcal/mol. It adds no preferences, so its scores stay comparable with the default's. Both
divergences are declared, including that the asymmetry threshold is this project's invention.

`exploratory` adds six optional rules to the default. It declares that its denominator is 13, so
its scores are **not comparable** with any other preset's, and that its ordering is a hypothesis
rather than a better estimate. `pos10-u` is deliberately excluded: shipping it alongside `pos10-g`
would guarantee every candidate loses one of the two points, which is noise dressed as scoring.

## Consequences

**Good:**

- The default remains exactly what the README claims it is.
- A user can test a hypothesis and see, in an export, which hypothesis they tested.
- The strongest optional rule is identified as such, rather than sitting in an undifferentiated
  list beside rules that contradict each other.
- The contradiction between `pos10-u` and `pos10-g` is surfaced as a property of the evidence
  instead of being resolved silently by whoever wrote the catalogue.

**Costs:**

- 26 rules instead of 19. The catalogue is now large enough that a UI must group it; `group` is
  set on every optional rule for exactly this reason.
- `exploratory`'s denominator of 13 is a real footgun. Mitigated by a declared divergence, but a
  user who compares "9 of 13" with "5 of 7" will still draw a wrong conclusion. A UI showing two
  presets side by side must suppress the comparison rather than render both numbers.
- Offering `pos3-a` and `no-g-at-13`, for which no mechanism exists in either pathway, is close to
  offering noise. Carried because Reynolds' set is widely used and a user will ask for it; each
  says plainly that no mechanism supports it.

## Alternatives rejected

**Leave them out entirely.** The honest-looking option, and the one that loses most. The user has
the thesis and the siRNA literature either way; without the catalogue they apply them by hand, with
no provenance, no export and no statement of cost. Offering a rule with its argument against is
more protective than not offering it.

**Add the well-supported ones to the default.** Rejected for `au-rich-guide-5-prime` even though
it is defensible, because the default's one job is to be Wang et al. 2015 Table 1, and a reader
comparing this tool's output to the paper must find no unexplained difference. It is one
selection away for anyone who wants it.

**Weight rules by evidence class.** Rejected: it reintroduces exactly the weighted sum ADR 0006
removed, and this time with weights derived from a judgement about evidence quality rather than
from any measurement. Counting criteria and reporting the count is still the only scheme that
survives a reader asking where a number came from.

**Report the thesis positions as a separate "hypothesis score".** Rejected as premature. A second
number next to the first invites a comparison the data cannot support. If the positional
hypotheses are ever tested prospectively (which needs one pre-registered position, not 76), a
score becomes worth building.
