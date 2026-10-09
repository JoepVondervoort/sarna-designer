# ADR 0002: Ship Wang 2015 conformance as the default; preserve the legacy Scanner as a pinned preset

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The legacy Scanner (ADR 0001) was written from **Wang J, Place RF, Portnoy V, Huang V, Kang MR,
Kosaka M, Ho MKC, Li L-C. *Inducing gene expression by targeting promoter sequences using small
activating RNAs.* J Biol Methods, 2015 (PMC4379447)**, whose Table 1 is the canonical rule set
and whose File S1 is an Excel implementation.

The Scanner diverges from Table 1 on five points:

| Wang 2015 Table 1 | Scanner | Nature of divergence |
|---|---|---|
| GC 40 – 60 % | 40 – **65 %** | more permissive |
| CpGs in target ≤ 1 | **zero** tolerance | more restrictive |
| CpG islands: avoid | `inisland` **commented out** | rule inert |
| Distance from TSS −1000 … −100 | **absent** | rule missing |
| Simple (di-/tri-nucleotide) repeats: avoid | **absent** | rule missing |
| 1st & 2nd nt preferably G or C | **absent** | rule missing |
| *(no such rule)* | **position 7 = `U`, weight 1000** | **unsourced**, yet a dominant score term |

Wang 2015 assigns **no weights to any criterion**. The Scanner's `1000 / 1000 / 1 / 0.01`
weighting is entirely local invention; the paper's own tool emits only an undisclosed "overall
quality score".

"Scanner is ground truth" (ADR 0001) and "Wang 2015 is the guideline" give *different answers*
on each row above. Choosing one wholesale either ships a tool that is non-conformant out of the
box, or discards the ability to reproduce historical designs.

## Decision

Separate **reproducibility** from **recommended practice** by shipping two presets.

**`legacy-scanner`** reproduces `legacy/sarna_promoter_scanner.py` bit-for-bit: GC ≤ 0.65,
zero-CpG tolerance, position-7 rule at weight 1000, islands rule off, no location rules, no
repeat rule, `i = 3` start offset included. This is what the regression suite asserts. It is
immutable: it is a historical record, not a tuning surface.

**`wang-2015`** is the **shipped default**: GC 40–60 %, CpGs ≤ 1, CpG islands avoided,
−1000 … −100 bp window enforced, di-/tri-nucleotide repeat and positions-1–2 G/C rules added,
and the position-7 rule **demoted to evidence class `exploratory` and off by default**.

Every divergence is therefore a visible, documented, switchable configuration difference. The
two presets can be diffed against the same promoter to show exactly which candidates the choice
costs or gains.

## Consequences

- Default output will differ from the legacy script. This is intentional and must be stated
  plainly in the UI and release notes.
- `legacy-scanner` carries known defects by design (notably the mis-paired nearest-neighbour
  parameters, B2). The preset description says so explicitly rather than hiding it.
- A third preset (`strict` / `exploratory`) is deferred: presets are cheap to add once the rule
  registry exists, and inventing them before the registry risks inventing biology.
- Expected success rates from Wang 2015, **15–20 % activate, ~70 % no effect, 5–10 % suppress**
 , are surfaced in the UI so a high rank is not read as a likely hit.
