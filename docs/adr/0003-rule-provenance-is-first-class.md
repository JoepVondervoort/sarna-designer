# ADR 0003 — Rule provenance and rationale are first-class, user-visible metadata

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The author's stated reason for valuing Wang et al. 2015 was that it **explained why the score
was the way it was**, and asked that the application do the same — including letting the user
see *where each option's information comes from* and *why the option could matter*.

Against this, the legacy Scanner buries its science in module-level constants and bare `if`
conditions. A reader cannot tell that GC ≤ 0.65 is a local deviation, that the position-7 rule
has no source, or that `1000` was chosen to create tiers rather than to express an effect size.

Investigation also showed the paper itself publishes **no weights at all**, so "explain the
score like the paper does" cannot mean "copy the paper's explanation" — there isn't one. The
requirement is therefore to do better than the source.

## Decision

A design rule is not a predicate. It is a record carrying, at minimum:

- `id` and human-readable `name`
- `kind`: `filter` or `score_term`
- `evidence_class`: `established` | `empirical` | `exploratory`
- `provenance`: a citation **with locator** (e.g. "Wang et al. 2015, Table 1"), or
  `legacy-script`, or `unsourced`
- `rationale`: why this rule is believed to matter, biologically
- `consequence`: what changing it does to the candidate set, in plain language
- `default`, `valid_range`, `unit` where meaningful

Consequences of this record being mandatory:

1. **No rule may be added without provenance.** `unsourced` is a permitted value and is
   rendered as a warning; silence is not permitted. This is what reclassifies the position-7
   rule honestly rather than deleting it.
2. **The UI renders provenance next to every setting**, in Advanced Settings and in each
   candidate's rule breakdown, so the user can see which rules are paper-backed, which are
   inherited from the original script, and which are speculative.
3. **Per-rule score contributions are shown per candidate**, answering "why is #1 above #7?"
   with arithmetic rather than a single opaque number.
4. **Exports carry the rule records**, not just the resolved values, so an exported analysis
   remains interpretable without this codebase.

## Consequences

- The scoring surface must expose contributions, which rules out implementing ranking as one
  fused arithmetic expression. The legacy `1000/1000/1/0.01` formula is reproduced as a
  weighted sum **over named, individually reported terms** — identical arithmetic, decomposed
  reporting.
- The legacy weighting is additionally flagged in the UI as producing **colliding tiers**
  (`pos18=2, pos7=0` ties `pos18=1, pos7=1` at 2000), which is a real artefact of the chosen
  weights and exactly the kind of thing this ADR exists to make visible.
- Whether ranking stays a weighted sum or becomes explicitly tiered was left open here; this ADR
  constrains only *reporting*, not the ranking model. **Closed by
  [ADR 0006](0006-criteria-count-ranking.md)**, which replaces the weighted sum with an
  equally-weighted count of criteria met.
