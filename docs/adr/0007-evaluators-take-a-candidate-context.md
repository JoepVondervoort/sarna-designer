# ADR 0007: Rule evaluators are pure functions over a single `CandidateContext`

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

[ADR 0003](0003-rule-provenance-is-first-class.md) made a rule *data*: a record carrying its
identity, provenance, rationale and parameters. The catalogue built on that
([`sarna.rules.catalog`](../../src/sarna/rules/catalog/)) can be imported and rendered in a UI
today without any evaluation code existing at all, which is a property worth keeping.

What does not exist yet is the other half: the code that decides whether a given 19-mer actually
*meets* a rule. That requires choosing how an evaluator receives its inputs, and the rules in
the catalogue do not all want the same inputs:

- `gc-content`, `no-homopolymer-run`, `cpg-count` and the position preferences need **only the
  target sequence**.
- `flank-20-au` … `flank-23-au` need the target **plus its 3' flank**, so they need the
  surrounding promoter, not just the window.
- `tss-window` and `tss-preferred-window` need a **strand-aware signed distance to a TSS**, and
  therefore assembly-tagged coordinates.
- `not-in-cpg-island` needs a **feature computed over the whole promoter**, not the window.
- `thermodynamic-asymmetry` and `asymmetry-margin` need **nearest-neighbour energies** and hence
  a chosen table.

So the inputs range from one short string to the full genomic context.

## Decision

Every rule evaluator is a **pure function** with one uniform signature:

```python
def evaluate(ctx: CandidateContext, params: ResolvedParameters) -> RuleResult
```

`CandidateContext` is an immutable record holding everything any catalogued rule could need: the
target window, its flanks, the promoter it came from, the signed TSS-relative distance, the
assembly-tagged `GenomicInterval` where one is known, promoter-level derived features such as
CpG-island spans, and the thermodynamic features computed from the configured table. A rule that
needs none of the genomic fields simply does not read them.

Three consequences are deliberate.

**The context is built once per candidate, not once per rule.** Promoter-level work (CpG-island
detection, which is O(n·200) in the legacy implementation (defect B7)) happens once for the
promoter, not once per window per rule.

**Genomic fields are optional, and a rule says so.** Workflow B (a pasted promoter with no
coordinates) produces a context whose `interval` and `tss_offset` are `None`. A coordinate-
dependent rule must then report `RuleOutcome.NOT_APPLICABLE` with the reason, and the UI shows
that it could not be evaluated. It must **not** silently pass, and it must not guess a position.
This is the mechanism that keeps "we could not check this" distinct from "this passed".

**Evaluation stays out of the catalogue.** Evaluators live beside the rules they implement but
are reached through a lookup keyed on rule id, so importing the catalogue still pulls in no
evaluation code, and a rule with no evaluator is a loud startup error rather than a rule that
quietly never fires.

## Alternatives considered

**Typed per-rule inputs**, each rule declaring exactly the arguments it needs. More honest
about dependencies and the smallest test setup per rule, but it needs a dispatch layer to call
heterogeneous signatures from a generic loop, and that layer is where the "which rules ran"
bookkeeping would get lost.

**An `evaluate` method on `Rule`**, metadata and logic in one place per rule. Rejected because
it couples the declarative catalogue to execution and gives up importing the catalogue without
the engine, which ADR 0003 bought on purpose.

## Consequences

- Every rule is testable from a hand-constructed context, with no network, no promoter file and
  no fixtures. This is what makes the strand tests cheap enough to be exhaustive.
- The context must be constructible without coordinates, so its genomic fields are `| None` and
  every evaluator reading them handles the `None` case explicitly.
- `RuleResult` needs a third outcome beyond pass/fail (`not applicable`) and the ranking must
  treat it as neither a met preference nor a failed gate, reported as unevaluated.
- Adding a field to the context touches one type rather than every evaluator signature; the risk
  in exchange is a context that accretes fields, so a field is added only when a catalogued rule
  needs it.
- Context construction is the single place where strand and coordinate arithmetic happens. That
  concentrates the defect-B1 class of bug into one well-tested module instead of scattering it.
