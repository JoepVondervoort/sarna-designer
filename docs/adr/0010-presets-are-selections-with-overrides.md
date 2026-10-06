# ADR 0010 — A preset is a named selection of rules with per-rule overrides, and `legacy-scanner` reproduces the original with two stated exceptions

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

[ADR 0002](0002-wang-2015-default-legacy-preset.md) committed to two presets: `wang-2015` as the
default, and `legacy-scanner` reproducing the original script. [ADR 0003](0003-rule-provenance-is-first-class.md)
required the resolved configuration to be inspectable and exportable. Neither said what a preset
*is* in the code.

The question is not cosmetic, because the two presets differ in three distinct ways and the model
has to express all of them:

1. **A parameter value** — the GC ceiling is 0.60 under the paper and 0.65 in the script.
2. **Whether a rule applies at all** — the original has no concept of a TSS, so the TSS window
   rules cannot meaningfully run under `legacy-scanner`; `no-simple-repeats` is absent from the
   script entirely; and `not-in-cpg-island` is present but commented out (defect B8).
3. **Whether a rule is enabled against its own default** — `pos7-u` is off by default under
   `wang-2015`, being unsourced, but the original script scored it as one of two dominant terms,
   so reproducing historical rankings requires it on.

## Decision

A preset is a **named selection of rule ids, each with optional parameter overrides**, resolved
into a `DesignConfig` holding one `ResolvedRule` per selected rule.

Three properties follow, and each is load-bearing:

**Selection, not mutation.** A preset never edits a rule. The catalogue stays the single immutable
statement of what rules exist, and a preset is a projection of it. Switching preset cannot
therefore leave residue behind, and two presets can be diffed by comparing selections and
overrides rather than by diffing program behaviour.

**A disabled rule is absent, not failed.** A rule a preset does not select is not evaluated and
produces no result at all — distinct from a rule that was evaluated and reported
`not-applicable` (ADR 0007), which means "this run could not answer this question". The
difference matters when reading a candidate: "`tss-window` was not part of this preset" and
"`tss-window` could not be checked because the promoter has no TSS" are different statements and
the UI must not merge them.

**The resolved configuration is the exportable artefact.** `DesignConfig` serialises to the full
list of rules with their resolved parameter values and provenance, so a stored analysis records
what actually ran rather than a preset name whose meaning may later change. A preset name alone is
not a reproducible record.

### `wang-2015`

Every catalogued rule that is enabled by default, with default parameters. `pos7-u` is excluded,
being unsourced.

### `legacy-scanner`

| Rule | Treatment | Why |
|---|---|---|
| `gc-content` | `maximum-gc-fraction` = 0.65 | The script's band |
| `cpg-count` | `maximum-cpg-count` = 0 | Script rejects any CpG |
| `thermodynamic-asymmetry` | table = `legacy-dna-kj` | Defect B2 preserved (ADR 0009) |
| `asymmetry-margin` | table = `legacy-dna-kj` | Order by the same parameters that admitted |
| `pos7-u` | **enabled** | One of the script's two dominant ranking terms |
| `tss-window` | **not selected** | The script has no TSS concept |
| `tss-preferred-window` | **not selected** | Likewise |
| `not-in-cpg-island` | **not selected** | Present but commented out (defect B8) |
| `no-simple-repeats` | **not selected** | Absent from the script |

## Exceptions to exact reproduction

`legacy-scanner` is a **close** reproduction, not bit-for-bit. Two stated exceptions:

**1. Position-1 G/C is a preference, not a gate.** Finding B16 established that the script's
`endbase` reduces exactly to "position 1 is G or C and position 19 is A or U", making a required
G/C at position 1 — which Wang Table 1 marks only as "preferably". The decision taken here is to
**follow `wang-2015` and treat it as a preference in both presets**, rather than add a
legacy-only gate.

The cost, stated so it is not discovered later: `legacy-scanner` will **admit candidates the
original script rejected** — specifically every candidate whose position 1 is A or T. Those
candidates will appear in its output, ranked, where the original would have discarded them. For
anyone re-deriving a historical design list, that is a visible difference in the candidate set,
not merely in ordering.

Reversing this is a small change: add a mandatory `pos1-gc-required` rule to the catalogue and
select it in `legacy-scanner` only. It is recorded here rather than left implicit so that the
choice is revisitable.

**2. The ranking formula.** `legacy-scanner` ranks by criteria count like every other preset
([ADR 0006](0006-criteria-count-ranking.md)), not by the script's
`0.01·Δstab + 1000·pos18 + 1000·pos7 + 1·flank`. ADR 0006 already settled this; it is repeated
here because together with the above it means the legacy preset reproduces the original's
**filtering** closely and its **ordering** not at all.

## Consequences

- No preset is a bit-for-bit reproduction of the original script. ADR 0002's wording — "reproduces
  the original script bit-for-bit" — is accordingly too strong and the README must say "closely,
  with two stated exceptions" instead.
- A regression test against the legacy script therefore compares the candidate **set** under
  stated exclusions, not the output file. That is weaker evidence than a byte comparison, and it
  is the strongest evidence available given the above.
- Adding a rule to the catalogue does not change any existing preset, because a preset names its
  rules explicitly rather than taking everything enabled. New rules must be added to `wang-2015`
  deliberately. This is the safe direction: a new rule cannot silently start excluding
  candidates from a stored configuration.
- Because a preset's meaning can change between versions, exports record resolved rules rather
  than the preset id alone.
