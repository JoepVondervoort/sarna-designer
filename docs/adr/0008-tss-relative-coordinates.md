# ADR 0008: TSS-relative positions use the biological no-zero convention, and the whole target must lie inside the window

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The legacy Scanner has no concept of a transcription start site. It scans whatever string is
pasted into it and reports a candidate's position as an index into that string, so there is no
existing behaviour to preserve here and no precedent to inherit. Two conventions therefore have
to be chosen outright, and both are load-bearing: the mandatory `tss-window` rule
(−1000 … −100 bp) and the `tss-preferred-window` preference (−500 … −200 bp) are both expressed
in these terms, so the choices decide which candidates exist and which are excluded.

Neither choice is reversible in a quiet way. Changing either later changes the candidate set for
every analysis and silently invalidates stored results, which is why they are recorded here
rather than left as an implementation detail.

## Decision

### 1. TSS-relative positions follow the biological convention, with no zero

`+1` is the first transcribed base. `−1` is the base immediately upstream of it. **There is no
position 0.** The distance from `−1` to `+1` is one base pair.

This is the convention of Wang et al. 2015 and of the promoter literature generally, and the one
the UI ruler already shows (`−1000 … −500 … −100 … TSS/+1 … +100`). A window written
`−1000 … −100` therefore means exactly what a reader of the paper expects, with no translation
between what the engine computes and what the interface displays.

The cost is a discontinuity at the TSS, which is a genuine source of off-by-one error. It is
contained by converting through a **continuous internal offset** in exactly one place:

| TSS position | internal offset |
|---|---|
| `−2` | `−2` |
| `−1` | `−1` |
| `+1` | `0` |
| `+2` | `1` |

All arithmetic (distances, shifts, window containment) happens on the offset. No other module
performs arithmetic on a TSS position directly.

### 2. A candidate satisfies a TSS window only if **every** base of it is inside

The mandatory window is a statement about where the target lies, so a target half of which sits
outside the window does not lie inside it. A 19-mer whose 5' end is at `−105` extends to `−87`,
overlapping the transcript, which is the outcome the upstream window exists to prevent.

Concretely, with the default `−1000 … −100`:

```
window     −1000 ──────────────────────── −100
candidate            [−347 … −329]              passes
candidate                      [−105 … −87]     fails: extends past −100
```

The same containment test applies to the `−500 … −200` preference.

## Alternatives considered

**Zero-based offsets throughout**, with `0` as the first transcribed base. Continuous arithmetic,
no special case, one fewer class of off-by-one. Rejected because every upstream number shown to a
researcher would then be one less than the cited paper's, so the interface would either display a
different convention from the one it computes in or contradict the literature it cites. The
discontinuity is real but it is confined to one conversion; a convention mismatch would be spread
across the UI, the exports and every conversation about a result.

**Anchoring on the 5' end only**, or **any overlap**. Both admit candidates reported as satisfying
a window while lying partly or mostly outside it. Rejected: a rule that is close to meaningless at
its own boundary is worse than a narrower rule.

## Consequences

- The effective window is 18 bp narrower than its stated bounds for a 19 nt target. The UI must
  therefore show the window **and** the candidate span, not just the bound, or the arithmetic will
  look wrong to a user who checks it.
- Candidate count falls slightly relative to a 5'-anchored test. This is intended.
- Changing the target length moves the effective window, because containment depends on length.
  That coupling is already noted in the `target-length-19` rule's `consequence` text.
- TSS positions are compared and shifted only through the dedicated value type. A raw `int`
  TSS position in any other module is a defect.
- A promoter with no known TSS has no TSS-relative positions at all, so both window rules report
  `not-applicable` rather than passing (ADR 0007).
