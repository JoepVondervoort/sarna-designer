# ADR 0009: Guide-strand prediction reports every table's verdict, each explaining itself

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The legacy Scanner gates candidates on two separate asymmetry criteria and reports neither:

```python
if endstability(alltargets[i]) < 0:   # nearest-neighbour, 4 terminal bases
    ...rejected
if endbase(alltargets[i]) == False:   # single terminal base
    ...rejected
```

`endstability` sums three nearest-neighbour stacks over each terminal 4 nt and returns
`E3prime − E5prime`. A positive value means the target's 5' end is the more stable, hence the
antisense strand's 5' end is the weaker, hence RISC is expected to load the antisense strand
(Schwarz et al. 2003; Khvorova et al. 2003). The `end = 4` window is sourced in the script's own
comment to the helicase footprint reported by Schwarz.

The algorithm is sound. Its **table is not**: defect B2 assigns three stacks to the wrong
parameter, so `AC` carries `CA`'s value, `AG` carries `GA`'s and `TC` carries `CT`'s, where the
nearest-neighbour model requires `AC ≡ GT`, `AG ≡ CT` and `TC ≡ GA`. The individual errors are
around 0.1 kJ/mol, but they enter every end-stability sum, so some candidates fall on the other
side of the threshold, and the threshold is a mandatory gate. The table also applies DNA
parameters in kJ/mol to an RNA:RNA duplex and omits the initiation and terminal penalties.

[ADR 0002](0002-wang-2015-default-legacy-preset.md) requires `legacy-scanner` to keep that table
so historical designs reproduce exactly, and `wang-2015` not to. The consequence is unavoidable:
**for some candidates the two presets predict different guide strands.** How that is surfaced is
this decision.

## Decision

Guide-strand prediction is **not** a single verdict. For a candidate, the engine computes an
assessment per selected nearest-neighbour table, and the UI shows them together, each labelled
with the table it came from:

| Table | 5' end | 3' end | Difference | Predicted guide |
|---|---|---|---|---|
| RNA (Xia et al. 1998) | −7.61 | −5.36 | +2.25 kcal/mol | antisense |
| Legacy script (defective) | −21.5 | −22.1 | −0.6 kJ/mol | **sense** |

Three requirements follow.

**Every assessment states its reasoning in words, not just a number.** Each carries the per-stack
breakdown it summed, the two end energies, the difference, and a plain-language explanation of
why that difference implies that guide, naming the strand whose 5' end is weaker and the
mechanism it rests on. A researcher must be able to answer "why does it think that?" without
reading the source.

**Disagreement is surfaced, never resolved silently.** Where tables disagree the UI says so and
shows both. It does not average them, pick the default, or hide the minority verdict. A
disagreement is a real statement about the candidate: it means the asymmetry is small enough that
the parameter set decides the answer, which is itself a reason for lower confidence.

**A third outcome exists.** Where the difference is within the comparison tolerance the
prediction is `ambiguous`, not a coin-flipped strand. The legacy script's `fudge = 1e-8` is
retained as the default tolerance so the legacy preset's boundary behaviour is reproduced
exactly, but it is a configurable parameter rather than a buried constant.

The mandatory gate itself still uses the **one** table the active preset selects, so a preset
remains a single reproducible configuration. The other tables' assessments are reported
alongside as context; they never silently change which candidates pass.

## Alternatives considered

**Report only the active preset's verdict.** Simpler, and arguably the purist position: a preset
is a configuration and the user chose it. Rejected because the most useful thing the project can
tell a researcher about a legacy design is precisely that the corrected parameters disagree about
its guide strand, and that is invisible if only the active table is shown.

**Silently correct the legacy table.** Rejected by ADR 0002: it is a reproduction target, and a
corrected `legacy-scanner` reproduces nothing.

**Average the tables, or vote.** Rejected. The parameter sets are not independent measurements of
the same quantity; two of them are wrong for an RNA duplex and one of those is additionally
mis-assigned. Averaging would manufacture a number with no source.

## Consequences

- A candidate carries several guide assessments, so the candidate model holds a collection rather
  than one `predicted_guide` field, and exports have a row per table.
- The results table needs a column that can show agreement at a glance, with the per-table detail
  on candidate selection, or the comparison will be too noisy to read.
- `legacy-scanner` and `wang-2015` will rank and filter some candidates differently for this
  reason alone. That is intended and is the point of being able to diff the presets.
- The explanation text is part of the engine's output, not the UI's. It is therefore covered by
  tests, and a change to it is a visible change to the engine.
- Defect B2's practical effect becomes measurable: running both tables over a real promoter shows
  how many candidates it actually moves. That number belongs in the scientific documentation once
  a promoter fixture exists.

## Addendum: a second criterion, newly identified

Working through `endbase` for this decision showed it is not an independent thermodynamic test at
all. It compares one terminal base per end, scoring `A`/`U` as `AAE/2` and `G`/`C` as `GGE/2`, and
returns true only when the 3' base is the less stable. With just two possible values per end, the
condition reduces exactly to:

> position 1 is G or C **and** position 19 is A or U

So the legacy script enforces **as a mandatory gate** what Wang 2015 Table 1 marks only as
"preferably", a G or C at position 1. Combined with the existing `pos19-a` gate, the second half
is automatically satisfied, leaving `endbase` equivalent to a required G/C at position 1.

This is recorded in the legacy inventory as finding B16. It is not a defect in the arithmetic; it
is a scientific divergence between the presets that was previously unnoticed, and it means
`legacy-scanner` needs a mandatory position-1 rule that `wang-2015` scores as a preference. The
catalogue does not yet carry it.
