# Why synthetic promoter fixtures must be CpG-depleted

**Found 2026-10-06** while running the first end-to-end design on a synthetic promoter.

## The observation

A 1200 bp promoter of uniform random sequence at 54 % GC produced **zero** passing candidates
under `wang-2015`. The cause was the island rule: the whole promoter was called as a single
974 bp CpG island, excluding 985 of 1182 candidates on its own.

The detector is right and the fixture is wrong.

| Sequence | GC | CpG obs/exp | Islands called |
|---|---|---|---|
| Uniform random (`ACGT`, 54 % GC) | 53.9 % | **1.00** | 1, covering 974 bp |
| CpG-depleted random | 42.2 % | **0.38** | 0 |

## Why

A CpG island is defined by CpG dinucleotides being *more frequent than base composition predicts*
— observed over expected at least 0.65 under Takai & Jones. Uniform random sequence has obs/exp
≈ 1.0 **by construction**, because expected is derived from the same composition that generated
it. Any correct island detector therefore calls random GC-balanced sequence an island.

Real genomic DNA does not look like that. Outside islands it is strongly CpG-depleted, obs/exp
typically around 0.2, because methylated cytosine deaminates to thymine and CpG has been lost to
TpG over evolutionary time. **That depletion is what makes an island detectable at all**: an
island is a region that escaped the depletion. A fixture without the depletion has no baseline for
an island to stand out from, so everything is an island.

## Consequence for this repository

A synthetic promoter fixture must either:

- be generated with CpG depletion — suppress most `CG` dinucleotides, targeting obs/exp around
  0.2; or
- be real genomic sequence, which is the better option for any fixture used to make a scientific
  claim.

Uniform `random.choices("ACGT")` is acceptable only for tests that do not involve the island rule.
A test fixture built that way and used with `not-in-cpg-island` selected will reject nearly every
candidate, and the failure will look like a bug in the design rules rather than in the fixture.

The same caution applies to interpreting any result: if a user pastes a sequence that is not CpG-
depleted — a synthetic construct, a shuffled control, a designed plasmid promoter — the island
rule will fire across it. That is the rule behaving correctly on input it was not designed for,
and the UI should say so rather than silently returning no candidates.

## Verified behaviour after correcting the fixture

The same promoter, CpG-depleted, at 42 % GC:

```
generated 1182   rejected 1141   passing 41   islands 0

pos19-a                  828   70.1 %
thermodynamic-asymmetry  587   49.7 %
gc-content               532   45.0 %
pos18-at                 499   42.2 %
tss-window               299   25.3 %
cpg-count                 42    3.6 %
no-homopolymer-run        15    1.3 %
```

The `pos19-a` rejection rate of 70 % independently corroborates the catalogue's own claim that it
"cuts roughly three quarters of candidates on its own", which was written from the paper rather
than from measurement.
