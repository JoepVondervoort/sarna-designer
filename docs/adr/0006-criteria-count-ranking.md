# ADR 0006 — Rank by counting criteria met, not by a weighted sum

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

[ADR 0003](0003-rule-provenance-is-first-class.md) constrained how ranking is *reported* but
explicitly left the ranking *model* open. This ADR closes that question.

The legacy Scanner ranks with

```
score = 0.01·Δstab + 1000·pos18 + 1000·pos7 + 1·flank
```

Three things are wrong with it. The `1000 / 1000 / 1 / 0.01` weights are **local invention** —
Wang et al. 2015 publishes no weights for any criterion. The weights **collide**: `pos18=2,
pos7=0` ties `pos18=1, pos7=1` at 2000, so the tiers the weights were presumably chosen to
create do not cleanly exist. And the single dominant term `pos7` is **unsourced**: position 7
appears nowhere in Wang 2015 (see the inventory), yet it carries the same weight as a
paper-backed criterion.

The requested direction is a ranking closer to the paper, on the grounds that it is more
logical and easier to understand. That is right, with one correction that must be stated
plainly: **Wang 2015 contains no ranking formula to copy.** Its tool emits only an undisclosed
"overall quality score". So "rank like Wang" cannot mean reproducing his arithmetic; there is
none to reproduce.

What Table 1 *does* supply is a structure: a list of criteria, separating hard requirements
("avoid…", fixed ranges) from preferences, which the paper marks with the word **"preferably"**.
Ranking faithful to that structure is **counting**, not weighting.

## Decision

Ranking has three stages, in order.

### 1. Mandatory gate

A candidate failing any mandatory criterion is **excluded**, and the UI names the specific rule
that excluded it. These are not scored; they are admission requirements.

| Criterion | Source |
|---|---|
| Target length 19 nt | Table 1 |
| GC content 40–60 % | Table 1 |
| No run of ≥ 5 identical nucleotides | Table 1 |
| CpGs within target ≤ 1 | Table 1 |
| Not within a CpG island | Table 1 |
| No di-/tri-nucleotide simple repeats | Table 1 |
| Distance from TSS within −1000 … −100 bp | Table 1 |
| Nucleotide 18 is A or T | Table 1 |
| Nucleotide 19 is A | Table 1 |
| Thermodynamic asymmetry satisfied | Table 1; Schwarz 2003, Khvorova 2003 |

### 2. Preference count — the rank

Each preference met contributes **exactly one point. All preferences are weighted equally.**
Range 0–7.

| Preference | Points | Source |
|---|---|---|
| Nucleotide 1 is G or C | 0–1 | Table 1, "preferably" |
| Nucleotide 2 is G or C | 0–1 | Table 1, "preferably" |
| Nucleotides 20, 21, 22, 23 each A or T | 0–4 | Table 1, "preferably" |
| TSS distance within the preferred −500 … −200 bp window | 0–1 | Wang 2015 **text**, not Table 1 |

Table 1 uses "preferably" of positions 1–2 and of the 20–23 flank. Those two groups are scored
per-position rather than all-or-nothing, because a candidate meeting three of the four flank
positions is genuinely closer to the stated preference than one meeting none, and collapsing
them to a single yes/no would leave almost every candidate tied.

The −500 … −200 preference is **not from Table 1** and is labelled accordingly in the UI; it
derives from the paper's narrative preference for that sub-region (see
[ADR 0001](0001-scanner-is-ground-truth.md)) and sits inside the mandatory −1000 … −100 window.

**Why equal weights:** the evidence provides no basis for ranking these preferences against one
another. Assigning unequal weights would be inventing an effect size the literature does not
report — precisely what produced the legacy `1000/1000/1/0.01` artefact. Equal weighting is a
deliberate, stated position, not a default.

### 3. Tie-breaking

Ties are broken by, in order:

1. **Thermodynamic asymmetry magnitude**, descending. Continuous and genuinely established
   (Schwarz 2003, Khvorova 2003), so it discriminates where the checklist cannot.
2. **Genomic position**, ascending, purely so that output is deterministic.

Tie-breaks are reported as tie-breaks, never folded into the score.

### Presentation

A candidate's rank is shown as **"meets 5 of 7 preferences"** with the checklist itemised, each
line carrying its provenance per ADR 0003. There is no composite number, because a composite
number is what made the legacy score uninterpretable.

The position-7 rule stays where [ADR 0002](0002-wang-2015-default-legacy-preset.md) put it:
evidence class `exploratory`, off by default in `wang-2015`. If enabled it appears as an
eighth preference worth one point, like every other — never a dominant term.

## Consequences

- Rankings become coarse by design: many candidates share a preference count. This is **honest**
  rather than a defect. Wang reports that only 15–20 % of designed duplexes activate, so fine
  discrimination between candidates would be false precision. The UI presents ties as ties.
- `legacy-scanner` is unaffected. It keeps the weighted sum verbatim, including the collisions,
  because it is a reproduction target and immutable (ADR 0002). The two presets will therefore
  rank differently, which is the point of being able to diff them.
- Scores are not comparable across presets, and the UI must not imply they are.
- Adding a preference later shifts every score, so the preference set is versioned with the
  preset and recorded in exports.
- The score is not a probability of activation and must never be rendered as one.
- Implementation stays decomposed per ADR 0003: each preference is a named, individually
  reported term. Here that falls out for free, since the score *is* the itemised list.
