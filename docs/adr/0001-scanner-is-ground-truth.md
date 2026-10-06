# ADR 0001 — `saRNA Promoter Scanner.py` is the scientific ground truth

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

Three legacy files were supplied. Two (`saRNA Promoter Scanner.py`,
`saRNA by Promoter Sequence.py`) are **byte-identical** — one implementation, not two. The
third, `EpiMarker checker .py`, is a separate pipeline that locates *already-designed* 19-mers
and analyses histone signal around them.

The two pipelines disagree on the design rules:

| | Scanner | EpiMarker |
|---|---|---|
| GC bounds | 0.40 – 0.65 | 40 – 60 % |
| Score | `0.01·Δstab + 1000·pos18 + 1000·pos7 + 1·flank` | `0.7·GCscore + 0.3·distanceScore`, ideal −500 bp |
| Candidate generation | yes, enumerates all 19-mers | no, string-searches a supplied list |
| TSS / coordinates | none at all | Ensembl canonical transcript, ±30 kb |

Regression tests cannot be written against two conflicting definitions, and Phase 2 of the
build is blocked until one is authoritative.

## Decision

**`saRNA Promoter Scanner.py` is the sole ground-truth implementation** for candidate
generation, filtering and ranking. It is preserved verbatim at
`legacy/sarna_promoter_scanner.py` and is the regression oracle for the `legacy-scanner`
preset.

`EpiMarker checker .py` is **not** a design-rule authority. Its `validate_dna` GC bounds and
its `compute_score` formula are retired rather than reconciled, because both were written to
screen an existing design list, not to generate one.

EpiMarker's genuinely valuable parts are retained but relocated:
- Ensembl TSS/transcript retrieval → genomics layer
- bigWig signal extraction → regulatory annotation layer
- POS_CTRL / NON_CTRL enrichment analysis → a validation tool, not part of design

## Consequences

- The regression suite pins Scanner behaviour exactly, including its known defects, under the
  `legacy-scanner` preset. See ADR 0002 for how corrections are surfaced without losing this.
- EpiMarker's `ideal = −500 bp` is **not** discarded as arbitrary: it traces to Wang et al.
  2015's preferred −500/−200 bp region. It returns as a location preference rule with proper
  provenance, not as EpiMarker's formula.
- The duplicate second script is retired. It is not preserved separately, being byte-identical.
- `Designs.xlsx`, EpiMarker's required input, is missing and cannot be supplied (upload
  limitation). Wang et al. 2015 **File S1** is proposed as a substitute validation oracle.
