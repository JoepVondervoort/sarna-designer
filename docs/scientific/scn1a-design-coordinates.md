# Where the thesis's SCN1A designs actually sit

**Finding.** All 15 SCN1A designs from chapters 2-4 map to GRCh38 unambiguously. They form two
clusters around two different transcription start sites, and **none of them targets the MANE Select
promoter**, which is what this tool selects by default. Five of the fifteen, including two of the
three activators, lie 5.5-6.0 kb from any annotated TSS — far outside the window Wang's criteria
permit.

Measured by `scripts/locate_scn1a_designs.py` against Ensembl GRCh38, 2026-10-06. Reproducible;
needs the network.

## Method

The chapter 4 appendix records target sequences without coordinates, which is why the location
rules could not be tested against measured activity. Recovering them needs only a search: each
19-mer is searched for in the SCN1A locus plus 25 kb of padding, in both orientations. SCN1A is on
the minus strand, so a target reading 5'→3' in transcript orientation appears reverse-complemented
in the plus-strand reference, and its 5' end is the **highest** coordinate of the match.

All 15 matched exactly once, on the minus strand, as expected for promoter-targeting designs
against a minus-strand gene. Each design's offset is then computed against **every** one of the
gene's 16 annotated TSSs rather than against an assumed one.

## The designs sit at two TSSs, neither of them MANE Select

SCN1A has 30 transcripts, 12 protein-coding, with 16 distinct TSSs spanning 75 kb.
MANE Select is `ENST00000674923`, TSS at chr2:166,128,019.

| Cluster | TSS | Designs | Offsets | In Wang's window | Activators |
|---|---|---|---|---|---|
| A | 166,149,213 | 5 | −124 to −744 | 5 / 5 | SCN1A-10 |
| B | 166,073,669 | 10 | −224 to −6,010 | 5 / 10 | SCN1A-12, SCN1A-14 |
| — | 166,128,019 (MANE Select) | **0** | — | — | — |

Cluster A is textbook: five designs from −124 to −744 bp, all inside the −1000 to −100 window,
targeting the TSS shared by `ENST00000409050` and three other transcripts.

Cluster B splits in two on inspection. SCN1A-1 to SCN1A-5 sit at −224 to −778 bp of TSS
166,073,669, which is again textbook. SCN1A-11 to SCN1A-15 sit at **−5,523 to −6,010 bp** of the
same TSS — the nearest of all sixteen — which is not a promoter position by any definition this
project uses.

## The distant cluster contains two of the three activators

This is the part that matters, because it bears on the mechanism rather than on bookkeeping.

| Design | Offset from nearest TSS | Activity |
|---|---|---|
| SCN1A-11 | −6,010 | no |
| SCN1A-12 | −5,669 | **activating** |
| SCN1A-13 | −5,659 | no |
| SCN1A-14 | −5,530 | **activating** |
| SCN1A-15 | −5,523 | no |

Two of the project's three SCN1A activators are ~5.5 kb upstream of the nearest annotated TSS. The
third, SCN1A-10, is at −124 bp of a different TSS. So the SCN1A activators do not share a promoter.

**Resolved by chapter 2.** The question was open only until chapter 2's supplementary material was
read. These designs are not mysterious and not misplaced: they target **experimentally mapped
promoters from the literature**, which chapter 2 names explicitly.

| Block | Chapter 2's attribution | Span | Nearest Ensembl TSS | Activators |
|---|---|---|---|---|
| SCN1A-1 to -5 | "promoter region described by **Dong (2014)**" | chr2:166,073,893-166,074,447 | 166,073,669, −778…−224 | none |
| SCN1A-6 to -10 | "promoter region described by **Long (2008)**" | chr2:166,149,337-166,149,957 | 166,149,213, −744…−124 | SCN1A-10 |
| SCN1A-11 to -15 | "promoter region by **Dong (2014)**" | chr2:166,079,192-166,079,679 | 166,073,669, **5.5 kb away** | SCN1A-12, SCN1A-14 |

Chapter 2 states that "the most proximal and distal promoters were selected as target regions of
the SCN1A gene", citing Dong et al. 2014 (*Mol Neurobiol* 50:438-448, SCN1A transcription repressed
by RACK1), Long et al. 2008 (*J Neurosci Res* 86:3375-3381, identification of the SCN1A promoter
region and 5′-untranslated exons) and Nakayama et al. 2010 (*Hum Mutat* 31:820-829, 5′ genomic
deletions with promoter activity).

Two of the three blocks line up with an Ensembl TSS. **The third does not, and it produced two of
the three SCN1A activators.** The most productive SCN1A promoter in this project is one Ensembl
annotates no transcript for.

That is the finding with the most bearing on the tool, and it generalises well past SCN1A: **a
promoter worth targeting is not always one the annotation knows about.** A tool that can only design
relative to Ensembl TSSs could not have produced SCN1A-12 or SCN1A-14 — two of this project's three
SCN1A hits — because `tss-window` excludes anything beyond −1000 bp of an annotated TSS. Supporting
literature-defined promoter regions as first-class design anchors is a requirement, not a
convenience.

One narrower question stays open: whether the Dong 2014 distal region acts as a promoter in the
cells used, as an enhancer, or through the SCN1A natural antisense transcript. Chapter 2's Figure 1
presents it as a promoter on the strength of the cited mapping studies.

## Ensembl's regulatory build annotates the MANE promoter and nothing else nearby

Queried at each TSS, 2026-10-06:

| Region | Features |
|---|---|
| MANE Select TSS (166,127,500-166,129,000) | `promoter` 166,128,005-166,128,231; `enhancer` 166,127,800-166,128,004; `CTCF_binding_site` 166,128,090-166,128,111 |
| TSS 166,149,213 ± 750 | none |
| TSS 166,073,669 ± 1,000 | none |
| Third cluster (166,078,000-166,081,000) | none |

The annotated promoter feature spans the MANE Select TSS, as it should. **Absence elsewhere is weak
evidence**, and should not be read as "these are not promoters": the Ensembl Regulatory Build is
built from epigenomic data aggregated over a panel of cell types that includes neither AF22 nor any
neural progenitor, and a promoter active only in neurons has every reason to be missing from it.
What the table does establish is that no *independent* source consulted here supports the two TSSs
the thesis designed against, while one does support the TSS it did not.

## Consequences for the tool

**Also established from chapter 2's supplementary material:** both duplex strands carry a **dTdT**
overhang, written `[dT][dT]`, not a ribo-UU. That settles the open question about overhang
chemistry in this project's own practice, and it is what an order-ready export should emit.

**Acted on:**

- ADR 0012 records why MANE Select is the default and why a tool offering only MANE Select would
  have been unable to reproduce this project's own earlier SCN1A work. The client returns every
  transcript, so a UI can offer the choice.
- The engine must show **all** TSSs on a promoter view, not only the selected one, so a user can
  see what they are not designing against. Noted as a UI requirement.

**Worth doing next:**

- Check an archived Ensembl release for a TSS near 166,079,200. One query, and it would settle
  reading 2 against readings 3 and 4.
- Feed these coordinates back into the validation harness. They are the promoter context the
  chapter 4 dataset lacks, and with them the TSS windows, CpG-island overlap and the 3' flank
  preferences become testable against measured activity for 15 of the 125 sequences — currently
  five of seven preferences cannot be evaluated at all
  (`docs/scientific/empirical-validation.md`).

## Designing against the declared anchors reproduces chapter 2 — and finds an edge effect

With the three promoters declared as anchors (`sarna.genome.anchors`), the engine can be pointed at
the same regions chapter 2 used. All 15 designs are then generated, and most pass:

| Anchor | Chapter 2's designs, as the engine judges them |
|---|---|
| `scn1a-long-2008` | SCN1A-7, -8, -9 and **-10 (activator)** pass. SCN1A-6 fails `no-homopolymer-run`. |
| `scn1a-dong-2014-proximal` | SCN1A-3, -4, -5 pass. SCN1A-1 fails `pos18-at`, SCN1A-2 fails `pos19-a`. |
| `scn1a-dong-2014-distal` | **SCN1A-12 (activator)** and SCN1A-13 pass. SCN1A-11 fails `pos19-a`. **SCN1A-14 (activator)** and SCN1A-15 fail `tss-window`. |

Before anchors existed, SCN1A-12 and SCN1A-14 were rejected for being 5.5 kb from any annotated
TSS. Anchored properly, SCN1A-12 passes and SCN1A-14 fails for a quite different and much narrower
reason.

### Total containment costs 2% of candidates, and one of them worked

SCN1A-14 sits at −107. Its 19-mer therefore spans **−107 to −89**, and ADR 0008 requires the
*whole* target inside the window, so the 7 bases past −100 exclude it. SCN1A-15 at −100 spans
−100 to −82 and goes the same way.

Measured over the distal anchor's 1,382 candidates: 901 have their start inside −1000…−100 and 883
have the whole target inside. **Total containment costs 18 candidates, 2.0%**, all of them
straddling one of the two edges — 36 bp of a 900 bp window, which is about what the arithmetic
predicts.

Two percent is a small price and ADR 0008's reasoning for it still stands: a target half outside
the window is half outside the promoter region the rule is about. But one of the 14 duplexes known
to activate its target falls in that 2%, which is worth knowing before trusting a rejection at the
window edge.

**Caveat on the number, not the mechanism.** SCN1A-14's −107 depends on this anchor's *inferred*
reference point, which chapter 2's designs constrain only to a ~410 bp range. Under a different
value within that range SCN1A-14 could sit comfortably inside the window. The edge effect itself is
a property of the rule and does not depend on the anchor: a 19-mer whose start is within 18 bases
of either edge is excluded, wherever the edges are.
