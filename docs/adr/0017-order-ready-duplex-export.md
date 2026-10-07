# ADR 0017 — The order-ready duplex is a separate export, computed in Python, carrying its caveats

Date: 2026-10-07. Status: accepted implementation of accepted constraints 4, 11 and 15.

## Context

A candidate carries a 19 nt RNA core, which is what every rule is evaluated against. What a
researcher buys is two 21-mers. `CLAUDE.md` named this the largest gap before bench use: the
chemistry question was already settled and only the export was missing.

Chapter 2's supplementary material settles the chemistry — a 2 nt **dTdT** 3' overhang on **both**
strands, written `[dT][dT]`, not the ribo-`UU` the legacy scripts assumed. Nothing in this ADR
reopens that; `docs/scientific/scn1a-design-coordinates.md` records the finding and `GLOSSARY.md`
the vocabulary. What needed deciding was where the chemistry is computed, what travels with it
and what the export may not imply.

## Decisions

**The ordered strand is not a `NucleicSequence`.** A 21-mer carrying ribonucleotides and two
deoxythymidines is chimeric: the RNA alphabet has no `T` and the DNA alphabet has no `U`, so
neither can hold it. Forcing one would fail validation or silently convert the overhang and record
a chemistry nobody chose — and `Nucleic` exists precisely to stop that class of error. `OrderStrand`
therefore keeps the core as an RNA `NucleicSequence` and the overhang as a separate, explicitly
deoxy value, joined only for display. The export emits core, overhang and joined form as separate
fields, because vendor notation for a chimeric oligo differs between suppliers and a researcher
must be able to reformat without re-deriving the chemistry.

**`export.order.csv` is separate from `export.csv`.** ADR 0015 made the core CSV deliberately a
core-sequence export and not order chemistry. Rather than change what that file means, the order
duplex gets its own route. The core CSV still answers "what did the engine evaluate"; the order CSV
answers "what do I buy".

**The sanity report travels in the same row as the sequences.** Five checks ship per duplex: core
length, strand complementarity, whether the named model expects the strand the duplex is designed
around, whether the parameter sets agree, and off-target screening. Putting them in a sibling file
would let a row be pasted into a vendor form without them. The cost is a wide, repetitive column;
the alternative is caveats that detach from the thing they qualify.

**Off-target screening reports `not-performed`, on every duplex, forever until it exists.** This
project implements no off-target analysis. Omitting the row would read as a clean screen, and the
engine's own three-outcome rule — missing data is not a negative assay — applies here as much as to
a rule verdict. Accepted constraint 15 asks for off-target findings in this report; this records
that the requirement is *not* met, rather than quietly shipping without it. Promoter-level
off-target analysis remains the prerequisite it always was.

**Findings report; they never withhold.** A concern never removes a duplex from the file. A
researcher who disagrees with a model's strand prediction may still order the duplex.

**Order rows are computed in Python and stored in the run record.** A frozen archive has no
backend, so its download must come from somewhere. Re-deriving the chemistry in JavaScript would
breach accepted constraint 11 and create two implementations that could name different molecules
for the same candidate — the ordering equivalent of the two-renderer risk in ADR 0015. Python
computes the rows once at run completion; JavaScript filters them by candidate id and joins
strings. A test asserts the stored rows and the live route's output are identical.

**Both CSV exports are written in rank order.** This changes existing behaviour. A snapshot's
candidates are in *position* order carrying their ranks, which is correct for the promoter map and
for a table the browser sorts for itself. A CSV is neither. Measured on cached SCN1A, `export.csv`
emitted ranks `8, 3, 7, 5, 13, 4, 2, 1, …`: a researcher taking the top three rows ordered the
wrong three molecules, with the rank column sitting right there contradicting them. Position
ordering remains correct in the record; only the flat exports sort.

**Rejected candidates are excluded by default from the order export, but orderable on request.**
A duplex that failed a gate is not something to buy. Naming its id explicitly still exports it,
because ADR 0007's distinction between "failed", "could not be checked" and "omitted" is the
researcher's to act on, not the exporter's to pre-empt.

## Costs and limits

The sanity report is verbose and repeats per row. The order CSV adds a second export whose
relationship to the first has to be explained in the interface. Storing order rows in the record
grows every snapshot slightly, and adds a field that a future schema change must carry.

**This export does not establish that a duplex works, is specific, or is safe.** It states the
molecule, the chemistry's source, and what was not checked. The strand-prediction checks report a
named model's expectation about duplex loading — not the identity of a promoter-associated target
RNA, and not activation. Ranking has not demonstrated performance above chance on the chapter 4
dataset, and nothing in an order file changes that.

Supported order chemistry remains 19 nt cores plus dTdT. A core of another length exports with a
concern naming the supported length rather than being silently refused or silently blessed.
