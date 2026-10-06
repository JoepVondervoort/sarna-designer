"""The project's first empirical dataset: saRNAs with measured activity.

Until this existed, "correct" in this repository meant "matches Wang et al. 2015 and matches the
legacy script". Nothing connected a design rule to a measured outcome. This module loads 125
target sequences whose activity was assayed, 14 of which activated their gene, and so it is the
first thing here that can make a rule *wrong* rather than merely unconventional.

What the set will and will not support is worth stating plainly, because it is easy to read more
into 125 labelled sequences than they carry.

**It cannot measure the mandatory gate.** Every one of the 125 was designed to satisfy Wang's
structural criteria, so the set is range-restricted by construction: it contains almost no
candidate the gate would reject. Asking "does the gate separate activators from non-activators?"
of this set asks about variation it does not contain. What it can test is whether anything
*beyond* the gate — ranking, preferences, thermodynamics — separates the two classes.

**"Not activating" is a measurement, not a property.** The 111 negatives failed to show
reproducible upregulation in HeLa or AF22 at 50 nM, 72 h, lipofected, by qPCR. A sequence that
does nothing under those conditions may still act in another cell type, at another dose, or on a
phenotype qPCR does not see. Chapter 4's own finding is that the *cell* changed the answer more
than the sequence did: AF22 produced no activators for any locus, including the P21 control that
works reliably in HeLa.

**The positives are not one experiment.** Six come from published work on six different genes,
cell types and assays; three are SCN1A duplexes from chapters 2-3; five are from chapter 4's own
HeLa screen, where several were reproducible but did not reach significance. The activation
threshold is >= 1.5-fold, which is the chapter's criterion, not a universal one.

**Sequences are target/sense strands**, 5'->3', so Wang's Table 1 positions apply directly: see
:func:`load_chapter4_dataset` for how that was established, since the chapter's own text calls
them guide strands and the distinction changes what every position means.
"""

import csv
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from enum import Enum
from importlib import resources
from typing import Any

from sarna.rules.provenance import Provenance
from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence

DATA_PACKAGE = "sarna.validation.data"
CHAPTER4_FILENAME = "vondervoort_2026_chapter4.csv"

ACTIVATION_THRESHOLD_FOLD = 1.5
"""Fold change over scrambled control at which chapter 4 called a design an activator."""


class Activity(Enum):
    """Whether a duplex raised its target's transcript level in the assay that tested it."""

    ACTIVATING = "activating"
    """Reproducible upregulation at or above 1.5-fold over scrambled control."""

    NOT_ACTIVATING = "not-activating"
    """No reproducible upregulation under the conditions tested. Not 'inert' in general."""

    def __str__(self) -> str:
        """Return the label, e.g. ``activating``."""
        return self.value


@dataclass(frozen=True, slots=True)
class ValidationSequence:
    """One assayed duplex: its target site and what the assay measured."""

    name: str
    """Identifier as published, e.g. ``CHD2-8``. Unique within a dataset."""

    gene: str
    """The gene whose promoter was targeted."""

    target: NucleicSequence
    """The target site as RNA, 5'->3'. Equivalently the duplex's sense strand."""

    activity: Activity
    evidence: str
    """Where the measurement comes from, as a citation a reader can check."""

    panel: int
    """Which screening panel the sequence came from, disambiguating reused names."""

    @property
    def is_activating(self) -> bool:
        """Return whether this duplex activated its target."""
        return self.activity is Activity.ACTIVATING

    @property
    def target_dna(self) -> NucleicSequence:
        """Return the target site as DNA, which is what the design rules are written about."""
        return self.target.as_dna()

    @property
    def guide(self) -> NucleicSequence:
        """Return the antisense strand, the strand intended for AGO2 loading."""
        return self.target.reverse_complement()


@dataclass(frozen=True, slots=True)
class ValidationDataset:
    """A labelled collection of assayed duplexes, with its provenance and its caveats."""

    id: str
    name: str
    sequences: tuple[ValidationSequence, ...]
    provenance: Provenance
    caveats: tuple[str, ...]
    """Why this set cannot answer certain questions. Carried so an export states its own limits."""

    def __len__(self) -> int:
        """Return how many sequences the set holds."""
        return len(self.sequences)

    def __iter__(self) -> Iterator[ValidationSequence]:
        """Iterate the sequences in file order."""
        return iter(self.sequences)

    @property
    def activating(self) -> tuple[ValidationSequence, ...]:
        """Return the duplexes that activated their target."""
        return tuple(s for s in self.sequences if s.is_activating)

    @property
    def not_activating(self) -> tuple[ValidationSequence, ...]:
        """Return the duplexes that did not."""
        return tuple(s for s in self.sequences if not s.is_activating)

    @property
    def genes(self) -> tuple[str, ...]:
        """Return the genes represented, in alphabetical order."""
        return tuple(sorted({s.gene for s in self.sequences}))

    def by_name(self, name: str) -> ValidationSequence:
        """Return the sequence called ``name``.

        Raises:
            KeyError: If the dataset holds no such sequence.
        """
        for sequence in self.sequences:
            if sequence.name == name:
                return sequence
        raise KeyError(f"dataset {self.id!r} has no sequence named {name!r}")

    def to_dict(self) -> dict[str, Any]:
        """Return the set as JSON-serialisable data, caveats included."""
        return {
            "id": self.id,
            "name": self.name,
            "provenance": str(self.provenance),
            "counts": {
                "total": len(self),
                "activating": len(self.activating),
                "not_activating": len(self.not_activating),
            },
            "genes": list(self.genes),
            "caveats": list(self.caveats),
        }


CHAPTER4_CAVEATS: tuple[str, ...] = (
    "All 125 sequences were designed to satisfy Wang et al. 2015 structural criteria, so the "
    "set is range-restricted and cannot measure how well the mandatory gate discriminates. It "
    "can only test whether anything beyond the gate separates the two classes.",
    "'Not activating' means no reproducible upregulation in HeLa or AF22 at 50 nM, 72 h, "
    "lipofected, measured by qPCR against GAPDH and GUSB. It is a negative result under those "
    "conditions, not evidence that a sequence is inert.",
    "Cell type dominated the outcome: AF22 neural stem cells produced no activators at any "
    "locus, including the P21 positive control. A sequence-only model cannot explain that, so "
    "no sequence rule fitted to this set should be expected to transfer across cell types.",
    "The 14 activators pool six published duplexes against six genes in six cell contexts, "
    "three SCN1A duplexes from chapters 2-3, and five from chapter 4's own HeLa screen, where "
    "several were reproducible but did not reach statistical significance.",
    "The appendix reuses the names GAD1-1 to GAD1-6 for two different panels of sequences. The "
    "second panel is suffixed 'b' here. Which gene its 12 sequences target is not stated in the "
    "chapter and has not been confirmed.",
    "GAD1-5b and GAD1-6b are the same sequence under two names, so the 111 negatives contain "
    "110 distinct sequences.",
)


def load_chapter4_dataset() -> ValidationDataset:
    """Load the chapter 4 screening dataset.

    The sequences are stored as **target/sense** strands, 5'->3'. The chapter's appendix labels
    the column only "Sequence", and its motif analysis calls these guide strands, so the strand
    was established from the data rather than assumed: 90% carry A at position 19 and 98% carry
    A or U at position 18, which are Wang Table 1 criteria written about the *target*. Read as
    guide strands the same criteria would appear reverse-complemented, and only 12% begin with U.
    External confirmation was attempted and not obtained: Li et al. 2006 place their dsRNA
    sequences in PNAS supplementary Table 1, which was not retrievable, so the strand assignment
    rests on the positional evidence above rather than on a published sequence match.

    This matters beyond bookkeeping: position 15 of a target is position 5 of its guide, so
    which strand the positions are counted on decides which end of which molecule a positional
    finding is about.

    Returns:
        The dataset, with its provenance and the caveats that bound what it can show.
    """
    source = resources.files(DATA_PACKAGE).joinpath(CHAPTER4_FILENAME)
    with resources.as_file(source) as path, path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return ValidationDataset(
        id="vondervoort-2026-chapter4",
        name="Promoter-targeting saRNA screen across seven genes (thesis chapter 4)",
        sequences=tuple(_sequence_from(row) for row in rows),
        provenance=Provenance.citation(
            "van de Vondervoort 2026",
            "Evaluating saRNA-Mediated Promoter Activation Beyond SCN1A, chapter 4, "
            "appendix tables S1-S2",
        ),
        caveats=CHAPTER4_CAVEATS,
    )


def _sequence_from(row: dict[str, str]) -> ValidationSequence:
    """Build one sequence from a CSV row, validating the sequence as it goes."""
    return ValidationSequence(
        name=row["name"],
        gene=row["gene"],
        target=NucleicSequence.parse(row["target_rna_5to3"], Nucleic.RNA),
        activity=Activity(row["activity"]),
        evidence=row["evidence"],
        panel=int(row["panel"]),
    )


def target_dna_sequences(dataset: Sequence[ValidationSequence]) -> tuple[NucleicSequence, ...]:
    """Return the target sites as DNA, which is the form the design rules are written about."""
    return tuple(sequence.target_dna for sequence in dataset)


__all__ = [
    "ACTIVATION_THRESHOLD_FOLD",
    "CHAPTER4_CAVEATS",
    "Activity",
    "ValidationDataset",
    "ValidationSequence",
    "load_chapter4_dataset",
    "target_dna_sequences",
]
