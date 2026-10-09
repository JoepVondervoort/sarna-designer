"""The order-ready duplex: what a synthesis vendor is actually asked to make.

A :class:`~sarna.design.candidate.Candidate` carries a 19 nt RNA core, which is what every rule
is evaluated against. What gets *ordered* is two 21-mers: each core plus a 2 nt **dTdT** 3'
overhang. Chapter 2's supplementary material settles that the overhang is deoxythymidine,
written ``[dT][dT]``, rather than a ribo-``UU``; the legacy inventory's reading of ``UU`` is
superseded. See `docs/scientific/scn1a-design-coordinates.md` and ``GLOSSARY.md``.

**The ordered strand is deliberately not a** :class:`~sarna.sequence.nucleic.NucleicSequence`.
A 21-mer carrying ribonucleotides and two deoxythymidines is chimeric: the RNA alphabet has no
``T`` and the DNA alphabet has no ``U``, so neither can hold it. Forcing one would either fail
validation or, worse, silently convert the overhang and record a chemistry nobody chose. The core
therefore stays an RNA :class:`~sarna.sequence.nucleic.NucleicSequence` and the overhang is a
separate, explicitly deoxy value; :attr:`OrderStrand.annotated` joins them only for display.

The overhang is a synthesis choice. It does not change the core, and it is never scanned or
scored, no rule sees it.

**This is not a screened order.** Accepted constraint 15 requires the order export to carry a
sanity report including off-target findings, and no off-target analysis exists in this project
yet. :data:`OrderOutcome.NOT_PERFORMED` says so on every duplex rather than omitting the row,
because an absent row reads as a clean screen. Findings report; they never withhold sequences.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any

from sarna.design.candidate import Candidate
from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.guide import GuideStrand

SUPPORTED_CORE_LENGTH_NT = 19
"""The only core length this project's order chemistry is defined for (accepted constraint 4)."""

OVERHANG_BASES = "TT"
"""The overhang's bases. Deoxythymidine, hence ``T`` rather than the core's ``U``."""

OVERHANG_LENGTH_NT = len(OVERHANG_BASES)

OVERHANG_ANNOTATED = "[dT][dT]"
"""How chapter 2's supplementary material writes the overhang, and how this project writes it.

The ``d`` is the point of the notation: it distinguishes the deoxythymidines from a ribo-``UU``,
which is what the legacy scripts assumed and what a bare ``TT`` would leave ambiguous beside an
RNA core.
"""

ORDERED_LENGTH_NT = SUPPORTED_CORE_LENGTH_NT + OVERHANG_LENGTH_NT
"""21 nt: what is actually synthesised, per strand."""

DESIGNED_GUIDE = GuideStrand.ANTISENSE
"""The duplex is designed so the antisense strand is loaded. A model may predict otherwise."""

DEFAULT_MODEL = "rna-xia-1998"
"""The table Wang's preset gates on, and so the default one to compare the design intent against."""


class OrderOutcome(Enum):
    """The three answers a sanity check can give, mirroring the engine's three rule outcomes.

    A check that cannot be run is distinct from one that found nothing wrong. Collapsing them
    would let an unscreened duplex read as a screened one.
    """

    OK = "ok"
    """Checked, and nothing to report."""

    CONCERN = "concern"
    """Checked, and worth reading before spending money. Never excludes the duplex."""

    NOT_PERFORMED = "not-performed"
    """Not checked at all. Not a clean result."""

    def __str__(self) -> str:
        """Return the wire value."""
        return self.value


@dataclass(frozen=True, slots=True)
class OrderStrand:
    """One strand as it will be synthesised: an RNA core plus a deoxy 3' overhang."""

    role: str
    """``sense`` or ``antisense``. Sequence identity, not predicted fate."""

    core: NucleicSequence
    """The RNA core, 5'->3'. What the rules were evaluated against."""

    @property
    def overhang_bases(self) -> str:
        """Return the overhang's bases, ``TT``."""
        return OVERHANG_BASES

    @property
    def overhang_alphabet(self) -> Nucleic:
        """Return :attr:`Nucleic.DNA`, because the overhang is deoxythymidine."""
        return Nucleic.DNA

    @property
    def length_nt(self) -> int:
        """Return the synthesised length: core plus overhang."""
        return len(self.core) + OVERHANG_LENGTH_NT

    @property
    def annotated(self) -> str:
        """Return the strand as this project writes it, 5'->3'.

        The overhang keeps its ``[dT]`` notation so the chemistry survives a copy-paste into a
        vendor form, where a bare ``TT`` beside an RNA core would be read as ribo-uridine.
        """
        return f"{self.core.bases}{OVERHANG_ANNOTATED}"

    def to_dict(self) -> dict[str, Any]:
        """Export the strand, keeping core and overhang separately addressable.

        Vendor notation for a chimeric oligo differs between suppliers, so a researcher needs
        the parts, not only this project's rendering of them.
        """
        return {
            "role": self.role,
            "core_rna_5to3": self.core.bases,
            "core_length_nt": len(self.core),
            "overhang": OVERHANG_ANNOTATED,
            "overhang_bases": self.overhang_bases,
            "overhang_alphabet": str(self.overhang_alphabet),
            "annotated": self.annotated,
            "length_nt": self.length_nt,
        }


@dataclass(frozen=True, slots=True)
class OrderCheck:
    """One verdict about a duplex, with the reasoning that produced it."""

    id: str
    outcome: OrderOutcome
    detail: str

    def to_dict(self) -> dict[str, Any]:
        """Export the check."""
        return {"id": self.id, "outcome": str(self.outcome), "detail": self.detail}


@dataclass(frozen=True, slots=True)
class DuplexOrder:
    """An order-ready duplex and everything checked about it before synthesis."""

    candidate_id: str
    rank: int | None
    target_site: NucleicSequence
    """The genomic locus, as DNA. Not a delivered molecule."""

    sense: OrderStrand
    antisense: OrderStrand
    checks: tuple[OrderCheck, ...]

    @property
    def concerns(self) -> tuple[OrderCheck, ...]:
        """Return the checks that found something, in report order."""
        return tuple(c for c in self.checks if c.outcome is OrderOutcome.CONCERN)

    def to_dict(self) -> dict[str, Any]:
        """Export the whole order."""
        return {
            "candidate": self.candidate_id,
            "rank": self.rank,
            "target_site_dna": self.target_site.bases,
            "sense": self.sense.to_dict(),
            "antisense": self.antisense.to_dict(),
            "checks": [c.to_dict() for c in self.checks],
            "concern_count": len(self.concerns),
        }


def _core_length_check(core: NucleicSequence) -> OrderCheck:
    """Report whether the core is the one length this chemistry is defined for."""
    if len(core) == SUPPORTED_CORE_LENGTH_NT:
        return OrderCheck(
            id="core-length",
            outcome=OrderOutcome.OK,
            detail=(
                f"The core is {SUPPORTED_CORE_LENGTH_NT} nt, the length this project's order "
                f"chemistry is defined for."
            ),
        )
    return OrderCheck(
        id="core-length",
        outcome=OrderOutcome.CONCERN,
        detail=(
            f"The core is {len(core)} nt. This project's order chemistry is defined only for "
            f"{SUPPORTED_CORE_LENGTH_NT} nt plus dTdT, so a {len(core)} nt core is designable "
            f"but not covered by that contract. Position-specific rules also lose their meaning "
            f"away from 19 nt: position 19 stops being the 3' terminal base."
        ),
    )


def _complementarity_check(sense: NucleicSequence, antisense: NucleicSequence) -> OrderCheck:
    """Report whether the two cores actually anneal."""
    if antisense.bases == sense.reverse_complement().bases:
        return OrderCheck(
            id="strand-complementarity",
            outcome=OrderOutcome.OK,
            detail="The two cores are exact reverse complements and anneal over all 19 bases.",
        )
    return OrderCheck(
        id="strand-complementarity",
        outcome=OrderOutcome.CONCERN,
        detail=(
            "The antisense core is not the exact reverse complement of the sense core, so the "
            "duplex would not anneal as designed. This indicates a defect upstream of ordering."
        ),
    )


@dataclass(frozen=True, slots=True)
class StrandPrediction:
    """One parameter set's verdict, reduced to what ordering needs.

    Exists so the checks below can be built identically from a live
    :class:`~sarna.design.candidate.Candidate` and from a frozen snapshot record, rather than the
    two growing separate implementations that could disagree about the same duplex.
    """

    model: str
    guide: GuideStrand
    explanation: str


def _intent_check(predictions: tuple[StrandPrediction, ...], model: str) -> OrderCheck:
    """Report whether the named model expects the strand the duplex was designed to load.

    Raises:
        KeyError: If ``model`` named no assessed table. Reporting nothing would read as
            agreement.
    """
    prediction = next((p for p in predictions if p.model == model), None)
    if prediction is None:
        assessed = ", ".join(p.model for p in predictions) or "none"
        raise KeyError(
            f"{model!r} is not among this candidate's assessed tables ({assessed}), so no "
            f"strand prediction can be reported for it"
        )
    if prediction.guide is DESIGNED_GUIDE:
        return OrderCheck(
            id="guide-matches-intent",
            outcome=OrderOutcome.OK,
            detail=(
                f"{model} predicts the antisense strand is loaded, which is the strand this "
                f"duplex is designed around. {prediction.explanation}"
            ),
        )
    return OrderCheck(
        id="guide-matches-intent",
        outcome=OrderOutcome.CONCERN,
        detail=(
            f"{model} predicts the {prediction.guide} strand is loaded, not the "
            f"antisense strand this duplex is designed around. The duplex may act through the "
            f"strand you did not intend. {prediction.explanation}"
        ),
    )


def _agreement_check(tables_agree: bool, summary: str) -> OrderCheck:
    """Report whether the parameter sets agree about which strand is loaded."""
    if tables_agree:
        return OrderCheck(
            id="model-agreement",
            outcome=OrderOutcome.OK,
            detail=f"Every assessed parameter set agrees. {summary}",
        )
    return OrderCheck(
        id="model-agreement",
        outcome=OrderOutcome.CONCERN,
        detail=(
            f"The parameter sets disagree about which strand is loaded, which means the "
            f"asymmetry is small enough for the choice of table to decide the answer, itself a "
            f"reason to trust the prediction less. This divides the tables for 13.72% of "
            f"possible inputs (docs/scientific/guide-prediction-divergence.md). {summary}"
        ),
    )


@dataclass(frozen=True, slots=True)
class ScreeningFinding:
    """One screening's result for one candidate, with the provenance that makes it readable.

    A screening is only meaningful with its database and its date: NCBI's holdings change week
    to week, so "nothing found" is a statement about one database on one day (ADR 0019). Both
    travel in the row rather than in a separate file, because a row pasted into a vendor form
    must not arrive without the caveats belonging to it.
    """

    scope: str
    database: str
    searched_at: str
    row: Mapping[str, Any]
    """The candidate's line from the specificity table: counts, never a combined score."""

    @property
    def is_concern(self) -> bool:
        """Whether anything here is worth reading before spending money.

        Deliberately a disjunction of named counts rather than a threshold on a score. ADR 0019
        refuses a combined specificity figure because weighing one perfect match elsewhere
        against forty weak partial ones is the researcher's judgement, not the tool's.
        """
        return any(
            int(self.row.get(key) or 0) > 0
            for key in ("perfect_hits", "near_perfect_hits", "in_run_window_hits")
        )

    def describe(self) -> str:
        """Return this screening's findings as one clause."""
        counts = ", ".join(
            f"{label} {int(self.row.get(key) or 0)}"
            for key, label in (
                ("hit_count", "alignments"),
                ("perfect_hits", "perfect"),
                ("near_perfect_hits", "near-perfect"),
                ("in_run_window_hits", "inside this run's window"),
            )
            if key != "in_run_window_hits" or self.scope == GENOMIC_SCOPE
        )
        return f"{self.scope} ({self.database}, {self.searched_at[:10]}): {counts}"


GENOMIC_SCOPE = "genomic"
"""Mirrors :data:`sarna.offtarget.models.GENOMIC_SCOPE` without importing it.

``sarna.design`` is the engine and must not depend on the screening package; the dependency
runs the other way. The string is the wire value both sides already agree on.
"""


def _screening_check(findings: Sequence[ScreeningFinding] = ()) -> OrderCheck:
    """Report what an off-target screening found, or that none was performed.

    Three outcomes, as everywhere else in this project. No screening is **not** a clean one,
    which is why the absence keeps its own outcome rather than collapsing into ``OK``.
    """
    if not findings:
        return OrderCheck(
            id="off-target-screening",
            outcome=OrderOutcome.NOT_PERFORMED,
            detail=(
                "No off-target screening was performed for this run. Nothing here establishes "
                "that these sequences are specific to the intended locus. Screen them before "
                "ordering, or run the screening in the workbench."
            ),
        )
    described = "; ".join(f.describe() for f in findings)
    if any(f.is_concern for f in findings):
        return OrderCheck(
            id="off-target-screening",
            outcome=OrderOutcome.CONCERN,
            detail=(
                f"Off-target screening found matches worth reading before ordering, "
                f"{described}. Sequence similarity is not measured binding, and nothing was "
                f"withheld: the duplex is above and the alignments are in the screening record."
            ),
        )
    return OrderCheck(
        id="off-target-screening",
        outcome=OrderOutcome.OK,
        detail=(
            f"Off-target screening found nothing of note, {described}. One search, of one "
            f"database, on one day, with parameters that can miss things: this does not "
            f"establish specificity."
        ),
    )


def _build(
    *,
    candidate_id: str,
    rank: int | None,
    target_site: NucleicSequence,
    sense_core: NucleicSequence,
    antisense_core: NucleicSequence,
    predictions: tuple[StrandPrediction, ...],
    tables_agree: bool,
    summary: str,
    model: str,
    screenings: Sequence[ScreeningFinding] = (),
) -> DuplexOrder:
    """Assemble a duplex order and its checks. The one place either entry point goes through."""
    return DuplexOrder(
        candidate_id=candidate_id,
        rank=rank,
        target_site=target_site,
        sense=OrderStrand(role="sense", core=sense_core),
        antisense=OrderStrand(role="antisense", core=antisense_core),
        checks=(
            _core_length_check(sense_core),
            _complementarity_check(sense_core, antisense_core),
            _intent_check(predictions, model),
            _agreement_check(tables_agree, summary),
            _screening_check(screenings),
        ),
    )


def order_for(
    candidate: Candidate, model: str = DEFAULT_MODEL,
    screenings: Sequence[ScreeningFinding] = (),
) -> DuplexOrder:
    """Build the order-ready duplex for ``candidate``, with its sanity report.

    Args:
        candidate: Any generated candidate. A rejected one still produces an order, because
            excluding it here would hide the reason it was rejected behind a missing row;
            its rank is ``None`` and its rule verdicts remain the record of why.
        model: The nearest-neighbour table whose strand prediction is compared against the
            design intent. Defaults to the RNA table, which is the one Wang's preset gates on.

        screenings: This candidate's findings from any completed off-target screening,
            one per scope. Empty means none was performed, which is not a clean
            result and is reported as its own outcome.

    Returns:
        The two 21-mers and every check made on them.

    Raises:
        KeyError: If ``model`` named no table assessed for this candidate.
    """
    comparison = candidate.guide_comparison
    return _build(
        candidate_id=candidate.id,
        rank=candidate.rank,
        target_site=candidate.target,
        sense_core=candidate.sense_strand,
        antisense_core=candidate.antisense_strand,
        predictions=tuple(
            StrandPrediction(
                model=a.table_id, guide=a.predicted_guide, explanation=a.explanation
            )
            for a in comparison.assessments
        ),
        tables_agree=comparison.tables_agree,
        summary=comparison.summary,
        model=model,
        screenings=screenings,
    )


def order_from_record(
    record: dict[str, Any], model: str = DEFAULT_MODEL,
    screenings: Sequence[ScreeningFinding] = (),
) -> DuplexOrder:
    """Build the order-ready duplex from an exported candidate record.

    A saved run is stored as records, not as :class:`~sarna.design.candidate.Candidate` objects,
    and reopening one must produce the same duplex and the same checks as the live run did. This
    goes through :func:`_build` for exactly that reason.

    Args:
        record: One entry of a run snapshot's ``candidates``.
        model: As for :func:`order_for`.
        screenings: As for :func:`order_for`.

    Returns:
        The two 21-mers and every check made on them.

    Raises:
        KeyError: If the record is missing a required field, or ``model`` named no assessed
            table.
    """
    prediction = record["guide_prediction"]
    return _build(
        candidate_id=record["id"],
        rank=record["rank"],
        target_site=NucleicSequence.parse(record["target_site_dna"], Nucleic.DNA),
        sense_core=NucleicSequence.parse(record["sense_strand_rna"], Nucleic.RNA),
        antisense_core=NucleicSequence.parse(record["antisense_strand_rna"], Nucleic.RNA),
        predictions=tuple(
            StrandPrediction(
                model=a["table"],
                guide=GuideStrand(a["predicted_guide"]),
                explanation=a["explanation"],
            )
            for a in prediction["assessments"]
        ),
        tables_agree=prediction["tables_agree"],
        summary=prediction["summary"],
        model=model,
        screenings=screenings,
    )


CSV_COLUMNS = (
    "run_id", "candidate", "rank", "target_site_dna_5to3",
    "sense_order_5to3", "antisense_order_5to3",
    "sense_core_rna_5to3", "antisense_core_rna_5to3", "overhang",
    "ordered_length_nt", "core_length_nt", "concerns", "sanity_report",
)
"""Columns of the order CSV. The ordered strands come first, then the parts they are built from."""


def csv_row(order: DuplexOrder, run_id: str) -> list[str]:
    """Return ``order`` as one CSV row, with its whole sanity report in the last column.

    The report travels with the sequences rather than in a separate file, so a row pasted into a
    vendor form cannot arrive without the caveats that belong to it.
    """
    report = " | ".join(f"{c.id} [{c.outcome}]: {c.detail}" for c in order.checks)
    return [
        run_id, order.candidate_id, "" if order.rank is None else str(order.rank),
        order.target_site.bases,
        order.sense.annotated, order.antisense.annotated,
        order.sense.core.bases, order.antisense.core.bases, OVERHANG_ANNOTATED,
        str(order.sense.length_nt), str(len(order.sense.core)),
        str(len(order.concerns)), report,
    ]
