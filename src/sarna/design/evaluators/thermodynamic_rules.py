"""Evaluators for the duplex asymmetry gate and its tie-break.

Both read the asymmetry between the target's two ends, computed by
:mod:`sarna.thermo.guide`. The gate asks whether the asymmetry is large enough and in the right
direction; the tie-break reports its magnitude for ordering candidates that are otherwise level
(ADR 0006).

**Units are converted, not assumed.** The threshold parameter is stated in kcal/mol, but the
legacy table's values are in kJ/mol. Comparing a kcal threshold against a kJ difference would be
wrong by a factor of 4.184 — invisible at the default threshold of 0.0, where only the sign
matters, and silently wrong for any other value. The threshold is therefore converted into
whichever unit the selected table uses.

The gate uses the **one** table the active preset selects, so a preset stays a single reproducible
configuration. Reporting several tables' verdicts side by side is a separate concern handled by
:func:`sarna.thermo.guide.compare_asymmetry`, and it never changes which candidates pass
(ADR 0009).
"""

from sarna.design.context import CandidateContext
from sarna.rules.parameters import ParameterValue
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult
from sarna.thermo.guide import GuideStrand, assess_asymmetry
from sarna.thermo.model import NearestNeighbourTable
from sarna.thermo.tables import KJ_PER_KCAL, TABLES

TABLE_PARAMETER = "nearest-neighbour-table"


def _selected_table(resolved: ResolvedRule) -> NearestNeighbourTable:
    """Return the table the rule's parameters select.

    Raises:
        KeyError: If the parameter names a table that does not ship. The parameter's own
            ``choices`` should prevent this, so it indicates the two have drifted apart.
    """
    table_id = resolved.text(TABLE_PARAMETER)
    try:
        return TABLES[table_id]
    except KeyError:
        raise KeyError(
            f"rule {resolved.rule.id!r} selects nearest-neighbour table {table_id!r}, which "
            f"is not one of {', '.join(sorted(TABLES))}"
        ) from None


def _threshold_in(unit: str, kcal_threshold: float) -> float:
    """Convert a kcal/mol threshold into ``unit``.

    Raises:
        ValueError: If the unit is not one this function knows how to convert to, rather than
            comparing two incompatible quantities.
    """
    if unit == "kcal/mol":
        return kcal_threshold
    if unit == "kJ/mol":
        return kcal_threshold * KJ_PER_KCAL
    raise ValueError(
        f"cannot compare a kcal/mol threshold against a table measured in {unit!r}"
    )


def thermodynamic_asymmetry(
    context: CandidateContext, resolved: ResolvedRule
) -> RuleResult:
    """Check the guide strand's 5' end is the less stable by at least the required margin."""
    table = _selected_table(resolved)
    end_length = resolved.integer("end-length-nt")
    assessment = assess_asymmetry(context.target, table, end_length=end_length)
    threshold = _threshold_in(table.unit, resolved.number("minimum-asymmetry-kcal"))
    measurements: dict[str, ParameterValue] = {
        "asymmetry": assessment.difference,
        "asymmetry-unit": table.unit,
        "five-prime-energy": assessment.five_prime.energy,
        "three-prime-energy": assessment.three_prime.energy,
        "predicted-guide": str(assessment.predicted_guide),
        "nearest-neighbour-table": table.id,
        "end-length-nt": end_length,
    }
    if assessment.predicted_guide is GuideStrand.ANTISENSE and (
        assessment.difference >= threshold
    ):
        return RuleResult.passed(resolved, assessment.explanation, measurements)
    if assessment.predicted_guide is GuideStrand.ANTISENSE:
        return RuleResult.failed(
            resolved,
            f"{assessment.explanation} The margin of "
            f"{assessment.difference:+.2f} {table.unit} is below the "
            f"{threshold:.2f} {table.unit} required, so the strand selection is not confident "
            f"enough.",
            measurements,
        )
    return RuleResult.failed(resolved, assessment.explanation, measurements)


def asymmetry_margin(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Report the asymmetry magnitude candidates are ordered by when otherwise tied.

    The table is this rule's own parameter rather than a constant, so the value a run ordered by
    is recorded. It should normally match the mandatory asymmetry rule's table: admitting
    candidates with one parameter set and ordering them with another would be incoherent.
    """
    table = _selected_table(resolved)
    end_length = resolved.integer("end-length-nt")
    assessment = assess_asymmetry(context.target, table, end_length=end_length)
    return RuleResult.tiebreak(
        resolved,
        assessment.difference,
        f"Asymmetry {assessment.difference:+.2f} {table.unit} "
        f"({table.name}). Larger is a more confident strand selection.",
        {
            "asymmetry": assessment.difference,
            "asymmetry-unit": table.unit,
            "nearest-neighbour-table": table.id,
            "end-length-nt": end_length,
        },
    )
