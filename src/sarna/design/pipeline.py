"""The design run: promoter plus configuration in, evaluated and ranked candidates out.

This is the engine's public entry point, and it is deliberately a plain function over plain
objects — no HTTP, no UI, no global state (ADR 0005), so it can be driven from a notebook, a CLI
or a test as easily as from the API.

The run has four stages and keeps all four populations:

1. **Generate** every window of the configured target length, from the first base. The original
   script starts at index 3 and discards three windows with no stated reason (defect B6).
2. **Evaluate** every selected rule against each candidate.
3. **Gate** on the mandatory rules. A rejected candidate is kept, with the rules that rejected
   it, because a sequence that vanishes without explanation is what this project replaces.
4. **Rank** the survivors by preference count, breaking ties on asymmetry and then position
   (ADR 0006).

Promoter-level work happens once: CpG islands are detected for the whole promoter before the scan,
and only if the configuration selects the island rule. Detecting them per candidate would make the
scan quadratic; detecting them when the rule is unselected would impose a definition the
configuration never chose.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from sarna.design.candidate import Candidate
from sarna.design.config import DesignConfig
from sarna.design.context import CandidateContext
from sarna.design.evaluators import evaluate
from sarna.genome.promoter import Promoter
from sarna.sequence.cpg import CpGIsland, find_cpg_islands
from sarna.thermo.guide import AsymmetryComparison, compare_asymmetry
from sarna.thermo.model import NearestNeighbourTable
from sarna.thermo.tables import DNA_SANTALUCIA_1998, LEGACY_DNA_KJ, RNA_XIA_1998

DEFAULT_COMPARISON_TABLES: tuple[NearestNeighbourTable, ...] = (
    RNA_XIA_1998,
    DNA_SANTALUCIA_1998,
    LEGACY_DNA_KJ,
)
"""Tables whose guide predictions are reported side by side for every candidate (ADR 0009).

All three by default, because the most useful thing the tool can say about a legacy design is
that the corrected parameters disagree with it. These are **reported**, not applied: the
mandatory asymmetry gate uses only the one table the configuration selects, so the comparison
never changes which candidates pass.
"""


@dataclass(frozen=True, slots=True)
class DesignResult:
    """Everything one design run produced, including what it rejected and why."""

    promoter: Promoter
    config: DesignConfig
    generated: tuple[Candidate, ...]
    """Every window scanned, in position order. Includes rejected candidates."""

    ranked: tuple[Candidate, ...]
    """The candidates that passed every gate, best first, each carrying its rank."""

    cpg_islands: tuple[CpGIsland, ...] | None
    """Islands detected across the promoter, or ``None`` if detection was not run."""

    @property
    def rejected(self) -> tuple[Candidate, ...]:
        """Return the candidates a mandatory rule excluded, in position order."""
        return tuple(c for c in self.generated if c.is_rejected)

    @property
    def passing(self) -> tuple[Candidate, ...]:
        """Return the candidates that failed no mandatory rule, in position order."""
        return tuple(c for c in self.generated if not c.is_rejected)

    @property
    def unevaluated_rule_ids(self) -> tuple[str, ...]:
        """Return rules that could not be evaluated for any candidate.

        Reported at the run level so a user is told once that, say, the TSS rules did not apply,
        rather than discovering it per candidate or not at all.
        """
        if not self.generated:
            return ()
        unevaluated = {
            result.rule_id
            for candidate in self.generated
            for result in candidate.unevaluated
        }
        return tuple(
            resolved.rule.id
            for resolved in self.config
            if resolved.rule.id in unevaluated
            and all(
                candidate.result(resolved.rule.id).is_unevaluated
                for candidate in self.generated
            )
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the whole run as JSON-serialisable data, for a reproducible export."""
        return {
            "configuration": self.config.to_dict(),
            "promoter": {
                "length": self.promoter.length,
                "tss_start": (
                    self.promoter.tss_start.value if self.promoter.tss_start else None
                ),
                "assembly": (
                    str(self.promoter.interval.assembly) if self.promoter.interval else None
                ),
                "strand": (
                    str(self.promoter.interval.strand) if self.promoter.interval else None
                ),
            },
            "counts": {
                "generated": len(self.generated),
                "rejected": len(self.rejected),
                "passing": len(self.passing),
            },
            "unevaluated_rules": list(self.unevaluated_rule_ids),
            "cpg_islands": (
                [
                    {
                        "start": island.start,
                        "end": island.end,
                        "gc_fraction": island.gc_fraction,
                        "observed_expected_ratio": island.observed_expected_ratio,
                    }
                    for island in self.cpg_islands
                ]
                if self.cpg_islands is not None
                else None
            ),
            "candidates": [candidate.to_dict() for candidate in self.generated],
        }


def design_sarnas(
    promoter: Promoter,
    config: DesignConfig,
    comparison_tables: Sequence[NearestNeighbourTable] = DEFAULT_COMPARISON_TABLES,
) -> DesignResult:
    """Design, evaluate and rank saRNA candidates across a promoter.

    Args:
        promoter: The region to scan, with whatever anchoring is known. A promoter with no
            coordinates and no TSS still yields candidates; the rules that need them report
            not-applicable rather than failing (ADR 0007).
        config: The resolved rules to apply, normally from a preset.
        comparison_tables: Nearest-neighbour tables whose guide predictions are reported beside
            each candidate. Reported only — the gate uses the configuration's table.

    Returns:
        The run's candidates, rankings and the configuration that produced them.

    Raises:
        KeyError: If the configuration does not define a target length, without which there is
            nothing to generate.
        ValueError: If the configuration is internally inconsistent.
    """
    target_length = config.target_length
    criteria = config.cpg_island_criteria
    islands = (
        find_cpg_islands(promoter.sequence, criteria) if criteria is not None else None
    )

    generated = tuple(
        _build(promoter, offset, target_length, config, islands, tuple(comparison_tables))
        for offset, _ in promoter.sequence.windows(target_length)
    )
    return DesignResult(
        promoter=promoter,
        config=config,
        generated=generated,
        ranked=_rank(tuple(c for c in generated if not c.is_rejected)),
        cpg_islands=islands,
    )


def _build(
    promoter: Promoter,
    offset: int,
    target_length: int,
    config: DesignConfig,
    islands: tuple[CpGIsland, ...] | None,
    comparison_tables: tuple[NearestNeighbourTable, ...],
) -> Candidate:
    """Evaluate one window and assemble its candidate."""
    context = CandidateContext.at(
        promoter, offset=offset, length=target_length, cpg_islands=islands
    )
    results = tuple(evaluate(context, resolved) for resolved in config)
    duplex = context.target.as_rna()
    return Candidate(
        # 1-based and zero-padded, so ids sort lexically in the order they were generated.
        id=f"saRNA-{offset + 1:03d}",
        offset=offset,
        target=context.target,
        sense_strand=duplex,
        antisense_strand=duplex.reverse_complement(),
        tss_span=context.tss_span,
        interval=context.interval,
        results=results,
        guide_comparison=_compare(context, comparison_tables, config),
        max_preference_points=config.max_preference_points,
    )


def _compare(
    context: CandidateContext, tables: tuple[NearestNeighbourTable, ...], config: DesignConfig,
) -> AsymmetryComparison:
    """Compare tables using the selected asymmetry rule's end length, or the default.

    The gate takes precedence over the optional ranking tiebreak. Comparing four bases while
    evaluating an explicitly chosen different end length could report the opposite strand.
    """
    selected = next((config.rule(key) for key in (
        "thermodynamic-asymmetry", "asymmetry-margin"
    ) if key in config), None)
    return compare_asymmetry(
        context.target, tables,
        end_length=selected.integer("end-length-nt") if selected else 4,
    )


def _rank(passing: tuple[Candidate, ...]) -> tuple[Candidate, ...]:
    """Order the survivors and assign ranks.

    Preference count descending, then tie-break value descending, then position ascending.
    The last term exists only so that output is deterministic; candidates reaching it are
    genuinely indistinguishable on the evidence and the interface should present them as tied
    (ADR 0006).

    A candidate with no tie-break value sorts below one that has it at the same preference
    count, rather than being treated as having a value of zero: an uncomputed asymmetry is not
    a neutral asymmetry.
    """
    ordered = sorted(
        passing,
        key=lambda c: (
            -c.preference_points,
            0 if c.tiebreak_value is not None else 1,
            -(c.tiebreak_value if c.tiebreak_value is not None else 0.0),
            c.offset,
        ),
    )
    return tuple(
        candidate.with_rank(rank) for rank, candidate in enumerate(ordered, start=1)
    )
