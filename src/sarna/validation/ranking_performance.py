"""How well each preset's ranking actually separates activity — measured, for the interface.

The interface shows a rating next to every candidate. A reader will take that as a prediction
unless told otherwise, so the tool states its own measured performance beside it: on the 125
assayed sequences of thesis chapter 4, **no preset's ranking beats chance**.

**The numbers here are computed, not transcribed.** `performance_for_preset` runs the real
benchmark over the bundled dataset and caches the result. A constant copied out of a document
goes stale the first time a threshold changes and nobody re-runs the measurement; a computed one
cannot. The full benchmark takes well under a second per preset, so there is no reason to copy
it.

What the figure does *not* mean is as important as what it does, and both travel with it in
:func:`performance_summary`:

- The dataset is **range-restricted**. Every sequence in it was built to pass Wang's gate, so it
  can show a gate that is too strict and never one that is too permissive.
- Four of seven preferences cannot be evaluated without promoter coordinates, so on this set two
  preferences are doing nearly all the work. A near-chance AUC says as much about the dataset as
  about the ranking.
- `exploratory` scores highest **because two of its rules were fitted to this very data**.
  Presenting its number without that would make the least-evidenced preset look the best.
- A figure belongs to a *preset*. A researcher who adjusts rules is no longer running the
  configuration that was measured, and :func:`performance_summary` is not theirs to claim.

See `docs/scientific/empirical-validation.md` for the method and
[ADR 0006](../../../docs/adr/0006-ranking-counts-criteria.md) for why ranking counts criteria
rather than weighting them.
"""

from functools import cache
from typing import Any

from sarna.design.presets import PRESETS
from sarna.validation.benchmark import UNINFORMATIVE_AUC, BenchmarkResult, benchmark
from sarna.validation.dataset import load_chapter4_dataset

CHANCE_AUC = UNINFORMATIVE_AUC
"""0.5 — an activator outranks a non-activator exactly half the time."""

CIRCULAR_PRESETS = frozenset({"exploratory"})
"""Presets containing rules derived from the dataset they are measured against.

`pos10-g` and `central-stability` are the two best-separating rules in the whole catalogue, and
both were derived from these 125 sequences. Their apparent performance is circular, and the
tests assert it so the circularity stays visible rather than becoming a selling point.
"""

CIRCULARITY_NOTE = (
    "This configuration scores highest because two of its rules — pos10-g and "
    "central-stability — were fitted to this very dataset. The figure is circular and is not "
    "evidence that the configuration predicts activity."
)

RANGE_RESTRICTION = (
    "The dataset is range-restricted: every sequence in it was designed to pass Wang's gate, so "
    "it can show that a gate is too strict but never that one is too permissive. Four of the "
    "seven preferences also cannot be evaluated without promoter coordinates, so on this set a "
    "near-chance result describes the dataset at least as much as the ranking."
)


@cache
def performance_for_preset(preset_id: str) -> BenchmarkResult | None:
    """Return the measured benchmark for ``preset_id``, or ``None`` if there is no such preset.

    Cached, because the interface asks for this on every run and the answer never changes within
    a process. ``None`` rather than a default: attaching another preset's measurement, or a
    neutral 0.5, would report a number nobody measured for this configuration.
    """
    preset = PRESETS.get(preset_id)
    if preset is None:
        return None
    # resolve() carries preset_id into the config, and benchmark() reads it from there.
    return benchmark(load_chapter4_dataset(), preset.resolve())


def performance_summary(preset_id: str) -> dict[str, Any] | None:
    """Return the measured performance as the interface displays it.

    Returns ``None`` for an unrecognised preset, so a caller cannot accidentally render a
    measurement that does not exist.
    """
    result = performance_for_preset(preset_id)
    if result is None:
        return None
    dataset = load_chapter4_dataset()
    activating = len(result.activating)
    circular = preset_id in CIRCULAR_PRESETS
    return {
        "preset": preset_id,
        "auc": round(result.auc, 3),
        "auc_p_value": round(result.auc_p_value, 3),
        "chance_auc": CHANCE_AUC,
        "beats_chance": result.discriminates,
        "headline": _headline(result),
        "detail": result.summary,
        "dataset": {
            "name": dataset.name,
            "assayed": len(result.scored),
            "activating": activating,
            "activators_excluded_by_a_gate": len(result.rejected_activators),
        },
        "caveat": RANGE_RESTRICTION,
        "circular": circular,
        "circularity_note": CIRCULARITY_NOTE if circular else "",
        "source": "docs/scientific/empirical-validation.md",
    }


def _headline(result: BenchmarkResult) -> str:
    """Return the one line shown beside a rating."""
    if result.discriminates:
        return (
            f"On {len(result.scored)} assayed sequences this ranking beat chance "
            f"(AUC {result.auc:.2f}, p = {result.auc_p_value:.2f})."
        )
    return (
        f"Measured on {len(result.scored)} assayed sequences, this ranking does not beat "
        f"chance (AUC {result.auc:.2f} against {CHANCE_AUC:.2f} for a coin flip, "
        f"p = {result.auc_p_value:.2f}). Treat it as a prioritisation aid, not a prediction."
    )
