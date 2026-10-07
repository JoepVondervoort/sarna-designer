"""Gene lookup and run assembly over the existing engine, with explicit evidence gaps."""

import logging
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from sarna import __version__, settings
from sarna.app.models import DesignRequest
from sarna.app.references import (
    anchor_record,
    custom_anchor,
    custom_anchor_record,
    gene_overview,
)
from sarna.app.store import RunStore
from sarna.design.candidate import Candidate
from sarna.design.order import CSV_COLUMNS, csv_row, order_for
from sarna.design.pipeline import design_sarnas
from sarna.genome.anchors import PromoterAnchor, load_anchors
from sarna.genome.ensembl import (
    EnsemblCacheMissError,
    EnsemblClient,
    EnsemblError,
    EnsemblUnavailableError,
    Gene,
    Transcript,
)
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.regulatory.chromatin import open_af22_reader
from sarna.regulatory.evidence import CELL_TYPES
from sarna.regulatory.features import fetch_regulatory_features
from sarna.regulatory.motifs import scan_all
from sarna.rules.catalog import CATALOGUE
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleOutcome
from sarna.sequence.nucleic import NucleicSequence
from sarna.validation.ranking_performance import performance_summary
from sarna.view.promoter_view import PromoterView, alternative_tss_spans, cell_panels_for

LOGGER = logging.getLogger(__name__)


def candidate_record(candidate: Candidate, model: str, end_length: int) -> dict[str, Any]:
    """Publish an itemised rule rating and one named model's guide, never a probability."""
    record = candidate.to_dict()
    record["guide_prediction"]["end_length_nt"] = end_length
    primary = next((a for a in candidate.guide_comparison.assessments if a.table_id == model), None)
    evaluated = candidate.evaluated_preference_points
    record["rating"] = {
        "met": candidate.preference_points, "evaluated": evaluated,
        "total": candidate.max_preference_points,
        "percent": round(100 * candidate.preference_points / evaluated) if evaluated else None,
        "label": f"{candidate.preference_points}/{evaluated} evaluated preferences",
        "unavailable": candidate.max_preference_points - evaluated,
        "unchecked_gates": [
            r.rule_id for r in candidate.results
            if str(r.kind) == "mandatory" and r.outcome is RuleOutcome.NOT_APPLICABLE
        ],
        "guide": str(primary.predicted_guide) if primary else None,
        "guide_model": model, "models_agree": candidate.guide_comparison.tables_agree,
    }
    return record


class DesignService:
    """Coordinate a local job while preserving the source inputs and engine results."""

    def __init__(
        self, data_dir: Path, client_factory: Callable[[bool], EnsemblClient] | None = None,
    ) -> None:
        """Use an injectable client for recorded-source tests and cached/offline requests."""
        self.store = RunStore(data_dir / "runs.sqlite3")
        self.client_factory = client_factory or (
            lambda offline: EnsemblClient(cache_dir=settings.cache_dir("ensembl"),
                                          offline=None if settings.is_offline() else offline,
                                          timeout=15, backoff_seconds=0.5)
        )

    @staticmethod
    def _gene(client: EnsemblClient, symbol: str) -> Gene:
        # Human symbols are normalised before both cache access and external lookup.
        return client.gene(symbol.strip().upper())

    def lookup(self, symbol: str, offline: bool = False) -> dict[str, Any]:
        """Resolve a human symbol or gene ID and expose every deliberate reference choice."""
        client = self.client_factory(offline)
        gene = self._gene(client, symbol)
        transcripts = client.transcripts(gene)
        return gene_overview(gene, transcripts)

    def run(self, run_id: str, request: DesignRequest) -> None:
        """Execute a submitted request, publishing progress and preserving errors."""
        try:
            self._run(run_id, request)
        except EnsemblCacheMissError as error:
            self.store.fail(run_id, str(error), "offline-cache-miss")
        except EnsemblUnavailableError as error:
            self.store.fail(run_id, str(error), "ensembl-unavailable")
        except Exception as error:
            LOGGER.exception("Design run %s failed", run_id)
            self.store.fail(run_id, str(error) or type(error).__name__)

    def _run(self, run_id: str, request: DesignRequest) -> None:
        config = request.resolve()
        client = self.client_factory(request.offline)
        self.store.progress(run_id, "Retrieving the promoter sequence")
        gene, transcript = None, None
        transcripts: tuple[Transcript, ...] = ()
        selected_anchor: dict[str, Any] | None = None
        if request.source == "gene":
            gene = self._gene(client, request.symbol)
            transcripts = client.transcripts(gene)
            if request.custom_reference is not None:
                anchor = custom_anchor(gene, request.custom_reference)
                selected_anchor = custom_anchor_record(anchor, request.custom_reference)
            elif request.anchor_id:
                anchor = load_anchors()[request.anchor_id]
                if anchor.gene.upper() != gene.symbol.upper():
                    raise ValueError("The chosen literature anchor belongs to another gene.")
            else:
                choices = [t for t in transcripts if (
                    t.ensembl_id == request.transcript_id if request.transcript_id
                    else t.is_mane_select
                )]
                if len(choices) != 1:
                    raise ValueError("Choose a transcript: no unique MANE Select or matching id.")
                transcript = choices[0]
                anchor = PromoterAnchor.from_transcript(
                    gene=gene.symbol, transcript_id=transcript.ensembl_id,
                    tss=transcript.tss_coordinate, interval=transcript.interval,
                )
            promoter = client.promoter_for_anchor(anchor, request.upstream, request.downstream)
            if selected_anchor is None:
                selected_anchor = anchor_record(anchor)
        else:
            lines = request.sequence.strip().splitlines()
            if lines[0].startswith(">"):
                if sum(line.startswith(">") for line in lines) != 1:
                    raise ValueError("Supply one promoter sequence, not multiple FASTA records.")
                lines = lines[1:]
            sequence = NucleicSequence.parse("".join(lines))
            if not 19 <= len(sequence) <= 20_000:
                raise ValueError("The supported pasted sequence length is 19 to 20,000 bp.")
            promoter = Promoter(sequence=sequence, tss_start=(
                TssPosition(request.tss_start) if request.tss_start is not None else None
            ))

        notices = [f"No gene-expression measurement is loaded for {request.cell_line}. "
                   "Cell choice does not establish promoter activity or change ranks."]
        self.store.progress(run_id, "Evaluating candidates with your settings")
        result = design_sarnas(promoter, config)
        regulatory = None
        if request.regulatory and promoter.interval is not None:
            self.store.progress(run_id, "Retrieving regulatory annotations")
            try:
                regulatory = fetch_regulatory_features(client, promoter.interval)
            except EnsemblError as error:
                notices.append(f"Regulatory annotations unavailable: {error}")
        else:
            notices.append("Regulatory annotations were not requested, "
                           "or genomic context is absent.")
        self.store.progress(run_id, "Annotating sequence motifs and cell context")
        scans = (scan_all(promoter.sequence, shuffles=request.motif_shuffles)
                 if request.motifs else ())
        if not request.motifs:
            notices.append("Motif scanning was switched off.")
        reader = None
        known_cell = next((c for c in CELL_TYPES.values()
                           if request.cell_line.lower() in (c.id, c.name.lower())), None)
        if request.chromatin and known_cell is not None and known_cell.id == "af22":
            try:
                reader = open_af22_reader()
            except (OSError, ImportError) as error:
                notices.append(f"AF22 tracks could not be loaded: {error}")
        if reader is None:
            notices.append("No histone measurements were loaded. "
                           "This does not establish closed chromatin.")
        view = PromoterView(
            result=result, gene=gene, transcript=transcript,
            alternative_tss=alternative_tss_spans(promoter, transcripts, transcript),
            regulatory=regulatory, motif_scans=scans,
            cell_panels=cell_panels_for(promoter=promoter, chromatin_reader=reader),
        )
        payload = view.to_dict(full_rules_for_rejected=False)
        # Keep complete verdicts in SQLite and lazy-load individual details in the UI.
        thermo_rule = next((config.rule(key) for key in (
            "thermodynamic-asymmetry", "asymmetry-margin"
        ) if key in config), None)
        model = thermo_rule.text("nearest-neighbour-table") if thermo_rule else "rna-xia-1998"
        end_length = thermo_rule.integer("end-length-nt") if thermo_rule else 4
        payload["candidates"] = [
            candidate_record(c, model, end_length) for c in view.ranked_candidates()
        ]
        # What this ranking was measured to achieve, computed from the bundled chapter 4 data
        # rather than transcribed. It belongs in the record so an archived run carries the same
        # caveat the live page shows. `applies` is false once a researcher adjusts the rules:
        # the measurement describes the preset, not whatever they changed it into.
        measured = performance_summary(request.preset_id)
        payload["ranking_performance"] = (
            None
            if measured is None
            else {
                **measured,
                "applies": not config.is_custom,
                "custom_note": (
                    "This run departs from the preset baseline, so the measurement below "
                    "describes the unmodified preset and not the configuration you ran. No "
                    "measurement exists for a custom configuration."
                    if config.is_custom
                    else ""
                ),
            }
        )
        # Order-ready duplexes are computed here, in Python, and carried in the record, so a
        # frozen archive can hand over the exact molecules without JavaScript re-deriving the
        # chemistry (accepted constraint 11). Passing candidates only, in rank order: a duplex
        # that failed a gate is not something to order.
        payload["order"] = {
            "columns": list(CSV_COLUMNS),
            "rows": [
                csv_row(order_for(c, model), run_id)
                for c in sorted(
                    (c for c in view.ranked_candidates() if c.rank is not None),
                    key=lambda c: c.rank or 0,
                )
            ],
        }
        if payload["promoter"]["tss_offset"] is None and promoter.tss_start is not None:
            # Relative-only pasted DNA has an offset but no genomic coordinate to invert.
            offset = promoter.tss_start.distance_to(TssPosition(1))
            if 0 <= offset < promoter.length:
                payload["promoter"]["tss_offset"] = offset
        payload["app"] = {
            "schema_version": 1, "engine_version": __version__, "run_id": run_id,
            "created": datetime.now(UTC).isoformat(), "request": request.model_dump(mode="json"),
            "anchor": selected_anchor, "cell_line": request.cell_line,
            "gene_overview": gene_overview(gene, transcripts) if gene is not None else None,
            "data_access": {"offline": client.offline},
            "notices": notices, "rating_note": "Criteria met, not a probability of activation.",
            "available_rules": [ResolvedRule(rule).to_dict() for rule in CATALOGUE],
        }
        self.store.progress(run_id, "Saving the complete local record")
        self.store.complete(run_id, payload)
