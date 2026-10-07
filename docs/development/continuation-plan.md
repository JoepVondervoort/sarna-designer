# Continuation plan — 7 October 2026

## Objective and authority

Build a human saRNA design workbench that other researchers can use to choose a promoter,
inspect cell-context evidence, design duplexes and export an auditable experimental record.
Scientific traceability takes priority over a persuasive score or interface.

The user's additional requirement is that researchers can adjust the application to their own
experiments. Defaults are explained starting points. Every supported adjustment must show what
it changes, why the default was chosen, where alternatives come from and what remains uncertain.
The configuration editor and its explanations are core functionality.

The user's supplied continuation discussion establishes 16 decisions. Their latest clarification
keeps Wang as the default and selects local installation; those supersede the earlier proposal
to change the default window role. Other papers inform optional recommendations, not universal
rules. These decisions supersede
the corresponding open questions in `CLAUDE.md`; implementation should not ask those questions
again. This document distinguishes implemented foundations from remaining application work.
Read the existing ADRs, glossary and scientific notes before changing their contracts.

## Verified starting point

- The Python engine, Ensembl adapter, rule catalogue, empirical benchmark and static HTML
  renderer exist. There is no interactive application backend or run store.
- Literature-defined anchors already exist in `sarna.genome.anchors`, including three SCN1A
  regions. The distal reference point is explicitly inferred and uncertain; it must not acquire
  a measured-TSS label simply because another annotation source becomes available.
- Cell lines already have a TOML registry and expression, chromatin and TF-occupancy axes.
  This is availability metadata, not an implemented expression-data adapter.
- The AF22 reader already implements enrichment over same-contig background and guarded
  GRCh38-to-hg19 conversion. Preserve those safeguards while generalising it.
- `target-length` already has a maximum of 19. Check the complete supported-length contract
  before describing the chemistry as strictly fixed at 19 nt; an upper bound alone is not that
  contract. Order-ready dTdT duplex export remains absent.
- `tss-window` remains mandatory in Wang's default, as the latest user instruction requires.
  Custom configurations can investigate a preference or disable it, with an explicit departure.
- An existing uncommitted change in `src/sarna/view/template.html` adds visible rule sources,
  reasoning, thresholds, alternatives and configuration-derived window bands. Preserve it,
  review it and exercise it in a browser before counting the transparency work as complete.
- Local verification: 1,279 tests passed, eight live-network tests deselected; strict mypy over
  `src tests scripts` and `ruff check .` passed. This does not establish coverage, live-source
  correctness, browser behavior or availability of the actual histone tracks.
- The subsequent custom-configuration slice implements effective roles and frozen baseline/choice
  metadata, with benchmark and browser checks recorded in
  [configuration verification](custom-configuration-verification.md). Full milestone 1 and the
  interactive application remain unfinished.
- The thesis PDF is present at `docs/191591_vandeVondervoort_BNW_def_DPR.pdf`. This review used
  the repository's extracted dataset and scientific notes; it did not independently re-extract
  the full thesis. `pdftotext` is unavailable in the environment.

README and CLAUDE.md now reflect settled dTdT chemistry, recovered coordinates, assembly handling
and the latest default/delivery decisions. Source metadata and scientific notes still need
reconciliation; historical passages sometimes precede later corrections.

The supplied discussion mentions “seven defects” without enumerating them. Do not invent their
identity or equate them with the 18 legacy defects. Establish an explicit current defect list
from code, tests and rendered output during the first milestone.

## Accepted constraints

1. Wang 2015 remains the default, including its window gate. Other settings and roles are
   optional research choices; literature conformance does not establish activation efficacy.
2. Distance is a sortable column with no points. A custom window preference must explain the
   overlap with the existing preferred-window criterion; do not silently give distance extra points.
3. `pos19-a` remains mandatory in the accepted default; methods must name the known CEBPA
   exception. Researchers may change its role in an advanced custom configuration, as accepted
   in [ADR 0013](../adr/0013-researchers-can-adjust-rule-roles.md).
4. Order chemistry is a 19 nt RNA core plus a 3′ dTdT overhang on both strands.
5. Cell-context assessment precedes candidates, is promoter-level and never changes ranking.
6. The assessment has three outcomes and names its sources. Missing data is not a negative assay.
7. The 5× enrichment threshold is an explicit, provenanced parameter and never blocks design.
8. Human-only support is a real, visible limitation.
9. The application runs locally, with accessible installation and launch paths on supported
   Windows, macOS and Linux systems. The library remains the core; exports are frozen records.
10. Transparency includes a truthful interface, per-run methods and live reranking with a diff.
11. Python is the sole design engine; JavaScript displays results rather than implementing rules.
12. Use the accepted MIT/Apache licensing direction. Audit the actual bundled JASPAR assets:
    official JASPAR documentation states CC BY 4.0, contradicting the repository's non-commercial
    label. Unbundling is therefore not established as a licensing necessity. Record exact asset
    terms and attribution, and the exact code licensing form before replacing LICENSE.
13. Adapter implementation order is archived EPDnew, DepMap, FANTOM5, ENCODE.
14. Saved runs live in one local SQLite file.
15. Order export includes CSV and a sanity report, including off-target findings. Findings report
    rather than exclude; build promoter-level off-target analysis first.
16. A short guided setup leads into a workbench with deliberate anchor and cell-line selection.

## Scientific interpretation to preserve

The AF22 measurements support low enrichment of the assayed active histone marks at the SCN1A
regions in that dataset. They do not directly measure baseline transcription, prove complete
silence or establish why every null result occurred. Keep a chromatin warning separate from an
expression observation and from a proposed biological explanation.

Gene-level expression cannot establish which alternative promoter is active. Display its scope.
HeLa expression remains unverified here; the DepMap adapter must first establish that the exact
line and relevant expression dataset are covered. A missing row is unknown, not zero expression.
Neither additional context nor the accepted rule changes establishes predictive accuracy: the
chapter 4 dataset is range-restricted and ranking has not demonstrated performance above chance.

## Researcher-adjustable settings and their explanations

### Controls to deliver

| Area | Researcher controls | Explanation beside the control |
|---|---|---|
| Starting configuration | Wang default, sourced alternatives and saved custom configurations | Which source supports each part, departures and why they were chosen |
| Sequence criteria | Supported GC bounds, CpG limits, repeat thresholds and related numeric/choice parameters | Definition, units, default origin, alternative evidence and direction of the filtering change |
| Location | Promoter span, anchor/reference, supported TSS window and distance sorting | Coordinate convention, total containment, anchor uncertainty and whether the control filters or scores |
| Rule participation | Set supported scientific rules to mandatory, preference or off in advanced custom configurations | Off means not evaluated; a preference can fail without rejection; a failed gate excludes |
| Chromatin assessment | Supported enrichment threshold, query window, marks and evidence source | Background definition, assay limitations and effects on the warning without effects on ranking |
| Evidence and context | Cell line, available samples/tracks and local evidence | Exact sample identity, assay scope, assembly and missing evidence |
| Results and ordering | Sort columns, candidate selection and saved-run comparison | Viewing/selection changes versus changes requiring a Python design run |

Keep the accepted equal-weight preference count. Adjustable scientific criteria do not imply
adding arbitrary score weights. Keep 19 nt cores with dTdT for the supported order workflow;
different lengths or chemistries require explicit engine/evaluator support, rather than an input
that merely accepts another number. Explain any supported limitation at its control.

Recommended first-release scope is the existing catalogue and local custom configurations;
researcher-authored executable rules and a general plugin system can be considered later.
Describe this as the proposed scope, not a user decision already made.

### What each setting must explain

Each control must expose the current value and role, the selected preset's value and role, and
the source rule's published setting where it exists. These are distinct baselines. Provide:

1. A short plain-language definition, units and the operation it affects.
2. Why this default was selected, with citation and precise table/section/page locator where
   available. Cite project policy as project policy and thesis findings as thesis findings.
3. Named alternatives with their own source and reason for choosing them. A source for the rule
   is not automatically a source for every value a user types.
4. The expected computational consequence of increasing/decreasing the value or changing the
   role, and the biological motivation or uncertainty. Separate these claims.
5. Measured consequences on this promoter after a preview/rerun: candidates admitted or excluded,
   preference denominator changes, rank movement and changed rejection reasons.
6. Known exceptions and limits of validation. A benchmark observation must state its dataset;
   an unmeasured custom setting must not acquire an efficacy claim.

The interface can show a brief explanation beside the setting and expand citations, methods
and caveats. The explanation must be available before applying the change. Unsupported evidence
for a custom value is an explicit answer, not a reason to forbid a computationally valid value.

For example, raising the maximum GC fraction from 0.60 to 0.65 relaxes that gate and allows
candidates above 60% and up to 65% GC to pass it. The additional survivors must be measured
after all other gates: some candidates may still fail elsewhere. The app must not imply that
those extra candidates have demonstrated activity, or attribute the custom value to Wang.

### Editing, comparison and reuse

- Start with an explained preset and offer an advanced editor using the same rule metadata.
  Keep the published Wang settings as the default baseline and clearly name custom departures.
- Label changed values and rules. Provide per-setting reset and whole-preset reset. Show which
  settings are unsupported or inapplicable in the current context and why.
- Keep draft settings separate from the displayed completed run. Apply/preview executes Python,
  then shows the resulting changes. If a calculation is pending or fails, identify which
  configuration the visible results still belong to.
- Show excluded candidates and their reasons; do not add up per-rule rejection counts as if they
  were disjoint. Compare settings with the same sequence, anchor and evidence snapshot where
  possible, and name any additional changes that prevent attributing the difference to one rule.
- Save a custom configuration with a name and version locally; support export/import of the
  configuration for collaborators. An optional research note can explain the user's rationale.
  Changing a saved configuration creates a new version and cannot alter a completed run.
- Record custom departures in run summaries, candidate explanations and methods. Changing the
  preference set changes the denominator; different scores are not directly comparable.

## Configuration foundation required by adjustability

Extend the Python contracts before exposing broad controls. The first implementation slice adds
effective roles, selected option identity, a frozen preset baseline, explicit disabled rules and
departures. It preserves catalogue records and shares validation across library entry points.
Configuration schema/versioning, import, user notes and the interactive editing path remain work
for the persistence and app milestones.

The planned configuration record contains:

- Versioned base preset, base values and base roles.
- Effective selected rules, supported roles and resolved values.
- Explicit choice origin: preset, named sourced option or researcher custom value. Record the
  option id when selected; a numerical match alone does not establish the user's source.
- Disabled rules and changed settings, so omissions are visible in methods as deliberate choices.
- Required explanatory metadata, units, supported ranges and source locators.
- Configuration schema version and optional user note; complete resolved data is sufficient to
  reopen a run without applying a changed preset definition.

Use one Python validation path for the library, import/export and eventual UI. Validate finite
numeric values, integer-only fields and relationships between settings, such as minimum GC not
exceeding maximum GC and window start not exceeding window end. Reject unknown and ineffective
overrides instead of silently ignoring them. Revalidate imported custom configurations.

Distinguish an invalid computation from an exploratory scientific choice. Bounds required for
algorithm correctness should be labelled as supported limits, and source recommendations as
recommendations. Every restriction needs a reason. Do not widen a range until the evaluator,
coordinate handling and methods description support it.

An observed item for milestone 1: `tss-window` offers a thesis option labelled “−2000 to +500 bp”
whose effective `window-end-bp` is 0. Its note admits the restriction, but the option label and
actual value still differ. Display the implemented range accurately and explain the source
range and restriction; allowing downstream targeting is a separate support decision.

## Milestone 1 — reconcile the current implementation

Review and preserve the existing template changes. Generate representative pages for both
strands, inferred and absent TSS references, and different presets/windows. Exercise candidate
selection, sorting, rule explanations and missing-data states in a browser. Inventory concrete
defects and fix only confirmed ones, with appropriate regression checks.

Update README and CLAUDE.md to distinguish current functionality, accepted choices and remaining
work. Record accepted scientific/architecture changes in ADRs, including measured costs and
superseded statements. Keep `legacy/` immutable.

**Exit:** page claims agree with visible behavior; the defect list has named dispositions;
documentation no longer contradicts settled chemistry, anchor support or assembly handling.

## Milestone 2 — implement the accepted rule policy

Support resolved rule roles without mutating global catalogue entries. Preserve Wang's default
rules and parameters. Label a changed run as custom and show its frozen baseline. A changed role
does not establish a better design recommendation. Preserve total-containment
semantics. Distance sorting must state what
position it uses and carry inferred-reference caveats.

Keep `pos19-a` mandatory in the default and pin 19 nt order chemistry. Rerun the chapter 4 benchmark using
coordinates, report exclusions by rule, candidate counts and preference denominators, and record
the costs. “Removing the window's rejections” must not be reported as “every activator passes
all gates”: CEBPA still fails `pos19-a`.

**Exit:** gates, scoring and exported configuration reflect the accepted policy; published and
historical configurations remain reproducible; independent coordinate expectations cover both
strands and window edges.

**Adjustability deliverable:** a validated custom configuration can select supported rules and
values and change scientific rule roles, preserve the chosen source and preset baseline, and
export every departure. Python validates this contract before any UI exposes it. A role change
must reach evaluation, gating, scoring, candidate explanations and exports consistently.

## Milestone 3 — generalise cell-context evidence

Extend the existing cell-line registry with explicit tracks: mark, local path or URL, assembly,
assay/sample identity, source accession and processing details. Separate track identity from a
reader tied to AF22. GSM2902406 is AF22 H3K4me3 and links both GSE108506 and GSE111327;
the accession itself is not shown to be wrong. Audit each track and correct the publication
attribution against its source record. Preserve liftover rejection reasons and local track use.

Make enrichment threshold, background sampling, seed, interval and statistic part of the run
record. Produce a three-valued promoter assessment with source, reason and missing-data details.
Do not reuse TF “bound/not bound” vocabulary for expression or histone measurements. Conflicting
evidence remains visible rather than disappearing into a summary label.

**Exit:** unknown, measured-negative and measured-positive evidence are distinguishable; warnings
appear before candidates, remain non-blocking and cannot affect rank.

## Milestone 4 — add provenance-aware sources

Build a source registry with release/version, retrieval date, assembly, identifiers, licence,
cache location and content checksum. Add EPD first, then DepMap, FANTOM5 and ENCODE. Verify each
source's current format, coverage and terms from its primary documentation during implementation.
Recorded fixtures must include absent, ambiguous and incompatible records; add targeted live
checks for external contracts.

EPD was retired; implement the final November 2023 EPDnew archive rather than depending on a
continuously updated API. Record archive DOI, version and checksums. EPDnew incorporates FANTOM5
and ENCODE data, so overlapping records are not independent corroboration.

EPD offers measured anchor alternatives without silently replacing an inferred reference.
DepMap adds gene expression in explicitly matched cell models. FANTOM5 adds supported TSS evidence;
ENCODE adds assay-specific cell-context evidence. Support researcher-supplied local evidence for
lines absent from public datasets.

Resolve JASPAR's incorrect licence metadata before distribution; retain attribution and choose
bundling versus cached retrieval based on verified terms and installation practicality.
Explain external retrieval in setup, including what gene/coordinate queries are sent. Keep
researchers' uploaded sequences and results local. “Local app” does not imply optional evidence
retrieval is network-free. Provide an explicit cached/offline mode.

Audit the September 2026 curated saRNA database and its linked primary experiments as optional
external validation material. Check sequence identity/orientation, endpoints, cell models, negative
results, duplicates and overlapping studies before drawing conclusions. Hold out genes/studies
where feasible. Database recommendations do not replace the Wang default or become truth by
being newer; rank performance and empirical exceptions remain claims requiring measurement.

**Exit:** every observation can be traced to an exact source record, and cached/offline runs do
not silently refresh their evidence.

## Milestone 5 — freeze runs before building the workbench

Define a versioned run record containing input sequence, assembly and coordinates, anchor and
uncertainty, cell context, resolved rules and kinds, software version, evidence snapshots and
source hashes, candidates, verdicts and ranks. Store snapshots in SQLite with schema migrations;
cached URLs alone are insufficient archival evidence.

Changed settings produce a derived run linked to its parent. Compare candidates by stable
sequence/location/orientation identity rather than rank or a display ID alone. Report admitted
and excluded candidates, rank movement, preference changes and changed reasons.

**Exit:** a saved run reopens offline without consulting today's preset definitions, and a diff
explains a deliberate parameter change.

**Adjustability deliverable:** saved custom configurations are versioned independently of run
snapshots. Import/export preserves choice origins and supported limitations; reopening a frozen
record requires no access to the current preset catalogue for interpretation. Recomputing it
also requires the recorded engine version or an explicit report of version differences.

## Milestone 6 — deliver the local workbench

Deliver a browser workbench backed by a thin local Python service. Bind to the loopback interface;
launch it through an application launcher. Keep science in Python and use HTML/CSS/JavaScript for
controls, display and progress. No hosted service or account is required. Store runs and caches
in a documented per-user data directory, separate from the installation, so updates preserve them.

Guided setup selects gene or pasted sequence, deliberate anchor, cell line and configuration.
The workbench shows promoter context, evidence assessment, passing and rejected candidates,
rule explanations, settings, reranking and saved-run comparisons. A pasted sequence can proceed
with explicit unknown genomic context. Exports remain frozen even when the workbench can rerun.

Use the researcher-settings contract above as the editor's acceptance criteria. Complete one
working path from changing a GC bound, through Python validation and rerun, to an explained diff
and saved configuration before adding all controls. Then add rule participation/roles within
the agreed scope, location controls, chromatin controls and the remaining catalogue fields.

**Exit:** a researcher can complete setup, inspect warnings, change a setting, compare and reopen
a run through the application without writing Python.

Also verify a sourced alternative, a custom value, an invalid combination, a reset, an imported
configuration and missing-context behavior. The researcher can distinguish the setting's
motivation, its measured effect on this run and the absence of validated efficacy evidence.

## Milestone 7 — produce auditable bench exports

Build promoter-level off-target reporting first, with explicit reference/search scope and
matching assumptions. Do not call it comprehensive genome-wide or transcriptome-wide safety.
Missing search results are labelled unavailable; off-target findings never exclude candidates.

Generate order CSV with both strands written 5′→3′, unambiguous RNA cores and `[dT][dT]`
overhang notation, core/total lengths, candidate identity and run identity. The sanity report
checks reverse complementarity of the cores, chemistry, lengths, duplicate orders, orientation
and available off-target evidence. Do not silently interpret mixed DNA/RNA ordering notation.

Generate per-run methods from the saved record: anchor derivation, cell-context evidence,
citations, resolved thresholds, gates/preferences, ranking/tie-breaks, CEBPA exception, empirical
validation limitations and source versions. Bundle frozen HTML, methods, CSV and sanity report.
Methods and chemistry serialization can be developed against milestone 5 before the UI is done.

**Exit:** an exported selection agrees with the saved run and can be reviewed independently of
the application. Claims about off-target coverage and activity remain bounded by actual evidence.

## Immediate next implementation slice

Complete milestone 1, then deliver milestone 2 as one reviewable scientific change with ADRs,
benchmark output and generated example pages. This resolves the accepted policy before building
data adapters or a user interface around stale semantics. No further confirmation of the accepted
decisions is needed. Milestone 1 browser checks and documentation cleanup remain unfinished.

Deliver milestone 2 in small steps:

1. Record the default-versus-custom policy and published baseline in ADRs.
2. Extend configuration resolution, choice origins and shared validation; exercise a sourced
   option and a hand-edited GC bound through the Python engine.
3. Preserve the Wang default and allow explicitly labelled custom roles and sourced alternatives.
4. Export the resolved values, effective roles, disabled rules and deviations; measure candidate
   changes and benchmark exclusions, then render example pages that explain them.

Run persistence and the first settings-to-diff workbench path depend on these contracts and the
source/evidence snapshot format. They do not need to wait for every public adapter to be complete.
Keep the EPD → DepMap → FANTOM5 → ENCODE adapter order while delivering the application
incrementally with explicit unavailable evidence.

The user has accepted advanced changes of supported scientific rules between mandatory,
preference and off. [ADR 0013](../adr/0013-researchers-can-adjust-rule-roles.md) records this choice,
its costs and the explanation/export requirements. Named preset defaults remain versioned
baselines. Engine and frozen-page support are implemented; the interactive editor, configuration
versioning and import/export remain pending.

## Local installation and release requirements

Start a packaging feasibility check before substantial app work: native genomics dependencies,
including bigWig access, can constrain supported operating systems and architectures. Test what
can actually be bundled; do not promise universal compatibility from Python alone.

The research-user path should install and launch the workbench without editing Python or requiring
Docker. Provide platform packages where feasible, a documented Python installation for developers
and an optional container route for labs. If a platform initially needs an alternative path,
state the support boundary and concrete steps rather than claiming a native installer exists.

Each released platform needs a clean-machine installation and upgrade check: launch, an offline
example design, an invalid custom setting, a saved run, export, restart and reopen. Verify that an
upgrade preserves the local database and provide backup/restore instructions. Include a concise
quick start, optional-source downloads and sizes, supported versions/architectures, troubleshooting,
and an example dataset. Large optional resources must not be mandatory just to start the tool.

Release packaging follows the first complete settings-to-saved-run app path. The feasibility
check starts earlier so distribution constraints can shape dependency choices in time.

## Primary sources for the corrected source strategy

- [EPD retirement notice and final archive](https://www.expasy.org/archives/epd)
- [EPDnew human data sources](https://epd.expasy.org/epd/human/human_database.php)
- [JASPAR documentation and licence](https://jaspar.elixir.no/about/)
- [AF22 H3K4me3 sample GSM2902406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2902406)
- [2026 curated saRNA database paper](https://doi.org/10.3390/molecules31183249)
- [Authors' database and source records](https://github.com/NIC-SBI/sarna-database)
