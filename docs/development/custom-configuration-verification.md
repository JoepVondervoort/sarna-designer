# Custom configuration verification: 7 October 2026

## Implemented scope

The engine supports run-specific mandatory/preference roles and explicit off selections,
without changing the rule catalogue or Wang default. Result factories use the resolved role so
gating, preference counts, ranking, candidate verdicts and exports agree. Continuous tie-breaks
cannot be converted into binary criteria; candidate-length generation cannot be disabled or scored.

`Preset.resolve` accepts `role_overrides` and `option_ids`, alongside existing numeric overrides.
A named option overrides conflicting preset values before manual overrides are applied. Exports
retain the effective roles, catalogue roles, frozen preset rules, disabled records, deliberate
option identity and custom departures. Bare matching values do not acquire option provenance.

Validation rejects unknown/ineffective preset changes, nonfinite values, fractional integer
parameters and inverted supported ranges. An overlapping pair of window preferences is allowed
but explained with a warning; users can disable one to retain one location criterion.

Frozen pages label custom configurations, show disabled rules and role changes, and distinguish
selected options from numerical matches. UTF-8 metadata fixes garbled scientific symbols when
opening a page locally. Existing uncommitted template improvements were preserved.

The interactive settings editor, versioned configuration import/export, local app service,
SQLite run persistence, methods/order exports and packaged installers remain unimplemented.
Existing parameter-origin labels distinguish catalogue defaults, named options and overrides;
full preset-versus-user editing history is still required for the later configuration schema.

## Independent expectations

Three overlapping 19 nt windows ending in A, T, A give a hand-computed test:

| Position-19 role | Passing offsets | Ranked offsets | Preference points in rank order |
|---|---|---|---|
| Mandatory | 0, 2 | 0, 2 | 0, 0 |
| Preference | 0, 1, 2 | 0, 2, 1 | 1, 1, 0 |
| Off | 0, 1, 2 | 0, 1, 2 | 0, 0, 0 |

Other regression checks cover preference-to-gate promotion, missing TSS evidence, immutable
defaults, explicit sourced choices, named-option/manual precedence, disabled-rule provenance,
required input restrictions, continuous tie-break restrictions and invalid custom values.

Verification: 1,300 tests passed with eight live-network tests deselected; strict mypy checked
101 source files successfully, Ruff passed and `git diff --check` found no whitespace errors.
The 21 new configuration cases are included in that total. Coverage was not remeasured.

## Chapter 4 measurement

Used `load_chapter4_dataset()`, `load_chapter4_locations()` and the existing `benchmark` harness:
125 sequences, 14 known activators; recovered locations for 123 sequences, 519 bp context.
Each configuration used 10,000 permutations and seed 20261006.

| Configuration | All rejected | Known activators rejected | Preference count | Score AUC |
|---|---:|---:|---:|---:|
| Wang default | 98 | 10 | 7 | 0.4575 |
| Window becomes a preference, preferred-window retained | 25 | 1 | 8 | 0.4550 |
| Window becomes a preference, preferred-window off | 25 | 1 | 7 | 0.4437 |
| Window off | 25 | 1 | 7 | 0.4575 |
| Window and position-19 off | 14 | 0 | 7 | 0.4575 |

CEBPA is the remaining rejected activator when only the window gate is relaxed. These measurements
explain admission changes, not better efficacy. The dataset is range-restricted, context cannot
support the default 500 bp island definition, and missing locations remain unknown. AUC here is
the existing harness's comparison of preference counts across all 125 assayed sequences, including
rejected sequences; it is not an evaluation restricted to the surviving shortlist. Scores with
different preference sets are not directly comparable. Retain these limitations in methods.

## Browser verification

Generated default plus-strand, custom minus-strand and custom unanchored example pages. A real
headless Chromium DOM check of the custom minus-strand page verified custom/baseline labels,
disabled records, the overlap warning, UTF-8 symbols, passing/all candidate views, changing the
selected candidate and visible explanations. External hostname resolution was disabled.
JavaScript syntax checking also passed.

This is a DOM/interaction smoke check, not visual or cross-platform certification. Chromium emitted
GPU initialization diagnostics in the headless environment; no screenshot-based verification was
performed. Broader browser/strand/context coverage and installation testing remain in the plan.
