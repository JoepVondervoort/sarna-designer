# ADR 0013 — Researchers can adjust scientific rules and their roles

- **Status**: Accepted; engine and frozen-page support implemented, interactive editor pending
- **Date**: 2026-10-07
- **Deciders**: Joep Vondervoort

## Context

The user requires an adjustable research application whose defaults, alternative settings and
consequences are explained. They explicitly accepted allowing advanced users to change supported
scientific rules between mandatory, preference and off, including rules mandatory in a preset.

[ADR 0003](0003-rule-provenance-is-first-class.md) requires explanations and provenance;
[ADR 0010](0010-presets-are-selections-with-overrides.md) requires immutable catalogue records
and resolved configurations. This decision extends custom configurations to include rule roles.
It retains the equal-weight preference count from
[ADR 0006](0006-criteria-count-ranking.md).

## Decision

Researchers may adjust supported parameter values and, in advanced custom configurations, set
supported scientific rules to:

- **Mandatory:** failure excludes the candidate; this rule contributes no preference points.
- **Preference:** passing earns one point; failure contributes no point and does not exclude.
- **Off:** the rule is not evaluated, contributes no point and cannot exclude a candidate.

Off remains distinct from an evaluated rule reporting `not-applicable`. The configuration and
methods record deliberately disabled rules even though candidates have no verdict for them.
Existing continuous tie-breaks retain their separately defined behavior; role editing does not
automatically turn a continuous value into a binary criterion.

Named presets are versioned starting points. User edits produce a custom configuration without
mutating either the preset or the rule catalogue. `pos19-a` remains mandatory in the accepted
default but may become a preference or be disabled in a custom configuration. This clarifies
the earlier default decision rather than reversing it.

Every editable setting must show the selected preset's baseline, its citation and rationale,
sourced alternatives, supported limits and consequences of editing. The UI distinguishes the
source of a rule from the origin of the user's chosen value or role. A custom change does not
inherit a claim of published support.

Applying a change runs the Python engine and explains changes to survivors, preference counts,
ranks and rejection reasons. Saved runs and methods include resolved values and roles, disabled
rules, choice origins and departures from the base preset. A reset can restore the baseline.

Algorithm validity and supported capabilities remain enforced: role editing cannot invalidate
sequence/coordinate types or remove required generation inputs. The supported order workflow
remains 19 nt RNA cores plus dTdT. Different lengths or chemistries require a separately supported
implementation and scientific description.

## Costs and alternatives

Keeping mandatory rules fixed would simplify comparison and preserve each preset's admission
policy, but would prevent researchers from investigating known exceptions. Threshold-only
editing cannot test whether a parameter-free rule such as `pos19-a` should be a preference.

The accepted flexibility allows configurations with less published support and increases the
surface for interpretation errors. The app must label departures before applying them, retain
rejected candidates and explain actual computational changes without claiming improved efficacy.
Changing the preference set changes the denominator, so raw scores across configurations are
not directly comparable.

The repository's [empirical validation](../scientific/empirical-validation.md) reports that
`pos19-a` rejects one of 14 measured activators, CEBPA. Relaxing it removes that rule's exclusion;
it does not establish that the candidate passes every other gate or that similar candidates
will activate. The candidate count and ranking cost of each supported role change must be
measured during implementation, using the same sequence and evidence snapshot for comparisons.

## Implementation requirements

Extend resolution with effective roles rather than modifying global rule kinds. Audit the full
path through evaluators, result construction, gating, scoring, views and serialization: current
code often derives behavior from the catalogue rule's kind. Verify both a gate-to-preference
change and disabling a rule with independently expected results.

Keep one Python validation contract for the library, UI and imported custom configurations.
Record why a role transition is unsupported where a rule cannot satisfy its semantics. Exports
must retain enough metadata to interpret custom choices without today's preset catalogue.
