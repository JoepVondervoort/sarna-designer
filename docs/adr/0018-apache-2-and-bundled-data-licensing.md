# ADR 0018 — Apache-2.0 for the code, and the bundled data's actual terms

Date: 2026-10-07. Status: accepted. Supersedes the placeholder `LICENSE` and the
"licensing undecided" entries in `README.md` and `CLAUDE.md`.

## Context

`LICENSE` was a deliberate placeholder reading *"all rights reserved... no permission to use,
copy, modify or distribute is granted"*, and `pyproject.toml` pointed at it. The package was
therefore undistributable: handing it to a colleague would have granted them no rights at all.
This blocked the whole delivery direction accepted in [ADR 0014](0014-local-delivery-and-wang-default.md).

Two questions had to be answered together: which licence the code carries, and what the bundled
data actually permit. The second was the riskier one, because the repository's own metadata was
wrong about it.

## Decision: Apache-2.0 for the software

MIT and Apache-2.0 were the candidates; both are OSI-approved and permissive, and the choice is
not between open and closed. Apache-2.0 was chosen for its **express patent grant** (section 3).
This tool designs molecules in a field with active therapeutic patents — MTL-CEBPA, built on a
duplex this engine's default configuration excludes, is in human trials. An express grant, and
the reciprocal termination clause for anyone who sues over it, protects both the author and
downstream users at no practical cost. MIT would have been defensible; it simply offers nothing
on that axis.

The accepted direction in `CLAUDE.md` was always "MIT/Apache"; this fixes the form.

Consequences: Apache-2.0 obliges anyone redistributing a modified version to state their changes
and to carry the `NOTICE` file (section 4). `NOTICE` is therefore part of the distribution, listed
in `license-files` and verified present in the built wheel.

## Decision: the bundled JASPAR subset is CC BY 4.0, not CC BY-NC-SA 4.0

`src/sarna/regulatory/data/jaspar_core_selection.json` recorded its own licence as
**CC BY-NC-SA 4.0**. [JASPAR's about page](https://jaspar.elixir.no/about/), checked on
2026-10-07, states that JASPAR "is licensed under the Creative Commons Attribution 4.0
International License".

The repository's label was wrong **in the restrictive direction**, which is the more insidious
error of the two: it would have led users to refuse uses that are in fact permitted, and it
suggested a share-alike obligation that does not exist. `CLAUDE.md` had already flagged the
contradiction and, correctly, warned against assuming unbundling was legally required. It is not.

CC BY 4.0 permits redistribution and commercial use with attribution and imposes no share-alike
term, so the subset ships as it is, with attribution in `NOTICE` and the corrected licence string
in the data file itself. The data file also records that it previously said otherwise, so anyone
who read the old label can see it was retracted rather than quietly edited.

## Decision: the chapter 4 dataset ships

The 125 assayed sequences and the recovered coordinates for 123 of them come from a publicly
defended thesis, so redistribution raises no embargo question. They ship with the thesis cited in
`NOTICE` and in `Provenance` inside `dataset.py`.

This matters beyond permission: it is the only data in the package that can show a design rule to
be **wrong** rather than merely unconventional, and removing it would leave the validation harness
unable to run. Its range restriction — every sequence was built to pass Wang's gate, so it can
reveal a gate that is too strict but never one too permissive — travels with it in
`CHAPTER4_CAVEATS`.

## Also removed: the unused `api` extra

`pyproject.toml` declared `api = ["fastapi", "uvicorn[standard]"]`. Nothing in `src`, `tests` or
`scripts` imports either: [ADR 0015](0015-local-workbench-and-frozen-ui.md) chose a stdlib HTTP
adapter precisely so no web framework would be required. A published package advertising an extra
that installs two unused dependencies invites a reader to conclude there is a web API. There is
not, and the extra is gone.

## Costs and limits

Apache-2.0 is longer than MIT and its `NOTICE` obligation is a small ongoing burden: new bundled
data must be added there, not only in a docstring. The patent grant is worth little if no patent
is ever asserted, which is the likely case.

This ADR records the terms of what is bundled **today**. Adapters planned for archived EPDnew,
DepMap, FANTOM5 and ENCODE each carry their own terms, and none has been audited here. Any
dataset that an optional adapter reads from local disk — the AF22 histone tracks today — is not
redistributed and keeps whatever terms it was obtained under; that distinction must survive as
more local sources are added.

A licence is not legal advice, and institutional policy on IP in doctoral work has not been
checked here. If the institution claims rights in this code, that supersedes this ADR.
