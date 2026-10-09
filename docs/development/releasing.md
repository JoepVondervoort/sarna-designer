# Publishing a new version, and how colleagues get it

Two questions, one answer each. "How do I ship an update" is the maintainer's, "how do I get the
update" is the colleague's, and they are coupled: the second command only works if the first one
put the work where it looks.

## The short version

| | Maintainer ships it by | Colleague gets it with |
|---|---|---|
| **Route B, from GitHub** | merging to the **default branch** (`main`) and pushing | `uv tool upgrade sarna-designer` |
| **Route A, a `.whl` file** | building a wheel and sending it | rerunning the install command on the new file |

Route B is the better setup, and it is the one to move to. It needs no attachment, no second
copy of anything, and the colleague's update is one command they can remember. Route A exists
because the repository is private and not everyone has access to it.

## What was actually measured

These behaviours were tested against `uv 0.12.23` on 8 October 2026, by building two releases
that **shared the same version string** and differed only in their code, installing one and then
trying to reach the other. Guessing at them is how the wrong advice gets written down.

| Situation | Command | Result |
|---|---|---|
| Installed from git, new commit pushed | `uv tool upgrade sarna-designer` | **Updates.** uv replaced `@66ab565` with `@e433276` |
| Installed from git, new commit pushed | `uv tool install "...@ git+https://..."` | Updates |
| Installed from a wheel, new wheel built | `uv tool install <new path>` | **Updates.** uv replaced the old file URL with the new one |
| Installed from a wheel | `uv tool upgrade sarna-designer` | **Does nothing**, and says `Nothing to upgrade` |

The finding that matters: **`uv tool upgrade` tracks the git commit, not the version number.**
An unchanged version does not stop a git-installed tool from updating. The version number is
needed for a different reason, below.

The second finding: a wheel install has no upstream, so `uv tool upgrade` is a no-op for it. The
message is correct rather than a failure, and the instructions now say so, because a colleague
who reads `Nothing to upgrade` otherwise concludes they are up to date when they are not.

`--reinstall` is not needed in either case. It rebuilds the environment from scratch, which is
worth reaching for only when an install looks damaged.

## Why the default branch decides everything

`uv tool install "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"`
clones the repository and installs **the default branch**, which is `main`. `uv tool upgrade`
later reinstalls from the same place. Neither command knows about any other branch.

So work that sits on a feature branch is invisible to every colleague, and `uv tool upgrade`
reports success while installing the same stale code it installed last time. That is the single
most important thing to get right about this setup, and it is a property of the branch, not of
the tooling.

## Releasing

1. **Decide the version** and set it in `pyproject.toml`. Bump the middle number for new
   capability, the last one for fixes. `tests/unit/test_release_metadata.py` fails if the wheel
   filename quoted in `README.md` or `INSTALL.md` no longer matches, so update those together.
2. **Run the gates.**

   ```bash
   .venv/bin/python -m pytest -q
   .venv/bin/python -m mypy src tests scripts
   .venv/bin/python -m ruff check .
   .venv/bin/python scripts/verify_workbench_browser.py --chrome <chromium> --gene SCN1A
   ```

3. **Merge to `main` and push.** Until this happens, nothing a colleague runs can reach the new
   work.
4. **Tag it**, so a version number points at an identifiable commit:

   ```bash
   git tag -a v0.2.0 -m "saRNA Designer 0.2.0"
   git push origin v0.2.0
   ```

5. **Build the wheel** for the colleagues who install by file:

   ```bash
   python -m build --wheel          # dist/sarna_designer-<version>-py3-none-any.whl
   ```

6. **Tell people which route they are on.** Route B colleagues need one command. Route A
   colleagues need the file.

## Why the version number still matters

It is not needed to make `uv tool upgrade` work. It is needed because the version is part of the
scientific record:

- `sarna-designer-check` prints it, which is the only way a colleague can tell you what they are
  running when something looks wrong.
- Every saved run stores it as `engine_version`, and every export carries it.

While the version was frozen at `0.1.0.dev0`, every run ever produced recorded the same string,
so no result could be traced to the code that produced it. Two runs that disagree would have
looked like they came from identical software. That is why releases now carry real numbers and
tags, and why the placeholder is gone.

## Worth doing next

- **GitHub Releases.** Attaching the `.whl` to a tagged release gives route A colleagues a
  stable link instead of an email attachment, and gives every build a dated, immutable record.
  It needs the colleague to have repository access, so it does not replace sending a file to
  someone who has none.
- **A release checklist run in CI.** The gates above are run by hand today. Nothing enforces
  that they passed before a tag was pushed.
