# Installing saRNA Designer

**Written for:** a colleague who wants to use the tool, with no assumption of programming
experience. If you are comfortable with Python packaging, [the short version](#the-short-version-for-people-who-know-python)
is at the bottom.

saRNA Designer runs **on your own computer**. It is not a website, there is no account, and
nothing is uploaded except the gene lookups described under [What leaves your computer](#what-leaves-your-computer).

You need about **ten minutes** and an internet connection for the setup.

---

## Before you start

You do **not** need to install Python yourself. The installer below fetches the exact version
the tool needs and keeps it separate from anything else on your machine, so it cannot break
other software you use.

You do **not** need WSL, Docker, or administrator rights on Windows.

---

## Windows

Most people will use this. Everything happens in a program called **PowerShell**, which is
already on your computer.

### Step 1 — Open PowerShell

Press the **Windows key**, type `powershell`, and press **Enter**. A blue window opens with a
blinking cursor. You type commands here and press Enter after each one.

> If a command seems to do nothing for a while, that is normal — it is downloading. Wait until
> the blinking cursor comes back before typing the next one.

### Step 2 — Install `uv`

Copy this whole line, paste it into PowerShell (right-click pastes), and press **Enter**:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

This installs a small program called `uv` that manages the installation for you.

**Now close PowerShell completely and open it again.** This is necessary — the new program is
only visible to a freshly opened window.

### Step 3 — Install saRNA Designer

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

This downloads the tool and the correct Python. It prints a list of packages and ends with
`Installed 2 executables: sarna-designer, sarna-designer-check`.

> **Given the tool as a folder or zip file instead?** Unzip it, then use the folder path:
> ```powershell
> uv tool install --python 3.11 "C:\path\to\sarna-designer[app]"
> ```

### Step 4 — Check that it works

```powershell
sarna-designer-check
```

You should see a list of checks and then **`RESULT: PASS`**. If any line says `FAIL`, the
check tells you what to do — see [If something goes wrong](#if-something-goes-wrong).

### Step 5 — Start it

```powershell
sarna-designer
```

It prints a web address like `http://localhost:8001/`. Open that in your browser. **Leave the
PowerShell window open** while you work — closing it stops the tool.

To stop it, click the PowerShell window and press **Ctrl+C**.

---

## macOS

### Step 1 — Open Terminal

Press **Cmd+Space**, type `terminal`, press **Enter**.

### Step 2 — Install `uv`

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then close Terminal and open it again.

### Step 3 — Install saRNA Designer

```bash
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

### Steps 4 and 5 — Check and start

```bash
sarna-designer-check
sarna-designer
```

Open the printed address in your browser, and leave Terminal open while you work.

---

## Linux

Identical to macOS. If `uv` is not found after installing, run `source $HOME/.local/bin/env`
or open a new terminal.

---

## Using it

Open the address the tool printed. Then:

1. Type a human gene symbol — `SCN1A`, `TP53`, `KCNQ2` — and click **Look up**.
2. Choose which promoter reference to design against. The tool picks the MANE Select transcript
   when there is exactly one; everything else is yours to choose, and the tool explains what
   each choice means.
3. Click **Design saRNAs**.
4. Inspect the candidates, the promoter map and the reasons behind every verdict.

The first lookup of any gene needs internet. After that the gene is cached and works offline.

**Read [what the tool can and cannot tell you](docs/development/local-workbench.md) before
acting on its output.** In short: it is a filter and an explainer with every decision traceable
to a published source — it is **not** a predictor of whether a design will work. Its ranking has
been measured against 125 assayed sequences and does not beat chance. The tool shows you that
number next to the ranking rather than hiding it.

---

## What leaves your computer

Transparency about this matters, so here it is in full.

| When | What is sent | Where |
|---|---|---|
| Looking up a gene for the first time | The gene symbol or Ensembl ID | Ensembl (rest.ensembl.org) |
| Fetching promoter DNA | Genomic coordinates | Ensembl |
| **Only if you click the off-target button** | The 19-nucleotide candidate sequences | NCBI BLAST |
| Everything else | Nothing | — |

Your designs, settings, saved runs and exports stay in a folder on your own machine
(`~/.sarna-designer` on macOS/Linux, `C:\Users\<you>\.sarna-designer` on Windows). The tool
listens only on `127.0.0.1`, which means your own computer and nothing else on the network can
reach it.

The sequences sent to NCBI are drawn from the public human genome, so no novel sequence is
disclosed — but running a search does tell NCBI which gene you are working on. That is why the
off-target search is a button you press, not something that happens automatically.

---

## If something goes wrong

**`uv` is not recognised** — you did not close and reopen the terminal after step 2. Do that.

**`sarna-designer` is not recognised** — step 3 did not finish. Run it again and read the last
few lines for an error.

**The check says FAIL** — each failing line prints a `what to do:` with it. If it says to report
the output, copy the whole thing and send it on; it is written to be useful to whoever maintains
the tool.

**Windows warns about a firewall** — allow it. The tool only listens on your own machine.

**The browser page is blank** — make sure the terminal window is still open and shows no error,
then reload.

**Nothing above helped** — run this and send the output:

```
sarna-designer-check --online
```

---

## Updating

```
uv tool upgrade sarna-designer
```

## Uninstalling

```
uv tool uninstall sarna-designer
```

This removes the tool and its private Python. Your saved runs stay in `.sarna-designer`; delete
that folder too if you want them gone.

---

## Why it is installed this way

Three choices worth explaining, since you are running this on your own machine.

**Why `uv` rather than `pip`.** The tool is pinned to Python 3.11. Windows and macOS usually
ship a different version, so `pip install` would simply refuse. `uv` fetches the right Python
itself, in a private folder, so you never have to touch or change your system Python — and
nothing this tool installs can affect other software.

**Why it runs locally rather than as a website.** Recorded in
[ADR 0014](docs/adr/0014-local-delivery-and-wang-default.md),
[ADR 0015](docs/adr/0015-local-workbench-and-frozen-ui.md) and, with the case for and against a
hosted service set out in full, [ADR 0020](docs/adr/0020-localhost-or-a-hosted-service.md): your
sequences and designs stay on your machine, the tool works offline once a gene is cached, and
there is no server to maintain or trust. The short version is that the gene you are working on
*is* the research question, and running locally means only Ensembl — and NCBI, if you ask for an
off-target screen — ever learns it.

**Why there are two commands.** `sarna-designer-check` makes a broken installation say so
clearly instead of failing later in a way that looks like a scientific problem. It is quick and
makes no network calls unless you ask.

---

## Sharing reference data in a group (optional)

Everyone who looks up SCN1A downloads the same gene, transcript and sequence records from
Ensembl. If several of you work on the same genes, you can point the cache at a shared folder so
it is fetched once:

- **Windows (PowerShell):** `setx SARNA_CACHE_DIR "\\\\server\\share\\sarna-cache"`
- **macOS / Linux:** add `export SARNA_CACHE_DIR=/Volumes/share/sarna-cache` to your shell profile

Then restart the app. Whoever looks a gene up first fills the cache; everyone else reads it, and
can work offline afterwards.

This is safe to share while somebody else is using it — cache entries are written to a temporary
name and renamed, so you never read a half-written file. It holds only public Ensembl records,
no designs: your runs stay in your own `--data-dir`. If the shared folder is unreachable the tool
fetches from Ensembl as usual rather than failing.

---

## The short version, for people who know Python

```bash
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
sarna-designer-check --online
sarna-designer --port 8001
```

- Requires Python 3.11 exactly; `[app]` adds `httpx` on top of `pydantic`. No compiler needed —
  both ship universal wheels.
- `[genomics]` additionally pulls `pyBigWig`, `numpy` and `pyliftover` for local AF22 histone
  tracks. You do not need it unless you have those bigWigs.
- State lives in `--data-dir` (default `~/.sarna-designer`), one SQLite file.
- Environment variables are documented in `src/sarna/settings.py`; `SARNA_CONTACT_EMAIL` is the
  address Ensembl asks you to send with queries, and `SARNA_OFFLINE=1` forbids all network use.
- Developing on it instead? `git clone`, then `.venv/bin/python -m pytest -q`. See
  [CLAUDE.md](CLAUDE.md).
