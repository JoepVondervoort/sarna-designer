# Installing saRNA Designer

**Written for:** a colleague who wants to use the tool, with no assumption of programming
experience. If you are comfortable with Python packaging, [the short version](#the-short-version-for-people-who-know-python)
is at the bottom.

saRNA Designer runs **on your own computer**. It is not a website, there is no account, and
nothing is uploaded except the gene lookups described under [What leaves your computer](#what-leaves-your-computer).

You need about **ten minutes** and an internet connection for the setup.

---

## Before you start

### What you need to install first: on Windows, nothing

This is the most common question, so here it is as a list.

| | Needed? | |
|---|---|---|
| **WSL** (Windows Subsystem for Linux) | **No** | Not used, with [one exception](#optional-af22-histone-tracks-the-one-thing-that-needs-wsl) almost nobody needs. |
| **Python** | **No** | Step 2 installs a tool that fetches the right Python itself. |
| **Administrator rights** | **No** | Everything installs under your own user account. |
| **Docker**, Anaconda, Visual Studio | **No** | None of them. |
| **Git** | **Only for Route B** | [Git for Windows](https://git-scm.com/download/win). Most people do not need it. |

The only thing you install is a small program called **`uv`**
([what it is](https://docs.astral.sh/uv/)), in step 2, and it is one line.

### Where you type the commands: PowerShell

Everything on Windows happens in **PowerShell**. Not Command Prompt, not a WSL/Ubuntu window,
not inside any other program. Step 1 shows how to open it.

> **If you already use WSL for something else:** do not install this inside it. The tool is
> installed on Windows and opens in your normal Windows browser. A copy installed inside WSL
> will not be on your Start menu and its saved runs will be somewhere you cannot easily find.

### How you will get the tool — read this first

There are two routes, and **which one applies to you is not a preference**:

| | **Route A — a file you were sent** | **Route B — from GitHub** |
|---|---|---|
| You need | the `.whl` file, nothing else | a GitHub account **with access to the repository**, and **Git installed** |
| Works for | anyone | collaborators on the repository only |
| Use it when | you were emailed or handed a file | you have been added to the repository |

**The repository is private.** A GitHub link will fail with *"repository not found"* for
anyone who has not been given access — that is what a private repository looks like from
outside, and it is not a sign that anything is broken. If you were not explicitly added to it,
you want **Route A**.

Route A is also the only one that needs no Git and no GitHub account at all, which is why it
is written first below.

> **If you are the person sending the file**, see
> [Sending it to a colleague](#sending-it-to-a-colleague) at the end — it tells you how to
> build the file and exactly where it appears.

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

This installs a small program called `uv` that manages the installation for you. If your
organisation blocks that command, `winget install --id=astral-sh.uv -e` does the same thing.

**Now close PowerShell completely and open it again.** This is necessary — the new program is
only visible to a freshly opened window.

### Step 3 — Install saRNA Designer

**Route A — from the file you were sent.** Save the `.whl` file somewhere you can find it,
for example your **Downloads** folder. Then type `uv tool install --python 3.11 "` , drag the
file from Explorer into the PowerShell window (this pastes its full path for you), and finish
the line with `[app]"` so it reads like this:

```powershell
uv tool install --python 3.11 "C:\Users\YourName\Downloads\sarna_designer-0.1.0.dev0-py3-none-any.whl[app]"
```

The `[app]` on the end is **not optional** — it is what installs the browser interface.

**Route B — from GitHub**, only if you have been added to the repository and have
[Git](https://git-scm.com/download/win) installed:

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

Either way, this downloads the tool and the correct Python, and ends with:

```
Installed 2 executables: sarna-designer, sarna-designer-check
```

> **This step is slow the first time, and that is normal.** It is fetching a complete private
> Python — about 90 MB unpacked across **3,862 files** — and on Windows the virus scanner
> inspects every one of those files as it is written. **Several minutes is expected**, longer
> on a slow or corporate network. It only happens once: installing again later takes under a
> second, because the download is cached.
>
> **Is it stuck or just slow?** Open Task Manager (<kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>Esc</kbd>)
> and look for **uv.exe** — if it is using CPU, disk or network, leave it alone. Worry only if
> it has sat completely idle for more than about fifteen minutes, or if the text stops at
> `Resolving` and never moves, which usually means a proxy or firewall is blocking the
> download.

**If you do not see that line, stop here** and read the last few lines of output — the rest of
these steps will not work.

### Step 3b — Close and reopen PowerShell, again

`uv` may print a warning that its folder *"is not on your PATH"*. Closing PowerShell and
opening it again fixes that. Do it even if you saw no warning; it costs five seconds and it is
the single most common reason step 4 fails.

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

From the file you were sent (Route A — drag the `.whl` into Terminal to paste its path):

```bash
uv tool install --python 3.11 "/Users/you/Downloads/sarna_designer-0.1.0.dev0-py3-none-any.whl[app]"
```

Or from GitHub (Route B), if you have access to the repository and `git` installed:

```bash
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

It should end with `Installed 2 executables`. Then close Terminal and open it again, so the
new commands are on your PATH.

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

**`sarna-designer` is not recognised** — almost always step 3b: close PowerShell or Terminal
completely and open a new one. If it still is not found, run `uv tool update-shell` and open
another new window. If that does not help either, step 3 did not finish — run it again and
read the last few lines for an error.

**`repository not found`, or Git asks for a password** — the repository is private and your
account has not been given access. Use Route A: ask for the `.whl` file.

**`git` is not recognised** — Route B needs [Git for Windows](https://git-scm.com/download/win).
Route A does not need it at all.

**The check says `httpx is missing` or `the local web interface is incomplete`** — the `[app]`
was left off the end of the command in step 3. This is the most common single mistake. Run
`uv tool uninstall sarna-designer` and do step 3 again, with `[app]` on the end.

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

If you installed from GitHub (Route B):

```
uv tool upgrade sarna-designer
```

If you installed from a file (Route A), there is nothing for `uv` to upgrade *from* — ask for
the new `.whl` and install it over the top with the same step 3 command. Your saved runs are
kept.

Either way, run `sarna-designer-check` afterwards.

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

## Optional: AF22 histone tracks, the one thing that needs WSL

Skip this unless you have been given the AF22 bigWig files and want the tool to read them.
**The design tool does not need any of it** — without the tracks, chromatin simply reports
"not loaded", which is an honest answer rather than a broken one.

The reader for those files, [pyBigWig](https://pypi.org/project/pyBigWig/), publishes **no
Windows wheel** and lists MacOS, POSIX and Unix as its supported platforms; building it needs
`libcurl` and `zlib` headers. So on Windows this extra cannot be installed natively, and WSL
is the way round it.

**1. Install WSL** from an **administrator** PowerShell — right-click PowerShell in the Start
menu and choose *Run as administrator*:

```powershell
wsl --install
```

Reboot when it asks. See [Microsoft's guide](https://learn.microsoft.com/windows/wsl/install)
if it does not complete.

**2. Open the Ubuntu window** it installed and follow the Linux steps, with `[app,genomics]`
in place of `[app]`:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# then, from a new terminal:
uv tool install --python 3.11 "/path/to/sarna_designer-0.1.0.dev0-py3-none-any.whl[app,genomics]"
```

**3. Point it at the tracks**, then start it:

```bash
export SARNA_AF22_DIR=/path/to/af22_data
sarna-designer
```

The browser still opens on Windows at the printed address. Note that a copy installed this way
keeps its saved runs inside WSL, not in your Windows user folder.

---

## Sending it to a colleague

*For whoever is distributing the tool, not for the person installing it.*

The repository is private, so the simplest thing to send is a single file: a **wheel** (`.whl`).
It is about 320 kB, it contains the whole tool, and the person receiving it needs no GitHub
account, no Git and no Python.

### 1. Build it

From a checkout, on any platform:

```bash
.venv/bin/python -m build --wheel
```

If that says `No module named build`, install it once with
`.venv/bin/python -m pip install build`.

It prints `Successfully built sarna_designer-0.1.0.dev0-py3-none-any.whl` and puts that one
file in the **`dist/`** folder next to the code.

### 2. Find it

```
<your checkout>/dist/sarna_designer-0.1.0.dev0-py3-none-any.whl
```

**If you develop inside WSL**, that path is inside Linux and Outlook cannot see it. Two ways
to reach it from Windows:

- Run `explorer.exe dist` from the project folder in your WSL terminal. A normal Windows
  Explorer window opens on the folder, and you can drag the file straight into an email.
- Or open Explorer and paste this into the address bar, substituting your distribution name
  and user:
  ```
  \\wsl$\Ubuntu\home\jvon\sarna-designer\dist
  ```
  (On newer Windows builds `\\wsl.localhost\Ubuntu\...` also works.)

### 3. Send it

Attach the `.whl` and link them to this file. A message that works:

> Here is the saRNA design tool. Save the attached file somewhere you can find it, such as your
> Downloads folder, then follow the Windows section of INSTALL.md — it is written for someone
> who has never installed anything like this.
>
> Three things worth knowing before you start: you do **not** need WSL, Python or
> administrator rights; everything happens in **PowerShell**; and at step 3 use **Route A**,
> the one that installs from the file I sent, not the GitHub one.
>
> When step 4 prints `RESULT: PASS` it is working. If anything says FAIL, send me the whole
> output and I will tell you what it means.

### If you would rather they install from GitHub

They need to be added as a collaborator on the repository **and** to install
[Git for Windows](https://git-scm.com/download/win) first. That is two extra things that can
go wrong, which is why the file is the default recommendation. Making the repository public
would remove both, and is a decision about publishing the work rather than a technical one.

---

## The short version, for people who know Python

```bash
# From a wheel you were sent (works without Git or repository access):
uv tool install --python 3.11 "./sarna_designer-0.1.0.dev0-py3-none-any.whl[app]"
# Or, with access to the private repository and git on PATH:
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"

sarna-designer-check --online
sarna-designer --port 8001
```

- The repository is **private**: the git URL resolves only for collaborators, and needs `git`
  installed. The wheel needs neither. Build one with `python -m build --wheel`.
- Requires Python 3.11 exactly; `[app]` adds `httpx` on top of `pydantic`. No compiler needed —
  both ship universal wheels.
- `[genomics]` additionally pulls `pyBigWig`, `numpy` and `pyliftover` for local AF22 histone
  tracks. **Not installable on Windows**: pyBigWig publishes no Windows wheel and lists only
  MacOS, POSIX and Unix, needing `libcurl` and `zlib` headers to build. Use WSL or a Mac for
  that extra. You do not need it unless you have those bigWigs.
- State lives in `--data-dir` (default `~/.sarna-designer`), one SQLite file.
- Environment variables are documented in `src/sarna/settings.py`; `SARNA_CONTACT_EMAIL` is the
  address Ensembl asks you to send with queries, and `SARNA_OFFLINE=1` forbids all network use.
- Developing on it instead? `git clone`, then `.venv/bin/python -m pytest -q`. See
  [CLAUDE.md](CLAUDE.md).
