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

### How you will get the tool: read this first

There are two routes, and **which one applies to you is not a preference**:

| | **Route A: a file you were sent** | **Route B: from GitHub** |
|---|---|---|
| You need | the `.whl` file, nothing else | a GitHub account **with access to the repository**, and **Git installed** |
| Works for | anyone | collaborators on the repository only |
| Use it when | you were emailed or handed a file | you have been added to the repository |

**The repository is private.** A GitHub link will fail with *"repository not found"* for
anyone who has not been given access. That is what a private repository looks like from
outside, and it is not a sign that anything is broken. If you were not explicitly added to it,
you want **Route A**.

Route A is also the only one that needs no Git and no GitHub account at all, which is why it
is written first below.

> **If you are the person sending the file**, see
> [Sending it to a colleague](#sending-it-to-a-colleague) at the end. It tells you how to
> build the file and exactly where it appears.

---

## Windows

Most people will use this. Everything happens in a program called **PowerShell**, which is
already on your computer.

### Step 1: Open PowerShell

Press the **Windows key**, type `powershell`, and press **Enter**. A blue window opens with a
blinking cursor. You type commands here and press Enter after each one.

Copy the **whole command** from each code box, paste it into PowerShell, then press
**Enter**. Do not copy a prompt such as `PS C:\...>` or the example output. Run one command
at a time, and wait for the `PS ...>` prompt to return before running the next one.

If you see `>>`, PowerShell is waiting for an unfinished command (often a missing quote).
Press **Ctrl+C**, then paste the complete command again. Do not type a command in pieces or
drag a file into a partially typed command.

### Step 2: Install `uv`

Copy this whole line, paste it into PowerShell (right-click pastes), and press **Enter**:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

This installs a small program called `uv` that manages the installation for you. If your
organisation blocks that command, `winget install --id=astral-sh.uv -e` does the same thing.

**Now close PowerShell completely and open it again.** This is necessary, because the new program is
only visible to a freshly opened window.

### Step 3: Install saRNA Designer

**Route A: from the file you were sent (recommended for colleagues).**

**A. Save the attachment.** You need a **file**, called
`sarna_designer-0.2.0-py3-none-any.whl`, in your **Downloads** folder.

- If you received the `.whl` directly, use your email's **Save attachment** option. Do not
  open or extract it.
- If you received a separate ZIP containing the wheel, save the ZIP, choose **Extract All**
  once, then copy the `.whl` file inside into Downloads. **Do not extract the inner `.whl`.**

A wheel is a ZIP archive internally. An archive icon is okay; an extracted folder is not.
In Explorer, enable **View → Show → File name extensions** (on older Windows, **View → File
name extensions**). If the actual wheel file was renamed to `.whl.zip`, remove only the final
`.zip`. Renaming an extracted folder to `.whl` does not turn it into a wheel.

**B. Check that the file is there.** Paste this whole line and press Enter:

```powershell
Test-Path -LiteralPath "$env:USERPROFILE\Downloads\sarna_designer-0.2.0-py3-none-any.whl" -PathType Leaf
```

It must print **`True`**. If it prints `False`, stop: the file is missing, is a folder, or has
another name. Check Downloads and the filename. Windows sometimes adds `(1)` to duplicate
attachments. The commands here assume the exact filename above; for a newer release, replace
only that filename in both commands with the one you received.

**C. Install it.** Copy this **entire line**, including `uv tool install --python 3.11`,
and press Enter:

```powershell
uv tool install --python 3.11 "$env:USERPROFILE\Downloads\sarna_designer-0.2.0-py3-none-any.whl[app]"
```

`$env:USERPROFILE` automatically uses your Windows user folder; leave it exactly as written.
The `[app]` inside the quotes installs the browser interface and must stay there. Pasting only
the quoted path will just print the path; it does **not** install anything.

Already have Python 3.14 or another version? Leave it installed. This command selects Python
3.11 separately, as required by this tool.

**Route B: from GitHub**, only if you have been added to the repository and have
[Git](https://git-scm.com/download/win) installed:

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

Either way, this downloads the tool and the correct Python, and ends with:

```
Installed 2 executables: sarna-designer, sarna-designer-check
```

> The first installation may take several minutes while it downloads Python and packages.
> Wait for the command prompt to return. If it appears stuck, press **Ctrl+C** and rerun the
> same complete install command with `-v` immediately after `uv` (`uv -v tool install ...`).
> Send the command and its output to the person who sent you the tool. Task Manager activity
> alone cannot tell you whether the command was entered correctly or why it is waiting.

**If you do not see that line, stop here** and read the last few lines of output. The rest of
these steps will not work.

### Step 3b: Close and reopen PowerShell, again

`uv` may print a warning that its folder *"is not on your PATH"*. Closing PowerShell and
opening it again fixes that. Do it even if you saw no warning; it costs five seconds and it is
the single most common reason step 4 fails.

### Step 4: Check that it works

```powershell
sarna-designer-check
```

You should see a list of checks and then **`RESULT: PASS`**. If any line says `FAIL`, the
check tells you what to do. See [If something goes wrong](#if-something-goes-wrong).

### Step 5: Start it

```powershell
sarna-designer
```

It prints a web address like `http://localhost:8001/`. Open that in your browser. **Leave the
PowerShell window open** while you work. Closing it stops the tool.

To stop it, click the PowerShell window and press **Ctrl+C**.

---

## macOS

### Step 1: Open Terminal

Press **Cmd+Space**, type `terminal`, press **Enter**.

### Step 2: Install `uv`

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then close Terminal and open it again.

### Step 3: Install saRNA Designer

From the file you were sent (Route A, and you can drag the `.whl` into Terminal to paste its path):

```bash
uv tool install --python 3.11 "/Users/you/Downloads/sarna_designer-0.2.0-py3-none-any.whl[app]"
```

Or from GitHub (Route B), if you have access to the repository and `git` installed:

```bash
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

It should end with `Installed 2 executables`. Then close Terminal and open it again, so the
new commands are on your PATH.

### Steps 4 and 5: Check and start

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

1. Type a human gene symbol such as `SCN1A`, `TP53` or `KCNQ2`, then click **Look up**.
2. Choose which promoter reference to design against. The tool picks the MANE Select transcript
   when there is exactly one; everything else is yours to choose, and the tool explains what
   each choice means.
3. Click **Design saRNAs**.
4. Inspect the candidates, the promoter map and the reasons behind every verdict.

The first lookup of any gene needs internet. After that the gene is cached and works offline.

**Read [what the tool can and cannot tell you](docs/development/local-workbench.md) before
acting on its output.** In short: it is a filter and an explainer with every decision traceable
to a published source. It is **not** a predictor of whether a design will work. Its ranking has
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
| Everything else | Nothing | n/a |

Your designs, settings, saved runs and exports stay in a folder on your own machine
(`~/.sarna-designer` on macOS/Linux, `C:\Users\<you>\.sarna-designer` on Windows). The tool
listens only on `127.0.0.1`, which means your own computer and nothing else on the network can
reach it.

The sequences sent to NCBI are drawn from the public human genome, so no novel sequence is
disclosed, but running a search does tell NCBI which gene you are working on. That is why the
off-target search is a button you press, not something that happens automatically.

---

## If something goes wrong

**`error: ... <PACKAGE>` is missing**. Only the beginning of the install command was
entered. Copy the whole line from Step 3C, including the quoted path and `[app]`.

**PowerShell prints the file path and nothing else**. Only the quoted path was entered.
Run the whole Step 3C command, starting with `uv tool install`.

**`neither pyproject.toml nor setup.py are present in the directory`**. Check whether the
`.whl` is actually an extracted folder. Run the Step 3B file check; it must print `True`.
Rename the extracted folder, then save the original wheel attachment to Downloads. Do not
create either of those missing files. If the check is `True` and the error persists, send the
full command and output to the sender.

**PowerShell shows `>>`**. Press **Ctrl+C** and paste the complete command again; see Step 1.

**`uv` is not recognised**. You did not close and reopen the terminal after step 2. Do that.

**`sarna-designer` is not recognised**. Almost always step 3b: close PowerShell or Terminal
completely and open a new one. If it still is not found, run `uv tool update-shell` and open
another new window. If that does not help either, step 3 did not finish. Run it again and
read the last few lines for an error.

**`repository not found`, or Git asks for a password**. The repository is private and your
account has not been given access. Use Route A: ask for the `.whl` file.

**`git` is not recognised.** Route B needs [Git for Windows](https://git-scm.com/download/win).
Route A does not need it at all.

**The check says `httpx is missing` or `the local web interface is incomplete`**. The `[app]`
was left off the end of the command in step 3. This is the most common single mistake. Run
`uv tool uninstall sarna-designer` and do step 3 again, with `[app]` on the end.

**The check says FAIL**. Each failing line prints a `what to do:` with it. If it says to report
the output, copy the whole thing and send it on; it is written to be useful to whoever maintains
the tool.

**Windows warns about a firewall**. Allow it. The tool only listens on your own machine.

**The browser page is blank**. Make sure the terminal window is still open and shows no error,
then reload.

**Nothing above helped**. Run this and send the output:

```
sarna-designer-check --online
```

---

## Updating to a newer version

The right command depends on which route you installed with in step 3. Using the wrong one is
not harmless: one of them reports success while leaving you on the old version.

### If you installed from GitHub (Route B)

```
uv tool upgrade sarna-designer
```

That is the whole update. `uv` remembers which repository you installed from and fetches
whatever is there now. It follows the **code**, not the version number, so this works even when
the version has not changed.

### If you installed from a file (Route A)

**`uv tool upgrade` will not work, and will not tell you that it did not work.** It prints
`Nothing to upgrade`, because there is no repository for it to look at. That message does not
mean you are up to date.

Instead, ask for the new `.whl` file, save it to **Downloads** as before, and run the step 3
install command again with the new filename:

```
uv tool install --python 3.11 "$env:USERPROFILE\Downloads\sarna_designer-0.2.0-py3-none-any.whl[app]"
```

On macOS or Linux, use the path to where you saved it:

```
uv tool install --python 3.11 "/Users/you/Downloads/sarna_designer-0.2.0-py3-none-any.whl[app]"
```

Installing over the top replaces the old version. You do not need to uninstall first.

> **Delete the old `.whl` from Downloads once the new one works.** Two files with different
> version numbers in the same folder is the easiest way to install the old one by accident.

### Afterwards, either way

Close the terminal, open a new one, and run:

```
sarna-designer-check
```

The last line names the version you now have. If it still shows the old number, the install did
not take: check that you used the new filename, then run the command again.

### Your saved work is not touched

Runs, snapshots and exports live in a `.sarna-designer` folder, separate from the program.
Updating, reinstalling and uninstalling all leave them alone.

### Which version am I on?

`sarna-designer-check` prints it, and so does the footer of the page in your browser. Every run
you save records the same version, so a result from today can be told apart from one produced by
an older build. Quote that number if you report a problem.

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
itself, in a private folder, so you never have to touch or change your system Python, and
nothing this tool installs can affect other software.

**Why it runs locally rather than as a website.** Recorded in
[ADR 0014](docs/adr/0014-local-delivery-and-wang-default.md),
[ADR 0015](docs/adr/0015-local-workbench-and-frozen-ui.md) and, with the case for and against a
hosted service set out in full, [ADR 0020](docs/adr/0020-localhost-or-a-hosted-service.md): your
sequences and designs stay on your machine, the tool works offline once a gene is cached, and
there is no server to maintain or trust. The short version is that the gene you are working on
*is* the research question. Running locally means only Ensembl, plus NCBI if you ask for an
off-target screen, ever learns it.

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

This is safe to share while somebody else is using it, because cache entries are written to a temporary
name and renamed, so you never read a half-written file. It holds only public Ensembl records,
no designs: your runs stay in your own `--data-dir`. If the shared folder is unreachable the tool
fetches from Ensembl as usual rather than failing.

---

## Optional: AF22 histone tracks, the one thing that needs WSL

Skip this unless you have been given the AF22 bigWig files and want the tool to read them.
**The design tool does not need any of it**. Without the tracks, chromatin simply reports
"not loaded", which is an honest answer rather than a broken one.

The reader for those files, [pyBigWig](https://pypi.org/project/pyBigWig/), publishes **no
Windows wheel** and lists MacOS, POSIX and Unix as its supported platforms; building it needs
`libcurl` and `zlib` headers. So on Windows this extra cannot be installed natively, and WSL
is the way round it.

**1. Install WSL** from an **administrator** PowerShell. Right-click PowerShell in the Start
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
uv tool install --python 3.11 "/path/to/sarna_designer-0.2.0-py3-none-any.whl[app,genomics]"
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

It prints `Successfully built sarna_designer-0.2.0-py3-none-any.whl` and puts that one
file in the **`dist/`** folder next to the code.

### 2. Find it

```
<your checkout>/dist/sarna_designer-0.2.0-py3-none-any.whl
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

For colleagues new to this, put the **actual wheel file** inside a separate ZIP, for example
`sarna-designer-install.zip`, and attach that ZIP together with this guide. This makes the
instruction explicit: extract the outer ZIP once, then leave the inner wheel intact. Do not
send a folder named `.whl`, and do not extract the wheel before packaging it.

From the project folder, you can create the outer ZIP with:

```bash
.venv/bin/python -m zipfile -c dist/sarna-designer-install.zip dist/sarna_designer-0.2.0-py3-none-any.whl
```

Open `dist` in Explorer as described above and attach `sarna-designer-install.zip`.
The wheel path is a file to attach, **not a command to execute** in your terminal.

A message you can send:

> Attached are the saRNA design tool and its installation guide. Save the ZIP and choose
> **Extract All** once. Copy the `.whl` file inside into your **Downloads** folder; do not
> open or extract that inner file.
>
> Follow **Windows → Route A** in INSTALL.md. You do not need Python, WSL, administrator
> rights, Git or a GitHub account. An existing Python installation can stay as it is.
> Paste each complete command into PowerShell and press Enter; do not type it in pieces.
>
> The file check must print `True`, the installation should report `Installed 2 executables`,
> and the final check should end with `RESULT: PASS`. If any step fails, stop and send me the
> command and its complete output.

### If you would rather they install from GitHub

Worth doing if they will keep receiving updates, because it turns every future update into one
command they run themselves instead of a file you build and send. The cost is two one-time
steps: they need a GitHub account added as a collaborator, and
[Git for Windows](https://git-scm.com/download/win) installed.

[Switching an existing install to GitHub](#switching-an-existing-install-to-github) is the
walkthrough to send them. Making the repository public would remove both steps, and is a
decision about publishing the work rather than a technical one.

---

## Switching an existing install to GitHub

If you already have this working from a `.whl` file and want updates to be one command,
[`docs/install-via-github.md`](docs/install-via-github.md) is the walkthrough. It is written to
be sent on its own.

In outline: accept the repository invitation, install
[Git for Windows](https://git-scm.com/download/win) accepting every default, then run the
Route B install command above once and sign in when the browser opens. It replaces your
current installation in place, and your saved runs are not touched. After that, every update
is `uv tool upgrade sarna-designer`.

You do not need WSL for any of it.

---

## The short version, for people who know Python

```bash
# From a wheel you were sent (works without Git or repository access):
uv tool install --python 3.11 "./sarna_designer-0.2.0-py3-none-any.whl[app]"
# Or, with access to the private repository and git on PATH:
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"

sarna-designer-check --online
sarna-designer --port 8001
```

- The repository is **private**: the git URL resolves only for collaborators, and needs `git`
  installed. The wheel needs neither. Build one with `python -m build --wheel`.
- Requires Python 3.11 exactly; `[app]` adds `httpx` on top of `pydantic`. No compiler needed,
  both ship universal wheels.
- `[genomics]` additionally pulls `pyBigWig`, `numpy` and `pyliftover` for local AF22 histone
  tracks. **Not installable on Windows**: pyBigWig publishes no Windows wheel and lists only
  MacOS, POSIX and Unix, needing `libcurl` and `zlib` headers to build. Use WSL or a Mac for
  that extra. You do not need it unless you have those bigWigs.
- State lives in `--data-dir` (default `~/.sarna-designer`), one SQLite file.
- Environment variables are documented in `src/sarna/settings.py`; `SARNA_CONTACT_EMAIL` is the
  address Ensembl asks you to send with queries, and `SARNA_OFFLINE=1` forbids all network use.
- Developing on it instead? `git clone`, then `.venv/bin/python -m pytest -q`. See
  [the development guide](docs/development/environment.md).
