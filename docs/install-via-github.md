g# Installing and updating saRNA Designer from GitHub on Windows

Follow this guide if you have access to the private saRNA Designer GitHub repository.
It works for a first installation and for replacing an installation from a `.whl` file.
Read access to the repository is enough.

The app runs on your computer. You open it in your normal browser after starting it in
PowerShell. Allow about fifteen minutes for the first installation and keep your internet
connection on. You do not need to install Python or WSL yourself.

## 1. Check your GitHub access

Sign in to GitHub with the account that received the invitation. Open:

[JoepVondervoort/sarna-designer](https://github.com/JoepVondervoort/sarna-designer)

You should see the repository files, including `README.md`. If you cannot see them, check
whether you still need to accept the invitation. Use that same account when signing in
later during installation.

## 2. Open PowerShell

Press the Windows key, type `powershell`, and press Enter.

Copy the whole command from each code box, paste it into PowerShell, and press Enter.
Run one command at a time. Wait until the prompt, usually starting with `PS`, returns
before entering the next command. Do not copy example output or the code-box heading.

If PowerShell displays `>>`, press Ctrl+C and paste the complete command again.

## 3. Install uv

`uv` installs the app and downloads the Python version it needs. If you already have uv,
skip the installation command and check its version below.

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

If your organisation blocks this command, try:

```powershell
winget install --id=astral-sh.uv -e
```

Close PowerShell completely and open it again. Then check:

```powershell
uv --version
```

It should print a version number. If it says the command is not recognised, the installation
has not succeeded or you need to open a fresh PowerShell window.

## 4. Install Git for Windows

If Git is already installed, check it in step 5 and skip this installation.

Download the installer from [Git for Windows](https://git-scm.com/download/win).
Run it and accept the default options. Keep Git Credential Manager enabled: it handles
signing in to the private GitHub repository.

You do not need to clone the repository manually. The next command downloads it for you.

## 5. Open a fresh PowerShell window and check both tools

Close PowerShell and open it again. Run these separately:

```powershell
git --version
```

```powershell
uv --version
```

Both should print version numbers. If either command is not recognised, resolve that
installation before continuing.

## 6. Install saRNA Designer from GitHub

Copy this whole line into PowerShell:

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

If a browser sign-in window opens, sign in with the GitHub account that can see the
repository, then return to PowerShell. If you are already signed in through Git Credential
Manager, it may continue without asking again.

The first installation may take several minutes. Wait for the PowerShell prompt to return.
A successful installation lists these two commands:

```text
sarna-designer
sarna-designer-check
```

If you previously installed from a `.whl`, this command switches the installation to
GitHub so future updates can use the upgrade command. Your saved runs are kept.

If installation reports an error, save the complete command and output to share with the
person helping you. Do not continue until installation succeeds.

## 7. Check the installation

Close PowerShell and open it again. Run:

```powershell
sarna-designer-check
```

The output should end with `RESULT: PASS`. It also prints the installed version.
The version number can stay the same when code changes, so an unchanged version alone
is not evidence that installation failed.

## 8. Start the app and open its page

In PowerShell, run:

```powershell
sarna-designer
```

Leave this window open. Open your normal browser and enter:

[http://localhost:8001/](http://localhost:8001/)

This address opens the app running on your own computer. It works while the PowerShell
process is running. You can bookmark it.

Enter a gene name or paste a promoter DNA sequence, choose your promoter reference and
settings, and start a design. Retrieving a new gene needs internet access; cache-only mode
works only when the required data has already been downloaded. If the app reports a cache
miss, use its online retrieval option and retry.

To stop the app, return to PowerShell and press Ctrl+C.

## Opening the app next time

You do not repeat the installation. Open PowerShell and run:

```powershell
sarna-designer
```

Then open [http://localhost:8001/](http://localhost:8001/) in your browser.
Keep PowerShell open while you work.

## Getting an update after new changes are pushed

Updates are manual. This installation follows the repository's default branch, `main`.
Changes must have been pushed to that branch before this command can retrieve them.

1. If the app is running, wait for your design or screening jobs to finish. Stop the app
   with Ctrl+C in its PowerShell window.
2. Run:

   ```powershell
   uv tool upgrade sarna-designer
   ```

3. Wait for it to finish, then check:

   ```powershell
   sarna-designer-check
   ```

4. If the check ends with `RESULT: PASS`, start the app again:

   ```powershell
   sarna-designer
   ```

5. Open the browser page or refresh the existing tab. If it still displays the old
   interface, press Ctrl+F5 to reload it.

The update follows Git commits even when the package version number has not changed.
Your saved runs remain in `%USERPROFILE%\.sarna-designer`, separate from the installed
program. Old runs retain their original results; start a new run to use the updated engine.

## If something goes wrong

| Message or problem | What to do |
|---|---|
| `repository not found` or access denied | Check step 1 and sign in with the invited GitHub account. |
| `git` or `uv` is not recognised | Reopen PowerShell. If it still fails, complete the corresponding installation. |
| Git asks for a password or authentication fails | Check that Git Credential Manager is enabled and use its GitHub sign-in flow. |
| `sarna-designer` is not recognised | Reopen PowerShell and check that step 6 succeeded. |
| The browser cannot connect to localhost | Start `sarna-designer` and keep its PowerShell window open. Read any startup error. |
| Port 8001 is already in use | An app may already be running. Alternatively run `sarna-designer --port 8002` and open `http://localhost:8002/`. |
| Upgrade says `Nothing to upgrade` after a new update | Check that the changes reached `main`. If your installation still came from a wheel, repeat step 6 to switch to GitHub. |
| Anything else | Share the complete error output and the output of `sarna-designer-check`. |

## Set up fast local off-target screening (once per computer)

After installing or updating the app, start it normally. No additional PowerShell command,
Python installation or WSL installation is needed for this setup on Windows x64.

1. Open **Local off-target screening · reference setup** near the top of the app.
2. Click **Download & prepare local references**. The app downloads about 1.2 GB from NCBI
   and the official Bowtie release, verifies the files, and builds search indexes locally.
   Have at least **20 GiB free disk space**. The first setup can take an hour or more on a modest computer;
   the time depends on your connection and computer. Keep the app process running.
3. Wait for **Local human references are ready**, then use **Local reference** in the
   off-target panel. Genome and transcript searches use separate indexes. Candidate sequences
   are not sent to NCBI in this mode.

You can click **Pause setup** and later **Resume / retry reference setup**. Interrupted
file downloads resume when the server supports byte ranges; already verified downloads are
reused. An interrupted index build is rebuilt safely. Closing the browser alone does not stop
setup while the app process continues running.

References are stored under `%USERPROFILE%\.sarna-designer\references`, outside the installed
app. Updating the app with `uv` keeps them and does **not** automatically download the genome
again. **Check reference status / updates** compares your installation with the reference
release supplied by your installed app. If a later app version supplies a new reference,
click **Download reference update** to install it explicitly. The previous active reference
remains usable until preparation of the new release finishes. It does not search online for
an unpinned “latest” dataset.

The first supplied reference is **GRCh38.p14, GCF_000001405.40, RefSeq annotation release 110**.
This is a pinned, reproducible release, not a claim to have the latest annotation. Source URLs,
checksums, tool version and search policy travel in the JSON screening record. Linux x64 is
also supported. Native Windows execution still needs validation on a Windows machine; the
installer uses the official Windows Bowtie 1.2.2 release, while Linux uses Bowtie 1.3.1.

Local screening searches **full-length 19-base cores with 0–3 mismatches on both strands**.
It does not assess gaps, bulges, shorter seed matches, personal genetic variants or biological
binding. All matching alignments are counted; up to 200 detail alignments per candidate are
stored for display. Very repetitive batches that exceed the resource budget fail explicitly
instead of producing a truncated “clean” result. NCBI BLAST remains an explicitly selected,
remote cross-check with different search behavior and potentially long queues.

### Build a shortlist after screening

The setup panel opens automatically when references are missing. It lists the genome
(972.9 MB), transcripts (132.9 MB), annotations (56.9 MB), chromosome map (79 kB) and the
platform's search tools (31.3 MB on Windows; 52.4 MB on Linux). The total is about 1.2 GB.
Downloading and building the indexes are separate steps; building can take an hour or more.

After a design run, choose **Screen all locally · both scopes** in **Your screened shortlist**.
The panel offers a top **10 or 15**, ordered by original design rank after screening exclusions.
Both local searches must be complete using the same reference. Older reports may need a refresh;
cached local alignments are reused. Missing evidence never counts as a negative finding.

The initial project policy excludes any candidate with an 18/19 or 19/19 match in either scope,
or five or more 17/19 matches. Genomic counts exclude the intended locus when its coordinates
are known. Counts concern reference alignments, including isoforms and alternate records.
**These are conservative project choices, not validated saRNA safety thresholds.** Read the
selection rules and linked papers in the panel. The shortlist may contain fewer than ten,
and is provisional while other passing candidates still lack evidence.

**Select these for export** uses only eligible candidates. Later adverse findings remove them
from that managed selection. Sequence exports cannot silently fall back to all designs when
it becomes empty. **Download selection record** preserves every decision, source search,
reference, counts and policy. Original scores and manual research notes remain available.
