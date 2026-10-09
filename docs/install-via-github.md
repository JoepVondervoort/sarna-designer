# Switching saRNA Designer to install from GitHub

**For a colleague who already has saRNA Designer working from a `.whl` file** and would like
updates to be one command instead of a file sent by email.

It replaces the installation in place. Nothing is uninstalled first, and your saved runs are
not touched. Budget about fifteen minutes, most of it waiting for two downloads.

You do **not** need WSL. You do **not** need to understand Python. You already have the one
tool that matters, `uv`, from the last installation.

---

## What changes

| | Today | After this |
|---|---|---|
| Getting a new version | wait for a `.whl`, save it, retype the install command with the new filename | `uv tool upgrade sarna-designer` |
| One-time setup | none | a GitHub account, an invitation, and Git for Windows |

`uv tool upgrade` does nothing useful for a file installation. It prints `Nothing to upgrade`,
which reads as "you are up to date" when you are not. That is the problem this fixes.

---

## 1. Accept the repository invitation

The person who sent you the tool adds you to the repository. You get an email from GitHub with
an **Accept invitation** button. Click it.

If you have no GitHub account yet, make a free one at [github.com](https://github.com) first,
then tell them the username you chose so the invitation goes to the right place.

**Nothing below works until this is done.** The repository is private, and without the
invitation step 4 fails with `repository not found`.

## 2. Install Git for Windows

Download it from [git-scm.com/download/win](https://git-scm.com/download/win). The download
starts on its own.

Run the installer and **click Next on every screen, changing nothing**. The defaults are the
ones you want. One of them, Git Credential Manager, is what lets you sign in to GitHub through
a browser window in step 4, instead of typing a password that GitHub no longer accepts for
this.

## 3. Close PowerShell and open it again

Press the **Windows key**, type `powershell`, press **Enter**. If a window was already open,
close it first: a program installed while a terminal is open is invisible to that terminal.

Check both tools are there, one line at a time:

```powershell
git --version
```

```powershell
uv --version
```

Each should print a version number. If `git` is not recognised, step 2 did not finish or you
did not reopen PowerShell.

## 4. Install from GitHub

Copy this **entire line**, paste it into PowerShell, press **Enter**:

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

**A browser window will open asking you to sign in to GitHub.** That is expected, and it
happens once. Sign in, authorise it, then go back to PowerShell.

The command downloads the tool and the right Python, and ends with:

```
Installed 2 executables: sarna-designer, sarna-designer-check
```

It can take several minutes the first time. **If you do not see that line, stop** and send the
whole output to the person who sent you the tool.

## 5. Close PowerShell, open it again, and check

```powershell
sarna-designer-check
```

The last lines name the version you now have and should end with `RESULT: PASS`. If the version
is the same one you had before, the install did not take: run step 4 again.

## 6. Start it as usual

```powershell
sarna-designer
```

Then open <http://localhost:8001/> in your browser. Leave the PowerShell window open while you
work; closing it stops the tool.

---

## From now on, updating is one line

```powershell
uv tool upgrade sarna-designer
```

That is the whole thing. No file to wait for, no filename to match. It follows the code rather
than the version number, so it picks up new work even when the version has not changed.

Run `sarna-designer-check` afterwards to see which version you are on. Quote that number if you
report a problem: every run you save records it, so a result from today can be told apart from
one produced by an older build.

Your saved runs live in a `.sarna-designer` folder, separate from the program. Updating,
reinstalling and uninstalling all leave them alone.

---

## If something goes wrong

**`repository not found`, or it asks for a password.** The invitation in step 1 has not been
accepted, or you signed in with a different GitHub account than the one invited.

**`git` is not recognised.** Step 2 did not finish, or PowerShell was not reopened after it.

**The browser sign-in never appeared and it failed instead.** Git Credential Manager was
switched off during the Git installation. Run the Git for Windows installer again and accept
the defaults.

**`sarna-designer` is not recognised.** Close PowerShell and open a new window. If that does
not help, step 4 did not finish.

**Anything else.** Run this and send the output:

```powershell
sarna-designer-check
```
