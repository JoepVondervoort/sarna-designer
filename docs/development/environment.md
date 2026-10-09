# Development environment (WSL)

Verified on this machine, 2026-10-06.

| Component | Version | Notes |
|---|---|---|
| OS | Linux 6.18.40.1-microsoft-standard-WSL2 | WSL2 on Windows |
| Python | 3.11.4 | miniconda at `~/miniconda3` |
| Node | v20.20.0 | via nvm at `~/.nvm` |
| git | 2.34.1 | |
| `gh` CLI | **not installed** | needed for the `to-tickets` / `pr` workflow |
| `uv` | **not installed** | recommended Python package manager |

## Filesystem placement

The repository lives at `~/sarna-designer`, i.e. on the **Linux** filesystem, deliberately
**not** under `/mnt/c`. Cross-OS filesystem access in WSL2 is slow and `inotify` file-watching
is unreliable there, which breaks test watchers and frontend dev servers. Keep it here.

The ~890 MB of AF22 bigWig files at `~/af22_data/` stay outside the repo and are gitignored.

## Caveat on Python

`python3` currently resolves to a miniconda base environment. The project will use an isolated
virtualenv (`.venv/`) so the base environment is never mutated.
