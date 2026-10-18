# Development environment

The package requires **Python 3.11**. Use an isolated virtual environment for development.
Git is required to clone the repository; Node.js is used by JavaScript checks in the test suite.
The application itself has no Node.js runtime requirement.

## Set up a checkout

```bash
git clone https://github.com/JoepVondervoort/sarna-designer.git
cd sarna-designer
python3.11 -m venv .venv
.venv/bin/python -m pip install -e ".[app,dev]"
```

The repository currently requires GitHub read access. On Windows, use
`.venv\Scripts\python.exe` in place of `.venv/bin/python`, or follow the
[application installation guide](../../INSTALL.md) if you only need to run the workbench.

## Verify changes

```bash
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m mypy src tests scripts
```

GitHub Actions runs these checks on Ubuntu with Python 3.11. The default suite excludes
live-network tests. Optional AF22 track tests skip when their data are unavailable.
Browser interactions have a separate [verification script](../../scripts/verify_workbench_browser.py)
that takes a Chromium executable through `--chrome`.

## Optional genomics data

Install `.[app,dev,genomics]` to read local AF22 bigWig tracks. The genomics extra adds
pyBigWig, NumPy and pyliftover; native Windows support for pyBigWig is limited, so use Linux,
WSL or macOS for these tracks. Keep downloaded genomic data outside the checkout and set
`SARNA_AF22_DIR` to their location. The core application works without them.

For WSL development, keep the checkout on the Linux filesystem to reduce cross-filesystem
I/O. Application runs and reference indexes are stored separately from the source tree;
see [reproducibility](reproducibility.md) and [local screening validation](local-screening-validation.md).
