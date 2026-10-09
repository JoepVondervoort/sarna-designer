"""saRNA Designer, design and rank promoter-targeting small activating RNA duplexes.

The engine is deliberately free of any web, HTTP or UI dependency so that the design rules and
their regression suite can be exercised in isolation. See ADR 0005.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("sarna-designer")
except PackageNotFoundError:  # pragma: no cover - only when run from an uninstalled tree
    __version__ = "0+unknown"

__all__ = ["__version__"]
