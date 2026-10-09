"""Every environment variable this package reads, in one place.

Scattered ``os.environ`` calls are how a documented variable comes to be read nowhere.
``.env.example`` named three (a contact address for Ensembl, a cache root and an offline switch)
and none was consulted; the client meanwhile sent a ``User-Agent`` carrying a dangling
``https://github.com/`` and no address, while the example file told the reader Ensembl asks for one.
A test now asserts that the set of variables documented and the set read are the same set, in both
directions.

Three more variables belong to the registries they extend and are read beside their loaders, where
their file formats are documented: ``SARNA_ANCHORS`` (:mod:`sarna.genome.anchors`),
``SARNA_CELL_LINES`` (:mod:`sarna.regulatory.evidence`) and ``SARNA_AF22_DIR``
(:mod:`sarna.regulatory.chromatin`). They are listed here so the inventory is complete.

**Nothing here reads the environment at import time.** Each value is fetched per call, so a caller
that sets a variable after importing is not surprised, and a test can change one without reloading
a module.
"""

import os
from pathlib import Path

from sarna import __version__

CONTACT_EMAIL_VARIABLE = "SARNA_CONTACT_EMAIL"
"""An address Ensembl can use to reach whoever is running the queries."""

CACHE_DIR_VARIABLE = "SARNA_CACHE_DIR"
"""Root for cached external data. Subdirectories are per source, e.g. ``ensembl``."""

OFFLINE_VARIABLE = "SARNA_OFFLINE"
"""Set truthy to refuse every external lookup, serving only what is already cached."""

DEFAULT_CACHE_ROOT = Path(".cache")
"""Where caches go when nothing says otherwise: beside the working directory, not in ``$HOME``."""

PLACEHOLDER_EMAILS = frozenset({"you@example.org", "you@example.com"})
"""Addresses ``.env.example`` ships. Copying the file unedited is the likeliest mistake, and
sending a placeholder is not identification, so these count as unset."""

TRUTHY = frozenset({"1", "true", "yes", "on"})
FALSY = frozenset({"0", "false", "no", "off", ""})


def contact_email() -> str | None:
    """Return the configured contact address, or ``None`` if there is none to send.

    Blank and placeholder values are ``None``: an exported-but-unfilled variable is the normal
    state of a copied ``.env``, and an unusable address is worse than no address because it
    looks like a working one.
    """
    raw = os.environ.get(CONTACT_EMAIL_VARIABLE, "").strip()
    if not raw or raw.lower() in PLACEHOLDER_EMAILS:
        return None
    return raw


def user_agent() -> str:
    """Return the ``User-Agent`` to send with external requests.

    Ensembl's usage guidance asks clients to identify themselves so they can be contacted about
    request patterns rather than simply blocked. Setting :data:`CONTACT_EMAIL_VARIABLE` extends
    that courtesy; leaving it unset still identifies the tool and its purpose, which is as much
    as can honestly be claimed without an address.
    """
    identity = f"sarna-designer/{__version__}"
    email = contact_email()
    if email is None:
        return f"{identity} (promoter-targeting saRNA design)"
    return f"{identity} (promoter-targeting saRNA design; {email})"


def cache_dir(source: str) -> Path:
    """Return the cache directory for one external source, without creating it.

    Args:
        source: A short name for the source, used as the subdirectory, e.g. ``"ensembl"``.

    Returns:
        The path. Resolving it is not a commitment to write there, whoever caches creates it,
        so asking where the cache would go never leaves an empty directory behind.
    """
    root = os.environ.get(CACHE_DIR_VARIABLE, "").strip()
    return (Path(root) if root else DEFAULT_CACHE_ROOT) / source


def is_offline() -> bool:
    """Return whether external lookups are refused.

    Raises:
        ValueError: If the variable is set to something that is neither truthy nor falsy.
            Guessing would silently decide whether the network is touched, and the two
            readings differ in whether a run is reproducible.
    """
    raw = os.environ.get(OFFLINE_VARIABLE, "").strip().lower()
    if raw in TRUTHY:
        return True
    if raw in FALSY:
        return False
    raise ValueError(
        f"{OFFLINE_VARIABLE}={raw!r} is neither on nor off; use one of "
        f"{sorted(TRUTHY)} or {sorted(FALSY)}"
    )


__all__ = [
    "CACHE_DIR_VARIABLE",
    "CONTACT_EMAIL_VARIABLE",
    "DEFAULT_CACHE_ROOT",
    "OFFLINE_VARIABLE",
    "PLACEHOLDER_EMAILS",
    "cache_dir",
    "contact_email",
    "is_offline",
    "user_agent",
]
