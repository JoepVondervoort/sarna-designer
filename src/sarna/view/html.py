"""Render a :class:`~sarna.view.promoter_view.PromoterView` as one self-contained HTML page.

Self-contained on purpose. A design run is a scientific record, and a record that depends on a
server, a network fetch or a sibling data file stops being reproducible the moment any of those
moves. The page embeds its own data, so it can be attached to an email, committed beside a
manuscript, or opened in five years.

The template is `template.html` in this package, which holds the whole interface and reads its data
from one embedded JSON block. Keeping the two apart means the page can be edited as a page, rather
than as strings inside Python.
"""

import json
from importlib import resources

from sarna.view.promoter_view import PromoterView

TEMPLATE_PACKAGE = "sarna.view"
TEMPLATE_FILENAME = "template.html"
DATA_PLACEHOLDER = "__VIEW_DATA__"
TITLE_PLACEHOLDER = "__VIEW_TITLE__"
FALLBACK_TITLE = "Promoter Design Browser"


def render(
    view: PromoterView,
    motif_minimum_enrichment: float = 1.0,
    full_rules_for_rejected: bool = False,
) -> str:
    """Return the view as a complete HTML page.

    Args:
        view: The assembled promoter view.
        motif_minimum_enrichment: Motifs below this enrichment over a composition-matched shuffle
            are left off the drawable track. Every scan still appears in the page's motif table,
            so a reader sees what was looked for and not found.
        full_rules_for_rejected: Whether rejected candidates carry every rule's measurements.
            False by default here, unlike the JSON export: a 2.5 kb promoter's full run is about
            21 MB, which no browser should be asked to parse to draw a diagram. Rejected
            candidates are still present with the rules that excluded them.

    Returns:
        The page's HTML.

    Raises:
        ValueError: If the template has lost its data placeholder, which would otherwise produce
            a page that renders an empty browser with no error.
    """
    source = resources.files(TEMPLATE_PACKAGE).joinpath(TEMPLATE_FILENAME)
    template = source.read_text(encoding="utf-8")
    if DATA_PLACEHOLDER not in template:
        raise ValueError(
            f"{TEMPLATE_FILENAME} no longer contains {DATA_PLACEHOLDER}, so the run's data has "
            f"nowhere to go and the page would render empty"
        )
    title = (
        f"{view.gene.symbol} Promoter Browser" if view.gene is not None else FALLBACK_TITLE
    )
    payload = view.to_dict(
        motif_minimum_enrichment=motif_minimum_enrichment,
        full_rules_for_rejected=full_rules_for_rejected,
    )
    # separators drop the whitespace json.dumps adds by default; the page is large enough already.
    # `</script` cannot appear inside the block or the browser closes it early, and a lone U+2028
    # or U+2029 is a line terminator in JavaScript but not in JSON.
    encoded = (
        json.dumps(payload, separators=(",", ":"))
        .replace("</", "<\\/")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )
    return template.replace(TITLE_PLACEHOLDER, title).replace(DATA_PLACEHOLDER, encoded)


__all__ = ["DATA_PLACEHOLDER", "FALLBACK_TITLE", "TITLE_PLACEHOLDER", "render"]
