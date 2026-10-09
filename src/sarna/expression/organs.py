"""Grouping GTEx's 54 sample groups into organs, and the symbols that stand for them.

GTEx samples at a finer grain than a reader thinks in. Thirteen of its groups are nervous
system and three are artery, so "is this transcript made in the brain" cannot be answered by
looking at one row. This module does the grouping, in the data layer rather than in the page,
so it can be tested against the real tissue list and so the page cannot invent a grouping of
its own.

**A symbol is a label, not a measurement.** The emoji here stand for the organ, never for the
value: an organ appears for a transcript only when GTEx measured a median above zero in at
least one of its groups, and the number is always shown or available beside it.

**Only symbols that mean what they show.** Most organs have no Unicode emoji, and inventing one
from a lookalike object would be both undignified in a scientific tool and ambiguous. Organs
without an honest symbol carry :data:`NO_SYMBOL` and are named in words instead. The thyroid
butterfly and the kidney bean are shape conventions in wide medical use, so they are kept.

**Two GTEx groups are not organs.** ``Cells - Cultured fibroblasts`` and ``Cells -
EBV-transformed lymphocytes`` are cultured from donors, not sampled from a body site. They are
grouped as cultured cells and carry a petri dish, because filing them under an organ would
claim a tissue measurement that was not made.
"""

from collections.abc import Mapping
from dataclasses import dataclass

NO_SYMBOL = ""
"""Carried by an organ with no honest emoji. The page prints the name instead."""


@dataclass(frozen=True, slots=True)
class Organ:
    """One organ or body system, as this module groups GTEx's sample groups."""

    id: str
    name: str
    symbol: str = NO_SYMBOL
    """An emoji that genuinely depicts the organ, or :data:`NO_SYMBOL`."""

    is_body_site: bool = True
    """False for cultured cells, which were never sampled from a body site."""

    short: str = ""
    """A one-word label for a narrow chip, or empty to use :attr:`name`.

    "Oesophagus, stomach and intestine" is the right name for the group and the wrong thing to
    print in a table cell beside three others.
    """

    @property
    def chip(self) -> str:
        """Return the shortest honest label for this organ."""
        return self.short or self.name


ORGANS: tuple[Organ, ...] = (
    Organ("nervous", "Brain and nerve", "\N{BRAIN}", short="brain"),
    Organ("heart", "Heart", "\N{ANATOMICAL HEART}"),
    Organ("lung", "Lung", "\N{LUNGS}"),
    Organ("blood", "Blood", "\N{DROP OF BLOOD}"),
    Organ("muscle", "Skeletal muscle", "\N{FLEXED BICEPS}", short="muscle"),
    Organ("thyroid", "Thyroid", "\N{BUTTERFLY}"),
    Organ("kidney", "Kidney", "\N{BEANS}"),
    Organ("liver", "Liver"),
    Organ("pancreas", "Pancreas"),
    Organ("spleen", "Spleen"),
    Organ("digestive", "Oesophagus, stomach and intestine", short="digestive"),
    Organ("salivary", "Salivary gland", short="salivary"),
    Organ("adrenal", "Adrenal gland", short="adrenal"),
    Organ("pituitary", "Pituitary"),
    Organ("artery", "Artery"),
    Organ("adipose", "Adipose tissue", short="adipose"),
    Organ("skin", "Skin"),
    Organ("breast", "Breast"),
    Organ("bladder", "Bladder"),
    Organ("reproductive", "Reproductive tract", short="repro."),
    Organ("cultured", "Cultured cells", "\N{PETRI DISH}", is_body_site=False,
          short="cultured"),
)
"""Every organ a GTEx group maps to. Ordered roughly by how often a reader asks about it."""

BY_ID: Mapping[str, Organ] = {organ.id: organ for organ in ORGANS}

TISSUE_ORGAN: Mapping[str, str] = {
    "Adipose_Subcutaneous": "adipose",
    "Adipose_Visceral_Omentum": "adipose",
    "Adrenal_Gland": "adrenal",
    "Artery_Aorta": "artery",
    "Artery_Coronary": "artery",
    "Artery_Tibial": "artery",
    "Bladder": "bladder",
    "Brain_Amygdala": "nervous",
    "Brain_Anterior_cingulate_cortex_BA24": "nervous",
    "Brain_Caudate_basal_ganglia": "nervous",
    "Brain_Cerebellar_Hemisphere": "nervous",
    "Brain_Cerebellum": "nervous",
    "Brain_Cortex": "nervous",
    "Brain_Frontal_Cortex_BA9": "nervous",
    "Brain_Hippocampus": "nervous",
    "Brain_Hypothalamus": "nervous",
    "Brain_Nucleus_accumbens_basal_ganglia": "nervous",
    "Brain_Putamen_basal_ganglia": "nervous",
    "Brain_Spinal_cord_cervical_c-1": "nervous",
    "Brain_Substantia_nigra": "nervous",
    "Nerve_Tibial": "nervous",
    "Breast_Mammary_Tissue": "breast",
    "Cells_Cultured_fibroblasts": "cultured",
    "Cells_EBV-transformed_lymphocytes": "cultured",
    "Cervix_Ectocervix": "reproductive",
    "Cervix_Endocervix": "reproductive",
    "Fallopian_Tube": "reproductive",
    "Ovary": "reproductive",
    "Prostate": "reproductive",
    "Testis": "reproductive",
    "Uterus": "reproductive",
    "Vagina": "reproductive",
    "Colon_Sigmoid": "digestive",
    "Colon_Transverse": "digestive",
    "Esophagus_Gastroesophageal_Junction": "digestive",
    "Esophagus_Mucosa": "digestive",
    "Esophagus_Muscularis": "digestive",
    "Small_Intestine_Terminal_Ileum": "digestive",
    "Stomach": "digestive",
    "Heart_Atrial_Appendage": "heart",
    "Heart_Left_Ventricle": "heart",
    "Kidney_Cortex": "kidney",
    "Kidney_Medulla": "kidney",
    "Liver": "liver",
    "Lung": "lung",
    "Minor_Salivary_Gland": "salivary",
    "Muscle_Skeletal": "muscle",
    "Pancreas": "pancreas",
    "Pituitary": "pituitary",
    "Skin_Not_Sun_Exposed_Suprapubic": "skin",
    "Skin_Sun_Exposed_Lower_leg": "skin",
    "Spleen": "spleen",
    "Thyroid": "thyroid",
    "Whole_Blood": "blood",
}
"""Every ``gtex_v10`` sample group, mapped by identifier.

Written out rather than inferred from the name prefix. ``Nerve_Tibial`` belongs with the
nervous system and does not start with ``Brain``, ``Minor_Salivary_Gland`` is not a muscle
despite ending in one of the words, and a prefix rule would put ``Cells_`` under whatever it
happened to match. An unmapped group is caught by a test rather than silently dropped.
"""


def organ_for(tissue_id: str) -> Organ | None:
    """Return the organ a sample group belongs to, or ``None`` if it is not mapped."""
    organ_id = TISSUE_ORGAN.get(tissue_id)
    return BY_ID[organ_id] if organ_id is not None else None


__all__ = [
    "BY_ID",
    "NO_SYMBOL",
    "ORGANS",
    "TISSUE_ORGAN",
    "Organ",
    "organ_for",
]
