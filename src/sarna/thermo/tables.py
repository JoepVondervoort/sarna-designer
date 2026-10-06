"""The nearest-neighbour stacking tables this project ships.

Three tables, all selectable, all carrying their source:

``rna-xia-1998``
    The recommended default. An saRNA duplex is RNA:RNA, so RNA parameters are the ones that
    describe it. Values confirmed against the Turner/Mathews Nearest Neighbor Database
    parameter file, and against the worked examples in the NNDB documentation.

``dna-santalucia-1998``
    DNA parameters in kcal/mol. Offered because much of the siRNA design literature computes
    asymmetry on the DNA-equivalent sequence, so this reproduces that convention.

``legacy-dna-kj``
    What the original script actually used: the SantaLucia DNA values converted to kJ/mol, with
    three stacks mis-assigned. Preserved unfixed so the ``legacy-scanner`` preset reproduces
    historical output exactly (ADR 0002). Not for new work.
"""

from collections.abc import Mapping

from sarna.rules.provenance import Provenance
from sarna.thermo.model import NearestNeighbourTable, Nucleic

KJ_PER_KCAL = 4.184
"""Thermochemical conversion: one kcal is 4.184 kJ.

Named for the direction it converts in. It was previously called ``KCAL_PER_KJ``, which read as
the reciprocal (0.239) and would have silently scaled a threshold by 17x if anyone had used it
that way. Nothing did.
"""


def _expand(unique: Mapping[str, float], alphabet: Nucleic) -> dict[str, float]:
    """Expand ten unique stacks into all sixteen by reverse-complement equivalence."""
    complement = alphabet.complement
    table: dict[str, float] = {}
    for step, energy in unique.items():
        table[step] = energy
        table[complement[step[1]] + complement[step[0]]] = energy
    return table


# --------------------------------------------------------------------------------------
# RNA — Xia et al. 1998. The recommended default.
# --------------------------------------------------------------------------------------

_XIA_1998_UNIQUE: Mapping[str, float] = {
    "AA": -0.93,
    "AU": -1.10,
    "UA": -1.33,
    "CU": -2.08,
    "CA": -2.11,
    "GU": -2.24,
    "GA": -2.35,
    "CG": -2.36,
    "GG": -3.26,
    "GC": -3.42,
}

RNA_XIA_1998 = NearestNeighbourTable(
    id="rna-xia-1998",
    name="RNA nearest-neighbour (Xia et al. 1998)",
    alphabet=Nucleic.RNA,
    unit="kcal/mol",
    provenance=Provenance.citation(
        "Xia et al. 1998",
        "Biochemistry 37(42):14719-14735, Watson-Crick nearest-neighbour parameters; "
        "values cross-checked against the Turner/Mathews NNDB parameter tables",
    ),
    stacks=_expand(_XIA_1998_UNIQUE, Nucleic.RNA),
    initiation=4.09,
    terminal_au_penalty=0.45,
    symmetry_correction=0.43,
    is_recommended=True,
    note=(
        "The correct choice for an saRNA duplex, which is RNA:RNA. Free energies at 37 degrees "
        "Celsius in 1 M NaCl."
    ),
)


# --------------------------------------------------------------------------------------
# DNA — SantaLucia 1998 unified parameters, in kcal/mol.
# --------------------------------------------------------------------------------------

_SANTALUCIA_1998_UNIQUE: Mapping[str, float] = {
    "AA": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "GT": -1.44,
    "CT": -1.28,
    "GA": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
}

DNA_SANTALUCIA_1998 = NearestNeighbourTable(
    id="dna-santalucia-1998",
    name="DNA nearest-neighbour (SantaLucia 1998)",
    alphabet=Nucleic.DNA,
    unit="kcal/mol",
    provenance=Provenance.citation(
        "SantaLucia 1998",
        "PNAS 95(4):1460-1465, unified DNA nearest-neighbour parameters",
    ),
    stacks=_expand(_SANTALUCIA_1998_UNIQUE, Nucleic.DNA),
    initiation=1.96,
    terminal_au_penalty=0.05,
    symmetry_correction=0.43,
    note=(
        "DNA parameters applied to the DNA-equivalent sequence. Offered because much of the "
        "siRNA design literature follows this convention, but the duplex that is actually "
        "loaded is RNA."
    ),
)


# --------------------------------------------------------------------------------------
# The legacy table, preserved with its defect.
# --------------------------------------------------------------------------------------

_LEGACY_VALUES: Mapping[str, float] = {
    "AA": -4.26,
    "AT": -3.67,
    "TA": -2.50,
    "CA": -6.12,
    "GT": -6.09,
    "CT": -5.40,
    "GA": -5.51,
    "CG": -9.07,
    "GC": -9.36,
    "GG": -7.66,
}

# The original script maps three stacks to the wrong parameter. Reproduced exactly: AC should
# take the GT value, AG the CT value, and TC the GA value. Correcting them here would break the
# legacy-scanner preset's reason for existing.
_LEGACY_STACKS: dict[str, float] = {
    **_LEGACY_VALUES,
    "TT": _LEGACY_VALUES["AA"],
    "AC": _LEGACY_VALUES["CA"],
    "AG": _LEGACY_VALUES["GA"],
    "TC": _LEGACY_VALUES["CT"],
    "TG": _LEGACY_VALUES["CA"],
    "CC": _LEGACY_VALUES["GG"],
}

LEGACY_DNA_KJ = NearestNeighbourTable(
    id="legacy-dna-kj",
    name="Legacy script table (DNA, kJ/mol, defective)",
    alphabet=Nucleic.DNA,
    unit="kJ/mol",
    provenance=Provenance.legacy_script(
        "sarna_promoter_scanner.py:23-32 — SantaLucia 1998 DNA values multiplied by 4.184"
    ),
    stacks=_LEGACY_STACKS,
    initiation=0.0,
    terminal_au_penalty=0.0,
    symmetry_correction=0.0,
    note=(
        "Carries defect B2: the stacks AC, AG and TC are assigned to the wrong parameters, so "
        "each differs from its reverse-complement twin by roughly 0.1 kJ/mol. The errors are "
        "small but real, and they bias every end-stability value and therefore the filter and "
        "the score. Preserved unfixed so historical designs remain reproducible; do not use it "
        "for new work. It also applies DNA parameters to an RNA duplex and omits the "
        "initiation and terminal penalties entirely."
    ),
)


TABLES: Mapping[str, NearestNeighbourTable] = {
    table.id: table
    for table in (RNA_XIA_1998, DNA_SANTALUCIA_1998, LEGACY_DNA_KJ)
}

__all__ = [
    "DNA_SANTALUCIA_1998",
    "KJ_PER_KCAL",
    "LEGACY_DNA_KJ",
    "RNA_XIA_1998",
    "TABLES",
    "NearestNeighbourTable",
    "Nucleic",
]
