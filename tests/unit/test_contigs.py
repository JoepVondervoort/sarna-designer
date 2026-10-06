"""The contig tables are reference data; these tests pin them against reality.

The GRCh37 table is additionally checked against the actual AF22 bigWig headers when those
files are present, because that comparison is what identified the tracks as GRCh37 (defect
B15) in the first place.
"""

import os
from pathlib import Path

import pytest

from sarna.genome.contigs import GRCH37_CONTIG_LENGTHS, GRCH38_CONTIG_LENGTHS

AF22_DIR = Path(os.environ.get("SARNA_AF22_DIR", Path.home() / "af22_data"))
NUCLEAR_CONTIGS = [str(n) for n in range(1, 23)] + ["X", "Y"]


@pytest.fixture(scope="module")
def af22_chroms() -> dict[str, int]:
    """Return the contig table of one AF22 track, skipping if the data is not present."""
    pybigwig = pytest.importorskip("pyBigWig")
    tracks = sorted(AF22_DIR.glob("*.bw")) if AF22_DIR.is_dir() else []
    if not tracks:
        pytest.skip(f"AF22 bigWigs not available at {AF22_DIR}")
    handle = pybigwig.open(str(tracks[0]))
    try:
        return dict(handle.chroms())
    finally:
        handle.close()


class TestTables:
    def test_both_assemblies_cover_the_same_contigs(self):
        assert set(GRCH37_CONTIG_LENGTHS) == set(GRCH38_CONTIG_LENGTHS)

    def test_primary_assembly_is_complete(self):
        assert set(GRCH37_CONTIG_LENGTHS) == {*NUCLEAR_CONTIGS, "MT"}

    def test_tables_are_immutable(self):
        with pytest.raises(TypeError):
            GRCH38_CONTIG_LENGTHS["1"] = 1  # type: ignore[index]

    @pytest.mark.parametrize("contig", NUCLEAR_CONTIGS)
    def test_the_assemblies_actually_differ_on_every_nuclear_contig(self, contig: str):
        # If these were ever equal, an assembly mix-up would be undetectable on that contig.
        assert GRCH37_CONTIG_LENGTHS[contig] != GRCH38_CONTIG_LENGTHS[contig]

    def test_mitochondrion_is_the_one_contig_the_assemblies_share(self):
        assert GRCH37_CONTIG_LENGTHS["MT"] == GRCH38_CONTIG_LENGTHS["MT"] == 16_569


class TestAgainstRealAf22Tracks:
    """Pins the GRCh37 table against the published GSE108506 bigWig headers."""

    def test_tracks_use_ucsc_contig_naming(self, af22_chroms: dict[str, int]):
        assert all(name.startswith("chr") for name in af22_chroms)

    @pytest.mark.parametrize("contig", NUCLEAR_CONTIGS)
    def test_every_nuclear_contig_matches_grch37_and_not_grch38(
        self, contig: str, af22_chroms: dict[str, int]
    ):
        observed = af22_chroms[f"chr{contig}"]
        assert observed == GRCH37_CONTIG_LENGTHS[contig]
        assert observed != GRCH38_CONTIG_LENGTHS[contig]

    def test_mitochondrion_matches_neither_because_hg19_is_not_grch37_there(
        self, af22_chroms: dict[str, int]
    ):
        # UCSC hg19 chrM is the superseded NC_001807 sequence; Ensembl GRCh37 MT is the rCRS.
        assert af22_chroms["chrM"] == 16_571
        assert af22_chroms["chrM"] != GRCH37_CONTIG_LENGTHS["MT"]
