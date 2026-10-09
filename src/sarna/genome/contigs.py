"""Primary-assembly contig lengths for GRCh37 and GRCh38.

These tables are reference data, not configuration. They exist so that an assembly can be
*asserted* rather than assumed: the GRCh37 values for ``1``, ``2`` and ``X`` below are the ones
read directly from the AF22 bigWig headers that identified those tracks as GRCh37
(see ``docs/scientific/legacy-inventory.md`` section 12).

Names follow Ensembl convention, bare, unprefixed, with the mitochondrion as ``MT``. The
UCSC ``chr`` prefix used by the bigWig files is a file-format detail handled at that adapter
and never appears in the domain model (ADR 0004).

**The mitochondrion is not interchangeable between Ensembl GRCh37 and UCSC hg19.** Checking
these tables against the AF22 headers found all 24 nuclear contigs identical and ``chrM``
differing: UCSC hg19 carries the superseded NC_001807 sequence (16,571 bp) while Ensembl
GRCh37 ``MT`` is the rCRS (16,569 bp). They are *different sequences*, so coordinates are not
transferable between them at all. The value below is Ensembl's, which is correct for this
model; the bigWig adapter rejects ``MT`` outright rather than pretending the two agree. This
costs nothing scientifically, mitochondrial DNA is not nucleosomal, so histone marks are not
meaningful there, and saRNA targets are nuclear promoters.
"""

from collections.abc import Mapping
from types import MappingProxyType
from typing import Final

GRCH37_CONTIG_LENGTHS: Final[Mapping[str, int]] = MappingProxyType(
    {
        "1": 249_250_621,
        "2": 243_199_373,
        "3": 198_022_430,
        "4": 191_154_276,
        "5": 180_915_260,
        "6": 171_115_067,
        "7": 159_138_663,
        "8": 146_364_022,
        "9": 141_213_431,
        "10": 135_534_747,
        "11": 135_006_516,
        "12": 133_851_895,
        "13": 115_169_878,
        "14": 107_349_540,
        "15": 102_531_392,
        "16": 90_354_753,
        "17": 81_195_210,
        "18": 78_077_248,
        "19": 59_128_983,
        "20": 63_025_520,
        "21": 48_129_895,
        "22": 51_304_566,
        "X": 155_270_560,
        "Y": 59_373_566,
        "MT": 16_569,
    }
)

GRCH38_CONTIG_LENGTHS: Final[Mapping[str, int]] = MappingProxyType(
    {
        "1": 248_956_422,
        "2": 242_193_529,
        "3": 198_295_559,
        "4": 190_214_555,
        "5": 181_538_259,
        "6": 170_805_979,
        "7": 159_345_973,
        "8": 145_138_636,
        "9": 138_394_717,
        "10": 133_797_422,
        "11": 135_086_622,
        "12": 133_275_309,
        "13": 114_364_328,
        "14": 107_043_718,
        "15": 101_991_189,
        "16": 90_338_345,
        "17": 83_257_441,
        "18": 80_373_285,
        "19": 58_617_616,
        "20": 64_444_167,
        "21": 46_709_983,
        "22": 50_818_468,
        "X": 156_040_895,
        "Y": 57_227_415,
        "MT": 16_569,
    }
)
