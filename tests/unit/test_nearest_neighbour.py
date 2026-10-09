"""Nearest-neighbour stacking tables.

The legacy script's table has defect B2: three stacks are mapped to the wrong parameter. Rather
than describe that in a comment, these tests make it a *computed property*, a table is
self-consistent when every stack equals its reverse complement, and the legacy table is not.
"""

import pytest
from pydantic import ValidationError

from sarna.thermo.tables import (
    LEGACY_DNA_KJ,
    RNA_XIA_1998,
    TABLES,
    NearestNeighbourTable,
    Nucleic,
)


class TestRnaTableMatchesThePublishedValues:
    """Values confirmed against the Turner/Mathews NNDB parameter file."""

    @pytest.mark.parametrize(
        ("stack", "expected"),
        [
            ("AA", -0.93),
            ("AU", -1.10),
            ("UA", -1.33),
            ("CU", -2.08),
            ("CA", -2.11),
            ("GU", -2.24),
            ("GA", -2.35),
            ("CG", -2.36),
            ("GG", -3.26),
            ("GC", -3.42),
        ],
    )
    def test_published_stack_values(self, stack: str, expected: float):
        assert RNA_XIA_1998.stack(stack) == pytest.approx(expected)

    def test_initiation_and_penalties(self):
        assert RNA_XIA_1998.initiation == pytest.approx(4.09)
        assert RNA_XIA_1998.terminal_au_penalty == pytest.approx(0.45)
        assert RNA_XIA_1998.symmetry_correction == pytest.approx(0.43)

    def test_is_rna_in_kcal_per_mole(self):
        assert RNA_XIA_1998.alphabet is Nucleic.RNA
        assert RNA_XIA_1998.unit == "kcal/mol"

    def test_cites_its_source_with_a_locator(self):
        assert "Xia" in str(RNA_XIA_1998.provenance)
        assert RNA_XIA_1998.provenance.locator


class TestEquivalentStacks:
    @pytest.mark.parametrize(
        ("stack", "equivalent"),
        [("AC", "GU"), ("AG", "CU"), ("UC", "GA"), ("UG", "CA"), ("UU", "AA"), ("CC", "GG")],
    )
    def test_reverse_complement_stacks_are_identified(self, stack: str, equivalent: str):
        assert RNA_XIA_1998.equivalent_stack(stack) == equivalent

    @pytest.mark.parametrize(
        ("stack", "equivalent"), [("AC", "GU"), ("AG", "CU"), ("UC", "GA"), ("UG", "CA")]
    )
    def test_equivalent_stacks_carry_equal_energies(self, stack: str, equivalent: str):
        assert RNA_XIA_1998.stack(stack) == RNA_XIA_1998.stack(equivalent)

    def test_a_correct_table_is_self_consistent(self):
        assert RNA_XIA_1998.is_self_consistent is True
        assert RNA_XIA_1998.inconsistent_stacks == ()


class TestLegacyTableDefectB2:
    def test_is_dna_in_kilojoules(self):
        # SantaLucia 1998 DNA values, converted to kJ/mol, used for an RNA duplex.
        assert LEGACY_DNA_KJ.alphabet is Nucleic.DNA
        assert LEGACY_DNA_KJ.unit == "kJ/mol"

    def test_is_not_self_consistent(self):
        # This *is* defect B2, asserted rather than described.
        assert LEGACY_DNA_KJ.is_self_consistent is False

    def test_names_exactly_the_three_mispaired_stacks(self):
        inconsistent = {frozenset(pair) for pair in LEGACY_DNA_KJ.inconsistent_stacks}
        assert inconsistent == {
            frozenset({"AC", "GT"}),
            frozenset({"AG", "CT"}),
            frozenset({"TC", "GA"}),
        }

    def test_the_mispairing_is_small_in_magnitude(self):
        # Honest accounting: B2 is real but each error is about 0.1 kJ/mol.
        for left, right in LEGACY_DNA_KJ.inconsistent_stacks:
            difference = abs(LEGACY_DNA_KJ.stack(left) - LEGACY_DNA_KJ.stack(right))
            assert 0 < difference < 0.2

    def test_its_provenance_points_at_the_legacy_script(self):
        assert "legacy" in str(LEGACY_DNA_KJ.provenance).lower()

    def test_carries_a_note_warning_that_it_is_defective(self):
        assert "B2" in (LEGACY_DNA_KJ.note or "")


class TestTableValidation:
    def test_every_stack_must_be_present(self):
        with pytest.raises(ValidationError, match="missing"):
            NearestNeighbourTable(
                id="broken",
                name="Broken",
                alphabet=Nucleic.RNA,
                unit="kcal/mol",
                provenance=RNA_XIA_1998.provenance,
                stacks={"AA": -0.93},
                initiation=0.0,
                terminal_au_penalty=0.0,
                symmetry_correction=0.0,
            )

    def test_a_stack_outside_the_alphabet_is_refused(self):
        stacks = dict(RNA_XIA_1998.stacks)
        stacks["TT"] = -1.0
        with pytest.raises(ValidationError, match="alphabet"):
            NearestNeighbourTable(
                id="broken",
                name="Broken",
                alphabet=Nucleic.RNA,
                unit="kcal/mol",
                provenance=RNA_XIA_1998.provenance,
                stacks=stacks,
                initiation=0.0,
                terminal_au_penalty=0.0,
                symmetry_correction=0.0,
            )

    def test_an_unknown_stack_lookup_raises(self):
        with pytest.raises(KeyError, match="XY"):
            RNA_XIA_1998.stack("XY")

    def test_tables_are_immutable(self):
        with pytest.raises(ValidationError):
            RNA_XIA_1998.initiation = 1.0


class TestRegistryOfTables:
    def test_every_table_is_reachable_by_id(self):
        assert TABLES[RNA_XIA_1998.id] is RNA_XIA_1998
        assert TABLES[LEGACY_DNA_KJ.id] is LEGACY_DNA_KJ

    def test_the_rna_table_is_the_default_choice(self):
        # The duplex is RNA, so RNA parameters are the scientifically correct default.
        assert RNA_XIA_1998.is_recommended is True
        assert LEGACY_DNA_KJ.is_recommended is False

    def test_every_table_carries_provenance_and_a_full_stack_set(self):
        for table in TABLES.values():
            assert str(table.provenance)
            assert len(table.stacks) == 16
