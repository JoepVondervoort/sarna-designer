"""Tests for the TSS-relative position value type (ADR 0008).

The no-zero convention is a real source of off-by-one error, deliberately confined to this one
type. These tests pin the discontinuity at the TSS from both sides.
"""

import pytest

from sarna.genome.tss import TssPosition


class TestConstruction:
    @pytest.mark.parametrize("value", [-1000, -1, 1, 200])
    def test_accepts_any_non_zero_position(self, value: int) -> None:
        assert TssPosition(value).value == value

    def test_rejects_zero(self) -> None:
        # ADR 0008: +1 is the first transcribed base and -1 the base before it. There is no 0,
        # so accepting one would mean accepting a position that does not exist.
        with pytest.raises(ValueError, match="no position 0"):
            TssPosition(0)

    def test_is_immutable(self) -> None:
        position = TssPosition(-347)
        with pytest.raises((AttributeError, TypeError)):
            position.value = 1  # type: ignore[misc]

    def test_is_hashable_and_compares_by_value(self) -> None:
        assert TssPosition(-347) == TssPosition(-347)
        assert len({TssPosition(-347), TssPosition(-347)}) == 1


class TestOffsetConversion:
    @pytest.mark.parametrize(
        ("position", "offset"),
        [
            (-1000, -1000),
            (-2, -2),
            (-1, -1),
            (1, 0),
            (2, 1),
            (200, 199),
        ],
    )
    def test_converts_to_a_continuous_internal_offset(
        self, position: int, offset: int
    ) -> None:
        assert TssPosition(position).offset == offset

    @pytest.mark.parametrize("position", [-1000, -2, -1, 1, 2, 200])
    def test_round_trips_through_the_offset(self, position: int) -> None:
        assert TssPosition.from_offset(TssPosition(position).offset).value == position

    @pytest.mark.parametrize(
        ("offset", "position"),
        [(-1, -1), (0, 1), (1, 2)],
    )
    def test_builds_from_an_offset(self, offset: int, position: int) -> None:
        assert TssPosition.from_offset(offset) == TssPosition(position)


class TestDistance:
    def test_the_tss_gap_is_one_base_pair(self) -> None:
        # The whole point of the convention: -1 and +1 are adjacent bases, not two apart.
        assert TssPosition(-1).distance_to(TssPosition(1)) == 1

    def test_distance_upstream_is_negative(self) -> None:
        assert TssPosition(1).distance_to(TssPosition(-1)) == -1

    def test_distance_within_the_upstream_region(self) -> None:
        assert TssPosition(-500).distance_to(TssPosition(-200)) == 300

    def test_distance_within_the_transcript(self) -> None:
        assert TssPosition(1).distance_to(TssPosition(200)) == 199

    def test_distance_to_itself_is_zero(self) -> None:
        assert TssPosition(-347).distance_to(TssPosition(-347)) == 0


class TestShifting:
    def test_shifting_upstream_across_the_tss_skips_zero(self) -> None:
        assert TssPosition(1).shifted_by(-1) == TssPosition(-1)

    def test_shifting_downstream_across_the_tss_skips_zero(self) -> None:
        assert TssPosition(-1).shifted_by(1) == TssPosition(1)

    def test_shifting_by_zero_is_the_same_position(self) -> None:
        assert TssPosition(-347).shifted_by(0) == TssPosition(-347)

    def test_a_nineteen_mer_spans_eighteen_base_pairs(self) -> None:
        # A 19 nt target starting at -347 ends at -329, not -328: the span is length - 1.
        assert TssPosition(-347).shifted_by(18) == TssPosition(-329)

    def test_shifting_far_downstream(self) -> None:
        assert TssPosition(-1000).shifted_by(1000) == TssPosition(1)


class TestRendering:
    @pytest.mark.parametrize(
        ("value", "rendered"),
        [(-347, "-347"), (1, "+1"), (200, "+200")],
    )
    def test_renders_with_an_explicit_sign(self, value: int, rendered: str) -> None:
        # "+1" rather than "1", because the sign is the information: it says which side of the
        # TSS the position is on.
        assert str(TssPosition(value)) == rendered

    def test_is_upstream(self) -> None:
        assert TssPosition(-1).is_upstream is True
        assert TssPosition(1).is_upstream is False
