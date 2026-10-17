"""Use virtual time for the shared BLAST scheduler; unit transports remain network-free."""

import pytest

from sarna.offtarget import ncbi


@pytest.fixture(autouse=True)
def virtual_blast_scheduler(monkeypatch: pytest.MonkeyPatch) -> None:
    """Exercise production pacing without spending real minutes on recorded responses."""
    current = [0.0]

    def advance(seconds: float) -> None:
        current[0] += seconds

    monkeypatch.setattr(ncbi, "_SCHEDULER", ncbi.RequestScheduler(lambda: current[0], advance))
