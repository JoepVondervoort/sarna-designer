"""``sarna-designer-check``: the first thing a colleague runs, and the last line of support.

Its value is entirely in whether the advice beside a failure is the advice that fixes it. A
colleague with no command line reads one sentence and acts on it; sending them to look at their
firewall when a dependency is missing costs them an afternoon and costs whoever maintains this
an email.
"""

import socket

import pytest

from sarna.app import check


def advice(capsys: pytest.CaptureFixture[str]) -> str:
    """Return everything printed, which is what the person actually sees."""
    return capsys.readouterr().out


class TestTheServerCheckDistinguishesItsFailures:
    """A missing dependency and a blocked port are different problems.

    They were reported identically, with the firewall advice, because both arrive as an
    exception from the same call. Leaving ``[app]`` off the install command is the likeliest
    single mistake in the whole procedure, and it produced "a firewall or security tool may be
    blocking local connections".
    """

    def test_a_missing_dependency_is_not_blamed_on_the_firewall(
        self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        def missing(_router: object, _port: int) -> None:
            raise ImportError("No module named 'httpx'")

        monkeypatch.setattr("sarna.app.web.make_server", missing)
        report = check.Report(colour=False)
        check.check_server(report)
        printed = advice(capsys)
        assert report.failures, "a missing dependency must fail the check"
        # The advice may mention a firewall in order to rule it out; what it must not do is
        # send the reader to go and look at one.
        assert "may be blocking" not in printed, printed
        assert "missing dependency" in printed
        assert "[app]" in printed, "the advice must name the thing that fixes it"

    def test_a_blocked_port_still_says_firewall(
        self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        def refused(_router: object, _port: int) -> None:
            raise OSError("address already in use")

        monkeypatch.setattr("sarna.app.web.make_server", refused)
        report = check.Report(colour=False)
        check.check_server(report)
        printed = advice(capsys)
        assert report.failures
        assert "firewall" in printed.lower(), printed
        assert "127.0.0.1" in printed, "say what the tool actually listens on"

    def test_an_unexpected_failure_asks_for_the_output_rather_than_guessing(
        self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        def odd(_router: object, _port: int) -> None:
            raise RuntimeError("something nobody predicted")

        monkeypatch.setattr("sarna.app.web.make_server", odd)
        report = check.Report(colour=False)
        check.check_server(report)
        printed = advice(capsys)
        assert report.failures
        assert "send it on" in printed
        assert "something nobody predicted" in printed, "the cause must survive into the report"

    def test_a_working_port_passes_and_names_it(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        report = check.Report(colour=False)
        check.check_server(report)
        assert not report.failures
        assert "Local web interface" in advice(capsys)


class TestTheAdviceNamesSomethingTheReaderCanDo:
    def test_the_missing_httpx_advice_points_at_the_install_command(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        # Previously "uv pip install", which is not how the tool is installed and would put
        # httpx somewhere the tool cannot see.
        _ = capsys
        import inspect

        body = inspect.getsource(check.check_imports)
        assert "uv tool install" in body, "the remedy must match how the tool is installed"
        assert "uv pip install" not in body

    def test_every_failure_path_supplies_a_fix(self) -> None:
        """``fail`` takes the remedy as a required argument, so none can be omitted."""
        import inspect

        signature = inspect.signature(check.Report.fail)
        assert list(signature.parameters) == ["self", "label", "fix"]
        assert signature.parameters["fix"].default is inspect.Parameter.empty


def test_the_check_runs_end_to_end_and_reports_a_result(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # The offline path: no network, so nothing here depends on Ensembl being reachable.
    report = check.Report(colour=False)
    for step in (check.check_python, check.check_imports, check.check_engine,
                 check.check_storage, check.check_text_encoding, check.check_server):
        step(report)
    printed = advice(capsys)
    assert "OK" in printed
    assert not report.failures, printed


def test_a_free_port_is_what_the_server_check_actually_proves() -> None:
    # Independent derivation: the check claims it can open a port, so opening one here must
    # also succeed. If this fails the environment is the problem, not the tool.
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        assert probe.getsockname()[1] > 0
