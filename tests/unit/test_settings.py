"""Tests for the environment contract.

The defect these exist to prevent: ``.env.example`` documented three variables — a contact email
for Ensembl, a cache root and an offline switch — that no code read. Two were merely inert. The
third was worse than inert: the file told the user Ensembl needs a contact address so they can be
reached about abusive usage, and the ``User-Agent`` the client actually sent carried a dangling
``https://github.com/`` and no address at all. Documenting a courtesy and then not extending it is
a claim about our own behaviour, not a missing feature.

The last test here is the guard: every variable named in ``.env.example`` must be read somewhere in
the package. It fails if a variable is documented and dead, in either direction.
"""

import re
from pathlib import Path

import pytest

from sarna import settings


class TestContactEmail:
    def test_an_unset_email_is_none_rather_than_a_placeholder(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv(settings.CONTACT_EMAIL_VARIABLE, raising=False)
        assert settings.contact_email() is None

    def test_a_blank_email_is_treated_as_unset(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # An exported-but-empty variable is the normal result of an unfilled .env, and sending
        # "()" as an identity would be worse than sending nothing.
        monkeypatch.setenv(settings.CONTACT_EMAIL_VARIABLE, "   ")
        assert settings.contact_email() is None

    def test_the_example_placeholder_is_rejected(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # .env.example ships `you@example.org`. Copying the file without editing it is the
        # likeliest mistake, and mailing a placeholder to Ensembl is not identification.
        monkeypatch.setenv(settings.CONTACT_EMAIL_VARIABLE, "you@example.org")
        assert settings.contact_email() is None

    def test_a_real_address_is_returned(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv(settings.CONTACT_EMAIL_VARIABLE, "lab@example.ac.uk")
        assert settings.contact_email() == "lab@example.ac.uk"


class TestUserAgent:
    def test_it_names_the_project_and_its_version(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv(settings.CONTACT_EMAIL_VARIABLE, raising=False)
        agent = settings.user_agent()
        assert agent.startswith("sarna-designer/")

    def test_it_carries_no_dangling_url_when_no_email_is_set(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv(settings.CONTACT_EMAIL_VARIABLE, raising=False)
        assert "https://github.com/;" not in settings.user_agent()

    def test_a_configured_email_reaches_the_header(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv(settings.CONTACT_EMAIL_VARIABLE, "lab@example.ac.uk")
        assert "lab@example.ac.uk" in settings.user_agent()


class TestCacheDir:
    def test_the_default_root_is_a_dot_cache_subdirectory(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv(settings.CACHE_DIR_VARIABLE, raising=False)
        assert settings.cache_dir("ensembl") == Path(".cache") / "ensembl"

    def test_the_environment_moves_the_root(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        monkeypatch.setenv(settings.CACHE_DIR_VARIABLE, str(tmp_path))
        assert settings.cache_dir("ensembl") == tmp_path / "ensembl"

    def test_it_does_not_create_the_directory(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        # Resolving a path is not the same as committing to write there; the client creates it.
        monkeypatch.setenv(settings.CACHE_DIR_VARIABLE, str(tmp_path))
        assert not settings.cache_dir("ensembl").exists()


class TestOfflineSwitch:
    @pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
    def test_truthy_spellings_enable_it(
        self, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, value)
        assert settings.is_offline() is True

    @pytest.mark.parametrize("value", ["0", "false", "no", "off", "", "   "])
    def test_falsy_spellings_leave_it_off(
        self, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, value)
        assert settings.is_offline() is False

    def test_unset_is_off(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv(settings.OFFLINE_VARIABLE, raising=False)
        assert settings.is_offline() is False

    def test_an_unrecognised_value_is_refused_rather_than_guessed(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # "maybe" is neither; quietly choosing one would decide whether the network is touched.
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, "maybe")
        with pytest.raises(ValueError, match=settings.OFFLINE_VARIABLE):
            settings.is_offline()


class TestNoDocumentedVariableIsDead:
    """Every variable ``.env.example`` names must be read by the package, and vice versa."""

    def example_file(self) -> Path:
        path = Path(__file__).resolve().parents[2] / ".env.example"
        if not path.is_file():
            pytest.skip(".env.example is not present, so there is nothing to check against")
        return path

    def documented(self) -> set[str]:
        # A commented-out assignment still counts. The registry variables should not be set in a
        # default .env, because they would point at files nobody has written yet, so showing them
        # commented is the correct way to document them.
        text = self.example_file().read_text()
        return set(re.findall(r"^#?\s*([A-Z][A-Z0-9_]*)=", text, flags=re.MULTILINE))

    def declared_in_code(self) -> set[str]:
        source = Path(__file__).resolve().parents[2] / "src" / "sarna"
        text = "\n".join(p.read_text() for p in source.rglob("*.py"))
        return set(re.findall(r'"(SARNA_[A-Z0-9_]+)"', text))

    def test_every_documented_variable_is_read_somewhere(self) -> None:
        dead = self.documented() - self.declared_in_code()
        assert not dead, f"documented in .env.example but read nowhere: {sorted(dead)}"

    def test_every_variable_the_package_reads_is_documented(self) -> None:
        undocumented = self.declared_in_code() - self.documented()
        assert not undocumented, (
            f"read by the package but absent from .env.example: {sorted(undocumented)}"
        )
