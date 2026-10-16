"""Embed build provenance in wheels and source archives without editing the source package."""

import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """Capture the source revision and content before packaging strips away Git history."""

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Add one generated provenance file to the artifact."""
        root = Path(self.root)
        spec = importlib.util.spec_from_file_location(
            '_sarna_reproducibility', root / 'src/sarna/reproducibility.py',
        )
        if spec is None or spec.loader is None:
            raise RuntimeError('Cannot load build identity helper')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        identity = module.software_identity(root / 'src/sarna')
        self.temporary = TemporaryDirectory(prefix='sarna-build-')
        stamp = Path(self.temporary.name) / '_build_identity.json'
        stamp.write_text(json.dumps(identity), encoding='utf-8')
        destination = ('sarna/_build_identity.json' if self.target_name == 'wheel'
                       else 'src/sarna/_build_identity.json')
        build_data.setdefault('force_include', {})[str(stamp)] = destination

    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        """Remove only the temporary file created by this build."""
        self.temporary.cleanup()
