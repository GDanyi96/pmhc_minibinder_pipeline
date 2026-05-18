"""Guard against absolute paths sneaking into committed fixture YAMLs.

Background: tests/fixtures/stage1/mock_target.yaml was once committed with
cleaned_pdb=/home/user/.../mock_clean.pdb because _make_fixtures.py used
str(absolute_path). Tests passed in the generator's sandbox but failed
the moment CI checked the repo out elsewhere (path didn't exist).

This test walks every YAML under tests/fixtures/ and asserts no string
value starts with "/". Cheap (sub-100ms) recurrence guard.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures"


def _iter_strings(node: Any) -> list[str]:
    out: list[str] = []
    if isinstance(node, str):
        out.append(node)
    elif isinstance(node, dict):
        for v in node.values():
            out.extend(_iter_strings(v))
    elif isinstance(node, list):
        for v in node:
            out.extend(_iter_strings(v))
    return out


def test_no_absolute_paths_in_committed_fixtures() -> None:
    yaml_files = sorted(FIXTURE_ROOT.rglob("*.yaml"))
    offenders: list[tuple[Path, str]] = []
    for path in yaml_files:
        try:
            doc = yaml.safe_load(path.read_text())
        except yaml.YAMLError:
            continue
        for value in _iter_strings(doc):
            if value.startswith("/"):
                offenders.append((path, value))
    assert not offenders, (
        "Absolute path leaked into a committed fixture YAML "
        "(see docs/known_traps.md trap #16): " + "; ".join(f"{p}: {v}" for p, v in offenders)
    )
