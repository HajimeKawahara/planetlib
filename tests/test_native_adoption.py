"""Materialized adoptions remain standalone contracts for native MCEG readers."""

from pathlib import Path
import subprocess
import sys

import pytest

from mceg import load_contract, mathematical_model_digest
from mceg.validator import read_document


ROOT = Path(__file__).resolve().parents[1]
ADOPTIONS = [
    pytest.param(lock, template, id=template.relative_to(ROOT).as_posix())
    for lock in sorted((ROOT / "examples").rglob("imports.lock.json"))
    for template in sorted(lock.parent.glob("*.yaml"))
]


@pytest.mark.parametrize("lock,template", ADOPTIONS)
def test_materialized_adoption_is_native_mceg(tmp_path, lock, template):
    output = tmp_path / "adoption.yaml"
    build = subprocess.run(
        [sys.executable, "-I", str(ROOT / "tools/planetlib.py"), "build", str(template),
         "--lock", str(lock), "--output", str(output)],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert build.returncode == 0, build.stdout + build.stderr

    # Isolated Python ignores PYTHONPATH and the repository's import path.
    # The exported file must work with installed MCEG alone.
    native = subprocess.run(
        [sys.executable, "-I", "-m", "mceg", str(output)],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert native.returncode == 0, native.stdout + native.stderr

    contract = load_contract(output)
    for imported in read_document(lock)["imports"]:
        model = imported["model"]
        assert mathematical_model_digest(contract, model["id"], model["revision"]) == model["sha256"]

    # Native record validation supplies defaults without reconstructing a
    # Planetlib-specific schema. No upstream evaluation or binding may appear.
    original = read_document(template)
    for collection in ("subjects", "contexts", "claims", "assessments", "bindings"):
        actual = getattr(contract, collection)
        expected = original.get(collection, [])
        assert len(actual) == len(expected), collection
        for record, source in zip(actual, expected):
            assert record == type(record).model_validate(source), collection
