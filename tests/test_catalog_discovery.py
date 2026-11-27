"""Native-readable contract formats cannot silently bypass catalog checks."""

import json
from pathlib import Path
import shutil

import pytest
import yaml

from mceg.validator import read_document
from tools.planetlib import ROOT, check


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(data) if path.suffix == ".json" else yaml.safe_dump(data)
    path.write_text(content, encoding="utf-8")


@pytest.mark.parametrize("suffix", [".yaml", ".yml", ".json"])
def test_invalid_catalog_format_is_checked(tmp_path, suffix):
    path = tmp_path / "catalog" / ("invalid" + suffix)
    write(path, {"schema_version": "0.6"})
    assert any(str(path.relative_to(tmp_path)) in failure for failure in check(tmp_path))


@pytest.mark.parametrize("suffix", [".yaml", ".yml", ".json"])
def test_invalid_adoption_format_is_checked(tmp_path, suffix):
    source = Path("catalog/solubility/sun_yao2026.yaml")
    (tmp_path / source.parent).mkdir(parents=True)
    shutil.copyfile(ROOT / source, tmp_path / source)
    folder = Path("examples/sun_yao")
    shutil.copytree(ROOT / folder, tmp_path / folder)
    assert check(tmp_path) == []  # The import lock itself is not an adopter.
    template = read_document(ROOT / folder / "provider.yaml")
    template["not_a_native_field"] = True
    path = tmp_path / folder / ("invalid" + suffix)
    write(path, template)
    assert any(str(path.relative_to(tmp_path)) in failure for failure in check(tmp_path))


def test_empty_catalog_is_not_reported_as_valid(tmp_path):
    (tmp_path / "catalog").mkdir()
    assert check(tmp_path) == ["No catalog contracts found under catalog/"]
