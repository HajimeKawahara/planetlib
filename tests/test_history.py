"""Check released definition preservation against actual Git snapshots."""

from copy import deepcopy
import json
import shutil
import subprocess
import sys

import pytest
import yaml

from mceg import load_contract
from mceg.validator import read_document
from tools.check_history import ROOT, check_history


SOURCE = ROOT / "catalog/solubility/sun_yao2026.yaml"


def git(root, *arguments):
    return subprocess.run(["git", "-C", str(root), *arguments], check=True,
                          capture_output=True, text=True).stdout.strip()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data, indent=2) if path.suffix == ".json" else yaml.safe_dump(data, sort_keys=False)
    path.write_text(text, encoding="utf-8")


def commit(root):
    git(root, "add", ".")
    git(root, "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null",
        "commit", "-qm", "Catalog snapshot")
    return git(root, "rev-parse", "HEAD")


@pytest.fixture
def repository(tmp_path, request):
    git(tmp_path, "init", "-q")
    git(tmp_path, "config", "user.name", "Planetlib test")
    git(tmp_path, "config", "user.email", "test@example.invalid")
    (tmp_path / "README.md").write_text("Test repository.\n")
    empty = commit(tmp_path)
    path = tmp_path / ("catalog/original" + getattr(request, "param", ".yaml"))
    write(path, read_document(SOURCE))
    base = commit(tmp_path)
    return tmp_path, path, base, empty


@pytest.mark.parametrize("repository", [".yaml", ".yml", ".json"], indirect=True)
@pytest.mark.parametrize("change,collection", [
    ("meaning", "quantities"), ("expression", "equations"),
    ("connection", "mathematical_models"),
])
def test_changed_published_definition_is_rejected(repository, change, collection):
    root, path, base, _ = repository
    data = read_document(path)
    if change == "meaning":
        data["quantities"][0]["description"] = "A different physical meaning."
    elif change == "expression":
        data["equations"][0]["sympy"]["expression"] = "Equality(Symbol('w'), Integer(0))"
    else:
        mapping = data["mathematical_models"][0]["relations"][0]["quantities"]
        mapping["P"], mapping["p_H2O"] = mapping["p_H2O"], mapping["P"]
    write(path, data)
    load_contract(path)  # Valid native records can still overwrite old identities.
    failures = check_history(root, base)
    assert any(f"Changed definition: {collection} " in failure for failure in failures)


@pytest.mark.parametrize("remove", ["model", "file", "catalog"])
def test_removing_published_definitions_is_rejected(repository, remove):
    root, path, base, _ = repository
    if remove == "model":
        data = read_document(path)
        data["mathematical_models"].pop()
        for claim in data["claims"]:
            claim["mathematical_models"] = [model for model in claim["mathematical_models"]
                                           if model["revision"] != "1"]
        # Unknown initial reviews contain no definition pins.
        write(path, data)
        load_contract(path)
    elif remove == "file":
        path.unlink()
    else:
        shutil.rmtree(root / "catalog")
    assert any("Removed definition:" in failure for failure in check_history(root, base))


def test_new_revisions_preserve_previous_records(repository):
    root, path, base, _ = repository
    data = read_document(path)
    equation = deepcopy(data["equations"][0])
    equation["revision"] = "3"
    equation["sympy"]["expression"] = "Equality(Symbol('w'), Integer(0))"
    data["equations"].append(equation)
    model = deepcopy(data["mathematical_models"][0])
    model["revision"] = "3"
    model["relations"][0]["equation"]["revision"] = "3"
    data["mathematical_models"].append(model)
    quantity = deepcopy(data["quantities"][0])
    quantity["revision"] = "2"
    data["quantities"].append(quantity)
    write(path, data)
    assert check_history(root, base) == []


@pytest.mark.parametrize("repository,destination", [
    (".yaml", ".json"), (".json", ".yml"), (".yml", ".yaml"),
], indirect=["repository"])
def test_native_normalization_and_file_moves_preserve_identity(repository, destination):
    root, path, base, _ = repository
    data = read_document(path)
    data["quantities"][0]["coordinates"] = []  # Native default.
    data["equations"][0]["latex"] = "Ignored display for a SymPy Equation."
    for collection in ("quantities", "equations", "mathematical_models"):
        data[collection].reverse()
    moved = root / ("catalog/reorganized/new name" + destination)
    write(moved, data)
    moved.write_text("\n" + moved.read_text())
    path.unlink()
    assert check_history(root, base) == []


def test_identical_duplicates_are_allowed_but_conflicts_fail(repository):
    root, path, base, _ = repository
    duplicate = root / "catalog/duplicate.yaml"
    data = read_document(path)
    write(duplicate, data)
    assert check_history(root, base) == []
    data["quantities"][0]["description"] = "Conflicting shared identity."
    write(duplicate, data)
    assert any("Conflicting definition" in failure for failure in check_history(root, base))


@pytest.mark.parametrize("base", ["missing-ref", "--help"])
def test_invalid_base_is_not_an_empty_snapshot(repository, base):
    root, _, _, _ = repository
    assert any("Git rev-parse failed" in failure for failure in check_history(root, base))


def test_valid_empty_base_allows_initial_catalog(repository):
    root, _, _, empty = repository
    assert check_history(root, empty) == []


@pytest.mark.parametrize("repository", [".yaml", ".yml", ".json"], indirect=True)
@pytest.mark.parametrize("side", ["current", "base"])
def test_native_validation_rejects_invalid_contracts(repository, side):
    root, path, base, _ = repository
    data = read_document(path)
    data["equations"][0]["sympy"]["expression"] = "Equality(Symbol('undeclared'), Integer(0))"
    write(path, data)
    if side == "base":
        base = commit(root)
        write(path, read_document(SOURCE))
    failures = check_history(root, base)
    assert any(f"{side}:{path.relative_to(root)}" in failure and "Undeclared" in failure
               for failure in failures)


@pytest.mark.parametrize("repository", [".json"], indirect=True)
@pytest.mark.parametrize("side", ["current", "base"])
def test_json_suffix_preserves_native_syntax_validation(repository, side):
    root, path, base, _ = repository
    # YAML syntax in a .json file must fail on either side of the comparison.
    path.write_text(yaml.safe_dump(read_document(SOURCE), sort_keys=False))
    with pytest.raises(ValueError, match="Expecting value"):
        load_contract(path)
    if side == "base":
        base = commit(root)
        write(path, read_document(SOURCE))
    assert any(f"{side}:{path.relative_to(root)}" in failure and "Expecting value" in failure
               for failure in check_history(root, base))


@pytest.mark.parametrize("side", ["current", "base"])
def test_catalog_symlinks_are_rejected(repository, side):
    root, path, base, _ = repository
    path.unlink()
    path.symlink_to(SOURCE)
    if side == "base":
        base = commit(root)
        path.unlink()
        shutil.copyfile(SOURCE, path)
    assert any("symlink" in failure for failure in check_history(root, base))


def test_cli_requires_explicit_base_and_fails_for_missing_commit():
    command = [sys.executable, str(ROOT / "tools/check_history.py")]
    missing = subprocess.run(command, capture_output=True, text=True)
    assert missing.returncode != 0 and "--base" in missing.stderr
    invalid = subprocess.run(command + ["--base", "missing-history-guard-test-ref"],
                             capture_output=True, text=True)
    assert invalid.returncode != 0 and "Git rev-parse failed" in invalid.stderr
