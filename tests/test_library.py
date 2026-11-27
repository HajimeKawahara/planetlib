"""Exercise identity and adoption boundaries, without running numerical providers."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import shutil

import pytest
import yaml

from mceg import load_contract, mathematical_model_digest
from mceg.validator import read_document
from mceg.schema_v03 import binding_digest, equation_digest
from tools.planetlib import ROOT, assemble, catalog_path, check


@pytest.fixture
def library(tmp_path):
    for name in ("catalog", "examples", "tests", "reviews"):
        shutil.copytree(ROOT / name, tmp_path / name, ignore=shutil.ignore_patterns("__pycache__"))
    return tmp_path


def adoption(root, name="provider"):
    folder = root / "examples/sun_yao"
    return assemble(folder / f"{name}.yaml", folder / "imports.lock.json", root)


def write_yaml(path, data):
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def test_catalog_and_examples():
    assert check() == []


def test_reviewed_check_code_cannot_change_silently(library):
    path = library / "tests/test_escape_relations.py"
    path.write_text(path.read_text() + "\n# Changed after the recorded review.\n")
    failures = check(library)
    assert any("tests/test_escape_relations.py: SHA-256 mismatch" in failure for failure in failures)


def test_shared_definition_preserves_separate_judgments():
    provider, consumer = adoption(ROOT), adoption(ROOT, "consumer")
    assert provider.schema_version == consumer.schema_version == "0.6"
    assert all(e.sympy is not None and e.latex is None for e in provider.equations)
    model = provider.mathematical_models[0]
    assert mathematical_model_digest(provider, model.id, model.revision) == (
        mathematical_model_digest(consumer, model.id, model.revision))
    for name, result in (("provider", provider), ("consumer", consumer)):
        template = read_document(ROOT / f"examples/sun_yao/{name}.yaml")
        for field in ("subjects", "contexts", "claims", "assessments"):
            assert [r.id for r in getattr(result, field)] == [r["id"] for r in template[field]]
        assert all(a.status == "unknown" for a in result.assessments)
    assert provider.bindings[0].subject == provider.subjects[0].id
    assert not consumer.bindings  # The coupled Gibbs scalar is not a water-output adapter.


@pytest.mark.parametrize("change", ["meaning", "connection"])
def test_changed_meaning_or_connection_invalidates_both_pins(library, change):
    path = library / "catalog/solubility/sun_yao2026.yaml"
    data = read_document(path)
    if change == "meaning":
        data["quantities"][-1]["description"] = "Water mass divided by dry mass only."
    else:
        mapping = data["mathematical_models"][0]["relations"][0]["quantities"]
        mapping["P"], mapping["p_H2O"] = mapping["p_H2O"], mapping["P"]
    write_yaml(path, data)
    with pytest.raises(ValueError, match="Changed catalog bytes"):
        adoption(library)

    lock_path = library / "examples/sun_yao/imports.lock.json"
    lock = read_document(lock_path)
    lock["imports"][0]["sha256"] = sha256(path.read_bytes()).hexdigest()
    lock_path.write_text(json.dumps(lock))
    with pytest.raises(ValueError, match="Changed model definition"):
        adoption(library)
    failures = check(library)
    assert any("provider.yaml" in failure for failure in failures)
    assert any("consumer.yaml" in failure for failure in failures)


def test_conflicting_local_definition_rejected(library):
    path = library / "examples/sun_yao/provider.yaml"
    data = read_document(path)
    quantity = adoption(library).quantities[0].model_dump(mode="json")
    quantity["description"] = "Different normalization under the same identity."
    data["quantities"] = [quantity]
    write_yaml(path, data)
    with pytest.raises(ValueError, match="Conflicting quantities"):
        adoption(library)


@pytest.mark.parametrize("collection", ["quantities", "equations", "mathematical_models"])
def test_native_equivalent_local_definition_is_accepted(library, collection):
    original = adoption(library)
    record = getattr(original, collection)[0].model_dump(mode="json")
    if collection == "quantities":
        record.pop("coordinates")  # Omitted native default means the same thing.
    elif collection == "equations":
        record["latex"] = "Ignored display; SymPy is authoritative."
    else:
        for relation in record["relations"]:
            relation.pop("kind")  # Native default is an equation relation.
    path = library / "examples/sun_yao/provider.yaml"
    data = read_document(path)
    data[collection] = [record]
    write_yaml(path, data)
    assert adoption(library).model_dump(mode="json") == original.model_dump(mode="json")


def test_duplicate_local_definitions_are_still_rejected(library):
    path = library / "examples/sun_yao/provider.yaml"
    data = read_document(path)
    quantity = adoption(library).quantities[0].model_dump(mode="json")
    data["quantities"] = [quantity, deepcopy(quantity)]
    write_yaml(path, data)
    with pytest.raises(ValueError, match="Conflicting quantities"):
        adoption(library)


def test_repeated_identical_import_is_idempotent(library):
    original = adoption(library).model_dump(mode="json")
    lock_path = library / "examples/sun_yao/imports.lock.json"
    lock = read_document(lock_path)
    lock["imports"].append(deepcopy(lock["imports"][0]))
    lock_path.write_text(json.dumps(lock))
    assert adoption(library).model_dump(mode="json") == original


@pytest.mark.parametrize("field", ["equation_pins", "mathematical_model_pins", "binding_pins"])
def test_import_pin_does_not_renew_an_assessment_pin(library, field):
    result = adoption(library)
    equation, model, binding = result.equations[0], result.mathematical_models[0], result.bindings[0]
    path = library / "examples/sun_yao/provider.yaml"
    data = read_document(path)
    assessment = data["assessments"][0]
    assessment["equation_pins"] = [{"id": equation.id, "revision": equation.revision,
                                    "sha256": equation_digest(equation)}]
    assessment["mathematical_model_pins"] = [{
        "id": model.id, "revision": model.revision,
        "sha256": mathematical_model_digest(result, model.id, model.revision)}]
    assessment["binding_pins"] = [{"id": binding.id, "sha256": binding_digest(binding)}]
    write_yaml(path, data)
    adoption(library)  # Fresh review pins are accepted, without changing the status.
    assessment[field][0]["sha256"] = "0" * 64
    write_yaml(path, data)
    with pytest.raises(ValueError, match="Stale Assessment"):
        adoption(library)


def test_scoped_evidence_cannot_be_detached_from_upstream_subject(library):
    path = library / "catalog/solubility/sun_yao2026.yaml"
    data = read_document(path)
    data["evidence"][0]["subject"] = data["subjects"][0]["id"]
    write_yaml(path, data)
    lock_path = library / "examples/sun_yao/imports.lock.json"
    lock = read_document(lock_path)
    lock["imports"][0]["sha256"] = sha256(path.read_bytes()).hexdigest()
    lock_path.write_text(json.dumps(lock))
    with pytest.raises(ValueError, match="scoped to an upstream evaluation"):
        adoption(library)


def test_hirschmann_variants_keep_distinct_concentration_identities():
    contract = load_contract(ROOT / "catalog/solubility/hirschmann2012.yaml")
    oxide, endmember = [m for m in contract.mathematical_models if m.revision == "2"]
    assert oxide.quantities["fraction"].quantity != endmember.quantities["fraction"].quantity
    assert contract.equations[0].sympy == contract.equations[1].sympy
    assert mathematical_model_digest(contract, oxide.id, oxide.revision) != (
        mathematical_model_digest(contract, endmember.id, endmember.revision))
    # Same formula text must not cause the other concentration variant to leak
    # into a selected definition closure.
    from mceg.schema_v05 import mathematical_model_records
    records = mathematical_model_records(contract, oxide.id, oxide.revision)
    assert not any(r["id"] == endmember.quantities["fraction"].quantity.id for r in records)


def test_line_coefficients_leave_population_closure_external():
    contract = load_contract(ROOT / "catalog/radiation/line_coefficients.yaml")
    model = contract.mathematical_models[0]
    for alias in ("upper", "lower", "profile"):
        assert model.quantities[alias].role == "input"
    extinction, emission = model.relations
    for symbol in ("nu_population", "phi", "frequency"):
        assert extinction.quantities[symbol] == emission.quantities[symbol]
    assert model.quantities["frequency"].role == "coordinate"


@pytest.mark.parametrize("path", ["/tmp/elsewhere.yaml", "catalog/../../outside.yaml",
                                  "catalog/./item.yaml", "examples/provider.yaml"])
def test_import_paths_stay_in_catalog(path):
    with pytest.raises(ValueError):
        catalog_path(ROOT, path)


def test_duplicate_lock_keys_are_rejected(library):
    path = library / "examples/sun_yao/imports.lock.json"
    path.write_text('{"format":"planetlib.imports.v1","format":"planetlib.imports.v1"}')
    with pytest.raises(ValueError, match="Duplicate mapping key"):
        adoption(library)
