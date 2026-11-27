"""Resolve pinned local definitions into ordinary MCEG 0.6 contracts."""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import sys
from typing import Literal

from pydantic import Field, TypeAdapter
import yaml

from mceg import ContractV06, load_contract, mathematical_model_digest, trace_provenance
from mceg.equation_checks import equation_checks
from mceg.integrity import file_integrity
from mceg.schema import Record, Text
from mceg.schema_v02 import Digest
from mceg.schema_v05 import MathematicalModelPin, mathematical_model_records
from mceg.validator import read_document


ROOT = Path(__file__).resolve().parents[1]
CATALOG_SUFFIXES = {".yaml", ".yml", ".json"}
COLLECTIONS = {"quantity": "quantities", "equation": "equations",
               "mathematical_model": "mathematical_models"}


class Import(Record):
    path: Text
    sha256: Digest
    model: MathematicalModelPin


class Lock(Record):
    format: Literal["planetlib.imports.v1"]
    imports: list[Import] = Field(min_length=1)


def catalog_path(root: Path, relative: str) -> Path:
    path = Path(relative)
    if (path.is_absolute() or ".." in path.parts or "\\" in relative
            or path.as_posix() != relative or not relative.startswith("catalog/")):
        raise ValueError(f"Expected a normalized catalog-relative path: {relative}")
    result = (root / path).resolve()
    if not result.is_relative_to((root / "catalog").resolve()):
        raise ValueError(f"Catalog path escapes its root: {relative}")
    return result


def pinned_contract(root: Path, pin: Import) -> ContractV06:
    path = catalog_path(root, pin.path)
    if sha256(path.read_bytes()).hexdigest() != pin.sha256:
        raise ValueError(f"Changed catalog bytes: {pin.path}; review before replacing the pin")
    contract = load_contract(path)
    if contract.schema_version != "0.6":
        raise ValueError("This prototype imports only MCEG schema 0.6")
    digest = mathematical_model_digest(contract, pin.model.id, pin.model.revision)
    if digest != pin.model.sha256:
        raise ValueError(f"Changed model definition: {pin.model.id}@{pin.model.revision}")
    return contract


def merge_record(document: dict, collection: str, record: dict) -> None:
    """Compare records using MCEG's defaults and mathematical-source priority."""
    adapter = TypeAdapter(ContractV06.model_fields[collection].annotation)
    record, = adapter.dump_python(adapter.validate_python([record]), mode="json")
    records = document.setdefault(collection, [])
    key = record["id"], record.get("revision")
    matches = [r for r in records if (r["id"], r.get("revision")) == key]
    if matches:
        matches = adapter.dump_python(adapter.validate_python(matches), mode="json")
        if len(matches) != 1 or matches[0] != record:
            raise ValueError(f"Conflicting {collection} definition: {key}")
    else:
        records.append(deepcopy(record))


def check_assessment_pins(contract: ContractV06, root: Path | None = None) -> None:
    """Use native MCEG freshness checks without renewing any judgment."""
    for assessment in contract.assessments:
        for result in equation_checks(contract, assessment):
            if result["status"] == "FAIL":
                raise ValueError(f"Stale Assessment {assessment.id}: " + "; ".join(result["details"]))
        if (root is not None and assessment.kind == "evidence_support"
                and assessment.status in {"supported", "contradicted"}):
            evidence_ids = {node["id"] for key in assessment.evidence
                            for node in trace_provenance(contract, key)["nodes"]}
            integrity = file_integrity(contract, assessment, evidence_ids, root)
            for item in integrity["selected_files"] + integrity["artifact_hashes"]:
                if item["status"] == "FAIL":
                    raise ValueError(f"Stale Assessment {assessment.id}: "
                                     f"{item['source']}: {item['message']}")


def assemble(template: Path, lock_path: Path, root: Path = ROOT) -> ContractV06:
    """Copy the selected definition closure; leave all adopter judgments local."""
    lock = Lock.model_validate(read_document(lock_path))
    document = deepcopy(read_document(template))
    for pin in lock.imports:
        source = pinned_contract(root, pin)
        records = mathematical_model_records(source, pin.model.id, pin.model.revision)
        for item in records:
            merge_record(document, COLLECTIONS[item["kind"]], item["record"])

        # Retain equation provenance, including its ancestry, without importing
        # upstream evaluations, Subjects, Contexts, Claims or Assessments.
        evidence = {e.id: e for e in source.evidence}
        pending = {e for r in records if r["kind"] == "equation"
                   for e in r["record"]["evidence"]}
        included = set()
        while pending:
            identifier = min(pending)
            pending.remove(identifier)
            if identifier in included:
                continue
            item = evidence[identifier]
            if item.subject is not None or item.context is not None:
                raise ValueError(f"Definition evidence is scoped to an upstream evaluation: {identifier}")
            merge_record(document, "evidence", item.model_dump(mode="json"))
            included.add(identifier)
            if item.provenance:
                pending.update(link.evidence for link in item.provenance.derived_from)
    contract = ContractV06.model_validate(document)
    check_assessment_pins(contract)
    return contract


def check(root: Path = ROOT) -> list[str]:
    """Check native catalog records and every local example adoption."""
    failures = []
    definitions = {}
    paths = sorted(path for path in (root / "catalog").rglob("*")
                   if path.is_file() and path.suffix in CATALOG_SUFFIXES)
    if not paths:
        failures.append("No catalog contracts found under catalog/")
    for path in paths:
        try:
            contract = load_contract(path)
            if contract.schema_version != "0.6":
                raise ValueError("The catalog uses schema 0.6")
            check_assessment_pins(contract, root)
            for collection in (*COLLECTIONS.values(), "evidence"):
                for item in getattr(contract, collection):
                    record = item.model_dump(mode="json")
                    if not item.id.startswith("planetlib."):
                        raise ValueError(f"Missing planetlib namespace: {item.id}")
                    key = collection, item.id, record.get("revision")
                    if key in definitions and definitions[key] != record:
                        raise ValueError(f"Conflicting catalog identity: {key}")
                    definitions[key] = record
        except (OSError, ValueError, yaml.YAMLError) as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    for lock in sorted((root / "examples").rglob("imports.lock.json")):
        templates = sorted(path for path in lock.parent.iterdir()
                           if path.is_file() and path != lock and path.suffix in CATALOG_SUFFIXES)
        for template in templates:
            try:
                assemble(template, lock, root)
            except (OSError, ValueError, yaml.YAMLError) as exc:
                failures.append(f"{template.relative_to(root)}: {exc}")
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("check", help="Check catalog and example pins without updating them")
    build = commands.add_parser("build", help="Materialize a standalone native MCEG contract")
    build.add_argument("template", type=Path)
    build.add_argument("--lock", required=True, type=Path)
    build.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "check":
            failures = check()
            if failures:
                print("\n".join(failures), file=sys.stderr)
                return 1
            print("Catalog and adoption pins valid (structure and identity only).")
        else:
            contract = assemble(args.template, args.lock)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(yaml.safe_dump(contract.model_dump(mode="json"),
                                                  sort_keys=False), encoding="utf-8")
            print(f"Wrote {args.output} (no support transferred).")
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
