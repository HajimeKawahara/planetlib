"""Preserve catalog definition revisions using MCEG-normalized records."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import yaml

from mceg import load_contract


ROOT = Path(__file__).resolve().parents[1]
COLLECTIONS = ("quantities", "equations", "mathematical_models")
CATALOG_SUFFIXES = {".yaml", ".yml", ".json"}


def git(root: Path, *arguments: str) -> bytes:
    result = subprocess.run(["git", "-C", str(root), *arguments], capture_output=True)
    if result.returncode:
        raise ValueError(f"Git {arguments[0]} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def collect(index: dict, path: Path, label: str) -> None:
    """Use MCEG's schema and serialization, including its source-priority rule."""
    try:
        contract = load_contract(path)
        for collection in COLLECTIONS:
            for item in getattr(contract, collection, []):
                key = collection, item.id, item.revision
                record = item.model_dump(mode="json")
                if key in index and index[key] != record:
                    raise ValueError(f"Conflicting definition: {collection} {item.id}@{item.revision}")
                index[key] = record
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise ValueError(f"{label}: {exc}") from exc


def check_history(root: Path, base: str) -> list[str]:
    """Reject changed/deleted definitions; allow additions and file reorganization."""
    root = root.resolve()
    before, after = {}, {}
    try:
        # Resolve to a commit before constructing object names. An absent ref is
        # an error, never an empty baseline or an implicitly selected branch.
        commit = git(root, "rev-parse", "--verify", "--end-of-options",
                     f"{base}^{{commit}}").decode().strip()
        entries = git(root, "ls-tree", "-rz", "--full-tree", commit, "--", "catalog")
        with TemporaryDirectory(prefix="planetlib-history-") as directory:
            for entry in entries.split(b"\0"):
                if not entry:
                    continue
                metadata, name = entry.split(b"\t", 1)
                mode, kind, object_id = metadata.split()
                label = os.fsdecode(name)
                if mode == b"120000":
                    raise ValueError(f"Base catalog contains a symlink: {label}")
                suffix = Path(label).suffix
                if suffix not in CATALOG_SUFFIXES:
                    continue
                if kind != b"blob":
                    raise ValueError(f"Base catalog entry is not a file: {label}")
                # Native loading uses the suffix to enforce strict JSON syntax.
                path = Path(directory) / ("contract" + suffix)
                path.write_bytes(git(root, "cat-file", "blob", object_id.decode()))
                collect(before, path, f"base:{label}")

        catalog = root / "catalog"
        if catalog.is_symlink():
            raise ValueError("Current catalog is a symlink")
        for path in sorted(catalog.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"Current catalog contains a symlink: {path.relative_to(root)}")
            if path.suffix in CATALOG_SUFFIXES and path.is_file():
                collect(after, path, f"current:{path.relative_to(root)}")
    except (OSError, ValueError) as exc:
        return [str(exc)]

    failures = []
    for key, record in sorted(before.items()):
        collection, identifier, revision = key
        if key not in after:
            failures.append(f"Removed definition: {collection} {identifier}@{revision}")
        elif after[key] != record:
            failures.append(f"Changed definition: {collection} {identifier}@{revision}; "
                            "retain the original and add a new revision")
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Existing Git commit/ref to preserve")
    args = parser.parse_args(argv)
    failures = check_history(ROOT, args.base)
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"Catalog definition revisions preserved against {args.base} (MCEG-normalized records).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
