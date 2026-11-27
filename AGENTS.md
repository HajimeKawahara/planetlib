# Development rules

Keep repository edits in English and implementations minimal. Planetlib supplies
scientific content; MCEG supplies its record types, decoding and native checks.

Before changing record structure or tooling, read the relevant specification at
the MCEG revision pinned in `requirements-dev.txt`. The current references are:

- [Symbolic source and optional editing workflows](https://github.com/HajimeKawahara/mceg/blob/4034c9259c9c2e7d55243d69f3ab2f22c9cc7d88/docs/en/symbolic.md)
- [Quantities and mathematical models](https://github.com/HajimeKawahara/mceg/blob/4034c9259c9c2e7d55243d69f3ab2f22c9cc7d88/docs/en/mathematical_models.md)
- [Definitions, bindings and review pins](https://github.com/HajimeKawahara/mceg/blob/4034c9259c9c2e7d55243d69f3ab2f22c9cc7d88/docs/en/equations.md)
- [Claims, evidence and assessment history](https://github.com/HajimeKawahara/mceg/blob/4034c9259c9c2e7d55243d69f3ab2f22c9cc7d88/docs/en/assessment_lifecycle.md)

Distinguish MCEG requirements, optional examples and Planetlib policies in the
PR description. Name the native APIs used and explain why any added local
processing is needed. Reuse native loading, normalization, symbolic codecs,
model digests and assessment checks; do not implement a parallel schema or
formula parser. When updating MCEG, inspect its changes and update these
references and the relevant compatibility checks.

Planetlib policy is to maintain each current equation only in the catalog YAML's
`Equation.sympy` field. Do not add permanent Python formula registries or
generators beside those records. Temporary `to_sympy` / `from_sympy` editing is
allowed. Tests read catalog expressions and may construct independent reference
values, consequences or limits; tests are not a second definition source.

Preserve published QuantityDefinition, Equation and MathematicalModel records
under their existing IDs and revisions. Append a new revision or variant when
the definition changes, retaining the old records. Preserve historical LaTeX
revisions, recorded assessments and execution artifacts. Identify affected
claims and reviews; changing a lock or hash does not renew a judgment or run.
Keep source correspondence, mathematical consistency and physical applicability
separate. Do not automatically infer support from schema validity or test success,
or transfer provider assessments into an adopter implicitly.

Run the checks appropriate to the change:

```bash
.venv/bin/python tools/planetlib.py check
.venv/bin/python tools/check_history.py --base origin/main
.venv/bin/python -m pytest -q
```

Use the actual target base for a PR; CI supplies its exact commit. Review the
single maintenance location and scientific meaning explicitly. CI checks record
structure, identity and selected invariants; it cannot generally detect duplicate
mathematics in arbitrary code or establish dimensional consistency, applicability
or physical truth.
