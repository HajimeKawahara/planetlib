# Planetlib

Planetlib is a shared library of scientific definitions for planetary science,
written in [MCEG](https://github.com/HajimeKawahara/mceg). MCEG supplies the record
types and checks; Planetlib supplies quantities, equations, assumptions, sources
and partial models. The relationship is analogous to mathlib and Lean at the
level of responsibilities. Inclusion does not prove a physical model or its
equivalence to numerical code.

This prototype contains shared relations extracted from local ExoFamily records
on 2026-11-27, together with symbolic consistency checks. It has no numerical
solver or runtime Python package. Numerical providers need neither Planetlib
nor MCEG at runtime.

## Initial catalog

All catalog files are ordinary MCEG schema `0.6` contracts. Each current Equation
is maintained only in its catalog record's `sympy` field as a serialized SymPy expression;
MCEG derives display LaTeX from those expressions and the declared symbol labels.
The SymPy backend is a development dependency for authoring and checking.

| File | Shared content | Deliberately open choices |
| --- | --- | --- |
| [Sun–Yao water](catalog/solubility/sun_yao2026.yaml) | Two pressure powers, ten coefficients, partial-pressure driving quantity, complete dry oxide predictors, binary H2O-equivalent mass fraction | Iron-redox convention, gas closure, other solutes, Gibbs continuation, empirical transfer |
| [Line coefficients](catalog/radiation/line_coefficients.yaml) | Signed extinction and spontaneous emission for one transition, ordinary frequency in Hz, volume populations, emissivity per steradian | LTE or non-LTE populations, line shape, broadening, transport and boundaries |
| [Hirschmann hydrogen](catalog/solubility/hirschmann2012.yaml) | Two separately identified variants of the basalt regression: reconstructed oxide concentration and ExoEOS endmember concentration | Exact original normalization, target-host validity, gas model, full Gibbs scalar |

The three initial families above use symbolic Equation and Model revision `2`.
Their historical LaTeX revisions `1` are also retained, with their original
model digests; their quantities remain at revision `1`.
Their Claims and unknown Assessments remain attached to the original definitions;
the symbolic revision has separate definition Subjects and unknown Assessments.
Each model declares its external inputs through `ModelQuantity.role`
and states unselected closures in `scope`. A graph records mathematical
incidence, including possible cycles; it does not specify execution order or
establish solvability.

For example, the two Hirschmann regressions have the same formula text but
different quantity meanings. One mole of Mg2SiO4 counts as three oxide component
moles and one MELTS endmember mole. Equal coefficients cannot justify merging
their concentration definitions. The oxide definition is the upstream review's
explicit reconstruction, with its remaining discrepancies preserved.

## Relations with symbolic consistency checks

The additional families start at revision `1`. Checks read the catalog expressions
and derive a second relation, eliminate a variable, or take a physical limit.

| Definition | Checked relationship | Scope retained |
| --- | --- | --- |
| [Energy exchange](catalog/escape/energy_budget.yaml) | Cancellation of material and radiation source terms; cell balance and telescoping | Steady spherical flow, the same source and quadrature, final coupling equal to one; escaping two-photon energy is counted once |
| [Parker wind](catalog/escape/isothermal_parker.yaml) | Eliminate density and pressure; derive sonic radius and accelerating slope | Isothermal, constant mean molecular weight, no tide or radiation force; no solver or branch-selection algorithm |
| [Atomic balance](catalog/atomic/rate_balance.yaml) | Closed two-level nuclear conservation; paired photoionization/Milne integrands and Kirchhoff identity | A common supplied cross section, Saha populations and Planck radiation; no full-atom LTE assertion or one-way two-photon inverse |
| [Coulomb recurrence](catalog/atomic/coulomb_recurrence.yaml) | Original and normalized coefficients, downward steps and channel-specific seeds | Explicit plus/minus index domains through n=20; exact algebra does not establish floating-point stability |
| [LTE and spectral coordinates](catalog/radiation/lte_coordinates.yaml) | Boltzmann populations, the existing line coefficients' LTE source function, and Hz-to-wavenumber transformations | Matched physical frequencies, line-center prefactors, per-steradian intensity; extinction receives no spectral-density Jacobian |
| [Ideal mixing](catalog/thermodynamics/ideal_mixing.yaml) | Chemical potentials from Gibbs energy, Euler identity and Hessian structure | Two positive components at fixed T/P; standards already include any appropriate pressure contribution |
| [Mass-fraction solute](catalog/thermodynamics/mass_fraction_solute.yaml) | Solute and host derivatives of the same extensive scalar | Two host components; solute standard independent of host composition; this does not supply the nonlinear Sun–Yao Gibbs continuation |
| [NASA9](catalog/thermodynamics/nasa9.yaml) | Enthalpy, entropy, Gibbs energy and heat-capacity derivative identities | One coefficient interval, positive temperature, explicit reference units; no coefficient dataset or interval stitching |
| [Molar Helmholtz responses](catalog/thermodynamics/helmholtz.yaml) | Maxwell relation, heat-capacity difference, sound speed and ideal-gas limit | Smooth single phase, fixed composition, molar density; stable response domain and no equilibrium sound speed |

The energy model keeps the integrated moment distinct from physical flux and
the spherical area factor. Its conserved total includes material excitation
and ionization energy. Collisional redistribution is internal to that energy;
it is not subtracted again as an external cooling term. The continuum and
finite-cell identities have separate scopes, and neither establishes numerical
grid convergence or atmospheric accuracy.

The [recorded symbolic run](reviews/symbolic_checks_20261127.json) passed 32 checks
for these nine families. Fourteen narrowly scoped mathematical Claims have
reviewed `supported` Assessments, with model digests and hashes of the test code
and run record. Their initial `unknown` Assessments are retained. These judgments
support the stated conditional identities; they do not establish provider
equivalence, empirical applicability or a formal proof certificate.

## Inspect and check

Use Python 3.10 or later. Development tooling pins the inspected MCEG commit
because its package version alone does not identify the required development
schema. Installation needs network access; checks and builds then work offline.

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements-dev.txt
.venv/bin/python tools/planetlib.py check
.venv/bin/python -m pytest -q
.venv/bin/python -m mceg catalog/solubility/sun_yao2026.yaml
```

With an existing local MCEG development environment that includes its
`symbolic` extra (SymPy 1.14.x), use that environment's Python instead.
The development dependency pins MCEG revision
`4034c9259c9c2e7d55243d69f3ab2f22c9cc7d88`, including the clarified authoring
guidance. Archived run records retain the MCEG revision actually used for that run.

MCEG can render the definitions without Planetlib-specific rendering code:

```bash
.venv/bin/python -m mceg model catalog/radiation/line_coefficients.yaml \
  --model planetlib.radiation.line.coefficients --revision 2 --format mermaid
.venv/bin/python -m mceg model catalog/solubility/sun_yao2026.yaml \
  --model planetlib.solubility.sun_yao2026.binary_water --revision 2 --format tex
```

The catalog is the sole maintained source of each equation. MCEG schema
`0.6` stores SymPy constructor expressions in `Equation.sympy`; it does not
require a companion Python definition module or a generation step.
For interactive editing, use `mceg.symbolic.to_sympy(equation.sympy)` to recover
the expression and `mceg.symbolic.from_sympy(expression).model_dump(mode="json")`
to serialize an edited SymPy object into the selected contract record. Keep the
resulting definition in the catalog. Tests read that same record, and MCEG
derives its displayed LaTeX. Review definition revisions and adoption pins
when changing an equation or its meaning.

Symbols represent numerical magnitudes in the units declared by each Equation:
Sun–Yao uses K and bar; Hirschmann uses GPa for melt pressure and bar for fugacity;
line coefficients use their stated cgs/Hz conventions. These unit normalizations
remain explicit assumptions, not symbolic dimensional algebra. Decimal fit
coefficients are constructed as exact rationals to preserve their published
digits; this adds no empirical precision. Migrating the mathematical source
changes definition identity, hence revision `2` for the three migrated initial
families. Newly added families start at revision `1`.

## Fixed definitions across provider and consumer

MCEG currently resolves definitions within a contract. The small
[local resolver](tools/planetlib.py) supplies an explicit authoring step for
external reuse. It reads a checked-out or vendored catalog; it does not fetch
files, choose a latest version, infer compatibility or renew review pins.

The [Sun–Yao example](examples/sun_yao/imports.lock.json) uses one shared lock
for a [provider template](examples/sun_yao/provider.yaml) and a
[consumer template](examples/sun_yao/consumer.yaml). These are proposed adoptions
stored in Planetlib; the external repositories have not been migrated. The
templates become native MCEG contracts only after materialization:

```bash
.venv/bin/python tools/planetlib.py build examples/sun_yao/provider.yaml \
  --lock examples/sun_yao/imports.lock.json --output build/provider.yaml
.venv/bin/python tools/planetlib.py build examples/sun_yao/consumer.yaml \
  --lock examples/sun_yao/imports.lock.json --output build/consumer.yaml
.venv/bin/python -m mceg build/provider.yaml build/consumer.yaml
```

The lock records the catalog-relative path, exact file SHA256, model ID,
explicit revision and MCEG model dependency digest. The digest covers the
selected model, equation definitions, quantity meanings and connection maps.
The byte hash additionally fixes provenance and the surrounding contract.
For distribution, retain the lock and vendor the referenced catalog files, or
check out an immutable Planetlib Git commit containing them. The standalone
output needs no external resolver when read by MCEG.

Only the selected model's quantities, equations and equation evidence ancestry
are copied. Evidence scoped to an upstream Subject or Context is rejected.
Required assumptions and applicability caveats therefore belong inside the
copied definitions or their evidence, rather than only in catalog-level prose.
Catalog Subjects, Contexts, Claims, Assessments, implementation bindings and
unselected models are not adopted implicitly. Conflicting definitions under
the same ID and revision fail; identical repeated imports are idempotent.
Comparison uses MCEG-normalized records, so omitted native defaults and ignored
LaTeX on a SymPy Equation do not create false conflicts.

The provider template declares an `EquationBinding` to a source-hashed
ExoGibbs function. The consumer template records the same partial model with
its own Subject, Context and separate correspondence/applicability claims.
Its Gibbs scalar does not directly return the water relation, so it needs a
separate correspondence review. Both examples retain `unknown` assessments.
Upstream source-formula support and target-host uncertainty remain traceable
through immutable evidence links, with their original Subjects and Contexts.

## Identity, revisions and reviews

Use `planetlib.<domain>.<family>.<item>` IDs and explicit string revisions.
A quantity's meaning includes its counted material, phase, denominator,
reference state and coordinate. Equal symbols or unit strings establish no
identity. Declare conversions as equations and materially different meanings
as separate IDs; do not silently alias them. The initial catalog keeps quantities
local to each family until a cross-family identity has actually been reviewed.

After adoption, preserve old definition revisions and their source snapshots.
Create a new revision for a changed equation, meaning, scope or connection;
keep historical files available through immutable Git commits or vendored copies.
Use a distinct variant ID when both meanings remain useful. New reviews belong
to the newly selected definitions and contexts. No automatic migration is provided.

To change an adoption:

1. Inspect the Git diff of the definitions, sources and connection maps.
2. Run `check`; stale pins report the affected adopter templates. Their Claims
   identify the selected models, and their Assessments identify the reviews.
3. Review the affected correspondence and applicability claims independently.
4. Explicitly choose the new revision, file hash and model digest in the lock,
   then rebuild. A replaced import pin does not renew any Assessment pin.

The resolver also invokes MCEG's native Equation, Model and Binding assessment
pin checks. A stale recorded review pin fails even when the import lock is current;
missing pins on an `unknown` assessment remain unknown. For supported or contradicted
catalog reviews, it also checks local Subject and Evidence artifact hashes. Editing
reviewed test code or its recorded result requires a new review; rerunning tests
does not renew any judgment. Provider source execution and source-file freshness
are separate checks in the provider's environment.

Native MCEG `mathematical_model_digest(contract, id, revision)` computes the
model pin; `hashlib.sha256(path.read_bytes()).hexdigest()` computes the file pin.
There is intentionally no command that blindly updates every lock.

Source-to-provider correspondence and physical applicability remain separate
from the explicitly checked mathematical identities. The catalog's required
MCEG registration metadata and Subjects identify definition
files, not numerical implementations. Null or absent self-file hashes avoid a
self-hash cycle; adopter locks fix the actual bytes. Scientific support remains
attached to a Claim, Subject and Context, rather than to library membership.

## Development guards

The [development rules](AGENTS.md) and [PR template](.github/pull_request_template.md)
separate MCEG's required record behavior, its optional examples and Planetlib's
own policies. For structural or tooling changes, identify the pinned specification
and native APIs first, then justify any glue that MCEG does not supply. In
Planetlib, each equation is maintained in its catalog record; temporary SymPy
editing operations do not become a parallel formula registry.

The guards use MCEG's loader, serializers, model digests and assessment checks.
Planetlib adds the Git comparison and adoption workflow, not a second scientific
schema or expression parser.

Catalog discovery covers `.yaml`, `.yml` and `.json` contracts, and rejects an
empty catalog. Each example lock applies to sibling contract files in these
formats; `imports.lock.json` itself is excluded. All formats use MCEG's native
loader, including its strict JSON handling.

| Guard | Enforcement |
| --- | --- |
| Native records and review pins | `python tools/planetlib.py check` uses the pinned MCEG implementation |
| Published definition history | `python tools/check_history.py --base origin/main` compares native-normalized Quantity, Equation and MathematicalModel records by collection, ID and revision |
| Standalone adoption | Tests run the build CLI, load the output through native MCEG outside the repository, and check model digests and adopter judgments |
| One maintained definition and scoped scientific claims | Development rules and PR review; arbitrary mathematical duplication or physical applicability cannot be decided by these checks |

The history guard rejects changed or removed definitions already present in the
chosen Git base. New IDs and revisions may be added while retaining old records;
formatting changes and file moves are allowed. This append-only catalog policy
belongs to Planetlib, not to the MCEG schema. It also catches definition changes
when import locks or assessment pins have been updated in the same PR. Rewriting
an old assessment or its evidence still requires review; history checking here
covers the three definition record types, not every kind of scientific record.

CI compares PRs with the exact target base commit, pushes to `main` with its
previous commit, and other branch pushes with `origin/main`. Unpublished additions
on a feature branch can therefore be revised before merging. An invalid or
unavailable base fails; a valid commit before the catalog existed allows the
initial catalog to be added. Updating
MCEG itself requires reviewing the upstream change and rerunning the guards; do
not replace native validation with a locally relaxed schema to make an upgrade pass.

## What the repository survey changed

The survey used local source files and committed reviews, not just package
names or paper titles. Catalog evidence links pin the inspected revisions.

| Repository and inspected HEAD | Findings used in the design |
| --- | --- |
| ExoGibbs `08d3d85a6680` | Sun–Yao returns binary mass fraction and distinguishes total pressure from water partial pressure. Older contract prose loosely mentions fugacity; the extracted definition follows the source and later paper review's partial-pressure interpretation. Selected source and contract files match HEAD; unrelated local changes were present. |
| ExoInventory `02540ee13cfb` | Historical provider snapshot and a later paper review are separate records. Static equation/coefficient agreement is supported in that review; empirical target-host transfer stays unknown. Consumer Gibbs continuation and live composition need their own evaluation. |
| ExoEOS `1032f58872bd` | Coefficients match the Hirschmann basalt regression, but oxide and MELTS endmember denominators differ. Follow-up supplementary-data review narrows the uncertainty without replacing old assessments. |
| exojax-escape `18064c9f743c` | Shared energy exchange, atomic inverse relations, current Coulomb rescaling and the isothermal wind are extracted with their distinct scopes. Historical small-atom comparisons and later extended-atom reviews remain upstream; their support is not transferred. |
| ExoJAX `19fe37b5509d` | No tracked MCEG contract was found. Inspected LTE strengths use wavenumber and per-absorber normalization. Their relation to Hz-based volume coefficients requires explicit density, population and Jacobian conventions. Untracked non-LTE work and unrelated local changes were excluded. |

The current ExoGibbs contract hash is
`c1a94e3e111e6e992b31c80eb6b36aad9558381bf56fa8ecdd2188b1081363b1`;
the historical ExoInventory provider snapshot hash is
`becaf73f32e226915d753186a4ac3ed0398fee2a8e37e0d7d45e007e582eab8f`.
Their difference is an evidence-path relocation, not a new solubility law.
Both retain distinct immutable references; no historical assessment is rewritten.

## Boundary of the prototype

Checks cover native structure, explicit identity, import integrity and scoped
record preservation. They do not perform dimensional algebra, interpret domain
prose as compatibility predicates, verify remote source availability, execute
providers, reproduce experiments or provide a general equivalence prover. The tests
exercise changed meanings and connections, conflicting IDs, scoped evidence,
separate adopter judgments and distinct concentration variants.
Symbolic checks also cover four frozen Sun–Yao reference samples, its zero-water
limit, the line-coefficient LTE source function and inversion sign, and the
Hirschmann regression at two pressures. These are bounded transcription checks,
not new physical validation.

Further work includes explicit provider mappings, multi-component and
study-level composition, and reviews of empirical applicability. A package
registry, solver, universal unit system and automatic support propagation are
outside this implementation.
