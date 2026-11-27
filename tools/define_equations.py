"""Author explicitly revisioned symbolic equations; --check never changes records."""

from __future__ import annotations

import argparse
from importlib import import_module
from pathlib import Path

import sympy as sp
import yaml

from mceg.symbolic import from_sympy
from mceg.validator import read_document


ROOT = Path(__file__).resolve().parents[1]


def initial_definitions() -> dict[str, sp.Equality]:
    """Use numerical magnitudes in the units declared by each Equation.

    Rational decimal coefficients preserve the published digits, not a claim
    of exact empirical knowledge. Units and applicability remain MCEG metadata.
    """
    T, P, p, na, al, w = sp.symbols("T P p_H2O x_Na2O x_Al2O3 w")
    a = [sp.Rational(value) for value in (
        "0.807805", "-0.700336", "2.949559e-6", "-74.546597", "7146.275452",
        "0.449387", "-0.360194", "-0.308259", "2.984703", "-1.693710")]
    molecular = p**a[0] * T**a[1] * sp.exp(
        a[2]*P**2/T + a[3]*sp.sqrt(P)*na/T + a[4]*na/T)
    hydroxyl = p**a[5] * T**a[6] * sp.exp((a[7] + a[8]*al + a[9]*na)*P/T)

    chi, j, c, h, A, nu0, gu, gl, nl, nu, phi = sp.symbols(
        "chi j c h A nu0 gu gl nl nu_population phi")
    X, f = sp.symbols("X f")
    regression = sp.Eq(sp.log(X/f), -sp.Rational("11.403") - sp.Rational("0.76")*P,
                       evaluate=False)
    return {
        "planetlib.solubility.sun_yao2026.water_mass_fraction":
            sp.Eq(w, (molecular + hydroxyl)/100, evaluate=False),
        "planetlib.radiation.line.extinction":
            sp.Eq(chi, c**2*A/(8*sp.pi*nu0**2)*(gu/gl*nl - nu)*phi, evaluate=False),
        "planetlib.radiation.line.emission":
            sp.Eq(j, h*nu0*A/(4*sp.pi)*nu*phi, evaluate=False),
        "planetlib.solubility.hirschmann2012.oxide_regression": regression,
        "planetlib.solubility.hirschmann2012.endmember_regression": regression,
    }


def definitions() -> dict[tuple[str, str], sp.Basic]:
    """Collect explicit revisions without selecting a newest definition implicitly."""
    result = {(identifier, "2"): expression
              for identifier, expression in initial_definitions().items()}
    for name in ("definitions_thermodynamics", "definitions_atomic", "definitions_escape"):
        module = import_module(f"{__package__}.{name}" if __package__ else name)
        additions = module.definitions()
        if result.keys() & additions.keys():
            raise ValueError(f"Repeated authoring definitions in {name}")
        result.update(additions)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    expected = {key: from_sympy(expr).model_dump() for key, expr in definitions().items()}
    seen, changed = set(), []
    for path in sorted((ROOT / "catalog").rglob("*.yaml")):
        document = read_document(path)
        modified = False
        for equation in document["equations"]:
            key = equation["id"], equation["revision"]
            if key not in expected and "sympy" not in equation:
                continue  # Retained historical LaTeX definitions are not rewritten.
            if key not in expected or key in seen:
                raise ValueError(f"Unexpected or repeated symbolic definition: {key}")
            seen.add(key)
            if equation.get("sympy") != expected[key] or "latex" in equation:
                changed.append(key)
                equation["sympy"] = expected[key]
                equation.pop("latex", None)
                modified = True
        if modified and not args.check:
            path.write_text(yaml.safe_dump(document, sort_keys=False, width=100), encoding="utf-8")
    if seen != expected.keys():
        raise ValueError(f"Missing symbolic definitions: {sorted(expected.keys() - seen)}")
    if args.check and changed:
        print("Stale symbolic definitions: " + ", ".join(f"{key}@{rev}" for key, rev in changed))
        return 1
    print("Symbolic definitions match the authoring script." if args.check
          else "Updated symbolic definitions; review revisions and pins before adoption.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
