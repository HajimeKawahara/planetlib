"""Derive thermodynamic consistency from the stored symbolic relation sets."""

import pytest
import sympy as sp
from sympy.core.function import AppliedUndef

from mceg import load_contract
from mceg.symbolic import to_sympy
from tools.planetlib import ROOT


def relations(family):
    contract = load_contract(ROOT / "catalog/thermodynamics" / f"{family}.yaml")
    return {item.id.rsplit(".", 1)[1]: to_sympy(item.sympy)
            for item in contract.equations if item.revision == "1"}


def symbols(equations):
    return {symbol.name: symbol for equation in equations.values()
            for symbol in equation.free_symbols}


def assert_zero(expression):
    assert sp.simplify(expression) == 0


@pytest.mark.parametrize("component", [1, 2])
def test_ideal_mixing_potentials_are_scalar_derivatives(component):
    eq = relations("ideal_mixing")
    variable = symbols(eq)[f"n{component}"]
    assert_zero(sp.diff(eq["gibbs"].rhs, variable) - eq[f"potential_{component}"].rhs)


def test_ideal_mixing_euler_and_hessian():
    eq = relations("ideal_mixing")
    s = symbols(eq)
    n = sp.Matrix([s["n1"], s["n2"]])
    potentials = {eq[f"potential_{i}"].lhs: eq[f"potential_{i}"].rhs for i in (1, 2)}
    assert_zero(eq["gibbs"].rhs - eq["euler"].rhs.subs(potentials))
    hessian = sp.hessian(eq["gibbs"].rhs, n)
    assert hessian == hessian.T
    for residual in hessian * n:
        assert_zero(residual)


def expanded_solute_relations():
    eq = relations("mass_fraction_solute")

    def expanded(expression):
        return expression.subs(eq["mass_fraction"].lhs, eq["mass_fraction"].rhs).subs(
            eq["equivalent_host"].lhs, eq["equivalent_host"].rhs)

    return eq, symbols(eq), expanded


@pytest.mark.parametrize("variable,potential", [
    ("ns", "solute_potential"), ("n1", "host_response_1"), ("n2", "host_response_2"),
])
def test_mass_fraction_host_and_solute_derivatives(variable, potential):
    eq, s, expanded = expanded_solute_relations()
    scalar = expanded(eq["gibbs_increment"].rhs)
    assert_zero(sp.diff(scalar, s[variable]) - expanded(eq[potential].rhs))


def test_mass_fraction_euler_and_cross_derivatives():
    eq, s, expanded = expanded_solute_relations()
    scalar = expanded(eq["gibbs_increment"].rhs)
    pairs = (("ns", "solute_potential"), ("n1", "host_response_1"), ("n2", "host_response_2"))
    assert_zero(sum(s[name] * expanded(eq[potential].rhs) for name, potential in pairs) - scalar)
    assert_zero(sp.diff(expanded(eq["host_response_1"].rhs), s["ns"])
                - sp.diff(expanded(eq["solute_potential"].rhs), s["n1"]))


def test_mass_fraction_zero_solute_is_an_energy_limit_not_a_finite_potential():
    eq, s, expanded = expanded_solute_relations()
    assert sp.limit(expanded(eq["gibbs_increment"].rhs), s["ns"], 0, dir="+") == 0
    assert sp.limit(expanded(eq["host_response_1"].rhs), s["ns"], 0, dir="+") == 0
    assert sp.limit(expanded(eq["solute_potential"].rhs), s["ns"], 0, dir="+") == -sp.oo


def test_nasa9_thermodynamic_derivatives_keep_reduced_normalizations():
    eq = relations("nasa9")
    T = symbols(eq)["T"]
    h, entropy, cp = (eq[key].rhs for key in ("enthalpy", "entropy", "heat_capacity"))
    g = eq["gibbs"].rhs.subs({eq["enthalpy"].lhs: h, eq["entropy"].lhs: entropy})
    assert_zero(sp.diff(T*h, T) - cp)
    assert_zero(T*sp.diff(entropy, T) - cp)
    assert_zero(-(g + T*sp.diff(g, T)) - entropy)


def helmholtz_relations():
    eq = relations("helmholtz")
    s = symbols(eq)
    P, entropy, cv = (eq[name].rhs for name in ("pressure", "entropy", "heat_capacity_cv"))
    substitutions = {eq["pressure"].lhs: P, eq["heat_capacity_cv"].lhs: cv}
    return eq, s, P, entropy, cv, substitutions


def test_helmholtz_maxwell_relation():
    _, s, P, entropy, _, _ = helmholtz_relations()
    assert_zero(sp.diff(entropy, s["rho"]) + sp.diff(P, s["T"])/s["rho"]**2)


def test_helmholtz_cp_follows_constant_pressure_derivative():
    eq, s, P, entropy, _, substitutions = helmholtz_relations()
    T, rho = s["T"], s["rho"]
    density_slope = -sp.diff(P, T)/sp.diff(P, rho)
    cp_from_entropy = T * (sp.diff(entropy, T) + sp.diff(entropy, rho)*density_slope)
    assert_zero(eq["heat_capacity_cp"].rhs.subs(substitutions).doit() - cp_from_entropy)


def test_helmholtz_sound_speed_follows_frozen_entropy_derivative():
    eq, s, P, entropy, _, substitutions = helmholtz_relations()
    T, rho = s["T"], s["rho"]
    temperature_slope = -sp.diff(entropy, rho)/sp.diff(entropy, T)
    sound_from_entropy = (sp.diff(P, rho) + sp.diff(P, T)*temperature_slope)/s["M"]
    assert_zero(eq["sound_speed"].rhs.subs(substitutions).doit() - sound_from_entropy)


def test_helmholtz_constant_cv_ideal_gas_limit():
    eq, s, P, _, cv, substitutions = helmholtz_relations()
    T, rho, M = s["T"], s["rho"], s["M"]
    a = next(field for field in P.atoms(AppliedUndef) if str(field.func) == "a")
    R, cv0, rho0, T0 = sp.symbols("R cv0 rho0 T0", positive=True)
    ideal = R*T*sp.log(rho/rho0) - cv0*T*sp.log(T/T0)

    def evaluate(expression):
        return expression.subs(substitutions).subs(a, ideal).doit()

    assert_zero(evaluate(P) - rho*R*T)
    assert_zero(evaluate(cv) - cv0)
    assert_zero(evaluate(eq["heat_capacity_cp"].rhs) - (cv0+R))
    assert_zero(evaluate(eq["sound_speed"].rhs) - (1+R/cv0)*R*T/M)
