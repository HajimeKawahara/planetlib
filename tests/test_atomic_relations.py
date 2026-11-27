"""Conditional identities derived from the stored atomic/radiation equations."""

from hashlib import sha256

import pytest
import sympy as sp

from mceg import load_contract
from mceg.symbolic import to_sympy
from tools.planetlib import ROOT


def records(path, revision="1"):
    contract = load_contract(ROOT / "catalog" / path)
    return {item.id.rsplit(".", 1)[-1]: to_sympy(item.sympy)
            for item in contract.equations if item.revision == revision}


def substitute(expression, **values):
    return expression.subs({symbol: values[symbol.name]
                            for symbol in expression.free_symbols if symbol.name in values},
                           simultaneous=True)


def test_closed_two_state_nuclear_conservation():
    equations = records("atomic/rate_balance.yaml")
    assert sp.expand(equations["lower_production"].rhs + equations["upper_production"].rhs) == 0


def test_milne_factor_uses_the_same_saha_state():
    equations = records("atomic/rate_balance.yaml")
    factor = substitute(equations["milne_factor"].rhs, S=equations["saha_factor"].rhs)
    symbols = {symbol.name: symbol for symbol in factor.free_symbols}
    h, k, T, f, binding = (symbols[key] for key in ("h", "k", "T", "f", "binding"))
    translational = equations["saha_factor"].rhs / sp.exp(binding/(k*T))
    expected = translational * sp.exp(-(h*f-binding)/(k*T))
    assert sp.simplify(factor - expected) == 0


def test_saha_planck_forward_inverse_integrands_cancel():
    equations = records("atomic/rate_balance.yaml")
    photo = equations["photoionization_integrand"].rhs
    recombination = equations["recombination_integrand"].rhs
    ni = equations["saha_population"].rhs
    ne, np = sp.symbols("ne np")
    residual = ni*photo-ne*np*recombination
    residual = substitute(residual, b=equations["milne_factor"].rhs,
                          eta=equations["planck_occupation"].rhs)
    assert sp.simplify(residual) == 0


def test_bound_free_kirchhoff_identity():
    equations = records("atomic/rate_balance.yaml")
    h, f, c = sp.symbols("h f c")
    planck = 2*h*f**3/c**2 * equations["planck_occupation"].rhs
    residual = equations["bound_free_emission"].rhs-equations["bound_free_extinction"].rhs*planck
    residual = substitute(residual, ni=equations["saha_population"].rhs,
                          b=equations["milne_factor"].rhs)
    assert sp.simplify(residual) == 0


@pytest.mark.parametrize("branch", [1, 2])
def test_coulomb_explicit_coefficients_rescale(branch):
    equations = records("atomic/coulomb_recurrence.yaml")
    energy = next(symbol for symbol in equations["old_a"].rhs.free_symbols if symbol.name == "E")
    for coefficient, power in (("a", 1), ("b", 2)):
        residual = equations["new_"+coefficient].rhs-equations["old_"+coefficient].rhs/energy**power
        assert sp.simplify(substitute(residual, branch=branch)) == 0


@pytest.mark.parametrize("branch", [1, 2])
def test_coulomb_entire_step_is_equivalent(branch):
    equations = records("atomic/coulomb_recurrence.yaml")
    old = equations["old_step"].lhs-equations["old_step"].rhs
    old = substitute(old, Aold=equations["old_a"].rhs, Bold=equations["old_b"].rhs,
                     G0=equations["rescale_0"].rhs, G1=equations["rescale_1"].rhs,
                     G2=equations["rescale_2"].rhs)
    scale = equations["rescale_0"].rhs / sp.Symbol("H0")
    new = equations["new_step"].lhs-equations["new_step"].rhs
    new = substitute(new, Anew=equations["new_a"].rhs, Bnew=equations["new_b"].rhs)
    # E is declared positive; n,q,branch are integers and Gstar is positive
    # in the model scope, so the common scale can be cancelled.
    assert sp.simplify(substitute(old/scale-new, branch=branch)) == 0


def test_coulomb_seed_normalizations():
    equations = records("atomic/coulomb_recurrence.yaml")
    n = next(symbol for symbol in equations["plus_next_seed"].rhs.free_symbols if symbol.name == "n")
    E = next(symbol for symbol in equations["minus_next_seed"].rhs.free_symbols if symbol.name == "E")
    Gstar = sp.Symbol("Gstar")
    # These original seeds are independently transcribed from the inspected
    # brocklehurst_recurrence source, valid only where the indices exist.
    old_plus_top = Gstar
    old_plus_next = n*(2*n-1)*E*Gstar
    old_minus_top = E/(2*n)*Gstar
    old_minus_next = (2*n-1)*(4+(n-1)*E)*old_minus_top
    assert sp.simplify(Gstar*equations["plus_top_seed"].rhs-old_plus_top) == 0
    assert sp.simplify(Gstar*E*equations["plus_next_seed"].rhs-old_plus_next) == 0
    assert sp.simplify(Gstar*E*equations["minus_top_seed"].rhs-old_minus_top) == 0
    normalized_minus_next = substitute(equations["minus_next_seed"].rhs,
                                      Hm1=equations["minus_top_seed"].rhs)
    assert sp.simplify(Gstar*E**2*normalized_minus_next-old_minus_next) == 0


def test_lte_source_function_of_existing_line_definition():
    lte = records("radiation/lte_coordinates.yaml")
    line = records("radiation/line_coefficients.yaml", "2")
    upper = sp.Symbol("nl") * lte["boltzmann_ratio"].rhs
    planck_at_center = substitute(lte["planck_frequency"].rhs, frequency=sp.Symbol("nu0"))
    residual = substitute(line["emission"].rhs - line["extinction"].rhs*planck_at_center,
                          nu_population=upper)
    # Avoid division by extinction: zero profile or zero A gives j = chi = 0.
    assert sp.simplify(residual) == 0


def test_lte_line_dependency_is_pinned():
    contract = load_contract(ROOT / "catalog/radiation/lte_coordinates.yaml")
    evidence = next(item for item in contract.evidence if item.id.endswith(".line_definition"))
    artifact, = evidence.artifacts
    assert sha256((ROOT / artifact.source).read_bytes()).hexdigest() == artifact.sha256


def test_spectral_density_jacobians_preserve_integrals():
    equations = records("radiation/lte_coordinates.yaml")
    jacobian = sp.diff(equations["frequency_coordinate"].rhs, sp.Symbol("wavenumber"))
    for relation, original in (("profile_conversion", "phif"),
                               ("emission_conversion", "jf"),
                               ("intensity_conversion", "Bf")):
        assert sp.simplify(equations[relation].rhs/jacobian-sp.Symbol(original)) == 0


def test_extinction_has_no_spectral_jacobian_in_transfer():
    equations = records("radiation/lte_coordinates.yaml")
    emission_k = equations["emission_conversion"].rhs
    extinction_k = equations["extinction_conversion"].rhs
    intensity_k = equations["intensity_conversion"].rhs
    c, jf, chif, Bf = sp.symbols("c jf chif Bf")
    assert sp.expand(emission_k-extinction_k*intensity_k-c*(jf-chif*Bf)) == 0


def test_planck_coordinate_and_hemisphere_conversion():
    equations = records("radiation/lte_coordinates.yaml")
    wavenumber_planck = substitute(equations["intensity_conversion"].rhs,
                                  Bf=equations["planck_frequency"].rhs)
    wavenumber_planck = substitute(wavenumber_planck, frequency=equations["frequency_coordinate"].rhs)
    h, c, v, k, T = sp.symbols("h c wavenumber k T")
    expected = 2*h*c**2*v**3/(sp.exp(h*c*v/(k*T))-1)
    assert sp.simplify(wavenumber_planck-expected) == 0
    angular_cosine = sp.Symbol("angular_cosine")
    hemisphere = 2*sp.pi*sp.integrate(angular_cosine, (angular_cosine, 0, 1))*sp.Symbol("Bk")
    assert sp.simplify(equations["hemisphere_planck"].rhs-hemisphere) == 0
