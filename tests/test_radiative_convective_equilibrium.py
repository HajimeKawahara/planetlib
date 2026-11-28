"""Check consequences of the stored RCE definitions without running a provider."""

import pytest
import sympy as sp

from mceg import load_contract, to_sympy
from tools.planetlib import ROOT


CONTRACT = load_contract(ROOT / "catalog/atmosphere/radiative_convective_equilibrium.yaml")
EQUATIONS = {item.id.rsplit(".", 1)[-1]: to_sympy(item.sympy) for item in CONTRACT.equations}
SYMBOLS = {symbol.name: symbol for equation in EQUATIONS.values() for symbol in equation.atoms(sp.Symbol)}


def substitute(expression, **values):
    return expression.subs({SYMBOLS[name]: value for name, value in values.items()})


def test_incident_beam_is_subtracted_once_from_the_flux_target():
    net = substitute(EQUATIONS["net_radiative_flux"].rhs, F_up=700, F_down=100, F_direct=500)
    assert net == 100
    balance = substitute(EQUATIONS["flux_balance"], F_rad=net, F_conv=200, F_int=300)
    assert balance is sp.true
    # Adding absorbed irradiation to the internal target would double-count it.
    assert substitute(EQUATIONS["flux_balance"], F_rad=net, F_conv=200, F_int=800) is sp.false


def test_distinct_top_boundary_recovers_radiative_balance():
    model = next(item for item in CONTRACT.mathematical_models if item.id.endswith(".balance"))
    top = next(item for item in model.relations if item.id.endswith(".top_balance"))
    assert top.quantities["F_rad"] == "F_rad_top"
    assert top.quantities["F_conv"] == "F_conv_top"
    mapping = {SYMBOLS[name]: SYMBOLS[alias] for name, alias in top.quantities.items()}
    balance = EQUATIONS["flux_balance"].xreplace(mapping)
    boundary = EQUATIONS["top_convection"]
    solved = sp.solve([balance, boundary], (SYMBOLS["F_rad_top"], SYMBOLS["F_conv_top"]))
    residual = EQUATIONS["top_residual"].rhs.subs(solved)
    assert sp.simplify(residual) == 0
    assert solved[SYMBOLS["F_rad_top"]] == SYMBOLS["F_int"]


@pytest.mark.parametrize("flux,excess,admitted", [
    (0, -2, True), (3, 0, True), (0, 0, True),
    (3, -2, False), (-1, 0, False), (0, 2, False),
])
def test_complementarity_keeps_both_signs_and_degenerate_case(flux, excess, admitted):
    assert bool(substitute(EQUATIONS["complementarity"], F_conv=flux, g=excess)) is admitted


@pytest.mark.parametrize("mask", [0, 1], ids=["radiative", "convective"])
def test_fixed_mask_root_and_recovery_obey_balance(mask):
    residual = substitute(EQUATIONS["fixed_mask_residual"].rhs, m=mask)
    convective = substitute(EQUATIONS["convective_flux_recovery"].rhs, m=mask)
    variable = SYMBOLS["F_rad"] if mask == 0 else SYMBOLS["g"]
    root, = sp.solve(residual, variable)
    balance = EQUATIONS["flux_balance"].lhs-EQUATIONS["flux_balance"].rhs
    assert sp.simplify(balance.subs(SYMBOLS["F_conv"], convective).subs(variable, root)) == 0
    if mask == 0:
        assert convective == 0
    else:
        assert root == 0


@pytest.mark.parametrize("mask,flux,excess", [(0, 10, 1), (1, 11, 0)])
def test_zero_reduced_residual_does_not_imply_admissible_equilibrium(mask, flux, excess):
    values = dict(m=mask, F_rad=flux, F_int=10, g=excess, flux_scale=3, stability_scale=2)
    assert substitute(EQUATIONS["fixed_mask_residual"].rhs, **values) == 0
    convective = substitute(EQUATIONS["convective_flux_recovery"].rhs, m=mask, F_rad=flux, F_int=10)
    assert substitute(EQUATIONS["complementarity"], F_conv=convective, g=excess) is sp.false


def test_gradient_power_law_and_pressure_unit_invariance():
    rhs = EQUATIONS["gradient_excess"].rhs
    # T increases by two when P increases by sixteen: the exact slope is 1/4.
    state = dict(T_upper=1000, T_lower=2000, P_upper=1, P_lower=16, nabla_neutral=sp.Rational(1, 4))
    assert sp.simplify(substitute(rhs, **state)) == 0
    converted = dict(state, P_upper=10**6, P_lower=16*10**6)
    assert sp.simplify(substitute(rhs, **converted)) == 0
    assert sp.simplify(substitute(rhs, **dict(state, T_lower=4000))) == sp.Rational(1, 4)


def test_entropy_offset_and_positive_normalization_preserve_the_criterion():
    rhs = EQUATIONS["entropy_excess"].rhs
    offset = sp.Symbol("offset", real=True)
    shifted = substitute(rhs, s_upper=SYMBOLS["s_upper"]+offset, s_lower=SYMBOLS["s_lower"]+offset)
    assert sp.simplify(shifted-rhs) == 0
    values = dict(s_upper=100, s_lower=120, s_ref=10, P_upper=1, P_lower=sp.E)
    assert sp.simplify(substitute(rhs, **values)) == 2
    assert sp.simplify(substitute(rhs, **dict(values, s_ref=100))) == sp.Rational(1, 5)
    assert sp.simplify(substitute(rhs, **dict(values, s_lower=80))) == -2


def test_constant_composition_ideal_gas_links_but_does_not_identify_the_variants():
    cp, gas_constant, T_upper, T_lower, P_upper, P_lower, s_ref = sp.symbols(
        "cp gas_constant Tu Td Pu Pd sref", positive=True)
    # Independent ideal-gas entropy consequence, only for constant cp and composition.
    entropy_difference = cp*sp.log(T_lower/T_upper)-gas_constant*sp.log(P_lower/P_upper)
    common = dict(T_upper=T_upper, T_lower=T_lower, P_upper=P_upper, P_lower=P_lower)
    gradient = substitute(EQUATIONS["gradient_excess"].rhs,
                          **common, nabla_neutral=gas_constant/cp)
    entropy = substitute(EQUATIONS["entropy_excess"].rhs,
                         P_upper=P_upper, P_lower=P_lower, s_ref=s_ref,
                         s_upper=0, s_lower=entropy_difference)
    assert sp.simplify(entropy-cp/s_ref*gradient) == 0
    extra_entropy = sp.Symbol("composition_entropy", real=True)
    changed = substitute(EQUATIONS["entropy_excess"].rhs,
                         P_upper=P_upper, P_lower=P_lower, s_ref=s_ref,
                         s_upper=0, s_lower=entropy_difference+extra_entropy)
    assert sp.simplify(changed-cp/s_ref*gradient) != 0


def tangent_component(matrix, response, parameter_partial, row):
    equation = substitute(EQUATIONS["implicit_tangent"], N=matrix.rows-1, i=row, a=0).doit()
    values = {sp.Function("A")(row, j): matrix[row, j] for j in range(matrix.cols)}
    values.update({sp.Function("D")(j, 0): response[j] for j in range(matrix.cols)})
    values[sp.Function("B")(row, 0)] = parameter_partial[row]
    return sp.simplify((equation.lhs-equation.rhs).subs(values))


def test_coupled_implicit_response_matches_independent_exact_root_derivative():
    parameter = sp.Symbol("parameter", real=True)
    u0, u1 = sp.symbols("u0 u1", real=True)
    state = sp.Matrix([u0, u1])
    # Synthetic coupled two-node residual with a known nonsingular smooth branch.
    residual = sp.Matrix([2*u0+u1-parameter, u0+3*u1-parameter**2])
    root = sp.solve(residual, (u0, u1))
    response = sp.Matrix([sp.diff(root[u0], parameter), sp.diff(root[u1], parameter)])
    matrix = residual.jacobian(state)
    partial = residual.diff(parameter)
    assert matrix.det() != 0
    assert all(tangent_component(matrix, response, partial, row) == 0 for row in range(2))
    assert any(tangent_component(matrix, -response, partial, row) != 0 for row in range(2))


def test_singular_jacobian_does_not_determine_a_unique_tangent():
    matrix = sp.Matrix([[1, 1], [2, 2]])
    partial = sp.Matrix([-1, -2])
    arbitrary = sp.Symbol("arbitrary", real=True)
    response = sp.Matrix([arbitrary, 1-arbitrary])
    assert matrix.det() == 0
    assert all(tangent_component(matrix, response, partial, row) == 0 for row in range(2))
    assert response.diff(arbitrary) != sp.zeros(2, 1)
