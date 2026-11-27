"""Check stored symbolic definitions against source samples and physical limits."""

import math

import pytest

from mceg import load_contract
from mceg.symbolic import to_sympy
from tools.planetlib import ROOT


def equation(path, identifier):
    contract = load_contract(ROOT / "catalog" / path)
    record = next(item for item in contract.equations
                  if item.id == identifier and item.revision == "2")
    return to_sympy(record.sympy)


def evaluate(expression, values):
    """Evaluate stored pointwise expressions in their declared numerical units."""
    replacements = {symbol: values[symbol.name] for symbol in expression.free_symbols}
    return float(expression.evalf(30, subs=replacements))


# Frozen scalar reference values, not laboratory measurements:
# https://github.com/HajimeKawahara/exogibbs/blob/08d3d85a6680eb81c4057bf584530d269de88253/tests/data/validation/mceg/20261021_volatile_sources/comparison.json
# The original review uses a relative tolerance of 2e-12.
@pytest.mark.parametrize("p,T,P,Na,Al,expected", [
    (1.0, 1373.0, 1000.0, 0.03, 0.1, 0.0007796489575224205),
    (100.0, 1673.0, 3000.0, 0.05, 0.15, 0.008569654899817519),
    (500.0, 1873.0, 10000.0, 0.01, 0.08, 0.0158568460565143),
    (2000.0, 2173.15, 20000.0, 0.002, 0.025, 0.03879584553140427),
])
def test_sun_yao_matches_frozen_source_samples(p, T, P, Na, Al, expected):
    expression = equation("solubility/sun_yao2026.yaml",
                          "planetlib.solubility.sun_yao2026.water_mass_fraction")
    actual = evaluate(expression.rhs, dict(p_H2O=p, T=T, P=P, x_Na2O=Na, x_Al2O3=Al))
    assert actual == pytest.approx(expected, rel=2e-12)


def test_sun_yao_zero_water_pressure_extension():
    expression = equation("solubility/sun_yao2026.yaml",
                          "planetlib.solubility.sun_yao2026.water_mass_fraction")
    assert evaluate(expression.rhs, dict(p_H2O=0, T=1673, P=3000,
                                         x_Na2O=0.05, x_Al2O3=0.15)) == 0


def line_expressions():
    return tuple(equation("radiation/line_coefficients.yaml", "planetlib.radiation.line." + name)
                 for name in ("extinction", "emission"))


def line_values():
    return dict(c=2.99792458e10, h=6.62607015e-27, A=6.265e8,
                nu0=2.466e15, frequency=2.466e15, gu=8, gl=2,
                nl=1e10, phi=1e-10)


def test_line_lte_source_function_matches_planck_at_line_center():
    extinction, emission = line_expressions()
    values = line_values()
    temperature, boltzmann = 10000.0, 1.380649e-16  # K and erg/K
    exponent = values["h"] * values["nu0"] / (boltzmann * temperature)
    values["nu_population"] = values["nl"] * values["gu"] / values["gl"] * math.exp(-exponent)
    source_function = evaluate(emission.rhs, values) / evaluate(extinction.rhs, values)
    planck = 2 * values["h"] * values["nu0"] ** 3 / values["c"] ** 2 / math.expm1(exponent)
    assert source_function == pytest.approx(planck, rel=2e-13)


def test_line_population_inversion_retains_negative_extinction():
    extinction, emission = line_expressions()
    values = dict(line_values(), nu_population=1e12)
    assert evaluate(extinction.rhs, values) < 0
    assert evaluate(emission.rhs, values) > 0


@pytest.mark.parametrize("variant", ["oxide", "endmember"])
@pytest.mark.parametrize("pressure,expected_log", [(1.0, -12.163), (3.0, -13.683)])
def test_hirschmann_basalt_log_regression(variant, pressure, expected_log):
    expression = equation("solubility/hirschmann2012.yaml",
                          f"planetlib.solubility.hirschmann2012.{variant}_regression")
    # The same coefficient check applies separately to the two concentration
    # definitions; it does not equate their oxide and endmember denominators.
    values = dict(P=pressure, f=1000.0, X=1000.0 * math.exp(expected_log))
    assert evaluate(expression.rhs, values) == pytest.approx(expected_log, abs=1e-12)
    assert evaluate(expression.lhs, values) == pytest.approx(expected_log, abs=1e-12)
