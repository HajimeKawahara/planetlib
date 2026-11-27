"""Derive spherical conservation and sonic identities from the stored relations."""

from functools import lru_cache

import sympy as sp

from mceg import load_contract
from mceg.symbolic import to_sympy
from tools.planetlib import ROOT


@lru_cache
def relations(family):
    filename = {"energy": "energy_budget", "parker": "isothermal_parker"}[family]
    contract = load_contract(ROOT / f"catalog/escape/{filename}.yaml")
    return {item.id.rsplit(".", 1)[-1]: to_sympy(item.sympy)
            for item in contract.equations if item.revision == "1"}


def symbol(expression, name):
    return next(item for item in expression.atoms(sp.Symbol) if item.name == name)


def residual(equation):
    return equation.lhs - equation.rhs


def test_continuum_energy_cancellation():
    equations = relations("energy")
    total = equations["total_balance"]
    r = symbol(total, "r")
    F, B, Q = (sp.Function(name)(r) for name in ("F", "B", "Q"))
    radiation = equations["radiation_balance"].doit()
    material = equations["material_balance"].subs(Q, equations["material_heating"].rhs)
    derivatives = {
        sp.diff(F, r): sp.solve(radiation, sp.diff(F, r))[0],
        sp.diff(B, r): sp.solve(material, sp.diff(B, r))[0],
    }
    assert sp.simplify(residual(total).doit().subs(derivatives)) == 0

    # Substituting F=4*pi*H must leave the two independent 4*pi factors:
    # radiation moment -> physical flux, and spherical area -> luminosity.
    H = sp.Function("H")(r)
    flux = equations["flux_moment"]
    luminosity = 4*sp.pi*r**2*flux.rhs
    assert sp.simplify(luminosity - 16*sp.pi**2*r**2*H) == 0


def test_cell_energy_cancellation():
    equations = relations("energy")
    total = equations["cell_total"]
    B_right = symbol(total, "B_right")
    F_right = symbol(total, "F_right")
    exchange = {
        B_right: sp.solve(equations["cell_material"], B_right)[0],
        F_right: sp.solve(equations["cell_radiation"], F_right)[0],
    }
    assert sp.simplify(residual(total).subs(exchange)) == 0


def test_adjacent_shells_telescope():
    equations = relations("energy")
    cell = residual(equations["cell_total"])
    volume = equations["shell_volume"]
    cell = cell.subs(volume.lhs, volume.rhs)
    symbols = {item.name: item for item in cell.atoms(sp.Symbol)}
    radii = sp.symbols("r0 r1 r2", positive=True)
    fluxes = sp.symbols("F0 F1 F2", real=True)
    energies = sp.symbols("B0 B1 B2", real=True)
    losses = sp.symbols("C0 C1", nonnegative=True)

    def shell(index):
        return cell.subs({
            symbols["r_left"]: radii[index], symbols["r_right"]: radii[index+1],
            symbols["F_left"]: fluxes[index], symbols["F_right"]: fluxes[index+1],
            symbols["B_left"]: energies[index], symbols["B_right"]: energies[index+1],
            symbols["C_cell"]: losses[index],
        }, simultaneous=True)

    summed = sp.expand(shell(0) + shell(1))
    expected = symbols["mdot"]*(energies[2]-energies[0])
    expected += 4*sp.pi*(radii[2]**2*fluxes[2]-radii[0]**2*fluxes[0])
    expected += 4*sp.pi*sum((radii[i+1]**3-radii[i]**3)*losses[i]/3 for i in range(2))
    assert sp.simplify(summed - expected) == 0
    assert fluxes[1] not in summed.free_symbols
    assert energies[1] not in summed.free_symbols


def test_parker_density_elimination():
    equations = relations("parker")
    reduced = equations["reduced_wind"]
    r = symbol(reduced, "r")
    rho, pressure = (sp.Function(name)(r) for name in ("rho", "P"))
    density_gradient = sp.solve(sp.diff(residual(equations["continuity"]), r),
                                sp.diff(rho, r))[0]
    pressure_gradient = sp.diff(equations["isothermal_pressure"].rhs, r)
    derived = residual(equations["momentum"]).subs(sp.diff(pressure, r), pressure_gradient)
    derived = derived.subs(sp.diff(rho, r), density_gradient)
    assert sp.simplify(derived - residual(reduced)) == 0


def test_parker_sonic_radius_and_accelerating_slope():
    equations = relations("parker")
    wind = residual(equations["reduced_wind"])
    r, a, gm = (symbol(wind, name) for name in ("r", "a", "gm"))
    sonic = equations["sonic_radius"].lhs
    slope = equations["accelerating_slope"].lhs
    v = sp.Function("v")(r)
    gravity_at_sonic = sp.solve(equations["sonic_radius"], gm)[0]
    at_sonic = {sp.diff(v, r): slope, v: a, r: sonic}
    assert sp.simplify(wind.subs(at_sonic, simultaneous=True).subs(gm, gravity_at_sonic)) == 0

    # At the regular critical point, the coefficient of v'' vanishes.
    derivative = sp.diff(wind, r)
    local = derivative.subs(at_sonic, simultaneous=True).subs(gm, gravity_at_sonic).doit()
    local = sp.simplify(local)
    roots = sp.solve(local, slope)
    assert set(roots) == {-a/sonic, a/sonic}
    selected = equations["accelerating_slope"].rhs
    assert selected.is_positive is True
    assert sp.simplify(local.subs(slope, selected)) == 0
