"""Symbolic source relations for spherical energy exchange and isothermal wind."""

import sympy as sp


def definitions() -> dict[tuple[str, str], sp.Basic]:
    """Use magnitudes in the declared cgs units; leave closures to the adopter."""
    r, mdot = sp.symbols("r mdot", positive=True)
    F, H, S, C, Q, B = (sp.Function(name)(r) for name in ("F", "H", "S", "C", "Q", "B"))
    left, right, V = sp.symbols("r_left r_right V", positive=True)
    F_left, F_right, B_left, B_right, S_cell, C_cell = sp.symbols(
        "F_left F_right B_left B_right S_cell C_cell", real=True)
    energy = {
        "flux_moment": sp.Eq(F, 4*sp.pi*H, evaluate=False),
        "material_heating": sp.Eq(Q, 4*sp.pi*S-C, evaluate=False),
        "radiation_balance": sp.Eq(sp.Derivative(r**2*F, r), -4*sp.pi*r**2*S, evaluate=False),
        "material_balance": sp.Eq(mdot*sp.Derivative(B, r), 4*sp.pi*r**2*Q, evaluate=False),
        "total_balance": sp.Eq(sp.Derivative(mdot*B+4*sp.pi*r**2*F, r),
                                -4*sp.pi*r**2*C, evaluate=False),
        "shell_volume": sp.Eq(V, (right**3-left**3)/3, evaluate=False),
        "cell_radiation": sp.Eq(right**2*F_right-left**2*F_left,
                                 -4*sp.pi*V*S_cell, evaluate=False),
        "cell_material": sp.Eq(mdot*(B_right-B_left)/(4*sp.pi),
                                V*(4*sp.pi*S_cell-C_cell), evaluate=False),
        "cell_total": sp.Eq(mdot*(B_right-B_left)
                            +4*sp.pi*(right**2*F_right-left**2*F_left)
                            +4*sp.pi*V*C_cell, 0, evaluate=False),
    }

    a, gm, sonic = sp.symbols("a gm r_s", positive=True)
    slope = sp.Symbol("slope", real=True)
    v, rho, pressure = (sp.Function(name)(r) for name in ("v", "rho", "P"))
    wind = {
        "continuity": sp.Eq(mdot, 4*sp.pi*r**2*rho*v, evaluate=False),
        "isothermal_pressure": sp.Eq(pressure, a**2*rho, evaluate=False),
        "momentum": sp.Eq(v*sp.Derivative(v, r),
                           -sp.Derivative(pressure, r)/rho-gm/r**2, evaluate=False),
        "reduced_wind": sp.Eq((v-a**2/v)*sp.Derivative(v, r),
                               2*a**2/r-gm/r**2, evaluate=False),
        "sonic_radius": sp.Eq(sonic, gm/(2*a**2), evaluate=False),
        "accelerating_slope": sp.Eq(slope, a/sonic, evaluate=False),
    }
    return {**{(f"planetlib.escape.energy.{name}", "1"): expr for name, expr in energy.items()},
            **{(f"planetlib.escape.parker.{name}", "1"): expr for name, expr in wind.items()}}
