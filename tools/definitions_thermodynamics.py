"""Author small thermodynamic relation sets using numerical SI magnitudes."""

import sympy as sp


def definitions() -> dict[tuple[str, str], sp.Basic]:
    equations = {}

    def add(family, name, lhs, rhs):
        equations[f"planetlib.thermodynamics.{family}.{name}", "1"] = sp.Eq(
            lhs, rhs, evaluate=False)

    n1, n2 = sp.symbols("n1 n2", positive=True)
    u01, u02, u1, u2, G = sp.symbols("u01 u02 u1 u2 G")
    total = n1 + n2
    add("ideal_mixing", "gibbs", G,
        n1*u01 + n2*u02 + n1*sp.log(n1/total) + n2*sp.log(n2/total))
    add("ideal_mixing", "potential_1", u1, u01 + sp.log(n1/total))
    add("ideal_mixing", "potential_2", u2, u02 + sp.log(n2/total))
    add("ideal_mixing", "euler", G, n1*u1 + n2*u2)

    M1, M2, Ms, ns, S = sp.symbols("M1 M2 Ms ns S", positive=True)
    w, u0, us, du1, du2, F = sp.symbols("w u0 us du1 du2 F")
    add("mass_fraction_solute", "equivalent_host", S, (M1*n1 + M2*n2)/Ms)
    add("mass_fraction_solute", "mass_fraction", w, ns/(S + ns))
    add("mass_fraction_solute", "gibbs_increment", F,
        ns*u0 + S*sp.log(1-w) + ns*sp.log(w))
    add("mass_fraction_solute", "solute_potential", us, u0 + sp.log(w))
    add("mass_fraction_solute", "host_response_1", du1, M1/Ms*sp.log(1-w))
    add("mass_fraction_solute", "host_response_2", du2, M2/Ms*sp.log(1-w))

    T = sp.Symbol("T", positive=True)
    a1, a2, a3, a4, a5, a6, a7, b1, b2 = sp.symbols("a1 a2 a3 a4 a5 a6 a7 b1 b2")
    cpbar, sbar, hbar, gbar = sp.symbols("cpbar sbar hbar gbar")
    add("nasa9", "heat_capacity", cpbar,
        a1/T**2 + a2/T + a3 + a4*T + a5*T**2 + a6*T**3 + a7*T**4)
    add("nasa9", "entropy", sbar,
        -a1/(2*T**2) - a2/T + a3*sp.log(T) + a4*T
        + a5*T**2/2 + a6*T**3/3 + a7*T**4/4 + b2)
    add("nasa9", "enthalpy", hbar,
        -a1/T**2 + a2*sp.log(T)/T + a3 + a4*T/2
        + a5*T**2/3 + a6*T**3/4 + a7*T**4/5 + b1/T)
    add("nasa9", "gibbs", gbar, hbar-sbar)

    rho, M = sp.symbols("rho M", positive=True)
    a, P, s, cv, cp, cs2 = (sp.Function(name)(T, rho)
                            for name in ("a", "P", "s", "cv", "cp", "cs2"))
    d = sp.Derivative
    add("helmholtz", "pressure", P, rho**2*d(a, rho))
    add("helmholtz", "entropy", s, -d(a, T))
    add("helmholtz", "heat_capacity_cv", cv, -T*d(a, T, T))
    add("helmholtz", "heat_capacity_cp", cp,
        cv + T*d(P, T)**2/(rho**2*d(P, rho)))
    add("helmholtz", "sound_speed", cs2,
        (d(P, rho) + T*d(P, T)**2/(rho**2*cv))/M)
    return equations
