"""Author bounded atomic identities and spectral-coordinate relations."""

import sympy as sp


def definitions():
    records = {}

    def add(family, name, lhs, rhs):
        records[(f"{family}.{name}", "1")] = sp.Eq(lhs, rhs, evaluate=False)

    family = "planetlib.atomic.balance"
    nl, nu, Plu, Pul, Ql, Qu = sp.symbols("nl nu Plu Pul Ql Qu")
    add(family, "lower_production", Ql, nu*Pul - nl*Plu)
    add(family, "upper_production", Qu, nl*Plu - nu*Pul)
    S, b, eta, T, f, h, k, c, me, g, binding = sp.symbols("S b eta T f h k c me g binding")
    sigma, ni, ne, np, gamma, alpha, chi, j = sp.symbols("sigma ni ne np gamma alpha chi j")
    add(family, "saha_factor", S, g/2*(h**2/(2*sp.pi*me*k*T))**sp.Rational(3, 2)*sp.exp(binding/(k*T)))
    add(family, "milne_factor", b, S*sp.exp(-h*f/(k*T)))
    add(family, "planck_occupation", eta, 1/(sp.exp(h*f/(k*T))-1))
    add(family, "saha_population", ni, ne*np*S)
    add(family, "photoionization_integrand", gamma, 8*sp.pi*f**2/c**2*sigma*eta)
    add(family, "recombination_integrand", alpha, 8*sp.pi*f**2/c**2*sigma*b*(1+eta))
    add(family, "bound_free_extinction", chi, sigma*(ni-ne*np*b))
    add(family, "bound_free_emission", j, 2*h*f**3/c**2*sigma*ne*np*b)

    family = "planetlib.atomic.coulomb"
    n, q, branch = sp.symbols("n q branch", integer=True)
    E = sp.Symbol("E", positive=True)
    momentum = sp.Symbol("momentum", nonnegative=True)
    Aold, Bold, Anew, Bnew, G0, G1, G2, H0, H1, H2, Gstar = sp.symbols(
        "Aold Bold Anew Bnew G0 G1 G2 H0 H1 H2 Gstar")
    add(family, "energy_momentum", E, 1+n**2*momentum**2)
    add(family, "old_a", Aold, 4*n**2-4*(q+branch)**2+(q+branch)*(2*q+3)*E)
    add(family, "old_b", Bold, 4*n**2*(n**2-(q+2)**2)*(1+(q+2*branch-1)**2*momentum**2))
    add(family, "new_a", Anew, (4*n**2-4*(q+branch)**2)/E+(q+branch)*(2*q+3))
    add(family, "new_b", Bnew, 4*n**2*(n**2-(q+2)**2)*(1+(q+2*branch-1)**2*momentum**2)/E**2)
    add(family, "old_step", G0, Aold*G1-Bold*G2)
    add(family, "new_step", H0, Anew*H1-Bnew*H2)
    for shift, original, scaled in ((0, G0, H0), (1, G1, H1), (2, G2, H2)):
        add(family, f"rescale_{shift}", original, Gstar*E**(n+1-branch-q-shift)*scaled)
    Hp1, Hp2, Hm1, Hm2 = sp.symbols("Hp1 Hp2 Hm1 Hm2")
    add(family, "plus_top_seed", Hp1, sp.Integer(1))
    add(family, "plus_next_seed", Hp2, n*(2*n-1))
    add(family, "minus_top_seed", Hm1, 1/(2*n))
    add(family, "minus_next_seed", Hm2, (2*n-1)*(4/E+n-1)*Hm1)

    family = "planetlib.radiation.lte"
    nl, nu, gu, gl, nu0, T, h, k, c = sp.symbols("nl nu gu gl nu0 T h k c")
    frequency, wavenumber, Bf, Bk, phif, phik, jf, jk, chif, chik, piB = sp.symbols(
        "frequency wavenumber Bf Bk phif phik jf jk chif chik piB")
    add(family, "boltzmann_ratio", nu/nl, gu/gl*sp.exp(-h*nu0/(k*T)))
    add(family, "planck_frequency", Bf, 2*h*frequency**3/c**2/(sp.exp(h*frequency/(k*T))-1))
    add(family, "frequency_coordinate", frequency, c*wavenumber)
    add(family, "profile_conversion", phik, c*phif)
    add(family, "emission_conversion", jk, c*jf)
    add(family, "intensity_conversion", Bk, c*Bf)
    add(family, "extinction_conversion", chik, chif)
    add(family, "hemisphere_planck", piB, sp.pi*Bk)
    return records
