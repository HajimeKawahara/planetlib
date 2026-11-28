"""Conditional consequences of the catalog expressions, without providers."""

from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

from mceg import load_contract, to_sympy


ROOT = Path(__file__).resolve().parents[1]
PATHS = {
    "lines": "catalog/spectroscopy/lte_lines.yaml",
    "columns": "catalog/atmosphere/columns.yaml",
    "emission": "catalog/radiation/pure_absorption.yaml",
    "transit": "catalog/radiation/spherical_transmission.yaml",
    "observation": "catalog/spectroscopy/observation.yaml",
}
EXPRESSIONS = {
    family: {equation.id.rsplit(".", 1)[-1]: to_sympy(equation.sympy)
             for equation in load_contract(ROOT / path).equations}
    for family, path in PATHS.items()
}


def equation(family, name):
    return EXPRESSIONS[family][name]


def symbols(expression):
    return {symbol.name: symbol for symbol in expression.atoms(sp.Symbol)}


def replace_values(expression, **values):
    return expression.subs({symbol: values[name]
                            for name, symbol in symbols(expression).items()
                            if name in values}, simultaneous=True)


def test_lte_reference_temperature_and_composable_temperature_changes():
    strength = equation("lines", "strength_temperature").rhs
    names = symbols(strength)
    assert sp.simplify(strength.subs({names["T"]: names["T_ref"], names["Q_ratio"]: 1})
                       - names["S_ref"]) == 0
    t0, t1, t2, q0, q1, q2 = sp.symbols("t0 t1 t2 q0 q1 q2", positive=True)
    first = replace_values(strength, T=t1, T_ref=t0, Q_ratio=q1/q0)
    composed = replace_values(strength, T=t2, T_ref=t1, Q_ratio=q2/q1, S_ref=first)
    direct = replace_values(strength, T=t2, T_ref=t0, Q_ratio=q2/q0)
    assert sp.simplify(composed-direct) == 0


@pytest.mark.parametrize("name,width", [("gaussian_profile", "sigma_D"), ("lorentz_profile", "gamma_L")])
def test_profile_normalization_and_evenness(name, width):
    profile = equation("lines", name).rhs
    x = symbols(profile)["x"]
    positive_width = sp.Symbol("positive_width", positive=True)
    profile = replace_values(profile, **{width: positive_width})
    assert sp.integrate(profile, (x, -sp.oo, sp.oo)) == 1
    assert sp.simplify(profile.subs(x, -x)-profile) == 0


@pytest.mark.parametrize("offset", [-2.0, 0.0, 3.0])
def test_folded_voigt_against_independent_complex_error_function(offset):
    # The independent reference evaluates the complex-error-function identity;
    # the catalog definition is the convolution of its normalized densities.
    relation = equation("lines", "voigt_profile")
    x, y = (symbols(relation)[name] for name in ("x", "y"))
    gaussian = equation("lines", "gaussian_profile").rhs
    lorentz = equation("lines", "lorentz_profile").rhs
    g = sp.Function("G")
    l = sp.Function("L")
    integrand = relation.rhs.function.subs({g(y): gaussian.subs(x, y),
                                           g(-y): gaussian.subs(x, -y),
                                           l(x-y): lorentz.subs(x, x-y),
                                           l(x+y): lorentz.subs(x, x+y)})
    integrand = replace_values(integrand, x=offset, sigma_D=0.7, gamma_L=0.3)
    with mp.workdps(35):
        actual = mp.quad(sp.lambdify(y, integrand, "mpmath"), [0, 1, 4, mp.inf])
        z = (mp.mpf(offset)+0.3j)/(mp.sqrt(2)*mp.mpf("0.7"))
        expected = mp.re(mp.exp(-z*z)*mp.erfc(-1j*z))/(mp.sqrt(2*mp.pi)*mp.mpf("0.7"))
        assert float(actual) == pytest.approx(float(expected), rel=3e-14)


def test_finite_line_sum_uses_independent_centers_and_strengths():
    relation = equation("lines", "line_sum")
    w, index = (symbols(relation)[name] for name in ("w", "i"))
    expanded = replace_values(relation.rhs, N_lines=2).doit()
    first = sp.Function("phi_line")(1, w-3)
    second = sp.Function("phi_line")(2, w-7)
    expanded = expanded.subs({sp.Function("S_line")(1): 2,
                             sp.Function("S_line")(2): 5,
                             sp.Function("w_line")(1): 3,
                             sp.Function("w_line")(2): 7})
    assert sp.expand(expanded-2*first-5*second) == 0
    assert index not in expanded.free_symbols


def test_voigt_convolution_preserves_profile_area_and_line_strength():
    relation = equation("lines", "voigt_profile")
    x, y = (symbols(relation)[name] for name in ("x", "y"))
    gaussian = replace_values(equation("lines", "gaussian_profile").rhs, sigma_D=sp.Rational(7, 10))
    lorentz = replace_values(equation("lines", "lorentz_profile").rhs, gamma_L=sp.Rational(3, 10))
    profiles = {"G": gaussian, "L": lorentz}
    integrand = relation.rhs.function.xreplace({
        call: profiles[call.func.__name__].subs(x, call.args[0])
        for call in relation.rhs.function.atoms(sp.Function)
    })
    # Nonnegative densities allow reversing the integrations. This checks the
    # folded convolution's area independently of the complex-error identity.
    inner_integral = 0
    for term in sp.Add.make_args(integrand):
        independent, dependent = term.as_independent(x)
        inner_integral += independent * sp.integrate(dependent, (x, -sp.oo, sp.oo))
    area = sp.integrate(inner_integral, (y, 0, sp.oo))
    assert area == 1
    cross_section = equation("lines", "line_cross_section").rhs
    integrated_cross_section = cross_section.subs(sp.Function("phi")(x), area)
    assert integrated_cross_section == symbols(cross_section)["S"]


def test_mass_and_volume_columns_agree_only_after_explicit_conversion():
    mass = equation("columns", "pressure_column_mass").rhs
    volume = equation("columns", "pressure_column_volume").rhs
    conversion = equation("columns", "mixing_conversion")
    assert sp.simplify(mass.subs(conversion.lhs, conversion.rhs)-volume) == 0
    values = dict(dP=2e6, x_abs=0.2, m_abs=4.0, m_mean=2.0, g=10.0)
    converted = float(replace_values(mass, q=0.4, **values))
    assert converted == float(replace_values(volume, **values))
    assert converted != float(replace_values(mass, q=0.2, **values))


def test_cia_pair_convention_and_zero_limits():
    expression = equation("columns", "cia_depth").rhs
    for name in ("n1", "n2", "dz"):
        assert replace_values(expression, **{name: 0}) == 0
    n = sp.Symbol("n", positive=True)
    identical = replace_values(expression, n1=n, n2=n)
    assert sp.simplify(identical / n**2) == symbols(expression)["k_cia"] * symbols(expression)["dz"]


def test_inverse_square_isothermal_hydrostatic_consequence():
    hydro = equation("columns", "hydrostatic")
    density = equation("columns", "ideal_gas_density")
    r = symbols(hydro)["r"]
    p0, t0, m0, g0, r0, kb = sp.symbols("p0 t0 m0 reference_gravity reference_radius kb", positive=True)
    profile = p0*sp.exp(m0*g0*r0**2/(kb*t0)*(1/r-1/r0))
    rho = density.rhs.subs({sp.Function("P")(r): profile, sp.Function("T")(r): t0,
                           sp.Function("m_mean_r")(r): m0})
    rho = replace_values(rho, k_B=kb)
    residual = hydro.lhs.subs(sp.Function("P")(r), profile).doit() - hydro.rhs.subs(sp.Function("rho")(r), rho)
    residual = replace_values(residual, g0=g0, r0=r0)
    assert sp.simplify(residual) == 0


def test_homogeneous_layer_limits_composition_and_equilibrium():
    solution = equation("emission", "homogeneous_layer").rhs
    n = symbols(solution)
    assert replace_values(solution, tau=0) == n["I_in"]
    assert sp.simplify(solution.subs(n["I_in"], n["S"])-n["S"]) == 0
    mu = sp.Symbol("mu_positive", positive=True)
    assert sp.limit(solution.subs(n["mu"], mu), n["tau"], sp.oo) == n["S"]
    a, b = sp.symbols("a b", nonnegative=True)
    first = replace_values(solution, tau=a)
    composed = replace_values(solution, tau=b, I_in=first)
    assert sp.simplify(composed-replace_values(solution, tau=a+b)) == 0


def test_formal_layer_solution_satisfies_the_ray_equation():
    ray = equation("emission", "transfer_equation")
    s = symbols(ray)["s"]
    a, source, incoming = sp.symbols("a source incoming", positive=True)
    solution = replace_values(equation("emission", "homogeneous_layer").rhs,
                              I_in=incoming, S=source, tau=a*s, mu=1)
    residual = (ray.lhs-ray.rhs).subs({sp.Function("I")(s): solution,
                                     sp.Function("alpha")(s): a,
                                     sp.Function("S_ray")(s): source}).doit()
    assert sp.simplify(residual) == 0


def test_isotropic_intensity_and_hemispheric_flux_have_one_pi():
    flux = equation("emission", "flux").rhs
    mu = symbols(flux)["mu"]
    intensity = sp.Symbol("intensity", positive=True)
    assert sp.simplify(flux.subs(sp.Function("I_mu")(mu), intensity).doit()-sp.pi*intensity) == 0


def test_constant_extinction_chord_recovers_geometric_path_length():
    relation = equation("transit", "chord_depth").rhs
    r = symbols(relation)["r"]
    actual = replace_values(relation.subs(sp.Function("alpha")(r), 2), b=3, R_top=5).doit()
    assert sp.simplify(actual-2*2*sp.sqrt(5**2-3**2)) == 0


def test_transit_area_limits_and_distinct_radius_normalizations():
    area = equation("transit", "blocked_area").rhs
    b = symbols(area)["b"]
    tau = sp.Symbol("tau", nonnegative=True)
    area = area.subs(sp.Function("tau_chord")(b), tau).doit()
    assert replace_values(area, tau=0) == symbols(area)["R0"]**2
    assert sp.simplify(sp.limit(area, tau, sp.oo)-symbols(area)["R_top"]**2) == 0
    normalized = replace_values(equation("transit", "normalized_radius").rhs, R_tr_sq=9, R0=2)
    depth = replace_values(equation("transit", "transit_depth").rhs, R_tr_sq=9, R_star=10)
    assert normalized == sp.Rational(9, 4)
    assert depth == sp.Rational(9, 100)


def test_observation_gaussian_normalization_and_rv_sign():
    kernel = equation("observation", "gaussian_kernel").rhs
    v = symbols(kernel)["v"]
    beta = sp.Symbol("beta_positive", positive=True)
    assert sp.integrate(replace_values(kernel, beta=beta), (v, -sp.oo, sp.oo)) == 1
    shift = equation("observation", "doppler_sampling").rhs
    assert replace_values(shift, RV=0) == symbols(shift)["w_obs"]
    assert float(replace_values(shift, RV=30, c=300000, w_obs=1000)) > 1000


def test_rotational_strip_integral_matches_quadratic_disk_intensity():
    raw = equation("observation", "rotation_raw").rhs
    x = symbols(raw)["x"]
    # The projected strip is integrated across the disk at a fixed x. Its
    # closed form is compared with independent numerical disk integration.
    intensity = equation("observation", "limb_darkening").rhs
    for offset in (0.0, 0.3, 0.8):
        expected = float(replace_values(raw, x=offset, u1=0.4, u2=0.2))
        with mp.workdps(30):
            f = sp.lambdify(symbols(intensity)["mu"], replace_values(intensity, u1=0.4, u2=0.2), "mpmath")
            half = mp.sqrt(1-mp.mpf(offset)**2)
            actual = mp.quad(lambda y: f(mp.sqrt(max(0, 1-mp.mpf(offset)**2-y*y))), [-half, half])
        assert expected == pytest.approx(float(actual), rel=2e-14)
    assert replace_values(raw, x=2) == 0
    assert replace_values(raw, x=-2) == 0


def test_continuous_rotation_normalization():
    raw = equation("observation", "rotation_raw").rhs
    x = symbols(raw)["x"]
    # On this integration domain Max(1-x**2, 0) equals 1-x**2.
    interior = raw.subs(sp.Max(1-x*x, 0), 1-x*x)
    area = sp.integrate(interior, (x, -1, 1))
    normalized = equation("observation", "rotation_normalized").rhs
    n = symbols(normalized)
    factor = normalized / sp.Function("K_raw")(n["v"]/n["vsini"])
    assert sp.simplify(area*n["vsini"]*factor) == 1


def test_folded_convolution_of_constant_spectrum_preserves_normalization():
    convolution = equation("observation", "convolution").rhs
    u, v = (symbols(convolution)[name] for name in ("u", "v"))
    gaussian = equation("observation", "gaussian_kernel").rhs
    beta = sp.Symbol("beta_positive", positive=True)
    gaussian = replace_values(gaussian, beta=beta)
    constant = sp.Symbol("constant", real=True)
    expression = convolution.xreplace({
        call: constant if call.func.__name__ == "F_in" else gaussian.subs(v, call.args[0])
        for call in convolution.function.atoms(sp.Function)
    })
    assert sp.simplify(expression.doit()-constant) == 0
