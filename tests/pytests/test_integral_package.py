"""INTEGRAL_PACKAGE falls back to Libint2 for whatever the chosen engine does not provide."""

import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

_H2O = """
0 1
O  0.000  0.000  0.117
H  0.000  0.757 -0.467
H  0.000 -0.757 -0.467
symmetry c1
no_reorient
no_com
"""


def _run(package, method, driver, options):
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(_H2O)
    psi4.set_options({"basis": "6-31g*", "e_convergence": 10, "d_convergence": 10, "integral_package": package})
    psi4.set_options(options)
    return driver(method)


@pytest.mark.parametrize(
    "method,driver,options",
    [
        pytest.param("hf", psi4.gradient, {"scf_type": "pk"}, id="pk-hf-gradient"),
        pytest.param("hf", psi4.gradient, {"scf_type": "df"}, id="df-hf-gradient"),
        pytest.param("hf", psi4.hessian, {"scf_type": "pk"}, id="pk-hf-hessian"),
        pytest.param("wb97x", psi4.energy, {"scf_type": "pk", "dft_radial_points": 75, "dft_spherical_points": 302},
                     id="pk-wb97x-energy"),
    ],
)
def test_simint_falls_back_to_libint2(method, driver, options):
    """SIMINT has no derivative or range-separated integrals (and may not be compiled in at all); those must come
    from Libint2 instead of aborting, and agree with an all-Libint2 calculation."""
    ref = _run("libint2", method, driver, options)
    res = _run("simint", method, driver, options)
    if driver is psi4.energy:
        assert psi4.compare_values(ref, res, 9, f"{method} energy, SIMINT vs LIBINT2")
    else:
        assert psi4.compare_values(ref, res, 8, f"{method} {driver.__name__}, SIMINT vs LIBINT2")
