"""INTEGRAL_PACKAGE validates the selected engine and falls back for unsupported integral types."""

import pytest
from addons import using, uusing

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


@pytest.mark.skipif(psi4.addons("simint"), reason="requires a build without SIMINT")
def test_unavailable_simint_fails():
    with pytest.raises(RuntimeError, match=r"INTEGRAL_PACKAGE SIMINT.*built without it.*ENABLE_simint=ON"):
        _run("simint", "hf", psi4.energy, {"scf_type": "pk"})


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
@uusing("simint")
def test_simint_falls_back_to_libint2(method, driver, options):
    """Integral types SIMINT does not provide must come from Libint2 and agree with an all-Libint2 calculation."""
    ref = _run("libint2", method, driver, options)
    res = _run("simint", method, driver, options)
    if driver is psi4.energy:
        assert psi4.compare_values(ref, res, 9, f"{method} energy, SIMINT vs LIBINT2")
    else:
        assert psi4.compare_values(ref, res, 8, f"{method} {driver.__name__}, SIMINT vs LIBINT2")


@pytest.mark.parametrize(
    "package,driver,options,expected",
    [
        pytest.param("libint2", psi4.energy, {"scf_type": "pk"}, "libint2", id="libint2-energy"),
        pytest.param("simint", psi4.energy, {"scf_type": "pk"}, "simint", marks=using("simint"), id="simint-energy"),
        pytest.param("simint", psi4.gradient, {"scf_type": "pk"}, "simint>libint2", marks=using("simint"),
                     id="simint-gradient"),
        pytest.param("libcint", psi4.energy, {"scf_type": "df"}, "libcint", marks=using("cint"), id="libcint-df-energy"),
        pytest.param("libcint", psi4.gradient, {"scf_type": "pk"}, "libcint", marks=using("cint"),
                     id="libcint-gradient"),
    ],
)
def test_coulomb_eri_integrals_role(package, driver, options, expected):
    """The wavefunction records which engines computed the plain ERIs, in order of precedence."""
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(_H2O)
    psi4.set_options({"basis": "6-31g*", "integral_package": package})
    psi4.set_options(options)
    _, wfn = driver("hf", return_wfn=True)
    assert wfn.module_role("coulomb_eri_integrals") == expected



@uusing("simint")
def test_coulomb_eri_integrals_role_sequential():
    """Each driver call starts a fresh record: a SIMINT gradient followed by a SIMINT energy reports the
    gradient's Libint2 derivative integrals on the first wavefunction only."""
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(_H2O)
    psi4.set_options({"basis": "6-31g*", "scf_type": "pk", "integral_package": "simint"})
    _, gwfn = psi4.gradient("hf", return_wfn=True)
    _, ewfn = psi4.energy("hf", return_wfn=True)
    assert gwfn.module_role("coulomb_eri_integrals") == "simint>libint2"
    assert ewfn.module_role("coulomb_eri_integrals") == "simint"

def test_coulomb_eri_integrals_role_in_schema():
    """The role reaches QCSchema provenance, and each calculation starts a fresh record."""
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.set_options({"integral_package": "libint2"})
    psi4.core.IntegralFactory.reset_coulomb_eri_engines()
    atin = {
        "molecule": {"symbols": ["O", "H", "H"], "geometry": [0, 0, 0.221, 0, 1.430, -0.883, 0, -1.430, -0.883]},
        "specification": {
            "driver": "energy",
            "model": {"method": "hf", "basis": "6-31g*"},
            "keywords": {"scf_type": "pk", "integral_package": "libint2"},
        },
    }
    ret = psi4.schema_wrapper.run_qcschema(atin)
    assert ret.success, ret.error
    assert ret.provenance.model_dump().get("module_roles", {}).get("coulomb_eri_integrals") == "libint2"
