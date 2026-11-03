"""INTEGRAL_PACKAGE validates the selected engine and falls back for unsupported integral types."""

import numpy as np
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
        pytest.param("auto", psi4.energy, {"scf_type": "pk"}, "libint2", id="auto-energy"),
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


# A fitting basis whose l=10 shell is past the angular momentum Libint2 is usually built for, and within libcint's
# l <= 15.
_AUX_L10 = """
    spherical
    ****
    H 0
    S 1 1.0
      2.0 1.0
    S 1 1.0
      0.5 1.0
    P 1 1.0
      1.0 1.0
    ****
    O 0
    S 1 1.0
      20.0 1.0
    S 1 1.0
      4.0 1.0
    S 1 1.0
      0.8 1.0
    P 1 1.0
      4.0 1.0
    P 1 1.0
      0.8 1.0
    D 1 1.0
      1.5 1.0
    L=10 1 1.0
      2.0 1.0
    ****
"""


def _setup_l10(package):
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.geometry(_H2O)
    psi4.basis_helper(_AUX_L10, name="aux_l10", key="df_basis_scf", set_option=True)
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "df", "e_convergence": 10, "d_convergence": 10,
                      "integral_package": package})


@uusing("cint")
def test_auto_falls_back_to_libcint(tmp_path):
    """Under AUTO, integrals past Libint2's angular momentum limit come from libcint, derivatives included, and the
    rest from Libint2."""
    _setup_l10("libcint")
    ref_grad, ref_wfn = psi4.gradient("hf", return_wfn=True)

    _setup_l10("auto")
    out = tmp_path / "auto.out"
    psi4.set_output_file(str(out), False)
    try:
        grad, wfn = psi4.gradient("hf", return_wfn=True)
    finally:
        psi4.set_output_file("output.dat", True)
    assert psi4.compare_values(ref_wfn.energy(), wfn.energy(), 10, "AUTO vs LIBCINT energy, l=10 fitting basis")
    assert psi4.compare_values(np.asarray(ref_grad), np.asarray(grad), 10, "AUTO vs LIBCINT gradient")
    assert wfn.module_role("coulomb_eri_integrals") == "libint2>libcint"
    text = out.read_text()
    for kind in ("3-center Coulomb integrals", "3-center Coulomb first-derivative integrals"):
        assert f"Libint2 was built for {kind} up to angular momentum" in text


def test_libint2_stops_past_its_limit():
    """LIBINT2 stops at integrals past Libint2's angular momentum limit and names the remedy."""
    _setup_l10("libint2")
    remedy = "INTEGRAL_PACKAGE AUTO or LIBCINT computes them" if psi4.addons("cint") else "ENABLE_cint=ON"
    try:
        with pytest.raises(RuntimeError, match=r"INTEGRAL_PACKAGE LIBINT2: Libint2 was built for [23]-center Coulomb "
                           r"integrals up to angular momentum \d+, but these need 10\. .*" + remedy):
            psi4.energy("hf")
    finally:
        psi4.core.clean_timers()  # the failed calculation leaves its outer timers on



@pytest.mark.parametrize("package,analytic", [
    pytest.param("auto", True, marks=using("cint"), id="auto"),
    pytest.param("libcint", True, marks=using("cint"), id="libcint"),
    pytest.param("libint2", False, id="libint2"),
])
def test_hessian_without_libint2_hessian_ints(package, analytic, monkeypatch):
    """With a Libint2 built without Hessian integrals (as on Windows), the driver keeps SCF Hessians analytic when
    libcint can compute those integrals, and otherwise falls back to finite differences of gradients."""
    from psi4.driver import driver_util, p4util

    real = p4util.libint2_configuration()
    monkeypatch.setattr(p4util, "libint2_configuration", lambda: dict(real, eri=[real["eri"][0], real["eri"][1], None]))
    psi4.core.clean_options()
    psi4.set_options({"integral_package": package})
    dertype, _ = driver_util.highest_analytic_derivative_available("hf", proc=None, managed_keywords={})
    assert dertype == (2 if analytic else 1)
