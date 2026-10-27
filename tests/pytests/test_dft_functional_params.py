"""Parameters in dictionary-defined DFT functionals must be honored or rejected, never silently dropped.

* Component keys the builder doesn't read (e.g., "omega" on an ``xc_functionals`` kernel) raise.
* A ``tweak`` that changes a full XC kernel's exact-exchange parameters (e.g., "_omega") must reach
  the SuperFunctional too, since that is what sets the HF exchange built by JK.
"""
import pytest

import psi4
from psi4.driver.procrouting.dft import dft_builder

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.dft]


@pytest.mark.quick
@pytest.mark.parametrize("part,comp,params", [
    pytest.param("xc_functionals", "HYB_GGA_XC_WB97X", {"omega": 0.4}, id="xc-omega"),
    pytest.param("xc_functionals", "HYB_GGA_XC_B3LYP", {"alpha": 0.3}, id="xc-alpha"),
    pytest.param("c_functionals", "GGA_C_PBE", {"omega": 0.4}, id="c-omega"),
    pytest.param("x_functionals", "GGA_X_PBE", {"beta": 0.5}, id="x-beta"),
])
def test_unread_component_parameter_raises(part, comp, params):
    func = {"name": "custom_unread_param", part: {comp: params}}
    if part == "x_functionals":
        func["c_functionals"] = {"GGA_C_PBE": {}}
    elif part == "c_functionals":
        func["x_functionals"] = {"GGA_X_PBE": {}}

    with pytest.raises(psi4.ValidationError, match="Unrecognized parameter"):
        dft_builder.check_consistency(func)


@pytest.mark.quick
def test_xc_build_tweak_updates_exact_exchange():
    default = psi4.core.SuperFunctional.XC_build("XC_HYB_GGA_XC_WB97X", True, {})
    tweaked = psi4.core.SuperFunctional.XC_build("XC_HYB_GGA_XC_WB97X", True, {"_omega": 0.4})

    assert psi4.compare_values(0.3, default.x_omega(), 12, "wB97X default omega")
    assert psi4.compare_values(0.4, tweaked.x_omega(), 12, "wB97X tweaked omega reaches SuperFunctional")
    assert psi4.compare_values(0.4, tweaked.c_functionals()[0].query_libxc("XC_HYB_CAM_COEF")["OMEGA"], 12,
                               "wB97X tweaked omega reaches LibXC")
    assert psi4.compare_values(default.x_alpha(), tweaked.x_alpha(), 12, "wB97X alpha unchanged")
    assert psi4.compare_values(default.x_beta(), tweaked.x_beta(), 12, "wB97X beta unchanged")


def test_xc_tweak_omega_matches_dft_omega():
    """Tweaking a full XC kernel's _omega is equivalent to the DFT_OMEGA option, which sets both halves."""
    psi4.geometry("""
        0 1
        O
        H 1 0.96
        H 1 0.96 2 104.5
        symmetry c1
    """)
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "df", "e_convergence": 10, "d_convergence": 8})

    wb97x_om = {"name": "wb97x_omega04", "xc_functionals": {"HYB_GGA_XC_WB97X": {"tweak": {"_omega": 0.4}}}}
    e_tweak = psi4.energy("scf", dft_functional=wb97x_om)

    psi4.set_options({"dft_omega": 0.4})
    e_option = psi4.energy("wb97x")

    psi4.set_options({"dft_omega": 0.3})
    e_default = psi4.energy("wb97x")

    assert psi4.compare_values(e_option, e_tweak, 8, "wB97X _omega tweak vs DFT_OMEGA")
    assert abs(e_tweak - e_default) > 1.e-4, "omega change must be visible in the energy"
