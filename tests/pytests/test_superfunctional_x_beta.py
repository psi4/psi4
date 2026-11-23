"""SuperFunctional x_beta default (0.0) and custom LRC validation tests.

SuperFunctional::common_init() initializes x_beta_ = 0.0 so non-LRC functionals
(HF, pure DFAs, global hybrids) do not carry a stale x_beta = 1.0. Custom LRC
functionals that set a nonzero omega must specify beta explicitly (both in
dictionary definitions via dft_builder.check_consistency and in callable /
direct SuperFunctional builders via SuperFunctional::allocate).
"""
import pytest

import psi4
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.procrouting.dft import dft_builder, superfunctionals

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.dft]


def test_blank_superfunctional_x_beta_defaults_to_zero():
    sup = psi4.core.SuperFunctional.blank()
    assert sup.x_alpha() == 0.0
    assert sup.x_beta() == 0.0
    assert sup.x_omega() == 0.0
    assert sup.is_x_hybrid() is False
    assert sup.is_x_lrc() is False


@pytest.mark.parametrize("name,exp_alpha", [
    ("hf", 1.0),
    ("pbe", 0.0),
    ("pbe0", 0.25),
    ("b3lyp", 0.20),
])
def test_non_lrc_builtin_functionals_have_zero_x_beta(name, exp_alpha):
    sup, _ = superfunctionals.build_superfunctional(name, restricted=True)
    assert sup.is_x_lrc() is False
    assert sup.x_omega() == 0.0
    assert sup.x_beta() == 0.0
    assert sup.x_alpha() == pytest.approx(exp_alpha)


@pytest.mark.parametrize("name,exp_alpha,exp_beta,exp_omega", [
    ("wpbe", 0.0, 1.0, 0.4),
    ("wpbe0", 0.25, 0.75, 0.3),
    ("wb97x", 0.157706, 0.842294, 0.3),
    ("cam-b3lyp", 0.19, 0.46, 0.33),
])
def test_lrc_builtin_functionals_preserve_x_beta(name, exp_alpha, exp_beta, exp_omega):
    sup, _ = superfunctionals.build_superfunctional(name, restricted=True)
    assert sup.is_x_lrc() is True
    assert sup.x_alpha() == pytest.approx(exp_alpha, abs=1e-6)
    assert sup.x_beta() == pytest.approx(exp_beta, abs=1e-6)
    assert sup.x_omega() == pytest.approx(exp_omega, abs=1e-6)


def test_all_registered_builtin_functionals_pass_consistency():
    for func_dict in dft_builder.dict_functionals.values():
        dft_builder.check_consistency(func_dict)


def test_custom_dict_nonzero_omega_without_beta_raises():
    bad = {
        "name": "my_missing_beta_lrc",
        "x_functionals": {"GGA_X_HJS_PBE": {"omega": 0.4}},
        "x_hf": {"omega": 0.4},
        "c_functionals": {"GGA_C_PBE": {}},
    }
    with pytest.raises(ValidationError, match="x_hf\\['beta'\\]"):
        dft_builder.build_superfunctional_from_dictionary(bad, 1, 1, True)


def test_custom_dict_nonzero_beta_without_omega_raises():
    bad = {
        "name": "my_missing_omega_lrc",
        "x_functionals": {"GGA_X_PBE": {"alpha": 0.75}},
        "x_hf": {"alpha": 0.25, "beta": 0.75},
        "c_functionals": {"GGA_C_PBE": {}},
    }
    with pytest.raises(ValidationError, match="x_hf\\['omega'\\]"):
        dft_builder.build_superfunctional_from_dictionary(bad, 1, 1, True)


def test_custom_dict_explicit_beta_with_nonzero_omega_succeeds():
    explicit_zero = {
        "name": "my_explicit_zero_beta_lrc",
        "x_functionals": {"GGA_X_HJS_PBE": {"omega": 0.4}},
        "x_hf": {"alpha": 0.25, "beta": 0.0, "omega": 0.4},
        "c_functionals": {"GGA_C_PBE": {}},
    }
    sup, _ = dft_builder.build_superfunctional_from_dictionary(explicit_zero, 1, 1, True)
    assert sup.x_omega() == pytest.approx(0.4)
    assert sup.x_beta() == 0.0
    assert sup.x_alpha() == pytest.approx(0.25)

    explicit_lrc = {
        "name": "my_explicit_lrc",
        "x_functionals": {"GGA_X_HJS_PBE": {"omega": 0.35, "alpha": 0.8}},
        "x_hf": {"alpha": 0.2, "beta": 0.8, "omega": 0.35},
        "c_functionals": {"GGA_C_PBE": {}},
    }
    sup2, _ = dft_builder.build_superfunctional_from_dictionary(explicit_lrc, 1, 1, True)
    assert sup2.x_omega() == pytest.approx(0.35)
    assert sup2.x_beta() == pytest.approx(0.8)
    assert sup2.x_alpha() == pytest.approx(0.2)


def test_custom_callable_nonzero_omega_without_set_x_beta_raises():
    def build_missing_beta(name, npoints, deriv, restricted):
        sup = psi4.core.SuperFunctional.blank()
        sup.set_name("custom_no_beta")
        x_func = psi4.core.LibXCFunctional("XC_GGA_X_HJS_PBE", restricted)
        x_func.set_omega(0.4)
        sup.add_x_functional(x_func)
        sup.add_c_functional(psi4.core.LibXCFunctional("XC_GGA_C_PBE", restricted))
        sup.set_x_omega(0.4)
        return sup

    with pytest.raises(RuntimeError, match="set_x_beta"):
        superfunctionals.build_superfunctional(build_missing_beta, restricted=True)


def test_custom_callable_with_explicit_set_x_beta_succeeds():
    def build_with_beta(name, npoints, deriv, restricted):
        sup = psi4.core.SuperFunctional.blank()
        sup.set_name("custom_with_beta")
        x_func = psi4.core.LibXCFunctional("XC_GGA_X_HJS_PBE", restricted)
        x_func.set_omega(0.4)
        sup.add_x_functional(x_func)
        sup.add_c_functional(psi4.core.LibXCFunctional("XC_GGA_C_PBE", restricted))
        sup.set_x_omega(0.4)
        sup.set_x_beta(1.0)
        return sup

    sup, _ = superfunctionals.build_superfunctional(build_with_beta, restricted=True)
    assert sup.x_omega() == pytest.approx(0.4)
    assert sup.x_beta() == pytest.approx(1.0)
    sup_uks, _ = superfunctionals.build_superfunctional(build_with_beta, restricted=False)
    assert sup_uks.x_omega() == pytest.approx(0.4)
    assert sup_uks.x_beta() == pytest.approx(1.0)


def test_scf_jk_zeroes_omega_beta_for_non_lrc_and_zero_dft_omega():
    psi4.core.clean()
    psi4.core.clean_options()
    mol = psi4.geometry("He 0 0 0\nsymmetry c1")
    psi4.set_options({
        "basis": "cc-pvdz",
        "scf_type": "df",
        "save_jk": True,
    })

    _, wfn_hf = psi4.energy("hf", molecule=mol, return_wfn=True)
    assert wfn_hf.functional().x_beta() == 0.0
    assert wfn_hf.jk().get_omega() == 0.0
    assert wfn_hf.jk().get_omega_beta() == 0.0

    # Overriding DFT_OMEGA to 0.0 on a built-in LRC functional (wPBE reduces
    # smoothly to PBE at omega=0) leaves functional.x_beta() == 1.0, so the
    # initialize_jk guard is still required to keep jk.get_omega_beta() == 0.0.
    psi4.set_options({"dft_omega": 0.0})
    _, wfn_zero = psi4.energy("wpbe", molecule=mol, return_wfn=True)
    assert wfn_zero.functional().is_x_lrc() is False
    assert wfn_zero.functional().x_omega() == 0.0
    assert wfn_zero.functional().x_beta() == pytest.approx(1.0)
    assert wfn_zero.jk().get_omega() == 0.0
    assert wfn_zero.jk().get_omega_beta() == 0.0
    psi4.core.clean()
    psi4.core.clean_options()
