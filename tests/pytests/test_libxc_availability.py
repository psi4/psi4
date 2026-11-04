"""Runtime LibXC-availability gating.

The DFT method table is built at import time by instantiating every functional
in the dictionary. Historically a functional whose underlying LibXC method was
absent from the linked LibXC build (e.g. renamed between LibXC versions) aborted
`import psi4` outright. These tests cover the guard that now checks each
functional's LibXC dependencies against the runtime library and skips /
diagnoses the unavailable ones instead of crashing.
"""
import pytest

import psi4
from psi4.driver.procrouting.dft import dft_builder

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

# A name that no LibXC version defines, used to simulate an unavailable method.
_FAKE = "XC_THIS_IS_NOT_A_REAL_LIBXC_FUNCTIONAL"


def test_available_true_for_real_functional():
    # LDA exchange is present in every LibXC build.
    assert psi4.core.LibXCFunctional.available("XC_LDA_X") is True


def test_available_false_for_missing_functional():
    assert psi4.core.LibXCFunctional.available(_FAKE) is False


def test_dependency_extraction():
    # A combined X+C definition depends on both LibXC pieces.
    d = {"name": "SVWN", "x_functionals": {"LDA_X": {}},
         "c_functionals": {"LDA_C_VWN_RPA": {}}}
    deps = dft_builder.libxc_functionals_in_dictionary(d)
    assert set(deps) == {"XC_LDA_X", "XC_LDA_C_VWN_RPA"}
    # The xc_functionals special case and an x_hf/use_libxc borrow are covered.
    d2 = {"name": "wB97X", "x_functionals": {"HYB_GGA_XC_WB97X": {"use_libxc": True}},
          "x_hf": {"use_libxc": "HYB_GGA_XC_WB97X"}}
    assert "XC_HYB_GGA_XC_WB97X" in dft_builder.libxc_functionals_in_dictionary(d2)


def test_functional_available_gate():
    real = {"name": "SVWN", "x_functionals": {"LDA_X": {}},
            "c_functionals": {"LDA_C_VWN_RPA": {}}}
    assert dft_builder.functional_available(real) is True

    bad = {"name": "BOGUS", "x_functionals": {"LDA_X": {}},
           "c_functionals": {_FAKE[3:]: {}}}  # strip the XC_ the builder re-adds
    assert dft_builder.functional_available(bad) is False
    assert dft_builder.unavailable_libxc_functionals(bad) == [_FAKE]


def test_build_raises_clear_error_when_unavailable():
    bad = {"name": "BOGUS", "x_functionals": {_FAKE[3:]: {}}, "c_functionals": {"LDA_C_VWN_RPA": {}}}
    with pytest.raises(Exception) as exc:
        dft_builder.build_superfunctional_from_dictionary(bad, 1, 1, True)
    msg = str(exc.value)
    assert _FAKE in msg and "LibXC" in msg


def test_th_fl_dual_spelling_registered():
    """TH-FL is GGA_XC_TH_FL in LibXC < 7.1 and LDA_XC_TH_FL in >= 7.1. The
    driver lists both spellings and registers whichever the linked LibXC
    provides, so TH-FL is a usable method on either version (and never aborts
    import on the other spelling)."""
    from psi4.driver.procrouting import proc_table
    has_th_fl = (psi4.core.LibXCFunctional.available("XC_GGA_XC_TH_FL")
                 or psi4.core.LibXCFunctional.available("XC_LDA_XC_TH_FL"))
    if not has_th_fl:
        pytest.skip("linked LibXC provides neither TH-FL spelling")
    assert "th-fl" in proc_table.procedures["energy"]


def test_available_functionals_still_registered():
    # Sanity: a staple functional that IS available remains a runnable method,
    # i.e. the gate did not over-prune.
    assert dft_builder.functional_available(dft_builder.functionals["pbe"]) is True
    psi4.core.clean()
    psi4.set_options({"scf_type": "pk"})
    mol = psi4.geometry("He 0 0 0\nsymmetry c1")
    _, wfn = psi4.energy("SVWN/cc-pvdz", return_wfn=True, molecule=mol)
    assert wfn is not None
    psi4.core.clean()
