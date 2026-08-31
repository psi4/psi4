import pytest
import numpy as np

import psi4
from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf, *using("einsums")]


def test_sad_guess_uhf_vs_cghf():
    """For any molecule, UHF and CGHF with GUESS=SAD must start from the same
    real, spin-restricted SAD density.

    The SAD guess is fundamentally real and spin-restricted (Da == Db).
    CGHF promotes it to a complex block-diagonal spinor density
       D_cghf = [[Da,  0 ],
                 [0 ,  Db]]
    with zero off-diagonal (alpha-beta coupling) blocks and zero imaginary
    part.  Because the underlying SADGuess object is shared, the UHF and
    CGHF initial densities are identical.

    Only the initial guess is tested: the SCF run is limited to a single
    iteration (maxiter=1) and allowed to fail on maxiter, so the captured
    density is exactly the SAD guess before any SCF relaxation.
    """

    # Capture the SAD guess density from inside the SCF driver before any
    # SCF iterations begin.
    guess_data = {}

    def capture_cghf_guess(wfn):
        guess_data["cghf_D"] = wfn.D().to_array().copy()

    def capture_uhf_guess(wfn):
        guess_data["uhf_Da"] = wfn.Da().to_array().copy()
        guess_data["uhf_Db"] = wfn.Db().to_array().copy()

    mol = psi4.geometry("""
        0 3
        O
        O 1 1.2
        symmetry c1
    """)

    common_options = {
        "basis": "6-31g",
        "scf_type": "direct",
        "df_scf_guess": False,
        "maxiter": 1,
        "fail_on_maxiter": False,
        "guess": "sad",
    }

    # ------------------------------------------------------------------
    # Run CGHF
    # ------------------------------------------------------------------
    psi4.set_options({**common_options, "reference": "cghf"})
    psi4.core.pre_scf_hook = capture_cghf_guess
    psi4.energy("scf", molecule=mol, return_wfn=True)

    # ------------------------------------------------------------------
    # Run UHF
    # ------------------------------------------------------------------
    psi4.set_options({**common_options, "reference": "uhf"})
    psi4.core.pre_scf_hook = capture_uhf_guess
    psi4.energy("scf", molecule=mol, return_wfn=True)

    # Make sure the hook does not leak to other tests.
    psi4.core.pre_scf_hook = None

    # ------------------------------------------------------------------
    # Extract densities
    # ------------------------------------------------------------------
    Da_uhf = guess_data["uhf_Da"]
    Db_uhf = guess_data["uhf_Db"]
    D_cghf = guess_data["cghf_D"]

    nbf = Da_uhf.shape[0]

    assert D_cghf.dtype == np.complex128, "CGHF density is not complex"
    assert D_cghf.shape == (2 * nbf, 2 * nbf), \
        f"CGHF density shape {D_cghf.shape} != expected {(2*nbf, 2*nbf)}"

    # ------------------------------------------------------------------
    # Decompose CGHF density into spin blocks
    # ------------------------------------------------------------------
    D_aa = D_cghf[:nbf, :nbf]    # alpha-alpha
    D_ab = D_cghf[:nbf, nbf:]    # alpha-beta
    D_ba = D_cghf[nbf:, :nbf]    # beta-alpha
    D_bb = D_cghf[nbf:, nbf:]    # beta-beta

    # Diagonal blocks must match UHF Da/Db (real part)
    np.testing.assert_allclose(D_aa.real, Da_uhf, atol=1e-12,
                               err_msg="CGHF D_aa != UHF Da")
    np.testing.assert_allclose(D_bb.real, Db_uhf, atol=1e-12,
                               err_msg="CGHF D_bb != UHF Db")

    # Imaginary parts of diagonal blocks must be ~zero
    np.testing.assert_allclose(D_aa.imag, 0.0, atol=1e-12,
                               err_msg="CGHF D_aa has non-zero imaginary part")
    np.testing.assert_allclose(D_bb.imag, 0.0, atol=1e-12,
                               err_msg="CGHF D_bb has non-zero imaginary part")

    np.testing.assert_allclose(np.abs(D_ab), 0.0, atol=1e-12,
                               err_msg="CGHF D alpha-beta block is non-zero")
    np.testing.assert_allclose(np.abs(D_ba), 0.0, atol=1e-12,
                               err_msg="CGHF D beta-alpha block is non-zero")

    # ------------------------------------------------------------------
    # Verify that if we build the expected SAD-style CGHF density from
    # UHF Da, it has the correct structure (block-diagonal, zero imag).
    # This models what CGHF::compute_SAD_guess() does internally.
    # ------------------------------------------------------------------
    D_sad_style = np.zeros((2 * nbf, 2 * nbf), dtype=np.complex128)
    D_sad_style[:nbf, :nbf] = Da_uhf
    D_sad_style[nbf:, nbf:] = Da_uhf   # SAD is spin-restricted: Da == Db

    # Imaginary part is identically zero
    np.testing.assert_allclose(D_sad_style.imag, 0.0, atol=1e-15)
    # Off-diagonal blocks are identically zero
    np.testing.assert_allclose(D_sad_style[:nbf, nbf:], 0.0, atol=1e-15)
    np.testing.assert_allclose(D_sad_style[nbf:, :nbf], 0.0, atol=1e-15)

    # Round-trip through psi4 ComplexMatrix to confirm API consistency
    D_cmplx = psi4.core.ComplexMatrix.from_array(D_sad_style, name="D_sad")
    np.testing.assert_allclose(D_cmplx.to_array(), D_sad_style, atol=1e-15)
