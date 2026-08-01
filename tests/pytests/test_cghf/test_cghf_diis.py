import pytest
import numpy as np

import psi4
from psi4.driver.procrouting.diis import DIIS, RemovalPolicy, StoragePolicy
from psi4.driver.p4util.exceptions import SCFConvergenceError

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf]


def _co_mol():
    return psi4.geometry(
        """
        0 1
        C
        O 1 1.110
        symmetry c1
        """
    )


def _manual_pulay_extrapolation(Fs, Es):
    """Reference Pulay/DIIS coefficients + extrapolated target, mirroring
    DIIS.diis_coefficients()/DIIS.extrapolate() but computed directly in NumPy with the
    Hermitian (conjugated) error metric appropriate for complex F/error matrices."""
    n = len(Fs)
    dim = n + 1
    B = np.zeros((dim, dim))
    for i in range(n):
        for j in range(n):
            B[i, j] = np.vdot(Es[i], Es[j]).real
    B[-1, :-1] = B[:-1, -1] = -1

    rhs = np.zeros(dim)
    rhs[-1] = -1

    diagonals = B.diagonal().copy()
    diagonals[-1] = 1
    if np.all(diagonals > 0):
        d = diagonals ** (-0.5)
        Bp = np.einsum("i,ij,j->ij", d, B, d)
        coeffs = np.linalg.lstsq(Bp, rhs, rcond=None)[0][:-1] * d[:-1]
    else:
        coeffs = np.linalg.lstsq(B, rhs, rcond=None)[0][:-1]

    return sum(c * F for c, F in zip(coeffs, Fs))


@pytest.mark.parametrize("storage_policy", [StoragePolicy.InCore, StoragePolicy.OnDisk])
def test_diis_extrapolation_with_complexmatrix_entries(storage_policy):
    """DIIS (the same class RHF/UHF/ROHF use) works unmodified with core.ComplexMatrix
    target/error entries, for both storage policies, thanks to the clone/axpy/vector_dot/
    zero/name/save/load pybind lambdas added to ComplexMatrix."""
    rng = np.random.default_rng(5)
    n = 4
    Fs = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    Es = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    F_ref = _manual_pulay_extrapolation(Fs, Es)

    dm = DIIS(6, "test DIIS vector", RemovalPolicy.LargestError, storage_policy, engines={"diis"})
    for F, E in zip(Fs, Es):
        dm.add_entry({
            "target": [psi4.core.ComplexMatrix.from_array(F)],
            "error": [psi4.core.ComplexMatrix.from_array(E)],
        })

    F_out = psi4.core.ComplexMatrix.from_array(np.zeros((n, n), dtype=complex))
    performed = dm.extrapolate(F_out, Dnorm=0.0)

    assert performed == {"DIIS"}
    np.testing.assert_allclose(F_out.to_array(), F_ref, atol=1e-10)


def test_cghf_diis_converges_co():
    """CO/cc-pVDZ is a classic core-guess DIIS torture test: with GUESS=CORE and no
    convergence acceleration it oscillates between two energies forever (see
    test_cghf_co_does_not_converge_without_diis below). With CGHF's DIIS implemented,
    it should converge cleanly to the closed-shell RHF energy (no spin-orbit coupling),
    matching the ~-112.750151 Eh reference quoted from other software.
    """
    mol = _co_mol()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "scf_initial_accelerator": "none",
    })
    e_cghf = psi4.energy("scf", molecule=mol)

    assert e_cghf == pytest.approx(-112.750151, abs=1e-5)

    psi4.core.clean()
    psi4.core.clean_options()
    mol_rhf = _co_mol()
    psi4.set_options({"basis": "cc-pVDZ", "reference": "rhf", "scf_type": "direct"})
    e_rhf = psi4.energy("scf", molecule=mol_rhf)

    assert e_cghf == pytest.approx(e_rhf, abs=1e-6)


def test_cghf_co_does_not_converge_without_diis():
    """Sanity check for the DIIS test above: CO's core-guess oscillation is real (not an
    artifact of some other change), so the previous test is actually exercising DIIS."""
    mol = _co_mol()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
        "maxiter": 20,
    })
    with pytest.raises(SCFConvergenceError):
        psi4.energy("scf", molecule=mol)
