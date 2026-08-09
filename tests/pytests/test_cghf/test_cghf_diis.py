import pytest
import numpy as np
from scipy.optimize import minimize

import psi4
from addons import using
from psi4.driver.procrouting.diis import DIIS, RemovalPolicy, StoragePolicy
from psi4.driver.p4util.exceptions import SCFConvergenceError

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf, *using("einsums")]


def _co_mol():
    return psi4.geometry(
        """
        0 1
        C
        O 1 1.110
        symmetry c1
        """
    )


def _hermitian_dot(A, B):
    """Re(Tr(A^H B)), matching ComplexMatrix.vector_dot."""
    return np.vdot(A, B).real


def _manual_pulay_extrapolation(Fs, Es):
    """Reference Pulay/DIIS coefficients + extrapolated target, mirroring
    DIIS.diis_coefficients()/DIIS.extrapolate() but computed directly in NumPy with the
    Hermitian (conjugated) error metric appropriate for complex F/error matrices."""
    n = len(Fs)
    dim = n + 1
    B = np.zeros((dim, dim))
    for i in range(n):
        for j in range(n):
            B[i, j] = _hermitian_dot(Es[i], Es[j])
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


def _manual_adiis_extrapolation(Fs, Ds):
    """Reference ADIIS extrapolation for closed_shell=False (CGHF convention)."""
    n = len(Fs)
    dD = [D - Ds[-1] for D in Ds]
    dF = [F - Fs[-1] for F in Fs]
    linear = np.array([_hermitian_dot(d, Fs[-1]) for d in dD])
    quadratic = np.array([[_hermitian_dot(dD[i], dF[j]) for j in range(n)] for i in range(n)])

    def energy(x):
        return np.dot(linear, x) + np.einsum("i,ij,j->", x, quadratic, x) / 2

    def gradient(x):
        return linear + np.einsum("i,ij->j", x, quadratic)

    result = minimize(
        energy, np.ones(n), method="SLSQP",
        bounds=tuple((0, 1) for _ in range(n)),
        constraints=[{"type": "eq", "fun": lambda x: sum(x) - 1, "jac": lambda x: np.ones_like(x)}],
        jac=gradient, tol=5e-6, options={"maxiter": 200},
    )
    assert result.success
    return sum(c * F for c, F in zip(result.x, Fs))


def _manual_ediis_extrapolation(Fs, Ds, energies):
    """Reference EDIIS extrapolation for closed_shell=False (CGHF convention)."""
    n = len(Fs)
    DF = np.array([[_hermitian_dot(Ds[i], Fs[j]) for j in range(n)] for i in range(n)])
    diag = np.diag(DF)
    quadratic = diag[:, None] + diag - DF - DF.T
    quadratic *= -1 / 2
    linear = np.asarray(energies, dtype=float)

    def energy(x):
        return np.dot(linear, x) + np.einsum("i,ij,j->", x, quadratic, x) / 2

    def gradient(x):
        return linear + np.einsum("i,ij->j", x, quadratic)

    result = minimize(
        energy, np.ones(n), method="SLSQP",
        bounds=tuple((0, 1) for _ in range(n)),
        constraints=[{"type": "eq", "fun": lambda x: sum(x) - 1, "jac": lambda x: np.ones_like(x)}],
        jac=gradient, tol=5e-6, options={"maxiter": 200},
    )
    assert result.success
    return sum(c * F for c, F in zip(result.x, Fs))


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


@pytest.mark.parametrize("storage_policy", [StoragePolicy.InCore, StoragePolicy.OnDisk])
def test_adiis_extrapolation_with_complexmatrix_entries(storage_policy):
    """ADIIS works with ComplexMatrix target/density entries (closed_shell=False), exercising
    ComplexMatrix.subtract used by adiis_populate."""
    rng = np.random.default_rng(7)
    n = 4
    Fs = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    Ds = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    F_ref = _manual_adiis_extrapolation(Fs, Ds)

    dm = DIIS(6, "test ADIIS vector", RemovalPolicy.OldestAdded, storage_policy, False, engines={"adiis"})
    for F, D in zip(Fs, Ds):
        dm.add_entry({
            "target": [psi4.core.ComplexMatrix.from_array(F)],
            "densities": [psi4.core.ComplexMatrix.from_array(D)],
        })

    F_out = psi4.core.ComplexMatrix.from_array(np.zeros((n, n), dtype=complex))
    # Pure ADIIS only extrapolates when Dnorm > SCF_INITIAL_FINISH_DIIS_TRANSITION.
    performed = dm.extrapolate(F_out, Dnorm=1.0)

    assert performed == {"ADIIS"}
    np.testing.assert_allclose(F_out.to_array(), F_ref, atol=1e-8)


@pytest.mark.parametrize("storage_policy", [StoragePolicy.InCore, StoragePolicy.OnDisk])
def test_ediis_extrapolation_with_complexmatrix_entries(storage_policy):
    """EDIIS works with ComplexMatrix target/density entries and scalar energies
    (closed_shell=False)."""
    rng = np.random.default_rng(11)
    n = 4
    Fs = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    Ds = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
    energies = [-10.0, -10.5, -11.0]
    F_ref = _manual_ediis_extrapolation(Fs, Ds, energies)

    dm = DIIS(6, "test EDIIS vector", RemovalPolicy.OldestAdded, storage_policy, False, engines={"ediis"})
    for F, D, E in zip(Fs, Ds, energies):
        dm.add_entry({
            "target": [psi4.core.ComplexMatrix.from_array(F)],
            "densities": [psi4.core.ComplexMatrix.from_array(D)],
            "energy": [E],
        })

    F_out = psi4.core.ComplexMatrix.from_array(np.zeros((n, n), dtype=complex))
    performed = dm.extrapolate(F_out, Dnorm=1.0)

    assert performed == {"EDIIS"}
    np.testing.assert_allclose(F_out.to_array(), F_ref, atol=1e-8)


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
        "orbital_optimizer_package": "internal",
    })
    e_cghf = psi4.energy("scf", molecule=mol)

    assert e_cghf == pytest.approx(-112.750151, abs=1e-5)

    psi4.core.clean()
    psi4.core.clean_options()
    mol_rhf = _co_mol()
    psi4.set_options({"basis": "cc-pVDZ", "reference": "rhf", "scf_type": "direct"})
    e_rhf = psi4.energy("scf", molecule=mol_rhf)

    assert e_cghf == pytest.approx(e_rhf, abs=1e-6)


@pytest.mark.parametrize("accelerator", ["adiis", "ediis"])
def test_cghf_aediis_converges_co(accelerator):
    """CO/cc-pVDZ with core guess also converges when ADIIS or EDIIS is the initial
    accelerator (blended with DIIS), matching the RHF energy."""
    mol = _co_mol()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "scf_initial_accelerator": accelerator,
        "orbital_optimizer_package": "internal",
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
        "orbital_optimizer_package": "internal",
        "maxiter": 20,
    })
    with pytest.raises(SCFConvergenceError):
        psi4.energy("scf", molecule=mol)
