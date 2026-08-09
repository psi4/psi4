"""
Test Psi4 INTERNAL DIIS for CGHF.

Note that the M-Intel runner forces orbital optimizer OOO by setting
ORBITAL_OPTIMIZER_PACKAGE to OOO (changing the default of INTERNAL). This does
not respect DIIS or SCF_INITIAL_ACCELERATOR options.

OOO-specific tests are in another file.
"""
import pytest
import numpy as np
from scipy.optimize import minimize

import psi4
from addons import using
from psi4.driver.procrouting.diis import DIIS, RemovalPolicy, StoragePolicy
from psi4.driver.p4util.exceptions import SCFConvergenceError

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf, *using("einsums")]


@pytest.fixture(autouse=True)
def reset_psi4_state():
    """Resets global Psi4 options and cleans up C++ core state before and after each test."""
    # Pre-test cleanup
    psi4.core.clean()
    psi4.core.clean_options()

    yield

    # Post-test cleanup
    psi4.core.clean()
    psi4.core.clean_options()


def _co_mol_str():
    return """
    0 1
    C
    O 1 1.110
    symmetry c1
    """

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
    """Reference ADIIS extrapolation for closed_shell=False."""
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
    """Reference EDIIS extrapolation for closed_shell=False."""
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
    zero/name/save/load methods added to ComplexMatrix."""
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


def test_cghf_co_does_not_converge_without_diis():
    """Sanity check for the test below: assert non-convergence without DIIS."""
    mol = psi4.geometry(_co_mol_str())
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "orbital_optimizer_package": "internal",
        "scf_initial_accelerator": "none",
        "maxiter": 20,
    })
    with pytest.raises(SCFConvergenceError):
        psi4.energy("scf", molecule=mol)


def test_cghf_diis_converges_co():
    """Confirms that the addition of DIIS converges CO/cc-pVDZ with GUESS=CORE.
    Otherwise oscillates between two energies forever. This fact is verified by
    test_cghf_co_does_not_converge_without_diis."""
    mol = psi4.geometry(_co_mol_str())
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "orbital_optimizer_package": "internal",
        "scf_initial_accelerator": "none",
    })
    e_cghf = psi4.energy("scf", molecule=mol)

    assert e_cghf == pytest.approx(-112.750151, abs=1e-5)


@pytest.mark.parametrize("accelerator", ["adiis", "ediis"])
def test_cghf_aediis_converges_co(accelerator):
    """CO/cc-pVDZ converges when ADIIS or EDIIS is the initial accelerator."""
    mol = psi4.geometry(_co_mol_str())
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "orbital_optimizer_package": "internal",
        "scf_initial_accelerator": accelerator,
    })
    e_cghf = psi4.energy("scf", molecule=mol)

    assert e_cghf == pytest.approx(-112.750151, abs=1e-5)
