import numpy as np
import pytest

import psi4
from psi4.driver.p4util.solvers import davidson_solver, hamiltonian_solver
from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api]

def _diag_dom_sym_mat(size, sep, scale, sym=1.0):
    M = np.zeros((size, size))
    # positive diag elements, sep factor determines how close
    v = 1.0
    for i in range(size):
        v += sep
        M[i, i] = v
    # scale factor determines how large the off diagonal elements are
    M = M + scale * np.random.randn(size, size)
    # sym factor determines how symmetric the thing is, 1.0=completely , 0.0 not at all
    M = (M.T * (sym) + M) / 2.0
    return M


class SimulateBase:
    @staticmethod
    def vector_dot(X, Y):
        return X.dot(Y)

    @staticmethod
    def vector_axpy(a, X, Y):
        Y += a * X
        return Y

    @staticmethod
    def vector_scale(a, X):
        return a * X

    @staticmethod
    def vector_copy(X):
        return X.copy()

    def new_vector(self):
        return np.zeros((self.size, ))


class DSProblemSimulate(SimulateBase):
    "Provide the interface of an engine, around an actual matrix stored in memory"

    def __init__(self, size, sep=1.0, scale=0.0001, sym=1.0):
        self.size = size
        self.A = _diag_dom_sym_mat(size, sep=sep, scale=scale, sym=sym)

    def compute_products(self, X):
        Xmat = np.column_stack(X)
        prods = self.A.dot(Xmat)
        return list(prods.T), len(X)

    def precondition(self, R, shift):
        return R / (shift - np.diag(self.A))


class HSProblemSimulate(SimulateBase):
    "Provide the interface of an engine, around an actual matrix stored in memory"

    def __init__(self, size, a_sep=1.0, b_sep=0.001, scale=0.0001, a_sym=1.0, b_sym=1.0):
        self.size = size
        self.A = _diag_dom_sym_mat(size, sep=a_sep, scale=scale, sym=a_sym)
        self.B = _diag_dom_sym_mat(size, sep=b_sep, scale=scale, sym=b_sym)

    def compute_products(self, X):
        Xmat = np.column_stack(X)
        ApBx_prd = np.dot(self.A + self.B, Xmat)
        AmBx_prd = np.dot(self.A - self.B, Xmat)
        H1x = list(ApBx_prd.T)
        H2x = list(AmBx_prd.T)
        return H1x, H2x, len(X)

    def precondition(self, R, shift):
        return R / (shift - np.diag(self.A))


@pytest.mark.unittest
@pytest.mark.solver
def test_davidson_solver_numpy():
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = DSProblemSimulate(BIGDIM)
    ret = davidson_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        #Don't let the ss grow to size of real thing
        max_ss_size=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=100)

    test_vals = ret["eigvals"]
    assert test_vals is not None, "Solver Failed to converge"

    ref_vals, ref_vectors = np.linalg.eigh(test_engine.A)
    idx = ref_vals.argsort()[:nroot]
    ref_vals = ref_vals[idx]
    ref_vectors = ref_vectors[:, idx]

    compare_arrays(ref_vals, test_vals, 6, "Davidson eigenvalues")
    # NOTE: The returned eigenvectors are a list of tuples of (whatever the engine used is using for a vector)
    # So in this case, we need to put the columns together into a matrix to compare directly to the np.LA.eigh result
    test_vectors = [x[0] for x in ret["eigvecs"]]
    compare_arrays(ref_vectors, np.column_stack(test_vectors), 8, "Davidson eigenvectors")


@pytest.mark.unittest
@pytest.mark.solver
def test_hamiltonian_solver():
    BIGDIM = 100
    nroot = 3
    guess = list(np.eye(BIGDIM)[:, :nroot * 2].T)
    test_engine = HSProblemSimulate(BIGDIM)
    ret = hamiltonian_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        # Don't let the ss grow to size of real thing
        max_ss_size=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=100)

    test_vals = ret["eigvals"]
    assert test_vals is not None, "The solver failed to converge"

    # compute the reference values, Use the partially reduced non-hermitian form of the problem, solve for RH-eigenvectors
    ref_H = np.dot(test_engine.A - test_engine.B, test_engine.A + test_engine.B)
    ref_vals, ref_rvecs = np.linalg.eig(ref_H)
    # Associated eigenvalues are the squares of the 2N dimensional problem
    ref_vals = np.sqrt(ref_vals)
    # sort the values/vectors
    sort_idx = ref_vals.argsort()
    ref_vals = ref_vals[sort_idx]
    ref_rvecs = ref_rvecs[:, sort_idx]

    # truncate to number of roots found
    ref_vals = ref_vals[:nroot]
    ref_rvecs = ref_rvecs[:, :nroot]

    # compare values
    compare_arrays(ref_vals, test_vals, 5, "Hamiltonian Eigenvalues")

    # compare roots
    # NOTE: The returned eigenvectors are a list of (whatever the engine used is using for a vector)
    # So in this case, we need to put the columns together into a matrix to compare directly to the np.LA.eig result
    test_rvecs = [x[0] for x in ret["eigvecs"]]
    compare_arrays(ref_rvecs, np.column_stack(test_rvecs), 6, "Hamiltonian Right Eigenvectors")

    #Can't compute LH eigenvectors with numpy but we have the RH, and the matrix
    # solver computed Vl^T * H
    test_lvecs = [x[1] for x in ret["eigvecs"]]
    vl_T_H = np.column_stack([np.dot(test_lvecs[i], ref_H) for i in range(nroot)])
    # value * right_vector (should not matter which one since we just checked them equal)
    w_Vr = np.column_stack([test_rvecs[i] * test_vals[i] for i in range(nroot)])
    #compare
    compare_arrays(vl_T_H, w_Vr, 6, "Hamiltonian Left Eigenvectors")


@pytest.mark.xfail(True, reason='Stress Testing', run=True)
@pytest.mark.solver
@pytest.mark.stress
@pytest.mark.parametrize("sparsity", [10**(-x) for x in range(4, -2, -1)])
def test_davidson_solver_stress_sparsity(sparsity):
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = DSProblemSimulate(BIGDIM, scale=sparsity)
    test_vals, test_vectors, _ = davidson_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        #Don't let the ss grow to size of real thing
        max_vecs_per_root=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=1000)
    assert test_vals is not None, "Solver Failed to converge"

    ref_vals, ref_vectors = np.linalg.eigh(test_engine.A)
    idx = ref_vals.argsort()[:nroot]
    ref_vals = ref_vals[idx]
    ref_vectors = ref_vectors[:, idx]

    compare_arrays(ref_vals, test_vals, 6, "Davidson eigenvalues")
    # NOTE: The returned eigenvectors are a list of (whatever the engine used is using for a vector)
    # So in this case, we need to put the columns together into a matrix to compare directly to the np.LA.eigh result
    compare_arrays(ref_vectors, np.column_stack(test_vectors), 8, "Davidson eigenvectors")


@pytest.mark.xfail(True, reason='Stress Testing', run=True)
@pytest.mark.solver
@pytest.mark.stress
@pytest.mark.parametrize("sep", [10**(-x) for x in range(0, 8)])
def test_davidson_solver_stress_eigval_sep(sep):
    """Solver can miss "skip" a root if they are close in the range we are looking at"""
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = DSProblemSimulate(BIGDIM, sep=sep)
    test_vals, test_vectors, _ = davidson_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        #Don't let the ss grow to size of real thing
        max_vecs_per_root=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=1000)

    assert test_vals is not None, "Solver Failed to converge"
    ref_vals, ref_vectors = np.linalg.eigh(test_engine.A)
    idx = ref_vals.argsort()[:nroot]
    ref_vals = ref_vals[idx]
    ref_vectors = ref_vectors[:, idx]

    compare_arrays(ref_vals, test_vals, 6, "Davidson eigenvalues")
    # NOTE: The returned eigenvectors are a list of (whatever the engine used is using for a vector)
    # So in this case, we need to put the columns together into a matrix to compare directly to the np.LA.eigh result
    compare_arrays(ref_vectors, np.column_stack(test_vectors), 8, "Davidson eigenvectors")


@pytest.mark.xfail(True, reason='Stress Testing', run=True)
@pytest.mark.solver
@pytest.mark.stress
@pytest.mark.parametrize("sparsity", [10**(-x) for x in range(4, -2, -1)])
def test_hamiltonian_solver_stress_sparsity(sparsity):
    """Algorithm is challenged as the off diagonal elements approach the magnitude of the diagonal elements"""
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = HSProblemSimulate(BIGDIM, scale=sparsity)
    test_vals, test_rvecs, test_lvecs, _ = hamiltonian_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        # Don't let the ss grow to size of real thing
        max_vecs_per_root=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=1000)

    assert test_vals is not None, "The solver failed to converge"

    # compute the reference values, Use the partially reduced non-hermitian form of the problem, solve for RH-eigenvectors
    ref_H = np.dot(test_engine.A - test_engine.B, test_engine.A + test_engine.B)
    ref_vals, ref_rvecs = np.linalg.eig(ref_H)
    # Associated eigenvalues are the squares of the 2N dimensional problem
    ref_vals = np.sqrt(ref_vals)
    # sort the values/vectors
    sort_idx = ref_vals.argsort()
    ref_vals = ref_vals[sort_idx]
    ref_rvecs = ref_rvecs[:, sort_idx]

    # truncate to number of roots found
    ref_vals = ref_vals[:nroot]
    ref_rvecs = ref_rvecs[:, :nroot]

    # compare values
    compare_arrays(ref_vals, test_vals, 5, "Hamiltonian Eigenvalues")


@pytest.mark.xfail(True, reason='Stress Testing', run=True)
@pytest.mark.solver
@pytest.mark.stress
@pytest.mark.parametrize("b_sep", [10**(-x) for x in range(8, -1, -1)])
def test_hamiltonian_solver_stress_b_diag_grows(b_sep):
    """Solver encounters problems when the B diagonal elements approach A"""
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = HSProblemSimulate(BIGDIM, b_sep=b_sep)
    test_vals, test_rvecs, test_lvecs, _ = hamiltonian_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        # Don't let the ss grow to size of real thing
        max_vecs_per_root=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=1000)

    assert test_vals is not None, "The solver failed to converge"

    # compute the reference values, Use the partially reduced non-hermitian form of the problem, solve for RH-eigenvectors
    ref_H = np.dot(test_engine.A - test_engine.B, test_engine.A + test_engine.B)
    ref_vals, ref_rvecs = np.linalg.eig(ref_H)
    # Associated eigenvalues are the squares of the 2N dimensional problem
    ref_vals = np.sqrt(ref_vals)
    # sort the values/vectors
    sort_idx = ref_vals.argsort()
    ref_vals = ref_vals[sort_idx]
    ref_rvecs = ref_rvecs[:, sort_idx]

    # truncate to number of roots found
    ref_vals = ref_vals[:nroot]
    ref_rvecs = ref_rvecs[:, :nroot]

    # compare values
    compare_arrays(ref_vals, test_vals, 5, "Hamiltonian Eigenvalues")


@pytest.mark.xfail(True, reason='Stress Testing', run=True)
@pytest.mark.solver
@pytest.mark.stress
@pytest.mark.parametrize("a_sep", [10**(-x) for x in range(0, 8)])
def test_hamiltonian_solver_stress_eigval_sep(a_sep):
    BIGDIM = 100
    nroot = 3
    guess = list(np.random.randn(BIGDIM, nroot).T)
    test_engine = HSProblemSimulate(BIGDIM, a_sep=a_sep, b_sep=a_sep / 10)
    test_vals, test_rvecs, test_lvecs, _ = hamiltonian_solver(
        engine=test_engine,
        guess=guess,
        nroot=nroot,
        # Don't let the ss grow to size of real thing
        max_vecs_per_root=20,
        # Don't assault stdout with logging
        verbose=0,
        maxiter=1000)

    assert test_vals is not None, "The solver failed to converge"

    # compute the reference values, Use the partially reduced non-hermitian form of the problem, solve for RH-eigenvectors
    ref_H = np.dot(test_engine.A - test_engine.B, test_engine.A + test_engine.B)
    ref_vals, ref_rvecs = np.linalg.eig(ref_H)
    # Associated eigenvalues are the squares of the 2N dimensional problem
    ref_vals = np.sqrt(ref_vals)
    # sort the values/vectors
    sort_idx = ref_vals.argsort()
    ref_vals = ref_vals[sort_idx]
    ref_rvecs = ref_rvecs[:, sort_idx]

    # truncate to number of roots found
    ref_vals = ref_vals[:nroot]
    ref_rvecs = ref_rvecs[:, :nroot]

    # compare values
    compare_arrays(ref_vals, test_vals, 5, "Hamiltonian Eigenvalues")
