import os
import numpy as np
import pytest
import time

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.mark.quick
def test_triplet_speedup():
    a = np.random.rand(10000, 10)
    b = np.random.rand(10, 10000)
    c = np.random.rand(10000, 10)

    A = psi4.core.Matrix.from_array(a)
    B = psi4.core.Matrix.from_array(b)
    C = psi4.core.Matrix.from_array(c)

    start = time.time()
    R = psi4.core.triplet(A, B, C, False, False, False)
    time1 = time.time() - start

    start = time.time()
    S = psi4.core.doublet(A, B, False, False)
    S = psi4.core.doublet(S, C, False, False)
    time2 = time.time() - start

    if int(os.environ.get("PYTEST_XDIST_WORKER_COUNT", 1)) == 1:
        assert time1 < time2
    assert psi4.compare_matrices(R, S, 8, "linalg::triplet speedup test")


@pytest.mark.quick
def test_triplet_normal_case():
    a = np.random.rand(500, 500)
    b = np.random.rand(500, 500)
    c = np.random.rand(500, 500)

    A = psi4.core.Matrix.from_array(a)
    B = psi4.core.Matrix.from_array(b)
    C = psi4.core.Matrix.from_array(c)

    R = psi4.core.triplet(A, B, C, False, False, False)

    S = psi4.core.doublet(A, B, False, False)
    S = psi4.core.doublet(S, C, False, False)

    assert psi4.compare_matrices(R, S, 8, "linalg::triplet avg_case test")


@pytest.mark.quick
def test_triplet_with_irreps():
    A = psi4.core.Matrix.from_array([np.random.rand(10000, 5), np.random.rand(5, 10000), np.random.rand(10000, 5)])
    B = psi4.core.Matrix.from_array([np.random.rand(5, 10000), np.random.rand(10000, 5), np.random.rand(5, 10000)])
    C = psi4.core.Matrix.from_array([np.random.rand(10000, 5), np.random.rand(5, 10000), np.random.rand(10000, 5)])

    start = time.time()
    R = psi4.core.triplet(A, B, C, False, False, False)
    time1 = time.time() - start

    start = time.time()
    S = psi4.core.doublet(A, B, False, False)
    S = psi4.core.doublet(S, C, False, False)
    time2 = time.time() - start

    if int(os.environ.get("PYTEST_XDIST_WORKER_COUNT", 1)) == 1:
        assert time1 < time2
    assert psi4.compare_matrices(R, S, 8, "linalg::triplet irrep test")


@pytest.mark.quick
def test_triplet_with_transpose():
    A = psi4.core.Matrix.from_array([np.random.rand(5, 10000), np.random.rand(10000, 5), np.random.rand(5, 10000)])
    B = psi4.core.Matrix.from_array([np.random.rand(5, 10000), np.random.rand(10000, 5), np.random.rand(5, 10000)])
    C = psi4.core.Matrix.from_array([np.random.rand(10000, 5), np.random.rand(5, 10000), np.random.rand(10000, 5)])

    start = time.time()
    R = psi4.core.triplet(A, B, C, True, False, False)
    time1 = time.time() - start

    start = time.time()
    S = psi4.core.doublet(A, B, True, False)
    S = psi4.core.doublet(S, C, False, False)
    time2 = time.time() - start

    if int(os.environ.get("PYTEST_XDIST_WORKER_COUNT", 1)) == 1:
        assert time1 < time2
    assert psi4.compare_matrices(R, S, 8, "linalg::triplet transpose test")
