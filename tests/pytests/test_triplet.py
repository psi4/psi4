import numpy as np
import pytest
import time

import psi4

@pytest.mark.quick
def test_triplet_speedup():
    a = np.random.rand(1000, 100);
    b = np.random.rand(100, 1000);
    c = np.random.rand(1000, 100);

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

    assert (time1 < time2)
    assert psi4.compare_matrices(R, S, 5, "linalg::triplet speedup test")

@pytest.mark.quick
def test_triplet_normal_case():
    a = np.random.rand(500, 500);
    b = np.random.rand(500, 500);
    c = np.random.rand(500, 500);

    A = psi4.core.Matrix.from_array(a)
    B = psi4.core.Matrix.from_array(b)
    C = psi4.core.Matrix.from_array(c)

    R = psi4.core.triplet(A, B, C, False, False, False)

    S = psi4.core.doublet(A, B, False, False)
    S = psi4.core.doublet(S, C, False, False)

    assert psi4.compare_matrices(R, S, 8, "linalg::triplet avg_case test")
