import pytest
import numpy as np

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]


def _random_complex(shape, rng):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _build_complex_jk():
    """Minimal ComplexDirectJK via ComplexJK.build_JK (SCF_TYPE=DIRECT)."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    return psi4.core.ComplexJK.build_JK(basis, basis)


def _run_compute_D(C_left, C_right=None):
    """Queue one or more C pairs and form D. C_* may be a matrix or a sequence of matrices."""
    if C_right is None:
        C_right = C_left

    lefts = C_left if isinstance(C_left, (list, tuple)) else [C_left]
    rights = C_right if isinstance(C_right, (list, tuple)) else [C_right]
    assert len(lefts) == len(rights)

    jk = _build_complex_jk()
    jk.C_clear()
    for Cl, Cr in zip(lefts, rights):
        jk.C_left_add(Cl)
        jk.C_right_add(Cr)
    jk.compute_D()
    assert len(jk.D()) == len(lefts)
    return jk.D()


def test_complexjk_compute_D_c1_rectangular():
    """C1 non-square C: D = C @ C.conj().T."""
    rng = np.random.default_rng(0)
    C_arr = _random_complex((5, 2), rng)
    C = psi4.core.ComplexMatrix.from_array(C_arr, name="C")

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), C_arr @ C_arr.conj().T)


def test_complexjk_compute_D_c1_asymmetric():
    """C1 non-square C_left != C_right: D = C_left @ C_right.conj().T."""
    rng = np.random.default_rng(2)
    Cl_arr = _random_complex((5, 2), rng)
    Cr_arr = _random_complex((5, 2), rng)
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), Cl_arr @ Cr_arr.conj().T)
    # Sanity: not the same as the symmetric products
    assert not np.allclose(D.to_array(), Cl_arr @ Cl_arr.conj().T)
    assert not np.allclose(D.to_array(), Cr_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_c1_rectangular_D():
    """C_left/C_right with different AO counts yield rectangular D = Cl @ Cr.H."""
    rng = np.random.default_rng(4)
    Cl_arr = _random_complex((5, 2), rng)  # 5 AOs × 2 occ
    Cr_arr = _random_complex((3, 2), rng)  # 3 AOs × 2 occ
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    D_arr = D.to_array()
    assert D_arr.shape == (5, 3)
    np.testing.assert_allclose(D_arr, Cl_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_multi_irrep_rectangular():
    """Multi-irrep non-square C: D_h = C_h @ C_h.conj().T per diagonal tile."""
    rng = np.random.default_rng(1)
    row_sizes = [4, 2, 3]
    col_sizes = [1, 2, 1]
    C = psi4.core.ComplexMatrix("C", row_sizes, col_sizes)

    refs = []
    for h, view in enumerate(C.array_interface()):
        C_h = _random_complex(view.shape, rng)
        view[:] = C_h
        refs.append(C_h @ C_h.conj().T)

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == len(row_sizes)
    D_blocks = D.to_array()
    assert isinstance(D_blocks, list)
    assert len(D_blocks) == len(row_sizes)
    for got, ref, nso in zip(D_blocks, refs, row_sizes):
        assert got.shape == (nso, nso)
        np.testing.assert_allclose(got, ref)


def test_complexjk_compute_D_multiple_C_pairs():
    """Multiple queued C_left/C_right pairs each produce their own D = Cl @ Cr.H."""
    rng = np.random.default_rng(3)
    shapes = [(4, 1), (5, 2), (3, 2)]
    lefts, rights, refs = [], [], []
    for i, shape in enumerate(shapes):
        Cl_arr = _random_complex(shape, rng)
        Cr_arr = _random_complex(shape, rng)
        lefts.append(psi4.core.ComplexMatrix.from_array(Cl_arr, name=f"Cl{i}"))
        rights.append(psi4.core.ComplexMatrix.from_array(Cr_arr, name=f"Cr{i}"))
        refs.append(Cl_arr @ Cr_arr.conj().T)

    Ds = _run_compute_D(lefts, rights)
    assert len(Ds) == len(shapes)
    for D, ref in zip(Ds, refs):
        assert D.num_blocks() == 1
        np.testing.assert_allclose(D.to_array(), ref)
