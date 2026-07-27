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


def _run_compute_D(C):
    jk = _build_complex_jk()
    jk.C_clear()
    jk.C_left_add(C)
    jk.C_right_add(C)
    jk.compute_D()
    assert len(jk.D()) == 1
    return jk.D()[0]


def test_complexjk_compute_D_c1_rectangular():
    """C1 non-square C: D = C @ C.conj().T."""
    rng = np.random.default_rng(0)
    C_arr = _random_complex((5, 2), rng)
    C = psi4.core.ComplexMatrix.from_array(C_arr, name="C")

    D = _run_compute_D(C)
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), C_arr @ C_arr.conj().T)


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

    D = _run_compute_D(C)
    assert D.num_blocks() == len(row_sizes)
    D_blocks = D.to_array()
    assert isinstance(D_blocks, list)
    assert len(D_blocks) == len(row_sizes)
    for got, ref, nso in zip(D_blocks, refs, row_sizes):
        assert got.shape == (nso, nso)
        np.testing.assert_allclose(got, ref)
