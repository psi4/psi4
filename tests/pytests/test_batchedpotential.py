"""Potential integrals of a shell pair over a batch of points (BatchedPotentialInt).

The reference is Libint2 through ExternalPotential: the potential matrix of a unit
point charge at C is V(C) = -(mu| 1/|r - C| |nu).
"""

import numpy as np
import pytest

from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def _water():
    import psi4

    return psi4.geometry("""
        O  0.000000  0.000000  0.117790
        H  0.000000  0.755453 -0.471160
        H  0.000000 -0.755453 -0.471160
        units angstrom
        symmetry c1
        no_com
        no_reorient
        """)


def _points(mol, n, seed=11):
    """n points [Bohr]: random ones around and far from the molecule, and one on each nucleus."""
    rng = np.random.default_rng(seed)
    xyz = np.asarray(mol.geometry())
    pts = rng.uniform(-4.0, 4.0, size=(n - xyz.shape[0] - 4, 3))
    far = rng.normal(size=(4, 3))
    far *= (12.0 / np.linalg.norm(far, axis=1))[:, None]
    return np.vstack([xyz, pts, far])


def _batch_tensor(bs1, bs2, points):
    """A[mu, nu, k] from BatchedPotentialInt, assembled shell pair by shell pair."""
    import psi4

    batch = psi4.core.BatchedPotentialInt(bs1, bs2)
    P = psi4.core.Matrix.from_array(points)
    A = np.zeros((bs1.nbf(), bs2.nbf(), points.shape[0]))
    for s1 in range(bs1.nshell()):
        o1, n1 = bs1.shell_to_basis_function(s1), bs1.shell(s1).nfunction
        for s2 in range(bs2.nshell()):
            o2, n2 = bs2.shell_to_basis_function(s2), bs2.shell(s2).nfunction
            block = batch.compute_shell(s1, s2, P).np
            A[o1:o1 + n1, o2:o2 + n2, :] = block.reshape(n1, n2, -1)
    return A


@pytest.mark.parametrize("basis,puream", [
    pytest.param("cc-pvdz", True, id="cc-pvdz"),
    pytest.param("6-31g*", False, id="6-31gs-cartesian"),
    pytest.param("cc-pvqz", True, id="cc-pvqz"),
    pytest.param("cc-pvtz", False, id="cc-pvtz-cartesian"),
])
def test_batchedpotential_vs_libint2(basis, puream):
    """Every integral at every point matches Libint2, including points on the nuclei and far away, and more
    points than one chunk of the kernel holds."""
    import psi4

    psi4.core.clean()
    mol = _water()
    bas = psi4.core.BasisSet.build(mol, "ORBITAL", basis, puream=puream)
    points = _points(mol, 150)

    A = _batch_tensor(bas, bas, points)
    assert np.abs(A).max() > 1.0
    for k, C in enumerate(points):
        ext = psi4.core.ExternalPotential()
        ext.addCharge(1.0, *C)
        ref = -np.asarray(ext.computePotentialMatrix(bas))
        assert compare_arrays(ref, A[:, :, k], 12, f"{basis} puream={puream} point {k}")


def test_batchedpotential_two_bases():
    """Bra and ket from different basis sets, spherical and cartesian: swapping them transposes the integrals."""
    import psi4

    psi4.core.clean()
    mol = _water()
    bs1 = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvtz", puream=True)
    bs2 = psi4.core.BasisSet.build(mol, "ORBITAL", "6-31g*", puream=False)
    points = _points(mol, 40)

    A12 = _batch_tensor(bs1, bs2, points)
    A21 = _batch_tensor(bs2, bs1, points)
    assert np.abs(A12).max() > 1.0
    assert compare_arrays(A12, A21.transpose(1, 0, 2), 12, "A(bs1, bs2) = A(bs2, bs1)^T")
