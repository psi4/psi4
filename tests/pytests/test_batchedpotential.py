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


@pytest.mark.parametrize("basis,puream", [
    pytest.param("cc-pvtz", True, id="cc-pvtz"),
    pytest.param("6-31g*", False, id="6-31gs-cartesian"),
])
def test_batchedpotential_contractions(basis, puream):
    """The density and charge contractions over 1100 points (more than one chunk of the drivers) match Libint2,
    for a nonsymmetric density."""
    import psi4

    psi4.core.clean()
    mol = _water()
    bas = psi4.core.BasisSet.build(mol, "ORBITAL", basis, puream=puream)
    rng = np.random.default_rng(5)
    points = _points(mol, 1100)
    D = rng.standard_normal((bas.nbf(), bas.nbf()))
    q = rng.standard_normal(points.shape[0])

    P = psi4.core.Matrix.from_array(points)
    V = np.asarray(psi4.core.batched_potential_contract_density(bas, bas, psi4.core.Matrix.from_array(D), P))
    F = np.asarray(psi4.core.batched_potential_contract_charges(bas, bas, psi4.core.Vector.from_array(q), P))

    # Libint2: the potential matrix of all the charges at once, and of each unit charge for the density
    ext = psi4.core.ExternalPotential()
    for qk, C in zip(q, points):
        ext.addCharge(qk, *C)
    F_ref = -np.asarray(ext.computePotentialMatrix(bas))
    V_ref = np.empty(points.shape[0])
    for k, C in enumerate(points[::10]):
        one = psi4.core.ExternalPotential()
        one.addCharge(1.0, *C)
        V_ref[10 * k] = -np.vdot(D, np.asarray(one.computePotentialMatrix(bas)))

    assert np.abs(F_ref).max() > 1.0
    assert compare_arrays(F_ref, F, 10, f"{basis} charges contraction")
    assert compare_arrays(V_ref[::10], V[::10], 10, f"{basis} density contraction")


def test_batchedpotential_zero_basis():
    """Against the zero basis (one constant function), as the density-fitted electrostatic potential of cube
    files uses them: the potential matrix of the nuclei, from charges Z_A at the nuclear positions, matches
    Libint2's nuclear attraction integrals between the two basis sets."""
    import psi4

    psi4.core.clean()
    mol = _water()
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "", "JKFIT", "cc-pvdz")
    zero = psi4.core.BasisSet.zero_ao_basis_set()
    nuclei = np.asarray(mol.geometry())
    Z = np.array([mol.Z(A) for A in range(mol.natom())])

    F = np.asarray(psi4.core.batched_potential_contract_charges(aux, zero, psi4.core.Vector.from_array(Z),
                                                                psi4.core.Matrix.from_array(nuclei)))
    ref = -np.asarray(psi4.core.MintsHelper(aux).ao_potential(aux, zero))
    assert F.shape == (aux.nbf(), 1)
    assert np.abs(ref).max() > 1.0
    assert compare_arrays(ref, F, 10, "aux x zero basis potential of the nuclei")


def test_electrostatic_potential_value():
    """MintsHelper.electrostatic_potential_value, which ddX uses: -q_k sum D (mu| 1/|r - C_k| |nu)."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": "cc-pvdz"})
    wfn = psi4.core.Wavefunction.build(mol, "cc-pvdz")
    bas = wfn.basisset()
    rng = np.random.default_rng(13)
    points = _points(mol, 25)
    D = rng.standard_normal((bas.nbf(), bas.nbf()))
    D = D + D.T
    q = rng.uniform(0.5, 2.0, points.shape[0])

    phi = np.asarray(psi4.core.MintsHelper(bas).electrostatic_potential_value(
        psi4.core.Vector.from_array(q), psi4.core.Matrix.from_array(points), psi4.core.Matrix.from_array(D)))
    ref = np.empty(points.shape[0])
    for k, C in enumerate(points):
        ext = psi4.core.ExternalPotential()
        ext.addCharge(q[k], *C)
        ref[k] = np.vdot(D, np.asarray(ext.computePotentialMatrix(bas)))
    assert compare_arrays(ref, phi, 10, "electrostatic_potential_value")
