import pytest
import numpy as np

import psi4
from addons import uusing
from os.path import isfile

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.fixture
def h2o_sto3g():
    mol = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)
    psi4.set_options({"basis": "sto-3g"})
    basis = psi4.core.BasisSet.build(mol, "BASIS", psi4.core.get_global_option("BASIS"))
    return mol, basis


@pytest.mark.parametrize("cls_str", [
    "Wavefunction",
    pytest.param("ComplexWavefunction", marks=pytest.mark.xfail)
])
def test_basewavefunction_constructors(cls_str, h2o_sto3g):
    """Tests constructors for wavefunctions derived from BaseWavefunction."""
    cls = getattr(psi4.core, cls_str)
    mol, basis = h2o_sto3g
    options = psi4.core.get_options()

    wfn_mol_basis = cls(mol, basis)
    assert isinstance(wfn_mol_basis, psi4.core.BaseWavefunction)
    assert isinstance(wfn_mol_basis, cls)

    wfn_mol_basis_opts = cls(mol, basis, options)
    assert isinstance(wfn_mol_basis_opts, psi4.core.BaseWavefunction)
    assert isinstance(wfn_mol_basis_opts, cls)

    assert wfn_mol_basis.molecule().natom() == mol.natom()
    assert wfn_mol_basis.basisset().nbf() == basis.nbf()
    assert wfn_mol_basis_opts.molecule().natom() == mol.natom()
    assert wfn_mol_basis_opts.basisset().nbf() == basis.nbf()


@uusing("einsums")
@pytest.mark.xfail(
    condition=(not hasattr(psi4.core, "ComplexWavefunction")),
    reason="You must have psi4.core.ComplexWavefunction to run ComplexWavefunction pytests."
)
def test_complexwavefunction_to_from_file(h2o_sto3g):
    """ComplexWavefunction.to_file / from_file round-trip including ComplexMatrix (C1)."""
    mol, basis = h2o_sto3g
    nbf = basis.nbf()

    rng = np.random.default_rng(0)
    C_ref = rng.normal(size=(nbf, nbf)) + 1.0j * rng.normal(size=(nbf, nbf))
    D_ref = rng.normal(size=(nbf, nbf)) + 1.0j * rng.normal(size=(nbf, nbf))
    # Hermitize density-like matrix for realism
    D_ref = 0.5 * (D_ref + D_ref.conj().T)
    F_ref = rng.normal(size=(nbf, nbf)) + 1.0j * rng.normal(size=(nbf, nbf))
    H_ref = rng.normal(size=(nbf, nbf))
    S_ref = np.eye(nbf, dtype=np.complex128)
    eps_ref = rng.normal(size=(nbf,))

    wfn = psi4.core.ComplexWavefunction(mol, basis)
    wfn.set_name("TEST-CWFN")
    wfn.set_module("pytest")
    wfn.set_energy(-1.23456789)
    wfn.set_variable("CURRENT ENERGY", -1.23456789)
    wfn.set_variable("MY SCALAR", 42.0)
    wfn.set_C(psi4.core.ComplexMatrix.from_array(C_ref, name="C"))
    wfn.set_D(psi4.core.ComplexMatrix.from_array(D_ref, name="D"))
    wfn.set_F(psi4.core.ComplexMatrix.from_array(F_ref, name="F"))
    wfn.set_H(psi4.core.ComplexMatrix.from_array(H_ref, name="H"))
    wfn.set_S(psi4.core.ComplexMatrix.from_array(S_ref, name="S"))
    wfn.set_epsilon(psi4.core.Vector.from_array(eps_ref, name="epsilon"))
    wfn.set_array_variable("TEST ARR", psi4.core.ComplexMatrix.from_array(D_ref, name="TEST ARR"))

    wfn.to_file("pytest_cwfn")
    assert isfile("pytest_cwfn.npy")

    wfn_dict = wfn.to_file()
    np.testing.assert_allclose(wfn_dict["matrix"]["C"], C_ref)
    np.testing.assert_allclose(wfn_dict["matrix"]["D"], D_ref)
    np.testing.assert_allclose(wfn_dict["matrix"]["F"], F_ref)
    np.testing.assert_allclose(wfn_dict["matrix"]["H"], H_ref)
    np.testing.assert_allclose(wfn_dict["matrix"]["S"], S_ref)
    np.testing.assert_allclose(wfn_dict["vector"]["epsilon"], eps_ref)
    np.testing.assert_allclose(wfn_dict["matrixarr"]["TEST ARR"], D_ref)

    wfn_disk = psi4.core.ComplexWavefunction.from_file("pytest_cwfn")
    wfn_mem = psi4.core.ComplexWavefunction.from_file(wfn_dict)

    for restored in (wfn_disk, wfn_mem):
        assert isinstance(restored, psi4.core.ComplexWavefunction)
        assert restored.name() == "TEST-CWFN"
        assert restored.module() == "pytest"
        assert restored.energy() == pytest.approx(-1.23456789)
        assert restored.scalar_variable("MY SCALAR") == pytest.approx(42.0)
        assert restored.molecule().natom() == mol.natom()
        assert restored.basisset().nbf() == nbf
        assert restored.nelec() == wfn.nelec()
        assert restored.nirrep() == 1
        np.testing.assert_allclose(restored.C().to_array(), C_ref)
        np.testing.assert_allclose(restored.D().to_array(), D_ref)
        np.testing.assert_allclose(restored.F().to_array(), F_ref)
        np.testing.assert_allclose(restored.H().to_array(), H_ref)
        np.testing.assert_allclose(restored.S().to_array(), S_ref)
        np.testing.assert_allclose(restored.epsilon().to_array(), eps_ref)
        np.testing.assert_allclose(restored.array_variable("TEST ARR").to_array(), D_ref)
