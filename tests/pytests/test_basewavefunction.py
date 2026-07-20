import pytest

import psi4

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


@pytest.mark.parametrize(
    "cls",
    [psi4.core.Wavefunction, psi4.core.ComplexWavefunction],
    ids=["Wavefunction", "ComplexWavefunction"],
)
def test_basewavefunction_constructors(cls, h2o_sto3g):
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
