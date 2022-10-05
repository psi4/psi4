"""
This file tests one-electron integrals from libmints computed with
the McMurchie-Davidson scheme
"""
import numpy as np
import pytest
import psi4
from pathlib import Path
import json


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.fixture
def reference_data():
    # Reference data generated using Psi4 1.5, still using OS86
    with open(Path(__file__).parent / "oei_reference_data.json") as f:
        reference_data = json.load(f)
    return reference_data


@pytest.fixture
def mol_h2o():
    return psi4.geometry("""
        units bohr
        O  0.000000000000  0.000000000000 -0.149544924691
        H  0.000000000000 -1.336238163149  1.186693238458
        H  0.000000000000  1.336238163149  1.186693238458
        no_com
        no_reorient
        symmetry c1
    """)


def matlist_to_ndarray(mats):
    """Converts a list of psi4.core.Matrix to np.ndarray"""
    return np.array([x.np for x in mats])


def cumulative_cart_dim(order):
    return (order + 1) * (order + 2) * (order + 3) // 6


def fdiff_multipole_integral(mol, basis_name, origin, order, step=1e-4):
    prefactors_5p = np.array([1.0, -8.0, 8.0, -1.0]) / 12.0
    multipliers_5p = [-2, -1, 1, 2]
    natoms = mol.natom()
    basis = psi4.core.BasisSet.build(mol, 'orbital', basis_name)
    nbf = basis.nbf()
    nmul = cumulative_cart_dim(order) - 1
    int_grad = np.zeros((natoms, 3, nmul, nbf, nbf))
    for i in range(natoms):
        for c in range(3):
            for f, p in zip(multipliers_5p, prefactors_5p):
                mol_p = mol.clone()
                coords_p = mol_p.geometry()
                coords_p.np[i, c] += f * step
                mol_p.set_geometry(coords_p)
                basis = psi4.core.BasisSet.build(mol_p, 'orbital', basis_name)
                mints = psi4.core.MintsHelper(basis)
                ints_pert = matlist_to_ndarray(mints.ao_multipoles(order, origin))
                int_grad[i, c, :, :, :] += p * ints_pert / step
    return int_grad
                

def test_mcmurchie_davidson_consistency_angmom(reference_data):
    assert reference_data['version'] == '1.5'
    refdata = reference_data['data']
    for molname in refdata:
        moldata = refdata[molname]
        mol = psi4.core.Molecule.from_string(moldata.pop('psi4string'))
        for basis_name in moldata:
            ref = moldata[basis_name]
            basis = psi4.core.BasisSet.build(mol, 'orbital', basis_name)
            mints = psi4.core.MintsHelper(basis)

            shape = ref[f'L_shape']
            integral_ref = np.array(ref["L"]).reshape(shape)
            ret = mints.ao_angular_momentum()
            ret_np = np.array([x.np for x in ret])
            np.testing.assert_allclose(ret_np, integral_ref, atol=1e-14)


def test_mcmurchie_davidson_multipoles(mol_h2o):
    basis = psi4.core.BasisSet.build(mol_h2o, 'orbital', 'cc-pvdz')
    mints = psi4.core.MintsHelper(basis)
    order = 6
    M = mints.ao_multipoles(order=order, origin=[0.0, 0.0, 0.0])

    with pytest.raises(Exception):
        # overlap integrals not accessible via multipole interface
        mints.ao_multipoles(order=0, origin=[0.0, 0.0, 0.0])

    with pytest.raises(Exception):
        # wrong origin specification
        mints.ao_multipoles(order=0, origin=[0.0, 0.0, 0.0, 0.0])

    Mnp = matlist_to_ndarray(M)
    # reference from l2
    dips = matlist_to_ndarray(mints.ao_dipole())
    quads = matlist_to_ndarray(mints.ao_quadrupole())
    
    np.testing.assert_allclose(dips, Mnp[:3], atol=1e-14)
    np.testing.assert_allclose(quads, Mnp[3:9], atol=1e-14)


def test_mcmurchie_davidson_multipoles_gradient(mol_h2o):
    psi4.set_options({'basis': 'cc-pvdz'})
    _, wfn = psi4.energy('HF', molecule=mol_h2o, return_wfn=True)

    order = 8
    ints_grad = fdiff_multipole_integral(mol_h2o, 'cc-pvdz', [1.0, 2.0, 3.0], order=order)
    nmul = cumulative_cart_dim(order) - 1
    grad_fdiff = np.einsum('ncmij,ij->ncm', ints_grad, wfn.Da().np).reshape(-1, nmul)

    mints = psi4.core.MintsHelper(wfn.basisset())
    # test against finite differences
    grad = mints.multipole_grad(D=wfn.Da(), order=order, origin=[1.0, 2.0, 3.0])
    np.testing.assert_allclose(grad_fdiff, grad.np, atol=1e-7)
    
    # test that we get the same result from the 'hard-wired' dipole_grad
    grad_dip = mints.dipole_grad(wfn.Da())
    grad = mints.multipole_grad(D=wfn.Da(), order=1, origin=[0.0, 0.0, 0.0])
    np.testing.assert_allclose(grad_dip.np, grad.np, atol=1e-14)
