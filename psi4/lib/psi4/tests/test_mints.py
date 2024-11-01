import pytest

import numpy as np
import psi4

from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_overlap_obs():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    psi4.set_options({'basis': 'aug-cc-pvdz'})

    conv = psi4.core.BasisSet.build(h2o,'BASIS', psi4.core.get_global_option('BASIS'))

    wfn = psi4.core.Wavefunction.build(h2o, psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    case1 = mints.ao_overlap()
    case2 = mints.ao_overlap(wfn.basisset(), wfn.basisset())

    assert psi4.compare_matrices(case1, case2, 10, "OVERLAP_TEST")  # TEST

def test_overlap_aux():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    psi4.set_options({'basis': 'aug-cc-pvdz',
                      'df_basis_mp2':'aug-cc-pvdz-ri'})

    conv = psi4.core.BasisSet.build(h2o,'BASIS', psi4.core.get_global_option('BASIS'))
    aux = psi4.core.BasisSet.build(h2o,'DF_BASIS_MP2',"", "RIFIT", psi4.core.get_global_option('DF_BASIS_MP2'))

    wfn = psi4.core.Wavefunction.build(h2o, psi4.core.get_global_option('BASIS'))
    mints = psi4.core.MintsHelper(wfn.basisset())

    tr = mints.ao_overlap(aux, aux).trace()

    assert psi4.compare_values(118, tr, 12, 'Test that diagonal elements of AO Overlap are 1.0')  # TEST

def test_export_ao_elec_dip_deriv():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    mints = psi4.core.MintsHelper(wfn.basisset())

    natoms = h2o.natom()
    cart = ['_X', '_Y', '_Z']

    D = wfn.Da()
    D.add(wfn.Db())
    D_np = np.asarray(D)

    deriv1_mat = {}
    deriv1_np = {}

    MU_Gradient = np.zeros((3 * natoms, 3))

    for atom in range(natoms):
        deriv1_mat["MU_" + str(atom)] = mints.ao_elec_dip_deriv1(atom)
        for mu_cart in range(3):
            for atom_cart in range(3):
                map_key = "MU" + cart[mu_cart] + "_" + str(atom) + cart[atom_cart]
                deriv1_np[map_key] = np.asarray(deriv1_mat["MU_" + str(atom)][3 * mu_cart + atom_cart])
                MU_Gradient[3 * atom + atom_cart, mu_cart] += np.einsum("uv,uv->", deriv1_np[map_key], D_np)

    PSI4_MU_Grad = mints.dipole_grad(D)
    G_python_MU_mat = psi4.core.Matrix.from_array(MU_Gradient)
    assert psi4.compare_matrices(PSI4_MU_Grad, G_python_MU_mat, 10, "DIPOLE_GRADIENT_TEST")  # TEST

def test_export_ao_overlap_half_deriv():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-PVDZ', molecule=h2o, return_wfn=True)
    C = wfn.Ca_subset("AO", "ALL")

    mints = psi4.core.MintsHelper(wfn.basisset())

    natoms = h2o.natom()
    cart = ['_X', '_Y', '_Z']

    deriv1_mat = {}
    deriv1_np = {}

    # Get total overlap derivative integrals along with both left and right half-derivative integrals
    for atom in range(natoms):
        deriv1_mat["S_LEFT_HALF_" + str(atom)] = mints.mo_overlap_half_deriv1("LEFT", atom, C, C)
        deriv1_mat["S_RIGHT_HALF_" + str(atom)] = mints.mo_overlap_half_deriv1("RIGHT", atom, C, C)
        deriv1_mat["S_" + str(atom)] = mints.mo_oei_deriv1("OVERLAP", atom, C, C)
        for atom_cart in range(3):
            map_key1 = "S_LEFT_HALF_" + str(atom) + cart[atom_cart]
            map_key2 = "S_RIGHT_HALF_" + str(atom) + cart[atom_cart]
            map_key3 = "S_" + str(atom) + cart[atom_cart]
            deriv1_np[map_key1] = np.asarray(deriv1_mat["S_LEFT_HALF_" + str(atom)][atom_cart])
            deriv1_np[map_key2] = np.asarray(deriv1_mat["S_RIGHT_HALF_" + str(atom)][atom_cart])
            deriv1_np[map_key3] = np.asarray(deriv1_mat["S_" + str(atom)][atom_cart])

            # Test (S_ii)^x = 2 * < i^x | i >
            assert compare_arrays(deriv1_np[map_key1].diagonal(), deriv1_np[map_key3].diagonal()/2)

            # Test (S_ii)^x = 2 * < i | i^x >
            assert compare_arrays(deriv1_np[map_key2].diagonal(), deriv1_np[map_key3].diagonal()/2)

            # Test (S_ij)^x = < i^x | j > + < j^x | i >
            assert compare_arrays(deriv1_np[map_key1] + deriv1_np[map_key1].transpose(), deriv1_np[map_key3])

            # Test (S_ij)^x = < i^x | j > + < i | j^x >
            assert compare_arrays(deriv1_np[map_key1] + deriv1_np[map_key2], deriv1_np[map_key3])
