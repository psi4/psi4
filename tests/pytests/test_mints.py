import pytest

import numpy as np
import psi4

from .utils import compare_arrays

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

def test_ang_mom_deriv1():
    h2o = psi4.geometry("""
        O      0.000000000000     0.000000000000     R
        H      0.000000000000    -0.866811832375     0.601435781623
        H      0.000000000000     0.866811832375     0.601435781623
    """)

    psi4.set_options({'basis': 'cc-pVDZ'})

    dz = 0.01
    o_z = -0.075791843897
    #z_list = [o_z - 2 * dz, o_z - dz, o_z, o_z + dz, o_z + 2 * dz]

    energies = dict()
    for l in [-2, -1, 0, 1, 2]:
        h2o.R = o_z + l * dz
        psi4.core.set_active_molecule(h2o)
        energies[l] = psi4.energy('SCF')
    print(energies)

