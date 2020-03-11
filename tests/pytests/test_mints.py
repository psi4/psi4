import pytest

import numpy as np
import psi4

from .utils import compare_values

def test_export_ao_elec_dip_deriv():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    psi4.core.set_active_molecule(h2o)
    psi4.set_options({'basis': 'cc-pVDZ'})

    rhf_e, wfn = psi4.energy('SCF', return_wfn=True)

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
    psi4.compare_matrices(PSI4_MU_Grad, G_python_MU_mat, 10, "DIPOLE_GRADIENT_TEST")  # TEST
