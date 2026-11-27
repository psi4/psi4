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

    # Matrix names from ao_elec_dip_deriv1_helper: must be ao_mu{X,Y,Z}_deriv1_{atom}{X,Y,Z},
    # not the accumulating names a previous bug produced ("ao_muX_deriv1_0X0Y0Z").
    for atom in range(natoms):
        ao_mats = deriv1_mat["MU_" + str(atom)]
        for mu_cart in range(3):
            for atom_cart in range(3):
                expected = "ao_mu" + ["X", "Y", "Z"][mu_cart] + "_deriv1_" + str(atom) + ["X", "Y", "Z"][atom_cart]
                actual = ao_mats[3 * mu_cart + atom_cart].name
                assert actual == expected, f"AO dipole-deriv matrix name mismatch: got {actual!r}, expected {expected!r}"

def test_export_mo_elec_dip_deriv():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    mints = psi4.core.MintsHelper(wfn.basisset())

    natoms = h2o.natom()
    cart = ["X", "Y", "Z"]

    # The square all x all case, and the rectangular occupied x virtual case that goes
    # through the reshaping branch of Matrix::transform.
    for space1, space2 in (("ALL", "ALL"), ("OCC", "VIR")):
        C1 = wfn.Ca_subset("AO", space1)
        C2 = wfn.Ca_subset("AO", space2)
        C1np = np.asarray(C1)
        C2np = np.asarray(C2)

        for atom in range(natoms):
            ao_mats = mints.ao_elec_dip_deriv1(atom)
            mo_mats = mints.mo_elec_dip_deriv1(atom, C1, C2)

            # mo_elec_dip_deriv1 must return 9 matrices indexed as 3*mu_cart + atom_cart.
            # A previous bug looped p in [0,9) and indexed cartcomp[p] (size 3), which is OOB.
            assert len(mo_mats) == 9, f"mo_elec_dip_deriv1 returned {len(mo_mats)} matrices, expected 9"
            for mu_cart in range(3):
                for atom_cart in range(3):
                    p = 3 * mu_cart + atom_cart

                    # Name: mo_mu{X,Y,Z}_deriv1_{atom}{X,Y,Z}
                    expected_name = "mo_mu" + cart[mu_cart] + "_deriv1_" + str(atom) + cart[atom_cart]
                    actual_name = mo_mats[p].name
                    assert actual_name == expected_name, \
                        f"MO dipole-deriv matrix name mismatch: got {actual_name!r}, expected {expected_name!r}"

                    # MO derivative must equal C1^T @ AO derivative @ C2.
                    ao = np.asarray(ao_mats[p])
                    mo = np.asarray(mo_mats[p])
                    assert mo.shape == (C1np.shape[1], C2np.shape[1])
                    label = f"mo_elec_dip_deriv1[{p}] {space1} x {space2}"
                    assert compare_arrays(C1np.T @ ao @ C2np, mo, 10, label), \
                        f"{label} != C1^T @ ao_elec_dip_deriv1[{p}] @ C2"

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


@pytest.mark.parametrize("puream,screening", [(True, "csam"), (False, "csam"), (True, "none")])
def test_contracted_df_deriv1(puream, screening):
    """ao_3center_deriv1_contracted / ao_metric_deriv1_contracted equal the per-atom derivatives contracted with G."""
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
        symmetry c1
    """)
    psi4.set_options({"basis": "cc-pvdz", "puream": puream, "screening": screening})
    orb = psi4.core.BasisSet.build(h2o, "ORBITAL", "cc-pvdz", puream=puream)
    aux = psi4.core.BasisSet.build(h2o, "DF_BASIS_SCF", "cc-pvdz-jkfit", puream=puream)
    mints = psi4.core.MintsHelper(orb)
    mints.set_basisset("DF_BASIS_SCF", aux)
    naux, nbf = aux.nbf(), orb.nbf()
    rng = np.random.default_rng(7)
    G3 = rng.standard_normal((naux, nbf * nbf))      # not symmetric in (m, n)
    G2 = rng.standard_normal((naux, naux))
    ref3 = np.array([[np.vdot(np.asarray(d).reshape(naux, -1), G3) for d in mints.ao_3center_deriv1(a, "DF_BASIS_SCF")]
                     for a in range(h2o.natom())])
    ref2 = np.array([[np.vdot(np.asarray(d), G2) for d in mints.ao_metric_deriv1(a, "DF_BASIS_SCF")]
                     for a in range(h2o.natom())])
    tst3 = mints.ao_3center_deriv1_contracted("DF_BASIS_SCF", psi4.core.Matrix.from_array(G3))
    tst2 = mints.ao_metric_deriv1_contracted("DF_BASIS_SCF", psi4.core.Matrix.from_array(G2))
    assert compare_arrays(ref3, np.asarray(tst3), 10, "contracted 3-center derivative")
    assert compare_arrays(ref2, np.asarray(tst2), 10, "contracted metric derivative")
    with pytest.raises(RuntimeError, match="must be naux x nbf"):
        mints.ao_3center_deriv1_contracted("DF_BASIS_SCF", psi4.core.Matrix.from_array(G2))
