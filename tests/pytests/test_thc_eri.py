import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.quick]

@pytest.mark.smoke
def test_ls_thc_df():
    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({'basis' : 'cc-pVDZ',
                    'ls_thc_df' : True,
                    'ls_thc_radial_points' : 10,
                    'ls_thc_spherical_points' : 50,
                    'ls_thc_basis_tolerance' : 1.0e-10,
                    'ls_thc_weights_tolerance' : 1.0e-12})

    e, wfn = psi4.energy('scf', return_wfn=True)

    primary = wfn.basisset()
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_MP2", "", "RIFIT", "cc-pVDZ")

    ls_thc_computer = psi4.core.LS_THC_Computer(mol, primary, aux)
    ls_thc_computer.compute_thc_factorization()

    # THC-Factored ERI
    Z_PQ = np.array(ls_thc_computer.get_Z())
    x1 = np.array(ls_thc_computer.get_x1())
    I_guess_thc = np.einsum('pu,pv,pq,qr,qt->uvrt', x1, x1, Z_PQ, x1, x1, optimize=True)

    mints = psi4.core.MintsHelper(primary)
    zero = psi4.core.BasisSet.zero_ao_basis_set()

    # AO 3-index ERIs
    Ppq = np.squeeze(mints.ao_eri(aux, zero, primary, primary))
    J_PQ = np.squeeze(mints.ao_eri(aux, zero, aux, zero))
    J_PQ_inv = np.linalg.pinv(J_PQ)
    I_guess_df = np.einsum('Ppq,PQ,Qrs->pqrs', Ppq, J_PQ_inv, Ppq, optimize=True)

    # AO 4-index ERIs
    I = np.array(mints.ao_eri(primary, primary, primary, primary))

    assert(np.sqrt(np.average(np.square(I_guess_thc-I))) < 3e-4)
    assert(np.sqrt(np.average(np.square(I_guess_df-I))) < 3e-4)

@pytest.mark.smoke
def test_ls_thc_exact():
    mol = psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)

    psi4.set_options({'basis' : 'cc-pVDZ',
                    'ls_thc_df' : False,
                    'ls_thc_radial_points' : 10,
                    'ls_thc_spherical_points' : 50,
                    'ls_thc_basis_tolerance' : 1.0e-10,
                    'ls_thc_weights_tolerance' : 1.0e-12})

    e, wfn = psi4.energy('scf', return_wfn=True)

    primary = wfn.basisset()
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "", "RIFIT", "cc-pVDZ")

    ls_thc_computer = psi4.core.LS_THC_Computer(mol, primary, None)
    ls_thc_computer.compute_thc_factorization()

    # THC-Factored ERI
    Z_PQ = np.array(ls_thc_computer.get_Z())
    x1 = np.array(ls_thc_computer.get_x1())
    I_guess_thc = np.einsum('pu,pv,pq,qr,qt->uvrt', x1, x1, Z_PQ, x1, x1, optimize=True)

    # AO 4-index ERIs
    mints = psi4.core.MintsHelper(primary)
    I = np.array(mints.ao_eri(primary, primary, primary, primary))

    assert(np.sqrt(np.average(np.square(I_guess_thc-I))) < 1e-4)