import numpy as np
import pytest

import psi4
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine,
                                                           TDUSCFEngine)
from utils import compare_arrays, compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

def build_RHF_AB_C1_singlet(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    V_iajb = mints.mo_eri(Co, Cv, Co, Cv).to_array()
    V_abij = mints.mo_eri(Cv, Cv, Co, Co).to_array()
    Fab = psi4.core.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    ni = Fij.shape[0]
    na = Fab.shape[0]
    nia = ni * na
    A_ref = np.einsum("ab,ij->iajb", Fab, np.eye(ni))
    A_ref -= np.einsum("ab,ij->iajb", np.eye(na), Fij)
    A_ref += 2 * V_iajb - np.einsum("abij->iajb", V_abij)
    B_ref = 2 * V_iajb - V_iajb.swapaxes(0, 2)
    return A_ref, B_ref


def build_RHF_AB_singlet(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    Fab = psi4.core.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    # mo_eri can't handle systems with symmetry. We need to work around this.
    ao_eri = mints.ao_eri()
    ao2so = wfn.aotoso()

    # The h'th irrep stores the block where ia has symmetry h.
    # The elements are indexed by ov pairs. Elements are in ascending order
    # of the occupied element of the pair.
    A_blocks = []
    B_blocks = []

    for hjb in range(wfn.nirrep()):
        hia = hjb
        A_block = []
        B_block = []
        for hi in range(wfn.nirrep()):
            A_block.append([])
            B_block.append([])
            ha = hia ^ hi
            Ca = np.matmul(ao2so.nph[ha], Cv.nph[ha])
            Ci = np.matmul(ao2so.nph[hi], Co.nph[hi])
            for hj in range(wfn.nirrep()):
                hb = hjb ^ hj
                Cb = np.matmul(ao2so.nph[hb], Cv.nph[hb])
                Cj = np.matmul(ao2so.nph[hj], Co.nph[hj])
                V_iajb = np.einsum("pqrs, pP, qQ, rR, sS -> PQRS", ao_eri, Ci, Ca, Cj, Cb, optimize=True)
                V_jaib = np.einsum("pqrs, pP, qQ, rR, sS -> PQRS", ao_eri, Cj, Ca, Ci, Cb, optimize=True)
                V_abij = np.einsum("pqrs, pP, qQ, rR, sS -> PQRS", ao_eri, Ca, Cb, Ci, Cj, optimize=True)
                A_ref =  2 * V_iajb - np.einsum("abij->iajb", V_abij)
                if ha == hb and hi == hj:
                    A_ref += np.einsum("ab,ij->iajb", Fab[ha], np.eye(Fij[hi].shape[0]), optimize=True)
                    A_ref -= np.einsum("ij,ab->iajb", Fij[hi], np.eye(Fab[ha].shape[0]), optimize=True)
                B_ref = 2 * V_iajb - V_jaib.swapaxes(0, 2)
                shape_tuple = (A_ref.shape[0] * A_ref.shape[1], A_ref.shape[2] * A_ref.shape[3])
                A_block[-1].append(A_ref.reshape(shape_tuple))
                B_block[-1].append(B_ref.reshape(shape_tuple))
        A_blocks.append(np.block(A_block))
        B_blocks.append(np.block(B_block))
    return A_blocks, B_blocks


def build_RHF_AB_C1_singlet_df(wfn):
    orb = wfn.get_basisset("ORBITAL")
    mints = psi4.core.MintsHelper(orb)
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
    aux = wfn.get_basisset("DF_BASIS_SCF")
    Ppq = np.squeeze(mints.ao_eri(aux, zero_bas, orb, orb))
    metric = mints.ao_eri(aux, zero_bas, aux, zero_bas)
    metric.power(-0.5, 1.e-14)
    metric = np.squeeze(metric)
    Qpq = np.einsum("QP,Ppq->Qpq", metric, Ppq, optimize=True)
    Qij = np.einsum("Qpq, pi, qj -> Qij", Qpq, Co, Co)
    Qab = np.einsum("Qpq, pa, qb -> Qab", Qpq, Cv, Cv)
    Qia = np.einsum("Qpq, pi, qa -> Qia", Qpq, Co, Cv)
    V_iajb = np.einsum("Qia, Qjb -> iajb", Qia, Qia)
    V_abij = np.einsum("Qab, Qij -> abij", Qab, Qij)
    Fab = psi4.core.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    ni = Fij.shape[0]
    na = Fab.shape[0]
    nia = ni * na
    A_ref = np.einsum("ab,ij->iajb", Fab, np.eye(ni))
    A_ref -= np.einsum("ab,ij->iajb", np.eye(na), Fij)
    A_ref += 2 * V_iajb - np.einsum("abij->iajb", V_abij)
    B_ref = 2 * V_iajb - V_iajb.swapaxes(0, 2)
    return A_ref, B_ref


def build_RHF_AB_singlet_df(wfn):
    orb = wfn.get_basisset("ORBITAL")
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    Fab = psi4.core.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
    aux = wfn.get_basisset("DF_BASIS_SCF")
    Ppq = np.squeeze(mints.ao_eri(aux, zero_bas, orb, orb))
    metric = mints.ao_eri(aux, zero_bas, aux, zero_bas)
    metric.power(-0.5, 1.e-14)
    metric = np.squeeze(metric)
    Qpq = np.einsum("QP,Ppq->Qpq", metric, Ppq, optimize=True)
    ao2so = wfn.aotoso()

    # The h'th irrep stores the block where ia has symmetry h.
    # The elements are indexed by ov pairs. Elements are in ascending order
    # of the occupied element of the pair.
    A_blocks = []
    B_blocks = []

    for hjb in range(wfn.nirrep()):
        hia = hjb
        A_block = []
        B_block = []
        for hi in range(wfn.nirrep()):
            A_block.append([])
            B_block.append([])
            ha = hia ^ hi
            Ca = np.matmul(ao2so.nph[ha], Cv.nph[ha])
            Ci = np.matmul(ao2so.nph[hi], Co.nph[hi])
            for hj in range(wfn.nirrep()):
                hb = hjb ^ hj
                Cb = np.matmul(ao2so.nph[hb], Cv.nph[hb])
                Cj = np.matmul(ao2so.nph[hj], Co.nph[hj])
                Qij = np.einsum("Ppq, pi, qj -> Pij", Qpq, Ci, Cj, optimize=True)
                Qab = np.einsum("Ppq, pa, qb -> Pab", Qpq, Ca, Cb, optimize=True)
                Qia = np.einsum("Ppq, pi, qa -> Pia", Qpq, Ci, Ca, optimize=True)
                Qjb = np.einsum("Ppq, pj, qb -> Pjb", Qpq, Cj, Cb, optimize=True)
                Qja = np.einsum("Ppq, pj, qa -> Pja", Qpq, Cj, Ca, optimize=True)
                Qib = np.einsum("Ppq, pi, qb -> Pib", Qpq, Ci, Cb, optimize=True)
                V_iajb = np.einsum("Pia, Pjb -> iajb", Qia, Qjb, optimize=True)
                V_jaib = np.einsum("Pia, Pjb -> iajb", Qja, Qib, optimize=True)
                V_abij = np.einsum("Pij, Pab -> abij", Qij, Qab, optimize=True)
                A_ref =  2 * V_iajb - np.einsum("abij->iajb", V_abij)
                if ha == hb and hi == hj:
                    A_ref += np.einsum("ab,ij->iajb", Fab[ha], np.eye(Fij[hi].shape[0]), optimize=True)
                    A_ref -= np.einsum("ij,ab->iajb", Fij[hi], np.eye(Fab[ha].shape[0]), optimize=True)
                B_ref = 2 * V_iajb - V_jaib.swapaxes(0, 2)
                shape_tuple = (A_ref.shape[0] * A_ref.shape[1], A_ref.shape[2] * A_ref.shape[3])
                A_block[-1].append(A_ref.reshape(shape_tuple))
                B_block[-1].append(B_ref.reshape(shape_tuple))
        A_blocks.append(np.block(A_block))
        B_blocks.append(np.block(B_block))
    return A_blocks, B_blocks


def build_RHF_AB_C1_triplet(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    V_iajb = mints.mo_eri(Co, Cv, Co, Cv).to_array()
    V_abij = mints.mo_eri(Cv, Cv, Co, Co).to_array()
    Fab = psi4.core.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    ni = Fij.shape[0]
    na = Fab.shape[0]
    nia = ni * na
    A_ref = np.einsum("ab,ij->iajb", Fab, np.eye(ni))
    A_ref -= np.einsum("ab,ij->iajb", np.eye(na), Fij)
    A_ref -= np.einsum("abij->iajb", V_abij)
    B_ref = -V_iajb.swapaxes(0, 2)
    return A_ref, B_ref


def build_UHF_AB_C1(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    CI = wfn.Ca_subset("SO", "OCC")
    CA = wfn.Ca_subset("SO", "VIR")
    V_IAJB = mints.mo_eri(CI, CA, CI, CA).to_array()
    V_ABIJ = mints.mo_eri(CA, CA, CI, CI).to_array()
    FAB = psi4.core.triplet(CA, wfn.Fa(), CA, True, False, False).to_array()
    FIJ = psi4.core.triplet(CI, wfn.Fa(), CI, True, False, False).to_array()
    nI = FIJ.shape[0]
    nA = FAB.shape[0]
    nIA = nI * nA

    A = {}
    B = {}

    A['IAJB'] = np.einsum("AB,IJ->IAJB", FAB, np.eye(nI))
    A['IAJB'] -= np.einsum("AB,IJ->IAJB", np.eye(nA), FIJ)
    A['IAJB'] += V_IAJB
    A['IAJB'] -= np.einsum("ABIJ->IAJB", V_ABIJ)

    B['IAJB'] = V_IAJB - V_IAJB.swapaxes(0, 2)

    Ci = wfn.Cb_subset("SO", "OCC")
    Ca = wfn.Cb_subset("SO", "VIR")
    V_iajb = mints.mo_eri(Ci, Ca, Ci, Ca).to_array()
    V_abij = mints.mo_eri(Ca, Ca, Ci, Ci).to_array()
    Fab = psi4.core.triplet(Ca, wfn.Fb(), Ca, True, False, False).to_array()
    Fij = psi4.core.triplet(Ci, wfn.Fb(), Ci, True, False, False).to_array()
    ni = Fij.shape[0]
    na = Fab.shape[0]
    nia = ni * na

    A['iajb'] = np.einsum("ab,ij->iajb", Fab, np.eye(ni))
    A['iajb'] -= np.einsum("ab,ij->iajb", np.eye(na), Fij)
    A['iajb'] += V_iajb
    A['iajb'] -= np.einsum('abij->iajb', V_abij)

    B['iajb'] = V_iajb - V_iajb.swapaxes(0, 2)

    V_IAjb = mints.mo_eri(CI, CA, Ci, Ca).to_array()
    V_iaJB = mints.mo_eri(Ci, Ca, CI, CA).to_array()
    A['IAjb'] = V_IAjb
    A['iaJB'] = V_iaJB

    B['IAjb'] = V_IAjb
    B['iaJB'] = V_iaJB
    return A, B


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_TDA_singlet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, _ = build_RHF_AB_C1_singlet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    A_test = np.column_stack([x.to_array().flatten() for x in eng.compute_products(ID)[0]])
    assert compare_arrays(A_ref, A_test, 8, "RHF Ax C1 products")

@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_TDA_singlet():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_blocks, B_blocks = build_RHF_AB_singlet(wfn)
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    vir_dim = wfn.nmopi() - wfn.doccpi()
    for hia, A_block in enumerate(A_blocks):
        ID = []
        # Construct a matrix for each (O, V) pair with hia symmetry.
        for hi in range(wfn.nirrep()):
            for i in range(wfn.Ca_subset("SO", "OCC").coldim()[hi]):
                for a in range(wfn.Ca_subset("SO", "VIR").coldim()[hi ^ hia]):
                    matrix = psi4.core.Matrix("Test Matrix", wfn.doccpi(), vir_dim, hia)
                    matrix.set(hi, i, a, 1)
                    ID.append(matrix)
        x = eng.compute_products(ID)[0][0]
        # Assemble the A values as a single (ia, jb) matrix, all possible ia and jb of symmetry hia.
        A_test = np.column_stack([np.concatenate([y.flatten() for y in x.to_array()]) for x in eng.compute_products(ID)[0]])
        assert compare_arrays(A_block, A_test, 8, "RHF Ax C2v products")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_TDA_singlet_df_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "df", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, _ = build_RHF_AB_C1_singlet_df(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    A_test = np.column_stack([x.to_array().flatten() for x in eng.compute_products(ID)[0]])
    assert compare_arrays(A_ref, A_test, 8, "DF-RHF Ax C1 products")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_TDA_singlet_df():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)
    psi4.set_options({"scf_type": "df", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_blocks, B_blocks = build_RHF_AB_singlet_df(wfn)
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    vir_dim = wfn.nmopi() - wfn.doccpi()
    for hia, A_block in enumerate(A_blocks):
        ID = []
        # Construct a matrix for each (O, V) pair with hia symmetry.
        for hi in range(wfn.nirrep()):
            for i in range(wfn.Ca_subset("SO", "OCC").coldim()[hi]):
                for a in range(wfn.Ca_subset("SO", "VIR").coldim()[hi ^ hia]):
                    matrix = psi4.core.Matrix("Test Matrix", wfn.doccpi(), vir_dim, hia)
                    matrix.set(hi, i, a, 1)
                    ID.append(matrix)
        x = eng.compute_products(ID)[0][0]
        # Assemble the A values as a single (ia, jb) matrix, all possible ia and jb of symmetry hia.
        A_test = np.column_stack([np.concatenate([y.flatten() for y in x.to_array()]) for x in eng.compute_products(ID)[0]])
        assert compare_arrays(A_block, A_test, 8, "DF-RHF Ax C2v products")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_TDA_triplet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, _ = build_RHF_AB_C1_triplet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=True)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    A_test = np.column_stack([x.to_array().flatten() for x in eng.compute_products(ID)[0]])
    assert compare_arrays(A_ref, A_test, 8, "RHF Ax C1 products")


@pytest.mark.unittest
@pytest.mark.tdscf
@pytest.mark.check_triplet
def test_RU_TDA_C1():
    h2o = psi4.geometry("""0 1
    O          0.000000    0.000000    0.135446
    H         -0.000000    0.866812   -0.541782
    H         -0.000000   -0.866812   -0.541782
    symmetry c1
    no_reorient
    no_com
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/sto-3g", molecule=h2o, return_wfn=True)
    A_ref, _ = build_UHF_AB_C1(wfn)
    ni, na, _, _ = A_ref['IAJB'].shape
    nia = ni * na
    A_sing_ref = A_ref['IAJB'] + A_ref['IAjb']
    A_sing_ref = A_sing_ref.reshape(nia, nia)
    A_trip_ref = A_ref['IAJB'] - A_ref['IAjb']
    A_trip_ref = A_trip_ref.reshape(nia, nia)
    sing_vals, _ = np.linalg.eigh(A_sing_ref)
    trip_vals, _ = np.linalg.eigh(A_trip_ref)

    trip_eng = TDRSCFEngine(wfn, ptype='tda', triplet=True)
    sing_eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    psi4.core.print_out("\nA sing:\n" + str(A_sing_ref) + "\n\n")
    psi4.core.print_out("\nA trip:\n" + str(A_trip_ref) + "\n\n")
    A_trip_test = np.column_stack([x.to_array().flatten() for x in trip_eng.compute_products(ID)[0]])
    assert compare_arrays(A_trip_ref, A_trip_test, 8, "Triplet Ax C1 products")
    A_sing_test = np.column_stack([x.to_array().flatten() for x in sing_eng.compute_products(ID)[0]])
    assert compare_arrays(A_sing_ref, A_sing_test, 8, "Singlet Ax C1 products")

    sing_vals_2, _ = np.linalg.eigh(A_sing_test)
    trip_vals_2, _ = np.linalg.eigh(A_trip_test)

    psi4.core.print_out("\n\n SINGLET EIGENVALUES\n")
    for x, y in zip(sing_vals, sing_vals_2):
        psi4.core.print_out("{:10.6f}  {:10.6f}\n".format(x, y))
        # assert compare_values(x, y, 4, "Singlet ROOT")
    psi4.core.print_out("\n\n Triplet EIGENVALUES\n")
    for x, y in zip(trip_vals, trip_vals_2):
        psi4.core.print_out("{:10.6f}  {:10.6f}\n".format(x, y))
        # assert compare_values(x, y, 4, "Triplet Root")

    for x, y in zip(sing_vals, sing_vals_2):
        assert compare_values(x, y, 4, "Singlet ROOT")
    for x, y in zip(trip_vals, trip_vals_2):
        assert compare_values(x, y, 4, "Triplet Root")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_RPA_singlet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, B_ref = build_RHF_AB_C1_singlet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    B_ref = B_ref.reshape((nia, nia))
    P_ref = A_ref + B_ref
    M_ref = A_ref - B_ref
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='rpa', triplet=False)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    Px, Mx = eng.compute_products(ID)[:-1]
    P_test = np.column_stack([x.to_array().flatten() for x in Px])
    assert compare_arrays(P_ref, P_test, 8, "RHF (A+B)x C1 products")
    M_test = np.column_stack([x.to_array().flatten() for x in Mx])
    assert compare_arrays(M_ref, M_test, 8, "RHF (A-B)x C1 products")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_restricted_RPA_triplet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, B_ref = build_RHF_AB_C1_triplet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    B_ref = B_ref.reshape((nia, nia))
    P_ref = A_ref + B_ref
    M_ref = A_ref - B_ref
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='rpa', triplet=True)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    Px, Mx = eng.compute_products(ID)[:-1]
    P_test = np.column_stack([x.to_array().flatten() for x in Px])
    assert compare_arrays(P_ref, P_test, 8, "RHF (A+B)x C1 products")
    M_test = np.column_stack([x.to_array().flatten() for x in Mx])
    assert compare_arrays(M_ref, M_test, 8, "RHF (A-B)x C1 products")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_unrestricted_TDA_C1():
    ch2 = psi4.geometry("""
    0 3
    c
    h 1 1.0
    h 1 1.0 2 125.0
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'reference': 'UHF', 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=ch2, return_wfn=True)
    A_ref, B_ref = build_UHF_AB_C1(wfn)
    nI, nA, _, _ = A_ref['IAJB'].shape
    nIA = nI * nA
    ni, na, _, _ = A_ref['iajb'].shape
    nia = ni * na

    eng = TDUSCFEngine(wfn, ptype='tda')
    X_jb = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    zero_jb = [psi4.core.Matrix(ni, na) for x in range(nIA)]
    X_JB = [psi4.core.Matrix.from_array(v.reshape((nI, nA))) for v in tuple(np.eye(nIA).T)]
    zero_JB = [psi4.core.Matrix(nI, nA) for x in range(nia)]
    # Guess Identity:
    #     X_I0          X_0I          =      X_I0      X_0I
    # [ I{nOV x nOV} | 0{nOV x nov}]  = [ X{KC,JB} | 0{KC, jb}]
    # [ 0{nov x nOV} | I{nov x nov}]    [ 0{kc,JB} | X{kc, jb}]

    # Products:
    # [ A{IA, KC}  A{IA, kc}] [ I{KC, JB} |  0{KC,jb}] = [A x X_I0] = [ A_IAJB, A_iaJB]
    # [ A{ia, KC}  A{ia, kc}] [ O{kc, JB} |  X{kc,jb}]   [A x X_0I] = [ A_IAjb, A_iajb]
    X_I0 = [[x, zero] for x, zero in zip(X_JB, zero_jb)]
    X_0I = [[zero, x] for zero, x in zip(zero_JB, X_jb)]
    Ax_I0 = eng.compute_products(X_I0)[0]
    Ax_0I = eng.compute_products(X_0I)[0]
    A_IAJB_test = np.column_stack([x[0].to_array().flatten() for x in Ax_I0])
    assert compare_arrays(A_ref['IAJB'].reshape(nIA, nIA), A_IAJB_test, 8, "A_IAJB")
    A_iaJB_test = np.column_stack([x[1].to_array().flatten() for x in Ax_I0])
    assert compare_arrays(A_ref['iaJB'].reshape(nia, nIA), A_iaJB_test, 8, "A_iaJB")
    A_IAjb_test = np.column_stack([x[0].to_array().flatten() for x in Ax_0I])
    assert compare_arrays(A_ref['IAjb'].reshape(nIA, nia), A_IAjb_test, 8, "A_IAjb")
    A_iajb_test = np.column_stack([x[1].to_array().flatten() for x in Ax_0I])
    assert compare_arrays(A_ref['iajb'].reshape(nia, nia), A_iajb_test, 8, "A_iajb")


@pytest.mark.unittest
@pytest.mark.tdscf
def test_unrestricted_RPA_C1():
    ch2 = psi4.geometry("""
    0 3
    c
    h 1 1.0
    h 1 1.0 2 125.0
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'reference': 'UHF', 'save_jk': True})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=ch2, return_wfn=True)
    A_ref, B_ref = build_UHF_AB_C1(wfn)
    nI, nA, _, _ = A_ref['IAJB'].shape
    nIA = nI * nA
    ni, na, _, _ = A_ref['iajb'].shape
    nia = ni * na

    P_ref = {k: A_ref[k] + B_ref[k] for k in A_ref.keys()}
    M_ref = {k: A_ref[k] - B_ref[k] for k in A_ref.keys()}

    eng = TDUSCFEngine(wfn, ptype='rpa')
    X_jb = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    zero_jb = [psi4.core.Matrix(ni, na) for x in range(nIA)]
    X_JB = [psi4.core.Matrix.from_array(v.reshape((nI, nA))) for v in tuple(np.eye(nIA).T)]
    zero_JB = [psi4.core.Matrix(nI, nA) for x in range(nia)]
    # Guess Identity:
    #     X_I0          X_0I          =      X_I0      X_0I
    # [ I{nOV x nOV} | 0{nOV x nov}]  = [ X{KC,JB} | 0{KC, jb}]
    # [ 0{nov x nOV} | I{nov x nov}]    [ 0{kc,JB} | X{kc, jb}]

    # Products:
    # [ A+/-B{IA, KC}  A+/-B{IA, kc}] [ I{KC, JB} |  0{KC,jb}] = [A+/-B x X_I0] = [ (A+/-B)_IAJB, (A+/-B)_iaJB]
    # [ A+/-B{ia, KC}  A+/-B{ia, kc}] [ O{kc, JB} |  X{kc,jb}]   [A+/-B x X_0I] = [ (A+/-B)_IAjb, (A+/-B)_iajb]
    X_I0 = [[x, zero] for x, zero in zip(X_JB, zero_jb)]
    X_0I = [[zero, x] for zero, x in zip(zero_JB, X_jb)]
    Px_I0, Mx_I0 = eng.compute_products(X_I0)[:-1]
    Px_0I, Mx_0I = eng.compute_products(X_0I)[:-1]

    P_IAJB_test = np.column_stack([x[0].to_array().flatten() for x in Px_I0])
    assert compare_arrays(P_ref['IAJB'].reshape(nIA, nIA), P_IAJB_test, 8, "A_IAJB")

    M_IAJB_test = np.column_stack([x[0].to_array().flatten() for x in Mx_I0])
    assert compare_arrays(M_ref['IAJB'].reshape(nIA, nIA), M_IAJB_test, 8, "A_IAJB")

    P_iaJB_test = np.column_stack([x[1].to_array().flatten() for x in Px_I0])
    assert compare_arrays(P_ref['iaJB'].reshape(nia, nIA), P_iaJB_test, 8, "P_iaJB")

    M_iaJB_test = np.column_stack([x[1].to_array().flatten() for x in Mx_I0])
    assert compare_arrays(M_ref['iaJB'].reshape(nia, nIA), M_iaJB_test, 8, "M_iaJB")

    P_IAjb_test = np.column_stack([x[0].to_array().flatten() for x in Px_0I])
    assert compare_arrays(P_ref['IAjb'].reshape(nIA, nia), P_IAjb_test, 8, "P_IAjb")

    M_IAjb_test = np.column_stack([x[0].to_array().flatten() for x in Mx_0I])
    assert compare_arrays(M_ref['IAjb'].reshape(nIA, nia), M_IAjb_test, 8, "M_IAjb")

    P_iajb_test = np.column_stack([x[1].to_array().flatten() for x in Px_0I])
    assert compare_arrays(P_ref['iajb'].reshape(nia, nia), P_iajb_test, 8, "P_iajb")

    M_iajb_test = np.column_stack([x[1].to_array().flatten() for x in Mx_0I])
    assert compare_arrays(M_ref['iajb'].reshape(nia, nia), M_iajb_test, 8, "M_iajb")
