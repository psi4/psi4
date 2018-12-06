import numpy as np
import pytest

import psi4
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine,
                                                           TDUSCFEngine)
from utils import compare_arrays


def build_RHF_AB_C1_singlet(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    V_iajb = mints.mo_eri(Co, Cv, Co, Cv).to_array()
    V_abij = mints.mo_eri(Cv, Cv, Co, Co).to_array()
    Fab = psi4.core.Matrix.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.Matrix.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
    ni = Fij.shape[0]
    na = Fab.shape[0]
    nia = ni * na
    A_ref = np.einsum("ab,ij->iajb", Fab, np.eye(ni))
    A_ref -= np.einsum("ab,ij->iajb", np.eye(na), Fij)
    A_ref += 2 * V_iajb - np.einsum("abij->iajb", V_abij)
    B_ref = 2 * V_iajb - V_iajb.swapaxes(0, 2)
    return A_ref, B_ref


def build_RHF_AB_C1_triplet(wfn):
    mints = psi4.core.MintsHelper(wfn.basisset())
    Co = wfn.Ca_subset("SO", "OCC")
    Cv = wfn.Ca_subset("SO", "VIR")
    V_iajb = mints.mo_eri(Co, Cv, Co, Cv).to_array()
    V_abij = mints.mo_eri(Cv, Cv, Co, Co).to_array()
    Fab = psi4.core.Matrix.triplet(Cv, wfn.Fa(), Cv, True, False, False).to_array()
    Fij = psi4.core.Matrix.triplet(Co, wfn.Fa(), Co, True, False, False).to_array()
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
    FAB = psi4.core.Matrix.triplet(CA, wfn.Fa(), CA, True, False, False).to_array()
    FIJ = psi4.core.Matrix.triplet(CI, wfn.Fa(), CI, True, False, False).to_array()
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
    Fab = psi4.core.Matrix.triplet(Ca, wfn.Fb(), Ca, True, False, False).to_array()
    Fij = psi4.core.Matrix.triplet(Ci, wfn.Fb(), Ci, True, False, False).to_array()
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


@pytest.mark.tddft
def test_restricted_TDA_singlet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk"})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, _ = build_RHF_AB_C1_singlet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=False)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    A_test = np.column_stack([x.to_array().flatten() for x in eng.compute_products(ID)])
    assert compare_arrays(A_ref, A_test, 8, "RHF Ax C1 products")


@pytest.mark.tddft
def test_restricted_TDA_triplet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk"})
    e, wfn = psi4.energy("hf/cc-pvdz", molecule=h2o, return_wfn=True)
    A_ref, _ = build_RHF_AB_C1_triplet(wfn)
    ni, na, _, _ = A_ref.shape
    nia = ni * na
    A_ref = A_ref.reshape((nia, nia))
    # Build engine
    eng = TDRSCFEngine(wfn, ptype='tda', triplet=True)
    # our "guess"" vectors
    ID = [psi4.core.Matrix.from_array(v.reshape((ni, na))) for v in tuple(np.eye(nia).T)]
    A_test = np.column_stack([x.to_array().flatten() for x in eng.compute_products(ID)])
    assert compare_arrays(A_ref, A_test, 8, "RHF Ax C1 products")


@pytest.mark.tddft
def test_restricted_RPA_singlet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk"})
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
    Px, Mx = eng.compute_products(ID)
    P_test = np.column_stack([x.to_array().flatten() for x in Px])
    assert compare_arrays(P_ref, P_test, 8, "RHF (A+B)x C1 products")
    M_test = np.column_stack([x.to_array().flatten() for x in Mx])
    assert compare_arrays(M_ref, M_test, 8, "RHF (A-B)x C1 products")


@pytest.mark.tddft
def test_restricted_RPA_triplet_c1():
    "Build out the full CIS/TDA hamiltonian (A) col by col with the product engine"
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk"})
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
    Px, Mx = eng.compute_products(ID)
    P_test = np.column_stack([x.to_array().flatten() for x in Px])
    assert compare_arrays(P_ref, P_test, 8, "RHF (A+B)x C1 products")
    M_test = np.column_stack([x.to_array().flatten() for x in Mx])
    assert compare_arrays(M_ref, M_test, 8, "RHF (A-B)x C1 products")


@pytest.mark.tddft
def test_unrestricted_TDA_C1():
    ch2 = psi4.geometry("""
    0 3
    c
    h 1 1.0
    h 1 1.0 2 125.0
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'reference': 'UHF'})
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
    X_I0 = list(zip(X_JB, zero_jb))
    X_0I = list(zip(zero_JB, X_jb))
    Ax_I0 = eng.compute_products(X_I0)
    Ax_0I = eng.compute_products(X_0I)
    A_IAJB_test = np.column_stack([x[0].to_array().flatten() for x in Ax_I0])
    assert compare_arrays(A_ref['IAJB'].reshape(nIA, nIA), A_IAJB_test, 8, "A_IAJB")
    A_iaJB_test = np.column_stack([x[1].to_array().flatten() for x in Ax_I0])
    assert compare_arrays(A_ref['iaJB'].reshape(nia, nIA), A_iaJB_test, 8, "A_iaJB")
    A_IAjb_test = np.column_stack([x[0].to_array().flatten() for x in Ax_0I])
    assert compare_arrays(A_ref['IAjb'].reshape(nIA, nia), A_IAjb_test, 8, "A_IAjb")
    A_iajb_test = np.column_stack([x[1].to_array().flatten() for x in Ax_0I])
    assert compare_arrays(A_ref['iajb'].reshape(nia, nia), A_iajb_test, 8, "A_iajb")


@pytest.mark.tddft
def test_unrestricted_RPA_C1():
    ch2 = psi4.geometry("""
    0 3
    c
    h 1 1.0
    h 1 1.0 2 125.0
    symmetry c1
    """)
    psi4.set_options({"scf_type": "pk", 'reference': 'UHF'})
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
    X_I0 = list(zip(X_JB, zero_jb))
    X_0I = list(zip(zero_JB, X_jb))
    Px_I0, Mx_I0 = eng.compute_products(X_I0)
    Px_0I, Mx_0I = eng.compute_products(X_0I)

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
