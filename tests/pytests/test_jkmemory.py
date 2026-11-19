"""
Tests for the memory estimators on JK objects
"""

import psi4
import pytest
from utils import *

pytestmark = [pytest.mark.psi, pytest.mark.api]

def _build_system(basis): 
    mol = psi4.geometry("""
    Ar 0 0  0
    Ar 0 0  5
    Ar 0 0  15
    Ar 0 0  25
    Ar 0 0  35
    """)
    
    #psi4.set_options({"INTS_TOLERANCE": 0.0})
    
    basis = psi4.core.BasisSet.build(mol, target=basis)
    aux = psi4.core.BasisSet.build(basis.molecule(), "DF_BASIS_SCF",
                                   psi4.core.get_option("SCF", "DF_BASIS_SCF"), "JKFIT",
                                   basis.name(), basis.has_puream())
    
    return basis, aux


@pytest.mark.parametrize("basis,jk_type,estimate,name",[
    # Note that MEM_DF and DISK_DF use nthreads in their mem estimate, so ref valid for -n1 only

    # Zero temps
    ["cc-pvdz", "DIRECT",  0, "DirectJK"],
    ["cc-pvdz", "OUT_OF_CORE",  0, "DiskJK"],

    # pvdz tests
    ["cc-pvdz", "MEM_DF",  1590520, "MemDFJK"],
    ["cc-pvdz", "DISK_DF", 1286244, "DiskDFJK"],
    ["cc-pvdz", "CD",      2916000, "CDJK"],
    ["cc-pvdz", "PK",      65610000, "PKJK"],

    # 5z tests
    ["cc-pv5z", "MEM_DF",  57020770, "MemDFJK"],
    ["cc-pv5z", "DISK_DF", 26984120, "DiskDFJK"],
]) # yapf: disable

def test_jk_memory_estimate(basis, jk_type, estimate, name):

    basis, aux = _build_system(basis)
    jk = psi4.core.JK.build(basis, aux=aux, jk_type=jk_type, do_wK=False, memory=1e9)

    assert compare(name, jk.name(), f"{jk_type} algorithm name")
    assert compare_integers(estimate, jk.memory_estimate(), "{} memory estimate".format(jk_type))



def test_diskdfjk_core_blocked_k():
    """In-core DiskDFJK whose budget is lowered after initialize builds K in many row blocks; J and K are unchanged."""
    import numpy as np

    basis, aux = _build_system("cc-pvdz")
    nbf, naux, nocc = basis.nbf(), aux.nbf(), 9
    C, _ = np.linalg.qr(np.random.default_rng(5).standard_normal((nbf, nocc)))
    C = psi4.core.Matrix.from_array(C)

    def jk_mats(rows_budget):
        jk = psi4.core.JK.build(basis, aux=aux, jk_type="DISK_DF", do_wK=False, memory=1e9)
        jk.initialize()                                   # in core: (Q|mn) resident
        if rows_budget is not None:                       # resident tensor + overheads + a few E rows
            jk.set_memory(int(jk.memory_estimate() - 2 * naux * naux + rows_budget * nocc * nbf))
        jk.C_left_add(C)
        jk.compute()
        return np.asarray(jk.J()[0]).copy(), np.asarray(jk.K()[0]).copy()

    J_ref, K_ref = jk_mats(None)
    J_blk, K_blk = jk_mats(3)
    assert compare_values(J_ref, J_blk, 10, "DiskDFJK blocked J")
    assert compare_values(K_ref, K_blk, 10, "DiskDFJK blocked K")
