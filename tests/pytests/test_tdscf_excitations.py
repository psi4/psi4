from pathlib import Path
import json

import numpy as np
import pytest

import psi4
from psi4.driver.p4util.solvers import davidson_solver, hamiltonian_solver
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine, TDUSCFEngine)
from .utils import *

## marks
# reference type
UHF = pytest.mark.unrestricted
RHF_singlet = pytest.mark.restricted_singlet
RHF_triplet = pytest.mark.restricted_triplet
# functional types
hf = pytest.mark.hf
lda = pytest.mark.lda
gga = pytest.mark.gga
hyb_gga = pytest.mark.hyb_gga
hyb_gga_lrc = pytest.mark.hyb_gga_lrc
# response type
RPA = pytest.mark.RPA
TDA = pytest.mark.TDA

# Reference data generated using Gaussian09
with open(Path(__file__).parent / "tdscf_reference_data.json") as f:
    reference_data = json.load(f)


@pytest.fixture
def tddft_systems():
    psi4.core.clean()

    # Canonical unrestricted system
    ch2 = psi4.geometry("""0 3
    C           0.000000    0.000000    0.159693
    H          -0.000000    0.895527   -0.479080
    H          -0.000000   -0.895527   -0.479080
    no_reorient
    no_com
    """)

    # Canonical restricted system
    h2o = psi4.geometry("""0 1
    O           0.000000    0.000000    0.135446
    H          -0.000000    0.866812   -0.541782
    H          -0.000000   -0.866812   -0.541782
    no_reorient
    no_com
    """)

    return {'UHF': ch2, 'RHF': h2o}


@pytest.fixture
def wfn_factory(tddft_systems):
    def _build_wfn(ref, func, basis, nosym):
        if ref.startswith('RHF'):
            mol = tddft_systems['RHF']
        else:
            mol = tddft_systems['UHF']
            psi4.set_options({'reference': 'UHF'})
        if nosym:
            mol.reset_point_group('c1')
        psi4.set_options({'scf_type': 'pk', 'e_convergence': 8, 'd_convergence': 8, 'save_jk': True})
        e, wfn = psi4.energy(f"{func}/{basis}", return_wfn=True, molecule=mol)
        return wfn

    return _build_wfn


@pytest.fixture
def solver_funcs():
    return {'TDA': davidson_solver, 'RPA': hamiltonian_solver}


@pytest.fixture
def engines():
    return {
        'RHF-1': lambda w, p: TDRSCFEngine(w, ptype=p.lower(), triplet=False),
        'RHF-3': lambda w, p: TDRSCFEngine(w, ptype=p.lower(), triplet=True),
        'UHF': lambda w, p: TDUSCFEngine(w, ptype=p.lower())
    }


@pytest.mark.tdscf
@pytest.mark.parametrize("ref,func,ptype,basis", [
    pytest.param(   'UHF',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',        'HF',  'RPA',  'cc-pvdz', marks=[hf, UHF, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param(   'UHF',        'HF',  'TDA',  'cc-pvdz', marks=[hf, UHF, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param(   'UHF',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]), # G09 rev E.01
]) # yapf: disable
def test_tdscf(ref, func, ptype, basis, wfn_factory, solver_funcs, engines):
    if (ref == 'RHF-1') or (func == "HF"):
        # RHF-singlet everything works and TDHF/CIS works for RHF-triplet, UHF
        pass
    elif (ref == 'RHF-3'):
        pytest.xfail("RKS Vx kernel only Spin Adapted for Singlet")
    elif (ref == 'UHF' and func != 'SVWN'):
        pytest.xfail("UKS Vx kernel bug for non-lda")

    # get wfn (don't use symmetry b/c too slow)
    wfn = wfn_factory(ref, func, basis, nosym=True)
    # select solver function (TDA->davidson/RPA->hamiltonian)
    solver = solver_funcs[ptype]

    # build engine
    engine = engines[ref](wfn, ptype)

    # skipping the entrypoint, just call the solver
    out = solver(engine=engine,
                 guess=engine.generate_guess(16),
                 max_vecs_per_root=10,
                 nroot=4,
                 verbose=1,
                 maxiter=30,
                 r_tol=1.0e-5,
                 schmidt_tol=1.0e-12)

    test_vals = out["eigvals"]
    stats = out["stats"]
    assert stats[-1]['done'], "Solver did not converge"

    mol = "CH2" if ref == "UHF" else "H2O"
    ref_v = reference_data[f"{mol}_{ref}_{func}_{ptype}"]

    for i, my_v in enumerate(test_vals):
        assert compare_values(ref_v[i]["EXCITATION ENERGY"], my_v, 4, f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1}")
