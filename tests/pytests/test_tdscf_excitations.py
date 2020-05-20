from pathlib import Path
import json

import numpy as np
import pytest

import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
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

# Reference data generated using G09.E01
with open(Path(__file__).parent / "tdscf_reference_data.json") as f:
    reference_data = json.load(f)

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


@pytest.mark.tdscf
@pytest.mark.parametrize("mol,ref,func,ptype,basis", [
    pytest.param( ch2,  'UHF',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, UHF, RPA]),
    pytest.param( ch2,  'UHF',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, UHF, TDA]),
    pytest.param( h2o,'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA]),
    pytest.param( h2o,'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA]),
    pytest.param( h2o,'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA]),
    pytest.param( h2o,'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA]),
    pytest.param( ch2,  'UHF',        'HF',  'RPA',  'cc-pvdz', marks=[hf, UHF, RPA, pytest.mark.quick]),
    pytest.param( ch2,  'UHF',        'HF',  'TDA',  'cc-pvdz', marks=[hf, UHF, TDA, pytest.mark.quick]),
    pytest.param( h2o,'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param( h2o,'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param( h2o,'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param( h2o,'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param( ch2,  'UHF',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, UHF, RPA]),
    pytest.param( ch2,  'UHF',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, UHF, TDA]),
    pytest.param( h2o,'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]),
    pytest.param( h2o,'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]),
    pytest.param( h2o,'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]),
    pytest.param( h2o,'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]),
    pytest.param( ch2,  'UHF',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, UHF, RPA]),
    pytest.param( ch2,  'UHF',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, UHF, TDA]),
    pytest.param( h2o,'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]),
    pytest.param( h2o,'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]),
    pytest.param( h2o,'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]),
    pytest.param( h2o,'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]),
    pytest.param( ch2,  'UHF',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, RPA]),
    pytest.param( ch2,  'UHF',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, TDA]),
    pytest.param( h2o,'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]),
    pytest.param( h2o,'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]),
    pytest.param( h2o,'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]),
    pytest.param( h2o,'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]),
]) # yapf: disable
def test_tdscf(mol, ref, func, ptype, basis):
    # expected failures
    if (ref == 'RHF-3'):
        pytest.xfail("RKS Vx kernel only Spin Adapted for Singlet")
    elif (ref == 'UHF' and func != 'SVWN'):
        pytest.xfail("UKS Vx kernel bug for non-LDA")

    psi4.core.clean()
    psi4.set_options({'scf_type': 'pk', 'e_convergence': 8, 'd_convergence': 8, 'save_jk': True})
    if ref == "UHF":
        psi4.set_options({'reference': 'UHF'})
    mol.reset_point_group('c1')
    _, wfn = psi4.energy(f"{func}/{basis}", return_wfn=True, molecule=mol)

    out = tdscf_excitations(wfn,
                            states=4,
                            maxiter=30,
                            r_tol=1.0e-5,
                            triplets="only" if ref == "RHF-3" else "none",
                            tda=True if ptype == "TDA" else False)

    mol = "CH2" if ref == "UHF" else "H2O"
    ref_v = reference_data[f"{mol}_{ref}_{func}_{ptype}"]

    for i, my_v in enumerate(out):
        assert compare_values(ref_v[i]["EXCITATION ENERGY"], my_v["EXCITATION ENERGY"], 4,
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1}")
