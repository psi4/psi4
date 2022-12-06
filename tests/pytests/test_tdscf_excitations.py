from pathlib import Path
import json

import numpy as np
import pytest

import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from utils import *

pytestmark = [pytest.mark.psi, pytest.mark.api]

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


@pytest.fixture
def reference_data():
    # Reference data generated using G09.E01
    with open(Path(__file__).parent / "tdscf_reference_data.json") as f:
        reference_data = json.load(f)
    return reference_data


@pytest.fixture
def molecules():
    smols = {
        # Canonical unrestricted system
        "CH2":
        """0 3
    C           0.000000    0.000000    0.159693
    H          -0.000000    0.895527   -0.479080
    H          -0.000000   -0.895527   -0.479080
    no_reorient
    no_com
    """,
        # Canonical restricted system
        "H2O":
        """0 1
    O           0.000000    0.000000    0.135446
    H          -0.000000    0.866812   -0.541782
    H          -0.000000   -0.866812   -0.541782
    no_reorient
    no_com
    """,
        # Canonical chiral system
        "H2O2":
        """0 1
    O        0.000000    0.695000   -0.092486
    O       -0.000000   -0.695000   -0.092486
    H       -0.388142    0.895249    0.739888
    H        0.388142   -0.895249    0.739888
    no_reorient
    no_com
    """,
        # Slightly larger chiral system
        "METHYLOXIRANE":
        """0 1
    C  0.152133 -0.035800  0.485797
    C -1.039475  0.615938 -0.061249
    C  1.507144  0.097806 -0.148460
    O -0.828215 -0.788248 -0.239431
    H  0.153725 -0.249258  1.552136
    H -1.863178  0.881921  0.593333
    H -0.949807  1.214210 -0.962771
    H  2.076806 -0.826189 -0.036671
    H  2.074465  0.901788  0.325106
    H  1.414895  0.315852 -1.212218
    no_reorient
    no_com
    """,
    }

    return {k: psi4.core.Molecule.from_string(v) for k, v in smols.items()}


def _oscillator_strength(e: float, tm: np.ndarray, gauge: str = "L") -> float:
    if gauge == "L":
        return ((2 * e) / 3) * np.sum(tm**2)
    else:
        return (2 / (3 * e)) * np.sum(tm**2)


def _rotatory_strength(e: float, etm: np.ndarray, mtm: np.ndarray, gauge: str = "L") -> float:
    """Compute the rotatory strength from the G09 reference values.

    Notes
    -----
    1. Signs are flipped with respect to the definition!
    2. The magnetic dipole moment is really the angular momentum, so we scale
    it by 1/2 (Bohr magneton) to get the magnetic dipole.
    """
    if gauge == "L":
        return -np.einsum("i,i", etm, 0.5 * mtm)
    else:
        return np.einsum("i,i", etm, 0.5 * mtm) / e


@pytest.mark.tdscf
@pytest.mark.parametrize("mol,ref,func,ptype,basis", [
    pytest.param(  "CH2",   'UHF',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, UHF, RPA, pytest.mark.quick]),
    pytest.param(  "CH2",   'UHF',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, UHF, TDA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param(  "CH2",   'UHF',        'HF',  'RPA',  'cc-pvdz', marks=[hf, UHF, RPA, pytest.mark.quick]),
    pytest.param(  "CH2",   'UHF',        'HF',  'TDA',  'cc-pvdz', marks=[hf, UHF, TDA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param(  "H2O", 'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA, pytest.mark.quick]),
    pytest.param( "H2O2", 'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA, pytest.mark.quick]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA]),
    pytest.param(  "CH2",   'UHF',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, UHF, RPA]),
    pytest.param(  "CH2",   'UHF',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, UHF, TDA]),
    pytest.param(  "H2O", 'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]),
    pytest.param(  "H2O", 'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]),
    pytest.param(  "H2O", 'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]),
    pytest.param(  "H2O", 'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]),
    pytest.param( "H2O2", 'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]),
    pytest.param( "H2O2", 'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]),
    pytest.param( "H2O2", 'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]),
    pytest.param( "H2O2", 'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]),
    pytest.param(  "CH2",   'UHF',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, UHF, RPA]),
    pytest.param(  "CH2",   'UHF',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, UHF, TDA]),
    pytest.param(  "H2O", 'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]),
    pytest.param(  "H2O", 'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]),
    pytest.param(  "H2O", 'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]),
    pytest.param(  "H2O", 'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]),
    pytest.param( "H2O2", 'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]),
    pytest.param( "H2O2", 'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]),
    pytest.param( "H2O2", 'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]),
    pytest.param( "H2O2", 'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]),
    pytest.param(  "CH2",   'UHF',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, RPA]),
    pytest.param(  "CH2",   'UHF',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, TDA]),
    pytest.param(  "H2O", 'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]),
    pytest.param(  "H2O", 'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]),
    pytest.param(  "H2O", 'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]),
    pytest.param(  "H2O", 'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]),
    pytest.param( "H2O2", 'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]),
    pytest.param( "H2O2", 'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]),
    pytest.param( "H2O2", 'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]),
    pytest.param( "H2O2", 'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]),
    pytest.param( "METHYLOXIRANE", 'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]),
]) # yapf: disable
def test_tdscf(mol, ref, func, ptype, basis, molecules, reference_data):
    # expected failures
    if (ref == 'RHF-3') and (func != "HF"):
        pytest.xfail("RKS Vx kernel only Spin Adapted for Singlet")

    molecule = molecules[mol]
    psi4.set_options({'scf_type': 'pk', 'e_convergence': 8, 'd_convergence': 8, 'save_jk': True})
    if ref == "UHF":
        psi4.set_options({'reference': 'UHF'})
    molecule.reset_point_group('c1')
    _, wfn = psi4.energy(f"{func}/{basis}", return_wfn=True, molecule=molecule)

    out = tdscf_excitations(wfn,
                            states=4,
                            maxiter=30,
                            r_convergence=1.0e-6,
                            triplets="ONLY" if ref == "RHF-3" else "NONE",
                            tda=True if ptype == "TDA" else False)

    ref_v = reference_data[f"{mol}_{ref}_{func}_{ptype}"]

    for i, my_v in enumerate(out):

        # compare excitation energies
        ref_e = ref_v[i]["EXCITATION ENERGY"]
        assert compare_values(ref_e,
                              my_v["EXCITATION ENERGY"],
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1} Excitation energy",
                              atol=2.0e-4)

        ref_edtm_L = np.array(ref_v[i]["LENGTH MU"])
        # compare length-gauge oscillator strength
        ref_f_L = _oscillator_strength(ref_e, ref_edtm_L, "L")
        assert compare_values(ref_f_L,
                              my_v["OSCILLATOR STRENGTH (LEN)"],
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1} Length-gauge oscillator strength",
                              atol=1.0e-3)

        ref_edtm_V = np.array(ref_v[i]["VELOCITY MU"])
        # compare velocity-gauge oscillator strengths
        ref_f_V = _oscillator_strength(ref_e, ref_edtm_V, "V")
        assert compare_values(ref_f_V,
                              my_v["OSCILLATOR STRENGTH (VEL)"],
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1} Velocity-gauge oscillator strength",
                              atol=1.0e-2)

        ref_mdtm = np.array(ref_v[i]["M"])
        # compare length-gauge rotatory strengths
        ref_R_L = _rotatory_strength(ref_e, ref_edtm_L, ref_mdtm, "L")
        assert compare_values(ref_R_L,
                              my_v["ROTATORY STRENGTH (LEN)"],
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1} Length-gauge rotatory strength",
                              atol=2.0e-3)

        # compare velocity-gauge rotatory strengths
        ref_R_V = _rotatory_strength(ref_e, ref_edtm_V, ref_mdtm, "V")
        assert compare_values(ref_R_V,
                              my_v["ROTATORY STRENGTH (VEL)"],
                              f"{mol}_{ref}_{func}_{ptype}-ROOT_{i+1} Velocity-gauge rotatory strength",
                              atol=2.0e-3)
