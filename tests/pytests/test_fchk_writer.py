import pytest
import psi4
import numpy as np
import re
from ast import literal_eval
import os
from shutil import copytree

pytestmark = [pytest.mark.psi, pytest.mark.api]

# checks for
# - correct HF density
# - principal execution 
# - comparison against reference file

@pytest.fixture
def datadir(tmpdir, request):
    """
    from: https://stackoverflow.com/a/29631801
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        copytree(test_dir, str(tmpdir), dirs_exist_ok=True)

    return tmpdir


def calcD(wfn):
    Ca_occ = wfn.Ca_subset("AO", "OCC").np
    Cb_occ = wfn.Cb_subset("AO", "OCC").np
    Da = np.dot(Ca_occ, Ca_occ.T)
    Db = np.dot(Cb_occ, Cb_occ.T)
    return Da + Db

@pytest.mark.parametrize('inp2', [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'} }, id='df-uhf'),
    pytest.param({'name': 'pbe', 'options': {'scf_type': 'df'} }, id='df-uhf-dft'), 
    pytest.param({'name': 'mp2', 'options': {'scf_type': 'df','qc_module': 'occ'} }, id='df-uhf-mp2'),
    pytest.param({'name': 'ccsd', 'options': {'scf_type': 'pk','cc_type':'conv'} }, id='conv-uhf-ccsd') 
    ])
def test_uhf_fchk(inp2, datadir):
    """  FCHK UHF """
    mol = psi4.geometry("""
  no_reorient
  0 2
  O
  O 1 1.46
  H 2 0.97 1 104.6
  """)
    psi4.set_options({
        "BASIS": "pcseg-0",
        'reference': 'uhf',
        'e_convergence': 1e-12,
        'd_convergence': 1e-12,
        'r_convergence': 1e-10,
        'pcg_convergence': 1e-10,
    })
    if psi4.core.get_option("scf", "orbital_optimizer_package") == "INTERNAL":
        dens_tol = 9
        fchk_tol = 1.e-8
    else:
        dens_tol = 3.e-8
        fchk_tol = 2.e-6
        psi4.set_options({"e_convergence": 9, "d_convergence": 2e-8})
    FCHK_file = f"uhf-{inp2['name']}.fchk"
    reference_file = datadir.join(f"uhf-{inp2['name']}.ref")
    psi4.set_options(inp2['options'])
    e, wfn = psi4.gradient(inp2['name'], return_wfn=True, molecule=mol)
    ret = psi4.driver.fchk(wfn, FCHK_file, debug=True)
    assert psi4.compare_arrays(ret["Total SCF Density"], calcD(wfn), dens_tol, "FCHK UHF Density")
    assert psi4.compare_fchkfiles(reference_file, FCHK_file, fchk_tol, f" File comparison: {FCHK_file}")

@pytest.mark.parametrize('inp', [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'} }, id='df-rhf)'),
    pytest.param({'name': 'pbe', 'options': {'scf_type': 'df'} }, id='df-rhf-dft)'),
    pytest.param({'name': 'mp2', 'options': {'scf_type': 'df','mp2_type':'df'} }, id='df-rhf-mp2'),
    pytest.param({'name': 'omp2', 'options': {'scf_type': 'df','mp2_type':'df'} }, id='df-rhf-omp2'),
    pytest.param({'name': 'cc2', 'options': {'scf_type': 'pk','cc_type':'conv'}}, id='conv-rhf-cc2'),
    pytest.param({'name': 'ccsd', 'options': {'scf_type': 'pk','cc_type':'conv'}}, id='conv-rhf-ccsd'),
    pytest.param({'name': 'dct', 'options': {'scf_type': 'pk','dct_type':'conv'}}, id='conv-rhf-dct'),
    pytest.param({'name': 'mp2', 'options': {'scf_type': 'pk','mp2_type':'conv','qc_module':'occ'}}, marks=pytest.mark.xfail(reason="OCC not allowed in FCHK"), id='conv-rhf-mp2(occ)'),
    ])
def test_rhf_fchk(inp, datadir):
    """  FCHK RHF """
    mol = psi4.geometry("""
  no_reorient
  0 1
  O
  H 1 1.01
  H 1 1.0 2 104.5
  """)
    psi4.set_options({
        "BASIS": "pcseg-0",
        'e_convergence': 1e-12,
        'd_convergence': 1e-12,
        'r_convergence': 1e-10,
        'pcg_convergence': 1e-10,
    })
    if psi4.core.get_option("scf", "orbital_optimizer_package") == "INTERNAL":
        dens_tol = 9
        fchk_tol = 1.e-8
    else:
        dens_tol = 3.e-8
        fchk_tol = 5.e-7
        psi4.set_options({"e_convergence": 9, "d_convergence": 5e-9})
    FCHK_file = f"rhf-{inp['name']}.fchk"
    reference_file = datadir.join(f"rhf-{inp['name']}.ref")
    psi4.set_options(inp['options'])
    e, wfn = psi4.gradient(inp['name'], return_wfn=True, molecule=mol)
    ret = psi4.driver.fchk(wfn, FCHK_file, debug=True)
    if inp['name'] in ['mp2', 'dct', 'omp2']:
        refwfn = wfn.reference_wavefunction()
        expected = calcD(refwfn)
    else:
        expected = calcD(wfn)
    assert psi4.compare_arrays(ret["Total SCF Density"], expected, dens_tol, "FCHK RHF Density")
    assert psi4.compare_fchkfiles(reference_file, FCHK_file, fchk_tol, f" File comparison: {FCHK_file}")

