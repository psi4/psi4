import pytest
import psi4
import numpy as np
import re
from ast import literal_eval

# checks for
# - correct HF density
# - principal execution 

def calcD(wfn):
  Ca_occ=wfn.Ca_subset("AO","OCC").np
  Cb_occ=wfn.Cb_subset("AO","OCC").np
  Da=np.dot(Ca_occ,Ca_occ.T)
  Db=np.dot(Cb_occ,Cb_occ.T)
  return Da+Db

@pytest.mark.parametrize('inp2', [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'} }, id='df-uhf'),
    pytest.param({'name': 'pbe', 'options': {'scf_type': 'df'} }, id='df-uhf-dft'), 
    pytest.param({'name': 'mp2', 'options': {'scf_type': 'df','qc_module': 'occ'} }, id='df-uhf-mp2'),
    pytest.param({'name': 'ccsd', 'options': {'scf_type': 'pk','cc_type':'conv'} }, id='conv-uhf-ccsd') 
    ])
def test_uhf_fchk(inp2):
  """  FCHK UHF """
  mol = psi4.geometry("""
  0 3
  O 0.0 0.0 0.0
  O 0.0 0.0 1.1
  """)
  psi4.set_options({
  "BASIS": "pcseg-0",
  'reference':'uhf',
  'd_convergence': 1e-10,
  })
  psi4.set_options(inp2['options'])
  e, wfn = psi4.gradient(inp2['name'], return_wfn=True,molecule=mol)
  ret = psi4.driver.fchk(wfn, f"uhf-{inp2['name']}.fchk", debug=True)
  psi4.compare_arrays(ret["Total SCF Density"],calcD(wfn),9)


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
def test_rhf_fchk(inp):
  """  FCHK RHF """
  mol = psi4.geometry("""
  0 1
  O
  H 1 1.0
  H 1 1.0 2 104.5
  # symmetry c1
  """)
  psi4.set_options({
  "BASIS": "pcseg-0",
  })
  psi4.set_options(inp['options'])
  e, wfn = psi4.gradient(inp['name'], return_wfn=True, molecule=mol)
  ret = psi4.driver.fchk(wfn, f"rhf-{inp['name']}.fchk", debug=True)
  if inp['name'] in ['mp2','dct','omp2']: 
      refwfn = wfn.reference_wavefunction()
      expected = calcD(refwfn)
  else:
      expected = calcD(wfn)
  psi4.compare_arrays(ret["Total SCF Density"],expected,9)
