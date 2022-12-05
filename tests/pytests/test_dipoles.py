import pytest

import qcelemental as qcel
import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

perturbation_strength = 0.001

@pytest.mark.slow
# TODO: That "true" needs to be a string is silly. Convert it to a boolean when you can do that without incurring a NaN energy.
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'omp2', 'options': {'mp2_type': 'df', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP2'}, id='df-omp2 ae'),
    pytest.param({'name': 'omp2', 'options': {'mp2_type': 'df', 'freeze_core': 'true', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP2'}, id='df-omp2 fc'),
    pytest.param({'name': 'omp3', 'options': {'mp_type': 'df', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP3'}, id='df-omp3 ae'),
    pytest.param({'name': 'omp3', 'options': {'mp_type': 'df', 'freeze_core': 'true', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP3'}, id='df-omp3 fc'),
    pytest.param({'name': 'omp2.5', 'options': {'mp_type': 'df', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP2.5'}, id='df-omp2.5 ae'),
    pytest.param({'name': 'omp2.5', 'options': {'mp_type': 'df', 'freeze_core': 'true', 'max_mograd_convergence': 6}, 'varname': 'DF-OMP2.5'}, id='df-omp2.5 fc'),
    pytest.param({'name': 'olccd', 'options': {'cc_type': 'df', 'max_mograd_convergence': 6}, 'varname': 'DF-OLCCD'}, id='df-olccd ae'),
    pytest.param({'name': 'olccd', 'options': {'cc_type': 'df', 'freeze_core': 'true', 'max_mograd_convergence': 6}, 'varname': 'DF-OLCCD'}, id='df-olccd fc'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df'}, 'varname': 'DCT'}, id='df-rdct'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df', 'reference': 'uhf'}, 'varname': 'DCT'}, id='df-udct'),
    pytest.param({'name': 'ccsd', 'options': {'opdm_relax': 'true'}, 'varname': 'CCSD'}, id='ccsd'),
    pytest.param({'name': 'ccsd', 'options': {'opdm_relax': 'true', 'reference': 'uhf'}, 'varname': 'CCSD'}, id='ccsd'),
    ]
)
def test_dipole(inp):
    h2o_singlet = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)
    h2o_doublet = psi4.geometry("""
        1 2
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    mol = h2o_singlet if inp["options"].get("reference", "rhf") == "rhf" else h2o_doublet
    psi4.set_options({'perturb_h': True, 'perturb_with': 'dipole', 'basis': 'cc-pvdz'})
    psi4.set_options(inp['options'])
    energies = dict()
    for l in [1, -1, 2, -2]:
        psi4.set_options({'perturb_dipole': [0, 0, l * perturbation_strength]})
        energies[l] = psi4.energy(inp['name'], molecule=mol)
    findif_dipole = [0, 0, (8 * energies[1] - 8 * energies[-1] - energies[2] + energies[-2]) / (12 * perturbation_strength)]

    psi4.set_options({'perturb_h': False})
    wfn = psi4.properties(inp['name'], properties=['dipole'], molecule=mol, return_wfn=True)[1]
    analytic_dipole = wfn.variable(inp['varname'] + " DIPOLE")

    assert compare_values(findif_dipole, analytic_dipole, 5, "findif vs. analytic dipole")

