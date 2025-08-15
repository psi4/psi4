import pytest
import numpy as np

import psi4
from utils import compare_values
from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]

data = {
     'hf': {
        'Ba2': {
            'grad': np.array([[0, 0,  2.08544743e-02],
                [ 0,  0, -2.08544743e-02]]),
            'energy': -50.12131931846057,
            'docc' : [ 3, 0, 1, 1, 0, 3, 1, 1 ]},
        'In2': {
            'grad': np.array([[0, 0,  -1.46964113e-02],
                [ 0,  0, 1.46964113e-02]]),
            'energy': -378.3986291158715,
            'docc' : [ 6, 1, 2, 2, 1, 5, 2, 2 ]},
        'Sr2': {
            'grad': np.array([[0, 0,  9.74576806e-03],
                [ 0,  0, -9.74576806e-03]]),
            'energy': -60.66907844700816,
            'docc': [ 3, 0, 1, 1, 0, 3, 1, 1 ]}
     },
     'mp2': {
        'Ba2': {
            'grad': np.array([[0, 0, 1.6421392746403734e-02],
                     [0, 0, -1.6421392746403734e-02]]),
            'energy': -50.34353063117147,
            'docc' : [ 3, 0, 1, 1, 0, 3, 1, 1 ]},
        'In2': {
            'grad': np.array([[0, 0, -1.4464867319381867e-02],
                     [0, 0, 1.4464867319381867e-02]]),
            'energy': -378.5493603721844,
            'docc' : [ 6, 1, 2, 2, 1, 5, 2, 2 ]},
        'Sr2': {
            'grad': np.array([[0, 0, 6.469022590311186e-03],
                     [0, 0, -6.469022590311186e-03]]),
            'energy': -60.83705034568311,
            'docc': [ 3, 0, 1, 1, 0, 3, 1, 1 ]}
     }
}

@pytest.fixture
def mols():
    smols = {
        "Ba2": """
Ba
Ba 1 8.0
units bohr
""",
        "In2": """
In
In 1 8.0
units bohr
""",
        "Sr2": """
Sr
Sr 1 8.0
units bohr
""",
    }

    return {k: psi4.core.Molecule.from_string(v) for k, v in smols.items()}

@pytest.mark.parametrize('inp', [
    pytest.param({'name': 'hf', 'method': 'hf', 'options': {'basis': 'def2-svp'}, 'ref': data['hf']}, id='rhf(df)'),
    pytest.param({'name': 'mp2', 'method': 'mp2', 'options': {'basis': 'def2-svp', 'mp2_type': 'df'}, 'ref': data['mp2']}, id='mp2(df)'),
])

@pytest.mark.parametrize('mol', [
    pytest.param('Ba2', id='Ba2'),
    pytest.param('In2', id='In2'),
    pytest.param('Sr2', id='Sr2')
])

@uusing("ecpint")
def test_large_atoms(inp, mol, mols):

    ref_e = inp['ref'][mol]['energy']
    ref_grad = inp['ref'][mol]['grad']

    psi4.set_options(inp['options'])
    psi4.set_options({'points': 5})
    # The atoms are so far apart SCF convergence can land on different
    # solutions, and this is why it is important to set the wanted
    # occupations
    psi4.set_options({'docc' : inp['ref'][mol]['docc']})

    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 5e-9})

    method = inp['method']
    analytic_grad, wfn = psi4.gradient(method, molecule=mols[mol], dertype=1, return_wfn=True)
    findif_grad = psi4.gradient(method, molecule=mols[mol], dertype=0)
    energy = wfn.energy()

    assert compare_values(ref_e, energy, 5, "calculated vs. reference energy")
    assert compare_values(findif_grad, analytic_grad, 5, "analytic vs. findif gradient")
    assert compare_values(ref_grad, analytic_grad.np, 5, "analytic vs. reference gradient")

