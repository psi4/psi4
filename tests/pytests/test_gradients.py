import numpy as np
import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Reference data generated from Psi's dfmp2 module
data = {
    "df-mp2 ae": np.array([[ 0, 0,  9.62190509e-03],
         [0,  5.49835030e-03, -4.81095255e-03],
         [0, -5.49835030e-03, -4.81095255e-03]]),
    "df-mp2 fc": np.array([[ 0, 0,  1.02432654e-02],
         [0,  5.88581965e-03, -5.12163268e-03],
         [0, -5.88581965e-03, -5.12163268e-03]]),
    "df-mp2 fv": np.array([[ 0, 0,  1.22918833e-02],
         [0,  5.52107556e-03, -6.14594166e-03],
         [0, -5.52107556e-03, -6.14594166e-03]]),
    "df-mp2 fc/fv": np.array([[ 0, 0,  1.25984187e-02],
         [0,  5.71563223e-03, -6.29920936e-03],
         [0, -5.71563223e-03, -6.29920936e-03]]),
    "df-dct": np.array([[0, 0, 0.008477558394],
        [0,  0.005148825942, -0.004238779197],
        [0, -0.005148825942, -0.004238779197]]),
    "df-cc2": np.array([[0, 0, 0.011903811700],
        [0,  0.006730035450, -0.005951905850],
        [0, -0.006730035450, -0.005951905850]])
    }

@pytest.mark.slow
# TODO: That "true" needs to be a string is silly. Convert it to a boolean when you can do that without incurring a NaN energy.
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df'}, 'ref': data["df-mp2 ae"]}, id='df-mp2 ae'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'freeze_core': 'true'}, 'ref': data["df-mp2 fc"]}, id='df-mp2 fc'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'num_frozen_uocc': 4}, 'ref': data["df-mp2 fv"]}, id='df-mp2 fv'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'freeze_core': 'true', 'num_frozen_uocc': 4}, 'ref': data["df-mp2 fc/fv"]}, id='df-omp2 fc/fv'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df'}, 'ref': data["df-dct"]}, id='df-rdct'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df', 'reference': 'uhf'}, 'ref': data["df-dct"]}, id='df-udct'),
    pytest.param({'name': 'cc2', 'options': {'cc_type': 'df'}, 'ref': data["df-cc2"]}, id='df-cc2')
    ]
)
def test_gradient(inp):
    h2o = psi4.geometry("""
        O
        H 1 0.958
        H 1 0.958 2 104.5
    """)

    psi4.set_options({'basis': 'aug-cc-pvdz', 'points': 5})
    psi4.set_options(inp['options'])
    
    analytic_gradient = psi4.gradient(inp['name'], dertype=1)
    print(analytic_gradient)
    findif_gradient = psi4.gradient(inp['name'], dertype=0)
    reference_gradient = inp["ref"]

    assert compare_values(findif_gradient, analytic_gradient, 5, "analytic vs. findif gradient")
    assert compare_values(reference_gradient, analytic_gradient.np, 5, "analytic vs. reference gradient")

def test_gradient_ref():
    h2o = psi4.geometry("""
        O
        H 1 0.958
        H 1 0.958 2 104.5
    """)

    psi4.set_options({"basis": "cc-pVDZ"})
    mp2_wfn = psi4.energy("mp2", return_wfn=True)[1]
    with pytest.raises(TypeError):
        psi4.gradient("scf", ref_wfn=mp2_wfn)
    scf_wfn = psi4.energy("scf", return_wfn=True)[1]
    psi4.gradient("scf", ref_wfn=scf_wfn)
