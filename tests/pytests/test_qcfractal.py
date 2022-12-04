import numpy as np

import pytest
from utils import *
from addons import uusing

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


@uusing("qcfractal")
@uusing("qcportal")
@pytest.mark.parametrize("parallel", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake"),
])
def test_qcf_cbs_mbe(parallel):
    if parallel:
        try:
            from qcfractal import FractalSnowflake
            qca_next_branch = False
        except ImportError:
            from qcfractal.snowflake import FractalSnowflake
            qca_next_branch = True

        if not qca_next_branch:
            snowflake = FractalSnowflake(logging=True, max_workers=4)
        else:
            snowflake = FractalSnowflake()
        client = snowflake.client()
    
    import psi4
    dimer = psi4.geometry("""
    He 2 0 0
    --
    He -2 0 0
    """)
    
    ref_grad = [[-3.987697668786454e-07, 0.0, 0.0], [3.987697668787284e-07, 0.0, 0.0]]

    if parallel:
        plan = psi4.gradient("HF/cc-pV[D,T]Z", bsse_type="CP", return_plan=True, return_total_data=True)
        plan.compute(client)
        
        snowflake.await_results()
        ret = plan.get_results(client)
        assert compare_values(ref_grad, ret.return_result)
    
    else:
        grad = psi4.gradient("HF/cc-pV[D,T]Z", bsse_type="CP", return_total_data=True)
        assert compare_values(ref_grad, grad)
    
    if parallel:
        print(f'Final energy   = {ret.extras["qcvars"]["CURRENT ENERGY"]} [E_h]')
        print(f'Final gradient = {ret.extras["qcvars"]["CURRENT GRADIENT"]} [E_h/a0]')
        print(f'Final gradient = {ret.return_result} [E_h/a0]')
        print(f'Final energy   = {ret.properties.return_energy} [E_h]')
    
        assert compare_values(-5.72527135604184, ret.extras["qcvars"]["CURRENT ENERGY"], atol=1e-5)
        assert compare_values(ref_grad, ret.extras["qcvars"]["CURRENT GRADIENT"], atol=1e-4)
    
        snowflake.stop()
