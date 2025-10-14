import numpy as np

import pytest
from utils import *
from addons import uusing

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


@uusing("qcfractal")
@uusing("qcportal")
@pytest.mark.cbs
@pytest.mark.smoke
@pytest.mark.quick
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake"),
])
def test_qcf_cbs_mbe(distributed, snowflake):
    
    import psi4
    dimer = psi4.geometry("""
    He 2 0 0
    --
    He -2 0 0
    """)
    
    # << pk >> (just to test options passing works)
    #psi4.set_options({"scf_type": "pk"}) #, "e_convergence": 10, "d_convergence": 9})
    #ref_ene = -5.725227008286358
    #ref_grad = [[-7.01014982e-07, 0.0, 0.0], [ 7.01014982e-07, 0.0, 0.0]]

    # << df >>
    ref_ene = -5.72527135604184
    ref_grad = [[-3.987697668786454e-07, 0.0, 0.0], [3.987697668787284e-07, 0.0, 0.0]]

    if distributed:
        client = snowflake.client()
        plan = psi4.gradient("HF/cc-pV[D,T]Z", bsse_type="CP", return_plan=True, return_total_data=True)
        plan.compute(client)
        snowflake.await_results()
        grad = plan.get_psi_results(client)

    else:
        grad = psi4.gradient("HF/cc-pV[D,T]Z", bsse_type="CP", return_total_data=True)

    assert compare_values(ref_ene, psi4.variable("CURRENT ENERGY"))
    assert compare_values(ref_grad, grad, atol=1.e-7)
    
    if distributed:
        # `get_results` is a closer-to-internals alternative to `get_psi_results`.
        #   It grabs the AtomicResult-compliant QCSchema model directly, rather
        #   than forming a dummy wfn and setting qcvars to it and to globals.
        #   As schema are formalized, be sure to favor official properties, and
        #   don't assume QCVariables or Psi4 specialties are available before get_psi_results.
        ret = plan.get_results(client)

        # print(f'Final energy   = {ret.extras["qcvars"]["CURRENT ENERGY"]} [E_h]')
        # print(f'Final gradient = {ret.extras["qcvars"]["CURRENT GRADIENT"]} [E_h/a0]')
        print(f'Final gradient = {ret.return_result} [E_h/a0]')
        print(f'Final energy   = {ret.properties.return_energy} [E_h]')
    
        # assert compare_values(ref_ene, ret.extras["qcvars"]["CURRENT ENERGY"], atol=1e-5)
        # assert compare_values(ref_grad, ret.extras["qcvars"]["CURRENT GRADIENT"], atol=1e-4)
        assert compare_values(ref_ene, ret.properties.return_energy)
        assert compare_values(ref_grad, ret.return_result)
    
        snowflake.stop()
