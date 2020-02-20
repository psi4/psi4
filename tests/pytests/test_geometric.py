import pytest

from .utils import *

import psi4
from qcengine.testing import using


@using('geometric')
@pytest.mark.parametrize('opt', [
    pytest.param({'engine' : 'optking'}, id='optking'),
    pytest.param({'engine' : 'geometric'}, id='geometric'),
]) #yapf : disable
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.027032783717, 'ref_nuc': 9.300794299874}, id='rhf(df)'),
    pytest.param({'name': 'hf', 'options': {'scf_type': 'pk'}, 'ref_ene' : -76.027053512764, 'ref_nuc': 9.300838770294}, id='rhf(pk)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df'}, 'ref_ene' : -76.230938589591, 'ref_nuc': 9.133271168193}, id='mp2(df)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'conv'}, 'ref_ene' : -76.230989373502, 'ref_nuc': 9.133125471291}, id='mp2(conv)'),
    pytest.param({'name': 'b3lyp', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.420645414834, 'ref_nuc': 9.090397129492}, id='b3lyp'),
])  # yapf: disable
def test_h2o(inp, opt):
    """Optimizer the square water molecule"""

    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 90.0
    """)

    psi4.set_options({'basis': 'cc-pvdz',
                      'g_convergence': 'gau_tight'
                     })
    psi4.set_options(inp['options'])

    e, wfn = psi4.optimize(inp['name'], return_wfn=True, engine=opt['engine'])
    assert compare_values(inp['ref_ene'], e, 6)
    assert compare_values(inp['ref_nuc'], h2o.nuclear_repulsion_energy(), 3)

