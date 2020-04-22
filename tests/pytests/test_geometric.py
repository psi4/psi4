import pytest

from .utils import *

import psi4
from qcengine.testing import using


@pytest.mark.parametrize('engine', [
    pytest.param('optking'),
    pytest.param('geometric', marks=using('geometric')),
])  # yapf: disable
@pytest.mark.parametrize('inp', [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.027032783717, 'ref_nuc': 9.300794299874}, id='rhf(df)'),
    pytest.param({'name': 'hf', 'options': {'scf_type': 'pk'}, 'ref_ene' : -76.027053512764, 'ref_nuc': 9.300838770294}, id='rhf(pk)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df'}, 'ref_ene' : -76.230938589591, 'ref_nuc': 9.133271168193}, id='mp2(df)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'conv'}, 'ref_ene' : -76.230989373502, 'ref_nuc': 9.133125471291}, id='mp2(conv)'),
    pytest.param({'name': 'b3lyp', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.420645414834, 'ref_nuc': 9.090397129492}, id='b3lyp'),
])  # yapf: disable
def test_h2o(inp, engine):
    """Optimization of the square water molecule"""

    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 90.0
    """)

    psi4.set_options({'basis': 'cc-pvdz',
                      'g_convergence': 'gau_tight'
                     })
    psi4.set_options(inp['options'])

    e, wfn = psi4.optimize(inp['name'], return_wfn=True, engine=engine)
    assert compare_values(inp['ref_ene'], e, 6)
    assert compare_values(inp['ref_nuc'], h2o.nuclear_repulsion_energy(), 3)


@using('geometric')
@pytest.mark.parametrize('inp', [
    pytest.param({'name': 'hf', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.02079629252714, 'ref_nuc': 9.265341708725257}, id='rhf(df)'),
    pytest.param({'name': 'hf', 'options': {'scf_type': 'pk'}, 'ref_ene' : -76.02082389228, 'ref_nuc': 9.26528625744628}, id='rhf(pk)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df'}, 'ref_ene' : -76.22711819393223, 'ref_nuc': 9.09137805747361}, id='mp2(df)'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'conv'}, 'ref_ene' : -76.2271678506303, 'ref_nuc': 9.091178486990861}, id='mp2(conv)'),
    pytest.param({'name': 'b3lyp', 'options': {'scf_type': 'df'}, 'ref_ene' : -76.41632755714534, 'ref_nuc': 9.04535641436914}, id='b3lyp'),
])  # yapf: disable
def test_h2o_constrained(inp):
    """Constrained optimization of the square water molecule"""

    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 90.0
    """)

    psi4.set_options({'basis': 'cc-pvdz',
                      'g_convergence': 'gau_tight'
                     })
    psi4.set_options(inp['options'])

    # geometric specific options 
    geometric_keywords = {
        'coordsys' : 'tric',
        'enforce' : 0.0,
        'constraints' : { 
            'set' : [{'type'    : 'angle',
                      'indices' : [1, 0, 2],
                      'value'   : 90.0 }]
        }
    }

    e, wfn = psi4.optimize(inp['name'], return_wfn=True, engine='geometric', optimizer_keywords=geometric_keywords)
    assert compare_values(inp['ref_ene'], e, 6)
    assert compare_values(inp['ref_nuc'], h2o.nuclear_repulsion_energy(), 3)

