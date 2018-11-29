import pytest

from utils import *

import psi4

#! conventional and density-fitting mp2 test of mp2 itself and setting scs-mp2

_ref_h2o_ccpvdz = {
    'df': {
        'HF TOTAL ENERGY':                      -76.0167614256151865,
        'MP2 SAME-SPIN CORRELATION ENERGY':      -0.0527406422061238,
        'MP2 OPPOSITE-SPIN CORRELATION ENERGY':  -0.1562926850310142,
        'MP2 CORRELATION ENERGY':                -0.2090333272371381,
        'MP2 TOTAL ENERGY':                     -76.2257947528523232,
        'SCS-MP2 CORRELATION ENERGY':            -0.2051314361059251,
        'SCS-MP2 TOTAL ENERGY':                 -76.2218928617211162,
    },
    'conv': {
        'HF TOTAL ENERGY':                      -76.01678947133706,
        'MP2 SAME-SPIN CORRELATION ENERGY':      -0.05268120425816,
        'MP2 OPPOSITE-SPIN CORRELATION ENERGY':  -0.15637564436589,
        'MP2 CORRELATION ENERGY':                -0.20905684862405,
        'MP2 TOTAL ENERGY':                     -76.22584631996111,
        'SCS-MP2 CORRELATION ENERGY':            -0.20521117465845,
        'SCS-MP2 TOTAL ENERGY':                 -76.22200064599551,
    },
}
for mp2type in ['df', 'conv']:
    _ref_h2o_ccpvdz[mp2type]['SCF TOTAL ENERGY'] = _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']
    _ref_h2o_ccpvdz[mp2type]['CURRENT REFERENCE ENERGY'] = _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']
    _ref_h2o_ccpvdz[mp2type]['5050SCS-MP2 CORRELATION ENERGY'] = (0.5 * (_ref_h2o_ccpvdz[mp2type]['MP2 SAME-SPIN CORRELATION ENERGY'] +
                                                                         _ref_h2o_ccpvdz[mp2type]['MP2 OPPOSITE-SPIN CORRELATION ENERGY']))
    _ref_h2o_ccpvdz[mp2type]['5050SCS-MP2 TOTAL ENERGY'] = _ref_h2o_ccpvdz[mp2type]['5050SCS-MP2 CORRELATION ENERGY'] + _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']


@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'df'}}, id='mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'conv'}}, id='mp2 (conv)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'df', 'mp2_os_scale': 1.2, 'mp2_ss_scale': 0.33333333}}, id='explicit scs mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'conv', 'mp2_os_scale': 1.2, 'mp2_ss_scale': 0.33333333}}, id='explicit scs mp2 (conv)'),
    pytest.param({'name': 'Mp2', 'custom': '5050SCS-MP2', 'options': {'mp2_type': 'df', 'mp2_os_scale': 0.5, 'mp2_ss_scale': 0.5}}, id='user-def scs mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': '5050SCS-MP2', 'options': {'mp2_type': 'conv', 'mp2_os_scale': 0.5, 'mp2_ss_scale': 0.5}}, id='user-def scs mp2 (conv)'),

    #pytest.param({'name': 'scs-mp2', 'options': {'mp2_type': 'df'}}, id='scs-mp2 (df)'),
    #pytest.param({'name': 'scs-mp2', 'options': {'mp2_type': 'conv'}}, id='scs-mp2 (conv)'),
])
def test_scsmp2(inp):
    """Formerly known as dfmp2-4"""

    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 90.0
    """)

    psi4.set_options({'basis': 'cc-pvdz'})
    psi4.set_options(inp['options'])

    ene, wfn = psi4.energy(inp['name'], return_wfn=True)

    ref_block = _ref_h2o_ccpvdz[inp['options']['mp2_type']]
    #ref_corl = ref_block[inp['pv'] + ' CORRELATION ENERGY']
    #ref_tot = ref_block[inp['pv'] + ' TOTAL ENERGY']
    ref_corl = ref_block['MP2 CORRELATION ENERGY']
    ref_tot = ref_block['MP2 TOTAL ENERGY']
    ref_custom_corl = ref_block[inp['custom'] + ' CORRELATION ENERGY']
    ref_custom_tot = ref_block[inp['custom'] + ' TOTAL ENERGY']

    for obj in [psi4.core, wfn]:
        for pv in [
            'HF TOTAL ENERGY',
            'SCF TOTAL ENERGY',
            'MP2 SAME-SPIN CORRELATION ENERGY',
            'MP2 OPPOSITE-SPIN CORRELATION ENERGY',
            'MP2 CORRELATION ENERGY',
            'MP2 TOTAL ENERGY',
            'SCS-MP2 CORRELATION ENERGY',
            'SCS-MP2 TOTAL ENERGY',
            'CURRENT REFERENCE ENERGY']:

            assert compare_values(ref_block[pv], obj.variable(pv), 5, pv)

        assert compare_values(ref_custom_corl, obj.variable('CUSTOM SCS-MP2 CORRELATION ENERGY'), 5, 'custom scsmp2 corl')
        assert compare_values(ref_custom_tot, obj.variable('CUSTOM SCS-MP2 TOTAL ENERGY'), 5, 'custom scsmp2 ')

        assert compare_values(ref_corl, obj.variable('CURRENT CORRELATION ENERGY'), 5, 'current corl')
        assert compare_values(ref_tot, obj.variable('CURRENT ENERGY'), 5, 'current')

    assert compare_values(ref_tot, ene, 5, 'return')
    assert compare_values(ref_tot, wfn.energy(), 5, 'wfn')



@pytest.fixture
def clsd_open_pmols():
    mols = {
        'hf': psi4.core.Molecule.from_string("""
                H
                F 1 0.917
              """),
        'bh_h2p': psi4.core.Molecule.from_string("""
                1 2
                B     0.10369114     0.00000000     0.00000000
                H    -1.13269886     0.00000000     0.00000000
                H     3.00000000     0.37149000     0.00000000
                H     3.00000000    -0.37149000     0.00000000
              """),
    }
    return mols

_ref_module = {ref: {alg: {} for alg in ['conv', 'df', 'cd']} for ref in ['rhf', 'uhf', 'rohf']}

_ref_module['rhf']['HF TOTAL ENERGY'] =                -100.019400605629
_ref_module['rhf']['conv']['MP2 CORRELATION ENERGY'] =   -0.201612517228
_ref_module['rhf']['conv']['MP2 TOTAL ENERGY'] =       -100.221013122857
_ref_module['rhf']['df']['MP2 CORRELATION ENERGY'] =     -0.201610660387
_ref_module['rhf']['df']['MP2 TOTAL ENERGY'] =         -100.221011266016
_ref_module['rhf']['cd']['MP2 CORRELATION ENERGY'] =     -0.201609396752
_ref_module['rhf']['cd']['MP2 TOTAL ENERGY'] =         -100.221010002381

_ref_module['uhf']['HF TOTAL ENERGY'] =                 -25.945130559147
_ref_module['uhf']['conv']['MP2 CORRELATION ENERGY'] =   -0.058421122206
_ref_module['uhf']['conv']['MP2 TOTAL ENERGY'] =        -26.003551681354
_ref_module['uhf']['df']['MP2 CORRELATION ENERGY'] =     -0.058390006825
_ref_module['uhf']['df']['MP2 TOTAL ENERGY'] =          -26.003520565972
_ref_module['uhf']['cd']['MP2 CORRELATION ENERGY'] =     -0.058409837177
_ref_module['uhf']['cd']['MP2 TOTAL ENERGY'] =          -26.003540396324

_ref_module['rohf']['HF TOTAL ENERGY'] =                -25.943606522029
_ref_module['rohf']['conv']['MP2 CORRELATION ENERGY'] =  -0.060939211739
_ref_module['rohf']['conv']['MP2 TOTAL ENERGY'] =       -26.004545733768
_ref_module['rohf']['df']['MP2 CORRELATION ENERGY'] =    -0.059372748391
_ref_module['rohf']['df']['MP2 TOTAL ENERGY'] =         -26.002979270420
_ref_module['rohf']['cd']['MP2 CORRELATION ENERGY'] =    -0.059393510962
_ref_module['rohf']['cd']['MP2 TOTAL ENERGY'] =         -26.003000032991


@pytest.mark.parametrize("inp", [
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'conv', 'qc_module': 'occ'}}, id='mp2 rhf conv: 3 occ*'),
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'conv', 'qc_module': 'fnocc'}}, id='mp2 rhf conv: 3 fnocc'),
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'conv', 'qc_module': 'detci'}}, id='mp2 rhf conv: 3 detci'),
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'df', 'qc_module': 'occ'}}, id='mp2 rhf df: 2 occ'),
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'df', 'qc_module': 'dfmp2'}}, id='mp2 rhf df: 2 dfmp2*'),
    pytest.param({'subject': 'hf', 'options': {'reference': 'rhf', 'mp2_type': 'cd', 'qc_module': 'occ'}}, id='mp2 rhf cd: 1 occ*'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'uhf', 'mp2_type': 'conv', 'qc_module': 'occ'}}, id='mp2 uhf conv: 1 occ*'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'uhf', 'mp2_type': 'df', 'qc_module': 'occ'}}, id='mp2 uhf df: 2 occ'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'uhf', 'mp2_type': 'df', 'qc_module': 'dfmp2'}}, id='mp2 uhf df: 2 dfmp2*'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'uhf', 'mp2_type': 'cd', 'qc_module': 'occ'}}, id='mp2 uhf cd: 1 occ*'),
    #pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ'}}, id='mp2 rohf conv: 2 occ*'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'detci'}}, id='mp2 rohf conv: 2 detci'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'df', 'qc_module': 'occ'}}, id='mp2 rohf df: 2 occ'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'df', 'qc_module': 'dfmp2'}}, id='mp2 rohf df: 2 dfmp2*'),
    pytest.param({'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'cd', 'qc_module': 'occ'}}, id='mp2 rohf cd: 1 occ*'),
])
def test_mp2_module(inp, clsd_open_pmols, request):
    tnm = request.node.name
    subject = clsd_open_pmols[inp['subject']]

    psi4.set_options({
        'basis': 'cc-pvdz',
        'scf_type': 'df',
        'guess': 'sad',
        'freeze_core': 'True',
        'e_convergence': 8,
        'd_convergence': 7,
    })

    psi4.set_options(inp['options'])
    ene, wfn = psi4.energy('mp2', molecule=subject, return_wfn=True)

    ref_ref = _ref_module[inp['options']['reference']]['HF TOTAL ENERGY']
    ref_corl = _ref_module[inp['options']['reference']][inp['options']['mp2_type']]['MP2 CORRELATION ENERGY']
    ref_tot = _ref_module[inp['options']['reference']][inp['options']['mp2_type']]['MP2 TOTAL ENERGY']

    for obj in [psi4.core, wfn]:
        for pv in ['HF TOTAL ENERGY',
                   'SCF TOTAL ENERGY',
                   'CURRENT REFERENCE ENERGY']:
            assert compare_values(ref_ref, obj.variable(pv), 6, tnm + ' ' + pv)

        for pv in ['MP2 CORRELATION ENERGY',
                   'CURRENT CORRELATION ENERGY',]:
            assert compare_values(ref_corl, obj.variable(pv), 6, tnm + ' ' + pv)

        for pv in ['MP2 TOTAL ENERGY',
                   'CURRENT ENERGY',]:
            assert compare_values(ref_tot, obj.variable(pv), 6, tnm + ' ' + pv)

    assert compare_values(ref_tot, ene, 6, tnm + ' return')
    assert compare_values(ref_tot, wfn.energy(), 6, tnm + ' wfn')
