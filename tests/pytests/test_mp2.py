import pytest

from utils import *

import psi4

pytestmark = [pytest.mark.quick, pytest.mark.mp2]


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
}  # yapf: disable
for mp2type in ['df', 'conv']:
    _ref_h2o_ccpvdz[mp2type]['SCF TOTAL ENERGY'] = _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']
    _ref_h2o_ccpvdz[mp2type]['CURRENT REFERENCE ENERGY'] = _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']
    _ref_h2o_ccpvdz[mp2type]['5050SCS-MP2 CORRELATION ENERGY'] = (
        0.5 * (_ref_h2o_ccpvdz[mp2type]['MP2 SAME-SPIN CORRELATION ENERGY'] +
               _ref_h2o_ccpvdz[mp2type]['MP2 OPPOSITE-SPIN CORRELATION ENERGY']))
    _ref_h2o_ccpvdz[mp2type]['5050SCS-MP2 TOTAL ENERGY'] = _ref_h2o_ccpvdz[mp2type][
        '5050SCS-MP2 CORRELATION ENERGY'] + _ref_h2o_ccpvdz[mp2type]['HF TOTAL ENERGY']


@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'df'}}, id='mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'conv'}}, id='mp2 (conv)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'df', 'mp2_os_scale': 1.2, 'mp2_ss_scale': 0.33333333}}, id='explicit scs mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': 'SCS-MP2', 'options': {'mp2_type': 'conv', 'mp2_os_scale': 1.2, 'mp2_ss_scale': 0.33333333}}, id='explicit scs mp2 (conv)'),
    pytest.param({'name': 'Mp2', 'custom': '5050SCS-MP2', 'options': {'mp2_type': 'df', 'mp2_os_scale': 0.5, 'mp2_ss_scale': 0.5}}, id='user-def scs mp2 (df)'),
    pytest.param({'name': 'Mp2', 'custom': '5050SCS-MP2', 'options': {'mp2_type': 'conv', 'mp2_os_scale': 0.5, 'mp2_ss_scale': 0.5}}, id='user-def scs mp2 (conv)'),
])  # yapf: disable
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
                'HF TOTAL ENERGY', 'SCF TOTAL ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY',
                'MP2 OPPOSITE-SPIN CORRELATION ENERGY', 'MP2 CORRELATION ENERGY', 'MP2 TOTAL ENERGY',
                'SCS-MP2 CORRELATION ENERGY', 'SCS-MP2 TOTAL ENERGY', 'CURRENT REFERENCE ENERGY'
        ]:

            assert compare_values(ref_block[pv], obj.variable(pv), 5, pv)

        assert compare_values(ref_custom_corl, obj.variable('CUSTOM SCS-MP2 CORRELATION ENERGY'), 5,
                              'custom scsmp2 corl')
        assert compare_values(ref_custom_tot, obj.variable('CUSTOM SCS-MP2 TOTAL ENERGY'), 5, 'custom scsmp2 ')

        assert compare_values(ref_corl, obj.variable('CURRENT CORRELATION ENERGY'), 5, 'current corl')
        assert compare_values(ref_tot, obj.variable('CURRENT ENERGY'), 5, 'current')

    assert compare_values(ref_tot, ene, 5, 'return')
    assert compare_values(ref_tot, wfn.energy(), 5, 'wfn')


@pytest.fixture
def clsd_open_pmols():
    mols = {
        'hf':
        psi4.core.Molecule.from_string("""
                H
                F 1 0.917
              """),
        'bh_h2p':
        psi4.core.Molecule.from_string("""
                1 2
                B     0.10369114     0.00000000     0.00000000
                H    -1.13269886     0.00000000     0.00000000
                H     3.00000000     0.37149000     0.00000000
                H     3.00000000    -0.37149000     0.00000000
              """),
    }
    return mols


# yapf: disable
_ref_module = {scftype: {ref: {frz: {mp2type: {} for mp2type in ['conv', 'df', 'cd']} for frz in ['true', 'false']} for ref in ['rhf', 'uhf', 'rohf']} for scftype in ['pk', 'df']}

_ref_module['df']['rhf']['HF TOTAL ENERGY'] =                        -100.019400605629
_ref_module['df']['uhf']['HF TOTAL ENERGY'] =                         -25.945130559147
_ref_module['df']['rohf']['HF TOTAL ENERGY'] =                        -25.943606522029
_ref_module['pk']['rhf']['HF TOTAL ENERGY'] =                        -100.01941126902270
_ref_module['pk']['uhf']['HF TOTAL ENERGY'] =                         -25.94513842869638
#_ref_module['pk']['rohf']['HF TOTAL ENERGY'] =

# <<<  scf DF, fc  >>>
_ref_module['df']['rhf']['true']['conv']['MP2 CORRELATION ENERGY'] =   -0.201612517228
_ref_module['df']['rhf']['true']['conv']['MP2 TOTAL ENERGY'] =       -100.221013122857
_ref_module['df']['rhf']['true']['df']['MP2 CORRELATION ENERGY'] =     -0.201610660387
_ref_module['df']['rhf']['true']['df']['MP2 TOTAL ENERGY'] =         -100.221011266016
_ref_module['df']['rhf']['true']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # dfmp2 findif-5 fc df+df
    [[    0.00000000000000,    0.00000000000000,    0.00314716362539],
     [    0.00000000000000,    0.00000000000000,   -0.00314716362539]])
_ref_module['df']['rhf']['true']['cd']['MP2 CORRELATION ENERGY'] =     -0.201609396752
_ref_module['df']['rhf']['true']['cd']['MP2 TOTAL ENERGY'] =         -100.221010002381

_ref_module['df']['uhf']['true']['conv']['MP2 CORRELATION ENERGY'] =   -0.058421122206
_ref_module['df']['uhf']['true']['conv']['MP2 TOTAL ENERGY'] =        -26.003551681354
_ref_module['df']['uhf']['true']['df']['MP2 CORRELATION ENERGY'] =     -0.058390006825
_ref_module['df']['uhf']['true']['df']['MP2 TOTAL ENERGY'] =          -26.003520565972
_ref_module['df']['uhf']['true']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # dfmp2 findif-5 fc df+df
    [[    0.00000000000000,     0.00000000000000,     0.01231996225662],
     [    0.00000000000000,     0.00000000000000,    -0.01186374280678],
     [    0.00000000000000,     0.01031743020277,    -0.00022810972492],
     [    0.00000000000000,    -0.01031743020277,    -0.00022810972492]])
_ref_module['df']['uhf']['true']['cd']['MP2 CORRELATION ENERGY'] =     -0.058409837177
_ref_module['df']['uhf']['true']['cd']['MP2 TOTAL ENERGY'] =          -26.003540396324

_ref_module['df']['rohf']['true']['conv']['MP2 CORRELATION ENERGY'] =  -0.060939211739
_ref_module['df']['rohf']['true']['conv']['MP2 TOTAL ENERGY'] =       -26.004545733768
_ref_module['df']['rohf']['true']['df']['MP2 CORRELATION ENERGY'] =    -0.059372748391
_ref_module['df']['rohf']['true']['df']['MP2 TOTAL ENERGY'] =         -26.002979270420
#_ref_module['df']['rohf']['true']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(
_ref_module['df']['rohf']['true']['cd']['MP2 CORRELATION ENERGY'] =    -0.059393510962
_ref_module['df']['rohf']['true']['cd']['MP2 TOTAL ENERGY'] =         -26.003000032991

# <<<  scf DF, nfc  >>>
#_ref_module['df']['rhf']['false']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['rhf']['false']['conv']['MP2 TOTAL ENERGY'] =
_ref_module['df']['rhf']['false']['df']['MP2 CORRELATION ENERGY'] =    -0.2037649370559149
_ref_module['df']['rhf']['false']['df']['MP2 TOTAL ENERGY'] =        -100.2231655426856776
_ref_module['df']['rhf']['false']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # dfmp2 findif-5 nfc df+df
    [[    0.00000000000000,     0.00000000000000,     0.00279211492833],
     [    0.00000000000000,     0.00000000000000,    -0.00279211492833]])
#_ref_module['df']['rhf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['rhf']['false']['cd']['MP2 TOTAL ENERGY'] =

#_ref_module['df']['uhf']['false']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['uhf']['false']['conv']['MP2 TOTAL ENERGY'] =
_ref_module['df']['uhf']['false']['df']['MP2 CORRELATION ENERGY'] =    -0.0594557966607590
_ref_module['df']['uhf']['false']['df']['MP2 TOTAL ENERGY'] =         -26.0045863558097601
_ref_module['df']['uhf']['false']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # dfmp2 findif-5 nfc df+df
    [[     0.00000000000000,     0.00000000000000,     0.01252024755551],
     [     0.00000000000000,     0.00000000000000,    -0.01207773525598],
     [     0.00000000000000,     0.01032204616770,    -0.00022125614977],
     [     0.00000000000000,    -0.01032204616770,    -0.00022125614977]])
#_ref_module['df']['uhf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['uhf']['false']['cd']['MP2 TOTAL ENERGY'] =

#_ref_module['df']['rohf']['false']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['rohf']['false']['conv']['MP2 TOTAL ENERGY'] =
_ref_module['df']['rohf']['false']['df']['MP2 CORRELATION ENERGY'] =   -0.0604436327328384
_ref_module['df']['rohf']['false']['df']['MP2 TOTAL ENERGY'] =        -26.0040501547626377
_ref_module['df']['rohf']['false']['df']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # occ findif-5 nfc df+df
    [[     0.00000000000000,     0.00000000000000,     0.01361287313486],
     [     0.00000000000000,     0.00000000000000,    -0.01314329502424],
     [     0.00000000000000,     0.01029838165151,    -0.00023478905531],
     [     0.00000000000000,    -0.01029838165151,    -0.00023478905531]])
#_ref_module['df']['rohf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['df']['rohf']['false']['cd']['MP2 TOTAL ENERGY'] =

# <<<  scf CONV, fc  >>>
_ref_module['pk']['rhf']['true']['conv']['MP2 CORRELATION ENERGY'] =   -0.201627516796
_ref_module['pk']['rhf']['true']['conv']['MP2 TOTAL ENERGY'] =       -100.221038785818
_ref_module['pk']['rhf']['true']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # fnocc findif-5 fc pk+conv
    [[     0.00000000000000,     0.00000000000000,     0.00317450456474],
     [     0.00000000000000,     0.00000000000000,    -0.00317450456474]])
#_ref_module['pk']['rhf']['true']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rhf']['true']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rhf']['true']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rhf']['true']['cd']['MP2 TOTAL ENERGY'] =

#_ref_module['pk']['uhf']['true']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['uhf']['true']['conv']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['uhf']['true']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(
#_ref_module['pk']['uhf']['true']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['uhf']['true']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['uhf']['true']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['uhf']['true']['cd']['MP2 TOTAL ENERGY'] =

#_ref_module['pk']['rohf']['true']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['true']['conv']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rohf']['true']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(
#_ref_module['pk']['rohf']['true']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['true']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rohf']['true']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['true']['cd']['MP2 TOTAL ENERGY'] =

# <<<  scf CONV, nfc  >>>
_ref_module['pk']['rhf']['false']['conv']['MP2 CORRELATION ENERGY'] =  -0.203781911950
_ref_module['pk']['rhf']['false']['conv']['MP2 TOTAL ENERGY'] =      -100.223193180973
_ref_module['pk']['rhf']['false']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # fnocc findif-5 nfc pk+conv
    [[   0.0000000000,       0.0000000000,       0.0028193375],
     [   0.0000000000,       0.0000000000,      -0.0028193375]])
#_ref_module['pk']['rhf']['false']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rhf']['false']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rhf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rhf']['false']['cd']['MP2 TOTAL ENERGY'] =

_ref_module['pk']['uhf']['false']['conv']['MP2 CORRELATION ENERGY'] =  -0.05948928003552
_ref_module['pk']['uhf']['false']['conv']['MP2 TOTAL ENERGY'] =       -26.00462770873190
_ref_module['pk']['uhf']['false']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(  # occ findif-5 nfc pk+conv
    [[     0.00000000000000,     0.00000000000000,     0.01250561195911],
     [     0.00000000000000,     0.00000000000000,    -0.01206536529299],
     [     0.00000000000000,     0.01033165380573,    -0.00022012333306],
     [     0.00000000000000,    -0.01033165380573,    -0.00022012333306]])
#_ref_module['pk']['uhf']['false']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['uhf']['false']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['uhf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['uhf']['false']['cd']['MP2 TOTAL ENERGY'] =

#_ref_module['pk']['rohf']['false']['conv']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['false']['conv']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rohf']['false']['conv']['MP2 TOTAL GRADIENT'] = psi4.core.Matrix.from_list(
#_ref_module['pk']['rohf']['false']['df']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['false']['df']['MP2 TOTAL ENERGY'] =
#_ref_module['pk']['rohf']['false']['cd']['MP2 CORRELATION ENERGY'] =
#_ref_module['pk']['rohf']['false']['cd']['MP2 TOTAL ENERGY'] =
# yapf: enable


@pytest.mark.parametrize("inp", [
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf conv: * occ'),
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'fnocc', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf conv:   fnocc'),
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'detci', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf conv:   detci'),
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf df:     occ'),
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'dfmp2', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf df:   * dfmp2'),
    pytest.param({'driver': 'energy',   'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'cd',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  rhf cd:   * occ'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  uhf conv: * occ'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  uhf df:     occ'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'df',   'qc_module': 'dfmp2', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  uhf df:   * dfmp2'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'cd',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2  uhf cd:   * occ'),
    #pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 rohf conv: * occ'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'detci', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 rohf conv:   detci'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 rohf df:     occ'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'df',   'qc_module': 'dfmp2', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 rohf df:   * dfmp2'),
    pytest.param({'driver': 'energy',   'subject': 'bh_h2p', 'options': {'reference': 'rohf', 'mp2_type': 'cd',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 rohf cd:   * occ'),

    # scf_type set for indexing, not for driver
    #pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'pk'}}, id='mp2 grad  rhf conv  fc: * occ'),
    pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'false', 'scf_type': 'pk'}}, id='mp2 grad  rhf conv nfc: * occ'),
    pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 grad  rhf df    fc:   occ'),
    pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'dfmp2', 'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 grad  rhf df    fc: * dfmp2'),
    pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'false', 'scf_type': 'df'}}, id='mp2 grad  rhf df   nfc:   occ'),
    pytest.param({'driver': 'gradient', 'subject': 'hf',     'options': {'reference': 'rhf',  'mp2_type': 'df',   'qc_module': 'dfmp2', 'freeze_core': 'false', 'scf_type': 'df'}}, id='mp2 grad  rhf df   nfc: * dfmp2'),
    pytest.param({'driver': 'gradient', 'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'false', 'scf_type': 'pk'}}, id='mp2 grad  uhf conv nfc: * occ'),
    pytest.param({'driver': 'gradient', 'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'df'}}, id='mp2 grad  uhf df    fc: * occ'),
    pytest.param({'driver': 'gradient', 'subject': 'bh_h2p', 'options': {'reference': 'uhf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'false', 'scf_type': 'df'}}, id='mp2 grad  uhf df   nfc: * occ'),
    # For gradients, this method would be found in the procedures table but would return a ManagedMethodError from proc.py. Normally, this would confuse the analytic-or-findif logic in gradient().
    #   This scenario is now managed, and this test case makes sure its routing stays managed.
    pytest.param({'driver': 'gradient', 'subject': 'bh_h2p', 'options': {'reference': 'rohf',  'mp2_type': 'df',   'qc_module': 'occ',   'freeze_core': 'false', 'scf_type': 'df', 'points': 5}}, id='mp2 grad rohf df   nfc: findif'),
])  # yapf: disable
def test_mp2_module(inp, clsd_open_pmols, request):
    tnm = request.node.name
    subject = clsd_open_pmols[inp['subject']]

    ref_block = _ref_module[inp['options']['scf_type']][inp['options']['reference']]
    ref_ref = ref_block['HF TOTAL ENERGY']
    ref_block = ref_block[inp['options']['freeze_core']][inp['options']['mp2_type']]
    ref_corl = ref_block['MP2 CORRELATION ENERGY']
    ref_tot = ref_block['MP2 TOTAL ENERGY']

    psi4.set_options({
        'basis': 'cc-pvdz',
        'guess': 'sad',
        'e_convergence': 8,
        'd_convergence': 7,
    })
    psi4.set_options(inp['options'])

    if inp['driver'] == 'energy':
        ene, wfn = psi4.energy('mp2', molecule=subject, return_wfn=True)

    elif inp['driver'] == 'gradient':
        grad, wfn = psi4.gradient('mp2', molecule=subject, return_wfn=True)
        ref_tot_grad = ref_block['MP2 TOTAL GRADIENT']

    for obj in [psi4.core, wfn]:
        for pv in ['HF TOTAL ENERGY', 'SCF TOTAL ENERGY', 'CURRENT REFERENCE ENERGY']:
            assert compare_values(ref_ref, obj.variable(pv), 6, tnm + ' ' + pv)

        for pv in [
                'MP2 CORRELATION ENERGY',
                'CURRENT CORRELATION ENERGY',
        ]:
            assert compare_values(ref_corl, obj.variable(pv), 6, tnm + ' ' + pv)

        for pv in [
                'MP2 TOTAL ENERGY',
                'CURRENT ENERGY',
        ]:
            assert compare_values(ref_tot, obj.variable(pv), 6, tnm + ' ' + pv)

    assert compare_values(ref_tot, wfn.energy(), 6, tnm + ' wfn')

    if inp['driver'] == 'energy':
        assert compare_values(ref_tot, ene, 6, tnm + ' return')

    elif inp['driver'] == 'gradient':
        assert compare_matrices(ref_tot_grad, wfn.gradient(), 6, tnm + ' grad wfn')
        assert compare_matrices(ref_tot_grad, grad, 6, tnm + ' grad return')
