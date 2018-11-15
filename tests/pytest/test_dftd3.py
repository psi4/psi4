import copy

import pytest
import qcelemental as qcel

from utils import *
from addons import *

import psi4
from psi4.driver import qcdb
from psi4.driver.qcdb import intf_dftd3

seneyne = """
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
"""

db3lypd3bj = {
    'dashlevel': 'd3bj',
    'dashparams': {
        's8': 1.9889,
        's6': 1.0,
        'a2': 4.4211,
        'a1': 0.3981
    },
    'dashparams_citation': '',
    'fctldash': 'b3lyp-d3(bj)'
}
db3lypd3bjcustom = copy.deepcopy(db3lypd3bj)
db3lypd3bjcustom['fctldash'] = ''
db3lypd3bjcustom['dashparams']['a2'] = 5.4211

dpbed3zero = {
    'dashlevel': 'd3zero',
    'dashparams': {
        's6': 1.0,
        's8': 0.722,
        'sr6': 1.217,
        'sr8': 1.0,
        'alpha6': 14.0
    },
    'dashparams_citation': '',
    'fctldash': 'pbe-d3'
}

@pytest.mark.parametrize("inp,expected", [
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3LYP', 'level_hint': 'D3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'param_tweaks': {'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'verbose': 3, 'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 5.4211}}, ''), db3lypd3bjcustom),
    (({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'pbe', 'level_hint': 'd3zero'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe', 'level_hint': 'd3'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe-d3'}, 'PBE-D3'), dpbed3zero),
])
def test_intf_dftd3_from_arrays(inp, expected):
    res = intf_dftd3.from_arrays(**inp[0])
    assert compare_dicts(expected, res, 4, tnm())
    assert compare_strings(inp[1], compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(expected, res, 4, tnm() + ' idempotent')


@pytest.mark.parametrize("inp", [
    ({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a3': 5.4211}}),
    ({'name_hint': 'fakeb3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211, 'zzz': 0.0}}),
    ({'name_hint': 'asdf-d4'}),
])
def test_intf_dftd3_from_arrays_error(inp):
    with pytest.raises(qcdb.ValidationError):
        intf_dftd3.from_arrays(**inp)



def test_intf_dftd3_from_arrays_supplement():
    ans = {'dashlevel': 'chg', 'dashparams': {'s6': 4.05}, 'fctldash': 'asdf-d4', 'dashparams_citation': '    mypaper\n'}
    supp = {'chg': {'definitions': {'asdf-d4': {'params': {'s6': 4.05}, 'citation': '    mypaper\n'}}}}

    res = intf_dftd3.from_arrays(name_hint='asdf-d4', level_hint='chg', dashcoeff_supplement=supp)
    print(res)
    assert compare_dicts(ans, res, 4, tnm())
    with pytest.raises(qcdb.ValidationError):
        intf_dftd3.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'],
        level_hint=res['dashlevel'],
        param_tweaks=res['dashparams'],
        dashcoeff_supplement=supp)
    assert compare_dicts(ans, res, 4, tnm() + ' idempotent')


def test_3():
    sys = qcel.molparse.from_string(seneyne)['qm']

    res = intf_dftd3.run_dftd3_from_arrays(molrec=sys, name_hint='b3lyp', level_hint='d3bj')
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res['options']), 'key')


def compute_key(pjrec):
    #return '-'.join([pjrec['functional'], pjrec['dashlevel']]).upper()
    return pjrec['fctldash'].upper()


# """dftd3/energy"""
# ! Exercises the various DFT-D corrections, both through python directly and through c++

ref_eneyne = {}
dmm = ['dimer', 'mA', 'mB']
ref_eneyne['B3LYP-D2'] = dict(zip(dmm, [-0.00390110, -0.00165271, -0.00058118]))
ref_eneyne['B3LYP-D3'] = dict(zip(dmm, [-0.00285088, -0.00084340, -0.00031923]))
ref_eneyne['B3LYP-D3(BJ)'] = dict(zip(dmm, [-0.00784595, -0.00394347, -0.00226683]))
ref_eneyne['PBE-D2'] = dict(zip(dmm, [-0.00278650, -0.00118051, -0.00041513]))
ref_eneyne['PBE-D3'] = dict(zip(dmm, [-0.00175474, -0.00045421, -0.00016839]))
ref_eneyne['PBE-D3(BJ)'] = dict(zip(dmm, [-0.00475937, -0.00235265, -0.00131239]))


def test_10_qmol():
    #eneyne = qcdb.set_molecule(seneyne)
    #eneyne.update_geometry()
    eneyne = qcdb.Molecule(seneyne)

    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    assert compare_values(ref_eneyne['B3LYP-D2']['dimer'], E, 7, 'Q: Ethene-Ethyne -D2')


@using_qcdb
def test_10_pmol():
    import psi4
    eneyne = psi4.set_molecule(seneyne)
    eneyne.update_geometry()

    E, G = eneyne.run_dftd3('b3lyp', 'd2gr')
    assert compare_values(ref_eneyne['B3LYP-D2']['dimer'], E, 7, 'P: Ethene-Ethyne -D2')


@using_qcdb
def test_11_energy():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, jrec = qcdb.energy('d3-b3lyp-d2', return_wfn=True)
    assert compare_values(ref_eneyne['B3LYP-D2']['dimer'], E, 7, 'P: Ethene-Ethyne -D2')
    assert compare_values(ref_eneyne['B3LYP-D2']['dimer'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(ref_eneyne['B3LYP-D2']['dimer'], jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm())


_dftd3_paramset = [
    ({'name': 'd3-b3lyp-d', 'subject': 'dimer', 'lbl': 'B3LYP-D2'}),
    ({'name': 'd3-b3lyp-d3bj', 'subject': 'mA', 'lbl': 'B3LYP-D3(BJ)'}),
    ({'name': 'd3-PBE-D3zero', 'subject': 'mB', 'lbl': 'PBE-D3'}),
]

@pytest.mark.parametrize("inp", _dftd3_paramset)
def test_11_dftd3_qmol(inp):
    eneyne = qcdb.Molecule(seneyne)
    subject = {'dimer': eneyne,
               'mA': eneyne.extract_subsets(1),
               'mB': eneyne.extract_subsets(2)}
    expected = ref_eneyne[inp['lbl']][inp['subject']]

    jrec = qcdb.intf_dftd3.run_dftd3(inp['name'], subject[inp['subject']], options={}, ptype='energy')
    assert compare_values(expected, jrec['qcvars']['CURRENT ENERGY'].data, 7, tnm())
    assert compare_values(expected, jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(expected, jrec['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'].data, 7, tnm())



@using_psi4
@pytest.mark.parametrize("inp", _dftd3_paramset)
def test_11_dftd3_pmol(inp):
    eneyne = psi4.core.Molecule.from_string(seneyne)
    subject = {'dimer': eneyne,
               'mA': eneyne.extract_subsets(1),
               'mB': eneyne.extract_subsets(2)}
    expected = ref_eneyne[inp['lbl']][inp['subject']]

    jrec = qcdb.intf_dftd3.run_dftd3(inp['name'], subject[inp['subject']], options={}, ptype='energy')
    assert compare_values(expected, jrec['qcvars']['CURRENT ENERGY'].data, 7, tnm())
    assert compare_values(expected, jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(expected, jrec['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'].data, 7, tnm())

    #assert compare_values(expected, psi4.get_variable('DISPERSION CORRECTION ENERGY'), 7, tnm())
    #assert compare_values(expected, psi4.get_variable(inp['lbl'] + ' DISPERSION CORRECTION ENERGY'), 7, tnm())


@using_qcdb
def test_11_b():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()
    mA = eneyne.extract_subsets(1)
    mB = eneyne.extract_subsets(2)

    E, jrec = qcdb.energy('d3-b3lyp-d3bj', return_wfn=True, molecule=mA)
    assert compare_values(ref_eneyne['B3LYP-D3(BJ)']['mA'], E, 7, tnm())
    assert compare_values(ref_eneyne['B3LYP-D3(BJ)']['mA'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(ref_eneyne['B3LYP-D3(BJ)']['mA'], jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm())
