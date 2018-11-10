import sys
import copy
import pprint

import pytest

import qcelemental as qcel

from utils import *

import qcdb
from qcdb import intf_dftd3

eneyne = """
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


def test_recon_1a():
    res = intf_dftd3.from_arrays(name_hint='b3lyp', level_hint='d3bj')
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1b():
    res = intf_dftd3.from_arrays(name_hint='b3LYP', level_hint='D3bj')
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1c():
    res = intf_dftd3.from_arrays(param_tweaks={'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, level_hint='d3bj')
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1d():
    res = intf_dftd3.from_arrays(name_hint='b3lyp', level_hint='d3bJ', param_tweaks={'a2': 4.4211})
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1e():
    ans = copy.deepcopy(db3lypd3bj)
    ans['fctldash'] = ''
    ans['dashparams']['a2'] = 5.4211

    res = intf_dftd3.from_arrays(verbose=3, name_hint='b3lyp', level_hint='d3bJ', param_tweaks={'a2': 5.4211})
    assert compare_dicts(ans, res, 4, sys._getframe().f_code.co_name)
    #assert compare_strings('-D3BJ', compute_key(res), 'key')
    assert compare_strings('', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        verbose=3, name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(ans, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1f():
    with pytest.raises(qcdb.ValidationError):
        intf_dftd3.from_arrays(name_hint='b3lyp', level_hint='d3bJ', param_tweaks={'a3': 5.4211})


def test_recon_1g():
    with pytest.raises(qcdb.ValidationError):
        intf_dftd3.from_arrays(name_hint='fakeb3lyp', level_hint='d3bJ', param_tweaks={'s6': 5.4211})


def test_recon_1h():
    with pytest.raises(qcdb.ValidationError):
        res = intf_dftd3.from_arrays(level_hint='d3bJ', param_tweaks={'s6': 5.4211})


def test_recon_1i():
    res = intf_dftd3.from_arrays(name_hint='b3lyp-d3bj', param_tweaks={'a2': 4.4211})
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_1j():
    with pytest.raises(qcdb.ValidationError):
        res = intf_dftd3.from_arrays(name_hint='b3lyp-d3bj', param_tweaks={'a2': 4.4211, 'zzz': 0.0})


def test_recon_2a():
    res = intf_dftd3.from_arrays(name_hint='pbe', level_hint='d3zero')
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name)
    #assert compare_strings('PBE-D3ZERO', compute_key(res), 'key')
    assert compare_strings('PBE-D3', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_2b():
    res = intf_dftd3.from_arrays(name_hint='pbe', level_hint='d3')
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name)
    assert compare_strings('PBE-D3', compute_key(res), 'key')
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_2c():
    res = intf_dftd3.from_arrays(name_hint='pbe-d3')
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name)
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_recon_2d():
    with pytest.raises(qcdb.ValidationError):
        res = intf_dftd3.from_arrays(name_hint='asdf-d4')


def test_recon_2e():
    ans = {'dashlevel': 'chg', 'dashparams': {'s6': 4.05}, 'fctldash': 'asdf-d4', 'dashparams_citation': '    mypaper\n'}
    supp = {'chg': {'definitions': {'asdf-d4': {'params': {'s6': 4.05}, 'citation': '    mypaper\n'}}}}

    res = intf_dftd3.from_arrays(name_hint='asdf-d4', level_hint='chg', dashcoeff_supplement=supp)
    print(res)
    assert compare_dicts(ans, res, 4, sys._getframe().f_code.co_name)
    with pytest.raises(qcdb.ValidationError):
        intf_dftd3.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    res = intf_dftd3.from_arrays(
        name_hint=res['fctldash'],
        level_hint=res['dashlevel'],
        param_tweaks=res['dashparams'],
        dashcoeff_supplement=supp)
    assert compare_dicts(ans, res, 4, sys._getframe().f_code.co_name + ' idempotent')


def test_3():
    sys = qcel.molparse.from_string(eneyne)['qm']

    res = intf_dftd3.run_dftd3_from_arrays(molrec=sys, name_hint='b3lyp', level_hint='d3bj')
    assert compare_strings('B3LYP-D3(BJ)', compute_key(res['options']), 'key')


def compute_key(pjrec):
    #return '-'.join([pjrec['functional'], pjrec['dashlevel']]).upper()
    return pjrec['fctldash'].upper()


# """dftd3/energy"""
# ! Exercises the various DFT-D corrections, both through python directly and through c++

ref_d2 = [-0.00390110, -0.00165271, -0.00058118]
ref_d3zero = [-0.00285088, -0.00084340, -0.00031923]
ref_d3bj = [-0.00784595, -0.00394347, -0.00226683]

ref_pbe_d2 = [-0.00278650, -0.00118051, -0.00041513]
ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]
ref_pbe_d3bj = [-0.00475937, -0.00235265, -0.00131239]

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


def test_10_qmol():
    #eneyne = qcdb.set_molecule(seneyne)
    #eneyne.update_geometry()
    eneyne = qcdb.Molecule(seneyne)

    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    assert compare_values(ref_d2[0], E, 7, 'Q: Ethene-Ethyne -D2')


def hide_test_10_pmol():
    import psi4
    eneyne = psi4.set_molecule(seneyne)
    eneyne.update_geometry()

    E, G = eneyne.run_dftd3('b3lyp', 'd2gr')
    assert compare_values(ref_d2[0], E, 7, 'P: Ethene-Ethyne -D2')


def hide_test_11_energy():
    tnm = sys._getframe().f_code.co_name
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, jrec = qcdb.energy('d3-b3lyp-d2', return_wfn=True)
    assert compare_values(ref_d2[0], E, 7, 'P: Ethene-Ethyne -D2')
    assert compare_values(ref_d2[0], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm)
    assert compare_values(ref_d2[0], jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm)


def test_11_a_energy():
    tnm = sys._getframe().f_code.co_name
    eneyne = qcdb.Molecule(seneyne)

    jrec = qcdb.intf_dftd3.run_dftd3('d3-b3lyp-d', eneyne, {}, ptype='energy')

    assert compare_values(ref_d2[0], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm)
    assert compare_values(ref_d2[0], jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm)


def hide_test_11_b():
    tnm = sys._getframe().f_code.co_name
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()
    mA = eneyne.extract_subsets(1)
    mB = eneyne.extract_subsets(2)

    E, jrec = qcdb.energy('d3-b3lyp-d3bj', return_wfn=True, molecule=mA)
    assert compare_values(ref_d3bj[1], E, 7, tnm)
    assert compare_values(ref_d3bj[1], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm)
    assert compare_values(ref_d3bj[1], jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm)


def test_11_c():
    tnm = sys._getframe().f_code.co_name
    eneyne = qcdb.Molecule(seneyne)
    mA = eneyne.extract_subsets(1)
    mB = eneyne.extract_subsets(2)

    jrec = qcdb.intf_dftd3.run_dftd3('d3-b3lyp-d3bj', molecule=mA, options={}, ptype='energy')
    assert compare_values(ref_d3bj[1], jrec['qcvars']['CURRENT ENERGY'].data, 7, tnm)
    assert compare_values(ref_d3bj[1], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm)
    assert compare_values(ref_d3bj[1], jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm)
