import copy

import pytest
from utils import *

import numpy as np
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


_vars_entered = {
    'VAR A': 4.0,
    'VaR B': -4.0,
    'MATVAR A': psi4.core.Matrix.from_array(np.arange(6).reshape(2, 3)),
    'MatvaR B': psi4.core.Matrix.from_array(np.arange(3).reshape(1, 3)),
    'NPVAR A': np.arange(8).reshape(2, 4),
    'NpvaR B': np.arange(4).reshape(1, 4),
}

_vars_stored = {
    'VAR A': 4.0,
    'VAR B': -4.0,
    'MATVAR A': psi4.core.Matrix.from_array(np.arange(6).reshape(2, 3)),
    'MATVAR B': psi4.core.Matrix.from_array(np.arange(3).reshape(1, 3)),
    'NPVAR A': psi4.core.Matrix.from_array(np.arange(8).reshape(2, 4)),
    'NPVAR B': psi4.core.Matrix.from_array(np.arange(4).reshape(1, 4)),
}


@pytest.fixture
def pe_wfn_qcvars():
    psi4.core.clean_variables()

    he = psi4.geometry('He')
    wfn = psi4.core.Wavefunction.build(he, 'cc-pvdz')

    for pv, pvv in _vars_entered.items():
        psi4.core.set_variable(pv, pvv)
        wfn.set_variable(pv, pvv)

    return wfn


# can't use compare_dicts with symmetry psi4.Matrix
def _compare_qcvars(ref, expected, decimal, label):
    assert set(ref.keys()) == set(expected.keys())
    for k, v in ref.items():
        if isinstance(v, psi4.core.Matrix):
            assert compare_matrices(v, expected[k], decimal, label)
        else:
            assert compare_values(v, expected[k], decimal, label)


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_variables(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    subject = obj.variables()
    _compare_qcvars(_vars_stored, subject, 8, '')

    obj.set_variable('npvar A', np.zeros(3).reshape(1, 3))
    _compare_qcvars(_vars_stored, subject, 8, '')


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_set_variable_overwrite(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    # fine to overwrite keys
    key = 'var D'
    val = 3.3
    val2 = 4.4
    obj.set_variable(key, val)
    assert compare_values(val, obj.variable(key), 8, tnm())
    obj.set_variable(key, val2)
    assert compare_values(val2, obj.variable(key), 8, tnm())

    # fine to overwrite array keys
    key = 'matvar D'
    mat = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    mat2 = psi4.core.Matrix.from_array(np.arange(6).reshape(3, 2))
    obj.set_variable(key, mat)
    assert compare_matrices(mat, obj.variable(key), 8, tnm())
    obj.set_variable(key, mat2)
    assert compare_matrices(mat2, obj.variable(key), 8, tnm())

    # not fine to shadow keys with both types
    with pytest.raises(psi4.ValidationError) as err:
        obj.set_variable('vAr D', mat)
    assert 'already a scalar variable' in str(err.value)

    with pytest.raises(psi4.ValidationError) as err:
        obj.set_variable('matvAr D', val)
    assert 'already an array variable' in str(err.value)


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_variable_none(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    with pytest.raises(KeyError):
        obj.variable('var f')


@pytest.mark.parametrize("mode,key", [
    pytest.param('globals', 'vAR B', id='globals scal'),
    pytest.param('globals', 'MatvAR B', id='globals mat'),
    pytest.param('globals', 'NpvAR B', id='globals np'),
    pytest.param('wfn', 'vAR B', id='wfn scal'),
    pytest.param('wfn', 'MatvAR B', id='wfn mat'),
    pytest.param('wfn', 'NpvAR B', id='wfn np'),
])
def test_variable(mode, key, pe_wfn_qcvars, request):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]
    if 'scal' in request.node.name:
        compare = compare_values
    else:
        compare = compare_matrices

    assert compare(_vars_stored[key.upper()], obj.variable(key), 8, tnm())


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_variable_mem_scal(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    key = 'VaR C'
    ref = 3.3
    val = 3.3
    obj.set_variable(key, val)

    assert compare_values(ref, val, 8, tnm())
    assert compare_values(ref, obj.variable(key), 8, tnm())

    val *= 2
    assert compare_values(ref, obj.variable(key), 8, tnm())

    accessed = obj.variable(key)
    accessed *= 3
    assert compare_values(ref, obj.variable(key), 8, tnm())


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_variable_mem_mat(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    key = 'MaTvAr C'
    ref = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    val = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    obj.set_variable(key, val)

    assert compare_matrices(ref, val, 8, tnm())
    assert compare_matrices(ref, obj.variable(key), 8, tnm())

    val.scale(2)
    assert compare_matrices(ref, obj.variable(key), 8, tnm())

    accessed = obj.variable(key)
    accessed.scale(3)
    assert compare_matrices(ref, obj.variable(key), 8, tnm())


@pytest.mark.parametrize("mode", [
    ("globals"),
    ("wfn"),
])
def test_variable_mem_np(mode, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    key = 'npVaR C'
    ref = np.arange(4).reshape(2, 2)
    val = np.arange(4).reshape(2, 2)
    obj.set_variable(key, val)

    assert compare_arrays(ref, val, 8, tnm())
    ref = psi4.core.Matrix.from_array(ref)
    assert compare_matrices(ref, obj.variable(key), 8, tnm())

    val *= 2
    assert compare_matrices(ref, obj.variable(key), 8, tnm())

    accessed = obj.variable(key)
    accessed.scale(3)
    assert compare_matrices(ref, obj.variable(key), 8, tnm())


@pytest.mark.parametrize("mode,tkey,fkey", [
    pytest.param('globals', 'var A', 'var C', id='globals scal'),
    pytest.param('globals', 'matvar A', 'var C', id='globals mat'),
    pytest.param('globals', 'npvar A', 'var C', id='globals np'),
    pytest.param('wfn', 'var A', 'var C', id='wfn scal'),
    pytest.param('wfn', 'matvar A', 'var C', id='wfn mat'),
    pytest.param('wfn', 'npvar A', 'var C', id='wfn np'),
])
def test_has_del_variable_scal(mode, tkey, fkey, pe_wfn_qcvars):
    obj = {'globals': psi4.core, 'wfn': pe_wfn_qcvars}[mode]

    assert obj.has_variable(tkey)
    assert not obj.has_variable(fkey)

    obj.del_variable(tkey)
    assert not obj.has_variable(tkey)

    obj.del_variable(fkey)


# <<<  TODO Deprecated! Delete in Psi4 v1.4  >>>


def test_deprecated_core_get_variable(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = psi4.core.get_variable('vAR B')

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_core_get_variables(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = psi4.core.get_variables()

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_core_get_array_variable(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = psi4.core.get_array_variable('MatvAR B')

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_core_get_array_variables(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = psi4.core.get_array_variables()

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_wfn_get_variable(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = pe_wfn_qcvars.get_variable('vAR B')

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_wfn_get_array(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = pe_wfn_qcvars.get_array('MatvAR B')

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_wfn_set_array(pe_wfn_qcvars):
    mat = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    with pytest.raises(psi4.UpgradeHelper) as err:
        pe_wfn_qcvars.set_array('matvar D', mat)

    assert 'is obsolete as of 1.9' in str(err.value)


def test_deprecated_wfn_arrays(pe_wfn_qcvars):
    with pytest.raises(psi4.UpgradeHelper) as err:
        subject = pe_wfn_qcvars.arrays()

    assert 'is obsolete as of 1.9' in str(err.value)
