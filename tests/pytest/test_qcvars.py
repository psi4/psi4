import copy

import pytest
from utils import *

import numpy as np
import psi4

_vars_entered = {
    'VAR A': 4.0,
    'VaR B': -4.0,
    'MATVAR A': psi4.core.Matrix.from_array(np.arange(6).reshape(2, 3)),
    'MatvaR B': psi4.core.Matrix.from_array(np.arange(3).reshape(1, 3)),
    'NPVAR A': np.arange(8).reshape(2, 4),
    'NpvaR B': np.arange(4).reshape(1, 4),
}

#_vars_stored = {k.upper(): v for k, v in _vars_entered.items()}
_vars_stored = {
    'VAR A': 4.0,
    'VAR B': -4.0,
    'MATVAR A': psi4.core.Matrix.from_array(np.arange(6).reshape(2, 3)),
    'MATVAR B': psi4.core.Matrix.from_array(np.arange(3).reshape(1, 3)),
    'NPVAR A': psi4.core.Matrix.from_array(np.arange(8).reshape(2, 4)),
    'NPVAR B': psi4.core.Matrix.from_array(np.arange(4).reshape(1, 4)),
}


@pytest.fixture
def pe_qcvars():
    psi4.core.clean_variables()
    for pv, pvv in _vars_entered.items():
        psi4.set_variable(pv, pvv)


@pytest.fixture
def wfn_qcvars():
    he = psi4.geometry('He')
    wfn = psi4.core.Wavefunction.build(he, 'cc-pvdz')

    for pv, pvv in _vars_entered.items():
        wfn.set_variable(pv, pvv)
    return wfn


# <<<  Deprecated/Replacement pair tests  >>>

## TODO delete in Psi4 v1.4
#def test_core_get_variables(pe_qcvars):
#    subject = psi4.core.get_variables()
#
#    assert compare_dicts(_vars_stored, subject, 8, tnm())


def test_core_variables(pe_qcvars):

    # can't use compare_dicts with symmetry psi4.Matrix
    def compare_qcvars(ref, expected, decimal, label):
        assert set(ref.keys()) == set(expected.keys())
        for k, v in ref.items():
            if isinstance(v, psi4.core.Matrix):
                compare_matrices(v, expected[k], decimal, label)
            else:
                compare_values(v, expected[k], decimal, label)

    subject = psi4.core.variables()
    compare_qcvars(_vars_stored, subject, 8, '')

    psi4.core.set_variable('npvar A', np.zeros(3).reshape(1, 3))
    compare_qcvars(_vars_stored, subject, 8, '')


## TODO delete in Psi4 v1.4
#def test_core_get_variable(pe_qcvars):
#    subject = psi4.core.get_variable('vAR B')
#
#    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


def test_core_set_variable_overwrite(pe_qcvars):
    # fine to overwrite keys
    key = 'var D'
    val = 3.3
    val2 = 4.4
    psi4.core.set_variable(key, val)
    assert compare_values(val, psi4.core.variable(key), 8, tnm())
    psi4.core.set_variable(key, val2)
    assert compare_values(val2, psi4.core.variable(key), 8, tnm())

    # fine to overwrite array keys
    key = 'matvar D'
    mat = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    mat2 = psi4.core.Matrix.from_array(np.arange(6).reshape(3, 2))
    psi4.core.set_variable(key, mat)
    assert compare_matrices(mat, psi4.core.variable(key), 8, tnm())
    psi4.core.set_variable(key, mat2)
    assert compare_matrices(mat2, psi4.core.variable(key), 8, tnm())

    # not fine to shadow keys with both types
    with pytest.raises(psi4.ValidationError) as err:
        psi4.core.set_variable('vAr D', mat)
    assert 'already a scalar variable' in str(err)

    with pytest.raises(psi4.ValidationError) as err:
        psi4.core.set_variable('matvAr D', val)
    assert 'already an array variable' in str(err)


def test_core_variable_none(pe_qcvars):
    with pytest.raises(KeyError):
        psi4.core.variable('var f')


def test_core_variable_scal(pe_qcvars):
    assert compare_values(_vars_stored['VAR B'], psi4.core.variable('vAR B'), 8, tnm())


def test_core_variable_mat(pe_qcvars):
    assert compare_matrices(_vars_stored['MATVAR B'], psi4.core.variable('MatvAR B'), 8, tnm())


def test_core_variable_np(pe_qcvars):
    assert compare_matrices(_vars_stored['NPVAR B'], psi4.core.variable('NpvAR B'), 8, tnm())


def test_core_variable_mem_scal(pe_qcvars):
    key = 'VaR C'
    ref = 3.3
    val = 3.3
    psi4.core.set_variable(key, val)

    assert compare_values(ref, val, 8, tnm())
    assert compare_values(ref, psi4.core.variable(key), 8, tnm())

    val *= 2
    assert compare_values(ref, psi4.core.variable(key), 8, tnm())

    accessed = psi4.core.variable(key)
    accessed *= 3
    assert compare_values(ref, psi4.core.variable(key), 8, tnm())


def test_core_variable_mem_mat(pe_qcvars):
    key = 'MaTvAr C'
    ref = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    val = psi4.core.Matrix.from_array(np.arange(4).reshape(2, 2))
    psi4.core.set_variable(key, val)

    assert compare_matrices(ref, val, 8, tnm())
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())

    val.scale(2)
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())

    accessed = psi4.core.variable(key)
    accessed.scale(3)
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())


def test_core_variable_mem_np(pe_qcvars):
    key = 'npVaR C'
    ref = np.arange(4).reshape(2, 2)
    val = np.arange(4).reshape(2, 2)
    psi4.core.set_variable(key, val)

    assert compare_arrays(ref, val, 8, tnm())
    ref = psi4.core.Matrix.from_array(ref)
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())

    val *= 2
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())

    accessed = psi4.core.variable(key)
    accessed.scale(3)
    assert compare_matrices(ref, psi4.core.variable(key), 8, tnm())


## TODO delete in Psi4 v1.4
#def test_wfn_get_variable(wfn_qcvars):
#    subject = wfn_qcvars.get_variable('vAR B')
#
#    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


#def test_wfn_variable(wfn_qcvars):
#    subject = wfn_qcvars.variable('vAR B')
#
#    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


# <<<  Ordinary tests  >>>

@pytest.mark.parametrize("tkey,fkey", [
    pytest.param('var A', 'var C', id='scal'),
    pytest.param('matvar A', 'var C', id='mat'),
    pytest.param('npvar A', 'var C', id='np'),
])
def test_core_has_del_variable_scal(tkey, fkey, pe_qcvars):
    assert psi4.core.has_variable(tkey)
    assert not psi4.core.has_variable(fkey)

    psi4.core.del_variable(tkey)
    assert not psi4.core.has_variable(tkey)

    psi4.core.del_variable(fkey)


#def test_wfn_has_del_variable(wfn_qcvars):
#    subject_t = wfn_qcvars.has_variable('VAR a')
#    subject_f = wfn_qcvars.has_variable('VAR c')
#
#    assert compare_integers(True, subject_t, tnm())
#    assert compare_integers(False, subject_f, tnm())
#
#    wfn_qcvars.del_variable('var A')
#    subject_d = wfn_qcvars.has_variable('VAR a')
#
#    assert compare_integers(False, subject_d, tnm())
#
#    wfn_qcvars.del_variable('var C')
