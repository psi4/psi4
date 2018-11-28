import pytest
from utils import *

import psi4

_vars_entered = {'VAR A': 4.0, 'VaR B': -4.0}

_vars_stored = {k.upper(): v for k, v in _vars_entered.items()}


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

# TODO delete in Psi4 v1.4
def test_core_get_variables(pe_qcvars):
    subject = psi4.core.get_variables()

    assert compare_dicts(_vars_stored, subject, 8, tnm())


def test_core_variables(pe_qcvars):
    subject = psi4.core.variables()

    assert compare_dicts(_vars_stored, subject, 8, tnm())


# TODO delete in Psi4 v1.4
def test_core_get_variable(pe_qcvars):
    subject = psi4.core.get_variable('vAR B')

    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


def test_core_variable(pe_qcvars):
    subject = psi4.core.variable('vAR B')

    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


# TODO delete in Psi4 v1.4
def test_wfn_get_variable(wfn_qcvars):
    subject = wfn_qcvars.get_variable('vAR B')

    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


def test_wfn_variable(wfn_qcvars):
    subject = wfn_qcvars.variable('vAR B')

    assert compare_values(_vars_stored['VAR B'], subject, 8, tnm())


# <<<  Ordinary tests  >>>

def test_core_has_del_variable(pe_qcvars):
    subject_t = psi4.core.has_variable('VAR a')
    subject_f = psi4.core.has_variable('VAR c')

    assert compare_integers(True, subject_t, tnm())
    assert compare_integers(False, subject_f, tnm())

    psi4.core.del_variable('var A')
    subject_d = psi4.core.has_variable('VAR a')

    assert compare_integers(False, subject_d, tnm())

    psi4.core.del_variable('var C')


def test_wfn_has_del_variable(wfn_qcvars):
    subject_t = wfn_qcvars.has_variable('VAR a')
    subject_f = wfn_qcvars.has_variable('VAR c')

    assert compare_integers(True, subject_t, tnm())
    assert compare_integers(False, subject_f, tnm())

    wfn_qcvars.del_variable('var A')
    subject_d = wfn_qcvars.has_variable('VAR a')

    assert compare_integers(False, subject_d, tnm())

    wfn_qcvars.del_variable('var C')
