import pytest
from utils import *

import psi4


@pytest.fixture
def peglob():
    dicary = {'VAR A': 4.0,
              'VAR B': -4.0}

    psi4.core.clean_variables()
    for pv, pvv in dicary.items():
        psi4.set_variable(pv, pvv)

    return dicary


# TODO delete in Psi4 v1.4
def test_core_get_variables(peglob):
    subject = psi4.core.get_variables()

    assert compare_dicts(peglob, subject, 8, tnm())


def test_core_variables(peglob):
    subject = psi4.core.variables()

    assert compare_dicts(peglob, subject, 8, tnm())
