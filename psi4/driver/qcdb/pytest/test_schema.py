import sys
import copy
import pprint

import pytest
import numpy as np

from utils import *
from addons import *

import qcdb
from qcdb.physconst import psi_bohr2angstroms


def test_1():
    ans = qcdb.Molecule.contiguize_from_fragment_pattern([[0], [1]],
                                                   geom=[0., 0., 0., 1., 0., 0.],
                                                   elbl=['O', 'H'])

    assert np.allclose(ans['geom'], np.array([0., 0., 0., 1., 0., 0.]))
    assert all(a == b for (a, b) in zip(ans['elbl'], ['O', 'H']))
    assert np.allclose(ans['fragment_separators'], np.array([]))


def test_2():
    ans = qcdb.Molecule.contiguize_from_fragment_pattern([[2, 0], [1]],
                                                   geom=np.array([[0., 0., 1.], [0., 0., 2.], [0., 0., 0.]]),
                                                   elem=np.array(['Li', 'H', 'He']))
    assert np.allclose(ans['geom'], np.array([0., 0., 0., 0., 0., 1., 0., 0., 2.]))
    assert all(a == b for (a, b) in zip(ans['elem'], ['He', 'Li', 'H']))
    assert np.allclose(ans['fragment_separators'], np.array([2]))


def test_3():
    ans = qcdb.Molecule.contiguize_from_fragment_pattern([[2, 0], [1]],
                                                         elez=[3, 1, 2])
    assert all(a == b for (a, b) in zip(ans['elez'], [2, 3, 1]))
    assert np.allclose(ans['fragment_separators'], np.array([2]))


def test_4():
    with pytest.raises(AssertionError):
        ans = qcdb.Molecule.contiguize_from_fragment_pattern([[2, 0], [1, 3]],
                                                              geom=np.array([[0., 0., 1.], [0., 0., 2.], [0., 0., 0.]]),
                                                              elem=np.array(['Li', 'H', 'He']))


def test_5():
    with pytest.raises(qcdb.ValidationError):
        ans = qcdb.Molecule.contiguize_from_fragment_pattern([[2, 0], [1, 4]])


def test_6():
    with pytest.raises(AssertionError):
        ans = qcdb.Molecule.contiguize_from_fragment_pattern([[2, 0], [1, 3]],
                                                              elem=np.array(['U', 'Li', 'H', 'He']),
                                                              elbl=np.array(['Li', 'H', 'He']))


#fullans1a = {'geom': np.array([ 0.,  0.,  0.,  1.,  0.,  0.]),
#             'elea': np.array([16,  1]),
#             'elez': np.array([8, 1]),
#             'elem': np.array(['O', 'H']),
#             'mass': np.array([ 15.99491462,   1.00782503]),
#             'real': np.array([ True,  True]),
#             'elbl': np.array(['', '']),
#             'units': 'Angstrom',
#             'fix_com': True,
#             'fix_orientation': False,
#             'fragment_separators': [],
#             'fragment_charges': [0.0],
#             'fragment_multiplicities': [2],
#             'molecular_charge': 0.0,
#             'molecular_multiplicity': 2,
#            }
