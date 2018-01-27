import pytest
from utils import *

import numpy as np

import qcdb

# system-shorthand   tot-chg, frag-chg, tot-mult, frag-mult         expected final tot/frag chg/mult


def test_validate_and_fill_chgmult_1():
    chgmult_tester(['He', 0, [0], 1, [1], (0, [0], 1, [1])])


def test_validate_and_fill_chgmult_2():
    chgmult_tester(['He', None, [None], None, [None], (0, [0], 1, [1])])


def test_validate_and_fill_chgmult_3():
    chgmult_tester(['He/He', None, [None, None], None, [None, None], (0, [0, 0], 1, [1, 1])])


def test_validate_and_fill_chgmult_4():
    chgmult_tester(['He/He', 2, [None, None], None, [None, None], (2, [2, 0], 1, [1, 1])])


def test_validate_and_fill_chgmult_5():
    chgmult_tester(['He/He', None, [2, None], None, [None, None], (2, [2, 0], 1, [1, 1])])


def test_validate_and_fill_chgmult_6():
    chgmult_tester(['He/He', 0, [2, None], None, [None, None], (0, [2, -2], 1, [1, 1])])


def test_validate_and_fill_chgmult_7():
    chgmult_tester(['Ne/He/He', -2, [None, 2, None], None, [None, None, None], (-2, [-4, 2, 0], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_8():
    chgmult_tester(['Ne/He/He', 2, [None, -2, None], None, [None, None, None], (2, [4, -2, 0], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_9():
    chgmult_tester(['He/He/Ne', 2, [None, -2, None], None, [None, None, None], (2, [0, -2, 4], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_10():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He/He/Ne', 2, [None, -2, 0], None, [None, None, None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_11():
    chgmult_tester(['He/He/Ne', 2, [2, -2, None], None, [None, None, None], (2, [2, -2, 2], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_12():
    chgmult_tester(['He/He', None, [-2, 2], None, [None, None], (0, [-2, 2], 1, [1, 1])])


def test_validate_and_fill_chgmult_13():
    chgmult_tester(['He/He', None, [None, -2], None, [None, None], (-2, [0, -2], 1, [1, 1])])


def test_validate_and_fill_chgmult_14():
    chgmult_tester(['Ne/Ne', 0, [None, 4], None, [None, None], (0, [-4, 4], 1, [1, 1])])


def test_validate_and_fill_chgmult_15():
    chgmult_tester(['He/He/He', 4, [2, None, None], None, [None, None, None], (4, [2, 2, 0], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_16():
    chgmult_tester(['He/He', 0, [-2, 2], None, [None, None], (0, [-2, 2], 1, [1, 1])])


def test_validate_and_fill_chgmult_17():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He/He', 0, [-2, -2], None, [None, None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_18():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He', None, [None], 0, [None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_19():
    chgmult_tester(['He', None, [None], None, [1], (0, [0], 1, [1])])


def test_validate_and_fill_chgmult_20():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He', None, [None], None, [2], 'Irreconcilable'])


def test_validate_and_fill_chgmult_21():
    chgmult_tester(['He', None, [None], None, [3], (0, [0], 3, [3])])


def test_validate_and_fill_chgmult_22():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He', None, [None], None, [5], 'Irreconcilable'])


def test_validate_and_fill_chgmult_23():
    chgmult_tester(['He', None, [-1], None, [2], (-1, [-1], 2, [2])])


def test_validate_and_fill_chgmult_24():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He', None, [-2], None, [2], 'Irreconcilable'])


def test_validate_and_fill_chgmult_25():
    chgmult_tester(['He/He', None, [None, None], None, [1, 1], (0, [0, 0], 1, [1, 1])])


def test_validate_and_fill_chgmult_26():
    chgmult_tester(['He/He', None, [None, None], None, [3, 1], (0, [0, 0], 3, [3, 1])])


def test_validate_and_fill_chgmult_27():
    chgmult_tester(['He/He', None, [None, None], None, [1, 3], (0, [0, 0], 3, [1, 3])])


def test_validate_and_fill_chgmult_28():
    chgmult_tester(['He/He', None, [None, None], None, [3, 3], (0, [0, 0], 5, [3, 3])])


def test_validate_and_fill_chgmult_29():
    chgmult_tester(['He/He', None, [None, None], 3, [3, 3], (0, [0, 0], 3, [3, 3])])


def test_validate_and_fill_chgmult_30():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['He/He', None, [None, None], 2, [3, 3], 'Irreconcilable'])


def test_validate_and_fill_chgmult_31():
    chgmult_tester(['H', None, [None], None, [None], (0, [0], 2, [2])])


def test_validate_and_fill_chgmult_32():
    chgmult_tester(['H', 1, [None], None, [None], (1, [1], 1, [1])])


def test_validate_and_fill_chgmult_33():
    chgmult_tester(['H', None, [-1], None, [None], (-1, [-1], 1, [1])])


def test_validate_and_fill_chgmult_34():
    chgmult_tester(['funnyH', None, [None], None, [None], (0, [0], 1, [1])])


def test_validate_and_fill_chgmult_35():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['funnierH', None, [None], None, [None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_36():
    chgmult_tester(['H/H', None, [None, None], None, [None, None], (0, [0, 0], 3, [2, 2])])


def test_validate_and_fill_chgmult_37():
    chgmult_tester(['H/He', None, [None, None], None, [None, None], (0, [0, 0], 2, [2, 1])])


def test_validate_and_fill_chgmult_38():
    chgmult_tester(['H/He', None, [1, 1], None, [None, None], (2, [1, 1], 2, [1, 2])])


def test_validate_and_fill_chgmult_39():
    chgmult_tester(['H/He', -2, [-1, None], None, [None, None], (-2, [-1, -1], 2, [1, 2])])


def test_validate_and_fill_chgmult_40():
    chgmult_tester(
        ['H/He/Na/Ne', None, [1, None, 1, None], None, [None, None, None, None], (2, [1, 0, 1, 0], 1, [1, 1, 1, 1])])


def test_validate_and_fill_chgmult_41():
    chgmult_tester(
        ['H/He/Na/Ne', None, [-1, None, 1, None], None, [None, None, None, None], (0, [-1, 0, 1, 0], 1, [1, 1, 1, 1])])


def test_validate_and_fill_chgmult_42():
    chgmult_tester(
        ['H/He/Na/Ne', 2, [None, None, 1, None], None, [None, None, None, None], (2, [1, 0, 1, 0], 1, [1, 1, 1, 1])])


def test_validate_and_fill_chgmult_43():
    chgmult_tester(
        ['H/He/Na/Ne', 3, [None, None, 1, None], None, [None, None, None, None], (3, [0, 2, 1, 0], 2, [2, 1, 1, 1])])


def test_validate_and_fill_chgmult_44():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['H/He', None, [1, None], None, [2, None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_45():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['H/He', None, [None, 0], None, [None, 2], 'Irreconcilable'])


def test_validate_and_fill_chgmult_46():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['H/He', None, [None, -1], None, [None, 3], 'Irreconcilable'])


def test_validate_and_fill_chgmult_47():
    chgmult_tester(
        ['H/He/Na/Ne', None, [None, 1, 0, 1], None, [None, None, None, None], (2, [0, 1, 0, 1], 5, [2, 2, 2, 2])])


def test_validate_and_fill_chgmult_48():
    chgmult_tester(
        ['H/He/Na/Ne', None, [None, 1, 0, None], None, [None, None, None, None], (1, [0, 1, 0, 0], 4, [2, 2, 2, 1])])


def test_validate_and_fill_chgmult_49():
    chgmult_tester(
        ['H/He/Na/Ne', None, [None, 1, 0, None], None, [None, None, 4, None], (1, [0, 1, 0, 0], 6, [2, 2, 4, 1])])


def test_validate_and_fill_chgmult_50():
    chgmult_tester(['He/He/He', 0, [None, None, 1], None, [1, None, 2], (0, [0, -1, 1], 3, [1, 2, 2])])


def test_validate_and_fill_chgmult_51():
    chgmult_tester(['N/N/N', None, [1, 1, 1], 3, [None, 3, None], (3, [1, 1, 1], 3, [1, 3, 1])])


def test_validate_and_fill_chgmult_52():
    chgmult_tester(['N/N/N', None, [1, 1, 1], 3, [None, None, None], (3, [1, 1, 1], 3, [3, 1, 1])])


def test_validate_and_fill_chgmult_53():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['N/N/N', None, [None, None, None], 3, [None, None, 2], 'Irreconcilable'])


def test_validate_and_fill_chgmult_54():
    chgmult_tester(['N/N/N', 1, [None, -1, None], 3, [None, None, 2], (1, [2, -1, 0], 3, [2, 1, 2])])


def test_validate_and_fill_chgmult_55():
    chgmult_tester(['N/Ne/N', 1, [None, None, None], 4, [None, 3, None], (1, [1, 0, 0], 4, [1, 3, 2])])


def test_validate_and_fill_chgmult_56():
    chgmult_tester(['N/Ne/N', None, [None, None, 1], 4, [None, 3, None], (1, [0, 0, 1], 4, [2, 3, 1])])


def test_validate_and_fill_chgmult_57():
    chgmult_tester(['He/He', None, [-1, 1], None, [None, None], (0, [-1, 1], 3, [2, 2])])


def test_validate_and_fill_chgmult_58():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['Gh', 1, [None], None, [None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_59():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['Gh', -1, [None], None, [None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_60():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['Gh', None, [None], 3, [None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_61():
    chgmult_tester(['He/Gh', None, [2, None], None, [None, None], (2, [2, 0], 1, [1, 1])])


def test_validate_and_fill_chgmult_62():
    with pytest.raises(qcdb.ValidationError):
        chgmult_tester(['Gh/He', None, [2, None], None, [None, None], 'Irreconcilable'])


def test_validate_and_fill_chgmult_63():
    chgmult_tester(['Gh/He/Ne', 2, [None, -2, None], None, [None, None, None], (2, [0, -2, 4], 1, [1, 1, 1])])


def test_validate_and_fill_chgmult_64():
    chgmult_tester(['Gh/He/Gh', 1, [None, None, None], None, [None, None, None], (1, [0, 1, 0], 2, [1, 2, 1])])


# Notes
#  9 - residual +4 distributes to first fragment able to wholly accept it (He+4 is no-go)
# 10 - residual +4 unsuited for only open fragment, He, so irreconcilable
# 11 - non-positive multiplicity
# 20 - doublet non consistent with closed-shell, neutral default charge
# 22 - insufficient electrons for pentuplet
# 24 - doublet not consistent with even charge
# 30 - bad parity btwn mult and total # electrons
# 35 - insufficient electrons
# 55 - both (1, (1, 0.0, 0.0), 4, (1, 3, 2)) and (1, (0.0, 0.0, 1), 4, (2, 3, 1)) plausible

_keys = ['molecular_charge', 'fragment_charges', 'molecular_multiplicity', 'fragment_multiplicities']
_systemtranslator = {
    'He': (np.array([2]), np.array([])),
    'He/He': (np.array([2, 2]), np.array([1])),
    'Ne/He/He': (np.array([10, 2, 2]), np.array([1, 2])),
    'He/He/Ne': (np.array([2, 2, 10]), np.array([1, 2])),
    'Ne/Ne': (np.array([10, 10]), np.array([1])),
    'He/He/He': (np.array([2, 2, 2]), np.array([1, 2])),
    'H': (np.array([1]), np.array([])),
    'funnyH': (np.array([0]), np.array([])),  # has no electrons
    'funnierH': (np.array([-1]), np.array([])),  # has positron
    'H/H': (np.array([1, 1]), np.array([1])),
    'H/He': (np.array([1, 2]), np.array([1])),
    'H/He/Na/Ne': (np.array([1, 2, 11, 10]), np.array([1, 2, 3])),
    'N/N/N': (np.array([7, 7, 7]), np.array([1, 2])),
    'N/Ne/N': (np.array([7, 10, 7]), np.array([1, 2])),
    'He/Gh': (np.array([2, 0]), np.array([1])),
    'Gh/He': (np.array([0, 2]), np.array([1])),
    'Gh': (np.array([0, 0]), np.array([])),
    'Gh/He/Ne': (np.array([0, 0, 2, 10]), np.array([2, 3])),
    'Gh/He/Gh': (np.array([0, 2, 0]), np.array([1, 2])),
}


def chgmult_tester(test):
    system = _systemtranslator[test[0]]
    ans = qcdb.molparse.validate_and_fill_chgmult(system[0], system[1], test[1], test[2], test[3], test[4], verbose=0)
    assert compare_integers(1, ans == dict(zip(_keys, test[5])), """{}: {}, {}, {}, {} --> {}""".format(*test))
