from utils import *
import qcdb


def test_parse_nucleus_label_1():
    label_tester('@ca_miNe', {'E': 'ca', 'Z': None, 'user': '_miNe', 'A': None, 'real': False, 'mass': None})


def test_parse_nucleus_label_2():
    label_tester('Gh(Ca_mine)', {'E': 'Ca', 'Z': None, 'user': '_mine', 'A': None, 'real': False, 'mass': None})


def test_parse_nucleus_label_3():
    label_tester('@Ca_mine@1.07', {'E': 'Ca', 'Z': None, 'user': '_mine', 'A': None, 'real': False, 'mass': 1.07})


def test_parse_nucleus_label_4():
    label_tester('Gh(cA_MINE@1.07)', {'E': 'cA', 'Z': None, 'user': '_MINE', 'A': None, 'real': False, 'mass': 1.07})


def test_parse_nucleus_label_5():
    label_tester('@40Ca_mine@1.07', {'E': 'Ca', 'Z': None, 'user': '_mine', 'A': 40, 'real': False, 'mass': 1.07})


def test_parse_nucleus_label_6():
    label_tester('Gh(40Ca_mine@1.07)', {'E': 'Ca', 'Z': None, 'user': '_mine', 'A': 40, 'real': False, 'mass': 1.07})


def test_parse_nucleus_label_7():
    label_tester('444lu333@4.0', {'E': 'lu', 'Z': None, 'user': '333', 'A': 444, 'real': True, 'mass': 4.0})


def test_parse_nucleus_label_8():
    label_tester('@444lu333@4.4', {'E': 'lu', 'Z': None, 'user': '333', 'A': 444, 'real': False, 'mass': 4.4})


def test_parse_nucleus_label_9():
    label_tester('8i', {'E': 'i', 'Z': None, 'user': None, 'A': 8, 'real': True, 'mass': None})


def test_parse_nucleus_label_10():
    label_tester('53_mI4', {'Z': 53, 'E': None, 'user': '_mI4', 'A': None, 'real': True, 'mass': None})


def test_parse_nucleus_label_11():
    label_tester('@5_MINEs3@4.4', {'Z': 5, 'E': None, 'user': '_MINEs3', 'A': None, 'real': False, 'mass': 4.4})


def test_parse_nucleus_label_12():
    label_tester('Gh(555_mines3@0.1)', {'Z': 555, 'E': None, 'user': '_mines3', 'A': None, 'real': False, 'mass': 0.1})


def label_tester(label, ans):
    lbl_A, lbl_Z, lbl_E, lbl_mass, lbl_real, lbl_user = qcdb.molparse.nucleus.parse_nucleus_label(label)

    assert compare_integers(ans['real'], lbl_real, label + " real")
    assert compare_integers(ans['A'], lbl_A, label + " A")
    assert compare_integers(ans['Z'], lbl_Z, label + " Z")
    assert compare_strings(ans['E'], lbl_E, label + " symbol")
    assert compare_strings(ans['user'], lbl_user, label + " user")
    assert compare_values(ans['mass'], lbl_mass, 6, label + " mass", passnone=True)
