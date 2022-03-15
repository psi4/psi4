import sys

import pytest
from utils import *
from addons import *

subject1 = """
3 au

Co 0 0 0
H  2 0 0
h_OTher -2 0 0
"""

ans1_au = """3 au
CoH2
Co                    0.000000000000     0.000000000000     0.000000000000
H                     2.000000000000     0.000000000000     0.000000000000
H                    -2.000000000000    -0.000000000000     0.000000000000
"""

ans1_ang = """3
CoH2
Co                    0.000000000000     0.000000000000     0.000000000000
H                     1.058354421340     0.000000000000     0.000000000000
H                    -1.058354421340    -0.000000000000     0.000000000000
"""

ans1c_ang = """3
CoH2
59Co                      0.00000000         0.00000000         0.00000000
1H                        1.05835442         0.00000000         0.00000000
1H_other                 -1.05835442        -0.00000000         0.00000000
"""

#subject2 = """
#Co 0 0 0
#units au
#no_reorient
#--
#@H  2 0 0
#h_OTher -2 0 0
#"""
#
#ans2_au = """3 au
#
#Co                    0.000000000000     0.000000000000     0.000000000000
#@H                    2.000000000000     0.000000000000     0.000000000000
#H                    -2.000000000000     0.000000000000     0.000000000000"""
#
#ans2_ang = """3
#
#Co                    0.000000000000     0.000000000000     0.000000000000
#Gh(1)                 1.058354417180     0.000000000000     0.000000000000
#H                    -1.058354417180     0.000000000000     0.000000000000"""
#
#ans2c_ang = """2
#
#Co                    0.000000000000     0.000000000000     0.000000000000
#H                    -1.058354417180     0.000000000000     0.000000000000"""


subject2 = """
Co 0 0 0
no_reorient
--
@H  1.05835442134 0 0
h_OTher -1.05835442134 0 0
"""

ans2_au = """3 au
CoH2
Co                    0.000000000000     0.000000000000     0.000000000000
@H                    2.000000000000     0.000000000000     0.000000000000
H                    -2.000000000000     0.000000000000     0.000000000000
"""

ans2_ang = """3
CoH2
Co                    0.000000000000     0.000000000000     0.000000000000
Gh(1)                 1.058354421340     0.000000000000     0.000000000000
H                    -1.058354421340     0.000000000000     0.000000000000
"""

ans2c_ang = """2
CoH2
Co                    0.000000000000     0.000000000000     0.000000000000
H                    -1.058354421340     0.000000000000     0.000000000000
"""


def test_toxyz_1a():
    subject = subject1
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', units='Bohr')

    assert compare_strings(ans1_au, xyz, sys._getframe().f_code.co_name)

def test_toxyz_1b():
    subject = subject1
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', units='Angstrom')

    assert compare_strings(ans1_ang, xyz, sys._getframe().f_code.co_name)

def test_toxyz_1c():
    subject = subject1
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', prec=8, atom_format='{elea}{elem}{elbl}')
    print(xyz)

    assert compare_strings(ans1c_ang, xyz, sys._getframe().f_code.co_name)

#def test_toxyz_2a():
#    subject = subject2
#    mol = qcdb.Molecule(subject)
#
#    xyz = mol.to_string(dtype='xyz', units='Bohr')
#
#    assert compare_strings(ans2_au, xyz, sys._getframe().f_code.co_name)
#
#def test_toxyz_2b():
#    subject = subject2
#    mol = qcdb.Molecule(subject)
#
#    xyz = mol.to_string(dtype='xyz', units='Angstrom', ghost_format='Gh({elez})')
#
#    assert compare_strings(ans2_ang, xyz, sys._getframe().f_code.co_name)
#
#def test_toxyz_2c():
#    subject = subject2
#    mol = qcdb.Molecule(subject)
#
#    xyz = mol.to_string(dtype='xyz', units='Angstrom', ghost_format='')
#
#    assert compare_strings(ans2c_ang, xyz, sys._getframe().f_code.co_name)

def test_toxyz_2a():
    subject = subject2
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', units='Bohr')

    assert compare_strings(ans2_au, xyz, sys._getframe().f_code.co_name)

def test_toxyz_2b():
    subject = subject2
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', units='Angstrom', ghost_format='Gh({elez})')

    assert compare_strings(ans2_ang, xyz, sys._getframe().f_code.co_name)

def test_toxyz_2c():
    subject = subject2
    mol = qcdb.Molecule(subject)

    xyz = mol.to_string(dtype='xyz', units='Angstrom', ghost_format='')

    assert compare_strings(ans2c_ang, xyz, sys._getframe().f_code.co_name)

@using_psi4_molrec
def test_toxyz_3a():
    import psi4

    subject = subject2
    mol = psi4.core.Molecule.from_string(subject)

    xyz = mol.to_string(dtype='xyz', units='Bohr')

    assert compare_strings(ans2_au, xyz, sys._getframe().f_code.co_name)

