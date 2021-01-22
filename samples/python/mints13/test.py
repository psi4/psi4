#! test fragment decomposition + to/from_dict

import numpy as np
import psi4
from psi4.driver import qcdb

psi4.set_output_file("output.dat", False)

def test_chgmult(expected, cgmpdict, label):
    rc, rfc, rm, rfm = expected
    qcdb.compare_integers(rc, cgmpdict['molecular_charge'], label + ': c')
    qcdb.compare_integers(rm, cgmpdict['molecular_multiplicity'], label + ': m')
    qcdb.compare_integers(True, np.allclose(cgmpdict['fragment_charges'], rfc), label + ': fc')
    qcdb.compare_integers(True, np.allclose(cgmpdict['fragment_multiplicities'], rfm), label + ': fm')


def test_dimer(mol, expected_cgmp, label, mtype):

    mol.update_geometry()
    dAB = mol.to_dict()
    test_chgmult(expected_cgmp['AB'], dAB, label + ' AB')
    mAB = mtype.from_dict(dAB)
    qcdb.compare_molrecs(dAB, mAB.to_dict(), label + ' AB roundtrip', atol=1.e-6)
    
    aB = mol.extract_subsets(2, 1)
    daB = aB.to_dict()
    test_chgmult(expected_cgmp['aB'], daB, label + ' aB')
    maB = mtype.from_dict(daB)
    qcdb.compare_molrecs(daB, maB.to_dict(), label + ' aB roundtrip', atol=1.e-6)
    
    Ab = mol.extract_subsets(1, 2)
    dAb = Ab.to_dict()
    test_chgmult(expected_cgmp['Ab'], dAb, label + ' Ab')
    mAb = mtype.from_dict(dAb)
    qcdb.compare_molrecs(dAb, mAb.to_dict(), label + ' Ab roundtrip', atol=1.e-6)
    
    A_ = mol.extract_subsets(1)
    dA_ = A_.to_dict()
    test_chgmult(expected_cgmp['A_'], dA_, label + ' A_')
    mA_ = mtype.from_dict(dA_)
    qcdb.compare_molrecs(dA_, mA_.to_dict(), label + ' A_ roundtrip', atol=1.e-6)
    
    _B = mol.extract_subsets(2)
    d_B = _B.to_dict()
    test_chgmult(expected_cgmp['_B'], d_B, label + ' _B')
    m_B = mtype.from_dict(d_B)
    qcdb.compare_molrecs(d_B, m_B.to_dict(), label + ' _B roundtrip', atol=1.e-6)

    qcdb.compare_integers(True, type(mol) == mtype, label + ': AB type')
    qcdb.compare_integers(True, type(Ab) == mtype, label + ': Ab type')



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

eneyne_cgmp = {
    'AB': (0, [0, 0], 1, [1, 1]),
    'aB': (0, [0, 0], 1, [1, 1]),
    'Ab': (0, [0, 0], 1, [1, 1]),
    'A_': (0, [0], 1, [1]),
    '_B': (0, [0], 1, [1]),
}

negpos = """
-1 1
O 0.0 0.0 0.0
H 0.0 0.0 1.0
--
1 1
O 2.0 2.0 2.0
H 3.0 2.0 2.0
H 2.0 3.0 2.0
H 2.0 2.0 3.0
"""

negpos_cgmp = {
    'AB': (0, [-1, 1], 1, [1, 1]),
    'A_': (-1, [-1], 1, [1]),
    '_B': (1, [1], 1, [1]),
    'Ab': (-1, [-1, 0], 1, [1, 1]),
    'aB': (1, [0, 1], 1, [1, 1]),
}

qeneyne = qcdb.Molecule(eneyne)
peneyne = psi4.geometry(eneyne)
qnegpos = qcdb.Molecule(negpos)
pnegpos = psi4.geometry(negpos)

test_dimer(qeneyne, eneyne_cgmp, 'Q: eneyne', qcdb.Molecule)
test_dimer(peneyne, eneyne_cgmp, 'P: eneyne', psi4.core.Molecule)
test_dimer(qnegpos, negpos_cgmp, 'Q: negpos', qcdb.Molecule)
test_dimer(pnegpos, negpos_cgmp, 'P: negpos', psi4.core.Molecule)

# Once user starts messing with cgmp other than in construction, user has
#   no way to mess with fragment cgmp, and Psi/QCDB Molecule classes don't do
#   much to set things in order. Upon to_dict, things get sorted into some
#   physical reality, but fragment charges in a complicated system like this
#   won't get sorted out to resemble thier initial state (could do more
#   try/catch, but that's really the class's job). So really all that can be
#   tested in the main dimer's total charge and total mult.

qnegpos.set_multiplicity(3)
qnegpos.set_molecular_charge(2)

qresetAB = qnegpos.to_dict()
qcdb.compare_integers(2, qresetAB['molecular_charge'], 'Q: reset-negpos: c')
qcdb.compare_integers(3, qresetAB['molecular_multiplicity'], 'Q: reset-negpos: m')

pnegpos.set_multiplicity(3)
pnegpos.set_molecular_charge(2)

presetAB = pnegpos.to_dict()
qcdb.compare_integers(2, presetAB['molecular_charge'], 'P: reset-negpos: c')
qcdb.compare_integers(3, presetAB['molecular_multiplicity'], 'P: reset-negpos: m')
