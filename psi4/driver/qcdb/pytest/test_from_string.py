import sys
import copy

import pytest
import numpy as np

from utils import *
import qcdb


subject1 = """O 0 0   0
no_com

H 1 ,, 0 \t  0 # stuff-n-nonsense"""

ans1 = {'geom': [0., 0., 0., 1., 0., 0.],
        'elbl': ['O', 'H'],
        'fix_com': True,
        'fragment_separators': [],
        'fragment_charges': [None],
        'fragment_multiplicities': [None],
       }

fullans1a = {'geom': np.array([ 0.,  0.,  0.,  1.,  0.,  0.]),
             'elea': np.array([16,  1]),
             'elez': np.array([8, 1]),
             'elem': np.array(['O', 'H']),
             'mass': np.array([ 15.99491462,   1.00782503]),
             'real': np.array([ True,  True]),
             'elbl': np.array(['', '']),
             'units': 'Angstrom', 
             'fix_com': True,
             'fix_orientation': False, 
             'fragment_separators': [],
             'fragment_charges': [0.0], 
             'fragment_multiplicities': [2],
             'molecular_charge': 0.0,
             'molecular_multiplicity': 2,
            }
fullans1c = copy.deepcopy(fullans1a)
fullans1c.update({'fragment_charges': [1.],
                  'fragment_multiplicities': [1],
                  'molecular_charge': 1.,
                  'molecular_multiplicity': 1})


def test_psi4_molstr1a():
    subject = subject1

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans1, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1a, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr1b():
    subject = '\n' + '\t' + subject1 + '\n\n'

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans1, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1a, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr1c():
    subject = '1 1\n  -- \n' + subject1
    ans = copy.deepcopy(ans1)
    ans.update({'molecular_charge': 1.,
                'molecular_multiplicity': 1})

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1c, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr1d():
    subject = subject1 + '\n1 1'
    ans = copy.deepcopy(ans1)
    ans.update({'fragment_charges': [1.],
                'fragment_multiplicities': [1]})

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1c, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr1e():
    """duplicate com"""
    subject = subject1 + '\n  nocom'

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)


subject2 = ["""
6Li 0.0 0.0 0.0
    units  a.u.
H_specIAL@2.014101  100 0 0""",
"""@Ne 2 4 6""",
"""h .0,1,2
Gh(he3) 0 1 3
noreorient"""]

ans2 = {'geom': [ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.],
        'elbl': ['6Li', 'H_specIAL@2.014101', '@Ne', 'h', 'Gh(he3)'],
        'units': 'Bohr',
        'fix_orientation': True,
        'fragment_separators': [2, 3],
        'fragment_charges': [None, None, None],
        'fragment_multiplicities': [None, None, None],
        }

fullans2 = {'geom': np.array([ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.]),
            'elea': np.array([6,  2, 20, 1, 4]),
            'elez': np.array([3, 1, 10, 1, 2]),
            'elem': np.array(['Li', 'H', 'Ne', 'H', 'He']),
            'mass': np.array([ 6.015122794, 2.014101, 19.99244017542, 1.00782503, 4.00260325415]),
            'real': np.array([ True,  True, False, True, False]),
            'elbl': np.array(['', '_special', '', '', '3']),
            'units': 'Bohr', 
            'fix_com': False,
            'fix_orientation': True, 
            'fragment_separators': [2, 3],
            }


def test_psi4_molstr2a():
    subject = '\n--\n'.join(subject2)
    fullans = copy.deepcopy(fullans2)
    fullans.update({'molecular_charge': 0.,
                    'molecular_multiplicity': 2,
                    'fragment_charges': [0., 0., 0.],
                    'fragment_multiplicities': [1, 1, 2]})

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans2, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr2b():
    subject = copy.deepcopy(subject2)    
    subject.insert(0, '1 3')
    subject = '\n--\n'.join(subject)
    ans = copy.deepcopy(ans2)
    ans.update({'molecular_charge': 1.,
                'molecular_multiplicity': 3})
    fullans = copy.deepcopy(fullans2)
    fullans.update({'molecular_charge': 1.,
                    'molecular_multiplicity': 3,
                    'fragment_charges': [1., 0., 0.],
                    'fragment_multiplicities': [2, 1, 2]})

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans, final, 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_molstr2c():
    """double overall chg/mult spec"""
    subject = copy.deepcopy(subject2)    
    subject.insert(0, '1 3\n1 3')
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

def test_psi4_molstr2d():
    """trailing comma"""
    subject = copy.deepcopy(subject2)    
    subject.insert(0, 'H 10,10,10,')
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

def test_psi4_molstr2e():
    """empty fragment"""
    subject = copy.deepcopy(subject2)    
    subject.insert(2, '\n')
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

def test_psi4_molstr2f():
    """double frag chgmult"""
    subject = copy.deepcopy(subject2)    
    subject[1] += '\n 1 2\n 5 6'
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

def test_psi4_molstr2g():
    """illegal chars in nucleus"""
    subject = copy.deepcopy(subject2)    
    subject[1] = """@Ne_{CN}_O 2 4 6"""
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

def test_psi4_molstr3():
    """psi4/psi4#731"""
    subject = """0 1
Mg 0 0"""

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)

subject4 = """pubchem:benzene"""

ans4 = {'geom': [
  -1.213100 , -0.688400 ,  0.000000 ,
  -1.202800 ,  0.706400 ,  0.000100 ,
  -0.010300 , -1.394800 ,  0.000000 ,
   0.010400 ,  1.394800 , -0.000100 ,
   1.202800 , -0.706300 ,  0.000000 ,
   1.213100 ,  0.688400 ,  0.000000 ,
  -2.157700 , -1.224400 ,  0.000000 ,
  -2.139300 ,  1.256400 ,  0.000100 ,
  -0.018400 , -2.480900 , -0.000100 ,
   0.018400 ,  2.480800 ,  0.000000 ,
   2.139400 , -1.256300 ,  0.000100 ,
   2.157700 ,  1.224500 ,  0.000000 ],
        'elbl': ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
        'units': 'Angstrom', 
        'fragment_separators': [],
        'fragment_charges': [None],
        'fragment_multiplicities': [None],
        }

fullans4 = {'geom': np.array([
  -1.213100 , -0.688400 ,  0.000000 ,
  -1.202800 ,  0.706400 ,  0.000100 ,
  -0.010300 , -1.394800 ,  0.000000 ,
   0.010400 ,  1.394800 , -0.000100 ,
   1.202800 , -0.706300 ,  0.000000 ,
   1.213100 ,  0.688400 ,  0.000000 ,
  -2.157700 , -1.224400 ,  0.000000 ,
  -2.139300 ,  1.256400 ,  0.000100 ,
  -0.018400 , -2.480900 , -0.000100 ,
   0.018400 ,  2.480800 ,  0.000000 ,
   2.139400 , -1.256300 ,  0.000100 ,
   2.157700 ,  1.224500 ,  0.000000 ]),
            'elbl': np.array(['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']),
            'elea': np.array([12, 12, 12, 12, 12, 12, 1, 1, 1, 1, 1, 1]),
            'elez': np.array([6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]),
            'elem': np.array(['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']),
            'mass': np.array([ 12., 12., 12., 12., 12., 12., 1.00782503,  1.00782503, 1.00782503, 1.00782503, 1.00782503, 1.00782503]),
            'real': np.array([ True, True, True, True, True, True, True, True, True, True, True, True]),
            'elbl': np.array(['', '', '', '', '', '', '', '', '', '', '', '']),
            'units': 'Angstrom', 
            'fix_com': False,
            'fix_orientation': False, 
            'fragment_separators': [],
            'molecular_charge': 0.,
            'molecular_multiplicity': 1,
            'fragment_charges': [0.],
            'fragment_multiplicities': [1],
            }

def test_pubchem_molstr4a():
    subject = subject4

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans4, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans4, final, 4, sys._getframe().f_code.co_name + ': full')

def test_pubchem_molstr4b():
    """user units potentially contradicting pubchem units"""
    subject = subject4 + '\nunits au'

    with pytest.raises(qcdb.ValidationError):
        final, intermed = qcdb.molparse.from_string(subject, return_processed=True)


def test_pubchem_molstr4a():
    subject = """
pubchem  : 241
"""

    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans4, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans4, final, 4, sys._getframe().f_code.co_name + ': full')

