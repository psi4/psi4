import sys
import copy
import pprint

import pytest
import numpy as np

from utils import *
from addons import *

import qcelemental as qcel

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
        'fragment_files': [],
        'geom_hints': [],
        'hint_types': [],
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


def test_psi4_qm_1a():
    subject = subject1

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans1, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1a, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_1b():
    subject = '\n' + '\t' + subject1 + '\n\n'

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans1, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1a, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_1c():
    subject = '1 1\n  -- \n' + subject1
    ans = copy.deepcopy(ans1)
    ans.update({'molecular_charge': 1.,
                'molecular_multiplicity': 1})

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1c, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_1d():
    subject = subject1 + '\n1 1'
    ans = copy.deepcopy(ans1)
    ans.update({'fragment_charges': [1.],
                'fragment_multiplicities': [1]})

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans1c, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_1e():
    """duplicate com"""
    subject = subject1 + '\n  nocom'

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)


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
        'fragment_files': [],
        'geom_hints': [],
        'hint_types': [],
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


def test_psi4_qm_2a():
    subject = '\n--\n'.join(subject2)
    fullans = copy.deepcopy(fullans2)
    fullans.update({'molecular_charge': 0.,
                    'molecular_multiplicity': 2,
                    'fragment_charges': [0., 0., 0.],
                    'fragment_multiplicities': [1, 1, 2]})

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans2, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_2b():
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

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_qm_2c():
    """double overall chg/mult spec"""
    subject = copy.deepcopy(subject2)
    subject.insert(0, '1 3\n1 3')
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

def test_psi4_qm_2d():
    """trailing comma"""
    subject = copy.deepcopy(subject2)
    subject.insert(0, 'H 10,10,10,')
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

#def test_psi4_qm_2e():
#    """empty fragment"""
#    subject = copy.deepcopy(subject2)
#    subject.insert(2, '\n')
#    subject = '\n--\n'.join(subject)
#
#    with pytest.raises(qcel.MoleculeFormatError):
#        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

def test_psi4_qm_2f():
    """double frag chgmult"""
    subject = copy.deepcopy(subject2)
    subject[1] += '\n 1 2\n 5 6'
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

def test_psi4_qm_2g():
    """illegal chars in nucleus"""
    subject = copy.deepcopy(subject2)
    subject[1] = """@Ne_{CN}_O 2 4 6"""
    subject = '\n--\n'.join(subject)

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

def test_psi4_qm_3():
    """psi4/psi4#731"""
    subject = """0 1
Mg 0 0"""

    with pytest.raises(qcel.ValidationError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)

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
        'fragment_files': [],
        'geom_hints': [],
        'hint_types': [],
'molecular_charge': 0.0,
'name': 'IUPAC benzene',

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
'name': 'IUPAC benzene',
            }

def test_psi4_pubchem_4a():
    subject = subject4

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_molrecs(ans4, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans4, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


def test_psi4_pubchem_4b():
    """user units potentially contradicting pubchem units"""
    subject = subject4 + '\nunits au'

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)


def test_psi4_pubchem_4c():
    subject = """
pubchem  : 241
"""

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_molrecs(ans4, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans4, final['qm'], 4, sys._getframe().f_code.co_name + ': full')


subject5 = """
efp C6H6 -0.30448173 -2.24210052 -0.29383131 -0.642499 7.817407 -0.568147  # second to last equiv to 1.534222
--
efp C6H6 -0.60075437  1.36443336  0.78647823  3.137879 1.557344 -2.568550
"""

ans5 = {
    'fragment_files': ['C6H6', 'C6H6'],
    'hint_types': ['xyzabc', 'xyzabc'],
    'geom_hints': [[-0.30448173, -2.24210052, -0.29383131, -0.642499, 7.817407, -0.568147],
                  [-0.60075437,  1.36443336,  0.78647823,  3.137879, 1.557344, -2.568550]],
    'geom': [],
    'elbl': [],
    'fragment_charges': [None],
    'fragment_multiplicities': [None],
    'fragment_separators': [],
       }

fullans5b = {'efp': {}}
fullans5b['efp']['hint_types'] = ans5['hint_types']
fullans5b['efp']['geom_hints'] = ans5['geom_hints']
fullans5b['efp']['units'] = 'Bohr'
fullans5b['efp']['fix_com'] =  True
fullans5b['efp']['fix_orientation'] = True
fullans5b['efp']['fix_symmetry'] = 'c1'
fullans5b['efp']['fragment_files'] = ['c6h6', 'c6h6']


def test_psi4_efp_5a():
    subject = subject5

    hintsans = [[(val / qcdb.constants.bohr2angstroms if i < 3 else val) for i, val in enumerate(ans5['geom_hints'][0])],
                [(val / qcdb.constants.bohr2angstroms if i < 3 else val) for i, val in enumerate(ans5['geom_hints'][1])]]
    hintsans[0][4] = 1.534222
    fullans = copy.deepcopy(fullans5b)
    fullans['efp']['units'] = 'Angstrom'

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans5, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': final efp')

    hintsstd = qcel.util.standardize_efp_angles_units('Angstrom', final['efp']['geom_hints'])
    final['efp']['geom_hints'] = hintsstd
    fullans['efp']['geom_hints'] = hintsans
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': final efp standardized')

def test_psi4_efp_5b():
    subject = subject5 + '\nunits bohr'

    ans = copy.deepcopy(ans5)
    ans['units'] = 'Bohr'

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans5b['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': final efp')


def test_psi4_efp_5c():
    """fix_orientation not mol kw"""
    subject = subject5 + '\nno_com\nfix_orientation\nsymmetry c1'

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)


def test_psi4_efp_5d():
    subject = subject5 + '\nno_com\nno_reorient\nsymmetry c1\nunits a.u.'

    ans = copy.deepcopy(ans5)
    ans['units'] = 'Bohr'
    ans['fix_com'] = True
    ans['fix_orientation'] =  True
    ans['fix_symmetry'] = 'c1'

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans5b['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': final')


def test_psi4_efp_5e():
    """symmetry w/efp"""
    subject = subject5 + 'symmetry cs\nunits a.u.'

    with pytest.raises(qcel.ValidationError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)


subject6 = """
    0 1
    O1    0         0     0.118720
    h2   -0.753299, 0.0, -0.474880

    H3    0.753299, 0.0, -0.474880
    
    --
    efp h2O -2.12417561  1.22597097 -0.95332054 -2.902133 1.734999 -1.953647
 --
efp ammoniA
     0.98792    1.87681    2.85174
units au
     1.68798    1.18856    3.09517

     1.45873    2.55904    2.27226

"""

ans6 = {'units': 'Bohr',
        'geom': [0., 0., 0.118720, -0.753299, 0.0, -0.474880, 0.753299, 0.0, -0.474880],
        'elbl': ['O1', 'h2', 'H3'],
        'fragment_charges': [0.],
        'fragment_multiplicities': [1],
        'fragment_separators': [],
        'fragment_files': ['h2O', 'ammoniA'],
        'geom_hints': [[-2.12417561,  1.22597097, -0.95332054, -2.902133, 1.734999, -1.953647],
                           [0.98792,    1.87681,    2.85174, 1.68798 ,   1.18856  ,  3.09517, 1.45873  ,  2.55904  ,  2.27226]],
        'hint_types': ['xyzabc', 'points'],
        }

fullans6 = {'qm': {'geom': np.array([0., 0., 0.118720, -0.753299, 0.0, -0.474880, 0.753299, 0.0, -0.474880]),
                   'elea': np.array([16, 1, 1]),
                   'elez': np.array([8, 1, 1]),
                   'elem': np.array(['O', 'H', 'H']),
                   'mass': np.array([ 15.99491462,   1.00782503,   1.00782503]),
                   'real': np.array([True, True, True]),
                   'elbl': np.array(['1', '2', '3']),
                   'units': 'Bohr',
                   'fix_com': True,
                   'fix_orientation': True,
                   'fix_symmetry': 'c1',
                   'fragment_charges': [0.],
                   'fragment_multiplicities': [1],
                   'fragment_separators': [],
                   'molecular_charge': 0.,
                   'molecular_multiplicity': 1},
           'efp': {'fragment_files': ['h2o', 'ammonia'],
                   'geom_hints': [[-2.12417561,  1.22597097, -0.95332054, -2.902133, 1.734999, -1.953647],
                                  [0.98792,    1.87681,    2.85174, 1.68798 ,   1.18856  ,  3.09517, 1.45873  ,  2.55904  ,  2.27226]],
                   'hint_types': ['xyzabc', 'points'],
                   'units': 'Bohr',
                   'fix_com': True,
                   'fix_orientation': True,
                   'fix_symmetry': 'c1',
        }}


def test_psi4_qmefp_6a():
    subject = subject6

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans6, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans6['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')
    assert compare_molrecs(fullans6['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')


def test_psi4_qmefp_6b():
    subject = subject6.replace('au', 'ang')

    ans = copy.deepcopy(ans6)
    ans['units'] = 'Angstrom'

    fullans = copy.deepcopy(fullans6)
    fullans['qm']['units'] = 'Angstrom'
    fullans['efp']['units'] = 'Angstrom'

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')


def test_psi4_qmefp_6c():
    """try to give chgmult to an efp"""

    subject = subject6.replace('    efp h2O', '0 1\n    efp h2O')

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True)


@using_pylibefp
def test_psi4_qmefp_6d():
    subject = subject6

    fullans = copy.deepcopy(fullans6)
    fullans['efp']['geom'] = np.array([-2.22978429,  1.19270015, -0.99721732, -1.85344873,  1.5734809 ,
        0.69660583, -0.71881655,  1.40649303, -1.90657336,  0.98792   ,
        1.87681   ,  2.85174   ,  2.31084386,  0.57620385,  3.31175679,
        1.87761143,  3.16604791,  1.75667803,  0.55253064,  2.78087794,
        4.47837555])
    fullans['efp']['elea'] = np.array([16, 1, 1, 14, 1, 1, 1])
    fullans['efp']['elez'] = np.array([8, 1, 1, 7, 1, 1, 1])
    fullans['efp']['elem'] = np.array(['O', 'H', 'H', 'N', 'H', 'H', 'H'])
    fullans['efp']['mass'] = np.array([15.99491462, 1.00782503, 1.00782503, 14.00307400478, 1.00782503, 1.00782503, 1.00782503])
    fullans['efp']['real'] = np.array([True, True, True, True, True, True, True])
    fullans['efp']['elbl'] = np.array(['_a01o1', '_a02h2', '_a03h3', '_a01n1', '_a02h2', '_a03h3', '_a04h4'])
    fullans['efp']['fragment_separators'] = [3]
    fullans['efp']['fragment_charges'] = [0., 0.]
    fullans['efp']['fragment_multiplicities'] = [1, 1]
    fullans['efp']['molecular_charge'] = 0.
    fullans['efp']['molecular_multiplicity'] = 1
    fullans['efp']['hint_types'] = ['xyzabc', 'xyzabc']
    fullans['efp']['geom_hints'][1] = [1.093116487139866, 1.9296501432128303, 2.9104336205167156, -1.1053108079381473, 2.0333070957565544, -1.488586877218809]

    final, intermed = qcel.molparse.from_string(subject, return_processed=True)

    import pylibefp
    efpobj = pylibefp.from_dict(final['efp'])
    efpfinal = efpobj.to_dict()
    efpfinal = qcel.molparse.from_arrays(speclabel=False, domain='efp', **efpfinal)

    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    assert compare_molrecs(fullans['efp'], efpfinal, 4, sys._getframe().f_code.co_name + ': full efp')


subject7 = """\
5
   stuffs 
6Li 0.0 0.0 0.0 
H_specIAL@2.014101  100 0 0
@Ne 2 4 6
h .0,1,2
Gh(he3) 0 1 3
"""

ans7 = {'geom': [ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.],
        'elbl': ['6Li', 'H_specIAL@2.014101', '@Ne', 'h', 'Gh(he3)'],
        'units': 'Angstrom',
        'geom_hints': [],  # shouldn't be needed
        }

fullans7 = {'geom': np.array([ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.]),
            'elea': np.array([6,  2, 20, 1, 4]),
            'elez': np.array([3, 1, 10, 1, 2]),
            'elem': np.array(['Li', 'H', 'Ne', 'H', 'He']),
            'mass': np.array([ 6.015122794, 2.014101, 19.99244017542, 1.00782503, 4.00260325415]),
            'real': np.array([ True,  True, False, True, False]),
            'elbl': np.array(['', '_special', '', '', '3']),
            'units': 'Angstrom',
            'fix_com': False,
            'fix_orientation': False,
            'fragment_separators': [],
            'fragment_charges': [0.],
            'fragment_multiplicities': [2],
            'molecular_charge': 0.,
            'molecular_multiplicity': 2,
            }


def test_xyzp_qm_7a():
    """XYZ doesn't fit into psi4 string"""
    subject = subject7

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True, dtype='psi4')


def test_xyzp_qm_7b():
    """XYZ doesn't fit into strict xyz string"""
    subject = subject7

    with pytest.raises(qcel.MoleculeFormatError):
        final, intermed = qcel.molparse.from_string(subject, return_processed=True, dtype='xyz')


def test_xyzp_qm_7c():
    subject = subject7

    final, intermed = qcel.molparse.from_string(subject, return_processed=True, dtype='xyz+')
    assert compare_dicts(ans7, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans7, final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')


def test_xyzp_qm_7d():
    subject = subject7.replace('5', '5 au ')
    subject = subject.replace('stuff', '-1 3 slkdjfl2 32#$^& ')

    ans = copy.deepcopy(ans7)
    ans['units'] = 'Bohr'
    ans['molecular_charge'] = -1.
    ans['molecular_multiplicity'] = 3

    fullans = copy.deepcopy(fullans7)
    fullans['units'] = 'Bohr'
    fullans['fragment_charges'] = [-1.]
    fullans['fragment_multiplicities'] = [3]
    fullans['molecular_charge'] = -1.
    fullans['molecular_multiplicity'] = 3

    final, intermed = qcel.molparse.from_string(subject, return_processed=True, dtype='xyz+')
    assert compare_dicts(ans, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans, final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')

subject8 = """\
3
   stuffs 
Li 0.0 0.0 0.0 
1  100 0 0
Ne 2 4 6
h .0,1,2
 2 0 1 3
"""

ans8 = {'geom': [ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.],
        'elbl': ['Li', '1', 'Ne', 'h', '2'],
        'units': 'Angstrom',
        'geom_hints': [],  # shouldn't be needed
        }

fullans8 = {'geom': np.array([ 0.,  0.,  0.,  100.,  0.,  0., 2., 4., 6., 0., 1., 2., 0., 1., 3.]),
            'elea': np.array([7,  1, 20, 1, 4]),
            'elez': np.array([3, 1, 10, 1, 2]),
            'elem': np.array(['Li', 'H', 'Ne', 'H', 'He']),
            'mass': np.array([ 7.016004548, 1.00782503, 19.99244017542, 1.00782503, 4.00260325415]),
            'real': np.array([ True,  True, True, True, True]),
            'elbl': np.array(['', '', '', '', '']),
            'units': 'Angstrom',
            'fix_com': False,
            'fix_orientation': False,
            'fragment_separators': [],
            'fragment_charges': [0.],
            'fragment_multiplicities': [2],
            'molecular_charge': 0.,
            'molecular_multiplicity': 2,
            }


def test_xyzp_qm_8a():
    subject = subject8

    final, intermed = qcel.molparse.from_string(subject, return_processed=True, dtype='xyz+')
    assert compare_dicts(ans8, intermed, 4, sys._getframe().f_code.co_name + ': intermediate')
    assert compare_molrecs(fullans8, final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')



fullans10qm = {'geom': np.array([0., 0., 0.]),
               'elea': np.array([12]),
               'elez': np.array([6]),
               'elem': np.array(['C']),
               'mass': np.array([12.]),
               'real': np.array([True]),
               'elbl': np.array(['']),
               'units': 'Angstrom',
               'fix_com': False,
               'fix_orientation': False,
               'fragment_separators': [],
               'fragment_charges': [0.],
               'fragment_multiplicities': [1],
               'molecular_charge': 0.,
               'molecular_multiplicity': 1}
fullans10efp = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0., 0., 0., 0., 0., 0.]],
               'units': 'Angstrom',
               'fix_com': True,
               'fix_orientation': True,
               'fix_symmetry': 'c1'}
blankqm =     {'geom': np.array([]),
               'elea': np.array([]),
               'elez': np.array([]),
               'elem': np.array([]),
               'mass': np.array([]),
               'real': np.array([]),
               'elbl': np.array([]),
               'units': 'Angstrom',
               'fix_com': False,
               'fix_orientation': False,
               'fragment_separators': [],
               'fragment_charges': [0.],
               'fragment_multiplicities': [1],
               'molecular_charge': 0.,
               'molecular_multiplicity': 1}
blankefp =    {'fragment_files': [],
               'hint_types': [],
               'geom_hints': [],
               'units': 'Angstrom',
               'fix_com': True,
               'fix_orientation': True,
               'fix_symmetry': 'c1'}


def test_arrays_10a():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': True}

    fullans = {'qm': copy.deepcopy(fullans10qm),
               'efp': copy.deepcopy(fullans10efp)}
    fullans['qm']['fix_com'] = True
    fullans['qm']['fix_orientation'] = True
    fullans['qm']['fix_symmetry'] = 'c1'

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_arrays_10b():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': False,
               'enable_efp': True}

    fullans = {'efp': fullans10efp}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_arrays_10c():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': False}

    fullans = {'qm': fullans10qm}

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    with pytest.raises(KeyError):
        final['efp']


def test_arrays_10d():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_efp': 'none'}

    fullans = {'qm': fullans10qm}

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    with pytest.raises(KeyError):
        final['efp']


def test_arrays_10e():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_efp': 'minimal'}

    fullans = {'qm': fullans10qm,
               'efp': blankefp}

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_arrays_10f():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_efp': 'error'}

    with pytest.raises(qcel.ValidationError):
       qcel.molparse.from_input_arrays(**subject)


def test_arrays_10g():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': False,
               'enable_efp': True,
               'missing_enabled_return_efp': 'none'}

    fullans = {}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    with pytest.raises(KeyError):
        final['efp']

def test_arrays_10h():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': False,
               'enable_efp': True,
               'missing_enabled_return_efp': 'minimal'}

    fullans = {'efp': blankefp}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')

def test_arrays_10i():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': False,
               'enable_efp': True,
               'missing_enabled_return_efp': 'error'}

    with pytest.raises(qcel.ValidationError):
        qcel.molparse.from_input_arrays(**subject)



def test_arrays_10j():
    subject = {'geom': [0, 0, 0],
               'elem': ['C'],
               'enable_qm': True,
               'enable_efp': False}

    fullans = {'qm': fullans10qm}

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    with pytest.raises(KeyError):
        final['efp']


def test_arrays_10k():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_qm': 'none'}

    fullans = {'efp': fullans10efp}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_arrays_10l():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_qm': 'minimal'}

    fullans = {'qm': copy.deepcopy(blankqm),
               'efp': fullans10efp}
    fullans['qm']['fix_com'] = True
    fullans['qm']['fix_orientation'] = True
    fullans['qm']['fix_symmetry'] = 'c1'

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')

def test_arrays_10m():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': True,
               'missing_enabled_return_qm': 'error'}

    with pytest.raises(qcel.ValidationError):
        qcel.molparse.from_input_arrays(**subject)


def test_arrays_10n():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': False,
               'enable_efp': True}

    fullans = {'efp': fullans10efp}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_arrays_10o():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': False,
               'missing_enabled_return_qm': 'none'}

    fullans = {}

    final = qcel.molparse.from_input_arrays(**subject)
    with pytest.raises(KeyError):
        final['qm']
    with pytest.raises(KeyError):
        final['efp']

def test_arrays_10p():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': False,
               'missing_enabled_return_qm': 'minimal'}

    fullans = {'qm': blankqm}

    final = qcel.molparse.from_input_arrays(**subject)
    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    with pytest.raises(KeyError):
        final['efp']

def test_arrays_10q():
    subject = {'fragment_files': ['cl2'],
               'hint_types': ['xyzabc'],
               'geom_hints': [[0, 0, 0, 0, 0, 0]],
               'enable_qm': True,
               'enable_efp': False,
               'missing_enabled_return_qm': 'error'}

    with pytest.raises(qcel.ValidationError):
        qcel.molparse.from_input_arrays(**subject)

def test_strings_10r():
    subject = ''

    final = qcel.molparse.from_string(subject, enable_qm=True,
                                               enable_efp=True,
                                               missing_enabled_return_qm='none',
                                               missing_enabled_return_efp='none')

    print('final', final)
    with pytest.raises(KeyError):
        final['qm']
    with pytest.raises(KeyError):
        final['efp']

def test_strings_10s():
    subject = ''

    final = qcel.molparse.from_string(subject, enable_qm=True,
                                               enable_efp=True,
                                               missing_enabled_return_qm='minimal',
                                               missing_enabled_return_efp='minimal')

    fullans = {'qm': blankqm,
               'efp': blankefp}

    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
    assert compare_molrecs(fullans['efp'], final['efp'], 4, sys._getframe().f_code.co_name + ': full efp')


def test_strings_10t():
    subject = ''

    with pytest.raises(qcel.ValidationError):

        qcel.molparse.from_string(subject, enable_qm=True,
                                           enable_efp=True,
                                           missing_enabled_return_qm='error',
                                           missing_enabled_return_efp='error')




def assess_mol_11(mol, label):
    dmol = mol.to_dict()
    assert compare_molrecs(fullans1a, dmol, 4, label, relative_geoms='align')
    assert compare_integers(2, mol.natom(), label)

def test_qmol_11a():
    asdf = qcdb.Molecule(fullans1a)
    assess_mol_11(asdf, '[1] qcdb.Molecule(dict)')

def test_qmol_11b():
    asdf = qcdb.Molecule(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True)
    assess_mol_11(asdf, '[2] qcdb.Molecule(geom, elez)')

def test_qmol_11c():
    asdf = qcdb.Molecule("""nocom\n8 0 0 0\n1 1 0 0""", dtype='psi4')
    assess_mol_11(asdf, '[3] qcdb.Molecule(str, dtype="psi4")')

def test_qmol_11d():
    asdf = qcdb.Molecule("""nocom\n8 0 0 0\n1 1 0 0""", dtype='psi4+')
    assess_mol_11(asdf, '[4] qcdb.Molecule(str, dtype="psi4+")')

def test_qmol_11e():
    asdf = qcdb.Molecule("""2\n\nO 0 0 0 \n1 1 0 0 """, dtype='xyz', fix_com=True)
    assess_mol_11(asdf, '[5] qcdb.Molecule(str, dtype="xyz")')

def test_qmol_11f():
    asdf = qcdb.Molecule.from_dict(fullans1a)
    assess_mol_11(asdf, '[6] qcdb.Molecule.from_dict(dict)')

def test_qmol_11g():
    asdf = qcdb.Molecule.from_arrays(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True)
    assess_mol_11(asdf, '[7] qcdb.Molecule.from_arrays(geom, elez)')

def test_qmol_11h():
    asdf = qcdb.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
    assess_mol_11(asdf, '[8] qcdb.Molecule.from_string(str, dtype="psi4")')

def test_qmol_11i():
    asdf = qcdb.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
    assess_mol_11(asdf, '[9] qcdb.Molecule.from_string(str, dtype="psi4+")')

def test_qmol_11j():
    asdf = qcdb.Molecule.from_string("""2\n\nO 0 0 0 \n1 1 0 0 """, fix_com=True)
    assess_mol_11(asdf, '[10] qcdb.Molecule.from_string(str, dtype="xyz")')


@using_psi4_molrec
def test_pmol_11k():
    import psi4
    asdf = psi4.core.Molecule.from_dict(fullans1a)
    assess_mol_11(asdf, '[16] psi4.core.Molecule.from_dict(dict)')

@using_psi4_molrec
def test_pmol_11l():
    import psi4
    asdf = psi4.core.Molecule.from_arrays(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True)
    assess_mol_11(asdf, '[17] psi4.core.Molecule.from_arrays(geom, elez)')

@using_psi4_molrec
def test_pmol_11m():
    import psi4
    asdf = psi4.core.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
    assess_mol_11(asdf, '[18] psi4.core.Molecule.from_string(str, dtype="psi4")')

@using_psi4_molrec
def test_pmol_11n():
    import psi4
    asdf = psi4.core.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
    assess_mol_11(asdf, '[19] psi4.core.Molecule.from_string(str, dtype="psi4+")')

@using_psi4_molrec
def test_pmol_11o():
    import psi4
    asdf = psi4.core.Molecule.from_string("""2\n\nO 0 0 0 \n1 1 0 0 """, fix_com=True)
    assess_mol_11(asdf, '[20] psi4.core.Molecule.from_string(str, dtype="xyz")')

def test_qmol_11p():
    asdf = qcdb.Molecule.from_arrays(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True, units='AngSTRom')
    assess_mol_11(asdf, '[7] qcdb.Molecule.from_arrays(geom, elez)')

def test_qmol_12():
    asdf = qcdb.Molecule(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True)
    assess_mol_11(asdf, 'qcdb.Molecule(geom, elez)')

    import json
    smol = json.dumps(asdf.to_dict(np_out=False))
    dmol = json.loads(smol)

    asdf2 = qcdb.Molecule(dmol)
    assess_mol_11(asdf, 'qcdb.Molecule(jsondict)')
