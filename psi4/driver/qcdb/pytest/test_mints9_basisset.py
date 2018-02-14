import pytest
from utils import *
from addons import *

import qcdb

#! A test of the basis specification.  Various basis sets are specified outright and in blocks, both
#! orbital and auxiliary. Constructs libmints BasisSet objects through the constructor that calls
#! qcdb.BasisSet infrastructure. Checks that the resulting bases are of the right size and checks
#! that symmetry of the Molecule observes the basis assignment to atoms.

#           cc-pvdz                 aug-cc-pvdz
# BASIS     H  5/ 5   C  14/15      H +4/ 4   C  +9/10
# RIFIT     H 14/15   C  56/66      H +9/10   C +16/20
# JKFIT     H 23/25   C  70/81      H +9/10   C +16/20


smol = """
C    0.0  0.0 0.0
O    1.4  0.0 0.0
H_r -0.5 -0.7 0.0
H_l -0.5  0.7 0.0
"""


BASIS = 'cc-pvdz'
verbose = 2

def test_1():
    """[1]    <<<  uniform cc-pVDZ  >>>"""

    qmol = qcdb.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'BASIS', BASIS, verbose=verbose, return_dict=True)

    compare_strings('CC-PVDZ', BASIS, 'name')
    compare_integers(38, wert.nbf(), 'nbf()')
    compare_integers(40, wert.nao(), 'nao()')  
    compare_strings('c2v', wert.molecule.schoenflies_symbol(), 'symm')  
    compare_strings('CC-PVDZ', dwert['name'], 'callby')  
    compare_strings('CC-PVDZ', dwert['blend'], 'blend')  

@using_psi4
def test_1b():
    """[1]    <<<  uniform cc-pVDZ  >>>"""
    import psi4

    qmol = qcdb.Molecule.from_string(smol)
    pmol = psi4.core.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'BASIS', BASIS, verbose=verbose, return_dict=True)

    pwert = psi4.core.BasisSet.construct_from_pydict(pmol, dwert, -1)
    compare_integers(38, pwert.nbf(), 'nbf()')
    compare_integers(40, pwert.nao(), 'nao()')  
    compare_strings('c2v', pwert.molecule().schoenflies_symbol(), 'symm')  
    compare_strings('CC-PVDZ', pwert.name(), 'callby')  
    compare_strings('CC-PVDZ', pwert.blend(), 'blend')  

def test_2():
    """[2]        <<<  RIFIT (default)  >>>"""

    qmol = qcdb.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'DF_BASIS_MP2', '', 'RIFIT', BASIS, verbose=verbose, return_dict=True)

    compare_integers(140, wert.nbf(), 'nbf()')
    compare_integers(162, wert.nao(), 'nao()')  
    compare_strings('c2v', wert.molecule.schoenflies_symbol(), 'symm')  
    compare_strings('(CC-PVDZ AUX)', dwert['name'], 'callby')  
    compare_strings('CC-PVDZ-RI', dwert['blend'], 'blend')  

@using_psi4
def test_2b():
    """[2]        <<<  RIFIT (default)  >>>"""
    import psi4

    qmol = qcdb.Molecule.from_string(smol)
    pmol = psi4.core.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'DF_BASIS_MP2', '', 'RIFIT', BASIS, verbose=verbose, return_dict=True)

    pwert = psi4.core.BasisSet.construct_from_pydict(pmol, dwert, -1)
    compare_integers(140, pwert.nbf(), 'nbf()')
    compare_integers(162, pwert.nao(), 'nao()')  
    compare_strings('c2v', pwert.molecule().schoenflies_symbol(), 'symm')  
    compare_strings('(CC-PVDZ AUX)', pwert.name(), 'callby')  
    compare_strings('CC-PVDZ-RI', pwert.blend(), 'blend')  

def test_3():
    """[3]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C  >>>"""

    def basisspec_psi4_yo__dz_plus(mol, role):
        basstrings = {}
        mol.set_basis_all_atoms("cc-pvdz", role=role)
        mol.set_basis_by_symbol("c", "aug-cc-pvdz", role=role)
        return basstrings
    qcdb.libmintsbasisset.basishorde['DZ_PLUS'] = basisspec_psi4_yo__dz_plus

    qmol = qcdb.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'BASIS', BASIS, verbose=verbose, return_dict=True)

    compare_integers(47, wert.nbf(), 'nbf()')
    compare_integers(50, wert.nao(), 'nao()')  
    compare_strings('c2v', wert.molecule.schoenflies_symbol(), 'symm')  
    compare_strings('DZ_PLUS', dwert['name'], 'callby')  
    compare_strings('AUG-CC-PVDZ + CC-PVDZ', dwert['blend'], 'blend')  

@using_psi4
def test_3b():
    """[3]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C  >>>"""
    import psi4

    psi4.basis_helper("""
    assign cc-pvdz
    assign c aug-cc-pvdz
""", name='dz_PLUS')

    qmol = qcdb.Molecule.from_string(smol)
    pmol = psi4.core.Molecule.from_string(smol)
    wert, dwert = qcdb.BasisSet.pyconstruct(qmol, 'BASIS', BASIS, verbose=verbose, return_dict=True)

    pwert = psi4.core.BasisSet.construct_from_pydict(pmol, dwert, -1)
    compare_integers(47, pwert.nbf(), 'nbf()')
    compare_integers(50, pwert.nao(), 'nao()')  
    compare_strings('c2v', pwert.molecule().schoenflies_symbol(), 'symm')  
    compare_strings('DZ_PLUS', pwert.name(), 'callby')  
    compare_strings('AUG-CC-PVDZ + CC-PVDZ', pwert.blend(), 'blend')  




#print('[4]        <<<  RIFIT (default)  >>>')
#wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))
#mymol.print_out()
#wert.print_out()
#psi4.compare_integers(156, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(182, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('(DZ_PLUS AUX)', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ-RI + CC-PVDZ-RI', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#print('[5]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H_R  >>>')
#psi4.basis_helper("""
#    assign cc-pvdz
#    assign c aug-cc-pvdz
#    assign h_r aug-cc-pvdz
#""",
#name='dz_PLUSplus',
#key='BASis')
#wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
#psi4.compare_strings('DZ_PLUSPLUS', psi4.core.get_global_option('BASIS'), 'name')  #TEST
#psi4.compare_integers(51, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(54, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('cs', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('DZ_PLUSPLUS', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#print('[6]    <<<  RIFIT (custom: force cc-pVDZ on H, default on C, O)  >>>')
#psi4.basis_helper("""
#    assign h cc-pvdz-ri
#""",
#name='dz_PLUSplusRI',
#key='df_basis_mp2')
#wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', psi4.core.get_global_option('DF_BASIS_MP2'), 'RIFIT', psi4.core.get_global_option('BASIS'))
#mymol.print_out()
#psi4.compare_integers(156, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(182, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('cs', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('DZ_PLUSPLUSRI', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ-RI + CC-PVDZ-RI', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#print('[7]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
#psi4.basis_helper("""
#    assign cc-pvdz
#    assign c aug-cc-pvdz
#    assign h aug-cc-pvdz
#""",
#name = 'dz_PLUSplusplus')
#wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(55, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(58, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('DZ_PLUSPLUSPLUS', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#print('[8]        <<<  JKFIT (default)  >>>')
#wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(220, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(252, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('(DZ_PLUSPLUSPLUS AUX)', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ-JKFIT + CC-PVDZ-JKFIT', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#psi4.set_options({'basis': 'aug-cc-pvdz'})
#
#print('[9]    <<<  aug-cc-pVDZ  >>>')
#wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(64, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(68, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#print('[10]       <<<  JKFIT (default)  >>>')
#wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(236, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(272, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('(AUG-CC-PVDZ AUX)', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ-JKFIT', wert.blend(), 'blend')  #TEST
#mymol.print_out()
#
#
#mymol2 = psi4.geometry("""
#0 2
#C    0.0  0.0 0.0
#O    1.4  0.0 0.0
#H_r -0.5 -0.6 0.3
#H_l -0.5  0.6 0.3
#H_c -0.5  0.0 0.7
#""")
#
#psi4.set_options({'basis': 'dz_plusplusplus'})
#
#print('[11]   <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
#wert = psi4.core.BasisSet.build(mymol2, 'BASIS', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(64, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(67, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('cs', mymol2.schoenflies_symbol(), 'symm')  #TEST
#psi4.compare_strings('DZ_PLUSPLUSPLUS', wert.name(), 'callby')  #TEST
#psi4.compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend')  #TEST
#mymol2.print_out()
#
#hene = psi4.geometry("""
#He
#Ne 1 2.0
#""")
#
#psi4.basis_helper("""
#    assign cc-pv5z
#""", name='disguised5z')
#
#psi4.core.set_global_option('DF_BASIS_MP2', '')  # clear df_basis_mp2 {...} to get autoaux below
#
#print('[12]   <<<  cc-pV5Z on HeNe  >>>')
#wert = psi4.core.BasisSet.build(hene, 'BASIS', psi4.core.get_global_option('BASIS'))
#hene.print_out()
#psi4.compare_integers(146, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(196, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('DISGUISED5Z', wert.name(), 'callby')  #TEST
#psi4.compare_strings('CC-PV5Z', wert.blend(), 'blend')  #TEST
#
#print('[13]   <<<  RI for cc-pV5Z on HeNe  >>>')
#wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))
#hene.print_out()
#psi4.compare_integers(284, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(413, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('(DISGUISED5Z AUX)', wert.name(), 'callby')  #TEST
#psi4.compare_strings('CC-PV5Z-RI', wert.blend(), 'blend')  #TEST
#
#print('[14]   <<<  impossible JK for cc-pV5Z on HeNe  >>>')
#error_tripped = 0
#try:
#    wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
#except qcdb.BasisSetNotFound:
#    error_tripped = 1
#psi4.compare_integers(1, error_tripped, 'squashed 4z aux for 5z orb')  #TEST
#
#psi4.basis_helper(key='df_basis_scf', name='uggh', block="""
#    assign he DEF2-QZVPP-JKFIT
#""")
#hene.print_out()
#
#print('[15]   <<<  forced JK for cc-pV5Z on HeNe  >>>')
#wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
#psi4.compare_integers(169, wert.nbf(), 'nbf()')  #TEST
#psi4.compare_integers(241, wert.nao(), 'nao()')  #TEST
#psi4.compare_strings('UGGH', wert.name(), 'callby')  #TEST
#psi4.compare_strings('CC-PV5Z-JKFIT + DEF2-QZVPP-JKFIT', wert.blend(), 'blend')  #TEST
#

#       print('[4]        <<<  RIFIT (default)  >>>')
#       wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', '', 'RIFIT', get_global_option('BASIS'))
#       compare_integers(156, wert.nbf(), 'nbf()') 
#       compare_integers(182, wert.nao(), 'nao()') 
#       compare_strings('c2v', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('(DZ_PLUS AUX)', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ-RI + CC-PVDZ-RI', wert.blend(), 'blend') 
#       mymol.print_out()
#       print('[5]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H_R  >>>')
#       def basisspec_psi4_yo__dz_plusplus(mol, role):
#           basstrings = {}
#           mol.set_basis_all_atoms("cc-pvdz", role=role)
#           mol.set_basis_by_symbol("c", "aug-cc-pvdz", role=role)
#           mol.set_basis_by_label("h_r", "aug-cc-pvdz", role=role)
#           return basstrings
#       qcdb.libmintsbasisset.basishorde['DZ_PLUSPLUS'] = basisspec_psi4_yo__dz_plusplus
#       core.set_global_option("BASIS", "dz_PLUSplus")
#       wert = psi4.core.BasisSet.build(mymol, 'BASIS', get_global_option('BASIS'))
#       compare_strings('DZ_PLUSPLUS', get_global_option('BASIS'), 'name') 
#       compare_integers(51, wert.nbf(), 'nbf()') 
#       compare_integers(54, wert.nao(), 'nao()') 
#       compare_strings('cs', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('DZ_PLUSPLUS', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend') 
#       mymol.print_out()
#       print('[6]    <<<  RIFIT (custom: force cc-pVDZ on H, default on C, O)  >>>')
#       def basisspec_psi4_yo__dz_plusplusri(mol, role):
#           basstrings = {}
#           mol.set_basis_by_symbol("h", "cc-pvdz-ri", role=role)
#           return basstrings
#       qcdb.libmintsbasisset.basishorde['DZ_PLUSPLUSRI'] = basisspec_psi4_yo__dz_plusplusri
#       core.set_global_option("DF_BASIS_MP2", "dz_PLUSplusRI")
#       wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', get_global_option('DF_BASIS_MP2'), 'RIFIT', get_global_option('BASIS'))
#       compare_integers(156, wert.nbf(), 'nbf()') 
#       compare_integers(182, wert.nao(), 'nao()') 
#       compare_strings('cs', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('DZ_PLUSPLUSRI', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ-RI + CC-PVDZ-RI', wert.blend(), 'blend') 
#       mymol.print_out()
#       print('[7]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
#       def basisspec_psi4_yo__dz_plusplusplus(mol, role):
#           basstrings = {}
#           mol.set_basis_all_atoms("cc-pvdz", role=role)
#           mol.set_basis_by_symbol("c", "aug-cc-pvdz", role=role)
#           mol.set_basis_by_symbol("h", "aug-cc-pvdz", role=role)
#           return basstrings
#       qcdb.libmintsbasisset.basishorde['DZ_PLUSPLUSPLUS'] = basisspec_psi4_yo__dz_plusplusplus
#       core.set_global_option("BASIS", "dz_PLUSplusplus")
#       wert = psi4.core.BasisSet.build(mymol, 'BASIS', get_global_option('BASIS'))
#       compare_integers(55, wert.nbf(), 'nbf()') 
#       compare_integers(58, wert.nao(), 'nao()') 
#       compare_strings('c2v', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('DZ_PLUSPLUSPLUS', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend') 
#       mymol.print_out()
#       print('[8]        <<<  JKFIT (default)  >>>')
#       wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', get_global_option('BASIS'))
#       compare_integers(220, wert.nbf(), 'nbf()') 
#       compare_integers(252, wert.nao(), 'nao()') 
#       compare_strings('c2v', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('(DZ_PLUSPLUSPLUS AUX)', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ-JKFIT + CC-PVDZ-JKFIT', wert.blend(), 'blend') 
#       mymol.print_out()
#       core.set_global_option("BASIS", "aug-cc-pvdz")
#       print('[9]    <<<  aug-cc-pVDZ  >>>')
#       wert = psi4.core.BasisSet.build(mymol, 'BASIS', get_global_option('BASIS'))
#       compare_integers(64, wert.nbf(), 'nbf()') 
#       compare_integers(68, wert.nao(), 'nao()') 
#       compare_strings('c2v', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('AUG-CC-PVDZ', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ', wert.blend(), 'blend') 
#       mymol.print_out()
#       print('[10]       <<<  JKFIT (default)  >>>')
#       wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', get_global_option('BASIS'))
#       compare_integers(236, wert.nbf(), 'nbf()') 
#       compare_integers(272, wert.nao(), 'nao()') 
#       compare_strings('c2v', mymol.schoenflies_symbol(), 'symm') 
#       compare_strings('(AUG-CC-PVDZ AUX)', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ-JKFIT', wert.blend(), 'blend') 
#       mymol.print_out()
#       mymol2 = geometry("""
#       C    0.0  0.0 0.0
#       O    1.4  0.0 0.0
#       H_r -0.5 -0.6 0.3
#       H_l -0.5  0.6 0.3
#       H_c -0.5  0.0 0.7
#       ""","mymol2")
#       core.IO.set_default_namespace("mymol2")
#       core.set_global_option("BASIS", "dz_plusplusplus")
#       print('[11]   <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
#       wert = psi4.core.BasisSet.build(mymol2, 'BASIS', get_global_option('BASIS'))
#       compare_integers(64, wert.nbf(), 'nbf()') 
#       compare_integers(67, wert.nao(), 'nao()') 
#       compare_strings('cs', mymol2.schoenflies_symbol(), 'symm') 
#       compare_strings('DZ_PLUSPLUSPLUS', wert.name(), 'callby') 
#       compare_strings('AUG-CC-PVDZ + CC-PVDZ', wert.blend(), 'blend') 
#       mymol2.print_out()
#       hene = geometry("""
#       He
#       Ne 1 2.0
#       ""","hene")
#       core.IO.set_default_namespace("hene")
#       def basisspec_psi4_yo__disguised5z(mol, role):
#           basstrings = {}
#           mol.set_basis_all_atoms("cc-pv5z", role=role)
#           return basstrings
#       qcdb.libmintsbasisset.basishorde['DISGUISED5Z'] = basisspec_psi4_yo__disguised5z
#       core.set_global_option("BASIS", "disguised5z")
#       set_global_option('DF_BASIS_MP2', '') 
#       print('[12]   <<<  cc-pV5Z on HeNe  >>>')
#       wert = psi4.core.BasisSet.build(hene, 'BASIS', get_global_option('BASIS'))
#       compare_integers(146, wert.nbf(), 'nbf()') 
#       compare_integers(196, wert.nao(), 'nao()') 
#       compare_strings('DISGUISED5Z', wert.name(), 'callby') 
#       compare_strings('CC-PV5Z', wert.blend(), 'blend') 
#       hene.print_out()
#       print('[13]   <<<  RI for cc-pV5Z on HeNe  >>>')
#       wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_MP2', '', 'RIFIT', get_global_option('BASIS'))
#       compare_integers(284, wert.nbf(), 'nbf()') 
#       compare_integers(413, wert.nao(), 'nao()') 
#       compare_strings('(DISGUISED5Z AUX)', wert.name(), 'callby') 
#       compare_strings('CC-PV5Z-RI', wert.blend(), 'blend') 
#       hene.print_out()
#       print('[14]   <<<  impossible JK for cc-pV5Z on HeNe  >>>')
#       error_tripped = 0
#       try:
#           wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', get_global_option('BASIS'))
#       except qcdb.BasisSetNotFound:
#           error_tripped = 1
#       compare_integers(1, error_tripped, 'squashed 4z aux for 5z orb') 
#       def basisspec_psi4_yo__uggh(mol, role):
#           basstrings = {}
#           mol.set_basis_by_symbol("he", "DEF2-QZVPP-JKFIT", role=role)
#           return basstrings
#       qcdb.libmintsbasisset.basishorde['UGGH'] = basisspec_psi4_yo__uggh
#       core.set_global_option("DF_BASIS_SCF", "uggh")
#       print('[15]   <<<  forced JK for cc-pV5Z on HeNe  >>>')
#       wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', get_global_option('BASIS'))
#       compare_integers(169, wert.nbf(), 'nbf()') 
#       compare_integers(241, wert.nao(), 'nao()') 
#       compare_strings('UGGH', wert.name(), 'callby') 
#       compare_strings('CC-PV5Z-JKFIT + DEF2-QZVPP-JKFIT', wert.blend(), 'blend') 
#       hene.print_out()

