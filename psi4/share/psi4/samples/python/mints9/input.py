from __future__ import print_function
import psi4
from psi4.driver import qcdb
#! A test of the basis specification.  Various basis sets are specified outright and in blocks, both
#! orbital and auxiliary. Constructs libmints BasisSet objects through the constructor that calls
#! qcdb.BasisSet infrastructure. Checks that the resulting bases are of the right size and checks
#! that symmetry of the Molecule observes the basis assignment to atoms.

#           cc-pvdz                 aug-cc-pvdz
# BASIS     H  5/ 5   C  14/15      H +4/ 4   C  +9/10
# RIFIT     H 14/15   C  56/66      H +9/10   C +16/20
# JKFIT     H 23/25   C  70/81      H +9/10   C +16/20

psi4.set_output_file("output.dat", False)

mymol = psi4.geometry("""
C    0.0  0.0 0.0
O    1.4  0.0 0.0
H_r -0.5 -0.7 0.0
H_l -0.5  0.7 0.0
""")

psi4.set_options({'basis': 'cc-pvdz'})

print('[1]    <<<  uniform cc-pVDZ  >>>')
wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol.print_out()


print('[2]        <<<  RIFIT (default)  >>>')
wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))
mymol.print_out()

print('[3]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C  >>>')
psi4.basis_helper("""
    assign cc-pvdz
    assign c aug-cc-pvdz
""", name='dz_PLUS')
wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol.print_out()


print('[4]        <<<  RIFIT (default)  >>>')
wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))
mymol.print_out()
wert.print_out()
mymol.print_out()


print('[5]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H_R  >>>')
psi4.basis_helper("""
    assign cc-pvdz
    assign c aug-cc-pvdz
    assign h_r aug-cc-pvdz
""",
name='dz_PLUSplus',
key='BASis')
wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol.print_out()


print('[6]    <<<  RIFIT (custom: force cc-pVDZ on H, default on C, O)  >>>')
psi4.basis_helper("""
    assign h cc-pvdz-ri
""",
name='dz_PLUSplusRI',
key='df_basis_mp2')
wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_MP2', psi4.core.get_global_option('DF_BASIS_MP2'), 'RIFIT', psi4.core.get_global_option('BASIS'))
mymol.print_out()
mymol.print_out()


print('[7]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
psi4.basis_helper("""
    assign cc-pvdz
    assign c aug-cc-pvdz
    assign h aug-cc-pvdz
""",
name = 'dz_PLUSplusplus')
wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol.print_out()


print('[8]        <<<  JKFIT (default)  >>>')
wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
mymol.print_out()

psi4.set_options({'basis': 'aug-cc-pvdz'})

print('[9]    <<<  aug-cc-pVDZ  >>>')
wert = psi4.core.BasisSet.build(mymol, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol.print_out()


print('[10]       <<<  JKFIT (default)  >>>')
wert = psi4.core.BasisSet.build(mymol, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
mymol.print_out()


mymol2 = psi4.geometry("""
0 2
C    0.0  0.0 0.0
O    1.4  0.0 0.0
H_r -0.5 -0.6 0.3
H_l -0.5  0.6 0.3
H_c -0.5  0.0 0.7
""")

psi4.set_options({'basis': 'dz_plusplusplus'})

print('[11]   <<<  cc-pVDZ w/ aug-cc-pVDZ on C, H  >>>')
wert = psi4.core.BasisSet.build(mymol2, 'BASIS', psi4.core.get_global_option('BASIS'))
mymol2.print_out()

hene = psi4.geometry("""
He
Ne 1 2.0
""")

psi4.basis_helper("""
    assign cc-pv5z
""", name='disguised5z')

psi4.core.set_global_option('DF_BASIS_MP2', '')  # clear df_basis_mp2 {...} to get autoaux below

print('[12]   <<<  cc-pV5Z on HeNe  >>>')
wert = psi4.core.BasisSet.build(hene, 'BASIS', psi4.core.get_global_option('BASIS'))
hene.print_out()

print('[13]   <<<  RI for cc-pV5Z on HeNe  >>>')
wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))
hene.print_out()

print('[14]   <<<  impossible JK for cc-pV5Z on HeNe  >>>')
error_tripped = 0
try:
    wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))
except qcdb.BasisSetNotFound:
    error_tripped = 1

psi4.basis_helper(key='df_basis_scf', name='uggh', block="""
    assign he DEF2-UNIVERSAL-JKFIT
""")
hene.print_out()

print('[15]   <<<  forced JK for cc-pV5Z on HeNe  >>>')
wert = psi4.core.BasisSet.build(hene, 'DF_BASIS_SCF', '', 'JKFIT', psi4.core.get_global_option('BASIS'))

