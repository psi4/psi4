import psi4

psi4.core.set_output_file('output.dat', False)
psi4.set_memory("2 GB")

# Reference ATM correction in kcal/mol from DFTD3
# >>> dftd3 geom.xyz -func b-p -abc -zero
ref_atm = -0.000110

mol = psi4.geometry("""
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
""")

mol.run_dftd3(f'pbe-atmgr')
psivars = psi4.get_variables()
test = psi4.constants.hartree2kcalmol * psivars['AXILROD-TELLER-MUTO 3-BODY DISPERSION ENERGY']
psi4.compare_values(ref_atm, test, 7, "Grimme ATM Correction")    # TEST
