import psi4
from psi4.driver.constants.physconst import hartree2ev

psi4.set_output_file("output.dat", False)

benz = psi4.geometry("""
    pubchem:benzene
""")

psi4.set_options({"REFERENCE"                : "RHF",
                  "MAX_ENERGY_G_CONVERGENCE" : 8,
                  "BASIS"                    : "STO-3G",
                  "DF_BASIS_SCF"             : "CC-PVDZ-RI"})

psi4.optimize('scf')

psi4.set_options({"REFERENCE"    : "RHF",
                  "BASIS"        : "CC-PVDZ",
                  "DF_BASIS_SCF" : "CC-PVDZ-JKFIT"})

e_sing_rhf = psi4.energy('scf')

benz.set_multiplicity(3)

psi4.set_options({"REFERENCE" : "ROHF"})
e_trip_rohf = psi4.energy('scf')
psi4.set_options({"REFERENCE" : "UHF"})
e_trip_uhf  = psi4.energy('scf')

vertical_uhf  = hartree2ev * (e_trip_uhf  - e_sing_rhf)
vertical_rohf = hartree2ev * (e_trip_rohf - e_sing_rhf)
psi4.core.print_out("\nSinglet-Triplet gap (vertical, UHF)   = %8.2f eV\n" % vertical_uhf)
psi4.core.print_out("\nSinglet-Triplet gap (vertical, ROHF)  = %8.2f eV\n" % vertical_rohf)

enuc  =  204.531600152395043                                                              #TEST
erhf  = -230.72190557842444                                                               #TEST
psi4.compare_values(enuc, benz.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy") #TEST 
psi4.compare_values(erhf,  e_sing_rhf,  6, "Singlet benzene RHF energy")                  #TEST 
