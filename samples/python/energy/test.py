#! PsiAPI energy example

import psi4

psi4.set_output_file("output.dat", False)

geom = psi4.geometry("""
C  # testing escaping comments
""")

psi4.set_options({"SCF_TYPE": "DF",
                  "BASIS": "cc-pVDZ"})

scf_e, scf_wfn = psi4.energy('SCF', return_wfn=True)
psi4.compare_values(-37.5959861862713893, scf_e, 6, 'SCF DF Energy')


psi4.core.set_local_option("SCF", "SCF_TYPE", "PK")
psi4.compare_values(-37.5959861862713893, scf_e, 6, 'SCF PK Energy')

