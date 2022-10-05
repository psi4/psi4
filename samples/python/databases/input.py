#! PsiAPI energy example

import psi4

psi4.set_output_file("output.dat", False)
psi4.set_options({"SCF_TYPE": "DF", "BASIS": "cc-pVDZ", "REFERENCE": "ROHF"})

# just check if it runs :)
psi4.wrapper_database.database('scf', 'O24by5mb', cp=1, subset=["4-0.9", "12-0.9", "20-0.9", "23-0.9"])
#psi4.wrapper_database.database('scf', 'O24by5mb', cp=0, subset=["4-0.9", "12-0.9", "20-0.9", "23-0.9"])
psi4.wrapper_database.database('scf', 'O24by5', cp=1, subset=["4-0.9", "12-0.9", "20-0.9", "23-0.9"])
psi4.wrapper_database.database('scf', 'O24by5', cp=0, subset=["4-0.9", "12-0.9", "20-0.9", "23-0.9"])
