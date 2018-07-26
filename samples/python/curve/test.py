#! PsiAPI scanning a potential energy curve

import psi4
import numpy as np

psi4.set_output_file("output.dat", False)

mol_string = """
He
He 1 **R**
"""

mols = [(d, mol_string.replace("**R**", "%5.3f" % d)) for d in np.linspace(2, 5, 4)] 
psi4.set_options({"SCF_TYPE": "DF",
                  "BASIS": "6-31G"})


out = {}
for d, mol in mols:
    geom = psi4.geometry(mol)
    out[d] = psi4.energy('SCF', molecule=geom)


bench_dict = {2.0: -5.70825982153,
              3.0: -5.71036140727,
              4.0: -5.71036306964,
              5.0: -5.71036203769}

for key in bench_dict.keys():
    val = out[key]
    bench = bench_dict[key]
    psi4.compare_values(val, bench, 6, 'SCF Energy @%5.3fR' % key)
