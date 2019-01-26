import psi4

import numpy as np
from psi4.driver.p4util.solvers import davidson_solver
from psi4.driver.procrouting.response.scf_products import SCFProducts

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")

psi4.set_options({"SAVE_JK": True})
psi4.set_options({"e_convergence": 1.e-1, "d_convergence": 1.e-1})
# psi4.set_options({"reference": "uhf"})
e, wfn = psi4.energy("HF/6-31G", return_wfn=True)

nmo = wfn.nmopi().sum()
ndocc = wfn.doccpi().sum()
nvir = nmo - ndocc
nrot = ndocc * nvir

wfn.form_D()

record = [(False, wfn.compute_E())]
for x in range(5):

    prod = SCFProducts(wfn)

    def func(vector):
        return -1 * prod.H1_product(vector)

    def precon(resid, i, A_w):
        return resid

    nvecs = 5
    guess = np.ones((prod.narot, nvecs))
    evals, evecs = davidson_solver(func, precon, guess, no_eigs=nvecs, e_conv=1.e-4)

    # If x == 0 we take a "bad" step, require rotation next step to get back on track
    if (x == 0) or (evals[0] < 0):
        stab = True
        evecs = -evecs[:, 0]
        rot = psi4.core.Matrix.from_array(evecs.reshape(ndocc, nvir))
        wfn.rotate_orbitals(wfn.Ca(), rot)
    else:
        stab = False
        wfn.form_C()

    wfn.form_D()
    wfn.form_G()
    wfn.form_F()
    print(wfn.compute_E())
    record.append((stab, wfn.compute_E()))

for stab, energy in record:
    print("% 8s % 14.10f" % (stab, energy))
