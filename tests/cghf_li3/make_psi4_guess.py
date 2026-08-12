#!/usr/bin/env python
"""
Psi4 portion of test generation.

Converts the PySCF-converged (with stability analysis) Li3 GHF occupied
orbitals -- produced by gen_reference_data.py, in PySCF's AO ordering -- into
Psi4's AO ordering, and writes them out as a plain-text guess file
(``li3_cocc_guess.dat``) that ``input.dat`` loads directly with
``numpy.loadtxt`` and feeds to ``energy(..., guess_C=...)``. No checkpoint
file and no manually-constructed Wavefunction are needed at test time; the
Wavefunction built below exists only to sanity-check the AO permutation.

Without stability:
>>> converged SCF energy = -21.2432554046251  <S^2> = 0.7579178  2S+1 = 2.0079022

With stability:
>>> converged SCF energy = -21.2568481970714  <S^2> = 2.3917645  2S+1 = 3.2507012

Rerunning gives energies accurate to 1e-6, and multiplicities up to 1e-4
"""

import psi4
import numpy as np
from scipy.linalg import block_diag
from psi4.driver.procrouting.dft import build_superfunctional

####### Load reference data

with np.load('pyscf_reference_data.npz') as data:
    mo_occ = data['C_occ'].astype(complex)
    ovlp = data['ovlp']
    hcore = data['hcore']

####### Convert PySCF AO basis to Psi4 AO basis ########

permutation_matrix = []
S_permute = np.asarray([1])
P_permute = np.asarray([
    [0,1,0],
    [0,0,1],
    [1,0,0]
])
for _ in range(6):
    permutation_matrix.append(S_permute)
    permutation_matrix.append(S_permute)
    permutation_matrix.append(P_permute)
perm = block_diag(*permutation_matrix)

S_pyscf = perm.T @ ovlp @ perm
H_pyscf = perm.T @ hcore @ perm

####### Sanity-check the permutation against Psi4's own AO integrals #######

mol = psi4.geometry("""
0 4
  Li 0 0 0
  Li 1 0 0
  Li .5 0 .866
symmetry C1
""")
psi4.set_options({
    "basis": "STO-3G",
    "reference": "cghf",
    "scf_type": "direct",
    "df_scf_guess": False,
})
basis = psi4.core.BasisSet.build(mol, "ORBITAL", "STO-3G")
ref_wfn = psi4.core.ComplexWavefunction(mol, basis)
superfunc, _ = build_superfunctional("hf", False)
wfn = psi4.core.CGHF(ref_wfn, superfunc)
wfn.form_H()
wfn.form_Shalf()

S = wfn.S().to_array().real
H = wfn.H().to_array().real
assert np.allclose(S, S_pyscf, atol=1e-5)
assert np.allclose(H, H_pyscf, atol=1e-5)

####### Permute the orbitals and write the plain-text guess file #######

mo_occ = perm.T @ mo_occ
np.savetxt("li3_cocc_guess.dat", mo_occ)
