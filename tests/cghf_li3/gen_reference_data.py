#!/usr/bin/env python
"""
Computes the reference data for Li3 test case.

Writes pyscf_reference_data.npz (gitignored, PySCF-AO-ordered). Follow up with
make_psi4_guess.py, which converts it to Psi4's AO ordering and writes the
plain-text occupied-orbital guess (li3_cocc_guess.dat) that input.dat uses.

Can be removed once stability analysis is implemented.
"""

from pyscf import gto, scf
import numpy as np
np.set_printoptions(suppress=True, precision=6, linewidth=90)

mol = gto.M(
    basis="sto-3g",
    atom="Li 0 0 0; Li 1 0 0; Li .5 0 .866",
    verbose=4,
    spin=3,
    symmetry=False
)
mf = scf.GHF(mol).set(conv_tol=1e-8)
sad_dm_rhf = scf.hf.SCF.init_guess_by_atom(mf, mol)
guess_dm = scf.ghf._from_rhf_init_dm(sad_dm_rhf, breaksym=True)
# guess_dm = scf.ghf._from_rhf_init_dm(sad_dm_rhf, breaksym=True).astype(complex)
# noise = np.random.random(guess_dm.shape).astype(complex)
# noise *= 2j * np.finfo(float).eps
# guess_dm += noise
e_tot = mf.kernel(guess_dm)

stable = False
iter = 0
while not stable:
    iter+=1
    if iter > 5:
        raise RuntimeError("Stability did not converge.")
    mo_coeff, stable = mf.stability(return_status=True)
    if stable:
        break
    dm = mf.make_rdm1(mo_coeff=mo_coeff)
    e_tot = mf.kernel(dm)
assert mf.converged

mo_occ = mf.mo_coeff[:,:mol.nelectron]

mo_occ_a = mo_occ[:mol.nao]
mo_occ_b = mo_occ[mol.nao:]

header = f" idx: {"alpha":>6}  {"beta":>6}"
print("\n"+header)
print('-'*len(header))
for i in range(mol.nelectron):
    occa = mo_occ[:mol.nao,i]
    occb = mo_occ[mol.nao:,i]
    a = np.dot(occa, occa.conj()).real
    b = np.dot(occb, occb.conj()).real
    print(f'{i+1:4g}  {a:6.4f}  {b:6.4f}')

dm = mf.make_rdm1()
d_aa = dm[:mol.nao,:mol.nao]
d_bb = dm[mol.nao:,mol.nao:]
d_ab = dm[:mol.nao,mol.nao:]
d_ba = dm[mol.nao:,:mol.nao]

mx = .5*(d_ba+d_ab)
my = -.5j*(d_ba-d_ab)
mz = .5*(d_aa-d_bb)
p = .5*(d_aa+d_bb)

print()
print(f"max P : {np.max(np.abs(p)):9.6f}")
print(f"max Mx: {np.max(np.abs(mx)):9.6f}")
print(f"max My: {np.max(np.abs(my)):9.6f}")
print(f"max Mz: {np.max(np.abs(mz)):9.6f}")

np.savez("pyscf_reference_data.npz", C_occ=mo_occ, ovlp=mf.get_ovlp(), hcore=mf.get_hcore())

## Without stability:
"""
converged SCF energy = -21.2432554046251  <S^2> = 0.7579178  2S+1 = 2.0079022
"""

## With stability:
"""
converged SCF energy = -21.2568481970714  <S^2> = 2.3917645  2S+1 = 3.2507012
"""

# Rerunning gives energies accurate to 1e-6, and multiplicities up to 1e-4

