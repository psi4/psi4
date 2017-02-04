# A simple Psi 4 input script to compute a SCF reference using second-order optimizations
# Requires scipy numpy 1.7.2+
#
# Created by: Daniel G. A. Smith
# Date: 2/27/15
# License: GPL v3.0
#

import time
import numpy as np
np.set_printoptions(precision=3, linewidth=200, suppress=True)
import psi4

# Grab a DIIS object, will be moved up later
from psi4.driver.procedures.mcscf.diis_helper import DIIS_helper

# Memory for Psi4 in GB
psi4.core.set_memory(int(2e9), False)
psi4.core.set_output_file('output.dat', False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 109
symmetry c1
""")

psi4.set_options({'scf_type':'pk',
                  'basis':'6-31G',
                  'e_convergence':1e1,
                  'd_convergence':1e1})

# Knobs
E_conv = 1.e-8
D_conv = 1.e-4
max_macro = 10
max_micro = 3
micro_conv = 1.e-3
micro_print = True

# Build objects

wfn = psi4.core.Wavefunction.build(mol)
ndocc = wfn.nalpha()
nvirt = wfn.nso() - ndocc

mints = psi4.core.MintsHelper(wfn.basisset())
eri = np.array(mints.ao_eri())

H = np.array(mints.ao_potential()) + np.array(mints.ao_kinetic())
S = np.array(mints.ao_overlap())

A = mints.ao_overlap()
A.power(-0.5, 1.e-14)
A = np.array(A)

# Build a DIIS object
diis = DIIS_helper(max_vec=6)


# Diagonalize routine
def build_orbitals(diag):
    Fp = A.dot(diag).dot(A)

    eigvals, Cp = np.linalg.eigh(Fp)
    C = A.dot(Cp)

    Cocc = C[:, :ndocc]
    D = np.dot(Cocc, Cocc.T)

    return C, Cocc, D

C, Cocc, D = build_orbitals(H)


print('\nStart SCF iterations:\n')
t = time.time()
E = 0.0
Eold = 0.0
iter_type = 'CORE'

def SCF_Hx(x, moF, Co, Cv):
    """
    Compute a hessian vector guess where x is a ov matrix of nonredundant operators.
    """
    F  = np.dot(moF[:ndocc, :ndocc], x)
    F -= np.dot(x, moF[ndocc:, ndocc:])

    #exit()

    # Build two electron part, M = -4 (4 G_{mnip} - g_{mpin} - g_{npim}) K_{ip}
    Dk = np.dot(Co, -x).dot(Cv.T)

    J = np.einsum('pqrs,rs->pq', eri, Dk) 
    K = np.einsum('prqs,rs->pq', eri, Dk)
    F += (Co.T).dot(4 * J - K.T - K).dot(Cv)
    F *= -4
    return F

for SCF_ITER in range(1, max_macro):

    # Build new fock matrix
    F = H + 2.0 * np.einsum('pqrs,rs->pq', eri, D) - np.einsum('prqs,rs->pq', eri, D)

    # DIIS error and update
    diis_e = F.dot(D).dot(S) - S.dot(D).dot(F)
    diis_e = (A).dot(diis_e).dot(A)
    diis.add(psi4.core.Matrix.from_array(F), psi4.core.Matrix.from_array(diis_e))

    # SCF energy and update
    scf_e = np.vdot(F + H, D) + mol.nuclear_repulsion_energy()
    dRMS = np.mean(diis_e ** 2) ** 0.5
    print 'SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E   %s' % (SCF_ITER, scf_e, (scf_e - Eold), dRMS, iter_type)
    if (abs(scf_e - Eold) < E_conv) and (dRMS < D_conv):
        break

    Eold = scf_e

    # Build MO fock ,matrix and gradient
    Co = C[:, :ndocc]
    Cv = C[:, ndocc:]
    moF = np.einsum('ui,vj,uv->ij', C, C, F)
    gradient = -4 * moF[:ndocc, ndocc:]
    grad_dot = np.vdot(gradient, gradient)

    if (np.max(np.abs(gradient)) > 0.2):
        F = diis.extrapolate()
        C, Cocc, D = build_orbitals(F)
        iter_type = 'DIIS'
    else:


        # DIIS helper and initial guess
        eps = np.diag(moF)
        precon = -4 * (eps[:ndocc].reshape(-1, 1) - eps[ndocc:])

        x = gradient / precon
        tmp = np.dot(Co.T, F).dot(Co).dot(x)
        tmp -= np.dot(x, Cv.T).dot(F).dot(Cv)
    #    print tmp

        Ax = SCF_Hx(x, moF, Co, Cv)
        r = gradient - Ax
        z = r / precon
        p = z.copy()
        rms = (np.vdot(r, r) / grad_dot) ** 0.5
        if micro_print:
            print('Micro Iteration Guess: Rel. RMS = %1.5e' %  (rms))

        # CG iterations
        for rot_iter in range(max_micro):
            rz_old = np.vdot(r, z)

            Ap = SCF_Hx(p, moF, Co, Cv)
            alpha = rz_old / np.vdot(Ap, p)

            x += alpha * p
            r -= alpha * Ap
            z = r / precon

            rms = (np.vdot(r, r) / grad_dot) ** 0.5

            if micro_print:
                print('Micro Iteration %5d: Rel. RMS = %1.5e' %  (rot_iter + 1, rms))
            if rms < micro_conv:
                break

            beta = np.vdot(r, z) / rz_old
            p = z + beta * p

        U = np.zeros_like(C)
        U[:ndocc, ndocc:] = x
        U[ndocc:, :ndocc] = -x.T

        U += 0.5 * np.dot(U, U)
        U[np.diag_indices_from(U)] += 1

        # Easy acess to shmidt orthogonalization
        U, r = np.linalg.qr(U.T)

        C = np.dot(C, U)
        Cocc = C[:, :ndocc]
        D = np.dot(Cocc, Cocc.T)

        iter_type = 'SOSCF, nmicro ' + str(rot_iter + 1)

print('Total time taken for SCF iterations: %.3f seconds \n' % (time.time()-t))

print('Final SCF energy:     %.8f hartree' % scf_e)

# Compute w/ Psi4 and compare
psi4.set_options({'e_convergence':1e-7,
                  'd_convergence':1e-7})

SCF_E_psi, scf_wfn = psi4.energy('SCF', return_wfn=True)
psi4.driver.p4util.compare_values(SCF_E_psi, scf_e, 6, 'SCF Energy')




