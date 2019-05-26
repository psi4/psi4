#! coded full CI vs. detci

import time
import numpy as np
import psi4
np.set_printoptions(precision=7, linewidth=200, threshold=2000, suppress=True)

psi4.set_output_file("output.dat", False)

# Set molecule to dimer
dimer = psi4.geometry("""
Be   0.000000000   0.000000000   0.000000000
""")

psi4.set_options({
   "scf_type"      :  "out_of_core",
   "basis"         :  "cc-pVTZ",
   "e_convergence" :  1e-8,
   "d_convergence" :  1e-8
})

# Build the SCF wavefunction
scf_energy, wfn = psi4.energy('SCF', return_wfn=True)

# Going FCI!
psi4.set_options({"ACTIVE": list(wfn.nmopi())})

# Python based options
num_eig = 1
ctol = 1.e-5
etol = 1.e-7

# Get the CI Wavefunction and compute required integrals
psi4.core.prepare_options_for_module("DETCI")
ciwfn = psi4.core.CIWavefunction(wfn)
print('The number of determinants is %d\n' % ciwfn.ndet())

# Transform integrals required for CI
ciwfn.transform_ci_integrals()

# Build trial vectors
guess_size = 10
print('Using %d determinants in the guess\n' % guess_size)
H = np.array(ciwfn.hamiltonian(guess_size))

gvecs = []
gevals, gevecs = np.linalg.eigh(H)
for x in range(num_eig):
    guess = np.zeros((ciwfn.ndet()))
    guess[:guess_size] = gevecs[:, x]
    gvecs.append(guess)

# Build a few CI vectors
max_guess = 20

Hd = ciwfn.Hd_vector(5)

cvecs = ciwfn.new_civector(max_guess, 200, True, True)
cvecs.set_nvec(max_guess)
cvecs.init_io_files(False)

swork_vec = max_guess
svecs = ciwfn.new_civector(max_guess + 1, 201, True, True)
svecs.set_nvec(max_guess)
svecs.init_io_files(False)

dwork_vec = num_eig
dvecs = ciwfn.new_civector(num_eig + 1, 202, True, True)
dvecs.init_io_files(False)
dvecs.set_nvec(num_eig + 1)
for x in range(num_eig + 1):
    dvecs.write(x, 0)

# Current number of vectors
num_vecs = num_eig

# Copy gvec data into in ci_gvecs
arr_cvecs = np.asarray(cvecs)
for x in range(num_eig):
    arr_cvecs[:] = gvecs[x]
    cvecs.write(x, 0)

# Start parameters
delta_c = 0.0
Eold = scf_energy

G = np.zeros((max_guess, max_guess))

for CI_ITER in range(max_guess - 1):

    # Subspace Matrix, Gij = < bi | H | bj >
    for i in range(num_vecs - num_eig, num_vecs):
        ciwfn.sigma(cvecs, svecs, i, i)
        for j in range(0, num_vecs):
            G[j,i] = G[i, j] = svecs.vdot(cvecs, i, j)

    evals, evecs = np.linalg.eigh(G[:num_vecs, :num_vecs])

    CI_E = evals[0]
    print('CI Iteration %3d: Energy = %4.16f   dE = % 1.5E   dC = %1.5E'
          % (CI_ITER, CI_E, (CI_E - Eold), delta_c))
    if (abs(CI_E - Eold) < etol) and (delta_c < ctol) and (CI_ITER > 3): 
        print('\nCI has converged!')
        break
    Eold = CI_E

    # Build new vectors as linear combinations of the subspace matrix, H
    for n in range(num_eig):

        # Build as linear combinations of previous vectors
        dvecs.zero()
        dvecs.write(dwork_vec, 0)
        for c in range(num_vecs):
            dvecs.axpy(evecs[c, n], cvecs, dwork_vec, c)

        # Build new vector new_vec = ((H * cvec) - evals[n] * cvec) / (evals[n] - Hd)
        ciwfn.sigma(dvecs, svecs, dwork_vec, swork_vec)
        svecs.axpy(-1 * evals[n], dvecs, swork_vec, dwork_vec)
        norm = svecs.dcalc(evals[n], Hd, swork_vec)
        svecs.symnormalize(1 / norm, swork_vec) 
        delta_c = norm

        # Build a new vector that is orthornormal to all previous vectors
        dvecs.copy(svecs, n, swork_vec)
        norm = dvecs.norm(n)
        dvecs.symnormalize(1 / norm, n)
        

        total_proj = 0
        for i in range(num_vecs):
            proj = svecs.vdot(cvecs, swork_vec, i)
            total_proj += proj
            dvecs.axpy(-proj, cvecs, n, i)

        norm = dvecs.norm(n)
        dvecs.symnormalize(1 / norm, n)

        # This *should* screen out contributions that are projected out by above
        if True:
            cvecs.write(num_vecs, 0)
            cvecs.copy(dvecs, num_vecs, n)
            num_vecs += 1


print("\nComparison to Psi4's DETCI module:")
psi_ci_energy, ciwfn = psi4.energy('DETCI', return_wfn=True)
psi4.compare_values(CI_E, psi_ci_energy, 6, 'CI Energy') 
