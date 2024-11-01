#! A simple Psi 4 input script to compute a SCF reference using Psi4's libJK

import time

import numpy as np
import psi4

psi4.set_output_file("output.dat", False)

# Benzene
mol = psi4.geometry("""
0 1
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")

psi4.set_options({"basis":         "aug-cc-pVDZ",
                  "scf_type":      "df",
                  "e_convergence": 1e-8
                 })

# Set tolerances
maxiter = 12
E_conv = 1.0E-6
D_conv = 1.0E-5

# Integral generation from Psi4's MintsHelper
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
mints = psi4.core.MintsHelper(wfn.basisset())
S = mints.ao_overlap()

# Get nbf and ndocc for closed shell molecules
nbf = wfn.nso()
ndocc = wfn.nalpha()
if wfn.nalpha() != wfn.nbeta():
    raise PsiException("Only valid for RHF wavefunctions!")

psi4.core.print_out('\nNumber of occupied orbitals: %d\n' % ndocc)
psi4.core.print_out('Number of basis functions:   %d\n\n' % nbf)

# Build H_core
V = mints.ao_potential()
T = mints.ao_kinetic()
H = T.clone()
H.add(V)

# Orthogonalizer A = S^(-1/2)
A = mints.ao_overlap()
A.power(-0.5, 1.e-16)

# Diagonalize routine
def build_orbitals(diag):
    Fp = psi4.core.triplet(A, diag, A, True, False, True)

    Cp = psi4.core.Matrix(nbf, nbf)
    eigvecs = psi4.core.Vector(nbf)
    Fp.diagonalize(Cp, eigvecs, psi4.core.DiagonalizeOrder.Ascending)

    C = psi4.core.doublet(A, Cp, False, False)

    Cocc = psi4.core.Matrix(nbf, ndocc)
    Cocc.np[:] = C.np[:, :ndocc]

    D = psi4.core.doublet(Cocc, Cocc, False, True)
    return C, Cocc, D

# Build core orbitals
C, Cocc, D = build_orbitals(H)

# Setup data for DIIS
t = time.time()
E = 0.0
Enuc = mol.nuclear_repulsion_energy()
Eold = 0.0

# Initialize the JK object
jk = psi4.core.JK.build(wfn.basisset())
jk.set_memory(int(1.25e8))  # 1GB
jk.initialize()
jk.print_header()

diis_obj = psi4.p4util.solvers.DIIS(max_vec=3, removal_policy="largest")

psi4.core.print_out('\nTotal time taken for setup: %.3f seconds\n' % (time.time() - t))

psi4.core.print_out('\nStart SCF iterations:\n\n')
t = time.time()

for SCF_ITER in range(1, maxiter + 1):

    # Compute JK
    jk.C_left_add(Cocc)
    jk.compute()
    jk.C_clear()

    # Build Fock matrix
    F = H.clone()
    F.axpy(2.0, jk.J()[0])
    F.axpy(-1.0, jk.K()[0])

    # DIIS error build and update
    diis_e = psi4.core.triplet(F, D, S, False, False, False)
    diis_e.subtract(psi4.core.triplet(S, D, F, False, False, False))
    diis_e = psi4.core.triplet(A, diis_e, A, False, False, False)

    diis_obj.add(F, diis_e)

    # SCF energy and update
    FH = F.clone()
    FH.add(H)
    SCF_E = FH.vector_dot(D) + Enuc

    dRMS = diis_e.rms()

    psi4.core.print_out('SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E\n'
          % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))
    if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
        break

    Eold = SCF_E

    # DIIS extrapolate
    F = diis_obj.extrapolate()

    # Diagonalize Fock matrix
    C, Cocc, D = build_orbitals(F)

    if SCF_ITER == maxiter:
        psi4.clean()
        raise Exception("Maximum number of SCF cycles exceeded.\n")

psi4.core.print_out('Total time for SCF iterations: %.3f seconds \n\n' % (time.time() - t))
#print(psi4.energy("SCF"))

psi4.core.print_out('Final SCF energy: %.8f hartree\n' % SCF_E)
psi4.compare_values(-76.0033389840197202, SCF_E, 6, 'SCF Energy')
