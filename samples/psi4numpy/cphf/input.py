#! Tests out the CG solver with CPHF Polarizabilities

import time
import numpy as np
import psi4

psi4.set_output_file("output.dat")

# Benzene
mol = psi4.geometry("""
    0 1
    O          0.000000000000     0.000000000000    -0.075791843589
    H          0.000000000000    -0.866811828967     0.601435779270
    H          0.000000000000     0.866811828967     0.601435779270
    symmetry c1
""")

psi4.set_options({"basis":         "aug-cc-pVDZ",
                  "scf_type":      "df",
                  "e_convergence": 1e-8,
                  "save_jk":       True,
                 })


scf_e, scf_wfn = psi4.energy("SCF", return_wfn=True)

# Orbitals
Co = scf_wfn.Ca_subset("AO", "OCC")
Cv = scf_wfn.Ca_subset("AO", "VIR")

# Mints object
mints = psi4.core.MintsHelper(scf_wfn.basisset())

# RHS Dipoles
dipoles_xyz = []
for dip in mints.ao_dipole():
    Fia = psi4.core.triplet(Co, dip, Cv, True, False, False)
    Fia.scale(-2.0)
    dipoles_xyz.append(Fia)

# Build up the preconditioner
precon = psi4.core.Matrix(Co.shape[1], Cv.shape[1])
occ = np.array(scf_wfn.epsilon_a_subset("AO", "OCC"))
vir = np.array(scf_wfn.epsilon_a_subset("AO", "VIR"))
precon.np[:] = (-occ.reshape(-1, 1) + vir)

# Build a preconditioner function
def precon_func(matrices, active_mask):
    ret = []
    for act, mat in zip(active_mask, matrices):
        if act:
            p = mat.clone()
            p.apply_denominator(precon)
            ret.append(p)
        else:
            ret.append(False)
    return ret 

def wrap_Hx(matrices, active_mask):
    x_vec = [mat for act, mat in zip(active_mask, matrices) if act]
    Hx_vec = scf_wfn.cphf_Hx(x_vec)

    ret = []
    cnt = 0
    for act, mat in zip(active_mask, matrices):
        if act:
            ret.append(Hx_vec[cnt])
            cnt += 1
        else:
            ret.append(False)

    return ret

# Solve
ret, resid = psi4.p4util.solvers.cg_solver(dipoles_xyz, wrap_Hx, precon_func, rcond=1.e-6)

polar = np.empty((3, 3))
for numx in range(3):
    for numf in range(3):
        polar[numx, numf] = -1 * ret[numx].vector_dot(dipoles_xyz[numf])

psi4.core.print_out("\n         " + "CPHF Dipole Polarizability:".center(44) + "\n")
tops = ("X", "Y", "Z")
psi4.core.print_out("       %12s %12s %12s\n" % tops)
for n, p in enumerate(tops):
    psi4.core.print_out("      %3s %12.4f %12.4f %12.4f\n" % (p, polar[n][0], polar[n][1], polar[n][2]))
psi4.core.print_out("\n")
    
