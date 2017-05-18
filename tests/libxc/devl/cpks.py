import numpy as np
import psi4
import time

psi4.core.set_output_file("output.dat", False)


mol = psi4.geometry("""
He 0 0 1
He 0 0  -1
no_reorient
symmetry c1
""")

psi4.core.set_global_option("SCF_TYPE", "PK")
e, wfn = psi4.energy("B3LYP/6-31G", molecule=mol, return_wfn=True)

# Yank some things from the wavefunction
V = wfn.V_potential()
D = np.asarray(wfn.Da())
C = wfn.Ca()
Co = wfn.Ca_subset("AO", "OCC")
Cv = wfn.Ca_subset("AO", "VIR")
nbf = wfn.nmo()
nocc = wfn.nalpha()
nvir = nbf - nocc

epsilon = np.asarray(wfn.epsilon_a())
points_func = V.properties()
superfunc = V.functional()
superfunc.set_deriv(2)
x_omega = superfunc.x_omega()

print("N blocks: %d" % V.nblocks())
# Lets do a LDA functional
Varr = np.zeros((nbf, nbf))
for x in range(V.nblocks()):
    grid = V.get_block(x)
    w = np.array(grid.w())
    npoints = w.shape[0]

    points_func.compute_points(grid)
    ret = superfunc.compute_functional(points_func.point_values(), -1)

    phi = np.array(points_func.basis_values()["PHI"])[:npoints]

    v_rho_a = np.array(ret["V_RHO_A"])[:npoints]
    
    tmp = np.zeros((npoints, nbf))

    # LDA
    tmpv  = 0.5 * np.einsum('pb,p,p,pa->ab', phi, v_rho_a, w, phi)

    Varr += tmpv + tmpv.T

# Build a V
print('\nV Allclose? %s' % np.allclose(Varr, np.array(wfn.Va())))

def build_XCderiv(k):
    # Lets do a LDA functional
    Varr = np.zeros((nbf, nbf))
    for x in range(V.nblocks()):
        grid = V.get_block(x)
        w = np.array(grid.w())
        npoints = w.shape[0]
    
        phi = np.array(points_func.basis_values()["PHI"])[:npoints]
        phi_x = np.array(points_func.basis_values()["PHI_X"])[:npoints] 
        phi_y = np.array(points_func.basis_values()["PHI_Y"])[:npoints] 
        phi_z = np.array(points_func.basis_values()["PHI_Z"])[:npoints] 
        rho = np.einsum('ab,bc,ca->a', phi, k, phi.T)

        inp = {
            'RHO_A' : psi4.core.Vector.from_array(rho),
            'RHO_B' : psi4.core.Vector.from_array(rho),
        }
    
        points_func.compute_points(grid)
        ret = superfunc.compute_functional(points_func.point_values(), -1)
        #ret = superfunc.compute_functional(inp, -1)
        # [u'V_RHO_B_RHO_B', u'V_RHO_A', u'V_RHO_B', u'V', u'V_RHO_A_RHO_B', u'V_RHO_A_RHO_A']
#        print ret.keys()
    
        v_rho_a = np.array(ret["V_RHO_A"])[:npoints]
        #v_rho_a = np.array(ret["V_RHO_A_RHO_A"])[:npoints]
        #v_rho_a = np.array(ret["V_RHO_A"])[:npoints]
        
        tmp = np.zeros((npoints, nbf))
   
        # LDA
        tmpv  = 0.5 * np.einsum('pb,p,p,pa->ab', phi, v_rho_a, w, phi)
        tmpv  = 0.5 * np.einsum('pb,p,p,pa->ab', phi, v_rho_a, w, phi)
        tmpv  = 0.5 * np.einsum('pb,p,p,pa->ab', phi, v_rho_a, w, phi)
    
        Varr += tmpv + tmpv.T
    print Varr
    print k * 1.e6
    exit()
    #print k
    return Varr * k
    


print('\nBuild CPKS objects')
maxiter = 3
conv = 1.e-6

# Grab perturbation tensors in MO basis
nCo = np.asarray(Co)
nCv = np.asarray(Cv)
mints = psi4.core.MintsHelper(wfn.basisset())
eri = mints.ao_eri()
tmp_dipoles = mints.so_dipole()
dipoles_xyz = []

for num in range(3):
    Fso = np.asarray(tmp_dipoles[num])
    Fia = (nCo.T).dot(Fso).dot(nCv)
    Fia *= -2
    dipoles_xyz.append(Fia)

# Build initial guess, previous vectors, diis object, and C_left updates
x = []
x_old = []
diis = []
ia_denom = - epsilon[:nocc].reshape(-1, 1) + epsilon[nocc:]
for xyz in range(3):
    x.append(dipoles_xyz[xyz] / ia_denom)
    x_old.append(np.zeros(ia_denom.shape))

# Convert Co and Cv to numpy arrays
mCo = Co
Co = np.asarray(Co)
Cv = np.asarray(Cv)
Va = np.array(wfn.Va())

print('\nStarting CPHF iterations:')
t = time.time()
for CPHF_ITER in range(1, maxiter + 1):

    # Update amplitudes
    for xyz in range(3):
        # Build J and K objects
        Kao = -np.dot(Co, x[xyz]).dot(Cv.T)
        J = np.einsum('pqrs,rs->pq', eri, Kao)
        K = np.einsum('prqs,rs->pq', eri, Kao)

        # Bulid new guess
        X = dipoles_xyz[xyz].copy()
        X -= x_omega * (Co.T).dot(4 * J - K.T - K).dot(Cv)
        #X += (Co.T).dot(Va * Kao).dot(Cv) 
        X += (Co.T).dot(build_XCderiv(Kao)).dot(Cv) 
        #print X
        X /= ia_denom
        #print Fia
        #exit() 

        # DIIS for good measure
        x[xyz] = X.copy()

    # Check for convergence
    rms = []
    for xyz in range(3):
        rms.append(np.max((x[xyz] - x_old[xyz]) ** 2))
        x_old[xyz] = x[xyz]

    avg_RMS = sum(rms) / 3
    max_RMS = max(rms)

    if max_RMS < conv:
        print('CPHF converged in %d iterations and %.2f seconds.' % (CPHF_ITER, time.time() - t))
        break

    print('CPHF Iteration %3d: Average RMS = %3.8f  Maximum RMS = %3.8f' %
            (CPHF_ITER, avg_RMS, max_RMS))


# Compute 3x3 polarizability tensor
polar = np.empty((3, 3))
for numx in range(3):
    for numf in range(3):
        polar[numx, numf] = -1 * np.einsum('ia,ia->', x[numx], dipoles_xyz[numf])
print('\nB3LYP Dipole Polarizability:')
print(np.around(polar, 5))
print("\nHF Dipole Polarizability:")
print("""[[-0.61569 -0.       0.     ]
 [-0.      -0.61569  0.     ]
 [ 0.       0.      -0.68648]]""")
