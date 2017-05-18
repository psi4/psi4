import numpy as np
import psi4

psi4.core.set_output_file("output.dat", False)


mol = psi4.geometry("""
He 0 0 1
He 0 0  -1
no_reorient
symmetry c1
""")

e, wfn = psi4.energy("PBE/cc-pVDZ", molecule=mol, return_wfn=True)

# Yank some things from the wavefunction
V = wfn.V_potential()
D = np.asarray(wfn.Da())
points_func = V.properties()
superfunc = V.functional()

print("N blocks: %d\n" % V.nblocks())

# Grab a random block for now
grid = V.get_block(20)
x = np.array(grid.x())
y = np.array(grid.y())
z = np.array(grid.z())
w = np.array(grid.w())
npoints = w.shape[0]
nlocal = len(grid.functions_local_to_global())

# Boost up deriv
points_func.set_deriv(2)
points_func.compute_points(grid)

# Grab phi's
phi = np.array(points_func.basis_values()["PHI"])
phi_x = np.array(points_func.basis_values()["PHI_X"])
phi_y = np.array(points_func.basis_values()["PHI_Y"])
phi_z = np.array(points_func.basis_values()["PHI_Z"])

phi_xx = np.array(points_func.basis_values()["PHI_XX"])
phi_yy = np.array(points_func.basis_values()["PHI_YY"])
phi_zz = np.array(points_func.basis_values()["PHI_ZZ"])
phi_xy = np.array(points_func.basis_values()["PHI_XY"])
phi_xz = np.array(points_func.basis_values()["PHI_XZ"])
phi_yz = np.array(points_func.basis_values()["PHI_YZ"])

tmp = np.dot(phi, D)

# Rho
rho = np.einsum('ab,bc,ca->a', phi, D, phi.T)
print('Rho close? %s' % np.allclose(rho, points_func.point_values()["RHO_A"]))

# Gradient of rho
rho_x = 2.0 * np.einsum('ab,bc,ca->a', phi, D, phi_x.T)
rho_y = 2.0 * np.einsum('ab,bc,ca->a', phi, D, phi_y.T)
rho_z = 2.0 * np.einsum('ab,bc,ca->a', phi, D, phi_z.T)
gamma = rho_x ** 2 + rho_y ** 2 + rho_z ** 2
print np.array(points_func.point_values()["GAMMA_AA"])[:20]
print gamma[:20]
print np.linalg.norm(gamma - points_func.point_values()["GAMMA_AA"])
print('Grad close? %s' % np.allclose(gamma, points_func.point_values()["GAMMA_AA"]))

# Build the laplacian of rho
rho_xx = np.einsum('ab,bc,ac->a', phi, D, phi_xx)
rho_xx += np.einsum('ab,bc,ac->a', phi_x, D, phi_x)
rho_xx *= 2

rho_yy = np.einsum('ab,bc,ac->a', phi, D, phi_yy)
rho_yy += np.einsum('ab,bc,ac->a', phi_y, D, phi_y)
rho_yy *= 2

rho_zz = np.einsum('ab,bc,ac->a', phi, D, phi_zz)
rho_zz += np.einsum('ab,bc,ac->a', phi_z, D, phi_z)
rho_zz *= 2

rho_xy = np.einsum('ab,bc,ac->a', phi, D, phi_xy)
rho_xy += np.einsum('ab,bc,ac->a', phi_x, D, phi_y)
rho_xy *= 2

rho_xz = np.einsum('ab,bc,ac->a', phi, D, phi_xz)
rho_xz += np.einsum('ab,bc,ac->a', phi_x, D, phi_z)
rho_xz *= 2

rho_yz = np.einsum('ab,bc,ac->a', phi, D, phi_yz)
rho_yz += np.einsum('ab,bc,ac->a', phi_y, D, phi_z)
rho_yz *= 2

rho_lapl  = rho_xx ** 2 
rho_lapl += rho_yy ** 2 
rho_lapl += rho_zz ** 2 
#rho_lapl += 2.0 * rho_xy ** 2 
#rho_lapl += 2.0 * rho_xz ** 2 
#rho_lapl += 2.0 * rho_yz ** 2


# Lets do a GGA functional
Varr = np.zeros((10, 10))
for x in range(V.nblocks()):
    grid = V.get_block(x)
    w = np.array(grid.w())
    npoints = w.shape[0]

    points_func.compute_points(grid)
    ret = superfunc.compute_functional(points_func.point_values(), -1)

    phi = np.array(points_func.basis_values()["PHI"])[:npoints]
    phi_x = np.array(points_func.basis_values()["PHI_X"])[:npoints]
    phi_y = np.array(points_func.basis_values()["PHI_Y"])[:npoints]
    phi_z = np.array(points_func.basis_values()["PHI_Z"])[:npoints]

    rho_x = np.array(points_func.point_values()["RHO_AX"])[:npoints]
    rho_y = np.array(points_func.point_values()["RHO_AY"])[:npoints]
    rho_z = np.array(points_func.point_values()["RHO_AZ"])[:npoints]

    v_rho_a = np.array(ret["V_RHO_A"])[:npoints]
    v_gamma_aa = np.array(ret["V_GAMMA_AA"])[:npoints]
    v_gamma_ab = np.array(ret["V_GAMMA_AB"])[:npoints]
    v_gamma = (2.0 * v_gamma_aa + v_gamma_ab)[:npoints]
    
    tmp = np.zeros((npoints, 10))

    # LDA
    tmpv  = 0.5 * np.einsum('pb,p,p,pa->ab', phi, v_rho_a, w, phi)

    # GGA
    tmpv += np.einsum('pb,p,p,p,pa->ab', phi, v_gamma, rho_x, w, phi_x)
    tmpv += np.einsum('pb,p,p,p,pa->ab', phi, v_gamma, rho_y, w, phi_y)
    tmpv += np.einsum('pb,p,p,p,pa->ab', phi, v_gamma, rho_z, w, phi_z)

    # Meta
    #tmpv += np.einsum('pb,p,p,p,pa->ab', phi, v_lapl, rho_xx, w, phi_xx)
    #tmpv += np.einsum('pb,p,p,p,pa->ab', phi_x, v_lapl, rho_xx, w, phi_x)
    # ...



    Varr += tmpv + tmpv.T

print('\nMy build:')
print Varr
print('\nPsi build:')
print np.array(wfn.Va())
print('\nAllclose? %s' % np.allclose(Varr, np.array(wfn.Va())))

