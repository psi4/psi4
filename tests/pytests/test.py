import pytest
import numpy as np
import psi4

# Set NumPy Options
np.set_printoptions(precision=15, linewidth=200, suppress=True)

# Specify Molecule
#hf = psi4.geometry("""
#    F       0.000000000000      0.000000000000      R
#    H       0.000000000000      0.000000000000      0.000000000000
#""")
#f_z = 1.0
h2o = psi4.geometry("""
    O            0.000000000000     0.000000000000     R
    H            0.000000000000    -0.866811832375     0.601435781623
    H            0.000000000000     0.866811832375     0.601435781623
""")
o_z = -0.075791843897

# Specify Basis Set
psi4.set_options({'basis': 'STO-3G'})

# Define h for finite difference method
h = 0.01

# Calculate Angular Momentum Integrals for various O_z coordinates
am_ints = dict()
np_am_ints = dict()
for l in [-2, -1, 0, 1, 2]: 
    h2o.R = o_z + l * h
    psi4.core.set_active_molecule(h2o)

    rhf_e, wfn = psi4.energy('SCF', return_wfn=True)

    mints = psi4.core.MintsHelper(wfn.basisset())
    am_ints[l] = mints.ao_angular_momentum()
    
    natoms = h2o.natom()
    cart = ['_X', '_Y', '_Z']

    np_am_ints[l] = dict()
    for am_cart in range(3):
        map_key = "L" + cart[am_cart]
        np_am_ints[l][map_key] = np.asarray(am_ints[l][am_cart])
# Print Angular Momentum Integrals for various O_z coordinates
#for l in [-2, -1, 0, 1, 2]: 
#    print("l = ", l, "\n")
#    for am_cart in range(3):
#        map_key = "L" + cart[am_cart]
#        print(map_key, ":\n", np_am_ints[l][map_key])
#    print("\n")
        
# Calculate Analytical Derivatives of Angular Momentum Integrals
analytic_deriv = dict()
np_analytic_deriv = dict()
for atom in range(natoms):
    analytic_deriv[str(atom)] = mints.ao_ang_mom_deriv1(atom)
    for am_cart in range(3):
        for atom_cart in range(3):
            map_key = "L" + cart[am_cart] + "_" + str(atom) + cart[atom_cart]
            np_analytic_deriv[map_key] = np.asarray(analytic_deriv[str(atom)][3 * am_cart + atom_cart])
# Print Analytical Derivatives of O_z coordinate
print("np_analytic_deriv[L_X_0_Z]:\n", np_analytic_deriv["L_X_0_Z"], "\n") # dL_x/dO_z
print("np_analytic_deriv[L_Y_0_Z]:\n", np_analytic_deriv["L_Y_0_Z"], "\n") # dL_y/dO_z
print("np_analytic_deriv[L_Z_0_Z]:\n", np_analytic_deriv["L_Z_0_Z"], "\n") # dL_z/dO_z

# Calculate 5-point Finite Difference of Angular Momentum Integrals
findif_deriv = dict()
for am_cart in range(3):
    map_key = "L" + cart[am_cart]
    findif_deriv[map_key] = ((8 * np_am_ints[1][map_key] - 8 * np_am_ints[-1][map_key] - np_am_ints[2][map_key] + np_am_ints[-2][map_key]) / (12 * h))
# Print Finite Difference Results
for i in findif_deriv:
    print("findif_deriv[", i, "]:\n", findif_deriv[i], "\n")

