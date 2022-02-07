#! rhf gradient code

"""
This script calculates nuclear gradients of RHF Wavefunction using
gradients of one and two electron integrals obtained from PSI4. 

Reference: "Derivative studies in Hartree--Fock and Moller--Plesset theories",
J. A. Pople, R. Krishnan, H. B. Schlegel and J. S. Binkley
DOI: 10.1002/qua.560160825
"""

__authors__ = "Ashutosh Kumar"
__credits__ = ["Ashutosh Kumar"]

__copyright__ = "(c) 2014-2017, The Psi4NumPy Developers"
__license__ = "BSD-3-Clause"
__date__ = "2017-12-17"

import time
import numpy as np
np.set_printoptions(precision=15, linewidth=200, suppress=True)
import psi4

psi4.set_output_file("output.dat", False)

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")

psi4.core.set_active_molecule(mol)

options = {'BASIS':'STO-3G', 'SCF_TYPE':'PK',
           'E_CONVERGENCE':1e-10,
           'D_CONVERGENCE':1e-10}

psi4.set_options(options)


rhf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Assuming C1 symmetry    
occ = wfn.doccpi()[0]
nmo = wfn.nmo()

C = wfn.Ca_subset("AO", "ALL")
npC = np.asarray(C)

mints = psi4.core.MintsHelper(wfn.basisset())
H_ao = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())

# Update H, transform to MO basis 
H = np.einsum('uj,vi,uv', npC, npC, H_ao)

# Integral generation from Psi4's MintsHelper
MO = np.asarray(mints.mo_eri(C, C, C, C))
# Physicist notation    
MO = MO.swapaxes(1,2)

F = H + 2.0 * np.einsum('pmqm->pq', MO[:, :occ, :, :occ])
F -= np.einsum('pmmq->pq', MO[:, :occ, :occ, :])
natoms = mol.natom()
cart = ['_X', '_Y', '_Z']
oei_dict = {"S" : "OVERLAP", "T" : "KINETIC", "V" : "POTENTIAL"}

deriv1_mat = {}
deriv1_np = {}

Gradient = {}

Gradient["N"] = np.zeros((natoms, 3))
Gradient["S"] = np.zeros((natoms, 3))
Gradient["S'"] = np.zeros((natoms, 3))
Gradient["V"] = np.zeros((natoms, 3))
Gradient["T"] = np.zeros((natoms, 3))
Gradient["J"] = np.zeros((natoms, 3))
Gradient["K"] = np.zeros((natoms, 3))
Gradient["Total"] = np.zeros((natoms, 3))

Gradient["N"] = np.asarray(mol.nuclear_repulsion_energy_deriv1([0,0,0]))

psi4.core.print_out("\n\n")
Mat = psi4.core.Matrix.from_array(Gradient["N"])
Mat.name = "NUCLEAR GRADIENT"
Mat.print_out()


# 1st Derivative of OEIs 

for atom in range(natoms):
    for key in  oei_dict:
        deriv1_mat[key + str(atom)] = mints.mo_oei_deriv1(oei_dict[key], atom, C, C)
        for p in range(3):
            map_key = key + str(atom) + cart[p]
            deriv1_np[map_key] = np.asarray(deriv1_mat[key + str(atom)][p])
            if key == "S":
                Gradient[key][atom, p] = -2.0 * np.einsum("ii,ii->", F[:occ,:occ], deriv1_np[map_key][:occ,:occ])
                Gradient["S'"][atom, p] = 2.0 * np.einsum("ii->", deriv1_np[map_key][:occ,:occ]) # For comparison with PSI4's overlap_grad
            else:
                Gradient[key][atom, p] = 2.0 * np.einsum("ii->", deriv1_np[map_key][:occ,:occ])

psi4.core.print_out("\n\n OEI Gradients\n\n")
for key in Gradient: 
    Mat = psi4.core.Matrix.from_array(Gradient[key])
    if key in oei_dict:
        Mat.name = oei_dict[key] + " GRADIENT"
        Mat.print_out()    
        psi4.core.print_out("\n")


Gradient["J"] = np.zeros((natoms, 3))
Gradient["K"] = np.zeros((natoms, 3))

# 1st Derivative of TEIs 

for atom in range(natoms):
    string = "TEI" + str(atom)
    deriv1_mat[string] = mints.mo_tei_deriv1(atom, C, C, C, C)
    for p in range(3):
        map_key = string + cart[p]
        deriv1_np[map_key] = np.asarray(deriv1_mat[string][p])
        Gradient["J"][atom, p] =  2.0 * np.einsum("iijj->", deriv1_np[map_key][:occ,:occ,:occ,:occ])
        Gradient["K"][atom, p] = -1.0 * np.einsum("ijij->", deriv1_np[map_key][:occ,:occ,:occ,:occ])

psi4.core.print_out("\n\n TEI Gradients\n\n")
JMat = psi4.core.Matrix.from_array(Gradient["J"])
KMat = psi4.core.Matrix.from_array(Gradient["K"])
JMat.name = " COULOMB  GRADIENT"
KMat.name = " EXCHANGE GRADIENT"
JMat.print_out()    
KMat.print_out()    

Gradient["OEI"] = Gradient["S"] + Gradient["V"] + Gradient["T"] 
Gradient["TEI"] = Gradient["J"] + Gradient["K"]
Gradient["Total"] = Gradient["OEI"] + Gradient["TEI"] + Gradient["N"]


# PIS4's overlap_grad, kinetic_grad and potential_grad

PSI4_Grad = {}
D = wfn.Da()
D.add(wfn.Db())

PSI4_Grad["S"] = mints.overlap_grad(D) 
PSI4_Grad["T"] = mints.kinetic_grad(D)
PSI4_Grad["V"] = mints.potential_grad(D)

#Convert np array into PSI4 Matrix 
G_python_S_mat = psi4.core.Matrix.from_array(Gradient["S'"])
G_python_T_mat = psi4.core.Matrix.from_array(Gradient["T"])
G_python_V_mat = psi4.core.Matrix.from_array(Gradient["V"])

# Test OEI gradients with that of PSI4

# PSI4's Total Gradient 
Total_G_psi4 = psi4.core.Matrix.from_list([                                     
             [ 0.000000000000,     0.00000000000000,    -0.09744143723018],
             [ 0.000000000000,    -0.08630009812231,     0.04872071861516],
             [ 0.000000000000,     0.08630009812231,     0.04872071861516],
       ])
G_python_Total_mat = psi4.core.Matrix.from_array(Gradient["Total"])
