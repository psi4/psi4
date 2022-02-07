#! RHF Hessian code

# -*- coding: utf-8 -*- 
"""
This script calculates nuclear Hessians of RHF Wavefunction using
derivatives of one and two electron integrals obtained from PSI4

Reference: "Derivative studies in hartree-fock and moller-plesset theories",
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

psi4.set_memory(int(1e9), False)
psi4.core.set_output_file('output.dat', False)
psi4.core.set_num_threads(4)

mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")

psi4.core.set_active_molecule(mol)

options = {'BASIS':'STO-3G', 'SCF_TYPE':'PK',
           'E_CONVERGENCE':1e-10,
           'D_CONVERGENCE':1e-10
           }

psi4.set_options(options)

rhf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Assuming C1 symmetry    
occ = wfn.doccpi()[0]
nmo = wfn.nmo()
vir = nmo - occ

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
deriv1 = {}

# 1st Derivative of OEIs 

for atom in range(natoms):
    for key in  oei_dict:
        deriv1_mat[key + str(atom)] = mints.mo_oei_deriv1(oei_dict[key], atom, C, C)
        for p in range(3):
            map_key = key + str(atom) + cart[p]
            deriv1[map_key] = np.asarray(deriv1_mat[key + str(atom)][p])

# 1st Derivative of TEIs 

for atom in range(natoms):
    string = "TEI" + str(atom)
    deriv1_mat[string] = mints.mo_tei_deriv1(atom, C, C, C, C)
    for p in range(3):
        map_key = string + cart[p]
        deriv1[map_key] = np.asarray(deriv1_mat[string][p])

Hes = {}
deriv2_mat = {}
deriv2 = {}

Hes["S"] = np.zeros((3 * natoms, 3 * natoms))
Hes["V"] = np.zeros((3 * natoms, 3 * natoms))
Hes["T"] = np.zeros((3 * natoms, 3 * natoms))
Hes["N"] = np.zeros((3 * natoms, 3 * natoms))
Hes["J"] = np.zeros((3 * natoms, 3 * natoms))
Hes["K"] = np.zeros((3 * natoms, 3 * natoms))
Hes["R"] = np.zeros((3 * natoms, 3 * natoms))
Hessian  = np.zeros((3 * natoms, 3 * natoms))

Hes["N"] = np.asarray(mol.nuclear_repulsion_energy_deriv2())

psi4.core.print_out("\n\n")
Mat = psi4.core.Matrix.from_array(Hes["N"])
Mat.name = "NUCLEAR HESSIAN"
Mat.print_out()

# 2nd Derivative of OEIs 

for atom1 in range(natoms):
    for atom2 in range(atom1 + 1):
        for key in  oei_dict:
            string = key + str(atom1) + str(atom2)
            deriv2_mat[string] = mints.mo_oei_deriv2(oei_dict[key], atom1, atom2, C, C)
            pq = 0
            for p in range(3):
                for q in range(3):
                    map_key = string + cart[p] + cart[q]
                    deriv2[map_key] = np.asarray(deriv2_mat[string][pq])
                    pq = pq+1
                    row = 3 * atom1 + p
                    col = 3 * atom2 + q
                    if key == "S":
                        Hes[key][row][col] = -2.0 * np.einsum("ii,ii->", F[:occ,:occ], deriv2[map_key][:occ,:occ])
                    else:
                        Hes[key][row][col] = 2.0 * np.einsum("ii->", deriv2[map_key][:occ,:occ])
                    Hes[key][col][row] = Hes[key][row][col]
                    Hes[key][col][row] = Hes[key][row][col]


for key in Hes: 
    Mat = psi4.core.Matrix.from_array(Hes[key])
    if key in oei_dict:
        Mat.name = oei_dict[key] + " HESSIAN"
        Mat.print_out()    
        psi4.core.print_out("\n")


# 2nd Derivative of TEIs

for atom1 in range(natoms):
    for atom2 in range(atom1 + 1):
        string = "TEI" + str(atom1) + str(atom2)
        deriv2_mat[string] = mints.mo_tei_deriv2(atom1, atom2, C, C, C, C)
        pq = 0
        for p in range(3):
            for q in range(3):
                map_key = string + cart[p] + cart[q]
                deriv2[map_key] = np.asarray(deriv2_mat[string][pq])
                pq = pq + 1
                row = 3 * atom1 + p
                col = 3 * atom2 + q
                Hes["J"][row][col] =  2.0 * np.einsum("iijj->", deriv2[map_key][:occ,:occ,:occ,:occ])
                Hes["K"][row][col] = -1.0 * np.einsum("ijij->", deriv2[map_key][:occ,:occ,:occ,:occ])

                Hes["J"][col][row] = Hes["J"][row][col]
                Hes["K"][col][row] = Hes["K"][row][col]

JMat = psi4.core.Matrix.from_array(Hes["J"])
KMat = psi4.core.Matrix.from_array(Hes["K"])
JMat.name = " COULOMB  HESSIAN"
KMat.name = " EXCHANGE HESSIAN"
JMat.print_out()    
KMat.print_out()    

# Solve the CPHF equations here,  G_aibj Ubj^x = Bai^x (Einstein summation),
# where G is the electronic hessian,
# G_aibj = delta_ij * delta_ab * epsilon_ij * epsilon_ab + 4 <ij|ab> - <ij|ba> - <ia|jb>, 
# where epsilon_ij = epsilon_i - epsilon_j, (epsilon -> orbital energies),
# x refers to the perturbation, Ubj^x are the corresponsing CPHF coefficients 
# and Bai^x = Sai^x * epsilon_ii - Fai^x + Smn^x  * (2<am|in> - <am|ni>),
# where, S^x =  del(S)/del(x), F^x =  del(F)/del(x).

I_occ = np.diag(np.ones(occ))
I_vir = np.diag(np.ones(vir))
epsilon = np.asarray(wfn.epsilon_a())
eps_diag = epsilon[occ:].reshape(-1, 1) - epsilon[:occ]

#  Build the electronic hessian G

G =  4 * MO[:occ, :occ, occ:, occ:]
G -= MO[:occ, :occ:, occ:, occ:].swapaxes(2,3)
G -= MO[:occ, occ:, :occ, occ:].swapaxes(1,2)
G = G.swapaxes(1,2)
G += np.einsum('ai,ij,ab->iajb', eps_diag, I_occ, I_vir)

# Inverse of G

Ginv = np.linalg.inv(G.reshape(occ * vir, -1))
Ginv = Ginv.reshape(occ,vir,occ,vir)

B = {}
F_grad = {}
U = {}

# Build Fpq^x now

for atom in range(natoms):
    for p in range(3):
        key = str(atom) + cart[p]
        F_grad[key] =  deriv1["T" + key]
        F_grad[key] += deriv1["V" + key]
        F_grad[key] += 2.0 * np.einsum('pqmm->pq', deriv1["TEI" + key][:,:,:occ,:occ])
        F_grad[key] -= 1.0 * np.einsum('pmmq->pq', deriv1["TEI" + key][:,:occ,:occ,:])


psi4.core.print_out("\n\n CPHF Coefficentsn:\n")

# Build Bai^x now

for atom in range(natoms):
    for p in range(3):
        key = str(atom) + cart[p]
        B[key] =  np.einsum("ai,ii->ai", deriv1["S" + key][occ:,:occ], F[:occ,:occ])
        B[key] -= F_grad[key][occ:,:occ]
        B[key] +=  2.0 * np.einsum("amin,mn->ai", MO[occ:,:occ,:occ,:occ], deriv1["S" + key][:occ,:occ])
        B[key] += -1.0 * np.einsum("amni,mn->ai", MO[occ:,:occ,:occ,:occ], deriv1["S" + key][:occ,:occ])

        # Compute U^x now: U_ai^x = G^(-1)_aibj * B_bj^x  

        U[key] = np.einsum("iajb,bj->ai", Ginv, B[key])
        psi4.core.print_out("\n")
        UMat = psi4.core.Matrix.from_array(U[key])
        UMat.name = key 
        UMat.print_out()    


# Build the response hessian now

for atom1 in range(natoms):
    for atom2 in range(atom1+1):
        for p in range(3):
            for q in range(3):
                key1  = str(atom1) + cart[p]
                key2  = str(atom2) + cart[q]
                key1S = "S" + key1
                key2S = "S" + key2
                r = 3 * atom1 + p
                c = 3 * atom2 + q

                Hes["R"][r][c] = -2.0 * np.einsum("ij,ij->", deriv1[key1S][:occ,:occ], F_grad[key2][:occ,:occ])
                Hes["R"][r][c] -= 2.0 * np.einsum("ij,ij->", deriv1[key2S][:occ,:occ], F_grad[key1][:occ,:occ])
                Hes["R"][r][c] += 4.0 * np.einsum("ii,mi,mi->", F[:occ,:occ], deriv1[key2S][:occ,:occ], deriv1[key1S][:occ,:occ])

                Hes["R"][r][c] += 4.0 * np.einsum("ij,mn,imjn->", deriv1[key1S][:occ,:occ], deriv1[key2S][:occ,:occ], MO[:occ,:occ,:occ,:occ])
                Hes["R"][r][c] -= 2.0 * np.einsum("ij,mn,imnj->", deriv1[key1S][:occ,:occ], deriv1[key2S][:occ,:occ], MO[:occ,:occ,:occ,:occ])

                Hes["R"][r][c] -= 4.0 * np.einsum("ai,ai->", U[key2], B[key1])
                Hes["R"][c][r] = Hes["R"][r][c]

Mat = psi4.core.Matrix.from_array(Hes["R"])
Mat.name = " RESPONSE HESSIAN"
Mat.print_out()    

for key in Hes:
    Hessian += Hes[key]

Mat = psi4.core.Matrix.from_array(Hessian)
Mat.name = " TOTAL HESSIAN"
Mat.print_out()    

H_psi4 = psi4.core.Matrix.from_list([
[ 0.07613952269361, -0.00000000000000, -0.00000000000011, -0.03806976134686,  0.00000000000009,  0.00000000000006, -0.03806976134685, -0.00000000000009,  0.00000000000006],
[-0.00000000000000,  0.48290537237134, -0.00000000000000,  0.00000000000005, -0.24145268618572,  0.15890015585447, -0.00000000000005, -0.24145268618572, -0.15890015585447],
[-0.00000000000011, -0.00000000000000,  0.43734495978407,  0.00000000000006,  0.07344233774235, -0.21867247989206,  0.00000000000006, -0.07344233774235, -0.21867247989205],
[-0.03806976134686,  0.00000000000005,  0.00000000000006,  0.04537741758844, -0.00000000000007, -0.00000000000004, -0.00730765624159,  0.00000000000002, -0.00000000000001],
[ 0.00000000000009, -0.24145268618572,  0.07344233774235, -0.00000000000007,  0.25786500659260, -0.11617124679841, -0.00000000000002, -0.01641232040687,  0.04272890905606],
[ 0.00000000000006,  0.15890015585447, -0.21867247989206, -0.00000000000004, -0.11617124679841,  0.19775198076992, -0.00000000000001, -0.04272890905606,  0.02092049912214],
[-0.03806976134685, -0.00000000000005,  0.00000000000006, -0.00730765624159, -0.00000000000002, -0.00000000000001,  0.04537741758844,  0.00000000000007, -0.00000000000004],
[-0.00000000000009, -0.24145268618572, -0.07344233774235,  0.00000000000002, -0.01641232040687, -0.04272890905606,  0.00000000000007,  0.25786500659260,  0.11617124679841],
[ 0.00000000000006, -0.15890015585447, -0.21867247989205, -0.00000000000001,  0.04272890905606,  0.02092049912214, -0.00000000000004,  0.11617124679841,  0.19775198076991],
])

H_python_mat = psi4.core.Matrix.from_array(Hessian)
psi4.compare_matrices(H_psi4, H_python_mat, 10, "RHF-HESSIAN-TEST") #TEST
