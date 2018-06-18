#! compare MemJK and DiskJK

import psi4
import numpy as np
import random

mol = psi4.geometry("""
O
H 1 1.00
H 1 1.00 2 103.1
""")

psi4.set_num_threads(6)
memory = 50000
primary = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ")
aux = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ-jkfit")

nbf = primary.nbf()
naux = aux.nbf()

# construct spaces
names = ['C1', 'C2', 'C3', 'C4', 'C5']
sizes = [16, 16, 20, 20, 30]
spaces = {names[ind]: psi4.core.Matrix.from_array(np.random.rand(nbf, size)) for ind, size in enumerate(sizes)}
space_pairs = [[0, 0], [0, 1], [1, 1], [2, 2], [3, 2], [3, 3], [4, 4]]

# space vectors
C_vectors  = [[spaces[names[left]], spaces[names[right]]] for left, right in space_pairs]

# now construct DiskJK and symm_JK objects

# DiskJK
psi4.set_options({"SCF_TYPE" : "DISK_DF"})
DiskJK = psi4.core.JK.build_JK(primary, aux)
DiskJK.initialize()
DiskJK.print_header()

# symm_JK
psi4.set_options({"SCF_TYPE" : "MEM_DF"})
MemJK = psi4.core.JK.build_JK(primary, aux)
MemJK.initialize()
MemJK.print_header()

# add C matrices
for Cleft, Cright in C_vectors:
    DiskJK.C_left_add(Cleft)
    MemJK.C_left_add(Cleft)
    DiskJK.C_right_add(Cright)
    MemJK.C_right_add(Cright)

# compute
DiskJK.compute()
MemJK.compute()

# get integrals
DiskJK_ints = [DiskJK.J(), DiskJK.K()]
MemJK_ints = [MemJK.J(), MemJK.K()]

# compare
for j, t in enumerate(['J', 'K']):
    for i in range(len(DiskJK_ints[0])):
        psi4.compare_arrays(np.asarray(DiskJK_ints[j][i]), 
                            np.asarray(MemJK_ints[j][i]), 9, t + str(i))

