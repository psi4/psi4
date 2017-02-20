import psi4
import numpy as np


mol = psi4.geometry("""
C
C 1 20
""")

primary = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ")
aux = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ-jkfit")

nbf = primary.nbf()
naux = aux.nbf()
print (nbf)

# form metric
mints = psi4.core.MintsHelper(primary)
zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
Jmetric = np.squeeze(mints.ao_eri(zero_bas, aux, zero_bas, aux))

# form inverse metric
Jmetric_inv = mints.ao_eri(zero_bas, aux, zero_bas, aux)
Jmetric_inv.power(-0.5, 1.e-12)
Jmetric_inv = np.squeeze(Jmetric_inv)

# form Qpq
Qpq = np.squeeze(mints.ao_eri(aux, zero_bas, primary, primary))
Qpq = Jmetric_inv.dot(Qpq.reshape(naux, -1))
Qpq = Qpq.reshape(naux, nbf, -1) 

# space sizes
c1 = 8
c2 = 16
c3 = 20
c4 = 24  

# space initiations
C1 = psi4.core.Matrix(nbf,c1) 
C2 = psi4.core.Matrix(nbf,c2) 
C3 = psi4.core.Matrix(nbf,c3)
C4 = psi4.core.Matrix(nbf,c4)

# get random numpy arrays
C1.np[:] = np.random.random()
C2.np[:] = np.random.random()
C3.np[:] = np.random.random()
C4.np[:] = np.random.random()

# set transformed integrals
Qmo = []
Qmo.append(np.zeros((naux, c1, c2)))
Qmo.append(np.zeros((naux, c3, c2)))
Qmo.append(np.zeros((naux, c4, c2)))
Qmo.append(np.zeros((naux, c1, c4)))
Qmo.append(np.zeros((naux, c3, c1)))
Qmo.append(np.zeros((naux, c3, c3)))

# transform
for i in range(naux):
    Qmo[0][i] = np.dot(C1.np.T, Qpq[i]).dot(C2)
    Qmo[1][i] = np.dot(C3.np.T, Qpq[i]).dot(C2)
    Qmo[2][i] = np.dot(C4.np.T, Qpq[i]).dot(C2)
    Qmo[3][i] = np.dot(C1.np.T, Qpq[i]).dot(C4)
    Qmo[4][i] = np.dot(C3.np.T, Qpq[i]).dot(C1)
    Qmo[5][i] = np.dot(C3.np.T, Qpq[i]).dot(C3)

# let's do two of each "Qpq, pQq, pqQ"
Qmo[2] = np.einsum("Qpq->pQq", Qmo[2])
Qmo[3] = np.einsum("Qpq->pQq", Qmo[3])

Qmo[4] = np.einsum("Qpq->pqQ", Qmo[4])
Qmo[5] = np.einsum("Qpq->pqQ", Qmo[5])

# now compare to DF_Helper
dfh = psi4.core.DF_Helper(primary, aux)

# tweak options
dfh.set_method("STORE")
dfh.set_memory(50000)
dfh.set_on_core(False)

# build
dfh.initialize()

# set spaces
dfh.add_space("C1", C1)
dfh.add_space("C2", C2)
dfh.add_space("C3", C3)
dfh.add_space("C4", C4)

# add transformations
dfh.add_transformation("Qmo1", "C1", "C2")      # best on left  (Q|bw)
dfh.add_transformation("Qmo2", "C3", "C2")      # best on right (Q|wb)
dfh.add_transformation("Qmo3", "C4", "C2")      # best on right (w|Qb)
dfh.add_transformation("Qmo4", "C1", "C4")      # best on left  (b|Qw)
dfh.add_transformation("Qmo5", "C3", "C1")      # best on right (wb|Q)
dfh.add_transformation("Qmo6", "C3", "C3")      # best on left  (bw|Q)

# invoke transformations
dfh.transform()

# tranpose necessary tensors
dfh.transpose("Qmo3", (1, 0, 2))
dfh.transpose("Qmo4", (1, 0, 2))
dfh.transpose("Qmo5", (1, 2, 0))
dfh.transpose("Qmo6", (1, 2, 0))

# grab transformed integrals
dfh_Qmo = []
dfh_Qmo.append(dfh.get_tensor("Qmo1"))
dfh_Qmo.append(dfh.get_tensor("Qmo2"))
dfh_Qmo.append(dfh.get_tensor("Qmo3"))
dfh_Qmo.append(dfh.get_tensor("Qmo4"))
dfh_Qmo.append(dfh.get_tensor("Qmo5"))
dfh_Qmo.append(dfh.get_tensor("Qmo6"))

# am i right?
for i in range(6):
    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "STORE with some tranposes - on disc" )     #TEST

# try again without tranposing --
# untranspose python side
Qmo[2] = np.einsum("pQq->Qpq", Qmo[2])
Qmo[3] = np.einsum("pQq->Qpq", Qmo[3])
Qmo[4] = np.einsum("pqQ->Qpq", Qmo[4])
Qmo[5] = np.einsum("pqQ->Qpq", Qmo[5])

# invoke transformations
dfh.transform()

# grab transformed integrals
dfh_Qmo[0] = (dfh.get_tensor("Qmo1"))
dfh_Qmo[1] = (dfh.get_tensor("Qmo2"))
dfh_Qmo[2] = (dfh.get_tensor("Qmo3"))
dfh_Qmo[3] = (dfh.get_tensor("Qmo4"))
dfh_Qmo[4] = (dfh.get_tensor("Qmo5"))
dfh_Qmo[5] = (dfh.get_tensor("Qmo6"))

# am i right?
for i in range(6):
    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "STORE w/o transposes- on disc" )     #TEST

# Let's use new orbital spaces
C1.np[:] = np.random.random()
C2.np[:] = np.random.random()
C3.np[:] = np.random.random()
C4.np[:] = np.random.random()

# adding a new space overwrites the previous one
dfh.add_space("C1", C1)
dfh.add_space("C2", C2)
dfh.add_space("C3", C3)
dfh.add_space("C4", C4)

# recompute python side
for i in range(naux):
    Qmo[0][i] = np.dot(C1.np.T, Qpq[i]).dot(C2)
    Qmo[1][i] = np.dot(C3.np.T, Qpq[i]).dot(C2)
    Qmo[2][i] = np.dot(C4.np.T, Qpq[i]).dot(C2)
    Qmo[3][i] = np.dot(C1.np.T, Qpq[i]).dot(C4)
    Qmo[4][i] = np.dot(C3.np.T, Qpq[i]).dot(C1)
    Qmo[5][i] = np.dot(C3.np.T, Qpq[i]).dot(C3)

# recompute using df_helper
dfh.transform()

# grab new transformed integrals
dfh_Qmo[0] = (dfh.get_tensor("Qmo1"))
dfh_Qmo[1] = (dfh.get_tensor("Qmo2"))
dfh_Qmo[2] = (dfh.get_tensor("Qmo3"))
dfh_Qmo[3] = (dfh.get_tensor("Qmo4"))
dfh_Qmo[4] = (dfh.get_tensor("Qmo5"))
dfh_Qmo[5] = (dfh.get_tensor("Qmo6"))

# am i right?
for i in range(6):
    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "STORE with new orbital spaces- on disc" )     #TEST

# okay, let's try rebuilding using direct ---------------------------------------------------------------------
# destoy the original instance!
del dfh

# redeclare
dfh = psi4.core.DF_Helper(primary, aux)

# tweak options
dfh.set_method("DIRECT")
dfh.set_memory(50000)
dfh.set_on_core(False)
#dfh.hold_met(True) #works

# build
dfh.initialize()

# set spaces
dfh.add_space("C1", C1)
dfh.add_space("C2", C2)
dfh.add_space("C3", C3)
dfh.add_space("C4", C4)

# add transformations
dfh.add_transformation("Qmo1", "C1", "C2")      # best on left  (Q|bw)
dfh.add_transformation("Qmo2", "C3", "C2")      # best on right (Q|wb)
dfh.add_transformation("Qmo3", "C4", "C2")      # best on right (w|Qb)
dfh.add_transformation("Qmo4", "C1", "C4")      # best on left  (b|Qw)
dfh.add_transformation("Qmo5", "C3", "C1")      # best on right (wb|Q)
dfh.add_transformation("Qmo6", "C3", "C3")      # best on left  (bw|Q)

# invoke transformations
dfh.transform()

# grab transformed integrals
dfh_Qmo[0] = (dfh.get_tensor("Qmo1"))
dfh_Qmo[1] = (dfh.get_tensor("Qmo2"))
dfh_Qmo[2] = (dfh.get_tensor("Qmo3"))
dfh_Qmo[3] = (dfh.get_tensor("Qmo4"))
dfh_Qmo[4] = (dfh.get_tensor("Qmo5"))
dfh_Qmo[5] = (dfh.get_tensor("Qmo6"))

# am i right?
for i in range(6):
    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "DIRECT - on disc" )     #TEST


# let's try clearing everything
dfh.clear()

# one more time with direct
C1.np[:] = np.random.random()
C2.np[:] = np.random.random()
C3.np[:] = np.random.random()
C4.np[:] = np.random.random()

# recompute python side
for i in range(naux):
    Qmo[0][i] = np.dot(C1.np.T, Qpq[i]).dot(C2)
    Qmo[1][i] = np.dot(C3.np.T, Qpq[i]).dot(C2)
    Qmo[2][i] = np.dot(C4.np.T, Qpq[i]).dot(C2)
    Qmo[3][i] = np.dot(C1.np.T, Qpq[i]).dot(C4)
    Qmo[4][i] = np.dot(C3.np.T, Qpq[i]).dot(C1)
    Qmo[5][i] = np.dot(C3.np.T, Qpq[i]).dot(C3)

# set spaces
dfh.add_space("C1", C1)
dfh.add_space("C2", C2)
dfh.add_space("C3", C3)
dfh.add_space("C4", C4)

# add transformations
dfh.add_transformation("Qmo1", "C1", "C2")      # best on left  (Q|bw)
dfh.add_transformation("Qmo2", "C3", "C2")      # best on right (Q|wb)
dfh.add_transformation("Qmo3", "C4", "C2")      # best on right (w|Qb)
dfh.add_transformation("Qmo4", "C1", "C4")      # best on left  (b|Qw)
dfh.add_transformation("Qmo5", "C3", "C1")      # best on right (wb|Q)
dfh.add_transformation("Qmo6", "C3", "C3")      # best on left  (bw|Q)

# invoke transformations
dfh.transform()

# grab transformed integrals
dfh_Qmo[0] = (dfh.get_tensor("Qmo1"))
dfh_Qmo[1] = (dfh.get_tensor("Qmo2"))
dfh_Qmo[2] = (dfh.get_tensor("Qmo3"))
dfh_Qmo[3] = (dfh.get_tensor("Qmo4"))
dfh_Qmo[4] = (dfh.get_tensor("Qmo5"))
dfh_Qmo[5] = (dfh.get_tensor("Qmo6"))

# am i right?
for i in range(6):
    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "DIRECT - on disc" )     #TEST

# switching to core ---------------------------------------------------------------------

#del dfh
#
## redeclare
#dfh = psi4.core.DF_Helper(primary, aux)
#
## tweak options
#dfh.set_method("STORE")
#dfh.set_memory(50000)
#dfh.set_on_core(True)
#
## build
#dfh.initialize()
#
## set spaces
#dfh.add_space("C1", C1)
#dfh.add_space("C2", C2)
#dfh.add_space("C3", C3)
#dfh.add_space("C4", C4)
#
## add transformations
#dfh.add_transformation("Qmo1", "C1", "C2")      # best on left  (Q|bw)
#dfh.add_transformation("Qmo2", "C3", "C2")      # best on right (Q|wb)
#dfh.add_transformation("Qmo3", "C4", "C2")      # best on right (w|Qb)
#dfh.add_transformation("Qmo4", "C1", "C4")      # best on left  (b|Qw)
#dfh.add_transformation("Qmo5", "C3", "C1")      # best on right (wb|Q)
#dfh.add_transformation("Qmo6", "C3", "C3")      # best on left  (bw|Q)
#
## invoke transformations
#dfh.transform()
#
## grab transformed integrals
#dfh_Qmo[0] = (dfh.get_tensor("Qmo1"))
#dfh_Qmo[1] = (dfh.get_tensor("Qmo2"))
#dfh_Qmo[2] = (dfh.get_tensor("Qmo3"))
#dfh_Qmo[3] = (dfh.get_tensor("Qmo4"))
#dfh_Qmo[4] = (dfh.get_tensor("Qmo5"))
#dfh_Qmo[5] = (dfh.get_tensor("Qmo6"))
#
## am i right?
#for i in range(6):
#    print (np.linalg.norm(np.asarray(dfh_Qmo[i])) - np.linalg.norm(Qmo[i]))     #TEST
#    
##    psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, "DIRECT - in core" )     #TEST

