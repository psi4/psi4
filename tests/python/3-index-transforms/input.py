#! examine JK packing forms

import psi4
import numpy as np
from collections import OrderedDict

psi4.set_output_file("output.dat", False)

mol = psi4.geometry("""
0 2
H 0.0  2.0 0.0    
H 0.0  4.0 0.0
H 0.0  6.0 0.0
H 0.0  8.0 0.0
H 0.0 10.0 0.0
""")

primary = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVTZ")
aux = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVTZ-jkfit")

# don't sieve the integrals, because we're checking to high precision
psi4.set_options({'ints_tolerance': 0.0,})

# setup ~
nbf = primary.nbf()
naux = aux.nbf()
mem = 200000
mem_bump = naux * naux
ntransforms = 6

# form metric
mints = psi4.core.MintsHelper(primary)
zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
Jmetric = np.squeeze(mints.ao_eri(aux, zero_bas, aux, zero_bas))

# form inverse metric
Jmetric_inv = mints.ao_eri(aux, zero_bas, aux, zero_bas)
Jmetric_inv.power(-0.5, 1.e-12)
Jmetric_inv = np.squeeze(Jmetric_inv)

# form Qpq
Qpq = np.squeeze(mints.ao_eri(aux, zero_bas, primary, primary))
Qpq = Jmetric_inv.dot(Qpq.reshape(naux, -1))
Qpq = Qpq.reshape(naux, nbf, -1) 

# construct spaces
names = ['C1', 'C2', 'C3', 'C4']
sizes = [_ * 5 for _ in range(1, 5)]
spaces = {names[ind]: psi4.core.Matrix.from_array(np.random.rand(nbf, _)) for ind, _ in enumerate(sizes)}

# set transformed integrals
Qmo = []
space_pairs = [[0, 1], [2, 1], [3, 1], [0, 3], [2, 0], [2, 2]]
for i in space_pairs:
    Qmo.append(np.zeros((naux, sizes[i[0]], sizes[i[1]])))

# transform
for ind, i in enumerate(space_pairs): 
    for Q in range(naux):
        C1 = spaces[names[space_pairs[ind][0]]]
        C2 = spaces[names[space_pairs[ind][1]]]
        Qmo[ind][Q] = np.dot(C1.np.T, Qpq[Q]).dot(C2)
    
# get other forms (pQq and pqQ)
Qmo_pQq = []    
Qmo_pqQ = []    
for i in range(ntransforms):
    Qmo_pQq.append(np.swapaxes(Qmo[i], 0, 1))    
    Qmo_pqQ.append(np.swapaxes(np.swapaxes(Qmo[i], 0, 2), 0, 1))    

# setup ~
methods = ['STORE', 'DIRECT', 'DIRECT_iaQ']
forms = ['Qpq', 'pQq', 'pqQ']
transformation_names = ['Qmo1', 'Qmo2', 'Qmo3', 'Qmo4', 'Qmo5', 'Qmo6']
transformations = OrderedDict({})
for ind, i in enumerate(transformation_names):
    transformations[i] = [names[space_pairs[ind][0]], names[space_pairs[ind][1]]] 

# somewhat exhuastive search on all options
for method in methods:
    for form in forms:
        if(form != 'pqQ' and method == 'DIRECT_iaQ'): continue
        for AO_core in [False, True]:
            for MO_core in [False, True]:
                for hold_met in [False, True]:
                            
                    # get object
                    dfh = psi4.core.DFHelper(primary, aux)
                    
                    # set test options
                    dfh.set_method(method)
                    memory = mem_bump if hold_met else 0
                    memory += 10*mem if AO_core else mem
                    dfh.set_memory(memory)
                    dfh.set_AO_core(AO_core)
                    dfh.set_MO_core(MO_core)
                    dfh.hold_met(hold_met)

                    # build
                    dfh.initialize()
                    dfh.print_header()                   
 
                    # add spaces
                    for i in spaces:
                        dfh.add_space(i, spaces[i])

                    # add transformations
                    for i in transformations:
                        j = transformations[i]
                        dfh.add_transformation(i, j[0], j[1], form) 

                    # invoke transformations
                    dfh.transform()

                    # grab transformed integrals
                    dfh_Qmo = []
                    if(form == 'pqQ'):    
                        for ind, i in enumerate(transformations):
                            j = space_pairs[ind]
                            dfh_Qmo.append(np.zeros((sizes[j[0]], sizes[j[1]], naux)))
                            for k in range(sizes[j[0]]):
                                dfh_Qmo[ind][k,:,:] = np.asarray(dfh.get_tensor(i, [k, k+1], [0, sizes[j[1]]], [0, naux]))
                    else:
                        for ind, i in enumerate(transformations):
                            dfh_Qmo.append(np.asarray(dfh.get_tensor(i)))

                    test_string = 'Alg: ' + method + ' + ' + form + ' core (AOs, MOs, met): [' 
                    test_string += str(AO_core) + ', ' + str(MO_core) + ', ' + str(hold_met) +  ']' 

                    print(test_string)
                    # am i right?
                    for i in range(ntransforms):
                        if(form == 'pqQ'):    
                            psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo_pqQ[i], 9, test_string)
                        elif(form == 'pQq'):
                            psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo_pQq[i], 9, test_string)
                        else:
                            psi4.compare_arrays(np.asarray(dfh_Qmo[i]), Qmo[i], 9, test_string)

                    del dfh

# TODO:
# test tensor slicing grabs
# test pQq and pqQ builds for store and direct0


