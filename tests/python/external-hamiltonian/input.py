#! Testing passing in an external hamiltonian
#! The external hamiltonian is taken from a simulation with an explicit
#! potential and passed into another one without an external potential
import numpy as np
import psi4.core
import psi4
coords = np.array([[ -0.778803000000 , 0.000000000000,  1.132683000000],
 [ -0.666682000000,  0.764099000000,  1.706291000000],
 [ -0.666682000000,  -0.764099000000 , 1.706290000000]]
                  )
elements = ["O","H","H"]
molecule = psi4.core.Molecule.from_arrays(geom=coords, elem=elements, fix_symmetry="c1", fix_com=True, fix_orientation=True)

b2a=0.529177249
external_potentials = [
    [-0.834, np.array([1.649232019048,0.0,-2.356023604706]) / b2a],
    [ 0.417, np.array([0.544757019107,0.0,-3.799961446760]) / b2a],
    [ 0.417, np.array([0.544757019107,0.0,-0.912085762652]) / b2a]
]

psi4.set_options( {
    "scf_type": "df",
    "d_convergence": 12,
    "basis": "STO-3G",
    "print": 1,
})

# This hamiltonian can be read out by setting print to 4 and
# then checking the output of form_H
equivalent_H = np.array([[-0.0299706 , -0.00709416,  0.00051994, -0.00022709,  0.        ,
        -0.00157738, -0.00157738],
       [-0.00709416, -0.0299706 ,  0.00656341, -0.00286662,  0.        ,
        -0.01201169, -0.01201171],
       [ 0.00051994,  0.00656341, -0.03114568,  0.00116943,  0.        ,
        -0.0032155 , -0.00321549],
       [-0.00022709, -0.00286662,  0.00116943, -0.02972556,  0.        ,
        -0.00210291, -0.00210291],
       [ 0.        ,  0.        ,  0.        ,  0.        , -0.02904055,
        -0.00754691,  0.00754692],
       [-0.00157738, -0.01201169, -0.0032155 , -0.00210291, -0.00754691,
        -0.01957676, -0.00515716],
       [-0.00157738, -0.01201171, -0.00321549, -0.00210291,  0.00754692,
        -0.00515716, -0.01957677]])

# This value is taken from the charged simulation without ext H. It is the expected 
# difference between the two simulations:
additional_nuclear_repulsion = 0.278918444763445

ene_from_charges = psi4.energy('scf', molecule=molecule, external_potentials=external_potentials)
ene_from_charges_without_nuclear_repulsion = ene_from_charges - additional_nuclear_repulsion

ene_from_external_hamiltonian = psi4.energy('scf', molecule=molecule, external_hamiltonian=equivalent_H)
psi4.compare_values(ene_from_charges_without_nuclear_repulsion, ene_from_external_hamiltonian, 6, "Checking identity of energy from explicit external field and additional hamiltonian contribution.")

