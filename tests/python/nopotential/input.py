#! Python equivalent of extern5 test:
#! External potential sanity check with 0 charge far away
#! Checks if all units behave the same and energy is same as no
#! potential
import numpy as np
import psi4.core
import psi4

b2a=0.529177249
# Coordinates added in angstrom
coords = np.array([[ -0.778803000000 , 0.000000000000,  1.132683000000],
 [ -0.666682000000,  0.764099000000,  1.706291000000],
 [ -0.666682000000,  -0.764099000000 , 1.706290000000]])
elements = ["O","H","H"]
molecule_ang  = psi4.core.Molecule.from_arrays(geom=coords,     elem=elements, fix_symmetry="c1", fix_com=True, fix_orientation=True)
molecule_bohr = psi4.core.Molecule.from_arrays(geom=coords/b2a, elem=elements, fix_symmetry="c1", fix_com=True, fix_orientation=True, units="Bohr")

external_potentials = [[0.00, np.array([10.0,10.0,10.0]) / b2a]]

psi4.set_options( {
    "scf_type": "df",
    "d_convergence": 12,
    "basis": "STO-3G",
    "print": 0,
    "debug": 0,
})


ene_bohr_charges = psi4.energy('scf', molecule=molecule_bohr, external_potentials=external_potentials)
ene_bohr_pure = psi4.energy('scf', molecule=molecule_bohr)
psi4.compare_values(ene_bohr_charges, ene_bohr_pure, 6, "Bohr geometry, charges vs no charges energy equality")

ene_ang_pure = psi4.energy('scf', molecule=molecule_ang)
psi4.compare_values(ene_ang_pure, ene_bohr_pure, 6, "No charges, Bohr vs Angstrom geometry energy equality")

ene_ang_charges = psi4.energy('scf', molecule=molecule_ang, external_potentials=external_potentials)
psi4.compare_values(ene_ang_charges, ene_ang_pure, 6, "Angstrom geometry, charges vs no charges energy equality")
