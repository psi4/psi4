import copy
import os
import p4util
import psi4
from p4const import *

def run_roa(name, **kwargs):

    # Generate input files
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()

    # Get list of displacements
    displacement_geoms = psi4.atomic_displacements()

    # Until we append the original geometry
    if not (6*natom) == len(displacement_geoms):
        raise Exception('The number of atoms and displacements should be the '
                        'same!')

    geometry = molecule.save_string_xyz()
    properties = kwargs['properties']

    for geom in displacement_geoms:
        molecule.set_geometry(geom)
        # save_string_xyz converts geometry to angstroms
#        print(molecule.save_string_xyz())

    coordinates =    ['x','y','z']
    step_direction = ['p','m']

    # List of displacements for iterating
    # Including original for convenience and to obtain other properties
    displacement_names = []
    for atom in range(1, natom+1):
        for coord in coordinates:
            for step in step_direction:
                displacement_names.append('{0}_{1}_{2}'.format(atom, coord, step))
    print(displacement_names)

    for n,entry in enumerate(displacement_names):
        if not os.path.exists(entry):
            os.makedirs(entry)
        # Set up 
        mol_open = 'molecule ' + molecule.name() + '_' + entry + ' {\n'
        mol_close = '}'
        molecule.set_geometry(displacement_geoms[n])

        # Writing
        inputfile = open('{0}/input.dat'.format(entry), 'w')
        inputfile.write("# This is a psi4 input file auto-generated for "
                        "computing Raman Optical Activity.\n\n")
        inputfile.write("{}{}{}"
                        .format(mol_open,molecule.save_string_xyz(),mol_close))
        inputfile.write('\n')
        inputfile.write(p4util.format_options_for_input())
        inputfile.write('\n')
        inputfile.write("property('{0}', properties=['roa_tensor'])".format(name))
        inputfile.close()

    # Check job status
    # Gather results


