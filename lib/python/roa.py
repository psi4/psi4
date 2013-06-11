
import copy
import os
import p4util
import psi4

def run_roa(name, **kwargs):

    # Generate input files
    molecule = psi4.get_active_molecule()
    geometry = molecule.save_string_xyz()
    natom = molecule.natom()
    properties = kwargs['properties']
    print(properties)

    coordinates = ['x','y','z']
    step_direction = ['p', 'm']
    displacements = []
    for atom in range(1, natom+1):
        for coord in coordinates:
            for step in step_direction:
                displacements.append('{0}_{1}_{2}'.format(atom, coord, step))
    print(displacements)
    atoms, xyz_original = convert_geometry(geometry)
    for displacement in displacements:
        if not os.path.exists(displacement):
            os.makedirs(displacement)
        xyz = displace_geometry(xyz_original, displacement)
        inputfile = open('{0}/input.dat'.format(displacement), 'w')
        inputfile.write('{0} {1}'.format(molecule.molecular_charge(),
                                         molecule.multiplicity()))
        inputfile.write('\n')
        for atom in range(natom):
            inputfile.write(atoms[atom])
            inputfile.write('   ')
            for i in range(3):
                inputfile.write(repr(xyz[3*atom+i]))
                inputfile.write('  ')
            inputfile.write('\n')
        inputfile.write(p4util.format_options_for_input())
        inputfile.write('\n')
        job_type = "property('{0}', properties='roa_tensor')".format(name)
        inputfile.write(job_type)
        inputfile.close()



    # Check job status
    # Gather results


# Takes a geometry as a string and returns it as a matrix
def convert_geometry(geometry):
    temp = []
    atom_array = []
    coord_matrix = []
    temp = geometry.split()
    temp.pop(0)
    temp.pop(0)
    for n,item in enumerate(temp):
        if not n % 4:
            atom_array.append(item)
        else:
            coord_matrix.append(float(item))

    return atom_array, coord_matrix

def displace_geometry(xyz, displacement):
    delta = 0.1
    new_xyz = copy.deepcopy(xyz)
    (atom, coord, step) = displacement.split('_')
    atom = int(atom)
    coord_dict = {'x': 0, 'y': 1, 'z': 2}
    if step == 'p':
        new_xyz[3*(atom-1)+coord_dict[coord]] = xyz[3*(atom-1)+coord_dict[coord]] + delta
    else:
        new_xyz[3*(atom-1)+coord_dict[coord]] = xyz[3*(atom-1)+coord_dict[coord]] - delta

    return new_xyz
