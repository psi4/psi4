import collections
import shelve
import copy
import os
import p4util
import psi4
from p4const import *

def run_roa(name, **kwargs):

    psi4.print_out('Running ROA computation. Subdirectories for each '
              'required displaced geometry have been created.')

    # Initialize database
    db = shelve.open('database',writeback=True)
    if not db.has_key('inputs_generated'):
        initialize_database(db)

    # Generate input files
    if not db['inputs_generated']:
        generate_inputs(name, db)
        db['inputs_generated'] = True
    
    # Check job status
    if db['inputs_generated'] and not db['jobs_complete']:
        print('Checking status')
        roa_stat(db)

    # Gather results
    if db['jobs_complete']:
        dip_polar_list = []
        opt_rot_list = []
        dip_quad_polar_list = []
        synthesize_dipole_polar(db,dip_polar_list)
        synthesize_opt_rot(db,opt_rot_list)
        synthesize_dip_quad_polar(db,dip_quad_polar_list)
        print(dip_quad_polar_list)
    db.close()
    # Run new function (scatter.cc)
    psi4.scatter()

def initialize_database(database):
    database['inputs_generated'] = False
    database['jobs_complete']    = False
    database['roa_computed']     = False
    database['job_status'] = collections.OrderedDict()
    # Populate job_status 
    molecule = psi4.get_active_molecule()
    natom    = molecule.natom()
    coordinates    = ['x','y','z']
    step_direction = ['p','m']

    for atom in range(1, natom+1): 
        for coord in coordinates:
            for step in step_direction:
                job_name = '{}_{}_{}'.format(atom,coord,step)
                database['job_status'].update({job_name: 'not_started'})

def generate_inputs(name,db):
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()

    # Get list of displacements
    displacement_geoms = psi4.atomic_displacements()

    # Sanity check!
    # Until we append the original geometry
    if not (6*natom) == len(displacement_geoms):
        raise Exception('The number displacements should be 6 times the number'
                        'of atoms!')

    # List of displacements for iterating
    displacement_names = db['job_status'].keys()

    for n,entry in enumerate(displacement_names):
        if not os.path.exists(entry):
            os.makedirs(entry)
        # Set up 
        mol_open = 'molecule ' + molecule.name() + '_' + entry + ' {\n'
        mol_close = '}'
        molecule.set_geometry(displacement_geoms[n])

        # Write input file
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

def roa_stat(db):
    # Number of jobs that are finished
    n_finished = 0
    for job,status in db['job_status'].items():
        if status == 'finished':
            n_finished += 1
        elif status in ('not_started','running'):
            try: 
                with open('{}/output.dat'.format(job)) as outfile:
                    outfile.seek(-150,2)
                    for line in outfile:
                        if 'PSI4 exiting successfully' in line:
                            db['job_status'][job] = 'finished'
                            n_finished += 1
                            break
                    else:
                        db['job_status'][job] = 'running'
            except:
                pass
    # Are all jobs finished?
    if n_finished == len(db['job_status'].keys()):
        db['jobs_complete'] = True

def synthesize_dipole_polar(db,dip_polar_list):
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            dip_polar_list.append(grab_psi4_matrix(outfile, 'Dipole '
                                                   'Polarizability', 3))

def synthesize_opt_rot(db, opt_rot_list):
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            opt_rot_list.append(grab_psi4_matrix(outfile, 'Optical Rotation Tensor (Modified Velocity Gauge)', 3))

def synthesize_dip_quad_polar(db, dip_quad_polar_list):
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            dip_quad_polar_list.append(grab_psi4_matrix(outfile,
                                       'Electric-Dipole/Quadrupole '
                                                        'Polarizability', 9))

def grab_psi4_matrix(outfile, matrix_name, row_tot):
    collect_matrix = False
    n_rows = 0
    n_tries = 0
    matrix_data = []
    for line in outfile:
        if matrix_name in line:
            collect_matrix = True
        if collect_matrix and (n_rows < row_tot):
            try:
                n_tries += 1
                if n_tries > (row_tot + 13):
                    raise Exception('{} matrix was unreadable. Scanned {} '
                                    'lines'.format(matrix_name, n_tries))
                else:
                    (index, x, y, z) = line.split()
                    print(x)
                    matrix_data.append(float(x))
                    matrix_data.append(float(y))
                    matrix_data.append(float(z))
                    n_rows += 1
            except:
                pass

        if (n_rows == row_tot) and (len(matrix_data) != 3*row_tot):
            raise Exception('Collecting matrix data failed!')

        if len(matrix_data) == 3*row_tot:
            return matrix_data




################################
###                          ###
###    DATABASE STRUCTURE    ###
###                          ###
################################

#Dict of dicts
#inputs_generated (boolean)
#job_status: job_list: status
#jobs_complete
#roa_computed
#
#data ? 
#data: dipole_polarizability
#    : optical_rotation
#    : dipole_quadrupole_polarizability


