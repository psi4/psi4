import collections
import shelve
import copy
import os
import p4util
import psi4
from p4const import *

def run_roa(name, **kwargs):

    # Get list of omega values -> Make sure we only have one wavelength
    # Catch this now before any real work gets done
    omega = psi4.get_option('CCRESPONSE','OMEGA')
    if len(omega) > 2:
        raise Exception('ROA scattering can only be performed for one wavelength.')
    else:
        pass

    psi4.print_out('Running ROA computation. Subdirectories for each '
              'required displaced geometry have been created.\n\n')

    ### Initialize database
    db = shelve.open('database',writeback=True)
    if not db.has_key('inputs_generated'):
        initialize_database(db)

    ### Generate input files
    if not db['inputs_generated']:
        generate_inputs(name, db)
        db['inputs_generated'] = True
    
    ### Check job status
    if db['inputs_generated'] and not db['jobs_complete']:
        print('Checking status')
        roa_stat(db)
        for job,status in db['job_status'].items():
            print("{} --> {}".format(job,status))

    ### Compute ROA Scattering
    if db['jobs_complete']:
#   SAVE this for when multiple wavelengths works
#        # Get list of omega values
#        omega = psi4.get_option('CCRESPONSE','OMEGA')
#        if len(omega) > 1:
#            units = copy.copy(omega[-1])
#            omega.pop()
#        else:
#            units = 'atomic'
#        wavelength = copy.copy(omega[0])
#        # Set up units for scatter.cc
#        if units == 'NM':
#            wavelength = (psi_c * psi_h * 1*(10**-9))/(wavelength * psi_hartree2J)
#        if units == 'HZ':
#            wavelength = wavelength * psi_h / psi_hartree2J
#        if units == 'EV':
#            wavelength = wavelength / psi_hartree2ev
#        if units == 'atomic':
#            pass
        # Initialize tensor lists
        dip_polar_list = []
        opt_rot_list = []
        dip_quad_polar_list = []
        gauge_list = []
        make_gauge_list(gauge_list)
        # Gather data
        synthesize_dipole_polar(db,dip_polar_list)
        synthesize_opt_rot(db,opt_rot_list)
        synthesize_dip_quad_polar(db,dip_quad_polar_list)
        # Compute Scattering
	    # Run new function (src/bin/ccresponse/scatter.cc)
        print('Running scatter function')
        step = psi4.get_local_option('FINDIF','DISP_SIZE')
        for gauge in opt_rot_list:
            g_idx = opt_rot_list.index(gauge)
            print('\n\n----------------------------------------------------------------------')
            print('\t%%%%%%%%%% {} %%%%%%%%%%'.format(gauge_list[g_idx]))
            print('----------------------------------------------------------------------\n\n')
            psi4.print_out('\n\n----------------------------------------------------------------------\n')
            psi4.print_out('\t%%%%%%%%%% {} %%%%%%%%%%\n'.format(gauge_list[g_idx]))
            psi4.print_out('----------------------------------------------------------------------\n\n')
            psi4.scatter(step, dip_polar_list, gauge, dip_quad_polar_list)

    db.close()

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
        molecule.fix_orientation(True)
        molecule.fix_com(True)

        # Write input file
        inputfile = open('{0}/input.dat'.format(entry), 'w')
        inputfile.write("# This is a psi4 input file auto-generated for "
                        "computing Raman Optical Activity.\n\n")
		#inputfile.write(basic_molecule_for_input(molecule))
        inputfile.write("{}{}{}"
						.format(mol_open,molecule.create_psi4_string_from_molecule(),mol_close))
                        #.format(mol_open,molecule.save_string_xyz(),mol_close))
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
    length = []
    velocity = []
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            mygauge = psi4.get_option('CCRESPONSE', 'GAUGE')
            if mygauge == 'LENGTH':
                length.append(grab_psi4_matrix(outfile, 'Optical Rotation Tensor (Length Gauge)', 3))
            elif mygauge == 'VELOCITY':
                velocity.append(grab_psi4_matrix(outfile, 'Optical Rotation Tensor (Modified Velocity Gauge)', 3))
            elif mygauge == 'BOTH':
                length.append(grab_psi4_matrix(outfile, 'Optical Rotation Tensor (Length Gauge)', 3))
                velocity.append(grab_psi4_matrix(outfile, 'Optical Rotation Tensor (Modified Velocity Gauge)', 3))
            else:
                print("There is no optical rotation tensor - something is wrong.")	
    if length:
        opt_rot_list.append(length)
    if velocity:
        opt_rot_list.append(velocity)
#    if length and not velocity:
#        opt_rot_list.append(length)
#        gauge_list.append('Length Gauge Results')
#    if velocity and not length:
#        opt_rot_list.append(velocity)
#        gauge_list.append('Modified Velocity Gauge Results')
#    if length and velocity:
#        gauge_list.append('Length Gauge Results')
#        gauge_list.append('Modified Velocity Gauge Results')
#        temp_list = []
#        temp_list.append(length)
#        temp_list.append(velocity)
#        opt_rot_list.append(temp_list)


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

# THIS IS UNECESSARY SINCE ONLY PURE VELOCITY GAUGE IS EVALUATED AT ZERO WAVELENGTH
# AND WHEN THE FUNCTION IS CALLED, IT IS PASSED MVG AS THE MATRIX NAME
#def grab_psi4_matrix(outfile, matrix_name, row_tot):
#    collect_matrix = False
#    collect_omega  = False
#    n_rows = 0
#    n_tries = 0
#    matrix_data = []
#    omega_zero = 'omega = 0.00'
#    #mygauge = psi4.get_option('CCRESPONSE','GAUGE')
#    for line in outfile:
#        # Quick gauge check

#            collect_omega = True
#        # Are we at the right kind of tensor?
#        if matrix_name in line:
#            collect_matrix = True
#        # Is this the zero wavelength data for MVG?
#        if omega_zero in line:
#            collect_omega  = True
#            collect_matrix = False
#        if collect_matrix and collect_omega and (n_rows < row_tot):
#            try:
#                n_tries += 1
#                if n_tries > (row_tot + 13):
#                    raise Exception('{} matrix was unreadable. Scanned {} '
#                                    'lines'.format(matrix_name, n_tries))
#                else:
#                    (index, x, y, z) = line.split()
#                    matrix_data.append(float(x))
#                    matrix_data.append(float(y))
#                    matrix_data.append(float(z))
#                    n_rows += 1
#            except:
#                pass
#
#        if (n_rows == row_tot) and (len(matrix_data) != 3*row_tot):
#            raise Exception('Collecting matrix data failed!')
#
#        if len(matrix_data) == 3*row_tot:
#            return matrix_data

def make_gauge_list(gauge_list):
    mygauge = psi4.get_option('CCRESPONSE', 'GAUGE')
    if mygauge == 'LENGTH':
        gauge_list.append('Length Gauge Results')
    elif mygauge == 'VELOCITY':
        gauge_list.append('Modified Velocity Gauge Results')
    elif mygauge == 'BOTH':
        gauge_list.append('Length Gauge Results')
        gauge_list.append('Modified Velocity Gauge Results')
    else:
        print("There is no optical rotation tensor - something is wrong.")
       
    



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


