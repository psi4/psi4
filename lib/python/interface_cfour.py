#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

from __future__ import print_function
"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function.

"""
import shutil
import os
import subprocess
import re
import psi4
import p4const
import p4util
import qcdb
from p4regex import *
#from extend_Molecule import *
from molutil import *
from functional import *
# never import driver, wrappers, or aliases into this file

# ATTN NEW ADDITIONS!
# consult http://sirius.chem.vt.edu/psi4manual/master/proc_py.html



def run_cfour(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Stanton and Gauss's CFOUR code.

    """
    # The parse_arbitrary_order method provides us the following information
    # We require that level be provided. level is a dictionary
    # of settings to be passed to psi4.cfour
    #if not('level' in kwargs):
    #    raise ValidationError('level parameter was not provided.')

    #level = kwargs['level']

    # Fullname is the string we need to search for in iface
    #fullname = level['fullname']

    # User can provide 'keep' to the method.
    # When provided, do not delete the CFOUR scratch directory.
    keep = False
    if 'keep' in kwargs:
        keep = kwargs['keep']

    # Save current directory location
    current_directory = os.getcwd()

    # Find environment by merging PSIPATH and PATH environment variables
    lenv = os.environ
    lenv['PATH'] = ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':')]) + ':' + lenv.get('PATH')

    # Need to move to the scratch directory, perferrably into a separate directory in that location
    psi_io = psi4.IOManager.shared_object()
    os.chdir(psi_io.get_default_path())

    # Make new directory specifically for cfour 
    cfour_tmpdir = 'cfour_' + str(os.getpid())
    if 'path' in kwargs:
        cfour_tmpdir = kwargs['path']

    # Check to see if directory already exists, if not, create.
    if os.path.exists(cfour_tmpdir) == False:
        os.mkdir(cfour_tmpdir)

    # Move into the new directory
    os.chdir(cfour_tmpdir)

    # Load the GENBAS file
    genbas_path = qcdb.search_file('GENBAS', lenv['PATH'])
    if genbas_path:
        shutil.copy2(genbas_path, psi_io.get_default_path() + cfour_tmpdir)
        psi4.print_out("\n  GENBAS loaded from %s\n" % (genbas_path))
        psi4.print_out("  CFOUR to be run from %s\n" % (psi_io.get_default_path() + cfour_tmpdir))
    else:
        psi4.print_out('\nGENBAS file for CFOUR interface not found\n\n')
        psi4.print_out('Search path that was tried:\n')
        temp = lenv['PATH']
        psi4.print_out(temp.replace(':', ', '))
        raise ValidationError("GENBAS file loading problem")

    # Generate integrals and input file (dumps files to the current directory)
    #psi4.cfour_generate_input(level)
    with open("ZMAT", "w") as cfour_infile:
        cfour_infile.write(psi4.get_global_option('CFOUR_LITERAL'))

    # Load the ZMAT file
    # and dump a copy into the outfile
    psi4.print_out('\n====== Begin ZMAT input for CFOUR ======\n')
    psi4.print_out(open('ZMAT', 'r').read())
    psi4.print_out('======= End ZMAT input for CFOUR =======\n\n')
    print('\n====== Begin ZMAT input for CFOUR ======\n', open('ZMAT', 'r').read(), '======= End ZMAT input for CFOUR =======\n\n')

    # Close output file and reopen
    psi4.close_outfile()
    p4out = open(current_directory + '/' + psi4.outfile_name(), 'a')

    # Modify the environment:
    #    PGI Fortan prints warning to screen if STOP is used
    #os.environ['NO_STOP_MESSAGE'] = '1'

    # Obtain user's OMP_NUM_THREADS so that we don't blow it away.
    omp_num_threads_found = 'OMP_NUM_THREADS' in os.environ
    if omp_num_threads_found == True:
        omp_num_threads_user = os.environ['OMP_NUM_THREADS']

    # If the user provided CFOUR_OMP_NUM_THREADS set the environ to it
    if psi4.has_option_changed('CFOUR', 'CFOUR_OMP_NUM_THREADS') == True:
        os.environ['OMP_NUM_THREADS'] = str(psi4.get_option('CFOUR', 'CFOUR_OMP_NUM_THREADS'))

    # Call xcfour, directing all screen output to the output file
    try:
        retcode = subprocess.Popen(['xcfour'], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        sys.stderr.write('Program xcfour not found in path or execution failed: %s\n' % (e.strerror))
        p4out.write('Program xcfour not found in path or execution failed: %s\n' % (e.strerror))
        sys.exit(1)

    c4out = ""
    while retcode.returncode == None:
        data = retcode.stdout.read(1)
        if psi4.outfile_name() == 'stdout':
            sys.stdout.write(data)
        else:
            p4out.write(data)
            p4out.flush()
        c4out += data
        retcode.poll()

    # Restore the OMP_NUM_THREADS that the user set.
    if omp_num_threads_found == True:
        if psi4.has_option_changed('CFOUR', 'CFOUR_OMP_NUM_THREADS') == True:
            os.environ['OMP_NUM_THREADS'] = omp_num_threads_user

    import qcprograms
    psivars = qcprograms.cfour.cfour_harvest(c4out)
    for key in psivars.keys():
        psi4.set_variable(key.upper(), psivars[key])


    # Scan iface file and grab the file energy.
    #e = 0.0
    #for line in file('iface'):
    #    fields = line.split()
    #    m = fields[1]
    #    try:
    #        e = float(fields[5])
    #        if m == "MP(2)":
    #            m = "MP2"
    #        psi4.set_variable(m + ' TOTAL ENERGY', e)
    #        psi4.set_variable(m + ' CORRELATION ENERGY', e - vscf)
    #    except ValueError:
    #        continue

    # The last 'e' in iface is the one the user requested.
    #psi4.set_variable('CURRENT ENERGY', e)
    #psi4.set_variable('CURRENT CORRELATION ENERGY', e - vscf)

    # Load the iface file
    #iface = open('iface', 'r')
    #iface_contents = iface.read()

    # Delete cfour tempdir
    os.chdir('..')
    try:
        # Delete unless we're told not to
        if (keep == False and not('path' in kwargs)):
            shutil.rmtree(cfour_tmpdir)
    except OSError as e:
        print('Unable to remove CFOUR temporary directory %s' % e, file=sys.stderr)
        exit(1)

    # Revert to previous current directory location
    os.chdir(current_directory)

    # Reopen output file
    psi4.reopen_outfile()
    psi4.print_variables()

    # If we're told to keep the files or the user provided a path, do nothing.
    if (yes.match(str(keep)) or ('path' in kwargs)):
        psi4.print_out('\nCFOUR scratch files have been kept.\n')
        psi4.print_out('They can be found in ' + psi_io.get_default_path() + cfour_tmpdir)

    # Dump iface contents to output
    #psi4.print_out('\n')
    #p4util.banner('Full results from MRCC')
    #psi4.print_out('\n')
    #psi4.print_out(iface_contents)

    #return e

def cfour_list():
    val = []
    val.append('cfour')
    val.append('c4-ccsd')
    return val

