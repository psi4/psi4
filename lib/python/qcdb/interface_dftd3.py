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

"""Module with functions that interface with Grimme's DFTD3 code."""
import os
import re
import sys
import math
import shutil
import socket
import random
import subprocess
#CUimport psi4
#CUimport p4const
#CUimport p4util
try:
    from p4xcpt import *
except ImportError:
    from exceptions import *
from dashparam import *


def run_dftd3(self, func=None, dashlvl=None, dashparam=None, dertype=None):
    """Function to call Grimme's dftd3 program (http://toc.uni-muenster.de/DFTD3/)
    to compute the -D correction of level *dashlvl* using parameters for
    the functional *func*. The dictionary *dashparam* can be used to supply
    a full set of dispersion parameters in the absense of *func* or to supply
    individual overrides in the presence of *func*. Returns energy if *dertype* is 0,
    gradient if *dertype* is 1, else tuple of energy and gradient if *dertype*
    unspecified. The dftd3 executable must be independently compiled and found in
    :envvar:`PATH` or :envvar:`PSIPATH`.

    """
    # Validate arguments
    if self is None:
        self = psi4.get_active_molecule()  # needed for C-side access

    dashlvl = dashlvl.lower()
    dashlvl = dash_alias['-' + dashlvl][1:] if ('-' + dashlvl) in dash_alias.keys() else dashlvl
    if dashlvl not in dashcoeff.keys():
        raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))

    if dertype is None:
        dertype = -1
    elif der0th.match(str(dertype)):
        dertype = 0
    elif der1st.match(str(dertype)):
        dertype = 1
    elif der2nd.match(str(dertype)):
        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))
    else:
        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))

    if func is None:
        if dashparam is None:
            # defunct case
            raise ValidationError("""Parameters for -D correction missing. Provide a func or a dashparam kwarg.""")
        else:
            # case where all param read from dashparam dict (which must have all correct keys)
            func = 'custom'
            dashcoeff[dashlvl][func] = {}
            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
            for key in dashcoeff[dashlvl]['b3lyp'].keys():
                if key in dashparam.keys():
                    dashcoeff[dashlvl][func][key] = dashparam[key]
                else:
                    raise ValidationError("""Parameter %s is missing from dashparam dict %s.""" % (key, dashparam))
    else:
        func = func.lower()
        if func not in dashcoeff[dashlvl].keys():
            raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))
        if dashparam is None:
            # (normal) case where all param taken from dashcoeff above
            pass
        else:
            # case where items in dashparam dict can override param taken from dashcoeff above
            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
            for key in dashcoeff[dashlvl]['b3lyp'].keys():
                if key in dashparam.keys():
                    dashcoeff[dashlvl][func][key] = dashparam[key]

    # Move ~/.dftd3par.<hostname> out of the way so it won't interfere
    defaultfile = os.path.expanduser('~') + '/.dftd3par.' + socket.gethostname()
    defmoved = False
    if os.path.isfile(defaultfile):
        os.rename(defaultfile, defaultfile + '_hide')
        defmoved = True

    # Find environment by merging PSIPATH and PATH environment variables
    lenv = os.environ
    lenv['PATH'] = ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':')]) + ':' + lenv.get('PATH')

    # Find out if running from Psi4 for scratch details and such
    try:
        psi4.version()
    except NameError:
        isP4regime = False
    else:
        isP4regime = True

    # Setup unique scratch directory and move in
    current_directory = os.getcwd()
    if isP4regime:
        psioh = psi4.IOManager.shared_object()
        psio = psi4.IO.shared_object()
        os.chdir(psioh.get_default_path())
        dftd3_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
            '.dftd3.' + str(random.randint(0, 99999))
    else:
        #dftd3_tmpdir = 'dftd3_' + str(random.randint(0, 99999))
        dftd3_tmpdir = os.environ['HOME'] + os.sep + 'dftd3_' + str(random.randint(0, 99999))
    if os.path.exists(dftd3_tmpdir) is False:
        os.mkdir(dftd3_tmpdir)
    os.chdir(dftd3_tmpdir)

    # Write dftd3_parameters file that governs dispersion calc
    paramfile = './dftd3_parameters'
    pfile = open(paramfile, 'w')
    pfile.write(dash_server(func, dashlvl, 'dftd3'))
    pfile.close()

    # Write dftd3_geometry file that supplies geometry to dispersion calc
    geomfile = './dftd3_geometry.xyz'
    gfile = open(geomfile, 'w')
    numAtoms = self.natom()
    geom = self.save_string_xyz()
    reals = []
    for line in geom.splitlines():
        lline = line.split()
        if len(lline) != 4:
            continue
        if lline[0] == 'Gh':
            numAtoms -= 1
        else:
            reals.append(line)

    geomtext = str(numAtoms) + '\n\n'
    for line in reals:
        geomtext += line.strip() + '\n'
    with open(geomfile, 'w') as handle:
        handle.write(geomtext)
    # TODO somehow the variations on save_string_xyz and
    #   whether natom and chgmult does or doesn't get written
    #   have gotten all tangled. I fear this doesn't work
    #   the same btwn libmints and qcdb or for ghosts

    # Call dftd3 program
    command = ['dftd3', geomfile]
    if dertype != 0:
        command.append('-grad')
    try:
        dashout = subprocess.Popen(command, stdout=subprocess.PIPE, env=lenv)
    except OSError:
        raise ValidationError('Program dftd3 not found in path.')
    out, err = dashout.communicate()

    # Parse output (could go further and break into E6, E8, E10 and Cn coeff)
    success = False
    for line in out.splitlines():
        if re.match(' Edisp /kcal,au', line):
            sline = line.split()
            dashd = float(sline[3])
        if re.match(' normal termination of dftd3', line):
            success = True

    if not success:
        raise ValidationError('Program dftd3 did not complete successfully.')

    # Parse grad output
    if dertype != 0:
        derivfile = './dftd3_gradient'
        dfile = open(derivfile, 'r')
        dashdderiv = []
        for line in geom.splitlines():
            lline = line.split()
            if len(lline) != 4:
                continue
            if lline[0] == 'Gh':
                dashdderiv.append([0.0, 0.0, 0.0])
            else:
                dashdderiv.append([float(x.replace('D', 'E')) for x in dfile.readline().split()])
        dfile.close()

        if len(dashdderiv) != self.natom():
            raise ValidationError('Program dftd3 gradient file has %d atoms- %d expected.' % \
                (len(dashdderiv), self.natom()))

    # Prepare results for Psi4
    if isP4regime:
        psi4.set_variable('DISPERSION CORRECTION ENERGY', dashd)
        psi_dashdderiv = psi4.Matrix(self.natom(), 3)
        psi_dashdderiv.set(dashdderiv)

    # Print program output to file if verbose
    verbose = psi4.get_option('SCF', 'PRINT') if isP4regime else 2  # TODO
    if verbose >= 3:

        text = '\n  ==> DFTD3 Output <==\n'
        text += out
        if dertype != 0:
            with open(derivfile, 'r') as handle:
                text += handle.read().replace('D', 'E')
            text += '\n'
        if isP4regime:
            psi4.print_out(text)  # TODO
        else:
            print text

    # Clean up files and remove scratch directory
    os.unlink(paramfile)
    os.unlink(geomfile)
    if dertype != 0:
        os.unlink(derivfile)
    if defmoved is True:
        os.rename(defaultfile + '_hide', defaultfile)

    os.chdir('..')
    try:
        shutil.rmtree(dftd3_tmpdir)
    except OSError as e:
        ValidationError('Unable to remove dftd3 temporary directory %s' % e, file=sys.stderr)
    os.chdir(current_directory)

    # return -D & d(-D)/dx
    if dertype == -1:
        return dashd, dashdderiv
    elif dertype == 0:
        return dashd
    elif dertype == 1:
        return psi_dashdderiv

try:
    # Attach method to libmints psi4.Molecule class
    psi4.Molecule.run_dftd3 = run_dftd3
except (NameError, AttributeError):
    # But don't worry if that doesn't work b/c
    #   it'll get attached to qcdb.Molecule class
    pass
