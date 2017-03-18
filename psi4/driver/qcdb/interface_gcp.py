#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

"""Module with functions that interface with Grimme's GCP code."""
from __future__ import absolute_import, print_function
import os
import re
import uuid
import socket
import subprocess
try:
    from psi4.driver.p4util.exceptions import *
    from psi4 import core
    isP4regime = True
except ImportError:
    from .exceptions import *
    isP4regime = False
from .p4regex import *
from .molecule import Molecule


def run_gcp(self, func=None, dertype=None, verbose=False):  # dashlvl=None, dashparam=None
    """Function to call Grimme's dftd3 program (http://toc.uni-muenster.de/DFTD3/)
    to compute the -D correction of level *dashlvl* using parameters for
    the functional *func*. The dictionary *dashparam* can be used to supply
    a full set of dispersion parameters in the absense of *func* or to supply
    individual overrides in the presence of *func*. Returns energy if *dertype* is 0,
    gradient if *dertype* is 1, else tuple of energy and gradient if *dertype*
    unspecified. The dftd3 executable must be independently compiled and found in
    :envvar:`PATH` or :envvar:`PSIPATH`.
    *self* may be either a qcdb.Molecule (sensibly) or a psi4.Molecule
    (works b/c psi4.Molecule has been extended by this method py-side and
    only public interface fns used) or a string that can be instantiated
    into a qcdb.Molecule.

    """
    # Create (if necessary) and update qcdb.Molecule
    if isinstance(self, Molecule):
        # called on a qcdb.Molecule
        pass
    elif isinstance(self, core.Molecule):
        # called on a python export of a psi4.core.Molecule (py-side through Psi4's driver)
        self.create_psi4_string_from_molecule()
    elif isinstance(self, basestring):
        # called on a string representation of a psi4.Molecule (c-side through psi4.Dispersion)
        self = Molecule(self)
    else:
        raise ValidationError("""Argument mol must be psi4string or qcdb.Molecule""")
    self.update_geometry()

#    # Validate arguments
#    dashlvl = dashlvl.lower()
#    dashlvl = dash_alias['-' + dashlvl][1:] if ('-' + dashlvl) in dash_alias.keys() else dashlvl
#    if dashlvl not in dashcoeff.keys():
#        raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))

    if dertype is None:
        dertype = -1
    elif der0th.match(str(dertype)):
        dertype = 0
    elif der1st.match(str(dertype)):
        dertype = 1
#    elif der2nd.match(str(dertype)):
#        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))
    else:
        raise ValidationError('Requested derivative level \'dertype\' %s not valid for run_dftd3.' % (dertype))

#    if func is None:
#        if dashparam is None:
#            # defunct case
#            raise ValidationError("""Parameters for -D correction missing. Provide a func or a dashparam kwarg.""")
#        else:
#            # case where all param read from dashparam dict (which must have all correct keys)
#            func = 'custom'
#            dashcoeff[dashlvl][func] = {}
#            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
#            for key in dashcoeff[dashlvl]['b3lyp'].keys():
#                if key in dashparam.keys():
#                    dashcoeff[dashlvl][func][key] = dashparam[key]
#                else:
#                    raise ValidationError("""Parameter %s is missing from dashparam dict %s.""" % (key, dashparam))
#    else:
#        func = func.lower()
#        if func not in dashcoeff[dashlvl].keys():
#            raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))
#        if dashparam is None:
#            # (normal) case where all param taken from dashcoeff above
#            pass
#        else:
#            # case where items in dashparam dict can override param taken from dashcoeff above
#            dashparam = dict((k.lower(), v) for k, v in dashparam.iteritems())
#            for key in dashcoeff[dashlvl]['b3lyp'].keys():
#                if key in dashparam.keys():
#                    dashcoeff[dashlvl][func][key] = dashparam[key]

    # TODO temp until figure out paramfile
    allowed_funcs = ['HF/MINIS', 'DFT/MINIS', 'HF/MINIX', 'DFT/MINIX',
        'HF/SV', 'DFT/SV', 'HF/def2-SV(P)', 'DFT/def2-SV(P)', 'HF/def2-SVP',
        'DFT/def2-SVP', 'HF/DZP', 'DFT/DZP', 'HF/def-TZVP', 'DFT/def-TZVP',
        'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/631Gd', 'DFT/631Gd',
        'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/cc-pVDZ', 'DFT/cc-pVDZ',
        'HF/aug-cc-pVDZ', 'DFT/aug-cc-pVDZ', 'DFT/SV(P/h,c)', 'DFT/LANL',
        'DFT/pobTZVP', 'TPSS/def2-SVP', 'PW6B95/def2-SVP',
        # specials
        'hf3c', 'pbeh3c']
    allowed_funcs = [f.lower() for f in allowed_funcs]
    if func.lower() not in allowed_funcs:
        raise Dftd3Error("""bad gCP func: %s. need one of: %r""" % (func, allowed_funcs))

    # Move ~/.dftd3par.<hostname> out of the way so it won't interfere
    defaultfile = os.path.expanduser('~') + '/.dftd3par.' + socket.gethostname()
    defmoved = False
    if os.path.isfile(defaultfile):
        os.rename(defaultfile, defaultfile + '_hide')
        defmoved = True

    # Find environment by merging PSIPATH and PATH environment variables
    lenv = {
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
                ':' + os.environ.get('PATH'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    #   Filter out None values as subprocess will fault on them
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # Find out if running from Psi4 for scratch details and such
    try:
        import psi4
    except ImportError as err:
        isP4regime = False
    else:
        isP4regime = True

    # Setup unique scratch directory and move in
    current_directory = os.getcwd()
    if isP4regime:
        psioh = core.IOManager.shared_object()
        psio = core.IO.shared_object()
        os.chdir(psioh.get_default_path())
        gcp_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
            '.gcp.' + str(uuid.uuid4())[:8]
    else:
        gcp_tmpdir = os.environ['HOME'] + os.sep + 'gcp_' + str(uuid.uuid4())[:8]
    if os.path.exists(gcp_tmpdir) is False:
        os.mkdir(gcp_tmpdir)
    os.chdir(gcp_tmpdir)

    # Write gcp_parameters file that governs cp correction
#    paramcontents = gcp_server(func, dashlvl, 'dftd3')
#    paramfile1 = 'dftd3_parameters'  # older patched name
#    with open(paramfile1, 'w') as handle:
#        handle.write(paramcontents)
#    paramfile2 = '.gcppar'
#    with open(paramfile2, 'w') as handle:
#        handle.write(paramcontents)

###Two kinds of parameter files can be read in: A short and an extended version. Both are read from
###$HOME/.gcppar.$HOSTNAME by default. If the option -local is specified the file is read in from
###the current working directory: .gcppar
###The short version reads in: basis-keywo

    # Write dftd3_geometry file that supplies geometry to dispersion calc
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
    geomfile = './gcp_geometry.xyz'
    with open(geomfile, 'w') as handle:
        handle.write(geomtext)
    # TODO somehow the variations on save_string_xyz and
    #   whether natom and chgmult does or doesn't get written
    #   have gotten all tangled. I fear this doesn't work
    #   the same btwn libmints and qcdb or for ghosts

    # Call gcp program
    command = ['gcp', geomfile]
    command.extend(['-level', func])
    if dertype != 0:
        command.append('-grad')
    try:
        #print('command', command)
        dashout = subprocess.Popen(command, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        raise ValidationError('Program gcp not found in path. %s' % e)
    out, err = dashout.communicate()

    # Parse output
    success = False
    for line in out.splitlines():
        line = line.decode('utf-8')
        if re.match('  Egcp:', line):
            sline = line.split()
            dashd = float(sline[1])
        if re.match('     normal termination of gCP', line):
            success = True

    if not success:
        os.chdir(current_directory)
        raise Dftd3Error("""Unsuccessful gCP run.""")

    # Parse grad output
    if dertype != 0:
        derivfile = './gcp_gradient'
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
            raise ValidationError('Program gcp gradient file has %d atoms- %d expected.' % \
                (len(dashdderiv), self.natom()))

    # Prepare results for Psi4
    if isP4regime and dertype != 0:
        core.set_variable('GCP CORRECTION ENERGY', dashd)
        psi_dashdderiv = core.Matrix(self.natom(), 3)
        psi_dashdderiv.set(dashdderiv)

    # Print program output to file if verbose
    if not verbose and isP4regime:
        verbose = True if core.get_option('SCF', 'PRINT') >= 3 else False
    if verbose:

        text = '\n  ==> GCP Output <==\n'
        text += out.decode('utf-8')
        if dertype != 0:
            with open(derivfile, 'r') as handle:
                text += handle.read().replace('D', 'E')
            text += '\n'
        if isP4regime:
            core.print_out(text)
        else:
            print(text)

#    # Clean up files and remove scratch directory
#    os.unlink(paramfile1)
#    os.unlink(paramfile2)
#    os.unlink(geomfile)
#    if dertype != 0:
#        os.unlink(derivfile)
#    if defmoved is True:
#        os.rename(defaultfile + '_hide', defaultfile)

    os.chdir('..')
#    try:
#        shutil.rmtree(dftd3_tmpdir)
#    except OSError as e:
#        ValidationError('Unable to remove dftd3 temporary directory %s' % e)
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
    core.Molecule.run_gcp = run_gcp
except (NameError, AttributeError):
    # But don't worry if that doesn't work b/c
    #   it'll get attached to qcdb.Molecule class
    pass
