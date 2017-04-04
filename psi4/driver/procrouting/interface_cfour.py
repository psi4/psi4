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

"""Module with functions for Psi4/Cfour interface. Portions that require
calls to Boost Python psi4 module are here, otherwise in qcdb module.
Also calls to qcdb module are here and not elsewhere in driver.
Organizationally, this module isolates qcdb code from psi4 code.

"""
from __future__ import print_function
from __future__ import absolute_import
import shutil
import os
import subprocess
import re
import sys
import uuid
import inspect

from psi4.driver import qcdb
from psi4.driver import p4util
from psi4.driver.molutil import *
from psi4.driver.p4util.exceptions import *
# never import driver, wrappers, or aliases into this file

P4C4_INFO = {}


def run_cfour(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Stanton and Gauss's CFOUR code.
    Also processes results back into Psi4 format.

    This function is not called directly but is instead called by
    :py:func:`~driver.energy` or :py:func:`~driver.optimize` when a Cfour
    method is requested (through *name* argument). In order to function
    correctly, the Cfour executable ``xcfour`` must be present in
    :envvar:`PATH` or :envvar:`PSIPATH`.

    .. hlist::
       :columns: 1

       * Many :ref:`PSI Variables <apdx:cfour_psivar>` extracted from the Cfour output
       * Python dictionary of associated file constants accessible as ``P4C4_INFO['zmat']``, ``P4C4_INFO['output']``, ``P4C4_INFO['grd']``, *etc.*


    :type name: string
    :param name: ``'c4-scf'`` || ``'c4-ccsd(t)'`` || ``'cfour'`` || etc.

        First argument, usually unlabeled. Indicates the computational
        method to be applied to the system.

    :type keep: :ref:`boolean <op_py_boolean>`
    :param keep: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether to delete the Cfour scratch directory upon
        completion of the Cfour job.

    :type path: string
    :param path:

        Indicates path to Cfour scratch directory (with respect to Psi4
        scratch directory). Otherwise, the default is a subdirectory
        within the Psi4 scratch directory.

        If specified, GENBAS and/or ZMAT within will be used.

    :type genbas: string
    :param genbas:

        Indicates that contents should be used for GENBAS file.

    GENBAS is a complicated topic. It is quite unnecessary if the
    molecule is from a molecule {...} block and basis is set through
    |Psifours| BASIS keyword. In that case, a GENBAS is written from
    LibMints and all is well. Otherwise, a GENBAS is looked for in
    the usual places: PSIPATH, PATH, PSIDATADIR/basis. If path kwarg is
    specified, also looks there preferentially for a GENBAS. Can
    also specify GENBAS within an input file through a string and
    setting the genbas kwarg. Note that due to the input parser's
    aggression, blank lines need to be replaced by the text blankline.

    """
    lowername = name.lower()
    internal_p4c4_info = {}
    return_wfn = kwargs.pop('return_wfn', False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    optstash = p4util.OptionsState(
            ['CFOUR', 'TRANSLATE_PSI4'])

    # Determine calling function and hence dertype
    calledby = inspect.stack()[1][3]
    dertype = ['energy', 'gradient', 'hessian'].index(calledby)
    #print('I am %s called by %s called by %s.\n' %
    #    (inspect.stack()[0][3], inspect.stack()[1][3], inspect.stack()[2][3]))

    # Save submission directory
    current_directory = os.getcwd()

    # Move into job scratch directory
    psioh = core.IOManager.shared_object()
    psio = core.IO.shared_object()
    os.chdir(psioh.get_default_path())

    # Construct and move into cfour subdirectory of job scratch directory
    cfour_tmpdir = kwargs['path'] if 'path' in kwargs else \
        'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
        '.cfour.' + str(uuid.uuid4())[:8]
    if not os.path.exists(cfour_tmpdir):
        os.mkdir(cfour_tmpdir)
    os.chdir(cfour_tmpdir)

    # Find environment by merging PSIPATH and PATH environment variables
    lenv = {
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
                ':' + os.environ.get('PATH') + \
                ':' + os.environ.get("PSIDATADIR") + '/basis',
        'GENBAS_PATH': os.environ.get("PSIDATADIR") + '/basis',
        'CFOUR_NUM_CORES': os.environ.get('CFOUR_NUM_CORES'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }

    if 'path' in kwargs:
        lenv['PATH'] = kwargs['path'] + ':' + lenv['PATH']
    #   Filter out None values as subprocess will fault on them
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # Load the GENBAS file
    genbas_path = qcdb.search_file('GENBAS', lenv['GENBAS_PATH'])
    if genbas_path:
        try:
            shutil.copy2(genbas_path, psioh.get_default_path() + cfour_tmpdir)
        except shutil.Error:  # should only fail if src and dest equivalent
            pass
        core.print_out("\n  GENBAS loaded from %s\n" % (genbas_path))
        core.print_out("  CFOUR to be run from %s\n" % (psioh.get_default_path() + cfour_tmpdir))
    else:
        message = """
  GENBAS file for CFOUR interface not found. Either:
  [1] Supply a GENBAS by placing it in PATH or PSIPATH
      [1a] Use cfour {} block with molecule and basis directives.
      [1b] Use molecule {} block and CFOUR_BASIS keyword.
  [2] Allow Psi4's internal basis sets to convert to GENBAS
      [2a] Use molecule {} block and BASIS keyword.

"""
        core.print_out(message)
        core.print_out('  Search path that was tried:\n')
        core.print_out(lenv['PATH'].replace(':', ', '))

    # Generate the ZMAT input file in scratch
    if 'path' in kwargs and os.path.isfile('ZMAT'):
        core.print_out("  ZMAT loaded from %s\n" % (psioh.get_default_path() + kwargs['path'] + '/ZMAT'))
    else:
        with open('ZMAT', 'w') as cfour_infile:
            cfour_infile.write(write_zmat(lowername, dertype, molecule))

    internal_p4c4_info['zmat'] = open('ZMAT', 'r').read()
    #core.print_out('\n====== Begin ZMAT input for CFOUR ======\n')
    #core.print_out(open('ZMAT', 'r').read())
    #core.print_out('======= End ZMAT input for CFOUR =======\n\n')
    #print('\n====== Begin ZMAT input for CFOUR ======')
    #print(open('ZMAT', 'r').read())
    #print('======= End ZMAT input for CFOUR =======\n')

    if 'genbas' in kwargs:
        with open('GENBAS', 'w') as cfour_basfile:
            cfour_basfile.write(kwargs['genbas'].replace('\nblankline\n', '\n\n'))
        core.print_out('  GENBAS loaded from kwargs string\n')

    # Close psi4 output file and reopen with filehandle
    print('output in', current_directory + '/' + core.outfile_name())
    pathfill = '' if os.path.isabs(core.outfile_name()) else current_directory + os.path.sep

    # Handle user's OMP_NUM_THREADS and CFOUR_OMP_NUM_THREADS
    omp_num_threads_found = 'OMP_NUM_THREADS' in os.environ
    if omp_num_threads_found == True:
        omp_num_threads_user = os.environ['OMP_NUM_THREADS']
    if core.has_option_changed('CFOUR', 'CFOUR_OMP_NUM_THREADS') == True:
        os.environ['OMP_NUM_THREADS'] = str(core.get_option('CFOUR', 'CFOUR_OMP_NUM_THREADS'))

    #print("""\n\n<<<<<  RUNNING CFOUR ...  >>>>>\n\n""")
    # Call executable xcfour, directing cfour output to the psi4 output file
    cfour_executable = kwargs['c4exec'] if 'c4exec' in kwargs else 'xcfour'
    try:
        retcode = subprocess.Popen([cfour_executable], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        sys.stderr.write('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        message = ('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        raise ValidationError(message)

    c4out = ''
    while True:
        data = retcode.stdout.readline()
        data = data.decode('utf-8')
        if not data:
            break
        core.print_out(data)
        c4out += data
    internal_p4c4_info['output'] = c4out

    # Restore user's OMP_NUM_THREADS
    if omp_num_threads_found == True:
        if core.has_option_changed('CFOUR', 'CFOUR_OMP_NUM_THREADS') == True:
            os.environ['OMP_NUM_THREADS'] = omp_num_threads_user

    c4files = {}
    core.print_out('\n')
    for item in ['GRD', 'FCMFINAL', 'DIPOL']:
        try:
            with open(psioh.get_default_path() + cfour_tmpdir + '/' + item, 'r') as handle:
                c4files[item] = handle.read()
                core.print_out('  CFOUR scratch file %s has been read\n' % (item))
                core.print_out('%s\n' % c4files[item])
                internal_p4c4_info[item.lower()] = c4files[item]
        except IOError:
            pass
    core.print_out('\n')

    if molecule.name() == 'blank_molecule_psi4_yo':
        qcdbmolecule = None
    else:
        molecule.update_geometry()
        qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
        qcdbmolecule.update_geometry()

    # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
    psivar, c4grad, c4mol = qcdb.cfour.harvest(qcdbmolecule, c4out, **c4files)

    # Absorb results into psi4 data structures
    for key in psivar.keys():
        core.set_variable(key.upper(), float(psivar[key]))

    if qcdbmolecule is None and c4mol is not None:
        molecule = geometry(c4mol.create_psi4_string_from_molecule(), name='blank_molecule_psi4_yo')
        molecule.update_geometry()
        # This case arises when no Molecule going into calc (cfour {} block) but want
        #   to know the orientation at which grad, properties, etc. are returned (c4mol).
        #   c4mol is dinky, w/o chg, mult, dummies and retains name
        #   blank_molecule_psi4_yo so as to not interfere with future cfour {} blocks

    if c4grad:
        mat = core.Matrix(len(c4grad), 3)
        mat.set(c4grad)
        core.set_gradient(mat)

        #print '    <<<   [3] C4-GRD-GRAD   >>>'
        #mat.print()

#    exit(1)

    # # Things needed core.so module to do
    # collect c4out string
    # read GRD
    # read FCMFINAL
    # see if theres an active molecule

    # # Things delegatable to qcdb
    # parsing c4out
    # reading GRD and FCMFINAL strings
    # reconciling p4 and c4 molecules (orient)
# reconciling c4out and GRD and FCMFINAL results
# transforming frame of results back to p4

    # # Things run_cfour needs to have back
    # psivar
# qcdb.Molecule of c4?
# coordinates?
# gradient in p4 frame





#    # Process the cfour output
#    psivar, c4coord, c4grad = qcdb.cfour.cfour_harvest(c4out)
#    for key in psivar.keys():
#        core.set_variable(key.upper(), float(psivar[key]))
#
#    # Awful Hack - Go Away TODO
#    if c4grad:
#        molecule = core.get_active_molecule()
#        molecule.update_geometry()
#
#        if molecule.name() == 'blank_molecule_psi4_yo':
#            p4grad = c4grad
#            p4coord = c4coord
#        else:
#            qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
#            #p4grad = qcdbmolecule.deorient_array_from_cfour(c4coord, c4grad)
#            #p4coord = qcdbmolecule.deorient_array_from_cfour(c4coord, c4coord)
#
#            with open(psioh.get_default_path() + cfour_tmpdir + '/GRD', 'r') as cfour_grdfile:
#                c4outgrd = cfour_grdfile.read()
#            print('GRD\n',c4outgrd)
#            c4coordGRD, c4gradGRD = qcdb.cfour.cfour_harvest_files(qcdbmolecule, grd=c4outgrd)
#
#        p4mat = core.Matrix(len(p4grad), 3)
#        p4mat.set(p4grad)
#        core.set_gradient(p4mat)

#    print('    <<<  P4 PSIVAR  >>>')
#    for item in psivar:
#        print('       %30s %16.8f' % (item, psivar[item]))
    #print('    <<<  P4 COORD   >>>')
    #for item in p4coord:
    #    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    print('    <<<   P4 GRAD   >>>')
#    for item in c4grad:
#        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))














    # Clean up cfour scratch directory unless user instructs otherwise
    keep = yes.match(str(kwargs['keep'])) if 'keep' in kwargs else False
    os.chdir('..')
    try:
        if keep or ('path' in kwargs):
            core.print_out('\n  CFOUR scratch files have been kept in %s\n' % (psioh.get_default_path() + cfour_tmpdir))
        else:
            shutil.rmtree(cfour_tmpdir)
    except OSError as e:
        print('Unable to remove CFOUR temporary directory %s' % e, file=sys.stderr)
        exit(1)

    # Return to submission directory and reopen output file
    os.chdir(current_directory)

    core.print_out('\n')
    p4util.banner(' Cfour %s %s Results ' % (name.lower(), calledby.capitalize()))
    core.print_variables()
    if c4grad:
        core.get_gradient().print_out()

    core.print_out('\n')
    p4util.banner(' Cfour %s %s Results ' % (name.lower(), calledby.capitalize()))
    core.print_variables()
    if c4grad:
        core.get_gradient().print_out()

    # Quit if Cfour threw error
    if core.get_variable('CFOUR ERROR CODE'):
        raise ValidationError("""Cfour exited abnormally.""")

    P4C4_INFO.clear()
    P4C4_INFO.update(internal_p4c4_info)

    optstash.restore()

    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    #   Feb 2017 hack. Could get proper basis in skel wfn even if not through p4 basis kw
    gobas = core.get_global_option('BASIS') if core.get_global_option('BASIS') else 'sto-3g'
    basis = core.BasisSet.build(molecule, "ORBITAL", gobas)
    wfn = core.Wavefunction(molecule, basis)

    optstash.restore()

    if dertype == 0:
        finalquantity = psivar['CURRENT ENERGY']
    elif dertype == 1:
        finalquantity = core.get_gradient()
        wfn.set_gradient(finalquantity)
        if finalquantity.rows(0) < 20:
            core.print_out('CURRENT GRADIENT')
            finalquantity.print_out()
    elif dertype == 2:
        pass
        #finalquantity = finalhessian
        #wfn.set_hessian(finalquantity)
        #if finalquantity.rows(0) < 20:
        #    core.print_out('CURRENT HESSIAN')
        #    finalquantity.print_out()

    return wfn


def cfour_list():
    """Form list of Cfour :py:func:`~driver.energy` arguments."""
    return qcdb.cfour.cfour_list()


def cfour_gradient_list():
    """Form list of Cfour analytic :py:func:`~driver.gradient` arguments."""
    return qcdb.cfour.cfour_gradient_list()


def cfour_psivar_list():
    """Form dictionary of :ref:`PSI Variables <apdx:cfour_psivar>` set by Cfour methods."""
    return qcdb.cfour.cfour_psivar_list()


def write_zmat(name, dertype, molecule):
    """Returns string with contents of Cfour ZMAT file as gathered from
    active molecule, current keyword settings, and cfour {...} block.

    """
    # Handle memory
    mem = int(0.000001 * core.get_memory())
    if mem == 256:
        memcmd, memkw = '', {}
    else:
        memcmd, memkw = qcdb.cfour.muster_memory(mem)

    # Handle molecule and basis set
    if molecule.name() == 'blank_molecule_psi4_yo':
        molcmd, molkw = '', {}
        bascmd, baskw = '', {}
        core.set_local_option('CFOUR', 'TRANSLATE_PSI4', False)
    else:
        molecule.update_geometry()
        #print(molecule.create_psi4_string_from_molecule())
        qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
        qcdbmolecule.tagline = molecule.name()
        molcmd, molkw = qcdbmolecule.format_molecule_for_cfour()

        if core.get_global_option('BASIS') == '':
            bascmd, baskw = '', {}
        else:
            user_pg = molecule.schoenflies_symbol()
            molecule.reset_point_group('c1')  # need basis printed for *every* atom
            qbs = core.BasisSet.build(molecule, "BASIS", core.get_global_option('BASIS'))
            with open('GENBAS', 'w') as cfour_basfile:
                cfour_basfile.write(qbs.genbas())
            core.print_out('  GENBAS loaded from Psi4 LibMints for basis %s\n' % (core.get_global_option('BASIS')))
            molecule.reset_point_group(user_pg)
            molecule.update_geometry()
            bascmd, baskw = qcdbmolecule.format_basis_for_cfour(qbs.has_puream())

    # Handle psi4 keywords implying cfour keyword values
    if core.get_option('CFOUR', 'TRANSLATE_PSI4'):
        psicmd, psikw = qcdb.cfour.muster_psi4options(p4util.prepare_options_for_modules(changedOnly=True))
    else:
        psicmd, psikw = '', {}

    # Handle calc type and quantum chemical method
    mdccmd, mdckw = qcdb.cfour.muster_modelchem(name, dertype)

    # Handle calc type and quantum chemical method
    mdccmd, mdckw = qcdb.cfour.muster_modelchem(name, dertype)

    # Handle driver vs input/default keyword reconciliation
    userkw = p4util.prepare_options_for_modules()
    userkw = qcdb.options.reconcile_options(userkw, memkw)
    userkw = qcdb.options.reconcile_options(userkw, molkw)
    userkw = qcdb.options.reconcile_options(userkw, baskw)
    userkw = qcdb.options.reconcile_options(userkw, psikw)
    userkw = qcdb.options.reconcile_options(userkw, mdckw)

    # Handle conversion of psi4 keyword structure into cfour format
    optcmd = qcdb.options.prepare_options_for_cfour(userkw)

    # Handle text to be passed untouched to cfour
    litcmd = core.get_global_option('LITERAL_CFOUR')

    # Assemble ZMAT pieces
    zmat = memcmd + molcmd + optcmd + mdccmd + psicmd + bascmd + litcmd

    if len(re.findall(r'^\*(ACES2|CFOUR|CRAPS)\(', zmat, re.MULTILINE)) != 1:
        core.print_out('\n  Faulty ZMAT constructed:\n%s' % (zmat))
        raise ValidationError("""
Multiple *CFOUR(...) blocks in input. This usually arises
because molecule or options are specified both the psi4 way through
molecule {...} and set ... and the cfour way through cfour {...}.""")

    return zmat
