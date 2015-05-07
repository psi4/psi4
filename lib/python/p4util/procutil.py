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

"""Module with utility functions used by several Python functions."""
import os
import sys
import pickle
import collections
import inspect
#CUimport psi4
#CUimport inputparser
from p4xcpt import *


if sys.version_info[0] > 2:
    basestring = str

def kwargs_lower(kwargs):
    """Function to rebuild and return *kwargs* dictionary
    with all keys made lowercase. Should be called by every
    function that could be called directly by the user.

    """
    caseless_kwargs = {}
    if sys.hexversion < 0x03000000:
        # Python 2; we have to explicitly use an iterator
        for key, value in kwargs.items():
            caseless_kwargs[key.lower()] = value
    else:
        # Python 3; an iterator is implicit
        for key, value in kwargs.items():
            caseless_kwargs[key.lower()] = value
    return caseless_kwargs


def get_psifile(fileno, pidspace=str(os.getpid())):
    """Function to return the full path and filename for psi file
    *fileno* (e.g., psi.32) in current namespace *pidspace*.

    """
    psioh = psi4.IOManager.shared_object()
    psio = psi4.IO.shared_object()
    filepath = psioh.get_file_path(fileno)
    namespace = psio.get_default_namespace()
    targetfile = filepath + 'psi' + '.' + pidspace + '.' + namespace + '.' + str(fileno)
    return targetfile


def format_molecule_for_input(mol, name=''):
    """Function to return a string of the output of
    :py:func:`inputparser.process_input` applied to the XYZ
    format of molecule, passed as either fragmented
    geometry string *mol* or molecule instance *mol*.
    Used to capture molecule information from database
    modules and for distributed (sow/reap) input files.
    For the reverse, see :py:func:`molutil.geometry`.

    """
    # when mol is already a string
    if isinstance(mol, basestring):
        mol_string = mol
        mol_name = name
    # when mol is psi4.Molecule or qcdb.Molecule object
    else:
        # save_string_for_psi4 is the more detailed choice as it includes fragment
        #   (and possibly no_com/no_reorient) info. but this is only available
        #   for qcdb Molecules. Since save_string_xyz was added to libmints just
        #   for the sow/reap purpose, may want to unify these fns sometime.
        # the time for unification is nigh
        mol_string = mol.create_psi4_string_from_molecule()
        #try:
        #    mol_string = mol.save_string_for_psi4()
        #except AttributeError:
        #    mol_string = mol.save_string_xyz()

        mol_name = mol.name()

    commands = """\nmolecule %s {\n%s\n}\n""" % (mol_name, mol_string)
    return commands
    #commands = 'inputparser.process_input("""\nmolecule %s {\n%s\n}\n""", 0)\n' % (mol_name, mol_string)
    #return eval(commands)


def format_options_for_input():
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Does not cover local (as opposed to global) options.

    """
    commands = ''
    commands += """\npsi4.set_memory(%s)\n\n""" % (psi4.get_memory())
    for chgdopt in psi4.get_global_option_list():
        if psi4.has_global_option_changed(chgdopt):
            chgdoptval = psi4.get_global_option(chgdopt)
            if isinstance(chgdoptval, basestring):
                commands += """psi4.set_global_option('%s', '%s')\n""" % (chgdopt, chgdoptval)
# Next four lines were conflict between master and roa branches (TDC, 10/29/2014)
            elif isinstance(chgdoptval, int) or isinstance(chgdoptval, float):
                commands += """psi4.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
            elif isinstance(chgdoptval, list):
                commands += """psi4.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
            else:
                commands += """psi4.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
    return commands


def format_kwargs_for_input(filename, lmode=1, **kwargs):
    """Function to pickle to file *filename* the options dictionary
    *kwargs*. Mode *lmode* =2 pickles appropriate settings for
    reap mode. Used to capture Python options information for
    distributed (sow/reap) input files.

    """
    if lmode == 2:
        kwargs['mode'] = 'reap'
        kwargs['linkage'] = os.getpid()
    filename.write('''\npickle_kw = ("""''')
    pickle.dump(kwargs, filename)
    filename.write('''""")\n''')
    filename.write("""\nkwargs = pickle.loads(pickle_kw)\n""")
    if lmode == 2:
        kwargs['mode'] = 'sow'
        del kwargs['linkage']


def drop_duplicates(seq):
    """Function that given an array *seq*, returns an array without any duplicate
    entries. There is no guarantee of which duplicate entry is dropped.

    """
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes


def all_casings(input_string):
    """Function to return a generator of all lettercase permutations
    of *input_string*.

    """
    if not input_string:
        yield ""
    else:
        first = input_string[:1]
        if first.lower() == first.upper():
            for sub_casing in all_casings(input_string[1:]):
                yield first + sub_casing
        else:
            for sub_casing in all_casings(input_string[1:]):
                yield first.lower() + sub_casing
                yield first.upper() + sub_casing


def getattr_ignorecase(module, attr):
    """Function to extract attribute *attr* from *module* if *attr*
    is available in any possible lettercase permutation. Returns
    attribute if available, None if not.

    """
    array = None
    for per in list(all_casings(attr)):
        try:
            getattr(module, per)
        except AttributeError:
            pass
        else:
            array = getattr(module, per)
            break

    return array


def import_ignorecase(module):
    """Function to import *module* in any possible lettercase
    permutation. Returns module object if available, None if not.

    """
    modobj = None
    for per in list(all_casings(module)):
        try:
            modobj = __import__(per)
        except ImportError:
            pass
        else:
            break

    return modobj


def extract_sowreap_from_output(sowout, quantity, sownum, linkage, allvital=False):
    """Function to examine file *sowout* from a sow/reap distributed job
    for formatted line with electronic energy information about index
    *sownum* to be used for construction of *quantity* computations as
    directed by master input file with *linkage* kwarg. When file *sowout*
    is missing or incomplete files, function will either return zero
    (*allvital* is ``False``) or terminate (*allvital* is ``True``) since
    some sow/reap procedures can produce meaningful results (database)
    from an incomplete set of sown files, while others cannot (gradient,
    hessian).

    """
    E = 0.0

    try:
        freagent = open('%s.out' % (sowout), 'r')
    except IOError:
        if allvital:
            raise ValidationError('Aborting upon output file \'%s.out\' not found.\n' % (sowout))
        else:
            ValidationError('Aborting upon output file \'%s.out\' not found.\n' % (sowout))
            return 0.0
    else:
        while True:
            line = freagent.readline()
            if not line:
                if E == 0.0:
                    if allvital:
                        raise ValidationError('Aborting upon output file \'%s.out\' has no %s RESULT line.\n' % (sowout, quantity))
                    else:
                        ValidationError('Aborting upon output file \'%s.out\' has no %s RESULT line.\n' % (sowout, quantity))
                break
            s = line.split()
            if (len(s) != 0) and (s[0:3] == [quantity, 'RESULT:', 'computation']):
                if int(s[3]) != linkage:
                    raise ValidationError('Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.'
                        % (sowout, str(s[3]), str(linkage)))
                if s[6] != str(sownum + 1):
                    raise ValidationError('Output file \'%s.out\' has nominal affiliation %s incompatible with item %s.'
                        % (sowout, s[6], str(sownum + 1)))
                if (s[8:10] == ['electronic', 'energy']):
                    E = float(s[10])
                    psi4.print_out('%s RESULT: electronic energy = %20.12f\n' % (quantity, E))
        freagent.close()
    return E

def prepare_options_for_modules(changedOnly=False, commandsInsteadDict=False):
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Need some option to get either all or changed

       - Need some option to either get dict or set string or psimod command list

       - command return doesn't revoke has_changed setting for unchanged with changedOnly=False

    """
    modules = [
        # PSI4 Modules
        "ADC", "CCENERGY", "CCEOM", "CCDENSITY", "CCLAMBDA", "CCHBAR",
        "CCRESPONSE", "CCSORT", "CCTRIPLES", "CLAG", "CPHF", "CIS",
        "DCFT", "DETCI", "DFMP2", "DFTSAPT", "FINDIF", "FNOCC", "LMP2",
        "MCSCF", "MINTS", "MRCC", "OCC", "OPTKING", "PSIMRCC", "RESPONSE",
        "SAPT", "SCF", "STABILITY", "THERMO", "TRANSQT", "TRANSQT2",
        # External Modules
        "CFOUR",
        ]

    options = {'GLOBALS': {}}
    commands = ''
    for opt in psi4.get_global_option_list():
        if psi4.has_global_option_changed(opt) or not changedOnly:
            val = psi4.get_global_option(opt)
            options['GLOBALS'][opt] = {'value': val,
                                       'has_changed': psi4.has_global_option_changed(opt)}
            if isinstance(val, basestring):
                commands += """psi4.set_global_option('%s', '%s')\n""" % (opt, val)
            else:
                commands += """psi4.set_global_option('%s', %s)\n""" % (opt, val)
            #if changedOnly:
            #    print('Appending module %s option %s value %s has_changed %s.' % \
            #        ('GLOBALS', opt, psi4.get_global_option(opt), psi4.has_global_option_changed(opt)))
        for module in modules:
            try:
                if psi4.has_option_changed(module, opt) or not changedOnly:
                    if not module in options:
                        options[module] = {}
                    val = psi4.get_option(module, opt)
                    options[module][opt] = {'value': val,
                                            'has_changed': psi4.has_option_changed(module, opt)}
                    if isinstance(val, basestring):
                        commands += """psi4.set_local_option('%s', '%s', '%s')\n""" % (module, opt, val)
                    else:
                        commands += """psi4.set_local_option('%s', '%s', %s)\n""" % (module, opt, val)
                    #if changedOnly:
                    #    print('Appending module %s option %s value %s has_changed %s.' % \
                    #        (module, opt, psi4.get_option(module, opt), psi4.has_option_changed(module, opt)))
            except RuntimeError:
                pass

    if commandsInsteadDict:
        return commands
    else:
        return options


def mat2arr(mat):
    """Function to convert psi4.Matrix *mat* to Python array of arrays.
    Expects psi4.Matrix to be flat with respect to symmetry.

    """
    if mat.rowdim().n() != 1:
        raise ValidationError('Cannot convert Matrix with symmetry.')
    arr = []
    for row in range(mat.rowdim()[0]):
        temp = []
        for col in range(mat.coldim()[0]):
            temp.append(mat.get(row, col))
        arr.append(temp)
    return arr


def format_currentstate_for_input(func, name, allButMol=False, **kwargs):
    """Function to return an input file in preprocessed psithon.
    Captures memory, molecule, options, function, method, and kwargs.
    Used to write distributed (sow/reap) input files.

    """
    commands = """\n# This is a psi4 input file auto-generated from the %s() wrapper.\n\n""" % (inspect.stack()[1][3])
    commands += """memory %d mb\n\n""" % (int(0.000001 * psi4.get_memory()))
    if not allButMol:
        molecule = psi4.get_active_molecule()
        molecule.update_geometry()
        commands += format_molecule_for_input(molecule)
        commands += '\n'
    commands += prepare_options_for_modules(changedOnly=True, commandsInsteadDict=True)
    commands += """\n%s('%s', """ % (func.__name__, name.lower())
    for key in kwargs.keys():
        commands += """%s=%r, """ % (key, kwargs[key])
    commands += ')\n\n'

    return commands
