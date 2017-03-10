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

"""Module with utility functions used by several Python functions."""
from __future__ import print_function
import os
import ast
import sys
import pickle
import collections
import inspect
from .exceptions import *
from . import p4regex


if sys.version_info[0] > 2:
    basestring = str

def kwargs_lower(kwargs):
    """Function to rebuild and return *kwargs* dictionary
    with all keys made lowercase. Should be called by every
    function that could be called directly by the user.
    Also turns boolean-like values into actual booleans.
    Also turns values lowercase if sensible.

    """
    caseless_kwargs = {}
    # items() inefficient on Py2 but this is small dict
    for key, value in kwargs.items():
        lkey = key.lower()
        if lkey in ['subset']:  # only kw for which case matters
            lvalue = value
        else:
            try:
                lvalue = value.lower()
            except (AttributeError, KeyError):
                lvalue = value

        if lkey in ['irrep', 'check_bsse', 'linkage', 'bsse_type']:
            caseless_kwargs[lkey] = lvalue

        elif 'dertype' in lkey:
            if p4regex.der0th.match(str(lvalue)):
                caseless_kwargs[lkey] = 0
            elif p4regex.der1st.match(str(lvalue)):
                caseless_kwargs[lkey] = 1
            elif p4regex.der2nd.match(str(lvalue)):
                caseless_kwargs[lkey] = 2
            else:
                raise KeyError('Derivative type key %s was not recognized' % str(key))

        elif p4regex.yes.match(str(lvalue)):
            caseless_kwargs[lkey] = True

        elif p4regex.no.match(str(lvalue)):
            caseless_kwargs[lkey] = False

        else:
            caseless_kwargs[lkey] = lvalue
    return caseless_kwargs


def get_psifile(fileno, pidspace=str(os.getpid())):
    """Function to return the full path and filename for psi file
    *fileno* (e.g., psi.32) in current namespace *pidspace*.

    """
    psioh = core.IOManager.shared_object()
    psio = core.IO.shared_object()
    filepath = psioh.get_file_path(fileno)
    namespace = psio.get_default_namespace()
    targetfile = filepath + 'psi' + '.' + pidspace + '.' + namespace + '.' + str(fileno)
    return targetfile


def format_molecule_for_input(mol, name='', forcexyz=False):
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
    # when mol is core.Molecule or qcdb.Molecule object
    else:
        # save_string_for_psi4 is the more detailed choice as it includes fragment
        #   (and possibly no_com/no_reorient) info. but this is only available
        #   for qcdb Molecules. Since save_string_xyz was added to libmints just
        #   for the sow/reap purpose, may want to unify these fns sometime.
        # the time for unification is nigh
        if forcexyz:
            mol_string = mol.save_string_xyz()
        else:
            mol_string = mol.create_psi4_string_from_molecule()
        mol_name = mol.name() if name == '' else name

    commands = """\nmolecule %s {\n%s%s\n}\n""" % (mol_name, mol_string,
               '\nno_com\nno_reorient' if forcexyz else '')
    return commands


def format_options_for_input(molecule=None, **kwargs):
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Does not cover local (as opposed to global) options.

    """
    if molecule is not None:
        symmetry = molecule.find_point_group(0.00001).symbol()
    commands = ''
    commands += """\ncore.set_memory_bytes(%s)\n\n""" % (core.get_memory())
    for chgdopt in core.get_global_option_list():
        if core.has_global_option_changed(chgdopt):
            chgdoptval = core.get_global_option(chgdopt)

            if molecule is not None:
                if chgdopt.lower() in kwargs:
                    if symmetry in kwargs[chgdopt.lower()]:
                        chgdoptval = kwargs[chgdopt.lower()][symmetry]

            if isinstance(chgdoptval, basestring):
                commands += """core.set_global_option('%s', '%s')\n""" % (chgdopt, chgdoptval)
# Next four lines were conflict between master and roa branches (TDC, 10/29/2014)
            elif isinstance(chgdoptval, int) or isinstance(chgdoptval, float):
                commands += """core.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
            elif isinstance(chgdoptval, list):
                commands += """core.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
            else:
                commands += """core.set_global_option('%s', %s)\n""" % (chgdopt, chgdoptval)
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
    filename.write('''\npickle_kw = ("""'''.encode('utf-8'))
    pickle.dump(kwargs, filename)
    filename.write('''""")\n'''.encode('utf-8'))
    filename.write("""\nkwargs = pickle.loads(pickle_kw)\n""".encode('utf-8'))
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


def extract_sowreap_from_output(sowout, quantity, sownum, linkage, allvital=False, label='electronic energy'):
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
            s = line.strip().split(None, 10)
            if (len(s) != 0) and (s[0:3] == [quantity, 'RESULT:', 'computation']):
                if int(s[3]) != linkage:
                    raise ValidationError('Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.'
                        % (sowout, str(s[3]), str(linkage)))
                if s[6] != str(sownum + 1):
                    raise ValidationError('Output file \'%s.out\' has nominal affiliation %s incompatible with item %s.'
                        % (sowout, s[6], str(sownum + 1)))
                if label == 'electronic energy' and s[8:10] == ['electronic', 'energy']:
                        E = float(s[10])
                        core.print_out('%s RESULT: electronic energy = %20.12f\n' % (quantity, E))
                if label == 'electronic gradient' and s[8:10] == ['electronic', 'gradient']:
                        E = ast.literal_eval(s[-1])
                        core.print_out('%s RESULT: electronic gradient = %r\n' % (quantity, E))
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
    for opt in core.get_global_option_list():
        if core.has_global_option_changed(opt) or not changedOnly:
            if opt in ['DFT_CUSTOM_FUNCTIONAL', 'EXTERN']:  # Feb 2017 hack
                continue
            val = core.get_global_option(opt)
            options['GLOBALS'][opt] = {'value': val,
                                       'has_changed': core.has_global_option_changed(opt)}
            if isinstance(val, basestring):
                commands += """core.set_global_option('%s', '%s')\n""" % (opt, val)
            else:
                commands += """core.set_global_option('%s', %s)\n""" % (opt, val)
            #if changedOnly:
            #    print('Appending module %s option %s value %s has_changed %s.' % \
            #        ('GLOBALS', opt, core.get_global_option(opt), core.has_global_option_changed(opt)))
        for module in modules:
            try:
                if core.has_option_changed(module, opt) or not changedOnly:
                    if not module in options:
                        options[module] = {}
                    val = core.get_option(module, opt)
                    options[module][opt] = {'value': val,
                                            'has_changed': core.has_option_changed(module, opt)}
                    if isinstance(val, basestring):
                        commands += """core.set_local_option('%s', '%s', '%s')\n""" % (module, opt, val)
                    else:
                        commands += """core.set_local_option('%s', '%s', %s)\n""" % (module, opt, val)
                    #if changedOnly:
                    #    print('Appending module %s option %s value %s has_changed %s.' % \
                    #        (module, opt, core.get_option(module, opt), core.has_option_changed(module, opt)))
            except RuntimeError:
                pass

    if commandsInsteadDict:
        return commands
    else:
        return options


def mat2arr(mat):
    """Function to convert core.Matrix *mat* to Python array of arrays.
    Expects core.Matrix to be flat with respect to symmetry.

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
    commands += """memory %d mb\n\n""" % (int(0.000001 * core.get_memory()))
    if not allButMol:
        molecule = core.get_active_molecule()
        molecule.update_geometry()
        commands += format_molecule_for_input(molecule)
        commands += '\n'
    commands += prepare_options_for_modules(changedOnly=True, commandsInsteadDict=True)
    commands += """\n%s('%s', """ % (func.__name__, name.lower())
    for key in kwargs.keys():
        commands += """%s=%r, """ % (key, kwargs[key])
    commands += ')\n\n'

    return commands


def expand_psivars(pvdefs):
    """Dictionary *pvdefs* has keys with names of PsiVariables to be
    created and values with dictionary of two keys: 'args', the
    PsiVariables that contribute to the key and 'func', a function (or
    lambda) to combine them. This function builds those PsiVariables if
    all the contributors are available. Helpful printing is available when
    PRINT > 2.

    """
    verbose = core.get_global_option('PRINT')

    for pvar, action in pvdefs.items():
        if verbose >= 2:
            print("""building %s %s""" % (pvar, '.' * (50 - len(pvar))), end='')

        psivars = core.get_variables()
        data_rich_args = []

        for pv in action['args']:
            if isinstance(pv, basestring):
                if pv in psivars:
                    data_rich_args.append(psivars[pv])
                else:
                    if verbose >= 2:
                        print("""EMPTY, missing {}""".format(pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action['func'](data_rich_args)
            core.set_variable(pvar, result)
            if verbose >= 2:
                print("""SUCCESS""")
