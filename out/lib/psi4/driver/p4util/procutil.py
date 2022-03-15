#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module with utility functions used by several Python functions."""
import os
import ast
import sys
import math
import pickle
import inspect
import warnings
import contextlib
import collections
from typing import List, Union

import numpy as np
import qcelemental as qcel

from psi4 import core
from psi4.metadata import __version__
from .exceptions import ValidationError
from . import p4regex


def kwargs_lower(kwargs):
    """Function to rebuild and return *kwargs* dictionary
    with all keys made lowercase. Should be called by every
    function that could be called directly by the user.
    Also turns boolean-like values into actual booleans.
    Also turns values lowercase if sensible.

    """
    caseless_kwargs = {}
    for key, value in kwargs.items():
        lkey = key.lower()
        if lkey in ['subset', 'banner', 'restart_file', 'write_orbitals']:  # only kw for which case matters
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
                raise KeyError(f'Derivative type key {key} was not recognized')

        elif lvalue is None:
            caseless_kwargs[lkey] = None

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
    if isinstance(mol, str):
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

    commands = """\nmolecule %s {\n%s%s\n}\n""" % (mol_name, mol_string, '\nno_com\nno_reorient' if forcexyz else '')
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

            if isinstance(chgdoptval, str):
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
    warnings.warn(
        "Using `psi4.driver.p4util.format_kwargs_for_input` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.get_legacy_gradient()

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
    """Return a generator of all lettercase permutations of `input_string`."""

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
    warnings.warn(
        "Using `psi4.driver.p4util.extract_sowreap_from_output` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    E = 0.0

    try:
        freagent = open('%s.out' % (sowout), 'r')
    except IOError:
        if allvital:
            raise ValidationError('Aborting upon output file \'%s.out\' not found.\n' % (sowout))
        else:
            return 0.0
    else:
        while True:
            line = freagent.readline()
            if not line:
                if E == 0.0:
                    raise ValidationError(
                        'Aborting upon output file \'%s.out\' has no %s RESULT line.\n' % (sowout, quantity))
                break
            s = line.strip().split(None, 10)
            if (len(s) != 0) and (s[0:3] == [quantity, 'RESULT:', 'computation']):
                if int(s[3]) != linkage:
                    raise ValidationError(
                        'Output file \'%s.out\' has linkage %s incompatible with master.in linkage %s.' %
                        (sowout, str(s[3]), str(linkage)))
                if s[6] != str(sownum + 1):
                    raise ValidationError(
                        'Output file \'%s.out\' has nominal affiliation %s incompatible with item %s.' %
                        (sowout, s[6], str(sownum + 1)))
                if label == 'electronic energy' and s[8:10] == ['electronic', 'energy']:
                    E = float(s[10])
                    core.print_out('%s RESULT: electronic energy = %20.12f\n' % (quantity, E))
                if label == 'electronic gradient' and s[8:10] == ['electronic', 'gradient']:
                    E = ast.literal_eval(s[-1])
                    core.print_out('%s RESULT: electronic gradient = %r\n' % (quantity, E))
        freagent.close()
    return E


_modules = [
    # Psi4 Modules
    "ADC",
    "CCENERGY",
    "CCEOM",
    "CCDENSITY",
    "CCLAMBDA",
    "CCHBAR",
    "CCRESPONSE",
    "CCTRANSORT",
    "CCTRIPLES",
    "CPHF",
    "DCT",
    "DETCI",
    "DFEP2",
    "DFMP2",
    "DFOCC",
    "DLPNO",
    "DMRG",
    "EFP",
    "FINDIF",
    "FISAPT",
    "FNOCC",
    "GDMA",
    "MCSCF",
    "MINTS",
    "MRCC",
    "OCC",
    "OPTKING",
    "PCM",
    "PE",
    "PSIMRCC",
    "RESPONSE",
    "SAPT",
    "SCF",
    "THERMO",
    # External Modules
    "CFOUR",
]


def reset_pe_options(pofm):
    """Acts on Process::environment.options to clear it, the set it to state encoded in `pofm`.

    Parameters
    ----------
    pofm : dict
        Result of psi4.driver.p4util.prepare_options_for_modules(changedOnly=True, commandsInsteadDict=False)

    Returns
    -------
    None

    """
    core.clean_options()

    for go, dgo in pofm['GLOBALS'].items():
        if dgo['has_changed']:
            core.set_global_option(go, dgo['value'])

    for module in _modules:
        for lo, dlo in pofm[module].items():
            if dlo['has_changed']:
                core.set_local_option(module, lo, dlo['value'])


def prepare_options_for_modules(changedOnly=False, commandsInsteadDict=False):
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Need some option to get either all or changed

       - Need some option to either get dict or set string or psimod command list

       - command return doesn't revoke has_changed setting for unchanged with changedOnly=False

    """
    options = collections.defaultdict(dict)
    commands = ''
    for opt in core.get_global_option_list():
        if core.has_global_option_changed(opt) or not changedOnly:
            if opt in ['DFT_CUSTOM_FUNCTIONAL', 'EXTERN']:  # Feb 2017 hack
                continue
            val = core.get_global_option(opt)
            options['GLOBALS'][opt] = {'value': val, 'has_changed': core.has_global_option_changed(opt)}
            if isinstance(val, str):
                commands += """core.set_global_option('%s', '%s')\n""" % (opt, val)
            else:
                commands += """core.set_global_option('%s', %s)\n""" % (opt, val)
            #if changedOnly:
            #    print('Appending module %s option %s value %s has_changed %s.' % \
            #        ('GLOBALS', opt, core.get_global_option(opt), core.has_global_option_changed(opt)))
        for module in _modules:
            if core.option_exists_in_module(module, opt):
                hoc = core.has_option_changed(module, opt)
                if hoc or not changedOnly:
                    val = core.get_option(module, opt)
                    options[module][opt] = {'value': val, 'has_changed': hoc}
                    if isinstance(val, str):
                        commands += """core.set_local_option('%s', '%s', '%s')\n""" % (module, opt, val)
                    else:
                        commands += """core.set_local_option('%s', '%s', %s)\n""" % (module, opt, val)
                    #if changedOnly:
                    #    print('Appending module %s option %s value %s has_changed %s.' % \
                    #        (module, opt, core.get_option(module, opt), hoc))

    if commandsInsteadDict:
        return commands
    else:
        return options


def prepare_options_for_set_options():
    """Capture current state of C++ psi4.core.Options information for reloading by `psi4.set_options()`.

    Returns
    -------
    dict
        Dictionary where keys are option names to be set globally or module__option
        mangled names to be set locally. Values are option values.

    """
    flat_options = {}
    has_changed_snapshot = {module: core.options_to_python(module) for module in _modules}

    for opt in core.get_global_option_list():

        handled_locally = False
        ghoc = core.has_global_option_changed(opt)
        opt_snapshot = {k: v[opt] for k, v in has_changed_snapshot.items() if opt in v}
        for module, (lhoc, ohoc) in opt_snapshot.items():
            if ohoc:
                if lhoc:
                    key = module + '__' + opt
                    val = core.get_local_option(module, opt)
                else:
                    key = opt
                    val = core.get_global_option(opt)
                    handled_locally = True
                flat_options[key] = val

        if ghoc and not handled_locally:
            # some options are globals section (not level) so not in any module
            flat_options[opt] = core.get_global_option(opt)

    return flat_options


def state_to_atomicinput(*, driver, method, basis=None, molecule=None, function_kwargs=None) -> "AtomicInput":
    """Form a QCSchema for job input from the current state of Psi4 settings."""

    if molecule is None:
        molecule = core.get_active_molecule()

    keywords = {k.lower(): v for k, v in prepare_options_for_set_options().items()}
    if function_kwargs is not None:
        keywords["function_kwargs"] = function_kwargs

    kw_basis = keywords.pop("basis", None)
    basis = basis or kw_basis

    resi = qcel.models.AtomicInput(
         **{
            "driver": driver,
            "extras": {
                "wfn_qcvars_only": True,
            },
            "model": {
                "method": method,
                "basis": basis,
            },
            "keywords": keywords,
            "molecule": molecule.to_schema(dtype=2),
            "provenance": provenance_stamp(__name__),
         })

    return resi


def mat2arr(mat):
    """Function to convert core.Matrix *mat* to Python array of arrays.
    Expects core.Matrix to be flat with respect to symmetry.

    """
    warnings.warn(
        "Using `psi4.driver.p4util.mat2arr` instead of `MatrixInstance.to_array().tolist()` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

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
    warnings.warn(
        "Using `psi4.driver.p4util.format_currentstate_for_input` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

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

        psivars = core.scalar_variables()
        data_rich_args = []

        for pv in action['args']:
            if isinstance(pv, str):
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


def provenance_stamp(routine):
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with Psi4's credentials for creator and version. The
    generating routine's name is passed in through `routine`.

    """
    return {'creator': 'Psi4', 'version': __version__, 'routine': routine}


def plump_qcvar(val: Union[float, str, List], shape_clue: str, ret: str = 'np') -> Union[float, np.ndarray, core.Matrix]:
    """Prepare serialized QCVariable for set_variable by convert flat arrays into shaped ones and floating strings.

    Parameters
    ----------
    val :
        flat (?, ) list or scalar or string, probably from JSON storage.
    shape_clue
        Label that includes (case insensitive) one of the following as
        a clue to the array's natural dimensions: 'gradient', 'hessian'
    ret
        {'np', 'psi4'}
        Whether to return `np.ndarray` or `psi4.core.Matrix`.

    Returns
    -------
    float or numpy.ndarray or Matrix
        Reshaped array of type `ret` with natural dimensions of `shape_clue`.

    """
    if isinstance(val, (np.ndarray, core.Matrix)):
        raise TypeError
    elif isinstance(val, list):
        tgt = np.asarray(val)
    else:
        # presumably scalar. may be string
        return float(val)
    # TODO choose float vs Decimal for return if string?

    if 'gradient' in shape_clue.lower():
        reshaper = (-1, 3)
    elif 'hessian' in shape_clue.lower():
        ndof = int(math.sqrt(len(tgt)))
        reshaper = (ndof, ndof)
    else:
        raise ValidationError(f'Uncertain how to reshape array: {shape_clue}')

    if ret == 'np':
        return tgt.reshape(reshaper)
    elif ret == 'psi4':
        return core.Matrix.from_array(tgt.reshape(reshaper))
    else:
        raise ValidationError(f'Return type not among [np, psi4]: {ret}')
