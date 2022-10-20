#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

__all__ = [
    "all_casings",
    "drop_duplicates",
    "expand_psivars",
    "format_molecule_for_input",
    "format_options_for_input",
    "get_psifile",
    "getattr_ignorecase",
    "hold_options_state",
    "import_ignorecase",
    "kwargs_lower",
    "mat2arr",
    "prepare_options_for_modules",
    "prepare_options_for_set_options",
    "provenance_stamp",
    "plump_qcvar",
    "state_to_atomicinput",
]

import os
import ast
import sys
import math
import pickle
import inspect
import warnings
from contextlib import contextmanager
import collections
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Union
from types import ModuleType

import numpy as np
from qcelemental.models import AtomicInput

from psi4 import core
from psi4.metadata import __version__
from .exceptions import ValidationError
from . import p4regex


def kwargs_lower(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Function to rebuild and return *kwargs* dictionary sanitized. Should be
    called by every function that could be called directly by the user.

    Parameters
    ----------
    kwargs
        Input kwargs for any user-facing function.

    Returns
    -------
    lowered : Dict[str, Any]
        Sanitized kwargs with all keys made lowercase. Also turns boolean-like
        values into actual booleans. Also turns values lowercase if sensible.

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


def get_psifile(fileno: int, pidspace: str = str(os.getpid())) -> str:
    """Form full path and filename for psi scratch file.

    Parameters
    ----------
    fileno
        Psi file, e.g., ``psi.32``.
    pidspace
        Current namespace. Defaults to ``os.getpid()``.

    Returns
    -------
    flpath : str
        Full path and filename for psi file.

    """
    psioh = core.IOManager.shared_object()
    psio = core.IO.shared_object()
    filepath = psioh.get_file_path(fileno)
    namespace = psio.get_default_namespace()
    targetfile = filepath + 'psi' + '.' + pidspace + '.' + namespace + '.' + str(fileno)
    return targetfile


def format_molecule_for_input(
    mol: Union[str, core.Molecule],
    name: str = '',
    forcexyz: bool = False) -> str:
    """Function to return a string of the output of
    :py:func:`~psi4.driver.inputparser.process_input` applied to the XYZ
    format of molecule, passed as either fragmented
    geometry string *mol* or molecule instance *mol*.
    Used to capture molecule information from database
    modules and for distributed (sow/reap) input files.
    For the reverse, see :py:func:`~psi4.driver.geometry`.

    Parameters
    ----------
    mol
        Fragmented geometry string or molecule instance.
    name
        Name to call the resulting molecule.
    forcexyz
        Use Cartesians, even for Z-Matrix molecules.

    """
    warnings.warn(
        "Using `psi4.driver.p4util.format_molecule_for_input` instead of `Molecule.to_string(dtype='psi4')` is deprecated, and as soon as 1.8 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

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


def format_options_for_input(molecule: Optional[core.Molecule] = None, **kwargs) -> str:
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Does not cover local (as opposed to global) options.

    """
    warnings.warn(
        "Using `psi4.driver.p4util.format_molecule_for_input` instead of `Molecule.to_string(dtype='psi4')` is deprecated, and as soon as 1.8 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

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


def drop_duplicates(seq: Iterable) -> List:
    """Return a copy of collection *seq* without any duplicate entries.

    Parameters
    ----------
    seq
        Collection to be de-duplicated. There is no guarantee of which
        duplicate entry is dropped.

    """
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes


def all_casings(input_string: str) -> Iterator[str]:
    """Return a generator of all lettercase permutations of `input_string`.

    Parameters
    ----------
    input_string
        String of which to permute the case.

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


def getattr_ignorecase(module: str, attr: str):
    """Extract attribute *attr* from *module* if *attr*
    is available in any possible lettercase permutation.

    Parameters
    ----------
    module
        Object on which to seek `attr`.
    attr
        Name of attribute with uncertain case.

    Returns
    -------
    attribute : Any
        Module attribute returned if available. None if not.

    """
    obj_attr = None
    for permutation in list(all_casings(attr)):
        try:
            getattr(module, permutation)
        except AttributeError:
            pass
        else:
            obj_attr = getattr(module, permutation)
            break

    return obj_attr


def import_ignorecase(module: str) -> ModuleType:
    """Import loader for *module* in any possible lettercase permutation.

    Parameters
    ----------
    module
        Name of module with uncertain case.

    Returns
    -------
    types.ModuleType
        Module object.

    """
    modobj = None
    for permutation in list(all_casings(module)):
        try:
            modobj = __import__(permutation)
        except ImportError:
            pass
        else:
            break

    return modobj


_modules = [
    # Psi4 Modules
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


@contextmanager
def hold_options_state() -> Iterator[None]:
    """Return a context manager that will collect the current state of
    ``Process:environment.options`` on entry to the with-statement and clear
    and restore the collected keywords state when exiting the with-statement.

    """
    pofm = prepare_options_for_modules(
        changedOnly=True, commandsInsteadDict=False, globalsOnly=False, stateInsteadMediated=True
    )
    yield
    _reset_pe_options(pofm)


def _reset_pe_options(pofm: Dict):
    """Acts on ``Process::environment.options`` to clear it, then set it to
    state encoded in `pofm`.

    Parameters
    ----------
    pofm
        Result of :py:func:`prepare_options_for_modules` with
        ``changedOnly=True``, ``commandsInsteadDict=False``, and
        ``stateInsteadMediated=True``.

    """
    core.clean_options()

    for go, dgo in pofm['GLOBALS'].items():
        if dgo['has_changed']:
            core.set_global_option(go, dgo['value'])

    for module in _modules:
        for lo, dlo in pofm[module].items():
            if dlo['has_changed']:
                core.set_local_option(module, lo, dlo['value'])

    ## this is a more defensive version if defaults may have changed
    # for go, dgo in pofm['GLOBALS'].items():
    #    core.set_global_option(go, dgo['value'])
    #    if not dgo['has_changed']:
    #        core.revoke_global_option_changed(go)
    # for module in _modules:
    #    for lo, dlo in pofm[module].items():
    #        core.set_local_option(module, lo, dlo['value'])
    #        if not dlo['has_changed']:
    #            core.revoke_local_option_changed(module, lo)


def prepare_options_for_modules(
    changedOnly: bool = False,
    commandsInsteadDict: bool = False,
    globalsOnly: bool = False,
    stateInsteadMediated: bool = False,
) -> Union[Dict, str]:
    """Capture current state of :py:class:`psi4.core.Options` information.

    Parameters
    ----------
    changedOnly
        Record info only for options that have been set (may still be default).
        When False, records values for every option.
    commandsInsteadDict
        Return string of commands to exec to reset options in current form.
        When False, return nested dictionary with globals in 'GLOBALS' subdictionary
        and locals in subdictionaries by module.
    globalsOnly
        Record only global options to save time querying the
        :py:class:`~psi4.core.Options` object.
        When False, record module-level options, too.
    stateInsteadMediated
        When ``True``, querying this function for options to be later *reset* into the same
        state -- the raw values and has_changed status at the global and local levels.
        When ``False``, querying this function for mediated options to be *used* -- the results
        of the globals/locals handshake as computed by the Options object itself. Here,
        ``dict[module][option][value]`` is the value to be used by module.

    Returns
    -------
    Dict
        When `commandsInsteadDict` is False.
    str
        When `commandsInsteadDict` is True.


    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - command return doesn't revoke has_changed setting for unchanged with changedOnly=False

       - not all kwargs are independent

    """
    has_changed_snapshot = {module: core.options_to_python(module) for module in _modules}
    options = collections.defaultdict(dict)
    commands = ''
    for opt in core.get_global_option_list():
        hoc = core.has_global_option_changed(opt)
        if hoc or not changedOnly:
            if opt in ['DFT_CUSTOM_FUNCTIONAL', 'EXTERN']:  # Feb 2017 hack
                continue
            val = core.get_global_option(opt)
            options['GLOBALS'][opt] = {'value': val, 'has_changed': hoc}
            if isinstance(val, str):
                commands += """core.set_global_option('%s', '%s')\n""" % (opt, val)
            else:
                commands += """core.set_global_option('%s', %s)\n""" % (opt, val)
        if globalsOnly:
            continue

        opt_snapshot = {k: v[opt] for k, v in has_changed_snapshot.items() if opt in v}
        for module, (lhoc, ohoc) in opt_snapshot.items():
            if stateInsteadMediated:
                hoc = lhoc
            else:
                hoc = ohoc
            if hoc or not changedOnly:
                if stateInsteadMediated:
                    val = core.get_local_option(module, opt)
                else:
                    val = core.get_option(module, opt)
                options[module][opt] = {'value': val, 'has_changed': hoc}
                if isinstance(val, str):
                    commands += """core.set_local_option('%s', '%s', '%s')\n""" % (module, opt, val)
                else:
                    commands += """core.set_local_option('%s', '%s', %s)\n""" % (module, opt, val)

    if commandsInsteadDict:
        return commands
    else:
        return options


def prepare_options_for_set_options() -> Dict[str, Any]:
    """Collect current state of :py:class:`psi4.core.Options` information for
    reloading by :py:func:`~psi4.driver.p4util.set_options`.

    Returns
    -------
    Dict[str, Any]
        Dictionary where keys are keyword names, either plain for those to be
        set globally or mangled "module__keyword" for those to be set locally,
        and values are keyword values.

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

    # The normal machinery to forward plugin options to Psi goes through 'plugin_load'.
    # Forte doesn't use this. Pending a larger options rewrite (move to a Python dictionary?),
    # we need the following dirty hack.

    try:
        import forte # Needed for Forte options to run.
    except ImportError:
        pass
    else:
        # Initialization tasks with Psi options
        psi_options = core.get_options()
        current_module = psi_options.get_current_module()
        # Get the current Forte options from Forte
        forte_options = forte.ForteOptions()
        forte.register_forte_options(forte_options)
        psi_options.set_current_module("FORTE")
        try:
            forte_options.get_options_from_psi4(psi_options)
        except RuntimeError:
            # If we're in this case, Forte hasn't pushed its options to Psi.
            pass
        else:
            # Load changed Forte options into `flat_options`
            for name, metadata in forte_options.dict().items():
                if metadata["value"] != metadata["default_value"]:
                    flat_options[f"forte__{name.lower()}"] = metadata["value"]
        finally:
            # Restore current module
            psi_options.set_current_module(current_module)
    return flat_options


def state_to_atomicinput(
    *,
    driver: str,
    method: str,
    basis: Optional[str] = None,
    molecule: Optional[core.Molecule] = None,
    function_kwargs: Optional[Dict[str, Any]] = None) -> AtomicInput:
    """Form a QCSchema for job input from the current state of |PSIfour| settings.

    Parameters
    ----------
    driver
        {'energy', 'gradient', 'hessian'}
        Target derivative level.
    method
        Level of theory for job.
    basis
        Basis set for job, if not to be extracted from :term:`BASIS <BASIS (MINTS)>`.
    molecule
        Molecule for job, if not the active one from
        :py:func:`~psi4.core.get_active_molecule`.
    function_kwargs
        Additional keyword arguments to pass to the driver function.

    Returns
    -------
    ~qcelemental.models.AtomicInput
        QCSchema instance including current keyword set and provenance.

    """
    if molecule is None:
        molecule = core.get_active_molecule()

    keywords = {k.lower(): v for k, v in prepare_options_for_set_options().items()}
    if function_kwargs is not None:
        keywords["function_kwargs"] = function_kwargs

    kw_basis = keywords.pop("basis", None)
    basis = basis or kw_basis

    resi = AtomicInput(
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


def mat2arr(mat: core.Matrix) -> List[List[float]]:
    """Convert Matrix to List.

    Parameters
    ----------
    mat
        |PSIfour| matrix. Should be flat with respect to symmetry.

    Returns
    -------
    List[List[float]]
        Pure Python representation of `mat`.

    """
    warnings.warn(
        "Using `psi4.driver.p4util.mat2arr` instead of `MatrixInstance.to_array().tolist()` is deprecated, and as soon as 1.4 it will stop working\n",
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


def expand_psivars(
    pvdefs: Dict[str, Dict[str, Union[List[str], Callable]]],
    verbose: Optional[int] = None):
    """From rules on building QCVariables from others, set new variables to
    P::e if all the contributors are available.

    Parameters
    ----------
    pvdefs
        Dictionary with keys with names of QCVariables to be created and values
        with dictionary of two keys: 'args', the QCVariables that contribute to
        the key and 'func', a function (or lambda) to combine them.
    verbose
        Control print level. If unspecified (None), value taken from
        :term:`PRINT <PRINT (GLOBALS)>`. Status printing when verbose > 2.

    Examples
    --------
    >>> pv1 = dict()
    >>> pv1['SAPT CCD DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4],
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP2(CCD) ENERGY',
                                     'SAPT DISP22(S)(CCD) ENERGY', 'SAPT EST.DISP22(T)(CCD) ENERGY']}
    >>> pv1['SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY']}
    >>> expand_psivars(pv1)

    """
    if verbose is None:
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


def provenance_stamp(routine: str, module: str = None) -> Dict[str, str]:
    """Prepare QCSchema Provenance with |PSIfour| credentials.

    Parameters
    ----------
    routine
        Name of driver function generating the QCSchema.
    module
        Primary contributing |PSIfour| library, like ``ccenergy`` or ``dfmp2``.

    Returns
    -------
    provenance : Dict[str, str]
        Dictionary satisfying QCSchema, with |PSIfour| credentials for creator
        and version.
        https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41

    """
    prov = {'creator': 'Psi4', 'version': __version__, 'routine': routine}
    if module:
        prov["module"] = module

    return prov


def plump_qcvar(
    val: Union[float, str, List],
    shape_clue: str,
    ret: str = 'np') -> Union[float, np.ndarray, core.Matrix]:
    """Prepare serialized QCVariable for :py:func:`~psi4.core.set_variable` by
    converting flat arrays into shaped ones and floating strings.

    Parameters
    ----------
    val
        flat (?, ) list or scalar or string, probably from JSON storage.
    shape_clue
        Label that includes (case insensitive) one of the following as
        a clue to the array's natural dimensions: 'gradient', 'hessian'
    ret
        {'np', 'psi4'}
        Whether for arrays to return :py:class:`numpy.ndarray` or
        :py:class:`psi4.core.Matrix`.

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
