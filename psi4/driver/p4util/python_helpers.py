#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

"""
Module with PsiAPI helpers for PSIthon `{...}` syntax.
Also, many Python extensions to core classes:

 - core (variable-related, gradient, python option),
 - Wavefunction (variable-related, freq, Lagrangian, constructor, scratch file, serialization),
 - Matrix (doublet, triplet),
 - BasisSet (constructor)
 - JK (constructor)
 - VBase (grid)
 - OEProp (avail prop)
 - ERISieve (constructor)
"""

__all__ = [
    "basis_helper",
    "pcm_helper",
    "plump_qcvar",
    "set_options",
    "set_module_options",
]


import os
import re
import sys
import math
import uuid
import warnings
from collections import Counter
from itertools import product
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

import qcelemental as qcel
from psi4 import core
from psi4 import extras
from psi4.driver import qcdb

from . import optproc
from .exceptions import TestComparisonError, ValidationError, UpgradeHelper

## Python basis helps

@staticmethod
def _pybuild_basis(
        mol: core.Molecule,
        key: Optional[str] = None,
        target: Optional[Union[str, Callable]] = None,
        fitrole: str = "ORBITAL",
        other: Optional[Union[str, Callable]] = None,
        puream: int = -1,
        return_atomlist: bool = False,
        *,
        quiet: bool = False,
    ) -> Union[core.BasisSet, List[core.BasisSet]]:
    """Build a primary or auxiliary basis set.

    Parameters
    ----------
    mol
        Molecule for which to build the basis set instance.
    key
        {'BASIS', 'ORBITAL', 'DF_BASIS_SCF', 'DF_BASIS_MP2', 'DF_BASIS_CC', 'BASIS_RELATIVISTIC', 'DF_BASIS_SAD'}
        Label (effectively Psi4 keyword) to append the basis on the molecule.
        The primary basis set is indicated by any of values None or
        ``"ORBITAL"`` or ``"BASIS"``.
    target
        Defines the basis set to be constructed. Can be a string (naming a
        basis file) or a callable (providing shells or multiple basis files).
        For auxiliary bases to be built entirely from primary default, can be
        an empty string. If None, value taken from `key` in global options. If
        a user-defined-basis callable is available at string `target`, `target`
        value will be set to it. In practice, setting this argument to a
        |PSIfour| keyword (e.g., ``core.get_option("SCF", "DF_BASIS_SCF")`` or
        ``core.get_global_option("BASIS")``) works to handle both simple and
        user-defined bases.
    fitrole
        {'ORBITAL', 'JKFIT', 'RIFIT', 'DECON'}
        Role for which to build basis. Only used when `key` indicates auxiliary
        (i.e., *is not* ``"BASIS"``) and auxiliary spec from processing `target`
        can't complete the `mol`. Then, primary spec from `other` can be used
        to complete the auxiliary basis by looking up suitable default basis
        according to `fitrole`.
    other
        Only used when building auxiliary basis sets. Defines the primary basis through a string or callable like `target`.
    puream
        Whether to override the native spherical/cartesian-ness of `target` for
        returned basis? Value ``1`` forces spherical, value ``0`` forces
        Cartesian, value ``-1`` (default) uses native puream. Note that
        explicitly setting :term:`PUREAM <PUREAM (GLOBALS)>` trumps both native
        puream and this `puream` argument.
    return_atomlist
        Build one-atom basis sets (e.g., for SAD) rather than one whole-`mol`
        basis set.
    quiet
        When True, do not print to the output file.

    Returns
    -------
    BasisSet or ~typing.List[BasisSet]
        Single basis for `mol`, unless `return_atomlist` is True.

    """
    if key == 'ORBITAL':
        key = 'BASIS'

    def _resolve_target(key, target):
        """Figure out exactly what basis set was intended by (key, target)
        """
        horde = qcdb.libmintsbasisset.basishorde
        if not target:
            if not key:
                key = 'BASIS'
            target = core.get_global_option(key)

        if target in horde:
            return horde[target]
        return target

    # Figure out what exactly was meant by 'target'.
    resolved_target = _resolve_target(key, target)

    # resolved_target needs to be either a string or function for pyconstuct.
    # if a string, they search for a gbs file with that name.
    # if a function, it needs to apply a basis to each atom.

    bs, basisdict = qcdb.BasisSet.pyconstruct(mol.to_dict(),
                                              key,
                                              resolved_target,
                                              fitrole,
                                              other,
                                              return_dict=True,
                                              return_atomlist=return_atomlist)

    if return_atomlist:
        atom_basis_list = []
        for atbs in basisdict:
            atommol = core.Molecule.from_dict(atbs['molecule'])
            lmbs = core.BasisSet.construct_from_pydict(atommol, atbs, puream)
            atom_basis_list.append(lmbs)
        return atom_basis_list
    if isinstance(resolved_target, str):
        basisdict['name'] = basisdict['name'].split('/')[-1].replace('.gbs', '')
    if callable(resolved_target):
        basisdict['name'] = resolved_target.__name__.replace('basisspec_psi4_yo__', '').upper()

    if not quiet:
        core.print_out(basisdict['message'])
        if 'ECP' in basisdict['message']:
            core.print_out('    !!!  WARNING: ECP capability is in beta. Please check occupations closely.  !!!\n\n')

    if basisdict['key'] is None:
        basisdict['key'] = 'BASIS'
    psibasis = core.BasisSet.construct_from_pydict(mol, basisdict, puream)
    return psibasis


core.BasisSet.build = _pybuild_basis

## Python wavefunction helps


@staticmethod
def _core_wavefunction_build(
        mol: core.Molecule,
        basis: Union[None, str, core.BasisSet] = None,
        *,
        quiet: bool = False,
    ) -> core.Wavefunction:
    """Build a wavefunction from minimal inputs, molecule and basis set.

    Parameters
    ----------
    mol
        Molecule for which to build the wavefunction instance.
    basis
        Basis set for which to build the wavefunction instance. If a
        :class:`BasisSet`, taken as-is. If a string, taken as a name for the
        primary basis. If None, name taken from :term:`BASIS <BASIS (MINTS)>`.
    quiet
        When True, do not print to the output file.

    """
    if basis is None:
        basis = core.BasisSet.build(mol, quiet=quiet)
    elif isinstance(basis, str):
        basis = core.BasisSet.build(mol, "ORBITAL", basis, quiet=quiet)

    wfn = core.Wavefunction(mol, basis)
    # Set basis for density-fitted calculations to the zero basis...
    # ...until the user explicitly provides a DF basis.
    wfn.set_basisset("DF_BASIS_SCF", core.BasisSet.zero_ao_basis_set())
    return wfn


core.Wavefunction.build = _core_wavefunction_build


def _core_wavefunction_get_scratch_filename(self: core.Wavefunction, filenumber: int) -> str:
    """Return canonical path to scratch file `filenumber` based on molecule on `self`.

    Parameters
    ----------
    self
        Wavefunction instance.
    filenumber
        Scratch file number from :source:`psi4/include/psi4/psifiles.h`.

    """
    fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(self.molecule().name())))[1]
    psi_scratch = core.IOManager.shared_object().get_default_path()
    return os.path.join(psi_scratch, fname + '.' + str(filenumber))


core.Wavefunction.get_scratch_filename = _core_wavefunction_get_scratch_filename


@staticmethod
def _core_wavefunction_from_file(wfn_data: Union[str, Dict, Path]) -> core.Wavefunction:
    r"""Build Wavefunction from data laid out like
    :meth:`~psi4.core.Wavefunction.to_file`.

    Parameters
    ----------
    wfn_data
        If a dict, use data directly. Otherwise, path-like passed to
        :py:func:`numpy.load` to read from disk.

    Returns
    -------
    Wavefunction
        A deserialized Wavefunction object

    """
    # load the wavefunction from file
    if isinstance(wfn_data, dict):
        pass
    elif isinstance(wfn_data, str):
        if not wfn_data.endswith(".npy"):
            wfn_data = wfn_data + ".npy"
        wfn_data = np.load(wfn_data, allow_pickle=True).item()
    else:
        # Could be path-like or file-like, let `np.load` handle it
        wfn_data = np.load(wfn_data, allow_pickle=True).item()

    # variable type specific dictionaries to be passed into C++ constructor
    wfn_matrix = wfn_data['matrix']
    wfn_vector = wfn_data['vector']
    wfn_dimension = wfn_data['dimension']
    wfn_int = wfn_data['int']
    wfn_string = wfn_data['string']
    wfn_boolean = wfn_data['boolean']
    wfn_float = wfn_data['float']
    wfn_floatvar = wfn_data['floatvar']
    wfn_matrixarr = wfn_data['matrixarr']

    # reconstruct molecule from dictionary representation
    wfn_molecule = wfn_data['molecule']
    molecule = core.Molecule.from_dict(wfn_molecule)

    # get basis set name and spherical harmonics boolean
    basis_name = wfn_string['basisname']
    if ".gbs" in basis_name:
        basis_name = basis_name.split('/')[-1].replace('.gbs', '')

    basis_puream = wfn_boolean['basispuream']
    basisset = core.BasisSet.build(molecule, 'ORBITAL', basis_name, puream=basis_puream)

    # change some variables to psi4 specific data types (Matrix, Vector, Dimension)
    for label in wfn_matrix:
        array = wfn_matrix[label]
        wfn_matrix[label] = core.Matrix.from_array(array, name=label) if array is not None else None

    for label in wfn_vector:
        array = wfn_vector[label]
        wfn_vector[label] = core.Vector.from_array(array, name=label) if array is not None else None

    for label in wfn_dimension:
        tup = wfn_dimension[label]
        wfn_dimension[label] = core.Dimension.from_list(tup, name=label) if tup is not None else None

    for label in wfn_matrixarr:
        array = wfn_matrixarr[label]
        wfn_matrixarr[label] = core.Matrix.from_array(array, name=label) if array is not None else None

    # make the wavefunction
    wfn = core.Wavefunction(molecule, basisset, wfn_matrix, wfn_vector, wfn_dimension, wfn_int, wfn_string,
                            wfn_boolean, wfn_float)

    # some of the wavefunction's variables can be changed directly
    for k, v in wfn_floatvar.items():
        wfn.set_variable(k, v)
    for k, v in wfn_matrixarr.items():
        wfn.set_variable(k, v)

    return wfn


core.Wavefunction.from_file = _core_wavefunction_from_file


def _core_wavefunction_to_file(wfn: core.Wavefunction, filename: str = None) -> Dict[str, Dict[str, Any]]:
    """Serialize a Wavefunction object. Opposite of
    :meth:`~psi4.core.Wavefunction.from_file`.

    Parameters
    ----------
    wfn
        Wavefunction or inherited class instance.
    filename
        An optional filename to which to write the data.

    Returns
    -------
    ~typing.Dict[str, ~typing.Dict[str, ~typing.Any]]
        A dictionary and NumPy representation of the Wavefunction.

    """

    # collect the wavefunction's variables in a dictionary indexed by varaible type
    # some of the data types have to be made numpy-friendly first
    if wfn.basisset().name().startswith("anonymous"):
        raise ValidationError("Cannot serialize wavefunction with custom basissets.")

    wfn_data = {
        'molecule': wfn.molecule().to_dict(),
        'matrix': {
            'Ca':       wfn.Ca().to_array()       if wfn.Ca()       else None,
            'Cb':       wfn.Cb().to_array()       if wfn.Cb()       else None,
            'Da':       wfn.Da().to_array()       if wfn.Da()       else None,
            'Db':       wfn.Db().to_array()       if wfn.Db()       else None,
            'Fa':       wfn.Fa().to_array()       if wfn.Fa()       else None,
            'Fb':       wfn.Fb().to_array()       if wfn.Fb()       else None,
            'H':        wfn.H().to_array()        if wfn.H()        else None,
            'S':        wfn.S().to_array()        if wfn.S()        else None,
            'X':        wfn.lagrangian().to_array() if wfn.lagrangian() else None,
            'aotoso':   wfn.aotoso().to_array()   if wfn.aotoso()   else None,
            'gradient': wfn.gradient().to_array() if wfn.gradient() else None,
            'hessian':  wfn.hessian().to_array()  if wfn.hessian()  else None
        },
        'vector': {
            'epsilon_a': wfn.epsilon_a().to_array() if wfn.epsilon_a() else None,
            'epsilon_b': wfn.epsilon_b().to_array() if wfn.epsilon_b() else None,
            'frequencies': wfn.frequencies().to_array() if wfn.frequencies() else None
        },
        'dimension': {
            'doccpi':   wfn.doccpi().to_tuple(),
            'frzcpi':   wfn.frzcpi().to_tuple(),
            'frzvpi':   wfn.frzvpi().to_tuple(),
            'nalphapi': wfn.nalphapi().to_tuple(),
            'nbetapi':  wfn.nbetapi().to_tuple(),
            'nmopi':    wfn.nmopi().to_tuple(),
            'nsopi':    wfn.nsopi().to_tuple(),
            'soccpi':   wfn.soccpi().to_tuple()
        },
        'int': {
            'nalpha': wfn.nalpha(),
            'nbeta':  wfn.nbeta(),
            'nfrzc':  wfn.nfrzc(),
            'nirrep': wfn.nirrep(),
            'nmo':    wfn.nmo(),
            'nso':    wfn.nso(),
            'print':  wfn.get_print(),
        },
        'string': {
            'name': wfn.name(),
            'module': wfn.module(),
            'basisname': wfn.basisset().name()
        },
        'boolean': {
            'PCM_enabled':    wfn.PCM_enabled(),
            'same_a_b_dens':  wfn.same_a_b_dens(),
            'same_a_b_orbs':  wfn.same_a_b_orbs(),
            'basispuream':    wfn.basisset().has_puream()
        },
        'float': {
            'energy': wfn.energy(),
            'efzc': wfn.efzc(),
            'dipole_field_x': wfn.get_dipole_field_strength()[0],
            'dipole_field_y': wfn.get_dipole_field_strength()[1],
            'dipole_field_z': wfn.get_dipole_field_strength()[2]
        },
        'floatvar': wfn.scalar_variables(),
        'matrixarr': {k: v.to_array() for k, v in wfn.array_variables().items()}
    }  # yapf: disable

    if filename is not None:
        if not filename.endswith('.npy'): filename += '.npy'
        np.save(filename, wfn_data, allow_pickle=True)

    return wfn_data


core.Wavefunction.to_file = _core_wavefunction_to_file

## Python JK helps


@staticmethod
def _core_jk_build(
        orbital_basis: core.BasisSet,
        aux: Optional[core.BasisSet] = None,
        jk_type: Optional[str] = None,
        do_wK: Optional[bool] = None,
        memory: Optional[int] = None,
    ) -> core.JK:
    """
    Constructs a Psi4 JK object from an input basis.

    Parameters
    ----------
    orbital_basis
        Orbital basis to use in the JK object.
    aux
        Optional auxiliary basis set for density-fitted tensors. Defaults
        to the :term:`DF_BASIS_SCF <DF_BASIS_SCF (SCF)>` if set, otherwise the corresponding JKFIT basis
        to the passed in `orbital_basis`.
    jk_type
        Type of JK object to build (DF, Direct, PK, etc). Defaults to the
        current :term:`SCF_TYPE <SCF_TYPE (GLOBALS)>` option.
    do_wK
        Set up JK to do omega K tasks. Set `do_wK` and `memory` together to
        activate either.
    memory
        Memory in doubles to use for JK. Set `do_wK` and `memory` together to
        activate either.

    Returns
    -------
    JK
        Uninitialized JK object.

    Example
    -------
    >>> jk = psi4.core.JK.build(bas)
    >>> jk.set_memory(int(5e8))  # 4GB of memory
    >>> jk.initialize()
    >>> ...
    >>> jk.C_left_add(matirx)
    >>> jk.compute()
    >>> jk.C_clear()
    >>> ...

    """

    optstash = optproc.OptionsState(["SCF_TYPE"])

    if jk_type is not None:
        core.set_global_option("SCF_TYPE", jk_type)

    if aux is None:
        if core.get_global_option("SCF_TYPE") == "DF":
            aux = core.BasisSet.build(orbital_basis.molecule(), "DF_BASIS_SCF", core.get_option("SCF", "DF_BASIS_SCF"),
                                      "JKFIT", orbital_basis.name(), orbital_basis.has_puream())
        else:
            aux = core.BasisSet.zero_ao_basis_set()

    if (do_wK is None) or (memory is None):
        jk = core.JK.build_JK(orbital_basis, aux)
    else:
        jk = core.JK.build_JK(orbital_basis, aux, bool(do_wK), int(memory))

    optstash.restore()
    return jk


core.JK.build = _core_jk_build

## Grid Helpers


def _core_vbase_get_np_xyzw(self: core.VBase) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns the x, y, z, and weights of a grid as a tuple of NumPy array objects.

    Parameters
    ----------
    self
        VBase instance.

    """
    x_list = []
    y_list = []
    z_list = []
    w_list = []

    # Loop over every block in the potenital
    for b in range(self.nblocks()):

        # Obtain the block
        block = self.get_block(b)

        # Obtain the x, y, and z coordinates along with the weight
        x_list.append(block.x())
        y_list.append(block.y())
        z_list.append(block.z())
        w_list.append(block.w())

    x = np.hstack(x_list)
    y = np.hstack(y_list)
    z = np.hstack(z_list)
    w = np.hstack(w_list)

    return (x, y, z, w)


core.VBase.get_np_xyzw = _core_vbase_get_np_xyzw

## Python other helps


def set_options(options_dict: Dict[str, Any], verbose: int = 1):
    """Sets Psi4 options from an input dictionary.

    Parameters
    ----------
    options_dict
        Dictionary where keys are case insensitive and values are the option value.

        - For global options, keys are ``"<option_name>"``.
        - For option local to "<module_name>", keys are ``"<module_name>__<option_name>"``
          (double underscore separation).
        - For contents that would be in ``pcm = {...}``, use ``"PCM__INPUT"`` key.
    verbose
        Control print volume.

    Returns
    -------
    None

    Examples
    --------
    >>> psi4.set_options({
            "basis": "cc-pvtz",
            "df_basis_scf": "cc-pvtz-jkfit",
            "scf__reference": "uhf",
            "print": 2})

    """
    optionre = re.compile(r'\A(?P<module>\w+__)?(?P<option>\w+)\Z', re.IGNORECASE)
    rejected = {}

    for k, v, in options_dict.items():

        mobj = optionre.match(k.strip())
        module = mobj.group('module').upper()[:-2] if mobj.group('module') else None
        option = mobj.group('option').upper()

        if module:
            if ((module, option, v) not in [('SCF', 'GUESS', 'READ')]) and ((module, option) not in [('PCM', 'INPUT')]):
                # TODO guess/read exception is for distributed driver. should be handled differently.
                try:
                    core.set_local_option(module, option, v)
                except RuntimeError as err:
                    rejected[k] = (v, err)
                if verbose > 1:
                    print('Setting: core.set_local_option', module, option, v)

            if (module, option) == ("PCM", "INPUT"):
                pcm_helper(v)

        else:
            try:
                core.set_global_option(option, v)
            except RuntimeError as err:
                rejected[k] = (v, err)
            if verbose > 1:
                print('Setting: core.set_global_option', option, v)

    if rejected:
        raise ValidationError(f'Error setting options: {rejected}')
        # TODO could subclass ValidationError and append rejected so that run_json could handle remanants.


def set_module_options(module: str, options_dict: Dict[str, Any]) -> None:
    """
    Sets Psi4 module options from a module specification and input dictionary.

    .. deprecated:: 1.5
       Use :py:func:`psi4.driver.p4util.set_options` instead.

    """
    warnings.warn(
        "Using `psi4.set_module_options(<module>, {<key>: <val>})` instead of `psi4.set_options({<module>__<key>: <val>})` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    for k, v, in options_dict.items():
        core.set_local_option(module.upper(), k.upper(), v)


## OEProp helpers


def pcm_helper(block: str):
    """Helper to specify the multiline PCMSolver syntax for PCM.
    Prefer to use :py:func:`set_options` with key ``"PCM__INPUT"``.

    Parameters
    ----------
    block
        Text that goes in a PSIthon ``pcm = {...}`` block.

    """
    import pcmsolver

    # delete=True works for Unix but not for Windows
    with NamedTemporaryFile(mode="w+t", delete=False) as fl:
        fl.write(block)
        fl.flush()
        parsed_pcm = pcmsolver.parse_pcm_input(fl.name)
        extras.register_scratch_file(fl.name)

    with NamedTemporaryFile(mode="w+t", delete=False) as fl:
        fl.write(parsed_pcm)
        core.set_local_option("PCM", "PCMSOLVER_PARSED_FNAME", fl.name)
        extras.register_scratch_file(fl.name)  # retain with -m (messy) option


def _basname(name: str) -> str:
    """Imitates :py:meth:`core.BasisSet.make_filename` without the gbs extension."""
    return name.lower().replace('+', 'p').replace('*', 's').replace('(', '_').replace(')', '_').replace(',', '_')


def basis_helper(block: str, name: str = '', key: str = 'BASIS', set_option: bool = True):
    """Helper to specify a custom basis set in PsiAPI mode.

    This function forms a basis specification function from *block*
    and associates it with keyword *key* under handle *name*. Registers
    the basis spec with Psi4 so that it can be applied again to future
    molecules. For usage, see :srcsample:`mints2`, :srcsample:`mints9`, and
    :srcsample:`cc54` test cases.

    Parameters
    ----------
    block
        Text that goes in a PSIthon ``basis {...}`` block.
    name
        Name label to associated with basis specified by `block`.
    key
        Basis keyword specified by `block`.
    set_option
        When True, execute the equivalent of ``set key name`` or ``set_option({key: name})``. When False, skip execution.

    """
    key = key.upper()
    name = ('anonymous' + str(uuid.uuid4())[:8]) if name == '' else name
    cleanbas = _basname(name).replace('-', '')  # further remove hyphens so can be function name
    block = qcel.util.filter_comments(block)
    command_lines = re.split('\n', block)

    symbol_re = re.compile(r'^\s*assign\s+(?P<symbol>[A-Z]{1,3})\s+(?P<basis>[-+*\(\)\w]+)\s*$', re.IGNORECASE)
    label_re = re.compile(
        r'^\s*assign\s+(?P<label>(?P<symbol>[A-Z]{1,3})(?:(_\w+)|(\d+))?)\s+(?P<basis>[-+*\(\)\w]+)\s*$',
        re.IGNORECASE)
    all_re = re.compile(r'^\s*assign\s+(?P<basis>[-+*\(\)\w]+)\s*$', re.IGNORECASE)
    basislabel = re.compile(r'\s*\[\s*([-*\(\)\w]+)\s*\]\s*')

    def anon(mol, role):
        basstrings = {}

        # Start by looking for assign lines, and remove them
        leftover_lines = []
        assignments = False
        for line in command_lines:
            if symbol_re.match(line):
                m = symbol_re.match(line)
                mol.set_basis_by_symbol(m.group('symbol'), m.group('basis'), role=role)
                assignments = True

            elif label_re.match(line):
                m = label_re.match(line)
                mol.set_basis_by_label(m.group('label'), m.group('basis'), role=role)
                assignments = True

            elif all_re.match(line):
                m = all_re.match(line)
                mol.set_basis_all_atoms(m.group('basis'), role=role)
                assignments = True

            else:
                # Ignore blank lines and accumulate remainder
                if line and not line.isspace():
                    leftover_lines.append(line.strip())

        # Now look for regular basis set definitions
        basblock = list(filter(None, basislabel.split('\n'.join(leftover_lines))))
        if len(basblock) == 1:
            if not assignments:
                # case with no [basname] markers where whole block is contents of gbs file
                mol.set_basis_all_atoms(name, role=role)
                basstrings[_basname(name)] = basblock[0]
            else:
                message = (
                    "Conflicting basis set specification: assign lines present but shells have no [basname] label."
                    "")
                raise TestComparisonError(message)
        else:
            # case with specs separated by [basname] markers
            for idx in range(0, len(basblock), 2):
                basstrings[_basname(basblock[idx])] = basblock[idx + 1]

        return basstrings

    anon.__name__ = 'basisspec_psi4_yo__' + cleanbas
    qcdb.libmintsbasisset.basishorde[name.upper()] = anon
    if set_option:
        core.set_global_option(key, name)


core.OEProp.valid_methods = [
    'DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'WIBERG_LOWDIN_INDICES', 'MAYER_INDICES',
    'MBIS_CHARGES','MBIS_VOLUME_RATIOS', 'MO_EXTENTS', 'GRID_FIELD', 'GRID_ESP', 'ESP_AT_NUCLEI', 'NO_OCCUPATIONS'
]

## Option helpers


def _core_set_global_option_python(key, EXTERN):
    """
    This is a fairly hacky way to get around EXTERN issues. Effectively we are routing this option Python side through attributes until the general Options overhaul.
    """
    if (key != "EXTERN"):
        raise ValidationError("Options: set_global_option_python does not recognize keyword %s" % key)

    if EXTERN is None:
        core.EXTERN = None
        core.set_global_option("EXTERN", False)
    elif isinstance(EXTERN, core.ExternalPotential):
        # Well this is probably the worst hack I have done, thats saying something
        core.EXTERN = EXTERN
        core.set_global_option("EXTERN", True)
    else:
        raise ValidationError("Options: set_global_option_python can either be a NULL or External Potential object")


core.set_global_option_python = _core_set_global_option_python

## QCvar helps

_qcvar_transitions = {
    # old: (replacement, release after next)
    "SCSN-MP2 CORRELATION ENERGY": ("SCS(N)-MP2 CORRELATION ENERGY", 1.5),
    "SCSN-MP2 TOTAL ENERGY": ("SCS(N)-MP2 TOTAL ENERGY", 1.5),
    "MAYER_INDICES": ("MAYER INDICES", 1.5),
    "WIBERG_LOWDIN_INDICES": ("WIBERG LOWDIN INDICES", 1.5),
    "LOWDIN_CHARGES": ("LOWDIN CHARGES", 1.5),
    "MULLIKEN_CHARGES": ("MULLIKEN CHARGES", 1.5),
    "(AT) CORRECTION ENERGY": ("A-(T) CORRECTION ENERGY", 1.5),
    "CCSD(AT) TOTAL ENERGY": ("A-CCSD(T) TOTAL ENERGY", 1.5),
    "CCSD(AT) CORRELATION ENERGY": ("A-CCSD(T) CORRELATION ENERGY", 1.5),
    "CP-CORRECTED 2-BODY INTERACTION ENERGY": ("CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY", 1.7),
    "CP-CORRECTED 3-BODY INTERACTION ENERGY": ("CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY", 1.7),
    "CP-CORRECTED 4-BODY INTERACTION ENERGY": ("CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY", 1.7),
    "CP-CORRECTED 5-BODY INTERACTION ENERGY": ("CP-CORRECTED INTERACTION ENERGY THROUGH 5-BODY", 1.7),
    "NOCP-CORRECTED 2-BODY INTERACTION ENERGY": ("NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY", 1.7),
    "NOCP-CORRECTED 3-BODY INTERACTION ENERGY": ("NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY", 1.7),
    "NOCP-CORRECTED 4-BODY INTERACTION ENERGY": ("NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY", 1.7),
    "NOCP-CORRECTED 5-BODY INTERACTION ENERGY": ("NOCP-CORRECTED INTERACTION ENERGY THROUGH 5-BODY", 1.7),
    "VMFC-CORRECTED 2-BODY INTERACTION ENERGY": ("VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY", 1.7),
    "VMFC-CORRECTED 3-BODY INTERACTION ENERGY": ("VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY", 1.7),
    "VMFC-CORRECTED 4-BODY INTERACTION ENERGY": ("VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY", 1.7),
    "VMFC-CORRECTED 5-BODY INTERACTION ENERGY": ("VMFC-CORRECTED INTERACTION ENERGY THROUGH 5-BODY", 1.7),
    "COUNTERPOISE CORRECTED TOTAL ENERGY": ("CP-CORRECTED TOTAL ENERGY", 1.7),
    "COUNTERPOISE CORRECTED INTERACTION ENERGY": ("CP-CORRECTED INTERACTION ENERGY", 1.7),
    "NON-COUNTERPOISE CORRECTED TOTAL ENERGY": ("NOCP-CORRECTED TOTAL ENERGY", 1.7),
    "NON-COUNTERPOISE CORRECTED INTERACTION ENERGY": ("NOCP-CORRECTED INTERACTION ENERGY", 1.7),
    "VALIRON-MAYER FUNCTION COUTERPOISE TOTAL ENERGY": ("VALIRON-MAYER FUNCTION COUNTERPOISE TOTAL ENERGY", 1.7),  # note misspelling
    "VALIRON-MAYER FUNCTION COUTERPOISE INTERACTION ENERGY": ("VMFC-CORRECTED INTERACTION ENERGY", 1.7),  # note misspelling
}

_qcvar_cancellations = {
    "SCSN-MP2 SAME-SPIN CORRELATION ENERGY": ["MP2 SAME-SPIN CORRELATION ENERGY"],
    "SCSN-MP2 OPPOSITE-SPIN CORRELATION ENERGY": ["MP2 OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCS-CCSD SAME-SPIN CORRELATION ENERGY": ["CCSD SAME-SPIN CORRELATION ENERGY"],
    "SCS-CCSD OPPOSITE-SPIN CORRELATION ENERGY": ["CCSD OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCS-MP2 SAME-SPIN CORRELATION ENERGY": ["MP2 SAME-SPIN CORRELATION ENERGY"],
    "SCS-MP2 OPPOSITE-SPIN CORRELATION ENERGY": ["MP2 OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCS(N)-OMP2 CORRELATION ENERGY": ["OMP2 SAME-SPIN CORRELATION ENERGY", "OMP2 OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCS(N)-OMP2 TOTAL ENERGY": ["OMP2 SAME-SPIN CORRELATION ENERGY", "OMP2 OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCSN-OMP2 CORRELATION ENERGY": ["OMP2 SAME-SPIN CORRELATION ENERGY", "OMP2 OPPOSITE-SPIN CORRELATION ENERGY"],
    "SCSN-OMP2 TOTAL ENERGY": ["OMP2 SAME-SPIN CORRELATION ENERGY", "OMP2 OPPOSITE-SPIN CORRELATION ENERGY"],
}


def _qcvar_warnings(key: str) -> str:
    """Intercept QCVariable keys to issue warnings or upgrade hints. Otherwise,
    pass through.

    """
    if any([key.upper().endswith(" DIPOLE " + cart) for cart in ["X", "Y", "Z"]]):
        raise UpgradeHelper(key.upper(), key.upper()[:-2], 1.6, " Note the Debye -> a.u. units change.")

    if any([key.upper().endswith(" QUADRUPOLE " + cart) for cart in ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]]):
        raise UpgradeHelper(key.upper(), key.upper()[:-3], 1.6, " Note the Debye -> a.u. units change.")

    if key.upper() in _qcvar_transitions:
        replacement, version = _qcvar_transitions[key.upper()]
        warnings.warn(
            f"Using QCVariable `{key.upper()}` instead of `{replacement}` is deprecated, and as soon as {version} it will stop working\n",
            category=FutureWarning,
            stacklevel=3)
        return replacement

    if key.upper() in _qcvar_cancellations:
        raise UpgradeHelper(key.upper(), "no direct replacement", 1.4, " Consult QCVariables " + ", ".join(_qcvar_cancellations[key.upper()]) + " to recompose the quantity.")

    return key


def plump_qcvar(
    key: str,
    val: Union[float, str, List]) -> Union[float, np.ndarray]:
    """Prepare serialized QCVariables for QCSchema AtomicResult.extras["qcvars"] by
    converting flat arrays into numpy, shaped ones and floating strings.
    Unlike _qcvar_reshape_get/set, multipoles aren't compressed or plumped, only reshaped.

    Parameters
    ----------
    key
        Shape clue (usually QCVariable key) that includes (case insensitive) an identifier like
        'gradient' as a clue to the array's natural dimensions.
    val
        flat (?, ) list or scalar or string, probably from JSON storage.

    Returns
    -------
    float or numpy.ndarray
        Reshaped array of `val` with natural dimensions of `key`.

    """
    if isinstance(val, (np.ndarray, core.Matrix)):
        raise TypeError
    elif isinstance(val, list):
        tgt = np.asarray(val)
    else:
        # presumably scalar. may be string
        return float(val)

    if key.upper().startswith("MBIS"):
        if key.upper().endswith("CHARGES"):
            reshaper = (-1, )
        elif key.upper().endswith("DIPOLES"):
            reshaper = (-1, 3)
        elif key.upper().endswith("QUADRUPOLES"):
            reshaper = (-1, 3, 3)
        elif key.upper().endswith("OCTUPOLES"):
            reshaper = (-1, 3, 3, 3)
    elif key.upper().endswith("DIPOLE") or "DIPOLE -" in key.upper():
        reshaper = (3, )
    elif "QUADRUPOLE POLARIZABILITY TENSOR" in key.upper():
        reshaper = (3, 3, 3)
    elif any((key.upper().endswith(p) or f"{p} -" in key.upper()) for p in _multipole_order):
        p = [p for p in _multipole_order if (key.upper().endswith(p) or f"{p} -" in key.upper())]
        reshaper = tuple([3] * _multipole_order.index(p[0]))
    elif key.upper() in ["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "MULLIKEN CHARGES", "LOWDIN CHARGES"]:
        reshaper = (-1, )
    elif "GRADIENT" in key.upper():
        reshaper = (-1, 3)
    elif "HESSIAN" in key.upper():
        ndof = int(math.sqrt(len(tgt)))
        reshaper = (ndof, ndof)
    else:
        raise ValidationError(f'Uncertain how to reshape array: {key}')

    return tgt.reshape(reshaper)


_multipole_order = ["dummy", "dummy", "QUADRUPOLE", "OCTUPOLE", "HEXADECAPOLE"]
for order in range(5, 10):
    _multipole_order.append(f"{int(2**order)}-POLE")


def _qcvar_reshape_set(key: str, val: np.ndarray) -> np.ndarray:
    """Reverse :py:func:`_qcvar_reshape_get` for internal
    :py:class:`psi4.core.Matrix` storage.

    """
    reshaper = None
    if key.upper().startswith("MBIS"):
        if key.upper().endswith("CHARGES"):
            return val
        elif key.upper().endswith("DIPOLES"):
            reshaper = (-1, 3)
            return val.reshape(reshaper)
        elif key.upper().endswith("QUADRUPOLES"):
            val = val.reshape(-1, 3, 3)
            val = np.array([_multipole_compressor(val[iat], 2) for iat in range(len(val))])
            return val
        elif key.upper().endswith("OCTUPOLES"):
            val = val.reshape(-1, 3, 3, 3)
            val = np.array([_multipole_compressor(val[iat], 3) for iat in range(len(val))])
            return val
    elif key.upper().endswith("DIPOLE") or "DIPOLE -" in key.upper():
        reshaper = (1, 3)
    elif "QUADRUPOLE POLARIZABILITY TENSOR" in key.upper():
        reshaper = (3, 3, 3)
    elif any((key.upper().endswith(p) or f"{p} -" in key.upper()) for p in _multipole_order):
        p = [p for p in _multipole_order if (key.upper().endswith(p) or f"{p} -" in key.upper())]
        val = _multipole_compressor(val, _multipole_order.index(p[0]))
        reshaper = (1, -1)
    elif key.upper() in ["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "MULLIKEN CHARGES", "LOWDIN CHARGES"]:
        reshaper = (1, -1)

    if reshaper:
        return val.reshape(reshaper)
    else:
        return val


def _qcvar_reshape_get(key: str, val: core.Matrix) -> Union[core.Matrix, np.ndarray]:
    """For QCVariables where the 2D :py:class:`psi4.core.Matrix` shape is
    unnatural, convert to natural shape in :class:`numpy.ndarray`.

    """
    reshaper = None
    if key.upper().startswith("MBIS"):
        if key.upper().endswith("CHARGES"):
            return val.np
        elif key.upper().endswith("DIPOLES"):
            reshaper = (-1, 3)
            return val.np.reshape(reshaper)
        elif key.upper().endswith("QUADRUPOLES"):
            val = val.np.reshape(-1, 6)
            val = np.array([_multipole_plumper(val[iat], 2) for iat in range(len(val))])
            return val
        elif key.upper().endswith("OCTUPOLES"):
            val = val.np.reshape(-1, 10)
            val = np.array([_multipole_plumper(val[iat], 3) for iat in range(len(val))])
            return val
    elif key.upper().endswith("DIPOLE") or "DIPOLE -" in key.upper():
        reshaper = (3, )
    elif "QUADRUPOLE POLARIZABILITY TENSOR" in key.upper():
        reshaper = (3, 3, 3)
    elif any((key.upper().endswith(p) or f"{p} -" in key.upper()) for p in _multipole_order):
        p = [p for p in _multipole_order if (key.upper().endswith(p) or f"{p} -" in key.upper())]
        return _multipole_plumper(val.np.reshape((-1, )), _multipole_order.index(p[0]))
    elif key.upper() in ["MULLIKEN_CHARGES", "LOWDIN_CHARGES", "MULLIKEN CHARGES", "LOWDIN CHARGES"]:
        reshaper = (-1, )

    if reshaper:
        return val.np.reshape(reshaper)
    else:
        return val


def _multipole_compressor(complete: np.ndarray, order: int) -> np.ndarray:
    """Form flat unique components multipole array from complete Cartesian array.

    Parameters
    ----------
    order
        Multipole order. e.g., 1 for dipole, 4 for hexadecapole.
    complete
        Multipole array, order-dimensional Cartesian array expanded to complete components.

    Returns
    -------
    compressed : numpy.ndarray
        Multipole array, length (order + 1) * (order + 2) / 2 compressed to unique components.

    """
    compressed = []
    for ii in range(order + 1):
        lx = order - ii
        for lz in range(ii + 1):
            ly = ii - lz

            np_index = []
            for xval in range(lx):
                np_index.append(0)
            for yval in range(ly):
                np_index.append(1)
            for zval in range(lz):
                np_index.append(2)
            compressed.append(complete[tuple(np_index)])

    assert len(compressed) == ((order + 1) * (order + 2) / 2)
    return np.array(compressed)


def _multipole_plumper(compressed: np.ndarray, order: int) -> np.ndarray:
    """Form multidimensional multipole array from unique components array.

    Parameters
    ----------
    order
        Multipole order. e.g., 1 for dipole, 4 for hexadecapole.
    compressed
        Multipole array, length (order + 1) * (order + 2) / 2 compressed to unique components.

    Returns
    -------
    complete : numpy.ndarray
        Multipole array, order-dimensional Cartesian array expanded to complete components.

    """
    shape = tuple([3] * order)
    complete = np.zeros(shape)

    def compound_index(counter):
        # thanks, https://www.pamoc.it/tpc_cart_mom.html Eqn 2.2!
        # jn = nz + (ny + nz)(ny + nz + 1) / 2
        return int(
            counter.get("2", 0) + (counter.get("1", 0) + counter.get("2", 0)) *
            (counter.get("1", 0) + counter.get("2", 0) + 1) / 2)

    for idx in product("012", repeat=order):
        xyz_counts = Counter(idx)  # "010" --> {"0": 2, "1": 1}
        np_index = tuple(int(x) for x in idx)  # ('0', '1') --> (0, 1)

        complete[np_index] = compressed[compound_index(xyz_counts)]

    return complete


def _core_has_variable(key: str) -> bool:
    """Whether scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key*
    has been set in global memory.

    Parameters
    ----------
    key
        Case-insensitive key to global double or :py:class:`~psi4.core.Matrix`
        storage maps.

    """
    return core.has_scalar_variable(key) or core.has_array_variable(key)


def _core_wavefunction_has_variable(self: core.Wavefunction, key: str) -> bool:
    """Whether scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key*
    has been set on *self*.

    Parameters
    ----------
    self
        Wavefunction instance.
    key
        Case-insensitive key to instance's double or
        :py:class:`~psi4.core.Matrix` storage maps.

    """
    return self.has_scalar_variable(key) or self.has_array_variable(key)


def _core_variable(key: str) -> Union[float, core.Matrix, np.ndarray]:
    """Return copy of scalar or array :ref:`QCVariable <sec:appendices:qcvars>`
    *key* from global memory.

    Parameters
    ----------
    key
        Case-insensitive key to global double or :py:class:`~psi4.core.Matrix`
        storage maps.

    Returns
    -------
    float or ~numpy.ndarray or Matrix
        Requested QCVariable from global memory.

        - Scalar variables are returned as floats.
        - Array variables not naturally 2D (like multipoles or per-atom charges)
          are returned as :class:`~numpy.ndarray` of natural dimensionality.
        - Other array variables are returned as :py:class:`~psi4.core.Matrix` and
          may have an extra dimension with symmetry information.

    Raises
    ------
    KeyError
        If `key` not set on `self`.

    Example
    -------
    >>> psi4.gradient("hf/cc-pvdz")
    >>> psi4.variable("CURRENT ENERGY")
    -100.00985995185668
    >>> psi4.variable("CURRENT DIPOLE")
    array([ 0.        ,  0.        , -0.83217802])
    >>> psi4.variable("CURRENT GRADIENT")
    <psi4.core.Matrix object at 0x12d884fc0>
    >>> psi4.variable("CURRENT GRADIENT").np
    array([[ 6.16297582e-33,  6.16297582e-33, -9.41037138e-02],
           [-6.16297582e-33, -6.16297582e-33,  9.41037138e-02]])

    """
    key = _qcvar_warnings(key)

    if core.has_scalar_variable(key):
        return core.scalar_variable(key)
    elif core.has_array_variable(key):
        return _qcvar_reshape_get(key, core.array_variable(key))
    else:
        raise KeyError(f"psi4.core.variable: Requested variable '{key}' was not set!\n")


def _core_wavefunction_variable(self: core.Wavefunction, key: str) -> Union[float, core.Matrix, np.ndarray]:
    """Return copy of scalar or array :ref:`QCVariable <sec:appendices:qcvars>`
    *key* from *self*.

    Parameters
    ----------
    self
        Wavefunction instance.
    key
        Case-insensitive key to instance's double or :py:class:`~psi4.core.Matrix`
        storage maps.

    Returns
    -------
    float or ~numpy.ndarray or Matrix
        Requested QCVariable from `self`.

        - Scalar variables are returned as floats.
        - Array variables not naturally 2D (like multipoles or per-atom charges)
          are returned as :class:`~numpy.ndarray` of natural dimensionality.
        - Other array variables are returned as :py:class:`~psi4.core.Matrix` and
          may have an extra dimension with symmetry information.

    Raises
    ------
    KeyError
        If `key` not set on `self`.

    Example
    -------
    >>> g, wfn = psi4.gradient("hf/cc-pvdz", return_wfn=True)
    >>> wfn.variable("CURRENT ENERGY")
    -100.00985995185668
    >>> wfn.variable("CURRENT DIPOLE")
    array([ 0.        ,  0.        , -0.83217802])
    >>> wfn.variable("CURRENT GRADIENT")
    <psi4.core.Matrix object at 0x12d884fc0>
    >>> wfn.variable("CURRENT GRADIENT").np
    array([[ 6.16297582e-33,  6.16297582e-33, -9.41037138e-02],
           [-6.16297582e-33, -6.16297582e-33,  9.41037138e-02]])

    """
    key = _qcvar_warnings(key)

    if self.has_scalar_variable(key):
        return self.scalar_variable(key)
    elif self.has_array_variable(key):
        return _qcvar_reshape_get(key, self.array_variable(key))
    else:
        raise KeyError(f"psi4.core.Wavefunction.variable: Requested variable '{key}' was not set!\n")


def _core_set_variable(key: str, val: Union[core.Matrix, np.ndarray, float]) -> None:
    """Sets scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key* to
    *val* in global memory.

    Parameters
    ----------
    key
        Case-insensitive key to global double or :py:class:`~psi4.core.Matrix`
        storage maps.
    val
        Scalar or array to be stored in `key`. If :class:`~numpy.ndarray` and
        data `key` does not naturally fit in 2D Matrix (often charge and
        multipole QCVariables), it will be reshaped, as all
        :class:`~numpy.ndarray` are stored as :class:`~psi4.core.Matrix`.

    Raises
    ------
    ValidationError
        If `val` is a scalar but `key` already exists as an array variable. Or
        if `val` is an array but `key` already exists as a scalar variable.

    """
    if isinstance(val, core.Matrix):
        if core.has_scalar_variable(key):
            raise ValidationError(f"psi4.core.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            core.set_array_variable(key, val)
    elif isinstance(val, np.ndarray):
        if core.has_scalar_variable(key):
            raise ValidationError(f"psi4.core.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            core.set_array_variable(key, core.Matrix.from_array(_qcvar_reshape_set(key, val)))
    else:
        if core.has_array_variable(key):
            raise ValidationError(f"psi4.core.set_variable: Target variable '{key}' already an array variable!")
        else:
            core.set_scalar_variable(key, val)

    # TODO _qcvar_warnings(key)


def _core_wavefunction_set_variable(self: core.Wavefunction, key: str, val: Union[core.Matrix, np.ndarray, float]) -> None:
    """Sets scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key* to
    *val* on *self*.

    Parameters
    ----------
    self
        Wavefunction instance.
    key
        Case-insensitive key to instance's double or :class:`~psi4.core.Matrix`
        storage maps.

        - If ``CURRENT ENERGY``, syncs with ``self.energy_``.
        - If ``CURRENT GRADIENT``, syncs with ``gradient_``.
        - If ``CURRENT HESSIAN``, syncs with ``self.hessian_``.
    val
        Scalar or array to be stored in `key`. If :class:`~numpy.ndarray` and
        data `key` does not naturally fit in 2D Matrix (often charge and
        multipole QCVariables), it will be reshaped, as all
        :class:`~numpy.ndarray` are stored as :class:`~psi4.core.Matrix`.

    Raises
    ------
    ~psi4.driver.ValidationError
        If `val` is a scalar but `key` already exists as an array variable. Or
        if `val` is an array but `key` already exists as a scalar variable.

    """
    if isinstance(val, core.Matrix):
        if self.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            self.set_array_variable(key, val)
    elif isinstance(val, np.ndarray):
        if self.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            self.set_array_variable(key, core.Matrix.from_array(_qcvar_reshape_set(key, val)))
    else:
        if self.has_array_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already an array variable!")
        else:
            self.set_scalar_variable(key, val)

    # TODO _qcvar_warnings(key)


def _core_del_variable(key: str) -> None:
    """Removes scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key*
    from global memory if present.

    Parameters
    ----------
    key
        Case-insensitive key to global double or :py:class:`~psi4.core.Matrix`
        storage maps.

    """
    if core.has_scalar_variable(key):
        core.del_scalar_variable(key)
    elif core.has_array_variable(key):
        core.del_array_variable(key)


def _core_wavefunction_del_variable(self: core.Wavefunction, key: str) -> None:
    """Removes scalar or array :ref:`QCVariable <sec:appendices:qcvars>` *key*
    from *self* if present.

    Parameters
    ----------
    self
        Wavefunction instance.
    key
        Case-insensitive key to instance's double or
        :py:class:`~psi4.core.Matrix` storage maps.

    """
    if self.has_scalar_variable(key):
        self.del_scalar_variable(key)
    elif self.has_array_variable(key):
        self.del_array_variable(key)


def _core_variables(include_deprecated_keys: bool = False) -> Dict[str, Union[float, core.Matrix, np.ndarray]]:
    """Return all scalar or array :ref:`QCVariables <sec:appendices:qcvars>`
    from global memory.

    Parameters
    ----------
    include_deprecated_keys
        Also return duplicate entries with keys that have been deprecated.

    Returns
    -------
    ~typing.Dict[str, ~typing.Union[float, ~numpy.ndarray, Matrix]
        Map of all QCVariables that have been set.

        - Scalar variables are returned as floats.
        - Array variables not naturally 2D (like multipoles or per-atom charges)
          are returned as :class:`~numpy.ndarray` of natural dimensionality.
        - Other array variables are returned as :py:class:`~psi4.core.Matrix` and
          may have an extra dimension with symmetry information.

    """
    dicary = {**core.scalar_variables(), **{k: _qcvar_reshape_get(k, v) for k, v in core.array_variables().items()}}

    if include_deprecated_keys:
        for old_key, (current_key, version) in _qcvar_transitions.items():
            if current_key in dicary:
                dicary[old_key] = dicary[current_key]

    return dicary


def _core_wavefunction_variables(self, include_deprecated_keys: bool = False) -> Dict[str, Union[float, core.Matrix, np.ndarray]]:
    """Return all scalar or array :ref:`QCVariables <sec:appendices:qcvars>`
    from *self*.

    Parameters
    ----------
    self
        Wavefunction instance.
    include_deprecated_keys
        Also return duplicate entries with keys that have been deprecated.

    Returns
    -------
    ~typing.Dict[str, ~typing.Union[float, ~numpy.ndarray, Matrix]
        Map of all QCVariables that have been set on `self`.

        - Scalar variables are returned as floats.
        - Array variables not naturally 2D (like multipoles or per-atom charges)
          are returned as :class:`~numpy.ndarray` of natural dimensionality.
        - Other array variables are returned as :py:class:`~psi4.core.Matrix` and
          may have an extra dimension with symmetry information.

    """
    dicary = {**self.scalar_variables(), **{k: _qcvar_reshape_get(k, v) for k, v in self.array_variables().items()}}

    if include_deprecated_keys:
        for old_key, (current_key, version) in _qcvar_transitions.items():
            if current_key in dicary:
                dicary[old_key] = dicary[current_key]

    return dicary


core.has_variable = _core_has_variable
core.variable = _core_variable
core.set_variable = _core_set_variable
core.del_variable = _core_del_variable
core.variables = _core_variables

core.Wavefunction.has_variable = _core_wavefunction_has_variable
core.Wavefunction.variable = _core_wavefunction_variable
core.Wavefunction.set_variable = _core_wavefunction_set_variable
core.Wavefunction.del_variable = _core_wavefunction_del_variable
core.Wavefunction.variables = _core_wavefunction_variables

## Psi4 v1.4 Export Deprecations


def _core_get_variable(key):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.variable` instead.

    """
    warnings.warn(
        "Using `psi4.core.get_variable` instead of `psi4.core.variable` (or `psi4.core.scalar_variable` for scalar variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.scalar_variable(key)


def _core_get_variables():
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.variables` instead.

    """
    warnings.warn(
        "Using `psi4.core.get_variables` instead of `psi4.core.variables` (or `psi4.core.scalar_variables` for scalar variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.scalar_variables()


def _core_get_array_variable(key):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.variable` instead.

    """
    warnings.warn(
        "Using `psi4.core.get_array_variable` instead of `psi4.core.variable` (or `psi4.core.array_variable` for array variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.array_variable(key)


def _core_get_array_variables():
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.variables` instead.

    """
    warnings.warn(
        "Using `psi4.core.get_array_variables` instead of `psi4.core.variables` (or `psi4.core.array_variables` for array variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.array_variables()


core.get_variable = _core_get_variable
core.get_variables = _core_get_variables
core.get_array_variable = _core_get_array_variable
core.get_array_variables = _core_get_array_variables


def _core_wavefunction_get_variable(cls, key):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.Wavefunction.variable` instead.

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.get_variable` instead of `psi4.core.Wavefunction.variable` (or `psi4.core.Wavefunction.scalar_variable` for scalar variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.scalar_variable(key)


def _core_wavefunction_get_array(cls, key):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.Wavefunction.variable` instead.

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.get_array` instead of `psi4.core.Wavefunction.variable` (or `psi4.core.Wavefunction.array_variable` for array variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.array_variable(key)


def _core_wavefunction_set_array(cls, key, val):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.Wavefunction.set_variable` instead.

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.set_array` instead of `psi4.core.Wavefunction.set_variable` (or `psi4.core.Wavefunction.set_array_variable` for array variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.set_array_variable(key, val)


def _core_wavefunction_arrays(cls):
    """
    .. deprecated:: 1.4
       Use :py:func:`psi4.core.Wavefunction.variables` instead.

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.arrays` instead of `psi4.core.Wavefunction.variables` (or `psi4.core.Wavefunction.array_variables` for array variables only) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.array_variables()


core.Wavefunction.get_variable = _core_wavefunction_get_variable
core.Wavefunction.get_array = _core_wavefunction_get_array
core.Wavefunction.set_array = _core_wavefunction_set_array
core.Wavefunction.arrays = _core_wavefunction_arrays


def _core_wavefunction_frequencies(self):
    """Returns the results of a frequency analysis.

    Parameters
    ----------
    self
        Wavefunction instance.

    Returns
    -------
    ~typing.Optional[~typing.Dict[str, ~numpy.ndarray]]
        A dictionary of vibrational information. See :py:func:`psi4.driver.qcdb.vib.harmonic_analysis`

    """
    if not hasattr(self, 'frequency_analysis'):
        return None

    vibinfo = self.frequency_analysis
    vibonly = qcdb.vib.filter_nonvib(vibinfo)
    return core.Vector.from_array(qcdb.vib.filter_omega_to_real(vibonly['omega'].data))


core.Wavefunction.frequencies = _core_wavefunction_frequencies


def _core_doublet(A, B, transA, transB):
    """Multiply two matrices together.

    .. deprecated:: 1.4
       Use :py:func:`psi4.core.doublet` instead.

    """
    warnings.warn(
        "Using `psi4.core.Matrix.doublet` instead of `psi4.core.doublet` is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.doublet(A, B, transA, transB)


def _core_triplet(A, B, C, transA, transB, transC):
    """Multiply three matrices together.

    .. deprecated:: 1.4
       Use :py:func:`psi4.core.triplet` instead.

    """
    warnings.warn(
        "Using `psi4.core.Matrix.triplet` instead of `psi4.core.triplet` is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.triplet(A, B, C, transA, transB, transC)


core.Matrix.doublet = staticmethod(_core_doublet)
core.Matrix.triplet = staticmethod(_core_triplet)


@staticmethod
def _core_erisieve_build(
        orbital_basis: core.BasisSet,
        cutoff: float = 0.0,
        do_csam: bool = False
    ) -> core.ERISieve:
    """
    This function previously constructed a Psi4 ERISieve object from an input basis set, with an optional cutoff threshold for
    ERI screening and an optional input to enable CSAM screening (over Schwarz screening).

    However, as the ERISieve class was removed from Psi4 in v1.9, the function now throws with an UpgradeHelper
    exception, and lets the user know to use TwoBodyAOInt instead.

    Parameters
    ----------
    orbital_basis
        Basis set to use in the ERISieve object.
    cutoff
        Integral cutoff threshold to use for Schwarz/CSAM screening. Defaults to 0.0, disabling screening entirely.
    do_csam
        Use CSAM screening? If True, CSAM screening is used; else, Schwarz screening is used. By default,
        Schwarz screening is utilized.

    Returns
    -------
    ERISieve
        Initialized ERISieve object.

    Example
    -------
    >>> sieve = psi4.core.ERISieve.build(bas, cutoff, csam)
    """

    raise UpgradeHelper("ERISieve", "TwoBodyAOInt", 1.8, " The ERISieve class has been removed and replaced with the TwoBodyAOInt class. ERISieve.build(orbital_basis, cutoff, do_csam) can be replaced with the command sequence factory = psi4.core.IntegralFactory(basis); factory.eri(0).")


core.ERISieve.build = _core_erisieve_build
