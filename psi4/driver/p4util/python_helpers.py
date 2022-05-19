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

import os
import re
import sys
import uuid
import warnings
from collections import Counter
from itertools import product
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, Union

import numpy as np

import qcelemental as qcel
from psi4 import core
from psi4.driver import qcdb

from . import optproc
from .exceptions import TestComparisonError, ValidationError, UpgradeHelper

## Python basis helps


@staticmethod
def _pybuild_basis(mol,
                   key=None,
                   target=None,
                   fitrole='ORBITAL',
                   other=None,
                   puream=-1,
                   return_atomlist=False,
                   *,
                   quiet=False):
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
def _core_wavefunction_build(mol, basis=None, *, quiet: bool = False):
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


def _core_wavefunction_get_scratch_filename(self, filenumber):
    """ Given a wavefunction and a scratch file number, canonicalizes the name
        so that files can be consistently written and read """
    fname = os.path.split(os.path.abspath(core.get_writer_file_prefix(self.molecule().name())))[1]
    psi_scratch = core.IOManager.shared_object().get_default_path()
    return os.path.join(psi_scratch, fname + '.' + str(filenumber))


core.Wavefunction.get_scratch_filename = _core_wavefunction_get_scratch_filename


@staticmethod
def _core_wavefunction_from_file(wfn_data: Union[str, Dict, Path]) -> core.Wavefunction:
    r"""Build Wavefunction from data.

    Parameters
    ----------
    wfn_data
        If a dict, use data directly. Otherwise, path-like passed to :py:func:`numpy.load`
        to read from disk.

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


def _core_wavefunction_to_file(wfn: core.Wavefunction, filename: str = None) -> Dict:
    """Converts a Wavefunction object to a base class

    Parameters
    ----------
    wfn
        A Wavefunction or inherited class
    filename
        An optional filename to write the data to

    Returns
    -------
    dict
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
            'density_fitted': wfn.density_fitted(),
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
def _core_jk_build(orbital_basis: core.BasisSet, aux: core.BasisSet = None, jk_type: str = None, do_wK: bool = None, memory: int = None) -> core.JK:
    """
    Constructs a Psi4 JK object from an input basis.

    Parameters
    ----------
    orbital_basis
        Orbital basis to use in the JK object.
    aux
        Optional auxiliary basis set for density-fitted tensors. Defaults
        to the DF_BASIS_SCF if set, otherwise the correspond JKFIT basis
        to the passed in `orbital_basis`.
    jk_type
        Type of JK object to build (DF, Direct, PK, etc). Defaults to the
        current global SCF_TYPE option.

    Returns
    -------
    JK
        Uninitialized JK object.

    Example
    -------

    jk = psi4.core.JK.build(bas)
    jk.set_memory(int(5e8)) # 4GB of memory
    jk.initialize()

    ...

    jk.C_left_add(matirx)
    jk.compute()
    jk.C_clear()

    ...

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


def _core_vbase_get_np_xyzw(Vpot):
    """
    Returns the x, y, z, and weights of a grid as a tuple of NumPy array objects.
    """
    x_list = []
    y_list = []
    z_list = []
    w_list = []

    # Loop over every block in the potenital
    for b in range(Vpot.nblocks()):

        # Obtain the block
        block = Vpot.get_block(b)

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


def set_options(options_dict: Dict[str, Any], verbose: int = 1) -> None:
    """Sets Psi4 options from an input dictionary.

    Parameters
    ----------
    options_dict
        Dictionary where keys are "option_name" for global options or
        "module_name__option_name" (double underscore separation) for
        option local to "module_name". Values are the option value. All
        are case insensitive.
    verbose
        Control print volume.

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
    """
    warnings.warn(
        "Using `psi4.set_module_options(<module>, {<key>: <val>})` instead of `psi4.set_options({<module>__<key>: <val>})` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)

    for k, v, in options_dict.items():
        core.set_local_option(module.upper(), k.upper(), v)


## OEProp helpers


def pcm_helper(block: str):
    """
    Passes multiline string *block* to PCMSolver parser.

    Parameters
    ----------
    block
        multiline string with PCM input in PCMSolver syntax.
    """
    import pcmsolver

    with NamedTemporaryFile(mode="w+t", delete=True) as fl:
        fl.write(block)
        fl.flush()
        parsed_pcm = pcmsolver.parse_pcm_input(fl.name)

    with NamedTemporaryFile(mode="w+t", delete=False) as fl:
        fl.write(parsed_pcm)
        core.set_local_option("PCM", "PCMSOLVER_PARSED_FNAME", fl.name)


def basname(name):
    """Imitates BasisSet.make_filename() without the gbs extension"""
    return name.lower().replace('+', 'p').replace('*', 's').replace('(', '_').replace(')', '_').replace(',', '_')


def temp_circular_import_blocker():
    pass


def basis_helper(block, name='', key='BASIS', set_option=True):
    """For PsiAPI mode, forms a basis specification function from *block*
    and associates it with keyword *key* under handle *name*. Registers
    the basis spec with Psi4 so that it can be applied again to future
    molecules. For usage, see mints2, mints9, and cc54 test cases. Unless
    *set_option* is False, *name* will be set as current active *key*,
    equivalent to `set key name` or `set_option({key: name})`.

    """
    key = key.upper()
    name = ('anonymous' + str(uuid.uuid4())[:8]) if name == '' else name
    cleanbas = basname(name).replace('-', '')  # further remove hyphens so can be function name
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
                basstrings[basname(name)] = basblock[0]
            else:
                message = (
                    "Conflicting basis set specification: assign lines present but shells have no [basname] label."
                    "")
                raise TestComparisonError(message)
        else:
            # case with specs separated by [basname] markers
            for idx in range(0, len(basblock), 2):
                basstrings[basname(basblock[idx])] = basblock[idx + 1]

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


_multipole_order = ["dummy", "dummy", "QUADRUPOLE", "OCTUPOLE", "HEXADECAPOLE"]
for order in range(5, 10):
    _multipole_order.append(f"{int(2**order)}-POLE")


def _qcvar_reshape_set(key, val):
    """Reverse `_qcvar_reshape_get` for internal psi4.core.Matrix storage."""

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


def _qcvar_reshape_get(key, val):
    """For QCVariables where the 2D psi4.core.Matrix shape is unnatural, convert to natural shape in ndarray."""

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

def _multipole_compressor(complete, order):
    """Form flat unique components multipole array from complete Cartesian array.

    Parameters
    ----------
    order : int
        Multipole order. e.g., 1 for dipole, 4 for hexadecapole.
    complete : ndarray
        Multipole array, order-dimensional Cartesian array expanded to complete components.

    Returns
    -------
    compressed : ndarray
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
    """Whether scalar or array QCVariable *key* has been set in global memory."""

    return core.has_scalar_variable(key) or core.has_array_variable(key)


def _core_wavefunction_has_variable(cls: core.Wavefunction, key: str) -> bool:
    """Whether scalar or array QCVariable *key* has been set on *self* :class:`psi4.core.Wavefunction`."""

    return cls.has_scalar_variable(key) or cls.has_array_variable(key)


def _core_variable(key: str) -> Union[float, core.Matrix, np.ndarray]:
    """Return copy of scalar or array QCVariable *key* from global memory.

    Returns
    -------
    float or numpy.ndarray or Matrix
        Scalar variables are returned as floats.
        Array variables not naturally 2D (like multipoles) are returned as :class:`numpy.ndarray` of natural dimensionality.
        Other array variables are returned as :py:class:`~psi4.core.Matrix` and may have an extra dimension with symmetry information.

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


def _core_wavefunction_variable(cls: core.Wavefunction, key: str) -> Union[float, core.Matrix, np.ndarray]:
    """Return copy of scalar or array QCVariable *key* from *self* :class:`psi4.core.Wavefunction`.

    Returns
    -------
    float or numpy.ndarray or Matrix
        Scalar variables are returned as floats.
        Array variables not naturally 2D (like multipoles) are returned as :class:`numpy.ndarray` of natural dimensionality.
        Other array variables are returned as :py:class:`~psi4.core.Matrix` and may have an extra dimension with symmetry information.

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

    if cls.has_scalar_variable(key):
        return cls.scalar_variable(key)
    elif cls.has_array_variable(key):
        return _qcvar_reshape_get(key, cls.array_variable(key))
    else:
        raise KeyError(f"psi4.core.Wavefunction.variable: Requested variable '{key}' was not set!\n")


def _core_set_variable(key: str, val: Union[core.Matrix, np.ndarray, float]) -> None:
    """Sets scalar or array QCVariable *key* to *val* in global memory."""

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


def _core_wavefunction_set_variable(cls: core.Wavefunction, key: str, val: Union[core.Matrix, np.ndarray, float]) -> None:
    """Sets scalar or array QCVariable *key* to *val* on *cls*."""

    if isinstance(val, core.Matrix):
        if cls.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            cls.set_array_variable(key, val)
    elif isinstance(val, np.ndarray):
        if cls.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already a scalar variable!")
        else:
            cls.set_array_variable(key, core.Matrix.from_array(_qcvar_reshape_set(key, val)))
    else:
        if cls.has_array_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable '{key}' already an array variable!")
        else:
            cls.set_scalar_variable(key, val)

    # TODO _qcvar_warnings(key)


def _core_del_variable(key: str) -> None:
    """Removes scalar or array QCVariable *key* from global memory if present."""

    if core.has_scalar_variable(key):
        core.del_scalar_variable(key)
    elif core.has_array_variable(key):
        core.del_array_variable(key)


def _core_wavefunction_del_variable(cls: core.Wavefunction, key: str) -> None:
    """Removes scalar or array QCVariable *key* from *cls* if present."""

    if cls.has_scalar_variable(key):
        cls.del_scalar_variable(key)
    elif cls.has_array_variable(key):
        cls.del_array_variable(key)


def _core_variables(include_deprecated_keys: bool = False) -> Dict[str, Union[float, core.Matrix, np.ndarray]]:
    """Return all scalar or array QCVariables from global memory."""

    dicary = {**core.scalar_variables(), **{k: _qcvar_reshape_get(k, v) for k, v in core.array_variables().items()}}

    if include_deprecated_keys:
        for old_key, current_key in _qcvar_transitions.items():
            if current_key in dicary:
                dicary[old_key] = dicary[current_key]

    return dicary


def _core_wavefunction_variables(cls, include_deprecated_keys: bool = False) -> Dict[str, Union[float, core.Matrix, np.ndarray]]:
    """Return all scalar or array QCVariables from *cls*."""

    dicary = {**cls.scalar_variables(), **{k: _qcvar_reshape_get(k, v) for k, v in cls.array_variables().items()}}

    if include_deprecated_keys:
        for old_key, current_key in _qcvar_transitions.items():
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
       Use :py:func:`psi4.variable` instead.

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
       Use :py:func:`psi4.variable` instead.

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


def _core_wavefunction_frequencies(cls):
    if not hasattr(cls, 'frequency_analysis'):
        return None

    vibinfo = cls.frequency_analysis
    vibonly = qcdb.vib.filter_nonvib(vibinfo)
    return core.Vector.from_array(qcdb.vib.filter_omega_to_real(vibonly['omega'].data))


def _core_wavefunction_legacy_frequencies(cls):
    """
    .. deprecated:: 1.4

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.legacy_frequencies` (accessing c-side member data) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.legacy_frequencies()


def _core_wavefunction_set_frequencies(cls, val):
    """
    .. deprecated:: 1.4

    """
    warnings.warn(
        "Using `psi4.core.Wavefunction.set_frequencies` (accessing c-side member data) instead of `psi4.core.Wavefunction.frequency_analysis` (py-side member data) is deprecated, and as soon as 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.set_legacy_frequencies(val)


core.Wavefunction.frequencies = _core_wavefunction_frequencies
core.Wavefunction.legacy_frequencies = _core_wavefunction_legacy_frequencies
core.Wavefunction.set_frequencies = _core_wavefunction_set_frequencies


def _core_wavefunction_X(cls):
    warnings.warn(
        "Using `psi4.core.Wavefunction.X` instead of `psi4.core.Wavefunction.lagrangian` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.lagrangian()


core.Wavefunction.X = _core_wavefunction_X

## Psi4 v1.3 Export Deprecations


def _core_get_gradient():
    """
    .. deprecated:: 1.2

    """
    warnings.warn(
        "Using `psi4.core.get_gradient` (only used internally for C++ optking; deprecated silently in 1.2) is deprecated, and as soon as 1.6 (or whenever Py optking is adopted) it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.get_legacy_gradient()


def _core_set_gradient(val):
    """
    .. deprecated:: 1.2

    """
    warnings.warn(
        "Using `psi4.core.set_gradient` (only used internally for C++ optking; deprecated silently in 1.2) is deprecated, and as soon as 1.6 (or whenever Py optking is adopted) it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.set_legacy_gradient(val)


core.get_gradient = _core_get_gradient
core.set_gradient = _core_set_gradient


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
