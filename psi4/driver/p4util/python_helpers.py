#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import numpy as np

import qcelemental as qcel
from psi4 import core
from psi4.driver import qcdb

from . import optproc
from .exceptions import TestComparisonError, ValidationError

## Python basis helps


@staticmethod
def _pybuild_basis(mol,
                   key=None,
                   target=None,
                   fitrole='ORBITAL',
                   other=None,
                   puream=-1,
                   return_atomlist=False,
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

    bs, basisdict = qcdb.BasisSet.pyconstruct(
        mol.to_dict(), key, resolved_target, fitrole, other, return_dict=True, return_atomlist=return_atomlist)

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
def _core_wavefunction_build(mol, basis=None):
    if basis is None:
        basis = core.BasisSet.build(mol)
    elif isinstance(basis, str):
        basis = core.BasisSet.build(mol, "ORBITAL", basis)

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
def _core_wavefunction_from_file(wfn_data):
    """Summary

    Parameters
    ----------
    wfn_data : str or dict
        If a str reads a Wavefunction from a disk otherwise, assumes the data
        is passed in.

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
        wfn_data = np.load(wfn_data).item()
    else:
        # Could be path-like or file-like, let `np.load` handle it
        wfn_data = np.load(wfn_data).item()

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


def _core_wavefunction_to_file(wfn, filename=None):
    """Converts a Wavefunction object to a base class

    Parameters
    ----------
    wfn : Wavefunction
        A Wavefunction or inherited class
    filename : None, optional
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
            'X':        wfn.X().to_array()        if wfn.X()        else None,
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
        np.save(filename, wfn_data)

    return wfn_data


core.Wavefunction.to_file = _core_wavefunction_to_file

## Python JK helps


@staticmethod
def _core_jk_build(orbital_basis, aux=None, jk_type=None, do_wK=None, memory=None):
    """
    Constructs a Psi4 JK object from an input basis.

    Parameters
    ----------
    orbital_basis : :py:class:`~psi4.core.BasisSet`
        Orbital basis to use in the JK object.
    aux : :py:class:`~psi4.core.BasisSet`, optional
        Optional auxiliary basis set for density-fitted tensors. Defaults
        to the DF_BASIS_SCF if set, otherwise the correspond JKFIT basis
        to the passed in `orbital_basis`.
    jk_type : str, optional
        Type of JK object to build (DF, Direct, PK, etc). Defaults to the
        current global SCF_TYPE option.

    Returns
    -------
    :py:class:`~psi4.core.JK`
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


def set_options(options_dict):
    """
    Sets Psi4 global options from an input dictionary.
    """

    for k, v, in options_dict.items():
        core.set_global_option(k.upper(), v)


def set_module_options(module, options_dict):
    """
    Sets Psi4 module options from a module specification and input dictionary.
    """

    for k, v, in options_dict.items():
        core.set_local_option(module.upper(), k.upper(), v)


## OEProp helpers


def pcm_helper(block):
    """
    Passes multiline string *block* to PCMSolver parser.

    Parameters
    ----------
    block: multiline string with PCM input in PCMSolver syntax.
    """

    suffix = str(os.getpid()) + '.' + str(uuid.uuid4())[:8]
    pcmsolver_fname = 'pcmsolver.' + suffix + '.inp'
    with open(pcmsolver_fname, 'w') as handle:
        handle.write(block)
    import pcmsolver
    parsed_pcm = pcmsolver.parse_pcm_input(pcmsolver_fname)
    os.remove(pcmsolver_fname)
    pcmsolver_parsed_fname = '@pcmsolver.' + suffix
    with open(pcmsolver_parsed_fname, 'w') as tmp:
        tmp.write(parsed_pcm)
    core.set_local_option('PCM', 'PCMSOLVER_PARSED_FNAME', '{}'.format(pcmsolver_parsed_fname))


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
    'MAYER_INDICES', 'MO_EXTENTS', 'GRID_FIELD', 'GRID_ESP', 'ESP_AT_NUCLEI', 'NO_OCCUPATIONS'
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


def _core_has_variable(key):
    return core.has_scalar_variable(key) or core.has_array_variable(key)


def _core_wavefunction_has_variable(cls, key):
    return cls.has_scalar_variable(key) or cls.has_array_variable(key)


def _core_variable(key):
    if core.has_scalar_variable(key):
        return core.scalar_variable(key)
    elif core.has_array_variable(key):
        return core.array_variable(key)
    else:
        raise KeyError("psi4.core.variable: Requested variable " + key + " was not set!\n")


def _core_wavefunction_variable(cls, key):
    if cls.has_scalar_variable(key):
        return cls.scalar_variable(key)
    elif cls.has_array_variable(key):
        return cls.array_variable(key)
    else:
        raise KeyError("psi4.core.Wavefunction.variable: Requested variable " + key + " was not set!\n")


def _core_set_variable(key, val):
    if isinstance(val, core.Matrix):
        if core.has_scalar_variable(key):
            raise ValidationError("psi4.core.set_variable: Target variable " + key + " already a scalar variable!")
        else:
            core.set_array_variable(key, val)
    elif isinstance(val, np.ndarray):
        if core.has_scalar_variable(key):
            raise ValidationError("psi4.core.set_variable: Target variable " + key + " already a scalar variable!")
        else:
            core.set_array_variable(key, core.Matrix.from_array(val))
    else:
        if core.has_array_variable(key):
            raise ValidationError("psi4.core.set_variable: Target variable " + key + " already an array variable!")
        else:
            core.set_scalar_variable(key, val)


def _core_wavefunction_set_variable(cls, key, val):
    if isinstance(val, core.Matrix):
        if cls.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable " + key +
                                  " already a scalar variable!")
        else:
            cls.set_array_variable(key, val)
    elif isinstance(val, np.ndarray):
        if cls.has_scalar_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable " + key +
                                  " already a scalar variable!")
        else:
            cls.set_array_variable(key, core.Matrix.from_array(val))
    else:
        if cls.has_array_variable(key):
            raise ValidationError("psi4.core.Wavefunction.set_variable: Target variable " + key +
                                  " already an array variable!")
        else:
            cls.set_scalar_variable(key, val)


def _core_del_variable(key):
    if core.has_scalar_variable(key):
        core.del_scalar_variable(key)
    elif core.has_array_variable(key):
        core.del_array_variable(key)


def _core_wavefunction_del_variable(cls, key):
    if cls.has_scalar_variable(key):
        cls.del_scalar_variable(key)
    elif cls.has_array_variable(key):
        cls.del_array_variable(key)


def _core_variables():
    return {**core.scalar_variables(), **core.array_variables()}


def _core_wavefunction_variables(cls):
    return {**cls.scalar_variables(), **cls.array_variables()}


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
    warnings.warn(
        "Using `psi4.core.get_variable` instead of `psi4.core.variable` (or `psi4.core.scalar_variable` for scalar variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.scalar_variable(key)


def _core_get_variables():
    warnings.warn(
        "Using `psi4.core.get_variables` instead of `psi4.core.variables` (or `psi4.core.scalar_variables` for scalar variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.scalar_variables()


def _core_get_array_variable(key):
    warnings.warn(
        "Using `psi4.core.get_array_variable` instead of `psi4.core.variable` (or `psi4.core.array_variable` for array variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.array_variable(key)


def _core_get_array_variables():
    warnings.warn(
        "Using `psi4.core.get_array_variables` instead of `psi4.core.variables` (or `psi4.core.array_variables` for array variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.array_variables()


core.get_variable = _core_get_variable
core.get_variables = _core_get_variables
core.get_array_variable = _core_get_array_variable
core.get_array_variables = _core_get_array_variables


def _core_wavefunction_get_variable(cls, key):
    warnings.warn(
        "Using `psi4.core.Wavefunction.get_variable` instead of `psi4.core.Wavefunction.variable` (or `psi4.core.Wavefunction.scalar_variable` for scalar variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.scalar_variable(key)


def _core_wavefunction_get_array(cls, key):
    warnings.warn(
        "Using `psi4.core.Wavefunction.get_array` instead of `psi4.core.Wavefunction.variable` (or `psi4.core.Wavefunction.array_variable` for array variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.array_variable(key)


def _core_wavefunction_set_array(cls, key, val):
    warnings.warn(
        "Using `psi4.core.Wavefunction.set_array` instead of `psi4.core.Wavefunction.set_variable` (or `psi4.core.Wavefunction.set_array_variable` for array variables only) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.set_array_variable(key, val)


def _core_wavefunction_arrays(cls):
    warnings.warn(
        "Using `psi4.core.Wavefunction.arrays` instead of `psi4.core.Wavefunction.variables` (or `psi4.core.Wavefunction.array_variables` for array variables only) is deprecated, and in 1.4 it will stop working\n",
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
    warnings.warn(
        "Using `psi4.core.Wavefunction.legacy_frequencies` (accessing c-side member data) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.legacy_frequencies()


def _core_wavefunction_set_frequencies(cls, val):
    warnings.warn(
        "Using `psi4.core.Wavefunction.set_frequencies` (accessing c-side member data) instead of `psi4.core.Wavefunction.frequency_analysis` (py-side member data) is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return cls.set_legacy_frequencies(val)


core.Wavefunction.frequencies = _core_wavefunction_frequencies
core.Wavefunction.legacy_frequencies = _core_wavefunction_legacy_frequencies
core.Wavefunction.set_frequencies = _core_wavefunction_set_frequencies

## Psi4 v1.3 Export Deprecations


def _core_get_gradient():
    warnings.warn(
        "Using `psi4.core.get_gradient` (only used internally for C++ optking; deprecated silently in 1.2) is deprecated, and in 1.4 (or whenever Py optking is adopted) it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.get_legacy_gradient()


def _core_set_gradient(val):
    warnings.warn(
        "Using `psi4.core.set_gradient` (only used internally for C++ optking; deprecated silently in 1.2) is deprecated, and in 1.4 (or whenever Py optking is adopted) it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.set_legacy_gradient(val)


core.get_gradient = _core_get_gradient
core.set_gradient = _core_set_gradient


def _core_doublet(A, B, transA, transB):
    warnings.warn(
        "Using `psi4.core.Matrix.doublet` instead of `psi4.core.doublet` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.doublet(A, B, transA, transB)


def _core_triplet(A, B, C, transA, transB, transC):
    warnings.warn(
        "Using `psi4.core.Matrix.triplet` instead of `psi4.core.triplet` is deprecated, and in 1.4 it will stop working\n",
        category=FutureWarning,
        stacklevel=2)
    return core.triplet(A, B, C, transA, transB, transC)


core.Matrix.doublet = staticmethod(_core_doublet)
core.Matrix.triplet = staticmethod(_core_triplet)
