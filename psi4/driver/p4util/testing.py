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
"""Module with comparison functions with output configured for Psi4."""

__all__ = [
    "compare",
    "compare_arrays",
    "compare_cubes",
    "compare_integers",
    "compare_matrices",
    "compare_molrecs",
    "compare_recursive",
    "compare_strings",
    "compare_values",
    "compare_vectors",
    "compare_wavefunctions",
]

import sys
from functools import partial

import numpy as np
import qcelemental as qcel

from psi4 import core
from psi4.driver import qcdb
from .exceptions import TestComparisonError

# TODO in multistage compare_* fns, we're potentially stopping the fn prematurely and not in the manner of the handling fn.


def _mergedapis_compare_cubes(expected, computed, *args, **kwargs):
    """Shim allowing Psi4-style or QCA-style testing interfaces for cube files."""

    qcdb.testing._merge_psi4_qcel_apis(args, kwargs)

    # Grab grid points. Skip the first nine lines and the last one
    expected = np.genfromtxt(expected, skip_header=9, skip_footer=1)
    computed = np.genfromtxt(computed, skip_header=9, skip_footer=1)

    return qcel.testing.compare_values(expected, computed, **kwargs)


def _mergedapis_compare_vectors(expected, computed, *args, **kwargs):
    """Shim allowing Psi4-style or QCA-style testing interfaces for :py:class:`psi4.core.Vector`."""

    qcdb.testing._merge_psi4_qcel_apis(args, kwargs)

    if kwargs.pop("check_name", False):
        compare(expected.name, computed.name, f'{expected.name} vs. {computed.name} name', quiet=True)

    compare(expected.nirrep(), computed.nirrep(), f'{expected.name} vs. {computed.name} irreps', quiet=True)
    for irrep in range(expected.nirrep()):
        compare(
            expected.dim(irrep),
            computed.dim(irrep),
            f'{expected.name} vs. {computed.name} irrep {irrep} dimensions ({expected.dim(irrep)} vs {computed.dim(irrep)})',
            quiet=True)

    expected = expected.to_array()
    computed = computed.to_array()

    return qcel.testing.compare_recursive(expected, computed, **kwargs)


def _mergedapis_compare_matrices(expected, computed, *args, **kwargs):
    """Shim allowing Psi4-style or QCA-style testing interfaces for :py:class:`psi4.core.Matrix`."""

    qcdb.testing._merge_psi4_qcel_apis(args, kwargs)

    if kwargs.pop("check_name", False):
        compare(expected.name, computed.name, f'{expected.name} vs. {computed.name} name', quiet=True)

    compare(expected.nirrep(), computed.nirrep(), f'{expected.name} vs. {computed.name} irreps', quiet=True)
    compare(expected.symmetry(), computed.symmetry(), f'{expected.name} vs. {computed.name} symmetry', quiet=True)
    for irrep in range(expected.nirrep()):
        compare(
            expected.rows(irrep),
            computed.rows(irrep),
            f'{expected.name} vs. {computed.name} irrep {irrep} rows ({expected.rows(irrep)} vs {computed.rows(irrep)})',
            quiet=True)
        compare(
            expected.cols(irrep ^ expected.symmetry()),
            computed.cols(irrep ^ expected.symmetry()),
            f'{expected.name} vs. {computed.name} irrep {irrep} columns ({expected.cols(irrep)} vs {computed.cols(irrep)})',
            quiet=True)

    expected = expected.to_array()
    computed = computed.to_array()

    return qcel.testing.compare_recursive(expected, computed, **kwargs)


def _mergedapis_compare_wavefunctions(expected, computed, *args, **kwargs):
    """Shim allowing Psi4-style or QCA-style testing interfaces for :py:class:`psi4.core.Wavefunction`."""

    qcdb.testing._merge_psi4_qcel_apis(args, kwargs)

    if kwargs['label'] == '<module>':
        kwargs['label'] = 'Wavefunctions equal'
    kwargscopy = kwargs.copy()
    kwargscopy.pop('label')
    atol = kwargscopy.pop('atol', 1.e-9)

    # yapf: disable
    if expected.Ca():          compare_matrices(expected.Ca(), computed.Ca(), 'compare Ca', atol=atol, **kwargscopy)
    if expected.Cb():          compare_matrices(expected.Cb(), computed.Cb(), 'compare Cb', atol=atol, **kwargscopy)
    if expected.Da():          compare_matrices(expected.Da(), computed.Da(), 'compare Da', atol=atol, **kwargscopy)
    if expected.Db():          compare_matrices(expected.Db(), computed.Db(), 'compare Db', atol=atol, **kwargscopy)
    if expected.Fa():          compare_matrices(expected.Fa(), computed.Fa(), 'compare Fa', atol=atol, **kwargscopy)
    if expected.Fb():          compare_matrices(expected.Fb(), computed.Fb(), 'compare Fb', atol=atol, **kwargscopy)
    if expected.H():           compare_matrices(expected.H(), computed.H(), 'compare H', atol=atol, **kwargscopy)
    if expected.S():           compare_matrices(expected.S(), computed.S(), 'compare S', atol=atol, **kwargscopy)
    if expected.lagrangian():  compare_matrices(expected.lagrangian(), computed.lagrangian(), 'compare lagrangian', atol=atol, **kwargscopy)
    if expected.aotoso():      compare_matrices(expected.aotoso(), computed.aotoso(), 'compare aotoso', atol=atol, **kwargscopy)
    if expected.gradient():    compare_matrices(expected.gradient(), computed.gradient(), 'compare gradient', atol=atol, **kwargscopy)
    if expected.hessian():     compare_matrices(expected.hessian(), computed.hessian(), 'compare hessian', atol=atol, **kwargscopy)
    if expected.epsilon_a():   compare_vectors(expected.epsilon_a(), computed.epsilon_a(), 'compare epsilon_a', atol=atol, **kwargscopy)
    if expected.epsilon_b():   compare_vectors(expected.epsilon_b(), computed.epsilon_b(), 'compare epsilon_b', atol=atol, **kwargscopy)
    if expected.frequencies(): compare_vectors(expected.frequencies(), computed.frequencies(), 'compare frequencies', atol=atol, **kwargscopy)
    # yapf: enable
    compare(expected.nalpha(), computed.nalpha(), 'compare nalpha', **kwargscopy)
    compare(expected.nbeta(), computed.nbeta(), 'compare nbeta', **kwargscopy)
    compare(expected.nfrzc(), computed.nfrzc(), 'compare nfrzc', **kwargscopy)
    compare(expected.nirrep(), computed.nirrep(), 'compare nirrep', **kwargscopy)
    compare(expected.nmo(), computed.nmo(), 'compare nmo', **kwargscopy)
    compare(expected.nso(), computed.nso(), 'compare nso', **kwargscopy)
    compare(expected.name(), computed.name(), 'compare name', **kwargscopy)
    compare(expected.module(), computed.module(), 'compare module', **kwargscopy)
    compare_values(expected.energy(), computed.energy(), 'compare energy', atol=atol, **kwargscopy)
    compare_values(expected.efzc(), computed.efzc(), 'compare frozen core energy', atol=atol, **kwargscopy)
    compare_values(expected.get_dipole_field_strength()[0],
                   computed.get_dipole_field_strength()[0],
                   'compare dipole field strength x',
                   atol=atol,
                   **kwargscopy)
    compare_values(expected.get_dipole_field_strength()[1],
                   computed.get_dipole_field_strength()[1],
                   'compare dipole field strength y',
                   atol=atol,
                   **kwargscopy)
    compare_values(expected.get_dipole_field_strength()[2],
                   computed.get_dipole_field_strength()[2],
                   'compare dipole field strength z',
                   atol=atol,
                   **kwargscopy)

    compare(expected.basisset().name(), computed.basisset().name(), 'compare basis set name', **kwargscopy)
    compare(expected.basisset().nbf(),
            computed.basisset().nbf(), 'compare number of basis functions in set', **kwargscopy)
    compare(expected.basisset().nprimitive(),
            computed.basisset().nprimitive(), 'compare total number of primitives in basis set', **kwargscopy)

    compare(expected.molecule().name(), computed.molecule().name(), 'compare molecule name', **kwargscopy)
    compare(expected.molecule().get_full_point_group(),
            computed.molecule().get_full_point_group(), 'compare molecule point group', **kwargscopy)
    compare_matrices(expected.molecule().geometry(),
                     computed.molecule().geometry(),
                     'compare molecule geometry',
                     atol=atol,
                     **kwargscopy)

    return qcel.testing.compare(True, True, **kwargs)


def _psi4_true_raise_handler(passfail, label, message, return_message=False, quiet=False):
    """Handle comparison result by printing to screen, printing to Psi output file, raising TestComparisonError, and (incidently) returning."""

    width = 86
    if passfail:
        if not quiet:
            core.print_out(f'    {label:.<{width}}PASSED\n')
            print(f'    {label:.<{width}}PASSED')
            sys.stdout.flush()
    else:
        core.print_out(f'    {label:.<{width}}FAILED')
        print(f'    {label:.<{width}}FAILED')
        sys.stdout.flush()
        raise TestComparisonError(message)

    return passfail


# ~POD
compare = partial(qcel.testing.compare, return_handler=_psi4_true_raise_handler)
compare_integers = partial(qcdb.testing._psi4_compare_integers, return_handler=_psi4_true_raise_handler)
compare_strings = compare_integers
compare_values = partial(qcdb.testing._mergedapis_compare_values, return_handler=_psi4_true_raise_handler)
compare_arrays = compare_values

# Psi4-only
compare_cubes = partial(_mergedapis_compare_cubes, return_handler=_psi4_true_raise_handler)
compare_vectors = partial(_mergedapis_compare_vectors, return_handler=_psi4_true_raise_handler)
compare_matrices = partial(_mergedapis_compare_matrices, return_handler=_psi4_true_raise_handler)
compare_wavefunctions = partial(_mergedapis_compare_wavefunctions, return_handler=_psi4_true_raise_handler)

# dict-like, QCEl API only
compare_recursive = partial(qcdb.testing._mergedapis_compare_recursive, return_handler=_psi4_true_raise_handler)
compare_molrecs = partial(qcdb.testing._mergedapis_compare_molrecs, return_handler=_psi4_true_raise_handler)

## docstrings

_psi4_style_doc = """
    Handles both Psi4-style signatures (``(expected, computed, atol_exponent, label)``; see **atol_exponent** parameter below) and QCA-style signatures (``(expected, computed, label)``).

    Parameters
    ----------
    atol_exponent : int or float
        Absolute tolerance (see formula in :py:func:`qcelemental.testing.compare_values` notes).
        Values less than one are taken literally; one or greater taken as decimal digits for comparison.
        So `1` means `atol=0.1` and `2` means `atol=0.01` but `0.04` means `atol=0.04`
        Note that the largest expressable processed atol will be `~0.99`.

"""

compare_values.__doc__ = """Comparison function for float or float array-like data structures.
    See :py:func:`qcelemental.testing.compare_values` for details.

    ``psi4.compare_arrays`` is an old comparison function for float NumPy arrays that is now an alias to this.

""" + _psi4_style_doc

compare_strings.__doc__ = """Comparison function for integers, strings, booleans, or integer array-like data structures.
    See :py:func:`qcelemental.testing.compare` for details.

    ``psi4.compare_strings`` is an alias to this.

"""
compare_integers.__doc__ = compare_strings.__doc__

compare_cubes.__doc__ = """Comparison function for volumetric data in cube file format.
    Compares only the volumetric data, not the voxel data or molecular geometry or other header contents.
    The volumetric data is passed to :py:func:`qcelemental.testing.compare_values`.

    Note only QCA-style signature (``(expected, computed, label)``) available.

    Parameters
    ----------
    expected
        Reference cube file against which `computed` is compared.
        Read by :func:`numpy.genfromtxt` so `expected` can be any of file, str,
        pathlib.Path, list of str, generator.
    computed
        Input cube file to compare against `expected`.
        Read by :func:`numpy.genfromtxt` so `computed` can be any of file, str,
        pathlib.Path, list of str, generator.

"""

compare_vectors.__doc__ = """Comparison function for :py:class:`psi4.core.Vector` objects.
    Compares Vector properties of ``name`` (optional through **check_name**), ``nirrep``, and dimension of each irrep.
    For comparing actual numerical contents, the vectors are serialized to NumPy array format and passed to :py:func:`qcelemental.testing.compare_recursive`.
""" + _psi4_style_doc

compare_matrices.__doc__ = """Comparison function for :py:class:`psi4.core.Matrix` objects.
    Compares Matrix properties of ``name`` (optional through **check_name**), ``nirrep``, ``symmetry``, and number of rows and columns for each irrep.
    For comparing actual numerical contents, the matrices are serialized to NumPy array format and passed to :py:func:`qcelemental.testing.compare_recursive`.
""" + _psi4_style_doc

compare_wavefunctions.__doc__ = """Comparison function for :py:class:`psi4.core.Wavefunction` objects.
    Compares over 30 Wavefunction properties, including ``nirrep``, ``nso``, molecule geometry, basis set ``nbf``, density matrices, gradient results, etc.

""" + _psi4_style_doc

compare_recursive.__doc__ = """Comparison function for recursively comparing mixed-type and nested structures such as dictionaries and lists.
    See :py:func:`qcelemental.testing.compare_recursive` for details.

"""

compare_molrecs.__doc__ = """Comparison function for :py:func:`psi4.core.Molecule.to_dict` objects.
    See :py:func:`qcelemental.testing.compare_molrecs` for details.

    Note only QCA-style signature (``(expected, computed, label)``) available.

"""
