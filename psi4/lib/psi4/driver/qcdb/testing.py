#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import sys
from typing import Callable
from functools import partial

import numpy as np
import qcelemental as qcel

from .exceptions import TestComparisonError, UpgradeHelper

__all__ = [
    'compare', 'compare_integers', 'compare_strings', 'compare_values', 'compare_arrays', 'compare_recursive',
    'compare_molrecs', 'compare_matrices', 'compare_dicts'
]


def _merge_psi4_qcel_apis(args, kwargs):
    """Outer shim to allow both Psi4-style and QCA-style testing interfaces through the same function.

    Notes
    -----
    `kwargs` modified (and returned) in-place

    """

    def process_digits(digits):
        if digits >= 1:
            return 10**-digits
        return digits

    if len(args) == 0:
        kwargs['label'] = sys._getframe().f_back.f_back.f_code.co_name

    elif len(args) == 1:
        if isinstance(args[0], str):
            kwargs['label'] = args[0]

        else:
            kwargs['atol'] = process_digits(args[0])
            kwargs['label'] = sys._getframe().f_back.f_back.f_code.co_name
            if 'verbose' in kwargs:
                kwargs['quiet'] = (kwargs.pop('verbose') < 1)

    elif len(args) == 2:
        kwargs['atol'] = process_digits(args[0])
        kwargs['label'] = args[1]
        if 'verbose' in kwargs:
            kwargs['quiet'] = (kwargs.pop('verbose') < 1)

    else:
        raise ValueError("""Not following either Psi4 or QCElemental API pattern for comparison.""")


def _psi4_compare_integers(expected, computed, label = None, *, verbose = 1,
                           return_handler = None):
    """Shim between Psi4-style and QCA-style testing interfaces for scalar ints, strings."""

    # uncomment to ferret out old function name
    #raise UpgradeHelper('qcdb.compare_integers', 'qcdb.compare', 'someday', ' Same API, just rename the function and convert `verbose` to `quiet`.')

    if label is None:
        label = sys._getframe().f_back.f_code.co_name

    return qcel.testing.compare(expected,
                                computed,
                                label=label,
                                quiet=(verbose == 0),
                                return_message=False,
                                return_handler=return_handler)


def _mergedapis_compare_values(expected, computed, *args, **kwargs):
    """Outer shim to allow both Psi4-style and QCA-style testing interfaces through the same function."""
    """Shim between Psi4-style and QCA-style testing interfaces for scalar and array floats."""

    _merge_psi4_qcel_apis(args, kwargs)
    return qcel.testing.compare_values(expected, computed, **kwargs)


def _mergedapis_compare_recursive(expected, computed, *args, **kwargs):

    if (len(args) > 0) and not isinstance(args[0], str):
        raise UpgradeHelper(
            'qcdb.compare_recursive', 'qcdb.compare_recursive', 1.4,
            ' Use the new `qcel.testing.compare_recursive` API, being sure to convert positional arg `digits` decimal places to keyword arg `atol` literal absolute tolerance.'
        )

    return qcel.testing.compare_recursive(expected, computed, *args, **kwargs)


def _mergedapis_compare_molrecs(expected, computed, *args, **kwargs):

    if (len(args) > 0) and not isinstance(args[0], str):
        raise UpgradeHelper(
            'qcdb.compare_molrecs', 'qcdb.compare_molrecs', 1.4,
            ' Use the new `qcel.testing.compare_molrecs` API, being sure to convert positional arg `digits` decimal places to keyword arg `atol` literal absolute tolerance.'
        )

    _merge_psi4_qcel_apis(args, kwargs)
    return qcel.testing.compare_molrecs(expected, computed, **kwargs)


def compare_matrices(expected, computed, *args, **kwargs):
    raise UpgradeHelper(
        'qcdb.compare_matrices', 'qcdb.compare_values', 1.4,
        ' Use the new qcel.testing.compare_values` API, being sure to convert `digits` decimal places to `atol` literal absolute tolerance.'
    )


def compare_dicts(expected, computed, *args, **kwargs):

    raise UpgradeHelper(
        'qcdb.compare_dicts', 'qcdb.compare_recursive', 1.4,
        ' Use the new `qcel.testing.compare_recursive` API, being sure to convert `tol` decimal places to `atol` literal absolute tolerance.'
    )


def _qcdb_true_raise_handler(passfail, label, message, return_message=False, quiet=False):
    """Handle comparison result by printing to screen and raising qcdb.TestComparisonError or returning True."""

    width = 66
    if passfail:
        if not quiet:
            print(f'    {label:.<{width}}PASSED')
            sys.stdout.flush()
    else:
        print(f'    {label:.<{width}}FAILED')
        sys.stdout.flush()
        raise TestComparisonError(message)

    return passfail


compare = partial(qcel.testing.compare, return_handler=_qcdb_true_raise_handler)
compare_integers = partial(_psi4_compare_integers, return_handler=_qcdb_true_raise_handler)
compare_strings = compare_integers
compare_values = partial(_mergedapis_compare_values, return_handler=_qcdb_true_raise_handler)
compare_arrays = compare_values

compare_recursive = partial(_mergedapis_compare_recursive, return_handler=_qcdb_true_raise_handler)
compare_molrecs = partial(_mergedapis_compare_molrecs, return_handler=_qcdb_true_raise_handler)

# Notes on testing fns migration
# PSI4
#   ADDED               def compare                                                                                         SINGLE
#   MERGED-APIs         def compare_integers(expected, computed, label, verbose=1):                                         SINGLE
#   MERGED-APIs         def compare_strings(expected, computed, label, verbose=1):                                          SINGLE
#   MERGED-APIs         def compare_values(expected, computed, digits, label, *, rtol=1.e-16, passnone=False, verbose=1):   SINGLE
#   MERGED-APIs         def compare_arrays(expected, computed, digits, label, rtol=1.e-16, verbose=1):                      SINGLE
#
#   ADDED-NEW-API       def compare_recursive                                                                               SINGLE
#   ADDED-NEW-API       def compare_molrecs                                                                                 SINGLE
#
#   MERGED-APIs         def compare_cubes(expected, computed, label, verbose=1):                                            SINGLE
#   MERGED-APIs         def compare_vectors(expected, computed, digits, label, *, rtol=1.e-16, verbose=1):                  MULTI
#   MERGED-APIs         def compare_matrices(expected, computed, digits, label, *, rtol=1.e-16, verbose=1):                 MULTI
#   MERGED-APIs         def compare_wavefunctions(expected, computed, digits=9, label='Wavefunctions equal'):               MULTI
#   def compare_fcidumps(expected, computed, label):                                                                        MULTI

# QCDB
#   ADDED               def compare                                                                                         SINGLE
#   MERGED-APIs-TRIVIAL def compare_integers(expected, computed, label, verbose=1):                                         SINGLE
#   MERGED-APIs-TRIVIAL def compare_strings(expected, computed, label, verbose=1):                                          SINGLE
#   MERGED-APIs         def compare_values(expected, computed, digits, label, passnone=False, verbose=1):                   SINGLE
#   MERGED-APIs         def compare_arrays(expected, computed, digits, label, verbose=1):                                   SINGLE
#
#   ADDED-NEW-API       def compare_recursive                                                                               SINGLE
#   STOPCONVERT/NEW-API def compare_molrecs(expected, computed, tol, label, forgive=None, verbose=1, relative_geoms='exact' SINGLE
#
#   STOPCONVERT/NEW-FN  def compare_matrices(expected, computed, digits, label, verbose=1):                                 ---
#   STOPCONVERT/NEW-FN  def compare_dicts(expected, computed, tol, label, forgive=None, verbose=1):                         ---
#   vib.py:def compare_vibinfos(expected, computed, tol, label, verbose=1, forgive=None, required=None, toldict=None):      SINGLE
