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
"""Module with utility functions for FCHK files."""

import numpy as np
from psi4.driver.p4util.testing import compare_strings, compare_arrays, compare_values, compare_integers
from psi4 import core
from .exceptions import ValidationError


def consume_fchk_section(input_list, index):
    """compare a float matrix section"""

    n = int(input_list[index].split()[-1])
    kind = input_list[index].split()[-3]

    if " R " in kind:
        dtype = np.float64
        format_counter = 5
    elif " I " in kind:
        dtype = np.float64
        format_counter = 6
    else:
        raise ValidationError('Unknow field type in FCHK reader\n')

    extra = 0 if n <= format_counter else n % format_counter
    lines = 1 if n <= format_counter else int(n / format_counter)
    offset = lines + 1 if extra > 0 else lines
    string = ''
    for j in range(lines):
        string += "".join(str(x) for x in input_list[index + 1 + j])
    if extra > 0:
        string += "".join(str(x) for x in input_list[index + 1 + lines])
    field = np.fromiter(string.split(), dtype=dtype)
    return offset + 1, field


def fchkfile_to_string(fname):
    """ Load FCHK file into a string"""
    with open(fname, 'r') as handle:
        fchk_string = handle.read()
    return fchk_string


def compare_fchkfiles(expected, computed, label):
    # """Function to compare two FCHK files.

    # :param expected: reference FCHK file name
    # :param computed: computed FCHK file name
    # :param label: string labelling the test
    # """

    high_accuracy = 9
    low_accuracy = 3

    # those need super high scf convergence (d_conv 1e-12). Unsure how it is
    # on different machines. Thus they get low_accuracy.
    # Rest seems okay with d_conv 1e-10.
    sensitive = ['Current cartesian coordinates', 'MO coefficients']

    fchk_ref = fchkfile_to_string('dct.ref').splitlines()
    fchk_calc = fchkfile_to_string('dct.fchk').splitlines()

    if len(fchk_ref) != len(fchk_calc):
        raise ValidationError('The two FCHK files to compare have a different file length! \n')

    i = 0
    max_length = len(fchk_calc)
    tests = []
    for start in range(max_length):
        if i >= max_length:
            break
        if "N=" in fchk_calc[i]:
            offset, calc = consume_fchk_section(fchk_calc, i)
            _, ref = consume_fchk_section(fchk_ref, i)
            if any(x in fchk_calc[i] for x in sensitive):
                test = compare_arrays(ref, calc, low_accuracy, f" matrix section: {fchk_calc[i]}")
            else:
                test = compare_arrays(ref, calc, high_accuracy, f" matrix section: {fchk_calc[i]}")
            i += offset
        elif " R " in fchk_calc[i] and not "N=" in fchk_calc[i]:
            calc = fchk_calc[i].split()[-1]
            ref = fchk_ref[i].split()[-1]
            test = compare_values(ref, calc, high_accuracy, f" float value: {fchk_calc[i]}")
            i += 1
        elif " I " in fchk_calc[i] and not "N=" in fchk_calc[i]:
            calc = fchk_calc[i].split()[-1]
            ref = fchk_ref[i].split()[-1]
            test = compare_integers(ref, calc, f" int value: {fchk_calc[i]}")
            i += 1
        else:
            test = compare_strings(fchk_calc[i], fchk_ref[i], f"FCK text line {i+1}.")
            i += 1
        tests.append(test)

    return compare_integers(True, all(tests), label)