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
"""Module with utility functions for FCHK files."""

import numpy as np
from psi4.driver.p4util.testing import compare_strings, compare_arrays, compare_values, compare_integers
from psi4 import core
from .exceptions import ValidationError

__all__ = ['fchkfile_to_string','compare_fchkfiles']

def _consume_fchk_section(input_list, index):
    """compare a float or integer matrix section"""

    n = int(input_list[index].split()[-1])
    kind = input_list[index].split()[-3]

    if "R" in kind:
        dtype = np.float64
        format_counter = 5
    elif "I" in kind:
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


def compare_fchkfiles(expected, computed, digits, label):
    # """Function to compare two FCHK files.

    # an older format description can be found here
    # http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/f_formchk.htm
    # It lists more fields (logical, character) that are not included in this
    # test function. They should be covered by the string comparison.
    # This function is only meant to work with PSI4's FCHK files.
    #
    # :param expected: reference FCHK file name
    # :param computed: computed FCHK file name
    # :param digits: tolerance for high accuracy fields -- 1.e-8 or 1.e-9 suitable
    # :param label: string labelling the test
    # """

    fchk_ref = fchkfile_to_string(expected).splitlines()
    fchk_calc = fchkfile_to_string(computed).splitlines()

    high_accuracy = digits
    low_accuracy = 3

    # Those listed below need super high scf convergence (d_conv 1e-12) and might
    # show machine dependence. They will be tested with low_accuracy.
    sensitive = ['Current cartesian coordinates', 'MO coefficients']

    if len(fchk_ref) != len(fchk_calc):
        raise ValidationError('The two FCHK files to compare have a different file length! \n')

    index = 0
    max_length = len(fchk_calc)
    tests = []
    for start in range(max_length):
        if index >= max_length:
            break
        line = fchk_calc[index]
        if "N=" in line:
            offset, calc = _consume_fchk_section(fchk_calc, index)
            _, ref = _consume_fchk_section(fchk_ref, index)
            if any(x in line for x in sensitive):
                test = compare_arrays(ref, calc, low_accuracy, f" matrix section: {line}")
            else:
                test = compare_arrays(ref, calc, high_accuracy, f" matrix section: {line}")
            index += offset
        elif " R " in line and not "N=" in line:
            calc = line.split()[-1]
            ref = fchk_ref[index].split()[-1]
            test = compare_values(ref, calc, high_accuracy, f" float value: {line}")
            index += 1
        elif " I " in line and not "N=" in line:
            calc = line.split()[-1]
            ref = fchk_ref[index].split()[-1]
            test = compare_integers(ref, calc, f" int value: {line}")
            index += 1
        else:
            test = compare_strings(line, fchk_ref[index], f"FCK text line {index+1}.")
            index += 1
        tests.append(test)

    return compare_integers(True, all(tests), label)
