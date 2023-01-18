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
"""Module with utility functions for FCHK files."""

import re
from typing import Union

import numpy as np
from psi4.driver.p4util.testing import compare_strings, compare_arrays, compare_values, compare_integers
from psi4 import core
from .exceptions import ValidationError

__all__ = [
    "compare_fchkfiles",
    "compare_moldenfiles",
]

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


def _fchkfile_to_string(fname: str) -> str:
    """ Load FCHK file into a string"""
    with open(fname, 'r') as handle:
        fchk_string = handle.read()
    return fchk_string


def compare_fchkfiles(expected: str, computed: str, atol_exponent: Union[int, float], label: str):
    """Comparison function for output data in FCHK (formatted checkpoint) file
    format. Compares many fields including number of electrons, highest angular
    momentum, basis set exponents, densities, final gradient.

    Note only Psi4-style signature (``(expected, computed, atol_exponent, label)``) available.

    An older format description can be found here
    http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/f_formchk.htm
    It lists more fields (logical, character) that are not included in this
    test function. They should be covered by the string comparison.
    This function is only meant to work with PSI4's FCHK files.

    Parameters
    ----------
    expected
        Path to reference FCHK file against which `computed` is compared.
    computed
        Path to input FCHK file to compare against `expected`.
    atol_exponent
        Absolute tolerance for high accuracy fields -- 1.e-8 or 1.e-9 is suitable.
        Values less than one are taken literally; one or greater taken as decimal digits for comparison.
        So `1` means `atol=0.1` and `2` means `atol=0.01` but `0.04` means `atol=0.04`
        Note that the largest expressable processed atol will be `~0.99`.
    label
        Label for passed and error messages.

    """
    fchk_ref = _fchkfile_to_string(expected).splitlines()
    fchk_calc = _fchkfile_to_string(computed).splitlines()

    high_accuracy = atol_exponent
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
        elif " R " in line and "N=" not in line:
            calc = line.split()[-1]
            ref = fchk_ref[index].split()[-1]
            test = compare_values(ref, calc, high_accuracy, f" float value: {line}")
            index += 1
        elif " I " in line and "N=" not in line:
            calc = line.split()[-1]
            ref = fchk_ref[index].split()[-1]
            test = compare_integers(ref, calc, f" int value: {line}")
            index += 1
        else:
            test = compare_strings(fchk_ref[index], line, f"FCK text line {index+1}.")
            index += 1
        tests.append(test)

    return compare_integers(True, all(tests), label)


def compare_moldenfiles(
    expected: str,
    computed: str,
    atol_exponent: Union[int, float] = 1.e-7,
    label: str = "Compare Molden"):
    """Comparison function for output data in Molden file format.
    Compares many fields including geometry, basis set, occupations, symmetries, energies.

    Note only Psi4-style signature (``(expected, computed, atol_exponent, label)``) available.

    A format description is found https://www3.cmbi.umcn.nl/molden/molden_format.html

    Parameters
    ----------
    expected
        Path to reference Molden file against which `computed` is compared.
    computed
        Path to input Molden file to compare against `expected`.
    atol_exponent
        Absolute tolerance for high accuracy fields -- 1.e-8 or 1.e-9 is suitable.
        Values less than one are taken literally; one or greater taken as decimal digits for comparison.
        So `1` means `atol=0.1` and `2` means `atol=0.01` but `0.04` means `atol=0.04`
        Note that the largest expressable processed atol will be `~0.99`.
    label
        Label for passed and error messages.

    """
    def moldenfile_to_string(fname):
        with open(fname, 'r') as fn:
            molden_string = fn.read()
        return molden_string

    ref = moldenfile_to_string(expected).splitlines()
    calc = moldenfile_to_string(computed).splitlines()
    if len(ref) != len(calc):
        raise ValidationError(f"These two molden files have different lengths...\n")

    high_accuracy = atol_exponent
    index = 0
    max_len = len(calc)
    tests = []
    section = 0

    geom_re = re.compile(r'^\s*(\w*)\s+(\d+)\s+(\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s*$')
    basis_header_re = re.compile(r'^\s*([s,p,d,f,g])\s*(\d*)\s*(\d*.\d*)\s*$')
    s1_re = re.compile(r'^\s*(\d+.?\d*)\s+(\d+.?\d*)$')
    s2_re = re.compile(r'^\s*(\d+)\s+(-?\d+.\d+[e,E][\+,-]\d+)\s*$')
    sym_re = re.compile(r'^\s*Sym\s*=\s*(\w*)\s*$')
    energy_re = re.compile(r'^\s*Ene\s*=\s*(-?\d*.?\d*[e,E]?\+?-?\d*)\s*$')
    spin_re = re.compile(r'^\s*Spin\s*=\s*(\w*)\s*$')
    occ_re = re.compile(r'^\s*Occup\s*=\s*(-?\d*.\d*[e,E]?-?\+?\d*)\s*$')

    for i in range(max_len):
        line = calc[i]

        if geom_re.match(line):
            c1, c2, c3, c4, c5, c6 = geom_re.match(line).groups()
            r1, r2, r3, r4, r5, r6 = geom_re.match(line).groups()
            test = compare_strings(r1, c1) and compare_integers(r2, c2) and compare_integers(r3, c3) and compare_values(r4, c4, high_accuracy) and compare_values(r5, c5, high_accuracy) and compare_values(r6, c6, high_accuracy)

        elif basis_header_re.match(line):
            c1, c2, c3 = basis_header_re.match(line).groups()
            r1, r2, r3 = basis_header_re.match(ref[i]).groups()
            test = compare_strings(r1,c1) and compare_integers(r2,c2) and compare_values(r3,c3,3)

        elif s1_re.match(line):
            c1, c2 = s1_re.match(line).groups()
            r1, r2 = s1_re.match(ref[i]).groups()
            test = compare_values(r1, c1, high_accuracy) and compare_values(r2, c2, high_accuracy)

        elif sym_re.match(line):
            c = sym_re.match(line).group(1)
            r = sym_re.match(ref[i]).group(1)
            test = compare_strings(r, c, f'text line: {line}')

        elif energy_re.match(line):
            c = energy_re.match(line).group(1)
            r = energy_re.match(ref[i]).group(1)
            test = compare_values(r, c, high_accuracy, f'float value: {line}')

        elif spin_re.match(line):
            c = spin_re.match(line).group(1)
            r = spin_re.match(ref[i]).group(1)
            test = compare_strings(r, c, f'text line: {line}')

        elif occ_re.match(line):
            c = occ_re.match(line).group(1)
            r = occ_re.match(ref[i]).group(1)
            test = compare_values(r, c, high_accuracy, f'float value: {line}')

        elif s2_re.match(line):
            c1, c2 = s2_re.match(line).groups()
            r1, r2 = s2_re.match(line).groups()
            test = compare_integers(r1, c1, f'int value: {line}') and compare_values(r2, c2, high_accuracy, f'float value: {line}')

        else:
            test = compare_strings(line, ref[i])

        tests.append(test)

    return compare_integers(True, all(tests), label)
