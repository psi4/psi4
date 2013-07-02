#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

r"""Stuff stolen from psi. Should import or not as necessary
or some better way. Apologies to the coders.

"""
import sys
import math
import re
from vecutil import *


def _success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    print('\t{0:.<66}PASSED'.format(label))
    sys.stdout.flush()


def compare_values(expected, computed, digits, label):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if abs(expected - computed) > 10 ** (-digits):
        print("\t%s: computed value (%f) does not match (%f) to %d digits." % (label, computed, expected, digits))
        sys.exit(1)
    if math.isnan(computed):
        print("\t%s: computed value (%f) does not match (%f) to %d digits.\n" % (label, computed, expected, digits))
        print("\tprobably because the computed value is nan.")
        sys.exit(1)
    _success(label)


def compare_matrices(expected, computed, digits, label):
    """Function to compare two matrices. Prints :py:func:`util.success`
    when elements of matrix *computed* match elements of matrix *expected* to
    number of *digits*. Performs a system exit on failure to match symmetry
    structure, dimensions, or element values. Used in input files in the test suite.

    """
    rows = len(expected)
    cols = len(expected[0])
    failed = 0
    for row in range(rows):
        for col in range(cols):
            if abs(expected[row][col] - computed[row][col]) > 10 ** (-digits):
                print("\t%s: computed value (%s) does not match (%s)." % (label, expected[row][col], computed[row][col]))
                failed = 1
                break

    if failed:
        print("The Failed Test Matrices\n")
        show(computed)
        print('\n')
        show(expected)
        sys.exit(1)
    _success(label)


def query_yes_no(question, default=True):
    """Ask a yes/no question via raw_input() and return their answer.

    *question* is a string that is presented to the user.
    *default* is the presumed answer if the user just hits <Enter>.
    It must be yes (the default), no or None (meaning
    an answer is required of the user).

    The return value is one of True or False.

    """

    yes = re.compile(r'^(y|yes|true|on|1)', re.IGNORECASE)
    no = re.compile(r'^(n|no|false|off|0)', re.IGNORECASE)

    if default == None:
        prompt = " [y/n] "
    elif default == True:
        prompt = " [Y/n] "
    elif default == False:
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().strip().lower()
        if default is not None and choice == "":
            return default
        elif yes.match(choice):
            return True
        elif no.match(choice):
            return False
        else:
            sys.stdout.write("    Please respond with 'yes' or 'no'.\n")
