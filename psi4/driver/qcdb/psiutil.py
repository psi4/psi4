#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

r"""Stuff stolen from psi. Should import or not as necessary
or some better way. Apologies to the coders.

"""
from __future__ import absolute_import
from __future__ import print_function
import sys
import math
import re
import os
import string
from .vecutil import *


def _success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    print('\t{0:.<66}PASSED'.format(label))
    sys.stdout.flush()


def compare_values(expected, computed, digits, label, exitonfail=True):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*
    (or to *digits* itself when *digits* > 1 e.g. digits=0.04). Performs
    a system exit on failure unless *exitonfail* False, in which case
    returns error message. Used in input files in the test suite.

    """
    thresh = 10 ** -digits if digits > 1 else digits
    if abs(expected - computed) > thresh:
        print("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
        if exitonfail:
            sys.exit(1)
        else:
            return
    if math.isnan(computed):
        print("\t%s: computed value (%f) does not match (%f)\n" % (label, computed, expected))
        print("\tprobably because the computed value is nan.")
        if exitonfail:
            sys.exit(1)
        else:
            return
    _success(label)


def compare_integers(expected, computed, label):
    """Function to compare two integers. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if (expected != computed):
        print("\t%s: computed value (%d) does not match (%d)." % (label, computed, expected))
        sys.exit(1)
    _success(label)


def compare_strings(expected, computed, label):
    """Function to compare two strings. Prints :py:func:`util.success`
    when string *computed* exactly matches string *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if(expected != computed):
        print("\t%s: computed value (%s) does not match (%s)." % (label, computed, expected))
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
        if default is not None and choice == '':
            return default
        elif yes.match(choice):
            return True
        elif no.match(choice):
            return False
        else:
            sys.stdout.write("    Please respond with 'yes' or 'no'.\n")


## {{{ http://code.activestate.com/recipes/52224/ (r1)
def search_file(filename, search_path):
    """Given an os.pathsep divided *search_path*, find first occurance of
    *filename*. Returns full path to file if found or None if unfound.

    """
    file_found = False
    paths = search_path.split(os.pathsep)
    #paths = string.split(search_path, os.pathsep)
    for path in paths:
        if os.path.exists(os.path.join(path, filename)):
            file_found = True
            break
    if file_found:
        return os.path.abspath(os.path.join(path, filename))
    else:
        return None
## end of http://code.activestate.com/recipes/52224/ }}}


def drop_duplicates(seq):
    """Function that given an array *seq*, returns an array without any duplicate
    entries. There is no guarantee of which duplicate entry is dropped.

    """
    #noDupes = []
    #[noDupes.append(i) for i in seq if not noDupes.count(i)]
    #return noDupes
    noDupes = []
    seq2 = sum(seq, [])
    [noDupes.append(i) for i in seq2 if not noDupes.count(i)]
    return noDupes


def all_casings(input_string):
    """Function to return a generator of all lettercase permutations
    of *input_string*.

    """
    if not input_string:
        yield ''
    else:
        first = input_string[:1]
        if first.lower() == first.upper():
            for sub_casing in all_casings(input_string[1:]):
                yield first + sub_casing
        else:
            for sub_casing in all_casings(input_string[1:]):
                yield first.lower() + sub_casing
                yield first.upper() + sub_casing


def getattr_ignorecase(module, attr):
    """Function to extract attribute *attr* from *module* if *attr*
    is available in any possible lettercase permutation. Returns
    attribute if available, None if not.

    """
    array = None
    for per in list(all_casings(attr)):
        try:
            getattr(module, per)
        except AttributeError:
            pass
        else:
            array = getattr(module, per)
            break

    return array


def import_ignorecase(module):
    """Function to import *module* in any possible lettercase
    permutation. Returns module object if available, None if not.

    """
    modobj = None
    for per in list(all_casings(module)):
        try:
            modobj = __import__(per)
        except ImportError:
            pass
        else:
            break

    return modobj

def findfile_ignorecase(fil, pre='', post=''):
    """Function to locate a file *pre* + *fil* + *post* in any possible 
    lettercase permutation of *fil*. Returns *pre* + *fil* + *post* if 
    available, None if not.

    """
    afil = None
    for per in list(all_casings(fil)):
        if os.path.isfile(pre + per + post):
            afil = pre + per + post
            break
        else:
            pass

    return afil

