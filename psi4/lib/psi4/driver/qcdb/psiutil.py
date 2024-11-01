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

r"""Stuff stolen from psi. Should import or not as necessary
or some better way. Apologies to the coders.

"""
import os
import re
import sys


def query_yes_no(question, default=True):
    """Ask a yes/no question via input() and return their answer.

    *question* is a string that is presented to the user.
    *default* is the presumed answer if the user just hits <Enter>.
    It must be yes (the default), no or None (meaning
    an answer is required of the user).

    The return value is one of True or False.

    """

    yes = re.compile(r'^(y|yes|true|on|1)', re.IGNORECASE)
    no = re.compile(r'^(n|no|false|off|0)', re.IGNORECASE)

    if default is None:
        prompt = " [y/n] "
    elif default is True:
        prompt = " [Y/n] "
    elif default is False:
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().strip().lower()
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
