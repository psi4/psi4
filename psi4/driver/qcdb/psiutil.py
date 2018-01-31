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


def compare_values(expected, computed, digits, label, passnone=False, verbose=1):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*
    (or to *digits* itself when *digits* > 1 e.g. digits=0.04).
    Used in input files in the test suite.

    Raises
    ------
    TestComparisonError
        If `computed` differs from `expected` by more than `digits`.

    """
    if passnone:
        if expected is None and computed is None:
            _success(label)
            return

    if digits > 1:
        thresh = 10 ** -digits
        message = ("\t%s: computed value (%.*f) does not match (%.*f) to %d digits." % (label, digits+1, computed, digits+1, expected, digits))
    else:
        thresh = digits
        message = ("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
    if abs(expected - computed) > thresh:
        raise TestComparisonError(message)
    if math.isnan(computed):
        message += "\tprobably because the computed value is nan."
        raise TestComparisonError(message)
    if verbose >= 1:
        _success(label)


def compare_integers(expected, computed, label, verbose=1):
    """Function to compare two integers. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if (expected != computed):
        message = "\t{}: computed value ({}) does not match ({}).".format(label, computed, expected)
        raise TestComparisonError(message)
    if verbose >= 1:
        _success(label)


def compare_strings(expected, computed, label, verbose=1):
    """Function to compare two strings. Prints :py:func:`util.success`
    when string *computed* exactly matches string *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if(expected != computed):
        message = "\t%s: computed value (%s) does not match (%s)." % (label, computed, expected)
        raise TestComparisonError(message)
    if verbose >= 1:
        _success(label)


def compare_matrices(expected, computed, digits, label, verbose=1):
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
        raise TestComparisonError('compare_matrices failed')
    if verbose >= 1:
        _success(label)


def compare_dicts(expected, computed, tol, label, forgive=None, verbose=1):
    """Compares dictionaries `computed` to `expected` using DeepDiff Float
    comparisons made to `tol` significant decimal places. Note that a clean
    DeepDiff returns {}, which evaluates to False, hence the compare_integers.
    Keys in `forgive` may change between `expected` and `computed` without
    triggering failure.

    """
    import collections
    import deepdiff
    import pprint

    if forgive is None:
        forgive = []
    forgiven = collections.defaultdict(dict)

    ans = deepdiff.DeepDiff(expected, computed, significant_digits=tol, verbose_level=2)

    for category in list(ans):
        for key in list(ans[category]):
            for fg in forgive:
                fgsig = "root['" + fg + "']"
                if key.startswith(fgsig):
                    forgiven[category][key] = ans[category].pop(key)
        if not ans[category]:
            del ans[category]

    clean = not bool(ans)
    if not clean:
        pprint.pprint(ans)
    if verbose >= 2:
        pprint.pprint(forgiven)
    return compare_integers(True, clean, label, verbose=verbose)


def compare_molrecs(expected, computed, tol, label, forgive=None, verbose=1):
    """Function to compare Molecule dictionaries. Prints
    :py:func:`util.success` when elements of `computed` match elements of
    `expected` to `tol` number of digits (for float arrays).

    """
    import copy
    import numpy as np
    from .align import B787

    thresh = 10 ** -tol if tol >= 1 else tol

    # Need to manipulate the dictionaries a bit, so hold values
    xptd = copy.deepcopy(expected)
    cptd = copy.deepcopy(computed)

    # deepdiff can't cope with np.int type
    #   https://github.com/seperman/deepdiff/issues/97
    xptd['elez'] = [int(z) for z in xptd['elez']]
    cptd['elez'] = [int(z) for z in cptd['elez']]
    xptd['elea'] = [int(a) for a in xptd['elea']]
    cptd['elea'] = [int(a) for a in cptd['elea']]
    # deepdiff w/py27 complains about unicode type and val errors
    xptd['elem'] = [str(e) for e in xptd['elem']]
    cptd['elem'] = [str(e) for e in cptd['elem']]
    xptd['elbl'] = [str(l) for l in xptd['elbl']]
    cptd['elbl'] = [str(l) for l in cptd['elbl']]
    if 'fix_symmetry' in xptd:
        xptd['fix_symmetry'] = str(xptd['fix_symmetry'])
    if 'fix_symmetry' in cptd:
        cptd['fix_symmetry'] = str(cptd['fix_symmetry'])
    xptd['units'] = str(xptd['units'])
    cptd['units'] = str(cptd['units'])
    # and about int vs long errors
    xptd['molecular_multiplicity'] = int(xptd['molecular_multiplicity'])
    cptd['molecular_multiplicity'] = int(cptd['molecular_multiplicity'])
    xptd['fragment_multiplicities'] = [(m if m is None else int(m)) for m in xptd['fragment_multiplicities']]
    cptd['fragment_multiplicities'] = [(m if m is None else int(m)) for m in cptd['fragment_multiplicities']]

    # can't just expect geometries to match, so we'll align them, check that
    #   they overlap and that the translation/rotation arrays jibe with
    #   fix_com/orientation, then attach the oriented geom to computed before the
    #   recursive dict comparison.
    cgeom = cptd['geom'].reshape((-1, 3))
    rmsd, mill = B787(rgeom=xptd['geom'].reshape((-1, 3)),
                      cgeom=cgeom,
                      runiq=None,
                      cuniq=None,
                      atoms_map=True,
                      mols_align=True,
                      run_mirror=False,
                      verbose=0)
    if cptd['fix_com']:
        compare_integers(1, np.allclose(np.zeros((3)), mill.shift, atol=thresh), 'null shift', verbose=verbose)
    if cptd['fix_orientation']:
        compare_integers(1, np.allclose(np.identity(3), mill.rotation, atol=thresh), 'null rotation', verbose=verbose)
    ageom = mill.align_coordinates(cgeom)
    cptd['geom'] = ageom.reshape((-1))

    compare_dicts(xptd, cptd, tol, label, forgive=forgive, verbose=verbose)


def compare_molrecs_simple(expected, computed, tol, label, forgive=None, verbose=1):
    """Function to compare Molecule dictionaries. Prints
    :py:func:`util.success` when elements of `computed` match elements of
    `expected` to `tol` number of digits (for float arrays).

    """
    import copy
    import numpy as np
    from .align import B787

    thresh = 10 ** -tol if tol >= 1 else tol

    # Need to manipulate the dictionaries a bit, so hold values
    hold_e_geom = np.copy(expected['geom'])
    hold_e_elez = np.copy(expected['elez'])
    hold_e_elea = np.copy(expected['elea'])
    hold_c_geom = np.copy(computed['geom'])
    hold_c_elez = np.copy(computed['elez'])
    hold_c_elea = np.copy(computed['elea'])

    # deepdiff can't cope with np.int type
    #   https://github.com/seperman/deepdiff/issues/97
    expected['elez'] = [int(z) for z in expected['elez']]
    computed['elez'] = [int(z) for z in computed['elez']]
    expected['elea'] = [int(a) for a in expected['elea']]
    computed['elea'] = [int(a) for a in computed['elea']]

    # can't just expect geometries to match, so we'll align them, check that
    #   they overlap and that the translation/rotation arrays jibe with
    #   fix_com/orientation, then attach the oriented geom to computed before the
    #   recursive dict comparison.
    cgeom = computed['geom'].reshape((-1, 3))
    rmsd, mill = B787(rgeom=expected['geom'].reshape((-1, 3)),
                      cgeom=cgeom,
                      runiq=None,
                      cuniq=None,
                      atoms_map=True,
                      mols_align=True,
                      run_mirror=False,
                      verbose=0)
    if computed['fix_com']:
        compare_integers(1, np.allclose(np.zeros((3)), mill.shift, atol=thresh), 'null shift', verbose=verbose)
    if computed['fix_orientation']:
        compare_integers(1, np.allclose(np.identity(3), mill.rotation, atol=thresh), 'null rotation', verbose=verbose)
    ageom = mill.align_coordinates(cgeom)
    computed['geom'] = ageom.reshape((-1))

    compare_dicts(expected, computed, tol, label, forgive=forgive, verbose=verbose)

    # Replace values so function is const
    expected['geom'] = hold_e_geom
    expected['elez'] = hold_e_elez
    expected['elea'] = hold_e_elea
    computed['geom'] = hold_c_geom
    computed['elez'] = hold_c_elez
    computed['elea'] = hold_c_elea


def compare_arrays(expected, computed, digits, label, verbose=1):
    """Function to compare two numpy arrays. Prints :py:func:`util.success`
    when elements of vector *computed* match elements of vector *expected* to
    number of *digits*. Performs a system exit on failure to match symmetry
    structure, dimension, or element values. Used in input files in the test suite.

    """
    import numpy as np

    try:
        shape1 = expected.shape
        shape2 = computed.shape
    except:
        raise TestComparisonError("Input objects do not have a shape attribute.")

    if shape1 != shape2:
        TestComparisonError("Input shapes do not match.")

    tol = 10 ** (-digits)
    if not np.allclose(expected, computed, atol=tol):
        message = "\tArray difference norm is %12.6f." % np.linalg.norm(expected - computed)
        raise TestComparisonError(message)
    if verbose >= 1:
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

