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

"""Module with utility functions for use in input files."""
from __future__ import division
import re
import sys
import os
import math
import numpy as np
from .exceptions import *


def oeprop(wfn, *args, **kwargs):
    """Evaluate one-electron properties.

    :returns: None

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to compute properties

    How to specify args, which are actually the most important

    :type title: string
    :param title: label prepended to all psivars computed

    :examples:

    >>> # [1] Moments with specific label
    >>> E, wfn = energy('hf', return_wfn=True)
    >>> oeprop(wfn, 'DIPOLE', 'QUADRUPOLE', title='H3O+ SCF')

    """
    oe = core.OEProp(wfn)
    if 'title' in kwargs:
        oe.set_title(kwargs['title'])
    for prop in args:
        oe.add(prop)
    oe.compute()


def cubeprop(wfn, **kwargs):
    """Evaluate properties on a grid and generate cube files.

    .. versionadded:: 0.5
       *wfn* parameter passed explicitly

    :returns: None

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, basis, orbitals from which to generate cube files

    :examples:

    >>> # [1] Cube files for all orbitals
    >>> E, wfn = energy('b3lyp', return_wfn=True)
    >>> cubeprop(wfn)

    >>> # [2] Cube files for density (alpha, beta, total, spin) and four orbitals
    >>> #     (two alpha, two beta)
    >>> set cubeprop_tasks ['orbitals', 'density']
    >>> set cubeprop_orbitals [5, 6, -5, -6]
    >>> E, wfn = energy('scf', return_wfn=True)
    >>> cubeprop(wfn)

    """
    # By default compute the orbitals
    if not core.has_global_option_changed('CUBEPROP_TASKS'):
        core.set_global_option('CUBEPROP_TASKS',['ORBITALS'])

    if ((core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and
        ('ESP' in core.get_global_option('CUBEPROP_TASKS'))):
        raise ValidationError('INTEGRAL_PACKAGE ERD does not play nicely with electrostatic potential, so stopping.')

    cp = core.CubeProperties(wfn)
    cp.compute_properties()


def set_memory(inputval, execute=True):
    """Function to reset the total memory allocation. Takes memory value
    *inputval* as type int, float, or str; int and float are taken literally
    as bytes to be set, string taken as a unit-containing value (e.g., 30 mb)
    which is case-insensitive. Set *execute* to False to interpret *inputval*
    without setting in Psi4 core.

    :returns: *memory_amount* (float) Number of bytes of memory set

    :raises: ValidationError when <500MiB or disallowed type or misformatted

    :examples:

    >>> # [1] Passing absolute number of bytes
    >>> psi4.set_memory(600000000)
    >>> psi4.get_memory()
    Out[1]: 600000000L

    >>> # [2] Passing memory value as string with units
    >>> psi4.set_memory('30 GB')
    >>> psi4.get_memory()
    Out[2]: 30000000000L

    :good examples:

    800000000         # 800000000
    2004088624.9      # 2004088624
    1.0e9             # 1000000000
    '600 mb'          # 600000000
    '600.0 MiB'       # 629145600
    '.6 Gb'           # 600000000
    ' 100000000kB '   # 100000000000
    '2 eb'            # 2000000000000000000

    :bad examples:

    {}         # odd type
    ''         # no info
    "8 dimms"  # unacceptable units
    "1e5 gb"   # string w/ exponent
    "5e5"      # string w/o units
    2000       # mem too small
    -5e5       # negative (and too small)

    """
    # Handle memory given in bytes directly (int or float)
    if isinstance(inputval, (int, float)):
        val = inputval
        units = ''
    # Handle memory given as a string
    elif isinstance(inputval, str):
        memory_string = re.compile(r'^\s*(\d*\.?\d+)\s*([KMGTPBE]i?B)\s*$', re.IGNORECASE)
        matchobj = re.search(memory_string, inputval)
        if matchobj:
            val = float(matchobj.group(1))
            units = matchobj.group(2)
        else:
            raise ValidationError("""Invalid memory specification: {}. Try 5e9 or '5 gb'.""".format(repr(inputval)))
    else:
        raise ValidationError("""Invalid type {} in memory specification: {}. Try 5e9 or '5 gb'.""".format(
            type(inputval), repr(inputval)))

    # Units decimal or binary?
    multiplier = 1000
    if "i" in units.lower():
        multiplier = 1024
        units = units.lower().replace("i", "").upper()

    # Build conversion factor, convert units
    unit_list = ["", "KB", "MB", "GB", "TB", "PB", "EB"]
    mult = 1
    for unit in unit_list:
        if units.upper() == unit:
            break
        mult *= multiplier

    memory_amount = int(val * mult)

    # Check minimum memory requirement
    min_mem_allowed = 262144000
    if memory_amount < min_mem_allowed:
        raise ValidationError("""set_memory(): Requested {:.3} MiB ({:.3} MB); minimum 250 MiB (263 MB). Please, sir, I want some more.""".format(
                memory_amount / 1024 ** 2, memory_amount / 1000 ** 2))

    if execute:
        core.set_memory_bytes(memory_amount)
    return memory_amount

def get_memory():
    """Function to return the total memory allocation."""
    return core.get_memory()


def success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    msg = '\t{0:.<66}PASSED'.format(label)
    print(msg)
    sys.stdout.flush()
    core.print_out(msg + '\n')


# Test functions
def compare_values(expected, computed, digits, label, exitonfail=True):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*
    (or to *digits* itself when *digits* < 1 e.g. digits=0.04). Performs
    a system exit on failure unless *exitonfail* False, in which case
    returns error message. Used in input files in the test suite.

    """
    if digits > 1:
        thresh = 10 ** -digits
        message = ("\t%s: computed value (%.*f) does not match (%.*f) to %d digits." % (label, digits+1, computed, digits+1, expected, digits))
    else:
        thresh = digits
        message = ("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
    if abs(expected - computed) > thresh:
        print(message)
        if exitonfail:
            raise TestComparisonError(message)
    if math.isnan(computed):
        print(message)
        print("\tprobably because the computed value is nan.")
        if exitonfail:
            raise TestComparisonError(message)
    success(label)
    return True


def compare_integers(expected, computed, label):
    """Function to compare two integers. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if (expected != computed):
        message = ("\t%s: computed value (%d) does not match (%d)." % (label, computed, expected))
        raise TestComparisonError(message)
    success(label)
    return True


def compare_strings(expected, computed, label):
    """Function to compare two strings. Prints :py:func:`util.success`
    when string *computed* exactly matches string *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if(expected != computed):
        message = ("\t%s: computed value (%s) does not match (%s)." % (label, computed, expected))
        raise TestComparisonError(message)
    success(label)
    return True


def compare_matrices(expected, computed, digits, label):
    """Function to compare two matrices. Prints :py:func:`util.success`
    when elements of matrix *computed* match elements of matrix *expected* to
    number of *digits*. Performs a system exit on failure to match symmetry
    structure, dimensions, or element values. Used in input files in the test suite.

    """
    if (expected.nirrep() != computed.nirrep()):
        message = ("\t%s has %d irreps, but %s has %d\n." % (expected.name(), expected.nirrep(), computed.name(), computed.nirrep()))
        raise TestComparisonError(message)
    if (expected.symmetry() != computed.symmetry()):
        message = ("\t%s has %d symmetry, but %s has %d\n." % (expected.name(), expected.symmetry(), computed.name(), computed.symmetry()))
        raise TestComparisonError(message)
    nirreps = expected.nirrep()
    symmetry = expected.symmetry()
    for irrep in range(nirreps):
        if(expected.rows(irrep) != computed.rows(irrep)):
            message = ("\t%s has %d rows in irrep %d, but %s has %d\n." % (expected.name(), expected.rows(irrep), irrep, computed.name(), computed.rows(irrep)))
            raise TestComparisonError(message)
        if(expected.cols(irrep ^ symmetry) != computed.cols(irrep ^ symmetry)):
            message = ("\t%s has %d columns in irrep, but %s has %d\n." % (expected.name(), expected.cols(irrep), irrep, computed.name(), computed.cols(irrep)))
            raise TestComparisonError(message)
        rows = expected.rows(irrep)
        cols = expected.cols(irrep ^ symmetry)
        failed = 0
        for row in range(rows):
            for col in range(cols):
                if(abs(expected.get(irrep, row, col) - computed.get(irrep, row, col)) > 10 ** (-digits)):
                    print("\t%s: computed value (%s) does not match (%s)." % (label, expected.get(irrep, row, col), computed.get(irrep, row, col)))
                    failed = 1
                    break

        if(failed):
            print("Check your output file for reporting of the matrices.")
            core.print_out("The Failed Test Matrices\n")
            core.print_out("Computed Matrix (2nd matrix passed in)\n")
            computed.print_out()
            core.print_out("Expected Matrix (1st matrix passed in)\n")
            expected.print_out()
            raise TestComparisonError("\n")
    success(label)
    return True


def compare_vectors(expected, computed, digits, label):
    """Function to compare two vectors. Prints :py:func:`util.success`
    when elements of vector *computed* match elements of vector *expected* to
    number of *digits*. Performs a system exit on failure to match symmetry
    structure, dimension, or element values. Used in input files in the test suite.

    """
    if (expected.nirrep() != computed.nirrep()):
        message = ("\t%s has %d irreps, but %s has %d\n." % (expected.name(), expected.nirrep(), computed.name(), computed.nirrep()))
        raise TestComparisonError(message)
    nirreps = expected.nirrep()
    for irrep in range(nirreps):
        if(expected.dim(irrep) != computed.dim(irrep)):
            message = ("\tThe reference has %d entries in irrep %d, but the computed vector has %d\n." % (expected.dim(irrep), irrep, computed.dim(irrep)))
            raise TestComparisonError(message)
        dim = expected.dim(irrep)
        failed = 0
        for entry in range(dim):
            if(abs(expected.get(irrep, entry) - computed.get(irrep, entry)) > 10 ** (-digits)):
                failed = 1
                break

        if(failed):
            core.print_out("The computed vector\n")
            computed.print_out()
            core.print_out("The reference vector\n")
            expected.print_out()
            message = ("\t%s: computed value (%s) does not match (%s)." % (label, computed.get(irrep, entry), expected.get(irrep, entry)))
            raise TestComparisonError(message)
    success(label)
    return True


def compare_arrays(expected, computed, digits, label):
    """Function to compare two numpy arrays. Prints :py:func:`util.success`
    when elements of vector *computed* match elements of vector *expected* to
    number of *digits*. Performs a system exit on failure to match symmetry
    structure, dimension, or element values. Used in input files in the test suite.

    """

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
    success(label)
    return True


def compare_cubes(expected, computed, label):
    """Function to compare two cube files. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    # Skip the first six elemets which are just labels
    evec = [float(k) for k in expected.split()[6:]]
    cvec = [float(k) for k in computed.split()[6:]]
    if len(evec) == len(cvec):
        for n in range(len(evec)):
            if (math.fabs(evec[n]-cvec[n]) > 1.0e-4):
                message = ("\t%s: computed cube file does not match expected cube file." % label)
                raise TestComparisonError(message)
    else:
        message = ("\t%s: computed cube file does not match expected cube file." % (label, computed, expected))
        raise TestComparisonError(message)
    success(label)
    return True


def copy_file_to_scratch(filename, prefix, namespace, unit, move = False):

    """Function to move file into scratch with correct naming
    convention.

    Arguments:

    @arg filename  full path to file
    @arg prefix    computation prefix, usually 'psi'
    @arg namespace context namespace, usually molecule name
    @arg unit      unit number, e.g. 32
    @arg move      copy or move? (default copy)

    Example:

    Assume PID is 12345 and SCRATCH is /scratch/parrish/

    copy_file_to_scratch('temp', 'psi', 'h2o', 32):
        -cp ./temp /scratch/parrish/psi.12345.h2o.32
    copy_file_to_scratch('/tmp/temp', 'psi', 'h2o', 32):
        -cp /tmp/temp /scratch/parrish/psi.12345.h2o.32
    copy_file_to_scratch('/tmp/temp', 'psi', '', 32):
        -cp /tmp/temp /scratch/parrish/psi.12345.32
    copy_file_to_scratch('/tmp/temp', 'psi', '', 32, True):
        -mv /tmp/temp /scratch/parrish/psi.12345.32

    """

    pid = str(os.getpid())
    scratch = core.IOManager.shared_object().get_file_path(int(unit))

    cp = '/bin/cp';
    if move:
        cp = '/bin/mv';

    unit = str(unit)

    target = ''
    target += prefix
    target += '.'
    target += pid
    if len(namespace):
        target += '.'
        target += namespace
    target += '.'
    target += unit

    command = ('%s %s %s/%s' % (cp, filename, scratch, target))

    os.system(command)
    #print command

def copy_file_from_scratch(filename, prefix, namespace, unit, move = False):

    """Function to move file out of scratch with correct naming
    convention.

    Arguments:

    @arg filename  full path to target file
    @arg prefix    computation prefix, usually 'psi'
    @arg namespace context namespace, usually molecule name
    @arg unit      unit number, e.g. 32
    @arg move      copy or move? (default copy)

    Example:

    Assume PID is 12345 and SCRATCH is /scratch/parrish/

    copy_file_to_scratch('temp', 'psi', 'h2o', 32):
        -cp /scratch/parrish/psi.12345.h2o.32 .temp
    copy_file_to_scratch('/tmp/temp', 'psi', 'h2o', 32):
        -cp /scratch/parrish/psi.12345.h2o.32 /tmp/temp
    copy_file_to_scratch('/tmp/temp', 'psi', '', 32):
        -cp /scratch/parrish/psi.12345.32 /tmp/temp
    copy_file_to_scratch('/tmp/temp', 'psi', '', 32, True):
        -mv /scratch/parrish/psi.12345.32 /tmp/temp

    """

    pid = str(os.getpid())
    scratch = core.IOManager.shared_object().get_file_path(int(unit))

    cp = '/bin/cp';
    if move:
        cp = '/bin/mv';

    unit = str(unit)

    target = ''
    target += prefix
    target += '.'
    target += pid
    if len(namespace):
        target += '.'
        target += namespace
    target += '.'
    target += unit

    command = ('%s %s/%s %s' % (cp, scratch, target, filename))

    os.system(command)


def xml2dict(filename=None):
    """Read XML *filename* into nested OrderedDict-s. *filename* defaults to
    active CSX file.

    """
    import xmltodict as xd
    if filename is None:
        csx = os.path.splitext(core.outfile_name())[0] + '.csx'
    else:
        csx = filename
    with open(csx, 'r') as handle:
        csxdict = xd.parse(handle)

    return csxdict


def getFromDict(dataDict, mapList):
    return reduce(lambda d, k: d[k], mapList, dataDict)


def csx2endict():
    """Grabs the CSX file as a dictionary, encodes translation of PSI variables
    to XML blocks, gathers all available energies from CSX file into returned
    dictionary.

    """
    blockprefix = ['chemicalSemantics', 'molecularCalculation', 'quantumMechanics', 'singleReferenceState', 'singleDeterminant']
    blockmidfix = ['energies', 'energy']
    prefix = 'cs:'

    pv2xml = {
        'MP2 CORRELATION ENERGY': [['mp2'], 'correlation'],
        'MP2 SAME-SPIN CORRELATION ENERGY': [['mp2'], 'sameSpin correlation'],
        'HF TOTAL ENERGY': [['abinitioScf'], 'electronic'],
        'NUCLEAR REPULSION ENERGY': [['abinitioScf'], 'nuclearRepulsion'],
        'DFT FUNCTIONAL TOTAL ENERGY': [['dft'], 'dftFunctional'],
        'DFT TOTAL ENERGY': [['dft'], 'electronic'],
        'DOUBLE-HYBRID CORRECTION ENERGY': [['dft'], 'doubleHybrid correction'],
        'DISPERSION CORRECTION ENERGY': [['dft'], 'dispersion correction'],
    }

    csxdict = xml2dict()
    enedict = {}
    for pv, lpv in pv2xml.items():
        address = blockprefix + lpv[0] + blockmidfix
        indices = [prefix + bit for bit in address]
        try:
            qwer = getFromDict(csxdict, indices)
        except KeyError:
            continue
        for v in qwer:
            vv = v.values()
            if vv[0] == prefix + lpv[1]:
                enedict[pv] = float(vv[1])

    return enedict


def compare_csx():
    """Function to validate energies in CSX files against PSIvariables. Only
    active if write_csx flag on.

    """
    if 'csx4psi' in sys.modules.keys():
        if core.get_global_option('WRITE_CSX'):
            enedict = csx2endict()
            compare_integers(len(enedict) >= 2, True, 'CSX harvested')
            for pv, en in enedict.items():
                compare_values(core.get_variable(pv), en, 6, 'CSX ' + pv + ' ' + str(round(en, 4)))
