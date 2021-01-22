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
"""Module with utility functions for use in input files."""

import os
import re
import sys
import warnings

from psi4 import core
from psi4.driver.procrouting import *
from .exceptions import ValidationError

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
        core.set_global_option('CUBEPROP_TASKS', ['ORBITALS'])

    if ((core.get_global_option('INTEGRAL_PACKAGE') == 'ERD') and ('ESP' in core.get_global_option('CUBEPROP_TASKS'))):
        raise ValidationError('INTEGRAL_PACKAGE ERD does not play nicely with electrostatic potential, so stopping.')

    cp = core.CubeProperties(wfn)
    cp.compute_properties()


def set_memory(inputval, execute=True, quiet=False):
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

    >>> # Good examples
    >>> psi4.set_memory(800000000)        # 800000000
    >>> psi4.set_memory(2004088624.9)     # 2004088624
    >>> psi4.set_memory(1.0e9)            # 1000000000
    >>> psi4.set_memory('600 mb')         # 600000000
    >>> psi4.set_memory('600.0 MiB')      # 629145600
    >>> psi4.set_memory('.6 Gb')          # 600000000
    >>> psi4.set_memory(' 100000000kB ')  # 100000000000
    >>> psi4.set_memory('2 eb')           # 2000000000000000000

    >>> # Bad examples
    >>> psi4.set_memory({})         # odd type
    >>> psi4.set_memory('')         # no info
    >>> psi4.set_memory("8 dimms")  # unacceptable units
    >>> psi4.set_memory("1e5 gb")   # string w/ exponent
    >>> psi4.set_memory("5e5")      # string w/o units
    >>> psi4.set_memory(2000)       # mem too small
    >>> psi4.set_memory(-5e5)       # negative (and too small)

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
        raise ValidationError(
            """set_memory(): Requested {:.3} MiB ({:.3} MB); minimum 250 MiB (263 MB). Please, sir, I want some more."""
            .format(memory_amount / 1024**2, memory_amount / 1000**2))

    if execute:
        core.set_memory_bytes(memory_amount, quiet)
    return memory_amount


def get_memory():
    """Function to return the total memory allocation."""
    return core.get_memory()


def copy_file_to_scratch(filename, prefix, namespace, unit, move=False):
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

    cp = '/bin/cp'
    if move:
        cp = '/bin/mv'

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


def copy_file_from_scratch(filename, prefix, namespace, unit, move=False):
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

    cp = '/bin/cp'
    if move:
        cp = '/bin/mv'

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
