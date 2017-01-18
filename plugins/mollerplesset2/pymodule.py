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

import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from molutil import *
import p4util
from p4util.exceptions import *


def run_plugin_mp2(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mollerplesset2 can be called via :py:func:`~driver.energy`.

    >>> energy('mollerplesset2')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_local_option('MOLLERPLESSET2', 'PRINT', 1)
    scf_wfn = scf_helper(lowername)

    # Need to semicanonicalize the ROHF orbitals
    if psi4.get_global_option('REFERENCE') == 'ROHF':
        scf_wfn.semicanonicalize()

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.get_option('SCF', 'SCF_TYPE'), scf_wfn)

    #psi4.set_legacy_wavefunction(scf_wfn)
    returnvalue = psi4.plugin('mollerplesset2.so', scf_wfn)

    return returnvalue


# Integration with driver routines
procedures['energy']['mollerplesset2'] = run_plugin_mp2


def exampleFN():
    # Your Python code goes here
    pass
