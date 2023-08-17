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

"""Catch NumPy based import errors and suggest solutions.
Define some prettyprint formats.
Collect physical constants for py-side from QCElemental at CODATA consistent with c-side header.
"""

__all__ = [
    "constants",
    "pp",
    "nppp",
    "nppp10",
]

# NumPy import
try:
    import numpy as np
except ImportError:
    msg = """
    NumPy is a runtime requirement for Psi4. Please install NumPy to proceed.

    NumPy installation with a package manager can be accomplished by the following lines:
        - conda install numpy
        - sudo yum install numpy
        - sudo apt-get install python-numpy
        - brew install numpy
    """
    raise ImportError(msg)

# printing and logging formatting niceties
import pprint
from functools import partial
import numpy as np
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
nppp = partial(np.array_str, max_line_width=120, precision=8, suppress_small=True)
nppp10 = partial(np.array_str, max_line_width=120, precision=10, suppress_small=True)
del np, partial, pprint

# ensure Psi4 py-side constants are fixed at CODATA 2014, regardless of qcel default
import qcelemental as qcel
constants = qcel.PhysicalConstantsContext("CODATA2014")
del qcel
