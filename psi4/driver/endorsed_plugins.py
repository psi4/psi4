#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

"""Import plugins eligible to be accessible in input files if detected."""

__all__ = []

import sys

try:
    import v2rdm_casscf
except ImportError:
    pass

try:
    import gpu_dfcc
except ImportError:
    pass

try:
    import forte
except ImportError:
    pass

try:
    import snsmp2
except ImportError as e:
    if 'scipy' in e.msg:
        raise ImportError("""Psi4 plugin 'snsmp2' available, but scipy missing. Try `conda install scipy` or `pip install scipy`.""")
    else:
        pass

try:
    import resp
except ImportError:
    pass

try:
    import adcc
except ImportError:
    pass

try:
    import psi4fockci
except ImportError:
    pass

try:
    import cct3
except ImportError:
    pass

