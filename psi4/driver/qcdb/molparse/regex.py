#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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

import re

NUCLEUS = r"""(?:
   (?P<gh1>@)|(?P<gh2>Gh\())?                # optional ghost: @stuff or Gh(stuff) ...
        (                                    # mandatory element: AEuser or Zuser
            (?P<label1>
                (?P<A>\d+)?                  # optional mass number, A
                (?P<E>[A-Z]{1,3})            # mandatory atomic symbol, E
                (?P<user1>(_\w+)|(\d+))?) |  # optional user label: Enumber or E_alphanum
            (?P<label2>
                (?P<Z>\d{1,3})               # mandatory atomic number, Z
                (?P<user2>(_\w+))?)          # optional user label: Z_alphanum
        )
        (?:@(?P<mass>\d+\.\d+))?             # optional mass value [u]
   (?(gh2)\)                                 # ... ghost
          )"""

NUMBER = r"""(
    (?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?) |   # .num with optional sign, exponent, wholenum
    (?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?) |   # num. with optional sign, exponent, decimals
    (?:[-+]?\d+(?:[DdEe][-+]?\d+)?)          # num with optional sign, exponent
         )"""

SEP = r"""[\t ,]+"""
ENDL = r"""[\t ,]*$"""

CHGMULT = r"""(?P<chg>""" + NUMBER + r')' + SEP + r"""(?P<mult>\d+)"""
CARTXYZ = r'(?P<x>' + NUMBER + r')' + SEP + r'(?P<y>' + NUMBER + r')' + SEP + r'(?P<z>' + NUMBER + r')'
