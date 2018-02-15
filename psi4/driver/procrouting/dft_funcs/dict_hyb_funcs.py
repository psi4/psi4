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
"""
List of hybrid functionals
"""

mn15 = {
    "name": "MN15",
    "x_functionals": {"HYB_MGGA_X_MN15": {}},
    "c_functionals": {"MGGA_C_MN15": {}},
}

sogga11_x = {
    "name": "SOGGA11-X",
    "x_functionals": {"HYB_GGA_X_SOGGA11_X": {}},
    "c_functionals": {"GGA_C_SOGGA11_X": {}},
}

mn12_sx = {
    "name": "MN12_SX",
    "x_functionals": {"HYB_MGGA_X_MN12_SX": {}},
    "c_functionals": {"MGGA_C_MN12_SX": {}},
}

functional_list = {
    "TEST-MN15": mn15,
    "TEST-SOGGA11-X": sogga11_x,
    "TEST-MN12-SX": mn12_sx,
}
