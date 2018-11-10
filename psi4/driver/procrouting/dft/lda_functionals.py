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
"""
List of LDA functionals
"""

import copy

funcs = []

funcs.append({
    "name": "SVWN",
    "x_functionals": {
        "LDA_X": {}
    },
    "c_functionals": {
        "LDA_C_VWN_RPA": {}
    },
    "citation": '    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n',
    "description": '    SVWN3 (RPA) LSDA Functional\n',
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
