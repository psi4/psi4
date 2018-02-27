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
List of meta-GGA functionals
"""

import copy

funcs = []

funcs.append({
    "name": "M06-L",
    "alias": ["M06L"],
    "description": '    M06-L Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101, 2006\n',
    "x_functionals": {
        "MGGA_X_M06_L": {}
    },
    "c_functionals": {
        "MGGA_C_M06_L": {}
    },
})

funcs.append({
    "name": "M11-L",
    "alias": ["M11L"],
    "description": '    M11-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Phys. Chem. Lett. 3, 117-124, 2012\n',
    "x_functionals": {
        "MGGA_X_M11_L": {}
    },
    "c_functionals": {
        "MGGA_C_M11_L": {}
    },
})

funcs.append({
    "name": "MN12-L",
    "alias": ["MN12L"],
    "description": '    MN12-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, Phys. Chem. Chem. Phys. 14, 13171-13174, 2012\n',
    "x_functionals": {
        "MGGA_X_MN12_L": {}
    },
    "c_functionals": {
        "MGGA_C_MN12_L": {}
    },
})

funcs.append({
    "name": "MN15-L",
    "alias": ["MN15L"],
    "description": '    MN15-L Meta-GGA XC Functional\n',
    "citation": '    H. S. Yu, X. He, and D. G. Truhlar, J. Chem. Theory Comput. 12, 1280-1293, 2016\n',
    "x_functionals": {
        "MGGA_X_MN15_L": {}
    },
    "c_functionals": {
        "MGGA_C_MN15_L": {}
    },
})

funcs.append({
    "name": "mGGA_MS0",
    "alias": ["MGGA-MS0"],
    "description": '    MGGA_MS0 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 137, 051101, 2012\n',
    "x_functionals": {
        "MGGA_X_MS0": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
})

funcs.append({
    "name": "mGGA_MS1",
    "alias": ["MGGA-MS1"],
    "description": '    MGGA_MS1 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
    "x_functionals": {
        "MGGA_X_MS1": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
})

funcs.append({
    "name": "mGGA_MS2",
    "alias": ["MGGA-MS2"],
    "description": '    MGGA_MS2 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
    "x_functionals": {
        "MGGA_X_MS2": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
})

funcs.append({
    "name":
    "TPSS",
    "description":
    '    TPSS Meta-GGA XC Functional\n',
    "citation":
    '   J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria, Phys. Rev. Lett., 91, 146401, 2003\n',
    "x_functionals": {
        "MGGA_X_TPSS": {}
    },
    "c_functionals": {
        "MGGA_C_TPSS": {}
    },
})

funcs.append({
    "name": "revTPSS",
    "description": '    revised TPSS Meta-GGA XC Functional\n',
    "citation": '   J. Sun  et. al., Phys. Rev. B, 84, 035117, 2011\n',
    "x_functionals": {
        "MGGA_X_REVTPSS": {}
    },
    "c_functionals": {
        "MGGA_C_REVTPSS": {}
    },
})

functional_list = {}
for functional in funcs:
    if "alias" in functional.keys():
        alias = functional.pop("alias")
        for a in alias:
            functional_list["TEST-" + a] = functional
    functional_list["TEST-" + functional["name"].upper()] = functional
