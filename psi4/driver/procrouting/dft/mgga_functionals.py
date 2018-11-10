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
List of meta-GGA functionals
"""

import copy

funcs = []

funcs.append({
    "name": "M06-L",
    "alias": ["M06L"],
    "x_functionals": {
        "MGGA_X_M06_L": {}
    },
    "c_functionals": {
        "MGGA_C_M06_L": {}
    },
    "description": '    M06-L Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101, 2006\n',
})

funcs.append({
    "name": "M11-L",
    "alias": ["M11L"],
    "x_functionals": {
        "MGGA_X_M11_L": {}
    },
    "c_functionals": {
        "MGGA_C_M11_L": {}
    },
    "description": '    M11-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Phys. Chem. Lett. 3, 117-124, 2012\n',
})

funcs.append({
    "name": "MN12-L",
    "alias": ["MN12L"],
    "x_functionals": {
        "MGGA_X_MN12_L": {}
    },
    "c_functionals": {
        "MGGA_C_MN12_L": {}
    },
    "description": '    MN12-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, Phys. Chem. Chem. Phys. 14, 13171-13174, 2012\n',
})

funcs.append({
    "name": "MN15-L",
    "alias": ["MN15L"],
    "x_functionals": {
        "MGGA_X_MN15_L": {}
    },
    "c_functionals": {
        "MGGA_C_MN15_L": {}
    },
    "description": '    MN15-L Meta-GGA XC Functional\n',
    "citation": '    H. S. Yu, X. He, and D. G. Truhlar, J. Chem. Theory Comput. 12, 1280-1293, 2016\n',
})

funcs.append({
    "name": "mGGA_MS0",
    "alias": ["MGGA-MS0"],
    "x_functionals": {
        "MGGA_X_MS0": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '    MGGA_MS0 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 137, 051101, 2012\n',
})

funcs.append({
    "name": "mGGA_MS1",
    "alias": ["MGGA-MS1"],
    "x_functionals": {
        "MGGA_X_MS1": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '    MGGA_MS1 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
})

funcs.append({
    "name": "mGGA_MS2",
    "alias": ["MGGA-MS2"],
    "x_functionals": {
        "MGGA_X_MS2": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '    MGGA_MS2 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
})

funcs.append({
    "name": "mGGA_MVS",
    "alias": ["MGGA-MVS"],
    "x_functionals": {
        "MGGA_X_MVS": {}
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '    MGGA_MVS Meta-GGA XC Functional\n',
    "citation": '    J. Sun, J.P. Perdew, A. Ruzsinszky, Proc. Natl. Acad. Sci. USA 112, 685, 2015\n',
})

funcs.append({
    "name": "TPSS",
    "x_functionals": {
        "MGGA_X_TPSS": {}
    },
    "c_functionals": {
        "MGGA_C_TPSS": {}
    },
    "description": '    TPSS Meta-GGA XC Functional\n',
    "citation": '    J. Tao, et al., Phys. Rev. Lett., 91, 146401, 2003\n',
})

funcs.append({
    "name": "revTPSS",
    "x_functionals": {
        "MGGA_X_REVTPSS": {}
    },
    "c_functionals": {
        "MGGA_C_REVTPSS": {}
    },
    "description": '    revised TPSS Meta-GGA XC Functional\n',
    "citation": '    J. Sun  et. al., Phys. Rev. B, 84, 035117, 2011\n',
})

funcs.append({
    "name": "PKZB",
    "x_functionals": {
        "MGGA_X_PKZB": {}
    },
    "c_functionals": {
        "MGGA_C_PKZB": {}
    },
    "description": '    PKZB Meta-GGA XC Functional\n',
    "citation": '    J.P. Perdew, S. Kurth, A. Zupan, P. Blaha, Phys. Rev. Lett. 82, 2544, 1999\n',
})

funcs.append({
    "name": "VSXC",
    "x_functionals": {
        "MGGA_X_GVT4": {}
    },
    "c_functionals": {
        "MGGA_C_VSXC": {}
    },
    "description": '    VSXC Meta-GGA XC Functional\n',
    "citation": '    T.V. Voorhis, G.E. Scuseria, J. Chem. Phys. 109, 400, 1998\n',
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
