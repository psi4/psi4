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
List of GGA functionals
"""

import copy

funcs = []

funcs.append({
    "name": "BLYP",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_LYP": {}},
    "citation": '    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n',
    "description": '    BLYP GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "BOP",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_OP_B88": {}},
    "citation": '    T. Tsuneda et. al., J. Chem. Phys. 110, 10664-10678, 1999\n',
    "description": '    BOP GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "B86BPBE",
    "x_functionals": {"GGA_X_B86_MGC": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n',
    "description": '    B86BPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PW86PBE",
    "x_functionals": {"GGA_X_PW86": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n',
    "description": '    PW86PBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PBE",
    "x_functionals": {"GGA_X_PBE": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n',
    "description": '    PBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "RPBE",
    "x_functionals": {"GGA_X_RPBE": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    B. Hammer, L. B. Hansen, and J. K. Nørskov, Phys. Rev. B, 59, 7413, 1999\n',
    "description": '    RPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "revPBE",
    "x_functionals": {"GGA_X_PBE_R": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    Y. Zhang and W. Yang, Phys. Rev. Lett. 80, 890, 1998\n',
    "description": '    revPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PW91",
    "x_functionals": {"GGA_X_PW91": {}},
    "c_functionals": {"GGA_C_Pw91": {}},
    "citation": '    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n',
    "description": '    PW91 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "mPWPW",
    "x_functionals": {"GGA_X_mPW91": {}},
    "c_functionals": {"GGA_C_Pw91": {}},
    "citation": '    C. Adamo, V. Barone, J. Chem. Phys., 108, 664, 1998\n',
    "description": '    mPWPW GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "BP86",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_P86": {}},
    "citation": '   A. D. Becke, Phys. Rev. A, 38, 3098-3100, 1988\n   J. P. Perdew, Phys. Rev. B, 33, 8822, 1986\n',
    "description": '    BP86 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "FT97",
    "x_functionals": {"GGA_X_FT97_B": {}},
    "c_functionals": {"GGA_C_FT97": {}},
    "citation": '    M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997\n',
    "description": '   FT97 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "SOGGA",
    "x_functionals": {"GGA_X_SOGGA": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    Y. Zhao, and D. G. Truhlar, J. Chem. Phys. 128, 184109, 2008\n' +
                '    J. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865-3868, 1996\n',
    "description": '   SOGGA Exchange + PBE Correlation Functional\n',
})

funcs.append({
    "name": "SOGGA11",
    "x_functionals": {"GGA_X_SOGGA11": {}},
    "c_functionals": {"GGA_C_SOGGA11": {}},
    "citation": '    R. Peverati, Y. Zhao, and D. G. Truhlar, J. Phys. Chem. Lett. 2, 1991-1997, 2011\n',
    "description": '   SOGGA11 Exchange-Correlation Functional\n',
})

functional_list = {}
for functional in funcs:
    if "alias" in functional.keys():
        alias = functional.pop("alias")
        for a in alias:
            functional_list["TEST-" + a] = functional
    functional_list["TEST-" + functional["name"].upper()] = functional
