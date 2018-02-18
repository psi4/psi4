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


dldf = {
    "name": "dlDF",
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
    "x_functionals": {"HYB_MGGA_X_DLDF": {"use_libxc": True}},
    "c_functionals": {"MGGA_C_DLDF": {}},
}


dldf_d09 = copy.deepcopy(dldf)
dldf_d09["name"] = "dlDF+d09"
dldf_d09["dispersion"] = {"type": "das2009", "params": {"s6": 1.0}}


dldf_d10 = copy.deepcopy(dldf)
dldf_d10["name"] = "dlDF+d10"
dldf_d10["dispersion"] = {"type": "das2010", "params": {"s6": 1.0}}


m06_l = {
    "name": "M06-L",
    "description": '    M06-L Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101, 2006\n',
    "x_functionals": {"MGGA_X_M06_L": {}},
    "c_functionals": {"MGGA_C_M06_L": {}},
}


m11_l = {
    "name": "M11-L",
    "description": '    M11-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Phys. Chem. Lett. 3, 117-124, 2012\n',
    "x_functionals": {"MGGA_X_M11_L": {}},
    "c_functionals": {"MGGA_C_M11_L": {}},
}


mn12_l = {
    "name": "MN12-L",
    "description": '    MN12-L Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, Phys. Chem. Chem. Phys. 14, 13171-13174, 2012\n',
    "x_functionals": {"MGGA_X_MN12_L": {}},
    "c_functionals": {"MGGA_C_MN12_L": {}},
}


mn15_l = {
    "name": "MN15-L",
    "description": '    MN15-L Meta-GGA XC Functional\n',
    "citation": '    H. S. Yu, X. He, and D. G. Truhlar, J. Chem. Theory Comput. 12, 1280-1293, 2016\n',
    "x_functionals": {"MGGA_X_MN15_L": {}},
    "c_functionals": {"MGGA_C_MN15_L": {}},
}


mgga_ms0 = {
    "name": "MGGA_MS0",
    "description": '    MGGA_MS0 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 137, 051101, 2012\n',
    "x_functionals": {"MGGA_X_MS0": {}},
    "c_functionals": {"GGA_C_REGTPSS": {}},
}


mgga_ms1 = {
    "name": "MGGA_MS1",
    "description": '    MGGA_MS1 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
    "x_functionals": {"MGGA_X_MS1": {}},
    "c_functionals": {"GGA_C_REGTPSS": {}},
}


mgga_ms2 = {
    "name": "MGGA_MS2",
    "description": '    MGGA_MS2 Meta-GGA XC Functional\n',
    "citation": '    J. Sun et. al., J. Chem. Phys. 138, 044113, 2013\n',
    "x_functionals": {"MGGA_X_MS2": {}},
    "c_functionals": {"GGA_C_REGTPSS": {}},
}


tpss = {
    "name": "TPSS",
    "description": '    TPSS Meta-GGA XC Functional\n',
    "citation": '   J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria, Phys. Rev. Lett., 91, 146401, 2003\n',
    "x_functionals": {"MGGA_X_TPSS": {}},
    "c_functionals": {"MGGA_C_TPSS": {}},
}


revtpss = {
    "name": "TPSS",
    "description": '    revised TPSS Meta-GGA XC Functional\n',
    "citation": '   J. Sun  et. al., Phys. Rev. B, 84, 035117, 2011\n',
    "x_functionals": {"MGGA_X_REVTPSS": {}},
    "c_functionals": {"MGGA_C_REVTPSS": {}},
}


functional_list = {
    "TEST-DLDF": dldf,
    "TEST-DLDF+D09": dldf_d09,
    "TEST-DLDF+D": dldf_d10,
    "TEST-M06-L": m06_l,
    "TEST-M11-L": m11_l,
    "TEST-MN12-L": mn12_l,
    "TEST-MN15-L": mn15_l,
    "TEST-MGGA_MS0": mgga_ms0,
    "TEST-MGGA_MS1": mgga_ms1,
    "TEST-MGGA_MS2": mgga_ms2,
    "TEST-TPSS": tpss,
    "TEST-REVTPSS": revtpss,
    "TEST-ZLP":       {"name": "ZLP",       "xc_functionals": {"MGGA_XC_ZLP": {}}},
    "TEST-OTPSS-D":   {"name": "OTPSS-D",   "xc_functionals": {"MGGA_XC_OTPSS_D": {}}, "dispersion": {"type": "d3zero", "params": {'s6': 1.0,  's8': 1.494, 'sr6': 1.128, 'alpha6': 14.0}}},
    "TEST-TPSSLYP1W": {"name": "TPSSLYP1W", "xc_functionals": {"MGGA_XC_TPSSLYP1W": {}}},
    "TEST-B97M-V":    {"name": "B97M-V",    "xc_functionals": {"MGGA_XC_B97M_V": {}}},
}
