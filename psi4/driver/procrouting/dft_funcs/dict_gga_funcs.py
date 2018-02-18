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


blyp = {
    "name": "BLYP",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_LYP": {}},
    "citation": '    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n',
    "description": '    BLYP GGA Exchange-Correlation Functional\n',
}

blyp_d2 = copy.deepcopy(blyp)
blyp_d2["name"] = "BLYP-D2"
blyp_d2["dispersion"] = {"type": "d2", "params": {"s6": 1.20}}

blyp_d3 = copy.deepcopy(blyp)
blyp_d3["name"] = "BLYP-D3"
blyp_d3["dispersion"] = {"type": "d3", "params": {'s6': 1.0,  's8': 1.682, 'sr6': 1.094, 'alpha6': 14.0}}

blyp_d3bj = copy.deepcopy(blyp)
blyp_d3bj["name"] = "BLYP-D3BJ"
blyp_d3bj["dispersion"] = {"type": "d3bj", "params": {'s6': 1.000, 's8':  2.6996, 'a1':  0.4298, 'a2': 4.2359}}

blyp_d3m = copy.deepcopy(blyp)
blyp_d3m["name"] = "BLYP-D3M"
blyp_d3m["dispersion"] = {"type": "d3m", "params": {'s6': 1.000, 's8':  1.841686, 'sr6': 1.279637, 'beta': 0.014370}}

blyp_d3mbj = copy.deepcopy(blyp)
blyp_d3mbj["name"] = "BLYP-D3MBJ"
blyp_d3mbj["dispersion"] = {"type": "d3mbj", "params": {'s6': 1.000, 's8': 1.875007, 'a1': 0.448486, 'a2': 3.610679}}




bop = {
    "name": "BOP",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_OP_B88": {}},
    "citation": '    T. Tsuneda et. al., J. Chem. Phys. 110, 10664-10678, 1999\n',
    "description": '    BOP GGA Exchange-Correlation Functional\n',
}

b86bpbe = {
    "name": "B86BPBE",
    "x_functionals": {"GGA_X_B86_MGC": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n',
    "description": '    B86BPBE GGA Exchange-Correlation Functional\n',
}

pw86pbe = {
    "name": "PW86PBE",
    "x_functionals": {"GGA_X_PW86": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n',
    "description": '    PW86PBE GGA Exchange-Correlation Functional\n',
}

pbe = {
    "name": "PBE",
    "x_functionals": {"GGA_X_PBE": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n',
    "description": '    PBE GGA Exchange-Correlation Functional\n',
}

pw91 = {
    "name": "PW91",
    "x_functionals": {"GGA_X_PW91": {}},
    "c_functionals": {"GGA_C_Pw91": {}},
    "citation": '    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n',
    "description": '    PW91 GGA Exchange-Correlation Functional\n',
}

mpwpw = {
    "name": "MPWPW",
    "x_functionals": {"GGA_X_mPW91": {}},
    "c_functionals": {"GGA_C_Pw91": {}},
    "citation": '    C. Adamo, V. Barone, J. Chem. Phys., 108, 664, 1998\n',
    "description": '    mPWPW GGA Exchange-Correlation Functional\n',
}

bp86 = {
    "name": "BP86",
    "x_functionals": {"GGA_X_B88": {}},
    "c_functionals": {"GGA_C_P86": {}},
    "citation": '   A. D. Becke, Phys. Rev. A, 38, 3098-3100, 1988\n   J. P. Perdew, Phys. Rev. B, 33, 8822, 1986\n',
    "description": '    BP86 GGA Exchange-Correlation Functional\n',
}

ft97 = {
    "name": "FT97",
    "x_functionals": {"GGA_X_FT97_B": {}},
    "c_functionals": {"GGA_C_FT97": {}},
    "citation": '    M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997\n',
    "description": '   FT97 GGA Exchange-Correlation Functional\n',
}

sogga = {
    "name": "SOGGA",
    "x_functionals": {"GGA_X_SOGGA": {}},
    "c_functionals": {"GGA_C_PBE": {}},
    "citation": '    Y. Zhao, and D. G. Truhlar, J. Chem. Phys. 128, 184109, 2008\n' +
                '    J. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865-3868, 1996\n',
    "description": '   SOGGA Exchange + PBE Correlation Functional\n',
}

sogga11 = {
    "name": "SOGGA11",
    "x_functionals": {"GGA_X_SOGGA11": {}},
    "c_functionals": {"GGA_C_SOGGA11": {}},
    "citation": '    R. Peverati, Y. Zhao, and D. G. Truhlar, J. Phys. Chem. Lett. 2, 1991-1997, 2011\n',
    "description": '   SOGGA11 Exchange-Correlation Functional\n',
}


functional_list = {
    "TEST-B86BPBE": b86bpbe,
    "TEST-BLYP": blyp,
    "TEST-BLYP-D2": blyp_d2,
    "TEST-BLYP-D3": blyp_d3,
    "TEST-BLYP-D3BJ": blyp_d3bj,
    "TEST-BLYP-D3M": blyp_d3m,
    "TEST-BLYP-D3MBJ": blyp_d3mbj,
    "TEST-PW86PBE": pw86pbe,
    "TEST-PBE": pbe,
    "TEST-BP86": bp86,
    "TEST-PW91": pw91,
    "TEST-FT97": ft97,
    "TEST-BOP": bop,
    "TEST-MPWPW": mpwpw,
    "TEST-SOGGA11": sogga11,
    "TEST-SOGGA": sogga,
    "TEST-OPBE-D"    : {"name": "OPBE-D"   , "xc_functionals": {"GGA_XC_OPBE_D"   : {}}, "dispersion": {"type": "d3zero", "params": {'s6': 1.0,  's8': 1.494, 'sr6': 1.128, 'alpha6': 14.0}}},
    #"TEST-OPWLYP_D"  : {"name": "OPWLYP-D" , "xc_functionals": {"GGA_XC_OPWLYP_D" : {}}}, # no dispersion
    #"TEST-OBLYP_D"   : {"name": "OBLYP-D"  , "xc_functionals": {"GGA_XC_OBLYP_D"  : {}}}, # no dispersion
    "TEST-HCTH407P" : {"name": "HCTH/407P", "xc_functionals": {"GGA_XC_HCTH_407P": {}}},
    "TEST-HCTHP76"  : {"name": "HCTH/P76" , "xc_functionals": {"GGA_XC_HCTH_P76" : {}}},
    "TEST-HCTHP14"  : {"name": "HCTH/P14" , "xc_functionals": {"GGA_XC_HCTH_P14" : {}}},
    "TEST-B97-GGA1"  : {"name": "B97-GGA1" , "xc_functionals": {"GGA_XC_B97_GGA1" : {}}},
    "TEST-KT2"       : {"name": "KT2"      , "xc_functionals": {"GGA_XC_KT2"      : {}}},
    "TEST-TH1"       : {"name": "TH1"      , "xc_functionals": {"GGA_XC_TH1"      : {}}},
    "TEST-TH2"       : {"name": "TH2"      , "xc_functionals": {"GGA_XC_TH2"      : {}}},
    "TEST-TH3"       : {"name": "TH3"      , "xc_functionals": {"GGA_XC_TH3"      : {}}},
    "TEST-TH4"       : {"name": "TH4"      , "xc_functionals": {"GGA_XC_TH4"      : {}}},
    "TEST-HCTH93"   : {"name": "HCTH/93"  , "xc_functionals": {"GGA_XC_HCTH_93"  : {}}},
    "TEST-HCTH120"  : {"name": "HCTH/120" , "xc_functionals": {"GGA_XC_HCTH_120" : {}}},
    "TEST-HCTH147"  : {"name": "HCTH/147" , "xc_functionals": {"GGA_XC_HCTH_147" : {}}},
    "TEST-HCTH407"  : {"name": "HCTH/407" , "xc_functionals": {"GGA_XC_HCTH_407" : {}}},
    "TEST-EDF1"      : {"name": "EDF1"     , "xc_functionals": {"GGA_XC_EDF1"     : {}}},
    "TEST-XLYP"      : {"name": "XLYP"     , "xc_functionals": {"GGA_XC_XLYP"     : {}}},
    "TEST-B97-D"     : {"name": "B97-D"    , "xc_functionals": {"GGA_XC_B97_D"    : {}}, "dispersion": {"type": "d2", "params": {'s6': 1.25}}},
    "TEST-PBE1W"     : {"name": "PBE1W"    , "xc_functionals": {"GGA_XC_PBE1W"    : {}}},
    "TEST-MPWLYP1W"  : {"name": "MPWLYP1W" , "xc_functionals": {"GGA_XC_MPWLYP1W" : {}}},
    "TEST-PBELYP1W"  : {"name": "PBELYP1W" , "xc_functionals": {"GGA_XC_PBELYP1W" : {}}},
    "TEST-MOHLYP"    : {"name": "MOHLYP"   , "xc_functionals": {"GGA_XC_MOHLYP"   : {}}},
    "TEST-MOHLYP2"   : {"name": "MOHLYP2"  , "xc_functionals": {"GGA_XC_MOHLYP2"  : {}}},
    "TEST-TH-FL"     : {"name": "TH-FL"    , "xc_functionals": {"GGA_XC_TH_FL"    : {}}},
    "TEST-TH-FC"     : {"name": "TH-FC"    , "xc_functionals": {"GGA_XC_TH_FC"    : {}}},
    "TEST-TH-FCFO"   : {"name": "TH-FCFO"  , "xc_functionals": {"GGA_XC_TH_FCFO"  : {}}},
    "TEST-TH-FCO"    : {"name": "TH-FCO"   , "xc_functionals": {"GGA_XC_TH_FCO"   : {}}},
    "TEST-VV10"      : {"name": "VV10"     , "xc_functionals": {"GGA_XC_VV10"     : {}}},
}
