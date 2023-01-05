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
"""
List of GGA functionals
"""

funcs = []

funcs.append({
    "name": "BLYP",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_LYP": {}
    },
    "citation":
    '    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n' +
    '    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n',
    "description":
    '    BLYP GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "BOP",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_OP_B88": {}
    },
    "citation": '    T. Tsuneda et. al., J. Chem. Phys. 110, 10664-10678, 1999\n',
    "description": '    BOP GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "B86BPBE",
    "x_functionals": {
        "GGA_X_B86_MGC": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation": '    A. D. Becke, J. Chem. Phys. 85:7184, 1986.\n',
    "description": '    B86BPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PW86PBE",
    "x_functionals": {
        "GGA_X_PW86": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation": '    J. P. Perdew and W. Yue, Phys. Rev. B 33:8800(R), 1986.\n',
    "description": '    PW86PBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PBE",
    "x_functionals": {
        "GGA_X_PBE": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation": '    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n',
    "description": '    PBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "RPBE",
    "x_functionals": {
        "GGA_X_RPBE": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation": '    B. Hammer, L. B. Hansen, and J. K. Norskov, Phys. Rev. B, 59, 7413, 1999\n',
    "description": '    RPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "revPBE",
    "x_functionals": {
        "GGA_X_PBE_R": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation": '    Y. Zhang and W. Yang, Phys. Rev. Lett. 80, 890, 1998\n',
    "description": '    revPBE GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "N12",
    "x_functionals": {
        "GGA_X_N12": {}
    },
    "c_functionals": {
        "GGA_C_N12": {}
    },
    "citation": '    R. Peverati, D.G. Truhlar, J. Chem. Theory Comput. 8, 2310, 2012\n',
    "description": '    N12 nonseparable GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PW91",
    "x_functionals": {
        "GGA_X_PW91": {}
    },
    "c_functionals": {
        "GGA_C_PW91": {}
    },
    "citation": '    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n',
    "description": '    PW91 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "mPWPW",
    "x_functionals": {
        "GGA_X_mPW91": {}
    },
    "c_functionals": {
        "GGA_C_PW91": {}
    },
    "citation": '    C. Adamo, V. Barone, J. Chem. Phys., 108, 664, 1998\n',
    "description": '    mPWPW GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "BP86",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_P86": {}
    },
    "citation":
    '    A. D. Becke, Phys. Rev. A, 38, 3098-3100, 1988\n' +
    '    J. P. Perdew, Phys. Rev. B, 33, 8822, 1986\n',
    "description":
    '    BP86 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "FT97",
    "x_functionals": {
        "GGA_X_FT97_B": {}
    },
    "c_functionals": {
        "GGA_C_FT97": {}
    },
    "citation": '    M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997\n',
    "description": '   FT97 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "SOGGA",
    "x_functionals": {
        "GGA_X_SOGGA": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation":
    '    Y. Zhao, and D. G. Truhlar, J. Chem. Phys. 128, 184109, 2008\n' +
    '    J. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865-3868, 1996\n',
    "description":
    '   SOGGA Exchange + PBE Correlation Functional\n',
})

funcs.append({
    "name": "SOGGA11",
    "x_functionals": {
        "GGA_X_SOGGA11": {}
    },
    "c_functionals": {
        "GGA_C_SOGGA11": {}
    },
    "citation": '    R. Peverati, Y. Zhao, and D. G. Truhlar, J. Phys. Chem. Lett. 2, 1991-1997, 2011\n',
    "description": '   SOGGA11 Exchange-Correlation Functional\n',
})

# empirical_dispersion_resources.dashcoeff now defines 'b97' (<-'b97-d') as Grimme's GGA
#   functional, so one could just define
#       funcs.append({
#            "name": "B97",
#           "xc_functionals": {
#               "GGA_XC_B97_D": {}
#           }})
#   and let the dict_builder fill out -D2, -D3(BJ), etc. Leaving the
#   explicit definitions below so plain 'b97' isn't further confused.

funcs.append({
    "name": "B97-D2",
    "xc_functionals": {
        "GGA_XC_B97_D": {}
    },
    "dispersion": {
        "type": "d2",
        "params": {
            's6': 1.25,
            'alpha6': 20.0,
            'sr6': 1.1
        }
    }
})

funcs.append({
    "name": "B97-D3ZERO",
    "xc_functionals": {
        "GGA_XC_B97_D": {}
    },
    "dispersion": {
        "type": "d3zero2b",
        "params": {
            's6': 1.0,
            's8': 0.909,
            'sr6': 0.892,
            'sr8': 1.0,
            'alpha6': 14.0
        }
    }
})

funcs.append({
    "name": "B97-D3BJ",
    "xc_functionals": {
        "GGA_XC_B97_D": {}
    },
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            's6': 1.000,
            's8': 2.2609,
            'a1': 0.5545,
            'a2': 3.2297
        }
    }
})

funcs.append({
    "name": "B97-D3MZERO",
    "xc_functionals": {
        "GGA_XC_B97_D": {}
    },
    "dispersion": {
        "type": "d3mzero2b",
        "params": {
            's6': 1.000,
            's8': 1.020078,
            'sr6': 1.151808,
            'beta': 0.035964
        }
    }
})

funcs.append({
    "name": "B97-D3MBJ",
    "xc_functionals": {
        "GGA_XC_B97_D": {}
    },
    "dispersion": {
        "type": "d3mbj2b",
        "params": {
            's6': 1.000,
            's8': 1.206988,
            'a1': 0.240184,
            'a2': 3.864426
        }
    }
})

funcs.append({
    "name": "GAM",
    "x_functionals": {
        "GGA_X_GAM": {}
    },
    "c_functionals": {
        "GGA_C_GAM": {}
    },
    "citation": '    H.S. Yu, et al., Phys. Chem. Chem. Phys. 17, 12146, 2015\n',
    "description": '   GAM GGA Minessota Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "OP-PBE",
    "alias": ["OPPBE", "PBEOP"],
    "x_functionals": {
        "GGA_X_PBE": {}
    },
    "c_functionals": {
        "GGA_C_OP_PBE": {}
    },
    "citation":
    '    T. Tsuneda, T. Suzumura, K. Hirao, J. Chem. Phys. 110, 10664, 1999\n' +
    '    T. Tsuneda, T. Suzumura, K. Hirao, J. Chem. Phys. 111, 5656, 1999\n',
    "description":
    '    BP86 GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PBE-SOL",
    "alias": ["PBESOL"],
    "x_functionals": {
        "GGA_X_PBE_SOL": {}
    },
    "c_functionals": {
        "GGA_C_PBE_SOL": {}
    },
    "description": '    Perdew, Burke & Ernzerhof exchange (solids)',
})

funcs.append({
    "name": "BP86-VWN",
    "alias": ["BP86VWN"],
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_P86VWN_FT": {}
    },
    "description":
    '    BP86 GGA XC Functional based on VWN5 corr. & more accurate ftilde value\n',
})

# B97-3c = modified B97(GGA) + D3BJ + SRB(through mctc-gcp)
funcs.append({
    "name": "B973c",
    "alias": ["B97-3c"],
    "xc_functionals": {
        # "GGA_XC_B97_3C": {} #for libxc >= v6.0
        "GGA_XC_B97_D": {
            "tweak": { # needed until libxc >= v6.0
            "_cx0": 1.076616,
            "_cx1":-0.469912,
            "_cx2":3.322442,
            "_css0":0.543788,
            "_css1":-1.444420,
            "_css2":1.637436,
            "_cos0":0.635047,
            "_cos1":5.532103,
            "_cos2":-15.301575,
            },
        },
    },
    "description": '    B97-3c GGA-based 3C composite method with a TZ basis set, D3 and short-range basis set correction.\n',
    "citation": '     J. G. Brandenburg, C.Bannwarth, A. Hansen, S. Grimme J. Chem. Phys. 148, 064104, 2018\n',
    "doi": "10.1063/1.5012601",
    "dispersion": {
        "type": "d3bjatm",
        "params": {
            's6': 1.000,
            's8': 1.500,
            'a1': 0.370,
            'a2': 4.100,
        },
    },
})


functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
