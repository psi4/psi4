#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
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
List of double-hybrid functionals
"""

funcs = []

funcs.append({
    "name": "MP2MP2",
    "description": "2nd-order MP perturbation theory",
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
    "c_mp2": {
        "alpha": 1.0
    },
})

funcs.append({
    "name": "MP2D",
    "description": "2nd-order MP perturbation theory plus dispersion",
    "alias": ["MP2-D"],
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
    "c_mp2": {
        "alpha": 1.0
    },
    "dispersion": {
        "type": "dmp2",
        "params": {
            "s8": 1.187,
            "a1": 0.944,
            "a2": 0.480,
            "rcut": 0.72,
            "w": 0.20,
        },
        "citation": "    Rezac, J.; Greenwell, C.; Beran, G. (2018), J. Chem. Theory Comput., 14: 4711-4721\n",
    },
})

funcs.append({
    "name": "B2PLYP",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.47
        }
    },
    "x_hf": {
        "alpha": 0.53
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.73
        }
    },
    "c_mp2": {
        "alpha": 0.27
    },
    "citation": '    S. Grimme, J. Chem. Phys., 124, 034108, 2006\n',
    "description": '    B2PLYP Double Hybrid Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "DSD-BLYP",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.25
        }
    },
    "x_hf": {
        "alpha": 0.75
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.53
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.60
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-BLYP SCS Double Hybrid XC Functional (not dispersion corrected)\n',
})

funcs.append({
    "name": "DSD-BLYP-D2",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.29
        }
    },
    "x_hf": {
        "alpha": 0.71
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.55
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.43
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-BLYP-D2 Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.35,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-BLYP-D3BJ",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.29
        }
    },
    "x_hf": {
        "alpha": 0.71
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.54
        }
    },
    "c_mp2": {
        "os": 0.47,
        "ss": 0.40,
    },
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            "s6": 0.57,
            "a2": 5.4,
            "a1": 0.0,
            "s8": 0.0
        },
        "citation": "    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n"
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-BLYP-D3BJ Dispersion-corrected SCS Double Hybrid XC Functional\n',
})

funcs.append({
    # note by H.Kruse: Uses the full-core parameters in the Yu paper. But my and L. Georigk's experience shows that it hardly matters.
    # Use FC is recommended by S. Grimme.
    # May this madness never end.
    "name": "DSD-BLYP-NL",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.29
        }
    },
    "x_hf": {
        "alpha": 0.71
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.54
        }
    },
    "c_mp2": {
        "alpha": 1.0,
        "os": 0.47,
        "ss": 0.40,
    },
    "dispersion": {
        "type": "nl",
        "params": {
            "b": 12.00,
            "c": 0.0093,
        },
        "citation": "    F. Yu J. Chem. Theory Comput. 10, 4400-4407, 2014\n"
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-BLYP-NL (D3BJ,FC parameters) VV10 SCS Double Hybrid XC Functional\n',
})

funcs.append({
    "name": "CORE-DSD-BLYP",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.31
        }
    },
    "x_hf": {
        "alpha": 0.69
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.54
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.37
    },
    "citation": '    S. Kozuch, Phys. Chem. Chem. Phys., 13, 20104, 2011\n',
    "description": '    DSD-BLYP SCS Double Hybrid XC Functional (full-core param.)\n'
})

funcs.append({
    "name": "PBE0-2",
    "alias": ["PBE02"],
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.206299
        }
    },
    "x_hf": {
        "alpha": 0.793701
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.5
        }
    },
    "c_mp2": {
        "alpha": 0.5
    },
    "citation": '    J. Chai, Chem. Phys. Lett., 538, 121-125, 2012\n',
    "description": '    PBE0-2 Double Hybrid Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "PBE0-DH",
    "alias": ["PBE0DH"],
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.5
        }
    },
    "x_hf": {
        "alpha": 0.5
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.875
        }
    },
    "c_mp2": {
        "alpha": 0.125
    },
    "citation": '    E. Bremond, C. Adamo, J. Chem. Phys., 135, 024106, 2011\n',
    "description": '    PBE0-DH Double Hybrid Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "DSD-PBEP86",
    "alias": ["DSDPBEP86"],
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.28
        }
    },
    "x_hf": {
        "alpha": 0.72
    },
    "c_functionals": {
        "GGA_C_P86": {
            "alpha": 0.44
        }
    },
    "c_mp2": {
        "os": 0.51,
        "ss": 0.36
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEP86 SCS Double Hybrid XC Functional (not dispersion corrected)\n',
})

funcs.append({
    "name": "DSD-PBEP86-D3BJ",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.31
        }
    },
    "x_hf": {
        "alpha": 0.69
    },
    "c_functionals": {
        "GGA_C_P86": {
            "alpha": 0.44
        }
    },
    "c_mp2": {
        "os": 0.52,
        "ss": 0.22
    },
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            "s6": 0.48,
            "a2": 5.6,
            "a1": 0.0,
            "s8": 0.0
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEP86-D3BJ Dispersion-corrected SCS Double Hybrid XC Functional\n',
})

funcs.append({
    # note: Using the D3BJ form for NL, which is sensible but not explicitly mentioned in the paper
    "name": "DSD-PBEP86-NL",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.31
        }
    },
    "x_hf": {
        "alpha": 0.69
    },
    "c_functionals": {
        "GGA_C_P86": {
            "alpha": 0.44
        }
    },
    "c_mp2": {
        "os": 0.52,
        "ss": 0.22
    },
    "dispersion": {
        "type": "nl",
        "params": {
            "b": 12.8,
            "c": 0.0093,
        },
        "citation": '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016 \n'
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEP86-NL (D3BJ parameters) VV10 SCS Double Hybrid XC Functional\n',
})

funcs.append({
    "name": "DSD-PBEP86-D2",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.32
        }
    },
    "x_hf": {
        "alpha": 0.68
    },
    "c_functionals": {
        "GGA_C_P86": {
            "alpha": 0.45
        }
    },
    "c_mp2": {
        "os": 0.51,
        "ss": 0.23
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEP86-D2 Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.29,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-PBEPBE",
    "alias": ["DSDPBEPBE"],
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.28
        }
    },
    "x_hf": {
        "alpha": 0.72
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.48
        }
    },
    "c_mp2": {
        "os": 0.54,
        "ss": 0.31
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEPBE SCS Double Hybrid XC Functional (not dispersion corrected)\n',
})

funcs.append({
    "name": "DSD-PBEPBE-D2",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.34
        }
    },
    "x_hf": {
        "alpha": 0.66
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.51
        }
    },
    "c_mp2": {
        "os": 0.53,
        "ss": 0.12
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEPBE-D2 Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.42,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-PBEPBE-D3BJ",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.32
        }
    },
    "x_hf": {
        "alpha": 0.68
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.49
        }
    },
    "c_mp2": {
        "os": 0.55,
        "ss": 0.13
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEPBE-D3BJ Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            "s6": 0.78,
            "a2": 6.1,
            "a1": 0.0,
            "s8": 0.0
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-PBEPBE-NL",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.32
        }
    },
    "x_hf": {
        "alpha": 0.68
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "alpha": 0.49
        }
    },
    "c_mp2": {
        "alpha": 1.0,
        "os": 0.55,
        "ss": 0.13
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEPBE-NL (D3BJ parameters) VV10 SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "nl",
        "params": {
            "b": 9.6,
            "c": 0.0093,
        },
        "citation": '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n'
    },
})

funcs.append({
    "name": "DSD-BP86-D2",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.33
        }
    },
    "x_hf": {
        "alpha": 0.67
    },
    "c_functionals": {
        "GGA_C_P86": {
            "alpha": 0.49
        }
    },
    "c_mp2": {
        "os": 0.49,
        "ss": 0.24
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-BP86-D2 Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.41,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-SVWN-D2",
    "x_functionals": {
        "LDA_X": {
            "alpha": 0.29
        }
    },
    "x_hf": {
        "alpha": 0.71
    },
    "c_functionals": {
        "LDA_C_VWN": {
            "alpha": 0.34
        }
    },
    "c_mp2": {
        "os": 0.58,
        "ss": 0.11
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-SVWN5-D2 Dispersion-corrected SCS Double Hybrid XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.28,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    },
})

funcs.append({
    "name": "B2GPPLYP",
    "x_functionals": {
        "GGA_X_B88": {
            "alpha": 0.35
        }
    },
    "x_hf": {
        "alpha": 0.65
    },
    "c_functionals": {
        "GGA_C_LYP": {
            "alpha": 0.64
        }
    },
    "c_mp2": {
        "alpha": 0.36
    },
    "citation": '    A. Karton, et al., J.Phys. Chem. A, 112, 12868-12886, 2008\n',
    "description": '    B2GPPLYP Double Hybrid Exchange-Correlation Functional\n',
})


def get_pwpb95_tweaks():
    X2S = 0.1282782438530421943003109254455883701296
    bt = 0.004440  # paper values
    c_pw = 0.32620  # paper values
    expo_pw6 = 3.7868  # paper values
    alpha_pw6 = c_pw / X2S / X2S
    return {"_bt": bt, "_alpha": alpha_pw6, "_expo": expo_pw6}


funcs.append({
    "name": "PWPB95",
    "x_functionals": {
        "GGA_X_MPW91": {  # only mpw91, not pw91, is tweakable
            "tweak": get_pwpb95_tweaks(),
            "alpha": 0.50
        }
    },
    "x_hf": {
        "alpha": 0.50
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "tweak": {
                "_css": 0.03241,
                "_copp": 0.00250,
            },
            "alpha": 0.731
        }
    },
    "c_mp2": {
        "ss": 0.0,
        "os": 0.269
    },
    "citation": '    L. Goerigk, S.Grimme, J.Chem. Theory Compt. 7, 291-309, 2011 \n',
    "description": '    PWPB95 SOS Double Hybrid XC Functional\n',
})

funcs.append({
    "name": "PTPSS",
    "x_functionals": {
        "MGGA_X_TPSS": {
            "tweak": {
                "_b": 0.15,
                "_c": 0.88491,
                "_e": 0.047,
                "_kappa": 0.872,
                "_mu": 0.16952,
            },
            "alpha": 0.50
        }
    },
    "x_hf": {
        "alpha": 0.50
    },
    "c_functionals": {
        "MGGA_C_TPSS": {
            "tweak": {
                "_beta": 0.06080,
                "_d": 6.3,
            },
            "alpha": 0.625
        }
    },
    "c_mp2": {
        "ss": 0.0,
        "os": 0.375
    },
    "citation": '    L. Goerigk, S.Grimme, J. Chem. Theory Comput., 7, 291-309, 2011 \n',
    "description": '    PTPSS SOS Double Hybrid XC Functional\n',
})

funcs.append({
    "name": "DSD-PBEB95",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.31
        }
    },
    "x_hf": {
        "alpha": 0.69
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "alpha": 0.54
        }
    },
    "c_mp2": {
        "os": 0.48,
        "ss": 0.22
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEB95 SCS Double Hybrid Meta-GGA XC Functional (not dispersion corrected)\n',
})

funcs.append({
    "name": "DSD-PBEB95-D2",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.35
        }
    },
    "x_hf": {
        "alpha": 0.65
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "alpha": 0.55
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.08
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEB95-D2 Dispersion-corrected SCS Double Hybrid Meta-GGA XC Functional\n',
    "dispersion": {
        "type": "d2",
        "params": {
            "s6": 0.32,
            "alpha6": 20.0,
            "sr6": 1.1
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    },
})

funcs.append({
    "name": "DSD-PBEB95-D3BJ",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.34
        }
    },
    "x_hf": {
        "alpha": 0.66
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "alpha": 0.55
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.09
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEB95-D3BJ Dispersion-corrected SCS Double Hybrid Meta-GGA XC Functional\n',
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            "s6": 0.61,
            "a2": 6.2,
            "a1": 0.0,
            "s8": 0.0
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n'
    },
})

funcs.append({
    "name": "DSD-PBEB95-NL",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.34
        }
    },
    "x_hf": {
        "alpha": 0.66
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "alpha": 0.55
        }
    },
    "c_mp2": {
        "os": 0.46,
        "ss": 0.09
    },
    "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
    "description": '    DSD-PBEB95-NL (D3BJ parameters) VV10 SCS Double Hybrid Meta-GGA XC Functional\n',
    "dispersion": {
        "type": "nl",
        "params": {
            "b": 12.50,
            "c": 0.0093,
        },
        "citation": '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n'
    },
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
