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
List of hybrid functionals
"""

funcs = []

funcs.append({
    # This should be replaced with a call to libxc's
    # HYB_GGA_XC_LC_WPBEH_WHS once it's possible
    "name": "wPBE0",
    "alias": ["LC-WPBE0"],
    "x_functionals": {
        "GGA_X_HJS_PBE": {
            "omega": 0.3,
            "alpha": 0.75
        }
    },
    "x_hf": {
        "alpha": 0.25,
        "beta": 0.75,
        "omega": 0.3
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "description":
    '    PBE0 SR-XC Functional (HJS Model)\n',
    "citation":
    '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n' + \
    '    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
})

funcs.append({
    # This should be replaced with a call to libxc's
    # HYB_GGA_XC_LC_WPBE_WHS once it's possible
    "name": "wPBE",
    "alias": ["LC-WPBE", "LCWPBE"],
    "x_functionals": {
        "GGA_X_HJS_PBE": {
            "omega": 0.4
        }
    },
    "x_hf": {
        "alpha": 0.0,
        "beta": 1.0,
        "omega": 0.4
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "description":
    '    PBE SR-XC Functional (HJS Model)\n',
    "citation":
    '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n' + \
    '    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
})

funcs.append({
    "name": "HF",
    "description": "Hartree--Fock theory",
    "alias": ["SCF"],
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
})

funcs.append({
    "name": "HF+D",
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
    "dispersion": {
        "type": "das2010",
        "params": {
            "s6": 1.0
        }
    }
})

funcs.append({
    "name": "HF3c",
    "alias": ["HF-3C"],
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            's6': 1.000,
            's8': 0.8777,
            'a1': 0.4171,
            'a2': 2.9149
        }
    },
    "description": '    Hartree Fock based 3C composite method with minimal basis set, gCP and D3(BJ)\n',
    "citation": '    Sure et al., J. Comput. Chem., 34, 1672-1685, 2013\n',
})

funcs.append({
    "name": "PBEh3c",
    "alias": ["PBEH-3C"],
    "x_functionals": {
        "GGA_X_PBE": {
            "tweak": {
                "_kappa": 1.0245,
                "_mu": 0.12345679,
            },
            "alpha": 0.58
        }
    },
    "x_hf": {
        "alpha": 0.42
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "tweak": {
                "_beta": 0.03,
            },
        }
    },
    "dispersion": {
        "type": "d3bj2b",
        "params": {
            's6': 1.000,
            's8': 0.0000,
            'a1': 0.4860,
            'a2': 4.5000
        }
    },
    "description": '    PBE Hybrid based 3C composite method with a small basis set, gCP and D3(BJ)\n',
    "citation": '    Grimme et. al., J. Chem. Phys., 143, 054107, 2015\n',
})

funcs.append({
    "name": "wB97X3c",
    "alias": ["wB97X-3C"],
    "xc_functionals": {
        "HYB_GGA_XC_WB97X_V": {}
    },
    "dispersion": {
        "type": "d4bjeeqatm",
        "nlc": False,
        "params": {
            'a1': 0.2464,
            'a2': 4.737,
            's6': 1.000,
            's8': 0.000,
            's9': 1.000,
        },
    },
    "description":
    '    wB97X basied 3C composite method with a small basis set, gCP and D4\n',
    "citation":
    '    M. Muller, A. Hansen, S. Grimme, J. Chem. Phys. 158, 014103 (2023)\n',
    "doi": "10.1063/5.0133026",
})

funcs.append({
    "name": "SOGGA11-X",
    "alias": ["SOGGA11X"],
    "x_functionals": {
        "HYB_GGA_X_SOGGA11_X": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_SOGGA11_X": {}
    },
    "description": '   SOGGA11-X Hybrid Exchange-Correlation Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Chem. Phys. 135, 191102, 2011\n',
})

funcs.append({
    "name": "MN12-SX",
    "alias": ["MN12SX"],
    "x_functionals": {
        "HYB_MGGA_X_MN12_SX": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_MN12_SX": {}
    },
    "description": '   MN12-SX Meta-GGA Hybrid Screened Exchange-Correlation Functional\n',
    "citation": '    R. Peverati, D. G. Truhlar, Phys. Chem. Chem. Phys 14, 16187, 2012\n',
})

funcs.append({
    "name": "MN15",
    "x_functionals": {
        "HYB_MGGA_X_MN15": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_MN15": {}
    },
    "description": '   MN15 Hybrid Meta-GGA Exchange-Correlation Functional\n',
    "citation": '    H. S. Yu, X. He, S. L. Li, and D. G. Truhlar, Chem. Sci. 7, 5032-5051, 2016\n',
})

funcs.append({
    "name": "N12-SX",
    "alias": ["N12SX"],
    "x_functionals": {
        "HYB_GGA_X_N12_SX": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_N12_SX": {}
    },
    "description": '   N12-SX Hybrid nonseparable GGA Exchange-Correlation Functional\n',
    "citation": '    R. Peverati, D.G. Truhlar, Phys. Chem. Chem. Phys. 14, 16187, 2012\n',
})

funcs.append({
    "name": "PBE50",
    "x_functionals": {
        "GGA_X_PBE": {
            "alpha": 0.5
        }
    },
    "x_hf": {
        "alpha": 0.5
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "description": '   PBE50 Hybrid GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "revPBE0",
    "x_functionals": {
        "GGA_X_PBE_R": {
            "alpha": 0.75
        }
    },
    "x_hf": {
        "alpha": 0.25
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "description": '   revPBE0 Hybrid GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "mPW1LYP",
    "x_functionals": {
        "GGA_X_MPW91": {
            "alpha": 0.75
        }
    },
    "x_hf": {
        "alpha": 0.25
    },
    "c_functionals": {
        "GGA_C_LYP": {}
    },
    "description": '   mPW1LYP Hybrid GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "mPW1PBE",
    "x_functionals": {
        "GGA_X_MPW91": {
            "alpha": 0.75
        }
    },
    "x_hf": {
        "alpha": 0.25
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "description": '   mPW1PBE Hybrid GGA Exchange-Correlation Functional\n',
})

funcs.append({
    "name": "MGGA_MS2h",
    "alias": ["MGGA-MS2H", "MS2H"],
    "x_functionals": {
        "HYB_MGGA_X_MS2H": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '   MGGA_MS2h Hybrid Meta-GGA XC Functional\n',
    "citation": '    J. Sun, et al., J. Chem. Phys. 138, 044113, 2013\n',
})

funcs.append({
    "name": "MGGA_MVSh",
    "alias": ["MGGA-MVSH", "MVSH"],
    "x_functionals": {
        "HYB_MGGA_X_MVSH": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_REGTPSS": {}
    },
    "description": '   MGGA_MV2h Hybrid Meta-GGA XC Functional\n',
    "citation": '    J. Sun, J.P. Perdew, A. Ruzsinsky, Proc. Natl. Acad. Sci. USA 112, 685, 2015\n',
})


def get_pw6b95_tweaks():
    X2S = 0.1282782438530421943003109254455883701296
    bt = 0.00538  # paper values
    c_pw = 1.7382  # paper values
    expo_pw6 = 3.8901  # paper values
    alpha_pw6 = c_pw / X2S / X2S
    return {"_bt": bt, "_alpha": alpha_pw6, "_expo": expo_pw6}


funcs.append({
    "name": "PW6B95",
    "x_functionals": {
        "GGA_X_MPW91": {
            "tweak": get_pw6b95_tweaks(),
            "alpha": 0.72
        }
    },
    "x_hf": {
        "alpha": 0.28
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "tweak": {
                "_css": 0.03668,
                "_copp": 0.00262,
            },
        }
    },
    "description": '    PW6B95 Hybrid Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. Truhlar, J. Phys. Chem. A., 109,5656-5667, 2005\n',
})

funcs.append({
    "name": "dlDF",
    "x_functionals": {
        "HYB_MGGA_X_DLDF": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_DLDF": {}
    },
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
})

funcs.append({
    "name": "dlDF+D09",
    "alias": ["DLDF-D09", "DLDF-DAS2009"],
    "x_functionals": {
        "HYB_MGGA_X_DLDF": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_DLDF": {}
    },
    "dispersion": {
        "type": "das2009",
        "params": {
            "s6": 1.0
        }
    },
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
})

funcs.append({
    "name": "dlDF+D10",
    "alias": ["DLDF-D10", "DLDF-DAS2010", "DLDF+D"],
    "x_functionals": {
        "HYB_MGGA_X_DLDF": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_DLDF": {}
    },
    "dispersion": {
        "type": "das2010",
        "params": {
            "s6": 1.0
        }
    },
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
})

funcs.append({
    "name": "M05",
    "x_functionals": {
        "HYB_MGGA_X_M05": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M05": {}
    },
    "description": '    M05 Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and N. E. Schultz and D. G. Truhlar, J. Chem. Phys. 123, 161103, 2005\n',
    "doi": '10.1063/1.2126975',
})

funcs.append({
    "name": "M05-2X",
    "alias": ["M052X"],
    "x_functionals": {
        "HYB_MGGA_X_M05_2X": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M05_2X": {}
    },
    "description": '    Heavily Parameterized Hybrid M05-2X Meta-GGA XC Functional\n',
    "citation": '    Zhao et. al., J. Chem. Theory Comput., 2, 364, 2006\n',
})

funcs.append({
    "name": "M06",
    "x_functionals": {
        "HYB_MGGA_X_M06": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M06": {}
    },
    "description": '    M06 Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar,  Theor. Chem. Acc., 120, 215, 2008\n',
    "doi": '10.1007/s00214-007-0310-x',
})

funcs.append({
    "name": "M06-2X",
    "alias": ["M062X"],
    "x_functionals": {
        "HYB_MGGA_X_M06_2X": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M06_2X": {}
    },
    "description": '    Hybrid M06-2X Meta-GGA XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar,  Theor. Chem. Acc., 120, 215, 2008\n',
    "doi": '10.1007/s00214-007-0310-x',
})

funcs.append({
    "name": "M06-HF",
    "alias": ["M06HF"],
    "x_functionals": {
        "HYB_MGGA_X_M06_HF": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M06_HF": {}
    },
    "description": '    Minnesota M06-HF Hybrid XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Phys. Chem. A, 110, 13126, 2006\n',
    "doi": '10.1021/jp066479k',
})

funcs.append({
    "name": "M08-HX",
    "alias": ["M08HX"],
    "x_functionals": {
        "HYB_MGGA_X_M08_HX": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M08_HX": {}
    },
    "description": '    Minnesota M08-HX Hybrid XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Chem. Theory Comput., 4, 1849, 2008\n',
    "doi": '10.1021/ct800246v',
})

funcs.append({
    "name": "M08-SO",
    "alias": ["M08SO"],
    "x_functionals": {
        "HYB_MGGA_X_M08_SO": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M08_SO": {}
    },
    "description": '    Minnesota M08-SO Hybrid XC Functional\n',
    "citation": '    Y. Zhao and D. G. Truhlar, J. Chem. Theory Comput., 4, 1849, 2008\n',
    "doi": '10.1021/ct800246v',
})

funcs.append({
    "name": "M11",
    "x_functionals": {
        "HYB_MGGA_X_M11": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_M11": {}
    },
    "description": '    M11 Meta-GGA XC Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Phys. Chem. Lett., 2, 2810, 2011\n',
    "doi": '10.1021/jz201170d',
})

funcs.append({
    "name": "wB97M-D3BJ",
    "xc_functionals": {
        "HYB_MGGA_XC_WB97M_V": {}
    },
    "dispersion": {
        "nlc": False,
        "type": "d3bj2b",
        "citation": '    A. Najib, L. Goerigk, J. Comput. Theory Chem.,14, 5725, 2018\n',
        "params": {
            's6': 1.000,
            's8': 0.3908,
            'a1': 0.5660,
            'a2': 3.1280
        },
    },
    "description":
    '    wB97M-V with D3(BJ) instead of VV10 dispersion \n',
    "citation":
    '    A. Najib, L. Goerigk, J. Comput. Theory Chem.,14, 5725, 2018\n' +
    '    N. Mardirossian, M. Head-Gordon, J. Chem. Phys. 144, 214110, 2016\n'
})

funcs.append({
    "name": "wB97X-D3BJ",
    "xc_functionals": {
        "HYB_GGA_XC_WB97X_V": {}
    },
    "dispersion": {
        "nlc": False,
        "type": "d3bj2b",
        "citation": '    A. Najib, L. Goerigk, J. Comput. Theory Chem.,14 5725, 2018\n',
        "params": {
            's6': 1.000,
            's8': 0.2641,
            'a1': 0.0000,
            'a2': 5.4959,
        },
    },
    "description":
    '    wB97X-V with D3(BJ) instead of VV10 dispersion \n',
    "citation":
    '    A. Najib, L. Goerigk, J. Comput. Theory Chem., 14 5725, 2018)\n' +
    '    N. Mardirossian, M. Head-Gordon, Phys. Chem. Chem. Phys, 16, 9904, 2014\n'
})

funcs.append({
    "name": "SCAN0",
    "x_functionals": {
        "HYB_MGGA_X_SCAN0": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_SCAN": {}
    },
    "description": '    SCAN0 Hybrid Meta-GGA XC Functional\n',
    "citation": '    K. Hui, J.-D. Chai, J. Chem. Phys. 144, 044114, 2016\n',
    "doi": '10.1063/1.4940734',
})

funcs.append({
    "name": "revSCAN0",
    "x_functionals": {
        "HYB_MGGA_X_REVSCAN0": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_REVSCAN": {}
    },
    "description": '    Revised SCAN0 Hybrid Meta-GGA XC Functional\n',
    "citation": '    P. D. Mezei, G. I. Csonka, M. Kallay J. Chem. Theory Comput. 14, 2469, 2018\n',
    "doi": "10.1021/acs.jctc.8b00072",
})

funcs.append({
    "name": "BMK",
    "x_functionals": {
        "HYB_MGGA_X_BMK": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_BMK": {}
    },
    "citation": '    A. D. Boese, J. M. L. Martin, J. Chem. Phys. 121, 3405, 2004\n',
    "description": '    BMK Hybrid Meta-GGA XC Functional for kinetics\n',
})

funcs.append({
    "name": "t-HCTHh",
    "alias": ["tHCTHh", "tau-HCTHh"],
    "x_functionals": {
        "HYB_MGGA_X_TAU_HCTH": {
            "use_libxc": True,
        },
    },
    "c_functionals": {
        "GGA_C_HYB_TAU_HCTH": {}
    },
    "description": '    Hybrid Tau HCTH Meta-GGA XC Functional\n',
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
