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

import copy

funcs = []

funcs.append({
    "name": "B5050LYP",
    "description": '    B5050LYP Hyb-GGA Exchange-Correlation Functional\n',
    "citation": '    Y. Shao et. al., J. Chem. Phys., 188, 4807-4818, 2003\n',
    "x_functionals": {
        "LDA_X": {
            "alpha": 0.08
        },
        "GGA_X_B88": {
            "alpha": 0.42
        }
    },
    "x_hf": {
        "alpha": 0.50
    },
    "c_functionals": {
        "LDA_C_VWN": {
            "alpha": 0.19
        },
        "GGA_C_LYP": {
            "alpha": 0.81
        }
    },
})

funcs.append({
    "name":
    "wPBE0",
    "alias": ["LC-WPBE0"],
    "description":
    '    PBE0 SR-XC Functional (HJS Model)\n',
    "citation":
    '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
    "x_functionals": {
        "GGA_X_HJS_PBE": {
            "omega": 0.3,
            "alpha": 0.75
        }
    },
    "x_hf": {
        "alpha": 0.25,
        "beta": 1.0,
        "omega": 0.3
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
})

funcs.append({
    "name":
    "wPBE",
    "alias": ["LC-WPBE"],
    "description":
    '    PBE SR-XC Functional (HJS Model)\n',
    "citation":
    '    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n    Weintraub, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754, 2009\n',
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
})

funcs.append({
    "name": "wB97X-D3",
    "description": '    Parameterized Hybrid LRC B97 GGA XC Functional with Dispersion\n',
    "citation": '    J.-D. Chai and M. Head-Gordon, Phys. Chem. Chem. Phys., 10, 6615-6620, 2008\n',
    "x_functionals": {
        "HYB_GGA_XC_WB97X": {
            "omega": 0.25
        }
    },
    "x_hf": {
        "alpha": 0.195728,
        "beta": 1.0,
        "omega": 0.25
    },
    "c_functionals": {},
    "dispersion": {
        "type": "d3zero",
        "params": {
            's6': 1.0,
            's8': 1.000,
            'sr6': 1.281,
            'sr8': 1.094,
            'alpha6': 14.0
        }
    }
})

funcs.append({
    "name": "HF",
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
    "name": "HF-3C",
    "alias": ["HF3C"],
    "description": '    Hartree Fock as Roothaan prescribed plus 3C\n',
    "citation": '    Sure et al., J. Comput. Chem., 34, 1672-1685, 2013\n',
    "x_hf": {
        "alpha": 1.0
    },
    "c_functionals": {},
    "dispersion": {
        "type": "d3bj",
        "params": {
            's6': 1.000,
            's8': 0.8777,
            'a1': 0.4171,
            'a2': 2.9149
        }
    },
})

funcs.append({
    "name": "PBEH-3C",
    "alias": ["PBEH3C"],
    "description": '    PBEH-3C Hybrid GGA Exchange-Correlation Functional plus 3C\n',
    "citation": '    Grimme et. al., J. Chem. Phys., 143, 054107, 2015\n',
    "x_functionals": {
        "GGA_X_PBE": {
            "tweak": [1.0245, 0.12345679],
            "alpha": 0.58
        }
    },
    "x_hf": {
        "alpha": 0.42
    },
    "c_functionals": {
        "GGA_C_PBE": {
            "tweak": [0.03]
        }
    },
    "dispersion": {
        "type": "d3bj",
        "params": {
            's6': 1.000,
            's8': 0.0000,
            'a1': 0.4860,
            'a2': 4.5000
        }
    },
})

funcs.append({
    "name": "SOGGA11-X",
    "alias": ["SOGGA11X"],
    "description": '   SOGGA11-X Hybrid Exchange-Correlation Functional\n',
    "citation": '    R. Peverati and D. G. Truhlar, J. Chem. Phys. 135, 191102, 2011\n',
    "x_functionals": {
        "HYB_GGA_X_SOGGA11_X": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "GGA_C_SOGGA11_X": {}
    },
})

funcs.append({
    "name": "MN12-SX",
    "alias": ["MN12SX"],
    "description": '   MN12-SX Meta-GGA Hybrid Screened Exchange-Correlation Functional\n',
    "citation": '    R. Peverati, D. G. Truhlar, Phys. Chem. Chem. Phys 14, 16187, 2012\n',
    "x_functionals": {
        "HYB_MGGA_X_MN12_SX": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_MN12_SX": {}
    },
})

funcs.append({
    "name": "MN15",
    "description": '   MN15 Hybrid Meta-GGA Exchange-Correlation Functional\n',
    "citation": '    H. S. Yu, X. He, S. L. Li, and D. G. Truhlar, Chem. Sci. 7, 5032-5051, 2016\n',
    "x_functionals": {
        "HYB_MGGA_X_MN15": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_MN15": {}
    },
})


def get_pw6b95_tweaks():
    beta = 0.0018903811666999256  # 5.0*(36.0*math.pi)**(-5.0/3.0)
    X2S = 0.1282782438530421943003109254455883701296
    X_FACTOR_C = 0.9305257363491000250020102180716672510262  #    /* 3/8*cur(3/pi)*4^(2/3) */
    bt = 0.00538  # paper values
    c_pw = 1.7382  # paper values
    expo_pw6 = 3.8901  # paperl values
    alpha_pw6 = c_pw / X2S / X2S
    a_pw6 = 6.0 * bt / X2S
    b_pw6 = 1.0 / X2S
    c_pw6 = bt / (X_FACTOR_C * X2S * X2S)
    d_pw6 = -(bt - beta) / (X_FACTOR_C * X2S * X2S)
    f_pw6 = 1.0e-6 / (X_FACTOR_C * X2S**expo_pw6)
    return ([a_pw6, b_pw6, c_pw6, d_pw6, f_pw6, alpha_pw6, expo_pw6])


funcs.append({
    "name": "PW6B95",
    "description": '    PW6B95 Hybrid Meta-GGA XC Functional\n',
    "citation": '  Y. Zhao and D. Truhlar, J. Phys. Chem. A., 109,5656-5667, 2005\n',
    "x_functionals": {
        "GGA_X_PW91": {
            "tweak": get_pw6b95_tweaks(),
            "alpha": 0.72
        }
    },
    "x_hf": {
        "alpha": 0.28
    },
    "c_functionals": {
        "MGGA_C_BC95": {
            "tweak": [0.03668, 0.00262]
        }
    },
})

funcs.append({
    "name": "dlDF",
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
    "x_functionals": {
        "HYB_MGGA_X_DLDF": {
            "use_libxc": True
        }
    },
    "c_functionals": {
        "MGGA_C_DLDF": {}
    },
})

funcs.append({
    "name": "dlDF+D09",
    "alias": ["DLDF-D09", "DLDF-DAS2009"],
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
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
    }
})

funcs.append({
    "name": "dlDF+D10",
    "alias": ["DLDF-D10", "DLDF-DAS2010", "DLDF+D"],
    "description": '    Dispersionless Hybrid Meta-GGA XC Functional\n',
    "citation": '    Pernal et. al., Phys. Rev. Lett., 103, 263201, 2009\n',
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
    }
})

functional_list = {}
for functional in funcs:
    if "alias" in functional.keys():
        alias = functional.pop("alias")
        for a in alias:
            functional_list["TEST-" + a] = functional
    functional_list["TEST-" + functional["name"].upper()] = functional
