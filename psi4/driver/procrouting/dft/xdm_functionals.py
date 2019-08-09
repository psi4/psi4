#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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
List of XDM-corrected functionals
"""

funcs = []

funcs.append({
    "name": "BLYP-XDM",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_LYP": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    BLYP Exchange-Correlation Functional plus XDM dispersion (GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.7647,
            'a2': 0.8457,
            'vol': 'blyp',
        }
    }
})
funcs.append({
    "name": "PBE-XDM",
    "x_functionals": {
        "GGA_X_PBE": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    PBE Exchange-Correlation Functional plus XDM dispersion (GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.4492,
            'a2': 2.5517,
            'vol': 'pbe',
        }
    }
})
funcs.append({
    "name": "PW86PBE-XDM",
    "x_functionals": {
        "GGA_X_PW86": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    PW86PBE Exchange-Correlation Functional plus XDM dispersion (GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.7564,
            'a2': 1.4545,
            'vol': 'pw86pbe',
        }
    }
})
funcs.append({
    "name": "B3LYP-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_B3LYP": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    B3LYP Exchange-Correlation Functional plus XDM dispersion (hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.6356,
            'a2': 1.5119,
            'vol': 'b3lyp',
        }
    }
})
funcs.append({
    "name": "B3PW91-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_B3PW91": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    B3PW91 Exchange-Correlation Functional plus XDM dispersion (hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.6002,
            'a2': 1.4043,
            'vol': '0.2',
        }
    }
})
funcs.append({
    "name": "PBE0-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_PBEH": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    PBE0 Exchange-Correlation Functional plus XDM dispersion (hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.4186,
            'a2': 2.6791,
            'vol': 'pbe0',
        }
    }
})
funcs.append({
    "name": "CAM-B3LYP-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_CAM_B3LYP": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    CAM-B3LYP Exchange-Correlation Functional plus XDM dispersion (rs-hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.3248,
            'a2': 2.8607,
            'vol': 'camb3lyp',
        }
    }
})
funcs.append({
    "name": "BHANDHLYP-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_BHANDHLYP": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    BHANDHLYP Exchange-Correlation Functional plus XDM dispersion (hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.5610,
            'a2': 1.9894,
            'vol': 'bhah',
        }
    }
})
funcs.append({
    "name": "TPSS-XDM",
    "x_functionals": {
        "MGGA_X_TPSS": {}
    },
    "c_functionals": {
        "MGGA_C_TPSS": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    TPSS Exchange-Correlation Functional plus XDM dispersion (meta-GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.6612,
            'a2': 1.5111,
            'vol': '0.0',
        }
    }
})
funcs.append({
    "name": "HSE06-XDM",
    "xc_functionals": {
        "HYB_GGA_XC_HSE06": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    HSE06 Exchange-Correlation Functional plus XDM dispersion (rs-hybrid).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.3691,
            'a2': 2.8793,
            'vol': '0.0',
        }
    }
})
funcs.append({
    "name": "BP86-XDM",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_P86": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    BP86 Exchange-Correlation Functional plus XDM dispersion (GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.8538,
            'a2': 0.6415,
            'vol': '0.0',
        }
    }
})
funcs.append({
    "name": "B86bPBE-XDM",
    "x_functionals": {
        "GGA_X_B86_MGC": {}
    },
    "c_functionals": {
        "GGA_C_PBE": {}
    },
    "citation":
    '    A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    B86bPBE Exchange-Correlation Functional plus XDM dispersion (GGA).\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            'a1': 0.7839,
            'a2': 1.2544,
            'vol': '0.0',
        }
    }
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional
