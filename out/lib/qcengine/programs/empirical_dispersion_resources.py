"""Collect empirical dispersion parameters."""

import collections
import copy

from ..exceptions import InputError

## ==> Dispersion Aliases and Parameters <== ##

# The dashcoeff dict below defines the -D parameters for most of the DFT methods. Some 'd2' are
#   taken from already defined functionals in psi4. Other parameters taken from the
#   references indicated in the "citation" parameter. The remainder of the parameters are
#   from http://toc.uni-muenster.de/DFTD3/ on September 25, 2012, with the dict keys
#   translated from Turbomole to Psi4 functional names.
#   Note: most DSD-functionals have their dispersion correction defined in dict_dh_funcs.py
dashcoeff = {
    "d2": {
        "formal": "D2",
        "alias": ["d"],
        "description": "    Grimme's -D2 Dispersion Correction",
        "citation": "    Grimme, S. (2006), J. Comp. Chem., 27: 1787-1799\n",
        "bibtex": "Grimme:2006:1787",
        "default": collections.OrderedDict([("s6", 1.0), ("alpha6", 20.0), ("sr6", 1.1)]),
        "definitions": {
            "blyp": {"params": {"s6": 1.2, "alpha6": 20.0, "sr6": 1.1}},
            "bp86": {"params": {"s6": 1.05, "alpha6": 20.0, "sr6": 1.1}},
            "b97": {"params": {"s6": 1.25, "alpha6": 20.0, "sr6": 1.1}},  # formerly b97-d
            "revpbe": {
                "params": {"s6": 1.25, "alpha6": 20.0, "sr6": 1.1},
                "citation": "    S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010\n",
            },
            "pbe": {"params": {"s6": 0.75, "alpha6": 20.0, "sr6": 1.1}},
            "tpss": {"params": {"s6": 1.0, "alpha6": 20.0, "sr6": 1.1}},
            "b3lyp": {"params": {"s6": 1.05, "alpha6": 20.0, "sr6": 1.1}},
            "pbe0": {"params": {"s6": 0.6, "alpha6": 20.0, "sr6": 1.1}},
            "pw6b95": {"params": {"s6": 0.5, "alpha6": 20.0, "sr6": 1.1}},
            "tpss0": {
                "params": {"s6": 0.85, "alpha6": 20.0, "sr6": 1.1},
                "citation": "    S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010\n",
            },
            "b2plyp": {"params": {"s6": 0.55, "alpha6": 20.0, "sr6": 1.1}},
            "b2gpplyp": {"params": {"s6": 0.4, "alpha6": 20.0, "sr6": 1.1}},
            "dsd-blyp": {"params": {"s6": 0.41, "alpha6": 60.0, "sr6": 1.1}},
            "core-dsd-blyp": {"params": {"s6": 0.41, "alpha6": 60.0, "sr6": 1.1}},
        },
    },
    "d3zero": {
        "formal": "D3",
        "alias": ["d3"],
        "description": "    Grimme's -D3 (zero-damping) Dispersion Correction",
        "citation": "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104\n",
        "bibtex": "Grimme:2010:154104",
        "default": collections.OrderedDict([("s6", 1.0), ("s8", 0.0), ("sr6", 1.0), ("alpha6", 14.0), ("sr8", 1.0)]),
        "definitions": {
            # S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010
            # 'b2plyp'      : {'s6': 0.5,  's8': 1.000, 'sr6': 1.332, 'alpha6': 14.0}, # superseded by a value below.
            "pw6b95": {"params": {"s6": 1.0, "s8": 0.862, "sr6": 1.532, "alpha6": 14.0, "sr8": 1.000}},
            "b97": {"params": {"s6": 1.0, "s8": 0.909, "sr6": 0.892, "alpha6": 14.0, "sr8": 1.000}},  # formerly b97-d
            "b3lyp": {"params": {"s6": 1.0, "s8": 1.703, "sr6": 1.261, "alpha6": 14.0, "sr8": 1.000}},
            "blyp": {"params": {"s6": 1.0, "s8": 1.682, "sr6": 1.094, "alpha6": 14.0, "sr8": 1.000}},
            "tpss0": {"params": {"s6": 1.0, "s8": 1.242, "sr6": 1.252, "alpha6": 14.0, "sr8": 1.000}},
            "pbe0": {"params": {"s6": 1.0, "s8": 0.928, "sr6": 1.287, "alpha6": 14.0, "sr8": 1.000}},
            "tpss": {"params": {"s6": 1.0, "s8": 1.105, "sr6": 1.166, "alpha6": 14.0, "sr8": 1.000}},
            "pbe": {"params": {"s6": 1.0, "s8": 0.722, "sr6": 1.217, "alpha6": 14.0, "sr8": 1.000}},
            "bp86": {"params": {"s6": 1.0, "s8": 1.683, "sr6": 1.139, "alpha6": 14.0, "sr8": 1.000}},
            # Later references
            "rpw86pbe": {
                "params": {"s6": 1.0, "s8": 0.901, "sr6": 1.224, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    S. Grimme, S. Ehrlich, L. Goerigk, J. Comput. Chem. 32, 1456-1465, 2011\n",
            },
            "b2plyp": {
                "params": {"s6": 0.64, "s8": 1.022, "sr6": 1.427, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n",
            },
            "b2gpplyp": {
                "params": {"s6": 0.56, "s8": 0.760, "sr6": 1.586, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n",
            },
            "pwpb95": {
                "params": {"s6": 0.82, "s8": 0.705, "sr6": 1.557, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n",
            },
            "dsd-blyp": {
                "params": {"s6": 0.50, "s8": 0.705, "sr6": 1.569, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n",
            },
            "bpbe": {
                "params": {"s6": 1.0, "s8": 2.033, "sr6": 1.087, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "opbe": {
                "params": {"s6": 1.0, "s8": 2.055, "sr6": 0.837, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "olyp": {
                "params": {"s6": 1.0, "s8": 1.764, "sr6": 0.806, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpwlyp": {
                "params": {"s6": 1.0, "s8": 1.098, "sr6": 1.239, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "bmk": {
                "params": {"s6": 1.0, "s8": 2.168, "sr6": 1.931, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "m05-2x": {
                "params": {"s6": 1.0, "s8": 0.00, "sr6": 1.417, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "m05": {
                "params": {"s6": 1.0, "s8": 0.595, "sr6": 1.373, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "m06-2x": {
                "params": {"s6": 1.0, "s8": 0.00, "sr6": 1.619, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "m06": {
                "params": {"s6": 1.0, "s8": 0.00, "sr6": 1.325, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "m06-l": {
                "params": {"s6": 1.0, "s8": 0.00, "sr6": 1.581, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b3pw91": {
                "params": {"s6": 1.0, "s8": 1.775, "sr6": 1.176, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "bhlyp": {
                "params": {"s6": 1.0, "s8": 1.442, "sr6": 1.370, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b1b95": {
                "params": {"s6": 1.0, "s8": 1.868, "sr6": 1.613, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpw1b95": {
                "params": {"s6": 1.0, "s8": 1.118, "sr6": 1.605, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpwb1k": {
                "params": {"s6": 1.0, "s8": 1.061, "sr6": 1.671, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "tpssh": {
                "params": {"s6": 1.0, "s8": 1.219, "sr6": 1.223, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "scan": {
                "params": {"s6": 1.0, "s8": 0.000, "sr6": 1.324, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    J.G. Brandenburg, J. E. Bates, J. Sun, J.P. Perdew, Phys. Rev. B 94, 115144, 2016\n",
            },
            "m11-l": {
                "params": {"s6": 1.0, "s8": 1.1129, "sr6": 2.3933, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "mn12-l": {
                "params": {"s6": 1.0, "s8": 0.9622, "sr6": 2.2329, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "n12": {
                "params": {"s6": 1.0, "s8": 2.3916, "sr6": 1.3493, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "sogga11-x": {
                "params": {"s6": 1.0, "s8": 1.8151, "sr6": 1.5431, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "m11": {
                "params": {"s6": 1.0, "s8": 0.6244, "sr6": 2.2300, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "mn12-sx": {
                "params": {"s6": 1.0, "s8": 0.8205, "sr6": 1.9572, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "n12-sx": {
                "params": {"s6": 1.0, "s8": 1.4713, "sr6": 1.4597, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "hse06": {
                "params": {"s6": 1.0, "s8": 0.1090, "sr6": 1.1290, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    J. Moellmann, S. Grimme, J. Chem. Phys. C 118, 7615-7621, 2014\n",
            },
            "wb97x": {
                "params": {"s6": 1.0, "s8": 1.0, "sr6": 1.281, "alpha6": 14.0, "sr8": 1.094},
                "citation": "    Y.-S. Lin, G.-D. Li, S.-P. Mao, J.-D. Chai, Chem. Theory Comput. 9, 263-272, 2013\n",
            },
            "pbehpbe": {
                "params": {"s6": 1.0, "s8": 1.4010, "sr6": 1.5703, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "rpbe": {
                "params": {"s6": 1.0, "s8": 0.5140, "sr6": 0.8720, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "xlyp": {
                "params": {"s6": 1.0, "s8": 0.7447, "sr6": 0.9384, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pw91p86": {
                "params": {"s6": 1.0, "s8": 0.8747, "sr6": 2.1040, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpwpw91": {
                "params": {"s6": 1.0, "s8": 1.9467, "sr6": 1.3725, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hcth407": {
                "params": {"s6": 1.0, "s8": 2.7694, "sr6": 4.0426, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pkzb": {
                "params": {"s6": 1.0, "s8": 0.0000, "sr6": 0.6327, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpss": {
                "params": {"s6": 1.0, "s8": 1.3666, "sr6": 1.3491, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "thcth": {
                "params": {"s6": 1.0, "s8": 0.5662, "sr6": 0.9320, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mn15-l": {
                "params": {"s6": 1.0, "s8": 0.0000, "sr6": 3.3388, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b3p86": {
                "params": {"s6": 1.0, "s8": 1.1961, "sr6": 1.1897, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b1p96": {
                "params": {"s6": 1.0, "s8": 1.1209, "sr6": 1.1815, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b1lyp": {
                "params": {"s6": 1.0, "s8": 1.9467, "sr6": 1.3725, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw1lyp": {
                "params": {"s6": 1.0, "s8": 1.9529, "sr6": 2.0512, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw1pw91": {
                "params": {"s6": 1.0, "s8": 1.4758, "sr6": 1.2892, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pw1pw": {
                "params": {"s6": 1.0, "s8": 1.1786, "sr6": 1.4968, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw1kcis": {
                "params": {"s6": 1.0, "s8": 2.2917, "sr6": 1.7231, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpwkcis1k": {
                "params": {"s6": 1.0, "s8": 1.7553, "sr6": 1.4853, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbeh1pbe": {
                "params": {"s6": 1.0, "s8": 1.0430, "sr6": 1.3719, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbe1kcis": {
                "params": {"s6": 1.0, "s8": 1.7934, "sr6": 3.6355, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "x3lyp": {
                "params": {"s6": 1.0, "s8": 0.2990, "sr6": 1.0000, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "o3lyp": {
                "params": {"s6": 1.0, "s8": 1.8058, "sr6": 1.4060, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b97-1": {
                "params": {"s6": 1.0, "s8": 1.6418, "sr6": 3.7924, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b97-2": {
                "params": {"s6": 1.0, "s8": 2.4661, "sr6": 1.7066, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b98": {
                "params": {"s6": 1.0, "s8": 1.9078, "sr6": 2.6895, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hiss": {
                "params": {"s6": 1.0, "s8": 0.7615, "sr6": 1.3338, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hse03": {
                "params": {"s6": 1.0, "s8": 1.0156, "sr6": 1.3944, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpssh": {
                "params": {"s6": 1.0, "s8": 1.2504, "sr6": 1.3224, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpss0": {
                "params": {"s6": 1.0, "s8": 1.0649, "sr6": 1.2881, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "tpss1kcis": {
                "params": {"s6": 1.0, "s8": 2.0902, "sr6": 1.7729, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "thcthhyb": {
                "params": {"s6": 1.0, "s8": 1.6302, "sr6": 1.5001, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "m08-hx": {
                "params": {"s6": 1.0, "s8": 0.0000, "sr6": 1.6247, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "lcwhpbe": {
                "params": {"s6": 1.0, "s8": 1.2797, "sr6": 1.3846, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw2plyp": {
                "params": {"s6": 0.66, "s8": 0.7529, "sr6": 1.5527, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbe0-dh": {
                "params": {"s6": 0.840, "s8": 0.748, "sr6": 1.394, "alpha6": 14.0, "sr8": 1.000},
                "citation": "    D. Bousquet, E. Bremond, J. C. Sancho-Garcia, I. Ciofini, C. Adamo, Theor. Chem. Acc. 134, 1602, 2015\n",
            },
            # https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/functionals
            "pwb6k": {"params": {"s6": 1.0, "s8": 0.550, "sr6": 1.660, "alpha6": 14.0, "sr8": 1.000}},
            "revpbe": {"params": {"s6": 1.0, "s8": 1.010, "sr6": 0.923, "alpha6": 14.0, "sr8": 1.000}},
            "bop": {"params": {"s6": 1.0, "s8": 1.975, "sr6": 0.929, "alpha6": 14.0, "sr8": 1.000}},
            "otpss": {"params": {"s6": 1.0, "s8": 1.494, "sr6": 1.128, "alpha6": 14.0, "sr8": 1.000}},
            "pbe38": {"params": {"s6": 1.0, "s8": 0.998, "sr6": 1.333, "alpha6": 14.0, "sr8": 1.000}},
            "pbesol": {"params": {"s6": 1.0, "s8": 0.612, "sr6": 1.345, "alpha6": 14.0, "sr8": 1.000}},
            "revssb": {"params": {"s6": 1.0, "s8": 0.560, "sr6": 1.221, "alpha6": 14.0, "sr8": 1.000}},
            "ssb": {"params": {"s6": 1.0, "s8": 0.663, "sr6": 1.215, "alpha6": 14.0, "sr8": 1.000}},
            "cam-b3lyp": {"params": {"s6": 1.0, "s8": 1.217, "sr6": 1.378, "alpha6": 14.0, "sr8": 1.000}},
            "wpbe": {"params": {"s6": 1.0, "s8": 1.279, "sr6": 1.355, "alpha6": 14.0, "sr8": 1.000}},  # formerly lcwpbe
            "m06-hf": {"params": {"s6": 1.0, "s8": 0.00, "sr6": 1.446, "alpha6": 14.0, "sr8": 1.000}},
            "hcth120": {"params": {"s6": 1.0, "s8": 1.206, "sr6": 1.221, "alpha6": 14.0, "sr8": 1.000}},
            "ptpss": {"params": {"s6": 0.75, "s8": 0.879, "sr6": 1.541, "alpha6": 14.0, "sr8": 1.000}},
            "revpbe0": {"params": {"s6": 1.0, "s8": 0.792, "sr6": 0.949, "alpha6": 14.0, "sr8": 1.000}},
            "revpbe38": {"params": {"s6": 1.0, "s8": 0.862, "sr6": 1.021, "alpha6": 14.0, "sr8": 1.000}},
            # unreferenced
            "hf": {"params": {"s6": 1.0, "s8": 1.746, "sr6": 1.158, "alpha6": 14.0, "sr8": 1.000}},
        },
    },
    "d3bj": {
        "formal": "D3(BJ)",
        "alias": [],
        "description": "    Grimme's -D3 (BJ-damping) Dispersion Correction",
        "citation": "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456\n",
        "bibtex": "Grimme:2011:1456",
        "default": collections.OrderedDict([("s6", 1.0), ("s8", 1.0), ("a1", 0.0), ("a2", 1.0)]),
        "definitions": {
            # S. Grimme, S. Ehrlich, L. Goerigk, J. Comput. Chem. 7, 3297-3305, 2011
            "pbe": {"params": {"s6": 1.000, "s8": 0.7875, "a1": 0.4289, "a2": 4.4407}},
            "revpbe": {"params": {"s6": 1.000, "s8": 2.3550, "a1": 0.5238, "a2": 3.5016}},
            "blyp": {"params": {"s6": 1.000, "s8": 2.6996, "a1": 0.4298, "a2": 4.2359}},
            "bp86": {"params": {"s6": 1.000, "s8": 3.2822, "a1": 0.3946, "a2": 4.8516}},
            "rpw86pbe": {"params": {"s6": 1.000, "s8": 1.3845, "a1": 0.4613, "a2": 4.5062}},
            "b97": {"params": {"s6": 1.000, "s8": 2.2609, "a1": 0.5545, "a2": 3.2297}},  # formerly b97-d
            "tpss": {"params": {"s6": 1.000, "s8": 1.9435, "a1": 0.4535, "a2": 4.4752}},
            "b3lyp": {"params": {"s6": 1.000, "s8": 1.9889, "a1": 0.3981, "a2": 4.4211}},
            "pw6b95": {"params": {"s6": 1.000, "s8": 0.7257, "a1": 0.2076, "a2": 6.3750}},
            "pbe0": {"params": {"s6": 1.000, "s8": 1.2177, "a1": 0.4145, "a2": 4.8593}},
            "tpss0": {"params": {"s6": 1.000, "s8": 1.2576, "a1": 0.3768, "a2": 4.5865}},
            # 'b2plyp'      : {'params': {'s6': 0.500, 's8':  1.0860, 'a1':  0.3451, 'a2': 4.7735}}, # superseded by a value below.
            "hf": {"params": {"s6": 1.000, "s8": 0.9171, "a1": 0.3385, "a2": 2.8830}},
            # Later references
            "bpbe": {
                "params": {"s6": 1.000, "s8": 4.0728, "a1": 0.4567, "a2": 4.3908},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "opbe": {
                "params": {"s6": 1.000, "s8": 3.3816, "a1": 0.5512, "a2": 2.9444},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "olyp": {
                "params": {"s6": 1.000, "s8": 2.6205, "a1": 0.5299, "a2": 2.8065},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpwlyp": {
                "params": {"s6": 1.000, "s8": 2.0077, "a1": 0.4831, "a2": 4.5323},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b3pw91": {
                "params": {"s6": 1.000, "s8": 2.8524, "a1": 0.4312, "a2": 4.4693},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "bhlyp": {
                "params": {"s6": 1.000, "s8": 1.0354, "a1": 0.2793, "a2": 4.9615},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b1b95": {
                "params": {"s6": 1.000, "s8": 1.4507, "a1": 0.2092, "a2": 5.5545},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpw1b95": {
                "params": {"s6": 1.000, "s8": 1.0508, "a1": 0.1955, "a2": 6.4177},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "mpwb1k": {
                "params": {"s6": 1.000, "s8": 0.9499, "a1": 0.1474, "a2": 6.6223},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "tpssh": {
                "params": {"s6": 1.000, "s8": 2.2382, "a1": 0.4529, "a2": 4.6550},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "bmk": {
                "params": {"s6": 1.000, "s8": 2.0860, "a1": 0.1940, "a2": 5.9197},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b2plyp": {
                "params": {"s6": 0.640, "s8": 0.9147, "a1": 0.3065, "a2": 5.0570},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "b2gpplyp": {
                "params": {"s6": 0.560, "s8": 0.2597, "a1": 0.0000, "a2": 6.3332},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "pwpb95": {
                "params": {"s6": 0.820, "s8": 0.2904, "a1": 0.0000, "a2": 7.3141},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "dsd-blyp": {
                "params": {"s6": 0.500, "s8": 0.2130, "a1": 0.0000, "a2": 6.0519},
                "citation": "    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n",
            },
            "hse06": {
                "params": {"s6": 1.000, "s8": 2.3100, "a1": 0.3830, "a2": 5.6850},
                "citation": "    J. Moellmann, S. Grimme, J. Chem. Phys. C 118, 7615-7621, 2014\n",
            },
            "pw91": {
                "params": {"s6": 1.000, "s8": 1.9598, "a1": 0.6319, "a2": 4.5718},
                "citation": "    J.R. Reimers et al., Proc. Natl. Acad. Sci. USA 112, E6101-E6110, 2015\n",
            },
            "m11-l": {
                "params": {"s6": 1.000, "s8": 0.4446, "a1": 0.0000, "a2": 7.2496},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "mn12-l": {
                "params": {"s6": 1.000, "s8": 2.2674, "a1": 0.0000, "a2": 9.1494},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "n12": {
                "params": {"s6": 1.000, "s8": 4.8491, "a1": 0.3842, "a2": 5.3545},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "sogga11-x": {
                "params": {"s6": 1.000, "s8": 1.1426, "a1": 0.1330, "a2": 5.7381},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "m11": {
                "params": {"s6": 1.000, "s8": 2.8112, "a1": 0.0000, "a2": 10.1389},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "mn12-sx": {
                "params": {"s6": 1.000, "s8": 1.1674, "a1": 0.0983, "a2": 8.0259},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "n12-sx": {
                "params": {"s6": 1.000, "s8": 2.4900, "a1": 0.3283, "a2": 5.7898},
                "citation": "    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n",
            },
            "scan": {
                "params": {"s6": 1.000, "s8": 0.0000, "a1": 0.5380, "a2": 5.4200},
                "citation": "    J.G. Brandenburg, J. E. Bates, J. Sun, J.P. Perdew, Phys. Rev. B 94, 115144, 2016\n",
            },
            "pbehpbe": {
                "params": {"s6": 1.000, "s8": 1.1152, "a1": 0.0000, "a2": 4.4407},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "rpbe": {
                "params": {"s6": 1.000, "s8": 0.8318, "a1": 0.1820, "a2": 4.0094},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "xlyp": {
                "params": {"s6": 1.000, "s8": 1.5669, "a1": 0.0809, "a2": 5.3166},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpwpw91": {
                "params": {"s6": 1.000, "s8": 0.3168, "a1": 0.3168, "a2": 4.7732},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hcth407": {
                "params": {"s6": 1.000, "s8": 0.6490, "a1": 0.0000, "a2": 4.8162},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpss": {
                "params": {"s6": 1.000, "s8": 1.4023, "a1": 0.4426, "a2": 4.4723},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "thcth": {
                "params": {"s6": 1.000, "s8": 1.2626, "a1": 0.0000, "a2": 5.6162},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b3p86": {
                "params": {"s6": 1.000, "s8": 3.3211, "a1": 0.4601, "a2": 4.9294},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b1p86": {
                "params": {"s6": 1.000, "s8": 3.5681, "a1": 0.4724, "a2": 4.9858},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b1lyp": {
                "params": {"s6": 1.000, "s8": 2.1167, "a1": 0.1986, "a2": 5.3875},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw1pw91": {
                "params": {"s6": 1.000, "s8": 1.8744, "a1": 0.3342, "a2": 4.9819},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw1kcis": {
                "params": {"s6": 1.000, "s8": 1.0893, "a1": 0.0576, "a2": 5.5314},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpwkcis1k": {
                "params": {"s6": 1.000, "s8": 1.2875, "a1": 0.0855, "a2": 5.8961},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbeh1pbe": {
                "params": {"s6": 1.000, "s8": 1.4877, "a1": 0.0000, "a2": 7.0385},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbe1kcis": {
                "params": {"s6": 1.000, "s8": 0.7688, "a1": 0.0000, "a2": 6.2794},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "x3lyp": {
                "params": {"s6": 1.000, "s8": 1.5744, "a1": 0.2022, "a2": 5.4184},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "o3lyp": {
                "params": {"s6": 1.000, "s8": 1.8171, "a1": 0.0963, "a2": 5.9940},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b97-1": {
                "params": {"s6": 1.000, "s8": 0.4814, "a1": 0.0000, "a2": 6.2279},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b97-2": {
                "params": {"s6": 1.000, "s8": 0.9448, "a1": 0.0000, "a2": 5.4603},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "b98": {
                "params": {"s6": 1.000, "s8": 0.7086, "a1": 0.0000, "a2": 6.0672},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hiss": {
                "params": {"s6": 1.000, "s8": 1.6112, "a1": 0.0000, "a2": 7.3539},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "hse03": {
                "params": {"s6": 1.000, "s8": 1.1243, "a1": 0.0000, "a2": 6.8889},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpssh": {
                "params": {"s6": 1.000, "s8": 1.4076, "a1": 0.2660, "a2": 5.3761},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "revtpss0": {
                "params": {"s6": 1.000, "s8": 1.6151, "a1": 0.2218, "a2": 5.7985},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "tpss1kcis": {
                "params": {"s6": 1.000, "s8": 1.0542, "a1": 0.0000, "a2": 6.0201},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "thcthhyb": {
                "params": {"s6": 1.000, "s8": 0.9585, "a1": 0.0000, "a2": 6.2303},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mn15": {
                "params": {"s6": 1.000, "s8": 2.0971, "a1": 0.7862, "a2": 7.5923},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "lcwhpbe": {
                "params": {"s6": 1.000, "s8": 1.1908, "a1": 0.2746, "a2": 5.3157},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "mpw2plyp": {
                "params": {"s6": 1.000, "s8": 0.6223, "a1": 0.4105, "a2": 5.0136},
                "citation": "    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n",
            },
            "pbe0-dh": {
                "params": {"s6": 0.840, "s8": 0.095, "a1": 0.000, "a2": 6.102},
                "citation": "    D. Bousquet, E. Bremond, J. C. Sancho-Garcia, I. Ciofini, C. Adamo, Theor. Chem. Acc. 134, 1602, 2015\n",
            },
            # https://www.chemie.u{'params': ni-bonn.de/pctc/mulliken-center/software/dft-d3/functiona}lsbj
            "bop": {"params": {"s6": 1.000, "s8": 3.295, "a1": 0.4870, "a2": 3.5043}},
            "cam-b3lyp": {"params": {"s6": 1.000, "s8": 2.0674, "a1": 0.3708, "a2": 5.4743}},
            "wpbe": {"params": {"s6": 1.000, "s8": 1.8541, "a1": 0.3919, "a2": 5.0897}},  # formerly lcwpbe
            "otpss": {"params": {"s6": 1.000, "s8": 2.7495, "a1": 0.4634, "a2": 4.3153}},
            "pbe38": {"params": {"s6": 1.000, "s8": 1.4623, "a1": 0.3995, "a2": 5.1405}},
            "pbesol": {"params": {"s6": 1.000, "s8": 2.9491, "a1": 0.4466, "a2": 6.1742}},
            "ptpss": {"params": {"s6": 0.750, "s8": 0.2804, "a1": 0.000, "a2": 6.5745}},
            "pwb6k": {"params": {"s6": 1.000, "s8": 0.9383, "a1": 0.1805, "a2": 7.7627}},
            "revssb": {"params": {"s6": 1.000, "s8": 0.4389, "a1": 0.4720, "a2": 4.0986}},
            "ssb": {"params": {"s6": 1.000, "s8": -0.1744, "a1": -0.0952, "a2": 5.2170}},
            "hcth120": {"params": {"s6": 1.000, "s8": 1.0821, "a1": 0.3563, "a2": 4.3359}},
            "revpbe0": {"params": {"s6": 1.000, "s8": 1.7588, "a1": 0.4679, "a2": 3.7619}},
            "revpbe38": {"params": {"s6": 1.000, "s8": 1.4760, "a1": 0.4309, "a2": 3.9446}},
            # unreferenced
            # special HF/DFT with eBSSE correction
            "hf/mixed": {"params": {"s6": 1.000, "s8": 3.9027, "a1": 0.5607, "a2": 4.5622}},
            "hf/sv": {"params": {"s6": 1.000, "s8": 2.1849, "a1": 0.4249, "a2": 4.2783}},
            "hf/minis": {"params": {"s6": 1.000, "s8": 0.9841, "a1": 0.1702, "a2": 3.8506}},
            "b3lyp/6-31gd": {"params": {"s6": 1.000, "s8": 4.0672, "a1": 0.5014, "a2": 4.8409}},
            # special HF-D3-gCP-SRB/MINIX parametrization
            "hf3c": {"params": {"s6": 1.000, "s8": 0.8777, "a1": 0.4171, "a2": 2.9149}},
            # special HF-D3-gCP-SRB2/ECP-2G parametrization
            "hf3cv": {"params": {"s6": 1.000, "s8": 0.5022, "a1": 0.3063, "a2": 3.9856}},
            # special PBEh-D3-gCP/def2-mSVP parametrization
            "pbeh3c": {"params": {"s6": 1.000, "s8": 0.0000, "a1": 0.4860, "a2": 4.5000}},
            "core-dsd-blyp": {"params": {"s6": 0.500, "s8": 0.2130, "a1": 0.0000, "a2": 6.0519}},
        },
    },
    "d3mzero": {
        "formal": "D3M",
        "alias": ["d3m"],
        "description": "    Grimme's -D3 (zero-damping, short-range refitted) Dispersion Correction",
        "citation": "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104\n"
        + "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197\n",
        "bibtex": "Grimme:2010:154104",
        "default": collections.OrderedDict([("s6", 1.0), ("s8", 1.0), ("sr6", 1.0), ("beta", 1.0)]),
        "definitions": {
            "b2plyp": {"params": {"s6": 0.640, "s8": 0.717543, "sr6": 1.313134, "beta": 0.016035}},
            "b3lyp": {"params": {"s6": 1.000, "s8": 1.532981, "sr6": 1.338153, "beta": 0.013988}},
            "b97": {"params": {"s6": 1.000, "s8": 1.020078, "sr6": 1.151808, "beta": 0.035964}},  # formerly b97-d
            "blyp": {"params": {"s6": 1.000, "s8": 1.841686, "sr6": 1.279637, "beta": 0.014370}},
            "bp86": {"params": {"s6": 1.000, "s8": 1.945174, "sr6": 1.233460, "beta": 0.000000}},
            "pbe": {"params": {"s6": 1.000, "s8": 0.000000, "sr6": 2.340218, "beta": 0.129434}},
            "pbe0": {"params": {"s6": 1.000, "s8": 0.000081, "sr6": 2.077949, "beta": 0.116755}},
            "wpbe": {"params": {"s6": 1.000, "s8": 1.280619, "sr6": 1.366361, "beta": 0.003160}},  # formerly lcwpbe
            "sapt0": {"params": {"s6": 1.000, "s8": 0.885517, "sr6": 1.383214, "beta": 0.075488}},  # JBS 01/2021
            "hf": {"params": {"s6": 1.000, "s8": 0.885517, "sr6": 1.383214, "beta": 0.075488}},  # JBS 01/2021
        },
    },
    "d3mbj": {
        "formal": "D3M(BJ)",
        "alias": [],
        "description": "    Grimme's -D3 (BJ-damping, short-range refitted) Dispersion Correction",
        "citation": "    Grimme S.; Ehrlich S.; Goerigk L. (2011), J. Comput. Chem., 32: 1456\n"
        + "    Smith, D. G. A.; Burns, L. A.; Patkowski, K.; Sherrill, C. D. (2016), J. Phys. Chem. Lett.; 7: 2197\n",
        "bibtex": "Grimme:2011:1456",
        "default": collections.OrderedDict([("s6", 1.0), ("s8", 1.0), ("a1", 1.0), ("a2", 1.0)]),
        "definitions": {
            "b2plyp": {"params": {"s6": 0.640, "s8": 0.672820, "a1": 0.486434, "a2": 3.656466}},
            "b3lyp": {"params": {"s6": 1.000, "s8": 1.466677, "a1": 0.278672, "a2": 4.606311}},
            "b97": {"params": {"s6": 1.000, "s8": 1.206988, "a1": 0.240184, "a2": 3.864426}},  # formerly b97-d
            "blyp": {"params": {"s6": 1.000, "s8": 1.875007, "a1": 0.448486, "a2": 3.610679}},
            "bp86": {"params": {"s6": 1.000, "s8": 3.140281, "a1": 0.821850, "a2": 2.728151}},
            "pbe": {"params": {"s6": 1.000, "s8": 0.358940, "a1": 0.012092, "a2": 5.938951}},
            "pbe0": {"params": {"s6": 1.000, "s8": 0.528823, "a1": 0.007912, "a2": 6.162326}},
            "wpbe": {"params": {"s6": 1.000, "s8": 0.906564, "a1": 0.563761, "a2": 3.593680}},  # formerly lcwpbe
            "sapt0": {"params": {"s6": 1.000, "s8": 0.713190, "a1": 0.079541, "a2": 3.627854}},  # JBS 01/2021
            "hf": {"params": {"s6": 1.000, "s8": 0.713190, "a1": 0.079541, "a2": 3.627854}},  # JBS 01/2021
        },
    },
    "nl": {
        "formal": "NL",
        "alias": [],
        "description": "    Grimme's -NL (DFT plus VV10 correlation) ",
        "citation": "    Hujo, W.; Grimme, S; (2011), J. Chem. Theory Comput.; 7:3866\n",
        "bibtex": "Hujo:2011:3866",
        "default": collections.OrderedDict([("b", 1.0), ("c", 0.0093)]),
        "definitions": {
            "blyp": {
                "params": {"b": 4.000, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "hf": {
                "params": {"b": 3.900, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },  # not implemented
            "revpbe": {
                "params": {"b": 3.700, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "revpbe38": {
                "params": {"b": 4.700, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "b3lyp": {
                "params": {"b": 4.800, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "b3pw91": {
                "params": {"b": 4.500, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "revpbe0": {
                "params": {"b": 4.300, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011\n",
            },
            "bp86": {
                "params": {"b": 4.400, "c": 0.0093},
                "citation": "    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n",
            },  # unclear if this is the real origin
            "pbe0": {
                "params": {"b": 6.900, "c": 0.0093},
                "citation": "    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n",
            },  # unclear if this is the real origin
            "pbe": {
                "params": {"b": 6.400, "c": 0.0093},
                "citation": "    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n",
            },  # unclear if this is the real origin
            "tpss0": {
                "params": {"b": 5.500, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme, J. Chem. Theory Comput. 9, 308-315, 2013\n",
            },
            "tpss": {
                "params": {"b": 5.000, "c": 0.0093},
                "citation": "    W. Hujo, S. Grimme, J. Chem. Theory Comput. 9, 308-315, 2013\n",
            },
            "b2gpplyp": {
                "params": {"b": 9.900, "c": 0.0093},
                "citation": "    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016\n",
            },
            "b2plyp": {
                "params": {"b": 7.800, "c": 0.0093},
                "citation": "    J. Calbo, E. Orti, J. C. Sancho-Garcia, J. Arago, J. Chem. Theory Comput. 11, 932-939, 1015\n",
            },
            "pwpb95": {
                "params": {"b": 11.100, "c": 0.0093},
                "citation": "    F. Yu J. Chem. Theory Comput. 10, 4400-4407, 2014\n",
            },
            "revtpss": {
                "params": {"b": 5.465, "c": 0.0093},
                "citation": "    H. Kruse, P. Banas, J. Sponer, JCTC 2018, 10.1021/acs.jctc.8b00643 \n",
            },  # "just accepted"
        },
    },
    "d1": {
        "formal": "D1",
        "alias": [],
        "description": "    Grimme's -D1 Dispersion Correction",
        "citation": "    Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473\n",
        "bibtex": "Grimme:2004:1463",
        "default": collections.OrderedDict([("s6", 1.0)]),
        "definitions": {"***": {"params": {"s6": 1.0}}},
    },
    "chg": {
        "formal": "CHG",
        "alias": [],
        "description": "    Chai and Head-Gordon Dispersion Correction",
        "citation": "    Chai, J.-D.; Head-Gordon, M. (2010), J. Chem. Phys., 132: 6615-6620\n",
        "bibtex": "Chai:2010:6615",
        "default": collections.OrderedDict([("s6", 1.0)]),
        "definitions": {"***": {"params": {"s6": 1.0}}},  # disp->d_ = 6.0;
    },
    "das2009": {
        "formal": "DAS2009",
        "alias": [],
        "description": "    Podeszwa and Szalewicz Dispersion Correction",
        "citation": "    Pernal, K.; Podeszwa, R.; Patkowski, K.; Szalewicz, K. (2009), Phys. Rev. Lett., 103: 263201\n",
        "bibtex": "Pernal:2009:263201",
        "default": collections.OrderedDict([("s6", 1.0)]),
        "definitions": {"***": {"params": {"s6": 1.0}}},
    },
    "das2010": {
        "formal": "DAS2010",
        "alias": [],
        "description": "    Podeszwa and Szalewicz Dispersion Correction",
        "citation": "    Podeszwa, R.; Pernal, K.; Patkowski, K.; Szalewicz, K. (2010), J. Phys. Chem. Lett., 1: 550\n",
        "bibtex": "Podeszwa:2010:550",
        "default": collections.OrderedDict([("s6", 1.0)]),
        "definitions": {"***": {"params": {"s6": 1.0}}},
    },
    "atmgr": {
        "formal": "ATM(GR)",
        "alias": [],
        "description": "    Grimme approximate Axilrod-Teller-Muto 3-body Dispersion Correction",
        "citation": "    Grimme S.; Antony J.; Ehrlich S.; Krieg H. (2010), J. Chem. Phys., 132: 154104\n",
        "bibtex": "Grimme:2010:154104",
        "default": collections.OrderedDict([("alpha6", 14.0)]),
        "definitions": {
            "***": {"params": {"alpha6": 14.0}, "citation": ""}  # alpha6 = 14 => damping parameter used is alp8 = 16
        },
    },
    "dmp2": {
        "formal": "DMP2",
        "alias": [],
        "description": "    Beran's MP2D Dispersion Correction for MP2",
        "citation": "    Rezac, J.; Greenwell, C.; Beran, G. (2018), J. Chem. Theory Comput., 14: 4711-4721\n",
        "bibtex": "Rezac:2018:4711",
        "doi": "10.1021/acs.jctc.8b00548",
        "default": collections.OrderedDict([("s8", 1.187), ("a1", 0.944), ("a2", 0.480), ("rcut", 0.72), ("w", 0.20)]),
        "definitions": {
            "mp2": {"params": {"s8": 1.187, "a1": 0.944, "a2": 0.480, "rcut": 0.72, "w": 0.20}}  # Rezac:2018:4711
        },
    },
    "d4bjeeqatm": {
        "formal": "D4(BJ,EEQ)ATM",
        "alias": ["d4", "d4bj", "d4(bj)"],
        "description": "    Grimme's -D4 (BJ-damping) Dispersion Correction with ATM",
        "citation": "    Caldeweyher, E.; Ehlert, S.; Hansen, A.; Neugebauer, H.; Spicher, S.; Bannwarth, C.; Grimmme, S., J. Chem. Phys. 150, 154122 (2019)\n",
        "bibtex": "Caldeweyher:2019:154122",
        "doi": "10.1063/1.5090222150",
        "default": collections.OrderedDict(
            [("a1", 1.0), ("a2", 1.0), ("alp", 16.0), ("s6", 1.0), ("s8", 1.0), ("s9", 1.0)]
        ),
        "definitions": {
            # D4 parameters loaded below from authoritative source below. Keep a couple for reference
            # "b3lyp": {"params": {"a1": 0.40868035, "a2": 4.53807137, "alp": 16.0, "s6": 1.0, "s8": 2.02929367, "s9": 1.0}},
            # "pbe":   {"params": {"a1": 0.38574991, "a2": 4.80688534, "alp": 16.0, "s6": 1.0, "s8": 0.95948085, "s9": 1.0}},
        },
    },
}


def _get_d4bj_definitions() -> dict:
    """DFTD4 provides access to damping parameters on per functional basis.
    But we want all of them.

    We let DFTD4 take care of finding and loading the parameter file,
    but to get all parameters, we implement the logic to read those
    parameters ourselves again.
    """

    try:
        from dftd4.parameters import get_data_file_name, load_data_base
    except ModuleNotFoundError:
        return {}

    def get_params(entry: dict, base: dict, defaults: list) -> dict:
        """Retrive the parameters from the data base, make sure the default
        values are applied correctly in the process. In case we have multiple
        defaults search for the first of the list defined for this method."""

        for default in defaults:
            try:
                params = base[default].copy()
                params.update(**entry[default])
                params.pop("mbd", None)
                params.pop("damping", None)
                return params
            except KeyError:
                continue

        raise KeyError("No entry for " + method + " in parameter data base")

    try:
        _data_base = load_data_base(get_data_file_name())
    except FileNotFoundError:
        return {}

    try:
        _defaults = _data_base["default"]["d4"]
        _base = _data_base["default"]["parameter"]["d4"]
        _parameters = _data_base["parameter"]
    except KeyError:
        return {}

    definitions = {}

    for method in _parameters:
        try:
            _entry = _parameters[method]["d4"]
            params = get_params(_entry, _base, _defaults)

            citation = params.pop("doi", None)
            # Make Psi4's citation style checker happy
            if citation is not None:
                citation = "    " + citation + "\n"
            definitions[method] = dict(
                params=params,
                citation=citation,
            )
        except KeyError:
            continue

    return definitions


dashcoeff["d4bjeeqatm"]["definitions"].update(_get_d4bj_definitions())


def get_dispersion_aliases():
    """Returns dict where keys are all (lowercased) strings that are
    interpretable as dashlevels, and values are the keys dealiased
    (d->d2) and deformalized (d3(bj)->d3bj) into valid dashcoeff keys.

    """
    alias = {}
    for dash, ddash in dashcoeff.items():
        alias[dash] = dash
        alias[ddash["formal"].lower()] = dash
        for al in ddash["alias"]:
            alias[al] = dash

    return alias


def from_arrays(name_hint=None, level_hint=None, param_tweaks=None, dashcoeff_supplement=None, verbose=1):
    """Use the three paths of empirical dispersion parameter information
    (DFT functional, dispersion correction level, and particular
    parameters) to populate the parameter array and validate a
    "functional-dispersion" label.

    Parameters
    ----------
    name_hint : str, optional
        Name of functional (func only, func & disp, or disp only) for
        which to compute dispersion (e.g., blyp, BLYP-D2, blyp-d3bj,
        blyp-d3(bj), hf+d). Any or all parameters initialized from
        `dashcoeff[dashlevel][functional-without-dashlevel]` or
        `dashcoeff_supplement[dashlevel][functional-with-dashlevel]
        can be overwritten via `param_tweaks`.
    level_hint : str, optional
        Name of dispersion correction to be applied (e.g., d, D2,
        d3(bj), das2010). Must be key in `dashcoeff` or "alias" or
        "formal" to one.
    param_tweaks : list or dict, optional
        Values for the same keys as `dashcoeff[dashlevel]['default']`
        (and same order if list) used to override any or all values
        initialized by `name_hint`.  Extra parameters will error.
    dashcoeff_supplement: dict, optional
        Dictionary of the same structure as `dashcoeff` that contains
        in "definitions" field full functional names, rather than
        fctl less dashlvl. Used to validate dict_builder fctls with
        dispersion or identify disp level from just `name_hint`.
    verbose : int, optional
        Amount of printing.

    Returns
    -------
    dict
        Metadata defining dispersion calculation.

        dashlevel : {'d1', 'd2', 'd3zero', 'd3bj', 'd3mzero', 'd3mbj', 'chg', 'das2009', 'das2010', 'nl', "d4bjeeqatm"}
            Name (de-aliased, de-formalized, lowercase) of dispersion
            correction -- atom data, dispersion model, damping functional
            form -- to be applied. Resolved from `name_hint` and/or
            `level_hint` into a key of `dashparam.dashcoeff`.
        dashparams : dict
            Complete (number and parameter names vary by `dashlevel`)
            set of parameter values defining the flexible parts
            of `dashlevel`. Resolved into a complete set (keys of
            dashcoeff[dashlevel]['default']) from `name_hint` and/or
            `dashcoeff_supplement` and/or user `param_tweaks`.
        fctldash : str
            If `dashparams` for `dashlevel` corresponds to a defined,
            named, untweaked "functional-dashlevel" set, then that
            label. Otherwise, empty string.
        dashparams_citation : str
            If `dashparams` for `dashlevel` corresponds to a defined,
            named, untweaked "functional-dashlevel" set, and that
            parameter set has a literature citation, then that
            citation. Otherwise, empty string.

    Notes
    -----
    * No parameter is required, but somewhere the dispersion level and
      intended parameters must be clear. Explicit contradictory info will
      cause an error.
    * Function intended to be idempotent.

    """
    if verbose > 1:
        print("dftd3.from_arrays HINTS:", name_hint, level_hint, param_tweaks, bool(dashcoeff_supplement))

    # << 0 >> prep
    if dashcoeff_supplement is not None:
        supplement_dashlevel_lookup = {}
        for disp, ddisp in dashcoeff_supplement.items():
            for func, params in ddisp["definitions"].items():
                if params["params"].keys() != dashcoeff[disp]["default"].keys():
                    if verbose > 2:
                        print(
                            "Warning: trouble in dict_builder def:",
                            func,
                            params["params"].keys(),
                            "!=",
                            dashcoeff[disp]["default"].keys(),
                        )
                else:
                    supplement_dashlevel_lookup[func] = disp

    # << 1 >> use name_hint and/or level_hint to determine intended dispersion level
    if name_hint is None and level_hint is None:
        raise InputError(
            """Can't guess -D level without name_hint ({}) or level_hint ({})""".format(name_hint, level_hint)
        )

    if level_hint is None:
        dashlevel_candidate_1 = None
    else:
        level_hint = level_hint.lower()
        try:
            dashlevel_candidate_1 = get_dispersion_aliases()[level_hint]
        except KeyError:
            raise InputError(
                """Requested -D correction level ({}) not among ({})""".format(level_hint, dashcoeff.keys())
            )

    if name_hint is None:
        dashlevel_candidate_2 = None
        name_key = None
        disp_params = {}
    else:
        name_hint = name_hint.lower()
        trial_split = name_hint.rsplit("-", 1)

        if name_hint in get_dispersion_aliases():
            dashlevel_candidate_2 = get_dispersion_aliases()[name_hint]
            if list(dashcoeff[dashlevel_candidate_2]["definitions"]) == ["***"]:
                # case disp-only fctl-indep: chg, atmgr
                name_key = "***"
                disp_params = dashcoeff[dashlevel_candidate_2]["definitions"][name_key]["params"]
            else:
                # case disp-only: d3, d3zero, d3(bj)
                name_key = None
                disp_params = {}
        elif (dashcoeff_supplement is not None) and name_hint in supplement_dashlevel_lookup:
            # case fctldisp: wb97x-d, hf+d
            dashlevel_candidate_2 = supplement_dashlevel_lookup[name_hint]
            name_key = name_hint
            disp_params = dashcoeff_supplement[dashlevel_candidate_2]["definitions"][name_hint]["params"]
        elif (
            len(trial_split) == 2
            and trial_split[1] in get_dispersion_aliases()
            and trial_split[0] in dashcoeff[get_dispersion_aliases()[trial_split[1]]]["definitions"]
        ):
            # case fctldisp: b3lyp-d3, b3lyp-d3zero, b3lyp-d3(bj)
            dashlevel_candidate_2 = get_dispersion_aliases()[trial_split[1]]
            name_key = trial_split[0]
            disp_params = dashcoeff[dashlevel_candidate_2]["definitions"][trial_split[0]]["params"]
        elif (
            len(trial_split) == 2
            and trial_split[1] in get_dispersion_aliases()
            and list(dashcoeff[get_dispersion_aliases()[trial_split[1]]]["definitions"]) == ["***"]
        ):
            # case fctldisp: asdf-chg, pbe-atmgr
            dashlevel_candidate_2 = get_dispersion_aliases()[trial_split[1]]
            name_key = "***"
            disp_params = dashcoeff[dashlevel_candidate_2]["definitions"][name_key]["params"]
        elif (level_hint is not None) and name_hint in dashcoeff[dashlevel_candidate_1]["definitions"]:
            # case fctl: b3lyp
            dashlevel_candidate_2 = None
            name_key = name_hint
            disp_params = dashcoeff[dashlevel_candidate_1]["definitions"][name_hint]["params"]
        elif name_hint == "":
            dashlevel_candidate_2 = None
            name_key = None
            disp_params = {}
        # elif ((level_hint is not None) and list(dashcoeff[get_dispersion_aliases()[level_hint]]['definitions']) == ['***']):
        #    dashlevel_candidate_2 = get_dispersion_aliases()[level_hint]
        #    name_key = '***'
        #    disp_params = dashcoeff[dashlevel_candidate_2]['definitions'][name_key]['params']
        else:
            # dashlevel_candidate_2 = None
            raise InputError("""Can't guess -D correction level from ({})""".format(name_hint))

    disp_params = copy.deepcopy(disp_params)

    if dashlevel_candidate_1 is None and dashlevel_candidate_2 is None:
        raise InputError(
            f"""Can't guess -D correction level from name_hint ({name_hint}) and level_hint ({level_hint})"""
        )
    elif dashlevel_candidate_1 is not None and dashlevel_candidate_2 is not None:
        if dashlevel_candidate_1 != dashlevel_candidate_2:
            raise InputError(
                f"""Inconsistent -D correction level ({dashlevel_candidate_2} != {dashlevel_candidate_1}) from name_hint ({name_hint}) and level_hint ({level_hint})"""
            )
    dashleveleff = dashlevel_candidate_1 or dashlevel_candidate_2

    allowed_params = dashcoeff[dashleveleff]["default"].keys()

    # << 2 >> use name_hint and/or param_tweaks to determine intended dispersion parameters
    if name_hint is None and param_tweaks is None:
        raise InputError(
            """Can't guess -D parameters without name_hint ({}) or param_tweaks ({})""".format(name_hint, param_tweaks)
        )

    if isinstance(param_tweaks, list):
        param_tweaks = dict(zip(allowed_params, param_tweaks))

    if param_tweaks is None:
        param_tweaks = {}

    if not set(param_tweaks.keys()).issubset(allowed_params):
        raise InputError(
            "Requested keys ({}) not among allowed ({}) for dispersion level ({})".format(
                list(param_tweaks.keys()), list(allowed_params), dashleveleff
            )
        )

    disp_params.update(param_tweaks)

    if disp_params.keys() != allowed_params:
        raise InputError(
            "Requested keys ({}) insufficient ({}) for dispersion level ({})".format(
                list(param_tweaks.keys()), list(allowed_params), dashleveleff
            )
        )

    # << 3 >> use final dashlevel and disp_params to determine if a defined "fctl-disp" label exists
    # * plucks any citation for the parameters from definition source
    # * if/elif chooses right label when some fctls have identical param sets
    if (
        (name_hint is not None)
        and (dashcoeff_supplement is not None)
        and (name_key in dashcoeff_supplement[dashleveleff]["definitions"])
        and (disp_params == dashcoeff_supplement[dashleveleff]["definitions"][name_key]["params"])
    ):
        citeff = dashcoeff_supplement[dashleveleff]["definitions"][name_key].get("citation", "")
        if name_key == "***":
            fctldasheff = ""
        else:
            fctldasheff = name_key
    elif name_hint not in [None, ""] and (disp_params == dashcoeff[dashleveleff]["definitions"][name_key]["params"]):
        citeff = dashcoeff[dashleveleff]["definitions"][name_key].get("citation", "")
        if name_key == "***":
            fctldasheff = dashcoeff[dashleveleff]["formal"].lower()
        else:
            fctldasheff = "-".join([name_key, dashcoeff[dashleveleff]["formal"].lower()])
    else:
        if dashcoeff_supplement is not None:
            for func, params in dashcoeff_supplement[dashleveleff]["definitions"].items():
                if disp_params == params["params"]:
                    fctldasheff = func
                    citeff = params.get("citation", "")
                    break

        for func, params in dashcoeff[dashleveleff]["definitions"].items():
            if disp_params == params["params"]:
                fctldasheff = "-".join([func, dashcoeff[dashleveleff]["formal"].lower()])
                citeff = params.get("citation", "")
                break
        else:
            fctldasheff = ""
            citeff = ""

    # TODO right now fctldasheff is empty if undefined. use '-dash' or 'custom dash' instead?
    # TODO right now citation is empty if undefined. remove key or use None or False instead?

    if verbose > 1:
        print(
            f"dftd3.from_arrays RESOLVED: dashlevel={dashleveleff}, dashparams={disp_params}, fctldash={fctldasheff}, dashparams_citation={citeff}"
        )

    return {
        "dashlevel": dashleveleff,
        "dashparams": disp_params,
        "fctldash": fctldasheff,
        "dashparams_citation": citeff,
    }
