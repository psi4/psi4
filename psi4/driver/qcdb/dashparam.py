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
Module to hold and distribute the -D dispersion correction parameters.
"""
from __future__ import absolute_import
from __future__ import print_function
try:
    from .p4regex import *
except ImportError:
    from .exceptions import *


## ==> Dispersion Aliases and Parameters <== ##

# This defines the -D aliases for all of psi4
def get_dispersion_aliases():
    # returns a dictionary consisting of "dispersion_name": "dispersion_type" entries
    # where "dispersion_type" is a valid dtype identifier in lowercase
    dispersion_names = {}
    dispersion_names["d2p4"] = "d2p4"
    dispersion_names["d2gr"] = "d2gr"
    dispersion_names["d3zero"] = "d3zero"
    dispersion_names["d3bj"] = "d3bj"
    dispersion_names["d3mzero"] = "d3mzero"
    dispersion_names["d3mbj"] = "d3mbj"
    dispersion_names["d"] = "d2p4"
    dispersion_names["d2"] = "d2p4"
    dispersion_names["d3"] = "d3zero"
    dispersion_names["d3m"] = "d3mzero"
    dispersion_names["nl"] = "nl"
    return(dispersion_names)

dash_alias = get_dispersion_aliases()

# The dashcoeff dict below defines the -D parameters for most of the DFT methods. 'd2p4' are
#   taken from already defined functionals in psi4. Other parameters taken from the
#   references indicated in the "citation" parameter. The remainder of the parameters are
#   from http://toc.uni-muenster.de/DFTD3/ on September 25, 2012, with the dict keys
#   translated from Turbomole to Psi4 functional names.
#   Note: most DSD-functionals have their dispersion correction defined in dict_dh_funcs.py
dashcoeff = {
    'd2p4': {
        'b97-d'       : {'s6': 1.25},
        'blyp'        : {'s6': 1.20},
        'b3lyp'       : {'s6': 1.05},
        'bp86'        : {'s6': 1.05},
        'tpss'        : {'s6': 1.00},
        'pbe'         : {'s6': 0.75},
        'revpbe'      : {'s6': 1.25, "citation": '    S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010\n'},
        'tpss0'       : {'s6': 0.86, "citation": '    S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010\n'},
        'b2plyp'      : {'s6': 0.55},
        'pbe0'        : {'s6': 0.60},

    },
    'd2gr': {
        'blyp'             : {'s6': 1.2,  'alpha6': 20.0, 'sr6': 1.1},
        'bp86'             : {'s6': 1.05, 'alpha6': 20.0, 'sr6': 1.1},
        'b97-d'            : {'s6': 1.25, 'alpha6': 20.0, 'sr6': 1.1},
        'revpbe'           : {'s6': 1.25, 'alpha6': 20.0, 'sr6': 1.1},
        'pbe'              : {'s6': 0.75, 'alpha6': 20.0, 'sr6': 1.1},
        'tpss'             : {'s6': 1.0,  'alpha6': 20.0, 'sr6': 1.1},
        'b3lyp'            : {'s6': 1.05, 'alpha6': 20.0, 'sr6': 1.1},
        'pbe0'             : {'s6': 0.6,  'alpha6': 20.0, 'sr6': 1.1},
        'pw6b95'           : {'s6': 0.5,  'alpha6': 20.0, 'sr6': 1.1},
        'tpss0'            : {'s6': 0.85, 'alpha6': 20.0, 'sr6': 1.1},
        'b2plyp'           : {'s6': 0.55, 'alpha6': 20.0, 'sr6': 1.1},
        'b2gpplyp'         : {'s6': 0.4,  'alpha6': 20.0, 'sr6': 1.1},
        'dsd-blyp'         : {'s6': 0.41, 'alpha6': 60.0, 'sr6': 1.1},
        'core-dsd-blyp'    : {'s6': 0.41, 'alpha6': 60.0, 'sr6': 1.1},
    },
    'd3zero': {
        # S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys 132, 154104, 2010
        # 'b2plyp'      : {'s6': 0.5,  's8': 1.000, 'sr6': 1.332, 'alpha6': 14.0}, # superseded by a value below.
        'pw6b95'      : {'s6': 1.0,  's8': 0.862,  'sr6': 1.532 , 'sr8': 1.000, 'alpha6': 14.0},
        'b97-d'       : {'s6': 1.0,  's8': 0.909,  'sr6': 0.892 , 'sr8': 1.000, 'alpha6': 14.0},
        'revpbe'      : {'s6': 1.0,  's8': 1.010,  'sr6': 0.923 , 'sr8': 1.000, 'alpha6': 14.0},
        'b3lyp'       : {'s6': 1.0,  's8': 1.703,  'sr6': 1.261 , 'sr8': 1.000, 'alpha6': 14.0},
        'blyp'        : {'s6': 1.0,  's8': 1.682,  'sr6': 1.094 , 'sr8': 1.000, 'alpha6': 14.0},
        'tpss0'       : {'s6': 1.0,  's8': 1.242,  'sr6': 1.252 , 'sr8': 1.000, 'alpha6': 14.0},
        'pbe0'        : {'s6': 1.0,  's8': 0.928,  'sr6': 1.287 , 'sr8': 1.000, 'alpha6': 14.0},
        'tpss'        : {'s6': 1.0,  's8': 1.105,  'sr6': 1.166 , 'sr8': 1.000, 'alpha6': 14.0},
        'pbe'         : {'s6': 1.0,  's8': 0.722,  'sr6': 1.217 , 'sr8': 1.000, 'alpha6': 14.0},
        'bp86'        : {'s6': 1.0,  's8': 1.683,  'sr6': 1.139 , 'sr8': 1.000, 'alpha6': 14.0},
    # Later references
        'rpw86pbe'    : {'s6': 1.0,  's8': 0.901,  'sr6': 1.224 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    S. Grimme, S. Ehrlich, L. Goerigk, J. Comput. Chem. 32, 1456-1465, 2011\n'},
        'b2plyp'      : {'s6': 0.64, 's8': 1.022,  'sr6': 1.427 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n'},
        'b2gpplyp'    : {'s6': 0.56, 's8': 0.760,  'sr6': 1.586 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n'},
        'pwpb95'      : {'s6': 0.82, 's8': 0.705,  'sr6': 1.557 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n'},
        'dsd-blyp'    : {'s6': 0.50, 's8': 0.705,  'sr6': 1.569 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, J. Chem. Theory Comput. 7, 291-309, 2011\n'},
        'bpbe'        : {'s6': 1.0,  's8': 2.033,  'sr6': 1.087 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'opbe'        : {'s6': 1.0,  's8': 2.055,  'sr6': 0.837 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'olyp'        : {'s6': 1.0,  's8': 1.764,  'sr6': 0.806 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpwlyp'      : {'s6': 1.0,  's8': 1.098,  'sr6': 1.239 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'bmk'         : {'s6': 1.0,  's8': 2.168,  'sr6': 1.931 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'm05-2x'      : {'s6': 1.0,  's8': 0.00 ,  'sr6': 1.417 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'm05'         : {'s6': 1.0,  's8': 0.595,  'sr6': 1.373 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'm06-2x'      : {'s6': 1.0,  's8': 0.00 ,  'sr6': 1.619 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'm06'         : {'s6': 1.0,  's8': 0.00 ,  'sr6': 1.325 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'm06-l'       : {'s6': 1.0,  's8': 0.00 ,  'sr6': 1.581 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b3pw91'      : {'s6': 1.0,  's8': 1.775,  'sr6': 1.176 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'bhlyp'       : {'s6': 1.0,  's8': 1.442,  'sr6': 1.370 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b1b95'       : {'s6': 1.0,  's8': 1.868,  'sr6': 1.613 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpw1b95'     : {'s6': 1.0,  's8': 1.118,  'sr6': 1.605 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpwb1k'      : {'s6': 1.0,  's8': 1.061,  'sr6': 1.671 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'tpssh'       : {'s6': 1.0,  's8': 1.219,  'sr6': 1.223 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'scan'        : {'s6': 1.0,  's8': 0.000,  'sr6': 1.324 , 'sr8': 1.000, 'alpha6': 14.0, "citation": '    J.G. Brandenburg, J. E. Bates, J. Sun, J.P. Perdew, Phys. Rev. B 94, 115144, 2016\n'},
        'm11-l'       : {'s6': 1.0,  's8': 1.1129, 'sr6': 2.3933, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'mn12-l'      : {'s6': 1.0,  's8': 0.9622, 'sr6': 2.2329, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'n12'         : {'s6': 1.0,  's8': 2.3916, 'sr6': 1.3493, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'sogga11-x'   : {'s6': 1.0,  's8': 1.8151, 'sr6': 1.5431, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'm11'         : {'s6': 1.0,  's8': 0.6244, 'sr6': 2.2300, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'mn12-sx'     : {'s6': 1.0,  's8': 0.8205, 'sr6': 1.9572, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'n12-sx'      : {'s6': 1.0,  's8': 1.4713, 'sr6': 1.4597, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'hse06'       : {'s6': 1.0,  's8': 0.1090, 'sr6': 1.1290, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    J. Moellmann, S. Grimme, J. Chem. Phys. C 118, 7615-7621, 2014\n'},
        'wb97x'       : {'s6': 1.0,  's8': 1.0,    'sr6': 1.281 , 'sr8': 1.094, 'alpha6': 14.0, "citation": '    Y.-S. Lin, G.-D. Li, S.-P. Mao, J.-D. Chai, Chem. Theory Comput. 9, 263-272, 2013\n'},
        'pbehpbe'     : {'s6': 1.0,  's8': 1.4010, 'sr6': 1.5703, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'rpbe'        : {'s6': 1.0,  's8': 0.5140, 'sr6': 0.8720, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'xlyp'        : {'s6': 1.0,  's8': 0.7447, 'sr6': 0.9384, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pw91p86'     : {'s6': 1.0,  's8': 0.8747, 'sr6': 2.1040, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpwpw91'     : {'s6': 1.0,  's8': 1.9467, 'sr6': 1.3725, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hcth407'     : {'s6': 1.0,  's8': 2.7694, 'sr6': 4.0426, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pkzb'        : {'s6': 1.0,  's8': 0.0000, 'sr6': 0.6327, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpss'     : {'s6': 1.0,  's8': 1.3666, 'sr6': 1.3491, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'thcth'       : {'s6': 1.0,  's8': 0.5662, 'sr6': 0.9320, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mn15-l'      : {'s6': 1.0,  's8': 0.0000, 'sr6': 3.3388, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b3p86'       : {'s6': 1.0,  's8': 1.1961, 'sr6': 1.1897, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b1p96'       : {'s6': 1.0,  's8': 1.1209, 'sr6': 1.1815, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b1lyp'       : {'s6': 1.0,  's8': 1.9467, 'sr6': 1.3725, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw1lyp'     : {'s6': 1.0,  's8': 1.9529, 'sr6': 2.0512, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw1pw91'    : {'s6': 1.0,  's8': 1.4758, 'sr6': 1.2892, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pw1pw'       : {'s6': 1.0,  's8': 1.1786, 'sr6': 1.4968, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw1kcis'    : {'s6': 1.0,  's8': 2.2917, 'sr6': 1.7231, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpwkcis1k'   : {'s6': 1.0,  's8': 1.7553, 'sr6': 1.4853, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pbeh1pbe'    : {'s6': 1.0,  's8': 1.0430, 'sr6': 1.3719, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pbe1kcis'    : {'s6': 1.0,  's8': 1.7934, 'sr6': 3.6355, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'x3lyp'       : {'s6': 1.0,  's8': 0.2990, 'sr6': 1.0000, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'o3lyp'       : {'s6': 1.0,  's8': 1.8058, 'sr6': 1.4060, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b97-1'       : {'s6': 1.0,  's8': 1.6418, 'sr6': 3.7924, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b97-2'       : {'s6': 1.0,  's8': 2.4661, 'sr6': 1.7066, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b98'         : {'s6': 1.0,  's8': 1.9078, 'sr6': 2.6895, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hiss'        : {'s6': 1.0,  's8': 0.7615, 'sr6': 1.3338, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hse03'       : {'s6': 1.0,  's8': 1.0156, 'sr6': 1.3944, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpssh'    : {'s6': 1.0,  's8': 1.2504, 'sr6': 1.3224, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpss0'    : {'s6': 1.0,  's8': 1.0649, 'sr6': 1.2881, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'tpss1kcis'   : {'s6': 1.0,  's8': 2.0902, 'sr6': 1.7729, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'thcthhyb'    : {'s6': 1.0,  's8': 1.6302, 'sr6': 1.5001, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'm08-hx'      : {'s6': 1.0,  's8': 0.0000, 'sr6': 1.6247, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'lcwhpbe'     : {'s6': 1.0,  's8': 1.2797, 'sr6': 1.3846, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw2plyp'    : {'s6': 0.66, 's8': 0.7529, 'sr6': 1.5527, 'sr8': 1.000, 'alpha6': 14.0, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
    # https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/functionals
        'pwb6k'       : {'s6': 1.0,  's8': 0.550,  'sr6': 1.660 , 'sr8': 1.000, 'alpha6': 14.0},
        'revpbe'      : {'s6': 1.0,  's8': 1.010,  'sr6': 0.923 , 'sr8': 1.000, 'alpha6': 14.0},
        'bop'         : {'s6': 1.0,  's8': 1.975,  'sr6': 0.929 , 'sr8': 1.000, 'alpha6': 14.0},
        'otpss'       : {'s6': 1.0,  's8': 1.494,  'sr6': 1.128 , 'sr8': 1.000, 'alpha6': 14.0},
        'pbe38'       : {'s6': 1.0,  's8': 0.998,  'sr6': 1.333 , 'sr8': 1.000, 'alpha6': 14.0},
        'pbesol'      : {'s6': 1.0,  's8': 0.612,  'sr6': 1.345 , 'sr8': 1.000, 'alpha6': 14.0},
        'revssb'      : {'s6': 1.0,  's8': 0.560,  'sr6': 1.221 , 'sr8': 1.000, 'alpha6': 14.0},
        'ssb'         : {'s6': 1.0,  's8': 0.663,  'sr6': 1.215 , 'sr8': 1.000, 'alpha6': 14.0},
        'cam-b3lyp'   : {'s6': 1.0,  's8': 1.217,  'sr6': 1.378 , 'sr8': 1.000, 'alpha6': 14.0},
        'lcwpbe'      : {'s6': 1.0,  's8': 1.279,  'sr6': 1.355 , 'sr8': 1.000, 'alpha6': 14.0},
        'm06-hf'      : {'s6': 1.0,  's8': 0.00 ,  'sr6': 1.446 , 'sr8': 1.000, 'alpha6': 14.0},
        'hcth120'     : {'s6': 1.0,  's8': 1.206,  'sr6': 1.221 , 'sr8': 1.000, 'alpha6': 14.0},
        'ptpss'       : {'s6': 0.75, 's8': 0.879,  'sr6': 1.541 , 'sr8': 1.000, 'alpha6': 14.0},
        'revpbe0'     : {'s6': 1.0,  's8': 0.792,  'sr6': 0.949 , 'sr8': 1.000, 'alpha6': 14.0},
        'revpbe38'    : {'s6': 1.0,  's8': 0.862,  'sr6': 1.021 , 'sr8': 1.000, 'alpha6': 14.0},
    # unreferenced
        'hf'          : {'s6': 1.0,  's8': 1.746,  'sr6': 1.158, 'sr8': 1.0, 'alpha6': 14.0},
    },
    'd3bj': {
    # S. Grimme, S. Ehrlich, L. Goerigk, J. Comput. Chem. 7, 3297-3305, 2011
        'pbe'         : {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407},
        'revpbe'      : {'s6': 1.000, 's8':  2.3550, 'a1':  0.5238, 'a2': 3.5016},
        'blyp'        : {'s6': 1.000, 's8':  2.6996, 'a1':  0.4298, 'a2': 4.2359},
        'bp86'        : {'s6': 1.000, 's8':  3.2822, 'a1':  0.3946, 'a2': 4.8516},
        'rpw86pbe'    : {'s6': 1.000, 's8':  1.3845, 'a1':  0.4613, 'a2': 4.5062},
        'b97-d'       : {'s6': 1.000, 's8':  2.2609, 'a1':  0.5545, 'a2': 3.2297},
        'tpss'        : {'s6': 1.000, 's8':  1.9435, 'a1':  0.4535, 'a2': 4.4752},
        'b3lyp'       : {'s6': 1.000, 's8':  1.9889, 'a1':  0.3981, 'a2': 4.4211},
        'pw6b95'      : {'s6': 1.000, 's8':  0.7257, 'a1':  0.2076, 'a2': 6.3750},
        'pbe0'        : {'s6': 1.000, 's8':  1.2177, 'a1':  0.4145, 'a2': 4.8593},
        'tpss0'       : {'s6': 1.000, 's8':  1.2576, 'a1':  0.3768, 'a2': 4.5865},
        # 'b2plyp'      : {'s6': 0.500, 's8':  1.0860, 'a1':  0.3451, 'a2': 4.7735}, # superseded by a value below.
        'hf'          : {'s6': 1.000, 's8':  0.9171, 'a1':  0.3385, 'a2': 2.8830},
    # Later references
        'bpbe'        : {'s6': 1.000, 's8':  4.0728, 'a1':  0.4567, 'a2': 4.3908, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'opbe'        : {'s6': 1.000, 's8':  3.3816, 'a1':  0.5512, 'a2': 2.9444, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'olyp'        : {'s6': 1.000, 's8':  2.6205, 'a1':  0.5299, 'a2': 2.8065, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpwlyp'      : {'s6': 1.000, 's8':  2.0077, 'a1':  0.4831, 'a2': 4.5323, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b3pw91'      : {'s6': 1.000, 's8':  2.8524, 'a1':  0.4312, 'a2': 4.4693, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'bhlyp'       : {'s6': 1.000, 's8':  1.0354, 'a1':  0.2793, 'a2': 4.9615, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b1b95'       : {'s6': 1.000, 's8':  1.4507, 'a1':  0.2092, 'a2': 5.5545, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpw1b95'     : {'s6': 1.000, 's8':  1.0508, 'a1':  0.1955, 'a2': 6.4177, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'mpwb1k'      : {'s6': 1.000, 's8':  0.9499, 'a1':  0.1474, 'a2': 6.6223, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'tpssh'       : {'s6': 1.000, 's8':  2.2382, 'a1':  0.4529, 'a2': 4.6550, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'bmk'         : {'s6': 1.000, 's8':  2.0860, 'a1':  0.1940, 'a2': 5.9197, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b2plyp'      : {'s6': 0.640, 's8':  0.9147, 'a1':  0.3065, 'a2': 5.0570, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'b2gpplyp'    : {'s6': 0.560, 's8':  0.2597, 'a1':  0.0000, 'a2': 6.3332, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'pwpb95'      : {'s6': 0.820, 's8':  0.2904, 'a1':  0.0000, 'a2': 7.3141, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'dsd-blyp'    : {'s6': 0.500, 's8':  0.2130, 'a1':  0.0000, 'a2': 6.0519, "citation": '    L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670-6688, 2011\n'},
        'hse06'       : {'s6': 1.000, 's8':  2.3100, 'a1':  0.3830, 'a2': 5.6850, "citation": '    J. Moellmann, S. Grimme, J. Chem. Phys. C 118, 7615-7621, 2014\n'},
        'pw91'        : {'s6': 1.000, 's8':  1.9598, 'a1':  0.6319, 'a2': 4.5718, "citation": '    J.R. Reimers et al., Proc. Natl. Acad. Sci. USA 112, E6101-E6110, 2015\n'},
        'm11-l'       : {'s6': 1.000, 's8':  0.4446, 'a1':  0.0000, 'a2': 7.2496, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'mn12-l'      : {'s6': 1.000, 's8':  2.2674, 'a1':  0.0000, 'a2': 9.1494, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'n12'         : {'s6': 1.000, 's8':  4.8491, 'a1':  0.3842, 'a2': 5.3545, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'sogga11-x'   : {'s6': 1.000, 's8':  1.1426, 'a1':  0.1330, 'a2': 5.7381, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'm11'         : {'s6': 1.000, 's8':  2.8112, 'a1':  0.0000, 'a2':10.1389, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'mn12-sx'     : {'s6': 1.000, 's8':  1.1674, 'a1':  0.0983, 'a2': 8.0259, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'n12-sx'      : {'s6': 1.000, 's8':  2.4900, 'a1':  0.3283, 'a2': 5.7898, "citation": '    L. Goerigk, J. Phys. Chem. Lett. 6, 3891-3896, 2015\n'},
        'scan'        : {'s6': 1.000, 's8':  0.0000, 'a1':  0.5380, 'a2': 5.4200, "citation": '    J.G. Brandenburg, J. E. Bates, J. Sun, J.P. Perdew, Phys. Rev. B 94, 115144, 2016\n'},
        'pbehpbe'     : {'s6': 1.000, 's8':  1.1152, 'a1':  0.0000, 'a2': 4.4407, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'rpbe'        : {'s6': 1.000, 's8':  0.8318, 'a1':  0.1820, 'a2': 4.0094, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'xlyp'        : {'s6': 1.000, 's8':  1.5669, 'a1':  0.0809, 'a2': 5.3166, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpwpw91'     : {'s6': 1.000, 's8':  0.3168, 'a1':  0.3168, 'a2': 4.7732, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hcth407'     : {'s6': 1.000, 's8':  0.6490, 'a1':  0.0000, 'a2': 4.8162, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpss'     : {'s6': 1.000, 's8':  1.4023, 'a1':  0.4426, 'a2': 4.4723, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'thcth'       : {'s6': 1.000, 's8':  1.2626, 'a1':  0.0000, 'a2': 5.6162, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b3p86'       : {'s6': 1.000, 's8':  3.3211, 'a1':  0.4601, 'a2': 4.9294, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b1p86'       : {'s6': 1.000, 's8':  3.5681, 'a1':  0.4724, 'a2': 4.9858, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b1lyp'       : {'s6': 1.000, 's8':  2.1167, 'a1':  0.1986, 'a2': 5.3875, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw1pw91'    : {'s6': 1.000, 's8':  1.8744, 'a1':  0.3342, 'a2': 4.9819, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw1kcis'    : {'s6': 1.000, 's8':  1.0893, 'a1':  0.0576, 'a2': 5.5314, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpwkcis1k'   : {'s6': 1.000, 's8':  1.2875, 'a1':  0.0855, 'a2': 5.8961, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pbeh1pbe'    : {'s6': 1.000, 's8':  1.4877, 'a1':  0.0000, 'a2': 7.0385, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'pbe1kcis'    : {'s6': 1.000, 's8':  0.7688, 'a1':  0.0000, 'a2': 6.2794, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'x3lyp'       : {'s6': 1.000, 's8':  1.5744, 'a1':  0.2022, 'a2': 5.4184, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'o3lyp'       : {'s6': 1.000, 's8':  1.8171, 'a1':  0.0963, 'a2': 5.9940, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b97-1'       : {'s6': 1.000, 's8':  0.4814, 'a1':  0.0000, 'a2': 6.2279, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b97-2'       : {'s6': 1.000, 's8':  0.9448, 'a1':  0.0000, 'a2': 5.4603, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'b98'         : {'s6': 1.000, 's8':  0.7086, 'a1':  0.0000, 'a2': 6.0672, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hiss'        : {'s6': 1.000, 's8':  1.6112, 'a1':  0.0000, 'a2': 7.3539, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'hse03'       : {'s6': 1.000, 's8':  1.1243, 'a1':  0.0000, 'a2': 6.8889, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpssh'    : {'s6': 1.000, 's8':  1.4076, 'a1':  0.2660, 'a2': 5.3761, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'revtpss0'    : {'s6': 1.000, 's8':  1.6151, 'a1':  0.2218, 'a2': 5.7985, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'tpss1kcis'   : {'s6': 1.000, 's8':  1.0542, 'a1':  0.0000, 'a2': 6.0201, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'thcthhyb'    : {'s6': 1.000, 's8':  0.9585, 'a1':  0.0000, 'a2': 6.2303, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mn15'        : {'s6': 1.000, 's8':  2.0971, 'a1':  0.7862, 'a2': 7.5923, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'lcwhpbe'     : {'s6': 1.000, 's8':  1.1908, 'a1':  0.2746, 'a2': 5.3157, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
        'mpw2plyp'    : {'s6': 1.000, 's8':  0.6223, 'a1':  0.4105, 'a2': 5.0136, "citation": '    L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, Phys. Chem. Chem. Phys. 19, 32184-32215, 2017\n'},
    # https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/functionalsbj
        'bop'         : {'s6': 1.000, 's8':  3.295,  'a1':  0.4870, 'a2': 3.5043},
        'cam-b3lyp'   : {'s6': 1.000, 's8':  2.0674, 'a1':  0.3708, 'a2': 5.4743},
        'lcwpbe'      : {'s6': 1.000, 's8':  1.8541, 'a1':  0.3919, 'a2': 5.0897},
        'otpss'       : {'s6': 1.000, 's8':  2.7495, 'a1':  0.4634, 'a2': 4.3153},
        'pbe38'       : {'s6': 1.000, 's8':  1.4623, 'a1':  0.3995, 'a2': 5.1405},
        'pbesol'      : {'s6': 1.000, 's8':  2.9491, 'a1':  0.4466, 'a2': 6.1742},
        'ptpss'       : {'s6': 0.750, 's8':  0.2804, 'a1':  0.000,  'a2': 6.5745},
        'pwb6k'       : {'s6': 1.000, 's8':  0.9383, 'a1':  0.1805, 'a2': 7.7627},
        'revssb'      : {'s6': 1.000, 's8':  0.4389, 'a1':  0.4720, 'a2': 4.0986},
        'ssb'         : {'s6': 1.000, 's8': -0.1744, 'a1': -0.0952, 'a2': 5.2170},
        'hcth120'     : {'s6': 1.000, 's8':  1.0821, 'a1':  0.3563, 'a2': 4.3359},
        'revpbe0'     : {'s6': 1.000, 's8':  1.7588, 'a1':  0.4679, 'a2': 3.7619},
        'revpbe38'    : {'s6': 1.000, 's8':  1.4760, 'a1':  0.4309, 'a2': 3.9446},
        # unreferenced
        # special HF/DFT with eBSSE correction
        'hf/mixed'    : {'s6': 1.000, 's8':  3.9027, 'a1':  0.5607, 'a2': 4.5622},
        'hf/sv'       : {'s6': 1.000, 's8':  2.1849, 'a1':  0.4249, 'a2': 4.2783},
        'hf/minis'    : {'s6': 1.000, 's8':  0.9841, 'a1':  0.1702, 'a2': 3.8506},
        'b3lyp/6-31gd': {'s6': 1.000, 's8':  4.0672, 'a1':  0.5014, 'a2': 4.8409},
        # special HF-D3-gCP-SRB/MINIX parametrization
        'hf3c'        : {'s6': 1.000, 's8':  0.8777, 'a1':  0.4171, 'a2': 2.9149},
        # special HF-D3-gCP-SRB2/ECP-2G parametrization
        'hf3cv'       : {'s6': 1.000, 's8':  0.5022, 'a1':  0.3063, 'a2': 3.9856},
        # special PBEh-D3-gCP/def2-mSVP parametrization
        'pbeh3c'      : {'s6': 1.000, 's8':  0.0000, 'a1':  0.4860, 'a2': 4.5000},
        'core-dsd-blyp'    : {'s6': 0.500, 's8':  0.2130, 'a1':  0.0000, 'a2': 6.0519},
    },
    'd3mzero': {
        'b2plyp'      : {'s6': 0.640, 's8':  0.717543, 'sr6': 1.313134, 'beta': 0.016035},
        'b3lyp'       : {'s6': 1.000, 's8':  1.532981, 'sr6': 1.338153, 'beta': 0.013988},
        'b97-d'       : {'s6': 1.000, 's8':  1.020078, 'sr6': 1.151808, 'beta': 0.035964},
        'blyp'        : {'s6': 1.000, 's8':  1.841686, 'sr6': 1.279637, 'beta': 0.014370},
        'bp86'        : {'s6': 1.000, 's8':  1.945174, 'sr6': 1.233460, 'beta': 0.000000},
        'pbe'         : {'s6': 1.000, 's8':  0.000000, 'sr6': 2.340218, 'beta': 0.129434},
        'pbe0'        : {'s6': 1.000, 's8':  0.000081, 'sr6': 2.077949, 'beta': 0.116755},
        'lcwpbe'      : {'s6': 1.000, 's8':  1.280619, 'sr6': 1.366361, 'beta': 0.003160},
    },
    'd3mbj': {
        'b2plyp'      : {'s6': 0.640, 's8': 0.672820, 'a1': 0.486434, 'a2': 3.656466},
        'b3lyp'       : {'s6': 1.000, 's8': 1.466677, 'a1': 0.278672, 'a2': 4.606311},
        'b97-d'       : {'s6': 1.000, 's8': 1.206988, 'a1': 0.240184, 'a2': 3.864426},
        'blyp'        : {'s6': 1.000, 's8': 1.875007, 'a1': 0.448486, 'a2': 3.610679},
        'bp86'        : {'s6': 1.000, 's8': 3.140281, 'a1': 0.821850, 'a2': 2.728151},
        'pbe'         : {'s6': 1.000, 's8': 0.358940, 'a1': 0.012092, 'a2': 5.938951},
        'pbe0'        : {'s6': 1.000, 's8': 0.528823, 'a1': 0.007912, 'a2': 6.162326},
        'lcwpbe'      : {'s6': 1.000, 's8': 0.906564, 'a1': 0.563761, 'a2': 3.593680},
    },
    'nl': { 
        'blyp'        : {'b':  4.000, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'}, 
        'hf'          : {'b':  3.900, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'}, # not implemented
        'revpbe'      : {'b':  3.700, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'}, 
        'revpbe38'    : {'b':  4.700, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'}, 
        'b3lyp'       : {'b':  4.800, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'},
        'b3pw91'      : {'b':  4.500, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'},
        'revpbe0'     : {'b':  4.300, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme J. Chem. Theory Comput. 7, 3866-3871, 2011 \n'},
        'bp86'        : {'b':  4.400, 'c': 0.0093, 'citation' : '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016 \n'}, # unclear if this is the real origin
        'pbe0'        : {'b':  6.900, 'c': 0.0093, 'citation' : '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016 \n'}, # unclear if this is the real origin
        'pbe'         : {'b':  6.400, 'c': 0.0093, 'citation' : '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016 \n'}, # unclear if this is the real origin
        'tpss0'       : {'b':  5.500, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme, J. Chem. Theory Comput. 9, 308-315, 2013 \n'},
        'tpss'        : {'b':  5.000, 'c': 0.0093, 'citation' : '    W. Hujo, S. Grimme, J. Chem. Theory Comput. 9, 308-315, 2013 \n'},
        'b2gpplyp'    : {'b':  9.900, 'c': 0.0093, 'citation' : '    M. K. Kesharwani, A. Karton, J.M. L. Martin, J. Chem. Theory Comput. 12, 444-454, 2016 \n'},
        'b2plyp'      : {'b':  7.800, 'c': 0.0093, 'citation' : '    J. Calbo, E. Orti, J. C. Sancho-Garcia, J. Arago, J. Chem. Theory Comput. 11, 932-939, 1015 \n'}, 
        'pwpb95'      : {'b': 11.100, 'c': 0.0093, 'citation' : '    F. Yu J. Chem. Theory Comput. 10, 4400-4407, 2014 \n'}, 
    },

} # yapf: disable


# Full list of all possible endings
# full_dash_keys = list(dashcoeff) + list(dash_alias)


def dash_server(func, dashlvl):
    """ Returns the dictionary of keys for default empirical parameters"""
    # Validate input arguments
    dashlvl = dashlvl.lower()
    dashlvleff = get_dispersion_aliases()[dashlvl] if dashlvl in get_dispersion_aliases() else dashlvl
    
    if dashlvleff not in dashcoeff.keys():
        raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl,
                                                                                                 dashcoeff.keys()))

    func = func.lower()
    if func not in dashcoeff[dashlvleff].keys():
        raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))

    # Return the values
    return dashcoeff[dashlvleff][func]


def dftd3_coeff_formatter(dashlvl, dashcoeff):
    # Return strings for dftd3 program parameter file
    #            s6 rs6 s18 rs8 alpha6 version
    #   d2p4:    s6 sr6=1.1 s8=0.0 a2=None alpha6=20.0 version=2
    #   d2gr:    s6 sr6=1.1 s8=0.0 a2=None alpha6 version=2
    #   d3zero:  s6 sr6 s8 a2=None alpha6 version=3
    #   d3bj:    s6 a1 s8 a2 alpha6=None version=4
    #   d3mzero: s6 sr6 s8 beta alpha6=14.0 version=5
    #   d3mbj:   s6 a1 s8 a2 alpha6=None version=6

    dashlvleff = get_dispersion_aliases()[dashlvl] if dashlvl in get_dispersion_aliases() else dashlvl
    dashformatter = """{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:6}\n"""

    if dashlvleff.lower() == 'd2p4':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             1.1, 0.0, 0.0, 20.0, 2)
    elif dashlvleff.lower() == 'd2gr':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             dashcoeff['sr6'],
             0.0, 0.0,
             dashcoeff['alpha6'],
             2)
    elif dashlvleff.lower() == 'd3zero':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             dashcoeff['sr6'],
             dashcoeff['s8'],
             dashcoeff['sr8'],
             dashcoeff['alpha6'],
             3)
    elif dashlvleff.lower() == 'd3bj':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             dashcoeff['a1'],
             dashcoeff['s8'],
             dashcoeff['a2'],
             0.0, 4)
    elif dashlvleff.lower() == 'd3mzero':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             dashcoeff['sr6'],
             dashcoeff['s8'],
             dashcoeff['beta'],
             14.0, 5)
    elif dashlvleff.lower() == 'd3mbj':
        returnstring = dashformatter.format(
             dashcoeff['s6'],
             dashcoeff['a1'],
             dashcoeff['s8'],
             dashcoeff['a2'],
             0.0, 6)
    else:
        raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl,
                                                                                                 dashcoeff.keys()))

    return returnstring


def get_default_dashparams(dtype):
    # returns a dictionary containing default dispersion parameters for a given method
    if dtype == 'd2p4':
        return ({"s6": 1.0})
    elif dtype == 'd2gr':
        return ({"s6": 1.0, "sr6": 1.1, "s8": 0.0, "alpha6": 20.0})
    elif dtype == 'd3zero':
        return ({"s6": 1.0, "sr6": 1.0, "s8": 0.0, "sr8": 1.0, "alpha6": 14.0})
    elif dtype == 'd3bj':
        return ({"s6": 1.0, "a1": 0.0, "s8": 1.0, "a2": 1.0})
    elif dtype == 'd3mzero':
        return ({"s6": 1.0, "sr6": 1.0, "s8": 1.0, "beta": 1.0, "alpha6": 14.0})
    elif dtype == 'd3mbj':
        return ({"s6": 1.0, "a1": 1.0, "s8": 1.0, "a2": 1.0})
    elif dtype == 'nl':
        return ({"b": 1.0, "c": 0.0093})
    else:
        return ({"s6": 1.0})
