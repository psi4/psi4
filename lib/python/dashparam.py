"""
Module to hold and distribute the -D dispersion correction parameters.
"""
from psiexceptions import *


## ==> Dispersion Aliases and Parameters <== ##

# This defines the -D aliases for all of psi4
dash_alias = {
    '-d': '-d2p4',     # means -D aliases to a -D2 level dispersion correction, as opposed to -D3
    '-d2': '-d2p4',    # means -D2 uses psi4's internal -D2 correction, as opposed to calling dftd3
    '-d3': '-d3zero'}  # means -D3 uses the original zero-damping fn, as opposed to bj-damping


# The dashcoeff dict below defines the -D parameters for all of psi4. 'd2p4' are
#   taken from already defined functionals in psi4. The remainder of the parameters are
#   from http://toc.uni-muenster.de/DFTD3/ on September 25, 2012, with the dict keys
#   translated from Turbomole to Psi4 functional names.
dashcoeff = {
    'd2p4': {
        'b97-d'       : {'s6': 1.25},  #IN
        'blyp'        : {'s6': 1.20},  #IN
        'b3lyp'       : {'s6': 1.05},  #IN
        'bp86'        : {'s6': 1.05},  #IN
        'pbe'         : {'s6': 0.75},  #IN
        'pbe0'        : {'s6': 0.60},  #IN
        'dsd-blyp'    : {'s6': 0.35},  # different btwn dftd3 and psi4
        'dsd-pbep86'  : {'s6': 0.29},  #IN
        'dsd-pbepbe'  : {'s6': 0.42},  #IN
        'b2plyp'      : {'s6': 0.55},  #IN
    },
    'd2gr': {
        'blyp'        : {'s6': 1.2,  'alpha6': 20.0},  # in psi4  #IN
        'bp86'        : {'s6': 1.05, 'alpha6': 20.0},  # in psi4  #IN
        'b97-d'       : {'s6': 1.25, 'alpha6': 20.0},  # in psi4  #IN
        'revpbe'      : {'s6': 1.25, 'alpha6': 20.0},
        'pbe'         : {'s6': 0.75, 'alpha6': 20.0},  # in psi4  #IN
        'tpss'        : {'s6': 1.0,  'alpha6': 20.0},
        'b3lyp'       : {'s6': 1.05, 'alpha6': 20.0},  # in psi4  #IN
        'pbe0'        : {'s6': 0.6,  'alpha6': 20.0},  # in psi4  #IN
        'pw6b95'      : {'s6': 0.5,  'alpha6': 20.0},
        'tpss0'       : {'s6': 0.85, 'alpha6': 20.0},
        'b2plyp'      : {'s6': 0.55, 'alpha6': 20.0},  # in psi4  #IN
        'b2gp-plyp'   : {'s6': 0.4,  'alpha6': 20.0},
        'dsd-blyp'    : {'s6': 0.41, 'alpha6': 60.0},  # in psi4
    },
    'd3zero': {
        'b1b95'       : {'s6': 1.0,  's8': 1.868, 'sr6': 1.613, 'alpha6': 14.0},
        'b2gpplyp'    : {'s6': 0.56, 's8': 0.760, 'sr6': 1.586, 'alpha6': 14.0},
        'b3lyp'       : {'s6': 1.0,  's8': 1.703, 'sr6': 1.261, 'alpha6': 14.0},  # in psi4  #IN
        'b97-d'       : {'s6': 1.0,  's8': 0.909, 'sr6': 0.892, 'alpha6': 14.0},  # in psi4  #IN
        'bhlyp'       : {'s6': 1.0,  's8': 1.442, 'sr6': 1.370, 'alpha6': 14.0},
        'blyp'        : {'s6': 1.0,  's8': 1.682, 'sr6': 1.094, 'alpha6': 14.0},  # in psi4  #IN
        'bp86'        : {'s6': 1.0,  's8': 1.683, 'sr6': 1.139, 'alpha6': 14.0},  # in psi4  #IN
        'bpbe'        : {'s6': 1.0,  's8': 2.033, 'sr6': 1.087, 'alpha6': 14.0},
        'mpwlyp'      : {'s6': 1.0,  's8': 1.098, 'sr6': 1.239, 'alpha6': 14.0},
        'pbe'         : {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0},  # in psi4  #IN
        'pbe0'        : {'s6': 1.0,  's8': 0.928, 'sr6': 1.287, 'alpha6': 14.0},  # in psi4  #IN
        'pw6b95'      : {'s6': 1.0,  's8': 0.862, 'sr6': 1.532, 'alpha6': 14.0},
        'pwb6k'       : {'s6': 1.0,  's8': 0.550, 'sr6': 1.660, 'alpha6': 14.0},
        'revpbe'      : {'s6': 1.0,  's8': 1.010, 'sr6': 0.923, 'alpha6': 14.0},
        'tpss'        : {'s6': 1.0,  's8': 1.105, 'sr6': 1.166, 'alpha6': 14.0},
        'tpss0'       : {'s6': 1.0,  's8': 1.242, 'sr6': 1.252, 'alpha6': 14.0},
        'tpssh'       : {'s6': 1.0,  's8': 1.219, 'sr6': 1.223, 'alpha6': 14.0},
        'bop'         : {'s6': 1.0,  's8': 1.975, 'sr6': 0.929, 'alpha6': 14.0},
        'mpw1b95'     : {'s6': 1.0,  's8': 1.118, 'sr6': 1.605, 'alpha6': 14.0},
        'mpwb1k'      : {'s6': 1.0,  's8': 1.061, 'sr6': 1.671, 'alpha6': 14.0},
        'olyp'        : {'s6': 1.0,  's8': 1.764, 'sr6': 0.806, 'alpha6': 14.0},
        'opbe'        : {'s6': 1.0,  's8': 2.055, 'sr6': 0.837, 'alpha6': 14.0},
        'otpss'       : {'s6': 1.0,  's8': 1.494, 'sr6': 1.128, 'alpha6': 14.0},
        'pbe38'       : {'s6': 1.0,  's8': 0.998, 'sr6': 1.333, 'alpha6': 14.0},
        'pbesol'      : {'s6': 1.0,  's8': 0.612, 'sr6': 1.345, 'alpha6': 14.0},
        'revssb'      : {'s6': 1.0,  's8': 0.560, 'sr6': 1.221, 'alpha6': 14.0},
        'ssb'         : {'s6': 1.0,  's8': 0.663, 'sr6': 1.215, 'alpha6': 14.0},
        'b3pw91'      : {'s6': 1.0,  's8': 1.775, 'sr6': 1.176, 'alpha6': 14.0},
        'bmk'         : {'s6': 1.0,  's8': 2.168, 'sr6': 1.931, 'alpha6': 14.0},
        'camb3lyp'    : {'s6': 1.0,  's8': 1.217, 'sr6': 1.378, 'alpha6': 14.0},
        'lcwpbe'      : {'s6': 1.0,  's8': 1.279, 'sr6': 1.355, 'alpha6': 14.0},
        'm05-2x'      : {'s6': 1.0,  's8': 0.00 , 'sr6': 1.417, 'alpha6': 14.0},  # in psi4  #IN
        'm05'         : {'s6': 1.0,  's8': 0.595, 'sr6': 1.373, 'alpha6': 14.0},  # in psi4  #IN
        'm062x'       : {'s6': 1.0,  's8': 0.00 , 'sr6': 1.619, 'alpha6': 14.0},
        'm06hf'       : {'s6': 1.0,  's8': 0.00 , 'sr6': 1.446, 'alpha6': 14.0},
        'm06l'        : {'s6': 1.0,  's8': 0.00 , 'sr6': 1.581, 'alpha6': 14.0},
        'm06'         : {'s6': 1.0,  's8': 0.00 , 'sr6': 1.325, 'alpha6': 14.0},
        'hcth120'     : {'s6': 1.0,  's8': 1.206, 'sr6': 1.221, 'alpha6': 14.0},  # in psi4  #IN
        'b2plyp'      : {'s6': 0.64, 's8': 1.022, 'sr6': 1.427, 'alpha6': 14.0},  # in psi4  #IN
        'dsd-blyp'    : {'s6': 0.50, 's8': 0.705, 'sr6': 1.569, 'alpha6': 14.0},  # in psi4
        'ptpss'       : {'s6': 0.75, 's8': 0.879, 'sr6': 1.541, 'alpha6': 14.0},
        'pwpb95'      : {'s6': 0.82, 's8': 0.705, 'sr6': 1.557, 'alpha6': 14.0},
        'revpbe0'     : {'s6': 1.0,  's8': 0.792, 'sr6': 0.949, 'alpha6': 14.0},
        'revpbe38'    : {'s6': 1.0,  's8': 0.862, 'sr6': 1.021, 'alpha6': 14.0},
        'rpw86pbe'    : {'s6': 1.0,  's8': 0.901, 'sr6': 1.224, 'alpha6': 14.0},
    },
    'd3bj': {
        'b1b95'       : {'s6': 1.000, 's8':  1.4507, 'a1':  0.2092, 'a2': 5.5545},
        'b2gpplyp'    : {'s6': 0.560, 's8':  0.2597, 'a1':  0.0000, 'a2': 6.3332},
        'b3pw91'      : {'s6': 1.000, 's8':  2.8524, 'a1':  0.4312, 'a2': 4.4693},
        'bhlyp'       : {'s6': 1.000, 's8':  1.0354, 'a1':  0.2793, 'a2': 4.9615},
        'bmk'         : {'s6': 1.000, 's8':  2.0860, 'a1':  0.1940, 'a2': 5.9197},
        'bop'         : {'s6': 1.000, 's8':  3.295,  'a1':  0.4870, 'a2': 3.5043},
        'bpbe'        : {'s6': 1.000, 's8':  4.0728, 'a1':  0.4567, 'a2': 4.3908},
        'camb3lyp'    : {'s6': 1.000, 's8':  2.0674, 'a1':  0.3708, 'a2': 5.4743},
        'lcwpbe'      : {'s6': 1.000, 's8':  1.8541, 'a1':  0.3919, 'a2': 5.0897},
        'mpw1b95'     : {'s6': 1.000, 's8':  1.0508, 'a1':  0.1955, 'a2': 6.4177},
        'mpwb1k'      : {'s6': 1.000, 's8':  0.9499, 'a1':  0.1474, 'a2': 6.6223},
        'mpwlyp'      : {'s6': 1.000, 's8':  2.0077, 'a1':  0.4831, 'a2': 4.5323},
        'olyp'        : {'s6': 1.000, 's8':  2.6205, 'a1':  0.5299, 'a2': 2.8065},
        'opbe'        : {'s6': 1.000, 's8':  3.3816, 'a1':  0.5512, 'a2': 2.9444},
        'otpss'       : {'s6': 1.000, 's8':  2.7495, 'a1':  0.4634, 'a2': 4.3153},
        'pbe38'       : {'s6': 1.000, 's8':  1.4623, 'a1':  0.3995, 'a2': 5.1405},
        'pbesol'      : {'s6': 1.000, 's8':  2.9491, 'a1':  0.4466, 'a2': 6.1742},
        'ptpss'       : {'s6': 0.750, 's8':  0.2804, 'a1':  0.000,  'a2': 6.5745},
        'pwb6k'       : {'s6': 1.000, 's8':  0.9383, 'a1':  0.1805, 'a2': 7.7627},
        'revssb'      : {'s6': 1.000, 's8':  0.4389, 'a1':  0.4720, 'a2': 4.0986},
        'ssb'         : {'s6': 1.000, 's8': -0.1744, 'a1': -0.0952, 'a2': 5.2170},
        'tpssh'       : {'s6': 1.000, 's8':  0.4243, 'a1':  0.0000, 'a2': 5.5253},
        'hcth120'     : {'s6': 1.000, 's8':  1.0821, 'a1':  0.3563, 'a2': 4.3359},  # in psi4  #IN
        'b2plyp'      : {'s6': 0.640, 's8':  0.9147, 'a1':  0.3065, 'a2': 5.0570},  # in psi4  #IN
        'b3lyp'       : {'s6': 1.000, 's8':  1.9889, 'a1':  0.3981, 'a2': 4.4211},  # in psi4  #IN
        'b97-d'       : {'s6': 1.000, 's8':  2.2609, 'a1':  0.5545, 'a2': 3.2297},  # in psi4  #IN
        'blyp'        : {'s6': 1.000, 's8':  2.6996, 'a1':  0.4298, 'a2': 4.2359},  # in psi4  #IN
        'bp86'        : {'s6': 1.000, 's8':  3.2822, 'a1':  0.3946, 'a2': 4.8516},  # in psi4  #IN
        'dsd-blyp'    : {'s6': 0.500, 's8':  0.2130, 'a1':  0.000,  'a2': 6.0519},  # in psi4
        'pbe0'        : {'s6': 1.000, 's8':  1.2177, 'a1':  0.4145, 'a2': 4.8593},  # in psi4  #IN
        'pbe'         : {'s6': 1.000, 's8':  0.7875, 'a1':  0.4289, 'a2': 4.4407},  # in psi4  #IN
        'pw6b95'      : {'s6': 1.000, 's8':  0.7257, 'a1':  0.2076, 'a2': 6.3750},
        'pwpb95'      : {'s6': 0.820, 's8':  0.2904, 'a1':  0.0000, 'a2': 7.3141},
        'revpbe0'     : {'s6': 1.000, 's8':  1.7588, 'a1':  0.4679, 'a2': 3.7619},
        'revpbe38'    : {'s6': 1.000, 's8':  1.4760, 'a1':  0.4309, 'a2': 3.9446},
        'revpbe'      : {'s6': 1.000, 's8':  2.3550, 'a1':  0.5238, 'a2': 3.5016},
        'rpw86pbe'    : {'s6': 1.000, 's8':  1.3845, 'a1':  0.4613, 'a2': 4.5062},
        'tpss0'       : {'s6': 1.000, 's8':  1.2576, 'a1':  0.3768, 'a2': 4.5865},
        'tpss'        : {'s6': 1.000, 's8':  1.9435, 'a1':  0.4535, 'a2': 4.4752},
    }
}


def dash_server(func, dashlvl, mode='psi4'):
    """Function to serve up dispersion correction parameters in whatever form needed.
    When *mode* is 'dftd3', returns a string suitable for writing to ./dftd3_parameters
    to calculuate the correction at *dashlvl* with the default parameters for functional
    *func*. When *mode* is 'psi4', returns a tuple of arguments suitable for building
    a Dispersion object with *dashlvl* parameters for functional *func*.

    There are four computational *dashlvl* choices. 'd2p4' calls the -D2 correction
    within psi4 (hence, faked for mode='dftd3'). The other three, 'd2gr', 'd3zero',
    and 'd3bj' call the three dftd3 modes of operation (corresponding to -old, -zero, -bj).
    Additionally, there are three aliased *dashlvl* choices since the aliases in dash_alias
    above are imposed.

    """
    # Validate input arguments
    dashlvl = dashlvl.lower()
    dashlvleff = dash_alias['-' + dashlvl][1:] if ('-' + dashlvl) in dash_alias.keys() else dashlvl

    func = func.lower()
    if func not in dashcoeff[dashlvleff].keys():
        raise ValidationError("""Functional %s is not available for -D level %s.""" % (func, dashlvl))

    # Return strings for dftd3 program parameter file
    #   d2p4:    s6 sr6=1.1 s8=0.0 a2=None alpha6=20.0 version=2
    #   d2gr:    s6 sr6=1.1 s8=0.0 a2=None alpha6 version=2
    #   d3zero:  s6 sr6 s8 a2=None alpha6 version=3
    #   d3bj:    s6 a1 s8 a2 alpha6=None version=4
    if mode.lower() == 'dftd3':
        if dashlvleff.lower() == 'd2p4':
            returnstring = '%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' % \
                (dashcoeff[dashlvleff][func]['s6'],
                 1.1, 0.0, 0.0, 20.0, 2)
        elif dashlvleff.lower() == 'd2gr':
            returnstring = '%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' % \
                (dashcoeff[dashlvleff][func]['s6'],
                 1.1, 0.0, 0.0,
                 dashcoeff[dashlvleff][func]['alpha6'],
                 2)
        elif dashlvleff.lower() == 'd3zero':
            returnstring = '%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' % \
                (dashcoeff[dashlvleff][func]['s6'],
                 dashcoeff[dashlvleff][func]['sr6'],
                 dashcoeff[dashlvleff][func]['s8'],
                 1.0,
                 dashcoeff[dashlvleff][func]['alpha6'],
                 3)
        elif dashlvleff.lower() == 'd3bj':
            returnstring = '%12.6f %12.6f %12.6f %12.6f %12.6f %6d\n' % \
                (dashcoeff[dashlvleff][func]['s6'],
                 dashcoeff[dashlvleff][func]['a1'],
                 dashcoeff[dashlvleff][func]['s8'],
                 dashcoeff[dashlvleff][func]['a2'],
                 0.0, 4)
        else:
            raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))
        return returnstring

    # Return parameter list for Dispersion.build function
    #   d2p4:    name=-D2     s6=s6 p1=None p2=None p3=None (damping parameter fixed in Dispersion class)
    #   d2gr:    name=-D2GR   s6=s6 p1=alpha6 p2=None p3=None
    #   d3zero:  name=-D3ZERO s6=s6 p1=sr6 p2=s8 p3=alpha6
    #   d3bj:    name=-D3BJ   s6=s6 p1=a1 p2=s8 p3=a2
    elif mode.lower() == 'psi4':
        if dashlvleff == 'd2p4':
            return '-D2', \
                dashcoeff[dashlvleff][func]['s6'], \
                0.0, 0.0, 0.0
        elif dashlvleff == 'd2gr':
            return '-D2GR', \
                dashcoeff[dashlvleff][func]['s6'], \
                dashcoeff[dashlvleff][func]['alpha6'], \
                0.0, 0.0
        elif dashlvleff == 'd3zero':
            return '-D3ZERO', \
                dashcoeff[dashlvleff][func]['s6'], \
                dashcoeff[dashlvleff][func]['s8'], \
                dashcoeff[dashlvleff][func]['sr6'], \
                dashcoeff[dashlvleff][func]['alpha6']
        elif dashlvleff == 'd3bj':
            return '-D3BJ', \
                dashcoeff[dashlvleff][func]['s6'], \
                dashcoeff[dashlvleff][func]['s8'], \
                dashcoeff[dashlvleff][func]['a1'], \
                dashcoeff[dashlvleff][func]['a2']
        else:
            raise ValidationError("""-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))

    else:
        raise ValidationError("""-D return format %s is not available. Choose 'psi4' or 'dftd3'.""")
