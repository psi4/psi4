#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

"""
| Database (Hobza) of interaction energies for bimolecular complexes.
| Geometries from Jurecka et al. PCCP 8 1985 (2006).
| First revision to interaction energies (S22A) from Takatani et al. JCP 132 144104 (2010).
| Second revision to interaction energies (S22B) from Marshall et al. JCP 135 194102 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **benchmark**

  - ``'S220'`` Jurecka et al. PCCP 8 1985 (2006).
  - ``'S22A'`` Takatani et al. JCP 132 144104 (2010).
  - |dl| ``'S22B'`` |dr| Marshall et al. JCP 135 194102 (2011).

- **subset**

  - ``'small'`` water dimer, methane dimer, ethene-ethine
  - ``'large'`` adenine-thymine
  - ``'HB'`` hydrogen-bonded systems
  - ``'MX'`` mixed-influence systems
  - ``'DD'`` dispersion-dominated systems
  - ``'S11'`` smaller systems in S22
  - ``'WATER'`` water dimer

"""
import qcdb

# <<< S22 Database Module >>>
dbse = 'S22'

# <<< Database Members >>>
HRXN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
HRXN_SM = [2, 8, 16]
HRXN_LG = [15]
HB = [1, 2, 3, 4, 5, 6, 7]
MX = [13, 15, 16, 17, 18, 19, 21, 22]
DD = [8, 9, 10, 11, 12, 14, 20]
S11 = [1, 2, 3, 4, 8, 9, 10, 16, 17, 18, 19]
WATER = [2]

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations
for rxn in HRXN:

    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }

    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-CP'   % (dbse, rxn),
                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]

    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-unCP' % (dbse, rxn),
                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values >>>
BIND = {}
# Original publication
BIND_S220 = {}
BIND_S220['%s-%s' % (dbse,  1)] =  -3.17
BIND_S220['%s-%s' % (dbse,  2)] =  -5.02
BIND_S220['%s-%s' % (dbse,  3)] = -18.61
BIND_S220['%s-%s' % (dbse,  4)] = -15.96
BIND_S220['%s-%s' % (dbse,  5)] = -20.65
BIND_S220['%s-%s' % (dbse,  6)] = -16.71
BIND_S220['%s-%s' % (dbse,  7)] = -16.37
BIND_S220['%s-%s' % (dbse,  8)] =  -0.53
BIND_S220['%s-%s' % (dbse,  9)] =  -1.51
BIND_S220['%s-%s' % (dbse, 10)] =  -1.50
BIND_S220['%s-%s' % (dbse, 11)] =  -2.73
BIND_S220['%s-%s' % (dbse, 12)] =  -4.42
BIND_S220['%s-%s' % (dbse, 13)] = -10.12
BIND_S220['%s-%s' % (dbse, 14)] =  -5.22
BIND_S220['%s-%s' % (dbse, 15)] = -12.23
BIND_S220['%s-%s' % (dbse, 16)] =  -1.53
BIND_S220['%s-%s' % (dbse, 17)] =  -3.28
BIND_S220['%s-%s' % (dbse, 18)] =  -2.35
BIND_S220['%s-%s' % (dbse, 19)] =  -4.46
BIND_S220['%s-%s' % (dbse, 20)] =  -2.74
BIND_S220['%s-%s' % (dbse, 21)] =  -5.73
BIND_S220['%s-%s' % (dbse, 22)] =  -7.05
# Revision
BIND_S22A = {}
BIND_S22A['%s-%s' % (dbse,  1)] =  -3.15
BIND_S22A['%s-%s' % (dbse,  2)] =  -5.07
BIND_S22A['%s-%s' % (dbse,  3)] = -18.81
BIND_S22A['%s-%s' % (dbse,  4)] = -16.11
BIND_S22A['%s-%s' % (dbse,  5)] = -20.69
BIND_S22A['%s-%s' % (dbse,  6)] = -17.00
BIND_S22A['%s-%s' % (dbse,  7)] = -16.74
BIND_S22A['%s-%s' % (dbse,  8)] =  -0.53
BIND_S22A['%s-%s' % (dbse,  9)] =  -1.48
BIND_S22A['%s-%s' % (dbse, 10)] =  -1.45
BIND_S22A['%s-%s' % (dbse, 11)] =  -2.62
BIND_S22A['%s-%s' % (dbse, 12)] =  -4.20
BIND_S22A['%s-%s' % (dbse, 13)] =  -9.74
BIND_S22A['%s-%s' % (dbse, 14)] =  -4.59
BIND_S22A['%s-%s' % (dbse, 15)] = -11.66
BIND_S22A['%s-%s' % (dbse, 16)] =  -1.50
BIND_S22A['%s-%s' % (dbse, 17)] =  -3.29
BIND_S22A['%s-%s' % (dbse, 18)] =  -2.32
BIND_S22A['%s-%s' % (dbse, 19)] =  -4.55
BIND_S22A['%s-%s' % (dbse, 20)] =  -2.71
BIND_S22A['%s-%s' % (dbse, 21)] =  -5.62
BIND_S22A['%s-%s' % (dbse, 22)] =  -7.09
# Current revision
BIND_S22B = {}
BIND_S22B['%s-%s' % (dbse,  1)] =  -3.133
BIND_S22B['%s-%s' % (dbse,  2)] =  -4.989
BIND_S22B['%s-%s' % (dbse,  3)] = -18.753
BIND_S22B['%s-%s' % (dbse,  4)] = -16.062
BIND_S22B['%s-%s' % (dbse,  5)] = -20.641
BIND_S22B['%s-%s' % (dbse,  6)] = -16.934
BIND_S22B['%s-%s' % (dbse,  7)] = -16.660
BIND_S22B['%s-%s' % (dbse,  8)] =  -0.527
BIND_S22B['%s-%s' % (dbse,  9)] =  -1.472
BIND_S22B['%s-%s' % (dbse, 10)] =  -1.448
BIND_S22B['%s-%s' % (dbse, 11)] =  -2.654
BIND_S22B['%s-%s' % (dbse, 12)] =  -4.255
BIND_S22B['%s-%s' % (dbse, 13)] =  -9.805
BIND_S22B['%s-%s' % (dbse, 14)] =  -4.524
BIND_S22B['%s-%s' % (dbse, 15)] = -11.730
BIND_S22B['%s-%s' % (dbse, 16)] =  -1.496
BIND_S22B['%s-%s' % (dbse, 17)] =  -3.275
BIND_S22B['%s-%s' % (dbse, 18)] =  -2.312
BIND_S22B['%s-%s' % (dbse, 19)] =  -4.541
BIND_S22B['%s-%s' % (dbse, 20)] =  -2.717
BIND_S22B['%s-%s' % (dbse, 21)] =  -5.627
BIND_S22B['%s-%s' % (dbse, 22)] =  -7.097
# Set default
BIND = BIND_S22B

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse,  1)] = 'HB-1 Ammonia Dimer, C2H'
TAGL['%s-%s-dimer'      % (dbse,  1)] =      'Ammonia Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  1)] =      'Ammonia from Ammonia Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  1)] =      'Ammonia from Ammonia Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  1)] =      'Ammonia from Ammonia Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  1)] =      'Ammonia from Ammonia Dimer'
TAGL['%s-%s'            % (dbse,  2)] = 'HB-2 Water Dimer, CS'
TAGL['%s-%s-dimer'      % (dbse,  2)] =      'Water Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  2)] =      'Water from Water Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  2)] =      'Water from Water Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  2)] =      'Water from Water Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  2)] =      'Water from Water Dimer'
TAGL['%s-%s'            % (dbse,  3)] = 'HB-3 Formic Acid Dimer, C2H'
TAGL['%s-%s-dimer'      % (dbse,  3)] =      'Formic Acid Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  3)] =      'Formic Acid from Formic Acid Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  3)] =      'Formic Acid from Formic Acid Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  3)] =      'Formic Acid from Formic Acid Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  3)] =      'Formic Acid from Formic Acid Dimer'
TAGL['%s-%s'            % (dbse,  4)] = 'HB-4 Formamide Dimer, C2H'
TAGL['%s-%s-dimer'      % (dbse,  4)] =      'Formamide Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  4)] =      'Formamide from Formamide Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  4)] =      'Formamide from Formamide Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  4)] =      'Formamide from Formamide Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  4)] =      'Formamide from Formamide Dimer'
TAGL['%s-%s'            % (dbse,  5)] = 'HB-5 Uracil Dimer HB, C2H'
TAGL['%s-%s-dimer'      % (dbse,  5)] =      'Uracil Dimer HB'
TAGL['%s-%s-monoA-CP'   % (dbse,  5)] =      'Uracil from Uracil Dimer HB'
TAGL['%s-%s-monoB-CP'   % (dbse,  5)] =      'Uracil from Uracil Dimer HB'
TAGL['%s-%s-monoA-unCP' % (dbse,  5)] =      'Uracil from Uracil Dimer HB'
TAGL['%s-%s-monoB-unCP' % (dbse,  5)] =      'Uracil from Uracil Dimer HB'
TAGL['%s-%s'            % (dbse,  6)] = 'HB-6 2-Pyridone-2-Aminopyridine Complex, C1'
TAGL['%s-%s-dimer'      % (dbse,  6)] =      '2-Pyridone-2-Aminopyridine Complex'
TAGL['%s-%s-monoA-CP'   % (dbse,  6)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex'
TAGL['%s-%s-monoB-CP'   % (dbse,  6)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex'
TAGL['%s-%s-monoA-unCP' % (dbse,  6)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex'
TAGL['%s-%s-monoB-unCP' % (dbse,  6)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex'
TAGL['%s-%s'            % (dbse,  7)] = 'HB-7 Adenine-Thymine Complex WC, C1'
TAGL['%s-%s-dimer'      % (dbse,  7)] =      'Adenine-Thymine Complex WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  7)] =      'Adenine from Adenine-Thymine Complex WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  7)] =      'Thymine from Adenine-Thymine Complex WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  7)] =      'Adenine from Adenine-Thymine Complex WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  7)] =      'Thymine from Adenine-Thymine Complex WC'
TAGL['%s-%s'            % (dbse,  8)] = 'DD-1 Methane Dimer, D3D'
TAGL['%s-%s-dimer'      % (dbse,  8)] =      'Methane Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  8)] =      'Methane from Methane Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  8)] =      'Methane from Methane Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  8)] =      'Methane from Methane Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  8)] =      'Methane from Methane Dimer'
TAGL['%s-%s'            % (dbse,  9)] = 'DD-2 Ethene Dimer, D2D'
TAGL['%s-%s-dimer'      % (dbse,  9)] =      'Ethene Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse,  9)] =      'Ethene from Ethene Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse,  9)] =      'Ethene from Ethene Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse,  9)] =      'Ethene from Ethene Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse,  9)] =      'Ethene from Ethene Dimer'
TAGL['%s-%s'            % (dbse, 10)] = 'DD-3 Benzene-Methane Complex, C3'
TAGL['%s-%s-dimer'      % (dbse, 10)] =      'Benzene-Methane Complex'
TAGL['%s-%s-monoA-CP'   % (dbse, 10)] =      'Benzene from Benzene-Methane Complex'
TAGL['%s-%s-monoB-CP'   % (dbse, 10)] =      'Methane from Benzene-Methane Complex'
TAGL['%s-%s-monoA-unCP' % (dbse, 10)] =      'Benzene from Benzene-Methane Complex'
TAGL['%s-%s-monoB-unCP' % (dbse, 10)] =      'Methane from Benzene-Methane Complex'
TAGL['%s-%s'            % (dbse, 11)] = 'DD-4 Benzene Dimer Parallel-Disp, C2H'
TAGL['%s-%s-dimer'      % (dbse, 11)] =      'Benzene Dimer PD'
TAGL['%s-%s-monoA-CP'   % (dbse, 11)] =      'Benzene from Benzene Dimer PD'
TAGL['%s-%s-monoB-CP'   % (dbse, 11)] =      'Benzene from Benzene Dimer PD'
TAGL['%s-%s-monoA-unCP' % (dbse, 11)] =      'Benzene from Benzene Dimer PD'
TAGL['%s-%s-monoB-unCP' % (dbse, 11)] =      'Benzene from Benzene Dimer PD'
TAGL['%s-%s'            % (dbse, 12)] = 'DD-6 Pyrazine Dimer, CS'
TAGL['%s-%s-dimer'      % (dbse, 12)] =      'Pyrazine Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse, 12)] =      'Pyrazine from Pyrazine Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse, 12)] =      'Pyrazine from Pyrazine Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse, 12)] =      'Pyrazine from Pyrazine Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse, 12)] =      'Pyrazine from Pyrazine Dimer'
TAGL['%s-%s'            % (dbse, 13)] = 'MX-5 Uracil Dimer Stack, C2'
TAGL['%s-%s-dimer'      % (dbse, 13)] =      'Uracil Dimer Stack'
TAGL['%s-%s-monoA-CP'   % (dbse, 13)] =      'Uracil from Uracil Dimer Stack'
TAGL['%s-%s-monoB-CP'   % (dbse, 13)] =      'Uracil from Uracil Dimer Stack'
TAGL['%s-%s-monoA-unCP' % (dbse, 13)] =      'Uracil from Uracil Dimer Stack'
TAGL['%s-%s-monoB-unCP' % (dbse, 13)] =      'Uracil from Uracil Dimer Stack'
TAGL['%s-%s'            % (dbse, 14)] = 'DD-7 Indole-Benzene Complex Stack, C1'
TAGL['%s-%s-dimer'      % (dbse, 14)] =      'Indole-Benzene Complex Stack'
TAGL['%s-%s-monoA-CP'   % (dbse, 14)] =      'Benzene from Indole-Benzene Complex Stack'
TAGL['%s-%s-monoB-CP'   % (dbse, 14)] =      'Indole from Indole-Benzene Complex Stack'
TAGL['%s-%s-monoA-unCP' % (dbse, 14)] =      'Benzene from Indole-Benzene Complex Stack'
TAGL['%s-%s-monoB-unCP' % (dbse, 14)] =      'Indole from Indole-Benzene Complex Stack'
TAGL['%s-%s'            % (dbse, 15)] = 'MX-8 Adenine-Thymine Complex Stack, C1'
TAGL['%s-%s-dimer'      % (dbse, 15)] =      'Adenine-Thymine Complex Stack'
TAGL['%s-%s-monoA-CP'   % (dbse, 15)] =      'Adenine from Adenine-Thymine Complex Stack'
TAGL['%s-%s-monoB-CP'   % (dbse, 15)] =      'Thymine from Adenine-Thymine Complex Stack'
TAGL['%s-%s-monoA-unCP' % (dbse, 15)] =      'Adenine from Adenine-Thymine Complex Stack'
TAGL['%s-%s-monoB-unCP' % (dbse, 15)] =      'Thymine from Adenine-Thymine Complex Stack'
TAGL['%s-%s'            % (dbse, 16)] = 'MX-1 Ethene-Ethine Complex, C2V'
TAGL['%s-%s-dimer'      % (dbse, 16)] =      'Ethene-Ethine Complex'
TAGL['%s-%s-monoA-CP'   % (dbse, 16)] =      'Ethene from Ethene-Ethine Complex'
TAGL['%s-%s-monoB-CP'   % (dbse, 16)] =      'Ethine from Ethene-Ethine Complex'
TAGL['%s-%s-monoA-unCP' % (dbse, 16)] =      'Ethene from Ethene-Ethine Complex'
TAGL['%s-%s-monoB-unCP' % (dbse, 16)] =      'Ethine from Ethene-Ethine Complex'
TAGL['%s-%s'            % (dbse, 17)] = 'MX-2 Benzene-Water Complex, CS'
TAGL['%s-%s-dimer'      % (dbse, 17)] =      'Benzene-Water Complex'
TAGL['%s-%s-monoA-CP'   % (dbse, 17)] =      'Benzene from Benzene-Water Complex'
TAGL['%s-%s-monoB-CP'   % (dbse, 17)] =      'Water from Benzene-Water Complex'
TAGL['%s-%s-monoA-unCP' % (dbse, 17)] =      'Benzene from Benzene-Water Complex'
TAGL['%s-%s-monoB-unCP' % (dbse, 17)] =      'Water from Benzene-Water Complex'
TAGL['%s-%s'            % (dbse, 18)] = 'MX-3 Benzene-Ammonia Complex, CS'
TAGL['%s-%s-dimer'      % (dbse, 18)] =      'Benzene-Ammonia Complex'
TAGL['%s-%s-monoA-CP'   % (dbse, 18)] =      'Benzene from Benzene-Ammonia Complex'
TAGL['%s-%s-monoB-CP'   % (dbse, 18)] =      'Ammonia from Benzene-Ammonia Complex'
TAGL['%s-%s-monoA-unCP' % (dbse, 18)] =      'Benzene from Benzene-Ammonia Complex'
TAGL['%s-%s-monoB-unCP' % (dbse, 18)] =      'Ammonia from Benzene-Ammonia Complex'
TAGL['%s-%s'            % (dbse, 19)] = 'MX-4 Benzene-HCN Complex, CS'
TAGL['%s-%s-dimer'      % (dbse, 19)] =      'Benzene-HCN Complex'
TAGL['%s-%s-monoA-CP'   % (dbse, 19)] =      'Benzene from Benzene-HCN Complex'
TAGL['%s-%s-monoB-CP'   % (dbse, 19)] =      'HCN from Benzene-HCN Complex'
TAGL['%s-%s-monoA-unCP' % (dbse, 19)] =      'Benzene from Benzene-HCN Complex'
TAGL['%s-%s-monoB-unCP' % (dbse, 19)] =      'HCN from Benzene-HCN Complex'
TAGL['%s-%s'            % (dbse, 20)] = 'DD-5 Benzene Dimer T-Shape, C2V'
TAGL['%s-%s-dimer'      % (dbse, 20)] =      'Benzene Dimer T-Shape'
TAGL['%s-%s-monoA-CP'   % (dbse, 20)] =      'Benzene from Benzene Dimer T-Shape'
TAGL['%s-%s-monoB-CP'   % (dbse, 20)] =      'Benzene from Benzene Dimer T-Shape'
TAGL['%s-%s-monoA-unCP' % (dbse, 20)] =      'Benzene from Benzene Dimer T-Shape'
TAGL['%s-%s-monoB-unCP' % (dbse, 20)] =      'Benzene from Benzene Dimer T-Shape'
TAGL['%s-%s'            % (dbse, 21)] = 'MX-6 Indole-Benzene Complex T-Shape, C1'
TAGL['%s-%s-dimer'      % (dbse, 21)] =      'Indole-Benzene Complex T-Shape'
TAGL['%s-%s-monoA-CP'   % (dbse, 21)] =      'Benzene from Indole-Benzene Complex T-Shape'
TAGL['%s-%s-monoB-CP'   % (dbse, 21)] =      'Indole from Indole-Benzene Complex T-Shape'
TAGL['%s-%s-monoA-unCP' % (dbse, 21)] =      'Benzene from Indole-Benzene Complex T-Shape'
TAGL['%s-%s-monoB-unCP' % (dbse, 21)] =      'Indole from Indole-Benzene Complex T-Shape'
TAGL['%s-%s'            % (dbse, 22)] = 'MX-7 Phenol Dimer, C1'
TAGL['%s-%s-dimer'      % (dbse, 22)] =      'Phenol Dimer'
TAGL['%s-%s-monoA-CP'   % (dbse, 22)] =      'Phenol from Phenol Dimer'
TAGL['%s-%s-monoB-CP'   % (dbse, 22)] =      'Phenol from Phenol Dimer'
TAGL['%s-%s-monoA-unCP' % (dbse, 22)] =      'Phenol from Phenol Dimer'
TAGL['%s-%s-monoB-unCP' % (dbse, 22)] =      'Phenol from Phenol Dimer'

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-dimer' % (dbse, '1')] = qcdb.Molecule("""
efp ammonia
--
efp ammonia
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2')] = qcdb.Molecule("""
efp water
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3')] = qcdb.Molecule("""
efp formicacid
--
efp formicacid
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4')] = qcdb.Molecule("""
efp formamide
--
efp formamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5')] = qcdb.Molecule("""
efp uracil
--
efp uracil
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6')] = qcdb.Molecule("""
efp pyridone
--
efp 2aminopyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7')] = qcdb.Molecule("""
efp adenine-wc
--
efp thymine-wc
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8')] = qcdb.Molecule("""
efp methane
--
efp methane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9')] = qcdb.Molecule("""
efp ethene
--
efp ethene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10')] = qcdb.Molecule("""
efp benzene
--
efp methane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11')] = qcdb.Molecule("""
efp benzene
--
efp benzene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12')] = qcdb.Molecule("""
efp pyrazine
--
efp pyrazine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13')] = qcdb.Molecule("""
efp uracil
--
efp uracil
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14')] = qcdb.Molecule("""
efp benzene
--
efp indole
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15')] = qcdb.Molecule("""
efp adenine-stack
--
efp thymine-stack
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16')] = qcdb.Molecule("""
efp ethene
--
efp ethyne
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17')] = qcdb.Molecule("""
efp benzene
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18')] = qcdb.Molecule("""
efp benzene
--
efp ammonia
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19')] = qcdb.Molecule("""
efp benzene
--
efp hydrogencyanide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20')] = qcdb.Molecule("""
efp benzene
--
efp benzene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21')] = qcdb.Molecule("""
efp benzene
--
efp indole
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22')] = qcdb.Molecule("""
efp phenol
--
efp phenol
units angstrom
""")

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)

#########################################################################

# <<< Supplementary Quantum Chemical Results >>>
DATA = {}

DATA['NUCLEAR REPULSION ENERGY'] = {}
DATA['NUCLEAR REPULSION ENERGY']['S22-1-dimer'] = 40.3142398391
DATA['NUCLEAR REPULSION ENERGY']['S22-2-dimer'] = 36.6628478528
DATA['NUCLEAR REPULSION ENERGY']['S22-3-dimer'] = 235.946620315
DATA['NUCLEAR REPULSION ENERGY']['S22-4-dimer'] = 230.794855209
DATA['NUCLEAR REPULSION ENERGY']['S22-5-dimer'] = 1032.28191174
DATA['NUCLEAR REPULSION ENERGY']['S22-6-dimer'] = 812.288526081
DATA['NUCLEAR REPULSION ENERGY']['S22-7-dimer'] = 1365.23227533
DATA['NUCLEAR REPULSION ENERGY']['S22-8-dimer'] = 41.0002637953
DATA['NUCLEAR REPULSION ENERGY']['S22-9-dimer'] = 102.165309277
DATA['NUCLEAR REPULSION ENERGY']['S22-10-dimer'] = 272.461820278
DATA['NUCLEAR REPULSION ENERGY']['S22-11-dimer'] = 628.972056837
DATA['NUCLEAR REPULSION ENERGY']['S22-12-dimer'] = 654.132022225
DATA['NUCLEAR REPULSION ENERGY']['S22-13-dimer'] = 1161.47069828
DATA['NUCLEAR REPULSION ENERGY']['S22-14-dimer'] = 935.530103761
DATA['NUCLEAR REPULSION ENERGY']['S22-15-dimer'] = 1542.1430487
DATA['NUCLEAR REPULSION ENERGY']['S22-16-dimer'] = 85.1890641964
DATA['NUCLEAR REPULSION ENERGY']['S22-17-dimer'] = 273.329424698
DATA['NUCLEAR REPULSION ENERGY']['S22-18-dimer'] = 273.279614381
DATA['NUCLEAR REPULSION ENERGY']['S22-19-dimer'] = 303.281397519
DATA['NUCLEAR REPULSION ENERGY']['S22-20-dimer'] = 592.416645285
DATA['NUCLEAR REPULSION ENERGY']['S22-21-dimer'] = 876.919230124
DATA['NUCLEAR REPULSION ENERGY']['S22-22-dimer'] = 805.117733746
DATA['NUCLEAR REPULSION ENERGY']['S22-1-monoA-unCP'] = 11.9474317239
DATA['NUCLEAR REPULSION ENERGY']['S22-1-monoB-unCP'] = 11.9474317239
DATA['NUCLEAR REPULSION ENERGY']['S22-2-monoA-unCP'] = 9.16383014597
DATA['NUCLEAR REPULSION ENERGY']['S22-2-monoB-unCP'] = 9.1780389049
DATA['NUCLEAR REPULSION ENERGY']['S22-3-monoA-unCP'] = 70.1157833033
DATA['NUCLEAR REPULSION ENERGY']['S22-3-monoB-unCP'] = 70.1157833033
DATA['NUCLEAR REPULSION ENERGY']['S22-4-monoA-unCP'] = 71.0728637475
DATA['NUCLEAR REPULSION ENERGY']['S22-4-monoB-unCP'] = 71.0728637475
DATA['NUCLEAR REPULSION ENERGY']['S22-5-monoA-unCP'] = 357.226773232
DATA['NUCLEAR REPULSION ENERGY']['S22-5-monoB-unCP'] = 357.226773232
DATA['NUCLEAR REPULSION ENERGY']['S22-6-monoA-unCP'] = 275.701873893
DATA['NUCLEAR REPULSION ENERGY']['S22-6-monoB-unCP'] = 275.671980226
DATA['NUCLEAR REPULSION ENERGY']['S22-7-monoA-unCP'] = 503.396306786
DATA['NUCLEAR REPULSION ENERGY']['S22-7-monoB-unCP'] = 440.301569251
DATA['NUCLEAR REPULSION ENERGY']['S22-8-monoA-unCP'] = 13.4480422656
DATA['NUCLEAR REPULSION ENERGY']['S22-8-monoB-unCP'] = 13.4480422656
DATA['NUCLEAR REPULSION ENERGY']['S22-9-monoA-unCP'] = 33.3602695815
DATA['NUCLEAR REPULSION ENERGY']['S22-9-monoB-unCP'] = 33.3602695815
DATA['NUCLEAR REPULSION ENERGY']['S22-10-monoA-unCP'] = 203.707991166
DATA['NUCLEAR REPULSION ENERGY']['S22-10-monoB-unCP'] = 13.4855266506
DATA['NUCLEAR REPULSION ENERGY']['S22-11-monoA-unCP'] = 203.71093056
DATA['NUCLEAR REPULSION ENERGY']['S22-11-monoB-unCP'] = 203.71093056
DATA['NUCLEAR REPULSION ENERGY']['S22-12-monoA-unCP'] = 208.639691163
DATA['NUCLEAR REPULSION ENERGY']['S22-12-monoB-unCP'] = 208.626286711
DATA['NUCLEAR REPULSION ENERGY']['S22-13-monoA-unCP'] = 357.160450068
DATA['NUCLEAR REPULSION ENERGY']['S22-13-monoB-unCP'] = 357.160450068
DATA['NUCLEAR REPULSION ENERGY']['S22-14-monoA-unCP'] = 203.669533561
DATA['NUCLEAR REPULSION ENERGY']['S22-14-monoB-unCP'] = 401.143592213
DATA['NUCLEAR REPULSION ENERGY']['S22-15-monoA-unCP'] = 503.365644851
DATA['NUCLEAR REPULSION ENERGY']['S22-15-monoB-unCP'] = 440.147006891
DATA['NUCLEAR REPULSION ENERGY']['S22-16-monoA-unCP'] = 33.3580720823
DATA['NUCLEAR REPULSION ENERGY']['S22-16-monoB-unCP'] = 24.6979461028
DATA['NUCLEAR REPULSION ENERGY']['S22-17-monoA-unCP'] = 203.633716029
DATA['NUCLEAR REPULSION ENERGY']['S22-17-monoB-unCP'] = 9.16734256253
DATA['NUCLEAR REPULSION ENERGY']['S22-18-monoA-unCP'] = 203.672752811
DATA['NUCLEAR REPULSION ENERGY']['S22-18-monoB-unCP'] = 11.9610533611
DATA['NUCLEAR REPULSION ENERGY']['S22-19-monoA-unCP'] = 203.595134421
DATA['NUCLEAR REPULSION ENERGY']['S22-19-monoB-unCP'] = 23.6698792311
DATA['NUCLEAR REPULSION ENERGY']['S22-20-monoA-unCP'] = 203.681438992
DATA['NUCLEAR REPULSION ENERGY']['S22-20-monoB-unCP'] = 203.664080154
DATA['NUCLEAR REPULSION ENERGY']['S22-21-monoA-unCP'] = 203.56582964
DATA['NUCLEAR REPULSION ENERGY']['S22-21-monoB-unCP'] = 401.056606452
DATA['NUCLEAR REPULSION ENERGY']['S22-22-monoA-unCP'] = 271.438700576
DATA['NUCLEAR REPULSION ENERGY']['S22-22-monoB-unCP'] = 271.346177694
DATA['NUCLEAR REPULSION ENERGY']['S22-1-monoA-CP'] = 11.9474317239
DATA['NUCLEAR REPULSION ENERGY']['S22-1-monoB-CP'] = 11.9474317239
DATA['NUCLEAR REPULSION ENERGY']['S22-2-monoA-CP'] = 9.16383014597
DATA['NUCLEAR REPULSION ENERGY']['S22-2-monoB-CP'] = 9.1780389049
DATA['NUCLEAR REPULSION ENERGY']['S22-3-monoA-CP'] = 70.1157833033
DATA['NUCLEAR REPULSION ENERGY']['S22-3-monoB-CP'] = 70.1157833033
DATA['NUCLEAR REPULSION ENERGY']['S22-4-monoA-CP'] = 71.0728637475
DATA['NUCLEAR REPULSION ENERGY']['S22-4-monoB-CP'] = 71.0728637475
DATA['NUCLEAR REPULSION ENERGY']['S22-5-monoA-CP'] = 357.226773232
DATA['NUCLEAR REPULSION ENERGY']['S22-5-monoB-CP'] = 357.226773232
DATA['NUCLEAR REPULSION ENERGY']['S22-6-monoA-CP'] = 275.701873893
DATA['NUCLEAR REPULSION ENERGY']['S22-6-monoB-CP'] = 275.671980226
DATA['NUCLEAR REPULSION ENERGY']['S22-7-monoA-CP'] = 503.396306786
DATA['NUCLEAR REPULSION ENERGY']['S22-7-monoB-CP'] = 440.301569251
DATA['NUCLEAR REPULSION ENERGY']['S22-8-monoA-CP'] = 13.4480422656
DATA['NUCLEAR REPULSION ENERGY']['S22-8-monoB-CP'] = 13.4480422656
DATA['NUCLEAR REPULSION ENERGY']['S22-9-monoA-CP'] = 33.3602695815
DATA['NUCLEAR REPULSION ENERGY']['S22-9-monoB-CP'] = 33.3602695815
DATA['NUCLEAR REPULSION ENERGY']['S22-10-monoA-CP'] = 203.707991166
DATA['NUCLEAR REPULSION ENERGY']['S22-10-monoB-CP'] = 13.4855266506
DATA['NUCLEAR REPULSION ENERGY']['S22-11-monoA-CP'] = 203.71093056
DATA['NUCLEAR REPULSION ENERGY']['S22-11-monoB-CP'] = 203.71093056
DATA['NUCLEAR REPULSION ENERGY']['S22-12-monoA-CP'] = 208.639691163
DATA['NUCLEAR REPULSION ENERGY']['S22-12-monoB-CP'] = 208.626286711
DATA['NUCLEAR REPULSION ENERGY']['S22-13-monoA-CP'] = 357.160450068
DATA['NUCLEAR REPULSION ENERGY']['S22-13-monoB-CP'] = 357.160450068
DATA['NUCLEAR REPULSION ENERGY']['S22-14-monoA-CP'] = 203.669533561
DATA['NUCLEAR REPULSION ENERGY']['S22-14-monoB-CP'] = 401.143592213
DATA['NUCLEAR REPULSION ENERGY']['S22-15-monoA-CP'] = 503.365644851
DATA['NUCLEAR REPULSION ENERGY']['S22-15-monoB-CP'] = 440.147006891
DATA['NUCLEAR REPULSION ENERGY']['S22-16-monoA-CP'] = 33.3580720823
DATA['NUCLEAR REPULSION ENERGY']['S22-16-monoB-CP'] = 24.6979461028
DATA['NUCLEAR REPULSION ENERGY']['S22-17-monoA-CP'] = 203.633716029
DATA['NUCLEAR REPULSION ENERGY']['S22-17-monoB-CP'] = 9.16734256253
DATA['NUCLEAR REPULSION ENERGY']['S22-18-monoA-CP'] = 203.672752811
DATA['NUCLEAR REPULSION ENERGY']['S22-18-monoB-CP'] = 11.9610533611
DATA['NUCLEAR REPULSION ENERGY']['S22-19-monoA-CP'] = 203.595134421
DATA['NUCLEAR REPULSION ENERGY']['S22-19-monoB-CP'] = 23.6698792311
DATA['NUCLEAR REPULSION ENERGY']['S22-20-monoA-CP'] = 203.681438992
DATA['NUCLEAR REPULSION ENERGY']['S22-20-monoB-CP'] = 203.664080154
DATA['NUCLEAR REPULSION ENERGY']['S22-21-monoA-CP'] = 203.56582964
DATA['NUCLEAR REPULSION ENERGY']['S22-21-monoB-CP'] = 401.056606452
DATA['NUCLEAR REPULSION ENERGY']['S22-22-monoA-CP'] = 271.438700576
DATA['NUCLEAR REPULSION ENERGY']['S22-22-monoB-CP'] = 271.346177694
