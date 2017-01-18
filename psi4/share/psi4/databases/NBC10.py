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
# @END LICENSE
#

"""
| Database (Sherrill) of interaction energies for dissociation curves of dispersion-bound bimolecular complexes.
| Geometries and Reference interaction energies from the following articles:
|   Benzene Dimers from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Hydrogen Sulfide from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Methane from Sherrill et al. JPCA 113 10146 (2009).
|   Methane Dimer from Takatani et al. PCCP 9 6106 (2007).
|   Pyridine Dimers from Hohenstein et al. JPCA 113 878 (2009).
|   Collection into NBC10 from Burns et al. JCP 134 084107 (2011).
|   Revised reference interaction energies (NBC10A) from Marshall et al. JCP 135 194102 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **benchmark**

  - ``'NBC100'`` Burns et al. JCP 134 084107 (2011).
  - |dl| ``'NBC10A'`` |dr| Marshall et al. JCP 135 194102 (2011).

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'BzBz_S'`` dissociation curve for benzene dimer, sandwich
  - ``'BzBz_T'`` dissociation curve for benzene dimer, t-shaped
  - ``'BzBz_PD34'`` dissociation curve for benzene dimer, parallel displaced by 3.4A
  - ``'BzH2S'`` dissociation curve for benzene-H2S
  - ``'BzMe'`` dissociation curve for benzene-methane
  - ``'MeMe'`` dissociation curve for methane dimer
  - ``'PyPy_S2'`` dissociation curve for pyridine dimer, sandwich
  - ``'PyPy_T3'`` dissociation curve for pyridine dimer, t-shaped
  - ``'BzBz_PD32'`` dissociation curve for benzene dimer, parallel displaced by 3.2A
  - ``'BzBz_PD36'`` dissociation curve for benzene dimer, parallel displaced by 3.6A

"""
import re
import qcdb

# <<< NBC10 Database Module >>>
dbse = 'NBC1'

# <<< Database Members >>>
BzBz_S = []
dist = [3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0]
for d in dist:
    BzBz_S.append('BzBz_S-' + str(d))

BzBz_T = []
dist = [4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 6.0, 6.5, 7.0, 7.5, 8.0]
for d in dist:
    BzBz_T.append('BzBz_T-' + str(d))

BzBz_PD34 = []
dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
for d in dist:
    BzBz_PD34.append('BzBz_PD34-' + str(d))

BzH2S = []
dist = [3.2, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 4.75, 5.0, 5.25, 5.5, 6.0, 6.5, 7.0, 7.5]
for d in dist:
    BzH2S.append('BzH2S-' + str(d))

BzMe = []
dist = [3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 6.0]
for d in dist:
    BzMe.append('BzMe-' + str(d))

MeMe = []
dist = [3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.6, 4.8, 5.0, 5.4, 5.8]
for d in dist:
    MeMe.append('MeMe-' + str(d))

PyPy_S2 = []
dist = [3.1, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.7, 5.0, 5.5, 6.0, 6.5, 7.0]
for d in dist:
    PyPy_S2.append('PyPy_S2-' + str(d))

PyPy_T3 = []
dist = [4.1, 4.3, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.7, 6.0, 6.5, 7.0, 8.0, 9.0]
for d in dist:
    PyPy_T3.append('PyPy_T3-' + str(d))

BzBz_PD32 = []
dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
for d in dist:
    BzBz_PD32.append('BzBz_PD32-' + str(d))

BzBz_PD36 = []
dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
for d in dist:
    BzBz_PD36.append('BzBz_PD36-' + str(d))

temp = [BzBz_S, BzBz_T, BzBz_PD34, BzH2S, BzMe, MeMe, PyPy_S2, PyPy_T3, BzBz_PD32, BzBz_PD36]
HRXN = sum(temp, [])

HRXN_SM = ['BzMe-6.0', 'MeMe-5.0']
HRXN_LG = ['BzBz_T-5.0']
HRXN_EQ = ['BzBz_S-3.9', 'BzBz_T-5.0', 'BzBz_PD34-1.8', 'BzH2S-3.8', 'BzMe-3.8',
           'MeMe-3.6', 'PyPy_S2-3.7', 'PyPy_T3-4.9', 'BzBz_PD32-1.9', 'BzBz_PD36-1.7']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations
for rxn in HRXN:

    if (rxn in BzBz_S) or (rxn in BzBz_PD34) or (rxn in BzBz_PD32) or (rxn in BzBz_PD36):
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse) ]

    elif rxn in BzBz_T:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse) ]

    elif rxn in BzH2S:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -1,
                                          '%s-H2S-mono-unCP' % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse),
                                          '%s-H2S-mono-unCP' % (dbse) ]

    elif rxn in BzMe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz2-mono-unCP' % (dbse)      : -1,
                                          '%s-Me-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz2-mono-unCP' % (dbse),
                                          '%s-Me-mono-unCP'  % (dbse) ]

    elif rxn in MeMe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Me-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Me-mono-unCP'  % (dbse) ]

    elif rxn in PyPy_S2:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Py-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Py-mono-unCP'  % (dbse) ]

    elif rxn in PyPy_T3:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Py-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Py-mono-unCP'  % (dbse) ]

# <<< Reference Values >>>
BIND = {}
# Original publication
BIND_NBC100 = {}
BIND_NBC100['%s-BzBz_S-3.2'  % (dbse)] =  3.522
BIND_NBC100['%s-BzBz_S-3.3'  % (dbse)] =  1.535
BIND_NBC100['%s-BzBz_S-3.4'  % (dbse)] =  0.189
BIND_NBC100['%s-BzBz_S-3.5'  % (dbse)] = -0.689
BIND_NBC100['%s-BzBz_S-3.6'  % (dbse)] = -1.231
BIND_NBC100['%s-BzBz_S-3.7'  % (dbse)] = -1.535
BIND_NBC100['%s-BzBz_S-3.8'  % (dbse)] = -1.674
BIND_NBC100['%s-BzBz_S-3.9'  % (dbse)] = -1.701  # BzBz_S minimum
BIND_NBC100['%s-BzBz_S-4.0'  % (dbse)] = -1.655
BIND_NBC100['%s-BzBz_S-4.1'  % (dbse)] = -1.565
BIND_NBC100['%s-BzBz_S-4.2'  % (dbse)] = -1.448
BIND_NBC100['%s-BzBz_S-4.5'  % (dbse)] = -1.058
BIND_NBC100['%s-BzBz_S-5.0'  % (dbse)] = -0.542
BIND_NBC100['%s-BzBz_S-5.5'  % (dbse)] = -0.248
BIND_NBC100['%s-BzBz_S-6.0'  % (dbse)] = -0.099
BIND_NBC100['%s-BzBz_S-6.5'  % (dbse)] = -0.028
BIND_NBC100['%s-BzBz_S-10.0' % (dbse)] =  0.018

BIND_NBC100['%s-BzBz_T-4.4' % (dbse)] =  0.626
BIND_NBC100['%s-BzBz_T-4.5' % (dbse)] = -0.760
BIND_NBC100['%s-BzBz_T-4.6' % (dbse)] = -1.673
BIND_NBC100['%s-BzBz_T-4.7' % (dbse)] = -2.239
BIND_NBC100['%s-BzBz_T-4.8' % (dbse)] = -2.552
BIND_NBC100['%s-BzBz_T-4.9' % (dbse)] = -2.687
BIND_NBC100['%s-BzBz_T-5.0' % (dbse)] = -2.698  # BzBz_T minimum
BIND_NBC100['%s-BzBz_T-5.1' % (dbse)] = -2.627
BIND_NBC100['%s-BzBz_T-5.2' % (dbse)] = -2.503
BIND_NBC100['%s-BzBz_T-5.3' % (dbse)] = -2.349
BIND_NBC100['%s-BzBz_T-5.4' % (dbse)] = -2.179
BIND_NBC100['%s-BzBz_T-5.5' % (dbse)] = -2.005
BIND_NBC100['%s-BzBz_T-5.6' % (dbse)] = -1.833
BIND_NBC100['%s-BzBz_T-6.0' % (dbse)] = -1.242
BIND_NBC100['%s-BzBz_T-6.5' % (dbse)] = -0.752
BIND_NBC100['%s-BzBz_T-7.0' % (dbse)] = -0.468
BIND_NBC100['%s-BzBz_T-7.5' % (dbse)] = -0.302
BIND_NBC100['%s-BzBz_T-8.0' % (dbse)] = -0.203

BIND_NBC100['%s-BzBz_PD34-0.2' % (dbse)] =  0.070
BIND_NBC100['%s-BzBz_PD34-0.4' % (dbse)] = -0.257
BIND_NBC100['%s-BzBz_PD34-0.6' % (dbse)] = -0.728
BIND_NBC100['%s-BzBz_PD34-0.8' % (dbse)] = -1.260
BIND_NBC100['%s-BzBz_PD34-1.0' % (dbse)] = -1.766
BIND_NBC100['%s-BzBz_PD34-1.2' % (dbse)] = -2.179
BIND_NBC100['%s-BzBz_PD34-1.4' % (dbse)] = -2.466
BIND_NBC100['%s-BzBz_PD34-1.5' % (dbse)] = -2.557
BIND_NBC100['%s-BzBz_PD34-1.6' % (dbse)] = -2.614
BIND_NBC100['%s-BzBz_PD34-1.7' % (dbse)] = -2.640
BIND_NBC100['%s-BzBz_PD34-1.8' % (dbse)] = -2.643  # BzBz_PD34 minimum
BIND_NBC100['%s-BzBz_PD34-1.9' % (dbse)] = -2.624
BIND_NBC100['%s-BzBz_PD34-2.0' % (dbse)] = -2.587
BIND_NBC100['%s-BzBz_PD34-2.2' % (dbse)] = -2.479
BIND_NBC100['%s-BzBz_PD34-2.4' % (dbse)] = -2.356
BIND_NBC100['%s-BzBz_PD34-2.6' % (dbse)] = -2.242
BIND_NBC100['%s-BzBz_PD34-2.8' % (dbse)] = -2.147
BIND_NBC100['%s-BzBz_PD34-3.0' % (dbse)] = -2.079

BIND_NBC100['%s-BzH2S-3.2'  % (dbse)] =  1.250
BIND_NBC100['%s-BzH2S-3.4'  % (dbse)] = -1.570
BIND_NBC100['%s-BzH2S-3.5'  % (dbse)] = -2.256
BIND_NBC100['%s-BzH2S-3.6'  % (dbse)] = -2.638
BIND_NBC100['%s-BzH2S-3.7'  % (dbse)] = -2.808
BIND_NBC100['%s-BzH2S-3.8'  % (dbse)] = -2.834  # BzH2S minimum
BIND_NBC100['%s-BzH2S-3.9'  % (dbse)] = -2.766
BIND_NBC100['%s-BzH2S-4.0'  % (dbse)] = -2.639
BIND_NBC100['%s-BzH2S-4.1'  % (dbse)] = -2.478
BIND_NBC100['%s-BzH2S-4.2'  % (dbse)] = -2.301
BIND_NBC100['%s-BzH2S-4.5'  % (dbse)] = -1.770
BIND_NBC100['%s-BzH2S-4.75' % (dbse)] = -1.393
BIND_NBC100['%s-BzH2S-5.0'  % (dbse)] = -1.093
BIND_NBC100['%s-BzH2S-5.25' % (dbse)] = -0.861
BIND_NBC100['%s-BzH2S-5.5'  % (dbse)] = -0.684
BIND_NBC100['%s-BzH2S-6.0'  % (dbse)] = -0.446
BIND_NBC100['%s-BzH2S-6.5'  % (dbse)] = -0.302
BIND_NBC100['%s-BzH2S-7.0'  % (dbse)] = -0.214
BIND_NBC100['%s-BzH2S-7.5'  % (dbse)] = -0.155

BIND_NBC100['%s-BzMe-3.2' % (dbse)] =  0.717
BIND_NBC100['%s-BzMe-3.3' % (dbse)] = -0.183
BIND_NBC100['%s-BzMe-3.4' % (dbse)] = -0.774
BIND_NBC100['%s-BzMe-3.5' % (dbse)] = -1.135
BIND_NBC100['%s-BzMe-3.6' % (dbse)] = -1.337
BIND_NBC100['%s-BzMe-3.7' % (dbse)] = -1.432
BIND_NBC100['%s-BzMe-3.8' % (dbse)] = -1.439  # BzMe minimum
BIND_NBC100['%s-BzMe-3.9' % (dbse)] = -1.414
BIND_NBC100['%s-BzMe-4.0' % (dbse)] = -1.327
BIND_NBC100['%s-BzMe-4.1' % (dbse)] = -1.232
BIND_NBC100['%s-BzMe-4.2' % (dbse)] = -1.138
BIND_NBC100['%s-BzMe-4.4' % (dbse)] = -0.950
BIND_NBC100['%s-BzMe-4.6' % (dbse)] = -0.760
BIND_NBC100['%s-BzMe-4.8' % (dbse)] = -0.606
BIND_NBC100['%s-BzMe-5.0' % (dbse)] = -0.475
BIND_NBC100['%s-BzMe-5.2' % (dbse)] = -0.370
BIND_NBC100['%s-BzMe-5.4' % (dbse)] = -0.286
BIND_NBC100['%s-BzMe-5.6' % (dbse)] = -0.230
BIND_NBC100['%s-BzMe-6.0' % (dbse)] = -0.141

BIND_NBC100['%s-MeMe-3.2' % (dbse)] =  0.069
BIND_NBC100['%s-MeMe-3.3' % (dbse)] = -0.239
BIND_NBC100['%s-MeMe-3.4' % (dbse)] = -0.417
BIND_NBC100['%s-MeMe-3.5' % (dbse)] = -0.508
BIND_NBC100['%s-MeMe-3.6' % (dbse)] = -0.541  # MeMe minimum
BIND_NBC100['%s-MeMe-3.7' % (dbse)] = -0.539
BIND_NBC100['%s-MeMe-3.8' % (dbse)] = -0.515
BIND_NBC100['%s-MeMe-3.9' % (dbse)] = -0.480
BIND_NBC100['%s-MeMe-4.0' % (dbse)] = -0.439
BIND_NBC100['%s-MeMe-4.1' % (dbse)] = -0.396
BIND_NBC100['%s-MeMe-4.2' % (dbse)] = -0.354
BIND_NBC100['%s-MeMe-4.3' % (dbse)] = -0.315
BIND_NBC100['%s-MeMe-4.4' % (dbse)] = -0.279
BIND_NBC100['%s-MeMe-4.6' % (dbse)] = -0.217
BIND_NBC100['%s-MeMe-4.8' % (dbse)] = -0.168
BIND_NBC100['%s-MeMe-5.0' % (dbse)] = -0.130
BIND_NBC100['%s-MeMe-5.4' % (dbse)] = -0.080
BIND_NBC100['%s-MeMe-5.8' % (dbse)] = -0.050

BIND_NBC100['%s-PyPy_S2-3.1' % (dbse)] =  2.442
BIND_NBC100['%s-PyPy_S2-3.3' % (dbse)] = -1.125
BIND_NBC100['%s-PyPy_S2-3.4' % (dbse)] = -2.016
BIND_NBC100['%s-PyPy_S2-3.5' % (dbse)] = -2.534
BIND_NBC100['%s-PyPy_S2-3.6' % (dbse)] = -2.791
BIND_NBC100['%s-PyPy_S2-3.7' % (dbse)] = -2.870  # PyPy_S2 minimum
BIND_NBC100['%s-PyPy_S2-3.8' % (dbse)] = -2.832
BIND_NBC100['%s-PyPy_S2-3.9' % (dbse)] = -2.719
BIND_NBC100['%s-PyPy_S2-4.0' % (dbse)] = -2.561
BIND_NBC100['%s-PyPy_S2-4.1' % (dbse)] = -2.381
BIND_NBC100['%s-PyPy_S2-4.2' % (dbse)] = -2.192
BIND_NBC100['%s-PyPy_S2-4.3' % (dbse)] = -2.005
BIND_NBC100['%s-PyPy_S2-4.4' % (dbse)] = -1.824
BIND_NBC100['%s-PyPy_S2-4.5' % (dbse)] = -1.655
BIND_NBC100['%s-PyPy_S2-4.7' % (dbse)] = -1.354
BIND_NBC100['%s-PyPy_S2-5.0' % (dbse)] = -0.999
BIND_NBC100['%s-PyPy_S2-5.5' % (dbse)] = -0.618
BIND_NBC100['%s-PyPy_S2-6.0' % (dbse)] = -0.402
BIND_NBC100['%s-PyPy_S2-6.5' % (dbse)] = -0.277
BIND_NBC100['%s-PyPy_S2-7.0' % (dbse)] = -0.200

BIND_NBC100['%s-PyPy_T3-4.1' % (dbse)] =  9.340
BIND_NBC100['%s-PyPy_T3-4.3' % (dbse)] =  1.991
BIND_NBC100['%s-PyPy_T3-4.5' % (dbse)] = -1.377
BIND_NBC100['%s-PyPy_T3-4.6' % (dbse)] = -2.203
BIND_NBC100['%s-PyPy_T3-4.7' % (dbse)] = -2.673
BIND_NBC100['%s-PyPy_T3-4.8' % (dbse)] = -2.897
BIND_NBC100['%s-PyPy_T3-4.9' % (dbse)] = -2.954  # PyPy_T3 minimum
BIND_NBC100['%s-PyPy_T3-5.0' % (dbse)] = -2.903
BIND_NBC100['%s-PyPy_T3-5.1' % (dbse)] = -2.784
BIND_NBC100['%s-PyPy_T3-5.2' % (dbse)] = -2.625
BIND_NBC100['%s-PyPy_T3-5.3' % (dbse)] = -2.447
BIND_NBC100['%s-PyPy_T3-5.4' % (dbse)] = -2.263
BIND_NBC100['%s-PyPy_T3-5.5' % (dbse)] = -2.080
BIND_NBC100['%s-PyPy_T3-5.7' % (dbse)] = -1.742
BIND_NBC100['%s-PyPy_T3-6.0' % (dbse)] = -1.324
BIND_NBC100['%s-PyPy_T3-6.5' % (dbse)] = -0.853
BIND_NBC100['%s-PyPy_T3-7.0' % (dbse)] = -0.574
BIND_NBC100['%s-PyPy_T3-8.0' % (dbse)] = -0.296
BIND_NBC100['%s-PyPy_T3-9.0' % (dbse)] = -0.175

BIND_NBC100['%s-BzBz_PD32-0.2' % (dbse)] =  3.301
BIND_NBC100['%s-BzBz_PD32-0.4' % (dbse)] =  2.678
BIND_NBC100['%s-BzBz_PD32-0.6' % (dbse)] =  1.783
BIND_NBC100['%s-BzBz_PD32-0.8' % (dbse)] =  0.781
BIND_NBC100['%s-BzBz_PD32-1.0' % (dbse)] = -0.171
BIND_NBC100['%s-BzBz_PD32-1.2' % (dbse)] = -0.954
BIND_NBC100['%s-BzBz_PD32-1.4' % (dbse)] = -1.508
BIND_NBC100['%s-BzBz_PD32-1.5' % (dbse)] = -1.695
BIND_NBC100['%s-BzBz_PD32-1.6' % (dbse)] = -1.827
BIND_NBC100['%s-BzBz_PD32-1.7' % (dbse)] = -1.911
BIND_NBC100['%s-BzBz_PD32-1.8' % (dbse)] = -1.950
BIND_NBC100['%s-BzBz_PD32-1.9' % (dbse)] = -1.957  # BzBz_PD32 minimum
BIND_NBC100['%s-BzBz_PD32-2.0' % (dbse)] = -1.937
BIND_NBC100['%s-BzBz_PD32-2.2' % (dbse)] = -1.860
BIND_NBC100['%s-BzBz_PD32-2.4' % (dbse)] = -1.767
BIND_NBC100['%s-BzBz_PD32-2.6' % (dbse)] = -1.702
BIND_NBC100['%s-BzBz_PD32-2.8' % (dbse)] = -1.680
BIND_NBC100['%s-BzBz_PD32-3.0' % (dbse)] = -1.705

BIND_NBC100['%s-BzBz_PD36-0.2' % (dbse)] = -1.293
BIND_NBC100['%s-BzBz_PD36-0.4' % (dbse)] = -1.462
BIND_NBC100['%s-BzBz_PD36-0.6' % (dbse)] = -1.708
BIND_NBC100['%s-BzBz_PD36-0.8' % (dbse)] = -1.984
BIND_NBC100['%s-BzBz_PD36-1.0' % (dbse)] = -2.248
BIND_NBC100['%s-BzBz_PD36-1.2' % (dbse)] = -2.458
BIND_NBC100['%s-BzBz_PD36-1.4' % (dbse)] = -2.597
BIND_NBC100['%s-BzBz_PD36-1.5' % (dbse)] = -2.635
BIND_NBC100['%s-BzBz_PD36-1.6' % (dbse)] = -2.652
BIND_NBC100['%s-BzBz_PD36-1.7' % (dbse)] = -2.654  # BzBz_PD36 minimum
BIND_NBC100['%s-BzBz_PD36-1.8' % (dbse)] = -2.642
BIND_NBC100['%s-BzBz_PD36-1.9' % (dbse)] = -2.615
BIND_NBC100['%s-BzBz_PD36-2.0' % (dbse)] = -2.575
BIND_NBC100['%s-BzBz_PD36-2.2' % (dbse)] = -2.473
BIND_NBC100['%s-BzBz_PD36-2.4' % (dbse)] = -2.356
BIND_NBC100['%s-BzBz_PD36-2.6' % (dbse)] = -2.240
BIND_NBC100['%s-BzBz_PD36-2.8' % (dbse)] = -2.130
BIND_NBC100['%s-BzBz_PD36-3.0' % (dbse)] = -2.035
# Current revision
BIND_NBC10A = {}
BIND_NBC10A['%s-BzBz_S-3.2'  % (dbse)] =  3.462
BIND_NBC10A['%s-BzBz_S-3.3'  % (dbse)] =  1.484
BIND_NBC10A['%s-BzBz_S-3.4'  % (dbse)] =  0.147
BIND_NBC10A['%s-BzBz_S-3.5'  % (dbse)] = -0.724
BIND_NBC10A['%s-BzBz_S-3.6'  % (dbse)] = -1.259
BIND_NBC10A['%s-BzBz_S-3.7'  % (dbse)] = -1.558
BIND_NBC10A['%s-BzBz_S-3.8'  % (dbse)] = -1.693
BIND_NBC10A['%s-BzBz_S-3.9'  % (dbse)] = -1.717  # BzBz_S minimum
BIND_NBC10A['%s-BzBz_S-4.0'  % (dbse)] = -1.669
BIND_NBC10A['%s-BzBz_S-4.1'  % (dbse)] = -1.577
BIND_NBC10A['%s-BzBz_S-4.2'  % (dbse)] = -1.459
BIND_NBC10A['%s-BzBz_S-4.5'  % (dbse)] = -1.066
BIND_NBC10A['%s-BzBz_S-5.0'  % (dbse)] = -0.546
BIND_NBC10A['%s-BzBz_S-5.5'  % (dbse)] = -0.251
BIND_NBC10A['%s-BzBz_S-6.0'  % (dbse)] = -0.101
BIND_NBC10A['%s-BzBz_S-6.5'  % (dbse)] = -0.029
BIND_NBC10A['%s-BzBz_S-10.0' % (dbse)] =  0.018

BIND_NBC10A['%s-BzBz_T-4.4' % (dbse)] =  0.617
BIND_NBC10A['%s-BzBz_T-4.5' % (dbse)] = -0.769
BIND_NBC10A['%s-BzBz_T-4.6' % (dbse)] = -1.682
BIND_NBC10A['%s-BzBz_T-4.7' % (dbse)] = -2.246
BIND_NBC10A['%s-BzBz_T-4.8' % (dbse)] = -2.559
BIND_NBC10A['%s-BzBz_T-4.9' % (dbse)] = -2.693
BIND_NBC10A['%s-BzBz_T-5.0' % (dbse)] = -2.703  # BzBz_T minimum
BIND_NBC10A['%s-BzBz_T-5.1' % (dbse)] = -2.630
BIND_NBC10A['%s-BzBz_T-5.2' % (dbse)] = -2.506
BIND_NBC10A['%s-BzBz_T-5.3' % (dbse)] = -2.351
BIND_NBC10A['%s-BzBz_T-5.4' % (dbse)] = -2.181
BIND_NBC10A['%s-BzBz_T-5.5' % (dbse)] = -2.006
BIND_NBC10A['%s-BzBz_T-5.6' % (dbse)] = -1.834
BIND_NBC10A['%s-BzBz_T-6.0' % (dbse)] = -1.242
BIND_NBC10A['%s-BzBz_T-6.5' % (dbse)] = -0.752
BIND_NBC10A['%s-BzBz_T-7.0' % (dbse)] = -0.468
BIND_NBC10A['%s-BzBz_T-7.5' % (dbse)] = -0.302
BIND_NBC10A['%s-BzBz_T-8.0' % (dbse)] = -0.203

BIND_NBC10A['%s-BzBz_PD34-0.2' % (dbse)] =  0.029
BIND_NBC10A['%s-BzBz_PD34-0.4' % (dbse)] = -0.298
BIND_NBC10A['%s-BzBz_PD34-0.6' % (dbse)] = -0.768
BIND_NBC10A['%s-BzBz_PD34-0.8' % (dbse)] = -1.298
BIND_NBC10A['%s-BzBz_PD34-1.0' % (dbse)] = -1.802
BIND_NBC10A['%s-BzBz_PD34-1.2' % (dbse)] = -2.213
BIND_NBC10A['%s-BzBz_PD34-1.4' % (dbse)] = -2.497
BIND_NBC10A['%s-BzBz_PD34-1.5' % (dbse)] = -2.586
BIND_NBC10A['%s-BzBz_PD34-1.6' % (dbse)] = -2.643
BIND_NBC10A['%s-BzBz_PD34-1.7' % (dbse)] = -2.668
BIND_NBC10A['%s-BzBz_PD34-1.8' % (dbse)] = -2.670  # BzBz_PD34 minimum
BIND_NBC10A['%s-BzBz_PD34-1.9' % (dbse)] = -2.649
BIND_NBC10A['%s-BzBz_PD34-2.0' % (dbse)] = -2.611
BIND_NBC10A['%s-BzBz_PD34-2.2' % (dbse)] = -2.501
BIND_NBC10A['%s-BzBz_PD34-2.4' % (dbse)] = -2.377
BIND_NBC10A['%s-BzBz_PD34-2.6' % (dbse)] = -2.260
BIND_NBC10A['%s-BzBz_PD34-2.8' % (dbse)] = -2.163
BIND_NBC10A['%s-BzBz_PD34-3.0' % (dbse)] = -2.093

BIND_NBC10A['%s-BzH2S-3.2'  % (dbse)] =  1.236
BIND_NBC10A['%s-BzH2S-3.4'  % (dbse)] = -1.584
BIND_NBC10A['%s-BzH2S-3.5'  % (dbse)] = -2.269
BIND_NBC10A['%s-BzH2S-3.6'  % (dbse)] = -2.649
BIND_NBC10A['%s-BzH2S-3.7'  % (dbse)] = -2.818
BIND_NBC10A['%s-BzH2S-3.8'  % (dbse)] = -2.843  # BzH2S minimum
BIND_NBC10A['%s-BzH2S-3.9'  % (dbse)] = -2.773
BIND_NBC10A['%s-BzH2S-4.0'  % (dbse)] = -2.645
BIND_NBC10A['%s-BzH2S-4.1'  % (dbse)] = -2.483
BIND_NBC10A['%s-BzH2S-4.2'  % (dbse)] = -2.305
BIND_NBC10A['%s-BzH2S-4.5'  % (dbse)] = -1.771
BIND_NBC10A['%s-BzH2S-4.75' % (dbse)] = -1.393
BIND_NBC10A['%s-BzH2S-5.0'  % (dbse)] = -1.092
BIND_NBC10A['%s-BzH2S-5.25' % (dbse)] = -0.859
BIND_NBC10A['%s-BzH2S-5.5'  % (dbse)] = -0.682
BIND_NBC10A['%s-BzH2S-6.0'  % (dbse)] = -0.444
BIND_NBC10A['%s-BzH2S-6.5'  % (dbse)] = -0.301
BIND_NBC10A['%s-BzH2S-7.0'  % (dbse)] = -0.212
BIND_NBC10A['%s-BzH2S-7.5'  % (dbse)] = -0.154

BIND_NBC10A['%s-BzMe-3.2' % (dbse)] =  0.686
BIND_NBC10A['%s-BzMe-3.3' % (dbse)] = -0.213
BIND_NBC10A['%s-BzMe-3.4' % (dbse)] = -0.805
BIND_NBC10A['%s-BzMe-3.5' % (dbse)] = -1.173
BIND_NBC10A['%s-BzMe-3.6' % (dbse)] = -1.378
BIND_NBC10A['%s-BzMe-3.7' % (dbse)] = -1.470
BIND_NBC10A['%s-BzMe-3.8' % (dbse)] = -1.484  # BzMe minimum
BIND_NBC10A['%s-BzMe-3.9' % (dbse)] = -1.445
BIND_NBC10A['%s-BzMe-4.0' % (dbse)] = -1.374
BIND_NBC10A['%s-BzMe-4.1' % (dbse)] = -1.284
BIND_NBC10A['%s-BzMe-4.2' % (dbse)] = -1.185
BIND_NBC10A['%s-BzMe-4.4' % (dbse)] = -0.984
BIND_NBC10A['%s-BzMe-4.6' % (dbse)] = -0.800
BIND_NBC10A['%s-BzMe-4.8' % (dbse)] = -0.643
BIND_NBC10A['%s-BzMe-5.0' % (dbse)] = -0.515
BIND_NBC10A['%s-BzMe-5.2' % (dbse)] = -0.413
BIND_NBC10A['%s-BzMe-5.4' % (dbse)] = -0.332
BIND_NBC10A['%s-BzMe-5.6' % (dbse)] = -0.268
BIND_NBC10A['%s-BzMe-6.0' % (dbse)] = -0.177

BIND_NBC10A['%s-MeMe-3.2' % (dbse)] =  0.069
BIND_NBC10A['%s-MeMe-3.3' % (dbse)] = -0.239
BIND_NBC10A['%s-MeMe-3.4' % (dbse)] = -0.417
BIND_NBC10A['%s-MeMe-3.5' % (dbse)] = -0.508
BIND_NBC10A['%s-MeMe-3.6' % (dbse)] = -0.541  # MeMe minimum
BIND_NBC10A['%s-MeMe-3.7' % (dbse)] = -0.539
BIND_NBC10A['%s-MeMe-3.8' % (dbse)] = -0.515
BIND_NBC10A['%s-MeMe-3.9' % (dbse)] = -0.480
BIND_NBC10A['%s-MeMe-4.0' % (dbse)] = -0.439
BIND_NBC10A['%s-MeMe-4.1' % (dbse)] = -0.396
BIND_NBC10A['%s-MeMe-4.2' % (dbse)] = -0.354
BIND_NBC10A['%s-MeMe-4.3' % (dbse)] = -0.315
BIND_NBC10A['%s-MeMe-4.4' % (dbse)] = -0.279
BIND_NBC10A['%s-MeMe-4.6' % (dbse)] = -0.217
BIND_NBC10A['%s-MeMe-4.8' % (dbse)] = -0.168
BIND_NBC10A['%s-MeMe-5.0' % (dbse)] = -0.130
BIND_NBC10A['%s-MeMe-5.4' % (dbse)] = -0.080
BIND_NBC10A['%s-MeMe-5.8' % (dbse)] = -0.050

BIND_NBC10A['%s-PyPy_S2-3.1' % (dbse)] =  2.387
BIND_NBC10A['%s-PyPy_S2-3.3' % (dbse)] = -1.165
BIND_NBC10A['%s-PyPy_S2-3.4' % (dbse)] = -2.050
BIND_NBC10A['%s-PyPy_S2-3.5' % (dbse)] = -2.562
BIND_NBC10A['%s-PyPy_S2-3.6' % (dbse)] = -2.815
BIND_NBC10A['%s-PyPy_S2-3.7' % (dbse)] = -2.890  # PyPy_S2 minimum
BIND_NBC10A['%s-PyPy_S2-3.8' % (dbse)] = -2.849
BIND_NBC10A['%s-PyPy_S2-3.9' % (dbse)] = -2.733
BIND_NBC10A['%s-PyPy_S2-4.0' % (dbse)] = -2.573
BIND_NBC10A['%s-PyPy_S2-4.1' % (dbse)] = -2.391
BIND_NBC10A['%s-PyPy_S2-4.2' % (dbse)] = -2.201
BIND_NBC10A['%s-PyPy_S2-4.3' % (dbse)] = -2.012
BIND_NBC10A['%s-PyPy_S2-4.4' % (dbse)] = -1.830
BIND_NBC10A['%s-PyPy_S2-4.5' % (dbse)] = -1.660
BIND_NBC10A['%s-PyPy_S2-4.7' % (dbse)] = -1.357
BIND_NBC10A['%s-PyPy_S2-5.0' % (dbse)] = -1.002
BIND_NBC10A['%s-PyPy_S2-5.5' % (dbse)] = -0.619
BIND_NBC10A['%s-PyPy_S2-6.0' % (dbse)] = -0.402
BIND_NBC10A['%s-PyPy_S2-6.5' % (dbse)] = -0.276
BIND_NBC10A['%s-PyPy_S2-7.0' % (dbse)] = -0.200

BIND_NBC10A['%s-PyPy_T3-4.1' % (dbse)] =  9.341
BIND_NBC10A['%s-PyPy_T3-4.3' % (dbse)] =  1.991
BIND_NBC10A['%s-PyPy_T3-4.5' % (dbse)] = -1.377
BIND_NBC10A['%s-PyPy_T3-4.6' % (dbse)] = -2.203
BIND_NBC10A['%s-PyPy_T3-4.7' % (dbse)] = -2.673
BIND_NBC10A['%s-PyPy_T3-4.8' % (dbse)] = -2.896
BIND_NBC10A['%s-PyPy_T3-4.9' % (dbse)] = -2.954  # PyPy_T3 minimum
BIND_NBC10A['%s-PyPy_T3-5.0' % (dbse)] = -2.903
BIND_NBC10A['%s-PyPy_T3-5.1' % (dbse)] = -2.783
BIND_NBC10A['%s-PyPy_T3-5.2' % (dbse)] = -2.625
BIND_NBC10A['%s-PyPy_T3-5.3' % (dbse)] = -2.447
BIND_NBC10A['%s-PyPy_T3-5.4' % (dbse)] = -2.262
BIND_NBC10A['%s-PyPy_T3-5.5' % (dbse)] = -2.080
BIND_NBC10A['%s-PyPy_T3-5.7' % (dbse)] = -1.741
BIND_NBC10A['%s-PyPy_T3-6.0' % (dbse)] = -1.323
BIND_NBC10A['%s-PyPy_T3-6.5' % (dbse)] = -0.852
BIND_NBC10A['%s-PyPy_T3-7.0' % (dbse)] = -0.573
BIND_NBC10A['%s-PyPy_T3-8.0' % (dbse)] = -0.296
BIND_NBC10A['%s-PyPy_T3-9.0' % (dbse)] = -0.174

BIND_NBC10A['%s-BzBz_PD32-0.2' % (dbse)] =  3.241
BIND_NBC10A['%s-BzBz_PD32-0.4' % (dbse)] =  2.619
BIND_NBC10A['%s-BzBz_PD32-0.6' % (dbse)] =  1.726
BIND_NBC10A['%s-BzBz_PD32-0.8' % (dbse)] =  0.726
BIND_NBC10A['%s-BzBz_PD32-1.0' % (dbse)] = -0.222
BIND_NBC10A['%s-BzBz_PD32-1.2' % (dbse)] = -1.002
BIND_NBC10A['%s-BzBz_PD32-1.4' % (dbse)] = -1.553
BIND_NBC10A['%s-BzBz_PD32-1.5' % (dbse)] = -1.738
BIND_NBC10A['%s-BzBz_PD32-1.6' % (dbse)] = -1.868
BIND_NBC10A['%s-BzBz_PD32-1.7' % (dbse)] = -1.949
BIND_NBC10A['%s-BzBz_PD32-1.8' % (dbse)] = -1.988
BIND_NBC10A['%s-BzBz_PD32-1.9' % (dbse)] = -1.992  # BzBz_PD32 minimum
BIND_NBC10A['%s-BzBz_PD32-2.0' % (dbse)] = -1.971
BIND_NBC10A['%s-BzBz_PD32-2.2' % (dbse)] = -1.891
BIND_NBC10A['%s-BzBz_PD32-2.4' % (dbse)] = -1.795
BIND_NBC10A['%s-BzBz_PD32-2.6' % (dbse)] = -1.727
BIND_NBC10A['%s-BzBz_PD32-2.8' % (dbse)] = -1.702
BIND_NBC10A['%s-BzBz_PD32-3.0' % (dbse)] = -1.725

BIND_NBC10A['%s-BzBz_PD36-0.2' % (dbse)] = -1.321
BIND_NBC10A['%s-BzBz_PD36-0.4' % (dbse)] = -1.490
BIND_NBC10A['%s-BzBz_PD36-0.6' % (dbse)] = -1.735
BIND_NBC10A['%s-BzBz_PD36-0.8' % (dbse)] = -2.011
BIND_NBC10A['%s-BzBz_PD36-1.0' % (dbse)] = -2.273
BIND_NBC10A['%s-BzBz_PD36-1.2' % (dbse)] = -2.482
BIND_NBC10A['%s-BzBz_PD36-1.4' % (dbse)] = -2.619
BIND_NBC10A['%s-BzBz_PD36-1.5' % (dbse)] = -2.657
BIND_NBC10A['%s-BzBz_PD36-1.6' % (dbse)] = -2.674
BIND_NBC10A['%s-BzBz_PD36-1.7' % (dbse)] = -2.675  # BzBz_PD36 minimum
BIND_NBC10A['%s-BzBz_PD36-1.8' % (dbse)] = -2.662
BIND_NBC10A['%s-BzBz_PD36-1.9' % (dbse)] = -2.633
BIND_NBC10A['%s-BzBz_PD36-2.0' % (dbse)] = -2.593
BIND_NBC10A['%s-BzBz_PD36-2.2' % (dbse)] = -2.489
BIND_NBC10A['%s-BzBz_PD36-2.4' % (dbse)] = -2.371
BIND_NBC10A['%s-BzBz_PD36-2.6' % (dbse)] = -2.253
BIND_NBC10A['%s-BzBz_PD36-2.8' % (dbse)] = -2.143
BIND_NBC10A['%s-BzBz_PD36-3.0' % (dbse)] = -2.046
# Set default
BIND = BIND_NBC10A

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')

for item in BzBz_S:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Sandwich Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Sandwich Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Sandwich Benzene Dimer at %s A' % (distance.group(2))

for item in BzBz_T:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Benzene from T-shaped Benzene Dimer at %s A' % (distance.group(2))

for item in BzBz_PD34:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))

for item in BzH2S:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Benzene-H2S at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Benzene-H2S at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Hydrogen Sulfide from Benzene-Methane at %s A' % (distance.group(2))

for item in BzMe:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Methane from Benzene-Methane at %s A' % (distance.group(2))

for item in MeMe:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Methane Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Methane Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Methane from Methane Dimer at %s A' % (distance.group(2))

for item in PyPy_S2:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))

for item in PyPy_T3:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))

for item in BzBz_PD32:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))

for item in BzBz_PD36:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))

TAGL['%s-Bz-mono-unCP'  % (dbse)] = 'Benzene'
TAGL['%s-H2S-mono-unCP' % (dbse)] = 'Hydrogen Sulfide'
TAGL['%s-Bz2-mono-unCP' % (dbse)] = 'Benzene (alt. geometry)'
TAGL['%s-Me-mono-unCP'  % (dbse)] = 'Methane'
TAGL['%s-Py-mono-unCP'  % (dbse)] = 'Pyridine'

#<<< Geometry Specification Strings >>>
GEOS = {}

for rxn in BzBz_S:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX   1  90.0
C  3  RCC   2  90.0  1  0.0
C  3  RCC   2  90.0  1  60.0
C  3  RCC   2  90.0  1  120.0
C  3  RCC   2  90.0  1  180.0
C  3  RCC   2  90.0  1  240.0
C  3  RCC   2  90.0  1  300.0
H  3  RCH   2  90.0  1  0.0
H  3  RCH   2  90.0  1  60.0
H  3  RCH   2  90.0  1  120.0
H  3  RCH   2  90.0  1  180.0
H  3  RCH   2  90.0  1  240.0
H  3  RCH   2  90.0  1  300.0
--
0 1
X  3  RXX   2  90.0  1  0.0
X  3  R     16 90.0  2  180.0
X  3  DRXX  16 90.0  2  180.0
X  18 RXX   17 90.0  16 0.0
C  17 RCC   18 90.0  19 0.0
C  17 RCC   18 90.0  19 60.0
C  17 RCC   18 90.0  19 120.0
C  17 RCC   18 90.0  19 180.0
C  17 RCC   18 90.0  19 240.0
C  17 RCC   18 90.0  19 300.0
H  17 RCH   18 90.0  19 0.0
H  17 RCH   18 90.0  19 60.0
H  17 RCH   18 90.0  19 120.0
H  17 RCH   18 90.0  19 180.0
H  17 RCH   18 90.0  19 240.0
H  17 RCH   18 90.0  19 300.0

RXX    = 1.0
DRXX   = 12.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in BzBz_T:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX  1  90.0
C  3  RCC  2  90.0  1  0.0
C  3  RCC  2  90.0  1  60.0
C  3  RCC  2  90.0  1  120.0
C  3  RCC  2  90.0  1  180.0
C  3  RCC  2  90.0  1  240.0
C  3  RCC  2  90.0  1  300.0
H  3  RCH  2  90.0  1  0.0
H  3  RCH  2  90.0  1  60.0
H  3  RCH  2  90.0  1  120.0
H  3  RCH  2  90.0  1  180.0
H  3  RCH  2  90.0  1  240.0
H  3  RCH  2  90.0  1  300.0
--
0 1
X  3  RXX  2  90.0  1  0.0
X  3  R    16 90.0  1  0.0
X  17 RXX  3  90.0  16 180.0
X  18 RXX  17 90.0  3  0.0
C  17 RCC  18 90.0  19 0.0
C  17 RCC  18 90.0  19 60.0
C  17 RCC  18 90.0  19 120.0
C  17 RCC  18 90.0  19 180.0
C  17 RCC  18 90.0  19 240.0
C  17 RCC  18 90.0  19 300.0
H  17 RCH  18 90.0  19 0.0
H  17 RCH  18 90.0  19 60.0
H  17 RCH  18 90.0  19 120.0
H  17 RCH  18 90.0  19 180.0
H  17 RCH  18 90.0  19 240.0
H  17 RCH  18 90.0  19 300.0

RXX    = 1.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in sum([BzBz_PD32, BzBz_PD34, BzBz_PD36], []):
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))
    if rxn in BzBz_PD32:
        R2val = 3.2
    elif rxn in BzBz_PD34:
        R2val = 3.4
    elif rxn in BzBz_PD36:
        R2val = 3.6

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  R2   1  90.0
C  3  RCC  2  90.0  1  0.0
C  3  RCC  2  90.0  1  60.0
C  3  RCC  2  90.0  1  120.0
C  3  RCC  2  90.0  1  180.0
C  3  RCC  2  90.0  1  240.0
C  3  RCC  2  90.0  1  300.0
H  3  RCH  2  90.0  1  0.0
H  3  RCH  2  90.0  1  60.0
H  3  RCH  2  90.0  1  120.0
H  3  RCH  2  90.0  1  180.0
H  3  RCH  2  90.0  1  240.0
H  3  RCH  2  90.0  1  300.0
--
0 1
X  3  RXX  2  90.0  1  0.0
X  2  R    3  90.0  16 90.0
X  17 RXX  2  90.0  1  90.0
X  18 RXX  17 90.0  2  90.0
C  17 RCC  18 90.0  19 0.0
C  17 RCC  18 90.0  19 60.0
C  17 RCC  18 90.0  19 120.0
C  17 RCC  18 90.0  19 180.0
C  17 RCC  18 90.0  19 240.0
C  17 RCC  18 90.0  19 300.0
H  17 RCH  18 90.0  19 0.0
H  17 RCH  18 90.0  19 60.0
H  17 RCH  18 90.0  19 120.0
H  17 RCH  18 90.0  19 180.0
H  17 RCH  18 90.0  19 240.0
H  17 RCH  18 90.0  19 300.0

RXX    = 1.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
R2     = %(R2val)s
units angstrom
""" % vars())

for rxn in BzH2S:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  BZCX  1 90.0
C  2  BZCX  1 90.0   3 60.0
C  2  BZCX  1 90.0   4 60.0
C  2  BZCX  1 90.0   5 60.0
C  2  BZCX  1 90.0   6 60.0
C  2  BZCX  1 90.0   7 60.0
X  3  1.0   2 90.0   1 0.0
H  3  BZHC  9 90.0   2 180.0
H  4  BZHC  3 120.0  2 180.0
H  5  BZHC  4 120.0  2 180.0
H  6  BZHC  5 120.0  2 180.0
H  7  BZHC  6 120.0  2 180.0
H  8  BZHC  7 120.0  2 180.0
--
0 1
S  2  R     3 90.0   4 90.0
H  16 HS    2 HHSH   9 180.0
H  16 HS    2 HHSH   9 0.0

BZCX   = 1.3915
BZHC   = 1.0800
HS     = 1.3356
HHSH   = 46.06
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in BzMe:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  CQ   1  90.0
C  3  CQ   2  60.0   1  90.0
C  4  CQ   2  60.0   1  90.0
C  5  CQ   2  60.0   1  90.0
C  6  CQ   2  60.0   1  90.0
C  7  CQ   2  60.0   1  90.0
X  3  1.0  2  90.0   1  0.0
H  3  CH1  9  90.0   2  180.0
H  4  CH1  3  120.0  2  180.0
H  5  CH1  4  120.0  2  180.0
H  6  CH1  5  120.0  2  180.0
H  7  CH1  6  120.0  2  180.0
H  8  CH1  7  120.0  2  180.0
--
0 1
C  2  R    3  90.0   9  0.0
H  16 CH2  2  0.0    3  0.0
H  16 CH2  2  HCH    3  0.0
H  16 CH2  17 HCH    18 120.0
H  16 CH2  17 HCH    18 240.0

CQ     = 1.405731
CH1    = 1.095210
CH2    = 1.099503
HCH    = 109.471209
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in MeMe:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
C
H  1 CH2
H  1 CH2  2 HCH
H  1 CH2  2 HCH    3  120.0
H  1 CH2  2 HCH    3  240.0
--
0 1
C  1 R    2 180.0  4  120.0
H  6 CH2  2 180.0  4  120.0
H  6 CH2  7 HCH    3  180.0
H  6 CH2  7 HCH    4  180.0
H  6 CH2  7 HCH    5  180.0

CH2    = 1.099503
HCH    = 109.471209
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in PyPy_S2:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  R
N  1  1.3980380  2  90.0
C  1  1.3371053  2  90.0    3  -58.504950
C  1  1.3822904  2  90.0    4  -61.640500
C  1  1.4067471  2  90.0    5  -59.854550
C  1  1.3822904  2  90.0    6  -59.854550
C  1  1.3371053  2  90.0    7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4  0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0
--
0 1
N  2  1.3980380  1  90.0    3  theta
C  2  1.3371053  1  90.0    14 -58.504950
C  2  1.3822904  1  90.0    15 -61.640500
C  2  1.4067471  1  90.0    16 -59.854550
C  2  1.3822904  1  90.0    17 -59.854550
C  2  1.3371053  1  90.0    18 -61.640500
H  15 1.08650    14 116.01  19 180.0
H  16 1.08260    15 120.12  14 180.0
H  17 1.08180    14 180.00  15 0.0
H  18 1.08260    19 120.12  14 180.0
H  19 1.08650    14 116.01  15 180.0

theta  = 180.0
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in PyPy_T3:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  R
N  1  1.3980380  2   90.0
C  1  1.3371053  2   90.0   3  -58.504950
C  1  1.3822904  2   90.0   4  -61.640500
C  1  1.4067471  2   90.0   5  -59.854550
C  1  1.3822904  2   90.0   6  -59.854550
C  1  1.3371053  2   90.0   7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4    0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0
--
0 1
X  2  2.0000000  1   90.0   3  theta
N  2  1.3980380  14  90.0   1  updown
C  2  1.3371053  14  90.0   15 -58.504950
C  2  1.3822904  14  90.0   16 -61.640500
C  2  1.4067471  14  90.0   17 -59.854550
C  2  1.3822904  14  90.0   18 -59.854550
C  2  1.3371053  14  90.0   19 -61.640500
H  16 1.08650    15 116.01  20 180.0
H  17 1.08260    16 120.12  15 180.0
H  18 1.08180    15 180.00  16 0.0
H  19 1.08260    20 120.12  15 180.0
H  20 1.08650    15 116.01  16 180.0

theta  = 90.0
updown = 270.0
R      = %(Rval)s
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Bz', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX 1  90.0
C  3  RCC 2  90.0 1  0.0
C  3  RCC 2  90.0 1  60.0
C  3  RCC 2  90.0 1  120.0
C  3  RCC 2  90.0 1  180.0
C  3  RCC 2  90.0 1  240.0
C  3  RCC 2  90.0 1  300.0
H  3  RCH 2  90.0 1  0.0
H  3  RCH 2  90.0 1  60.0
H  3  RCH 2  90.0 1  120.0
H  3  RCH 2  90.0 1  180.0
H  3  RCH 2  90.0 1  240.0
H  3  RCH 2  90.0 1  300.0

RXX = 1.0
RCC = 1.3915
RCH = 2.4715
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'H2S', 'mono-unCP')] = qcdb.Molecule("""
0 1
S
H  1 HS
H  1 HS  2 HSH

HS  = 1.3356
HSH = 92.12
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Bz2', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  CQ   1  90.0
C  3  CQ   2  60.0  1  90.0
C  4  CQ   2  60.0  1  90.0
C  5  CQ   2  60.0  1  90.0
C  6  CQ   2  60.0  1  90.0
C  7  CQ   2  60.0  1  90.0
X  3  1.0  2  90.0  1  0.0
H  3  CH1  9  90.0  2  180.0
H  4  CH1  3  120.0 2  180.0
H  5  CH1  4  120.0 2  180.0
H  6  CH1  5  120.0 2  180.0
H  7  CH1  6  120.0 2  180.0
H  8  CH1  7  120.0 2  180.0

CQ  = 1.405731
CH1 = 1.095210
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Me', 'mono-unCP')] = qcdb.Molecule("""
0 1
C
H  1  CH2
H  1  CH2  2  HCH
H  1  CH2  2  HCH   3  120.0
H  1  CH2  2  HCH   3  240.0

CH2 = 1.099503
HCH = 109.471209
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Py', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  RXX
N  1  1.3980380  2  90.0
C  1  1.3371053  2  90.0    3  -58.504950
C  1  1.3822904  2  90.0    4  -61.640500
C  1  1.4067471  2  90.0    5  -59.854550
C  1  1.3822904  2  90.0    6  -59.854550
C  1  1.3371053  2  90.0    7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4  0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0

RXX = 1.0
units angstrom
""" % vars())

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)

#########################################################################

# <<< Supplementary Quantum Chemical Results >>>
DATA = {}

DATA['NUCLEAR REPULSION ENERGY'] = {}
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-dimer'          ] =     652.58240326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Bz-mono-unCP'              ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-dimer'          ] =     647.08083072
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-dimer'          ] =     641.79881504
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-dimer'          ] =     636.72435401
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-dimer'          ] =     631.84627841
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-dimer'          ] =     627.15417831
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-dimer'          ] =     622.63833806
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-dimer'          ] =     618.28967853
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-dimer'          ] =     614.09970566
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-dimer'          ] =     610.06046424
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-dimer'          ] =     606.16449631
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-dimer'          ] =     595.26834684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-dimer'          ] =     579.39688238
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-dimer'          ] =     565.87021271
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-dimer'          ] =     554.22625379
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-dimer'          ] =     544.11253672
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-dimer'         ] =     499.16037479
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-dimer'          ] =     613.04854518
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-dimer'          ] =     608.81636557
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-dimer'          ] =     604.74550671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-dimer'          ] =     600.82787505
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-dimer'          ] =     597.05577907
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-dimer'          ] =     593.42192782
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-dimer'          ] =     589.91942332
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-dimer'          ] =     586.54174882
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-dimer'          ] =     583.28275414
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-dimer'          ] =     580.13663931
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-dimer'          ] =     577.09793714
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-dimer'          ] =     574.16149552
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-dimer'          ] =     571.32245963
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-dimer'          ] =     560.85272572
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-dimer'          ] =     549.47925556
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-dimer'          ] =     539.65622514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-dimer'          ] =     531.09189940
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-dimer'          ] =     523.56205991
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-dimer'       ] =     641.59153721
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-dimer'       ] =     640.97218086
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-dimer'       ] =     639.94808010
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-dimer'       ] =     638.53114770
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-dimer'       ] =     636.73745247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-dimer'       ] =     634.58670201
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-dimer'       ] =     632.10168144
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-dimer'       ] =     630.74164257
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-dimer'       ] =     629.30768985
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-dimer'       ] =     627.80329032
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-dimer'       ] =     626.23200316
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-dimer'       ] =     624.59746513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-dimer'       ] =     622.90337667
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-dimer'       ] =     619.35158842
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-dimer'       ] =     615.60701452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-dimer'       ] =     611.70022314
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-dimer'       ] =     607.66157487
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-dimer'       ] =     603.52082284
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-dimer'           ] =     332.50866690
DATA['NUCLEAR REPULSION ENERGY']['NBC1-H2S-mono-unCP'             ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-dimer'           ] =     326.76493049
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-dimer'           ] =     324.08312886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-dimer'           ] =     321.51823084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-dimer'           ] =     319.06348175
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-dimer'           ] =     316.71257239
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-dimer'           ] =     314.45961051
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-dimer'           ] =     312.29909326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-dimer'           ] =     310.22588084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-dimer'           ] =     308.23517159
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-dimer'           ] =     302.71463310
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-dimer'          ] =     298.57449040
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-dimer'           ] =     294.79763877
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-dimer'          ] =     291.34045574
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-dimer'           ] =     288.16568982
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-dimer'           ] =     282.54011405
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-dimer'           ] =     277.71464354
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-dimer'           ] =     273.53417452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-dimer'           ] =     269.88029141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-dimer'            ] =     277.70122037
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Bz2-mono-unCP'             ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Me-mono-unCP'              ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-dimer'            ] =     276.14505886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-dimer'            ] =     274.65657480
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-dimer'            ] =     273.23211647
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-dimer'            ] =     271.86820659
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-dimer'            ] =     270.56154682
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-dimer'            ] =     269.30901798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-dimer'            ] =     268.10767718
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-dimer'            ] =     266.95475267
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-dimer'            ] =     265.84763738
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-dimer'            ] =     264.78388141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-dimer'            ] =     262.77738579
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-dimer'            ] =     260.91850385
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-dimer'            ] =     259.19247204
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-dimer'            ] =     257.58628148
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-dimer'            ] =     256.08845607
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-dimer'            ] =     254.68885527
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-dimer'            ] =     253.37850109
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-dimer'            ] =     250.99455064
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-dimer'            ] =      42.94051671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-dimer'            ] =      42.46449704
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-dimer'            ] =      42.01471911
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-dimer'            ] =      41.58914043
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-dimer'            ] =      41.18591734
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-dimer'            ] =      40.80338247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-dimer'            ] =      40.44002498
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-dimer'            ] =      40.09447330
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-dimer'            ] =      39.76547998
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-dimer'            ] =      39.45190844
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-dimer'            ] =      39.15272123
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-dimer'            ] =      38.86696980
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-dimer'            ] =      38.59378540
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-dimer'            ] =      38.08199453
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-dimer'            ] =      37.61171219
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-dimer'            ] =      37.17815187
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-dimer'            ] =      36.40542136
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-dimer'            ] =      35.73746090
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-dimer'         ] =     664.74968142
DATA['NUCLEAR REPULSION ENERGY']['NBC1-Py-mono-unCP'              ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-dimer'         ] =     653.28897360
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-dimer'         ] =     647.90584891
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-dimer'         ] =     642.73711461
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-dimer'         ] =     637.77107423
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-dimer'         ] =     632.99683541
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-dimer'         ] =     628.40424073
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-dimer'         ] =     623.98380628
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-dimer'         ] =     619.72666684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-dimer'         ] =     615.62452662
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-dimer'         ] =     611.66961499
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-dimer'         ] =     607.85464633
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-dimer'         ] =     604.17278378
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-dimer'         ] =     600.61760611
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-dimer'         ] =     593.86352067
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-dimer'         ] =     584.54275675
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-dimer'         ] =     570.86466240
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-dimer'         ] =     559.10620798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-dimer'         ] =     548.90465922
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-dimer'         ] =     539.98032943
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-dimer'         ] =     631.74018099
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-dimer'         ] =     622.28221702
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-dimer'         ] =     613.57422251
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-dimer'         ] =     609.47520868
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-dimer'         ] =     605.53368830
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-dimer'         ] =     601.74111111
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-dimer'         ] =     598.08951503
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-dimer'         ] =     594.57147649
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-dimer'         ] =     591.18006603
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-dimer'         ] =     587.90880856
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-dimer'         ] =     584.75164753
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-dimer'         ] =     581.70291245
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-dimer'         ] =     578.75728949
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-dimer'         ] =     573.15574951
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-dimer'         ] =     565.41165299
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-dimer'         ] =     554.01089095
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-dimer'         ] =     544.16644693
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-dimer'         ] =     528.04095562
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-dimer'         ] =     515.40150653
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-dimer'       ] =     652.35026383
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-dimer'       ] =     651.65685475
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-dimer'       ] =     650.51106101
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-dimer'       ] =     648.92723975
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-dimer'       ] =     646.92462020
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-dimer'       ] =     644.52659143
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-dimer'       ] =     641.75995892
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-dimer'       ] =     640.24755050
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-dimer'       ] =     638.65423207
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-dimer'       ] =     636.98400901
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-dimer'       ] =     635.24097954
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-dimer'       ] =     633.42931896
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-dimer'       ] =     631.55326486
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-dimer'       ] =     627.62515488
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-dimer'       ] =     623.49127864
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-dimer'       ] =     619.18640729
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-dimer'       ] =     614.74502815
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-dimer'       ] =     610.20089775
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-dimer'       ] =     631.66053374
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-dimer'       ] =     631.10536715
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-dimer'       ] =     630.18691177
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-dimer'       ] =     628.91516711
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-dimer'       ] =     627.30369102
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-dimer'       ] =     625.36921338
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-dimer'       ] =     623.13120361
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-dimer'       ] =     621.90509666
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-dimer'       ] =     620.61142042
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-dimer'       ] =     619.25317914
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-dimer'       ] =     617.83346514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-dimer'       ] =     616.35544587
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-dimer'       ] =     614.82235130
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-dimer'       ] =     611.60409513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-dimer'       ] =     608.20532569
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-dimer'       ] =     604.65291019
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-dimer'       ] =     600.97358989
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-dimer'       ] =     597.19362514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-dimer'          ] =     652.58240326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-dimer'          ] =     647.08083072
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.3-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-dimer'          ] =     641.79881504
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-dimer'          ] =     636.72435401
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-dimer'          ] =     631.84627841
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-dimer'          ] =     627.15417831
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.7-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-dimer'          ] =     622.63833806
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.8-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-dimer'          ] =     618.28967853
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-3.9-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-dimer'          ] =     614.09970566
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-dimer'          ] =     610.06046424
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.1-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-dimer'          ] =     606.16449631
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-dimer'          ] =     595.26834684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-4.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-dimer'          ] =     579.39688238
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-dimer'          ] =     565.87021271
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-5.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-dimer'          ] =     554.22625379
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-dimer'          ] =     544.11253672
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-6.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-dimer'         ] =     499.16037479
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_S-10.0-monoA-CP'      ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-dimer'          ] =     613.04854518
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.4-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-dimer'          ] =     608.81636557
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-dimer'          ] =     604.74550671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.6-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-dimer'          ] =     600.82787505
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.7-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-dimer'          ] =     597.05577907
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.8-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-dimer'          ] =     593.42192782
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-4.9-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-dimer'          ] =     589.91942332
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-dimer'          ] =     586.54174882
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.1-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-dimer'          ] =     583.28275414
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.2-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-dimer'          ] =     580.13663931
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.3-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-dimer'          ] =     577.09793714
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.4-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-dimer'          ] =     574.16149552
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-dimer'          ] =     571.32245963
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-5.6-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-dimer'          ] =     560.85272572
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-dimer'          ] =     549.47925556
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-6.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-dimer'          ] =     539.65622514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-dimer'          ] =     531.09189940
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-7.5-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-dimer'          ] =     523.56205991
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_T-8.0-monoB-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-dimer'       ] =     641.59153721
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-dimer'       ] =     640.97218086
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-dimer'       ] =     639.94808010
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-dimer'       ] =     638.53114770
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-dimer'       ] =     636.73745247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-dimer'       ] =     634.58670201
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-dimer'       ] =     632.10168144
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-dimer'       ] =     630.74164257
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-dimer'       ] =     629.30768985
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-dimer'       ] =     627.80329032
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-dimer'       ] =     626.23200316
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-dimer'       ] =     624.59746513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-dimer'       ] =     622.90337667
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-dimer'       ] =     619.35158842
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-dimer'       ] =     615.60701452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-dimer'       ] =     611.70022314
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-dimer'       ] =     607.66157487
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-dimer'       ] =     603.52082284
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD34-3.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-dimer'           ] =     332.50866690
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.2-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-dimer'           ] =     326.76493049
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.4-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-dimer'           ] =     324.08312886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-dimer'           ] =     321.51823084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.6-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-dimer'           ] =     319.06348175
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.7-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-dimer'           ] =     316.71257239
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.8-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-dimer'           ] =     314.45961051
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-3.9-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-dimer'           ] =     312.29909326
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-dimer'           ] =     310.22588084
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.1-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-dimer'           ] =     308.23517159
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.2-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-dimer'           ] =     302.71463310
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-dimer'          ] =     298.57449040
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-4.75-monoB-CP'       ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-dimer'           ] =     294.79763877
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-dimer'          ] =     291.34045574
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-monoA-CP'       ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.25-monoB-CP'       ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-dimer'           ] =     288.16568982
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-5.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-dimer'           ] =     282.54011405
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-dimer'           ] =     277.71464354
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-6.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-dimer'           ] =     273.53417452
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.0-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-dimer'           ] =     269.88029141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-monoA-CP'        ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzH2S-7.5-monoB-CP'        ] =      12.95382185
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-dimer'            ] =     277.70122037
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-dimer'            ] =     276.14505886
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.3-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-dimer'            ] =     274.65657480
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-dimer'            ] =     273.23211647
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.5-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-dimer'            ] =     271.86820659
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-dimer'            ] =     270.56154682
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.7-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-dimer'            ] =     269.30901798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.8-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-dimer'            ] =     268.10767718
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-3.9-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-dimer'            ] =     266.95475267
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-dimer'            ] =     265.84763738
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.1-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-dimer'            ] =     264.78388141
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-dimer'            ] =     262.77738579
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-dimer'            ] =     260.91850385
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-dimer'            ] =     259.19247204
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-4.8-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-dimer'            ] =     257.58628148
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-dimer'            ] =     256.08845607
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.2-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-dimer'            ] =     254.68885527
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.4-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-dimer'            ] =     253.37850109
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-5.6-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-dimer'            ] =     250.99455064
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-monoA-CP'         ] =     201.83853774
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzMe-6.0-monoB-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-dimer'            ] =      42.94051671
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.2-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-dimer'            ] =      42.46449704
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.3-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-dimer'            ] =      42.01471911
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-dimer'            ] =      41.58914043
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.5-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-dimer'            ] =      41.18591734
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.6-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-dimer'            ] =      40.80338247
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.7-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-dimer'            ] =      40.44002498
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-dimer'            ] =      40.09447330
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-3.9-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-dimer'            ] =      39.76547998
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.0-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-dimer'            ] =      39.45190844
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.1-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-dimer'            ] =      39.15272123
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.2-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-dimer'            ] =      38.86696980
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.3-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-dimer'            ] =      38.59378540
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-dimer'            ] =      38.08199453
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.6-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-dimer'            ] =      37.61171219
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-4.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-dimer'            ] =      37.17815187
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.0-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-dimer'            ] =      36.40542136
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.4-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-dimer'            ] =      35.73746090
DATA['NUCLEAR REPULSION ENERGY']['NBC1-MeMe-5.8-monoA-CP'         ] =      13.31926457
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-dimer'         ] =     664.74968142
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-dimer'         ] =     653.28897360
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-dimer'         ] =     647.90584891
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-dimer'         ] =     642.73711461
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-dimer'         ] =     637.77107423
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.6-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-dimer'         ] =     632.99683541
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-dimer'         ] =     628.40424073
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.8-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-dimer'         ] =     623.98380628
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-3.9-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-dimer'         ] =     619.72666684
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-dimer'         ] =     615.62452662
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-dimer'         ] =     611.66961499
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.2-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-dimer'         ] =     607.85464633
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-dimer'         ] =     604.17278378
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-dimer'         ] =     600.61760611
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-dimer'         ] =     593.86352067
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-4.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-dimer'         ] =     584.54275675
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-dimer'         ] =     570.86466240
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-5.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-dimer'         ] =     559.10620798
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-dimer'         ] =     548.90465922
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-6.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-dimer'         ] =     539.98032943
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_S2-7.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-dimer'         ] =     631.74018099
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.1-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-dimer'         ] =     622.28221702
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.3-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-dimer'         ] =     613.57422251
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-dimer'         ] =     609.47520868
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.6-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-dimer'         ] =     605.53368830
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.7-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-dimer'         ] =     601.74111111
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.8-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-dimer'         ] =     598.08951503
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-4.9-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-dimer'         ] =     594.57147649
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-dimer'         ] =     591.18006603
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.1-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-dimer'         ] =     587.90880856
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.2-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-dimer'         ] =     584.75164753
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.3-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-dimer'         ] =     581.70291245
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.4-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-dimer'         ] =     578.75728949
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-dimer'         ] =     573.15574951
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-5.7-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-dimer'         ] =     565.41165299
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-dimer'         ] =     554.01089095
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-6.5-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-dimer'         ] =     544.16644693
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-7.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-dimer'         ] =     528.04095562
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-8.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-dimer'         ] =     515.40150653
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-monoA-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-PyPy_T3-9.0-monoB-CP'      ] =     206.21910131
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-dimer'       ] =     652.35026383
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-dimer'       ] =     651.65685475
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-dimer'       ] =     650.51106101
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-dimer'       ] =     648.92723975
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-dimer'       ] =     646.92462020
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-dimer'       ] =     644.52659143
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-dimer'       ] =     641.75995892
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-dimer'       ] =     640.24755050
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-dimer'       ] =     638.65423207
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-dimer'       ] =     636.98400901
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-dimer'       ] =     635.24097954
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-dimer'       ] =     633.42931896
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-dimer'       ] =     631.55326486
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-dimer'       ] =     627.62515488
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-dimer'       ] =     623.49127864
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-dimer'       ] =     619.18640729
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-dimer'       ] =     614.74502815
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-dimer'       ] =     610.20089775
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD32-3.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-dimer'       ] =     631.66053374
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-dimer'       ] =     631.10536715
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-dimer'       ] =     630.18691177
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-dimer'       ] =     628.91516711
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-0.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-dimer'       ] =     627.30369102
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-dimer'       ] =     625.36921338
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-dimer'       ] =     623.13120361
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-dimer'       ] =     621.90509666
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.5-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-dimer'       ] =     620.61142042
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-dimer'       ] =     619.25317914
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.7-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-dimer'       ] =     617.83346514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-dimer'       ] =     616.35544587
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-1.9-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-dimer'       ] =     614.82235130
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.0-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-dimer'       ] =     611.60409513
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.2-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-dimer'       ] =     608.20532569
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.4-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-dimer'       ] =     604.65291019
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.6-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-dimer'       ] =     600.97358989
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-2.8-monoA-CP'    ] =     204.01997321
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-dimer'       ] =     597.19362514
DATA['NUCLEAR REPULSION ENERGY']['NBC1-BzBz_PD36-3.0-monoA-CP'    ] =     204.01997321
