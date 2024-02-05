#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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
| Database (Sherrill) of interaction energies for dissociation curves of rare-gas biatomic complexes.
| Geometries and reference interaction energies from Tang et al. JCP 118 4976 (2003).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'HeHe'`` 18-point dissociation curve for helium dimer
  - ``'HeNe'`` 18-point dissociation curve for helium-neon complex
  - ``'HeAr'`` 18-point dissociation curve for helium-argon complex
  - ``'HeKr'`` 18-point dissociation curve for helium-krypton complex
  - ``'NeNe'`` 18-point dissociation curve for neon dimer
  - ``'NeAr'`` 18-point dissociation curve for neon-argon complex
  - ``'NeKr'`` 18-point dissociation curve for neon-krypton complex
  - ``'ArAr'`` 18-point dissociation curve for argon dimer
  - ``'ArKr'`` 18-point dissociation curve for argon-krypton complex
  - ``'KrKr'`` 18-point dissociation curve for krypton dimer

"""
import re

import qcdb

# <<< RGC10 Database Module >>>
dbse = 'RGC1'

# <<< Database Members >>>
HeHe = []
HeNe = []
HeAr = []
HeKr = []
NeNe = []
NeAr = []
NeKr = []
ArAr = []
ArKr = []
KrKr = []
dist = [0.85, 0.9, 0.95, 0.975, 1.0, 1.025, 1.05, 1.1, 1.15,
        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2]
for d in dist:
    HeHe.append('HeHe-' + str(d))
    HeNe.append('HeNe-' + str(d))
    HeAr.append('HeAr-' + str(d))
    HeKr.append('HeKr-' + str(d))
    NeNe.append('NeNe-' + str(d))
    NeAr.append('NeAr-' + str(d))
    NeKr.append('NeKr-' + str(d))
    ArAr.append('ArAr-' + str(d))
    ArKr.append('ArKr-' + str(d))
    KrKr.append('KrKr-' + str(d))

temp = [HeHe, HeNe, HeAr, HeKr, NeNe, NeAr, NeKr, ArAr, ArKr, KrKr]
HRXN = sum(temp, [])

HRXN_SM = ['NeNe-1.0', 'NeNe-1.1', 'NeAr-0.85']
HRXN_LG = ['KrKr-0.85']
HRXN_EQ = ['HeHe-1.0', 'HeNe-1.0', 'HeAr-1.0', 'HeKr-1.0', 'NeNe-1.0',
           'NeAr-1.0', 'NeKr-1.0', 'ArAr-1.0', 'ArKr-1.0', 'KrKr-1.0']

Req = {}
Req['HeHe'] = 2.98
Req['HeNe'] = 3.05
Req['HeAr'] = 3.50
Req['HeKr'] = 3.70
Req['NeNe'] = 3.09
Req['NeAr'] = 3.48
Req['NeKr'] = 3.65
Req['ArAr'] = 3.75
Req['ArKr'] = 3.89
Req['KrKr'] = 4.01

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations
for rxn in HRXN:

    if rxn in HeHe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-He-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-He-mono-unCP'  % (dbse) ]

    elif rxn in HeNe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-He-mono-unCP'  % (dbse)      : -1,
                                          '%s-Ne-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-He-mono-unCP'  % (dbse),
                                          '%s-Ne-mono-unCP'  % (dbse) ]

    elif rxn in HeAr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-He-mono-unCP'  % (dbse)      : -1,
                                          '%s-Ar-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-He-mono-unCP'  % (dbse),
                                          '%s-Ar-mono-unCP'  % (dbse) ]

    elif rxn in HeKr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-He-mono-unCP'  % (dbse)      : -1,
                                          '%s-Kr-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-He-mono-unCP'  % (dbse),
                                          '%s-Kr-mono-unCP'  % (dbse) ]

    elif rxn in NeNe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Ne-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Ne-mono-unCP'  % (dbse) ]

    elif rxn in NeAr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Ne-mono-unCP'  % (dbse)      : -1,
                                          '%s-Ar-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Ne-mono-unCP'  % (dbse),
                                          '%s-Ar-mono-unCP'  % (dbse) ]

    elif rxn in NeKr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Ne-mono-unCP'  % (dbse)      : -1,
                                          '%s-Kr-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Ne-mono-unCP'  % (dbse),
                                          '%s-Kr-mono-unCP'  % (dbse) ]

    elif rxn in ArAr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Ar-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Ar-mono-unCP'  % (dbse) ]

    elif rxn in ArKr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Ar-mono-unCP'  % (dbse)      : -1,
                                          '%s-Kr-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Ar-mono-unCP'  % (dbse),
                                          '%s-Kr-mono-unCP'  % (dbse) ]

    elif rxn in KrKr:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Kr-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Kr-mono-unCP'  % (dbse) ]

# <<< Reference Values >>>
BIND = {}
BIND['%s-HeHe-0.85'  % (dbse)] =   0.03759
BIND['%s-HeHe-0.9'   % (dbse)] =  -0.00449
BIND['%s-HeHe-0.95'  % (dbse)] =  -0.01905
BIND['%s-HeHe-0.975' % (dbse)] =  -0.02135
BIND['%s-HeHe-1.0'   % (dbse)] =  -0.02188  # HeHe minimum
BIND['%s-HeHe-1.025' % (dbse)] =  -0.02133
BIND['%s-HeHe-1.05'  % (dbse)] =  -0.02017
BIND['%s-HeHe-1.1'   % (dbse)] =  -0.01708
BIND['%s-HeHe-1.15'  % (dbse)] =  -0.01392
BIND['%s-HeHe-1.2'   % (dbse)] =  -0.01113
BIND['%s-HeHe-1.3'   % (dbse)] =  -0.00702
BIND['%s-HeHe-1.4'   % (dbse)] =  -0.00447
BIND['%s-HeHe-1.5'   % (dbse)] =  -0.00291
BIND['%s-HeHe-1.6'   % (dbse)] =  -0.00195
BIND['%s-HeHe-1.7'   % (dbse)] =  -0.00133
BIND['%s-HeHe-1.8'   % (dbse)] =  -0.00093
BIND['%s-HeHe-2.0'   % (dbse)] =  -0.00049
BIND['%s-HeHe-2.2'   % (dbse)] =  -0.00027

BIND['%s-HeNe-0.85'  % (dbse)] =   0.08105
BIND['%s-HeNe-0.9'   % (dbse)] =  -0.00535
BIND['%s-HeNe-0.95'  % (dbse)] =  -0.03530
BIND['%s-HeNe-0.975' % (dbse)] =  -0.04012
BIND['%s-HeNe-1.0'   % (dbse)] =  -0.04136  # HeNe minimum
BIND['%s-HeNe-1.025' % (dbse)] =  -0.04043
BIND['%s-HeNe-1.05'  % (dbse)] =  -0.03825
BIND['%s-HeNe-1.1'   % (dbse)] =  -0.03236
BIND['%s-HeNe-1.15'  % (dbse)] =  -0.02629
BIND['%s-HeNe-1.2'   % (dbse)] =  -0.02097
BIND['%s-HeNe-1.3'   % (dbse)] =  -0.01315
BIND['%s-HeNe-1.4'   % (dbse)] =  -0.00832
BIND['%s-HeNe-1.5'   % (dbse)] =  -0.00540
BIND['%s-HeNe-1.6'   % (dbse)] =  -0.00359
BIND['%s-HeNe-1.7'   % (dbse)] =  -0.00246
BIND['%s-HeNe-1.8'   % (dbse)] =  -0.00172
BIND['%s-HeNe-2.0'   % (dbse)] =  -0.00089
BIND['%s-HeNe-2.2'   % (dbse)] =  -0.00049

BIND['%s-HeAr-0.85'  % (dbse)] =   0.11196
BIND['%s-HeAr-0.9'   % (dbse)] =  -0.00862
BIND['%s-HeAr-0.95'  % (dbse)] =  -0.05048
BIND['%s-HeAr-0.975' % (dbse)] =  -0.05720
BIND['%s-HeAr-1.0'   % (dbse)] =  -0.05889  # HeAr minimum
BIND['%s-HeAr-1.025' % (dbse)] =  -0.05752
BIND['%s-HeAr-1.05'  % (dbse)] =  -0.05440
BIND['%s-HeAr-1.1'   % (dbse)] =  -0.04600
BIND['%s-HeAr-1.15'  % (dbse)] =  -0.03735
BIND['%s-HeAr-1.2'   % (dbse)] =  -0.02977
BIND['%s-HeAr-1.3'   % (dbse)] =  -0.01862
BIND['%s-HeAr-1.4'   % (dbse)] =  -0.01176
BIND['%s-HeAr-1.5'   % (dbse)] =  -0.00760
BIND['%s-HeAr-1.6'   % (dbse)] =  -0.00505
BIND['%s-HeAr-1.7'   % (dbse)] =  -0.00344
BIND['%s-HeAr-1.8'   % (dbse)] =  -0.00240
BIND['%s-HeAr-2.0'   % (dbse)] =  -0.00124
BIND['%s-HeAr-2.2'   % (dbse)] =  -0.00069

BIND['%s-HeKr-0.85'  % (dbse)] =   0.11043
BIND['%s-HeKr-0.9'   % (dbse)] =  -0.01063
BIND['%s-HeKr-0.95'  % (dbse)] =  -0.05251
BIND['%s-HeKr-0.975' % (dbse)] =  -0.05914
BIND['%s-HeKr-1.0'   % (dbse)] =  -0.06071  # HeKr minimum
BIND['%s-HeKr-1.025' % (dbse)] =  -0.05919
BIND['%s-HeKr-1.05'  % (dbse)] =  -0.05592
BIND['%s-HeKr-1.1'   % (dbse)] =  -0.04721
BIND['%s-HeKr-1.15'  % (dbse)] =  -0.03830
BIND['%s-HeKr-1.2'   % (dbse)] =  -0.03050
BIND['%s-HeKr-1.3'   % (dbse)] =  -0.01904
BIND['%s-HeKr-1.4'   % (dbse)] =  -0.01201
BIND['%s-HeKr-1.5'   % (dbse)] =  -0.00775
BIND['%s-HeKr-1.6'   % (dbse)] =  -0.00514
BIND['%s-HeKr-1.7'   % (dbse)] =  -0.00350
BIND['%s-HeKr-1.8'   % (dbse)] =  -0.00244
BIND['%s-HeKr-2.0'   % (dbse)] =  -0.00126
BIND['%s-HeKr-2.2'   % (dbse)] =  -0.00070

BIND['%s-NeNe-0.85'  % (dbse)] =   0.16931
BIND['%s-NeNe-0.9'   % (dbse)] =  -0.00949
BIND['%s-NeNe-0.95'  % (dbse)] =  -0.07154
BIND['%s-NeNe-0.975' % (dbse)] =  -0.08158
BIND['%s-NeNe-1.0'   % (dbse)] =  -0.08420  # NeNe minimum
BIND['%s-NeNe-1.025' % (dbse)] =  -0.08233
BIND['%s-NeNe-1.05'  % (dbse)] =  -0.07789
BIND['%s-NeNe-1.1'   % (dbse)] =  -0.06582
BIND['%s-NeNe-1.15'  % (dbse)] =  -0.05340
BIND['%s-NeNe-1.2'   % (dbse)] =  -0.04251
BIND['%s-NeNe-1.3'   % (dbse)] =  -0.02653
BIND['%s-NeNe-1.4'   % (dbse)] =  -0.01673
BIND['%s-NeNe-1.5'   % (dbse)] =  -0.01081
BIND['%s-NeNe-1.6'   % (dbse)] =  -0.00718
BIND['%s-NeNe-1.7'   % (dbse)] =  -0.00489
BIND['%s-NeNe-1.8'   % (dbse)] =  -0.00341
BIND['%s-NeNe-2.0'   % (dbse)] =  -0.00176
BIND['%s-NeNe-2.2'   % (dbse)] =  -0.00098

BIND['%s-NeAr-0.85'  % (dbse)] =   0.26334
BIND['%s-NeAr-0.9'   % (dbse)] =  -0.01713
BIND['%s-NeAr-0.95'  % (dbse)] =  -0.11367
BIND['%s-NeAr-0.975' % (dbse)] =  -0.12900
BIND['%s-NeAr-1.0'   % (dbse)] =  -0.13273  # NeAr minimum
BIND['%s-NeAr-1.025' % (dbse)] =  -0.12947
BIND['%s-NeAr-1.05'  % (dbse)] =  -0.12224
BIND['%s-NeAr-1.1'   % (dbse)] =  -0.10295
BIND['%s-NeAr-1.15'  % (dbse)] =  -0.08326
BIND['%s-NeAr-1.2'   % (dbse)] =  -0.06610
BIND['%s-NeAr-1.3'   % (dbse)] =  -0.04105
BIND['%s-NeAr-1.4'   % (dbse)] =  -0.02577
BIND['%s-NeAr-1.5'   % (dbse)] =  -0.01659
BIND['%s-NeAr-1.6'   % (dbse)] =  -0.01098
BIND['%s-NeAr-1.7'   % (dbse)] =  -0.00746
BIND['%s-NeAr-1.8'   % (dbse)] =  -0.00519
BIND['%s-NeAr-2.0'   % (dbse)] =  -0.00267
BIND['%s-NeAr-2.2'   % (dbse)] =  -0.00148

BIND['%s-NeKr-0.85'  % (dbse)] =   0.26707
BIND['%s-NeKr-0.9'   % (dbse)] =  -0.02063
BIND['%s-NeKr-0.95'  % (dbse)] =  -0.12057
BIND['%s-NeKr-0.975' % (dbse)] =  -0.13659
BIND['%s-NeKr-1.0'   % (dbse)] =  -0.14056  # NeKr minimum
BIND['%s-NeKr-1.025' % (dbse)] =  -0.13722
BIND['%s-NeKr-1.05'  % (dbse)] =  -0.12969
BIND['%s-NeKr-1.1'   % (dbse)] =  -0.10946
BIND['%s-NeKr-1.15'  % (dbse)] =  -0.08868
BIND['%s-NeKr-1.2'   % (dbse)] =  -0.07049
BIND['%s-NeKr-1.3'   % (dbse)] =  -0.04382
BIND['%s-NeKr-1.4'   % (dbse)] =  -0.02751
BIND['%s-NeKr-1.5'   % (dbse)] =  -0.01769
BIND['%s-NeKr-1.6'   % (dbse)] =  -0.01169
BIND['%s-NeKr-1.7'   % (dbse)] =  -0.00793
BIND['%s-NeKr-1.8'   % (dbse)] =  -0.00551
BIND['%s-NeKr-2.0'   % (dbse)] =  -0.00284
BIND['%s-NeKr-2.2'   % (dbse)] =  -0.00156

BIND['%s-ArAr-0.85'  % (dbse)] =   0.63637
BIND['%s-ArAr-0.9'   % (dbse)] =  -0.01138
BIND['%s-ArAr-0.95'  % (dbse)] =  -0.23729
BIND['%s-ArAr-0.975' % (dbse)] =  -0.27458
BIND['%s-ArAr-1.0'   % (dbse)] =  -0.28517  # ArAr minimum
BIND['%s-ArAr-1.025' % (dbse)] =  -0.27957
BIND['%s-ArAr-1.05'  % (dbse)] =  -0.26471
BIND['%s-ArAr-1.1'   % (dbse)] =  -0.22350
BIND['%s-ArAr-1.15'  % (dbse)] =  -0.18080
BIND['%s-ArAr-1.2'   % (dbse)] =  -0.14343
BIND['%s-ArAr-1.3'   % (dbse)] =  -0.08883
BIND['%s-ArAr-1.4'   % (dbse)] =  -0.05560
BIND['%s-ArAr-1.5'   % (dbse)] =  -0.03568
BIND['%s-ArAr-1.6'   % (dbse)] =  -0.02355
BIND['%s-ArAr-1.7'   % (dbse)] =  -0.01596
BIND['%s-ArAr-1.8'   % (dbse)] =  -0.01109
BIND['%s-ArAr-2.0'   % (dbse)] =  -0.00570
BIND['%s-ArAr-2.2'   % (dbse)] =  -0.00314

BIND['%s-ArKr-0.85'  % (dbse)] =   0.69499
BIND['%s-ArKr-0.9'   % (dbse)] =  -0.02873
BIND['%s-ArKr-0.95'  % (dbse)] =  -0.28040
BIND['%s-ArKr-0.975' % (dbse)] =  -0.32132
BIND['%s-ArKr-1.0'   % (dbse)] =  -0.33222  # ArKr minimum
BIND['%s-ArKr-1.025' % (dbse)] =  -0.32494
BIND['%s-ArKr-1.05'  % (dbse)] =  -0.30726
BIND['%s-ArKr-1.1'   % (dbse)] =  -0.25907
BIND['%s-ArKr-1.15'  % (dbse)] =  -0.20945
BIND['%s-ArKr-1.2'   % (dbse)] =  -0.16607
BIND['%s-ArKr-1.3'   % (dbse)] =  -0.10275
BIND['%s-ArKr-1.4'   % (dbse)] =  -0.06422
BIND['%s-ArKr-1.5'   % (dbse)] =  -0.04115
BIND['%s-ArKr-1.6'   % (dbse)] =  -0.02711
BIND['%s-ArKr-1.7'   % (dbse)] =  -0.01836
BIND['%s-ArKr-1.8'   % (dbse)] =  -0.01274
BIND['%s-ArKr-2.0'   % (dbse)] =  -0.00653
BIND['%s-ArKr-2.2'   % (dbse)] =  -0.00359

BIND['%s-KrKr-0.85'  % (dbse)] =   0.81252
BIND['%s-KrKr-0.9'   % (dbse)] =  -0.04237
BIND['%s-KrKr-0.95'  % (dbse)] =  -0.34024
BIND['%s-KrKr-0.975' % (dbse)] =  -0.38858
BIND['%s-KrKr-1.0'   % (dbse)] =  -0.40124  # KrKr minimum
BIND['%s-KrKr-1.025' % (dbse)] =  -0.39225
BIND['%s-KrKr-1.05'  % (dbse)] =  -0.37085
BIND['%s-KrKr-1.1'   % (dbse)] =  -0.31271
BIND['%s-KrKr-1.15'  % (dbse)] =  -0.25283
BIND['%s-KrKr-1.2'   % (dbse)] =  -0.20046
BIND['%s-KrKr-1.3'   % (dbse)] =  -0.12394
BIND['%s-KrKr-1.4'   % (dbse)] =  -0.07736
BIND['%s-KrKr-1.5'   % (dbse)] =  -0.04949
BIND['%s-KrKr-1.6'   % (dbse)] =  -0.03256
BIND['%s-KrKr-1.7'   % (dbse)] =  -0.02201
BIND['%s-KrKr-1.8'   % (dbse)] =  -0.01525
BIND['%s-KrKr-2.0'   % (dbse)] =  -0.00781
BIND['%s-KrKr-2.2'   % (dbse)] =  -0.00429

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')
for item in HeHe:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Helium Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Helium Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Helium from Helium Dimer at %s Req' % (molname.group(2))

for item in HeNe:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Helium-Neon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Helium-Neon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Helium from Helium-Neon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Neon from Helium-Neon Complex at %s Req' % (molname.group(2))

for item in HeAr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Helium-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Helium-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Helium from Helium-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Argon from Helium-Argon Complex at %s Req' % (molname.group(2))

for item in HeKr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Helium-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Helium-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Helium from Helium-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Krypton from Helium-Krypton Complex at %s Req' % (molname.group(2))

for item in NeNe:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Neon Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Neon Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Neon from Neon Dimer at %s Req' % (molname.group(2))

for item in NeAr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Neon-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Neon-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Neon from Neon-Argon Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Argon from Neon-Argon Complex at %s Req' % (molname.group(2))

for item in NeKr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Neon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Neon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Neon from Neon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Krypton from Neon-Krypton Complex at %s Req' % (molname.group(2))

for item in ArAr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Argon Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Argon Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Argon from Argon Dimer at %s Req' % (molname.group(2))

for item in ArKr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Argon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Argon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Argon from Argon-Krypton Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Krypton from Argon-Krypton Complex at %s Req' % (molname.group(2))

for item in KrKr:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'Krypton Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] = 'Krypton Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Krypton from Krypton Dimer at %s Req' % (molname.group(2))

TAGL['%s-He-mono-unCP' % (dbse)] = 'Helium Atom'
TAGL['%s-Ne-mono-unCP' % (dbse)] = 'Neon Atom'
TAGL['%s-Ar-mono-unCP' % (dbse)] = 'Argon Atom'
TAGL['%s-Kr-mono-unCP' % (dbse)] = 'Krypton Atom'

#<<< Geometry Specification Strings >>>
GEOS = {}
rxnpattern2 = re.compile(r'^(..)(..)-(.+)$')

for rxn in HRXN:
    molname = rxnpattern2.match(rxn)
    m1 = molname.group(1)
    m2 = molname.group(2)
    dm = m1 + m2
    Rscal = molname.group(3)
    Rval = float(Rscal) * Req[dm]

    GEOS['%s-%s-%s-%s' % (dbse, dm, Rscal, 'dimer')] = qcdb.Molecule("""
0 1
%(m1)s 0.0 0.0 0.0
--
0 1
%(m2)s 0.0 0.0 R

R = %(Rval)s
units angstrom
""" % vars())

for item in ['He', 'Ne', 'Ar', 'Kr']:
    GEOS['%s-%s-%s' % (dbse, item, 'mono-unCP')] = qcdb.Molecule("""
0 1
%(item)s 0.0 0.0 0.0
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
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.85-dimer'           ] =       0.83565292
DATA['NUCLEAR REPULSION ENERGY']['RGC1-He-mono-unCP'              ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.9-dimer'            ] =       0.78922775
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.95-dimer'           ] =       0.74768945
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.975-dimer'          ] =       0.72851793
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.0-dimer'            ] =       0.71030498
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.025-dimer'          ] =       0.69298047
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.05-dimer'           ] =       0.67648093
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.1-dimer'            ] =       0.64573180
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.15-dimer'           ] =       0.61765650
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.2-dimer'            ] =       0.59192081
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.3-dimer'            ] =       0.54638844
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.4-dimer'            ] =       0.50736070
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.5-dimer'            ] =       0.47353665
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.6-dimer'            ] =       0.44394061
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.7-dimer'            ] =       0.41782646
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.8-dimer'            ] =       0.39461388
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.0-dimer'            ] =       0.35515249
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.2-dimer'            ] =       0.32286590
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.85-dimer'           ] =       4.08236998
DATA['NUCLEAR REPULSION ENERGY']['RGC1-Ne-mono-unCP'              ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.9-dimer'            ] =       3.85557165
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.95-dimer'           ] =       3.65264682
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.975-dimer'          ] =       3.55898921
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.0-dimer'            ] =       3.47001448
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.025-dimer'          ] =       3.38537998
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.05-dimer'           ] =       3.30477570
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.1-dimer'            ] =       3.15455862
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.15-dimer'           ] =       3.01740390
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.2-dimer'            ] =       2.89167874
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.3-dimer'            ] =       2.66924191
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.4-dimer'            ] =       2.47858177
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.5-dimer'            ] =       2.31334299
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.6-dimer'            ] =       2.16875905
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.7-dimer'            ] =       2.04118499
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.8-dimer'            ] =       1.92778582
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.0-dimer'            ] =       1.73500724
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.2-dimer'            ] =       1.57727931
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.85-dimer'           ] =       6.40348891
DATA['NUCLEAR REPULSION ENERGY']['RGC1-Ar-mono-unCP'              ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.9-dimer'            ] =       6.04773953
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.95-dimer'           ] =       5.72943745
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.975-dimer'          ] =       5.58252879
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.0-dimer'            ] =       5.44296557
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.025-dimer'          ] =       5.31021032
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.05-dimer'           ] =       5.18377674
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.1-dimer'            ] =       4.94815052
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.15-dimer'           ] =       4.73301354
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.2-dimer'            ] =       4.53580465
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.3-dimer'            ] =       4.18689660
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.4-dimer'            ] =       3.88783255
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.5-dimer'            ] =       3.62864372
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.6-dimer'            ] =       3.40185348
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.7-dimer'            ] =       3.20174446
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.8-dimer'            ] =       3.02386976
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.0-dimer'            ] =       2.72148279
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.2-dimer'            ] =       2.47407526
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.85-dimer'           ] =      12.11470875
DATA['NUCLEAR REPULSION ENERGY']['RGC1-Kr-mono-unCP'              ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.9-dimer'            ] =      11.44166937
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.95-dimer'           ] =      10.83947625
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.975-dimer'          ] =      10.56154096
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.0-dimer'            ] =      10.29750244
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.025-dimer'          ] =      10.04634384
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.05-dimer'           ] =       9.80714518
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.1-dimer'            ] =       9.36136585
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.15-dimer'           ] =       8.95434995
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.2-dimer'            ] =       8.58125203
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.3-dimer'            ] =       7.92115572
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.4-dimer'            ] =       7.35535888
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.5-dimer'            ] =       6.86500162
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.6-dimer'            ] =       6.43593902
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.7-dimer'            ] =       6.05735437
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.8-dimer'            ] =       5.72083469
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.0-dimer'            ] =       5.14875122
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.2-dimer'            ] =       4.68068293
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.85-dimer'           ] =      20.14761883
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.9-dimer'            ] =      19.02830667
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.95-dimer'           ] =      18.02681685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.975-dimer'          ] =      17.56459078
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.0-dimer'            ] =      17.12547601
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.025-dimer'          ] =      16.70778147
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.05-dimer'           ] =      16.30997715
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.1-dimer'            ] =      15.56861455
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.15-dimer'           ] =      14.89171827
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.2-dimer'            ] =      14.27123001
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.3-dimer'            ] =      13.17344308
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.4-dimer'            ] =      12.23248286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.5-dimer'            ] =      11.41698400
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.6-dimer'            ] =      10.70342250
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.7-dimer'            ] =      10.07380942
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.8-dimer'            ] =       9.51415334
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.0-dimer'            ] =       8.56273800
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.2-dimer'            ] =       7.78430728
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.85-dimer'           ] =      32.20145286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.9-dimer'            ] =      30.41248325
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.95-dimer'           ] =      28.81182624
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.975-dimer'          ] =      28.07306146
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.0-dimer'            ] =      27.37123493
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.025-dimer'          ] =      26.70364383
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.05-dimer'           ] =      26.06784279
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.1-dimer'            ] =      24.88294084
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.15-dimer'           ] =      23.80107385
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.2-dimer'            ] =      22.80936244
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.3-dimer'            ] =      21.05479610
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.4-dimer'            ] =      19.55088209
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.5-dimer'            ] =      18.24748995
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.6-dimer'            ] =      17.10702183
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.7-dimer'            ] =      16.10072643
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.8-dimer'            ] =      15.20624163
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.0-dimer'            ] =      13.68561746
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.2-dimer'            ] =      12.44147042
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.85-dimer'           ] =      61.40331832
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.9-dimer'            ] =      57.99202286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.95-dimer'           ] =      54.93981113
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.975-dimer'          ] =      53.53109802
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.0-dimer'            ] =      52.19282057
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.025-dimer'          ] =      50.91982495
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.05-dimer'           ] =      49.70744817
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.1-dimer'            ] =      47.44801870
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.15-dimer'           ] =      45.38506137
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.2-dimer'            ] =      43.49401714
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.3-dimer'            ] =      40.14832352
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.4-dimer'            ] =      37.28058612
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.5-dimer'            ] =      34.79521372
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.6-dimer'            ] =      32.62051286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.7-dimer'            ] =      30.70165916
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.8-dimer'            ] =      28.99601143
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.0-dimer'            ] =      26.09641029
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.2-dimer'            ] =      23.72400935
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.85-dimer'           ] =      53.78930685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.9-dimer'            ] =      50.80101202
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.95-dimer'           ] =      48.12727455
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.975-dimer'          ] =      46.89324187
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.0-dimer'            ] =      45.72091082
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.025-dimer'          ] =      44.60576666
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.05-dimer'           ] =      43.54372459
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.1-dimer'            ] =      41.56446438
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.15-dimer'           ] =      39.75731376
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.2-dimer'            ] =      38.10075902
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.3-dimer'            ] =      35.16993140
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.4-dimer'            ] =      32.65779344
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.5-dimer'            ] =      30.48060721
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.6-dimer'            ] =      28.57556926
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.7-dimer'            ] =      26.89465342
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.8-dimer'            ] =      25.40050601
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.0-dimer'            ] =      22.86045541
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.2-dimer'            ] =      20.78223219
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.85-dimer'           ] =     103.70688981
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.9-dimer'            ] =      97.94539593
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.95-dimer'           ] =      92.79037510
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.975-dimer'          ] =      90.41113471
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.0-dimer'            ] =      88.15085634
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.025-dimer'          ] =      86.00083545
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.05-dimer'           ] =      83.95319652
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.1-dimer'            ] =      80.13714213
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.15-dimer'           ] =      76.65291856
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.2-dimer'            ] =      73.45904695
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.3-dimer'            ] =      67.80835103
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.4-dimer'            ] =      62.96489739
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.5-dimer'            ] =      58.76723756
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.6-dimer'            ] =      55.09428521
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.7-dimer'            ] =      51.85344491
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.8-dimer'            ] =      48.97269797
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.0-dimer'            ] =      44.07542817
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.2-dimer'            ] =      40.06857106
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.85-dimer'           ] =     201.20688348
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.9-dimer'            ] =     190.02872328
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.95-dimer'           ] =     180.02721153
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.975-dimer'          ] =     175.41112919
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.0-dimer'            ] =     171.02585096
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.025-dimer'          ] =     166.85448874
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.05-dimer'           ] =     162.88176282
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.1-dimer'            ] =     155.47804632
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.15-dimer'           ] =     148.71813127
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.2-dimer'            ] =     142.52154246
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.3-dimer'            ] =     131.55834689
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.4-dimer'            ] =     122.16132211
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.5-dimer'            ] =     114.01723397
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.6-dimer'            ] =     106.89115685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.7-dimer'            ] =     100.60344174
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.8-dimer'            ] =      95.01436164
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.0-dimer'            ] =      85.51292548
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.2-dimer'            ] =      77.73902316
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.85-dimer'           ] =       0.83565292
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.9-dimer'            ] =       0.78922775
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.95-dimer'           ] =       0.74768945
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.975-dimer'          ] =       0.72851793
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.0-dimer'            ] =       0.71030498
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.025-dimer'          ] =       0.69298047
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.05-dimer'           ] =       0.67648093
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.1-dimer'            ] =       0.64573180
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.15-dimer'           ] =       0.61765650
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.2-dimer'            ] =       0.59192081
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.3-dimer'            ] =       0.54638844
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.4-dimer'            ] =       0.50736070
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.5-dimer'            ] =       0.47353665
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.6-dimer'            ] =       0.44394061
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.7-dimer'            ] =       0.41782646
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.8-dimer'            ] =       0.39461388
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.0-dimer'            ] =       0.35515249
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.2-dimer'            ] =       0.32286590
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeHe-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.85-dimer'           ] =       4.08236998
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.9-dimer'            ] =       3.85557165
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.95-dimer'           ] =       3.65264682
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.975-dimer'          ] =       3.55898921
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.0-dimer'            ] =       3.47001448
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.025-dimer'          ] =       3.38537998
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.05-dimer'           ] =       3.30477570
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.1-dimer'            ] =       3.15455862
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.15-dimer'           ] =       3.01740390
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.2-dimer'            ] =       2.89167874
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.3-dimer'            ] =       2.66924191
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.4-dimer'            ] =       2.47858177
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.5-dimer'            ] =       2.31334299
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.6-dimer'            ] =       2.16875905
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.7-dimer'            ] =       2.04118499
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.8-dimer'            ] =       1.92778582
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.0-dimer'            ] =       1.73500724
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.2-dimer'            ] =       1.57727931
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeNe-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.85-dimer'           ] =       6.40348891
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.9-dimer'            ] =       6.04773953
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.95-dimer'           ] =       5.72943745
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.975-dimer'          ] =       5.58252879
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.0-dimer'            ] =       5.44296557
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.025-dimer'          ] =       5.31021032
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.05-dimer'           ] =       5.18377674
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.1-dimer'            ] =       4.94815052
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.15-dimer'           ] =       4.73301354
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.2-dimer'            ] =       4.53580465
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.3-dimer'            ] =       4.18689660
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.4-dimer'            ] =       3.88783255
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.5-dimer'            ] =       3.62864372
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.6-dimer'            ] =       3.40185348
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.7-dimer'            ] =       3.20174446
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.8-dimer'            ] =       3.02386976
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.0-dimer'            ] =       2.72148279
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.2-dimer'            ] =       2.47407526
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeAr-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.85-dimer'           ] =      12.11470875
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.9-dimer'            ] =      11.44166937
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.95-dimer'           ] =      10.83947625
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.975-dimer'          ] =      10.56154096
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.0-dimer'            ] =      10.29750244
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.025-dimer'          ] =      10.04634384
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.05-dimer'           ] =       9.80714518
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.1-dimer'            ] =       9.36136585
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.15-dimer'           ] =       8.95434995
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.2-dimer'            ] =       8.58125203
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.3-dimer'            ] =       7.92115572
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.4-dimer'            ] =       7.35535888
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.5-dimer'            ] =       6.86500162
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.6-dimer'            ] =       6.43593902
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.7-dimer'            ] =       6.05735437
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.8-dimer'            ] =       5.72083469
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.0-dimer'            ] =       5.14875122
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.2-dimer'            ] =       4.68068293
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-HeKr-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.85-dimer'           ] =      20.14761883
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.9-dimer'            ] =      19.02830667
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.95-dimer'           ] =      18.02681685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.975-dimer'          ] =      17.56459078
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.0-dimer'            ] =      17.12547601
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.025-dimer'          ] =      16.70778147
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.05-dimer'           ] =      16.30997715
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.1-dimer'            ] =      15.56861455
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.15-dimer'           ] =      14.89171827
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.2-dimer'            ] =      14.27123001
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.3-dimer'            ] =      13.17344308
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.4-dimer'            ] =      12.23248286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.5-dimer'            ] =      11.41698400
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.6-dimer'            ] =      10.70342250
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.7-dimer'            ] =      10.07380942
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.8-dimer'            ] =       9.51415334
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.0-dimer'            ] =       8.56273800
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.2-dimer'            ] =       7.78430728
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeNe-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.85-dimer'           ] =      32.20145286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.9-dimer'            ] =      30.41248325
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.95-dimer'           ] =      28.81182624
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.975-dimer'          ] =      28.07306146
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.0-dimer'            ] =      27.37123493
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.025-dimer'          ] =      26.70364383
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.05-dimer'           ] =      26.06784279
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.1-dimer'            ] =      24.88294084
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.15-dimer'           ] =      23.80107385
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.2-dimer'            ] =      22.80936244
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.3-dimer'            ] =      21.05479610
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.4-dimer'            ] =      19.55088209
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.5-dimer'            ] =      18.24748995
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.6-dimer'            ] =      17.10702183
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.7-dimer'            ] =      16.10072643
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.8-dimer'            ] =      15.20624163
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.0-dimer'            ] =      13.68561746
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.2-dimer'            ] =      12.44147042
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeAr-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.85-dimer'           ] =      61.40331832
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.9-dimer'            ] =      57.99202286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.95-dimer'           ] =      54.93981113
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.975-dimer'          ] =      53.53109802
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.0-dimer'            ] =      52.19282057
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.025-dimer'          ] =      50.91982495
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.05-dimer'           ] =      49.70744817
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.1-dimer'            ] =      47.44801870
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.15-dimer'           ] =      45.38506137
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.2-dimer'            ] =      43.49401714
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.3-dimer'            ] =      40.14832352
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.4-dimer'            ] =      37.28058612
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.5-dimer'            ] =      34.79521372
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.6-dimer'            ] =      32.62051286
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.7-dimer'            ] =      30.70165916
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.8-dimer'            ] =      28.99601143
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.0-dimer'            ] =      26.09641029
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.2-dimer'            ] =      23.72400935
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-NeKr-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.85-dimer'           ] =      53.78930685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.9-dimer'            ] =      50.80101202
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.95-dimer'           ] =      48.12727455
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.975-dimer'          ] =      46.89324187
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.0-dimer'            ] =      45.72091082
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.025-dimer'          ] =      44.60576666
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.05-dimer'           ] =      43.54372459
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.1-dimer'            ] =      41.56446438
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.15-dimer'           ] =      39.75731376
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.2-dimer'            ] =      38.10075902
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.3-dimer'            ] =      35.16993140
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.4-dimer'            ] =      32.65779344
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.5-dimer'            ] =      30.48060721
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.6-dimer'            ] =      28.57556926
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.7-dimer'            ] =      26.89465342
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.8-dimer'            ] =      25.40050601
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.0-dimer'            ] =      22.86045541
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.2-dimer'            ] =      20.78223219
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArAr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.85-dimer'           ] =     103.70688981
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.85-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.9-dimer'            ] =      97.94539593
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.9-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.95-dimer'           ] =      92.79037510
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.95-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.975-dimer'          ] =      90.41113471
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-0.975-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.0-dimer'            ] =      88.15085634
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.025-dimer'          ] =      86.00083545
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.025-monoB-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.05-dimer'           ] =      83.95319652
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.05-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.1-dimer'            ] =      80.13714213
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.1-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.15-dimer'           ] =      76.65291856
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.15-monoB-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.2-dimer'            ] =      73.45904695
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.3-dimer'            ] =      67.80835103
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.3-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.4-dimer'            ] =      62.96489739
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.4-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.5-dimer'            ] =      58.76723756
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.5-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.6-dimer'            ] =      55.09428521
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.6-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.7-dimer'            ] =      51.85344491
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.7-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.8-dimer'            ] =      48.97269797
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-1.8-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.0-dimer'            ] =      44.07542817
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.0-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.2-dimer'            ] =      40.06857106
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-ArKr-2.2-monoB-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.85-dimer'           ] =     201.20688348
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.85-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.9-dimer'            ] =     190.02872328
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.9-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.95-dimer'           ] =     180.02721153
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.95-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.975-dimer'          ] =     175.41112919
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-0.975-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.0-dimer'            ] =     171.02585096
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.025-dimer'          ] =     166.85448874
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.025-monoA-CP'       ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.05-dimer'           ] =     162.88176282
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.05-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.1-dimer'            ] =     155.47804632
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.1-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.15-dimer'           ] =     148.71813127
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.15-monoA-CP'        ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.2-dimer'            ] =     142.52154246
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.2-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.3-dimer'            ] =     131.55834689
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.3-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.4-dimer'            ] =     122.16132211
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.4-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.5-dimer'            ] =     114.01723397
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.5-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.6-dimer'            ] =     106.89115685
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.6-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.7-dimer'            ] =     100.60344174
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.7-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.8-dimer'            ] =      95.01436164
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-1.8-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.0-dimer'            ] =      85.51292548
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.0-monoA-CP'         ] =       0.00000000
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.2-dimer'            ] =      77.73902316
DATA['NUCLEAR REPULSION ENERGY']['RGC1-KrKr-2.2-monoA-CP'         ] =       0.00000000
