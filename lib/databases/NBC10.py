"""
| Database (Sherrill) of interaction energies for dissociation curves of dispersion-bound bimolecular complexes.
| Geometries and Reference interaction energies from the following articles:
|   Benzene Dimers from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Hydrogen Sulfide from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Methane from Sherrill et al. JPCA 113 10146 (2009).
|   Methane Dimer from Takatani et al. PCCP 9 6106 (2007).
|   Pyridine Dimers from Hohenstein et al. JPCA 113 878 (2009).
|   Collection into NBC10 from Burns et al. JCP 134 084107 (2011).
|   Reference from Marshall et al. JCP 135 194102 (2011).

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
import input

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

# <<< Molecule Specifications >>>
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

NBC1_BzBz_S_3p2 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.200000     0.000000    -0.391500
C       4.200000    -1.205074     0.304250
C       4.200000    -1.205074     1.695750
C       4.200000     0.000000     2.391500
C       4.200000     1.205074     1.695750
C       4.200000     1.205074     0.304250
H       4.200000     0.000000    -1.471500
H       4.200000    -2.140382    -0.235750
H       4.200000    -2.140382     2.235750
H       4.200000     0.000000     3.471500
H       4.200000     2.140382     2.235750
H       4.200000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p3 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.300000     0.000000    -0.391500
C       4.300000    -1.205074     0.304250
C       4.300000    -1.205074     1.695750
C       4.300000     0.000000     2.391500
C       4.300000     1.205074     1.695750
C       4.300000     1.205074     0.304250
H       4.300000     0.000000    -1.471500
H       4.300000    -2.140382    -0.235750
H       4.300000    -2.140382     2.235750
H       4.300000     0.000000     3.471500
H       4.300000     2.140382     2.235750
H       4.300000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p4 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.400000     0.000000    -0.391500
C       4.400000    -1.205074     0.304250
C       4.400000    -1.205074     1.695750
C       4.400000     0.000000     2.391500
C       4.400000     1.205074     1.695750
C       4.400000     1.205074     0.304250
H       4.400000     0.000000    -1.471500
H       4.400000    -2.140382    -0.235750
H       4.400000    -2.140382     2.235750
H       4.400000     0.000000     3.471500
H       4.400000     2.140382     2.235750
H       4.400000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.500000     0.000000    -0.391500
C       4.500000    -1.205074     0.304250
C       4.500000    -1.205074     1.695750
C       4.500000     0.000000     2.391500
C       4.500000     1.205074     1.695750
C       4.500000     1.205074     0.304250
H       4.500000     0.000000    -1.471500
H       4.500000    -2.140382    -0.235750
H       4.500000    -2.140382     2.235750
H       4.500000     0.000000     3.471500
H       4.500000     2.140382     2.235750
H       4.500000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p6 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.600000     0.000000    -0.391500
C       4.600000    -1.205074     0.304250
C       4.600000    -1.205074     1.695750
C       4.600000     0.000000     2.391500
C       4.600000     1.205074     1.695750
C       4.600000     1.205074     0.304250
H       4.600000     0.000000    -1.471500
H       4.600000    -2.140382    -0.235750
H       4.600000    -2.140382     2.235750
H       4.600000     0.000000     3.471500
H       4.600000     2.140382     2.235750
H       4.600000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p7 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.700000     0.000000    -0.391500
C       4.700000    -1.205074     0.304250
C       4.700000    -1.205074     1.695750
C       4.700000     0.000000     2.391500
C       4.700000     1.205074     1.695750
C       4.700000     1.205074     0.304250
H       4.700000     0.000000    -1.471500
H       4.700000    -2.140382    -0.235750
H       4.700000    -2.140382     2.235750
H       4.700000     0.000000     3.471500
H       4.700000     2.140382     2.235750
H       4.700000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p8 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.800000     0.000000    -0.391500
C       4.800000    -1.205074     0.304250
C       4.800000    -1.205074     1.695750
C       4.800000     0.000000     2.391500
C       4.800000     1.205074     1.695750
C       4.800000     1.205074     0.304250
H       4.800000     0.000000    -1.471500
H       4.800000    -2.140382    -0.235750
H       4.800000    -2.140382     2.235750
H       4.800000     0.000000     3.471500
H       4.800000     2.140382     2.235750
H       4.800000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_3p9 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       4.900000     0.000000    -0.391500
C       4.900000    -1.205074     0.304250
C       4.900000    -1.205074     1.695750
C       4.900000     0.000000     2.391500
C       4.900000     1.205074     1.695750
C       4.900000     1.205074     0.304250
H       4.900000     0.000000    -1.471500
H       4.900000    -2.140382    -0.235750
H       4.900000    -2.140382     2.235750
H       4.900000     0.000000     3.471500
H       4.900000     2.140382     2.235750
H       4.900000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_4p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       5.000000     0.000000    -0.391500
C       5.000000    -1.205074     0.304250
C       5.000000    -1.205074     1.695750
C       5.000000     0.000000     2.391500
C       5.000000     1.205074     1.695750
C       5.000000     1.205074     0.304250
H       5.000000     0.000000    -1.471500
H       5.000000    -2.140382    -0.235750
H       5.000000    -2.140382     2.235750
H       5.000000     0.000000     3.471500
H       5.000000     2.140382     2.235750
H       5.000000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_4p1 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       5.100000     0.000000    -0.391500
C       5.100000    -1.205074     0.304250
C       5.100000    -1.205074     1.695750
C       5.100000     0.000000     2.391500
C       5.100000     1.205074     1.695750
C       5.100000     1.205074     0.304250
H       5.100000     0.000000    -1.471500
H       5.100000    -2.140382    -0.235750
H       5.100000    -2.140382     2.235750
H       5.100000     0.000000     3.471500
H       5.100000     2.140382     2.235750
H       5.100000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_4p2 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       5.200000     0.000000    -0.391500
C       5.200000    -1.205074     0.304250
C       5.200000    -1.205074     1.695750
C       5.200000     0.000000     2.391500
C       5.200000     1.205074     1.695750
C       5.200000     1.205074     0.304250
H       5.200000     0.000000    -1.471500
H       5.200000    -2.140382    -0.235750
H       5.200000    -2.140382     2.235750
H       5.200000     0.000000     3.471500
H       5.200000     2.140382     2.235750
H       5.200000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_4p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       5.500000     0.000000    -0.391500
C       5.500000    -1.205074     0.304250
C       5.500000    -1.205074     1.695750
C       5.500000     0.000000     2.391500
C       5.500000     1.205074     1.695750
C       5.500000     1.205074     0.304250
H       5.500000     0.000000    -1.471500
H       5.500000    -2.140382    -0.235750
H       5.500000    -2.140382     2.235750
H       5.500000     0.000000     3.471500
H       5.500000     2.140382     2.235750
H       5.500000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_5p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       6.000000     0.000000    -0.391500
C       6.000000    -1.205074     0.304250
C       6.000000    -1.205074     1.695750
C       6.000000     0.000000     2.391500
C       6.000000     1.205074     1.695750
C       6.000000     1.205074     0.304250
H       6.000000     0.000000    -1.471500
H       6.000000    -2.140382    -0.235750
H       6.000000    -2.140382     2.235750
H       6.000000     0.000000     3.471500
H       6.000000     2.140382     2.235750
H       6.000000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_5p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       6.500000     0.000000    -0.391500
C       6.500000    -1.205074     0.304250
C       6.500000    -1.205074     1.695750
C       6.500000     0.000000     2.391500
C       6.500000     1.205074     1.695750
C       6.500000     1.205074     0.304250
H       6.500000     0.000000    -1.471500
H       6.500000    -2.140382    -0.235750
H       6.500000    -2.140382     2.235750
H       6.500000     0.000000     3.471500
H       6.500000     2.140382     2.235750
H       6.500000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_6p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       7.000000     0.000000    -0.391500
C       7.000000    -1.205074     0.304250
C       7.000000    -1.205074     1.695750
C       7.000000     0.000000     2.391500
C       7.000000     1.205074     1.695750
C       7.000000     1.205074     0.304250
H       7.000000     0.000000    -1.471500
H       7.000000    -2.140382    -0.235750
H       7.000000    -2.140382     2.235750
H       7.000000     0.000000     3.471500
H       7.000000     2.140382     2.235750
H       7.000000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_6p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C       7.500000     0.000000    -0.391500
C       7.500000    -1.205074     0.304250
C       7.500000    -1.205074     1.695750
C       7.500000     0.000000     2.391500
C       7.500000     1.205074     1.695750
C       7.500000     1.205074     0.304250
H       7.500000     0.000000    -1.471500
H       7.500000    -2.140382    -0.235750
H       7.500000    -2.140382     2.235750
H       7.500000     0.000000     3.471500
H       7.500000     2.140382     2.235750
H       7.500000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_S_10p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      11.000000     0.000000    -0.391500
C      11.000000    -1.205074     0.304250
C      11.000000    -1.205074     1.695750
C      11.000000     0.000000     2.391500
C      11.000000     1.205074     1.695750
C      11.000000     1.205074     0.304250
H      11.000000     0.000000    -1.471500
H      11.000000    -2.140382    -0.235750
H      11.000000    -2.140382     2.235750
H      11.000000     0.000000     3.471500
H      11.000000     2.140382     2.235750
H      11.000000     2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_BzBz_T_4p4 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.008500     0.000000     1.000000
C      -2.704250    -1.205074     1.000000
C      -4.095750    -1.205074     1.000000
C      -4.791500     0.000000     1.000000
C      -4.095750     1.205074     1.000000
C      -2.704250     1.205074     1.000000
H      -0.928500     0.000000     1.000000
H      -2.164250    -2.140382     1.000000
H      -4.635750    -2.140382     1.000000
H      -5.871500     0.000000     1.000000
H      -4.635750     2.140382     1.000000
H      -2.164250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_4p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.108500     0.000000     1.000000
C      -2.804250    -1.205074     1.000000
C      -4.195750    -1.205074     1.000000
C      -4.891500     0.000000     1.000000
C      -4.195750     1.205074     1.000000
C      -2.804250     1.205074     1.000000
H      -1.028500     0.000000     1.000000
H      -2.264250    -2.140382     1.000000
H      -4.735750    -2.140382     1.000000
H      -5.971500     0.000000     1.000000
H      -4.735750     2.140382     1.000000
H      -2.264250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_4p6 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.208500     0.000000     1.000000
C      -2.904250    -1.205074     1.000000
C      -4.295750    -1.205074     1.000000
C      -4.991500     0.000000     1.000000
C      -4.295750     1.205074     1.000000
C      -2.904250     1.205074     1.000000
H      -1.128500     0.000000     1.000000
H      -2.364250    -2.140382     1.000000
H      -4.835750    -2.140382     1.000000
H      -6.071500     0.000000     1.000000
H      -4.835750     2.140382     1.000000
H      -2.364250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_4p7 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.308500     0.000000     1.000000
C      -3.004250    -1.205074     1.000000
C      -4.395750    -1.205074     1.000000
C      -5.091500     0.000000     1.000000
C      -4.395750     1.205074     1.000000
C      -3.004250     1.205074     1.000000
H      -1.228500     0.000000     1.000000
H      -2.464250    -2.140382     1.000000
H      -4.935750    -2.140382     1.000000
H      -6.171500     0.000000     1.000000
H      -4.935750     2.140382     1.000000
H      -2.464250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_4p8 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.408500     0.000000     1.000000
C      -3.104250    -1.205074     1.000000
C      -4.495750    -1.205074     1.000000
C      -5.191500     0.000000     1.000000
C      -4.495750     1.205074     1.000000
C      -3.104250     1.205074     1.000000
H      -1.328500     0.000000     1.000000
H      -2.564250    -2.140382     1.000000
H      -5.035750    -2.140382     1.000000
H      -6.271500     0.000000     1.000000
H      -5.035750     2.140382     1.000000
H      -2.564250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_4p9 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.508500     0.000000     1.000000
C      -3.204250    -1.205074     1.000000
C      -4.595750    -1.205074     1.000000
C      -5.291500     0.000000     1.000000
C      -4.595750     1.205074     1.000000
C      -3.204250     1.205074     1.000000
H      -1.428500     0.000000     1.000000
H      -2.664250    -2.140382     1.000000
H      -5.135750    -2.140382     1.000000
H      -6.371500     0.000000     1.000000
H      -5.135750     2.140382     1.000000
H      -2.664250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.608500     0.000000     1.000000
C      -3.304250    -1.205074     1.000000
C      -4.695750    -1.205074     1.000000
C      -5.391500     0.000000     1.000000
C      -4.695750     1.205074     1.000000
C      -3.304250     1.205074     1.000000
H      -1.528500     0.000000     1.000000
H      -2.764250    -2.140382     1.000000
H      -5.235750    -2.140382     1.000000
H      -6.471500     0.000000     1.000000
H      -5.235750     2.140382     1.000000
H      -2.764250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p1 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.708500     0.000000     1.000000
C      -3.404250    -1.205074     1.000000
C      -4.795750    -1.205074     1.000000
C      -5.491500     0.000000     1.000000
C      -4.795750     1.205074     1.000000
C      -3.404250     1.205074     1.000000
H      -1.628500     0.000000     1.000000
H      -2.864250    -2.140382     1.000000
H      -5.335750    -2.140382     1.000000
H      -6.571500     0.000000     1.000000
H      -5.335750     2.140382     1.000000
H      -2.864250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p2 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.808500     0.000000     1.000000
C      -3.504250    -1.205074     1.000000
C      -4.895750    -1.205074     1.000000
C      -5.591500     0.000000     1.000000
C      -4.895750     1.205074     1.000000
C      -3.504250     1.205074     1.000000
H      -1.728500     0.000000     1.000000
H      -2.964250    -2.140382     1.000000
H      -5.435750    -2.140382     1.000000
H      -6.671500     0.000000     1.000000
H      -5.435750     2.140382     1.000000
H      -2.964250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p3 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -2.908500     0.000000     1.000000
C      -3.604250    -1.205074     1.000000
C      -4.995750    -1.205074     1.000000
C      -5.691500     0.000000     1.000000
C      -4.995750     1.205074     1.000000
C      -3.604250     1.205074     1.000000
H      -1.828500     0.000000     1.000000
H      -3.064250    -2.140382     1.000000
H      -5.535750    -2.140382     1.000000
H      -6.771500     0.000000     1.000000
H      -5.535750     2.140382     1.000000
H      -3.064250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p4 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -3.008500     0.000000     1.000000
C      -3.704250    -1.205074     1.000000
C      -5.095750    -1.205074     1.000000
C      -5.791500     0.000000     1.000000
C      -5.095750     1.205074     1.000000
C      -3.704250     1.205074     1.000000
H      -1.928500     0.000000     1.000000
H      -3.164250    -2.140382     1.000000
H      -5.635750    -2.140382     1.000000
H      -6.871500     0.000000     1.000000
H      -5.635750     2.140382     1.000000
H      -3.164250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -3.108500     0.000000     1.000000
C      -3.804250    -1.205074     1.000000
C      -5.195750    -1.205074     1.000000
C      -5.891500     0.000000     1.000000
C      -5.195750     1.205074     1.000000
C      -3.804250     1.205074     1.000000
H      -2.028500     0.000000     1.000000
H      -3.264250    -2.140382     1.000000
H      -5.735750    -2.140382     1.000000
H      -6.971500     0.000000     1.000000
H      -5.735750     2.140382     1.000000
H      -3.264250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_5p6 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -3.208500     0.000000     1.000000
C      -3.904250    -1.205074     1.000000
C      -5.295750    -1.205074     1.000000
C      -5.991500     0.000000     1.000000
C      -5.295750     1.205074     1.000000
C      -3.904250     1.205074     1.000000
H      -2.128500     0.000000     1.000000
H      -3.364250    -2.140382     1.000000
H      -5.835750    -2.140382     1.000000
H      -7.071500     0.000000     1.000000
H      -5.835750     2.140382     1.000000
H      -3.364250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_6p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -3.608500     0.000000     1.000000
C      -4.304250    -1.205074     1.000000
C      -5.695750    -1.205074     1.000000
C      -6.391500     0.000000     1.000000
C      -5.695750     1.205074     1.000000
C      -4.304250     1.205074     1.000000
H      -2.528500     0.000000     1.000000
H      -3.764250    -2.140382     1.000000
H      -6.235750    -2.140382     1.000000
H      -7.471500     0.000000     1.000000
H      -6.235750     2.140382     1.000000
H      -3.764250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_6p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -4.108500     0.000000     1.000000
C      -4.804250    -1.205074     1.000000
C      -6.195750    -1.205074     1.000000
C      -6.891500     0.000000     1.000000
C      -6.195750     1.205074     1.000000
C      -4.804250     1.205074     1.000000
H      -3.028500     0.000000     1.000000
H      -4.264250    -2.140382     1.000000
H      -6.735750    -2.140382     1.000000
H      -7.971500     0.000000     1.000000
H      -6.735750     2.140382     1.000000
H      -4.264250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_7p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -4.608500     0.000000     1.000000
C      -5.304250    -1.205074     1.000000
C      -6.695750    -1.205074     1.000000
C      -7.391500     0.000000     1.000000
C      -6.695750     1.205074     1.000000
C      -5.304250     1.205074     1.000000
H      -3.528500     0.000000     1.000000
H      -4.764250    -2.140382     1.000000
H      -7.235750    -2.140382     1.000000
H      -8.471500     0.000000     1.000000
H      -7.235750     2.140382     1.000000
H      -4.764250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_7p5 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -5.108500     0.000000     1.000000
C      -5.804250    -1.205074     1.000000
C      -7.195750    -1.205074     1.000000
C      -7.891500     0.000000     1.000000
C      -7.195750     1.205074     1.000000
C      -5.804250     1.205074     1.000000
H      -4.028500     0.000000     1.000000
H      -5.264250    -2.140382     1.000000
H      -7.735750    -2.140382     1.000000
H      -8.971500     0.000000     1.000000
H      -7.735750     2.140382     1.000000
H      -5.264250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_T_8p0 = input.process_input("""
molecule dimer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
--
0 1
C      -5.608500     0.000000     1.000000
C      -6.304250    -1.205074     1.000000
C      -7.695750    -1.205074     1.000000
C      -8.391500     0.000000     1.000000
C      -7.695750     1.205074     1.000000
C      -6.304250     1.205074     1.000000
H      -4.528500     0.000000     1.000000
H      -5.764250    -2.140382     1.000000
H      -8.235750    -2.140382     1.000000
H      -9.471500     0.000000     1.000000
H      -8.235750     2.140382     1.000000
H      -5.764250     2.140382     1.000000
units angstrom
}
""", 0)

NBC1_BzBz_PD34_0p2 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.200000     2.391500
C       0.000000     1.005074     1.695750
C       0.000000     1.005074     0.304250
C       0.000000    -0.200000    -0.391500
C       0.000000    -1.405074     0.304250
C       0.000000    -1.405074     1.695750
H       0.000000    -0.200000     3.471500
H       0.000000     1.940382     2.235750
H       0.000000     1.940382    -0.235750
H       0.000000    -0.200000    -1.471500
H       0.000000    -2.340382    -0.235750
H       0.000000    -2.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_0p4 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.400000     2.391500
C       0.000000     0.805074     1.695750
C       0.000000     0.805074     0.304250
C       0.000000    -0.400000    -0.391500
C       0.000000    -1.605074     0.304250
C       0.000000    -1.605074     1.695750
H       0.000000    -0.400000     3.471500
H       0.000000     1.740382     2.235750
H       0.000000     1.740382    -0.235750
H       0.000000    -0.400000    -1.471500
H       0.000000    -2.540382    -0.235750
H       0.000000    -2.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_0p6 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.600000     2.391500
C       0.000000     0.605074     1.695750
C       0.000000     0.605074     0.304250
C       0.000000    -0.600000    -0.391500
C       0.000000    -1.805074     0.304250
C       0.000000    -1.805074     1.695750
H       0.000000    -0.600000     3.471500
H       0.000000     1.540382     2.235750
H       0.000000     1.540382    -0.235750
H       0.000000    -0.600000    -1.471500
H       0.000000    -2.740382    -0.235750
H       0.000000    -2.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_0p8 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.800000     2.391500
C       0.000000     0.405074     1.695750
C       0.000000     0.405074     0.304250
C       0.000000    -0.800000    -0.391500
C       0.000000    -2.005074     0.304250
C       0.000000    -2.005074     1.695750
H       0.000000    -0.800000     3.471500
H       0.000000     1.340382     2.235750
H       0.000000     1.340382    -0.235750
H       0.000000    -0.800000    -1.471500
H       0.000000    -2.940382    -0.235750
H       0.000000    -2.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p0 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.000000     2.391500
C       0.000000     0.205074     1.695750
C       0.000000     0.205074     0.304250
C       0.000000    -1.000000    -0.391500
C       0.000000    -2.205074     0.304250
C       0.000000    -2.205074     1.695750
H       0.000000    -1.000000     3.471500
H       0.000000     1.140382     2.235750
H       0.000000     1.140382    -0.235750
H       0.000000    -1.000000    -1.471500
H       0.000000    -3.140382    -0.235750
H       0.000000    -3.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p2 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.200000     2.391500
C       0.000000     0.005074     1.695750
C       0.000000     0.005074     0.304250
C       0.000000    -1.200000    -0.391500
C       0.000000    -2.405074     0.304250
C       0.000000    -2.405074     1.695750
H       0.000000    -1.200000     3.471500
H       0.000000     0.940382     2.235750
H       0.000000     0.940382    -0.235750
H       0.000000    -1.200000    -1.471500
H       0.000000    -3.340382    -0.235750
H       0.000000    -3.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p4 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.400000     2.391500
C       0.000000    -0.194926     1.695750
C       0.000000    -0.194926     0.304250
C       0.000000    -1.400000    -0.391500
C       0.000000    -2.605074     0.304250
C       0.000000    -2.605074     1.695750
H       0.000000    -1.400000     3.471500
H       0.000000     0.740382     2.235750
H       0.000000     0.740382    -0.235750
H       0.000000    -1.400000    -1.471500
H       0.000000    -3.540382    -0.235750
H       0.000000    -3.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p5 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.500000     2.391500
C       0.000000    -0.294926     1.695750
C       0.000000    -0.294926     0.304250
C       0.000000    -1.500000    -0.391500
C       0.000000    -2.705074     0.304250
C       0.000000    -2.705074     1.695750
H       0.000000    -1.500000     3.471500
H       0.000000     0.640382     2.235750
H       0.000000     0.640382    -0.235750
H       0.000000    -1.500000    -1.471500
H       0.000000    -3.640382    -0.235750
H       0.000000    -3.640382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p6 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.600000     2.391500
C       0.000000    -0.394926     1.695750
C       0.000000    -0.394926     0.304250
C       0.000000    -1.600000    -0.391500
C       0.000000    -2.805074     0.304250
C       0.000000    -2.805074     1.695750
H       0.000000    -1.600000     3.471500
H       0.000000     0.540382     2.235750
H       0.000000     0.540382    -0.235750
H       0.000000    -1.600000    -1.471500
H       0.000000    -3.740382    -0.235750
H       0.000000    -3.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p7 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.700000     2.391500
C       0.000000    -0.494926     1.695750
C       0.000000    -0.494926     0.304250
C       0.000000    -1.700000    -0.391500
C       0.000000    -2.905074     0.304250
C       0.000000    -2.905074     1.695750
H       0.000000    -1.700000     3.471500
H       0.000000     0.440382     2.235750
H       0.000000     0.440382    -0.235750
H       0.000000    -1.700000    -1.471500
H       0.000000    -3.840382    -0.235750
H       0.000000    -3.840382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p8 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.800000     2.391500
C       0.000000    -0.594926     1.695750
C       0.000000    -0.594926     0.304250
C       0.000000    -1.800000    -0.391500
C       0.000000    -3.005074     0.304250
C       0.000000    -3.005074     1.695750
H       0.000000    -1.800000     3.471500
H       0.000000     0.340382     2.235750
H       0.000000     0.340382    -0.235750
H       0.000000    -1.800000    -1.471500
H       0.000000    -3.940382    -0.235750
H       0.000000    -3.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_1p9 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.900000     2.391500
C       0.000000    -0.694926     1.695750
C       0.000000    -0.694926     0.304250
C       0.000000    -1.900000    -0.391500
C       0.000000    -3.105074     0.304250
C       0.000000    -3.105074     1.695750
H       0.000000    -1.900000     3.471500
H       0.000000     0.240382     2.235750
H       0.000000     0.240382    -0.235750
H       0.000000    -1.900000    -1.471500
H       0.000000    -4.040382    -0.235750
H       0.000000    -4.040382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_2p0 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.000000     2.391500
C       0.000000    -0.794926     1.695750
C       0.000000    -0.794926     0.304250
C       0.000000    -2.000000    -0.391500
C       0.000000    -3.205074     0.304250
C       0.000000    -3.205074     1.695750
H       0.000000    -2.000000     3.471500
H       0.000000     0.140382     2.235750
H       0.000000     0.140382    -0.235750
H       0.000000    -2.000000    -1.471500
H       0.000000    -4.140382    -0.235750
H       0.000000    -4.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_2p2 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.200000     2.391500
C       0.000000    -0.994926     1.695750
C       0.000000    -0.994926     0.304250
C       0.000000    -2.200000    -0.391500
C       0.000000    -3.405074     0.304250
C       0.000000    -3.405074     1.695750
H       0.000000    -2.200000     3.471500
H       0.000000    -0.059618     2.235750
H       0.000000    -0.059618    -0.235750
H       0.000000    -2.200000    -1.471500
H       0.000000    -4.340382    -0.235750
H       0.000000    -4.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_2p4 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.400000     2.391500
C       0.000000    -1.194926     1.695750
C       0.000000    -1.194926     0.304250
C       0.000000    -2.400000    -0.391500
C       0.000000    -3.605074     0.304250
C       0.000000    -3.605074     1.695750
H       0.000000    -2.400000     3.471500
H       0.000000    -0.259618     2.235750
H       0.000000    -0.259618    -0.235750
H       0.000000    -2.400000    -1.471500
H       0.000000    -4.540382    -0.235750
H       0.000000    -4.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_2p6 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.600000     2.391500
C       0.000000    -1.394926     1.695750
C       0.000000    -1.394926     0.304250
C       0.000000    -2.600000    -0.391500
C       0.000000    -3.805074     0.304250
C       0.000000    -3.805074     1.695750
H       0.000000    -2.600000     3.471500
H       0.000000    -0.459618     2.235750
H       0.000000    -0.459618    -0.235750
H       0.000000    -2.600000    -1.471500
H       0.000000    -4.740382    -0.235750
H       0.000000    -4.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_2p8 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.800000     2.391500
C       0.000000    -1.594926     1.695750
C       0.000000    -1.594926     0.304250
C       0.000000    -2.800000    -0.391500
C       0.000000    -4.005074     0.304250
C       0.000000    -4.005074     1.695750
H       0.000000    -2.800000     3.471500
H       0.000000    -0.659618     2.235750
H       0.000000    -0.659618    -0.235750
H       0.000000    -2.800000    -1.471500
H       0.000000    -4.940382    -0.235750
H       0.000000    -4.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD34_3p0 = input.process_input("""
molecule dimer {
0 1
C       3.400000     0.000000    -0.391500
C       3.400000     1.205074     0.304250
C       3.400000     1.205074     1.695750
C       3.400000     0.000000     2.391500
C       3.400000    -1.205074     1.695750
C       3.400000    -1.205074     0.304250
H       3.400000     0.000000    -1.471500
H       3.400000     2.140382    -0.235750
H       3.400000     2.140382     2.235750
H       3.400000     0.000000     3.471500
H       3.400000    -2.140382     2.235750
H       3.400000    -2.140382    -0.235750
--
0 1
C       0.000000    -3.000000     2.391500
C       0.000000    -1.794926     1.695750
C       0.000000    -1.794926     0.304250
C       0.000000    -3.000000    -0.391500
C       0.000000    -4.205074     0.304250
C       0.000000    -4.205074     1.695750
H       0.000000    -3.000000     3.471500
H       0.000000    -0.859618     2.235750
H       0.000000    -0.859618    -0.235750
H       0.000000    -3.000000    -1.471500
H       0.000000    -5.140382    -0.235750
H       0.000000    -5.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzH2S_3p2 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.200000
H      -0.961721     0.000000    -1.273221
H       0.961721     0.000000    -1.273221
units angstrom
}
""", 0)

NBC1_BzH2S_3p4 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.400000
H      -0.961721     0.000000    -1.473221
H       0.961721     0.000000    -1.473221
units angstrom
}
""", 0)

NBC1_BzH2S_3p5 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.500000
H      -0.961721     0.000000    -1.573221
H       0.961721     0.000000    -1.573221
units angstrom
}
""", 0)

NBC1_BzH2S_3p6 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.600000
H      -0.961721     0.000000    -1.673221
H       0.961721     0.000000    -1.673221
units angstrom
}
""", 0)

NBC1_BzH2S_3p7 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.700000
H      -0.961721     0.000000    -1.773221
H       0.961721     0.000000    -1.773221
units angstrom
}
""", 0)

NBC1_BzH2S_3p8 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.800000
H      -0.961721     0.000000    -1.873221
H       0.961721     0.000000    -1.873221
units angstrom
}
""", 0)

NBC1_BzH2S_3p9 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -2.900000
H      -0.961721     0.000000    -1.973221
H       0.961721     0.000000    -1.973221
units angstrom
}
""", 0)

NBC1_BzH2S_4p0 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -3.000000
H      -0.961721     0.000000    -2.073221
H       0.961721     0.000000    -2.073221
units angstrom
}
""", 0)

NBC1_BzH2S_4p1 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -3.100000
H      -0.961721     0.000000    -2.173221
H       0.961721     0.000000    -2.173221
units angstrom
}
""", 0)

NBC1_BzH2S_4p2 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -3.200000
H      -0.961721     0.000000    -2.273221
H       0.961721     0.000000    -2.273221
units angstrom
}
""", 0)

NBC1_BzH2S_4p5 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -3.500000
H      -0.961721     0.000000    -2.573221
H       0.961721     0.000000    -2.573221
units angstrom
}
""", 0)

NBC1_BzH2S_4p75 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -3.750000
H      -0.961721     0.000000    -2.823221
H       0.961721     0.000000    -2.823221
units angstrom
}
""", 0)

NBC1_BzH2S_5p0 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -4.000000
H      -0.961721     0.000000    -3.073221
H       0.961721     0.000000    -3.073221
units angstrom
}
""", 0)

NBC1_BzH2S_5p25 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -4.250000
H      -0.961721     0.000000    -3.323221
H       0.961721     0.000000    -3.323221
units angstrom
}
""", 0)

NBC1_BzH2S_5p5 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -4.500000
H      -0.961721     0.000000    -3.573221
H       0.961721     0.000000    -3.573221
units angstrom
}
""", 0)

NBC1_BzH2S_6p0 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -5.000000
H      -0.961721     0.000000    -4.073221
H       0.961721     0.000000    -4.073221
units angstrom
}
""", 0)

NBC1_BzH2S_6p5 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -5.500000
H      -0.961721     0.000000    -4.573221
H       0.961721     0.000000    -4.573221
units angstrom
}
""", 0)

NBC1_BzH2S_7p0 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -6.000000
H      -0.961721     0.000000    -5.073221
H       0.961721     0.000000    -5.073221
units angstrom
}
""", 0)

NBC1_BzH2S_7p5 = input.process_input("""
molecule dimer {
0 1
C       1.391500     0.000000     1.000000
C       0.695750     1.205074     1.000000
C      -0.695750     1.205074     1.000000
C      -1.391500     0.000000     1.000000
C      -0.695750    -1.205074     1.000000
C       0.695750    -1.205074     1.000000
H       2.471500     0.000000     1.000000
H       1.235750     2.140382     1.000000
H      -1.235750     2.140382     1.000000
H      -2.471500     0.000000     1.000000
H      -1.235750    -2.140382     1.000000
H       1.235750    -2.140382     1.000000
--
0 1
S       0.000000     0.000000    -6.500000
H      -0.961721     0.000000    -5.573221
H       0.961721     0.000000    -5.573221
units angstrom
}
""", 0)

NBC1_BzMe_3p2 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.200000
H       0.000000     0.000000    -1.100497
H       1.036621     0.000000    -2.566501
H      -0.518311    -0.897741    -2.566501
H      -0.518311     0.897741    -2.566501
units angstrom
}
""", 0)

NBC1_BzMe_3p3 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.300000
H       0.000000     0.000000    -1.200497
H       1.036621     0.000000    -2.666501
H      -0.518311    -0.897741    -2.666501
H      -0.518311     0.897741    -2.666501
units angstrom
}
""", 0)

NBC1_BzMe_3p4 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.400000
H       0.000000     0.000000    -1.300497
H       1.036621     0.000000    -2.766501
H      -0.518311    -0.897741    -2.766501
H      -0.518311     0.897741    -2.766501
units angstrom
}
""", 0)

NBC1_BzMe_3p5 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.500000
H       0.000000     0.000000    -1.400497
H       1.036621     0.000000    -2.866501
H      -0.518311    -0.897741    -2.866501
H      -0.518311     0.897741    -2.866501
units angstrom
}
""", 0)

NBC1_BzMe_3p6 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.600000
H       0.000000     0.000000    -1.500497
H       1.036621     0.000000    -2.966501
H      -0.518311    -0.897741    -2.966501
H      -0.518311     0.897741    -2.966501
units angstrom
}
""", 0)

NBC1_BzMe_3p7 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.700000
H       0.000000     0.000000    -1.600497
H       1.036621     0.000000    -3.066501
H      -0.518311    -0.897741    -3.066501
H      -0.518311     0.897741    -3.066501
units angstrom
}
""", 0)

NBC1_BzMe_3p8 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.800000
H       0.000000     0.000000    -1.700497
H       1.036621     0.000000    -3.166501
H      -0.518311    -0.897741    -3.166501
H      -0.518311     0.897741    -3.166501
units angstrom
}
""", 0)

NBC1_BzMe_3p9 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -2.900000
H       0.000000     0.000000    -1.800497
H       1.036621     0.000000    -3.266501
H      -0.518311    -0.897741    -3.266501
H      -0.518311     0.897741    -3.266501
units angstrom
}
""", 0)

NBC1_BzMe_4p0 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.000000
H       0.000000     0.000000    -1.900497
H       1.036621     0.000000    -3.366501
H      -0.518311    -0.897741    -3.366501
H      -0.518311     0.897741    -3.366501
units angstrom
}
""", 0)

NBC1_BzMe_4p1 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.100000
H       0.000000     0.000000    -2.000497
H       1.036621     0.000000    -3.466501
H      -0.518311    -0.897741    -3.466501
H      -0.518311     0.897741    -3.466501
units angstrom
}
""", 0)

NBC1_BzMe_4p2 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.200000
H       0.000000     0.000000    -2.100497
H       1.036621     0.000000    -3.566501
H      -0.518311    -0.897741    -3.566501
H      -0.518311     0.897741    -3.566501
units angstrom
}
""", 0)

NBC1_BzMe_4p4 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.400000
H       0.000000     0.000000    -2.300497
H       1.036621     0.000000    -3.766501
H      -0.518311    -0.897741    -3.766501
H      -0.518311     0.897741    -3.766501
units angstrom
}
""", 0)

NBC1_BzMe_4p6 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.600000
H       0.000000     0.000000    -2.500497
H       1.036621     0.000000    -3.966501
H      -0.518311    -0.897741    -3.966501
H      -0.518311     0.897741    -3.966501
units angstrom
}
""", 0)

NBC1_BzMe_4p8 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -3.800000
H       0.000000     0.000000    -2.700497
H       1.036621     0.000000    -4.166501
H      -0.518311    -0.897741    -4.166501
H      -0.518311     0.897741    -4.166501
units angstrom
}
""", 0)

NBC1_BzMe_5p0 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -4.000000
H       0.000000     0.000000    -2.900497
H       1.036621     0.000000    -4.366501
H      -0.518311    -0.897741    -4.366501
H      -0.518311     0.897741    -4.366501
units angstrom
}
""", 0)

NBC1_BzMe_5p2 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -4.200000
H       0.000000     0.000000    -3.100497
H       1.036621     0.000000    -4.566501
H      -0.518311    -0.897741    -4.566501
H      -0.518311     0.897741    -4.566501
units angstrom
}
""", 0)

NBC1_BzMe_5p4 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -4.400000
H       0.000000     0.000000    -3.300497
H       1.036621     0.000000    -4.766501
H      -0.518311    -0.897741    -4.766501
H      -0.518311     0.897741    -4.766501
units angstrom
}
""", 0)

NBC1_BzMe_5p6 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -4.600000
H       0.000000     0.000000    -3.500497
H       1.036621     0.000000    -4.966501
H      -0.518311    -0.897741    -4.966501
H      -0.518311     0.897741    -4.966501
units angstrom
}
""", 0)

NBC1_BzMe_6p0 = input.process_input("""
molecule dimer {
0 1
C       1.405731     0.000000     1.000000
C       0.702865     1.217399     1.000000
C      -0.702865     1.217399     1.000000
C      -1.405731     0.000000     1.000000
C      -0.702865    -1.217399     1.000000
C       0.702866    -1.217399     1.000000
H       2.500941     0.000000     1.000000
H       1.250471     2.165878     1.000000
H      -1.250470     2.165878     1.000000
H      -2.500941     0.000000     1.000000
H      -1.250470    -2.165878     1.000000
H       1.250471    -2.165878     1.000000
--
0 1
C       0.000000     0.000000    -5.000000
H       0.000000     0.000000    -3.900497
H       1.036621     0.000000    -5.366501
H      -0.518311    -0.897741    -5.366501
H      -0.518311     0.897741    -5.366501
units angstrom
}
""", 0)

NBC1_MeMe_3p2 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.200000
H       0.000000     0.000000    -4.299503
H      -1.036621     0.000000    -2.833499
H       0.518311     0.897741    -2.833499
H       0.518311    -0.897741    -2.833499
units angstrom
}
""", 0)

NBC1_MeMe_3p3 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.300000
H       0.000000     0.000000    -4.399503
H      -1.036621     0.000000    -2.933499
H       0.518311     0.897741    -2.933499
H       0.518311    -0.897741    -2.933499
units angstrom
}
""", 0)

NBC1_MeMe_3p4 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.400000
H       0.000000     0.000000    -4.499503
H      -1.036621     0.000000    -3.033499
H       0.518311     0.897741    -3.033499
H       0.518311    -0.897741    -3.033499
units angstrom
}
""", 0)

NBC1_MeMe_3p5 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.500000
H       0.000000     0.000000    -4.599503
H      -1.036621     0.000000    -3.133499
H       0.518311     0.897741    -3.133499
H       0.518311    -0.897741    -3.133499
units angstrom
}
""", 0)

NBC1_MeMe_3p6 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.600000
H       0.000000     0.000000    -4.699503
H      -1.036621     0.000000    -3.233499
H       0.518311     0.897741    -3.233499
H       0.518311    -0.897741    -3.233499
units angstrom
}
""", 0)

NBC1_MeMe_3p7 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.700000
H       0.000000     0.000000    -4.799503
H      -1.036621     0.000000    -3.333499
H       0.518311     0.897741    -3.333499
H       0.518311    -0.897741    -3.333499
units angstrom
}
""", 0)

NBC1_MeMe_3p8 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.800000
H       0.000000     0.000000    -4.899503
H      -1.036621     0.000000    -3.433499
H       0.518311     0.897741    -3.433499
H       0.518311    -0.897741    -3.433499
units angstrom
}
""", 0)

NBC1_MeMe_3p9 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -3.900000
H       0.000000     0.000000    -4.999503
H      -1.036621     0.000000    -3.533499
H       0.518311     0.897741    -3.533499
H       0.518311    -0.897741    -3.533499
units angstrom
}
""", 0)

NBC1_MeMe_4p0 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.000000
H       0.000000     0.000000    -5.099503
H      -1.036621     0.000000    -3.633499
H       0.518311     0.897741    -3.633499
H       0.518311    -0.897741    -3.633499
units angstrom
}
""", 0)

NBC1_MeMe_4p1 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.100000
H       0.000000     0.000000    -5.199503
H      -1.036621     0.000000    -3.733499
H       0.518311     0.897741    -3.733499
H       0.518311    -0.897741    -3.733499
units angstrom
}
""", 0)

NBC1_MeMe_4p2 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.200000
H       0.000000     0.000000    -5.299503
H      -1.036621     0.000000    -3.833499
H       0.518311     0.897741    -3.833499
H       0.518311    -0.897741    -3.833499
units angstrom
}
""", 0)

NBC1_MeMe_4p3 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.300000
H       0.000000     0.000000    -5.399503
H      -1.036621     0.000000    -3.933499
H       0.518311     0.897741    -3.933499
H       0.518311    -0.897741    -3.933499
units angstrom
}
""", 0)

NBC1_MeMe_4p4 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.400000
H       0.000000     0.000000    -5.499503
H      -1.036621     0.000000    -4.033499
H       0.518311     0.897741    -4.033499
H       0.518311    -0.897741    -4.033499
units angstrom
}
""", 0)

NBC1_MeMe_4p6 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.600000
H       0.000000     0.000000    -5.699503
H      -1.036621     0.000000    -4.233499
H       0.518311     0.897741    -4.233499
H       0.518311    -0.897741    -4.233499
units angstrom
}
""", 0)

NBC1_MeMe_4p8 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -4.800000
H       0.000000     0.000000    -5.899503
H      -1.036621     0.000000    -4.433499
H       0.518311     0.897741    -4.433499
H       0.518311    -0.897741    -4.433499
units angstrom
}
""", 0)

NBC1_MeMe_5p0 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -5.000000
H       0.000000     0.000000    -6.099503
H      -1.036621     0.000000    -4.633499
H       0.518311     0.897741    -4.633499
H       0.518311    -0.897741    -4.633499
units angstrom
}
""", 0)

NBC1_MeMe_5p4 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -5.400000
H       0.000000     0.000000    -6.499503
H      -1.036621     0.000000    -5.033499
H       0.518311     0.897741    -5.033499
H       0.518311    -0.897741    -5.033499
units angstrom
}
""", 0)

NBC1_MeMe_5p8 = input.process_input("""
molecule dimer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
--
0 1
C       0.000000     0.000000    -5.800000
H       0.000000     0.000000    -6.899503
H      -1.036621     0.000000    -5.433499
H       0.518311     0.897741    -5.433499
H       0.518311    -0.897741    -5.433499
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p1 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.100000
C      -0.698537     1.140130     3.100000
C       0.694182     1.195340     3.100000
C       1.406747     0.000000     3.100000
C       0.694182    -1.195340     3.100000
C      -0.698537    -1.140130     3.100000
H      -1.281669     2.056885     3.100000
H       1.199925     2.152548     3.100000
H       2.488547     0.000000     3.100000
H       1.199925    -2.152548     3.100000
H      -1.281669    -2.056885     3.100000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p3 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.300000
C      -0.698537     1.140130     3.300000
C       0.694182     1.195340     3.300000
C       1.406747     0.000000     3.300000
C       0.694182    -1.195340     3.300000
C      -0.698537    -1.140130     3.300000
H      -1.281669     2.056885     3.300000
H       1.199925     2.152548     3.300000
H       2.488547     0.000000     3.300000
H       1.199925    -2.152548     3.300000
H      -1.281669    -2.056885     3.300000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p4 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.400000
C      -0.698537     1.140130     3.400000
C       0.694182     1.195340     3.400000
C       1.406747     0.000000     3.400000
C       0.694182    -1.195340     3.400000
C      -0.698537    -1.140130     3.400000
H      -1.281669     2.056885     3.400000
H       1.199925     2.152548     3.400000
H       2.488547     0.000000     3.400000
H       1.199925    -2.152548     3.400000
H      -1.281669    -2.056885     3.400000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.500000
C      -0.698537     1.140130     3.500000
C       0.694182     1.195340     3.500000
C       1.406747     0.000000     3.500000
C       0.694182    -1.195340     3.500000
C      -0.698537    -1.140130     3.500000
H      -1.281669     2.056885     3.500000
H       1.199925     2.152548     3.500000
H       2.488547     0.000000     3.500000
H       1.199925    -2.152548     3.500000
H      -1.281669    -2.056885     3.500000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p6 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.600000
C      -0.698537     1.140130     3.600000
C       0.694182     1.195340     3.600000
C       1.406747     0.000000     3.600000
C       0.694182    -1.195340     3.600000
C      -0.698537    -1.140130     3.600000
H      -1.281669     2.056885     3.600000
H       1.199925     2.152548     3.600000
H       2.488547     0.000000     3.600000
H       1.199925    -2.152548     3.600000
H      -1.281669    -2.056885     3.600000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p7 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.700000
C      -0.698537     1.140130     3.700000
C       0.694182     1.195340     3.700000
C       1.406747     0.000000     3.700000
C       0.694182    -1.195340     3.700000
C      -0.698537    -1.140130     3.700000
H      -1.281669     2.056885     3.700000
H       1.199925     2.152548     3.700000
H       2.488547     0.000000     3.700000
H       1.199925    -2.152548     3.700000
H      -1.281669    -2.056885     3.700000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p8 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.800000
C      -0.698537     1.140130     3.800000
C       0.694182     1.195340     3.800000
C       1.406747     0.000000     3.800000
C       0.694182    -1.195340     3.800000
C      -0.698537    -1.140130     3.800000
H      -1.281669     2.056885     3.800000
H       1.199925     2.152548     3.800000
H       2.488547     0.000000     3.800000
H       1.199925    -2.152548     3.800000
H      -1.281669    -2.056885     3.800000
units angstrom
}
""", 0)

NBC1_PyPy_S2_3p9 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     3.900000
C      -0.698537     1.140130     3.900000
C       0.694182     1.195340     3.900000
C       1.406747     0.000000     3.900000
C       0.694182    -1.195340     3.900000
C      -0.698537    -1.140130     3.900000
H      -1.281669     2.056885     3.900000
H       1.199925     2.152548     3.900000
H       2.488547     0.000000     3.900000
H       1.199925    -2.152548     3.900000
H      -1.281669    -2.056885     3.900000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.000000
C      -0.698537     1.140130     4.000000
C       0.694182     1.195340     4.000000
C       1.406747     0.000000     4.000000
C       0.694182    -1.195340     4.000000
C      -0.698537    -1.140130     4.000000
H      -1.281669     2.056885     4.000000
H       1.199925     2.152548     4.000000
H       2.488547     0.000000     4.000000
H       1.199925    -2.152548     4.000000
H      -1.281669    -2.056885     4.000000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p1 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.100000
C      -0.698537     1.140130     4.100000
C       0.694182     1.195340     4.100000
C       1.406747     0.000000     4.100000
C       0.694182    -1.195340     4.100000
C      -0.698537    -1.140130     4.100000
H      -1.281669     2.056885     4.100000
H       1.199925     2.152548     4.100000
H       2.488547     0.000000     4.100000
H       1.199925    -2.152548     4.100000
H      -1.281669    -2.056885     4.100000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p2 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.200000
C      -0.698537     1.140130     4.200000
C       0.694182     1.195340     4.200000
C       1.406747     0.000000     4.200000
C       0.694182    -1.195340     4.200000
C      -0.698537    -1.140130     4.200000
H      -1.281669     2.056885     4.200000
H       1.199925     2.152548     4.200000
H       2.488547     0.000000     4.200000
H       1.199925    -2.152548     4.200000
H      -1.281669    -2.056885     4.200000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p3 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.300000
C      -0.698537     1.140130     4.300000
C       0.694182     1.195340     4.300000
C       1.406747     0.000000     4.300000
C       0.694182    -1.195340     4.300000
C      -0.698537    -1.140130     4.300000
H      -1.281669     2.056885     4.300000
H       1.199925     2.152548     4.300000
H       2.488547     0.000000     4.300000
H       1.199925    -2.152548     4.300000
H      -1.281669    -2.056885     4.300000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p4 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.400000
C      -0.698537     1.140130     4.400000
C       0.694182     1.195340     4.400000
C       1.406747     0.000000     4.400000
C       0.694182    -1.195340     4.400000
C      -0.698537    -1.140130     4.400000
H      -1.281669     2.056885     4.400000
H       1.199925     2.152548     4.400000
H       2.488547     0.000000     4.400000
H       1.199925    -2.152548     4.400000
H      -1.281669    -2.056885     4.400000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.500000
C      -0.698537     1.140130     4.500000
C       0.694182     1.195340     4.500000
C       1.406747     0.000000     4.500000
C       0.694182    -1.195340     4.500000
C      -0.698537    -1.140130     4.500000
H      -1.281669     2.056885     4.500000
H       1.199925     2.152548     4.500000
H       2.488547     0.000000     4.500000
H       1.199925    -2.152548     4.500000
H      -1.281669    -2.056885     4.500000
units angstrom
}
""", 0)

NBC1_PyPy_S2_4p7 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.700000
C      -0.698537     1.140130     4.700000
C       0.694182     1.195340     4.700000
C       1.406747     0.000000     4.700000
C       0.694182    -1.195340     4.700000
C      -0.698537    -1.140130     4.700000
H      -1.281669     2.056885     4.700000
H       1.199925     2.152548     4.700000
H       2.488547     0.000000     4.700000
H       1.199925    -2.152548     4.700000
H      -1.281669    -2.056885     4.700000
units angstrom
}
""", 0)

NBC1_PyPy_S2_5p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.000000
C      -0.698537     1.140130     5.000000
C       0.694182     1.195340     5.000000
C       1.406747     0.000000     5.000000
C       0.694182    -1.195340     5.000000
C      -0.698537    -1.140130     5.000000
H      -1.281669     2.056885     5.000000
H       1.199925     2.152548     5.000000
H       2.488547     0.000000     5.000000
H       1.199925    -2.152548     5.000000
H      -1.281669    -2.056885     5.000000
units angstrom
}
""", 0)

NBC1_PyPy_S2_5p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.500000
C      -0.698537     1.140130     5.500000
C       0.694182     1.195340     5.500000
C       1.406747     0.000000     5.500000
C       0.694182    -1.195340     5.500000
C      -0.698537    -1.140130     5.500000
H      -1.281669     2.056885     5.500000
H       1.199925     2.152548     5.500000
H       2.488547     0.000000     5.500000
H       1.199925    -2.152548     5.500000
H      -1.281669    -2.056885     5.500000
units angstrom
}
""", 0)

NBC1_PyPy_S2_6p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     6.000000
C      -0.698537     1.140130     6.000000
C       0.694182     1.195340     6.000000
C       1.406747     0.000000     6.000000
C       0.694182    -1.195340     6.000000
C      -0.698537    -1.140130     6.000000
H      -1.281669     2.056885     6.000000
H       1.199925     2.152548     6.000000
H       2.488547     0.000000     6.000000
H       1.199925    -2.152548     6.000000
H      -1.281669    -2.056885     6.000000
units angstrom
}
""", 0)

NBC1_PyPy_S2_6p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     6.500000
C      -0.698537     1.140130     6.500000
C       0.694182     1.195340     6.500000
C       1.406747     0.000000     6.500000
C       0.694182    -1.195340     6.500000
C      -0.698537    -1.140130     6.500000
H      -1.281669     2.056885     6.500000
H       1.199925     2.152548     6.500000
H       2.488547     0.000000     6.500000
H       1.199925    -2.152548     6.500000
H      -1.281669    -2.056885     6.500000
units angstrom
}
""", 0)

NBC1_PyPy_S2_7p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     7.000000
C      -0.698537     1.140130     7.000000
C       0.694182     1.195340     7.000000
C       1.406747     0.000000     7.000000
C       0.694182    -1.195340     7.000000
C      -0.698537    -1.140130     7.000000
H      -1.281669     2.056885     7.000000
H       1.199925     2.152548     7.000000
H       2.488547     0.000000     7.000000
H       1.199925    -2.152548     7.000000
H      -1.281669    -2.056885     7.000000
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p1 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.100000
C      -0.698537     0.000000     5.240130
C       0.694182     0.000000     5.295340
C       1.406747     0.000000     4.100000
C       0.694182     0.000000     2.904660
C      -0.698537     0.000000     2.959870
H      -1.281669     0.000000     6.156885
H       1.199925     0.000000     6.252548
H       2.488547     0.000000     4.100000
H       1.199925     0.000000     1.947452
H      -1.281669     0.000000     2.043115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p3 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.300000
C      -0.698537     0.000000     5.440130
C       0.694182     0.000000     5.495340
C       1.406747     0.000000     4.300000
C       0.694182     0.000000     3.104660
C      -0.698537     0.000000     3.159870
H      -1.281669     0.000000     6.356885
H       1.199925     0.000000     6.452548
H       2.488547     0.000000     4.300000
H       1.199925     0.000000     2.147452
H      -1.281669     0.000000     2.243115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.500000
C      -0.698537     0.000000     5.640130
C       0.694182     0.000000     5.695340
C       1.406747     0.000000     4.500000
C       0.694182     0.000000     3.304660
C      -0.698537     0.000000     3.359870
H      -1.281669     0.000000     6.556885
H       1.199925     0.000000     6.652548
H       2.488547     0.000000     4.500000
H       1.199925     0.000000     2.347452
H      -1.281669     0.000000     2.443115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p6 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.600000
C      -0.698537     0.000000     5.740130
C       0.694182     0.000000     5.795340
C       1.406747     0.000000     4.600000
C       0.694182     0.000000     3.404660
C      -0.698537     0.000000     3.459870
H      -1.281669     0.000000     6.656885
H       1.199925     0.000000     6.752548
H       2.488547     0.000000     4.600000
H       1.199925     0.000000     2.447452
H      -1.281669     0.000000     2.543115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p7 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.700000
C      -0.698537     0.000000     5.840130
C       0.694182     0.000000     5.895340
C       1.406747     0.000000     4.700000
C       0.694182     0.000000     3.504660
C      -0.698537     0.000000     3.559870
H      -1.281669     0.000000     6.756885
H       1.199925     0.000000     6.852548
H       2.488547     0.000000     4.700000
H       1.199925     0.000000     2.547452
H      -1.281669     0.000000     2.643115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p8 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.800000
C      -0.698537     0.000000     5.940130
C       0.694182     0.000000     5.995340
C       1.406747     0.000000     4.800000
C       0.694182     0.000000     3.604660
C      -0.698537     0.000000     3.659870
H      -1.281669     0.000000     6.856885
H       1.199925     0.000000     6.952548
H       2.488547     0.000000     4.800000
H       1.199925     0.000000     2.647452
H      -1.281669     0.000000     2.743115
units angstrom
}
""", 0)

NBC1_PyPy_T3_4p9 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     4.900000
C      -0.698537     0.000000     6.040130
C       0.694182     0.000000     6.095340
C       1.406747     0.000000     4.900000
C       0.694182     0.000000     3.704660
C      -0.698537     0.000000     3.759870
H      -1.281669     0.000000     6.956885
H       1.199925     0.000000     7.052548
H       2.488547     0.000000     4.900000
H       1.199925     0.000000     2.747452
H      -1.281669     0.000000     2.843115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.000000
C      -0.698537     0.000000     6.140130
C       0.694182     0.000000     6.195340
C       1.406747     0.000000     5.000000
C       0.694182     0.000000     3.804660
C      -0.698537     0.000000     3.859870
H      -1.281669     0.000000     7.056885
H       1.199925     0.000000     7.152548
H       2.488547     0.000000     5.000000
H       1.199925     0.000000     2.847452
H      -1.281669     0.000000     2.943115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p1 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.100000
C      -0.698537     0.000000     6.240130
C       0.694182     0.000000     6.295340
C       1.406747     0.000000     5.100000
C       0.694182     0.000000     3.904660
C      -0.698537     0.000000     3.959870
H      -1.281669     0.000000     7.156885
H       1.199925     0.000000     7.252548
H       2.488547     0.000000     5.100000
H       1.199925     0.000000     2.947452
H      -1.281669     0.000000     3.043115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p2 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.200000
C      -0.698537     0.000000     6.340130
C       0.694182     0.000000     6.395340
C       1.406747     0.000000     5.200000
C       0.694182     0.000000     4.004660
C      -0.698537     0.000000     4.059870
H      -1.281669     0.000000     7.256885
H       1.199925     0.000000     7.352548
H       2.488547     0.000000     5.200000
H       1.199925     0.000000     3.047452
H      -1.281669     0.000000     3.143115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p3 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.300000
C      -0.698537     0.000000     6.440130
C       0.694182     0.000000     6.495340
C       1.406747     0.000000     5.300000
C       0.694182     0.000000     4.104660
C      -0.698537     0.000000     4.159870
H      -1.281669     0.000000     7.356885
H       1.199925     0.000000     7.452548
H       2.488547     0.000000     5.300000
H       1.199925     0.000000     3.147452
H      -1.281669     0.000000     3.243115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p4 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.400000
C      -0.698537     0.000000     6.540130
C       0.694182     0.000000     6.595340
C       1.406747     0.000000     5.400000
C       0.694182     0.000000     4.204660
C      -0.698537     0.000000     4.259870
H      -1.281669     0.000000     7.456885
H       1.199925     0.000000     7.552548
H       2.488547     0.000000     5.400000
H       1.199925     0.000000     3.247452
H      -1.281669     0.000000     3.343115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.500000
C      -0.698537     0.000000     6.640130
C       0.694182     0.000000     6.695340
C       1.406747     0.000000     5.500000
C       0.694182     0.000000     4.304660
C      -0.698537     0.000000     4.359870
H      -1.281669     0.000000     7.556885
H       1.199925     0.000000     7.652548
H       2.488547     0.000000     5.500000
H       1.199925     0.000000     3.347452
H      -1.281669     0.000000     3.443115
units angstrom
}
""", 0)

NBC1_PyPy_T3_5p7 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     5.700000
C      -0.698537     0.000000     6.840130
C       0.694182     0.000000     6.895340
C       1.406747     0.000000     5.700000
C       0.694182     0.000000     4.504660
C      -0.698537     0.000000     4.559870
H      -1.281669     0.000000     7.756885
H       1.199925     0.000000     7.852548
H       2.488547     0.000000     5.700000
H       1.199925     0.000000     3.547452
H      -1.281669     0.000000     3.643115
units angstrom
}
""", 0)

NBC1_PyPy_T3_6p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     6.000000
C      -0.698537     0.000000     7.140130
C       0.694182     0.000000     7.195340
C       1.406747     0.000000     6.000000
C       0.694182     0.000000     4.804660
C      -0.698537     0.000000     4.859870
H      -1.281669     0.000000     8.056885
H       1.199925     0.000000     8.152548
H       2.488547     0.000000     6.000000
H       1.199925     0.000000     3.847452
H      -1.281669     0.000000     3.943115
units angstrom
}
""", 0)

NBC1_PyPy_T3_6p5 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     6.500000
C      -0.698537     0.000000     7.640130
C       0.694182     0.000000     7.695340
C       1.406747     0.000000     6.500000
C       0.694182     0.000000     5.304660
C      -0.698537     0.000000     5.359870
H      -1.281669     0.000000     8.556885
H       1.199925     0.000000     8.652548
H       2.488547     0.000000     6.500000
H       1.199925     0.000000     4.347452
H      -1.281669     0.000000     4.443115
units angstrom
}
""", 0)

NBC1_PyPy_T3_7p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     7.000000
C      -0.698537     0.000000     8.140130
C       0.694182     0.000000     8.195340
C       1.406747     0.000000     7.000000
C       0.694182     0.000000     5.804660
C      -0.698537     0.000000     5.859870
H      -1.281669     0.000000     9.056885
H       1.199925     0.000000     9.152548
H       2.488547     0.000000     7.000000
H       1.199925     0.000000     4.847452
H      -1.281669     0.000000     4.943115
units angstrom
}
""", 0)

NBC1_PyPy_T3_8p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     8.000000
C      -0.698537     0.000000     9.140130
C       0.694182     0.000000     9.195340
C       1.406747     0.000000     8.000000
C       0.694182     0.000000     6.804660
C      -0.698537     0.000000     6.859870
H      -1.281669     0.000000    10.056885
H       1.199925     0.000000    10.152548
H       2.488547     0.000000     8.000000
H       1.199925     0.000000     5.847452
H      -1.281669     0.000000     5.943115
units angstrom
}
""", 0)

NBC1_PyPy_T3_9p0 = input.process_input("""
molecule dimer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
--
0 1
N      -1.398038     0.000000     9.000000
C      -0.698537     0.000000    10.140130
C       0.694182     0.000000    10.195340
C       1.406747     0.000000     9.000000
C       0.694182     0.000000     7.804660
C      -0.698537     0.000000     7.859870
H      -1.281669     0.000000    11.056885
H       1.199925     0.000000    11.152548
H       2.488547     0.000000     9.000000
H       1.199925     0.000000     6.847452
H      -1.281669     0.000000     6.943115
units angstrom
}
""", 0)

NBC1_BzBz_PD32_0p2 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.200000     2.391500
C       0.000000     1.005074     1.695750
C       0.000000     1.005074     0.304250
C       0.000000    -0.200000    -0.391500
C       0.000000    -1.405074     0.304250
C       0.000000    -1.405074     1.695750
H       0.000000    -0.200000     3.471500
H       0.000000     1.940382     2.235750
H       0.000000     1.940382    -0.235750
H       0.000000    -0.200000    -1.471500
H       0.000000    -2.340382    -0.235750
H       0.000000    -2.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_0p4 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.400000     2.391500
C       0.000000     0.805074     1.695750
C       0.000000     0.805074     0.304250
C       0.000000    -0.400000    -0.391500
C       0.000000    -1.605074     0.304250
C       0.000000    -1.605074     1.695750
H       0.000000    -0.400000     3.471500
H       0.000000     1.740382     2.235750
H       0.000000     1.740382    -0.235750
H       0.000000    -0.400000    -1.471500
H       0.000000    -2.540382    -0.235750
H       0.000000    -2.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_0p6 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.600000     2.391500
C       0.000000     0.605074     1.695750
C       0.000000     0.605074     0.304250
C       0.000000    -0.600000    -0.391500
C       0.000000    -1.805074     0.304250
C       0.000000    -1.805074     1.695750
H       0.000000    -0.600000     3.471500
H       0.000000     1.540382     2.235750
H       0.000000     1.540382    -0.235750
H       0.000000    -0.600000    -1.471500
H       0.000000    -2.740382    -0.235750
H       0.000000    -2.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_0p8 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.800000     2.391500
C       0.000000     0.405074     1.695750
C       0.000000     0.405074     0.304250
C       0.000000    -0.800000    -0.391500
C       0.000000    -2.005074     0.304250
C       0.000000    -2.005074     1.695750
H       0.000000    -0.800000     3.471500
H       0.000000     1.340382     2.235750
H       0.000000     1.340382    -0.235750
H       0.000000    -0.800000    -1.471500
H       0.000000    -2.940382    -0.235750
H       0.000000    -2.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p0 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.000000     2.391500
C       0.000000     0.205074     1.695750
C       0.000000     0.205074     0.304250
C       0.000000    -1.000000    -0.391500
C       0.000000    -2.205074     0.304250
C       0.000000    -2.205074     1.695750
H       0.000000    -1.000000     3.471500
H       0.000000     1.140382     2.235750
H       0.000000     1.140382    -0.235750
H       0.000000    -1.000000    -1.471500
H       0.000000    -3.140382    -0.235750
H       0.000000    -3.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p2 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.200000     2.391500
C       0.000000     0.005074     1.695750
C       0.000000     0.005074     0.304250
C       0.000000    -1.200000    -0.391500
C       0.000000    -2.405074     0.304250
C       0.000000    -2.405074     1.695750
H       0.000000    -1.200000     3.471500
H       0.000000     0.940382     2.235750
H       0.000000     0.940382    -0.235750
H       0.000000    -1.200000    -1.471500
H       0.000000    -3.340382    -0.235750
H       0.000000    -3.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p4 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.400000     2.391500
C       0.000000    -0.194926     1.695750
C       0.000000    -0.194926     0.304250
C       0.000000    -1.400000    -0.391500
C       0.000000    -2.605074     0.304250
C       0.000000    -2.605074     1.695750
H       0.000000    -1.400000     3.471500
H       0.000000     0.740382     2.235750
H       0.000000     0.740382    -0.235750
H       0.000000    -1.400000    -1.471500
H       0.000000    -3.540382    -0.235750
H       0.000000    -3.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p5 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.500000     2.391500
C       0.000000    -0.294926     1.695750
C       0.000000    -0.294926     0.304250
C       0.000000    -1.500000    -0.391500
C       0.000000    -2.705074     0.304250
C       0.000000    -2.705074     1.695750
H       0.000000    -1.500000     3.471500
H       0.000000     0.640382     2.235750
H       0.000000     0.640382    -0.235750
H       0.000000    -1.500000    -1.471500
H       0.000000    -3.640382    -0.235750
H       0.000000    -3.640382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p6 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.600000     2.391500
C       0.000000    -0.394926     1.695750
C       0.000000    -0.394926     0.304250
C       0.000000    -1.600000    -0.391500
C       0.000000    -2.805074     0.304250
C       0.000000    -2.805074     1.695750
H       0.000000    -1.600000     3.471500
H       0.000000     0.540382     2.235750
H       0.000000     0.540382    -0.235750
H       0.000000    -1.600000    -1.471500
H       0.000000    -3.740382    -0.235750
H       0.000000    -3.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p7 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.700000     2.391500
C       0.000000    -0.494926     1.695750
C       0.000000    -0.494926     0.304250
C       0.000000    -1.700000    -0.391500
C       0.000000    -2.905074     0.304250
C       0.000000    -2.905074     1.695750
H       0.000000    -1.700000     3.471500
H       0.000000     0.440382     2.235750
H       0.000000     0.440382    -0.235750
H       0.000000    -1.700000    -1.471500
H       0.000000    -3.840382    -0.235750
H       0.000000    -3.840382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p8 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.800000     2.391500
C       0.000000    -0.594926     1.695750
C       0.000000    -0.594926     0.304250
C       0.000000    -1.800000    -0.391500
C       0.000000    -3.005074     0.304250
C       0.000000    -3.005074     1.695750
H       0.000000    -1.800000     3.471500
H       0.000000     0.340382     2.235750
H       0.000000     0.340382    -0.235750
H       0.000000    -1.800000    -1.471500
H       0.000000    -3.940382    -0.235750
H       0.000000    -3.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_1p9 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.900000     2.391500
C       0.000000    -0.694926     1.695750
C       0.000000    -0.694926     0.304250
C       0.000000    -1.900000    -0.391500
C       0.000000    -3.105074     0.304250
C       0.000000    -3.105074     1.695750
H       0.000000    -1.900000     3.471500
H       0.000000     0.240382     2.235750
H       0.000000     0.240382    -0.235750
H       0.000000    -1.900000    -1.471500
H       0.000000    -4.040382    -0.235750
H       0.000000    -4.040382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_2p0 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.000000     2.391500
C       0.000000    -0.794926     1.695750
C       0.000000    -0.794926     0.304250
C       0.000000    -2.000000    -0.391500
C       0.000000    -3.205074     0.304250
C       0.000000    -3.205074     1.695750
H       0.000000    -2.000000     3.471500
H       0.000000     0.140382     2.235750
H       0.000000     0.140382    -0.235750
H       0.000000    -2.000000    -1.471500
H       0.000000    -4.140382    -0.235750
H       0.000000    -4.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_2p2 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.200000     2.391500
C       0.000000    -0.994926     1.695750
C       0.000000    -0.994926     0.304250
C       0.000000    -2.200000    -0.391500
C       0.000000    -3.405074     0.304250
C       0.000000    -3.405074     1.695750
H       0.000000    -2.200000     3.471500
H       0.000000    -0.059618     2.235750
H       0.000000    -0.059618    -0.235750
H       0.000000    -2.200000    -1.471500
H       0.000000    -4.340382    -0.235750
H       0.000000    -4.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_2p4 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.400000     2.391500
C       0.000000    -1.194926     1.695750
C       0.000000    -1.194926     0.304250
C       0.000000    -2.400000    -0.391500
C       0.000000    -3.605074     0.304250
C       0.000000    -3.605074     1.695750
H       0.000000    -2.400000     3.471500
H       0.000000    -0.259618     2.235750
H       0.000000    -0.259618    -0.235750
H       0.000000    -2.400000    -1.471500
H       0.000000    -4.540382    -0.235750
H       0.000000    -4.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_2p6 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.600000     2.391500
C       0.000000    -1.394926     1.695750
C       0.000000    -1.394926     0.304250
C       0.000000    -2.600000    -0.391500
C       0.000000    -3.805074     0.304250
C       0.000000    -3.805074     1.695750
H       0.000000    -2.600000     3.471500
H       0.000000    -0.459618     2.235750
H       0.000000    -0.459618    -0.235750
H       0.000000    -2.600000    -1.471500
H       0.000000    -4.740382    -0.235750
H       0.000000    -4.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_2p8 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.800000     2.391500
C       0.000000    -1.594926     1.695750
C       0.000000    -1.594926     0.304250
C       0.000000    -2.800000    -0.391500
C       0.000000    -4.005074     0.304250
C       0.000000    -4.005074     1.695750
H       0.000000    -2.800000     3.471500
H       0.000000    -0.659618     2.235750
H       0.000000    -0.659618    -0.235750
H       0.000000    -2.800000    -1.471500
H       0.000000    -4.940382    -0.235750
H       0.000000    -4.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD32_3p0 = input.process_input("""
molecule dimer {
0 1
C       3.200000     0.000000    -0.391500
C       3.200000     1.205074     0.304250
C       3.200000     1.205074     1.695750
C       3.200000     0.000000     2.391500
C       3.200000    -1.205074     1.695750
C       3.200000    -1.205074     0.304250
H       3.200000     0.000000    -1.471500
H       3.200000     2.140382    -0.235750
H       3.200000     2.140382     2.235750
H       3.200000     0.000000     3.471500
H       3.200000    -2.140382     2.235750
H       3.200000    -2.140382    -0.235750
--
0 1
C       0.000000    -3.000000     2.391500
C       0.000000    -1.794926     1.695750
C       0.000000    -1.794926     0.304250
C       0.000000    -3.000000    -0.391500
C       0.000000    -4.205074     0.304250
C       0.000000    -4.205074     1.695750
H       0.000000    -3.000000     3.471500
H       0.000000    -0.859618     2.235750
H       0.000000    -0.859618    -0.235750
H       0.000000    -3.000000    -1.471500
H       0.000000    -5.140382    -0.235750
H       0.000000    -5.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_0p2 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.200000     2.391500
C       0.000000     1.005074     1.695750
C       0.000000     1.005074     0.304250
C       0.000000    -0.200000    -0.391500
C       0.000000    -1.405074     0.304250
C       0.000000    -1.405074     1.695750
H       0.000000    -0.200000     3.471500
H       0.000000     1.940382     2.235750
H       0.000000     1.940382    -0.235750
H       0.000000    -0.200000    -1.471500
H       0.000000    -2.340382    -0.235750
H       0.000000    -2.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_0p4 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.400000     2.391500
C       0.000000     0.805074     1.695750
C       0.000000     0.805074     0.304250
C       0.000000    -0.400000    -0.391500
C       0.000000    -1.605074     0.304250
C       0.000000    -1.605074     1.695750
H       0.000000    -0.400000     3.471500
H       0.000000     1.740382     2.235750
H       0.000000     1.740382    -0.235750
H       0.000000    -0.400000    -1.471500
H       0.000000    -2.540382    -0.235750
H       0.000000    -2.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_0p6 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.600000     2.391500
C       0.000000     0.605074     1.695750
C       0.000000     0.605074     0.304250
C       0.000000    -0.600000    -0.391500
C       0.000000    -1.805074     0.304250
C       0.000000    -1.805074     1.695750
H       0.000000    -0.600000     3.471500
H       0.000000     1.540382     2.235750
H       0.000000     1.540382    -0.235750
H       0.000000    -0.600000    -1.471500
H       0.000000    -2.740382    -0.235750
H       0.000000    -2.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_0p8 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -0.800000     2.391500
C       0.000000     0.405074     1.695750
C       0.000000     0.405074     0.304250
C       0.000000    -0.800000    -0.391500
C       0.000000    -2.005074     0.304250
C       0.000000    -2.005074     1.695750
H       0.000000    -0.800000     3.471500
H       0.000000     1.340382     2.235750
H       0.000000     1.340382    -0.235750
H       0.000000    -0.800000    -1.471500
H       0.000000    -2.940382    -0.235750
H       0.000000    -2.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p0 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.000000     2.391500
C       0.000000     0.205074     1.695750
C       0.000000     0.205074     0.304250
C       0.000000    -1.000000    -0.391500
C       0.000000    -2.205074     0.304250
C       0.000000    -2.205074     1.695750
H       0.000000    -1.000000     3.471500
H       0.000000     1.140382     2.235750
H       0.000000     1.140382    -0.235750
H       0.000000    -1.000000    -1.471500
H       0.000000    -3.140382    -0.235750
H       0.000000    -3.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p2 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.200000     2.391500
C       0.000000     0.005074     1.695750
C       0.000000     0.005074     0.304250
C       0.000000    -1.200000    -0.391500
C       0.000000    -2.405074     0.304250
C       0.000000    -2.405074     1.695750
H       0.000000    -1.200000     3.471500
H       0.000000     0.940382     2.235750
H       0.000000     0.940382    -0.235750
H       0.000000    -1.200000    -1.471500
H       0.000000    -3.340382    -0.235750
H       0.000000    -3.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p4 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.400000     2.391500
C       0.000000    -0.194926     1.695750
C       0.000000    -0.194926     0.304250
C       0.000000    -1.400000    -0.391500
C       0.000000    -2.605074     0.304250
C       0.000000    -2.605074     1.695750
H       0.000000    -1.400000     3.471500
H       0.000000     0.740382     2.235750
H       0.000000     0.740382    -0.235750
H       0.000000    -1.400000    -1.471500
H       0.000000    -3.540382    -0.235750
H       0.000000    -3.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p5 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.500000     2.391500
C       0.000000    -0.294926     1.695750
C       0.000000    -0.294926     0.304250
C       0.000000    -1.500000    -0.391500
C       0.000000    -2.705074     0.304250
C       0.000000    -2.705074     1.695750
H       0.000000    -1.500000     3.471500
H       0.000000     0.640382     2.235750
H       0.000000     0.640382    -0.235750
H       0.000000    -1.500000    -1.471500
H       0.000000    -3.640382    -0.235750
H       0.000000    -3.640382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p6 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.600000     2.391500
C       0.000000    -0.394926     1.695750
C       0.000000    -0.394926     0.304250
C       0.000000    -1.600000    -0.391500
C       0.000000    -2.805074     0.304250
C       0.000000    -2.805074     1.695750
H       0.000000    -1.600000     3.471500
H       0.000000     0.540382     2.235750
H       0.000000     0.540382    -0.235750
H       0.000000    -1.600000    -1.471500
H       0.000000    -3.740382    -0.235750
H       0.000000    -3.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p7 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.700000     2.391500
C       0.000000    -0.494926     1.695750
C       0.000000    -0.494926     0.304250
C       0.000000    -1.700000    -0.391500
C       0.000000    -2.905074     0.304250
C       0.000000    -2.905074     1.695750
H       0.000000    -1.700000     3.471500
H       0.000000     0.440382     2.235750
H       0.000000     0.440382    -0.235750
H       0.000000    -1.700000    -1.471500
H       0.000000    -3.840382    -0.235750
H       0.000000    -3.840382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p8 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.800000     2.391500
C       0.000000    -0.594926     1.695750
C       0.000000    -0.594926     0.304250
C       0.000000    -1.800000    -0.391500
C       0.000000    -3.005074     0.304250
C       0.000000    -3.005074     1.695750
H       0.000000    -1.800000     3.471500
H       0.000000     0.340382     2.235750
H       0.000000     0.340382    -0.235750
H       0.000000    -1.800000    -1.471500
H       0.000000    -3.940382    -0.235750
H       0.000000    -3.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_1p9 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -1.900000     2.391500
C       0.000000    -0.694926     1.695750
C       0.000000    -0.694926     0.304250
C       0.000000    -1.900000    -0.391500
C       0.000000    -3.105074     0.304250
C       0.000000    -3.105074     1.695750
H       0.000000    -1.900000     3.471500
H       0.000000     0.240382     2.235750
H       0.000000     0.240382    -0.235750
H       0.000000    -1.900000    -1.471500
H       0.000000    -4.040382    -0.235750
H       0.000000    -4.040382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_2p0 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.000000     2.391500
C       0.000000    -0.794926     1.695750
C       0.000000    -0.794926     0.304250
C       0.000000    -2.000000    -0.391500
C       0.000000    -3.205074     0.304250
C       0.000000    -3.205074     1.695750
H       0.000000    -2.000000     3.471500
H       0.000000     0.140382     2.235750
H       0.000000     0.140382    -0.235750
H       0.000000    -2.000000    -1.471500
H       0.000000    -4.140382    -0.235750
H       0.000000    -4.140382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_2p2 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.200000     2.391500
C       0.000000    -0.994926     1.695750
C       0.000000    -0.994926     0.304250
C       0.000000    -2.200000    -0.391500
C       0.000000    -3.405074     0.304250
C       0.000000    -3.405074     1.695750
H       0.000000    -2.200000     3.471500
H       0.000000    -0.059618     2.235750
H       0.000000    -0.059618    -0.235750
H       0.000000    -2.200000    -1.471500
H       0.000000    -4.340382    -0.235750
H       0.000000    -4.340382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_2p4 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.400000     2.391500
C       0.000000    -1.194926     1.695750
C       0.000000    -1.194926     0.304250
C       0.000000    -2.400000    -0.391500
C       0.000000    -3.605074     0.304250
C       0.000000    -3.605074     1.695750
H       0.000000    -2.400000     3.471500
H       0.000000    -0.259618     2.235750
H       0.000000    -0.259618    -0.235750
H       0.000000    -2.400000    -1.471500
H       0.000000    -4.540382    -0.235750
H       0.000000    -4.540382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_2p6 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.600000     2.391500
C       0.000000    -1.394926     1.695750
C       0.000000    -1.394926     0.304250
C       0.000000    -2.600000    -0.391500
C       0.000000    -3.805074     0.304250
C       0.000000    -3.805074     1.695750
H       0.000000    -2.600000     3.471500
H       0.000000    -0.459618     2.235750
H       0.000000    -0.459618    -0.235750
H       0.000000    -2.600000    -1.471500
H       0.000000    -4.740382    -0.235750
H       0.000000    -4.740382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_2p8 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -2.800000     2.391500
C       0.000000    -1.594926     1.695750
C       0.000000    -1.594926     0.304250
C       0.000000    -2.800000    -0.391500
C       0.000000    -4.005074     0.304250
C       0.000000    -4.005074     1.695750
H       0.000000    -2.800000     3.471500
H       0.000000    -0.659618     2.235750
H       0.000000    -0.659618    -0.235750
H       0.000000    -2.800000    -1.471500
H       0.000000    -4.940382    -0.235750
H       0.000000    -4.940382     2.235750
units angstrom
}
""", 0)

NBC1_BzBz_PD36_3p0 = input.process_input("""
molecule dimer {
0 1
C       3.600000     0.000000    -0.391500
C       3.600000     1.205074     0.304250
C       3.600000     1.205074     1.695750
C       3.600000     0.000000     2.391500
C       3.600000    -1.205074     1.695750
C       3.600000    -1.205074     0.304250
H       3.600000     0.000000    -1.471500
H       3.600000     2.140382    -0.235750
H       3.600000     2.140382     2.235750
H       3.600000     0.000000     3.471500
H       3.600000    -2.140382     2.235750
H       3.600000    -2.140382    -0.235750
--
0 1
C       0.000000    -3.000000     2.391500
C       0.000000    -1.794926     1.695750
C       0.000000    -1.794926     0.304250
C       0.000000    -3.000000    -0.391500
C       0.000000    -4.205074     0.304250
C       0.000000    -4.205074     1.695750
H       0.000000    -3.000000     3.471500
H       0.000000    -0.859618     2.235750
H       0.000000    -0.859618    -0.235750
H       0.000000    -3.000000    -1.471500
H       0.000000    -5.140382    -0.235750
H       0.000000    -5.140382     2.235750
units angstrom
}
""", 0)

NBC1_Bz_monomer = input.process_input("""
molecule monomer {
0 1
C       1.000000     0.000000    -0.391500
C       1.000000     1.205074     0.304250
C       1.000000     1.205074     1.695750
C       1.000000     0.000000     2.391500
C       1.000000    -1.205074     1.695750
C       1.000000    -1.205074     0.304250
H       1.000000     0.000000    -1.471500
H       1.000000     2.140382    -0.235750
H       1.000000     2.140382     2.235750
H       1.000000     0.000000     3.471500
H       1.000000    -2.140382     2.235750
H       1.000000    -2.140382    -0.235750
units angstrom
}
""", 0)

NBC1_H2S_monomer = input.process_input("""
molecule monomer {
0 1
S       0.000000     0.000000     0.000000
H      -0.961721     0.000000     0.926779
H       0.961721     0.000000     0.926779
units angstrom
}
""", 0)

NBC1_Bz2_monomer = input.process_input("""
molecule monomer {
0 1
C       1.405731     0.000000     0.000000
C       0.702865     1.217399     0.000000
C      -0.702865     1.217399     0.000000
C      -1.405731     0.000000     0.000000
C      -0.702865    -1.217399     0.000000
C       0.702866    -1.217399     0.000000
H       2.500941     0.000000     0.000000
H       1.250471     2.165878     0.000000
H      -1.250470     2.165878     0.000000
H      -2.500941     0.000000     0.000000
H      -1.250470    -2.165878     0.000000
H       1.250471    -2.165878     0.000000
units angstrom
}
""", 0)

NBC1_Me_monomer = input.process_input("""
molecule monomer {
0 1
C       0.000000     0.000000     0.000000
H       0.000000     0.000000     1.099503
H       1.036621     0.000000    -0.366501
H      -0.518311    -0.897741    -0.366501
H      -0.518311     0.897741    -0.366501
units angstrom
}
""", 0)

NBC1_Py_monomer = input.process_input("""
molecule monomer {
0 1
N       1.398038     0.000000     0.000000
C       0.698537     1.140130     0.000000
C      -0.694182     1.195340     0.000000
C      -1.406747     0.000000     0.000000
C      -0.694182    -1.195340     0.000000
C       0.698537    -1.140130     0.000000
H       1.281669     2.056885     0.000000
H      -1.199925     2.152548     0.000000
H      -2.488547     0.000000     0.000000
H      -1.199925    -2.152548     0.000000
H       1.281669    -2.056885     0.000000
units angstrom
}
""", 0)

#<<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:
    distance = rxnpattern.match(rxn)

    GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))
    GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoA_CP
    GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoB_CP

GEOS['%s-Bz-mono-unCP'  % (dbse)] = eval('%s_Bz_monomer'  % (dbse))
GEOS['%s-H2S-mono-unCP' % (dbse)] = eval('%s_H2S_monomer' % (dbse))
GEOS['%s-Bz2-mono-unCP' % (dbse)] = eval('%s_Bz2_monomer' % (dbse))
GEOS['%s-Me-mono-unCP'  % (dbse)] = eval('%s_Me_monomer'  % (dbse))
GEOS['%s-Py-mono-unCP'  % (dbse)] = eval('%s_Py_monomer'  % (dbse))
