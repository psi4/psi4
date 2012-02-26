"""
**HBC6**

| Database (Sherrill) of interaction energies for dissociation curves of doubly hydrogen-bonded bimolecular complexes.
| Geometries from Thanthiriwatte et al. JCTC 7 88 (2011).
| Reference interaction energies from Marshall et al. JCP 135 194102 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'`` || ``'on'``

- **benchmark**

  - ``'HBC60'`` Thanthiriwatte et al. JCTC 7 88 (2011).
  - |dl| ``'HBC6A'`` |dr| Marshall et al. JCP 135 194102 (2011).
  - ``'HBC6ARLX'`` Sherrill group, unpublished.

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``

----

"""
import re
import input

# <<< HBC6 Database Module >>>
dbse = 'HBC1'

# <<< Database Members >>>
FaOOFaOO = []
FaONFaON = []
FaNNFaNN = []
FaOOFaON = []
FaONFaNN = []
FaOOFaNN = []
dist = [3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.6,4.8,5.0,5.4,5.8,6.4,7.0,8.0,10.0]
for d in dist: FaOOFaOO.append('FaOOFaOO-' + str(d))
for d in dist: FaONFaON.append('FaONFaON-' + str(d))
for d in dist: FaNNFaNN.append('FaNNFaNN-' + str(d))
for d in dist: FaOOFaON.append('FaOOFaON-' + str(d))
for d in dist: FaONFaNN.append('FaONFaNN-' + str(d))
dist = [3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.6,4.8,5.0,5.4,5.8,6.4,7.0,8.0,10.0]
for d in dist: FaOOFaNN.append('FaOOFaNN-' + str(d))

temp = [FaOOFaOO, FaONFaON, FaNNFaNN, FaOOFaON, FaONFaNN, FaOOFaNN]
HRXN = sum(temp, [])

HRXN_SM = ['FaOOFaOO-8.0','FaOOFaON-5.0']
HRXN_LG = ['FaNNFaNN-3.6']
HRXN_EQ = ['FaOOFaOO-3.6', 'FaONFaON-4.0', 'FaNNFaNN-4.1', 'FaOOFaON-3.8', 'FaONFaNN-4.0', 'FaOOFaNN-3.6']

# <<< Chemical Systems Involved >>>
RXNM = {}        # reaction matrix of reagent contributions per reaction
RXNM_CPRLX = {}  # reaction matrix of reagent contributions per reaction for counterpoise- and deformation-corrected
ACTV = {}        # order of active reagents per reaction
ACTV_CP = {}     # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}     # order of active reagents for non-supramolecular calculations
ACTV_RLX = {}    # order of active reagents for deformation-corrected reaction
ACTV_CPRLX = {}  # order of active reagents for counterpoise- and deformation-corrected reaction
monopattern = re.compile(r'^(....)(....)-(.+)$')
for rxn in HRXN:
   molname = monopattern.match(rxn)

   if (rxn in FaOOFaOO) or (rxn in FaONFaON) or (rxn in FaNNFaNN):
      RXNM[      '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : -2,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) : -2 }

      RXNM_CPRLX['%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : +2,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) : -2 }

      ACTV_SA[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

      ACTV_CP[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn) ]

      ACTV_RLX[  '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) ]

      ACTV_CPRLX['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) ]

      ACTV[      '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn) ]

   elif (rxn in FaOOFaON) or (rxn in FaONFaNN) or (rxn in FaOOFaNN):
      RXNM[      '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                           '%s-%s-monoB-unCP' % (dbse, rxn) : -1,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) : -1,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(2)) : -1 }

      RXNM_CPRLX['%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : +1,
                                           '%s-%s-monoB-unCP' % (dbse, rxn) : +1,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)) : -1,
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(2)) : -1 }

      ACTV_SA[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

      ACTV_CP[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn),
                                           '%s-%s-monoB-CP'   % (dbse, rxn) ]

      ACTV_RLX[  '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(2)) ]

      ACTV_CPRLX['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn),
                                           '%s-%s-monoB-CP'   % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn),
                                           '%s-%s-monoB-unCP' % (dbse, rxn),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(1)),
                                           '%s-%s-mono-RLX'   % (dbse, molname.group(2)) ]

      ACTV[      '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn),
                                           '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values >>>
BIND = {}
# Original publication
BIND_HBC60 = {}
BIND_HBC60['%s-FaOOFaOO-3.4'  % (dbse)] = -19.834
BIND_HBC60['%s-FaOOFaOO-3.5'  % (dbse)] = -20.027
BIND_HBC60['%s-FaOOFaOO-3.6'  % (dbse)] = -20.060  # FaOOFaOO minimum
BIND_HBC60['%s-FaOOFaOO-3.7'  % (dbse)] = -19.776
BIND_HBC60['%s-FaOOFaOO-3.8'  % (dbse)] = -19.132
BIND_HBC60['%s-FaOOFaOO-3.9'  % (dbse)] = -18.161
BIND_HBC60['%s-FaOOFaOO-4.0'  % (dbse)] = -16.943
BIND_HBC60['%s-FaOOFaOO-4.1'  % (dbse)] = -15.574
BIND_HBC60['%s-FaOOFaOO-4.2'  % (dbse)] = -14.148
BIND_HBC60['%s-FaOOFaOO-4.3'  % (dbse)] = -12.736
BIND_HBC60['%s-FaOOFaOO-4.4'  % (dbse)] = -11.392
BIND_HBC60['%s-FaOOFaOO-4.6'  % (dbse)] =  -9.014
BIND_HBC60['%s-FaOOFaOO-4.8'  % (dbse)] =  -7.091
BIND_HBC60['%s-FaOOFaOO-5.0'  % (dbse)] =  -5.590
BIND_HBC60['%s-FaOOFaOO-5.4'  % (dbse)] =  -3.548
BIND_HBC60['%s-FaOOFaOO-5.8'  % (dbse)] =  -2.325
BIND_HBC60['%s-FaOOFaOO-6.4'  % (dbse)] =  -1.320
BIND_HBC60['%s-FaOOFaOO-7.0'  % (dbse)] =  -0.801
BIND_HBC60['%s-FaOOFaOO-8.0'  % (dbse)] =  -0.394
BIND_HBC60['%s-FaOOFaOO-10.0' % (dbse)] =  -0.132

BIND_HBC60['%s-FaONFaON-3.4'  % (dbse)] =  -6.726
BIND_HBC60['%s-FaONFaON-3.5'  % (dbse)] = -10.191
BIND_HBC60['%s-FaONFaON-3.6'  % (dbse)] = -12.781
BIND_HBC60['%s-FaONFaON-3.7'  % (dbse)] = -14.667
BIND_HBC60['%s-FaONFaON-3.8'  % (dbse)] = -15.919
BIND_HBC60['%s-FaONFaON-3.9'  % (dbse)] = -16.582
BIND_HBC60['%s-FaONFaON-4.0'  % (dbse)] = -16.714  # FaONFaON minimum
BIND_HBC60['%s-FaONFaON-4.1'  % (dbse)] = -16.391
BIND_HBC60['%s-FaONFaON-4.2'  % (dbse)] = -15.713
BIND_HBC60['%s-FaONFaON-4.3'  % (dbse)] = -14.790
BIND_HBC60['%s-FaONFaON-4.4'  % (dbse)] = -13.723
BIND_HBC60['%s-FaONFaON-4.6'  % (dbse)] = -11.480
BIND_HBC60['%s-FaONFaON-4.8'  % (dbse)] =  -9.401
BIND_HBC60['%s-FaONFaON-5.0'  % (dbse)] =  -7.642
BIND_HBC60['%s-FaONFaON-5.4'  % (dbse)] =  -5.108
BIND_HBC60['%s-FaONFaON-5.8'  % (dbse)] =  -3.537
BIND_HBC60['%s-FaONFaON-6.4'  % (dbse)] =  -2.187
BIND_HBC60['%s-FaONFaON-7.0'  % (dbse)] =  -1.448
BIND_HBC60['%s-FaONFaON-8.0'  % (dbse)] =  -0.816
BIND_HBC60['%s-FaONFaON-10.0' % (dbse)] =  -0.340

BIND_HBC60['%s-FaNNFaNN-3.4'  % (dbse)] =  -8.987
BIND_HBC60['%s-FaNNFaNN-3.5'  % (dbse)] = -10.969
BIND_HBC60['%s-FaNNFaNN-3.6'  % (dbse)] = -12.693
BIND_HBC60['%s-FaNNFaNN-3.7'  % (dbse)] = -14.144
BIND_HBC60['%s-FaNNFaNN-3.8'  % (dbse)] = -15.287
BIND_HBC60['%s-FaNNFaNN-3.9'  % (dbse)] = -16.118
BIND_HBC60['%s-FaNNFaNN-4.0'  % (dbse)] = -16.587
BIND_HBC60['%s-FaNNFaNN-4.1'  % (dbse)] = -16.702  # FaNNFaNN minimum
BIND_HBC60['%s-FaNNFaNN-4.2'  % (dbse)] = -16.452
BIND_HBC60['%s-FaNNFaNN-4.3'  % (dbse)] = -15.901
BIND_HBC60['%s-FaNNFaNN-4.4'  % (dbse)] = -15.102
BIND_HBC60['%s-FaNNFaNN-4.6'  % (dbse)] = -13.047
BIND_HBC60['%s-FaNNFaNN-4.8'  % (dbse)] = -10.810
BIND_HBC60['%s-FaNNFaNN-5.0'  % (dbse)] =  -8.733
BIND_HBC60['%s-FaNNFaNN-5.4'  % (dbse)] =  -5.539
BIND_HBC60['%s-FaNNFaNN-5.8'  % (dbse)] =  -3.521
BIND_HBC60['%s-FaNNFaNN-6.4'  % (dbse)] =  -1.861
BIND_HBC60['%s-FaNNFaNN-7.0'  % (dbse)] =  -1.050
BIND_HBC60['%s-FaNNFaNN-8.0'  % (dbse)] =  -0.463
BIND_HBC60['%s-FaNNFaNN-10.0' % (dbse)] =  -0.123

BIND_HBC60['%s-FaOOFaON-3.4'  % (dbse)] = -14.356
BIND_HBC60['%s-FaOOFaON-3.5'  % (dbse)] = -16.486
BIND_HBC60['%s-FaOOFaON-3.6'  % (dbse)] = -17.833
BIND_HBC60['%s-FaOOFaON-3.7'  % (dbse)] = -18.543
BIND_HBC60['%s-FaOOFaON-3.8'  % (dbse)] = -18.692  # FaOOFaON minimum
BIND_HBC60['%s-FaOOFaON-3.9'  % (dbse)] = -18.347
BIND_HBC60['%s-FaOOFaON-4.0'  % (dbse)] = -17.592
BIND_HBC60['%s-FaOOFaON-4.1'  % (dbse)] = -16.537
BIND_HBC60['%s-FaOOFaON-4.2'  % (dbse)] = -15.300
BIND_HBC60['%s-FaOOFaON-4.3'  % (dbse)] = -13.989
BIND_HBC60['%s-FaOOFaON-4.4'  % (dbse)] = -12.684
BIND_HBC60['%s-FaOOFaON-4.6'  % (dbse)] = -10.274
BIND_HBC60['%s-FaOOFaON-4.8'  % (dbse)] =  -8.245
BIND_HBC60['%s-FaOOFaON-5.0'  % (dbse)] =  -6.613
BIND_HBC60['%s-FaOOFaON-5.4'  % (dbse)] =  -4.330
BIND_HBC60['%s-FaOOFaON-5.8'  % (dbse)] =  -2.935
BIND_HBC60['%s-FaOOFaON-6.4'  % (dbse)] =  -1.753
BIND_HBC60['%s-FaOOFaON-7.0'  % (dbse)] =  -1.121
BIND_HBC60['%s-FaOOFaON-8.0'  % (dbse)] =  -0.598
BIND_HBC60['%s-FaOOFaON-10.0' % (dbse)] =  -0.227

BIND_HBC60['%s-FaONFaNN-3.4'  % (dbse)] =  -8.239
BIND_HBC60['%s-FaONFaNN-3.5'  % (dbse)] = -10.918
BIND_HBC60['%s-FaONFaNN-3.6'  % (dbse)] = -13.055
BIND_HBC60['%s-FaONFaNN-3.7'  % (dbse)] = -14.717
BIND_HBC60['%s-FaONFaNN-3.8'  % (dbse)] = -15.921
BIND_HBC60['%s-FaONFaNN-3.9'  % (dbse)] = -16.672
BIND_HBC60['%s-FaONFaNN-4.0'  % (dbse)] = -16.977  # FaONFaNN minimum
BIND_HBC60['%s-FaONFaNN-4.1'  % (dbse)] = -16.865
BIND_HBC60['%s-FaONFaNN-4.2'  % (dbse)] = -16.390
BIND_HBC60['%s-FaONFaNN-4.3'  % (dbse)] = -15.631
BIND_HBC60['%s-FaONFaNN-4.4'  % (dbse)] = -14.676
BIND_HBC60['%s-FaONFaNN-4.6'  % (dbse)] = -12.490
BIND_HBC60['%s-FaONFaNN-4.8'  % (dbse)] = -10.304
BIND_HBC60['%s-FaONFaNN-5.0'  % (dbse)] =  -8.362
BIND_HBC60['%s-FaONFaNN-5.4'  % (dbse)] =  -5.445
BIND_HBC60['%s-FaONFaNN-5.8'  % (dbse)] =  -3.617
BIND_HBC60['%s-FaONFaNN-6.4'  % (dbse)] =  -2.087
BIND_HBC60['%s-FaONFaNN-7.0'  % (dbse)] =  -1.295
BIND_HBC60['%s-FaONFaNN-8.0'  % (dbse)] =  -0.663
BIND_HBC60['%s-FaONFaNN-10.0' % (dbse)] =  -0.237

BIND_HBC60['%s-FaOOFaNN-3.6'  % (dbse)] = -26.289  # FaNNFaNN minimum
BIND_HBC60['%s-FaOOFaNN-3.7'  % (dbse)] = -24.035
BIND_HBC60['%s-FaOOFaNN-3.8'  % (dbse)] = -23.017
BIND_HBC60['%s-FaOOFaNN-3.9'  % (dbse)] = -22.133
BIND_HBC60['%s-FaOOFaNN-4.0'  % (dbse)] = -21.122
BIND_HBC60['%s-FaOOFaNN-4.1'  % (dbse)] = -19.920
BIND_HBC60['%s-FaOOFaNN-4.2'  % (dbse)] = -18.544
BIND_HBC60['%s-FaOOFaNN-4.3'  % (dbse)] = -17.056
BIND_HBC60['%s-FaOOFaNN-4.4'  % (dbse)] = -15.526
BIND_HBC60['%s-FaOOFaNN-4.6'  % (dbse)] = -12.583
BIND_HBC60['%s-FaOOFaNN-4.8'  % (dbse)] = -10.031
BIND_HBC60['%s-FaOOFaNN-5.0'  % (dbse)] =  -7.960
BIND_HBC60['%s-FaOOFaNN-5.4'  % (dbse)] =  -5.069
BIND_HBC60['%s-FaOOFaNN-5.8'  % (dbse)] =  -3.336
BIND_HBC60['%s-FaOOFaNN-6.4'  % (dbse)] =  -1.906
BIND_HBC60['%s-FaOOFaNN-7.0'  % (dbse)] =  -1.170
BIND_HBC60['%s-FaOOFaNN-8.0'  % (dbse)] =  -0.587
BIND_HBC60['%s-FaOOFaNN-10.0' % (dbse)] =  -0.202
# Current revision
BIND_HBC6A = {}
BIND_HBC6A['%s-FaOOFaOO-3.4'  % (dbse)] = -19.627
BIND_HBC6A['%s-FaOOFaOO-3.5'  % (dbse)] = -19.850
BIND_HBC6A['%s-FaOOFaOO-3.6'  % (dbse)] = -19.910  # FaOOFaOO minimum
BIND_HBC6A['%s-FaOOFaOO-3.7'  % (dbse)] = -19.650
BIND_HBC6A['%s-FaOOFaOO-3.8'  % (dbse)] = -19.027
BIND_HBC6A['%s-FaOOFaOO-3.9'  % (dbse)] = -18.075
BIND_HBC6A['%s-FaOOFaOO-4.0'  % (dbse)] = -16.873
BIND_HBC6A['%s-FaOOFaOO-4.1'  % (dbse)] = -15.517
BIND_HBC6A['%s-FaOOFaOO-4.2'  % (dbse)] = -14.100
BIND_HBC6A['%s-FaOOFaOO-4.3'  % (dbse)] = -12.697
BIND_HBC6A['%s-FaOOFaOO-4.4'  % (dbse)] = -11.360
BIND_HBC6A['%s-FaOOFaOO-4.6'  % (dbse)] =  -8.990
BIND_HBC6A['%s-FaOOFaOO-4.8'  % (dbse)] =  -7.074
BIND_HBC6A['%s-FaOOFaOO-5.0'  % (dbse)] =  -5.577
BIND_HBC6A['%s-FaOOFaOO-5.4'  % (dbse)] =  -3.539
BIND_HBC6A['%s-FaOOFaOO-5.8'  % (dbse)] =  -2.323
BIND_HBC6A['%s-FaOOFaOO-6.4'  % (dbse)] =  -1.320
BIND_HBC6A['%s-FaOOFaOO-7.0'  % (dbse)] =  -0.802
BIND_HBC6A['%s-FaOOFaOO-8.0'  % (dbse)] =  -0.397
BIND_HBC6A['%s-FaOOFaOO-10.0' % (dbse)] =  -0.135

BIND_HBC6A['%s-FaONFaON-3.4'  % (dbse)] =  -6.556
BIND_HBC6A['%s-FaONFaON-3.5'  % (dbse)] = -10.027
BIND_HBC6A['%s-FaONFaON-3.6'  % (dbse)] = -12.628
BIND_HBC6A['%s-FaONFaON-3.7'  % (dbse)] = -14.529
BIND_HBC6A['%s-FaONFaON-3.8'  % (dbse)] = -15.796
BIND_HBC6A['%s-FaONFaON-3.9'  % (dbse)] = -16.475
BIND_HBC6A['%s-FaONFaON-4.0'  % (dbse)] = -16.622  # FaONFaON minimum
BIND_HBC6A['%s-FaONFaON-4.1'  % (dbse)] = -16.313
BIND_HBC6A['%s-FaONFaON-4.2'  % (dbse)] = -15.647
BIND_HBC6A['%s-FaONFaON-4.3'  % (dbse)] = -14.735
BIND_HBC6A['%s-FaONFaON-4.4'  % (dbse)] = -13.678
BIND_HBC6A['%s-FaONFaON-4.6'  % (dbse)] = -11.448
BIND_HBC6A['%s-FaONFaON-4.8'  % (dbse)] =  -9.379
BIND_HBC6A['%s-FaONFaON-5.0'  % (dbse)] =  -7.626
BIND_HBC6A['%s-FaONFaON-5.4'  % (dbse)] =  -5.097
BIND_HBC6A['%s-FaONFaON-5.8'  % (dbse)] =  -3.528
BIND_HBC6A['%s-FaONFaON-6.4'  % (dbse)] =  -2.181
BIND_HBC6A['%s-FaONFaON-7.0'  % (dbse)] =  -1.443
BIND_HBC6A['%s-FaONFaON-8.0'  % (dbse)] =  -0.813
BIND_HBC6A['%s-FaONFaON-10.0' % (dbse)] =  -0.337

BIND_HBC6A['%s-FaNNFaNN-3.4'  % (dbse)] =  -8.730
BIND_HBC6A['%s-FaNNFaNN-3.5'  % (dbse)] = -10.725
BIND_HBC6A['%s-FaNNFaNN-3.6'  % (dbse)] = -12.463
BIND_HBC6A['%s-FaNNFaNN-3.7'  % (dbse)] = -13.932
BIND_HBC6A['%s-FaNNFaNN-3.8'  % (dbse)] = -15.106
BIND_HBC6A['%s-FaNNFaNN-3.9'  % (dbse)] = -15.950
BIND_HBC6A['%s-FaNNFaNN-4.0'  % (dbse)] = -16.440
BIND_HBC6A['%s-FaNNFaNN-4.1'  % (dbse)] = -16.575  # FaNNFaNN minimum
BIND_HBC6A['%s-FaNNFaNN-4.2'  % (dbse)] = -16.344
BIND_HBC6A['%s-FaNNFaNN-4.3'  % (dbse)] = -15.811
BIND_HBC6A['%s-FaNNFaNN-4.4'  % (dbse)] = -15.028
BIND_HBC6A['%s-FaNNFaNN-4.6'  % (dbse)] = -12.999
BIND_HBC6A['%s-FaNNFaNN-4.8'  % (dbse)] = -10.780
BIND_HBC6A['%s-FaNNFaNN-5.0'  % (dbse)] =  -8.715
BIND_HBC6A['%s-FaNNFaNN-5.4'  % (dbse)] =  -5.532
BIND_HBC6A['%s-FaNNFaNN-5.8'  % (dbse)] =  -3.517
BIND_HBC6A['%s-FaNNFaNN-6.4'  % (dbse)] =  -1.861
BIND_HBC6A['%s-FaNNFaNN-7.0'  % (dbse)] =  -1.051
BIND_HBC6A['%s-FaNNFaNN-8.0'  % (dbse)] =  -0.466
BIND_HBC6A['%s-FaNNFaNN-10.0' % (dbse)] =  -0.127

BIND_HBC6A['%s-FaOOFaON-3.4'  % (dbse)] = -14.164
BIND_HBC6A['%s-FaOOFaON-3.5'  % (dbse)] = -16.312
BIND_HBC6A['%s-FaOOFaON-3.6'  % (dbse)] = -17.679
BIND_HBC6A['%s-FaOOFaON-3.7'  % (dbse)] = -18.409
BIND_HBC6A['%s-FaOOFaON-3.8'  % (dbse)] = -18.578  # FaOOFaON minimum
BIND_HBC6A['%s-FaOOFaON-3.9'  % (dbse)] = -18.250
BIND_HBC6A['%s-FaOOFaON-4.0'  % (dbse)] = -17.512
BIND_HBC6A['%s-FaOOFaON-4.1'  % (dbse)] = -16.471
BIND_HBC6A['%s-FaOOFaON-4.2'  % (dbse)] = -15.245
BIND_HBC6A['%s-FaOOFaON-4.3'  % (dbse)] = -13.944
BIND_HBC6A['%s-FaOOFaON-4.4'  % (dbse)] = -12.647
BIND_HBC6A['%s-FaOOFaON-4.6'  % (dbse)] = -10.248
BIND_HBC6A['%s-FaOOFaON-4.8'  % (dbse)] =  -8.227
BIND_HBC6A['%s-FaOOFaON-5.0'  % (dbse)] =  -6.597
BIND_HBC6A['%s-FaOOFaON-5.4'  % (dbse)] =  -4.321
BIND_HBC6A['%s-FaOOFaON-5.8'  % (dbse)] =  -2.931
BIND_HBC6A['%s-FaOOFaON-6.4'  % (dbse)] =  -1.751
BIND_HBC6A['%s-FaOOFaON-7.0'  % (dbse)] =  -1.119
BIND_HBC6A['%s-FaOOFaON-8.0'  % (dbse)] =  -0.597
BIND_HBC6A['%s-FaOOFaON-10.0' % (dbse)] =  -0.228

BIND_HBC6A['%s-FaONFaNN-3.4'  % (dbse)] =  -8.021
BIND_HBC6A['%s-FaONFaNN-3.5'  % (dbse)] = -10.711
BIND_HBC6A['%s-FaONFaNN-3.6'  % (dbse)] = -12.862
BIND_HBC6A['%s-FaONFaNN-3.7'  % (dbse)] = -14.539
BIND_HBC6A['%s-FaONFaNN-3.8'  % (dbse)] = -15.763
BIND_HBC6A['%s-FaONFaNN-3.9'  % (dbse)] = -16.532
BIND_HBC6A['%s-FaONFaNN-4.0'  % (dbse)] = -16.856  # FaONFaNN minimum
BIND_HBC6A['%s-FaONFaNN-4.1'  % (dbse)] = -16.760
BIND_HBC6A['%s-FaONFaNN-4.2'  % (dbse)] = -16.301
BIND_HBC6A['%s-FaONFaNN-4.3'  % (dbse)] = -15.557
BIND_HBC6A['%s-FaONFaNN-4.4'  % (dbse)] = -14.614
BIND_HBC6A['%s-FaONFaNN-4.6'  % (dbse)] = -12.448
BIND_HBC6A['%s-FaONFaNN-4.8'  % (dbse)] = -10.277
BIND_HBC6A['%s-FaONFaNN-5.0'  % (dbse)] =  -8.341
BIND_HBC6A['%s-FaONFaNN-5.4'  % (dbse)] =  -5.434
BIND_HBC6A['%s-FaONFaNN-5.8'  % (dbse)] =  -3.609
BIND_HBC6A['%s-FaONFaNN-6.4'  % (dbse)] =  -2.082
BIND_HBC6A['%s-FaONFaNN-7.0'  % (dbse)] =  -1.292
BIND_HBC6A['%s-FaONFaNN-8.0'  % (dbse)] =  -0.661
BIND_HBC6A['%s-FaONFaNN-10.0' % (dbse)] =  -0.237

BIND_HBC6A['%s-FaOOFaNN-3.6'  % (dbse)] = -26.064  # FaOOFaNN minimum
BIND_HBC6A['%s-FaOOFaNN-3.7'  % (dbse)] = -23.841
BIND_HBC6A['%s-FaOOFaNN-3.8'  % (dbse)] = -22.850
BIND_HBC6A['%s-FaOOFaNN-3.9'  % (dbse)] = -21.990
BIND_HBC6A['%s-FaOOFaNN-4.0'  % (dbse)] = -21.002
BIND_HBC6A['%s-FaOOFaNN-4.1'  % (dbse)] = -19.819
BIND_HBC6A['%s-FaOOFaNN-4.2'  % (dbse)] = -18.461
BIND_HBC6A['%s-FaOOFaNN-4.3'  % (dbse)] = -16.988
BIND_HBC6A['%s-FaOOFaNN-4.4'  % (dbse)] = -15.471
BIND_HBC6A['%s-FaOOFaNN-4.6'  % (dbse)] = -12.546
BIND_HBC6A['%s-FaOOFaNN-4.8'  % (dbse)] = -10.006
BIND_HBC6A['%s-FaOOFaNN-5.0'  % (dbse)] =  -7.942
BIND_HBC6A['%s-FaOOFaNN-5.4'  % (dbse)] =  -5.058
BIND_HBC6A['%s-FaOOFaNN-5.8'  % (dbse)] =  -3.328
BIND_HBC6A['%s-FaOOFaNN-6.4'  % (dbse)] =  -1.900
BIND_HBC6A['%s-FaOOFaNN-7.0'  % (dbse)] =  -1.166
BIND_HBC6A['%s-FaOOFaNN-8.0'  % (dbse)] =  -0.584
BIND_HBC6A['%s-FaOOFaNN-10.0' % (dbse)] =  -0.200
# Current revision level with deformation correction
BIND_HBC6ARLX = {}
BIND_HBC6ARLX['%s-FaOOFaOO-3.4'  % (dbse)] =  -7.072
BIND_HBC6ARLX['%s-FaOOFaOO-3.5'  % (dbse)] = -11.415
BIND_HBC6ARLX['%s-FaOOFaOO-3.6'  % (dbse)] = -14.186
BIND_HBC6ARLX['%s-FaOOFaOO-3.7'  % (dbse)] = -15.667
BIND_HBC6ARLX['%s-FaOOFaOO-3.8'  % (dbse)] = -16.146  # FaOOFaOO minimum
BIND_HBC6ARLX['%s-FaOOFaOO-3.9'  % (dbse)] = -15.900
BIND_HBC6ARLX['%s-FaOOFaOO-4.0'  % (dbse)] = -15.171
BIND_HBC6ARLX['%s-FaOOFaOO-4.1'  % (dbse)] = -14.153
BIND_HBC6ARLX['%s-FaOOFaOO-4.2'  % (dbse)] = -12.993
BIND_HBC6ARLX['%s-FaOOFaOO-4.3'  % (dbse)] = -11.792
BIND_HBC6ARLX['%s-FaOOFaOO-4.4'  % (dbse)] = -10.617
BIND_HBC6ARLX['%s-FaOOFaOO-4.6'  % (dbse)] =  -8.486
BIND_HBC6ARLX['%s-FaOOFaOO-4.8'  % (dbse)] =  -6.731
BIND_HBC6ARLX['%s-FaOOFaOO-5.0'  % (dbse)] =  -5.341
BIND_HBC6ARLX['%s-FaOOFaOO-5.4'  % (dbse)] =  -3.416
BIND_HBC6ARLX['%s-FaOOFaOO-5.8'  % (dbse)] =  -2.251
BIND_HBC6ARLX['%s-FaOOFaOO-6.4'  % (dbse)] =  -1.284
BIND_HBC6ARLX['%s-FaOOFaOO-7.0'  % (dbse)] =  -0.784
BIND_HBC6ARLX['%s-FaOOFaOO-8.0'  % (dbse)] =  -0.389
BIND_HBC6ARLX['%s-FaOOFaOO-10.0' % (dbse)] =  -0.133

BIND_HBC6ARLX['%s-FaONFaON-3.4'  % (dbse)] =   4.943
BIND_HBC6ARLX['%s-FaONFaON-3.5'  % (dbse)] =  -1.431
BIND_HBC6ARLX['%s-FaONFaON-3.6'  % (dbse)] =  -6.432
BIND_HBC6ARLX['%s-FaONFaON-3.7'  % (dbse)] = -10.102
BIND_HBC6ARLX['%s-FaONFaON-3.8'  % (dbse)] = -12.566
BIND_HBC6ARLX['%s-FaONFaON-3.9'  % (dbse)] = -14.000
BIND_HBC6ARLX['%s-FaONFaON-4.0'  % (dbse)] = -14.603  # FaONFaON minimum
BIND_HBC6ARLX['%s-FaONFaON-4.1'  % (dbse)] = -14.579
BIND_HBC6ARLX['%s-FaONFaON-4.2'  % (dbse)] = -14.112
BIND_HBC6ARLX['%s-FaONFaON-4.3'  % (dbse)] = -13.361
BIND_HBC6ARLX['%s-FaONFaON-4.4'  % (dbse)] = -12.451
BIND_HBC6ARLX['%s-FaONFaON-4.6'  % (dbse)] = -10.489
BIND_HBC6ARLX['%s-FaONFaON-4.8'  % (dbse)] =  -8.655
BIND_HBC6ARLX['%s-FaONFaON-5.0'  % (dbse)] =  -7.101
BIND_HBC6ARLX['%s-FaONFaON-5.4'  % (dbse)] =  -4.830
BIND_HBC6ARLX['%s-FaONFaON-5.8'  % (dbse)] =  -3.380
BIND_HBC6ARLX['%s-FaONFaON-6.4'  % (dbse)] =  -2.110
BIND_HBC6ARLX['%s-FaONFaON-7.0'  % (dbse)] =  -1.403
BIND_HBC6ARLX['%s-FaONFaON-8.0'  % (dbse)] =  -0.794
BIND_HBC6ARLX['%s-FaONFaON-10.0' % (dbse)] =  -0.331

BIND_HBC6ARLX['%s-FaNNFaNN-3.4'  % (dbse)] =  14.652
BIND_HBC6ARLX['%s-FaNNFaNN-3.5'  % (dbse)] =   6.948
BIND_HBC6ARLX['%s-FaNNFaNN-3.6'  % (dbse)] =   0.563
BIND_HBC6ARLX['%s-FaNNFaNN-3.7'  % (dbse)] =  -4.544
BIND_HBC6ARLX['%s-FaNNFaNN-3.8'  % (dbse)] =  -8.441
BIND_HBC6ARLX['%s-FaNNFaNN-3.9'  % (dbse)] = -11.223
BIND_HBC6ARLX['%s-FaNNFaNN-4.0'  % (dbse)] = -13.021
BIND_HBC6ARLX['%s-FaNNFaNN-4.1'  % (dbse)] = -13.996
BIND_HBC6ARLX['%s-FaNNFaNN-4.2'  % (dbse)] = -14.285  # FaNNFaNN minimum
BIND_HBC6ARLX['%s-FaNNFaNN-4.3'  % (dbse)] = -14.074
BIND_HBC6ARLX['%s-FaNNFaNN-4.4'  % (dbse)] = -13.501
BIND_HBC6ARLX['%s-FaNNFaNN-4.6'  % (dbse)] = -11.755
BIND_HBC6ARLX['%s-FaNNFaNN-4.8'  % (dbse)] =  -9.767
BIND_HBC6ARLX['%s-FaNNFaNN-5.0'  % (dbse)] =  -7.915
BIND_HBC6ARLX['%s-FaNNFaNN-5.4'  % (dbse)] =  -5.073
BIND_HBC6ARLX['%s-FaNNFaNN-5.8'  % (dbse)] =  -3.259
BIND_HBC6ARLX['%s-FaNNFaNN-6.4'  % (dbse)] =  -1.742
BIND_HBC6ARLX['%s-FaNNFaNN-7.0'  % (dbse)] =  -0.990
BIND_HBC6ARLX['%s-FaNNFaNN-8.0'  % (dbse)] =  -0.441
BIND_HBC6ARLX['%s-FaNNFaNN-10.0' % (dbse)] =  -0.121

BIND_HBC6ARLX['%s-FaOOFaON-3.4'  % (dbse)] =  -2.134
BIND_HBC6ARLX['%s-FaOOFaON-3.5'  % (dbse)] =  -7.505
BIND_HBC6ARLX['%s-FaOOFaON-3.6'  % (dbse)] = -11.323
BIND_HBC6ARLX['%s-FaOOFaON-3.7'  % (dbse)] = -13.775
BIND_HBC6ARLX['%s-FaOOFaON-3.8'  % (dbse)] = -15.093
BIND_HBC6ARLX['%s-FaOOFaON-3.9'  % (dbse)] = -15.524  # FaOOFaON minimum
BIND_HBC6ARLX['%s-FaOOFaON-4.0'  % (dbse)] = -15.308
BIND_HBC6ARLX['%s-FaOOFaON-4.1'  % (dbse)] = -14.658
BIND_HBC6ARLX['%s-FaOOFaON-4.2'  % (dbse)] = -13.747
BIND_HBC6ARLX['%s-FaOOFaON-4.3'  % (dbse)] = -12.705
BIND_HBC6ARLX['%s-FaOOFaON-4.4'  % (dbse)] = -11.620
BIND_HBC6ARLX['%s-FaOOFaON-4.6'  % (dbse)] =  -9.536
BIND_HBC6ARLX['%s-FaOOFaON-4.8'  % (dbse)] =  -7.730
BIND_HBC6ARLX['%s-FaOOFaON-5.0'  % (dbse)] =  -6.252
BIND_HBC6ARLX['%s-FaOOFaON-5.4'  % (dbse)] =  -4.148
BIND_HBC6ARLX['%s-FaOOFaON-5.8'  % (dbse)] =  -2.834
BIND_HBC6ARLX['%s-FaOOFaON-6.4'  % (dbse)] =  -1.704
BIND_HBC6ARLX['%s-FaOOFaON-7.0'  % (dbse)] =  -1.094
BIND_HBC6ARLX['%s-FaOOFaON-8.0'  % (dbse)] =  -0.587
BIND_HBC6ARLX['%s-FaOOFaON-10.0' % (dbse)] =  -0.226

BIND_HBC6ARLX['%s-FaONFaNN-3.4'  % (dbse)] =   9.365
BIND_HBC6ARLX['%s-FaONFaNN-3.5'  % (dbse)] =   2.303
BIND_HBC6ARLX['%s-FaONFaNN-3.6'  % (dbse)] =  -3.396
BIND_HBC6ARLX['%s-FaONFaNN-3.7'  % (dbse)] =  -7.780
BIND_HBC6ARLX['%s-FaONFaNN-3.8'  % (dbse)] = -10.944
BIND_HBC6ARLX['%s-FaONFaNN-3.9'  % (dbse)] = -13.026
BIND_HBC6ARLX['%s-FaONFaNN-4.0'  % (dbse)] = -14.191
BIND_HBC6ARLX['%s-FaONFaNN-4.1'  % (dbse)] = -14.622  # FaONFaNN minimum
BIND_HBC6ARLX['%s-FaONFaNN-4.2'  % (dbse)] = -14.499
BIND_HBC6ARLX['%s-FaONFaNN-4.3'  % (dbse)] = -13.984
BIND_HBC6ARLX['%s-FaONFaNN-4.4'  % (dbse)] = -13.216
BIND_HBC6ARLX['%s-FaONFaNN-4.6'  % (dbse)] = -11.325
BIND_HBC6ARLX['%s-FaONFaNN-4.8'  % (dbse)] =  -9.389
BIND_HBC6ARLX['%s-FaONFaNN-5.0'  % (dbse)] =  -7.664
BIND_HBC6ARLX['%s-FaONFaNN-5.4'  % (dbse)] =  -5.069
BIND_HBC6ARLX['%s-FaONFaNN-5.8'  % (dbse)] =  -3.412
BIND_HBC6ARLX['%s-FaONFaNN-6.4'  % (dbse)] =  -1.993
BIND_HBC6ARLX['%s-FaONFaNN-7.0'  % (dbse)] =  -1.245
BIND_HBC6ARLX['%s-FaONFaNN-8.0'  % (dbse)] =  -0.642
BIND_HBC6ARLX['%s-FaONFaNN-10.0' % (dbse)] =  -0.232

BIND_HBC6ARLX['%s-FaOOFaNN-3.6'  % (dbse)] = -12.415
BIND_HBC6ARLX['%s-FaOOFaNN-3.7'  % (dbse)] = -15.329
BIND_HBC6ARLX['%s-FaOOFaNN-3.8'  % (dbse)] = -17.085
BIND_HBC6ARLX['%s-FaOOFaNN-3.9'  % (dbse)] = -17.872
BIND_HBC6ARLX['%s-FaOOFaNN-4.0'  % (dbse)] = -17.895  # FaOOFaNN minimum
BIND_HBC6ARLX['%s-FaOOFaNN-4.1'  % (dbse)] = -17.356
BIND_HBC6ARLX['%s-FaOOFaNN-4.2'  % (dbse)] = -16.438
BIND_HBC6ARLX['%s-FaOOFaNN-4.3'  % (dbse)] = -15.294
BIND_HBC6ARLX['%s-FaOOFaNN-4.4'  % (dbse)] = -14.044
BIND_HBC6ARLX['%s-FaOOFaNN-4.6'  % (dbse)] = -11.535
BIND_HBC6ARLX['%s-FaOOFaNN-4.8'  % (dbse)] =  -9.301
BIND_HBC6ARLX['%s-FaOOFaNN-5.0'  % (dbse)] =  -7.458
BIND_HBC6ARLX['%s-FaOOFaNN-5.4'  % (dbse)] =  -4.830
BIND_HBC6ARLX['%s-FaOOFaNN-5.8'  % (dbse)] =  -3.212
BIND_HBC6ARLX['%s-FaOOFaNN-6.4'  % (dbse)] =  -1.850
BIND_HBC6ARLX['%s-FaOOFaNN-7.0'  % (dbse)] =  -1.140
BIND_HBC6ARLX['%s-FaOOFaNN-8.0'  % (dbse)] =  -0.575
BIND_HBC6ARLX['%s-FaOOFaNN-10.0' % (dbse)] =  -0.197
# Set default
BIND = BIND_HBC6A

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')
for item in FaOOFaOO:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formic Acid Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formic Acid Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formic Acid from Formic Acid Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formic Acid from Formic Acid Dimer at %s A' % (molname.group(2))

for item in FaONFaON:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formamide Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formamide Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formamide from Formamide Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formamide from Formamide Dimer at %s A' % (molname.group(2))

for item in FaNNFaNN:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formamidine Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formamidine Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formamidine from Formamidine Dimer at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formamidine from Formamidine Dimer at %s A' % (molname.group(2))

for item in FaOOFaON:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formic Acid-Formamide Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formic Acid-Formamide Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formic Acid from Formic Acid-Formamide Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Formamide from Formic Acid-Formamide Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formic Acid from Formic Acid-Formamide Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] = 'Formamide from Formic Acid-Formamide Complex at %s A' % (molname.group(2))

for item in FaONFaNN:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formamide-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formamide-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formamide from Formamide-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Formamidine from Formamide-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formamide from Formamide-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] = 'Formamidine from Formamide-Formamidine Complex at %s A' % (molname.group(2))

for item in FaOOFaNN:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'Formic Acid-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] = 'Formic Acid-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] = 'Formic Acid from Formic Acid-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] = 'Formamidine from Formic Acid-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] = 'Formic Acid from Formic Acid-Formamidine Complex at %s A' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] = 'Formamidine from Formic Acid-Formamidine Complex at %s A' % (molname.group(2))

TAGL['%s-FaOO-mono-RLX'  % (dbse)]  = 'Formic Acid Relaxed Monomer'
TAGL['%s-FaON-mono-RLX'  % (dbse)]  = 'Formamide Relaxed Monomer'
TAGL['%s-FaNN-mono-RLX'  % (dbse)]  = 'Formamidine Relaxed Monomer'

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

HBC1_FaOOFaOO_3p4 = input.process_input("""
molecule dimer {
0 1
C        1.69147262      -0.17006280       0.00000000
H        2.79500199      -0.28101305       0.00000000
O        1.02814129      -1.21720864       0.00000000
O        1.36966587       1.08860681       0.00000000
H        0.34380745       1.18798183       0.00000000
--
0 1
C       -1.69147262       0.17006280       0.00000000
H       -2.79500199       0.28101305       0.00000000
O       -1.02814129       1.21720864       0.00000000
O       -1.36966587      -1.08860681       0.00000000
H       -0.34380745      -1.18798183       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_3p5 = input.process_input("""
molecule dimer {
0 1
C        1.74073379      -0.17985247       0.00000000
H        2.84248921      -0.29368574       0.00000000
O        1.04839226      -1.20544675       0.00000000
O        1.40587723       1.08303481       0.00000000
H        0.38948927       1.16733829       0.00000000
--
0 1
C       -1.74073379       0.17985247       0.00000000
H       -2.84248921       0.29368574       0.00000000
O       -1.04839226       1.20544675       0.00000000
O       -1.40587723      -1.08303481       0.00000000
H       -0.38948927      -1.16733829       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_3p6 = input.process_input("""
molecule dimer {
0 1
C        1.79035823      -0.18606050       0.00000000
H        2.89087214      -0.30042988       0.00000000
O        1.07568931      -1.19425943       0.00000000
O        1.44185816       1.08049605       0.00000000
H        0.43274661       1.15045330       0.00000000
--
0 1
C       -1.79035823       0.18606050       0.00000000
H       -2.89087214       0.30042988       0.00000000
O       -1.07568931       1.19425943       0.00000000
O       -1.44185816      -1.08049605       0.00000000
H       -0.43274661      -1.15045330       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_3p7 = input.process_input("""
molecule dimer {
0 1
C        1.84016492      -0.19051039       0.00000000
H        2.93982222      -0.30435679       0.00000000
O        1.10803623      -1.18439540       0.00000000
O        1.47971186       1.07967254       0.00000000
H        0.47644336       1.13716323       0.00000000
--
0 1
C       -1.84016492       0.19051039       0.00000000
H       -2.93982222       0.30435679       0.00000000
O       -1.10803623       1.18439540       0.00000000
O       -1.47971186      -1.07967254       0.00000000
H       -0.47644336      -1.13716323       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_3p8 = input.process_input("""
molecule dimer {
0 1
C        1.89005169      -0.19417983       0.00000000
H        2.98915191      -0.30709901       0.00000000
O        1.14427894      -1.17618187       0.00000000
O        1.52059687       1.07982181       0.00000000
H        0.52216282       1.12736965       0.00000000
--
0 1
C       -1.89005169       0.19417983       0.00000000
H       -2.98915191       0.30709901       0.00000000
O       -1.14427894       1.17618187       0.00000000
O       -1.52059687      -1.07982181       0.00000000
H       -0.52216282      -1.12736965       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_3p9 = input.process_input("""
molecule dimer {
0 1
C        1.93997588      -0.19747119       0.00000000
H        3.03876443      -0.30931746       0.00000000
O        1.18372180      -1.16964417       0.00000000
O        1.56485002       1.08053745       0.00000000
H        0.57046128       1.12083453       0.00000000
--
0 1
C       -1.93997588       0.19747119       0.00000000
H       -3.03876443       0.30931746       0.00000000
O       -1.18372180       1.16964417       0.00000000
O       -1.56485002      -1.08053745       0.00000000
H       -0.57046128      -1.12083453       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p0 = input.process_input("""
molecule dimer {
0 1
C        1.98993209      -0.20042861       0.00000000
H        3.08861310      -0.31108922       0.00000000
O        1.22585974      -1.16459704       0.00000000
O        1.61210992       1.08161540       0.00000000
H        0.62108112       1.11708149       0.00000000
--
0 1
C       -1.98993209       0.20042861       0.00000000
H       -3.08861310       0.31108922       0.00000000
O       -1.22585974       1.16459704       0.00000000
O       -1.61210992      -1.08161540       0.00000000
H       -0.62108112      -1.11708149       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p1 = input.process_input("""
molecule dimer {
0 1
C        2.03992872      -0.20295867       0.00000000
H        3.13866374      -0.31227513       0.00000000
O        1.27022640      -1.16073661       0.00000000
O        1.66162256       1.08294418       0.00000000
H        0.67334166       1.11546236       0.00000000
--
0 1
C       -2.03992872       0.20295867       0.00000000
H       -3.13866374       0.31227513       0.00000000
O       -1.27022640       1.16073661       0.00000000
O       -1.66162256      -1.08294418       0.00000000
H       -0.67334166      -1.11546236       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p2 = input.process_input("""
molecule dimer {
0 1
C        2.08997101      -0.20499428       0.00000000
H        3.18887965      -0.31278045       0.00000000
O        1.31635032      -1.15775239       0.00000000
O        1.71259791       1.08444718       0.00000000
H        0.72653811       1.11533705       0.00000000
--
0 1
C       -2.08997101       0.20499428       0.00000000
H       -3.18887965       0.31278045       0.00000000
O       -1.31635032       1.15775239       0.00000000
O       -1.71259791      -1.08444718       0.00000000
H       -0.72653811      -1.11533705       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p3 = input.process_input("""
molecule dimer {
0 1
C        2.14005789      -0.20652795       0.00000000
H        3.23921879      -0.31260333       0.00000000
O        1.36379229      -1.15537668       0.00000000
O        1.76438167       1.08605925       0.00000000
H        0.78011615       1.11617462       0.00000000
--
0 1
C       -2.14005789       0.20652795       0.00000000
H       -3.23921879       0.31260333       0.00000000
O       -1.36379229       1.15537668       0.00000000
O       -1.76438167      -1.08605925       0.00000000
H       -0.78011615      -1.11617462       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p4 = input.process_input("""
molecule dimer {
0 1
C        2.19018321      -0.20760327       0.00000000
H        3.28964362      -0.31181902       0.00000000
O        1.41218943      -1.15341758       0.00000000
O        1.81652565       1.08773186       0.00000000
H        0.83371879       1.11760094       0.00000000
--
0 1
C       -2.19018321       0.20760327       0.00000000
H       -3.28964362       0.31181902       0.00000000
O       -1.41218943       1.15341758       0.00000000
O       -1.81652565      -1.08773186       0.00000000
H       -0.83371879      -1.11760094       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p6 = input.process_input("""
molecule dimer {
0 1
C        2.29051747      -0.20864206       0.00000000
H        3.39063528      -0.30885122       0.00000000
O        1.51085984      -1.15029332       0.00000000
O        1.92095251       1.09113093       0.00000000
H        0.94036724       1.12135540       0.00000000
--
0 1
C       -2.29051747       0.20864206       0.00000000
H       -3.39063528       0.30885122       0.00000000
O       -1.51085984       1.15029332       0.00000000
O       -1.92095251      -1.09113093       0.00000000
H       -0.94036724      -1.12135540       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_4p8 = input.process_input("""
molecule dimer {
0 1
C        2.39092365      -0.20853307       0.00000000
H        3.49171140      -0.30454226       0.00000000
O        1.61118653      -1.14785646       0.00000000
O        2.02488292       1.09452849       0.00000000
H        1.04594759       1.12557184       0.00000000
--
0 1
C       -2.39092365       0.20853307       0.00000000
H       -3.49171140       0.30454226       0.00000000
O       -1.61118653       1.14785646       0.00000000
O       -2.02488292      -1.09452849       0.00000000
H       -1.04594759      -1.12557184       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_5p0 = input.process_input("""
molecule dimer {
0 1
C        2.49138346      -0.20738986       0.00000000
H        3.59280864      -0.29907564       0.00000000
O        1.71273443      -1.14581025       0.00000000
O        2.12758466       1.09795049       0.00000000
H        1.14991023       1.12938889       0.00000000
--
0 1
C       -2.49138346       0.20738986       0.00000000
H       -3.59280864       0.29907564       0.00000000
O       -1.71273443       1.14581025       0.00000000
O       -2.12758466      -1.09795049       0.00000000
H       -1.14991023      -1.12938889       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_5p4 = input.process_input("""
molecule dimer {
0 1
C        2.69240165      -0.20242478       0.00000000
H        3.79490134      -0.28531481       0.00000000
O        1.91833290      -1.14219531       0.00000000
O        2.32770522       1.10482267       0.00000000
H        1.35155974       1.13337052       0.00000000
--
0 1
C       -2.69240165       0.20242478       0.00000000
H       -3.79490134       0.28531481       0.00000000
O       -1.91833290       1.14219531       0.00000000
O       -2.32770522      -1.10482267       0.00000000
H       -1.35155974      -1.13337052       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_5p8 = input.process_input("""
molecule dimer {
0 1
C        2.89338469      -0.19577496       0.00000000
H        3.99667367      -0.27042675       0.00000000
O        2.12532927      -1.13908975       0.00000000
O        2.52303889       1.11116394       0.00000000
H        1.54758478       1.13325561       0.00000000
--
0 1
C       -2.89338469       0.19577496       0.00000000
H       -3.99667367       0.27042675       0.00000000
O       -2.12532927       1.13908975       0.00000000
O       -2.52303889      -1.11116394       0.00000000
H       -1.54758478      -1.13325561       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_6p4 = input.process_input("""
molecule dimer {
0 1
C        3.19458871      -0.18602798       0.00000000
H        4.29867458      -0.25032135       0.00000000
O        2.43545364      -1.13545523       0.00000000
O        2.81355286       1.11890079       0.00000000
H        1.83855220       1.13009890       0.00000000
--
0 1
C       -3.19458871       0.18602798       0.00000000
H       -4.29867458       0.25032135       0.00000000
O       -2.43545364       1.13545523       0.00000000
O       -2.81355286      -1.11890079       0.00000000
H       -1.83855220      -1.13009890       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_7p0 = input.process_input("""
molecule dimer {
0 1
C        3.49547568      -0.17791547       0.00000000
H        4.60006941      -0.23413795       0.00000000
O        2.74388512      -1.13278853       0.00000000
O        3.10454835       1.12465267       0.00000000
H        2.12977769       1.12626632       0.00000000
--
0 1
C       -3.49547568       0.17791547       0.00000000
H       -4.60006941       0.23413795       0.00000000
O       -2.74388512       1.13278853       0.00000000
O       -3.10454835      -1.12465267       0.00000000
H       -2.12977769      -1.12626632       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_8p0 = input.process_input("""
molecule dimer {
0 1
C        3.99646944      -0.16803989       0.00000000
H        5.10156877      -0.21450609       0.00000000
O        3.25414999      -1.12972545       0.00000000
O        3.59284622       1.13117851       0.00000000
H        2.61830479       1.12095879       0.00000000
--
0 1
C       -3.99646944       0.16803989       0.00000000
H       -5.10156877       0.21450609       0.00000000
O       -3.25414999       1.12972545       0.00000000
O       -3.59284622      -1.13117851       0.00000000
H       -2.61830479      -1.12095879       0.00000000
units angstrom
}
""")

HBC1_FaOOFaOO_10p0 = input.process_input("""
molecule dimer {
0 1
C        4.99755344      -0.15642268       0.00000000
H        6.10311728      -0.19102667       0.00000000
O        4.26634092      -1.12629100       0.00000000
O        4.57854479       1.13834246       0.00000000
H        3.60431482       1.11461219       0.00000000
--
0 1
C       -4.99755344       0.15642268       0.00000000
H       -6.10311728       0.19102667       0.00000000
O       -4.26634092       1.12629100       0.00000000
O       -4.57854479      -1.13834246       0.00000000
H       -3.60431482      -1.11461219       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p4 = input.process_input("""
molecule dimer {
0 1
C        1.68040472      -0.25737318       0.00000000
H        2.78876519      -0.42713125       0.00000000
O        0.98387212      -1.27944113       0.00000000
N        1.47197100       1.06623396       0.00000000
H        0.51175066       1.43614881       0.00000000
H        2.30581639       1.63815514       0.00000000
--
0 1
C       -1.68040472       0.25737318       0.00000000
H       -2.78876519       0.42713125       0.00000000
O       -0.98387212       1.27944113       0.00000000
N       -1.47197100      -1.06623396       0.00000000
H       -0.51175066      -1.43614881       0.00000000
H       -2.30581639      -1.63815514       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p5 = input.process_input("""
molecule dimer {
0 1
C        1.73857937      -0.19960660       0.00000000
H        2.84972331      -0.32717722       0.00000000
O        1.07441392      -1.24725582       0.00000000
N        1.42842301       1.10192478       0.00000000
H        0.43552514       1.39016691       0.00000000
H        2.20654586       1.74794356       0.00000000
--
0 1
C       -1.73857937       0.19960660       0.00000000
H       -2.84972331       0.32717722       0.00000000
O       -1.07441392       1.24725582       0.00000000
N       -1.42842301      -1.10192478       0.00000000
H       -0.43552514      -1.39016691       0.00000000
H       -2.20654586      -1.74794356       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p6 = input.process_input("""
molecule dimer {
0 1
C        1.79214890      -0.16793868       0.00000000
H        2.90351314      -0.27208239       0.00000000
O        1.13654074      -1.22317842       0.00000000
N        1.41655781       1.11709296       0.00000000
H        0.40800979       1.35225476       0.00000000
H        2.15197143       1.81156642       0.00000000
--
0 1
C       -1.79214890       0.16793868       0.00000000
H       -2.90351314       0.27208239       0.00000000
O       -1.13654074       1.22317842       0.00000000
N       -1.41655781      -1.11709296       0.00000000
H       -0.40800979      -1.35225476       0.00000000
H       -2.15197143      -1.81156642       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p7 = input.process_input("""
molecule dimer {
0 1
C        1.84396206      -0.14934877       0.00000000
H        2.95481237      -0.23932033       0.00000000
O        1.18668392      -1.20479449       0.00000000
N        1.42230702       1.12368101       0.00000000
H        0.40573498       1.32114549       0.00000000
H        2.12347975       1.85285206       0.00000000
--
0 1
C       -1.84396206       0.14934877       0.00000000
H       -2.95481237       0.23932033       0.00000000
O       -1.18668392       1.20479449       0.00000000
N       -1.42230702      -1.12368101       0.00000000
H       -0.40573498      -1.32114549       0.00000000
H       -2.12347975      -1.85285206       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p8 = input.process_input("""
molecule dimer {
0 1
C        1.89495225      -0.13840942       0.00000000
H        3.00512682      -0.21949780       0.00000000
O        1.23153916      -1.19064307       0.00000000
N        1.43998729       1.12629331       0.00000000
H        0.41948596       1.29629750       0.00000000
H        2.11367741       1.88100998       0.00000000
--
0 1
C       -1.89495225       0.13840942       0.00000000
H       -3.00512682       0.21949780       0.00000000
O       -1.23153916       1.19064307       0.00000000
N       -1.43998729      -1.12629331       0.00000000
H       -0.41948596      -1.29629750       0.00000000
H       -2.11367741      -1.88100998       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_3p9 = input.process_input("""
molecule dimer {
0 1
C        1.94550597      -0.13231698       0.00000000
H        3.05506108      -0.20777961       0.00000000
O        1.27444696      -1.17984091       0.00000000
N        1.46671040       1.12708677       0.00000000
H        0.44465999       1.27731960       0.00000000
H        2.11884877       1.90053586       0.00000000
--
0 1
C       -1.94550597       0.13231698       0.00000000
H       -3.05506108       0.20777961       0.00000000
O       -1.27444696       1.17984091       0.00000000
N       -1.46671040      -1.12708677       0.00000000
H       -0.44465999      -1.27731960       0.00000000
H       -2.11884877      -1.90053586       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p0 = input.process_input("""
molecule dimer {
0 1
C        1.99581166      -0.12937222       0.00000000
H        3.10488983      -0.20126473       0.00000000
O        1.31722472      -1.17176432       0.00000000
N        1.50061758       1.12716306       0.00000000
H        0.47842639       1.26369880       0.00000000
H        2.13649227       1.91403520       0.00000000
--
0 1
C       -1.99581166       0.12937222       0.00000000
H       -3.10488983       0.20126473       0.00000000
O       -1.31722472       1.17176432       0.00000000
N       -1.50061758      -1.12716306       0.00000000
H       -0.47842639      -1.26369880       0.00000000
H       -2.13649227      -1.91403520       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p1 = input.process_input("""
molecule dimer {
0 1
C        2.04597250      -0.12844429       0.00000000
H        3.15473976      -0.19805169       0.00000000
O        1.36080054      -1.16589088       0.00000000
N        1.54023054       1.12707489       0.00000000
H        0.51870599       1.25472290       0.00000000
H        2.16442795       1.92321328       0.00000000
--
0 1
C       -2.04597250       0.12844429       0.00000000
H       -3.15473976       0.19805169       0.00000000
O       -1.36080054       1.16589088       0.00000000
N       -1.54023054      -1.12707489       0.00000000
H       -0.51870599      -1.25472290       0.00000000
H       -2.16442795      -1.92321328       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p2 = input.process_input("""
molecule dimer {
0 1
C        2.09604836      -0.12877408       0.00000000
H        3.20466046      -0.19688344       0.00000000
O        1.40550526      -1.16174638       0.00000000
N        1.58426128       1.12705799       0.00000000
H        0.56383358       1.24955108       0.00000000
H        2.20059240       1.92925606       0.00000000
--
0 1
C       -2.09604836       0.12877408       0.00000000
H       -3.20466046       0.19688344       0.00000000
O       -1.40550526       1.16174638       0.00000000
N       -1.58426128      -1.12705799       0.00000000
H       -0.56383358      -1.24955108       0.00000000
H       -2.20059240      -1.92925606       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p3 = input.process_input("""
molecule dimer {
0 1
C        2.14607510      -0.12985845       0.00000000
H        3.25466017      -0.19693865       0.00000000
O        1.45130851      -1.15890002       0.00000000
N        1.63157574       1.12717703       0.00000000
H        0.61244489       1.24732130       0.00000000
H        2.24305687       1.93302516       0.00000000
--
0 1
C       -2.14607510       0.12985845       0.00000000
H       -3.25466017       0.19693865       0.00000000
O       -1.45130851       1.15890002       0.00000000
N       -1.63157574      -1.12717703       0.00000000
H       -0.61244489      -1.24732130       0.00000000
H       -2.24305687      -1.93302516       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p4 = input.process_input("""
molecule dimer {
0 1
C        2.19607517      -0.13136009       0.00000000
H        3.30472871      -0.19767514       0.00000000
O        1.49802227      -1.15698038       0.00000000
N        1.68120383       1.12742404       0.00000000
H        0.66343356       1.24723827       0.00000000
H        2.29010119       1.93517379       0.00000000
--
0 1
C       -2.19607517       0.13136009       0.00000000
H       -3.30472871       0.19767514       0.00000000
O       -1.49802227       1.15698038       0.00000000
N       -1.68120383      -1.12742404       0.00000000
H       -0.66343356      -1.24723827       0.00000000
H       -2.29010119      -1.93517379       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p6 = input.process_input("""
molecule dimer {
0 1
C        2.29604573      -0.13481753       0.00000000
H        3.40501312      -0.19993307       0.00000000
O        1.59337012      -1.15484644       0.00000000
N        1.78457816       1.12819246       0.00000000
H        0.76947608       1.25110815       0.00000000
H        2.39273653       1.93641581       0.00000000
--
0 1
C       -2.29604573       0.13481753       0.00000000
H       -3.40501312       0.19993307       0.00000000
O       -1.59337012       1.15484644       0.00000000
N       -1.78457816      -1.12819246       0.00000000
H       -0.76947608      -1.25110815       0.00000000
H       -2.39273653      -1.93641581       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_4p8 = input.process_input("""
molecule dimer {
0 1
C        2.39601143      -0.13831545       0.00000000
H        3.50541998      -0.20235870       0.00000000
O        1.69032936      -1.15398956       0.00000000
N        1.89087102       1.12918383       0.00000000
H        0.87833293       1.25816233       0.00000000
H        2.50182267       1.93524951       0.00000000
--
0 1
C       -2.39601143       0.13831545       0.00000000
H       -3.50541998       0.20235870       0.00000000
O       -1.69032936       1.15398956       0.00000000
N       -1.89087102      -1.12918383       0.00000000
H       -0.87833293      -1.25816233       0.00000000
H       -2.50182267      -1.93524951       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_5p0 = input.process_input("""
molecule dimer {
0 1
C        2.49599331      -0.14149026       0.00000000
H        3.60589309      -0.20440710       0.00000000
O        1.78830561      -1.15373352       0.00000000
N        1.99827465       1.13031983       0.00000000
H        0.98820948       1.26676542       0.00000000
H        2.61374032       1.93291817       0.00000000
--
0 1
C       -2.49599331       0.14149026       0.00000000
H       -3.60589309       0.20440710       0.00000000
O       -1.78830561       1.15373352       0.00000000
N       -1.99827465      -1.13031983       0.00000000
H       -0.98820948      -1.26676542       0.00000000
H       -2.61374032      -1.93291817       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_5p4 = input.process_input("""
molecule dimer {
0 1
C        2.69608420      -0.14536998       0.00000000
H        3.80688469      -0.20526316       0.00000000
O        1.98650872      -1.15320220       0.00000000
N        2.21035574       1.13294283       0.00000000
H        1.20419348       1.28204746       0.00000000
H        2.83403804       1.92915087       0.00000000
--
0 1
C       -2.69608420       0.14536998       0.00000000
H       -3.80688469       0.20526316       0.00000000
O       -1.98650872       1.15320220       0.00000000
N       -2.21035574      -1.13294283       0.00000000
H       -1.20419348      -1.28204746       0.00000000
H       -2.83403804      -1.92915087       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_5p8 = input.process_input("""
molecule dimer {
0 1
C        2.89633015      -0.14585783       0.00000000
H        4.00779084      -0.20183046       0.00000000
O        2.18702642      -1.15195636       0.00000000
N        2.41611039       1.13579344       0.00000000
H        1.41208234       1.29109009       0.00000000
H        3.04382363       1.92881894       0.00000000
--
0 1
C       -2.89633015       0.14585783       0.00000000
H       -4.00779084       0.20183046       0.00000000
O       -2.18702642       1.15195636       0.00000000
N       -2.41611039      -1.13579344       0.00000000
H       -1.41208234      -1.29109009       0.00000000
H       -3.04382363      -1.92881894       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_6p4 = input.process_input("""
molecule dimer {
0 1
C        3.19679060      -0.14329444       0.00000000
H        4.30890404      -0.19314433       0.00000000
O        2.49028966      -1.14973005       0.00000000
N        2.71746327       1.13986108       0.00000000
H        1.71465256       1.29686270       0.00000000
H        3.34615668       1.93211051       0.00000000
--
0 1
C       -3.19679060       0.14329444       0.00000000
H       -4.30890404       0.19314433       0.00000000
O       -2.49028966       1.14973005       0.00000000
N       -2.71746327      -1.13986108       0.00000000
H       -1.71465256      -1.29686270       0.00000000
H       -3.34615668      -1.93211051       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_7p0 = input.process_input("""
molecule dimer {
0 1
C        3.49721450      -0.13962406       0.00000000
H        4.60973564      -0.18404076       0.00000000
O        2.79431063      -1.14771035       0.00000000
N        3.01520025       1.14321645       0.00000000
H        2.01252678       1.29808475       0.00000000
H        3.64229745       1.93673374       0.00000000
--
0 1
C       -3.49721450       0.13962406       0.00000000
H       -4.60973564       0.18404076       0.00000000
O       -2.79431063       1.14771035       0.00000000
N       -3.01520025      -1.14321645       0.00000000
H       -2.01252678      -1.29808475       0.00000000
H       -3.64229745      -1.93673374       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_8p0 = input.process_input("""
molecule dimer {
0 1
C        3.99775075      -0.13414290       0.00000000
H        5.11067288      -0.17148655       0.00000000
O        3.30011621      -1.14514769       0.00000000
N        3.51015001       1.14724727       0.00000000
H        2.50709338       1.29726882       0.00000000
H        4.13374763       1.94351670       0.00000000
--
0 1
C       -3.99775075       0.13414290       0.00000000
H       -5.11067288       0.17148655       0.00000000
O       -3.30011621       1.14514769       0.00000000
N       -3.51015001      -1.14724727       0.00000000
H       -2.50709338      -1.29726882       0.00000000
H       -4.13374763      -1.94351670       0.00000000
units angstrom
}
""")

HBC1_FaONFaON_10p0 = input.process_input("""
molecule dimer {
0 1
C        4.99839541      -0.12669527       0.00000000
H        6.11168302      -0.15491398       0.00000000
O        4.30803664      -1.14211866       0.00000000
N        4.50192844       1.15187623       0.00000000
H        3.49800610       1.29436065       0.00000000
H        5.11998480       1.95244688       0.00000000
--
0 1
C       -4.99839541       0.12669527       0.00000000
H       -6.11168302       0.15491398       0.00000000
O       -4.30803664       1.14211866       0.00000000
N       -4.50192844      -1.15187623       0.00000000
H       -3.49800610      -1.29436065       0.00000000
H       -5.11998480      -1.95244688       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p4 = input.process_input("""
molecule dimer {
0 1
C        1.69498253      -0.13051889       0.00000000
H        2.80633596      -0.21609653       0.00000000
N        1.09905882      -1.29491251       0.00000000
H        1.82089470      -2.01478416       0.00000000
N        1.37286317       1.16800831       0.00000000
H        0.35212595       1.46532985       0.00000000
H        2.16581591       1.79503126       0.00000000
--
0 1
C       -1.69498253       0.13051889       0.00000000
H       -2.80633596       0.21609653       0.00000000
N       -1.09905882       1.29491251       0.00000000
H       -1.82089470       2.01478416       0.00000000
N       -1.37286317      -1.16800831       0.00000000
H       -0.35212595      -1.46532985       0.00000000
H       -2.16581591      -1.79503126       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p5 = input.process_input("""
molecule dimer {
0 1
C        1.74573454      -0.12211451       0.00000000
H        2.85588849      -0.19977002       0.00000000
N        1.13008345      -1.27542245       0.00000000
H        1.83665594      -2.01109010       0.00000000
N        1.38209118       1.16632433       0.00000000
H        0.35550066       1.43042391       0.00000000
H        2.14642150       1.82776097       0.00000000
--
0 1
C       -1.74573454       0.12211451       0.00000000
H       -2.85588849       0.19977002       0.00000000
N       -1.13008345       1.27542245       0.00000000
H       -1.83665594       2.01109010       0.00000000
N       -1.38209118      -1.16632433       0.00000000
H       -0.35550066      -1.43042391       0.00000000
H       -2.14642150      -1.82776097       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p6 = input.process_input("""
molecule dimer {
0 1
C        1.79626801      -0.11585476       0.00000000
H        2.90515430      -0.18737513       0.00000000
N        1.16061610      -1.25790391       0.00000000
H        1.85182448      -2.00884269       0.00000000
N        1.39570164       1.16335736       0.00000000
H        0.36558676       1.39719019       0.00000000
H        2.13213780       1.85550007       0.00000000
--
0 1
C       -1.79626801       0.11585476       0.00000000
H       -2.90515430       0.18737513       0.00000000
N       -1.16061610       1.25790391       0.00000000
H       -1.85182448       2.00884269       0.00000000
N       -1.39570164      -1.16335736       0.00000000
H       -0.36558676      -1.39719019       0.00000000
H       -2.13213780      -1.85550007       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p7 = input.process_input("""
molecule dimer {
0 1
C        1.84665575      -0.11119213       0.00000000
H        2.95429208      -0.17788590       0.00000000
N        1.19156721      -1.24220748       0.00000000
H        1.86789891      -2.00737666       0.00000000
N        1.41335406       1.15979090       0.00000000
H        0.38141976       1.36618041       0.00000000
H        2.12308262       1.87905403       0.00000000
--
0 1
C       -1.84665575       0.11119213       0.00000000
H       -2.95429208       0.17788590       0.00000000
N       -1.19156721       1.24220748       0.00000000
H       -1.86789891       2.00737666       0.00000000
N       -1.41335406      -1.15979090       0.00000000
H       -0.38141976      -1.36618041       0.00000000
H       -2.12308262      -1.87905403       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p8 = input.process_input("""
molecule dimer {
0 1
C        1.89693820      -0.10782700       0.00000000
H        3.00339990      -0.17072121       0.00000000
N        1.22354936      -1.22831838       0.00000000
H        1.88594286      -2.00636060       0.00000000
N        1.43499889       1.15607525       0.00000000
H        0.40250444       1.33793562       0.00000000
H        2.11977084       1.89888692       0.00000000
--
0 1
C       -1.89693820       0.10782700       0.00000000
H       -3.00339990       0.17072121       0.00000000
N       -1.22354936       1.22831838       0.00000000
H       -1.88594286       2.00636060       0.00000000
N       -1.43499889      -1.15607525       0.00000000
H       -0.40250444      -1.33793562       0.00000000
H       -2.11977084      -1.89888692       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_3p9 = input.process_input("""
molecule dimer {
0 1
C        1.94713954      -0.10558828       0.00000000
H        3.05254284      -0.16553141       0.00000000
N        1.25700152      -1.21627302       0.00000000
H        1.90673736      -2.00565055       0.00000000
N        1.46072682       1.15251669       0.00000000
H        0.42859974       1.31294465       0.00000000
H        2.12286880       1.91531762       0.00000000
--
0 1
C       -1.94713954       0.10558828       0.00000000
H       -3.05254284       0.16553141       0.00000000
N       -1.25700152       1.21627302       0.00000000
H       -1.90673736       2.00565055       0.00000000
N       -1.46072682      -1.15251669       0.00000000
H       -0.42859974      -1.31294465       0.00000000
H       -2.12286880      -1.91531762       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p0 = input.process_input("""
molecule dimer {
0 1
C        1.99727808      -0.10431489       0.00000000
H        3.10176903      -0.16200082       0.00000000
N        1.29228783      -1.20609151       0.00000000
H        1.93090747      -2.00517101       0.00000000
N        1.49059653       1.14934169       0.00000000
H        0.45949690       1.29156428       0.00000000
H        2.13289668       1.92867591       0.00000000
--
0 1
C       -1.99727808       0.10431489       0.00000000
H       -3.10176903       0.16200082       0.00000000
N       -1.29228783       1.20609151       0.00000000
H       -1.93090747       2.00517101       0.00000000
N       -1.49059653      -1.14934169       0.00000000
H       -0.45949690      -1.29156428       0.00000000
H       -2.13289668      -1.92867591       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p1 = input.process_input("""
molecule dimer {
0 1
C        2.04736978      -0.10381896       0.00000000
H        3.15111193      -0.15978802       0.00000000
N        1.32968607      -1.19773394       0.00000000
H        1.95888022      -2.00487461       0.00000000
N        1.52453108       1.14670134       0.00000000
H        0.49489283       1.27394294       0.00000000
H        2.15004136       1.93933564       0.00000000
--
0 1
C       -2.04736978       0.10381896       0.00000000
H       -3.15111193       0.15978802       0.00000000
N       -1.32968607       1.19773394       0.00000000
H       -1.95888022       2.00487461       0.00000000
N       -1.52453108      -1.14670134       0.00000000
H       -0.49489283      -1.27394294       0.00000000
H       -2.15004136      -1.93933564       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p2 = input.process_input("""
molecule dimer {
0 1
C        2.09742845      -0.10390077       0.00000000
H        3.20059131      -0.15854840       0.00000000
N        1.36932506      -1.19108785       0.00000000
H        1.99079936      -2.00473399       0.00000000
N        1.56227621       1.14466550       0.00000000
H        0.53433608       1.25999374       0.00000000
H        2.17408498       1.94769642       0.00000000
--
0 1
C       -2.09742845       0.10390077       0.00000000
H       -3.20059131       0.15854840       0.00000000
N       -1.36932506       1.19108785       0.00000000
H       -1.99079936       2.00473399       0.00000000
N       -1.56227621      -1.14466550       0.00000000
H       -0.53433608      -1.25999374       0.00000000
H       -2.17408498      -1.94769642       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p3 = input.process_input("""
molecule dimer {
0 1
C        2.14746429      -0.10439674       0.00000000
H        3.25020808      -0.15800548       0.00000000
N        1.41111901      -1.18596705       0.00000000
H        2.02643555      -2.00474792       0.00000000
N        1.60342559       1.14321301       0.00000000
H        0.57725829       1.24942013       0.00000000
H        2.20446212       1.95412762       0.00000000
--
0 1
C       -2.14746429       0.10439674       0.00000000
H       -3.25020808       0.15800548       0.00000000
N       -1.41111901       1.18596705       0.00000000
H       -2.02643555       2.00474792       0.00000000
N       -1.60342559      -1.14321301       0.00000000
H       -0.57725829      -1.24942013       0.00000000
H       -2.20446212      -1.95412762       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p4 = input.process_input("""
molecule dimer {
0 1
C        2.19748298      -0.10521494       0.00000000
H        3.29994728      -0.15800065       0.00000000
N        1.45477374      -1.18214133       0.00000000
H        2.06523281      -2.00494056       0.00000000
N        1.64750275       1.14225070       0.00000000
H        0.62306965       1.24181216       0.00000000
H        2.24041883       1.95893145       0.00000000
--
0 1
C       -2.19748298       0.10521494       0.00000000
H       -3.29994728       0.15800065       0.00000000
N       -1.45477374       1.18214133       0.00000000
H       -2.06523281       2.00494056       0.00000000
N       -1.64750275      -1.14225070       0.00000000
H       -0.62306965      -1.24181216       0.00000000
H       -2.24041883      -1.95893145       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p6 = input.process_input("""
molecule dimer {
0 1
C        2.29747604      -0.10772948       0.00000000
H        3.39968456      -0.15941244       0.00000000
N        1.54604150      -1.17743388       0.00000000
H        2.14951038      -2.00594097       0.00000000
N        1.74263302       1.14129835       0.00000000
H        0.72137666       1.23384217       0.00000000
H        2.32599570       1.96451573       0.00000000
--
0 1
C       -2.29747604       0.10772948       0.00000000
H       -3.39968456       0.15941244       0.00000000
N       -1.54604150       1.17743388       0.00000000
H       -2.14951038       2.00594097       0.00000000
N       -1.74263302      -1.14129835       0.00000000
H       -0.72137666      -1.23384217       0.00000000
H       -2.32599570      -1.96451573       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_4p8 = input.process_input("""
molecule dimer {
0 1
C        2.39742249      -0.11120861       0.00000000
H        3.49964090      -0.16233693       0.00000000
N        1.64041315      -1.17535654       0.00000000
H        2.23917184      -2.00764182       0.00000000
N        1.84436251       1.14107385       0.00000000
H        0.82590336       1.23313494       0.00000000
H        2.42493207       1.96604670       0.00000000
--
0 1
C       -2.39742249       0.11120861       0.00000000
H       -3.49964090       0.16233693       0.00000000
N       -1.64041315       1.17535654       0.00000000
H       -2.23917184       2.00764182       0.00000000
N       -1.84436251      -1.14107385       0.00000000
H       -0.82590336      -1.23313494       0.00000000
H       -2.42493207      -1.96604670       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_5p0 = input.process_input("""
molecule dimer {
0 1
C        2.49735133      -0.11505842       0.00000000
H        3.59972532      -0.16584719       0.00000000
N        1.73648339      -1.17469470       0.00000000
H        2.33188443      -2.00962361       0.00000000
N        1.94973034       1.14125016       0.00000000
H        0.93377159       1.23689923       0.00000000
H        2.53164153       1.96514450       0.00000000
--
0 1
C       -2.49735133       0.11505842       0.00000000
H       -3.59972532       0.16584719       0.00000000
N       -1.73648339       1.17469470       0.00000000
H       -2.33188443       2.00962361       0.00000000
N       -1.94973034      -1.14125016       0.00000000
H       -0.93377159      -1.23689923       0.00000000
H       -2.53164153      -1.96514450       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_5p4 = input.process_input("""
molecule dimer {
0 1
C        2.69724643      -0.12191870       0.00000000
H        3.80009014      -0.17176852       0.00000000
N        1.93139799      -1.17474481       0.00000000
H        2.52207952      -2.01322139       0.00000000
N        2.16339519       1.14223517       0.00000000
H        1.15164366       1.24943138       0.00000000
H        2.75270049       1.96073403       0.00000000
--
0 1
C       -2.69724643       0.12191870       0.00000000
H       -3.80009014       0.17176852       0.00000000
N       -1.93139799       1.17474481       0.00000000
H       -2.52207952       2.01322139       0.00000000
N       -2.16339519      -1.14223517       0.00000000
H       -1.15164366      -1.24943138       0.00000000
H       -2.75270049      -1.96073403       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_5p8 = input.process_input("""
molecule dimer {
0 1
C        2.89724776      -0.12632587       0.00000000
H        4.00054647      -0.17443193       0.00000000
N        2.12821587      -1.17451446       0.00000000
H        2.71529330      -2.01553546       0.00000000
N        2.37424699       1.14355972       0.00000000
H        1.36542121       1.26070273       0.00000000
H        2.97042339       1.95703238       0.00000000
--
0 1
C       -2.89724776       0.12632587       0.00000000
H       -4.00054647       0.17443193       0.00000000
N       -2.12821587       1.17451446       0.00000000
H       -2.71529330       2.01553546       0.00000000
N       -2.37424699      -1.14355972       0.00000000
H       -1.36542121      -1.26070273       0.00000000
H       -2.97042339      -1.95703238       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_6p4 = input.process_input("""
molecule dimer {
0 1
C        3.19740673      -0.12881636       0.00000000
H        4.30120367      -0.17328586       0.00000000
N        2.42581025      -1.17305230       0.00000000
H        3.00905003      -2.01668226       0.00000000
N        2.68267636       1.14568753       0.00000000
H        1.67620388       1.27071277       0.00000000
H        3.28436882       1.95508270       0.00000000
--
0 1
C       -3.19740673       0.12881636       0.00000000
H       -4.30120367       0.17328586       0.00000000
N       -2.42581025       1.17305230       0.00000000
H       -3.00905003       2.01668226       0.00000000
N       -2.68267636      -1.14568753       0.00000000
H       -1.67620388      -1.27071277       0.00000000
H       -3.28436882      -1.95508270       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_7p0 = input.process_input("""
molecule dimer {
0 1
C        3.49764421      -0.12840997       0.00000000
H        4.60176692      -0.16894593       0.00000000
N        2.72569862      -1.17126633       0.00000000
H        3.30724390      -2.01600726       0.00000000
N        2.98518869       1.14779235       0.00000000
H        1.97966414       1.27507136       0.00000000
H        3.58843631       1.95603470       0.00000000
--
0 1
C       -3.49764421       0.12840997       0.00000000
H       -4.60176692       0.16894593       0.00000000
N       -2.72569862       1.17126633       0.00000000
H       -3.30724390       2.01600726       0.00000000
N       -2.98518869      -1.14779235       0.00000000
H       -1.97966414      -1.27507136       0.00000000
H       -3.58843631      -1.95603470       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_8p0 = input.process_input("""
molecule dimer {
0 1
C        3.99802556      -0.12568606       0.00000000
H        5.10246875      -0.16040648       0.00000000
N        3.22763840      -1.16876965       0.00000000
H        3.80888196      -2.01363519       0.00000000
N        3.48433462       1.15076117       0.00000000
H        2.47917010       1.27706322       0.00000000
H        4.08677787       1.95961623       0.00000000
--
0 1
C       -3.99802556       0.12568606       0.00000000
H       -5.10246875       0.16040648       0.00000000
N       -3.22763840       1.16876965       0.00000000
H       -3.80888196       2.01363519       0.00000000
N       -3.48433462      -1.15076117       0.00000000
H       -2.47917010      -1.27706322       0.00000000
H       -4.08677787      -1.95961623       0.00000000
units angstrom
}
""")

HBC1_FaNNFaNN_10p0 = input.process_input("""
molecule dimer {
0 1
C        4.99855171      -0.12037120       0.00000000
H        6.10329053      -0.14697466       0.00000000
N        4.23217198      -1.16574424       0.00000000
H        4.81516518      -2.00932380       0.00000000
N        4.47966407       1.15461725       0.00000000
H        3.47428778       1.27673150       0.00000000
H        5.07887010       1.96589319       0.00000000
--
0 1
C       -4.99855171       0.12037120       0.00000000
H       -6.10329053       0.14697466       0.00000000
N       -4.23217198       1.16574424       0.00000000
H       -4.81516518       2.00932380       0.00000000
N       -4.47966407      -1.15461725       0.00000000
H       -3.47428778      -1.27673150       0.00000000
H       -5.07887010      -1.96589319       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p4 = input.process_input("""
molecule dimer {
0 1
C       -1.68442643       0.20364693       0.00000000
H       -2.78917270       0.34190070       0.00000000
O       -1.02504588       1.24486390       0.00000000
O       -1.38766062      -1.06424798       0.00000000
H       -0.36064758      -1.17759060       0.00000000
--
0 1
C        1.68925857      -0.21855377       0.00000000
H        2.79602984      -0.35706095       0.00000000
O        0.99264772      -1.25039598       0.00000000
N        1.44623488       1.09491586       0.00000000
H        0.47728536       1.44199769       0.00000000
H        2.26183892       1.69328776       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p5 = input.process_input("""
molecule dimer {
0 1
C       -1.73348582       0.18223355       0.00000000
H       -2.83797096       0.30364160       0.00000000
O       -1.07653282       1.22594575       0.00000000
O       -1.38752757      -1.07701453       0.00000000
H       -0.36204548      -1.15331467       0.00000000
--
0 1
C        1.74555921      -0.20019268       0.00000000
H        2.85267994      -0.32189045       0.00000000
O        1.03958692      -1.22627875       0.00000000
N        1.44809759       1.10348563       0.00000000
H        0.46246915       1.40389863       0.00000000
H        2.22816418       1.74757255       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p6 = input.process_input("""
molecule dimer {
0 1
C       -1.78217824       0.17067463       0.00000000
H       -2.88586235       0.28192476       0.00000000
O       -1.11959419       1.21087120       0.00000000
O       -1.40116395      -1.08369245       0.00000000
H       -0.38043759      -1.13241294       0.00000000
--
0 1
C        1.79967181      -0.19037187       0.00000000
H        2.90662488      -0.30195151       0.00000000
O        1.08054052      -1.20751884       0.00000000
N        1.46296486       1.10656646       0.00000000
H        0.46794114       1.37317096       0.00000000
H        2.21379368       1.78455973       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p7 = input.process_input("""
molecule dimer {
0 1
C       -1.83088576       0.16446589       0.00000000
H       -2.93370803       0.26905294       0.00000000
O       -1.15951512       1.19856960       0.00000000
O       -1.42420415      -1.08772422       0.00000000
H       -0.40924981      -1.11532142       0.00000000
--
0 1
C        1.85258766      -0.18485933       0.00000000
H        2.95926374      -0.28981186       0.00000000
O        1.12048297      -1.19270075       0.00000000
N        1.48618925       1.10754437       0.00000000
H        0.48561914       1.34851310       0.00000000
H        2.21266192       1.81160865       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p8 = input.process_input("""
molecule dimer {
0 1
C       -1.87970385       0.16163715       0.00000000
H       -2.98177538       0.26162942       0.00000000
O       -1.19849611       1.18849760       0.00000000
O       -1.45487899      -1.09052331       0.00000000
H       -0.44572390      -1.10219203       0.00000000
--
0 1
C        1.90475161      -0.18173103       0.00000000
H        3.01119328      -0.28211981       0.00000000
O        1.16164142      -1.18119246       0.00000000
N        1.51565227       1.10780268       0.00000000
H        0.51184265       1.32923289       0.00000000
H        2.22211656       1.83196397       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_3p9 = input.process_input("""
molecule dimer {
0 1
C       -1.92862050       0.16119369       0.00000000
H       -3.03011644       0.25790262       0.00000000
O       -1.23743130       1.18028070       0.00000000
O       -1.49236998      -1.09274175       0.00000000
H       -0.48861382      -1.09313129       0.00000000
--
0 1
C        1.95643502      -0.17990567       0.00000000
H        3.06274786      -0.27703752       0.00000000
O        1.20521341      -1.17264505       0.00000000
N        1.55000583       1.10795923       0.00000000
H        0.54440542       1.31463815       0.00000000
H        2.24038620       1.84747610       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p0 = input.process_input("""
molecule dimer {
0 1
C       -1.97758328       0.16248610       0.00000000
H       -3.07869694       0.25672827       0.00000000
O       -1.27661752       1.17355076       0.00000000
O       -1.53599418      -1.09470560       0.00000000
H       -0.53703710      -1.08809981       0.00000000
--
0 1
C        2.00784676      -0.17861908       0.00000000
H        3.11415420      -0.27330577       0.00000000
O        1.25185521      -1.16673528       0.00000000
N        1.58801714       1.10829528       0.00000000
H        0.58149375       1.30380344       0.00000000
H        2.26561575       1.85952935       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p1 = input.process_input("""
molecule dimer {
0 1
C       -2.02654503       0.16495904       0.00000000
H       -3.12746053       0.25719799       0.00000000
O       -1.31615809       1.16794329       0.00000000
O       -1.58476587      -1.09656137       0.00000000
H       -0.58990682      -1.08671689       0.00000000
--
0 1
C        2.05914052      -0.17735549       0.00000000
H        3.16555405      -0.27005509       0.00000000
O        1.30164029      -1.16302967       0.00000000
N        1.62857259       1.10890330       0.00000000
H        0.62159334       1.29571943       0.00000000
H        2.29594551       1.86920816       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p2 = input.process_input("""
molecule dimer {
0 1
C       -2.07549477       0.16808404       0.00000000
H       -3.17636801       0.25855097       0.00000000
O       -1.35622177       1.16316834       0.00000000
O       -1.63731502      -1.09834467       0.00000000
H       -0.64583681      -1.08819660       0.00000000
--
0 1
C        2.11039579      -0.17590164       0.00000000
H        3.21699962      -0.26683949       0.00000000
O        1.35401341      -1.16097028       0.00000000
N        1.67091982       1.10975954       0.00000000
H        0.66372245       1.28959945       0.00000000
H        2.33000288       1.87722287       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p3 = input.process_input("""
molecule dimer {
0 1
C       -2.12445982       0.17139548       0.00000000
H       -3.22540810       0.26019471       0.00000000
O       -1.39708971       1.15905225       0.00000000
O       -1.69213335      -1.10003833       0.00000000
H       -0.70338392      -1.09153044       0.00000000
--
0 1
C        2.16162183      -0.17430715       0.00000000
H        3.26847086      -0.26358230       0.00000000
O        1.40801417      -1.15997383       0.00000000
N        1.71474383       1.11078804       0.00000000
H        0.70747623       1.28502211       0.00000000
H        2.36711721       1.88391930       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p4 = input.process_input("""
molecule dimer {
0 1
C       -2.17348283       0.17454732       0.00000000
H       -3.27458485       0.26172858       0.00000000
O       -1.43904397       1.15551433       0.00000000
O       -1.74791159      -1.10161980       0.00000000
H       -0.76134836      -1.09577404       0.00000000
--
0 1
C        2.21279083      -0.17274192       0.00000000
H        3.31991768      -0.26040021       0.00000000
O        1.46265433      -1.15955428       0.00000000
N        1.75999690       1.11191304       0.00000000
H        0.75279213       1.28183397       0.00000000
H        2.40708847       1.88942849       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p6 = input.process_input("""
molecule dimer {
0 1
C       -2.27180865       0.17975090       0.00000000
H       -3.37333986       0.26380020       0.00000000
O       -1.52660441       1.14998794       0.00000000
O       -1.85932830      -1.10443098       0.00000000
H       -0.87597638      -1.10468982       0.00000000
--
0 1
C        2.31485963      -0.17022217       0.00000000
H        3.42259050      -0.25474451       0.00000000
O        1.57138429      -1.15932816       0.00000000
N        1.85465707       1.11426632       0.00000000
H        0.84799242       1.27930626       0.00000000
H        2.49507804       1.89722164       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_4p8 = input.process_input("""
molecule dimer {
0 1
C       -2.37053947       0.18335144       0.00000000
H       -3.47255698       0.26439714       0.00000000
O       -1.61858011       1.14618338       0.00000000
O       -1.96882669      -1.10687094       0.00000000
H       -0.98773609      -1.11271132       0.00000000
--
0 1
C        2.41653317      -0.16870436       0.00000000
H        3.52489949      -0.25021697       0.00000000
O        1.67809174      -1.15922822       0.00000000
N        1.95399445       1.11666331       0.00000000
H        0.94841269       1.28120182       0.00000000
H        2.59207415       1.90149404       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_5p0 = input.process_input("""
molecule dimer {
0 1
C       -2.46965740       0.18532900       0.00000000
H       -3.57215472       0.26340833       0.00000000
O       -1.71425854       1.14365130       0.00000000
O       -2.07552238      -1.10913287       0.00000000
H       -1.09609784      -1.11892804       0.00000000
--
0 1
C        2.51785158      -0.16788849       0.00000000
H        3.62685662      -0.24642870       0.00000000
O        1.78274154      -1.15890784       0.00000000
N        2.05645986       1.11912237       0.00000000
H        1.05236659       1.28627334       0.00000000
H        2.69518733       1.90341149       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_5p4 = input.process_input("""
molecule dimer {
0 1
C       -2.66850366       0.18514825       0.00000000
H       -3.77181125       0.25702537       0.00000000
O       -1.91234427       1.14028296       0.00000000
O       -2.28050679      -1.11367376       0.00000000
H       -1.30314992      -1.12550085       0.00000000
--
0 1
C        2.72007445      -0.16590118       0.00000000
H        3.83020333      -0.23822269       0.00000000
O        1.98935069      -1.15730079       0.00000000
N        2.26253222       1.12413593       0.00000000
H        1.26110058       1.29739823       0.00000000
H        2.90434733       1.90588628       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_5p8 = input.process_input("""
molecule dimer {
0 1
C       -2.86735875       0.18206412       0.00000000
H       -3.97126829       0.24767419       0.00000000
O       -2.11347561       1.13744992       0.00000000
O       -2.47888370      -1.11819686       0.00000000
H       -1.50260023      -1.12749792       0.00000000
--
0 1
C        2.92242523      -0.16204751       0.00000000
H        4.03336881      -0.22807565       0.00000000
O        2.19639756      -1.15519802       0.00000000
N        2.46432498       1.12882796       0.00000000
H        1.46400695       1.30306213       0.00000000
H        3.10628291       1.91044267       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_6p4 = input.process_input("""
molecule dimer {
0 1
C       -3.16535897       0.17586134       0.00000000
H       -4.26988030       0.23307508       0.00000000
O       -2.41665650       1.13403709       0.00000000
O       -2.77135326      -1.12414764       0.00000000
H       -1.79587783      -1.12655724       0.00000000
--
0 1
C        3.22607312      -0.15521212       0.00000000
H        4.33781868      -0.21280006       0.00000000
O        2.50723536      -1.15220261       0.00000000
N        2.76274851       1.13469399       0.00000000
H        1.76253556       1.30519593       0.00000000
H        3.40173066       1.91872273       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_7p0 = input.process_input("""
molecule dimer {
0 1
C       -3.46303289       0.16994236       0.00000000
H       -4.56795025       0.22032609       0.00000000
O       -2.71952203       1.13148929       0.00000000
O       -3.06229197      -1.12882028       0.00000000
H       -2.08723785      -1.12403398       0.00000000
--
0 1
C        3.52970200      -0.14892316       0.00000000
H        4.64195129      -0.19964122       0.00000000
O        2.81721464      -1.14975890       0.00000000
N        3.05992779       1.13914187       0.00000000
H        2.05915506       1.30417766       0.00000000
H        3.69489864       1.92641219       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_8p0 = input.process_input("""
molecule dimer {
0 1
C       -3.95866293       0.16221564       0.00000000
H       -5.06397335       0.20411777       0.00000000
O       -3.22210995       1.12855078       0.00000000
O       -3.54823170      -1.13429286       0.00000000
H       -2.57357553      -1.11984997       0.00000000
--
0 1
C        4.03559598      -0.14084533       0.00000000
H        5.14833239      -0.18302898       0.00000000
O        3.33114210      -1.14677626       0.00000000
N        3.55651747       1.14422533       0.00000000
H        2.55468792       1.30114823       0.00000000
H        4.18559310       1.93619532       0.00000000
units angstrom
}
""")

HBC1_FaOOFaON_10p0 = input.process_input("""
molecule dimer {
0 1
C       -4.94899199       0.15261986       0.00000000
H       -6.05466697       0.18398552       0.00000000
O       -4.22139968       1.12529943       0.00000000
O       -4.52579681      -1.14045626       0.00000000
H       -3.55156372      -1.11418323       0.00000000
--
0 1
C        5.04698843      -0.13094494       0.00000000
H        6.16016561      -0.16252342       0.00000000
O        4.35243631      -1.14332963       0.00000000
N        4.55565401       1.14984957       0.00000000
H        3.55239498       1.29639409       0.00000000
H        5.17702398       1.94785550       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p4 = input.process_input("""
molecule dimer {
0 1
C       -1.68224527       0.20007686       0.00000000
H       -2.79851036       0.32731687       0.00000000
O       -1.03337695       1.25685782       0.00000000
N       -1.42371554      -1.11002338       0.00000000
H       -0.41476039      -1.44699423       0.00000000
H       -2.24807926      -1.69769344       0.00000000
--
0 1
C        1.69587997      -0.18498638       0.00000000
H        2.80117416      -0.31097584       0.00000000
N        1.05127931      -1.32038822       0.00000000
H        1.74163210      -2.07111657       0.00000000
N        1.41625166       1.12650330       0.00000000
H        0.43768191       1.44670293       0.00000000
H        2.21582070       1.74279349       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p5 = input.process_input("""
molecule dimer {
0 1
C       -1.73161493       0.16885424       0.00000000
H       -2.84774407       0.27229826       0.00000000
O       -1.09203207       1.23345283       0.00000000
N       -1.40935402      -1.12703884       0.00000000
H       -0.38663445      -1.41522920       0.00000000
H       -2.19618815      -1.76421756       0.00000000
--
0 1
C        1.75344970      -0.15414523       0.00000000
H        2.85932238      -0.25663867       0.00000000
N        1.10888334      -1.28995398       0.00000000
H        1.79513417      -2.04501619       0.00000000
N        1.40087130       1.13935068       0.00000000
H        0.40196256       1.40065024       0.00000000
H        2.15615427       1.80923145       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p6 = input.process_input("""
molecule dimer {
0 1
C       -1.78022870       0.14987298       0.00000000
H       -2.89545364       0.23836736       0.00000000
O       -1.13881825       1.21461091       0.00000000
N       -1.41027191      -1.13460552       0.00000000
H       -0.38077495      -1.38453958       0.00000000
H       -2.16426193      -1.81034042       0.00000000
--
0 1
C        1.80849126      -0.13489609       0.00000000
H        2.91397761      -0.22261770       0.00000000
N        1.15496921      -1.26599407       0.00000000
H        1.83228961      -2.02978225       0.00000000
N        1.40124345       1.14388415       0.00000000
H        0.39074481       1.36097240       0.00000000
H        2.11873085       1.85410086       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p7 = input.process_input("""
molecule dimer {
0 1
C       -1.82868510       0.13792020       0.00000000
H       -2.94280654       0.21653844       0.00000000
O       -1.18029155       1.19916421       0.00000000
N       -1.42121847      -1.13750418       0.00000000
H       -0.38876840      -1.35619673       0.00000000
H       -2.14603498      -1.84437855       0.00000000
--
0 1
C        1.86213779      -0.12252360       0.00000000
H        2.96689030      -0.20048072       0.00000000
N        1.19568826      -1.24658561       0.00000000
H        1.86217213      -2.02059464       0.00000000
N        1.41233657       1.14451863       0.00000000
H        0.39516324       1.32705396       0.00000000
H        2.09724900       1.88616308       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p8 = input.process_input("""
molecule dimer {
0 1
C       -1.87720296       0.13037316       0.00000000
H       -2.99024880       0.20231026       0.00000000
O       -1.21984382       1.18650518       0.00000000
N       -1.43964187      -1.13803078       0.00000000
H       -0.40644181      -1.33086253       0.00000000
H       -2.13868738      -1.87029717       0.00000000
--
0 1
C        1.91488588      -0.11471276       0.00000000
H        3.01882141      -0.18606105       0.00000000
N        1.23404120      -1.23071467       0.00000000
H        1.88923247      -2.01504741       0.00000000
N        1.43187244       1.14338561       0.00000000
H        0.41107420       1.29878706       0.00000000
H        2.08912295       1.90958573       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_3p9 = input.process_input("""
molecule dimer {
0 1
C       -1.92586862       0.12572004       0.00000000
H       -3.03796709       0.19311602       0.00000000
O       -1.25948828       1.17632094       0.00000000
N       -1.46411088      -1.13741671       0.00000000
H       -0.43144516      -1.30889080       0.00000000
H       -2.14076353      -1.89029703       0.00000000
--
0 1
C        1.96699000      -0.11019704       0.00000000
H        3.07015648      -0.17705172       0.00000000
N        1.27163617      -1.21778291       0.00000000
H        1.91584742      -2.01186499       0.00000000
N        1.45874144       1.14162823       0.00000000
H        0.43630966       1.27622061       0.00000000
H        2.09327641       1.92665702       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p0 = input.process_input("""
molecule dimer {
0 1
C       -1.97471721       0.12295853       0.00000000
H       -3.08604013       0.18727896       0.00000000
O       -1.30050144       1.16840525       0.00000000
N       -1.49363393      -1.13636823       0.00000000
H       -0.46222530      -1.29037121       0.00000000
H       -2.15122398      -1.90581462       0.00000000
--
0 1
C        2.01860067      -0.10816415       0.00000000
H        3.12111663      -0.17197486       0.00000000
N        1.30940050      -1.20735856       0.00000000
H        1.94334630      -2.01030940       0.00000000
N        1.49223980       1.13988904       0.00000000
H        0.46948873       1.25930607       0.00000000
H        2.10902806       1.93883512       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p1 = input.process_input("""
molecule dimer {
0 1
C       -2.02376044       0.12137130       0.00000000
H       -3.13449129       0.18360767       0.00000000
O       -1.34361035       1.16254229       0.00000000
N       -1.52737885      -1.13528332       0.00000000
H       -0.49758591      -1.27514722       0.00000000
H       -2.16905903      -1.91790148       0.00000000
--
0 1
C        2.06981923      -0.10799978       0.00000000
H        3.17183809      -0.16974799       0.00000000
N        1.34787373      -1.19906629       0.00000000
H        1.97244050      -2.00991317       0.00000000
N        1.53167706       1.13850851       0.00000000
H        0.50945141       1.24775624       0.00000000
H        2.13551288       1.94717878       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p2 = input.process_input("""
molecule dimer {
0 1
C       -2.07299095       0.12047764       0.00000000
H       -3.18330109       0.18128485       0.00000000
O       -1.38900859       1.15844919       0.00000000
N       -1.56461787      -1.13434964       0.00000000
H       -0.53656198      -1.26288353       0.00000000
H       -2.19324500      -1.92734321       0.00000000
--
0 1
C        2.12072533      -0.10919529       0.00000000
H        3.22240596      -0.16952990       0.00000000
N        1.38733366      -1.19254323       0.00000000
H        2.00342544      -2.01034739       0.00000000
N        1.57620126       1.13760159       0.00000000
H        0.55501008       1.24096826       0.00000000
H        2.17143215       1.95251421       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p3 = input.process_input("""
molecule dimer {
0 1
C       -2.12238102       0.12001003       0.00000000
H       -3.23241727       0.17981127       0.00000000
O       -1.43641896       1.15577795       0.00000000
N       -1.60476162      -1.13362638       0.00000000
H       -0.57842142      -1.25321622       0.00000000
H       -2.22287071      -1.93470173       0.00000000
--
0 1
C        2.17139325      -0.11130946       0.00000000
H        3.27288112      -0.17065017       0.00000000
N        1.42791432      -1.18745819       0.00000000
H        2.03640187      -2.01135424       0.00000000
N        1.62476864       1.13714001       0.00000000
H        0.60488448       1.23808048       0.00000000
H        2.21504656       1.95554011       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p4 = input.process_input("""
molecule dimer {
0 1
C       -2.17188140       0.11988045       0.00000000
H       -3.28175966       0.17894887       0.00000000
O       -1.48523553       1.15415106       0.00000000
N       -1.64740473      -1.13310658       0.00000000
H       -0.62269709      -1.24587569       0.00000000
H       -2.25726874      -1.94034051       0.00000000
--
0 1
C        2.22190118      -0.11395941       0.00000000
H        3.32331433      -0.17257731       0.00000000
N        1.46969819      -1.18352562       0.00000000
H        2.07143556      -2.01270689       0.00000000
N        1.67623295       1.13700945       0.00000000
H        0.65776551       1.23805882       0.00000000
H        2.26434801       1.95686864       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p6 = input.process_input("""
molecule dimer {
0 1
C       -2.27098481       0.12071738       0.00000000
H       -3.38079959       0.17880004       0.00000000
O       -1.58441741       1.15276122       0.00000000
N       -1.73927197      -1.13261218       0.00000000
H       -0.71758673      -1.23769913       0.00000000
H       -2.33880523      -1.94731780       0.00000000
--
0 1
C        2.32272918      -0.11969671       0.00000000
H        3.42422935      -0.17734422       0.00000000
N        1.55712717      -1.17831430       0.00000000
H        2.14807243      -2.01572191       0.00000000
N        1.78379346       1.13726371       0.00000000
H        0.76811221       1.24281216       0.00000000
H        2.37265274       1.95644249       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_4p8 = input.process_input("""
molecule dimer {
0 1
C       -2.37001060       0.12298207       0.00000000
H       -3.47998730       0.18042379       0.00000000
O       -1.68320651       1.15259030       0.00000000
N       -1.83827892      -1.13269684       0.00000000
H       -0.81940832      -1.23705071       0.00000000
H       -2.43438758      -1.94974931       0.00000000
--
0 1
C        2.42357564      -0.12508788       0.00000000
H        3.52532822      -0.18210400       0.00000000
N        1.64921165      -1.17564041       0.00000000
H        2.23277319      -2.01853948       0.00000000
N        1.89367698       1.13774591       0.00000000
H        0.88057502       1.25034358       0.00000000
H        2.48658153       1.95391122       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_5p0 = input.process_input("""
molecule dimer {
0 1
C       -2.46891156       0.12607812       0.00000000
H       -3.57918367       0.18293782       0.00000000
O       -1.78100213       1.15283979       0.00000000
N       -1.94185120      -1.13321981       0.00000000
H       -0.92564680      -1.24163694       0.00000000
H       -2.53922912      -1.94924812       0.00000000
--
0 1
C        2.52454538      -0.12964879       0.00000000
H        3.62661434      -0.18608839       0.00000000
N        1.74484282      -1.17452020       0.00000000
H        2.32376341      -2.02083242       0.00000000
N        2.00367207       1.13832154       0.00000000
H        0.99293782       1.25860595       0.00000000
H        2.60169997       1.95068341       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_5p4 = input.process_input("""
molecule dimer {
0 1
C       -2.66627924       0.13246574       0.00000000
H       -3.77728386       0.18752690       0.00000000
O       -1.97376452       1.15272174       0.00000000
N       -2.15448827      -1.13525942       0.00000000
H       -1.14315421      -1.25838806       0.00000000
H       -2.76095466      -1.94449331       0.00000000
--
0 1
C        2.72710220      -0.13482909       0.00000000
H        3.82978184      -0.18947766       0.00000000
N        1.94426866      -1.17439124       0.00000000
H        2.52008027      -2.02303547       0.00000000
N        2.21808889       1.13969445       0.00000000
H        1.21077019       1.27082675       0.00000000
H        2.82358807       1.94643158       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_5p8 = input.process_input("""
molecule dimer {
0 1
C       -2.86354863       0.13635439       0.00000000
H       -3.97519119       0.18862454       0.00000000
O       -2.16694926       1.15157613       0.00000000
N       -2.36429186      -1.13776075       0.00000000
H       -1.35648460      -1.27362001       0.00000000
H       -2.97920858      -1.94060039       0.00000000
--
0 1
C        2.93005123      -0.13606441       0.00000000
H        4.03321099      -0.18793570       0.00000000
N        2.14816261      -1.17431528       0.00000000
H        2.72399657      -2.02298552       0.00000000
N        2.42557457       1.14138666       0.00000000
H        1.42000011       1.27687473       0.00000000
H        3.03409167       1.94580891       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_6p4 = input.process_input("""
molecule dimer {
0 1
C       -3.15990010       0.13729708       0.00000000
H       -4.27217260       0.18465048       0.00000000
O       -2.46166577       1.14950598       0.00000000
N       -2.66867708      -1.14125954       0.00000000
H       -1.66354000      -1.28595015       0.00000000
H       -3.28934347      -1.93968912       0.00000000
--
0 1
C        3.23430882      -0.13492712       0.00000000
H        4.33798775      -0.18191466       0.00000000
N        2.45446506      -1.17291481       0.00000000
H        3.03058764      -2.02136402       0.00000000
N        2.73117559       1.14418616       0.00000000
H        1.72670616       1.28092608       0.00000000
H        3.34065139       1.94785431       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_7p0 = input.process_input("""
molecule dimer {
0 1
C       -3.45655227       0.13564152       0.00000000
H       -4.56921064       0.17826784       0.00000000
O       -2.75931575       1.14757333       0.00000000
N       -2.96668321      -1.14421239       0.00000000
H       -1.96251864      -1.29094605       0.00000000
H       -3.58850495      -1.94176941       0.00000000
--
0 1
C        3.53831767      -0.13233432       0.00000000
H        4.64233959      -0.17462977       0.00000000
N        2.76094801      -1.17118270       0.00000000
H        3.33794153      -2.01900317       0.00000000
N        3.03335964       1.14677845       0.00000000
H        2.02910179       1.28169924       0.00000000
H        3.64162361       1.95135046       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_8p0 = input.process_input("""
molecule dimer {
0 1
C       -3.95119391       0.13180726       0.00000000
H       -5.06421494       0.16795264       0.00000000
O       -3.25700879       1.14506295       0.00000000
N       -3.45900266      -1.14786697       0.00000000
H       -2.45511963      -1.29313383       0.00000000
H       -4.07950006      -1.94648494       0.00000000
--
0 1
C        4.04459226      -0.12785598       0.00000000
H        5.14896015      -0.16372034       0.00000000
N        3.27114001      -1.16876693       0.00000000
H        3.84998346      -2.01526327       0.00000000
N        3.53511346       1.15017556       0.00000000
H        2.53063281       1.28082664       0.00000000
H        4.14039777       1.95697653       0.00000000
units angstrom
}
""")

HBC1_FaONFaNN_10p0 = input.process_input("""
molecule dimer {
0 1
C       -4.94051459       0.12561969       0.00000000
H       -6.05385655       0.15311771       0.00000000
O       -4.25179826       1.14206138       0.00000000
N       -4.44195445      -1.15218600       0.00000000
H       -3.43769145      -1.29256009       0.00000000
H       -5.05860183      -1.95380737       0.00000000
--
0 1
C        5.05643837      -0.12129136       0.00000000
H        6.16113003      -0.14857573       0.00000000
N        4.28876097      -1.16580235       0.00000000
H        4.87088844      -2.00997853       0.00000000
N        4.53935952       1.15434078       0.00000000
H        3.53425488       1.27830661       0.00000000
H        5.13978851       1.96474473       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_3p6 = input.process_input("""
molecule dimer {
0 1
C       -1.76494924       0.15235779       0.00000000
H       -2.87327871       0.25118115       0.00000000
O       -1.11608332       1.20862564       0.00000000
O       -1.38646683      -1.09127823       0.00000000
H       -0.29340506      -1.17844056       0.00000000
--
0 1
C        1.82082566      -0.16736510       0.00000000
H        2.92191627      -0.26554301       0.00000000
N        1.07144779      -1.23580930       0.00000000
H        1.63171355      -2.08650902       0.00000000
N        1.45046411       1.12181032       0.00000000
H        0.44433636       1.35982341       0.00000000
H        2.18036640       1.81973064       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_3p7 = input.process_input("""
molecule dimer {
0 1
C       -1.81131260       0.15832706       0.00000000
H       -2.91739607       0.25779849       0.00000000
O       -1.14367167       1.19940056       0.00000000
O       -1.41977566      -1.09000937       0.00000000
H       -0.35849776      -1.15413164       0.00000000
--
0 1
C        1.87381606      -0.17308101       0.00000000
H        2.97427347      -0.27204647       0.00000000
N        1.09844368      -1.22316593       0.00000000
H        1.64649981      -2.08331650       0.00000000
N        1.48635545       1.11647699       0.00000000
H        0.48139283       1.34214247       0.00000000
H        2.19917472       1.83148239       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_3p8 = input.process_input("""
molecule dimer {
0 1
C       -1.85851203       0.15849692       0.00000000
H       -2.96318957       0.25525336       0.00000000
O       -1.17920577       1.19001178       0.00000000
O       -1.44984324      -1.09154980       0.00000000
H       -0.40704328      -1.13449846       0.00000000
--
0 1
C        1.92699581      -0.17306795       0.00000000
H        3.02703743      -0.26941835       0.00000000
N        1.13373794      -1.21065729       0.00000000
H        1.67232278      -2.07790126       0.00000000
N        1.51877490       1.11463167       0.00000000
H        0.51255074       1.32398473       0.00000000
H        2.21287037       1.84762526       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_3p9 = input.process_input("""
molecule dimer {
0 1
C       -1.90625728       0.15789816       0.00000000
H       -3.00992908       0.25165586       0.00000000
O       -1.21798842       1.18183546       0.00000000
O       -1.48225867      -1.09353763       0.00000000
H       -0.45272461      -1.11834425       0.00000000
--
0 1
C        1.97974667      -0.17222063       0.00000000
H        3.07948948      -0.26564455       0.00000000
N        1.17265485      -1.20000472       0.00000000
H        1.70375525      -2.07277955       0.00000000
N        1.55373649       1.11388209       0.00000000
H        0.54641308       1.30898417       0.00000000
H        2.23141500       1.86193091       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p0 = input.process_input("""
molecule dimer {
0 1
C       -1.95441797       0.15753953       0.00000000
H       -3.05736547       0.24858075       0.00000000
O       -1.25897479       1.17508156       0.00000000
O       -1.51819340      -1.09554164       0.00000000
H       -0.49885308      -1.10531756       0.00000000
--
0 1
C        2.03202506      -0.17151568       0.00000000
H        3.13158097      -0.26227695       0.00000000
N        1.21422382      -1.19144395       0.00000000
H        1.73970852      -2.06839079       0.00000000
N        1.59233540       1.11373305       0.00000000
H        0.58444574       1.29777092       0.00000000
H        2.25678662       1.87343146       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p1 = input.process_input("""
molecule dimer {
0 1
C       -2.00288568       0.15766229       0.00000000
H       -3.10533408       0.24633212       0.00000000
O       -1.30161674       1.16963551       0.00000000
O       -1.55791984      -1.09749465       0.00000000
H       -0.54668214      -1.09529505       0.00000000
--
0 1
C        2.08391761      -0.17103897       0.00000000
H        3.18340265      -0.25947047       0.00000000
N        1.25831080      -1.18492468       0.00000000
H        1.78008989      -2.06474177       0.00000000
N        1.63434437       1.11404135       0.00000000
H        0.62635840       1.29004572       0.00000000
H        2.28868862       1.88235333       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p2 = input.process_input("""
molecule dimer {
0 1
C       -2.05155512       0.15829850       0.00000000
H       -3.15369321       0.24486622       0.00000000
O       -1.34536667       1.16523819       0.00000000
O       -1.60135034      -1.09941548       0.00000000
H       -0.59665904      -1.08819693       0.00000000
--
0 1
C        2.13554958      -0.17057869       0.00000000
H        3.23507103      -0.25694088       0.00000000
N        1.30498819      -1.18028525       0.00000000
H        1.82500469      -2.06170111       0.00000000
N        1.67894465       1.11473347       0.00000000
H        0.67115531       1.28502625       0.00000000
H        2.32574278       1.88931006       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p3 = input.process_input("""
molecule dimer {
0 1
C       -2.10033250       0.15939803       0.00000000
H       -3.20231249       0.24404256       0.00000000
O       -1.38972180       1.16159856       0.00000000
O       -1.64810100      -1.10130943       0.00000000
H       -0.64870385      -1.08383587       0.00000000
--
0 1
C        2.18703912      -0.16992061       0.00000000
H        3.28668154      -0.25438558       0.00000000
N        1.35416505      -1.17725815       0.00000000
H        1.87421236      -2.05911458       0.00000000
N        1.72522255       1.11572240       0.00000000
H        0.71775960       1.28182144       0.00000000
H        2.36632423       1.89492965       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p4 = input.process_input("""
molecule dimer {
0 1
C       -2.14914999       0.16086986       0.00000000
H       -3.25109115       0.24371672       0.00000000
O       -1.43435230       1.15849662       0.00000000
O       -1.69755858      -1.10315341       0.00000000
H       -0.70240660      -1.08186738       0.00000000
--
0 1
C        2.23846787      -0.16900290       0.00000000
H        3.33829123      -0.25169053       0.00000000
N        1.40543473      -1.17550354       0.00000000
H        1.92694842      -2.05687020       0.00000000
N        1.77253439       1.11690493       0.00000000
H        0.76542789       1.27975581       0.00000000
H        2.40921970       1.89964705       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p6 = input.process_input("""
molecule dimer {
0 1
C       -2.24680225       0.16448460       0.00000000
H       -3.34889590       0.24405206       0.00000000
O       -1.52408934       1.15343926       0.00000000
O       -1.80181709      -1.10653912       0.00000000
H       -0.81274434      -1.08325157       0.00000000
--
0 1
C        2.34125676      -0.16675788       0.00000000
H        3.44152405      -0.24619349       0.00000000
N        1.51181334      -1.17441280       0.00000000
H        2.03873187      -2.05311586       0.00000000
N        1.86930785       1.11945644       0.00000000
H        0.86293152       1.27779567       0.00000000
H        2.49986583       1.90702451       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_4p8 = input.process_input("""
molecule dimer {
0 1
C       -2.34452726       0.16814176       0.00000000
H       -3.44693630       0.24477572       0.00000000
O       -1.61510986       1.14962012       0.00000000
O       -1.90907877      -1.10941692       0.00000000
H       -0.92394094      -1.08850288       0.00000000
--
0 1
C        2.44391781      -0.16472698       0.00000000
H        3.54466133      -0.24124517       0.00000000
N        1.61976719      -1.17478374       0.00000000
H        2.15288994      -2.05011167       0.00000000
N        1.96885825       1.12195665       0.00000000
H        0.96338922       1.27845400       0.00000000
H        2.59632255       1.91191489       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_5p0 = input.process_input("""
molecule dimer {
0 1
C       -2.44239623       0.17109139       0.00000000
H       -3.54516867       0.24497624       0.00000000
O       -1.70804752       1.14675140       0.00000000
O       -2.01618143      -1.11187675       0.00000000
H       -1.03372502      -1.09479015       0.00000000
--
0 1
C        2.54642003      -0.16315517       0.00000000
H        3.64763409      -0.23693560       0.00000000
N        1.72722352      -1.17543768       0.00000000
H        2.26584227      -2.04765017       0.00000000
N        2.07067818       1.12435634       0.00000000
H        1.06632369       1.28126044       0.00000000
H        2.69733971       1.91489993       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_5p4 = input.process_input("""
molecule dimer {
0 1
C       -2.63833181       0.17421223       0.00000000
H       -3.74182732       0.24259469       0.00000000
O       -1.89813960       1.14236018       0.00000000
O       -2.22587715      -1.11630144       0.00000000
H       -1.24674157      -1.10613581       0.00000000
--
0 1
C        2.75133046      -0.15977949       0.00000000
H        3.85338801      -0.22807284       0.00000000
N        1.94079641      -1.17624150       0.00000000
H        2.48828614      -2.04318335       0.00000000
N        2.27517829       1.12893018       0.00000000
H        1.27261881       1.28724241       0.00000000
H        2.90172978       1.91949574       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_5p8 = input.process_input("""
molecule dimer {
0 1
C       -2.83430026       0.17412588       0.00000000
H       -3.93836628       0.23696684       0.00000000
O       -2.09196716       1.13877607       0.00000000
O       -2.42840583      -1.12038969       0.00000000
H       -1.45107051      -1.11332814       0.00000000
--
0 1
C        2.95632855      -0.15546371       0.00000000
H        4.05902283      -0.21822659       0.00000000
N        2.15295059      -1.17586911       0.00000000
H        2.70690065      -2.03878687       0.00000000
N        2.47731363       1.13303523       0.00000000
H        1.47540835       1.28966180       0.00000000
H        3.10222174       1.92484945       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_6p4 = input.process_input("""
molecule dimer {
0 1
C       -3.12831510       0.17084992       0.00000000
H       -4.23295676       0.22611350       0.00000000
O       -2.38703135       1.13479845       0.00000000
O       -2.72349375      -1.12568269       0.00000000
H       -1.74750630      -1.11756461       0.00000000
--
0 1
C        3.26369185      -0.14893273       0.00000000
H        4.36705371      -0.20413229       0.00000000
N        2.46830003      -1.17413526       0.00000000
H        3.02869067      -2.03292199       0.00000000
N        2.77871570       1.13816091       0.00000000
H        1.77677840       1.28984953       0.00000000
H        3.39990002       1.93285794       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_7p0 = input.process_input("""
molecule dimer {
0 1
C       -3.42232127       0.16653527       0.00000000
H       -4.52734265       0.21545352       0.00000000
O       -2.68394725       1.13198626       0.00000000
O       -3.01424928      -1.12993738       0.00000000
H       -2.03892166      -1.11784311       0.00000000
--
0 1
C        3.57083085      -0.14304496       0.00000000
H        4.67461743      -0.19190855       0.00000000
N        2.78161679      -1.17225021       0.00000000
H        3.34675787      -2.02791543       0.00000000
N        3.07956461       1.14215114       0.00000000
H        2.07713948       1.28822281       0.00000000
H        3.69667809       1.93998571       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_8p0 = input.process_input("""
molecule dimer {
0 1
C       -3.91213194       0.16018637       0.00000000
H       -5.01752318       0.20105168       0.00000000
O       -3.17868840       1.12881552       0.00000000
O       -3.49736776      -1.13502496       0.00000000
H       -2.52261475      -1.11615900       0.00000000
--
0 1
C        4.08240815      -0.13536451       0.00000000
H        5.18661439      -0.17618601       0.00000000
N        3.30082536      -1.16974476       0.00000000
H        3.87167108      -2.02158353       0.00000000
N        3.58224397       1.14684614       0.00000000
H        2.57892509       1.28493460       0.00000000
H        4.19355750       1.94910244       0.00000000
units angstrom
}
""")

HBC1_FaOOFaNN_10p0 = input.process_input("""
molecule dimer {
0 1
C       -4.89119351       0.15168886       0.00000000
H       -5.99692859       0.18239314       0.00000000
O       -4.16502480       1.12540358       0.00000000
O       -4.46601853      -1.14087529       0.00000000
H       -3.49177332      -1.11255947       0.00000000
--
0 1
C        5.10495503      -0.12588632       0.00000000
H        6.20955485      -0.15655909       0.00000000
N        4.33254137      -1.16667991       0.00000000
H        4.91019302      -2.01389087       0.00000000
N        4.59328914       1.15218338       0.00000000
H        3.58879137       1.28039339       0.00000000
H        5.19724654       1.95996355       0.00000000
units angstrom
}
""")

HBC1_FaOO_monomer_RLX = input.process_input("""
molecule dimer {
0 1
C       -0.10067338      -0.41790840       0.00000000
H       -0.02994601      -1.52175017       0.00000000
O       -1.13575810       0.21734849       0.00000000
O        1.14827203       0.12334561       0.00000000
H        1.03004150       1.09065136       0.00000000
units angstrom
}
""")

HBC1_FaON_monomer_RLX = input.process_input("""
molecule dimer {
0 1
C       -0.08832415      -0.41231959       0.00000000
H       -0.04284111      -1.52506057       0.00000000
O       -1.14609563       0.21052951       0.00000000
N        1.15566522       0.16652753       0.00000000
H        1.23188425       1.17751832       0.00000000
H        1.99476943      -0.39808693       0.00000000
units angstrom
}
""")

HBC1_FaNN_monomer_RLX = input.process_input("""
molecule dimer {
0 1
C       -0.08463791      -0.42651676       0.00000000
H       -0.04431604      -1.53084724       0.00000000
N       -1.17217255       0.27814579       0.00000000
H       -1.98138857      -0.35160983       0.00000000
N        1.15840794       0.16588754       0.00000000
H        1.22157301       1.17651874       0.00000000
H        2.00315093      -0.38515422       0.00000000
units angstrom
}
""")

#<<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:
   molname = rxnpattern.match(rxn)

   GEOS['%s-%s-dimer'      % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, molname.group(1), re.sub(r'\.', 'p', molname.group(2) )))
   GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, molname.group(1), re.sub(r'\.', 'p', molname.group(2) )))) + monoA_CP
   GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, molname.group(1), re.sub(r'\.', 'p', molname.group(2) )))) + monoB_CP
   GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, molname.group(1), re.sub(r'\.', 'p', molname.group(2) )))) + monoA_unCP
   GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, molname.group(1), re.sub(r'\.', 'p', molname.group(2) )))) + monoB_unCP

GEOS['%s-FaOO-mono-RLX' % (dbse)] = eval('%s_FaOO_monomer_RLX' % (dbse))
GEOS['%s-FaON-mono-RLX' % (dbse)] = eval('%s_FaON_monomer_RLX' % (dbse))
GEOS['%s-FaNN-mono-RLX' % (dbse)] = eval('%s_FaNN_monomer_RLX' % (dbse))

