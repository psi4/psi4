import re
import input

# <<< RGC10 Database Module >>>
# Geometries and Reference interaction energies from
#   Tang et al. JCP 118 4976 (2003).
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
dist = [0.85,0.9,0.95,0.975,1.0,1.025,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.0,2.2]
for d in dist: HeHe.append('HeHe-' + str(d))
for d in dist: HeNe.append('HeNe-' + str(d))
for d in dist: HeAr.append('HeAr-' + str(d))
for d in dist: HeKr.append('HeKr-' + str(d))
for d in dist: NeNe.append('NeNe-' + str(d))
for d in dist: NeAr.append('NeAr-' + str(d))
for d in dist: NeKr.append('NeKr-' + str(d))
for d in dist: ArAr.append('ArAr-' + str(d))
for d in dist: ArKr.append('ArKr-' + str(d))
for d in dist: KrKr.append('KrKr-' + str(d))

temp = [HeHe, HeNe, HeAr, HeKr, NeNe, NeAr, NeKr, ArAr, ArKr, KrKr]
HRXN = sum(temp, [])

HRXN_SM = ['NeNe-1.0','NeNe-1.1','NeAr-0.85']
HRXN_LG = ['KrKr-0.85']
HRXN_EQ = ['HeHe-1.0','HeNe-1.0','HeAr-1.0','HeKr-1.0','NeNe-1.0','NeAr-1.0','NeKr-1.0','ArAr-1.0','ArKr-1.0','KrKr-1.0']

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

# <<< Molecule Specifications >>>
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

RGC1_HeHe_0p85 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  2.533000
units angstrom
}
""", 0)

RGC1_HeHe_0p9 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  2.682000
units angstrom
}
""", 0)

RGC1_HeHe_0p95 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  2.831000
units angstrom
}
""", 0)

RGC1_HeHe_0p975 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  2.905500
units angstrom
}
""", 0)

RGC1_HeHe_1p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  2.980000
units angstrom
}
""", 0)

RGC1_HeHe_1p025 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.054500
units angstrom
}
""", 0)

RGC1_HeHe_1p05 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.129000
units angstrom
}
""", 0)

RGC1_HeHe_1p1 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.278000
units angstrom
}
""", 0)

RGC1_HeHe_1p15 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.427000
units angstrom
}
""", 0)

RGC1_HeHe_1p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.576000
units angstrom
}
""", 0)

RGC1_HeHe_1p3 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  3.874000
units angstrom
}
""", 0)

RGC1_HeHe_1p4 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  4.172000
units angstrom
}
""", 0)

RGC1_HeHe_1p5 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  4.470000
units angstrom
}
""", 0)

RGC1_HeHe_1p6 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  4.768000
units angstrom
}
""", 0)

RGC1_HeHe_1p7 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  5.066000
units angstrom
}
""", 0)

RGC1_HeHe_1p8 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  5.364000
units angstrom
}
""", 0)

RGC1_HeHe_2p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  5.960000
units angstrom
}
""", 0)

RGC1_HeHe_2p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
He 1  R
R =  6.556000
units angstrom
}
""", 0)

RGC1_HeNe_0p85 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  2.592500
units angstrom
}
""", 0)

RGC1_HeNe_0p9 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  2.745000
units angstrom
}
""", 0)

RGC1_HeNe_0p95 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  2.897500
units angstrom
}
""", 0)

RGC1_HeNe_0p975 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  2.973750
units angstrom
}
""", 0)

RGC1_HeNe_1p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.050000
units angstrom
}
""", 0)

RGC1_HeNe_1p025 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.126250
units angstrom
}
""", 0)

RGC1_HeNe_1p05 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.202500
units angstrom
}
""", 0)

RGC1_HeNe_1p1 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.355000
units angstrom
}
""", 0)

RGC1_HeNe_1p15 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.507500
units angstrom
}
""", 0)

RGC1_HeNe_1p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.660000
units angstrom
}
""", 0)

RGC1_HeNe_1p3 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  3.965000
units angstrom
}
""", 0)

RGC1_HeNe_1p4 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  4.270000
units angstrom
}
""", 0)

RGC1_HeNe_1p5 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  4.575000
units angstrom
}
""", 0)

RGC1_HeNe_1p6 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  4.880000
units angstrom
}
""", 0)

RGC1_HeNe_1p7 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  5.185000
units angstrom
}
""", 0)

RGC1_HeNe_1p8 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  5.490000
units angstrom
}
""", 0)

RGC1_HeNe_2p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  6.100000
units angstrom
}
""", 0)

RGC1_HeNe_2p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ne 1  R
R =  6.710000
units angstrom
}
""", 0)

RGC1_HeAr_0p85 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  2.975000
units angstrom
}
""", 0)

RGC1_HeAr_0p9 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.150000
units angstrom
}
""", 0)

RGC1_HeAr_0p95 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.325000
units angstrom
}
""", 0)

RGC1_HeAr_0p975 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.412500
units angstrom
}
""", 0)

RGC1_HeAr_1p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.500000
units angstrom
}
""", 0)

RGC1_HeAr_1p025 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.587500
units angstrom
}
""", 0)

RGC1_HeAr_1p05 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.675000
units angstrom
}
""", 0)

RGC1_HeAr_1p1 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  3.850000
units angstrom
}
""", 0)

RGC1_HeAr_1p15 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  4.025000
units angstrom
}
""", 0)

RGC1_HeAr_1p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  4.200000
units angstrom
}
""", 0)

RGC1_HeAr_1p3 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  4.550000
units angstrom
}
""", 0)

RGC1_HeAr_1p4 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  4.900000
units angstrom
}
""", 0)

RGC1_HeAr_1p5 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  5.250000
units angstrom
}
""", 0)

RGC1_HeAr_1p6 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  5.600000
units angstrom
}
""", 0)

RGC1_HeAr_1p7 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  5.950000
units angstrom
}
""", 0)

RGC1_HeAr_1p8 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  6.300000
units angstrom
}
""", 0)

RGC1_HeAr_2p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  7.000000
units angstrom
}
""", 0)

RGC1_HeAr_2p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Ar 1  R
R =  7.700000
units angstrom
}
""", 0)

RGC1_HeKr_0p85 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.145000
units angstrom
}
""", 0)

RGC1_HeKr_0p9 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.330000
units angstrom
}
""", 0)

RGC1_HeKr_0p95 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.515000
units angstrom
}
""", 0)

RGC1_HeKr_0p975 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.607500
units angstrom
}
""", 0)

RGC1_HeKr_1p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.700000
units angstrom
}
""", 0)

RGC1_HeKr_1p025 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.792500
units angstrom
}
""", 0)

RGC1_HeKr_1p05 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  3.885000
units angstrom
}
""", 0)

RGC1_HeKr_1p1 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  4.070000
units angstrom
}
""", 0)

RGC1_HeKr_1p15 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  4.255000
units angstrom
}
""", 0)

RGC1_HeKr_1p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  4.440000
units angstrom
}
""", 0)

RGC1_HeKr_1p3 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  4.810000
units angstrom
}
""", 0)

RGC1_HeKr_1p4 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  5.180000
units angstrom
}
""", 0)

RGC1_HeKr_1p5 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  5.550000
units angstrom
}
""", 0)

RGC1_HeKr_1p6 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  5.920000
units angstrom
}
""", 0)

RGC1_HeKr_1p7 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  6.290000
units angstrom
}
""", 0)

RGC1_HeKr_1p8 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  6.660000
units angstrom
}
""", 0)

RGC1_HeKr_2p0 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  7.400000
units angstrom
}
""", 0)

RGC1_HeKr_2p2 = input.process_input("""
molecule dimer {
0 1
He
--
0 1
Kr 1  R
R =  8.140000
units angstrom
}
""", 0)

RGC1_NeNe_0p85 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  2.626500
units angstrom
}
""", 0)

RGC1_NeNe_0p9 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  2.781000
units angstrom
}
""", 0)

RGC1_NeNe_0p95 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  2.935500
units angstrom
}
""", 0)

RGC1_NeNe_0p975 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.012750
units angstrom
}
""", 0)

RGC1_NeNe_1p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.090000
units angstrom
}
""", 0)

RGC1_NeNe_1p025 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.167250
units angstrom
}
""", 0)

RGC1_NeNe_1p05 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.244500
units angstrom
}
""", 0)

RGC1_NeNe_1p1 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.399000
units angstrom
}
""", 0)

RGC1_NeNe_1p15 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.553500
units angstrom
}
""", 0)

RGC1_NeNe_1p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  3.708000
units angstrom
}
""", 0)

RGC1_NeNe_1p3 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  4.017000
units angstrom
}
""", 0)

RGC1_NeNe_1p4 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  4.326000
units angstrom
}
""", 0)

RGC1_NeNe_1p5 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  4.635000
units angstrom
}
""", 0)

RGC1_NeNe_1p6 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  4.944000
units angstrom
}
""", 0)

RGC1_NeNe_1p7 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  5.253000
units angstrom
}
""", 0)

RGC1_NeNe_1p8 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  5.562000
units angstrom
}
""", 0)

RGC1_NeNe_2p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  6.180000
units angstrom
}
""", 0)

RGC1_NeNe_2p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ne 1  R
R =  6.798000
units angstrom
}
""", 0)

RGC1_NeAr_0p85 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  2.958000
units angstrom
}
""", 0)

RGC1_NeAr_0p9 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.132000
units angstrom
}
""", 0)

RGC1_NeAr_0p95 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.306000
units angstrom
}
""", 0)

RGC1_NeAr_0p975 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.393000
units angstrom
}
""", 0)

RGC1_NeAr_1p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.480000
units angstrom
}
""", 0)

RGC1_NeAr_1p025 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.567000
units angstrom
}
""", 0)

RGC1_NeAr_1p05 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.654000
units angstrom
}
""", 0)

RGC1_NeAr_1p1 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  3.828000
units angstrom
}
""", 0)

RGC1_NeAr_1p15 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  4.002000
units angstrom
}
""", 0)

RGC1_NeAr_1p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  4.176000
units angstrom
}
""", 0)

RGC1_NeAr_1p3 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  4.524000
units angstrom
}
""", 0)

RGC1_NeAr_1p4 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  4.872000
units angstrom
}
""", 0)

RGC1_NeAr_1p5 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  5.220000
units angstrom
}
""", 0)

RGC1_NeAr_1p6 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  5.568000
units angstrom
}
""", 0)

RGC1_NeAr_1p7 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  5.916000
units angstrom
}
""", 0)

RGC1_NeAr_1p8 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  6.264000
units angstrom
}
""", 0)

RGC1_NeAr_2p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  6.960000
units angstrom
}
""", 0)

RGC1_NeAr_2p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Ar 1  R
R =  7.656000
units angstrom
}
""", 0)

RGC1_NeKr_0p85 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.102500
units angstrom
}
""", 0)

RGC1_NeKr_0p9 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.285000
units angstrom
}
""", 0)

RGC1_NeKr_0p95 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.467500
units angstrom
}
""", 0)

RGC1_NeKr_0p975 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.558750
units angstrom
}
""", 0)

RGC1_NeKr_1p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.650000
units angstrom
}
""", 0)

RGC1_NeKr_1p025 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.741250
units angstrom
}
""", 0)

RGC1_NeKr_1p05 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  3.832500
units angstrom
}
""", 0)

RGC1_NeKr_1p1 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  4.015000
units angstrom
}
""", 0)

RGC1_NeKr_1p15 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  4.197500
units angstrom
}
""", 0)

RGC1_NeKr_1p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  4.380000
units angstrom
}
""", 0)

RGC1_NeKr_1p3 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  4.745000
units angstrom
}
""", 0)

RGC1_NeKr_1p4 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  5.110000
units angstrom
}
""", 0)

RGC1_NeKr_1p5 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  5.475000
units angstrom
}
""", 0)

RGC1_NeKr_1p6 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  5.840000
units angstrom
}
""", 0)

RGC1_NeKr_1p7 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  6.205000
units angstrom
}
""", 0)

RGC1_NeKr_1p8 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  6.570000
units angstrom
}
""", 0)

RGC1_NeKr_2p0 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  7.300000
units angstrom
}
""", 0)

RGC1_NeKr_2p2 = input.process_input("""
molecule dimer {
0 1
Ne
--
0 1
Kr 1  R
R =  8.030000
units angstrom
}
""", 0)

RGC1_ArAr_0p85 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.187500
units angstrom
}
""", 0)

RGC1_ArAr_0p9 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.375000
units angstrom
}
""", 0)

RGC1_ArAr_0p95 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.562500
units angstrom
}
""", 0)

RGC1_ArAr_0p975 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.656250
units angstrom
}
""", 0)

RGC1_ArAr_1p0 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.750000
units angstrom
}
""", 0)

RGC1_ArAr_1p025 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.843750
units angstrom
}
""", 0)

RGC1_ArAr_1p05 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  3.937500
units angstrom
}
""", 0)

RGC1_ArAr_1p1 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  4.125000
units angstrom
}
""", 0)

RGC1_ArAr_1p15 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  4.312500
units angstrom
}
""", 0)

RGC1_ArAr_1p2 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  4.500000
units angstrom
}
""", 0)

RGC1_ArAr_1p3 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  4.875000
units angstrom
}
""", 0)

RGC1_ArAr_1p4 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  5.250000
units angstrom
}
""", 0)

RGC1_ArAr_1p5 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  5.625000
units angstrom
}
""", 0)

RGC1_ArAr_1p6 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  6.000000
units angstrom
}
""", 0)

RGC1_ArAr_1p7 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  6.375000
units angstrom
}
""", 0)

RGC1_ArAr_1p8 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  6.750000
units angstrom
}
""", 0)

RGC1_ArAr_2p0 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  7.500000
units angstrom
}
""", 0)

RGC1_ArAr_2p2 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Ar 1  R
R =  8.250000
units angstrom
}
""", 0)

RGC1_ArKr_0p85 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.306500
units angstrom
}
""", 0)

RGC1_ArKr_0p9 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.501000
units angstrom
}
""", 0)

RGC1_ArKr_0p95 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.695500
units angstrom
}
""", 0)

RGC1_ArKr_0p975 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.792750
units angstrom
}
""", 0)

RGC1_ArKr_1p0 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.890000
units angstrom
}
""", 0)

RGC1_ArKr_1p025 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  3.987250
units angstrom
}
""", 0)

RGC1_ArKr_1p05 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  4.084500
units angstrom
}
""", 0)

RGC1_ArKr_1p1 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  4.279000
units angstrom
}
""", 0)

RGC1_ArKr_1p15 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  4.473500
units angstrom
}
""", 0)

RGC1_ArKr_1p2 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  4.668000
units angstrom
}
""", 0)

RGC1_ArKr_1p3 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  5.057000
units angstrom
}
""", 0)

RGC1_ArKr_1p4 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  5.446000
units angstrom
}
""", 0)

RGC1_ArKr_1p5 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  5.835000
units angstrom
}
""", 0)

RGC1_ArKr_1p6 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  6.224000
units angstrom
}
""", 0)

RGC1_ArKr_1p7 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  6.613000
units angstrom
}
""", 0)

RGC1_ArKr_1p8 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  7.002000
units angstrom
}
""", 0)

RGC1_ArKr_2p0 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  7.780000
units angstrom
}
""", 0)

RGC1_ArKr_2p2 = input.process_input("""
molecule dimer {
0 1
Ar
--
0 1
Kr 1  R
R =  8.558000
units angstrom
}
""", 0)

RGC1_KrKr_0p85 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  3.408500
units angstrom
}
""", 0)

RGC1_KrKr_0p9 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  3.609000
units angstrom
}
""", 0)

RGC1_KrKr_0p95 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  3.809500
units angstrom
}
""", 0)

RGC1_KrKr_0p975 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  3.909750
units angstrom
}
""", 0)

RGC1_KrKr_1p0 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.010000
units angstrom
}
""", 0)

RGC1_KrKr_1p025 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.110250
units angstrom
}
""", 0)

RGC1_KrKr_1p05 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.210500
units angstrom
}
""", 0)

RGC1_KrKr_1p1 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.411000
units angstrom
}
""", 0)

RGC1_KrKr_1p15 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.611500
units angstrom
}
""", 0)

RGC1_KrKr_1p2 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  4.812000
units angstrom
}
""", 0)

RGC1_KrKr_1p3 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  5.213000
units angstrom
}
""", 0)

RGC1_KrKr_1p4 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  5.614000
units angstrom
}
""", 0)

RGC1_KrKr_1p5 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  6.015000
units angstrom
}
""", 0)

RGC1_KrKr_1p6 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  6.416000
units angstrom
}
""", 0)

RGC1_KrKr_1p7 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  6.817000
units angstrom
}
""", 0)

RGC1_KrKr_1p8 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  7.218000
units angstrom
}
""", 0)

RGC1_KrKr_2p0 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  8.020000
units angstrom
}
""", 0)

RGC1_KrKr_2p2 = input.process_input("""
molecule dimer {
0 1
Kr
--
0 1
Kr 1  R
R =  8.822000
units angstrom
}
""", 0)

RGC1_He_monomer = input.process_input("""
molecule monomer {
0 1
He
units angstrom
}
""", 0)

RGC1_Ne_monomer = input.process_input("""
molecule monomer {
0 1
Ne
units angstrom
}
""", 0)

RGC1_Ar_monomer = input.process_input("""
molecule monomer {
0 1
Ar
units angstrom
}
""", 0)

RGC1_Kr_monomer = input.process_input("""
molecule monomer {
0 1
Kr
units angstrom
}
""", 0)

#<<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:
   distance = rxnpattern.match(rxn)

   GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))
   GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoA_CP
   GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoB_CP

GEOS['%s-He-mono-unCP' % (dbse)] = eval('%s_He_monomer' % (dbse))
GEOS['%s-Ne-mono-unCP' % (dbse)] = eval('%s_Ne_monomer' % (dbse))
GEOS['%s-Ar-mono-unCP' % (dbse)] = eval('%s_Ar_monomer' % (dbse))
GEOS['%s-Kr-mono-unCP' % (dbse)] = eval('%s_Kr_monomer' % (dbse))

