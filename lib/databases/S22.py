import input

# <<< S22 Database Module >>>
# Geometries from Jurecka et al. PCCP 8 1985 (2006).
# Reference interaction energies from the following articles:
#   S220: Jurecka et al. PCCP 8 1985 (2006).
#   S22A: Takatani et al. JCP 132 144104 (2010).
#   S22B: Marshall et al. JCP 135 194102 (2011).  *** DEFAULT ***
dbse = 'S22'

# <<< Database Members >>>
HRXN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
HRXN_SM = [2, 8, 16]
HRXN_LG = [15]

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

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

S22_1 = input.process_input("""
molecule dimer {
0 1
N  -1.578718  -0.046611   0.000000
H  -2.158621   0.136396  -0.809565
H  -2.158621   0.136396   0.809565
H  -0.849471   0.658193   0.000000
--
0 1
N   1.578718   0.046611   0.000000
H   2.158621  -0.136396  -0.809565
H   0.849471  -0.658193   0.000000
H   2.158621  -0.136396   0.809565
units angstrom
}
""", 0)

S22_2 = input.process_input("""
molecule dimer {
0 1
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
--
0 1
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
units angstrom
}
""", 0)

S22_3 = input.process_input("""
molecule dimer {
0 1
C  -1.888896  -0.179692   0.000000
O  -1.493280   1.073689   0.000000
O  -1.170435  -1.166590   0.000000
H  -2.979488  -0.258829   0.000000
H  -0.498833   1.107195   0.000000
--
0 1
C   1.888896   0.179692   0.000000
O   1.493280  -1.073689   0.000000
O   1.170435   1.166590   0.000000
H   2.979488   0.258829   0.000000
H   0.498833  -1.107195   0.000000
units angstrom
}
""", 0)

S22_4 = input.process_input("""
molecule dimer {
0 1
C  -2.018649   0.052883   0.000000
O  -1.452200   1.143634   0.000000
N  -1.407770  -1.142484   0.000000
H  -1.964596  -1.977036   0.000000
H  -0.387244  -1.207782   0.000000
H  -3.117061  -0.013701   0.000000
--
0 1
C   2.018649  -0.052883   0.000000
O   1.452200  -1.143634   0.000000
N   1.407770   1.142484   0.000000
H   1.964596   1.977036   0.000000
H   0.387244   1.207782   0.000000
H   3.117061   0.013701   0.000000
units angstrom
}
""", 0)

S22_5 = input.process_input("""
molecule dimer {
0 1
O    -1.4663316    1.0121693    0.0000000
C    -0.6281464    1.9142678    0.0000000
N     0.7205093    1.6882688    0.0000000
C     1.6367290    2.7052764    0.0000000
C     1.2769036    4.0061763    0.0000000
C    -0.1286005    4.3621549    0.0000000
N    -0.9777230    3.2396433    0.0000000
O    -0.5972229    5.4864066    0.0000000
H     2.0103504    4.7938642    0.0000000
H     1.0232515    0.7061820    0.0000000
H    -1.9700268    3.4323850    0.0000000
H     2.6690620    2.3883417    0.0000000
--
0 1
O     1.4663316   -1.0121693    0.0000000
C     0.6281464   -1.9142678    0.0000000
N    -0.7205093   -1.6882688    0.0000000
C    -1.6367290   -2.7052764    0.0000000
C    -1.2769036   -4.0061763    0.0000000
C     0.1286005   -4.3621549    0.0000000
N     0.9777230   -3.2396433    0.0000000
O     0.5972229   -5.4864066    0.0000000
H    -2.0103504   -4.7938642    0.0000000
H    -1.0232515   -0.7061820    0.0000000
H     1.9700268   -3.4323850    0.0000000
H    -2.6690620   -2.3883417    0.0000000
units angstrom
}
""", 0)

S22_6 = input.process_input("""
molecule dimer {
0 1
O    -1.3976213   -1.8858368   -0.3673061
N    -1.4642550    0.3641828    0.0192301
C    -4.1857398    0.3696669    0.0360960
C    -3.4832598    1.5783111    0.2500752
C    -2.1179502    1.5307048    0.2338383
C    -2.0773833   -0.8637492   -0.1899414
C    -3.5156032   -0.8051950   -0.1757585
H    -5.2678045    0.3707428    0.0411419
H    -3.9920334    2.5127560    0.4214414
H    -1.4929196    2.3984096    0.3885018
H    -4.0401226   -1.7348452   -0.3379269
H    -0.4265266    0.3612127    0.0073538
--
0 1
N     1.4327616    0.3639703   -0.0159508
C     2.1154200   -0.7803450    0.1681099
C     3.5237586   -0.8016096    0.1545027
C     4.2185897    0.3735783   -0.0525929
C     3.5099708    1.5615014   -0.2449763
C     2.1280138    1.4953324   -0.2175374
H     4.0459206   -1.7361356    0.3076883
H     5.2999426    0.3666009   -0.0663349
H     4.0110923    2.5024313   -0.4130052
H     1.5339878    2.3893837   -0.3670565
N     1.3883123   -1.9083038    0.4198149
H     1.8694714   -2.7812773    0.2940385
H     0.4089067   -1.9079942    0.1300860
units angstrom
}
""", 0)

S22_7 = input.process_input("""
molecule dimer {
0 1
N     0.9350155   -0.0279801   -0.3788916
C     1.6739638   -0.0357766    0.7424316
C     3.0747955   -0.0094480    0.5994562
C     3.5646109    0.0195446   -0.7059872
N     2.8531510    0.0258031   -1.8409596
C     1.5490760    0.0012569   -1.5808009
N     4.0885824   -0.0054429    1.5289786
C     5.1829921    0.0253971    0.7872176
N     4.9294871    0.0412404   -0.5567274
N     1.0716177   -0.0765366    1.9391390
H     0.8794435    0.0050260   -2.4315709
H     6.1882591    0.0375542    1.1738824
H     5.6035368    0.0648755   -1.3036811
H     0.0586915   -0.0423765    2.0039181
H     1.6443796   -0.0347395    2.7619159
--
0 1
N    -3.9211729   -0.0009646   -1.5163659
C    -4.6136833    0.0169051   -0.3336520
C    -3.9917387    0.0219348    0.8663338
C    -2.5361367    0.0074651    0.8766724
N    -1.9256484   -0.0110593   -0.3638948
C    -2.5395897   -0.0149474   -1.5962357
C    -4.7106131    0.0413373    2.1738637
O    -1.8674730    0.0112093    1.9120833
O    -1.9416783   -0.0291878   -2.6573783
H    -4.4017172   -0.0036078   -2.4004924
H    -0.8838255   -0.0216168   -0.3784269
H    -5.6909220    0.0269347   -0.4227183
H    -4.4439282   -0.8302573    2.7695655
H    -4.4267056    0.9186178    2.7530256
H    -5.7883971    0.0505530    2.0247280
units angstrom
}
""", 0)

S22_8 = input.process_input("""
molecule dimer {
0 1
C   0.000000  -0.000140   1.859161
H  -0.888551   0.513060   1.494685
H   0.888551   0.513060   1.494685
H   0.000000  -1.026339   1.494868
H   0.000000   0.000089   2.948284
--
0 1
C   0.000000   0.000140  -1.859161
H   0.000000  -0.000089  -2.948284
H  -0.888551  -0.513060  -1.494685
H   0.888551  -0.513060  -1.494685
H   0.000000   1.026339  -1.494868
units angstrom
}
""", 0)

S22_9 = input.process_input("""
molecule dimer {
0 1
C  -0.471925  -0.471925  -1.859111
C   0.471925   0.471925  -1.859111
H  -0.872422  -0.872422  -0.936125
H   0.872422   0.872422  -0.936125
H  -0.870464  -0.870464  -2.783308
H   0.870464   0.870464  -2.783308
--
0 1
C  -0.471925   0.471925   1.859111
C   0.471925  -0.471925   1.859111
H  -0.872422   0.872422   0.936125
H   0.872422  -0.872422   0.936125
H  -0.870464   0.870464   2.783308
H   0.870464  -0.870464   2.783308
units angstrom
}
""", 0)

S22_10 = input.process_input("""
molecule dimer {
0 1
C     1.3932178    0.0362913   -0.6332803
C     0.7280364   -1.1884015   -0.6333017
C    -0.6651797   -1.2247077   -0.6332803
C    -1.3932041   -0.0362972   -0.6333017
C    -0.7280381    1.1884163   -0.6332803
C     0.6651677    1.2246987   -0.6333017
H     2.4742737    0.0644484   -0.6317240
H     1.2929588   -2.1105409   -0.6317401
H    -1.1813229   -2.1750081   -0.6317240
H    -2.4742614   -0.0644647   -0.6317401
H    -1.2929508    2.1105596   -0.6317240
H     1.1813026    2.1750056   -0.6317401
--
0 1
C     0.0000000    0.0000000    3.0826195
H     0.5868776    0.8381742    3.4463772
H    -1.0193189    0.0891638    3.4463772
H     0.0000000    0.0000000    1.9966697
H     0.4324413   -0.9273380    3.4463772
units angstrom
}
""", 0)

S22_11 = input.process_input("""
molecule dimer {
0 1
C    -1.0478252   -1.4216736    0.0000000
C    -1.4545034   -0.8554459    1.2062048
C    -1.4545034   -0.8554459   -1.2062048
C    -2.2667970    0.2771610    1.2069539
C    -2.6714781    0.8450211    0.0000000
C    -2.2667970    0.2771610   -1.2069539
H    -1.1338534   -1.2920593   -2.1423150
H    -2.5824943    0.7163066   -2.1437977
H    -3.3030422    1.7232700    0.0000000
H    -2.5824943    0.7163066    2.1437977
H    -1.1338534   -1.2920593    2.1423150
H    -0.4060253   -2.2919049    0.0000000
--
0 1
C     1.0478252    1.4216736    0.0000000
C     1.4545034    0.8554459   -1.2062048
C     1.4545034    0.8554459    1.2062048
C     2.2667970   -0.2771610   -1.2069539
C     2.6714781   -0.8450211    0.0000000
C     2.2667970   -0.2771610    1.2069539
H     0.4060253    2.2919049    0.0000000
H     1.1338534    1.2920593    2.1423150
H     2.5824943   -0.7163066    2.1437977
H     3.3030422   -1.7232700    0.0000000
H     2.5824943   -0.7163066   -2.1437977
H     1.1338534    1.2920593   -2.1423150
units angstrom
}
""", 0)

S22_12 = input.process_input("""
molecule dimer {
0 1
C    -1.2471894   -1.1718212   -0.6961388
C    -1.2471894   -1.1718212    0.6961388
N    -0.2589510   -1.7235771    1.4144796
C     0.7315327   -2.2652221    0.6967288
C     0.7315327   -2.2652221   -0.6967288
N    -0.2589510   -1.7235771   -1.4144796
H    -2.0634363   -0.7223199   -1.2472797
H    -2.0634363   -0.7223199    1.2472797
H     1.5488004   -2.7128282    1.2475604
H     1.5488004   -2.7128282   -1.2475604
--
0 1
C    -0.3380031    2.0800608    1.1300452
C     0.8540254    1.3593471    1.1306308
N     1.4701787    0.9907598    0.0000000
C     0.8540254    1.3593471   -1.1306308
C    -0.3380031    2.0800608   -1.1300452
N    -0.9523059    2.4528836    0.0000000
H    -0.8103758    2.3643033    2.0618643
H     1.3208583    1.0670610    2.0623986
H     1.3208583    1.0670610   -2.0623986
H    -0.8103758    2.3643033   -2.0618643
units angstrom
}
""", 0)

S22_13 = input.process_input("""
molecule dimer {
0 1
N     2.0113587   -1.2132073   -0.0980673
C     2.0257076   -0.6971797   -1.3644029
H     2.2975208   -1.3910592   -2.1456459
C     1.7145226    0.5919651   -1.6124892
H     1.7272873    0.9908466   -2.6120050
C     1.3089605    1.4575340   -0.5205890
O     0.9205926    2.6110864   -0.6260457
N     1.3768885    0.8397454    0.7346356
H     1.0518040    1.3862229    1.5233710
C     1.6459909   -0.4852113    1.0187267
O     1.5611090   -0.9718061    2.1298059
H     2.1294635   -2.2015046    0.0568134
--
0 1
N    -2.0113587    1.2132073   -0.0980673
C    -2.0257076    0.6971797   -1.3644029
H    -2.2975208    1.3910592   -2.1456459
C    -1.7145226   -0.5919651   -1.6124892
H    -1.7272873   -0.9908466   -2.6120050
C    -1.3089605   -1.4575340   -0.5205890
O    -0.9205926   -2.6110864   -0.6260457
N    -1.3768885   -0.8397454    0.7346356
H    -1.0518040   -1.3862229    1.5233710
C    -1.6459909    0.4852113    1.0187267
O    -1.5611090    0.9718061    2.1298059
H    -2.1294635    2.2015046    0.0568134
units angstrom
}
""", 0)

S22_14 = input.process_input("""
molecule dimer {
0 1
C    -0.0210742    1.5318615   -1.3639345
C    -1.2746794    0.9741030   -1.6074097
C    -1.3783055   -0.2256981   -2.3084154
C    -0.2289426   -0.8664053   -2.7687944
C     1.0247882   -0.3035171   -2.5312410
C     1.1289996    0.8966787   -1.8299830
H     0.0600740    2.4565627   -0.8093957
H    -2.1651002    1.4654521   -1.2405676
H    -2.3509735   -0.6616122   -2.4926698
H    -0.3103419   -1.7955762   -3.3172704
H     1.9165847   -0.7940845   -2.8993942
H     2.1000347    1.3326757   -1.6400420
--
0 1
H    -2.9417647    0.8953834    2.2239054
C    -2.0220674    0.4258540    1.9013549
C    -0.8149418    1.0740453    2.1066982
H    -0.7851529    2.0443812    2.5856086
C     0.3704286    0.4492852    1.6847458
C     1.7508619    0.8038935    1.7194004
H     2.1870108    1.6998281    2.1275903
C     2.4451359   -0.2310742    1.1353313
N     1.5646462   -1.2137812    0.7555384
C     0.2861214   -0.8269486    1.0618752
C    -0.9284667   -1.4853121    0.8606937
H    -0.9729200   -2.4554847    0.3834013
C    -2.0792848   -0.8417668    1.2876443
H    -3.0389974   -1.3203846    1.1468400
H     1.8075741   -2.0366963    0.2333038
H     3.5028794   -0.3485344    0.9695233
units angstrom
}
""", 0)

S22_15 = input.process_input("""
molecule dimer {
0 1
N     0.2793014    2.4068393   -0.6057517
C    -1.0848570    2.4457461   -0.5511608
H    -1.6594403    3.0230294   -1.2560905
N    -1.5977117    1.7179877    0.4287543
C    -0.4897255    1.1714358    1.0301910
C    -0.3461366    0.2914710    2.1172343
N    -1.4187090   -0.1677767    2.8101441
H    -1.2388750   -0.9594802    3.4047578
H    -2.2918734   -0.1788223    2.3073619
N     0.8857630   -0.0700763    2.4919494
C     1.9352348    0.4072878    1.7968022
H     2.9060330    0.0788414    2.1458181
N     1.9409775    1.2242019    0.7402202
C     0.6952186    1.5779858    0.4063984
H     0.8610073    2.8298045   -1.3104502
--
0 1
N     1.2754606   -0.6478993   -1.9779104
C     1.4130533   -1.5536850   -0.9550667
H     2.4258769   -1.8670780   -0.7468778
C     0.3575976   -2.0239499   -0.2530575
C     0.4821292   -3.0179494    0.8521221
H     0.1757705   -2.5756065    1.7986281
H    -0.1601691   -3.8770412    0.6639498
H     1.5112443   -3.3572767    0.9513659
C    -0.9684711   -1.5298112   -0.5939792
O    -2.0029280   -1.8396957   -0.0199453
N    -0.9956916   -0.6383870   -1.6720420
H    -1.9014057   -0.2501720   -1.8985760
C     0.0684702   -0.1191762   -2.3763759
O    -0.0397875    0.7227006   -3.2531083
H     2.0853289   -0.2760176   -2.4454577
units angstrom
}
""", 0)

S22_16 = input.process_input("""
molecule dimer {
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
}
""", 0)

S22_17 = input.process_input("""
molecule dimer {
0 1
C     0.7806117   -0.6098875   -1.2075426
C     0.4784039    0.7510406   -1.2079040
C     0.3276592    1.4318573    0.0000000
C     0.4784039    0.7510406    1.2079040
C     0.7806117   -0.6098875    1.2075426
C     0.9321510   -1.2899614    0.0000000
H     0.8966688   -1.1376051   -2.1441482
H     0.3573895    1.2782091   -2.1440546
H     0.0918593    2.4871407    0.0000000
H     0.3573895    1.2782091    2.1440546
H     0.8966688   -1.1376051    2.1441482
H     1.1690064   -2.3451668    0.0000000
--
0 1
O    -2.7885270   -0.2744854    0.0000000
H    -2.6229114   -1.2190831    0.0000000
H    -1.9015103    0.0979110    0.0000000
units angstrom
}
""", 0)

S22_18 = input.process_input("""
molecule dimer {
0 1
C    -0.7392810    0.5158785   -1.2071079
C    -1.4261442    0.3965455    0.0000000
C    -0.7392810    0.5158785    1.2071079
C     0.6342269    0.7546398    1.2070735
C     1.3210434    0.8737566    0.0000000
C     0.6342269    0.7546398   -1.2070735
H    -1.2719495    0.4206316   -2.1432894
H    -2.4902205    0.2052381    0.0000000
H    -1.2719495    0.4206316    2.1432894
H     1.1668005    0.8474885    2.1436950
H     2.3863585    1.0596312    0.0000000
H     1.1668005    0.8474885   -2.1436950
--
0 1
N     0.1803930   -2.9491231    0.0000000
H     0.7595495   -3.1459477   -0.8060729
H     0.7595495   -3.1459477    0.8060729
H     0.0444167   -1.9449399    0.0000000
units angstrom
}
""", 0)

S22_19 = input.process_input("""
molecule dimer {
0 1
C    -0.7097741   -0.9904230    1.2077018
C    -1.4065340   -0.9653529    0.0000000
C    -0.7097741   -0.9904230   -1.2077018
C     0.6839651   -1.0405105   -1.2078652
C     1.3809779   -1.0655522    0.0000000
C     0.6839651   -1.0405105    1.2078652
H    -1.2499482   -0.9686280    2.1440507
H    -2.4869197   -0.9237060    0.0000000
H    -1.2499482   -0.9686280   -2.1440507
H     1.2242882   -1.0580753   -2.1442563
H     2.4615886   -1.1029818    0.0000000
H     1.2242882   -1.0580753    2.1442563
--
0 1
N    -0.0034118    3.5353926    0.0000000
C     0.0751963    2.3707040    0.0000000
H     0.1476295    1.3052847    0.0000000
units angstrom
}
""", 0)

S22_20 = input.process_input("""
molecule dimer {
0 1
C     0.0000000    0.0000000    1.0590353
C     0.0000000   -1.2060084    1.7576742
C     0.0000000   -1.2071767    3.1515905
C     0.0000000    0.0000000    3.8485751
C     0.0000000    1.2071767    3.1515905
C     0.0000000    1.2060084    1.7576742
H     0.0000000    0.0000000   -0.0215805
H     0.0000000   -2.1416387    1.2144217
H     0.0000000   -2.1435657    3.6929953
H     0.0000000    0.0000000    4.9301499
H     0.0000000    2.1435657    3.6929953
H     0.0000000    2.1416387    1.2144217
--
0 1
C    -1.3940633    0.0000000   -2.4541524
C    -0.6970468    1.2072378   -2.4546277
C     0.6970468    1.2072378   -2.4546277
C     1.3940633    0.0000000   -2.4541524
C     0.6970468   -1.2072378   -2.4546277
C    -0.6970468   -1.2072378   -2.4546277
H    -2.4753995    0.0000000   -2.4503221
H    -1.2382321    2.1435655   -2.4536764
H     1.2382321    2.1435655   -2.4536764
H     2.4753995    0.0000000   -2.4503221
H     1.2382321   -2.1435655   -2.4536764
H    -1.2382321   -2.1435655   -2.4536764
units angstrom
}
""", 0)

S22_21 = input.process_input("""
molecule dimer {
0 1
C     2.5118997    1.6250148    0.0000000
C     2.7130094    0.9578537   -1.2082918
C     3.1177821   -0.3767436   -1.2083647
C     3.3213848   -1.0437307    0.0000000
C     3.1177821   -0.3767436    1.2083647
C     2.7130094    0.9578537    1.2082918
H     2.2024038    2.6611358    0.0000000
H     2.5511760    1.4736908   -2.1445900
H     3.2702999   -0.8951406   -2.1448379
H     3.6368139   -2.0781521    0.0000000
H     3.2702999   -0.8951406    2.1448379
H     2.5511760    1.4736908    2.1445900
--
0 1
H     0.8065245   -0.4358866    0.0000000
N    -0.1442408   -0.7686927    0.0000000
C    -0.5161122   -2.0893220    0.0000000
C    -1.8898755   -2.1814495    0.0000000
C    -2.3932317   -0.8470830    0.0000000
C    -1.2640653    0.0195887    0.0000000
C    -1.3896004    1.4117668    0.0000000
C    -2.6726501    1.9366450    0.0000000
C    -3.8054511    1.0974790    0.0000000
C    -3.6798167   -0.2817209    0.0000000
H     0.2310024   -2.8653173    0.0000000
H    -2.4585759   -3.0956052    0.0000000
H    -0.5188733    2.0539520    0.0000000
H    -2.8077570    3.0097859    0.0000000
H    -4.7905991    1.5439372    0.0000000
H    -4.5580187   -0.9142916    0.0000000
units angstrom
}
""", 0)

S22_22 = input.process_input("""
molecule dimer {
0 1
C    -2.0071056    0.7638459   -0.1083509
O    -1.3885044    1.9298523   -0.4431206
H    -0.5238121    1.9646519   -0.0064609
C    -1.4630807   -0.1519120    0.7949930
C    -2.1475789   -1.3295094    1.0883677
C    -3.3743208   -1.6031427    0.4895864
C    -3.9143727   -0.6838545   -0.4091028
C    -3.2370496    0.4929609   -0.7096126
H    -0.5106510    0.0566569    1.2642563
H    -1.7151135   -2.0321452    1.7878417
H    -3.9024664   -2.5173865    0.7197947
H    -4.8670730   -0.8822939   -0.8811319
H    -3.6431662    1.2134345   -1.4057590
--
0 1
O     1.3531168    1.9382724    0.4723133
C     2.0369747    0.7865043    0.1495491
H     1.7842846    2.3487495    1.2297110
C     1.5904026    0.0696860   -0.9574153
C     2.2417367   -1.1069765   -1.3128110
C     3.3315674   -1.5665603   -0.5748636
C     3.7696838   -0.8396901    0.5286439
C     3.1224836    0.3383498    0.8960491
H     0.7445512    0.4367983   -1.5218583
H     1.8921463   -1.6649726   -2.1701843
H     3.8330227   -2.4811537   -0.8566666
H     4.6137632   -1.1850101    1.1092635
H     3.4598854    0.9030376    1.7569489
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:

    GEOS["%s-%s-dimer"      % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn))
    GEOS["%s-%s-monoA-CP"   % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn)) + monoA_CP
    GEOS["%s-%s-monoB-CP"   % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn)) + monoB_CP
    GEOS["%s-%s-monoA-unCP" % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn)) + monoA_unCP
    GEOS["%s-%s-monoB-unCP" % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn)) + monoB_unCP
