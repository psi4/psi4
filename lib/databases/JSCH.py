"""
**JSCH**

| Database (Hobza) of interaction energies for nucelobase pairs.
| Geometries and reference interaction energies from Jurecka et al. PCCP 8 1985 (2006).
| Corrections implemented from footnote 92 of Burns et al., JCP 134 084107 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``

----

"""
import input

# <<< JSCH Database Module >>>
dbse = 'JSCH'

# <<< Database Members >>>
HRXN = range(1, 125)
HRXN_SM = [9, 97]
HRXN_LG = [63]

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
BIND['%s-%s' % (dbse,   1)] = -32.06
BIND['%s-%s' % (dbse,   2)] = -31.59
BIND['%s-%s' % (dbse,   3)] = -16.86
BIND['%s-%s' % (dbse,   4)] = -18.16
BIND['%s-%s' % (dbse,   5)] = -33.30
BIND['%s-%s' % (dbse,   6)] = -24.90
BIND['%s-%s' % (dbse,   7)] = -19.10
BIND['%s-%s' % (dbse,   8)] = -51.40
BIND['%s-%s' % (dbse,   9)] = -10.30
BIND['%s-%s' % (dbse,  10)] = -13.70
BIND['%s-%s' % (dbse,  11)] = -29.50
BIND['%s-%s' % (dbse,  12)] = -14.20
BIND['%s-%s' % (dbse,  13)] = -19.50
BIND['%s-%s' % (dbse,  14)] = -19.70
BIND['%s-%s' % (dbse,  15)] =  -5.20
BIND['%s-%s' % (dbse,  16)] = -17.80
BIND['%s-%s' % (dbse,  17)] = -16.60
BIND['%s-%s' % (dbse,  18)] = -17.60
BIND['%s-%s' % (dbse,  19)] = -21.30
BIND['%s-%s' % (dbse,  20)] = -21.80
BIND['%s-%s' % (dbse,  21)] = -22.70
BIND['%s-%s' % (dbse,  22)] = -19.40
BIND['%s-%s' % (dbse,  23)] = -18.90
BIND['%s-%s' % (dbse,  24)] = -14.40
BIND['%s-%s' % (dbse,  25)] = -12.80
BIND['%s-%s' % (dbse,  26)] = -18.80
BIND['%s-%s' % (dbse,  27)] = -13.50
BIND['%s-%s' % (dbse,  28)] = -14.50
BIND['%s-%s' % (dbse,  29)] = -13.70
BIND['%s-%s' % (dbse,  30)] = -12.20
BIND['%s-%s' % (dbse,  31)] = -22.80
BIND['%s-%s' % (dbse,  32)] = -12.60
BIND['%s-%s' % (dbse,  33)] = -16.40
BIND['%s-%s' % (dbse,  34)] = -35.80
BIND['%s-%s' % (dbse,  35)] = -18.40
BIND['%s-%s' % (dbse,  36)] = -11.30
BIND['%s-%s' % (dbse,  37)] = -30.70
BIND['%s-%s' % (dbse,  38)] = -31.40
BIND['%s-%s' % (dbse,  39)] =  -3.68
BIND['%s-%s' % (dbse,  40)] =  -4.82
BIND['%s-%s' % (dbse,  41)] =  -2.34
BIND['%s-%s' % (dbse,  42)] =  -2.16
BIND['%s-%s' % (dbse,  43)] =   3.09
BIND['%s-%s' % (dbse,  44)] =   1.93
BIND['%s-%s' % (dbse,  45)] =  -3.91
BIND['%s-%s' % (dbse,  46)] =   1.24
BIND['%s-%s' % (dbse,  47)] =  -0.31
BIND['%s-%s' % (dbse,  48)] =   0.58
BIND['%s-%s' % (dbse,  49)] =  -0.47
BIND['%s-%s' % (dbse,  50)] =  -0.18
BIND['%s-%s' % (dbse,  51)] =  -4.22
BIND['%s-%s' % (dbse,  52)] =  -1.15
BIND['%s-%s' % (dbse,  53)] =   0.30
BIND['%s-%s' % (dbse,  54)] =  -4.06
BIND['%s-%s' % (dbse,  55)] =   0.88
BIND['%s-%s' % (dbse,  56)] =  -0.92
BIND['%s-%s' % (dbse,  57)] =  -1.55
BIND['%s-%s' % (dbse,  58)] =   0.70
BIND['%s-%s' % (dbse,  59)] =  -1.71
BIND['%s-%s' % (dbse,  60)] =  -1.30
BIND['%s-%s' % (dbse,  61)] =  -0.70
BIND['%s-%s' % (dbse,  62)] =   1.00
BIND['%s-%s' % (dbse,  63)] =  -4.50
BIND['%s-%s' % (dbse,  64)] =   1.40
BIND['%s-%s' % (dbse,  65)] =  -4.80
BIND['%s-%s' % (dbse,  66)] =  -0.10
BIND['%s-%s' % (dbse,  67)] =  -3.00
BIND['%s-%s' % (dbse,  68)] =  -5.20
BIND['%s-%s' % (dbse,  69)] =   0.80
BIND['%s-%s' % (dbse,  70)] =   3.10
BIND['%s-%s' % (dbse,  71)] = -19.02
BIND['%s-%s' % (dbse,  72)] = -20.35
BIND['%s-%s' % (dbse,  73)] = -12.30
BIND['%s-%s' % (dbse,  74)] = -14.57
BIND['%s-%s' % (dbse,  75)] =   2.45
BIND['%s-%s' % (dbse,  76)] =  -3.85
BIND['%s-%s' % (dbse,  77)] =  -8.88
BIND['%s-%s' % (dbse,  78)] =  -9.92
BIND['%s-%s' % (dbse,  79)] =   0.32
BIND['%s-%s' % (dbse,  80)] =   0.64
BIND['%s-%s' % (dbse,  81)] =  -0.98
BIND['%s-%s' % (dbse,  82)] =  -9.10
BIND['%s-%s' % (dbse,  83)] =  -9.11
BIND['%s-%s' % (dbse,  84)] =  -8.27
BIND['%s-%s' % (dbse,  85)] =  -9.43
BIND['%s-%s' % (dbse,  86)] =  -7.43
BIND['%s-%s' % (dbse,  87)] =  -8.80
BIND['%s-%s' % (dbse,  88)] =  -9.11
BIND['%s-%s' % (dbse,  89)] =  -8.58
BIND['%s-%s' % (dbse,  90)] = -12.67
BIND['%s-%s' % (dbse,  91)] = -10.22
BIND['%s-%s' % (dbse,  92)] = -11.38
BIND['%s-%s' % (dbse,  93)] = -10.02
BIND['%s-%s' % (dbse,  94)] =  -9.79
BIND['%s-%s' % (dbse,  95)] = -10.60
BIND['%s-%s' % (dbse,  96)] = -10.42
BIND['%s-%s' % (dbse,  97)] =  -7.46
BIND['%s-%s' % (dbse,  98)] = -12.09
BIND['%s-%s' % (dbse,  99)] =  -3.54
BIND['%s-%s' % (dbse, 100)] =  -1.62
BIND['%s-%s' % (dbse, 101)] =  -6.06
BIND['%s-%s' % (dbse, 102)] =  -4.18
BIND['%s-%s' % (dbse, 103)] = -10.80
BIND['%s-%s' % (dbse, 104)] =  -7.88
BIND['%s-%s' % (dbse, 105)] =  -9.14
BIND['%s-%s' % (dbse, 106)] =  -4.69
BIND['%s-%s' % (dbse, 107)] =  -7.58
BIND['%s-%s' % (dbse, 108)] =  -6.07
BIND['%s-%s' % (dbse, 109)] =  -5.67
BIND['%s-%s' % (dbse, 110)] =  -4.96
BIND['%s-%s' % (dbse, 111)] =  -4.96
BIND['%s-%s' % (dbse, 112)] =  -5.44
BIND['%s-%s' % (dbse, 113)] =  -6.64
BIND['%s-%s' % (dbse, 114)] =  -6.07
BIND['%s-%s' % (dbse, 115)] =  -6.25
BIND['%s-%s' % (dbse, 116)] =  -3.86
BIND['%s-%s' % (dbse, 117)] =  -8.10
BIND['%s-%s' % (dbse, 118)] =  -7.90
BIND['%s-%s' % (dbse, 119)] =  -6.70
BIND['%s-%s' % (dbse, 120)] =  -6.20
BIND['%s-%s' % (dbse, 121)] =  -7.70
BIND['%s-%s' % (dbse, 122)] =  -6.50
BIND['%s-%s' % (dbse, 123)] = -12.40
BIND['%s-%s' % (dbse, 124)] = -11.60

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse,   1)] = 'HB-01 G...C WC'
TAGL['%s-%s-dimer'      % (dbse,   1)] =       'G...C WC'
TAGL['%s-%s-monoA-CP'   % (dbse,   1)] =       'Cytosine from G...C WC'
TAGL['%s-%s-monoB-CP'   % (dbse,   1)] =       'Guanine from G...C WC'
TAGL['%s-%s-monoA-unCP' % (dbse,   1)] =       'Cytosine from G...C WC'
TAGL['%s-%s-monoB-unCP' % (dbse,   1)] =       'Guanine from G...C WC'
TAGL['%s-%s'            % (dbse,   2)] = 'HB-02 mG...mC WC'
TAGL['%s-%s-dimer'      % (dbse,   2)] =       'mG...mC WC'
TAGL['%s-%s-monoA-CP'   % (dbse,   2)] =       'methyl-Cytosine from mG...mC WC'
TAGL['%s-%s-monoB-CP'   % (dbse,   2)] =       'methyl-Guanine from mG...mC WC'
TAGL['%s-%s-monoA-unCP' % (dbse,   2)] =       'methyl-Cytosine from mG...mC WC'
TAGL['%s-%s-monoB-unCP' % (dbse,   2)] =       'methyl-Guanine from mG...mC WC'
TAGL['%s-%s'            % (dbse,   3)] = 'HB-03 A...T WC'
TAGL['%s-%s-dimer'      % (dbse,   3)] =       'A...T WC'
TAGL['%s-%s-monoA-CP'   % (dbse,   3)] =       'Adenine from A...T WC'
TAGL['%s-%s-monoB-CP'   % (dbse,   3)] =       'Thymine from A...T WC'
TAGL['%s-%s-monoA-unCP' % (dbse,   3)] =       'Adenine from A...T WC'
TAGL['%s-%s-monoB-unCP' % (dbse,   3)] =       'Thymine from A...T WC'
TAGL['%s-%s'            % (dbse,   4)] = 'HB-04 mA...mT H'
TAGL['%s-%s-dimer'      % (dbse,   4)] =       'mA...mT H'
TAGL['%s-%s-monoA-CP'   % (dbse,   4)] =       'methyl-Adenine from mA...mT H'
TAGL['%s-%s-monoB-CP'   % (dbse,   4)] =       'methyl-Thymine from mA...mT H'
TAGL['%s-%s-monoA-unCP' % (dbse,   4)] =       'methyl-Adenine from mA...mT H'
TAGL['%s-%s-monoB-unCP' % (dbse,   4)] =       'methyl-Thymine from mA...mT H'
TAGL['%s-%s'            % (dbse,   5)] = 'HB-05 8oG...C WC pl'
TAGL['%s-%s-dimer'      % (dbse,   5)] =       '8oG...C WC pl'
TAGL['%s-%s-monoA-CP'   % (dbse,   5)] =       '8-oxo-Guanine from 8oG...C WC pl'
TAGL['%s-%s-monoB-CP'   % (dbse,   5)] =       'Cytosine from 8oG...C WC pl'
TAGL['%s-%s-monoA-unCP' % (dbse,   5)] =       '8-oxo-Guanine from 8oG...C WC pl'
TAGL['%s-%s-monoB-unCP' % (dbse,   5)] =       'Cytosine from 8oG...C WC pl'
TAGL['%s-%s'            % (dbse,   6)] = 'HB-06 I...C WC pl'
TAGL['%s-%s-dimer'      % (dbse,   6)] =       'I...C WC pl'
TAGL['%s-%s-monoA-CP'   % (dbse,   6)] =       'Cytosine from I...C WC pl'
TAGL['%s-%s-monoB-CP'   % (dbse,   6)] =       'Inosine from I...C WC pl'
TAGL['%s-%s-monoA-unCP' % (dbse,   6)] =       'Cytosine from I...C WC pl'
TAGL['%s-%s-monoB-unCP' % (dbse,   6)] =       'Inosine from I...C WC pl'
TAGL['%s-%s'            % (dbse,   7)] = 'HB-07 G...U wobble'
TAGL['%s-%s-dimer'      % (dbse,   7)] =       'G...U wobble'
TAGL['%s-%s-monoA-CP'   % (dbse,   7)] =       'Guanine from G...U wobble'
TAGL['%s-%s-monoB-CP'   % (dbse,   7)] =       'Uracil from G...U wobble'
TAGL['%s-%s-monoA-unCP' % (dbse,   7)] =       'Guanine from G...U wobble'
TAGL['%s-%s-monoB-unCP' % (dbse,   7)] =       'Uracil from G...U wobble'
TAGL['%s-%s'            % (dbse,   8)] = 'HB-08 CCH+'
TAGL['%s-%s-dimer'      % (dbse,   8)] =       'CCH+'
TAGL['%s-%s-monoA-CP'   % (dbse,   8)] =       'Cytosine from CCH+'
TAGL['%s-%s-monoB-CP'   % (dbse,   8)] =       'protonated-Cytosine from CCH+'
TAGL['%s-%s-monoA-unCP' % (dbse,   8)] =       'Cytosine from CCH+'
TAGL['%s-%s-monoB-unCP' % (dbse,   8)] =       'protonated-Cytosine from CCH+'
TAGL['%s-%s'            % (dbse,   9)] = 'HB-09 U...U Calcutta pl'
TAGL['%s-%s-dimer'      % (dbse,   9)] =       'U...U Calcutta pl'
TAGL['%s-%s-monoA-CP'   % (dbse,   9)] =       'Uracil from U...U Calcutta pl'
TAGL['%s-%s-monoB-CP'   % (dbse,   9)] =       'Uracil from U...U Calcutta pl'
TAGL['%s-%s-monoA-unCP' % (dbse,   9)] =       'Uracil from U...U Calcutta pl'
TAGL['%s-%s-monoB-unCP' % (dbse,   9)] =       'Uracil from U...U Calcutta pl'
TAGL['%s-%s'            % (dbse,  10)] = 'HB-10 U...U pl'
TAGL['%s-%s-dimer'      % (dbse,  10)] =       'U...U pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  10)] =       'Uracil from U...U pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  10)] =       'Uracil from U...U pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  10)] =       'Uracil from U...U pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  10)] =       'Uracil from U...U pl'
TAGL['%s-%s'            % (dbse,  11)] = 'HB-11 6tG...C WC pl'
TAGL['%s-%s-dimer'      % (dbse,  11)] =       '6tG...C WC pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  11)] =       'Cytosine from 6tG...C WC pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  11)] =       '6-thio-Guanine from 6tG...C WC pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  11)] =       'Cytosine from 6tG...C WC pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  11)] =       '6-thio-Guanine from 6tG...C WC pl'
TAGL['%s-%s'            % (dbse,  12)] = 'HB-12 A...4tU WC'
TAGL['%s-%s-dimer'      % (dbse,  12)] =       'A...4tU WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  12)] =       'Adenine from A...4tU WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  12)] =       '4-thio-Uracil from A...4tU WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  12)] =       'Adenine from A...4tU WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  12)] =       '4-thio-Uracil from A...4tU WC'
TAGL['%s-%s'            % (dbse,  13)] = 'HB-13 2-aminoA...T'
TAGL['%s-%s-dimer'      % (dbse,  13)] =       '2-aminoA...T'
TAGL['%s-%s-monoA-CP'   % (dbse,  13)] =       '2-amino-Adenine from 2-aminoA...T'
TAGL['%s-%s-monoB-CP'   % (dbse,  13)] =       'Thymine from 2-aminoA...T'
TAGL['%s-%s-monoA-unCP' % (dbse,  13)] =       '2-amino-Adenine from 2-aminoA...T'
TAGL['%s-%s-monoB-unCP' % (dbse,  13)] =       'Thymine from 2-aminoA...T'
TAGL['%s-%s'            % (dbse,  14)] = 'HB-14 2-aminoA...T pl'
TAGL['%s-%s-dimer'      % (dbse,  14)] =       '2-aminoA...T pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  14)] =       '2-amino-Adenine from 2-aminoA...T pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  14)] =       'Thymine from 2-aminoA...T pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  14)] =       '2-amino-Adenine from 2-aminoA...T pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  14)] =       'Thymine from 2-aminoA...T pl'
TAGL['%s-%s'            % (dbse,  15)] = 'HB-15 A...F'
TAGL['%s-%s-dimer'      % (dbse,  15)] =       'A...F'
TAGL['%s-%s-monoA-CP'   % (dbse,  15)] =       'Adenine from A...F'
TAGL['%s-%s-monoB-CP'   % (dbse,  15)] =       'difluorotoluene from A...F'
TAGL['%s-%s-monoA-unCP' % (dbse,  15)] =       'Adenine from A...F'
TAGL['%s-%s-monoB-unCP' % (dbse,  15)] =       'difluorotoluene from A...F'
TAGL['%s-%s'            % (dbse,  16)] = 'HB-16 G...4tU'
TAGL['%s-%s-dimer'      % (dbse,  16)] =       'G...4tU'
TAGL['%s-%s-monoA-CP'   % (dbse,  16)] =       'Guanine from G...4tU'
TAGL['%s-%s-monoB-CP'   % (dbse,  16)] =       '4-thio-Uracil from G...4tU'
TAGL['%s-%s-monoA-unCP' % (dbse,  16)] =       'Guanine from G...4tU'
TAGL['%s-%s-monoB-unCP' % (dbse,  16)] =       '4-thio-Uracil from G...4tU'
TAGL['%s-%s'            % (dbse,  17)] = 'HB-17 G...2tU'
TAGL['%s-%s-dimer'      % (dbse,  17)] =       'G...2tU'
TAGL['%s-%s-monoA-CP'   % (dbse,  17)] =       'Guanine from G...2tU'
TAGL['%s-%s-monoB-CP'   % (dbse,  17)] =       '2-thio-Uracil from G...2tU'
TAGL['%s-%s-monoA-unCP' % (dbse,  17)] =       'Guanine from G...2tU'
TAGL['%s-%s-monoB-unCP' % (dbse,  17)] =       '2-thio-Uracil from G...2tU'
TAGL['%s-%s'            % (dbse,  18)] = 'HB-18 A...C pl'
TAGL['%s-%s-dimer'      % (dbse,  18)] =       'A...C pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  18)] =       'Cytosine from A...C pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  18)] =       'Adenine from A...C pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  18)] =       'Cytosine from A...C pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  18)] =       'Adenine from A...C pl'
TAGL['%s-%s'            % (dbse,  19)] = 'HB-19 G...G pl'
TAGL['%s-%s-dimer'      % (dbse,  19)] =       'G...G pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  19)] =       'Guanine from G...G pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  19)] =       'Guanine from G...G pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  19)] =       'Guanine from G...G pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  19)] =       'Guanine from G...G pl'
TAGL['%s-%s'            % (dbse,  20)] = 'HB-20 G...6tG pl'
TAGL['%s-%s-dimer'      % (dbse,  20)] =       'G...6tG pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  20)] =       'Guanine from G...6tG pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  20)] =       '6-thio-Guanine from G...6tG pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  20)] =       'Guanine from G...6tG pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  20)] =       '6-thio-Guanine from G...6tG pl'
TAGL['%s-%s'            % (dbse,  21)] = 'HB-21 6tG...G pl'
TAGL['%s-%s-dimer'      % (dbse,  21)] =       '6tG...G pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  21)] =       '6-thio-Guanine from 6tG...G pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  21)] =       'Guanine from 6tG...G pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  21)] =       '6-thio-Guanine from 6tG...G pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  21)] =       'Guanine from 6tG...G pl'
TAGL['%s-%s'            % (dbse,  22)] = 'HB-22 G...A 1'
TAGL['%s-%s-dimer'      % (dbse,  22)] =       'G...A 1'
TAGL['%s-%s-monoA-CP'   % (dbse,  22)] =       'Guanine from G...A 1'
TAGL['%s-%s-monoB-CP'   % (dbse,  22)] =       'Adenine from G...A 1'
TAGL['%s-%s-monoA-unCP' % (dbse,  22)] =       'Guanine from G...A 1'
TAGL['%s-%s-monoB-unCP' % (dbse,  22)] =       'Adenine from G...A 1'
TAGL['%s-%s'            % (dbse,  23)] = 'HB-23 G...A 1 pl'
TAGL['%s-%s-dimer'      % (dbse,  23)] =       'G...A 1 pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  23)] =       'Adenine from G...A 1 pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  23)] =       'Guanine from G...A 1 pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  23)] =       'Adenine from G...A 1 pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  23)] =       'Guanine from G...A 1 pl'
TAGL['%s-%s'            % (dbse,  24)] = 'HB-24 G...A 2'
TAGL['%s-%s-dimer'      % (dbse,  24)] =       'G...A 2'
TAGL['%s-%s-monoA-CP'   % (dbse,  24)] =       'Guanine from G...A 2'
TAGL['%s-%s-monoB-CP'   % (dbse,  24)] =       'Adenine from G...A 2'
TAGL['%s-%s-monoA-unCP' % (dbse,  24)] =       'Guanine from G...A 2'
TAGL['%s-%s-monoB-unCP' % (dbse,  24)] =       'Adenine from G...A 2'
TAGL['%s-%s'            % (dbse,  25)] = 'HB-25 G...A 2 pl'
TAGL['%s-%s-dimer'      % (dbse,  25)] =       'G...A 2 pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  25)] =       'Guanine from G...A 2 pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  25)] =       'Adenine from G...A 2 pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  25)] =       'Guanine from G...A 2 pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  25)] =       'Adenine from G...A 2 pl'
TAGL['%s-%s'            % (dbse,  26)] = 'HB-26 G...A 3'
TAGL['%s-%s-dimer'      % (dbse,  26)] =       'G...A 3'
TAGL['%s-%s-monoA-CP'   % (dbse,  26)] =       'Guanine from G...A 3'
TAGL['%s-%s-monoB-CP'   % (dbse,  26)] =       'Adenine from G...A 3'
TAGL['%s-%s-monoA-unCP' % (dbse,  26)] =       'Guanine from G...A 3'
TAGL['%s-%s-monoB-unCP' % (dbse,  26)] =       'Adenine from G...A 3'
TAGL['%s-%s'            % (dbse,  27)] = 'HB-27 G...A 4'
TAGL['%s-%s-dimer'      % (dbse,  27)] =       'G...A 4'
TAGL['%s-%s-monoA-CP'   % (dbse,  27)] =       'Guanine from G...A 4'
TAGL['%s-%s-monoB-CP'   % (dbse,  27)] =       'Adenine from G...A 4'
TAGL['%s-%s-monoA-unCP' % (dbse,  27)] =       'Guanine from G...A 4'
TAGL['%s-%s-monoB-unCP' % (dbse,  27)] =       'Adenine from G...A 4'
TAGL['%s-%s'            % (dbse,  28)] = 'HB-28 A...A 1 pl'
TAGL['%s-%s-dimer'      % (dbse,  28)] =       'A...A 1 pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  28)] =       'Adenine from A...A 1 pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  28)] =       'Adenine from A...A 1 pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  28)] =       'Adenine from A...A 1 pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  28)] =       'Adenine from A...A 1 pl'
TAGL['%s-%s'            % (dbse,  29)] = 'HB-29 A...A 2 pl'
TAGL['%s-%s-dimer'      % (dbse,  29)] =       'A...A 2 pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  29)] =       'Adenine from A...A 2 pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  29)] =       'Adenine from A...A 2 pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  29)] =       'Adenine from A...A 2 pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  29)] =       'Adenine from A...A 2 pl'
TAGL['%s-%s'            % (dbse,  30)] = 'HB-30 A...A 3 pl'
TAGL['%s-%s-dimer'      % (dbse,  30)] =       'A...A 3 pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  30)] =       'Adenine from A...A 3 pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  30)] =       'Adenine from A...A 3 pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  30)] =       'Adenine from A...A 3 pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  30)] =       'Adenine from A...A 3 pl'
TAGL['%s-%s'            % (dbse,  31)] = 'HB-31 8oG...G'
TAGL['%s-%s-dimer'      % (dbse,  31)] =       '8oG...G'
TAGL['%s-%s-monoA-CP'   % (dbse,  31)] =       'Guanine from 8oG...G'
TAGL['%s-%s-monoB-CP'   % (dbse,  31)] =       '8-oxo-Guanine from 8oG...G'
TAGL['%s-%s-monoA-unCP' % (dbse,  31)] =       'Guanine from 8oG...G'
TAGL['%s-%s-monoB-unCP' % (dbse,  31)] =       '8-oxo-Guanine from 8oG...G'
TAGL['%s-%s'            % (dbse,  32)] = 'HB-32 2tU....2tU pl'
TAGL['%s-%s-dimer'      % (dbse,  32)] =       '2tU....2tU pl'
TAGL['%s-%s-monoA-CP'   % (dbse,  32)] =       '2-thio-Uracil from 2tU....2tU pl'
TAGL['%s-%s-monoB-CP'   % (dbse,  32)] =       '2-thio-Uracil from 2tU....2tU pl'
TAGL['%s-%s-monoA-unCP' % (dbse,  32)] =       '2-thio-Uracil from 2tU....2tU pl'
TAGL['%s-%s-monoB-unCP' % (dbse,  32)] =       '2-thio-Uracil from 2tU....2tU pl'
TAGL['%s-%s'            % (dbse,  33)] = 'HB-33 A...T WC'
TAGL['%s-%s-dimer'      % (dbse,  33)] =       'A...T WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  33)] =       'methyl-Adenine from A...T WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  33)] =       'methyl-Thymine from A...T WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  33)] =       'methyl-Adenine from A...T WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  33)] =       'methyl-Thymine from A...T WC'
TAGL['%s-%s'            % (dbse,  34)] = 'HB-34 G...C WC'
TAGL['%s-%s-dimer'      % (dbse,  34)] =       'G...C WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  34)] =       'methyl-Cytosine from G...C WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  34)] =       'methyl-Guanine from G...C WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  34)] =       'methyl-Cytosine from G...C WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  34)] =       'methyl-Guanine from G...C WC'
TAGL['%s-%s'            % (dbse,  35)] = 'HB-35 A...T WC'
TAGL['%s-%s-dimer'      % (dbse,  35)] =       'A...T WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  35)] =       'methyl-Adenine from A...T WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  35)] =       'methyl-Thymine from A...T WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  35)] =       'methyl-Adenine from A...T WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  35)] =       'methyl-Thymine from A...T WC'
TAGL['%s-%s'            % (dbse,  36)] = 'HB-36 G...A HB'
TAGL['%s-%s-dimer'      % (dbse,  36)] =       'G...A HB'
TAGL['%s-%s-monoA-CP'   % (dbse,  36)] =       'Guanine from G...A HB'
TAGL['%s-%s-monoB-CP'   % (dbse,  36)] =       'Adenine from G...A HB'
TAGL['%s-%s-monoA-unCP' % (dbse,  36)] =       'Guanine from G...A HB'
TAGL['%s-%s-monoB-unCP' % (dbse,  36)] =       'Adenine from G...A HB'
TAGL['%s-%s'            % (dbse,  37)] = 'HB-37 C...G WC'
TAGL['%s-%s-dimer'      % (dbse,  37)] =       'C...G WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  37)] =       'Cytosine from C...G WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  37)] =       'Guanine from C...G WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  37)] =       'Cytosine from C...G WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  37)] =       'Guanine from C...G WC'
TAGL['%s-%s'            % (dbse,  38)] = 'HB-38 G...C WC'
TAGL['%s-%s-dimer'      % (dbse,  38)] =       'G...C WC'
TAGL['%s-%s-monoA-CP'   % (dbse,  38)] =       'Guanine from G...C WC'
TAGL['%s-%s-monoB-CP'   % (dbse,  38)] =       'Cytosine from G...C WC'
TAGL['%s-%s-monoA-unCP' % (dbse,  38)] =       'Guanine from G...C WC'
TAGL['%s-%s-monoB-unCP' % (dbse,  38)] =       'Cytosine from G...C WC'
TAGL['%s-%s'            % (dbse,  39)] = 'IS-01 GG0/3.36 CGis036'
TAGL['%s-%s-dimer'      % (dbse,  39)] =       'GG0/3.36 CGis036'
TAGL['%s-%s-monoA-CP'   % (dbse,  39)] =       'Guanine from GG0/3.36 CGis036'
TAGL['%s-%s-monoB-CP'   % (dbse,  39)] =       'Cytosine from GG0/3.36 CGis036'
TAGL['%s-%s-monoA-unCP' % (dbse,  39)] =       'Guanine from GG0/3.36 CGis036'
TAGL['%s-%s-monoB-unCP' % (dbse,  39)] =       'Cytosine from GG0/3.36 CGis036'
TAGL['%s-%s'            % (dbse,  40)] = 'IS-02 GG0/3.36 GCis036'
TAGL['%s-%s-dimer'      % (dbse,  40)] =       'GG0/3.36 GCis036'
TAGL['%s-%s-monoA-CP'   % (dbse,  40)] =       'Cytosine from GG0/3.36 GCis036'
TAGL['%s-%s-monoB-CP'   % (dbse,  40)] =       'Guanine from GG0/3.36 GCis036'
TAGL['%s-%s-monoA-unCP' % (dbse,  40)] =       'Cytosine from GG0/3.36 GCis036'
TAGL['%s-%s-monoB-unCP' % (dbse,  40)] =       'Guanine from GG0/3.36 GCis036'
TAGL['%s-%s'            % (dbse,  41)] = 'IS-03 AA20/3.05 ATis2005'
TAGL['%s-%s-dimer'      % (dbse,  41)] =       'AA20/3.05 ATis2005'
TAGL['%s-%s-monoA-CP'   % (dbse,  41)] =       'Adenine from AA20/3.05 ATis2005'
TAGL['%s-%s-monoB-CP'   % (dbse,  41)] =       'Thymine from AA20/3.05 ATis2005'
TAGL['%s-%s-monoA-unCP' % (dbse,  41)] =       'Adenine from AA20/3.05 ATis2005'
TAGL['%s-%s-monoB-unCP' % (dbse,  41)] =       'Thymine from AA20/3.05 ATis2005'
TAGL['%s-%s'            % (dbse,  42)] = 'IS-04 AA20/3.05 TAis2005'
TAGL['%s-%s-dimer'      % (dbse,  42)] =       'AA20/3.05 TAis2005'
TAGL['%s-%s-monoA-CP'   % (dbse,  42)] =       'Thymine from AA20/3.05 TAis2005'
TAGL['%s-%s-monoB-CP'   % (dbse,  42)] =       'Adenine from AA20/3.05 TAis2005'
TAGL['%s-%s-monoA-unCP' % (dbse,  42)] =       'Thymine from AA20/3.05 TAis2005'
TAGL['%s-%s-monoB-unCP' % (dbse,  42)] =       'Adenine from AA20/3.05 TAis2005'
TAGL['%s-%s'            % (dbse,  43)] = 'IS-05 GC0/3.25 C//Cis'
TAGL['%s-%s-dimer'      % (dbse,  43)] =       'GC0/3.25 C//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  43)] =       'Cytosine from GC0/3.25 C//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  43)] =       'Cytosine from GC0/3.25 C//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  43)] =       'Cytosine from GC0/3.25 C//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  43)] =       'Cytosine from GC0/3.25 C//Cis'
TAGL['%s-%s'            % (dbse,  44)] = 'IS-06 GC0/3.25 G//Gis'
TAGL['%s-%s-dimer'      % (dbse,  44)] =       'GC0/3.25 G//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  44)] =       'Guanine from GC0/3.25 G//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  44)] =       'Guanine from GC0/3.25 G//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  44)] =       'Guanine from GC0/3.25 G//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  44)] =       'Guanine from GC0/3.25 G//Gis'
TAGL['%s-%s'            % (dbse,  45)] = 'IS-07 CG0/3.19 G//Gis'
TAGL['%s-%s-dimer'      % (dbse,  45)] =       'CG0/3.19 G//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  45)] =       'Guanine from CG0/3.19 G//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  45)] =       'Guanine from CG0/3.19 G//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  45)] =       'Guanine from CG0/3.19 G//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  45)] =       'Guanine from CG0/3.19 G//Gis'
TAGL['%s-%s'            % (dbse,  46)] = 'IS-08 CG0/3.19 C//Cis'
TAGL['%s-%s-dimer'      % (dbse,  46)] =       'CG0/3.19 C//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  46)] =       'Cytosine from CG0/3.19 C//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  46)] =       'Cytosine from CG0/3.19 C//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  46)] =       'Cytosine from CG0/3.19 C//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  46)] =       'Cytosine from CG0/3.19 C//Cis'
TAGL['%s-%s'            % (dbse,  47)] = 'IS-09 GA10/3.15 A//Cis'
TAGL['%s-%s-dimer'      % (dbse,  47)] =       'GA10/3.15 A//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  47)] =       'Adenine from GA10/3.15 A//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  47)] =       'Cytosine from GA10/3.15 A//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  47)] =       'Adenine from GA10/3.15 A//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  47)] =       'Cytosine from GA10/3.15 A//Cis'
TAGL['%s-%s'            % (dbse,  48)] = 'IS-10 GA10/3.15 T//Gis'
TAGL['%s-%s-dimer'      % (dbse,  48)] =       'GA10/3.15 T//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  48)] =       'Thymine from GA10/3.15 T//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  48)] =       'Guanine from GA10/3.15 T//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  48)] =       'Thymine from GA10/3.15 T//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  48)] =       'Guanine from GA10/3.15 T//Gis'
TAGL['%s-%s'            % (dbse,  49)] = 'IS-11 AG08/3.19 T//Gis'
TAGL['%s-%s-dimer'      % (dbse,  49)] =       'AG08/3.19 T//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  49)] =       'Guanine from AG08/3.19 T//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  49)] =       'Thymine from AG08/3.19 T//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  49)] =       'Guanine from AG08/3.19 T//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  49)] =       'Thymine from AG08/3.19 T//Gis'
TAGL['%s-%s'            % (dbse,  50)] = 'IS-12 AG08/3.19 A//Cis'
TAGL['%s-%s-dimer'      % (dbse,  50)] =       'AG08/3.19 A//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  50)] =       'Adenine from AG08/3.19 A//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  50)] =       'Cytosine from AG08/3.19 A//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  50)] =       'Adenine from AG08/3.19 A//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  50)] =       'Cytosine from AG08/3.19 A//Cis'
TAGL['%s-%s'            % (dbse,  51)] = 'IS-13 TG03.19 A//Gis'
TAGL['%s-%s-dimer'      % (dbse,  51)] =       'TG03.19 A//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  51)] =       'Adenine from TG03.19 A//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  51)] =       'Guanine from TG03.19 A//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  51)] =       'Adenine from TG03.19 A//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  51)] =       'Guanine from TG03.19 A//Gis'
TAGL['%s-%s'            % (dbse,  52)] = 'IS-14 TG03.19 T//Cis'
TAGL['%s-%s-dimer'      % (dbse,  52)] =       'TG03.19 T//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  52)] =       'Thymine from TG03.19 T//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  52)] =       'Cytosine from TG03.19 T//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  52)] =       'Thymine from TG03.19 T//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  52)] =       'Cytosine from TG03.19 T//Cis'
TAGL['%s-%s'            % (dbse,  53)] = 'IS-15 GT10/3.15 T//Cis'
TAGL['%s-%s-dimer'      % (dbse,  53)] =       'GT10/3.15 T//Cis'
TAGL['%s-%s-monoA-CP'   % (dbse,  53)] =       'Thymine from GT10/3.15 T//Cis'
TAGL['%s-%s-monoB-CP'   % (dbse,  53)] =       'Cytosine from GT10/3.15 T//Cis'
TAGL['%s-%s-monoA-unCP' % (dbse,  53)] =       'Thymine from GT10/3.15 T//Cis'
TAGL['%s-%s-monoB-unCP' % (dbse,  53)] =       'Cytosine from GT10/3.15 T//Cis'
TAGL['%s-%s'            % (dbse,  54)] = 'IS-16 GT10/3.15 A//Gis'
TAGL['%s-%s-dimer'      % (dbse,  54)] =       'GT10/3.15 A//Gis'
TAGL['%s-%s-monoA-CP'   % (dbse,  54)] =       'Adenine from GT10/3.15 A//Gis'
TAGL['%s-%s-monoB-CP'   % (dbse,  54)] =       'Guanine from GT10/3.15 A//Gis'
TAGL['%s-%s-monoA-unCP' % (dbse,  54)] =       'Adenine from GT10/3.15 A//Gis'
TAGL['%s-%s-monoB-unCP' % (dbse,  54)] =       'Guanine from GT10/3.15 A//Gis'
TAGL['%s-%s'            % (dbse,  55)] = 'IS-17 AT10/3.26 T//Tis'
TAGL['%s-%s-dimer'      % (dbse,  55)] =       'AT10/3.26 T//Tis'
TAGL['%s-%s-monoA-CP'   % (dbse,  55)] =       'Thymine from AT10/3.26 T//Tis'
TAGL['%s-%s-monoB-CP'   % (dbse,  55)] =       'Thymine from AT10/3.26 T//Tis'
TAGL['%s-%s-monoA-unCP' % (dbse,  55)] =       'Thymine from AT10/3.26 T//Tis'
TAGL['%s-%s-monoB-unCP' % (dbse,  55)] =       'Thymine from AT10/3.26 T//Tis'
TAGL['%s-%s'            % (dbse,  56)] = 'IS-18 AT10/3.26 A//Ais'
TAGL['%s-%s-dimer'      % (dbse,  56)] =       'AT10/3.26 A//Ais'
TAGL['%s-%s-monoA-CP'   % (dbse,  56)] =       'Adenine from AT10/3.26 A//Ais'
TAGL['%s-%s-monoB-CP'   % (dbse,  56)] =       'Adenine from AT10/3.26 A//Ais'
TAGL['%s-%s-monoA-unCP' % (dbse,  56)] =       'Adenine from AT10/3.26 A//Ais'
TAGL['%s-%s-monoB-unCP' % (dbse,  56)] =       'Adenine from AT10/3.26 A//Ais'
TAGL['%s-%s'            % (dbse,  57)] = 'IS-19 TA08/3.16 A//Ais'
TAGL['%s-%s-dimer'      % (dbse,  57)] =       'TA08/3.16 A//Ais'
TAGL['%s-%s-monoA-CP'   % (dbse,  57)] =       'Adenine from TA08/3.16 A//Ais'
TAGL['%s-%s-monoB-CP'   % (dbse,  57)] =       'Adenine from TA08/3.16 A//Ais'
TAGL['%s-%s-monoA-unCP' % (dbse,  57)] =       'Adenine from TA08/3.16 A//Ais'
TAGL['%s-%s-monoB-unCP' % (dbse,  57)] =       'Adenine from TA08/3.16 A//Ais'
TAGL['%s-%s'            % (dbse,  58)] = 'IS-20 TA08/3.16 T//Tis'
TAGL['%s-%s-dimer'      % (dbse,  58)] =       'TA08/3.16 T//Tis'
TAGL['%s-%s-monoA-CP'   % (dbse,  58)] =       'Thymine from TA08/3.16 T//Tis'
TAGL['%s-%s-monoB-CP'   % (dbse,  58)] =       'Thymine from TA08/3.16 T//Tis'
TAGL['%s-%s-monoA-unCP' % (dbse,  58)] =       'Thymine from TA08/3.16 T//Tis'
TAGL['%s-%s-monoB-unCP' % (dbse,  58)] =       'Thymine from TA08/3.16 T//Tis'
TAGL['%s-%s'            % (dbse,  59)] = 'IS-21 AA0/3.24 A//Tis'
TAGL['%s-%s-dimer'      % (dbse,  59)] =       'AA0/3.24 A//Tis'
TAGL['%s-%s-monoA-CP'   % (dbse,  59)] =       'Adenine from AA0/3.24 A//Tis'
TAGL['%s-%s-monoB-CP'   % (dbse,  59)] =       'Thymine from AA0/3.24 A//Tis'
TAGL['%s-%s-monoA-unCP' % (dbse,  59)] =       'Adenine from AA0/3.24 A//Tis'
TAGL['%s-%s-monoB-unCP' % (dbse,  59)] =       'Thymine from AA0/3.24 A//Tis'
TAGL['%s-%s'            % (dbse,  60)] = 'IS-22 AA0/3.24 T//Ais'
TAGL['%s-%s-dimer'      % (dbse,  60)] =       'AA0/3.24 T//Ais'
TAGL['%s-%s-monoA-CP'   % (dbse,  60)] =       'Adenine from AA0/3.24 T//Ais'
TAGL['%s-%s-monoB-CP'   % (dbse,  60)] =       'Thymine from AA0/3.24 T//Ais'
TAGL['%s-%s-monoA-unCP' % (dbse,  60)] =       'Adenine from AA0/3.24 T//Ais'
TAGL['%s-%s-monoB-unCP' % (dbse,  60)] =       'Thymine from AA0/3.24 T//Ais'
TAGL['%s-%s'            % (dbse,  61)] = 'IS-23 A...A IS'
TAGL['%s-%s-dimer'      % (dbse,  61)] =       'A...A IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  61)] =       'methyl-Adenine from A...A IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  61)] =       'methyl-Adenine from A...A IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  61)] =       'methyl-Adenine from A...A IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  61)] =       'methyl-Adenine from A...A IS'
TAGL['%s-%s'            % (dbse,  62)] = 'IS-24 T...T IS'
TAGL['%s-%s-dimer'      % (dbse,  62)] =       'T...T IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  62)] =       'methyl-Thymine from T...T IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  62)] =       'methyl-Thymine from T...T IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  62)] =       'methyl-Thymine from T...T IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  62)] =       'methyl-Thymine from T...T IS'
TAGL['%s-%s'            % (dbse,  63)] = 'IS-25 G...G IS'
TAGL['%s-%s-dimer'      % (dbse,  63)] =       'G...G IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  63)] =       'methyl-Guanine from G...G IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  63)] =       'methyl-Guanine from G...G IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  63)] =       'methyl-Guanine from G...G IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  63)] =       'methyl-Guanine from G...G IS'
TAGL['%s-%s'            % (dbse,  64)] = 'IS-26 C...C IS'
TAGL['%s-%s-dimer'      % (dbse,  64)] =       'C...C IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  64)] =       'methyl-Cytosine from C...C IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  64)] =       'methyl-Cytosine from C...C IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  64)] =       'methyl-Cytosine from C...C IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  64)] =       'methyl-Cytosine from C...C IS'
TAGL['%s-%s'            % (dbse,  65)] = 'IS-27 A...G IS'
TAGL['%s-%s-dimer'      % (dbse,  65)] =       'A...G IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  65)] =       'methyl-Adenine from A...G IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  65)] =       'methyl-Guanine from A...G IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  65)] =       'methyl-Adenine from A...G IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  65)] =       'methyl-Guanine from A...G IS'
TAGL['%s-%s'            % (dbse,  66)] = 'IS-28 T...C IS'
TAGL['%s-%s-dimer'      % (dbse,  66)] =       'T...C IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  66)] =       'methyl-Cytosine from T...C IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  66)] =       'methyl-Thymine from T...C IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  66)] =       'methyl-Cytosine from T...C IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  66)] =       'methyl-Thymine from T...C IS'
TAGL['%s-%s'            % (dbse,  67)] = 'IS-29 C...A IS'
TAGL['%s-%s-dimer'      % (dbse,  67)] =       'C...A IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  67)] =       'Cytosine from C...A IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  67)] =       'Adenine from C...A IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  67)] =       'Cytosine from C...A IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  67)] =       'Adenine from C...A IS'
TAGL['%s-%s'            % (dbse,  68)] = 'IS-30 G...G IS'
TAGL['%s-%s-dimer'      % (dbse,  68)] =       'G...G IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  68)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  68)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  68)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  68)] =       'Guanine from G...G IS'
TAGL['%s-%s'            % (dbse,  69)] = 'IS-31 G...G IS'
TAGL['%s-%s-dimer'      % (dbse,  69)] =       'G...G IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  69)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  69)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  69)] =       'Guanine from G...G IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  69)] =       'Guanine from G...G IS'
TAGL['%s-%s'            % (dbse,  70)] = 'IS-32 C...C IS'
TAGL['%s-%s-dimer'      % (dbse,  70)] =       'C...C IS'
TAGL['%s-%s-monoA-CP'   % (dbse,  70)] =       'Cytosine from C...C IS'
TAGL['%s-%s-monoB-CP'   % (dbse,  70)] =       'Cytosine from C...C IS'
TAGL['%s-%s-monoA-unCP' % (dbse,  70)] =       'Cytosine from C...C IS'
TAGL['%s-%s-monoB-unCP' % (dbse,  70)] =       'Cytosine from C...C IS'
TAGL['%s-%s'            % (dbse,  71)] = 'ST-01 G...C S'
TAGL['%s-%s-dimer'      % (dbse,  71)] =       'G...C S'
TAGL['%s-%s-monoA-CP'   % (dbse,  71)] =       'Guanine from G...C S'
TAGL['%s-%s-monoB-CP'   % (dbse,  71)] =       'Cytosine from G...C S'
TAGL['%s-%s-monoA-unCP' % (dbse,  71)] =       'Guanine from G...C S'
TAGL['%s-%s-monoB-unCP' % (dbse,  71)] =       'Cytosine from G...C S'
TAGL['%s-%s'            % (dbse,  72)] = 'ST-02 mG...mC S'
TAGL['%s-%s-dimer'      % (dbse,  72)] =       'mG...mC S'
TAGL['%s-%s-monoA-CP'   % (dbse,  72)] =       'methyl-Guanine from mG...mC S'
TAGL['%s-%s-monoB-CP'   % (dbse,  72)] =       'methyl-Cytosine from mG...mC S'
TAGL['%s-%s-monoA-unCP' % (dbse,  72)] =       'methyl-Guanine from mG...mC S'
TAGL['%s-%s-monoB-unCP' % (dbse,  72)] =       'methyl-Cytosine from mG...mC S'
TAGL['%s-%s'            % (dbse,  73)] = 'ST-03 A...T S'
TAGL['%s-%s-dimer'      % (dbse,  73)] =       'A...T S'
TAGL['%s-%s-monoA-CP'   % (dbse,  73)] =       'Adenine from A...T S'
TAGL['%s-%s-monoB-CP'   % (dbse,  73)] =       'Thymine from A...T S'
TAGL['%s-%s-monoA-unCP' % (dbse,  73)] =       'Adenine from A...T S'
TAGL['%s-%s-monoB-unCP' % (dbse,  73)] =       'Thymine from A...T S'
TAGL['%s-%s'            % (dbse,  74)] = 'ST-04 mA...mT S'
TAGL['%s-%s-dimer'      % (dbse,  74)] =       'mA...mT S'
TAGL['%s-%s-monoA-CP'   % (dbse,  74)] =       'methyl-Adenine from mA...mT S'
TAGL['%s-%s-monoB-CP'   % (dbse,  74)] =       'methyl-Thymine from mA...mT S'
TAGL['%s-%s-monoA-unCP' % (dbse,  74)] =       'methyl-Adenine from mA...mT S'
TAGL['%s-%s-monoB-unCP' % (dbse,  74)] =       'methyl-Thymine from mA...mT S'
TAGL['%s-%s'            % (dbse,  75)] = 'ST-05 CC1'
TAGL['%s-%s-dimer'      % (dbse,  75)] =       'CC1'
TAGL['%s-%s-monoA-CP'   % (dbse,  75)] =       'Cytosine from CC1'
TAGL['%s-%s-monoB-CP'   % (dbse,  75)] =       'Cytosine from CC1'
TAGL['%s-%s-monoA-unCP' % (dbse,  75)] =       'Cytosine from CC1'
TAGL['%s-%s-monoB-unCP' % (dbse,  75)] =       'Cytosine from CC1'
TAGL['%s-%s'            % (dbse,  76)] = 'ST-06 CC2'
TAGL['%s-%s-dimer'      % (dbse,  76)] =       'CC2'
TAGL['%s-%s-monoA-CP'   % (dbse,  76)] =       'Cytosine from CC2'
TAGL['%s-%s-monoB-CP'   % (dbse,  76)] =       'Cytosine from CC2'
TAGL['%s-%s-monoA-unCP' % (dbse,  76)] =       'Cytosine from CC2'
TAGL['%s-%s-monoB-unCP' % (dbse,  76)] =       'Cytosine from CC2'
TAGL['%s-%s'            % (dbse,  77)] = 'ST-07 CC3'
TAGL['%s-%s-dimer'      % (dbse,  77)] =       'CC3'
TAGL['%s-%s-monoA-CP'   % (dbse,  77)] =       'Cytosine from CC3'
TAGL['%s-%s-monoB-CP'   % (dbse,  77)] =       'Cytosine from CC3'
TAGL['%s-%s-monoA-unCP' % (dbse,  77)] =       'Cytosine from CC3'
TAGL['%s-%s-monoB-unCP' % (dbse,  77)] =       'Cytosine from CC3'
TAGL['%s-%s'            % (dbse,  78)] = 'ST-08 CC4'
TAGL['%s-%s-dimer'      % (dbse,  78)] =       'CC4'
TAGL['%s-%s-monoA-CP'   % (dbse,  78)] =       'Cytosine from CC4'
TAGL['%s-%s-monoB-CP'   % (dbse,  78)] =       'Cytosine from CC4'
TAGL['%s-%s-monoA-unCP' % (dbse,  78)] =       'Cytosine from CC4'
TAGL['%s-%s-monoB-unCP' % (dbse,  78)] =       'Cytosine from CC4'
TAGL['%s-%s'            % (dbse,  79)] = 'ST-09 CC5'
TAGL['%s-%s-dimer'      % (dbse,  79)] =       'CC5'
TAGL['%s-%s-monoA-CP'   % (dbse,  79)] =       'Cytosine from CC5'
TAGL['%s-%s-monoB-CP'   % (dbse,  79)] =       'Cytosine from CC5'
TAGL['%s-%s-monoA-unCP' % (dbse,  79)] =       'Cytosine from CC5'
TAGL['%s-%s-monoB-unCP' % (dbse,  79)] =       'Cytosine from CC5'
TAGL['%s-%s'            % (dbse,  80)] = 'ST-10 CC6'
TAGL['%s-%s-dimer'      % (dbse,  80)] =       'CC6'
TAGL['%s-%s-monoA-CP'   % (dbse,  80)] =       'Cytosine from CC6'
TAGL['%s-%s-monoB-CP'   % (dbse,  80)] =       'Cytosine from CC6'
TAGL['%s-%s-monoA-unCP' % (dbse,  80)] =       'Cytosine from CC6'
TAGL['%s-%s-monoB-unCP' % (dbse,  80)] =       'Cytosine from CC6'
TAGL['%s-%s'            % (dbse,  81)] = 'ST-11 CC7'
TAGL['%s-%s-dimer'      % (dbse,  81)] =       'CC7'
TAGL['%s-%s-monoA-CP'   % (dbse,  81)] =       'Cytosine from CC7'
TAGL['%s-%s-monoB-CP'   % (dbse,  81)] =       'Cytosine from CC7'
TAGL['%s-%s-monoA-unCP' % (dbse,  81)] =       'Cytosine from CC7'
TAGL['%s-%s-monoB-unCP' % (dbse,  81)] =       'Cytosine from CC7'
TAGL['%s-%s'            % (dbse,  82)] = 'ST-12 CC8'
TAGL['%s-%s-dimer'      % (dbse,  82)] =       'CC8'
TAGL['%s-%s-monoA-CP'   % (dbse,  82)] =       'Cytosine from CC8'
TAGL['%s-%s-monoB-CP'   % (dbse,  82)] =       'Cytosine from CC8'
TAGL['%s-%s-monoA-unCP' % (dbse,  82)] =       'Cytosine from CC8'
TAGL['%s-%s-monoB-unCP' % (dbse,  82)] =       'Cytosine from CC8'
TAGL['%s-%s'            % (dbse,  83)] = 'ST-13 CC9'
TAGL['%s-%s-dimer'      % (dbse,  83)] =       'CC9'
TAGL['%s-%s-monoA-CP'   % (dbse,  83)] =       'Cytosine from CC9'
TAGL['%s-%s-monoB-CP'   % (dbse,  83)] =       'Cytosine from CC9'
TAGL['%s-%s-monoA-unCP' % (dbse,  83)] =       'Cytosine from CC9'
TAGL['%s-%s-monoB-unCP' % (dbse,  83)] =       'Cytosine from CC9'
TAGL['%s-%s'            % (dbse,  84)] = 'ST-14 CC10'
TAGL['%s-%s-dimer'      % (dbse,  84)] =       'CC10'
TAGL['%s-%s-monoA-CP'   % (dbse,  84)] =       'Cytosine from CC10'
TAGL['%s-%s-monoB-CP'   % (dbse,  84)] =       'Cytosine from CC10'
TAGL['%s-%s-monoA-unCP' % (dbse,  84)] =       'Cytosine from CC10'
TAGL['%s-%s-monoB-unCP' % (dbse,  84)] =       'Cytosine from CC10'
TAGL['%s-%s'            % (dbse,  85)] = 'ST-15 CC11'
TAGL['%s-%s-dimer'      % (dbse,  85)] =       'CC11'
TAGL['%s-%s-monoA-CP'   % (dbse,  85)] =       'Cytosine from CC11'
TAGL['%s-%s-monoB-CP'   % (dbse,  85)] =       'Cytosine from CC11'
TAGL['%s-%s-monoA-unCP' % (dbse,  85)] =       'Cytosine from CC11'
TAGL['%s-%s-monoB-unCP' % (dbse,  85)] =       'Cytosine from CC11'
TAGL['%s-%s'            % (dbse,  86)] = 'ST-16 CC12'
TAGL['%s-%s-dimer'      % (dbse,  86)] =       'CC12'
TAGL['%s-%s-monoA-CP'   % (dbse,  86)] =       'Cytosine from CC12'
TAGL['%s-%s-monoB-CP'   % (dbse,  86)] =       'Cytosine from CC12'
TAGL['%s-%s-monoA-unCP' % (dbse,  86)] =       'Cytosine from CC12'
TAGL['%s-%s-monoB-unCP' % (dbse,  86)] =       'Cytosine from CC12'
TAGL['%s-%s'            % (dbse,  87)] = 'ST-17 CC13'
TAGL['%s-%s-dimer'      % (dbse,  87)] =       'CC13'
TAGL['%s-%s-monoA-CP'   % (dbse,  87)] =       'Cytosine from CC13'
TAGL['%s-%s-monoB-CP'   % (dbse,  87)] =       'Cytosine from CC13'
TAGL['%s-%s-monoA-unCP' % (dbse,  87)] =       'Cytosine from CC13'
TAGL['%s-%s-monoB-unCP' % (dbse,  87)] =       'Cytosine from CC13'
TAGL['%s-%s'            % (dbse,  88)] = 'ST-18 CC14'
TAGL['%s-%s-dimer'      % (dbse,  88)] =       'CC14'
TAGL['%s-%s-monoA-CP'   % (dbse,  88)] =       'Cytosine from CC14'
TAGL['%s-%s-monoB-CP'   % (dbse,  88)] =       'Cytosine from CC14'
TAGL['%s-%s-monoA-unCP' % (dbse,  88)] =       'Cytosine from CC14'
TAGL['%s-%s-monoB-unCP' % (dbse,  88)] =       'Cytosine from CC14'
TAGL['%s-%s'            % (dbse,  89)] = 'ST-19 AAst'
TAGL['%s-%s-dimer'      % (dbse,  89)] =       'AAst'
TAGL['%s-%s-monoA-CP'   % (dbse,  89)] =       'Adenine from AAst'
TAGL['%s-%s-monoB-CP'   % (dbse,  89)] =       'Adenine from AAst'
TAGL['%s-%s-monoA-unCP' % (dbse,  89)] =       'Adenine from AAst'
TAGL['%s-%s-monoB-unCP' % (dbse,  89)] =       'Adenine from AAst'
TAGL['%s-%s'            % (dbse,  90)] = 'ST-20 GGst'
TAGL['%s-%s-dimer'      % (dbse,  90)] =       'GGst'
TAGL['%s-%s-monoA-CP'   % (dbse,  90)] =       'Guanine from GGst'
TAGL['%s-%s-monoB-CP'   % (dbse,  90)] =       'Guanine from GGst'
TAGL['%s-%s-monoA-unCP' % (dbse,  90)] =       'Guanine from GGst'
TAGL['%s-%s-monoB-unCP' % (dbse,  90)] =       'Guanine from GGst'
TAGL['%s-%s'            % (dbse,  91)] = 'ST-21 ACst'
TAGL['%s-%s-dimer'      % (dbse,  91)] =       'ACst'
TAGL['%s-%s-monoA-CP'   % (dbse,  91)] =       'Adenine from ACst'
TAGL['%s-%s-monoB-CP'   % (dbse,  91)] =       'Cytosine from ACst'
TAGL['%s-%s-monoA-unCP' % (dbse,  91)] =       'Adenine from ACst'
TAGL['%s-%s-monoB-unCP' % (dbse,  91)] =       'Cytosine from ACst'
TAGL['%s-%s'            % (dbse,  92)] = 'ST-22 GAst'
TAGL['%s-%s-dimer'      % (dbse,  92)] =       'GAst'
TAGL['%s-%s-monoA-CP'   % (dbse,  92)] =       'Guanine from GAst'
TAGL['%s-%s-monoB-CP'   % (dbse,  92)] =       'Adenine from GAst'
TAGL['%s-%s-monoA-unCP' % (dbse,  92)] =       'Guanine from GAst'
TAGL['%s-%s-monoB-unCP' % (dbse,  92)] =       'Adenine from GAst'
TAGL['%s-%s'            % (dbse,  93)] = 'ST-23 CCst'
TAGL['%s-%s-dimer'      % (dbse,  93)] =       'CCst'
TAGL['%s-%s-monoA-CP'   % (dbse,  93)] =       'Cytosine from CCst'
TAGL['%s-%s-monoB-CP'   % (dbse,  93)] =       'Cytosine from CCst'
TAGL['%s-%s-monoA-unCP' % (dbse,  93)] =       'Cytosine from CCst'
TAGL['%s-%s-monoB-unCP' % (dbse,  93)] =       'Cytosine from CCst'
TAGL['%s-%s'            % (dbse,  94)] = 'ST-24 AUst'
TAGL['%s-%s-dimer'      % (dbse,  94)] =       'AUst'
TAGL['%s-%s-monoA-CP'   % (dbse,  94)] =       'Adenine from AUst'
TAGL['%s-%s-monoB-CP'   % (dbse,  94)] =       'Uracil from AUst'
TAGL['%s-%s-monoA-unCP' % (dbse,  94)] =       'Adenine from AUst'
TAGL['%s-%s-monoB-unCP' % (dbse,  94)] =       'Uracil from AUst'
TAGL['%s-%s'            % (dbse,  95)] = 'ST-25 GCst'
TAGL['%s-%s-dimer'      % (dbse,  95)] =       'GCst'
TAGL['%s-%s-monoA-CP'   % (dbse,  95)] =       'Guanine from GCst'
TAGL['%s-%s-monoB-CP'   % (dbse,  95)] =       'Cytosine from GCst'
TAGL['%s-%s-monoA-unCP' % (dbse,  95)] =       'Guanine from GCst'
TAGL['%s-%s-monoB-unCP' % (dbse,  95)] =       'Cytosine from GCst'
TAGL['%s-%s'            % (dbse,  96)] = 'ST-26 CUst'
TAGL['%s-%s-dimer'      % (dbse,  96)] =       'CUst'
TAGL['%s-%s-monoA-CP'   % (dbse,  96)] =       'Cytosine from CUst'
TAGL['%s-%s-monoB-CP'   % (dbse,  96)] =       'Uracil from CUst'
TAGL['%s-%s-monoA-unCP' % (dbse,  96)] =       'Cytosine from CUst'
TAGL['%s-%s-monoB-unCP' % (dbse,  96)] =       'Uracil from CUst'
TAGL['%s-%s'            % (dbse,  97)] = 'ST-27 UUst'
TAGL['%s-%s-dimer'      % (dbse,  97)] =       'UUst'
TAGL['%s-%s-monoA-CP'   % (dbse,  97)] =       'Uracil from UUst'
TAGL['%s-%s-monoB-CP'   % (dbse,  97)] =       'Uracil from UUst'
TAGL['%s-%s-monoA-unCP' % (dbse,  97)] =       'Uracil from UUst'
TAGL['%s-%s-monoB-unCP' % (dbse,  97)] =       'Uracil from UUst'
TAGL['%s-%s'            % (dbse,  98)] = 'ST-28 GUst'
TAGL['%s-%s-dimer'      % (dbse,  98)] =       'GUst'
TAGL['%s-%s-monoA-CP'   % (dbse,  98)] =       'Guanine from GUst'
TAGL['%s-%s-monoB-CP'   % (dbse,  98)] =       'Uracil from GUst'
TAGL['%s-%s-monoA-unCP' % (dbse,  98)] =       'Guanine from GUst'
TAGL['%s-%s-monoB-unCP' % (dbse,  98)] =       'Uracil from GUst'
TAGL['%s-%s'            % (dbse,  99)] = 'ST-29 GG0/3.36 GGs036'
TAGL['%s-%s-dimer'      % (dbse,  99)] =       'GGs036'
TAGL['%s-%s-monoA-CP'   % (dbse,  99)] =       'Guanine from GGs036'
TAGL['%s-%s-monoB-CP'   % (dbse,  99)] =       'Guanine from GGs036'
TAGL['%s-%s-monoA-unCP' % (dbse,  99)] =       'Guanine from GGs036'
TAGL['%s-%s-monoB-unCP' % (dbse,  99)] =       'Guanine from GGs036'
TAGL['%s-%s'            % (dbse, 100)] = 'ST-30 GG0/3.36 CCs036'
TAGL['%s-%s-dimer'      % (dbse, 100)] =       'CCs036'
TAGL['%s-%s-monoA-CP'   % (dbse, 100)] =       'Cytosine from CCs036'
TAGL['%s-%s-monoB-CP'   % (dbse, 100)] =       'Cytosine from CCs036'
TAGL['%s-%s-monoA-unCP' % (dbse, 100)] =       'Cytosine from CCs036'
TAGL['%s-%s-monoB-unCP' % (dbse, 100)] =       'Cytosine from CCs036'
TAGL['%s-%s'            % (dbse, 101)] = 'ST-31 AA20/3.05 AAs2005'
TAGL['%s-%s-dimer'      % (dbse, 101)] =       'AAs2005'
TAGL['%s-%s-monoA-CP'   % (dbse, 101)] =       'Adenine from AAs2005'
TAGL['%s-%s-monoB-CP'   % (dbse, 101)] =       'Adenine from AAs2005'
TAGL['%s-%s-monoA-unCP' % (dbse, 101)] =       'Adenine from AAs2005'
TAGL['%s-%s-monoB-unCP' % (dbse, 101)] =       'Adenine from AAs2005'
TAGL['%s-%s'            % (dbse, 102)] = 'ST-32 AA20/3.05 TTs2005'
TAGL['%s-%s-dimer'      % (dbse, 102)] =       'TTs2005'
TAGL['%s-%s-monoA-CP'   % (dbse, 102)] =       'Thymine from TTs2005'
TAGL['%s-%s-monoB-CP'   % (dbse, 102)] =       'Thymine from TTs2005'
TAGL['%s-%s-monoA-unCP' % (dbse, 102)] =       'Thymine from TTs2005'
TAGL['%s-%s-monoB-unCP' % (dbse, 102)] =       'Thymine from TTs2005'
TAGL['%s-%s'            % (dbse, 103)] = 'ST-33 GC0/3.25 G//Cs'
TAGL['%s-%s-dimer'      % (dbse, 103)] =       'GC0/3.25 G//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 103)] =       'Cytosine from GC0/3.25 G//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 103)] =       'Guanine from GC0/3.25 G//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 103)] =       'Cytosine from GC0/3.25 G//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 103)] =       'Guanine from GC0/3.25 G//Cs'
TAGL['%s-%s'            % (dbse, 104)] = 'ST-34 CG0/3.19 G//Cs'
TAGL['%s-%s-dimer'      % (dbse, 104)] =       'CG0/3.19 G//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 104)] =       'Cytosine from CG0/3.19 G//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 104)] =       'Guanine from CG0/3.19 G//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 104)] =       'Cytosine from CG0/3.19 G//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 104)] =       'Guanine from CG0/3.19 G//Cs'
TAGL['%s-%s'            % (dbse, 105)] = 'ST-35 GA10/3.15 A//Gs'
TAGL['%s-%s-dimer'      % (dbse, 105)] =       'GA10/3.15 A//Gs'
TAGL['%s-%s-monoA-CP'   % (dbse, 105)] =       'Adenine from GA10/3.15 A//Gs'
TAGL['%s-%s-monoB-CP'   % (dbse, 105)] =       'Guanine from GA10/3.15 A//Gs'
TAGL['%s-%s-monoA-unCP' % (dbse, 105)] =       'Adenine from GA10/3.15 A//Gs'
TAGL['%s-%s-monoB-unCP' % (dbse, 105)] =       'Guanine from GA10/3.15 A//Gs'
TAGL['%s-%s'            % (dbse, 106)] = 'ST-36 GA10/3.15 T//Cs'
TAGL['%s-%s-dimer'      % (dbse, 106)] =       'GA10/3.15 T//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 106)] =       'Thymine from GA10/3.15 T//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 106)] =       'Cytosine from GA10/3.15 T//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 106)] =       'Thymine from GA10/3.15 T//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 106)] =       'Cytosine from GA10/3.15 T//Cs'
TAGL['%s-%s'            % (dbse, 107)] = 'ST-37 AG08/3.19 A//Gs'
TAGL['%s-%s-dimer'      % (dbse, 107)] =       'AG08/3.19 A//Gs'
TAGL['%s-%s-monoA-CP'   % (dbse, 107)] =       'Adenine from AG08/3.19 A//Gs'
TAGL['%s-%s-monoB-CP'   % (dbse, 107)] =       'Guanine from AG08/3.19 A//Gs'
TAGL['%s-%s-monoA-unCP' % (dbse, 107)] =       'Adenine from AG08/3.19 A//Gs'
TAGL['%s-%s-monoB-unCP' % (dbse, 107)] =       'Guanine from AG08/3.19 A//Gs'
TAGL['%s-%s'            % (dbse, 108)] = 'ST-38 AG08/3.19 T//Cs'
TAGL['%s-%s-dimer'      % (dbse, 108)] =       'AG08/3.19 T//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 108)] =       'Thymine from AG08/3.19 T//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 108)] =       'Cytosine from AG08/3.19 T//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 108)] =       'Thymine from AG08/3.19 T//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 108)] =       'Cytosine from AG08/3.19 T//Cs'
TAGL['%s-%s'            % (dbse, 109)] = 'ST-39 TG03.19 T//Gs'
TAGL['%s-%s-dimer'      % (dbse, 109)] =       'TG03.19 T//Gs'
TAGL['%s-%s-monoA-CP'   % (dbse, 109)] =       'Thymine from TG03.19 T//Gs'
TAGL['%s-%s-monoB-CP'   % (dbse, 109)] =       'Guanine from TG03.19 T//Gs'
TAGL['%s-%s-monoA-unCP' % (dbse, 109)] =       'Thymine from TG03.19 T//Gs'
TAGL['%s-%s-monoB-unCP' % (dbse, 109)] =       'Guanine from TG03.19 T//Gs'
TAGL['%s-%s'            % (dbse, 110)] = 'ST-40 TG03.19 A//Cs'
TAGL['%s-%s-dimer'      % (dbse, 110)] =       'TG03.19 A//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 110)] =       'Adenine from TG03.19 A//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 110)] =       'Cytosine from TG03.19 A//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 110)] =       'Adenine from TG03.19 A//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 110)] =       'Cytosine from TG03.19 A//Cs'
TAGL['%s-%s'            % (dbse, 111)] = 'ST-41 GT10/3.15 T//Gs'
TAGL['%s-%s-dimer'      % (dbse, 111)] =       'GT10/3.15 T//Gs'
TAGL['%s-%s-monoA-CP'   % (dbse, 111)] =       'Thymine from GT10/3.15 T//Gs'
TAGL['%s-%s-monoB-CP'   % (dbse, 111)] =       'Guanine from GT10/3.15 T//Gs'
TAGL['%s-%s-monoA-unCP' % (dbse, 111)] =       'Thymine from GT10/3.15 T//Gs'
TAGL['%s-%s-monoB-unCP' % (dbse, 111)] =       'Guanine from GT10/3.15 T//Gs'
TAGL['%s-%s'            % (dbse, 112)] = 'ST-42 GT10/3.15 A//Cs'
TAGL['%s-%s-dimer'      % (dbse, 112)] =       'GT10/3.15 A//Cs'
TAGL['%s-%s-monoA-CP'   % (dbse, 112)] =       'Adenine from GT10/3.15 A//Cs'
TAGL['%s-%s-monoB-CP'   % (dbse, 112)] =       'Cytosine from GT10/3.15 A//Cs'
TAGL['%s-%s-monoA-unCP' % (dbse, 112)] =       'Adenine from GT10/3.15 A//Cs'
TAGL['%s-%s-monoB-unCP' % (dbse, 112)] =       'Cytosine from GT10/3.15 A//Cs'
TAGL['%s-%s'            % (dbse, 113)] = 'ST-43 AT10/3.26 A//Ts'
TAGL['%s-%s-dimer'      % (dbse, 113)] =       'AT10/3.26 A//Ts'
TAGL['%s-%s-monoA-CP'   % (dbse, 113)] =       'Adenine from AT10/3.26 A//Ts'
TAGL['%s-%s-monoB-CP'   % (dbse, 113)] =       'Thymine from AT10/3.26 A//Ts'
TAGL['%s-%s-monoA-unCP' % (dbse, 113)] =       'Adenine from AT10/3.26 A//Ts'
TAGL['%s-%s-monoB-unCP' % (dbse, 113)] =       'Thymine from AT10/3.26 A//Ts'
TAGL['%s-%s'            % (dbse, 114)] = 'ST-44 TA08/3.16 A//Ts'
TAGL['%s-%s-dimer'      % (dbse, 114)] =       'TA08/3.16 A//Ts'
TAGL['%s-%s-monoA-CP'   % (dbse, 114)] =       'Adenine from TA08/3.16 A//Ts'
TAGL['%s-%s-monoB-CP'   % (dbse, 114)] =       'Thymine from TA08/3.16 A//Ts'
TAGL['%s-%s-monoA-unCP' % (dbse, 114)] =       'Adenine from TA08/3.16 A//Ts'
TAGL['%s-%s-monoB-unCP' % (dbse, 114)] =       'Thymine from TA08/3.16 A//Ts'
TAGL['%s-%s'            % (dbse, 115)] = 'ST-45 AA0/3.24 A//As'
TAGL['%s-%s-dimer'      % (dbse, 115)] =       'AA0/3.24 A//As'
TAGL['%s-%s-monoA-CP'   % (dbse, 115)] =       'Adenine from AA0/3.24 A//As'
TAGL['%s-%s-monoB-CP'   % (dbse, 115)] =       'Adenine from AA0/3.24 A//As'
TAGL['%s-%s-monoA-unCP' % (dbse, 115)] =       'Adenine from AA0/3.24 A//As'
TAGL['%s-%s-monoB-unCP' % (dbse, 115)] =       'Adenine from AA0/3.24 A//As'
TAGL['%s-%s'            % (dbse, 116)] = 'ST-46 AA0/3.24 T//Ts'
TAGL['%s-%s-dimer'      % (dbse, 116)] =       'AA0/3.24 T//Ts'
TAGL['%s-%s-monoA-CP'   % (dbse, 116)] =       'Thymine from AA0/3.24 T//Ts'
TAGL['%s-%s-monoB-CP'   % (dbse, 116)] =       'Thymine from AA0/3.24 T//Ts'
TAGL['%s-%s-monoA-unCP' % (dbse, 116)] =       'Thymine from AA0/3.24 T//Ts'
TAGL['%s-%s-monoB-unCP' % (dbse, 116)] =       'Thymine from AA0/3.24 T//Ts'
TAGL['%s-%s'            % (dbse, 117)] = 'ST-47 A...T S'
TAGL['%s-%s-dimer'      % (dbse, 117)] =       'A...T S'
TAGL['%s-%s-monoA-CP'   % (dbse, 117)] =       'methyl-Adenine from A...T S'
TAGL['%s-%s-monoB-CP'   % (dbse, 117)] =       'methyl-Thymine from A...T S'
TAGL['%s-%s-monoA-unCP' % (dbse, 117)] =       'methyl-Adenine from A...T S'
TAGL['%s-%s-monoB-unCP' % (dbse, 117)] =       'methyl-Thymine from A...T S'
TAGL['%s-%s'            % (dbse, 118)] = 'ST-48 G...C S'
TAGL['%s-%s-dimer'      % (dbse, 118)] =       'G...C S'
TAGL['%s-%s-monoA-CP'   % (dbse, 118)] =       'methyl-Cytosine from G...C S'
TAGL['%s-%s-monoB-CP'   % (dbse, 118)] =       'methyl-Guanine from G...C S'
TAGL['%s-%s-monoA-unCP' % (dbse, 118)] =       'methyl-Cytosine from G...C S'
TAGL['%s-%s-monoB-unCP' % (dbse, 118)] =       'methyl-Guanine from G...C S'
TAGL['%s-%s'            % (dbse, 119)] = 'ST-49 A...C S'
TAGL['%s-%s-dimer'      % (dbse, 119)] =       'A...C S'
TAGL['%s-%s-monoA-CP'   % (dbse, 119)] =       'methyl-Adenine from A...C S'
TAGL['%s-%s-monoB-CP'   % (dbse, 119)] =       'methyl-Cytosine from A...C S'
TAGL['%s-%s-monoA-unCP' % (dbse, 119)] =       'methyl-Adenine from A...C S'
TAGL['%s-%s-monoB-unCP' % (dbse, 119)] =       'methyl-Cytosine from A...C S'
TAGL['%s-%s'            % (dbse, 120)] = 'ST-50 T...G S'
TAGL['%s-%s-dimer'      % (dbse, 120)] =       'T...G S'
TAGL['%s-%s-monoA-CP'   % (dbse, 120)] =       'methyl-Thymine from T...G S'
TAGL['%s-%s-monoB-CP'   % (dbse, 120)] =       'methyl-Guanine from T...G S'
TAGL['%s-%s-monoA-unCP' % (dbse, 120)] =       'methyl-Thymine from T...G S'
TAGL['%s-%s-monoB-unCP' % (dbse, 120)] =       'methyl-Guanine from T...G S'
TAGL['%s-%s'            % (dbse, 121)] = 'ST-51 G...C S'
TAGL['%s-%s-dimer'      % (dbse, 121)] =       'G...C S'
TAGL['%s-%s-monoA-CP'   % (dbse, 121)] =       'Cytosine from G...C S'
TAGL['%s-%s-monoB-CP'   % (dbse, 121)] =       'Guanine from G...C S'
TAGL['%s-%s-monoA-unCP' % (dbse, 121)] =       'Cytosine from G...C S'
TAGL['%s-%s-monoB-unCP' % (dbse, 121)] =       'Guanine from G...C S'
TAGL['%s-%s'            % (dbse, 122)] = 'ST-52 A...G S'
TAGL['%s-%s-dimer'      % (dbse, 122)] =       'A...G S'
TAGL['%s-%s-monoA-CP'   % (dbse, 122)] =       'Adenine from A...G S'
TAGL['%s-%s-monoB-CP'   % (dbse, 122)] =       'Guanine from A...G S'
TAGL['%s-%s-monoA-unCP' % (dbse, 122)] =       'Adenine from A...G S'
TAGL['%s-%s-monoB-unCP' % (dbse, 122)] =       'Guanine from A...G S'
TAGL['%s-%s'            % (dbse, 123)] = 'ST-53 C...G S'
TAGL['%s-%s-dimer'      % (dbse, 123)] =       'C...G S'
TAGL['%s-%s-monoA-CP'   % (dbse, 123)] =       'Guanine from C...G S'
TAGL['%s-%s-monoB-CP'   % (dbse, 123)] =       'Cytosine from C...G S'
TAGL['%s-%s-monoA-unCP' % (dbse, 123)] =       'Guanine from C...G S'
TAGL['%s-%s-monoB-unCP' % (dbse, 123)] =       'Cytosine from C...G S'
TAGL['%s-%s'            % (dbse, 124)] = 'ST-54 G...C S'
TAGL['%s-%s-dimer'      % (dbse, 124)] =       'G...C S'
TAGL['%s-%s-monoA-CP'   % (dbse, 124)] =       'Guanine from G...C S'
TAGL['%s-%s-monoB-CP'   % (dbse, 124)] =       'Cytosine from G...C S'
TAGL['%s-%s-monoA-unCP' % (dbse, 124)] =       'Guanine from G...C S'
TAGL['%s-%s-monoB-unCP' % (dbse, 124)] =       'Cytosine from G...C S'

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

JSCH_1 = input.process_input("""
molecule dimer {
0 1
C    -1.0398599    -0.0950435     2.9628987
N    -0.8760506    -0.1198953     4.3522101
C     0.3372729    -0.0573522     4.9526643
C     1.4603152     0.0294729     4.2021231
C     1.2876371     0.0522766     2.7771415
N     0.0866353    -0.0006919     2.2061593
O    -2.1779850    -0.1592983     2.4996990
N     2.3517978     0.1313296     1.9777210
H    -1.7254816    -0.1869061     4.8897274
H     0.3482118    -0.0833071     6.0321432
H     2.4345221     0.0778275     4.6597911
H     3.2714721     0.1534551     2.3764404
H     2.2350290     0.1077513     0.9551229
--
0 1
O     2.0171439     0.0263963    -0.7905108
C     0.9445057     0.0313388    -1.4013109
N    -0.2671137     0.0963439    -0.7051367
C    -1.5207327     0.1136461    -1.2552546
N    -1.7528129     0.0544172    -2.5494108
C    -0.6040129    -0.0113445    -3.2574879
C     0.7161247    -0.0244271    -2.8113172
N     1.6041685    -0.0981114    -3.8601422
C     0.8295480    -0.1292217    -4.9265187
N    -0.5075993    -0.0802063    -4.6198760
N    -2.5513427     0.2447649    -0.3850923
H    -0.1820496     0.1041077     0.3219703
H     1.1760819    -0.1871623    -5.9443460
H    -1.2844954    -0.0872596    -5.2590531
H    -3.4573855     0.0691895    -0.7801319
H    -2.4169221     0.0545062     0.6045745
units angstrom
}
""")

JSCH_2 = input.process_input("""
molecule dimer {
0 1
C    -0.8133331    -0.0866715     2.9789275
N    -0.6596824    -0.0799211     4.3738426
C     0.5690989     0.0065483     4.9361002
C     1.6895413     0.0878038     4.1758957
C     1.5121802     0.0788250     2.7574552
N     0.3010370     0.0009434     2.2093349
O    -1.9566795    -0.1721648     2.5197622
N     2.5630702     0.1536301     1.9395171
C    -1.8703132    -0.1697062     5.1759469
H     0.5992246     0.0046040     6.0165438
H     2.6640151     0.1533452     4.6318609
H     3.4878053     0.1855364     2.3258249
H     2.4345725     0.0922431     0.9177005
H    -2.4003892    -1.0888680     4.9431219
H    -2.5250547     0.6677091     4.9519253
H    -1.5896839    -0.1549564     6.2247983
--
0 1
O     2.2125806    -0.0495287    -0.7919145
C     1.1295844    -0.0091268    -1.3874888
N    -0.0672022     0.0855128    -0.6713507
C    -1.3285435     0.1453481    -1.2005422
N    -1.5812601     0.1027260    -2.4911624
C    -0.4460576     0.0088365    -3.2186883
C     0.8814078    -0.0466784    -2.7933188
N     1.7434051    -0.1366731    -3.8585399
C     0.9419745    -0.1354113    -4.9081560
N    -0.3886567    -0.0498363    -4.5825778
N    -2.3394192     0.3029765    -0.3108597
H     0.0331014     0.0898348     0.3565931
H     1.2604461    -0.1939844    -5.9363793
C    -1.5258946    -0.0225882    -5.4753265
H    -3.2566158     0.1569838    -0.6925288
H    -2.1947491     0.0836423     0.6730029
H    -2.1741480    -0.8737721    -5.2827842
H    -2.0939825     0.8923400    -5.3277349
H    -1.1628640    -0.0654072    -6.4980769
units angstrom
}
""")

JSCH_3 = input.process_input("""
molecule dimer {
0 1
N     0.9350155    -0.0279801    -0.3788916
C     1.6739638    -0.0357766     0.7424316
C     3.0747955    -0.0094480     0.5994562
C     3.5646109     0.0195446    -0.7059872
N     2.8531510     0.0258031    -1.8409596
C     1.5490760     0.0012569    -1.5808009
N     4.0885824    -0.0054429     1.5289786
C     5.1829921     0.0253971     0.7872176
N     4.9294871     0.0412404    -0.5567274
N     1.0716177    -0.0765366     1.9391390
H     0.8794435     0.0050260    -2.4315709
H     6.1882591     0.0375542     1.1738824
H     5.6035368     0.0648755    -1.3036811
H     0.0586915    -0.0423765     2.0039181
H     1.6443796    -0.0347395     2.7619159
--
0 1
N    -3.9211729    -0.0009646    -1.5163659
C    -4.6136833     0.0169051    -0.3336520
C    -3.9917387     0.0219348     0.8663338
C    -2.5361367     0.0074651     0.8766724
N    -1.9256484    -0.0110593    -0.3638948
C    -2.5395897    -0.0149474    -1.5962357
C    -4.7106131     0.0413373     2.1738637
O    -1.8674730     0.0112093     1.9120833
O    -1.9416783    -0.0291878    -2.6573783
H    -4.4017172    -0.0036078    -2.4004924
H    -0.8838255    -0.0216168    -0.3784269
H    -5.6909220     0.0269347    -0.4227183
H    -4.4439282    -0.8302573     2.7695655
H    -4.4267056     0.9186178     2.7530256
H    -5.7883971     0.0505530     2.0247280
units angstrom
}
""")

JSCH_4 = input.process_input("""
molecule dimer {
0 1
N     1.4233678    -2.5755572    -0.0177928
C     2.4164068    -1.6737862    -0.0069340
N     3.6894208    -2.0970069     0.0011813
C     4.6820785    -1.1903949     0.0099205
N     4.6008957     0.1417597     0.0124822
C     3.3223927     0.5333463     0.0042869
C     2.1894368    -0.2816944    -0.0057960
N     1.0473244     0.4790858    -0.0114405
C     1.4850348     1.7325334    -0.0051890
N     2.8426160     1.8225456     0.0046153
H     5.6822446    -1.6040945     0.0161209
H     0.8406886     2.5978945    -0.0067705
C     3.6543123     3.0217765     0.0114721
H     1.6817907    -3.5456829    -0.0088365
H     0.4430037    -2.3144424    -0.0119382
H     4.2913482     3.0300392     0.8917344
H     4.2797691     3.0498799    -0.8767515
H     2.9957930     3.8849549     0.0253861
--
0 1
N    -1.7137952     0.0896100    -0.0160334
C    -2.4110219     1.2713542    -0.0102985
N    -3.7875706     1.1256389    -0.0070257
C    -4.3748768    -0.1132637     0.0070053
C    -3.6761626    -1.2738156     0.0090266
C    -2.2277642    -1.1903377    -0.0021447
O    -1.8742013     2.3720388    -0.0087627
O    -1.4812933    -2.1730005    -0.0012821
C    -4.3178351    -2.6216898     0.0206209
C    -4.5749913     2.3488054     0.0099294
H    -0.6740657     0.1873910    -0.0185734
H    -5.4567167    -0.1063625     0.0145634
H    -4.0169219    -3.1954743    -0.8545578
H    -5.4025389    -2.5351168     0.0301255
H    -4.0008440    -3.1876298     0.8951907
H    -4.2676272     2.9948502    -0.8071476
H    -4.4269734     2.8838422     0.9449871
H    -5.6215286     2.0832557    -0.1007218
units angstrom
}
""")

JSCH_5 = input.process_input("""
molecule dimer {
0 1
N    -0.7878100     0.0020606    -4.2304138
C    -0.8092899     0.0015891    -2.8583291
C     0.5077793    -0.0002173    -2.4445540
N     1.3042354    -0.0007671    -3.5776764
C     0.5251204     0.0003646    -4.7183358
C     0.8109948    -0.0022184    -1.0723889
N    -0.3541156    -0.0012118    -0.3048371
C    -1.6266470     0.0012022    -0.7959383
N    -1.9085135     0.0021412    -2.0882229
O     1.9375199    -0.0047453    -0.5373376
N    -2.6234345     0.0038170     0.1093233
O     0.8790764    -0.0004039    -5.8828291
H    -1.5908258     0.0034199    -4.8351827
H    -0.2139163    -0.0018761     0.7197062
H    -3.5584869     0.0004853    -0.2509512
H    -2.4522376    -0.0015362     1.1120222
H     2.3082381    -0.0028321    -3.6013832
--
0 1
N     2.3837158     0.0024084     2.1880115
C     1.3635083     0.0013168     3.0441674
C     1.6066986     0.0037303     4.4583838
C     0.5209494     0.0030055     5.2667859
N    -0.7245613    -0.0002446     4.7321333
C    -0.9573175    -0.0028235     3.3551883
N     0.1301240    -0.0019248     2.5406741
O    -2.1208999    -0.0056022     2.9504908
H    -1.5467235    -0.0008598     5.3145350
H     0.5868645     0.0048893     6.3446020
H     2.6043649     0.0062410     4.8655851
H     3.3229417     0.0045233     2.5396202
H     2.2148236     0.0000552     1.1689293
units angstrom
}
""")

JSCH_6 = input.process_input("""
molecule dimer {
0 1
C    -0.0399020     0.0000000    -0.0353727
N    -0.0114814     0.0000000     1.3676751
C     1.1387066     0.0000000     2.0831816
C     2.3367544     0.0000000     1.4511865
C     2.3093819     0.0000000     0.0167889
N     1.1708246     0.0000000    -0.6653863
O    -1.1150036     0.0000000    -0.6203815
N     3.4490015     0.0000000    -0.6790500
H    -0.9108314     0.0000000     1.8214404
H     1.0410884     0.0000000     3.1587270
H     3.2607830     0.0000000     2.0056829
H     4.3278873     0.0000000    -0.1979709
H     3.4179786     0.0000000    -1.7058513
--
0 1
O     3.2403090     0.0000000    -3.4870523
C     2.1818608     0.0000000    -4.1224126
N     0.9469308     0.0000000    -3.4715093
C    -0.2749703     0.0000000    -4.0613039
N    -0.4833351     0.0000000    -5.3532242
C     0.6850532     0.0000000    -6.0363493
C     1.9918547     0.0000000    -5.5441185
N     2.9083577     0.0000000    -6.5651625
C     2.1671699     0.0000000    -7.6581930
N     0.8245995     0.0000000    -7.3970665
H    -1.1112956     0.0000000    -3.3753120
H     0.9890303     0.0000000    -2.4307715
H     2.5507838     0.0000000    -8.6645924
H     0.0674277     0.0000000    -8.0602299
units angstrom
}
""")

JSCH_7 = input.process_input("""
molecule dimer {
0 1
O    -1.3445145    -0.0017812     0.2109785
C    -0.5827723    -0.0011739     1.1802676
N     0.8024859    -0.0013909     0.9807253
C     1.7623754    -0.0004370     1.9497163
N     1.5138727     0.0005531     3.2372784
C     0.1802572     0.0006138     3.4773368
C    -0.8864943    -0.0001193     2.5787139
N    -2.0960454     0.0003342     3.2323438
C    -1.7654415     0.0013112     4.5084219
N    -0.4076744     0.0015312     4.7105883
N     3.0436067    -0.0006786     1.5010046
H     1.1059736    -0.0017584    -0.0014454
H    -2.4611653     0.0018946     5.3301928
H     0.0813329     0.0022447     5.5900062
H     3.7807912     0.2008871     2.1781729
H     3.2486220    -0.0002735     0.5197891
--
0 1
O     1.7201107    -0.0009041    -1.6341114
C     0.9160823    -0.0002542    -2.5726633
N    -0.4364135    -0.0009531    -2.4309422
C    -1.3901935    -0.0001633    -3.4689726
C    -0.8112768     0.0014209    -4.8053842
C     0.5255837     0.0020253    -4.9588359
N     1.3652416     0.0012542    -3.8724327
O    -2.5802179    -0.0008229    -3.2195788
H    -1.4757079     0.0020523    -5.6523281
H    -0.8069002    -0.0016737    -1.4589893
H     1.0090483     0.0031967    -5.9240562
H     2.3651020     0.0016683    -3.9825177
units angstrom
}
""")

JSCH_8 = input.process_input("""
molecule dimer {
0 1
C    -0.0546262    -0.0000666    -0.0449366
N    -0.0199900    -0.0000391     1.3370041
C     1.1383570     0.0000076     2.0419758
C     2.3258577     0.0000418     1.3914589
C     2.2850892     0.0000244    -0.0362838
N     1.1291861    -0.0000259    -0.7116098
O    -1.1572546    -0.0001260    -0.6121938
N     3.4163218     0.0000650    -0.7389180
H    -0.9156158    -0.0000626     1.8021973
H     1.0513420     0.0000150     3.1177694
H     3.2570460     0.0000679     1.9330954
H     4.3021685     0.0000746    -0.2670886
H     3.3885124     0.0000089    -1.7541744
--
1 1
N     0.8209043     0.0000137    -3.4723115
C    -0.4088139     0.0000345    -4.0385840
C    -0.4948596     0.0000126    -5.4617196
C     0.6604919    -0.0000154    -6.1675372
N     1.8676936    -0.0000144    -5.5424117
C     2.0133448    -0.0000297    -4.1696950
O     3.1022731    -0.0000760    -3.6244111
N    -1.4596996     0.0000728    -3.2533202
H     2.7287536    -0.0000519    -6.0699661
H     0.6811915    -0.0000377    -7.2470448
H    -1.4500232     0.0000004    -5.9585035
H    -2.3755461     0.0000545    -3.6707993
H    -1.3585129     0.0000190    -2.2095314
H     0.9279995     0.0000010    -2.4136765
units angstrom
}
""")

JSCH_9 = input.process_input("""
molecule dimer {
0 1
O     3.5986069     0.3187715    -0.0000425
C     3.1043656    -0.7907721     0.0000494
N     3.8766606    -1.9463148     0.0000974
C     3.3523333    -3.2081585     0.0000314
C     2.0204180    -3.4182033     0.0000005
C     1.1157326    -2.2823273    -0.0000024
N     1.7481024    -1.0416244     0.0001189
O    -0.1074745    -2.3680037    -0.0001125
H     1.6059877    -4.4112397    -0.0000642
H     4.8710378    -1.7927175    -0.0000238
H     1.1448708    -0.2099855     0.0000657
H     4.0702664    -4.0151053    -0.0000042
--
0 1
O     0.0832848     1.3018469    -0.0000275
C    -1.1439944     1.3419495     0.0000322
N    -1.7902509     2.5769847     0.0001252
C    -3.1492468     2.8315090     0.0000032
N    -3.9059521     1.6745104     0.0001278
C    -3.3673000     0.4141979     0.0000554
C    -2.0348331     0.2020916     0.0000036
O    -3.6310511     3.9463823    -0.0001798
H    -1.6060027    -0.7883399    -0.0000828
H    -4.9024679     1.8155815    -0.0000087
H    -1.1951712     3.3944157     0.0000479
H    -4.0815013    -0.3955360     0.0000308
units angstrom
}
""")

JSCH_10 = input.process_input("""
molecule dimer {
0 1
O     3.0139530    -0.0000663    -2.3714607
C     3.0149686    -0.0000275    -1.1572038
N     4.1880062    -0.0001806    -0.4139005
C     4.2197761     0.0000190     0.9517028
C     3.0868691     0.0000352     1.6837875
C     1.8067676    -0.0000168     1.0072050
N     1.8766809     0.0001650    -0.3776725
O     0.7209062     0.0000058     1.5860958
H     3.1056982    -0.0000395     2.7597452
H     5.0361993     0.0000319    -0.9554725
H     0.9797808    -0.0001787    -0.8806223
H     5.2026538     0.0000504     1.3991588
--
0 1
O    -0.6997337     0.0001361    -1.5583592
C    -1.7610439     0.0000592    -0.9401093
N    -2.9691203     0.0001496    -1.6022795
C    -4.1792891     0.0000396    -0.9575726
C    -4.2572226    -0.0001289     0.3871173
C    -3.0393068    -0.0001349     1.1840453
N    -1.8606585    -0.0000544     0.4209370
O    -2.9932145    -0.0002449     2.4009987
H    -5.2057276     0.0004029     0.8960881
H    -2.9068168     0.0002610    -2.6064637
H    -0.9648574    -0.0000753     0.9279392
H    -5.0495438     0.0000916    -1.5966575
units angstrom
}
""")

JSCH_11 = input.process_input("""
molecule dimer {
0 1
C    -0.0268479     0.0001243    -0.0484048
N     0.0107678    -0.0000464     1.3475182
C     1.1662769    -0.0000863     2.0541365
C     2.3502411    -0.0000220     1.4003156
C     2.3025646     0.0000460    -0.0345481
N     1.1568381     0.0001109    -0.7170898
O    -1.1319705     0.0002807    -0.5958040
N     3.4410032     0.0000681    -0.7244676
H    -0.8870732    -0.0000629     1.8046190
H     1.0808205    -0.0001565     3.1305463
H     3.2850793    -0.0000335     1.9362836
H     4.3196967    -0.0000384    -0.2408021
H     3.4310344     0.0000688    -1.7533114
--
0 1
S     3.6261481     0.0000331    -3.9489529
C     2.0961273    -0.0000288    -4.6037680
N     0.9509789    -0.0001104    -3.8229476
C    -0.3400955    -0.0001574    -4.2939721
N    -0.6487077    -0.0000023    -5.5786621
C     0.4418258     0.0000453    -6.3611558
C     1.7866043    -0.0000027    -5.9848867
N     2.6124570     0.0000014    -7.0854825
C     1.7805999     0.0000520    -8.1059933
N     0.4593064     0.0000986    -7.7272032
N    -1.3218629    -0.0005782    -3.3778714
H     1.0832701    -0.0001308    -2.8047028
H     2.0715160     0.0000504    -9.1429108
H    -0.3508790     0.0000799    -8.3240783
H    -2.2569895     0.0000670    -3.7397613
H    -1.1666925     0.0000488    -2.3691590
units angstrom
}
""")

JSCH_12 = input.process_input("""
molecule dimer {
0 1
N    -2.3081868     0.7091068     0.0000285
C    -1.2051249     1.4749768    -0.0000074
N     0.0068404     0.8955116    -0.0000239
C     1.1080909     1.6831490    -0.0000001
N     1.1718597     3.0108449     0.0000248
C    -0.0467899     3.5587068     0.0000144
C    -1.2685607     2.8863732    -0.0000096
N    -2.3312449     3.7631581    -0.0000201
C    -1.7532251     4.9525322     0.0000002
N    -0.3875112     4.8900390     0.0000301
H     2.0473735     1.1474715    -0.0000027
H    -2.2800394     5.8934739    -0.0000005
H     0.2491479     5.6682121    -0.0000014
H    -2.2351122    -0.2997556    -0.0000060
H    -3.2108389     1.1496550    -0.0000063
--
0 1
S    -1.7524155    -2.8468749    -0.0000015
C    -0.1157899    -3.0959877     0.0000038
N     0.7757259    -2.0472741     0.0000119
C     2.1570215    -2.1260472     0.0000127
N     2.6399748    -3.4231595     0.0000214
C     1.8308010    -4.5235818     0.0000026
C     0.4852439    -4.4036601    -0.0000120
O     2.8928976    -1.1571548    -0.0000022
H    -0.1424437    -5.2790356    -0.0000032
H     3.6448728    -3.5071844    -0.0000023
H     0.4048510    -1.0816361    -0.0000008
H     2.3325010    -5.4794294    -0.0000001
units angstrom
}
""")

JSCH_13 = input.process_input("""
molecule dimer {
0 1
N    -0.9044942     0.3053428    -1.9849463
C    -1.5722006     0.1028596    -0.8342896
N    -0.8984868    -0.0828082     0.3086249
C    -1.5939985    -0.2540305     1.4676013
N    -2.9198181    -0.2548010     1.6283568
C    -3.5512123    -0.0734507     0.4647187
C    -2.9785368     0.1103594    -0.7891549
N    -3.9279869     0.2613244    -1.7761223
C    -5.0691734     0.1716827    -1.1179520
N    -4.9024796    -0.0280544     0.2301645
N    -0.8371284    -0.4982415     2.5786136
H    -6.0473812     0.2435326    -1.5626408
H    -5.6242357    -0.1304918     0.9234162
H     0.0816120     0.0641263    -2.0404453
H    -1.4560458     0.2977892    -2.8251132
H    -1.3448761    -0.4008730     3.4406879
H     0.1151602    -0.1525559     2.5800463
--
0 1
O     1.9075808    -0.3384204    -1.9978152
C     2.5680509    -0.1510537    -0.9753123
N     1.9481026     0.0211081     0.2519027
C     2.5671978     0.2530499     1.4535237
N     3.9422564     0.3000243     1.3759119
C     4.6421414     0.1305464     0.2085191
C     4.0227571    -0.0934717    -0.9708747
O     1.9735147     0.4072231     2.5136021
C     4.7417280    -0.2795724    -2.2650127
H     4.4155554     0.4660468     2.2485693
H     0.8994704    -0.0257608     0.2762592
H     5.7177824     0.1899609     0.2946839
H     4.4976921    -1.2470914    -2.7006324
H     4.4340213     0.4777637    -2.9842064
H     5.8186182    -0.2165735    -2.1235955
units angstrom
}
""")

JSCH_14 = input.process_input("""
molecule dimer {
0 1
N     0.3803518     5.2710590     0.0000000
C    -0.9827443     5.4382698     0.0000000
N    -1.6445695     4.2964704     0.0000000
C    -0.6462489     3.3468246     0.0000000
C     0.6196111     3.9194590     0.0000000
N     1.7959710     3.2884069     0.0000000
C     1.6381971     1.9602261     0.0000000
N     0.4661224     1.2614481     0.0000000
C    -0.6897222     1.9395705     0.0000000
N     2.7687654     1.2167374     0.0000000
N    -1.8599123     1.2852151     0.0000000
H     1.0814580     5.9922927     0.0000000
H    -1.4336046     6.4162469     0.0000000
H    -1.9018421     0.2719784     0.0000000
H    -2.7003404     1.8337683     0.0000000
H     3.6406811     1.7091587     0.0000000
H     2.7459099     0.2072280     0.0000000
--
0 1
C     1.6184078    -2.2819447     0.0000000
N     0.4001260    -1.6517356     0.0000000
C    -0.8434818    -2.2665916     0.0000000
C    -0.8382446    -3.7219904     0.0000000
C     0.3574635    -4.3499319     0.0000000
N     1.5408193    -3.6563635     0.0000000
C    -2.1496426    -4.4344119     0.0000000
O    -1.8782774    -1.6012490     0.0000000
O     2.6941241    -1.6972351     0.0000000
H     0.4220912    -0.6062649     0.0000000
H     2.4263852    -4.1347427     0.0000000
H     0.4460260    -5.4270144     0.0000000
H    -2.7348414    -4.1561542    -0.8748178
H    -2.0041973    -5.5125757     0.0000000
H    -2.7348414    -4.1561542     0.8748178
units angstrom
}
""")

JSCH_15 = input.process_input("""
molecule dimer {
0 1
N    -5.1985541     0.4936739     0.0318901
C    -5.3820983    -0.8616291     0.0378780
H    -6.3663575    -1.2989623     0.0563254
N    -4.2515351    -1.5472672     0.0187374
C    -3.2863426    -0.5669089    -0.0011454
C    -1.8812597    -0.6327273    -0.0265724
N    -1.2256807    -1.8085660    -0.0449899
H    -1.7524474    -2.6605390     0.0028082
H    -0.2208316    -1.8257421    -0.0311968
N    -1.1915263     0.5144581    -0.0390690
C    -1.8701238     1.6793288    -0.0286861
H    -1.2534253     2.5689235    -0.0396068
N    -3.1871042     1.8787826    -0.0064054
C    -3.8427150     0.7125075     0.0068719
H    -5.9109500     1.2044355     0.0439248
--
0 1
C     4.4082682     1.3958429     0.0182886
C     4.9187035     0.0992764     0.0212789
H     5.9905108    -0.0483384     0.0314539
C     4.0880886    -1.0223564     0.0114675
C     4.6130388    -2.4267390     0.0143734
H     4.2620014    -2.9754574     0.8873371
H     4.2783956    -2.9729351    -0.8665709
H     5.7002762    -2.4237810     0.0245680
C     2.7198280    -0.7726204    -0.0014366
F     1.8841246    -1.8434019    -0.0116101
C     2.1541812     0.4899806    -0.0052011
H     1.0780578     0.6238951    -0.0168588
C     3.0326742     1.5626996     0.0050555
F     2.5236686     2.8066583     0.0017348
H     5.0549009     2.2603064     0.0258240
units angstrom
}
""")

JSCH_16 = input.process_input("""
molecule dimer {
0 1
O     0.3144345    -1.1442948     0.0144949
C     1.3535742    -0.4837792     0.0219615
N     1.2957095     0.9167643     0.0619980
C     2.3569978     1.7719209     0.0564464
N     3.6092442     1.3935257     0.0165258
C     3.7175450     0.0431271    -0.0038819
C     2.7154197    -0.9260868    -0.0053966
N     3.2450349    -2.1934744    -0.0364308
C     4.5477687    -1.9938501    -0.0526937
N     4.8854175    -0.6639760    -0.0358567
N     2.0456853     3.1046346     0.1569813
H     0.3457459     1.3051947     0.0833829
H     5.2946011    -2.7690789    -0.0763507
H     5.8090734    -0.2651897    -0.0402909
H     2.8020253     3.7137889    -0.1056218
H     1.1380549     3.3808140    -0.1794301
--
0 1
O    -1.3169188     1.9540889    -0.0350694
C    -2.3291669     1.2492795    -0.0287951
N    -2.3117178    -0.1155253    -0.0025392
C    -3.4195015    -0.9596499     0.0079063
C    -4.6908528    -0.2770302    -0.0148404
C    -4.7376559     1.0709532    -0.0414746
N    -3.5816109     1.8128549    -0.0480085
S    -3.2558282    -2.5902318     0.0425862
H    -5.5930477    -0.8648219    -0.0095256
H    -1.3692770    -0.5503586     0.0090372
H    -5.6601199     1.6316629    -0.0586834
H    -3.6021712     2.8189896    -0.0637914
units angstrom
}
""")

JSCH_17 = input.process_input("""
molecule dimer {
0 1
O     0.2958263    -1.2112383     0.3424002
C     1.3078216    -0.5286166     0.1794526
N     1.2037967     0.8617133     0.0406924
C     2.2285795     1.7353917    -0.1726765
N     3.4850711     1.3842589    -0.2697628
C     3.6408598     0.0475045    -0.1126641
C     2.6766347    -0.9371114     0.0982060
N     3.2448896    -2.1840776     0.1912254
C     4.5343218    -1.9574061     0.0391817
N     4.8260700    -0.6298861    -0.1485794
N     1.8780889     3.0609691    -0.2292539
H     0.2531805     1.2355886     0.1433595
H     5.3034685    -2.7105809     0.0557498
H     5.7326070    -0.2147977    -0.2854525
H     2.6033664     3.6462550    -0.6086135
H     0.9511283     3.2665112    -0.5646975
--
0 1
S    -1.8220636     2.0964300     0.3299922
C    -2.7962507     0.7575832     0.0980480
N    -2.3768437    -0.5243435     0.0586503
C    -3.1865402    -1.6725337    -0.1260894
C    -4.6014876    -1.3770894    -0.2724329
C    -5.0258914    -0.0988346    -0.2322480
N    -4.1419110     0.9324971    -0.0553359
O    -2.6885486    -2.7802641    -0.1499990
H    -5.2868250    -2.1961027    -0.4114699
H    -1.3618003    -0.7128334     0.1812758
H    -6.0623003     0.1861824    -0.3346801
H    -4.4551541     1.8895240    -0.0237988
units angstrom
}
""")

JSCH_18 = input.process_input("""
molecule dimer {
0 1
C    -1.2382495     0.0003068     3.2761967
N    -0.8377699    -0.0002822     4.6262520
C     0.4580599    -0.0008039     5.0168008
C     1.4462017    -0.0007290     4.0905459
C     1.0380467    -0.0000946     2.7139311
N    -0.2347225     0.0003919     2.3461594
O    -2.4294469     0.0006718     3.0053097
N     1.9638405     0.0000201     1.7458948
H    -1.5882710    -0.0003035     5.2983303
H     0.6465072    -0.0012211     6.0804324
H     2.4837567    -0.0011342     4.3810193
H     2.9358093    -0.0004389     1.9899669
H     1.6736506     0.0004167     0.7589718
--
0 1
N    -1.1590741     0.0004019    -0.4138632
C    -0.2319446     0.0003452    -1.3716397
N     1.0782989     0.0006222    -1.0483779
C     1.9971055     0.0005759    -2.0347184
N     1.8153528     0.0002525    -3.3521684
C     0.5065246    -0.0000627    -3.6438316
C    -0.5584383    -0.0000492    -2.7449616
N    -1.7730910    -0.0004726    -3.3901944
C    -1.4412382    -0.0006880    -4.6700413
N    -0.0894486    -0.0004698    -4.8806889
H     3.0276349     0.0008204    -1.7005582
H    -2.1424713    -0.0010329    -5.4876572
H     0.3894075    -0.0006068    -5.7657319
H    -0.8924030     0.0005975     0.5753416
H    -2.1264189     0.0000818    -0.6840313
units angstrom
}
""")

JSCH_19 = input.process_input("""
molecule dimer {
0 1
O     1.7709955     2.3306811     0.0000007
C     0.5807567     2.6278573    -0.0000017
N    -0.3963320     1.6134139     0.0000053
C    -1.7496702     1.7814670    -0.0000008
N    -2.3431748     2.9556940     0.0000044
C    -1.4426427     3.9663046    -0.0000011
C    -0.0499495     3.9187042     0.0000001
N     0.4958359     5.1819761    -0.0000084
C    -0.5511865     5.9829108     0.0000000
N    -1.7428313     5.3021768     0.0000079
N    -2.4981258     0.6511289    -0.0000042
H    -0.0110714     0.6591960    -0.0000028
H    -0.5111672     7.0591526     0.0000000
H    -2.6703203     5.6919780    -0.0000047
H    -3.4931951     0.7676564     0.0000001
H    -2.1057181    -0.2785664     0.0000054
--
0 1
O    -1.7956163    -2.4184989     0.0000035
C    -0.6750612    -2.9133012    -0.0000022
N    -0.5075308    -4.3172987     0.0000000
C     0.6872157    -4.9931784     0.0000035
N     1.8564548    -4.4063354     0.0000002
C     1.7518235    -3.0547138    -0.0000036
C     0.6020117    -2.2761182     0.0000000
N     0.9066335    -0.9422663    -0.0000018
C     2.2267133    -0.8962763     0.0000026
N     2.7793888    -2.1508756     0.0000002
N     0.6219231    -6.3495411    -0.0000038
H    -1.3767427    -4.8336984     0.0000000
H     2.7978322     0.0182460     0.0000000
H     3.7606420    -2.3771663    -0.0000001
H     1.4878769    -6.8540474     0.0000006
H    -0.2426849    -6.8513669     0.0000027
units angstrom
}
""")

JSCH_20 = input.process_input("""
molecule dimer {
0 1
O    -2.1042101     2.1109877    -0.0000208
C    -0.9861671     2.6187107    -0.0000392
N     0.1630602     1.8075078     0.0000102
C     1.4605433     2.2338425     0.0001555
N     1.8191828     3.4984777    -0.0000542
C     0.7428289     4.3179348    -0.0001135
C    -0.6139590     4.0028703    -0.0000745
N    -1.3931127     5.1362751    -0.0000141
C    -0.5202985     6.1236314    -0.0000063
N     0.7801628     5.6851901    -0.0000822
N     2.4058015     1.2657395     0.0011086
H    -0.0159870     0.7952237     0.0001524
H    -0.7678388     7.1715226     0.0000486
H     1.6159456     6.2449111    -0.0000097
H     3.3631616     1.5602679    -0.0004150
H     2.1837654     0.2818725    -0.0005076
--
0 1
S     2.4306485    -2.2874888    -0.0000855
C     0.9168812    -2.9468359    -0.0000744
N     0.7467080    -4.3320487    -0.0000158
C    -0.4440553    -5.0077630     0.0002108
N    -1.6149750    -4.4178761    -0.0000758
C    -1.5042335    -3.0717080    -0.0001193
C    -0.3421815    -2.2986807    -0.0000424
N    -0.6475784    -0.9636606     0.0000146
C    -1.9693802    -0.9175554    -0.0000159
N    -2.5269735    -2.1675435    -0.0001085
N    -0.3741558    -6.3629805     0.0018230
H     1.6159815    -4.8501428     0.0002176
H    -2.5316778     0.0048746     0.0000297
H    -3.5101913    -2.3858700    -0.0000954
H    -1.2351647    -6.8756289    -0.0006400
H     0.4955304    -6.8577215    -0.0011689
units angstrom
}
""")

JSCH_21 = input.process_input("""
molecule dimer {
0 1
S    -1.8166246     2.6821898    -0.0001323
C    -0.1580956     2.6810732    -0.0000609
N     0.5668809     1.4952672     0.0000612
C     1.9289258     1.3802774     0.0002313
N     2.7529724     2.4104948     0.0000355
C     2.0872014     3.5809186    -0.0000521
C     0.7090604     3.8070941    -0.0000816
N     0.4280062     5.1539248    -0.0000845
C     1.6128743     5.7291564    -0.0000496
N     2.6486647     4.8271187    -0.0000590
N     2.4376947     0.1291202     0.0011300
H     0.0139020     0.6266213     0.0001583
H     1.7867621     6.7921840    -0.0000082
H     3.6347397     5.0250186     0.0000158
H     3.4373976     0.0573910    -0.0001455
H     1.8889030    -0.7194971    -0.0004252
--
0 1
O     1.5845436    -2.6967539    -0.0001530
C     0.4227905    -3.0880579    -0.0000036
N     0.1350039    -4.4687763     0.0000593
C    -1.1124849    -5.0407117    -0.0002399
N    -2.2264983    -4.3563305    -0.0000043
C    -2.0034350    -3.0192765     0.0001096
C    -0.7909434    -2.3397387     0.0000690
N    -0.9816420    -0.9826389     0.0000103
C    -2.2944101    -0.8306632     0.0000066
N    -2.9493710    -2.0333946     0.0000657
N    -1.1503939    -6.3976637    -0.0016596
H     0.9553753    -5.0587791    -0.0001303
H    -2.7848264     0.1304484    -0.0000198
H    -3.9468229    -2.1737881     0.0000098
H    -2.0477574    -6.8436736     0.0002445
H    -0.3206311    -6.9565723     0.0005285
units angstrom
}
""")

JSCH_22 = input.process_input("""
molecule dimer {
0 1
O    -1.3058058     0.3353432    -1.9024452
C    -1.9900049     0.1709800    -0.8920389
N    -1.3797575     0.0147348     0.3602195
C    -2.0188483    -0.2185512     1.5398296
N    -3.3172819    -0.2995486     1.6835149
C    -3.9535777    -0.1111832     0.5011380
C    -3.4163316     0.1070768    -0.7658722
N    -4.4021933     0.2259169    -1.7159019
C    -5.5203687     0.0828308    -1.0307622
N    -5.3061034    -0.1241201     0.3081537
N    -1.1984729    -0.3198503     2.6442461
H    -0.3465494     0.1232217     0.3842921
H    -6.5125630     0.1194762    -1.4473490
H    -5.9980490    -0.2607961     1.0259833
H    -1.6811020    -0.7086455     3.4387722
H    -0.3023253    -0.7473222     2.4686403
--
0 1
N     1.4487686    -0.3061821    -1.8063482
C     2.1291340    -0.0639194    -0.6809352
N     1.4721887     0.2311616     0.4586152
C     2.1772393     0.4822275     1.5830193
N     3.4937302     0.4667260     1.7609342
C     4.1213725     0.1541959     0.6180841
C     3.5365208    -0.1183154    -0.6181769
N     4.4783743    -0.3859600    -1.5825801
C     5.6237470    -0.2779583    -0.9307987
N     5.4697684     0.0435735     0.3900585
H     1.5797459     0.7376023     2.4499547
H     6.5972124    -0.4228222    -1.3684582
H     6.1971892     0.1827155     1.0717315
H     0.4522459    -0.0785075    -1.8709034
H     1.9849455    -0.4365244    -2.6464500
units angstrom
}
""")

JSCH_23 = input.process_input("""
molecule dimer {
0 1
N     0.5317472    -1.5315785     0.0000000
C     1.6654518    -2.2639451     0.0000000
N     1.8178303    -3.5826276     0.0000000
C     0.6256673    -4.1952397     0.0000000
C    -0.6270275    -3.5869658     0.0000000
C    -0.6548527    -2.1760608     0.0000000
N     0.3553662    -5.5402802     0.0000000
C    -1.0067488    -5.6700992     0.0000000
N    -1.6430990    -4.5118780     0.0000000
H     1.0374873    -6.2804981     0.0000000
H    -1.4834180    -6.6357628     0.0000000
H     2.5904811    -1.6988711     0.0000000
N    -1.8018223    -1.4963325     0.0000000
H    -2.6555758    -2.0258909     0.0000000
H    -1.8291436    -0.4726173     0.0000000
--
0 1
C     1.5820983     2.1166821     0.0000000
N     0.3982589     1.4363592     0.0000000
C    -0.8811915     2.0133126     0.0000000
C    -0.7922275     3.4409053     0.0000000
C     0.4748735     4.0188928     0.0000000
N     1.6897724     3.4229374     0.0000000
N     0.2422898     5.3650391     0.0000000
C    -1.1196961     5.5363613     0.0000000
N    -1.7806879     4.3962812     0.0000000
N     2.7046701     1.3533161     0.0000000
H     2.6551620     0.3544770     0.0000000
H     3.5911246     1.8188331     0.0000000
H     0.4300601     0.4044282     0.0000000
O    -1.8846384     1.3003960     0.0000000
H    -1.5681188     6.5152131     0.0000000
H     0.9500064     6.0803260     0.0000000
units angstrom
}
""")

JSCH_24 = input.process_input("""
molecule dimer {
0 1
O    -1.3082180    -0.2837400    -5.2857830
C    -1.2591040    -0.1735270    -4.0737910
N    -0.0039450    -0.3880190    -3.4022040
C     0.2384190    -0.3204430    -2.0528570
N    -0.7115810    -0.0321780    -1.1774330
C    -1.9176870     0.2005520    -1.7572050
C    -2.2714950     0.1569770    -3.1054550
N    -3.6087050     0.4637830    -3.2773380
C    -4.0524060     0.6893350    -2.0695760
N    -3.0715720     0.5411370    -1.1003800
N     1.4926750    -0.5945950    -1.6152550
H     0.7445330    -0.6355620    -4.0377220
H    -5.0658850     0.9619290    -1.8095140
H    -3.1606740     0.6975300    -0.1076790
H     1.7292200    -0.3123020    -0.6555620
H     2.2464200    -0.6005410    -2.2841000
--
0 1
N    -0.6357410    -0.7643850     1.7470590
C     0.2971930    -0.4605580     2.6687480
N    -0.0180070    -0.6125000     3.9730480
C     0.8974790    -0.3069680     4.9024040
N     2.1418480     0.1550240     4.7356170
C     2.4316560     0.2831700     3.4365580
C     1.5994220    -0.0008880     2.3483350
N     2.2633010     0.2340730     1.1506720
C     3.4544220     0.6547710     1.5098070
N     3.6165670     0.7091770     2.8719470
H     0.5799100    -0.4518430     5.9329470
H     4.2477100     0.9349840     0.8305350
H     4.4400710     0.9945670     3.3802140
H    -1.5203500    -1.0951610     2.1001930
H    -0.5364890    -0.5388150     0.7550260
units angstrom
}
""")

JSCH_25 = input.process_input("""
molecule dimer {
0 1
C     1.2592909     1.6400416     0.0000000
N    -0.0329944     1.3939366     0.0000000
C    -0.7724457     2.5351119     0.0000000
C    -0.3568209     3.8625834     0.0000000
C     1.0524521     4.1391914     0.0000000
N     1.7707185     2.9096522     0.0000000
N     2.1434672     0.6224336     0.0000000
N    -1.4258984     4.7253576     0.0000000
C    -2.4770874     3.9323264     0.0000000
N    -2.1383285     2.6014184     0.0000000
H    -2.7704070     1.8193281     0.0000000
H    -3.5031441     4.2568795     0.0000000
O     1.6606912     5.1923670     0.0000000
H     2.7730000     3.0373368     0.0000000
H     1.8138435    -0.3438647     0.0000000
H     3.1276914     0.8060391     0.0000000
--
0 1
C     2.2859985    -3.1747071     0.0000000
N     1.3685098    -2.2195054     0.0000000
C     0.1720555    -2.9042803     0.0000000
N     1.7524294    -4.4267217     0.0000000
C     0.3848788    -4.2845108     0.0000000
C    -1.1754152    -2.4860287     0.0000000
N    -2.1251928    -3.4313144     0.0000000
C    -1.7646506    -4.7253100     0.0000000
N    -0.5383069    -5.2487516     0.0000000
H     3.3496602    -3.0075287     0.0000000
H     2.2521180    -5.3008143     0.0000000
H    -2.5835604    -5.4328271     0.0000000
N    -1.5512019    -1.1969440     0.0000000
H    -0.8988350    -0.4160731     0.0000000
H    -2.5417242    -1.0304237     0.0000000
units angstrom
}
""")

JSCH_26 = input.process_input("""
molecule dimer {
0 1
O     1.0272885    -1.7927509    -0.4061508
C     1.7699883    -0.8358415    -0.1967480
N     1.2319553     0.4331963     0.0618047
C     1.9364691     1.5620604     0.3460907
N     3.2411330     1.6380980     0.3913752
C     3.8110938     0.4396217     0.1071531
C     3.2032408    -0.7831141    -0.1709853
N     4.1349400    -1.7700995    -0.3829089
C     5.2907910    -1.1513051    -0.2354273
N     5.1516924     0.1802511     0.0617335
N     1.1709724     2.6978601     0.5477231
H     0.2039336     0.5207236    -0.0280886
H     6.2587500    -1.6129704    -0.3315009
H     5.8835828     0.8525276     0.2204394
H     1.7146455     3.4340823     0.9714033
H     0.3029369     2.5238305     1.0323857
--
0 1
N    -1.6634540    -2.1503266     0.4844345
C    -2.7337243    -1.3645844     0.2851342
N    -3.9617362    -1.8916536     0.4117484
C    -5.0361049    -1.1108105     0.2142857
N    -5.0856134     0.1814610    -0.1169411
C    -3.8521385     0.6813689    -0.2227871
C    -2.6434417     0.0104508    -0.0286100
N    -1.5757669     0.8650951    -0.2045833
C    -2.1325255     2.0285854    -0.5143055
N    -3.4907538     1.9716246    -0.5383952
H    -5.9917639    -1.6042053     0.3341775
H    -1.5835904     2.9296030    -0.7337374
H    -4.1262380     2.7236162    -0.7489234
H    -1.8639094    -3.1302901     0.5932195
H    -0.7358064    -1.8946489     0.1475504
units angstrom
}
""")

JSCH_27 = input.process_input("""
molecule dimer {
0 1
O     5.3545637    -1.5839084     0.1643820
C     4.3016967    -0.9781139     0.1123080
N     3.0713646    -1.6809667     0.2659334
C     1.8042754    -1.1676035     0.2494339
N     1.5592029     0.1129485     0.0745442
C     2.6978871     0.8367062    -0.0926982
C     4.0240233     0.4175244    -0.0912729
N     4.8841613     1.4681703    -0.3008504
C     4.0879153     2.5107452    -0.4295656
N     2.7596551     2.1849030    -0.3099959
N     0.7798744    -2.0309505     0.4761070
H     3.1965102    -2.6688068     0.4391670
H     4.4098825     3.5222321    -0.6091205
H     1.9722362     2.8053438    -0.3963660
H    -0.1527293    -1.7022269     0.2002243
H     0.9606535    -3.0046040     0.3028696
--
0 1
N    -1.2242031     1.0428672     0.4916050
C    -2.2040220     0.1989638     0.1394948
N    -1.9060140    -1.0547654    -0.2435106
C    -2.9084418    -1.8851915    -0.5992935
N    -4.2169079    -1.6502337    -0.6282927
C    -4.4819036    -0.3976671    -0.2338853
C    -3.5616520     0.5734472     0.1583798
N    -4.1759166     1.7562758     0.4961983
C    -5.4596333     1.4978989     0.3100336
N    -5.7008222     0.2244483    -0.1272061
H    -2.5970464    -2.8772539    -0.9009149
H    -6.2592754     2.1999750     0.4769792
H    -6.5949340    -0.1885904    -0.3349675
H    -0.2461099     0.7667993     0.3993852
H    -1.4852074     1.9635263     0.7951137
units angstrom
}
""")

JSCH_28 = input.process_input("""
molecule dimer {
0 1
N    -1.2744921    -0.0017953    -1.3782659
C    -2.1152770    -0.0015179    -0.3381875
N    -1.6375936    -0.0033248     0.9198720
C    -2.5074052    -0.0030670     1.9505655
N    -3.8371522    -0.0012513     1.9239415
C    -4.2828181     0.0005849     0.6601536
C    -3.5155952     0.0006434    -0.5034237
N    -4.2984407     0.0028211    -1.6341051
C    -5.5305121     0.0040931    -1.1536448
N    -5.5811183     0.0027919     0.2133605
H    -2.0528877    -0.0045852     2.9334633
H    -6.4253358     0.0059818    -1.7533163
H    -6.4045390     0.0033962     0.7920186
H    -0.2596925    -0.0029765    -1.2406098
H    -1.6728767    -0.0000588    -2.3000471
--
0 1
N     1.2734087    -0.0017991     1.3765409
C     2.1149782    -0.0015194     0.3371102
N     1.6382658    -0.0033244    -0.9213222
C     2.5089058    -0.0030642    -1.9513246
N     3.8386381    -0.0012477    -1.9236383
C     4.2833025     0.0005866    -0.6594937
C     3.5151554     0.0006426     0.5034653
N     4.2970793     0.0028188     1.6347785
C     5.5295364     0.0040924     1.1553161
N     5.5812428     0.0027936    -0.2116508
H     2.0551790    -0.0045813    -2.9345889
H     6.4238744     0.0059806     1.7557114
H     6.4051319     0.0033994    -0.7896412
H     0.2587211    -0.0029792     1.2389462
H     1.6714569    -0.0000645     2.2984633
units angstrom
}
""")

JSCH_29 = input.process_input("""
molecule dimer {
0 1
N    -1.1366363    -0.0666108    -1.4922977
C    -1.9363007    -0.0408578    -0.4194776
N    -1.4071800    -0.0593819     0.8158785
C    -2.2352343    -0.0402428     1.8802964
N    -3.5645999    -0.0046868     1.9074668
C    -4.0608616     0.0128417     0.6631804
C    -3.3407627    -0.0019209    -0.5302183
N    -4.1677258     0.0274614    -1.6288352
C    -5.3791670     0.0592559    -1.0991816
N    -5.3753601     0.0520242     0.2686301
H    -1.7428602    -0.0568025     2.8448773
H    -6.2968580     0.0879838    -1.6624329
H    -6.1747025     0.0724119     0.8798069
H    -0.1211878    -0.0517077    -1.3879919
H    -1.5677584    -0.0159058    -2.3974751
--
0 1
N     1.8123343    -0.0245408    -1.2588738
C     2.7152804    -0.0039456    -0.2161678
C     2.5630991    -0.0007999     1.1886196
N     3.6772242     0.0214724     1.9383141
C     4.8749026     0.0394397     1.3318824
N     5.1551234     0.0392981     0.0273056
C     4.0284465     0.0171443    -0.6906539
N     3.9051515     0.0090860    -2.0602851
C     2.5735625    -0.0163877    -2.3426881
N     1.3736145    -0.0168867     1.8043279
H     5.7275827     0.0569289     1.9982511
H     2.2043398    -0.0288835    -3.3545151
H     4.6669159     0.0187011    -2.7180078
H     1.3935693    -0.0168974     2.8088363
H     0.4799202    -0.0417563     1.3125324
units angstrom
}
""")

JSCH_30 = input.process_input("""
molecule dimer {
0 1
N     1.9051383    -0.1221668     1.3018987
C     1.2306207    -0.0724496     2.4599521
N     1.9188011    -0.0911776     3.6114745
C     1.2468662    -0.0458415     4.7741848
N    -0.0695564     0.0217460     4.9796590
C    -0.7242768     0.0390915     3.8148784
C    -0.1771642    -0.0054840     2.5309668
N    -1.1622374     0.0280993     1.5676491
C    -2.2864196     0.0935989     2.2653565
N    -2.0826685     0.1041226     3.6120321
H     1.8637884    -0.0652578     5.6631315
H    -3.2746103     0.1375165     1.8384429
H    -2.7831685     0.1478311     4.3334811
H     2.9067312    -0.1495396     1.3694272
H     1.4544353    -0.0668897     0.3920587
--
0 1
N    -1.9061558    -0.0641839    -1.3006378
C    -1.2309020    -0.0492565    -2.4591003
N    -1.9184496    -0.0947841    -3.6103362
C    -1.2461054    -0.0720254    -4.7734664
N     0.0702056    -0.0055843    -4.9798973
C     0.7244044     0.0360165    -3.8154103
C     0.1769597     0.0162410    -2.5310350
N     1.1617828     0.0675370    -1.5681922
C     2.2861580     0.1193292    -2.2667167
N     2.0827966     0.1037780    -3.6134111
H    -1.8625819    -0.1122096    -5.6620263
H     3.2742511     0.1708248    -1.8404275
H     2.7834317     0.1345490    -4.3353905
H    -2.9057988    -0.1329805    -1.3671674
H    -1.4533881    -0.0384046    -0.3902022
units angstrom
}
""")

JSCH_31 = input.process_input("""
molecule dimer {
0 1
C    -5.2998476     2.1696769    -0.1527418
N    -4.0012332     2.3280352     0.0113365
C    -3.5133825     1.0430020     0.0349501
N    -5.6731934     0.8519995    -0.2385098
C    -4.5348158     0.1058659    -0.1189924
N    -2.1562547    -0.8394339     0.1532486
C    -3.2316395    -1.6611847    -0.0119881
N    -4.4661255    -1.2464671    -0.1599515
N    -2.9527518    -3.0009490     0.0413727
H    -6.0192845     2.9682219    -0.2165198
H    -2.0059608    -3.2718138    -0.1741253
H    -3.6766313    -3.5832294    -0.3425488
H    -6.6016641     0.4841764    -0.3616632
H    -1.2330572    -1.2973531     0.2284474
C    -2.1744618     0.5575368     0.1742023
O    -1.1240037     1.1950526     0.2956190
--
0 1
N     1.2638556    -0.1422038     0.1193086
C     1.2024891    -1.5065889     0.0798953
N     2.5150609    -1.9393418    -0.0155197
C     3.3711631    -0.8652191    -0.0318660
C     2.5823266     0.2608198     0.0542710
O     0.2004187    -2.2363595     0.1229392
C     3.1762467     1.5568980     0.0478905
N     4.5883289     1.4099390    -0.0474268
C     5.2845121     0.2404436    -0.1157151
N     4.7218815    -0.9391943    -0.1032764
O     2.6645381     2.6620253     0.1017547
N     6.6513076     0.3472133    -0.2677366
H     7.0896806     1.0869563     0.2575269
H     7.1055907    -0.5431980    -0.1372108
H     2.7741742    -2.9100741    -0.0564635
H     0.4151138     0.4528548     0.1972152
H     5.0828171     2.2880907    -0.1283728
units angstrom
}
""")

JSCH_32 = input.process_input("""
molecule dimer {
0 1
S    -0.2983354    -0.0000513     0.0606545
C    -0.2090863    -0.0000888     1.7027085
N     0.9916329    -0.0001803     2.3727915
C     1.1063707     0.0000024     3.7325291
C     0.0163475     0.0001662     4.5304894
C    -1.2953719     0.0000721     3.9249070
N    -1.2941374    -0.0000268     2.5353745
O    -2.3533146     0.0000820     4.5510728
H     0.0943097     0.0003348     5.6041112
H     1.8067829    -0.0002554     1.7802008
H    -2.2194333     0.0000304     2.0853365
H     2.1154718     0.0000158     4.1172107
--
0 1
S    -4.3480040     0.0005221     1.2455679
C    -5.4129697     0.0002108     2.5234518
N    -6.7626348     0.0001784     2.2970286
C    -7.6987363     0.0000487     3.2957366
C    -7.3354268    -0.0001000     4.5945383
C    -5.9267360    -0.0001870     4.9466675
N    -5.0628029     0.0000198     3.8318987
O    -5.4752372    -0.0004103     6.0770807
H    -8.0659143    -0.0001724     5.3856633
H    -7.0337303     0.0003183     1.3267291
H    -4.0549036     0.0000246     4.0487228
H    -8.7287725     0.0000774     2.9713072
units angstrom
}
""")

JSCH_33 = input.process_input("""
molecule dimer {
0 1
C    12.1619966    21.5469940    -0.5249999
N    12.0019966    20.1249944    -0.3349999
C    12.9959964    19.1989946    -0.1290000
N    12.5899965    17.9429950    -0.1260000
C    11.2289969    18.0629949    -0.3469999
C    10.2259971    17.0909952    -0.4599999
N    10.4079971    15.7719956    -0.3739999
N     8.9619975    17.5199951    -0.6819998
C     8.7349976    18.8509947    -0.7899998
N     9.6049973    19.8469944    -0.7019998
C    10.8559970    19.3909946    -0.4999999
H    12.8450824    21.9515608     0.2257099
H    12.5490085    21.7744749    -1.5236356
H    11.1843859    22.0177918    -0.4120399
H    14.0220821    19.5129525     0.0161520
H    11.3436468    15.4109067    -0.2800629
H     9.6382753    15.1406078    -0.5991948
H     7.6909448    19.1156876    -0.9420537
--
0 1
C     3.0629991    16.2869954    -0.5529998
N     4.3679988    15.6949956    -0.7379998
C     5.4889985    16.5069954    -0.6549998
O     5.3979985    17.7169950    -0.4679999
N     6.6749981    15.8589956    -0.7949998
C     6.8699981    14.5069959    -0.9999997
O     8.0199978    14.0679961    -1.0789997
C     5.6559984    13.7139962    -1.1019997
C     5.7709984    12.2569966    -1.4029996
C     4.4739987    14.3319960    -0.9639997
H     7.5313379    16.4637704    -0.7443448
H     6.3741672    11.7424167    -0.6472968
H     4.7881707    11.7797217    -1.4448876
H     6.2751442    12.0930036    -2.3618343
H     3.5293140    13.8026561    -1.0289747
H     2.3790703    15.9479585    -1.3364316
H     2.6423583    16.0249025     0.4245489
H     3.1730521    17.3682771    -0.6086068
units angstrom
}
""")

JSCH_34 = input.process_input("""
molecule dimer {
0 1
N    10.3469971    14.4959959     8.8169975
C    11.5789968    13.8469961     8.7069976
O    11.6019967    12.6419965     8.4119976
N    12.6939964    14.5549959     8.8809975
C    12.6739964    15.9259955     9.1859974
N    13.8309961    16.5099954     9.3349974
C    11.4219968    16.5639954     9.2669974
C    10.3209971    15.8539956     9.0929975
H     9.3699974    16.4009954     9.1789974
H    11.3019968    17.6379951     9.4699973
H    14.6739959    15.9769955     9.2609974
H    13.8749961    17.4909951     9.5239973
C     9.1059774    13.7460371     8.6280336
H     9.4001314    12.7260934     8.3864956
H     8.5051816    13.7537151     9.5428113
H     8.5206636    14.1698120     7.8064238
--
0 1
C    18.8919947     9.6579973     9.7709973
N    18.5279948    11.0699969     9.5879973
C    19.3769946    12.1419966     9.6129973
N    18.7759947    13.3089963     9.4319974
C    17.4529951    12.9639964     9.3169974
C    16.2779954    13.7529961     9.1209974
O    16.2219955    14.9839958     9.0219975
N    15.1359958    13.0409963     9.0449975
C    15.0849958    11.6719967     9.1349974
N    13.8449961    11.1639969     9.0359975
N    16.1359955    10.8809970     9.3169974
C    17.2759952    11.5909968     9.3939974
H    14.2561290    13.5779002     8.9264415
H    13.0353973    11.7259537     8.7509445
H    13.7773141    10.1594092     9.0213535
H    17.9866610     9.0649795     9.6385253
H    19.2909706     9.4904943    10.7753660
H    19.6360815     9.3587324     9.0282525
H    20.4431063    12.0114766     9.7460263
units angstrom
}
""")

JSCH_35 = input.process_input("""
molecule dimer {
0 1
N    10.9240000    16.7550000     5.5620000
C    11.6470000    17.8510000     5.8140000
N    12.9490000    17.6590000     5.9790000
C    13.0500000    16.2780000     5.7950000
C    14.1950000    15.4230000     5.8560000
N    15.4060000    15.8590000     6.0610000
N    13.9020000    14.1180000     5.6250000
C    12.6770000    13.6430000     5.3990000
N    11.5490000    14.4040000     5.3300000
C    11.8450000    15.6910000     5.5460000
H    11.1804230    18.8265530     5.8822870
H    12.5884030    12.5696370     5.2620740
H    16.1977530    15.2199420     5.9750360
H    15.5570940    16.8510580     6.1500010
C     9.4931860    16.6413650     5.3399050
H     9.0446590    17.6337380     5.4112840
H     9.2947180    16.2234190     4.3499330
H     9.0442270    15.9854440     6.0897950
--
0 1
N    16.2460000     9.7810000     5.9650000
C    17.5950000    10.0510000     5.9930000
C    18.0920000    11.2690000     5.9020000
C    17.1390000    12.3410000     5.7640000
O    17.4920000    13.5330000     5.6630000
N    15.8280000    12.0550000     5.7130000
C    15.3100000    10.7970000     5.7960000
O    14.1120000    10.5770000     5.7580000
H    18.2280000     9.1744860     6.1031120
C    19.5529600    11.6051630     5.9357380
H    20.1631860    10.7042230     6.0438290
H    19.7760320    12.2828240     6.7658180
H    19.8526100    12.1260780     5.0209680
H    15.1383860    12.8499570     5.6472680
C    15.7717470     8.4029560     6.0779300
H    14.6864640     8.4223240     6.0045990
H    16.1825380     7.7884380     5.2708940
H    16.0652090     7.9755790     7.0417370
units angstrom
}
""")

JSCH_36 = input.process_input("""
molecule dimer {
0 1
H     0.0112670     4.2441280     0.3057270
N    -0.1600000     4.2010000     1.2990000
C     0.1490000     5.1520000     2.2350000
H     0.8336150     5.9557770     2.0023890
N    -0.3040000     4.9000000     3.4380000
C    -1.1470000     3.7970000     3.2290000
C    -2.0790000     3.1160000     4.0900000
O    -2.3440000     3.3110000     5.2740000
N    -2.7730000     2.0930000     3.4630000
H    -3.4444620     1.6202680     4.0533010
C    -2.5700000     1.7190000     2.1650000
N    -3.2200000     0.6740000     1.7040000
H    -3.7884800     0.1079360     2.3113460
H    -3.0424470     0.3264300     0.7529310
N    -1.7100000     2.3160000     1.3470000
C    -1.0480000     3.3630000     1.9240000
--
0 1
H    -3.4958570    -1.4150050    -3.9137580
N    -3.0510000    -1.0010000    -3.1090000
C    -3.5590000    -0.8800000    -1.8360000
H    -4.5790060    -1.1582720    -1.6128580
N    -2.7220000    -0.3740000    -0.9680000
C    -1.5590000    -0.1810000    -1.7250000
C    -0.2720000     0.3480000    -1.4650000
N     0.1070000     0.8840000    -0.3230000
H     1.0433330     1.2579620    -0.3065570
H    -0.5751070     1.2407790     0.3499520
N     0.6670000     0.3750000    -2.4130000
C     0.3480000    -0.0810000    -3.6160000
H     1.1321870    -0.0417550    -4.3673920
N    -0.8160000    -0.5790000    -4.0190000
C    -1.7380000    -0.6050000    -3.0150000
units angstrom
}
""")

JSCH_37 = input.process_input("""
molecule dimer {
0 1
H     3.1762460     2.3738070     2.9634160
N     2.3770000     1.8470000     3.2830000
C     1.6370000     2.2160000     4.3790000
H     1.9902970     3.0843050     4.9210710
C     0.5610000     1.4930000     4.7730000
H    -0.0085000     1.7736330     5.6470440
C     0.1830000     0.3990000     3.9430000
N    -0.8510000    -0.3400000     4.2540000
H    -1.1799330    -1.0651510     3.5908230
H    -1.4362750    -0.1022370     5.0377650
N     0.8500000     0.0580000     2.8540000
C     1.9550000     0.7640000     2.4990000
O     2.5580000     0.4150000     1.4830000
--
0 1
H    -1.2611710    -4.7286740    -2.6257100
N    -1.6090000    -4.2940000    -1.7860000
C    -2.7550000    -4.5990000    -1.0690000
H    -3.5136190    -5.2427470    -1.4922410
N    -2.8650000    -3.9860000     0.0730000
C    -1.6740000    -3.2820000     0.1910000
C    -1.1780000    -2.4570000     1.2560000
O    -1.7150000    -2.1460000     2.3170000
N     0.0980000    -1.9830000     1.0200000
H     0.4562670    -1.3045040     1.7132710
C     0.8280000    -2.2730000    -0.0890000
N     2.0180000    -1.7250000    -0.1770000
H     2.3044660    -0.9690820     0.4476800
H     2.5064670    -1.8555350    -1.0472790
N     0.3920000    -3.0250000    -1.1030000
C    -0.8790000    -3.5010000    -0.9150000
units angstrom
}
""")

JSCH_38 = input.process_input("""
molecule dimer {
0 1
H     4.0780890     0.2050200     6.5267380
N     3.3380000    -0.4520000     6.3380000
C     2.1440000    -0.6140000     7.0100000
H     1.9445960    -0.0744500     7.9251340
N     1.3390000    -1.4880000     6.4770000
C     2.0190000    -1.9110000     5.3320000
C     1.6500000    -2.8430000     4.3020000
O     0.6370000    -3.5330000     4.1980000
N     2.5960000    -2.9520000     3.3010000
H     2.3705000    -3.6388980     2.5623150
C     3.7610000    -2.2490000     3.2730000
N     4.5620000    -2.4690000     2.2580000
H     4.3528370    -3.1696290     1.5459440
H     5.4428290    -1.9835850     2.2550440
N     4.1450000    -1.3880000     4.2160000
C     3.2280000    -1.2560000     5.2240000
--
0 1
H     3.2823840    -6.1134940    -1.3105350
N     2.5530000    -6.0070000    -0.6210000
C     1.3990000    -6.7620000    -0.6490000
H     1.3017290    -7.4646550    -1.4662410
C     0.4550000    -6.5890000     0.3070000
H    -0.4593850    -7.1648600     0.2947650
C     0.7210000    -5.6290000     1.3280000
N    -0.1590000    -5.3940000     2.2700000
H    -1.0266130    -5.9017830     2.3125200
H     0.0709100    -4.7127400     3.0149280
N     1.8460000    -4.9310000     1.3860000
C     2.7800000    -5.0940000     0.4140000
O     3.8210000    -4.4400000     0.4780000
units angstrom
}
""")

JSCH_39 = input.process_input("""
molecule dimer {
0 1
O     0.9601320     1.3436400     0.0000000
C     1.5166980     0.2684520     0.0000000
N     0.7573320    -0.9011610     0.0000000
C     1.2481620    -2.1702510     0.0000000
N     2.5209460    -2.4496950     0.0000000
C     3.2915230    -1.3476830     0.0000000
C     2.9121790    -0.0279190     0.0000000
N     4.0200060     0.7969640     0.0000000
C     5.0170310     0.0003310     0.0000000
N     4.6446780    -1.3255770     0.0000000
N     0.3459700    -3.1553460     0.0000000
H    -0.2412520    -0.7659240     0.0000000
H     6.0483360     0.2895830     0.0000000
H     5.2362800    -2.1226110     0.0000000
H     0.6928700    -4.0838600     0.0000000
H    -0.6408270    -2.9885130     0.0000000
--
0 1
C    -1.5982280    -2.9490360     3.3600000
N    -2.8308990    -3.5868360     3.3600000
C    -4.0005400    -2.9065270     3.3600000
C    -4.0107280    -1.5698660     3.3600000
C    -2.7192980    -0.9187180     3.3600000
N    -1.5949260    -1.5998660     3.3600000
O    -0.5980710    -3.6295230     3.3600000
N    -2.6531990     0.4024280     3.3600000
H    -2.8066410    -4.5810390     3.3600000
H    -4.8972920    -3.4971900     3.3600000
H    -4.9235800    -1.0089750     3.3600000
H    -3.4794940     0.9500750     3.3600000
H    -1.7581040     0.8646590     3.3600000
units angstrom
}
""")

JSCH_40 = input.process_input("""
molecule dimer {
0 1
C    -3.0263940    -1.4464050     0.0000000
N    -4.3985350    -1.2378500     0.0000000
C    -4.9449180     0.0000290     0.0000000
C    -4.1674910     1.0873990     0.0000000
C    -2.7399670     0.8551050     0.0000000
N    -2.2307000    -0.3568440     0.0000000
O    -2.6172300    -2.5848080     0.0000000
N    -1.9099420     1.8850830     0.0000000
H    -4.9632880    -2.0564360     0.0000000
H    -6.0175890     0.0492700     0.0000000
H    -4.5763200     2.0777300     0.0000000
H    -2.2565290     2.8138220     0.0000000
H    -0.9141020     1.7329110     0.0000000
--
0 1
O    -0.0130090     1.6513790     3.3600000
C     1.0692420     1.1086750     3.3600000
N     1.1423840    -0.2839060     3.3600000
C     2.2854260    -1.0221180     3.3600000
N     3.4793830    -0.5000700     3.3600000
C     3.4550460     0.8444100     3.3600000
C     2.3724120     1.6891490     3.3600000
N     2.7838090     3.0076580     3.3600000
C     4.0586690     2.9492050     3.3600000
N     4.5367780     1.6576590     3.3600000
N     2.1345620    -2.3493720     3.3600000
H     0.2550220    -0.7614500     3.3600000
H     4.7229940     3.7894010     3.3600000
H     5.4838790     1.3605800     3.3600000
H     2.9609760    -2.8966530     3.3600000
H     1.2381640    -2.7944260     3.3600000
units angstrom
}
""")

JSCH_41 = input.process_input("""
molecule dimer {
0 1
N    -1.3923840    -1.5825730    -0.2790500
C    -1.8533500    -0.3518640    -0.0620430
N    -0.9943890     0.6521290     0.1149880
C    -1.4604570     1.8814980     0.3317590
N    -2.7070820     2.2763020     0.4013740
C    -3.5527210     1.2640760     0.2228910
C    -3.2236500    -0.0504790    -0.0089010
N    -4.3580740    -0.8272780    -0.1458710
C    -5.3247240    -0.0009840    -0.0001730
N    -4.9130980     1.2870000     0.2269330
H    -0.7040060     2.6348130     0.4645890
H    -6.3651290    -0.2529400    -0.0446000
H    -5.4840420     2.0871050     0.3680130
H    -0.4093220    -1.7576030    -0.3099130
H    -2.0356960    -2.3259680    -0.4101310
--
0 1
O     2.4555320    -0.5209070     3.3788050
C     2.5333330     0.6704300     3.2169230
N     1.4067200     1.4246400     2.9925690
C     1.3497150     2.7756270     2.7939400
N     2.5708460     3.3948650     2.8321660
C     3.7496420     2.7230470     3.0501750
C     3.8036770     1.4072660     3.2434740
O     0.3307920     3.3742620     2.6029400
C     5.0685490     0.6342580     3.4848390
H     2.5588020     4.3783590     2.6906210
H     0.5200190     0.9342720     2.9706210
H     4.6288470     3.3398640     3.0533080
H     5.0316430     0.1233820     4.4405330
H     5.2089370    -0.1230850     2.7216660
H     5.9296670     1.2931580     3.4791940
units angstrom
}
""")

JSCH_42 = input.process_input("""
molecule dimer {
0 1
O     1.6803850    -1.8647480     0.3288050
C     2.4435780    -0.9466670     0.1669230
N     1.9754430     0.3257080    -0.0574310
C     2.7234150     1.4521870    -0.2560600
N     4.0753100     1.2353980    -0.2178340
C     4.6340910    -0.0009940     0.0001750
C     3.9044100    -1.0972430     0.1934740
O     2.2509580     2.5354000    -0.4470600
C     4.4733490    -2.4660930     0.4348390
H     4.6436490     2.0381400    -0.3593790
H     0.9698550     0.4501820    -0.0793790
H     5.7079390    -0.0187620     0.0033080
H     4.1432070    -2.8577080     1.3905330
H     4.1417710    -3.1613140    -0.3283340
H     5.5573000    -2.4391850     0.4291940
--
0 1
N    -0.1962490    -2.0987510     2.7709500
C    -1.2925710    -1.3740360     2.9879570
N    -1.1877890    -0.0569040     3.1649880
C    -2.2874510     0.6637280     3.3817590
N    -3.5280520     0.2503840     3.4513740
C    -3.6172170    -1.0655780     3.2728910
C    -2.5783160    -1.9356530     3.0410990
N    -3.0394940    -3.2308940     2.9041290
C    -4.3072140    -3.1305910     3.0498260
N    -4.7312590    -1.8466420     3.2769330
H    -2.1182570     1.7178040     3.5145890
H    -5.0008230    -3.9459620     3.0054000
H    -5.6634530    -1.5349360     3.4180130
H     0.7019450    -1.6625240     2.7400870
H    -0.2797430    -3.0783000     2.6398690
units angstrom
}
""")

JSCH_43 = input.process_input("""
molecule dimer {
0 1
C     2.4313070     1.6249990    -1.4530130
N     3.8007370     1.6249990    -1.6786800
C     4.7029040     1.6249990    -0.6702290
C     4.2995430     1.6249990     0.6041590
C     2.8701050     1.6249990     0.8243640
N     2.0112500     1.6249990    -0.1708960
O     1.6903830     1.6249990    -2.4092590
N     2.3989850     1.6249990     2.0604220
H     4.0848920     1.6249990    -2.6317200
H     5.7382910     1.6249990    -0.9548720
H     4.9943920     1.6249990     1.4196840
H     3.0156050     1.6249990     2.8366040
H     1.4048620     1.6249990     2.2234300
--
0 1
C    -2.4313070    -1.6249990    -1.4530130
N    -3.8007370    -1.6249990    -1.6786800
C    -4.7029040    -1.6249990    -0.6702290
C    -4.2995430    -1.6249990     0.6041590
C    -2.8701050    -1.6249990     0.8243640
N    -2.0112500    -1.6249990    -0.1708960
O    -1.6903830    -1.6249990    -2.4092590
N    -2.3989850    -1.6249990     2.0604220
H    -4.0848920    -1.6249990    -2.6317200
H    -5.7382910    -1.6249990    -0.9548720
H    -4.9943920    -1.6249990     1.4196840
H    -3.0156050    -1.6249990     2.8366040
H    -1.4048620    -1.6249990     2.2234300
units angstrom
}
""")

JSCH_44 = input.process_input("""
molecule dimer {
0 1
O    -0.4979320     1.6249990     1.9422390
C    -1.3595090     1.6249990     1.0916630
N    -0.9987390     1.6249990    -0.2553610
C    -1.8577170     1.6249990    -1.3106620
N    -3.1545590     1.6249990    -1.1831180
C    -3.5468800     1.6249990     0.1030790
C    -2.7782730     1.6249990     1.2410250
N    -3.5769760     1.6249990     2.3678730
C    -4.7713760     1.6249990     1.9183280
N    -4.8269760     1.6249990     0.5422510
N    -1.3040920     1.6249990    -2.5263360
H    -0.0072390     1.6249990    -0.4353230
H    -5.6628210     1.6249990     2.5121130
H    -5.6359190     1.6249990    -0.0329580
H    -1.9209400     1.6249990    -3.3022070
H    -0.3140390     1.6249990    -2.6726050
--
0 1
O     0.4979320    -1.6249990     1.9422390
C     1.3595090    -1.6249990     1.0916630
N     0.9987390    -1.6249990    -0.2553610
C     1.8577170    -1.6249990    -1.3106620
N     3.1545590    -1.6249990    -1.1831180
C     3.5468800    -1.6249990     0.1030790
C     2.7782730    -1.6249990     1.2410250
N     3.5769760    -1.6249990     2.3678730
C     4.7713760    -1.6249990     1.9183280
N     4.8269760    -1.6249990     0.5422510
N     1.3040920    -1.6249990    -2.5263360
H     0.0072390    -1.6249990    -0.4353230
H     5.6628210    -1.6249990     2.5121130
H     5.6359190    -1.6249990    -0.0329580
H     1.9209400    -1.6249990    -3.3022070
H     0.3140390    -1.6249990    -2.6726050
units angstrom
}
""")

JSCH_45 = input.process_input("""
molecule dimer {
0 1
O     0.9601320     1.3436400     0.0000000
C     1.5166980     0.2684520     0.0000000
N     0.7573320    -0.9011610     0.0000000
C     1.2481620    -2.1702510     0.0000000
N     2.5209460    -2.4496950     0.0000000
C     3.2915230    -1.3476830     0.0000000
C     2.9121790    -0.0279190     0.0000000
N     4.0200060     0.7969640     0.0000000
C     5.0170310     0.0003310     0.0000000
N     4.6446780    -1.3255770     0.0000000
N     0.3459700    -3.1553460     0.0000000
H    -0.2412520    -0.7659240     0.0000000
H     6.0483360     0.2895830     0.0000000
H     5.2362800    -2.1226110     0.0000000
H     0.6928700    -4.0838600     0.0000000
H    -0.6408270    -2.9885130     0.0000000
--
0 1
O    -1.5665350     0.5226760     3.1900000
C    -1.3848270    -0.6743110     3.1900000
N    -0.0830050    -1.1742030     3.1900000
C     0.2658570    -2.4894210     3.1900000
N    -0.5995930    -3.4636200     3.1900000
C    -1.8707500    -3.0250070     3.1900000
C    -2.3395920    -1.7343230     3.1900000
N    -3.7206970    -1.7181430     3.1900000
C    -4.0590580    -2.9486690     3.1900000
N    -2.9784690    -3.8024880     3.1900000
N     1.5747700    -2.7560840     3.1900000
H     0.6453760    -0.4778410     3.1900000
H    -5.0634190    -3.3208450     3.1900000
H    -2.9886000    -4.7950370     3.1900000
H     1.8398890    -3.7111710     3.1900000
H     2.2750440    -2.0410890     3.1900000
units angstrom
}
""")

JSCH_46 = input.process_input("""
molecule dimer {
0 1
C    -3.0263940    -1.4464050     0.0000000
N    -4.3985350    -1.2378500     0.0000000
C    -4.9449180     0.0000290     0.0000000
C    -4.1674910     1.0873990     0.0000000
C    -2.7399670     0.8551050     0.0000000
N    -2.2307000    -0.3568440     0.0000000
O    -2.6172300    -2.5848080     0.0000000
N    -1.9099420     1.8850830     0.0000000
H    -4.9632880    -2.0564360     0.0000000
H    -6.0175890     0.0492700     0.0000000
H    -4.5763200     2.0777300     0.0000000
H    -2.2565290     2.8138220     0.0000000
H    -0.9141020     1.7329110     0.0000000
--
0 1
C     3.2985790     0.6087040     3.1900000
N     4.2860790     1.5839520     3.1900000
C     4.0005050     2.9065740     3.1900000
C     2.7324140     3.3293140     3.1900000
C     1.7140620     2.3023070     3.1900000
N     2.0144220     1.0224800     3.1900000
O     3.6366950    -0.5527840     3.1900000
N     0.4371510     2.6477000     3.1900000
H     5.2241270     1.2536560     3.1900000
H     4.8393710     3.5769110     3.1900000
H     2.4810610     4.3708130     3.1900000
H     0.1716470     3.6027840     3.1900000
H    -0.2790560     1.9392500     3.1900000
units angstrom
}
""")

JSCH_47 = input.process_input("""
molecule dimer {
0 1
N     1.0423840    -1.6008720     0.1400580
C     1.5033500    -0.3559320     0.0311400
N     0.6443890     0.6596690    -0.0577140
C     1.1104570     1.9032530    -0.1665130
N     2.3570820     2.3026230    -0.2014530
C     3.2027210     1.2786930    -0.1118710
C     2.8736500    -0.0510630     0.0044670
N     4.0080740    -0.8368430     0.0732140
C     4.9747240    -0.0009950     0.0000870
N     4.5630980     1.3018810    -0.1139000
H     0.3540060     2.6652780    -0.2331820
H     6.0151290    -0.2558650     0.0223850
H     5.1340420     2.1112380    -0.1847090
H     0.0593220    -1.7779260     0.1555480
H     1.6856960    -2.3528630     0.2058490
--
0 1
C    -1.6419140     2.9739730    -3.0239370
N    -2.8741190     3.6124140    -3.0421140
C    -4.0409900     2.9359160    -3.1500030
C    -4.0487470     1.6026030    -3.2447730
C    -2.7578360     0.9507400    -3.2245270
N    -1.6361750     1.6281560    -3.1188990
O    -0.6443040     3.6509540    -2.9247190
N    -2.6894340    -0.3672360    -3.3142960
H    -2.8516920     4.6040980    -2.9707700
H    -4.9376330     3.5267310    -3.1542940
H    -4.9593830     1.0447610    -3.3310860
H    -3.5136510    -0.9120230    -3.3952410
H    -1.7946790    -0.8299350    -3.3010330
units angstrom
}
""")

JSCH_48 = input.process_input("""
molecule dimer {
0 1
O    -2.0303850    -1.8863100    -0.1650310
C    -2.7935780    -0.9576130    -0.0837800
N    -2.3254430     0.3294740     0.0288250
C    -3.0734150     1.4689780     0.1285190
N    -4.4253100     1.2496820     0.1093330
C    -4.9840910    -0.0010050    -0.0000880
C    -4.2544100    -1.1099300    -0.0971060
O    -2.6009580     2.5647160     0.2243840
C    -4.8233490    -2.4946080    -0.2182500
H    -4.9936490     2.0617070     0.1803760
H    -1.3198550     0.4553880     0.0398410
H    -6.0579390    -0.0189790    -0.0016600
H    -4.4932070    -3.1202310     0.6035230
H    -4.4917710    -2.9686160    -1.1353540
H    -5.9073000    -2.4671550    -0.2167380
--
0 1
O    -0.0504540    -1.6178530    -3.0328940
C     1.0293920    -1.0784590    -3.1266030
N     1.0999170     0.3105210    -3.2285410
C     2.2401210     1.0448270    -3.3391500
N     3.4334530     0.5219170    -3.3635050
C     3.4115810    -0.8191700    -3.2674580
C     2.3318990    -1.6598460    -3.1524330
N     2.7451410    -2.9758150    -3.0805400
C     4.0182190    -2.9198150    -3.1499710
N     4.4933620    -1.6323510    -3.2655320
N     2.0870530     2.3690480    -3.4250070
H     0.2128580     0.7884810    -3.2167550
H     4.6831900    -3.7591200    -3.1247610
H     5.4386800    -1.3377250    -3.3349980
H     2.9113910     2.9134700    -3.5059320
H     1.1910290     2.8146150    -3.4104660
units angstrom
}
""")

JSCH_49 = input.process_input("""
molecule dimer {
0 1
O     1.5241600    -0.5494170     3.3837280
C     1.3439910     0.6454500     3.3087260
N     0.0438450     1.1430380     3.2271380
C    -0.3032010     2.4557550     3.1386110
N     0.5626500     3.4294030     3.1191180
C     1.8322280     2.9929620     3.1959900
C     2.2991810     1.7048790     3.2880520
N     3.6791050     1.6903250     3.3455930
C     4.0186060     2.9192810     3.2900230
N     2.9399160     3.7704860     3.1975320
N    -1.6107030     2.7204770     3.0698940
H    -0.6847300     0.4469420     3.2365720
H     5.0225530     3.2920270     3.3102000
H     2.9511880     4.7614640     3.1419340
H    -1.8744930     3.6737330     3.0051240
H    -2.3112160     2.0058100     3.0815320
--
0 1
O    -2.0303850    -1.8889030    -0.1320850
C    -2.7935780    -0.9589290    -0.0670550
N    -2.3254430     0.3299270     0.0230710
C    -3.0734150     1.4709970     0.1028620
N    -4.4253100     1.2514000     0.0875060
C    -4.9840910    -0.0010070    -0.0000700
C    -4.2544100    -1.1114560    -0.0777210
O    -2.6009580     2.5682420     0.1795890
C    -4.8233490    -2.4980370    -0.1746800
H    -4.9936490     2.0645410     0.1443670
H    -1.3198550     0.4560130     0.0318880
H    -6.0579390    -0.0190050    -0.0013290
H    -4.4932070    -3.1092230     0.6578860
H    -4.4917710    -2.9879780    -1.0833720
H    -5.9073000    -2.4705620    -0.1736470
units angstrom
}
""")

JSCH_50 = input.process_input("""
molecule dimer {
0 1
N     1.0423840    -1.6030730     0.1120980
C     1.5033500    -0.3564220     0.0249230
N     0.6443890     0.6605760    -0.0461920
C     1.1104570     1.9058690    -0.1332710
N     2.3570820     2.3057880    -0.1612360
C     3.2027210     1.2804500    -0.0895380
C     2.8736500    -0.0511330     0.0035760
N     4.0080740    -0.8379940     0.0585980
C     4.9747240    -0.0009970     0.0000700
N     4.5630980     1.3036710    -0.0911620
H     0.3540060     2.6689420    -0.1866310
H     6.0151290    -0.2562160     0.0179160
H     5.1340420     2.1141400    -0.1478350
H     0.0593220    -1.7803700     0.1244960
H     1.6856960    -2.3560970     0.1647540
--
0 1
C    -3.3369590    -0.6409430     3.3908960
N    -4.3247580    -1.6157810     3.3763480
C    -4.0409560    -2.9359630     3.2899980
C    -2.7744220    -3.3565600     3.2141470
C    -1.7557370    -2.3300110     3.2303510
N    -2.0543620    -1.0525720     3.3148920
O    -3.6734450     0.5183010     3.4703070
N    -0.4803010    -2.6733740     3.1585030
H    -5.2616330    -1.2870980     3.4334500
H    -4.8798930    -3.6062030     3.2865630
H    -2.5244870    -4.3961080     3.1450650
H    -0.2161270    -3.6266280     3.0937180
H     0.2361240    -1.9652240     3.1691180
units angstrom
}
""")

JSCH_51 = input.process_input("""
molecule dimer {
0 1
N    -1.0423840    -1.6069870     0.0000000
C    -1.5033500    -0.3572920     0.0000000
N    -0.6443890     0.6621890     0.0000000
C    -1.1104570     1.9105230     0.0000000
N    -2.3570820     2.3114180     0.0000000
C    -3.2027210     1.2835770     0.0000000
C    -2.8736500    -0.0512580     0.0000000
N    -4.0080740    -0.8400400     0.0000000
C    -4.9747240    -0.0009990     0.0000000
N    -4.5630980     1.3068540     0.0000000
H    -0.3540060     2.6754590     0.0000000
H    -6.0151290    -0.2568420     0.0000000
H    -5.1340420     2.1193020     0.0000000
H    -0.0593220    -1.7847170     0.0000000
H    -1.6856960    -2.3618500     0.0000000
--
0 1
O     1.5260840    -0.5520650     3.1800000
C     1.3443760     0.6449210     3.1800000
N     0.0425540     1.1448140     3.1800000
C    -0.3063080     2.4600320     3.1800000
N     0.5591430     3.4342300     3.1800000
C     1.8302990     2.9956180     3.1800000
C     2.2991410     1.7049340     3.1800000
N     3.6802460     1.6887540     3.1800000
C     4.0186070     2.9192800     3.1800000
N     2.9380180     3.7730980     3.1800000
N    -1.6152210     2.7266950     3.1800000
H    -0.6858270     0.4484520     3.1800000
H     5.0229680     3.2914560     3.1800000
H     2.9481490     4.7656470     3.1800000
H    -1.8803400     3.6817820     3.1800000
H    -2.3154950     2.0117000     3.1800000
units angstrom
}
""")

JSCH_52 = input.process_input("""
molecule dimer {
0 1
O     2.0303850    -1.8935150     0.0000000
C     2.7935780    -0.9612710     0.0000000
N     2.3254430     0.3307330     0.0000000
C     3.0734150     1.4745890     0.0000000
N     4.4253100     1.2544560     0.0000000
C     4.9840910    -0.0010090     0.0000000
C     4.2544100    -1.1141700     0.0000000
O     2.6009580     2.5745130     0.0000000
C     4.8233490    -2.5041370     0.0000000
H     4.9936490     2.0695820     0.0000000
H     1.3198550     0.4571270     0.0000000
H     6.0579390    -0.0190510     0.0000000
H     4.4932070    -3.0557570     0.8731720
H     4.4917710    -3.0562720    -0.8723020
H     5.9073000    -2.4766570    -0.0008860
--
0 1
C    -3.3390300    -0.6380930     3.1800000
N    -4.3265300    -1.6133420     3.1800000
C    -4.0409560    -2.9359630     3.1800000
C    -2.7728650    -3.3587040     3.1800000
C    -1.7545120    -2.3316960     3.1800000
N    -2.0548730    -1.0518690     3.1800000
O    -3.6771460     0.5233950     3.1800000
N    -0.4776020    -2.6770890     3.1800000
H    -5.2645780    -1.2830450     3.1800000
H    -4.8798220    -3.6063000     3.1800000
H    -2.5215120    -4.4002020     3.1800000
H    -0.2120980    -3.6321740     3.1800000
H     0.2386050    -1.9686390     3.1800000
units angstrom
}
""")

JSCH_53 = input.process_input("""
molecule dimer {
0 1
O     2.0303850    -1.8863100     0.1650310
C     2.7935780    -0.9576130     0.0837800
N     2.3254430     0.3294740    -0.0288250
C     3.0734150     1.4689780    -0.1285190
N     4.4253100     1.2496820    -0.1093330
C     4.9840910    -0.0010050     0.0000880
C     4.2544100    -1.1099300     0.0971060
O     2.6009580     2.5647160    -0.2243840
C     4.8233490    -2.4946080     0.2182500
H     4.9936490     2.0617070    -0.1803760
H     1.3198550     0.4553880    -0.0398410
H     6.0579390    -0.0189790     0.0016600
H     4.4932070    -2.9680270     1.1361760
H     4.4917710    -3.1206680    -0.6026110
H     5.9073000    -2.4673100     0.2149720
--
0 1
C    -1.6419140     2.9739730    -3.0239370
N    -2.8741190     3.6124140    -3.0421140
C    -4.0409900     2.9359160    -3.1500030
C    -4.0487470     1.6026030    -3.2447730
C    -2.7578360     0.9507400    -3.2245270
N    -1.6361750     1.6281560    -3.1188990
O    -0.6443040     3.6509540    -2.9247190
N    -2.6894340    -0.3672360    -3.3142960
H    -2.8516920     4.6040980    -2.9707700
H    -4.9376330     3.5267310    -3.1542940
H    -4.9593830     1.0447610    -3.3310860
H    -3.5136510    -0.9120230    -3.3952410
H    -1.7946790    -0.8299350    -3.3010330
units angstrom
}
""")

JSCH_54 = input.process_input("""
molecule dimer {
0 1
N    -1.0423840    -1.6008720    -0.1400580
C    -1.5033500    -0.3559320    -0.0311400
N    -0.6443890     0.6596690     0.0577140
C    -1.1104570     1.9032530     0.1665130
N    -2.3570820     2.3026230     0.2014530
C    -3.2027210     1.2786930     0.1118710
C    -2.8736500    -0.0510630    -0.0044670
N    -4.0080740    -0.8368430    -0.0732140
C    -4.9747240    -0.0009950    -0.0000870
N    -4.5630980     1.3018810     0.1139000
H    -0.3540060     2.6652780     0.2331820
H    -6.0151290    -0.2558650    -0.0223850
H    -5.1340420     2.1112380     0.1847090
H    -0.0593220    -1.7779260    -0.1555480
H    -1.6856960    -2.3528630    -0.2058490
--
0 1
O    -0.0504540    -1.6178530    -3.0328940
C     1.0293920    -1.0784590    -3.1266030
N     1.0999170     0.3105210    -3.2285410
C     2.2401210     1.0448270    -3.3391500
N     3.4334530     0.5219170    -3.3635050
C     3.4115810    -0.8191700    -3.2674580
C     2.3318990    -1.6598460    -3.1524330
N     2.7451410    -2.9758150    -3.0805400
C     4.0182190    -2.9198150    -3.1499710
N     4.4933620    -1.6323510    -3.2655320
N     2.0870530     2.3690480    -3.4250070
H     0.2128580     0.7884810    -3.2167550
H     4.6831900    -3.7591200    -3.1247610
H     5.4386800    -1.3377250    -3.3349980
H     2.9113910     2.9134700    -3.5059320
H     1.1910290     2.8146150    -3.4104660
units angstrom
}
""")

JSCH_55 = input.process_input("""
molecule dimer {
0 1
O    -1.6803850    -1.8863100    -0.1650310
C    -2.4435780    -0.9576130    -0.0837800
N    -1.9754430     0.3294740     0.0288250
C    -2.7234150     1.4689780     0.1285190
N    -4.0753100     1.2496820     0.1093330
C    -4.6340910    -0.0010050    -0.0000880
C    -3.9044100    -1.1099300    -0.0971060
O    -2.2509580     2.5647160     0.2243840
C    -4.4733490    -2.4946080    -0.2182500
H    -4.6436490     2.0617070     0.1803760
H    -0.9698550     0.4553880     0.0398410
H    -5.7079390    -0.0189790    -0.0016600
H    -4.1432070    -3.1202310     0.6035230
H    -4.1417710    -2.9686160    -1.1353540
H    -5.5573000    -2.4671550    -0.2167380
--
0 1
O     2.4682050    -0.5383510     3.4250310
C     2.5397670     0.6615740     3.3437800
N     1.4045070     1.4276870     3.2311750
C     1.3398450     2.7892110     3.1314810
N     2.5624500     3.4064220     3.1506670
C     3.7496490     2.7230370     3.2600880
C     3.8111350     1.3970020     3.3571060
O     0.3135610     3.3979790     3.0356160
C     5.0853090     0.6111890     3.4782500
H     2.5449500     4.3974240     3.0796240
H     0.5169590     0.9384830     3.2201590
H     4.6289750     3.3396890     3.2616600
H     5.0964880     0.0341320     4.3961760
H     5.1850460    -0.0902010     2.6573890
H     5.9461980     1.2704040     3.4749720
units angstrom
}
""")

JSCH_56 = input.process_input("""
molecule dimer {
0 1
N     1.3923840    -1.6008720     0.1400580
C     1.8533500    -0.3559320     0.0311400
N     0.9943890     0.6596690    -0.0577140
C     1.4604570     1.9032530    -0.1665130
N     2.7070820     2.3026230    -0.2014530
C     3.5527210     1.2786930    -0.1118710
C     3.2236500    -0.0510630     0.0044670
N     4.3580740    -0.8368430     0.0732140
C     5.3247240    -0.0009950     0.0000870
N     4.9130980     1.3018810    -0.1139000
H     0.7040060     2.6652780    -0.2331820
H     6.3651290    -0.2558650     0.0223850
H     5.4840420     2.1112380    -0.1847090
H     0.4093220    -1.7779260     0.1555480
H     2.0356960    -2.3528630     0.2058490
--
0 1
N    -0.1854930    -2.1135550     3.1199420
C    -1.2901800    -1.3773270     3.2288600
N    -1.1922210    -0.0508040     3.3177130
C    -2.3002380     0.6813290     3.4265130
N    -3.5435230     0.2716780     3.4614530
C    -3.6258080    -1.0537530     3.3718710
C    -2.5779730    -1.9361250     3.2555330
N    -3.0338720    -3.2386320     3.1867860
C    -4.3072070    -3.1306000     3.2599130
N    -4.7400060    -1.8346030     3.3739000
H    -2.1361640     1.7424510     3.4931820
H    -4.9991040    -3.9483280     3.2376150
H    -5.6776380    -1.5154120     3.4447090
H     0.7138900    -1.6789650     3.1044520
H    -0.2639350    -3.1000580     3.0541510
units angstrom
}
""")

JSCH_57 = input.process_input("""
molecule dimer {
0 1
N    -1.4867430     1.6920980    -2.3336600
C    -1.5399110     1.6049230    -1.0055780
N    -0.4087210     1.5338080    -0.3037890
C    -0.4671620     1.4467290     1.0245780
N    -1.5291910     1.4187640     1.7901520
C    -2.6502880     1.4904620     1.0763140
C    -2.7488050     1.5835760    -0.2917850
N    -4.0708590     1.6385980    -0.6895780
C    -4.7315520     1.5800700     0.4051650
N    -3.9369080     1.4888380     1.5187780
H     0.4880690     1.3933690     1.5165470
H    -5.7999030     1.5979160     0.4839400
H    -4.2294590     1.4321650     2.4660110
H    -0.6065830     1.7044960    -2.8060620
H    -2.3312660     1.7447540    -2.8510340
--
0 1
N     1.4867430    -1.6920980    -2.3336600
C     1.5399110    -1.6049230    -1.0055780
N     0.4087210    -1.5338080    -0.3037890
C     0.4671620    -1.4467290     1.0245780
N     1.5291910    -1.4187640     1.7901520
C     2.6502880    -1.4904620     1.0763140
C     2.7488050    -1.5835760    -0.2917850
N     4.0708590    -1.6385980    -0.6895780
C     4.7315520    -1.5800700     0.4051650
N     3.9369080    -1.4888380     1.5187780
H    -0.4880690    -1.3933690     1.5165470
H     5.7999030    -1.5979160     0.4839400
H     4.2294590    -1.4321650     2.4660110
H     0.6065830    -1.7044960    -2.8060620
H     2.3312660    -1.7447540    -2.8510340
units angstrom
}
""")

JSCH_58 = input.process_input("""
molecule dimer {
0 1
O     1.3473090     1.4479140    -0.7794320
C     2.3605260     1.5129430    -0.1308130
N     2.3135810     1.6030670     1.2396240
C     3.3775550     1.6828570     2.0937100
N     4.5954240     1.6675010     1.4671030
C     4.7398420     1.5799270     0.1033200
C     3.7027260     1.5022770    -0.7272960
O     3.2672880     1.7595820     3.2832490
C     3.8153430     1.4053200    -2.2218250
H     5.3872210     1.7243610     2.0648200
H     1.3961730     1.6118840     1.6702820
H     5.7555700     1.5786680    -0.2456340
H     3.3501360     0.4957970    -2.5852230
H     3.3109870     2.2369830    -2.7010640
H     4.8546930     1.4081210    -2.5307720
--
0 1
O    -1.3473090    -1.4479140    -0.7794320
C    -2.3605260    -1.5129430    -0.1308130
N    -2.3135810    -1.6030670     1.2396240
C    -3.3775550    -1.6828570     2.0937100
N    -4.5954240    -1.6675010     1.4671030
C    -4.7398420    -1.5799270     0.1033200
C    -3.7027260    -1.5022770    -0.7272960
O    -3.2672880    -1.7595820     3.2832490
C    -3.8153430    -1.4053200    -2.2218250
H    -5.3872210    -1.7243610     2.0648200
H    -1.3961730    -1.6118840     1.6702820
H    -5.7555700    -1.5786680    -0.2456340
H    -3.3109870    -2.2369830    -2.7010640
H    -3.3501360    -0.4957970    -2.5852230
H    -4.8546930    -1.4081210    -2.5307720
units angstrom
}
""")

JSCH_59 = input.process_input("""
molecule dimer {
0 1
N    -1.3923840    -1.6069870     0.0000000
C    -1.8533500    -0.3572920     0.0000000
N    -0.9943890     0.6621890     0.0000000
C    -1.4604570     1.9105230     0.0000000
N    -2.7070820     2.3114180     0.0000000
C    -3.5527210     1.2835770     0.0000000
C    -3.2236500    -0.0512580     0.0000000
N    -4.3580740    -0.8400400     0.0000000
C    -5.3247240    -0.0009990     0.0000000
N    -4.9130980     1.3068540     0.0000000
H    -0.7040060     2.6754590     0.0000000
H    -6.3651290    -0.2568420     0.0000000
H    -5.4840420     2.1193020     0.0000000
H    -0.4093220    -1.7847170     0.0000000
H    -2.0356960    -2.3618500     0.0000000
--
0 1
O     2.4724400    -0.5441800     3.2400000
C     2.5419170     0.6586150     3.2400000
N     1.4037670     1.4287050     3.2400000
C     1.3365470     2.7937510     3.2400000
N     2.5596440     3.4102840     3.2400000
C     3.7496510     2.7230340     3.2400000
C     3.8136270     1.3935720     3.2400000
O     0.3078020     3.4059050     3.2400000
C     5.0909100     0.6034800     3.2400000
H     2.5403210     4.4037960     3.2400000
H     0.5159370     0.9398900     3.2400000
H     4.6290170     3.3396300     3.2400000
H     5.1480540    -0.0368430     4.1131720
H     5.1471940    -0.0381040     2.3676980
H     5.9516930     1.2628420     3.2391140
units angstrom
}
""")

JSCH_60 = input.process_input("""
molecule dimer {
0 1
N    -0.1818990    -2.1185030     3.2400000
C    -1.2893810    -1.3784270     3.2400000
N    -1.1937020    -0.0487650     3.2400000
C    -2.3045120     0.6872100     3.2400000
N    -3.5486930     0.2787930     3.2400000
C    -3.6286790    -1.0498020     3.2400000
C    -2.5778590    -1.9362830     3.2400000
N    -3.0319930    -3.2412190     3.2400000
C    -4.3072050    -3.1306030     3.2400000
N    -4.7429290    -1.8305800     3.2400000
H    -2.1421480     1.7506870     3.2400000
H    -4.9985290    -3.9491190     3.2400000
H    -5.6823780    -1.5088880     3.2400000
H     0.7178820    -1.6844600     3.2400000
H    -0.2586520    -3.1073290     3.2400000
--
0 1
O     1.6803850    -1.8935150     0.0000000
C     2.4435780    -0.9612710     0.0000000
N     1.9754430     0.3307330     0.0000000
C     2.7234150     1.4745890     0.0000000
N     4.0753100     1.2544560     0.0000000
C     4.6340910    -0.0010090     0.0000000
C     3.9044100    -1.1141700     0.0000000
O     2.2509580     2.5745130     0.0000000
C     4.4733490    -2.5041370     0.0000000
H     4.6436490     2.0695820     0.0000000
H     0.9698550     0.4571270     0.0000000
H     5.7079390    -0.0190510     0.0000000
H     4.1432070    -3.0557570     0.8731720
H     4.1417710    -3.0562720    -0.8723020
H     5.5573000    -2.4766570    -0.0008860
units angstrom
}
""")

JSCH_61 = input.process_input("""
molecule dimer {
0 1
C    12.1619966    21.5469940    -0.5249999
N    12.0019966    20.1249944    -0.3349999
C    12.9959964    19.1989946    -0.1290000
N    12.5899965    17.9429950    -0.1260000
C    11.2289969    18.0629949    -0.3469999
C    10.2259971    17.0909952    -0.4599999
N    10.4079971    15.7719956    -0.3739999
N     8.9619975    17.5199951    -0.6819998
C     8.7349976    18.8509947    -0.7899998
N     9.6049973    19.8469944    -0.7019998
C    10.8559970    19.3909946    -0.4999999
H    12.8450824    21.9515608     0.2257099
H    12.5490085    21.7744749    -1.5236356
H    11.1843859    22.0177918    -0.4120399
H    14.0220821    19.5129525     0.0161520
H    11.3436468    15.4109067    -0.2800629
H     9.6382753    15.1406078    -0.5991948
H     7.6909448    19.1156876    -0.9420537
--
0 1
C     3.5239990    12.7489964     2.4389993
N     4.9449986    12.8539964     2.2449994
C     5.8529984    11.8509967     2.0569994
N     7.1019980    12.2539966     2.0409994
C     6.9979980    13.6219962     2.2459994
C     7.9829978    14.6269959     2.3449993
N     9.3019974    14.3749960     2.2649994
N     7.5379979    15.8889955     2.5409993
C     6.2229983    16.1279955     2.6329993
N     5.2169985    15.2499957     2.5399993
C     5.6739984    14.0109961     2.3699993
H     9.6079353    13.4170922     2.2138804
H     9.9620862    15.1183578     2.4869203
H     5.5326604    10.8241690     1.9326585
H     5.9571083    17.1738952     2.7655592
H     3.0968081    13.7487911     2.3499173
H     3.0789261    12.1004316     1.6796125
H     3.2840151    12.3521085     3.4311880
units angstrom
}
""")

JSCH_62 = input.process_input("""
molecule dimer {
0 1
C     3.0629991    16.2869954    -0.5529998
N     4.3679988    15.6949956    -0.7379998
C     5.4889985    16.5069954    -0.6549998
O     5.3979985    17.7169950    -0.4679999
N     6.6749981    15.8589956    -0.7949998
C     6.8699981    14.5069959    -0.9999997
O     8.0199978    14.0679961    -1.0789997
C     5.6559984    13.7139962    -1.1019997
C     5.7709984    12.2569966    -1.4029996
C     4.4739987    14.3319960    -0.9639997
H     7.5313379    16.4637704    -0.7443448
H     6.3741672    11.7424167    -0.6472968
H     4.7881707    11.7797217    -1.4448876
H     6.2751442    12.0930036    -2.3618343
H     3.5293140    13.8026561    -1.0289747
H     2.3790703    15.9479585    -1.3364316
H     2.6423583    16.0249025     0.4245489
H     3.1730521    17.3682771    -0.6086068
--
0 1
C     8.5479976    21.7979939     2.3959993
N     9.1919974    20.5259942     2.6589993
C     8.4229976    19.3799946     2.5429993
O     7.2269980    19.3959946     2.3429993
N     9.0979975    18.2049949     2.7069992
C    10.4579971    18.0869949     2.9379992
O    10.9519969    16.9699952     3.0289992
C    11.2079969    19.3189946     3.0599991
C    12.6759964    19.2659946     3.3619991
C    10.5419970    20.4719943     2.8979992
H     7.4741299    21.6651819     2.5133333
H     8.9049615    22.5495287     3.1049871
H     8.7503455    22.1445498     1.3760436
H    11.0339909    21.4374260     2.9618352
H    13.2133913    18.6878638     2.6029743
H    13.1061963    20.2701373     3.4050200
H    12.8619664    18.7673097     4.3193848
H     8.5371916    17.3217571     2.6353613
units angstrom
}
""")

JSCH_63 = input.process_input("""
molecule dimer {
0 1
C    10.7049970     9.6579973    11.8009967
N    11.0689969    11.0699969    11.9839966
C    10.2199971    12.1419966    11.9589966
N    10.8209970    13.3089963    12.1399966
C    12.1439966    12.9639964    12.2549966
C    13.3189963    13.7529961    12.4509965
O    13.3749963    14.9839958    12.5499965
N    14.4609959    13.0409963    12.5269965
C    14.5119959    11.6719967    12.4369965
N    15.7519956    11.1639969    12.5359965
N    13.4609962    10.8809970    12.2549966
C    12.3209965    11.5909968    12.1779966
H    11.6087247     9.0642815    11.9411017
H    10.3130781     9.4887283    10.7941210
H     9.9552752     9.3611644    12.5389945
H    15.3408647    13.5779012    12.6455145
H     9.1538724    12.0114576    11.8260867
H    15.8197976    10.1594152    12.5501065
H    16.5616854    11.7259467    12.8207994
--
0 1
C    18.8919947     9.6579973     9.7709973
N    18.5279948    11.0699969     9.5879973
C    19.3769946    12.1419966     9.6129973
N    18.7759947    13.3089963     9.4319974
C    17.4529951    12.9639964     9.3169974
C    16.2779954    13.7529961     9.1209974
O    16.2219955    14.9839958     9.0219975
N    15.1359958    13.0409963     9.0449975
C    15.0849958    11.6719967     9.1349974
N    13.8449961    11.1639969     9.0359975
N    16.1359955    10.8809970     9.3169974
C    17.2759952    11.5909968     9.3939974
H    14.2561290    13.5779002     8.9264415
H    13.0353973    11.7259537     8.7509445
H    13.7773141    10.1594092     9.0213535
H    17.9866610     9.0649795     9.6385253
H    19.2909706     9.4904943    10.7753660
H    19.6360815     9.3587324     9.0282525
H    20.4431063    12.0114766     9.7460263
units angstrom
}
""")

JSCH_64 = input.process_input("""
molecule dimer {
0 1
N    10.3469971    14.4959959     8.8169975
C    11.5789968    13.8469961     8.7069976
O    11.6019967    12.6419965     8.4119976
N    12.6939964    14.5549959     8.8809975
C    12.6739964    15.9259955     9.1859974
N    13.8309961    16.5099954     9.3349974
C    11.4219968    16.5639954     9.2669974
C    10.3209971    15.8539956     9.0929975
H     9.3699974    16.4009954     9.1789974
H    11.3019968    17.6379951     9.4699973
H    14.6739959    15.9769955     9.2609974
H    13.8749961    17.4909951     9.5239973
C     9.1059774    13.7460371     8.6280336
H     9.4001314    12.7260934     8.3864956
H     8.5051816    13.7537151     9.5428113
H     8.5206636    14.1698120     7.8064238
--
0 1
N    19.2499946    14.4959959    12.7549964
C    18.0179950    13.8469961    12.8649964
O    17.9949950    12.6419965    13.1599963
N    16.9029953    14.5549959    12.6909964
C    16.9229953    15.9259955    12.3859965
N    15.7659956    16.5099954    12.2369966
C    18.1749949    16.5639954    12.3049966
C    19.2759946    15.8539956    12.4789965
H    20.2269943    16.4009954    12.3929965
H    18.2949949    17.6379951    12.1019966
H    14.9229958    15.9769955    12.3109965
H    15.7219956    17.4909951    12.0479966
C    20.4910143    13.7460371    12.9439604
H    20.1968603    12.7260934    13.1854983
H    21.0918101    13.7537151    12.0291826
H    21.0763281    14.1698120    13.7655701
units angstrom
}
""")

JSCH_65 = input.process_input("""
molecule dimer {
0 1
N    10.9240000    16.7550000     5.5620000
C    11.6470000    17.8510000     5.8140000
N    12.9490000    17.6590000     5.9790000
C    13.0500000    16.2780000     5.7950000
C    14.1950000    15.4230000     5.8560000
N    15.4060000    15.8590000     6.0610000
N    13.9020000    14.1180000     5.6250000
C    12.6770000    13.6430000     5.3990000
N    11.5490000    14.4040000     5.3300000
C    11.8450000    15.6910000     5.5460000
H    11.1804230    18.8265530     5.8822870
H    12.5884030    12.5696370     5.2620740
H    16.1977530    15.2199420     5.9750360
H    15.5570940    16.8510580     6.1500010
C     9.4931860    16.6413650     5.3399050
H     9.0446590    17.6337380     5.4112840
H     9.2947180    16.2234190     4.3499330
H     9.0442270    15.9854440     6.0897950
--
0 1
C    18.8920000     9.6580000     9.7710000
N    18.5280000    11.0700000     9.5880000
C    19.3770000    12.1420000     9.6130000
N    18.7760000    13.3090000     9.4320000
C    17.4530000    12.9640000     9.3170000
C    16.2780000    13.7530000     9.1210000
O    16.2220000    14.9840000     9.0220000
N    15.1360000    13.0410000     9.0450000
C    15.0850000    11.6720000     9.1350000
N    13.8450000    11.1640000     9.0360000
N    16.1360000    10.8810000     9.3170000
C    17.2760000    11.5910000     9.3940000
H    14.2561290    13.5779040     8.9264920
H    13.0354310    11.7259420     8.7508330
H    13.7773690    10.1594100     9.0211800
H    17.9880060     9.0643740     9.6322420
H    19.2851540     9.4890700    10.7774400
H    19.6407390     9.3607520     9.0321720
H    20.4431070    12.0114850     9.7460660
units angstrom
}
""")

JSCH_66 = input.process_input("""
molecule dimer {
0 1
C     9.1690000    13.6920000     8.6010000
N    10.3470000    14.4960000     8.8170000
C    11.5790000    13.8470000     8.7070000
O    11.6020000    12.6420000     8.4120000
N    12.6940000    14.5550000     8.8810000
C    12.6740000    15.9260000     9.1860000
N    13.8310000    16.5100000     9.3350000
C    11.4220000    16.5640000     9.2670000
C    10.3210000    15.8540000     9.0930000
H     9.1403680    12.8642760     9.3131620
H     8.2785600    14.3117950     8.7260530
H     9.1795130    13.2651190     7.5953140
H    11.3501160    17.6252970     9.4808030
H     9.3300790    16.2918180     9.1491660
H    14.7113690    15.9651740     9.2135180
H    13.8876420    17.4962710     9.5342540
--
0 1
N    16.2460000     9.7810000     5.9650000
C    17.5950000    10.0510000     5.9930000
C    18.0920000    11.2690000     5.9020000
C    17.1390000    12.3410000     5.7640000
O    17.4920000    13.5330000     5.6630000
N    15.8280000    12.0550000     5.7130000
C    15.3100000    10.7970000     5.7960000
O    14.1120000    10.5770000     5.7580000
H    18.2280000     9.1744860     6.1031120
C    19.5529600    11.6051630     5.9357380
H    20.1631860    10.7042230     6.0438290
H    19.7760320    12.2828240     6.7658180
H    19.8526100    12.1260780     5.0209680
H    15.1383860    12.8499570     5.6472680
C    15.7717470     8.4029560     6.0779300
H    14.6864640     8.4223240     6.0045990
H    16.1825380     7.7884380     5.2708940
H    16.0652090     7.9755790     7.0417370
units angstrom
}
""")

JSCH_67 = input.process_input("""
molecule dimer {
0 1
H     3.1762460     2.3738070     2.9634160
N     2.3770000     1.8470000     3.2830000
C     1.6370000     2.2160000     4.3790000
H     1.9902970     3.0843050     4.9210710
C     0.5610000     1.4930000     4.7730000
H    -0.0085000     1.7736330     5.6470440
C     0.1830000     0.3990000     3.9430000
N    -0.8510000    -0.3400000     4.2540000
H    -1.1799330    -1.0651510     3.5908230
H    -1.4362750    -0.1022370     5.0377650
N     0.8500000     0.0580000     2.8540000
C     1.9550000     0.7640000     2.4990000
O     2.5580000     0.4150000     1.4830000
--
0 1
H    -3.4958570    -1.4150050    -3.9137580
N    -3.0510000    -1.0010000    -3.1090000
C    -3.5590000    -0.8800000    -1.8360000
H    -4.5790060    -1.1582720    -1.6128580
N    -2.7220000    -0.3740000    -0.9680000
C    -1.5590000    -0.1810000    -1.7250000
C    -0.2720000     0.3480000    -1.4650000
N     0.1070000     0.8840000    -0.3230000
H     1.0433330     1.2579620    -0.3065570
H    -0.5751070     1.2407790     0.3499520
N     0.6670000     0.3750000    -2.4130000
C     0.3480000    -0.0810000    -3.6160000
H     1.1321870    -0.0417550    -4.3673920
N    -0.8160000    -0.5790000    -4.0190000
C    -1.7380000    -0.6050000    -3.0150000
units angstrom
}
""")

JSCH_68 = input.process_input("""
molecule dimer {
0 1
H     0.0112670     4.2441280     0.3057270
N    -0.1600000     4.2010000     1.2990000
C     0.1490000     5.1520000     2.2350000
H     0.8336150     5.9557770     2.0023890
N    -0.3040000     4.9000000     3.4380000
C    -1.1470000     3.7970000     3.2290000
C    -2.0790000     3.1160000     4.0900000
O    -2.3440000     3.3110000     5.2740000
N    -2.7730000     2.0930000     3.4630000
H    -3.4444620     1.6202680     4.0533010
C    -2.5700000     1.7190000     2.1650000
N    -3.2200000     0.6740000     1.7040000
H    -3.7884800     0.1079360     2.3113460
H    -3.0424470     0.3264300     0.7529310
N    -1.7100000     2.3160000     1.3470000
C    -1.0480000     3.3630000     1.9240000
--
0 1
H    -1.2611710    -4.7286740    -2.6257100
N    -1.6090000    -4.2940000    -1.7860000
C    -2.7550000    -4.5990000    -1.0690000
H    -3.5136190    -5.2427470    -1.4922410
N    -2.8650000    -3.9860000     0.0730000
C    -1.6740000    -3.2820000     0.1910000
C    -1.1780000    -2.4570000     1.2560000
O    -1.7150000    -2.1460000     2.3170000
N     0.0980000    -1.9830000     1.0200000
H     0.4562670    -1.3045040     1.7132710
C     0.8280000    -2.2730000    -0.0890000
N     2.0180000    -1.7250000    -0.1770000
H     2.3044660    -0.9690820     0.4476800
H     2.5064670    -1.8555350    -1.0472790
N     0.3920000    -3.0250000    -1.1030000
C    -0.8790000    -3.5010000    -0.9150000
units angstrom
}
""")

JSCH_69 = input.process_input("""
molecule dimer {
0 1
H     4.0780890     0.2050200     6.5267380
N     3.3380000    -0.4520000     6.3380000
C     2.1440000    -0.6140000     7.0100000
H     1.9445960    -0.0744500     7.9251340
N     1.3390000    -1.4880000     6.4770000
C     2.0190000    -1.9110000     5.3320000
C     1.6500000    -2.8430000     4.3020000
O     0.6370000    -3.5330000     4.1980000
N     2.5960000    -2.9520000     3.3010000
H     2.3705000    -3.6388980     2.5623150
C     3.7610000    -2.2490000     3.2730000
N     4.5620000    -2.4690000     2.2580000
H     4.3528370    -3.1696290     1.5459440
H     5.4428290    -1.9835850     2.2550440
N     4.1450000    -1.3880000     4.2160000
C     3.2280000    -1.2560000     5.2240000
--
0 1
H    -1.2611710    -4.7286740    -2.6257100
N    -1.6090000    -4.2940000    -1.7860000
C    -2.7550000    -4.5990000    -1.0690000
H    -3.5136190    -5.2427470    -1.4922410
N    -2.8650000    -3.9860000     0.0730000
C    -1.6740000    -3.2820000     0.1910000
C    -1.1780000    -2.4570000     1.2560000
O    -1.7150000    -2.1460000     2.3170000
N     0.0980000    -1.9830000     1.0200000
H     0.4562670    -1.3045040     1.7132710
C     0.8280000    -2.2730000    -0.0890000
N     2.0180000    -1.7250000    -0.1770000
H     2.3044660    -0.9690820     0.4476800
H     2.5064670    -1.8555350    -1.0472790
N     0.3920000    -3.0250000    -1.1030000
C    -0.8790000    -3.5010000    -0.9150000
units angstrom
}
""")

JSCH_70 = input.process_input("""
molecule dimer {
0 1
H     3.1762460     2.3738070     2.9634160
N     2.3770000     1.8470000     3.2830000
C     1.6370000     2.2160000     4.3790000
H     1.9902970     3.0843050     4.9210710
C     0.5610000     1.4930000     4.7730000
H    -0.0085000     1.7736330     5.6470440
C     0.1830000     0.3990000     3.9430000
N    -0.8510000    -0.3400000     4.2540000
H    -1.1799330    -1.0651510     3.5908230
H    -1.4362750    -0.1022370     5.0377650
N     0.8500000     0.0580000     2.8540000
C     1.9550000     0.7640000     2.4990000
O     2.5580000     0.4150000     1.4830000
--
0 1
H     3.2823840    -6.1134940    -1.3105350
N     2.5530000    -6.0070000    -0.6210000
C     1.3990000    -6.7620000    -0.6490000
H     1.3017290    -7.4646550    -1.4662410
C     0.4550000    -6.5890000     0.3070000
H    -0.4593850    -7.1648600     0.2947650
C     0.7210000    -5.6290000     1.3280000
N    -0.1590000    -5.3940000     2.2700000
H    -1.0266130    -5.9017830     2.3125200
H     0.0709100    -4.7127400     3.0149280
N     1.8460000    -4.9310000     1.3860000
C     2.7800000    -5.0940000     0.4140000
O     3.8210000    -4.4400000     0.4780000
units angstrom
}
""")

JSCH_71 = input.process_input("""
molecule dimer {
0 1
O    -1.2390176    -2.5490521     0.6548924
C    -1.0284571    -1.3714583     0.9008651
N    -0.0318511    -0.9949528     1.8248233
C     0.3841646     0.2706806     2.1182164
N    -0.1910285     1.3513281     1.6527710
C    -1.2092305     1.0513624     0.8089237
C    -1.6565083    -0.1915101     0.3706051
N    -2.6541580    -0.0639048    -0.5661534
C    -2.8177333     1.2431899    -0.6803818
N    -1.9753657     1.9574414     0.1290579
N     1.4525454     0.3558875     2.9872621
H     0.4866119    -1.7695272     2.2174674
H    -3.5338415     1.7253425    -1.3240899
H    -1.9138820     2.9580997     0.2181746
H     1.7298659     1.3225221     3.0797421
H     2.2376547    -0.1901480     2.6476325
--
0 1
C     2.2123373    -0.0590839    -0.4645529
N     2.1205577     1.1822577    -1.1169007
C     1.2003987     1.4553092    -2.0711004
C     0.3300220     0.4917324    -2.4615962
C     0.4626198    -0.7818118    -1.8195186
N     1.3658705    -1.0412675    -0.8919664
O     3.0203933    -0.1851683     0.4516286
N    -0.3645719    -1.7922584    -2.1870353
H     2.7522574     1.8928231    -0.7832658
H     1.2077958     2.4531651    -2.4849122
H    -0.4090987     0.6756472    -3.2236499
H    -1.2619684    -1.5144746    -2.5505972
H    -0.4171470    -2.5417743    -1.5096178
units angstrom
}
""")

JSCH_72 = input.process_input("""
molecule dimer {
0 1
O    -1.6144948    -2.7570519    -0.2060980
C    -1.1842160    -1.7852952     0.3968248
N    -0.0652525    -1.9057030     1.2481671
C     0.5967440    -0.9002665     1.8887926
N     0.1973220     0.3471640     1.8866901
C    -0.9263931     0.5128083     1.1439041
C    -1.6366545    -0.4213094     0.3921297
N    -2.6684418     0.1720265    -0.2870095
C    -2.5827038     1.4477282     0.0619406
N    -1.5537529     1.7066351     0.9265244
N     1.7211367    -1.2836170     2.5904485
H     0.3310406    -2.8355291     1.2850314
H    -3.2456406     2.2307336    -0.2693001
C    -1.1552104     2.9806914     1.4836593
H     2.1791458    -0.4655292     2.9651066
H     2.3669253    -1.7840635     1.9899163
H    -0.1579041     3.2445863     1.1378850
H    -1.8652414     3.7371235     1.1620610
H    -1.1470859     2.9258678     2.5690805
--
0 1
C     1.9196368    -0.2893692    -0.7963336
N     1.8412129     1.1192423    -0.8402021
C     0.8544988     1.7477483    -1.5184980
C    -0.1003666     1.0437295    -2.1815479
C    -0.0071603    -0.3764662    -2.1036673
N     0.9713285    -0.9972682    -1.4703627
O     2.8208711    -0.8014248    -0.1303355
N    -0.9502154    -1.1428384    -2.7108142
C     2.8076706     1.8576425    -0.0455040
H     0.8766575     2.8292157    -1.5014613
H    -0.8834995     1.5417865    -2.7295416
H    -1.8580411    -0.7155517    -2.8051204
H    -0.9899772    -2.0895195    -2.3567972
H     3.8087243     1.5028551    -0.2699237
H     2.6081184     1.6978177     1.0127736
H     2.7199672     2.9140199    -0.2854865
units angstrom
}
""")

JSCH_73 = input.process_input("""
molecule dimer {
0 1
N     0.2793014     2.4068393    -0.6057517
C    -1.0848570     2.4457461    -0.5511608
H    -1.6594403     3.0230294    -1.2560905
N    -1.5977117     1.7179877     0.4287543
C    -0.4897255     1.1714358     1.0301910
C    -0.3461366     0.2914710     2.1172343
N    -1.4187090    -0.1677767     2.8101441
H    -1.2388750    -0.9594802     3.4047578
H    -2.2918734    -0.1788223     2.3073619
N     0.8857630    -0.0700763     2.4919494
C     1.9352348     0.4072878     1.7968022
H     2.9060330     0.0788414     2.1458181
N     1.9409775     1.2242019     0.7402202
C     0.6952186     1.5779858     0.4063984
H     0.8610073     2.8298045    -1.3104502
--
0 1
N     1.2754606    -0.6478993    -1.9779104
C     1.4130533    -1.5536850    -0.9550667
H     2.4258769    -1.8670780    -0.7468778
C     0.3575976    -2.0239499    -0.2530575
C     0.4821292    -3.0179494     0.8521221
H     0.1757705    -2.5756065     1.7986281
H    -0.1601691    -3.8770412     0.6639498
H     1.5112443    -3.3572767     0.9513659
C    -0.9684711    -1.5298112    -0.5939792
O    -2.0029280    -1.8396957    -0.0199453
N    -0.9956916    -0.6383870    -1.6720420
H    -1.9014057    -0.2501720    -1.8985760
C     0.0684702    -0.1191762    -2.3763759
O    -0.0397875     0.7227006    -3.2531083
H     2.0853289    -0.2760176    -2.4454577
units angstrom
}
""")

JSCH_74 = input.process_input("""
molecule dimer {
0 1
N    -0.3455004     1.7703632     1.4950792
C    -1.6474050     1.3634505     1.5386766
H    -2.4523693     2.0803127     1.5703490
N    -1.8053639     0.0450392     1.5375118
C    -0.5193842    -0.4240596     1.4834056
C     0.0152186    -1.7249725     1.4821754
N    -0.7782381    -2.8218524     1.5417158
H    -0.3281681    -3.6995564     1.3432557
H    -1.7192874    -2.7111068     1.1983318
N     1.3452903    -1.8718583     1.4651757
C     2.1159101    -0.7701212     1.4213994
H     3.1830548    -0.9527061     1.4028830
N     1.7419114     0.5131994     1.4043323
C     0.4096081     0.6245403     1.4501833
C     0.1512980     3.1326941     1.4689984
H    -0.0424219     3.5749692     0.4946347
H    -0.3347704     3.7141185     2.2479916
H     1.2201020     3.0964449     1.6609900
--
0 1
N     0.8076098     1.0547322    -1.6591556
C     1.2548662    -0.2426109    -1.7103022
H     2.3275169    -0.3452707    -1.8079765
C     0.4450062    -1.3265062    -1.6516166
C     0.9521849    -2.7269269    -1.7314581
H     0.7400336    -3.2617591    -0.8079097
H     0.4633129    -3.2624261    -2.5442266
H     2.0282910    -2.7371647    -1.8922437
C    -0.9813923    -1.1031917    -1.5070970
O    -1.8286775    -1.9792278    -1.3834794
N    -1.3482304     0.2425066    -1.5277277
H    -2.3301840     0.4328271    -1.3774912
C    -0.5338001     1.3531361    -1.5364698
O    -0.9719491     2.4936714    -1.4469221
C     1.7769330     2.1404377    -1.6201082
H     2.2553696     2.1734198    -0.6418794
H     2.5269559     1.9765736    -2.3901046
H     1.2518503     3.0690500    -1.8139519
units angstrom
}
""")

JSCH_75 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     0.9181960    -0.9215090     3.4000000
N    -0.3693690    -1.5141310     3.4000000
C    -1.5252510    -0.8082010     3.4000000
C    -1.4858600     0.5568310     3.4000000
C    -0.1723650     1.1455600     3.4000000
N     0.9526960     0.4540270     3.4000000
O     1.9020460    -1.6508420     3.4000000
N    -0.0596430     2.5018840     3.4000000
H    -0.3693760    -2.5309310     3.4000000
H    -2.4533460    -1.3813360     3.4000000
H    -2.3977890     1.1506030     3.4000000
H    -0.8684590     3.1052900     3.4000000
H     0.8696630     2.9040660     3.4000000
units angstrom
}
""")

JSCH_76 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     1.2571480     0.3344260     3.3000000
N     1.1265920    -1.0769480     3.3000000
C    -0.0627020    -1.7250070     3.3000000
C    -1.2251600    -1.0083780     3.3000000
C    -1.0782670     0.4235070     3.3000000
N     0.0831490     1.0520730     3.3000000
O     2.3806940     0.8218000     3.3000000
N    -2.1965170     1.1992900     3.3000000
H     2.0071630    -1.5853550     3.3000000
H    -0.0304010    -2.8153280     3.3000000
H    -2.1953460    -1.5012440     3.3000000
H    -3.1234900     0.8005370     3.3000000
H    -2.0801630     2.2051830     3.3000000
units angstrom
}
""")

JSCH_77 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     0.3389520     1.2559350     3.3000000
N     1.4959600     0.4371830     3.3000000
C     1.4625490    -0.9168050     3.3000000
C     0.2607010    -1.5652080     3.3000000
C    -0.9059010    -0.7220530     3.3000000
N    -0.8695470     0.5980460     3.3000000
O     0.4786470     2.4726420     3.3000000
N    -2.1368730    -1.3025950     3.3000000
H     2.3765390     0.9455770     3.3000000
H     2.4229450    -1.4339920     3.3000000
H     0.2024430    -2.6518470     3.3000000
H    -2.2550290    -2.3047530     3.3000000
H    -2.9498260    -0.6988840     3.3000000
units angstrom
}
""")

JSCH_78 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -0.9181960     0.9215090     3.3000000
N     0.3693690     1.5141310     3.3000000
C     1.5252510     0.8082010     3.3000000
C     1.4858600    -0.5568310     3.3000000
C     0.1723650    -1.1455600     3.3000000
N    -0.9526960    -0.4540270     3.3000000
O    -1.9020460     1.6508420     3.3000000
N     0.0596430    -2.5018840     3.3000000
H     0.3693760     2.5309310     3.3000000
H     2.4533460     1.3813360     3.3000000
H     2.3977890    -1.1506030     3.3000000
H     0.8684590    -3.1052900     3.3000000
H    -0.8696630    -2.9040660     3.3000000
units angstrom
}
""")

JSCH_79 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     1.9181950    -0.9215090     3.3000000
N     0.6306310    -1.5141310     3.3000000
C    -0.5252510    -0.8082010     3.3000000
C    -0.4858600     0.5568310     3.3000000
C     0.8276350     1.1455600     3.3000000
N     1.9526960     0.4540270     3.3000000
O     2.9020460    -1.6508420     3.3000000
N     0.9403570     2.5018840     3.3000000
H     0.6306240    -2.5309310     3.3000000
H    -1.4533460    -1.3813360     3.3000000
H    -1.3977890     1.1506030     3.3000000
H     0.1315410     3.1052900     3.3000000
H     1.8696630     2.9040660     3.3000000
units angstrom
}
""")

JSCH_80 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     0.9181960     0.0784910     3.3000000
N    -0.3693690    -0.5141310     3.3000000
C    -1.5252510     0.1917990     3.3000000
C    -1.4858600     1.5568310     3.3000000
C    -0.1723650     2.1455600     3.3000000
N     0.9526960     1.4540270     3.3000000
O     1.9020460    -0.6508420     3.3000000
N    -0.0596430     3.5018840     3.3000000
H    -0.3693760    -1.5309310     3.3000000
H    -2.4533460    -0.3813360     3.3000000
H    -2.3977890     2.1506030     3.3000000
H    -0.8684590     4.1052890     3.3000000
H     0.8696630     3.9040660     3.3000000
units angstrom
}
""")

JSCH_81 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     2.9181950    -2.9215090     3.3000000
N     1.6306310    -3.5141310     3.3000000
C     0.4747490    -2.8082010     3.3000000
C     0.5141400    -1.4431690     3.3000000
C     1.8276350    -0.8544400     3.3000000
N     2.9526960    -1.5459730     3.3000000
O     3.9020460    -3.6508420     3.3000000
N     1.9403570     0.5018840     3.3000000
H     1.6306240    -4.5309310     3.3000000
H    -0.4533460    -3.3813360     3.3000000
H    -0.3977890    -0.8493970     3.3000000
H     1.1315410     1.1052900     3.3000000
H     2.8696630     0.9040660     3.3000000
units angstrom
}
""")

JSCH_82 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     0.0818040     0.9215090     3.3000000
N     1.3693690     1.5141310     3.3000000
C     2.5252510     0.8082010     3.3000000
C     2.4858600    -0.5568310     3.3000000
C     1.1723650    -1.1455600     3.3000000
N     0.0473040    -0.4540270     3.3000000
O    -0.9020460     1.6508420     3.3000000
N     1.0596430    -2.5018840     3.3000000
H     1.3693760     2.5309310     3.3000000
H     3.4533460     1.3813360     3.3000000
H     3.3977890    -1.1506030     3.3000000
H     1.8684590    -3.1052900     3.3000000
H     0.1303370    -2.9040660     3.3000000
units angstrom
}
""")

JSCH_83 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -1.9181960     0.9215090     3.3000000
N    -0.6306310     1.5141310     3.3000000
C     0.5252510     0.8082010     3.3000000
C     0.4858600    -0.5568310     3.3000000
C    -0.8276350    -1.1455600     3.3000000
N    -1.9526960    -0.4540270     3.3000000
O    -2.9020460     1.6508420     3.3000000
N    -0.9403570    -2.5018840     3.3000000
H    -0.6306240     2.5309320     3.3000000
H     1.4533460     1.3813360     3.3000000
H     1.3977890    -1.1506030     3.3000000
H    -0.1315410    -3.1052900     3.3000000
H    -1.8696630    -2.9040660     3.3000000
units angstrom
}
""")

JSCH_84 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -0.9181960     1.9215090     3.3000000
N     0.3693690     2.5141310     3.3000000
C     1.5252510     1.8082010     3.3000000
C     1.4858600     0.4431690     3.3000000
C     0.1723650    -0.1455600     3.3000000
N    -0.9526960     0.5459730     3.3000000
O    -1.9020460     2.6508420     3.3000000
N     0.0596430    -1.5018840     3.3000000
H     0.3693760     3.5309310     3.3000000
H     2.4533460     2.3813360     3.3000000
H     2.3977890    -0.1506030     3.3000000
H     0.8684590    -2.1052900     3.3000000
H    -0.8696630    -1.9040660     3.3000000
units angstrom
}
""")

JSCH_85 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -0.9181960    -0.0784910     3.3000000
N     0.3693690     0.5141310     3.3000000
C     1.5252510    -0.1917990     3.3000000
C     1.4858600    -1.5568310     3.3000000
C     0.1723650    -2.1455600     3.3000000
N    -0.9526960    -1.4540270     3.3000000
O    -1.9020460     0.6508420     3.3000000
N     0.0596430    -3.5018840     3.3000000
H     0.3693760     1.5309310     3.3000000
H     2.4533460     0.3813360     3.3000000
H     2.3977890    -2.1506030     3.3000000
H     0.8684590    -4.1052900     3.3000000
H    -0.8696630    -3.9040660     3.3000000
units angstrom
}
""")

JSCH_86 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C     1.0818040    -1.0784910     3.3000000
N     2.3693690    -0.4858690     3.3000000
C     3.5252510    -1.1917990     3.3000000
C     3.4858600    -2.5568310     3.3000000
C     2.1723650    -3.1455600     3.3000000
N     1.0473050    -2.4540270     3.3000000
O     0.0979540    -0.3491580     3.3000000
N     2.0596430    -4.5018840     3.3000000
H     2.3693760     0.5309310     3.3000000
H     4.4533460    -0.6186640     3.3000000
H     4.3977890    -3.1506030     3.3000000
H     2.8684590    -5.1052900     3.3000000
H     1.1303370    -4.9040660     3.3000000
units angstrom
}
""")

JSCH_87 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -0.9181960     2.9215090     3.3000000
N     0.3693690     3.5141310     3.3000000
C     1.5252510     2.8082010     3.3000000
C     1.4858600     1.4431690     3.3000000
C     0.1723660     0.8544400     3.3000000
N    -0.9526960     1.5459730     3.3000000
O    -1.9020460     3.6508420     3.3000000
N     0.0596430    -0.5018840     3.3000000
H     0.3693760     4.5309310     3.3000000
H     2.4533460     3.3813360     3.3000000
H     2.3977890     0.8493970     3.3000000
H     0.8684590    -1.1052900     3.3000000
H    -0.8696630    -0.9040660     3.3000000
units angstrom
}
""")

JSCH_88 = input.process_input("""
molecule dimer {
0 1
C     0.9181960    -0.9215090     0.0000000
N    -0.3693690    -1.5141310     0.0000000
C    -1.5252510    -0.8082010     0.0000000
C    -1.4858600     0.5568310     0.0000000
C    -0.1723650     1.1455600     0.0000000
N     0.9526960     0.4540270     0.0000000
O     1.9020460    -1.6508420     0.0000000
N    -0.0596430     2.5018840     0.0000000
H    -0.3693760    -2.5309310     0.0000000
H    -2.4533460    -1.3813360     0.0000000
H    -2.3977890     1.1506030     0.0000000
H    -0.8684590     3.1052900     0.0000000
H     0.8696630     2.9040660     0.0000000
--
0 1
C    -0.9181960    -0.9215090     3.3000000
N     0.3693690    -1.5141310     3.3000000
C     1.5252510    -0.8082010     3.3000000
C     1.4858600     0.5568310     3.3000000
C     0.1723650     1.1455600     3.3000000
N    -0.9526960     0.4540270     3.3000000
O    -1.9020460    -1.6508420     3.3000000
N     0.0596430     2.5018840     3.3000000
H     0.3693760    -2.5309310     3.3000000
H     2.4533460    -1.3813360     3.3000000
H     2.3977890     1.1506030     3.3000000
H     0.8684590     3.1052900     3.3000000
H    -0.8696630     2.9040660     3.3000000
units angstrom
}
""")

JSCH_89 = input.process_input("""
molecule dimer {
0 1
N    -1.9000000    -0.3579200     0.0000000
C    -1.9000000     0.9808800     0.0000000
N    -0.8497640     1.8329300     0.0000000
C     0.3886960     1.3139590     0.0000000
C     0.5403660    -0.0879600     0.0000000
C    -0.6427490    -0.8323270     0.0000000
N    -0.2184070    -2.1434690     0.0000000
C     1.1532860    -2.1196680     0.0000000
N     1.6612330    -0.8947060     0.0000000
N     1.4545080     2.1469300     0.0000000
H    -2.8785740     1.4564190     0.0000000
H     1.7387770    -3.0306410     0.0000000
H    -0.8201870    -2.9587230     0.0000000
H     1.2990760     3.1441900     0.0000000
H     2.3918120     1.7733760     0.0000000
--
0 1
N     2.0400330     1.0244080     3.3000000
C     3.1994670     0.3550080     3.3000000
N     3.4122460    -0.9805480     3.3000000
C     2.3435740    -1.7936000     3.3000000
C     1.0536410    -1.2239910     3.3000000
C     1.0005580     0.1728010     3.3000000
N    -0.3470950     0.4608810     3.3000000
C    -1.0123290    -0.7389400     3.3000000
N    -0.2054550    -1.7913170     3.3000000
N     2.5320410    -3.1331060     3.3000000
H     4.1005830     0.9647080     3.3000000
H    -2.0940000    -0.7905040     3.3000000
H    -0.7522350     1.3896650     3.3000000
H     3.4734100    -3.4971280     3.3000000
H     1.7398820    -3.7580580     3.3000000
units angstrom
}
""")

JSCH_90 = input.process_input("""
molecule dimer {
0 1
O     0.2392880    -2.6920590     0.0000000
C     0.2392880    -1.4664590     0.0000000
N     1.4831650    -0.7585720     0.0000000
C     1.6585390     0.6049970     0.0000000
N     0.6694090     1.4698410     0.0000000
C    -0.5424070     0.8439610     0.0000000
C    -0.8433090    -0.5171760     0.0000000
N    -2.2044110    -0.7367510     0.0000000
C    -2.7203200     0.4816210     0.0000000
N    -1.7621040     1.4687240     0.0000000
N     2.9429720     1.0619610     0.0000000
H     2.2894380    -1.3782560     0.0000000
H    -3.7780480     0.7114080     0.0000000
H    -1.9055140     2.4718250     0.0000000
H     3.0816010     2.0608880     0.0000000
H     3.7442410     0.4525240     0.0000000
--
0 1
O    -3.4927120     1.0318190     3.3000000
C    -2.2857320     0.8189950     3.3000000
N    -1.8045960    -0.5289080     3.3000000
C    -0.4921960    -0.9383990     3.3000000
N     0.5312690    -0.1144750     3.3000000
C     0.1253280     1.1876140     3.3000000
C    -1.1628790     1.7203040     3.3000000
N    -1.1427660     3.0988560     3.3000000
C     0.1466840     3.3953590     3.3000000
N     0.9523970     2.2802920     3.3000000
N    -0.2652140    -2.2826700     3.3000000
H    -2.5548740    -1.2153250     3.3000000
H     0.5566520     4.3971160     3.3000000
H     1.9651610     2.2473370     3.3000000
H     0.6944640    -2.5926540     3.3000000
H    -1.0045320    -2.9659380     3.3000000
units angstrom
}
""")

JSCH_91 = input.process_input("""
molecule dimer {
0 1
N    -1.9000000    -0.3579200     0.0000000
C    -1.9000000     0.9808800     0.0000000
N    -0.8497640     1.8329300     0.0000000
C     0.3886960     1.3139590     0.0000000
C     0.5403660    -0.0879600     0.0000000
C    -0.6427490    -0.8323270     0.0000000
N    -0.2184070    -2.1434690     0.0000000
C     1.1532860    -2.1196680     0.0000000
N     1.6612330    -0.8947060     0.0000000
N     1.4545080     2.1469300     0.0000000
H    -2.8785740     1.4564190     0.0000000
H     1.7387770    -3.0306410     0.0000000
H    -0.8201870    -2.9587230     0.0000000
H     1.2990760     3.1441900     0.0000000
H     2.3918120     1.7733760     0.0000000
--
0 1
C     0.2481770    -2.1847000     3.3000000
N     0.9375100    -0.9462160     3.3000000
C     0.3220220     0.2602560     3.3000000
C    -1.0420250     0.3253500     3.3000000
C    -1.7294600    -0.9392870     3.3000000
N    -1.1259700    -2.1139300     3.3000000
O     0.9001510    -3.2214350     3.3000000
N    -3.0904320    -0.9479780     3.3000000
H     1.9513330    -1.0239520     3.3000000
H     0.9644400     1.1418140     3.3000000
H    -1.5643350     1.2800070     3.3000000
H    -3.6302290    -0.0953930     3.3000000
H    -3.5624890    -1.8438140     3.3000000
units angstrom
}
""")

JSCH_92 = input.process_input("""
molecule dimer {
0 1
O     0.2392880    -2.6920580     0.0000000
C     0.2392880    -1.4664580     0.0000000
N     1.4831650    -0.7585710     0.0000000
C     1.6585390     0.6049980     0.0000000
N     0.6694090     1.4698420     0.0000000
C    -0.5424070     0.8439620     0.0000000
C    -0.8433090    -0.5171750     0.0000000
N    -2.2044110    -0.7367500     0.0000000
C    -2.7203200     0.4816220     0.0000000
N    -1.7621040     1.4687250     0.0000000
N     2.9429720     1.0619620     0.0000000
H     2.2894380    -1.3782550     0.0000000
H    -3.7780480     0.7114090     0.0000000
H    -1.9055140     2.4718250     0.0000000
H     3.0816010     2.0608890     0.0000000
H     3.7442410     0.4525250     0.0000000
--
0 1
N    -2.9334410     0.9976940    -3.3000000
C    -2.6871670     2.3136480    -3.3000000
N    -1.4981170     2.9579650    -3.3000000
C    -0.3762570     2.2200360    -3.3000000
C    -0.4850590     0.8141390    -3.3000000
C    -1.7849120     0.3001090    -3.3000000
N    -1.6089980    -1.0667170    -3.3000000
C    -0.2563330    -1.2956460    -3.3000000
N     0.4682790    -0.1850240    -3.3000000
N     0.8245940     2.8427340    -3.3000000
H    -3.5615660     2.9610820    -3.3000000
H     0.1515920    -2.2987750    -3.3000000
H    -2.3504750    -1.7573600    -3.3000000
H     0.8552610     3.8515680    -3.3000000
H     1.6771880     2.3031370    -3.3000000
units angstrom
}
""")

JSCH_93 = input.process_input("""
molecule dimer {
0 1
C    -1.2210000    -0.4488000     0.0000000
N    -1.2210000     0.9686000     0.0000000
C    -0.0964540     1.7234490     0.0000000
C     1.1270700     1.1169400     0.0000000
C     1.1126920    -0.3223880     0.0000000
N     0.0141100    -1.0552590     0.0000000
O    -2.2948780    -1.0375910     0.0000000
N     2.2976450    -0.9918710     0.0000000
H    -2.1446570     1.3937360     0.0000000
H    -0.2290470     2.8061600     0.0000000
H     2.0477340     1.6970750     0.0000000
H     3.1839480    -0.5094300     0.0000000
H     2.2744390    -2.0042050     0.0000000
--
0 1
C     0.8210000     1.4488000     3.3000000
N     0.8210000     0.0314000     3.3000000
C    -0.3035460    -0.7234490     3.3000000
C    -1.5270700    -0.1169400     3.3000000
C    -1.5126920     1.3223880     3.3000000
N    -0.4141100     2.0552590     3.3000000
O     1.8948780     2.0375920     3.3000000
N    -2.6976450     1.9918700     3.3000000
H     1.7446570    -0.3937350     3.3000000
H    -0.1709520    -1.8061600     3.3000000
H    -2.4477340    -0.6970760     3.3000000
H    -3.5839490     1.5094290     3.3000000
H    -2.6744400     3.0042040     3.3000000
units angstrom
}
""")

JSCH_94 = input.process_input("""
molecule dimer {
0 1
N    -1.9000000    -0.3579200     0.0000000
C    -1.9000000     0.9808800     0.0000000
N    -0.8497640     1.8329300     0.0000000
C     0.3886960     1.3139590     0.0000000
C     0.5403660    -0.0879600     0.0000000
C    -0.6427490    -0.8323270     0.0000000
N    -0.2184070    -2.1434690     0.0000000
C     1.1532860    -2.1196680     0.0000000
N     1.6612330    -0.8947060     0.0000000
N     1.4545080     2.1469300     0.0000000
H    -2.8785740     1.4564190     0.0000000
H     1.7387770    -3.0306410     0.0000000
H    -0.8201870    -2.9587230     0.0000000
H     1.2990760     3.1441900     0.0000000
H     2.3918120     1.7733760     0.0000000
--
0 1
O     0.8290540    -2.7349420     3.3000000
C     0.6180080    -1.5257000     3.3000000
N     1.6882380    -0.6117020     3.3000000
C     1.6473780     0.7731900     3.3000000
N     0.3483000     1.2653660     3.3000000
C    -0.7746550     0.4682380     3.3000000
C    -0.6885800    -0.8807150     3.3000000
O     2.6370680     1.4932450     3.3000000
H    -1.5708660    -1.5073280     3.3000000
H     2.6203010    -1.0186580     3.3000000
H     0.2745460     2.2757330     3.3000000
H    -1.7216600     0.9981100     3.3000000
units angstrom
}
""")

JSCH_95 = input.process_input("""
molecule dimer {
0 1
O     0.2392880    -2.6920580     0.0000000
C     0.2392880    -1.4664580     0.0000000
N     1.4831650    -0.7585710     0.0000000
C     1.6585390     0.6049980     0.0000000
N     0.6694090     1.4698420     0.0000000
C    -0.5424070     0.8439620     0.0000000
C    -0.8433090    -0.5171750     0.0000000
N    -2.2044110    -0.7367500     0.0000000
C    -2.7203200     0.4816220     0.0000000
N    -1.7621040     1.4687250     0.0000000
N     2.9429720     1.0619620     0.0000000
H     2.2894380    -1.3782550     0.0000000
H    -3.7780480     0.7114090     0.0000000
H    -1.9055140     2.4718250     0.0000000
H     3.0816010     2.0608890     0.0000000
H     3.7442410     0.4525250     0.0000000
--
0 1
C     2.5140740     0.5212730    -3.3000000
N     1.2263530     1.1135940    -3.3000000
C     0.0706410     0.4073870    -3.3000000
C     0.1103570    -0.9576360    -3.3000000
C     1.4239940    -1.5460500    -3.3000000
N     2.5488900    -0.8542470    -3.3000000
O     3.4977500     1.2508430    -3.3000000
N     1.5370420    -2.9023480    -3.3000000
H     1.2261020     2.1303940    -3.3000000
H    -0.8575930     0.9802990    -3.3000000
H    -0.8014280    -1.5516270    -3.3000000
H     0.7283690    -3.5059470    -3.3000000
H     2.4664440    -3.3043070    -3.3000000
units angstrom
}
""")

JSCH_96 = input.process_input("""
molecule dimer {
0 1
C    -1.2210000    -0.4488000     0.0000000
N    -1.2210000     0.9686000     0.0000000
C    -0.0964540     1.7234490     0.0000000
C     1.1270700     1.1169400     0.0000000
C     1.1126920    -0.3223880     0.0000000
N     0.0141100    -1.0552590     0.0000000
O    -2.2948780    -1.0375910     0.0000000
N     2.2976450    -0.9918710     0.0000000
H    -2.1446570     1.3937360     0.0000000
H    -0.2290470     2.8061600     0.0000000
H     2.0477340     1.6970750     0.0000000
H     3.1839480    -0.5094300     0.0000000
H     2.2744390    -2.0042050     0.0000000
--
0 1
O    -0.2290330     2.1280320     3.3000000
C    -0.3146620     0.9035020     3.3000000
N     0.8438480     0.1043590     3.3000000
C     0.9455600    -1.2773930     3.3000000
N    -0.2960490    -1.9004980     3.3000000
C    -1.4949900    -1.2230250     3.3000000
C    -1.5480330     0.1276340     3.3000000
O     2.0040190    -1.8919000     3.3000000
H    -2.4900550     0.6602360     3.3000000
H     1.7291390     0.6049670     3.3000000
H    -0.2655530    -2.9130890     3.3000000
H    -2.3825080    -1.8474300     3.3000000
units angstrom
}
""")

JSCH_97 = input.process_input("""
molecule dimer {
0 1
O    -0.4072070    -2.5021900     0.0000000
C    -0.4072070    -1.2746690     0.0000000
N     0.8042290    -0.5582850     0.0000000
C     1.0020800     0.8130100     0.0000000
N    -0.1930340     1.5212070     0.0000000
C    -1.4363170     0.9290170     0.0000000
C    -1.5834480    -0.4146490     0.0000000
O     2.1008310     1.3521860     0.0000000
H    -2.5603280    -0.8802410     0.0000000
H     1.6524450    -1.1194300     0.0000000
H    -0.0919780     2.5292090     0.0000000
H    -2.2781210     1.6138160     0.0000000
--
0 1
O    -1.1927920     2.1021900     3.3000000
C    -1.1927930     0.8746690     3.3000000
N    -2.4042290     0.1582850     3.3000000
C    -2.6020800    -1.2130100     3.3000000
N    -1.4069670    -1.9212070     3.3000000
C    -0.1636830    -1.3290170     3.3000000
C    -0.0165520     0.0146480     3.3000000
O    -3.7008310    -1.7521850     3.3000000
H     0.9603280     0.4802400     3.3000000
H    -3.2524450     0.7194310     3.3000000
H    -1.5080230    -2.9292090     3.3000000
H     0.6781200    -2.0138170     3.3000000
units angstrom
}
""")

JSCH_98 = input.process_input("""
molecule dimer {
0 1
O     0.2392880    -2.6920580     0.0000000
C     0.2392880    -1.4664580     0.0000000
N     1.4831650    -0.7585710     0.0000000
C     1.6585390     0.6049980     0.0000000
N     0.6694090     1.4698420     0.0000000
C    -0.5424070     0.8439620     0.0000000
C    -0.8433090    -0.5171750     0.0000000
N    -2.2044110    -0.7367500     0.0000000
C    -2.7203200     0.4816220     0.0000000
N    -1.7621040     1.4687250     0.0000000
N     2.9429720     1.0619620     0.0000000
H     2.2894380    -1.3782550     0.0000000
H    -3.7780480     0.7114090     0.0000000
H    -1.9055140     2.4718250     0.0000000
H     3.0816010     2.0608890     0.0000000
H     3.7442410     0.4525250     0.0000000
--
0 1
O     2.7274930     0.0284280    -3.3000000
C     1.5380200    -0.2749280    -3.3000000
N     0.5444810     0.7218990    -3.3000000
C    -0.8331760     0.5747330    -3.3000000
N    -1.2240650    -0.7583250    -3.3000000
C    -0.3429990    -1.8167020    -3.3000000
C     0.9953510    -1.6272180    -3.3000000
O    -1.6271560     1.5061620    -3.3000000
H     1.6879130    -2.4587410    -3.3000000
H     0.8786070     1.6824790    -3.3000000
H    -2.2257760    -0.9095050    -3.3000000
H    -0.7985270    -2.8016270    -3.3000000
units angstrom
}
""")

JSCH_99 = input.process_input("""
molecule dimer {
0 1
O     0.9601320     1.3436400     0.0000000
C     1.5166980     0.2684520     0.0000000
N     0.7573320    -0.9011610     0.0000000
C     1.2481620    -2.1702510     0.0000000
N     2.5209460    -2.4496950     0.0000000
C     3.2915230    -1.3476830     0.0000000
C     2.9121790    -0.0279190     0.0000000
N     4.0200060     0.7969640     0.0000000
C     5.0170310     0.0003310     0.0000000
N     4.6446780    -1.3255770     0.0000000
N     0.3459700    -3.1553460     0.0000000
H    -0.2412520    -0.7659240     0.0000000
H     6.0483360     0.2895830     0.0000000
H     5.2362800    -2.1226110     0.0000000
H     0.6928700    -4.0838600     0.0000000
H    -0.6408270    -2.9885130     0.0000000
--
0 1
O    -0.0130090     1.6513790     3.3600000
C     1.0692420     1.1086750     3.3600000
N     1.1423840    -0.2839060     3.3600000
C     2.2854260    -1.0221180     3.3600000
N     3.4793830    -0.5000700     3.3600000
C     3.4550460     0.8444100     3.3600000
C     2.3724120     1.6891490     3.3600000
N     2.7838090     3.0076580     3.3600000
C     4.0586690     2.9492050     3.3600000
N     4.5367780     1.6576590     3.3600000
N     2.1345620    -2.3493720     3.3600000
H     0.2550220    -0.7614500     3.3600000
H     4.7229940     3.7894010     3.3600000
H     5.4838790     1.3605800     3.3600000
H     2.9609760    -2.8966530     3.3600000
H     1.2381640    -2.7944260     3.3600000
units angstrom
}
""")

JSCH_100 = input.process_input("""
molecule dimer {
0 1
C    -3.0263940    -1.4464050     0.0000000
N    -4.3985350    -1.2378500     0.0000000
C    -4.9449180     0.0000290     0.0000000
C    -4.1674910     1.0873990     0.0000000
C    -2.7399670     0.8551050     0.0000000
N    -2.2307000    -0.3568440     0.0000000
O    -2.6172300    -2.5848080     0.0000000
N    -1.9099420     1.8850830     0.0000000
H    -4.9632880    -2.0564360     0.0000000
H    -6.0175890     0.0492700     0.0000000
H    -4.5763200     2.0777300     0.0000000
H    -2.2565290     2.8138220     0.0000000
H    -0.9141020     1.7329110     0.0000000
--
0 1
C    -1.5982280    -2.9490360     3.3600000
N    -2.8308990    -3.5868360     3.3600000
C    -4.0005400    -2.9065270     3.3600000
C    -4.0107280    -1.5698660     3.3600000
C    -2.7192980    -0.9187180     3.3600000
N    -1.5949260    -1.5998660     3.3600000
O    -0.5980710    -3.6295230     3.3600000
N    -2.6531990     0.4024280     3.3600000
H    -2.8066410    -4.5810390     3.3600000
H    -4.8972920    -3.4971900     3.3600000
H    -4.9235800    -1.0089750     3.3600000
H    -3.4794940     0.9500750     3.3600000
H    -1.7581040     0.8646590     3.3600000
units angstrom
}
""")

JSCH_101 = input.process_input("""
molecule dimer {
0 1
N    -1.3923840    -1.5825730    -0.2790500
C    -1.8533500    -0.3518640    -0.0620430
N    -0.9943890     0.6521290     0.1149880
C    -1.4604570     1.8814980     0.3317590
N    -2.7070820     2.2763020     0.4013740
C    -3.5527210     1.2640760     0.2228910
C    -3.2236500    -0.0504790    -0.0089010
N    -4.3580740    -0.8272780    -0.1458710
C    -5.3247240    -0.0009840    -0.0001730
N    -4.9130980     1.2870000     0.2269330
H    -0.7040060     2.6348130     0.4645890
H    -6.3651290    -0.2529400    -0.0446000
H    -5.4840420     2.0871050     0.3680130
H    -0.4093220    -1.7576030    -0.3099130
H    -2.0356960    -2.3259680    -0.4101310
--
0 1
N    -0.1962490    -2.0987510     2.7709500
C    -1.2925710    -1.3740360     2.9879570
N    -1.1877890    -0.0569040     3.1649880
C    -2.2874510     0.6637280     3.3817590
N    -3.5280520     0.2503840     3.4513740
C    -3.6172170    -1.0655780     3.2728910
C    -2.5783160    -1.9356530     3.0410990
N    -3.0394940    -3.2308940     2.9041290
C    -4.3072140    -3.1305910     3.0498260
N    -4.7312590    -1.8466420     3.2769330
H    -2.1182570     1.7178040     3.5145890
H    -5.0008230    -3.9459620     3.0054000
H    -5.6634530    -1.5349360     3.4180130
H     0.7019450    -1.6625240     2.7400870
H    -0.2797430    -3.0783000     2.6398690
units angstrom
}
""")

JSCH_102 = input.process_input("""
molecule dimer {
0 1
O     1.6803850    -1.8647480     0.3288050
C     2.4435780    -0.9466670     0.1669230
N     1.9754430     0.3257080    -0.0574310
C     2.7234150     1.4521870    -0.2560600
N     4.0753100     1.2353980    -0.2178340
C     4.6340910    -0.0009940     0.0001750
C     3.9044100    -1.0972430     0.1934740
O     2.2509580     2.5354000    -0.4470600
C     4.4733490    -2.4660930     0.4348390
H     4.6436490     2.0381400    -0.3593790
H     0.9698550     0.4501820    -0.0793790
H     5.7079390    -0.0187620     0.0033080
H     4.1432070    -2.8577080     1.3905330
H     4.1417710    -3.1613140    -0.3283340
H     5.5573000    -2.4391850     0.4291940
--
0 1
O     2.4555320    -0.5209070     3.3788050
C     2.5333330     0.6704300     3.2169230
N     1.4067200     1.4246400     2.9925690
C     1.3497150     2.7756270     2.7939400
N     2.5708460     3.3948650     2.8321660
C     3.7496420     2.7230470     3.0501750
C     3.8036770     1.4072660     3.2434740
O     0.3307920     3.3742620     2.6029400
C     5.0685490     0.6342580     3.4848390
H     2.5588020     4.3783590     2.6906210
H     0.5200190     0.9342720     2.9706210
H     4.6288470     3.3398640     3.0533080
H     5.0316430     0.1233820     4.4405330
H     5.2089370    -0.1230850     2.7216660
H     5.9296670     1.2931580     3.4791940
units angstrom
}
""")

JSCH_103 = input.process_input("""
molecule dimer {
0 1
C     2.4313070     1.6249990    -1.4530130
N     3.8007370     1.6249990    -1.6786800
C     4.7029040     1.6249990    -0.6702290
C     4.2995430     1.6249990     0.6041590
C     2.8701050     1.6249990     0.8243640
N     2.0112500     1.6249990    -0.1708960
O     1.6903830     1.6249990    -2.4092590
N     2.3989850     1.6249990     2.0604220
H     4.0848920     1.6249990    -2.6317200
H     5.7382910     1.6249990    -0.9548720
H     4.9943920     1.6249990     1.4196840
H     3.0156050     1.6249990     2.8366040
H     1.4048620     1.6249990     2.2234300
--
0 1
O     0.4979320    -1.6249990     1.9422390
C     1.3595090    -1.6249990     1.0916630
N     0.9987390    -1.6249990    -0.2553610
C     1.8577170    -1.6249990    -1.3106620
N     3.1545590    -1.6249990    -1.1831180
C     3.5468800    -1.6249990     0.1030790
C     2.7782730    -1.6249990     1.2410250
N     3.5769760    -1.6249990     2.3678730
C     4.7713760    -1.6249990     1.9183280
N     4.8269760    -1.6249990     0.5422510
N     1.3040920    -1.6249990    -2.5263360
H     0.0072390    -1.6249990    -0.4353230
H     5.6628210    -1.6249990     2.5121130
H     5.6359190    -1.6249990    -0.0329580
H     1.9209400    -1.6249990    -3.3022070
H     0.3140390    -1.6249990    -2.6726050
units angstrom
}
""")

JSCH_104 = input.process_input("""
molecule dimer {
0 1
C    -3.0263940    -1.4464050     0.0000000
N    -4.3985350    -1.2378500     0.0000000
C    -4.9449180     0.0000290     0.0000000
C    -4.1674910     1.0873990     0.0000000
C    -2.7399670     0.8551050     0.0000000
N    -2.2307000    -0.3568440     0.0000000
O    -2.6172300    -2.5848080     0.0000000
N    -1.9099420     1.8850830     0.0000000
H    -4.9632880    -2.0564360     0.0000000
H    -6.0175890     0.0492700     0.0000000
H    -4.5763200     2.0777300     0.0000000
H    -2.2565290     2.8138220     0.0000000
H    -0.9141020     1.7329110     0.0000000
--
0 1
O    -1.5665350     0.5226760     3.1900000
C    -1.3848270    -0.6743110     3.1900000
N    -0.0830050    -1.1742030     3.1900000
C     0.2658570    -2.4894210     3.1900000
N    -0.5995930    -3.4636200     3.1900000
C    -1.8707500    -3.0250070     3.1900000
C    -2.3395920    -1.7343230     3.1900000
N    -3.7206970    -1.7181430     3.1900000
C    -4.0590580    -2.9486690     3.1900000
N    -2.9784690    -3.8024880     3.1900000
N     1.5747700    -2.7560840     3.1900000
H     0.6453760    -0.4778410     3.1900000
H    -5.0634190    -3.3208450     3.1900000
H    -2.9886000    -4.7950370     3.1900000
H     1.8398890    -3.7111710     3.1900000
H     2.2750440    -2.0410890     3.1900000
units angstrom
}
""")

JSCH_105 = input.process_input("""
molecule dimer {
0 1
N     1.0423840    -1.6008720     0.1400580
C     1.5033500    -0.3559320     0.0311400
N     0.6443890     0.6596690    -0.0577140
C     1.1104570     1.9032530    -0.1665130
N     2.3570820     2.3026230    -0.2014530
C     3.2027210     1.2786930    -0.1118710
C     2.8736500    -0.0510630     0.0044670
N     4.0080740    -0.8368430     0.0732140
C     4.9747240    -0.0009950     0.0000870
N     4.5630980     1.3018810    -0.1139000
H     0.3540060     2.6652780    -0.2331820
H     6.0151290    -0.2558650     0.0223850
H     5.1340420     2.1112380    -0.1847090
H     0.0593220    -1.7779260     0.1555480
H     1.6856960    -2.3528630     0.2058490
--
0 1
O    -0.0504540    -1.6178530    -3.0328940
C     1.0293920    -1.0784590    -3.1266030
N     1.0999170     0.3105210    -3.2285410
C     2.2401210     1.0448270    -3.3391500
N     3.4334530     0.5219170    -3.3635050
C     3.4115810    -0.8191700    -3.2674580
C     2.3318990    -1.6598460    -3.1524330
N     2.7451410    -2.9758150    -3.0805400
C     4.0182190    -2.9198150    -3.1499710
N     4.4933620    -1.6323510    -3.2655320
N     2.0870530     2.3690480    -3.4250070
H     0.2128580     0.7884810    -3.2167550
H     4.6831900    -3.7591200    -3.1247610
H     5.4386800    -1.3377250    -3.3349980
H     2.9113910     2.9134700    -3.5059320
H     1.1910290     2.8146150    -3.4104660
units angstrom
}
""")

JSCH_106 = input.process_input("""
molecule dimer {
0 1
O    -2.0303850    -1.8863100    -0.1650310
C    -2.7935780    -0.9576130    -0.0837800
N    -2.3254430     0.3294740     0.0288250
C    -3.0734150     1.4689780     0.1285190
N    -4.4253100     1.2496820     0.1093330
C    -4.9840910    -0.0010050    -0.0000880
C    -4.2544100    -1.1099300    -0.0971060
O    -2.6009580     2.5647160     0.2243840
C    -4.8233490    -2.4946080    -0.2182500
H    -4.9936490     2.0617070     0.1803760
H    -1.3198550     0.4553880     0.0398410
H    -6.0579390    -0.0189790    -0.0016600
H    -4.4932070    -3.1202310     0.6035230
H    -4.4917710    -2.9686160    -1.1353540
H    -5.9073000    -2.4671550    -0.2167380
--
0 1
C    -1.6419140     2.9739730    -3.0239370
N    -2.8741190     3.6124140    -3.0421140
C    -4.0409900     2.9359160    -3.1500030
C    -4.0487470     1.6026030    -3.2447730
C    -2.7578360     0.9507400    -3.2245270
N    -1.6361750     1.6281560    -3.1188990
O    -0.6443040     3.6509540    -2.9247190
N    -2.6894340    -0.3672360    -3.3142960
H    -2.8516920     4.6040980    -2.9707700
H    -4.9376330     3.5267310    -3.1542940
H    -4.9593830     1.0447610    -3.3310860
H    -3.5136510    -0.9120230    -3.3952410
H    -1.7946790    -0.8299350    -3.3010330
units angstrom
}
""")

JSCH_107 = input.process_input("""
molecule dimer {
0 1
N     1.0423840    -1.6030730     0.1120980
C     1.5033500    -0.3564220     0.0249230
N     0.6443890     0.6605760    -0.0461920
C     1.1104570     1.9058690    -0.1332710
N     2.3570820     2.3057880    -0.1612360
C     3.2027210     1.2804500    -0.0895380
C     2.8736500    -0.0511330     0.0035760
N     4.0080740    -0.8379940     0.0585980
C     4.9747240    -0.0009970     0.0000700
N     4.5630980     1.3036710    -0.0911620
H     0.3540060     2.6689420    -0.1866310
H     6.0151290    -0.2562160     0.0179160
H     5.1340420     2.1141400    -0.1478350
H     0.0593220    -1.7803700     0.1244960
H     1.6856960    -2.3560970     0.1647540
--
0 1
O     1.5241600    -0.5494170     3.3837280
C     1.3439910     0.6454500     3.3087260
N     0.0438450     1.1430380     3.2271380
C    -0.3032010     2.4557550     3.1386110
N     0.5626500     3.4294030     3.1191180
C     1.8322280     2.9929620     3.1959900
C     2.2991810     1.7048790     3.2880520
N     3.6791050     1.6903250     3.3455930
C     4.0186060     2.9192810     3.2900230
N     2.9399160     3.7704860     3.1975320
N    -1.6107030     2.7204770     3.0698940
H    -0.6847300     0.4469420     3.2365720
H     5.0225530     3.2920270     3.3102000
H     2.9511880     4.7614640     3.1419340
H    -1.8744930     3.6737330     3.0051240
H    -2.3112160     2.0058100     3.0815320
units angstrom
}
""")

JSCH_108 = input.process_input("""
molecule dimer {
0 1
O    -2.0303850    -1.8889030    -0.1320850
C    -2.7935780    -0.9589290    -0.0670550
N    -2.3254430     0.3299270     0.0230710
C    -3.0734150     1.4709970     0.1028620
N    -4.4253100     1.2514000     0.0875060
C    -4.9840910    -0.0010070    -0.0000700
C    -4.2544100    -1.1114560    -0.0777210
O    -2.6009580     2.5682420     0.1795890
C    -4.8233490    -2.4980370    -0.1746800
H    -4.9936490     2.0645410     0.1443670
H    -1.3198550     0.4560130     0.0318880
H    -6.0579390    -0.0190050    -0.0013290
H    -4.4932070    -3.1092230     0.6578860
H    -4.4917710    -2.9879780    -1.0833720
H    -5.9073000    -2.4705620    -0.1736470
--
0 1
C    -3.3369590    -0.6409430     3.3908960
N    -4.3247580    -1.6157810     3.3763480
C    -4.0409560    -2.9359630     3.2899980
C    -2.7744220    -3.3565600     3.2141470
C    -1.7557370    -2.3300110     3.2303510
N    -2.0543620    -1.0525720     3.3148920
O    -3.6734450     0.5183010     3.4703070
N    -0.4803010    -2.6733740     3.1585030
H    -5.2616330    -1.2870980     3.4334500
H    -4.8798930    -3.6062030     3.2865630
H    -2.5244870    -4.3961080     3.1450650
H    -0.2161270    -3.6266280     3.0937180
H     0.2361240    -1.9652240     3.1691180
units angstrom
}
""")

JSCH_109 = input.process_input("""
molecule dimer {
0 1
O     2.0303850    -1.8935150     0.0000000
C     2.7935780    -0.9612710     0.0000000
N     2.3254430     0.3307330     0.0000000
C     3.0734150     1.4745890     0.0000000
N     4.4253100     1.2544560     0.0000000
C     4.9840910    -0.0010090     0.0000000
C     4.2544100    -1.1141700     0.0000000
O     2.6009580     2.5745130     0.0000000
C     4.8233490    -2.5041370     0.0000000
H     4.9936490     2.0695820     0.0000000
H     1.3198550     0.4571270     0.0000000
H     6.0579390    -0.0190510     0.0000000
H     4.4932070    -3.0557570     0.8731720
H     4.4917710    -3.0562720    -0.8723020
H     5.9073000    -2.4766570    -0.0008860
--
0 1
O     1.5260840    -0.5520650     3.1800000
C     1.3443760     0.6449210     3.1800000
N     0.0425540     1.1448140     3.1800000
C    -0.3063080     2.4600320     3.1800000
N     0.5591430     3.4342300     3.1800000
C     1.8302990     2.9956180     3.1800000
C     2.2991410     1.7049340     3.1800000
N     3.6802460     1.6887540     3.1800000
C     4.0186070     2.9192800     3.1800000
N     2.9380180     3.7730980     3.1800000
N    -1.6152210     2.7266950     3.1800000
H    -0.6858270     0.4484520     3.1800000
H     5.0229680     3.2914560     3.1800000
H     2.9481490     4.7656470     3.1800000
H    -1.8803400     3.6817820     3.1800000
H    -2.3154950     2.0117000     3.1800000
units angstrom
}
""")

JSCH_110 = input.process_input("""
molecule dimer {
0 1
N    -1.0423840    -1.6069870     0.0000000
C    -1.5033500    -0.3572920     0.0000000
N    -0.6443890     0.6621890     0.0000000
C    -1.1104570     1.9105230     0.0000000
N    -2.3570820     2.3114180     0.0000000
C    -3.2027210     1.2835770     0.0000000
C    -2.8736500    -0.0512580     0.0000000
N    -4.0080740    -0.8400400     0.0000000
C    -4.9747240    -0.0009990     0.0000000
N    -4.5630980     1.3068540     0.0000000
H    -0.3540060     2.6754590     0.0000000
H    -6.0151290    -0.2568420     0.0000000
H    -5.1340420     2.1193020     0.0000000
H    -0.0593220    -1.7847170     0.0000000
H    -1.6856960    -2.3618500     0.0000000
--
0 1
C    -3.3390300    -0.6380930     3.1800000
N    -4.3265300    -1.6133420     3.1800000
C    -4.0409560    -2.9359630     3.1800000
C    -2.7728650    -3.3587040     3.1800000
C    -1.7545120    -2.3316960     3.1800000
N    -2.0548730    -1.0518690     3.1800000
O    -3.6771460     0.5233950     3.1800000
N    -0.4776020    -2.6770890     3.1800000
H    -5.2645780    -1.2830450     3.1800000
H    -4.8798220    -3.6063000     3.1800000
H    -2.5215120    -4.4002020     3.1800000
H    -0.2120980    -3.6321740     3.1800000
H     0.2386050    -1.9686390     3.1800000
units angstrom
}
""")

JSCH_111 = input.process_input("""
molecule dimer {
0 1
O     2.0303850    -1.8863100     0.1650310
C     2.7935780    -0.9576130     0.0837800
N     2.3254430     0.3294740    -0.0288250
C     3.0734150     1.4689780    -0.1285190
N     4.4253100     1.2496820    -0.1093330
C     4.9840910    -0.0010050     0.0000880
C     4.2544100    -1.1099300     0.0971060
O     2.6009580     2.5647160    -0.2243840
C     4.8233490    -2.4946080     0.2182500
H     4.9936490     2.0617070    -0.1803760
H     1.3198550     0.4553880    -0.0398410
H     6.0579390    -0.0189790     0.0016600
H     4.4932070    -2.9680270     1.1361760
H     4.4917710    -3.1206680    -0.6026110
H     5.9073000    -2.4673100     0.2149720
--
0 1
O    -0.0504540    -1.6178530    -3.0328940
C     1.0293920    -1.0784590    -3.1266030
N     1.0999170     0.3105210    -3.2285410
C     2.2401210     1.0448270    -3.3391500
N     3.4334530     0.5219170    -3.3635050
C     3.4115810    -0.8191700    -3.2674580
C     2.3318990    -1.6598460    -3.1524330
N     2.7451410    -2.9758150    -3.0805400
C     4.0182190    -2.9198150    -3.1499710
N     4.4933620    -1.6323510    -3.2655320
N     2.0870530     2.3690480    -3.4250070
H     0.2128580     0.7884810    -3.2167550
H     4.6831900    -3.7591200    -3.1247610
H     5.4386800    -1.3377250    -3.3349980
H     2.9113910     2.9134700    -3.5059320
H     1.1910290     2.8146150    -3.4104660
units angstrom
}
""")

JSCH_112 = input.process_input("""
molecule dimer {
0 1
N    -1.0423840    -1.6008720    -0.1400580
C    -1.5033500    -0.3559320    -0.0311400
N    -0.6443890     0.6596690     0.0577140
C    -1.1104570     1.9032530     0.1665130
N    -2.3570820     2.3026230     0.2014530
C    -3.2027210     1.2786930     0.1118710
C    -2.8736500    -0.0510630    -0.0044670
N    -4.0080740    -0.8368430    -0.0732140
C    -4.9747240    -0.0009950    -0.0000870
N    -4.5630980     1.3018810     0.1139000
H    -0.3540060     2.6652780     0.2331820
H    -6.0151290    -0.2558650    -0.0223850
H    -5.1340420     2.1112380     0.1847090
H    -0.0593220    -1.7779260    -0.1555480
H    -1.6856960    -2.3528630    -0.2058490
--
0 1
C    -1.6419140     2.9739730    -3.0239370
N    -2.8741190     3.6124140    -3.0421140
C    -4.0409900     2.9359160    -3.1500030
C    -4.0487470     1.6026030    -3.2447730
C    -2.7578360     0.9507400    -3.2245270
N    -1.6361750     1.6281560    -3.1188990
O    -0.6443040     3.6509540    -2.9247190
N    -2.6894340    -0.3672360    -3.3142960
H    -2.8516920     4.6040980    -2.9707700
H    -4.9376330     3.5267310    -3.1542940
H    -4.9593830     1.0447610    -3.3310860
H    -3.5136510    -0.9120230    -3.3952410
H    -1.7946790    -0.8299350    -3.3010330
units angstrom
}
""")

JSCH_113 = input.process_input("""
molecule dimer {
0 1
N     1.3923840    -1.6008720     0.1400580
C     1.8533500    -0.3559320     0.0311400
N     0.9943890     0.6596690    -0.0577140
C     1.4604570     1.9032530    -0.1665130
N     2.7070820     2.3026230    -0.2014530
C     3.5527210     1.2786930    -0.1118710
C     3.2236500    -0.0510630     0.0044670
N     4.3580740    -0.8368430     0.0732140
C     5.3247240    -0.0009950     0.0000870
N     4.9130980     1.3018810    -0.1139000
H     0.7040060     2.6652780    -0.2331820
H     6.3651290    -0.2558650     0.0223850
H     5.4840420     2.1112380    -0.1847090
H     0.4093220    -1.7779260     0.1555480
H     2.0356960    -2.3528630     0.2058490
--
0 1
O     2.4682050    -0.5383510     3.4250310
C     2.5397670     0.6615740     3.3437800
N     1.4045070     1.4276870     3.2311750
C     1.3398450     2.7892110     3.1314810
N     2.5624500     3.4064220     3.1506670
C     3.7496490     2.7230370     3.2600880
C     3.8111350     1.3970020     3.3571060
O     0.3135610     3.3979790     3.0356160
C     5.0853090     0.6111890     3.4782500
H     2.5449500     4.3974240     3.0796240
H     0.5169590     0.9384830     3.2201590
H     4.6289750     3.3396890     3.2616600
H     5.0964880     0.0341320     4.3961760
H     5.1850460    -0.0902010     2.6573890
H     5.9461980     1.2704040     3.4749720
units angstrom
}
""")

JSCH_114 = input.process_input("""
molecule dimer {
0 1
N    -1.4867430     1.6920980    -2.3336600
C    -1.5399110     1.6049230    -1.0055780
N    -0.4087210     1.5338080    -0.3037890
C    -0.4671620     1.4467290     1.0245780
N    -1.5291910     1.4187640     1.7901520
C    -2.6502880     1.4904620     1.0763140
C    -2.7488050     1.5835760    -0.2917850
N    -4.0708590     1.6385980    -0.6895780
C    -4.7315520     1.5800700     0.4051650
N    -3.9369080     1.4888380     1.5187780
H     0.4880690     1.3933690     1.5165470
H    -5.7999030     1.5979160     0.4839400
H    -4.2294590     1.4321650     2.4660110
H    -0.6065830     1.7044960    -2.8060620
H    -2.3312660     1.7447540    -2.8510340
--
0 1
O    -1.3473090    -1.4479140    -0.7794320
C    -2.3605260    -1.5129430    -0.1308130
N    -2.3135810    -1.6030670     1.2396240
C    -3.3775550    -1.6828570     2.0937100
N    -4.5954240    -1.6675010     1.4671030
C    -4.7398420    -1.5799270     0.1033200
C    -3.7027260    -1.5022770    -0.7272960
O    -3.2672880    -1.7595820     3.2832490
C    -3.8153430    -1.4053200    -2.2218250
H    -5.3872210    -1.7243610     2.0648200
H    -1.3961730    -1.6118840     1.6702820
H    -5.7555700    -1.5786680    -0.2456340
H    -3.3109870    -2.2369830    -2.7010640
H    -3.3501360    -0.4957970    -2.5852230
H    -4.8546930    -1.4081210    -2.5307720
units angstrom
}
""")

JSCH_115 = input.process_input("""
molecule dimer {
0 1
N    -1.3923840    -1.6069870     0.0000000
C    -1.8533500    -0.3572920     0.0000000
N    -0.9943890     0.6621890     0.0000000
C    -1.4604570     1.9105230     0.0000000
N    -2.7070820     2.3114180     0.0000000
C    -3.5527210     1.2835770     0.0000000
C    -3.2236500    -0.0512580     0.0000000
N    -4.3580740    -0.8400400     0.0000000
C    -5.3247240    -0.0009990     0.0000000
N    -4.9130980     1.3068540     0.0000000
H    -0.7040060     2.6754590     0.0000000
H    -6.3651290    -0.2568420     0.0000000
H    -5.4840420     2.1193020     0.0000000
H    -0.4093220    -1.7847170     0.0000000
H    -2.0356960    -2.3618500     0.0000000
--
0 1
N    -0.1818990    -2.1185030     3.2400000
C    -1.2893810    -1.3784270     3.2400000
N    -1.1937020    -0.0487650     3.2400000
C    -2.3045120     0.6872100     3.2400000
N    -3.5486930     0.2787930     3.2400000
C    -3.6286790    -1.0498020     3.2400000
C    -2.5778590    -1.9362830     3.2400000
N    -3.0319930    -3.2412190     3.2400000
C    -4.3072050    -3.1306030     3.2400000
N    -4.7429290    -1.8305800     3.2400000
H    -2.1421480     1.7506870     3.2400000
H    -4.9985290    -3.9491190     3.2400000
H    -5.6823780    -1.5088880     3.2400000
H     0.7178820    -1.6844600     3.2400000
H    -0.2586520    -3.1073290     3.2400000
units angstrom
}
""")

JSCH_116 = input.process_input("""
molecule dimer {
0 1
O     1.6803850    -1.8935150     0.0000000
C     2.4435780    -0.9612710     0.0000000
N     1.9754430     0.3307330     0.0000000
C     2.7234150     1.4745890     0.0000000
N     4.0753100     1.2544560     0.0000000
C     4.6340910    -0.0010090     0.0000000
C     3.9044100    -1.1141700     0.0000000
O     2.2509580     2.5745130     0.0000000
C     4.4733490    -2.5041370     0.0000000
H     4.6436490     2.0695820     0.0000000
H     0.9698550     0.4571270     0.0000000
H     5.7079390    -0.0190510     0.0000000
H     4.1432070    -3.0557570     0.8731720
H     4.1417710    -3.0562720    -0.8723020
H     5.5573000    -2.4766570    -0.0008860
--
0 1
O     2.4724400    -0.5441800     3.2400000
C     2.5419170     0.6586150     3.2400000
N     1.4037670     1.4287050     3.2400000
C     1.3365470     2.7937510     3.2400000
N     2.5596440     3.4102840     3.2400000
C     3.7496510     2.7230340     3.2400000
C     3.8136270     1.3935720     3.2400000
O     0.3078020     3.4059050     3.2400000
C     5.0909100     0.6034800     3.2400000
H     2.5403210     4.4037960     3.2400000
H     0.5159370     0.9398900     3.2400000
H     4.6290170     3.3396300     3.2400000
H     5.1480540    -0.0368430     4.1131720
H     5.1471940    -0.0381040     2.3676980
H     5.9516930     1.2628420     3.2391140
units angstrom
}
""")

JSCH_117 = input.process_input("""
molecule dimer {
0 1
C    12.1619966    21.5469940    -0.5249999
N    12.0019966    20.1249944    -0.3349999
C    12.9959964    19.1989946    -0.1290000
N    12.5899965    17.9429950    -0.1260000
C    11.2289969    18.0629949    -0.3469999
C    10.2259971    17.0909952    -0.4599999
N    10.4079971    15.7719956    -0.3739999
N     8.9619975    17.5199951    -0.6819998
C     8.7349976    18.8509947    -0.7899998
N     9.6049973    19.8469944    -0.7019998
C    10.8559970    19.3909946    -0.4999999
H    12.8450824    21.9515608     0.2257099
H    12.5490085    21.7744749    -1.5236356
H    11.1843859    22.0177918    -0.4120399
H    14.0220821    19.5129525     0.0161520
H    11.3436468    15.4109067    -0.2800629
H     9.6382753    15.1406078    -0.5991948
H     7.6909448    19.1156876    -0.9420537
--
0 1
C     8.5479976    21.7979939     2.3959993
N     9.1919974    20.5259942     2.6589993
C     8.4229976    19.3799946     2.5429993
O     7.2269980    19.3959946     2.3429993
N     9.0979975    18.2049949     2.7069992
C    10.4579971    18.0869949     2.9379992
O    10.9519969    16.9699952     3.0289992
C    11.2079969    19.3189946     3.0599991
C    12.6759964    19.2659946     3.3619991
C    10.5419970    20.4719943     2.8979992
H     7.4741299    21.6651819     2.5133333
H     8.9049615    22.5495287     3.1049871
H     8.7503455    22.1445498     1.3760436
H    11.0339909    21.4374260     2.9618352
H    13.2133913    18.6878638     2.6029743
H    13.1061963    20.2701373     3.4050200
H    12.8619664    18.7673097     4.3193848
H     8.5371916    17.3217571     2.6353613
units angstrom
}
""")

JSCH_118 = input.process_input("""
molecule dimer {
0 1
N    10.3469971    14.4959959     8.8169975
C    11.5789968    13.8469961     8.7069976
O    11.6019967    12.6419965     8.4119976
N    12.6939964    14.5549959     8.8809975
C    12.6739964    15.9259955     9.1859974
N    13.8309961    16.5099954     9.3349974
C    11.4219968    16.5639954     9.2669974
C    10.3209971    15.8539956     9.0929975
H     9.3699974    16.4009954     9.1789974
H    11.3019968    17.6379951     9.4699973
H    14.6739959    15.9769955     9.2609974
H    13.8749961    17.4909951     9.5239973
C     9.1059774    13.7460371     8.6280336
H     9.4001314    12.7260934     8.3864956
H     8.5051816    13.7537151     9.5428113
H     8.5206636    14.1698120     7.8064238
--
0 1
C    10.7049970     9.6579973    11.8009967
N    11.0689969    11.0699969    11.9839966
C    10.2199971    12.1419966    11.9589966
N    10.8209970    13.3089963    12.1399966
C    12.1439966    12.9639964    12.2549966
C    13.3189963    13.7529961    12.4509965
O    13.3749963    14.9839958    12.5499965
N    14.4609959    13.0409963    12.5269965
C    14.5119959    11.6719967    12.4369965
N    15.7519956    11.1639969    12.5359965
N    13.4609962    10.8809970    12.2549966
C    12.3209965    11.5909968    12.1779966
H    11.6087247     9.0642815    11.9411017
H    10.3130781     9.4887283    10.7941210
H     9.9552752     9.3611644    12.5389945
H    15.3408647    13.5779012    12.6455145
H     9.1538724    12.0114576    11.8260867
H    15.8197976    10.1594152    12.5501065
H    16.5616854    11.7259467    12.8207994
units angstrom
}
""")

JSCH_119 = input.process_input("""
molecule dimer {
0 1
N    10.9240000    16.7550000     5.5620000
C    11.6470000    17.8510000     5.8140000
N    12.9490000    17.6590000     5.9790000
C    13.0500000    16.2780000     5.7950000
C    14.1950000    15.4230000     5.8560000
N    15.4060000    15.8590000     6.0610000
N    13.9020000    14.1180000     5.6250000
C    12.6770000    13.6430000     5.3990000
N    11.5490000    14.4040000     5.3300000
C    11.8450000    15.6910000     5.5460000
H    11.1804230    18.8265530     5.8822870
H    12.5884030    12.5696370     5.2620740
H    16.1977530    15.2199420     5.9750360
H    15.5570940    16.8510580     6.1500010
C     9.4931860    16.6413650     5.3399050
H     9.0446590    17.6337380     5.4112840
H     9.2947180    16.2234190     4.3499330
H     9.0442270    15.9854440     6.0897950
--
0 1
C     9.1690000    13.6920000     8.6010000
N    10.3470000    14.4960000     8.8170000
C    11.5790000    13.8470000     8.7070000
O    11.6020000    12.6420000     8.4120000
N    12.6940000    14.5550000     8.8810000
C    12.6740000    15.9260000     9.1860000
N    13.8310000    16.5100000     9.3350000
C    11.4220000    16.5640000     9.2670000
C    10.3210000    15.8540000     9.0930000
H     9.1403680    12.8642760     9.3131620
H     8.2785600    14.3117950     8.7260530
H     9.1795130    13.2651190     7.5953140
H    11.3501160    17.6252970     9.4808030
H     9.3300790    16.2918180     9.1491660
H    14.7113690    15.9651740     9.2135180
H    13.8876420    17.4962710     9.5342540
units angstrom
}
""")

JSCH_120 = input.process_input("""
molecule dimer {
0 1
N    16.2460000     9.7810000     5.9650000
C    17.5950000    10.0510000     5.9930000
C    18.0920000    11.2690000     5.9020000
C    17.1390000    12.3410000     5.7640000
O    17.4920000    13.5330000     5.6630000
N    15.8280000    12.0550000     5.7130000
C    15.3100000    10.7970000     5.7960000
O    14.1120000    10.5770000     5.7580000
H    18.2280000     9.1744860     6.1031120
C    19.5529600    11.6051630     5.9357380
H    20.1631860    10.7042230     6.0438290
H    19.7760320    12.2828240     6.7658180
H    19.8526100    12.1260780     5.0209680
H    15.1383860    12.8499570     5.6472680
C    15.7717470     8.4029560     6.0779300
H    14.6864640     8.4223240     6.0045990
H    16.1825380     7.7884380     5.2708940
H    16.0652090     7.9755790     7.0417370
--
0 1
C    18.8920000     9.6580000     9.7710000
N    18.5280000    11.0700000     9.5880000
C    19.3770000    12.1420000     9.6130000
N    18.7760000    13.3090000     9.4320000
C    17.4530000    12.9640000     9.3170000
C    16.2780000    13.7530000     9.1210000
O    16.2220000    14.9840000     9.0220000
N    15.1360000    13.0410000     9.0450000
C    15.0850000    11.6720000     9.1350000
N    13.8450000    11.1640000     9.0360000
N    16.1360000    10.8810000     9.3170000
C    17.2760000    11.5910000     9.3940000
H    14.2561290    13.5779040     8.9264920
H    13.0354310    11.7259420     8.7508330
H    13.7773690    10.1594100     9.0211800
H    17.9880060     9.0643740     9.6322420
H    19.2851540     9.4890700    10.7774400
H    19.6407390     9.3607520     9.0321720
H    20.4431070    12.0114850     9.7460660
units angstrom
}
""")

JSCH_121 = input.process_input("""
molecule dimer {
0 1
H     3.1762460     2.3738070     2.9634160
N     2.3770000     1.8470000     3.2830000
C     1.6370000     2.2160000     4.3790000
H     1.9902970     3.0843050     4.9210710
C     0.5610000     1.4930000     4.7730000
H    -0.0085000     1.7736330     5.6470440
C     0.1830000     0.3990000     3.9430000
N    -0.8510000    -0.3400000     4.2540000
H    -1.1799330    -1.0651510     3.5908230
H    -1.4362750    -0.1022370     5.0377650
N     0.8500000     0.0580000     2.8540000
C     1.9550000     0.7640000     2.4990000
O     2.5580000     0.4150000     1.4830000
--
0 1
H     0.0112670     4.2441280     0.3057270
N    -0.1600000     4.2010000     1.2990000
C     0.1490000     5.1520000     2.2350000
H     0.8336150     5.9557770     2.0023890
N    -0.3040000     4.9000000     3.4380000
C    -1.1470000     3.7970000     3.2290000
C    -2.0790000     3.1160000     4.0900000
O    -2.3440000     3.3110000     5.2740000
N    -2.7730000     2.0930000     3.4630000
H    -3.4444620     1.6202680     4.0533010
C    -2.5700000     1.7190000     2.1650000
N    -3.2200000     0.6740000     1.7040000
H    -3.7884800     0.1079360     2.3113460
H    -3.0424470     0.3264300     0.7529310
N    -1.7100000     2.3160000     1.3470000
C    -1.0480000     3.3630000     1.9240000
units angstrom
}
""")

JSCH_122 = input.process_input("""
molecule dimer {
0 1
H    -3.4958570    -1.4150050    -3.9137580
N    -3.0510000    -1.0010000    -3.1090000
C    -3.5590000    -0.8800000    -1.8360000
H    -4.5790060    -1.1582720    -1.6128580
N    -2.7220000    -0.3740000    -0.9680000
C    -1.5590000    -0.1810000    -1.7250000
C    -0.2720000     0.3480000    -1.4650000
N     0.1070000     0.8840000    -0.3230000
H     1.0433330     1.2579620    -0.3065570
H    -0.5751070     1.2407790     0.3499520
N     0.6670000     0.3750000    -2.4130000
C     0.3480000    -0.0810000    -3.6160000
H     1.1321870    -0.0417550    -4.3673920
N    -0.8160000    -0.5790000    -4.0190000
C    -1.7380000    -0.6050000    -3.0150000
--
0 1
H    -1.2611710    -4.7286740    -2.6257100
N    -1.6090000    -4.2940000    -1.7860000
C    -2.7550000    -4.5990000    -1.0690000
H    -3.5136190    -5.2427470    -1.4922410
N    -2.8650000    -3.9860000     0.0730000
C    -1.6740000    -3.2820000     0.1910000
C    -1.1780000    -2.4570000     1.2560000
O    -1.7150000    -2.1460000     2.3170000
N     0.0980000    -1.9830000     1.0200000
H     0.4562670    -1.3045040     1.7132710
C     0.8280000    -2.2730000    -0.0890000
N     2.0180000    -1.7250000    -0.1770000
H     2.3044660    -0.9690820     0.4476800
H     2.5064670    -1.8555350    -1.0472790
N     0.3920000    -3.0250000    -1.1030000
C    -0.8790000    -3.5010000    -0.9150000
units angstrom
}
""")

JSCH_123 = input.process_input("""
molecule dimer {
0 1
H     4.0780890     0.2050200     6.5267380
N     3.3380000    -0.4520000     6.3380000
C     2.1440000    -0.6140000     7.0100000
H     1.9445960    -0.0744500     7.9251340
N     1.3390000    -1.4880000     6.4770000
C     2.0190000    -1.9110000     5.3320000
C     1.6500000    -2.8430000     4.3020000
O     0.6370000    -3.5330000     4.1980000
N     2.5960000    -2.9520000     3.3010000
H     2.3705000    -3.6388980     2.5623150
C     3.7610000    -2.2490000     3.2730000
N     4.5620000    -2.4690000     2.2580000
H     4.3528370    -3.1696290     1.5459440
H     5.4428290    -1.9835850     2.2550440
N     4.1450000    -1.3880000     4.2160000
C     3.2280000    -1.2560000     5.2240000
--
0 1
H     3.1762460     2.3738070     2.9634160
N     2.3770000     1.8470000     3.2830000
C     1.6370000     2.2160000     4.3790000
H     1.9902970     3.0843050     4.9210710
C     0.5610000     1.4930000     4.7730000
H    -0.0085000     1.7736330     5.6470440
C     0.1830000     0.3990000     3.9430000
N    -0.8510000    -0.3400000     4.2540000
H    -1.1799330    -1.0651510     3.5908230
H    -1.4362750    -0.1022370     5.0377650
N     0.8500000     0.0580000     2.8540000
C     1.9550000     0.7640000     2.4990000
O     2.5580000     0.4150000     1.4830000
units angstrom
}
""")

JSCH_124 = input.process_input("""
molecule dimer {
0 1
H    -1.2611710    -4.7286740    -2.6257100
N    -1.6090000    -4.2940000    -1.7860000
C    -2.7550000    -4.5990000    -1.0690000
H    -3.5136190    -5.2427470    -1.4922410
N    -2.8650000    -3.9860000     0.0730000
C    -1.6740000    -3.2820000     0.1910000
C    -1.1780000    -2.4570000     1.2560000
O    -1.7150000    -2.1460000     2.3170000
N     0.0980000    -1.9830000     1.0200000
H     0.4562670    -1.3045040     1.7132710
C     0.8280000    -2.2730000    -0.0890000
N     2.0180000    -1.7250000    -0.1770000
H     2.3044660    -0.9690820     0.4476800
H     2.5064670    -1.8555350    -1.0472790
N     0.3920000    -3.0250000    -1.1030000
C    -0.8790000    -3.5010000    -0.9150000
--
0 1
H     3.2823840    -6.1134940    -1.3105350
N     2.5530000    -6.0070000    -0.6210000
C     1.3990000    -6.7620000    -0.6490000
H     1.3017290    -7.4646550    -1.4662410
C     0.4550000    -6.5890000     0.3070000
H    -0.4593850    -7.1648600     0.2947650
C     0.7210000    -5.6290000     1.3280000
N    -0.1590000    -5.3940000     2.2700000
H    -1.0266130    -5.9017830     2.3125200
H     0.0709100    -4.7127400     3.0149280
N     1.8460000    -4.9310000     1.3860000
C     2.7800000    -5.0940000     0.4140000
O     3.8210000    -4.4400000     0.4780000
units angstrom
}
""")

# <<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:

    GEOS["%s-%s-dimer"      % (dbse, rxn)] = eval("%s_%s" % (dbse, rxn))
    GEOS["%s-%s-monoA-CP"   % (dbse, rxn)] = str(eval("%s_%s" % (dbse, rxn))) + monoA_CP
    GEOS["%s-%s-monoB-CP"   % (dbse, rxn)] = str(eval("%s_%s" % (dbse, rxn))) + monoB_CP
    GEOS["%s-%s-monoA-unCP" % (dbse, rxn)] = str(eval("%s_%s" % (dbse, rxn))) + monoA_unCP
    GEOS["%s-%s-monoB-unCP" % (dbse, rxn)] = str(eval("%s_%s" % (dbse, rxn))) + monoB_unCP
