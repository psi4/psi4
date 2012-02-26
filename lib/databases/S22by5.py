"""
**S22by5**

| Database (Hobza) of interaction energies for dissociation curves of bimolecular complexes.
| Geometries and reference interaction energies from Grafova et al. JCTC 6 2365 (2010).
| Note that the S22by5-N-1.0 members are essentially the same geometries as S22-N (there's trivial round-off error) but the reference interaction energies for S22by5 are of lower quality than those of S22.

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'mol1'`` five-point (0.9, 1.0, 1.2, 1.5, 2.0) :math:`\\times R_{eq}` dissociation curve for molecule 1
  - ...
  - ``'mol22'`` five-point (0.9, 1.0, 1.2, 1.5, 2.0) :math:`\\times R_{eq}` dissociation curve for molecule 22

----

"""
import re
import input

# <<< S22by5 Database Module >>>
dbse = 'S22by5'

# <<< Database Members >>>
mol1  = []
mol2  = []
mol3  = []
mol4  = []
mol5  = []
mol6  = []
mol7  = []
mol8  = []
mol9  = []
mol10 = []
mol11 = []
mol12 = []
mol13 = []
mol14 = []
mol15 = []
mol16 = []
mol17 = []
mol18 = []
mol19 = []
mol20 = []
mol21 = []
mol22 = []
dist = [0.9,1.0,1.2,1.5,2.0]
for d in dist:  mol1.append('1-'  + str(d))
for d in dist:  mol2.append('2-'  + str(d))
for d in dist:  mol3.append('3-'  + str(d))
for d in dist:  mol4.append('4-'  + str(d))
for d in dist:  mol5.append('5-'  + str(d))
for d in dist:  mol6.append('6-'  + str(d))
for d in dist:  mol7.append('7-'  + str(d))
for d in dist:  mol8.append('8-'  + str(d))
for d in dist:  mol9.append('9-'  + str(d))
for d in dist: mol10.append('10-' + str(d))
for d in dist: mol11.append('11-' + str(d))
for d in dist: mol12.append('12-' + str(d))
for d in dist: mol13.append('13-' + str(d))
for d in dist: mol14.append('14-' + str(d))
for d in dist: mol15.append('15-' + str(d))
for d in dist: mol16.append('16-' + str(d))
for d in dist: mol17.append('17-' + str(d))
for d in dist: mol18.append('18-' + str(d))
for d in dist: mol19.append('19-' + str(d))
for d in dist: mol20.append('20-' + str(d))
for d in dist: mol21.append('21-' + str(d))
for d in dist: mol22.append('22-' + str(d))

temp = [mol1, mol2, mol3, mol4,  mol5, mol6, mol7, mol8, mol9, mol10, mol11,
        mol12, mol13, mol14, mol15,  mol16, mol17, mol18, mol19, mol20, mol21, mol22]
HRXN = sum(temp, [])

HRXN_SM = ['1-0.9', '2-1.0', '8-1.5', '16-2.0']
HRXN_LG = ['15-0.9']
HRXN_EQ = []
for m in range(1,23): HRXN_EQ.append(str(m) + '-1.0')

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
BIND['%s-1-0.9'  % (dbse)] =  -2.41
BIND['%s-1-1.0'  % (dbse)] =  -3.14
BIND['%s-1-1.2'  % (dbse)] =  -2.36
BIND['%s-1-1.5'  % (dbse)] =  -1.11
BIND['%s-1-2.0'  % (dbse)] =  -0.36
BIND['%s-2-0.9'  % (dbse)] =  -4.32
BIND['%s-2-1.0'  % (dbse)] =  -4.97
BIND['%s-2-1.2'  % (dbse)] =  -4.04
BIND['%s-2-1.5'  % (dbse)] =  -2.29
BIND['%s-2-2.0'  % (dbse)] =  -0.96
BIND['%s-3-0.9'  % (dbse)] = -16.34
BIND['%s-3-1.0'  % (dbse)] = -18.59
BIND['%s-3-1.2'  % (dbse)] = -15.62
BIND['%s-3-1.5'  % (dbse)] =  -9.24
BIND['%s-3-2.0'  % (dbse)] =  -3.63
BIND['%s-4-0.9'  % (dbse)] = -14.14
BIND['%s-4-1.0'  % (dbse)] = -15.95
BIND['%s-4-1.2'  % (dbse)] = -13.40
BIND['%s-4-1.5'  % (dbse)] =  -8.10
BIND['%s-4-2.0'  % (dbse)] =  -3.51
BIND['%s-5-0.9'  % (dbse)] = -18.73
BIND['%s-5-1.0'  % (dbse)] = -20.46
BIND['%s-5-1.2'  % (dbse)] = -17.16
BIND['%s-5-1.5'  % (dbse)] = -10.46
BIND['%s-5-2.0'  % (dbse)] =  -4.58
BIND['%s-6-0.9'  % (dbse)] = -15.13
BIND['%s-6-1.0'  % (dbse)] = -16.70
BIND['%s-6-1.2'  % (dbse)] = -13.93
BIND['%s-6-1.5'  % (dbse)] =  -8.18
BIND['%s-6-2.0'  % (dbse)] =  -3.26
BIND['%s-7-0.9'  % (dbse)] = -15.02
BIND['%s-7-1.0'  % (dbse)] = -16.37
BIND['%s-7-1.2'  % (dbse)] = -13.30
BIND['%s-7-1.5'  % (dbse)] =  -7.43
BIND['%s-7-2.0'  % (dbse)] =  -2.59
BIND['%s-8-0.9'  % (dbse)] =  -0.34
BIND['%s-8-1.0'  % (dbse)] =  -0.53
BIND['%s-8-1.2'  % (dbse)] =  -0.25
BIND['%s-8-1.5'  % (dbse)] =  -0.06
BIND['%s-8-2.0'  % (dbse)] =  -0.01
BIND['%s-9-0.9'  % (dbse)] =  -0.68
BIND['%s-9-1.0'  % (dbse)] =  -1.48
BIND['%s-9-1.2'  % (dbse)] =  -0.81
BIND['%s-9-1.5'  % (dbse)] =  -0.20
BIND['%s-9-2.0'  % (dbse)] =  -0.03
BIND['%s-10-0.9' % (dbse)] =  -1.09
BIND['%s-10-1.0' % (dbse)] =  -1.50
BIND['%s-10-1.2' % (dbse)] =  -1.13
BIND['%s-10-1.5' % (dbse)] =  -0.48
BIND['%s-10-2.0' % (dbse)] =  -0.12
BIND['%s-11-0.9' % (dbse)] =  -0.15
BIND['%s-11-1.0' % (dbse)] =  -2.81
BIND['%s-11-1.2' % (dbse)] =  -1.92
BIND['%s-11-1.5' % (dbse)] =  -0.53
BIND['%s-11-2.0' % (dbse)] =  -0.07
BIND['%s-12-0.9' % (dbse)] =  -1.69
BIND['%s-12-1.0' % (dbse)] =  -4.51
BIND['%s-12-1.2' % (dbse)] =  -3.02
BIND['%s-12-1.5' % (dbse)] =  -0.98
BIND['%s-12-2.0' % (dbse)] =  -0.19
BIND['%s-13-0.9' % (dbse)] =  -6.76
BIND['%s-13-1.0' % (dbse)] =  -9.87
BIND['%s-13-1.2' % (dbse)] =  -6.26
BIND['%s-13-1.5' % (dbse)] =  -2.42
BIND['%s-13-2.0' % (dbse)] =  -0.69
BIND['%s-14-0.9' % (dbse)] =  -2.13
BIND['%s-14-1.0' % (dbse)] =  -5.18
BIND['%s-14-1.2' % (dbse)] =  -3.61
BIND['%s-14-1.5' % (dbse)] =  -1.08
BIND['%s-14-2.0' % (dbse)] =  -0.10
BIND['%s-15-0.9' % (dbse)] =  -7.99
BIND['%s-15-1.0' % (dbse)] = -12.22
BIND['%s-15-1.2' % (dbse)] =  -8.23
BIND['%s-15-1.5' % (dbse)] =  -3.25
BIND['%s-15-2.0' % (dbse)] =  -0.92
BIND['%s-16-0.9' % (dbse)] =  -1.17
BIND['%s-16-1.0' % (dbse)] =  -1.49
BIND['%s-16-1.2' % (dbse)] =  -1.08
BIND['%s-16-1.5' % (dbse)] =  -0.49
BIND['%s-16-2.0' % (dbse)] =  -0.15
BIND['%s-17-0.9' % (dbse)] =  -3.01
BIND['%s-17-1.0' % (dbse)] =  -3.27
BIND['%s-17-1.2' % (dbse)] =  -2.47
BIND['%s-17-1.5' % (dbse)] =  -1.30
BIND['%s-17-2.0' % (dbse)] =  -0.49
BIND['%s-18-0.9' % (dbse)] =  -2.04
BIND['%s-18-1.0' % (dbse)] =  -2.35
BIND['%s-18-1.2' % (dbse)] =  -1.75
BIND['%s-18-1.5' % (dbse)] =  -0.85
BIND['%s-18-2.0' % (dbse)] =  -0.28
BIND['%s-19-0.9' % (dbse)] =  -4.02
BIND['%s-19-1.0' % (dbse)] =  -4.52
BIND['%s-19-1.2' % (dbse)] =  -3.68
BIND['%s-19-1.5' % (dbse)] =  -2.09
BIND['%s-19-2.0' % (dbse)] =  -0.85
BIND['%s-20-0.9' % (dbse)] =  -2.20
BIND['%s-20-1.0' % (dbse)] =  -2.80
BIND['%s-20-1.2' % (dbse)] =  -2.25
BIND['%s-20-1.5' % (dbse)] =  -1.12
BIND['%s-20-2.0' % (dbse)] =  -0.35
BIND['%s-21-0.9' % (dbse)] =  -4.99
BIND['%s-21-1.0' % (dbse)] =  -5.74
BIND['%s-21-1.2' % (dbse)] =  -4.88
BIND['%s-21-1.5' % (dbse)] =  -2.80
BIND['%s-21-2.0' % (dbse)] =  -1.10
BIND['%s-22-0.9' % (dbse)] =  -6.42
BIND['%s-22-1.0' % (dbse)] =  -7.05
BIND['%s-22-1.2' % (dbse)] =  -5.79
BIND['%s-22-1.5' % (dbse)] =  -3.41
BIND['%s-22-2.0' % (dbse)] =  -1.38

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')
for item in mol1:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-1 Ammonia Dimer at %s Req, C2H' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ammonia Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))

for item in mol2:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-2 Water Dimer at %s Req, CS' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Water Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))

for item in mol3:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-3 Formic Acid Dimer at %s Req, C2H' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Formic Acid Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))

for item in mol4:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-4 Formamide Dimer at %s Req, C2H' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Formamide Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))

for item in mol5:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-5 Uracil Dimer HB at %s Req, C2H' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Uracil Dimer HB at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))

for item in mol6:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-6 2-Pyridone-2-Aminopyridine Complex at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      '2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))

for item in mol7:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'HB-7 Adenine-Thymine Complex WC at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Adenine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Thymine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Adenine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Thymine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))

for item in mol8:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-1 Methane Dimer at %s Req, D3D' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Methane Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))

for item in mol9:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-2 Ethene Dimer at %s Req, D2D' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ethene Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))

for item in mol10:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-3 Benzene-Methane Complex at %s Req, C3' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Methane Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Methane Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Methane from Benzene-Methane Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Methane Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Methane from Benzene-Methane Complex at %s Req' % (molname.group(2))

for item in mol11:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-4 Benzene Dimer Parallel-Disp at %s Req, C2H' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene Dimer PD at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))

for item in mol12:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-6 Pyrazine Dimer at %s Req, CS' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Pyrazine Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))

for item in mol13:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-5 Uracil Dimer Stack at %s Req, C2' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Uracil Dimer Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))

for item in mol14:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-7 Indole-Benzene Complex Stack at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Indole from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Indole from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))

for item in mol15:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-8 Adenine-Thymine Complex Stack at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Adenine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Thymine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Adenine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Thymine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))

for item in mol16:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-1 Ethene-Ethine Complex at %s Req, C2V' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ethene-Ethine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ethene from Ethene-Ethine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ethine from Ethene-Ethine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ethene from Ethene-Ethine Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ethine from Ethene-Ethine Complex at %s Req' % (molname.group(2))

for item in mol17:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-2 Benzene-Water Complex at %s Req, CS' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Water Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Water Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Water from Benzene-Water Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Water Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Water from Benzene-Water Complex at %s Req' % (molname.group(2))

for item in mol18:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-3 Benzene-Ammonia Complex at %s Req, CS' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Ammonia Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ammonia from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ammonia from Benzene-Ammonia Complex at %s Req' % (molname.group(2))

for item in mol19:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-4 Benzene-HCN Complex at %s Req, CS' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-HCN Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-HCN Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'HCN from Benzene-HCN Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-HCN Complex at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'HCN from Benzene-HCN Complex at %s Req' % (molname.group(2))

for item in mol20:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'DD-5 Benzene Dimer T-Shape at %s Req, C2V' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene Dimer T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))

for item in mol21:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-6 Indole-Benzene Complex T-Shape at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Indole from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Indole from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))

for item in mol22:
   molname = rxnpattern.match(item)
   TAGL['%s-%s'            % (dbse, item)] = 'MX-7 Phenol Dimer at %s Req, C1' % (molname.group(2))
   TAGL['%s-%s-dimer'      % (dbse, item)] =      'Phenol Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
   TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

S22by5_1_0p9 = input.process_input("""
molecule dimer {
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.253621272   0.000000000   0.000000000
H    2.860659421  -0.035829274  -0.809565000
H    1.718600721  -0.861570006   0.000000000
H    2.860659421  -0.035829274   0.809565000
units angstrom
}
""")

S22by5_1_1p0 = input.process_input("""
molecule dimer {
0 1
N   -1.578718000  -0.046611000   0.000000000
H   -2.158621000   0.136396000  -0.809565000
H   -2.158621000   0.136396000   0.809565000
H   -0.849471000   0.658193000   0.000000000
--
0 1
N    1.578718000   0.046611000   0.000000000
H    2.158621000  -0.136396000  -0.809565000
H    0.849471000  -0.658193000   0.000000000
H    2.158621000  -0.136396000   0.809565000
units angstrom
}
""")

S22by5_1_1p2 = input.process_input("""
molecule dimer {
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.004828362   0.000000000   0.000000000
H    3.611866511  -0.035829274  -0.809565000
H    2.469807811  -0.861570006   0.000000000
H    3.611866511  -0.035829274   0.809565000
units angstrom
}
""")

S22by5_1_1p5 = input.process_input("""
molecule dimer {
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.756035452   0.000000000   0.000000000
H    4.363073601  -0.035829274  -0.809565000
H    3.221014901  -0.861570006   0.000000000
H    4.363073601  -0.035829274   0.809565000
units angstrom
}
""")

S22by5_1_2p0 = input.process_input("""
molecule dimer {
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    5.008047270   0.000000000   0.000000000
H    5.615085419  -0.035829274  -0.809565000
H    4.473026719  -0.861570006   0.000000000
H    5.615085419  -0.035829274   0.809565000
units angstrom
}
""")

S22by5_2_0p9 = input.process_input("""
molecule dimer {
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    1.756426600   0.000000000   0.000000000
H    2.068390928  -0.496847294  -0.758561000
H    2.068390928  -0.496847294   0.758561000
units angstrom
}
""")

S22by5_2_1p0 = input.process_input("""
molecule dimer {
0 1
O   -1.551007000  -0.114520000   0.000000000
H   -1.934259000   0.762503000   0.000000000
H   -0.599677000   0.040712000   0.000000000
--
0 1
O    1.350625000   0.111469000   0.000000000
H    1.680398000  -0.373741000  -0.758561000
H    1.680398000  -0.373741000   0.758561000
units angstrom
}
""")

S22by5_2_1p2 = input.process_input("""
molecule dimer {
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    2.341902133   0.000000000   0.000000000
H    2.653866461  -0.496847294  -0.758561000
H    2.653866461  -0.496847294   0.758561000
units angstrom
}
""")

S22by5_2_1p5 = input.process_input("""
molecule dimer {
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    2.927377666   0.000000000   0.000000000
H    3.239341994  -0.496847294  -0.758561000
H    3.239341994  -0.496847294   0.758561000
units angstrom
}
""")

S22by5_2_2p0 = input.process_input("""
molecule dimer {
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    3.903170222   0.000000000   0.000000000
H    4.215134550  -0.496847294  -0.758561000
H    4.215134550  -0.496847294   0.758561000
units angstrom
}
""")

S22by5_3_0p9 = input.process_input("""
molecule dimer {
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    2.186205474  -1.011821594   0.000000000
O    1.746270742  -2.250342236   0.000000000
O    1.503291911   0.000000000   0.000000000
H    3.278921791  -0.971514961   0.000000000
H    0.751261211  -2.248465543   0.000000000
units angstrom
}
""")

S22by5_3_1p0 = input.process_input("""
molecule dimer {
0 1
C   -1.888896000  -0.179692000   0.000000000
O   -1.493280000   1.073689000   0.000000000
O   -1.170435000  -1.166590000   0.000000000
H   -2.979488000  -0.258829000   0.000000000
H   -0.498833000   1.107195000   0.000000000
--
0 1
C    1.888896000   0.179692000   0.000000000
O    1.493280000  -1.073689000   0.000000000
O    1.170435000   1.166590000   0.000000000
H    2.979488000   0.258829000   0.000000000
H    0.498833000  -1.107195000   0.000000000
units angstrom
}
""")

S22by5_3_1p2 = input.process_input("""
molecule dimer {
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    2.687302778  -1.011821594   0.000000000
O    2.247368046  -2.250342236   0.000000000
O    2.004389215   0.000000000   0.000000000
H    3.780019095  -0.971514961   0.000000000
H    1.252358515  -2.248465543   0.000000000
units angstrom
}
""")

S22by5_3_1p5 = input.process_input("""
molecule dimer {
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    3.188400082  -1.011821594   0.000000000
O    2.748465350  -2.250342236   0.000000000
O    2.505486519   0.000000000   0.000000000
H    4.281116399  -0.971514961   0.000000000
H    1.753455819  -2.248465543   0.000000000
units angstrom
}
""")

S22by5_3_2p0 = input.process_input("""
molecule dimer {
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    4.023562255  -1.011821594   0.000000000
O    3.583627523  -2.250342236   0.000000000
O    3.340648692   0.000000000   0.000000000
H    5.116278572  -0.971514961   0.000000000
H    2.588617992  -2.248465543   0.000000000
units angstrom
}
""")

S22by5_4_0p9 = input.process_input("""
molecule dimer {
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    3.242982655  -1.316757480   0.000000000
O    2.638862505  -2.387103713   0.000000000
N    2.674136184  -0.100826104   0.000000000
H    3.259710032   0.713812161   0.000000000
H    1.656505975   0.000000000   0.000000000
H    4.343047949  -1.288496220   0.000000000
units angstrom
}
""")

S22by5_4_1p0 = input.process_input("""
molecule dimer {
0 1
C   -2.018649000   0.052883000   0.000000000
O   -1.452200000   1.143634000   0.000000000
N   -1.407770000  -1.142484000   0.000000000
H   -1.964596000  -1.977036000   0.000000000
H   -0.387244000  -1.207782000   0.000000000
H   -3.117061000  -0.013701000   0.000000000
--
0 1
C    2.018649000  -0.052883000   0.000000000
O    1.452200000  -1.143634000   0.000000000
N    1.407770000   1.142484000   0.000000000
H    1.964596000   1.977036000   0.000000000
H    0.387244000   1.207782000   0.000000000
H    3.117061000   0.013701000   0.000000000
units angstrom
}
""")

S22by5_4_1p2 = input.process_input("""
molecule dimer {
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    3.795151314  -1.316757480   0.000000000
O    3.191031164  -2.387103713   0.000000000
N    3.226304843  -0.100826104   0.000000000
H    3.811878691   0.713812161   0.000000000
H    2.208674634   0.000000000   0.000000000
H    4.895216608  -1.288496220   0.000000000
units angstrom
}
""")

S22by5_4_1p5 = input.process_input("""
molecule dimer {
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    4.347319973  -1.316757480   0.000000000
O    3.743199823  -2.387103713   0.000000000
N    3.778473502  -0.100826104   0.000000000
H    4.364047350   0.713812161   0.000000000
H    2.760843293   0.000000000   0.000000000
H    5.447385267  -1.288496220   0.000000000
units angstrom
}
""")

S22by5_4_2p0 = input.process_input("""
molecule dimer {
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    5.267601070  -1.316757480   0.000000000
O    4.663480920  -2.387103713   0.000000000
N    4.698754599  -0.100826104   0.000000000
H    5.284328447   0.713812161   0.000000000
H    3.681124390   0.000000000   0.000000000
H    6.367666364  -1.288496220   0.000000000
units angstrom
}
""")

S22by5_5_0p9 = input.process_input("""
molecule dimer {
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    2.515009084   2.334329905   0.000000000
C    3.179253022   1.297450757   0.000000000
N    2.623672521   0.047940387   0.000000000
C    3.379701020  -1.093192048   0.000000000
C    4.729240680  -1.069579627   0.000000000
C    5.424878943   0.202526014   0.000000000
N    4.549933708   1.305028711   0.000000000
O    6.630530608   0.375595946   0.000000000
H    5.308849416  -1.976469441   0.000000000
H    1.597100890   0.000000000   0.000000000
H    4.984334888   2.217778579   0.000000000
H    2.815046715  -2.013694138   0.000000000
units angstrom
}
""")

S22by5_5_1p0 = input.process_input("""
molecule dimer {
0 1
O   -1.466332000   1.012169000   0.000000000
C   -0.628146000   1.914268000   0.000000000
N    0.720509000   1.688269000   0.000000000
C    1.636729000   2.705276000   0.000000000
C    1.276904000   4.006176000   0.000000000
C   -0.128601000   4.362155000   0.000000000
N   -0.977723000   3.239643000   0.000000000
O   -0.597223000   5.486407000   0.000000000
H    2.010350000   4.793864000   0.000000000
H    1.023251000   0.706182000   0.000000000
H   -1.970027000   3.432385000   0.000000000
H    2.669062000   2.388342000   0.000000000
--
0 1
O    1.466332000  -1.012169000   0.000000000
C    0.628146000  -1.914268000   0.000000000
N   -0.720509000  -1.688269000   0.000000000
C   -1.636729000  -2.705276000   0.000000000
C   -1.276904000  -4.006176000   0.000000000
C    0.128601000  -4.362155000   0.000000000
N    0.977723000  -3.239643000   0.000000000
O    0.597223000  -5.486407000   0.000000000
H   -2.010350000  -4.793864000   0.000000000
H   -1.023251000  -0.706182000   0.000000000
H    1.970027000  -3.432385000   0.000000000
H   -2.669062000  -2.388342000   0.000000000
units angstrom
}
""")

S22by5_5_1p2 = input.process_input("""
molecule dimer {
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    3.047376048   2.334329905   0.000000000
C    3.711619986   1.297450757   0.000000000
N    3.156039485   0.047940387   0.000000000
C    3.912067984  -1.093192048   0.000000000
C    5.261607644  -1.069579627   0.000000000
C    5.957245907   0.202526014   0.000000000
N    5.082300672   1.305028711   0.000000000
O    7.162897572   0.375595946   0.000000000
H    5.841216380  -1.976469441   0.000000000
H    2.129467854   0.000000000   0.000000000
H    5.516701852   2.217778579   0.000000000
H    3.347413679  -2.013694138   0.000000000
units angstrom
}
""")

S22by5_5_1p5 = input.process_input("""
molecule dimer {
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    3.579743012   2.334329905   0.000000000
C    4.243986950   1.297450757   0.000000000
N    3.688406449   0.047940387   0.000000000
C    4.444434948  -1.093192048   0.000000000
C    5.793974608  -1.069579627   0.000000000
C    6.489612871   0.202526014   0.000000000
N    5.614667636   1.305028711   0.000000000
O    7.695264536   0.375595946   0.000000000
H    6.373583344  -1.976469441   0.000000000
H    2.661834818   0.000000000   0.000000000
H    6.049068816   2.217778579   0.000000000
H    3.879780643  -2.013694138   0.000000000
units angstrom
}
""")

S22by5_5_2p0 = input.process_input("""
molecule dimer {
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    4.467021284   2.334329905   0.000000000
C    5.131265222   1.297450757   0.000000000
N    4.575684721   0.047940387   0.000000000
C    5.331713220  -1.093192048   0.000000000
C    6.681252880  -1.069579627   0.000000000
C    7.376891143   0.202526014   0.000000000
N    6.501945908   1.305028711   0.000000000
O    8.582542808   0.375595946   0.000000000
H    7.260861616  -1.976469441   0.000000000
H    3.549113090   0.000000000   0.000000000
H    6.936347088   2.217778579   0.000000000
H    4.767058915  -2.013694138   0.000000000
units angstrom
}
""")

S22by5_6_0p9 = input.process_input("""
molecule dimer {
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    1.673493386   0.000000000   0.000000000
C    2.352093429  -1.145324213   0.192591910
C    3.760459273  -1.168677470   0.196637005
C    4.459573002   0.005477083  -0.001723239
C    3.755182987   1.194447664  -0.202961469
C    2.372894041   1.130328028  -0.192845808
H    4.279274134  -2.103975233   0.356345736
H    5.541001766  -0.003103367  -0.001911235
H    4.259765167   2.134632052  -0.364687797
H    1.782114958   2.025258423  -0.349790900
N    1.620216197  -2.272201547   0.435153550
H    2.101618920  -3.145888174   0.315408858
H    0.644520940  -2.270442069   0.133172072
units angstrom
}
""")

S22by5_6_1p0 = input.process_input("""
molecule dimer {
0 1
O   -1.397621000  -1.885837000  -0.367306000
N   -1.464255000   0.364183000   0.019230000
C   -4.185740000   0.369667000   0.036096000
C   -3.483260000   1.578311000   0.250075000
C   -2.117950000   1.530705000   0.233838000
C   -2.077383000  -0.863749000  -0.189941000
C   -3.515603000  -0.805195000  -0.175759000
H   -5.267804000   0.370743000   0.041142000
H   -3.992033000   2.512756000   0.421441000
H   -1.492920000   2.398410000   0.388502000
H   -4.040123000  -1.734845000  -0.337927000
H   -0.426527000   0.361213000   0.007354000
--
0 1
N    1.432762000   0.363970000  -0.015951000
C    2.115420000  -0.780345000   0.168110000
C    3.523759000  -0.801610000   0.154503000
C    4.218590000   0.373578000  -0.052593000
C    3.509971000   1.561501000  -0.244976000
C    2.128014000   1.495332000  -0.217537000
H    4.045921000  -1.736136000   0.307688000
H    5.299943000   0.366601000  -0.066335000
H    4.011092000   2.502431000  -0.413005000
H    1.533988000   2.389384000  -0.367057000
N    1.388312000  -1.908304000   0.419815000
H    1.869471000  -2.781277000   0.294038000
H    0.408907000  -1.907994000   0.130086000
units angstrom
}
""")

S22by5_6_1p2 = input.process_input("""
molecule dimer {
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.231324514   0.000000000   0.000000000
C    2.909924557  -1.145324213   0.192591910
C    4.318290401  -1.168677470   0.196637005
C    5.017404130   0.005477083  -0.001723239
C    4.313014115   1.194447664  -0.202961469
C    2.930725169   1.130328028  -0.192845808
H    4.837105262  -2.103975233   0.356345736
H    6.098832894  -0.003103367  -0.001911235
H    4.817596295   2.134632052  -0.364687797
H    2.339946086   2.025258423  -0.349790900
N    2.178047325  -2.272201547   0.435153550
H    2.659450048  -3.145888174   0.315408858
H    1.202352068  -2.270442069   0.133172072
units angstrom
}
""")

S22by5_6_1p5 = input.process_input("""
molecule dimer {
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.789155642   0.000000000   0.000000000
C    3.467755685  -1.145324213   0.192591910
C    4.876121529  -1.168677470   0.196637005
C    5.575235258   0.005477083  -0.001723239
C    4.870845243   1.194447664  -0.202961469
C    3.488556297   1.130328028  -0.192845808
H    5.394936390  -2.103975233   0.356345736
H    6.656664022  -0.003103367  -0.001911235
H    5.375427423   2.134632052  -0.364687797
H    2.897777214   2.025258423  -0.349790900
N    2.735878453  -2.272201547   0.435153550
H    3.217281176  -3.145888174   0.315408858
H    1.760183196  -2.270442069   0.133172072
units angstrom
}
""")

S22by5_6_2p0 = input.process_input("""
molecule dimer {
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.718874190   0.000000000   0.000000000
C    4.397474233  -1.145324213   0.192591910
C    5.805840077  -1.168677470   0.196637005
C    6.504953806   0.005477083  -0.001723239
C    5.800563791   1.194447664  -0.202961469
C    4.418274845   1.130328028  -0.192845808
H    6.324654938  -2.103975233   0.356345736
H    7.586382570  -0.003103367  -0.001911235
H    6.305145971   2.134632052  -0.364687797
H    3.827495762   2.025258423  -0.349790900
N    3.665597001  -2.272201547   0.435153550
H    4.146999724  -3.145888174   0.315408858
H    2.689901744  -2.270442069   0.133172072
units angstrom
}
""")

S22by5_7_0p9 = input.process_input("""
molecule dimer {
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    4.674076612   0.155627547  -1.128075158
C    5.366947235  -0.031573530   0.039652507
C    4.745331442  -0.213180550   1.225999310
C    3.289690418  -0.205459536   1.237959001
N    2.678823212  -0.008913767   0.013109028
C    3.292432779   0.176239188  -1.205417098
C    5.464603172  -0.419950938   2.517000917
O    2.621308338  -0.362031655   2.261654302
O    2.694203350   0.342506569  -2.253367774
H    5.154382378   0.288458351  -2.002300903
H    1.636966971   0.000000000   0.000000000
H    6.444191927  -0.024779868  -0.049650000
H    5.195022957   0.354841198   3.233018736
H    5.183915029  -1.373098243   2.962397530
H    6.542374655  -0.403617008   2.368385087
units angstrom
}
""")

S22by5_7_1p0 = input.process_input("""
molecule dimer {
0 1
N    0.935015000  -0.027980000  -0.378892000
C    1.673964000  -0.035777000   0.742432000
C    3.074796000  -0.009448000   0.599456000
C    3.564611000   0.019545000  -0.705987000
N    2.853151000   0.025803000  -1.840960000
C    1.549076000   0.001257000  -1.580801000
N    4.088582000  -0.005443000   1.528979000
C    5.182992000   0.025397000   0.787218000
N    4.929487000   0.041240000  -0.556727000
N    1.071618000  -0.076537000   1.939139000
H    0.879444000   0.005026000  -2.431571000
H    6.188259000   0.037554000   1.173882000
H    5.603537000   0.064876000  -1.303681000
H    0.058692000  -0.042376000   2.003918000
H    1.644380000  -0.034739000   2.761916000
--
0 1
N   -3.921173000  -0.000965000  -1.516366000
C   -4.613683000   0.016905000  -0.333652000
C   -3.991739000   0.021935000   0.866334000
C   -2.536137000   0.007465000   0.876672000
N   -1.925648000  -0.011059000  -0.363895000
C   -2.539590000  -0.014947000  -1.596236000
C   -4.710613000   0.041337000   2.173864000
O   -1.867473000   0.011209000   1.912083000
O   -1.941678000  -0.029188000  -2.657378000
H   -4.401717000  -0.003608000  -2.400492000
H   -0.883826000  -0.021617000  -0.378427000
H   -5.690922000   0.026935000  -0.422718000
H   -4.443928000  -0.830257000   2.769566000
H   -4.426706000   0.918618000   2.753026000
H   -5.788397000   0.050553000   2.024728000
units angstrom
}
""")

S22by5_7_1p2 = input.process_input("""
molecule dimer {
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    5.219732269   0.155627547  -1.128075158
C    5.912602892  -0.031573530   0.039652507
C    5.290987099  -0.213180550   1.225999310
C    3.835346075  -0.205459536   1.237959001
N    3.224478869  -0.008913767   0.013109028
C    3.838088436   0.176239188  -1.205417098
C    6.010258829  -0.419950938   2.517000917
O    3.166963995  -0.362031655   2.261654302
O    3.239859007   0.342506569  -2.253367774
H    5.700038035   0.288458351  -2.002300903
H    2.182622628   0.000000000   0.000000000
H    6.989847584  -0.024779868  -0.049650000
H    5.740678614   0.354841198   3.233018736
H    5.729570686  -1.373098243   2.962397530
H    7.088030312  -0.403617008   2.368385087
units angstrom
}
""")

S22by5_7_1p5 = input.process_input("""
molecule dimer {
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    5.765387926   0.155627547  -1.128075158
C    6.458258549  -0.031573530   0.039652507
C    5.836642756  -0.213180550   1.225999310
C    4.381001732  -0.205459536   1.237959001
N    3.770134526  -0.008913767   0.013109028
C    4.383744093   0.176239188  -1.205417098
C    6.555914486  -0.419950938   2.517000917
O    3.712619652  -0.362031655   2.261654302
O    3.785514664   0.342506569  -2.253367774
H    6.245693692   0.288458351  -2.002300903
H    2.728278285   0.000000000   0.000000000
H    7.535503241  -0.024779868  -0.049650000
H    6.286334271   0.354841198   3.233018736
H    6.275226343  -1.373098243   2.962397530
H    7.633685969  -0.403617008   2.368385087
units angstrom
}
""")

S22by5_7_2p0 = input.process_input("""
molecule dimer {
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    6.674814021   0.155627547  -1.128075158
C    7.367684644  -0.031573530   0.039652507
C    6.746068851  -0.213180550   1.225999310
C    5.290427827  -0.205459536   1.237959001
N    4.679560621  -0.008913767   0.013109028
C    5.293170188   0.176239188  -1.205417098
C    7.465340581  -0.419950938   2.517000917
O    4.622045747  -0.362031655   2.261654302
O    4.694940759   0.342506569  -2.253367774
H    7.155119787   0.288458351  -2.002300903
H    3.637704380   0.000000000   0.000000000
H    8.444929336  -0.024779868  -0.049650000
H    7.195760366   0.354841198   3.233018736
H    7.184652438  -1.373098243   2.962397530
H    8.543112064  -0.403617008   2.368385087
units angstrom
}
""")

S22by5_8_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    3.346489810   0.000000000   0.000000000
H    4.435612789  -0.000311014  -0.000000023
H    2.981975165  -0.513105641  -0.888589641
H    2.981975165  -0.513239461   0.888512354
H    2.982274086   1.026226426   0.000077278
units angstrom
}
""")

S22by5_8_1p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.000140000   1.859161000
H   -0.888551000   0.513060000   1.494685000
H    0.888551000   0.513060000   1.494685000
H    0.000000000  -1.026339000   1.494868000
H    0.000000000   0.000089000   2.948284000
--
0 1
C    0.000000000   0.000140000  -1.859161000
H    0.000000000  -0.000089000  -2.948284000
H   -0.888551000  -0.513060000  -1.494685000
H    0.888551000  -0.513060000  -1.494685000
H    0.000000000   1.026339000  -1.494868000
units angstrom
}
""")

S22by5_8_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    4.461986413   0.000000000   0.000000000
H    5.551109392  -0.000311014  -0.000000023
H    4.097471768  -0.513105641  -0.888589641
H    4.097471768  -0.513239461   0.888512354
H    4.097770689   1.026226426   0.000077278
units angstrom
}
""")

S22by5_8_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    5.577483016   0.000000000   0.000000000
H    6.666605995  -0.000311014  -0.000000023
H    5.212968371  -0.513105641  -0.888589641
H    5.212968371  -0.513239461   0.888512354
H    5.213267292   1.026226426   0.000077278
units angstrom
}
""")

S22by5_8_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    7.436644022   0.000000000   0.000000000
H    8.525767001  -0.000311014  -0.000000023
H    7.072129377  -0.513105641  -0.888589641
H    7.072129377  -0.513239461   0.888512354
H    7.072428298   1.026226426   0.000077278
units angstrom
}
""")

S22by5_9_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    3.346399800   0.471925000   0.471925000
C    3.346399800  -0.471925000  -0.471925000
H    2.423413800   0.872422000   0.872422000
H    2.423413800  -0.872422000  -0.872422000
H    4.270596800   0.870464000   0.870464000
H    4.270596800  -0.870464000  -0.870464000
units angstrom
}
""")

S22by5_9_1p0 = input.process_input("""
molecule dimer {
0 1
C   -0.471925000  -0.471925000  -1.859111000
C    0.471925000   0.471925000  -1.859111000
H   -0.872422000  -0.872422000  -0.936125000
H    0.872422000   0.872422000  -0.936125000
H   -0.870464000  -0.870464000  -2.783308000
H    0.870464000   0.870464000  -2.783308000
--
0 1
C   -0.471925000   0.471925000   1.859111000
C    0.471925000  -0.471925000   1.859111000
H   -0.872422000   0.872422000   0.936125000
H    0.872422000  -0.872422000   0.936125000
H   -0.870464000   0.870464000   2.783308000
H    0.870464000  -0.870464000   2.783308000
units angstrom
}
""")

S22by5_9_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    4.461866400   0.471925000   0.471925000
C    4.461866400  -0.471925000  -0.471925000
H    3.538880400   0.872422000   0.872422000
H    3.538880400  -0.872422000  -0.872422000
H    5.386063400   0.870464000   0.870464000
H    5.386063400  -0.870464000  -0.870464000
units angstrom
}
""")

S22by5_9_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    5.577333000   0.471925000   0.471925000
C    5.577333000  -0.471925000  -0.471925000
H    4.654347000   0.872422000   0.872422000
H    4.654347000  -0.872422000  -0.872422000
H    6.501530000   0.870464000   0.870464000
H    6.501530000  -0.870464000  -0.870464000
units angstrom
}
""")

S22by5_9_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    7.436444000   0.471925000   0.471925000
C    7.436444000  -0.471925000  -0.471925000
H    6.513458000   0.872422000   0.872422000
H    6.513458000  -0.872422000  -0.872422000
H    8.360641000   0.870464000   0.870464000
H    8.360641000  -0.870464000  -0.870464000
units angstrom
}
""")

S22by5_10_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    3.452913900  -0.000000069   0.000000000
H    3.816671953   0.838173871  -0.586878053
H    3.816671906   0.089163973   1.019318994
H    2.366964900   0.000000000   0.000000000
H    3.816671841  -0.927338119  -0.432440941
units angstrom
}
""")

S22by5_10_1p0 = input.process_input("""
molecule dimer {
0 1
C    1.393218000   0.036291000  -0.633280000
C    0.728036000  -1.188402000  -0.633302000
C   -0.665180000  -1.224708000  -0.633280000
C   -1.393204000  -0.036297000  -0.633302000
C   -0.728038000   1.188416000  -0.633280000
C    0.665168000   1.224699000  -0.633302000
H    2.474274000   0.064448000  -0.631724000
H    1.292959000  -2.110541000  -0.631740000
H   -1.181323000  -2.175008000  -0.631724000
H   -2.474261000  -0.064465000  -0.631740000
H   -1.292951000   2.110560000  -0.631724000
H    1.181303000   2.175006000  -0.631740000
--
0 1
C    0.000000000   0.000000000   3.082619000
H    0.586878000   0.838174000   3.446377000
H   -1.019319000   0.089164000   3.446377000
H    0.000000000   0.000000000   1.996670000
H    0.432441000  -0.927338000   3.446377000
units angstrom
}
""")

S22by5_10_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    4.241902200  -0.000000069   0.000000000
H    4.605660253   0.838173871  -0.586878053
H    4.605660206   0.089163973   1.019318994
H    3.155953200   0.000000000   0.000000000
H    4.605660141  -0.927338119  -0.432440941
units angstrom
}
""")

S22by5_10_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    5.030890500  -0.000000069   0.000000000
H    5.394648553   0.838173871  -0.586878053
H    5.394648506   0.089163973   1.019318994
H    3.944941500   0.000000000   0.000000000
H    5.394648441  -0.927338119  -0.432440941
units angstrom
}
""")

S22by5_10_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    6.345871000  -0.000000069   0.000000000
H    6.709629053   0.838173871  -0.586878053
H    6.709629006   0.089163973   1.019318994
H    5.259922000   0.000000000   0.000000000
H    6.709628941  -0.927338119  -0.432440941
units angstrom
}
""")

S22by5_11_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    2.759649224   1.244058476   0.000000000
C    3.074628440   0.622134657  -1.206205000
C    3.074628440   0.622134657   1.206205000
C    3.703514278  -0.621699240  -1.206954000
C    4.016269727  -1.244929310   0.000000000
C    3.703514278  -0.621699240   1.206954000
H    2.258979020   2.202462660   0.000000000
H    2.824770156   1.102778154   2.142315000
H    3.948089550  -1.104085746   2.143798000
H    4.505594855  -2.209685917   0.000000000
H    3.948089550  -1.104085746  -2.143798000
H    2.824770156   1.102778154  -2.142315000
units angstrom
}
""")

S22by5_11_1p0 = input.process_input("""
molecule dimer {
0 1
C   -1.047825000  -1.421674000   0.000000000
C   -1.454503000  -0.855446000   1.206205000
C   -1.454503000  -0.855446000  -1.206205000
C   -2.266797000   0.277161000   1.206954000
C   -2.671478000   0.845021000   0.000000000
C   -2.266797000   0.277161000  -1.206954000
H   -1.133853000  -1.292059000  -2.142315000
H   -2.582494000   0.716307000  -2.143798000
H   -3.303042000   1.723270000   0.000000000
H   -2.582494000   0.716307000   2.143798000
H   -1.133853000  -1.292059000   2.142315000
H   -0.406025000  -2.291905000   0.000000000
--
0 1
C    1.047825000   1.421674000   0.000000000
C    1.454503000   0.855446000  -1.206205000
C    1.454503000   0.855446000   1.206205000
C    2.266797000  -0.277161000  -1.206954000
C    2.671478000  -0.845021000   0.000000000
C    2.266797000  -0.277161000   1.206954000
H    0.406025000   2.291905000   0.000000000
H    1.133853000   1.292059000   2.142315000
H    2.582494000  -0.716307000   2.143798000
H    3.303042000  -1.723270000   0.000000000
H    2.582494000  -0.716307000  -2.143798000
H    1.133853000   1.292059000  -2.142315000
units angstrom
}
""")

S22by5_11_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    3.889216135   1.244058476   0.000000000
C    4.204195351   0.622134657  -1.206205000
C    4.204195351   0.622134657   1.206205000
C    4.833081189  -0.621699240  -1.206954000
C    5.145836638  -1.244929310   0.000000000
C    4.833081189  -0.621699240   1.206954000
H    3.388545931   2.202462660   0.000000000
H    3.954337067   1.102778154   2.142315000
H    5.077656461  -1.104085746   2.143798000
H    5.635161766  -2.209685917   0.000000000
H    5.077656461  -1.104085746  -2.143798000
H    3.954337067   1.102778154  -2.142315000
units angstrom
}
""")

S22by5_11_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    5.018783046   1.244058476   0.000000000
C    5.333762262   0.622134657  -1.206205000
C    5.333762262   0.622134657   1.206205000
C    5.962648100  -0.621699240  -1.206954000
C    6.275403549  -1.244929310   0.000000000
C    5.962648100  -0.621699240   1.206954000
H    4.518112842   2.202462660   0.000000000
H    5.083903978   1.102778154   2.142315000
H    6.207223372  -1.104085746   2.143798000
H    6.764728677  -2.209685917   0.000000000
H    6.207223372  -1.104085746  -2.143798000
H    5.083903978   1.102778154  -2.142315000
units angstrom
}
""")

S22by5_11_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    6.901394563   1.244058476   0.000000000
C    7.216373779   0.622134657  -1.206205000
C    7.216373779   0.622134657   1.206205000
C    7.845259617  -0.621699240  -1.206954000
C    8.158015066  -1.244929310   0.000000000
C    7.845259617  -0.621699240   1.206954000
H    6.400724359   2.202462660   0.000000000
H    6.966515495   1.102778154   2.142315000
H    8.089834889  -1.104085746   2.143798000
H    8.647340194  -2.209685917   0.000000000
H    8.089834889  -1.104085746  -2.143798000
H    6.966515495   1.102778154  -2.142315000
units angstrom
}
""")

S22by5_12_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    3.398538200   0.643131999   1.130045000
C    2.862793235  -0.642689433   1.130631000
N    2.589772167  -1.306738847   0.000000000
C    2.862793235  -0.642689433  -1.130631000
C    3.398538200   0.643131999  -1.130045000
N    3.676023139   1.305979850   0.000000000
H    3.609496345   1.152471205   2.061864000
H    2.643057716  -1.147744338   2.062399000
H    2.643057716  -1.147744338  -2.062399000
H    3.609496345   1.152471205  -2.061864000
units angstrom
}
""")

S22by5_12_1p0 = input.process_input("""
molecule dimer {
0 1
C   -1.247189000  -1.171821000  -0.696139000
C   -1.247189000  -1.171821000   0.696139000
N   -0.258951000  -1.723577000   1.414480000
C    0.731533000  -2.265222000   0.696729000
C    0.731533000  -2.265222000  -0.696729000
N   -0.258951000  -1.723577000  -1.414480000
H   -2.063436000  -0.722320000  -1.247280000
H   -2.063436000  -0.722320000   1.247280000
H    1.548800000  -2.712828000   1.247560000
H    1.548800000  -2.712828000  -1.247560000
--
0 1
C   -0.338003000   2.080061000   1.130045000
C    0.854025000   1.359347000   1.130631000
N    1.470179000   0.990760000   0.000000000
C    0.854025000   1.359347000  -1.130631000
C   -0.338003000   2.080061000  -1.130045000
N   -0.952306000   2.452884000   0.000000000
H   -0.810376000   2.364303000   2.061864000
H    1.320858000   1.067061000   2.062399000
H    1.320858000   1.067061000  -2.062399000
H   -0.810376000   2.364303000  -2.061864000
units angstrom
}
""")

S22by5_12_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    4.442367465   0.643131999   1.130045000
C    3.906622500  -0.642689433   1.130631000
N    3.633601432  -1.306738847   0.000000000
C    3.906622500  -0.642689433  -1.130631000
C    4.442367465   0.643131999  -1.130045000
N    4.719852404   1.305979850   0.000000000
H    4.653325610   1.152471205   2.061864000
H    3.686886981  -1.147744338   2.062399000
H    3.686886981  -1.147744338  -2.062399000
H    4.653325610   1.152471205  -2.061864000
units angstrom
}
""")

S22by5_12_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    5.486196730   0.643131999   1.130045000
C    4.950451765  -0.642689433   1.130631000
N    4.677430697  -1.306738847   0.000000000
C    4.950451765  -0.642689433  -1.130631000
C    5.486196730   0.643131999  -1.130045000
N    5.763681669   1.305979850   0.000000000
H    5.697154875   1.152471205   2.061864000
H    4.730716246  -1.147744338   2.062399000
H    4.730716246  -1.147744338  -2.062399000
H    5.697154875   1.152471205  -2.061864000
units angstrom
}
""")

S22by5_12_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    7.225912172   0.643131999   1.130045000
C    6.690167207  -0.642689433   1.130631000
N    6.417146139  -1.306738847   0.000000000
C    6.690167207  -0.642689433  -1.130631000
C    7.225912172   0.643131999  -1.130045000
N    7.503397111   1.305979850   0.000000000
H    7.436870317   1.152471205   2.061864000
H    6.470431688  -1.147744338   2.062399000
H    6.470431688  -1.147744338  -2.062399000
H    7.436870317   1.152471205  -2.061864000
units angstrom
}
""")

S22by5_13_0p9 = input.process_input("""
molecule dimer {
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    3.306699593  -1.293679543   0.176141970
C    3.341937987  -0.778657200  -1.090194030
H    3.585423040  -1.482976305  -1.871437030
C    3.083223911   0.522034140  -1.338280030
H    3.112133763   0.920071815  -2.337796030
C    2.713052753   1.403319766  -0.246380030
O    2.371727953   2.571655559  -0.351837030
N    2.755902070   0.783286382   1.008844970
H    2.453219399   1.342483138   1.797579970
C    2.971118189  -0.551482081   1.292935970
O    2.866596791  -1.034239706   2.404014970
H    3.384676629  -2.285950208   0.331021970
units angstrom
}
""")

S22by5_13_1p0 = input.process_input("""
molecule dimer {
0 1
N    2.011359000  -1.213207000  -0.098067000
C    2.025708000  -0.697180000  -1.364403000
H    2.297521000  -1.391059000  -2.145646000
C    1.714523000   0.591965000  -1.612489000
H    1.727287000   0.990847000  -2.612005000
C    1.308960000   1.457534000  -0.520589000
O    0.920593000   2.611086000  -0.626046000
N    1.376888000   0.839745000   0.734636000
H    1.051804000   1.386223000   1.523371000
C    1.645991000  -0.485211000   1.018727000
O    1.561109000  -0.971806000   2.129806000
H    2.129463000  -2.201505000   0.056813000
--
0 1
N   -2.011359000   1.213207000  -0.098067000
C   -2.025708000   0.697180000  -1.364403000
H   -2.297521000   1.391059000  -2.145646000
C   -1.714523000  -0.591965000  -1.612489000
H   -1.727287000  -0.990847000  -2.612005000
C   -1.308960000  -1.457534000  -0.520589000
O   -0.920593000  -2.611086000  -0.626046000
N   -1.376888000  -0.839745000   0.734636000
H   -1.051804000  -1.386223000   1.523371000
C   -1.645991000   0.485211000   1.018727000
O   -1.561109000   0.971806000   2.129806000
H   -2.129463000   2.201505000   0.056813000
units angstrom
}
""")

S22by5_13_1p2 = input.process_input("""
molecule dimer {
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    4.316297789  -1.293679543   0.176141970
C    4.351536183  -0.778657200  -1.090194030
H    4.595021236  -1.482976305  -1.871437030
C    4.092822107   0.522034140  -1.338280030
H    4.121731959   0.920071815  -2.337796030
C    3.722650949   1.403319766  -0.246380030
O    3.381326149   2.571655559  -0.351837030
N    3.765500266   0.783286382   1.008844970
H    3.462817595   1.342483138   1.797579970
C    3.980716385  -0.551482081   1.292935970
O    3.876194987  -1.034239706   2.404014970
H    4.394274825  -2.285950208   0.331021970
units angstrom
}
""")

S22by5_13_1p5 = input.process_input("""
molecule dimer {
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    5.325895984  -1.293679543   0.176141970
C    5.361134378  -0.778657200  -1.090194030
H    5.604619431  -1.482976305  -1.871437030
C    5.102420302   0.522034140  -1.338280030
H    5.131330154   0.920071815  -2.337796030
C    4.732249144   1.403319766  -0.246380030
O    4.390924344   2.571655559  -0.351837030
N    4.775098461   0.783286382   1.008844970
H    4.472415790   1.342483138   1.797579970
C    4.990314580  -0.551482081   1.292935970
O    4.885793182  -1.034239706   2.404014970
H    5.403873020  -2.285950208   0.331021970
units angstrom
}
""")

S22by5_13_2p0 = input.process_input("""
molecule dimer {
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    7.008559644  -1.293679543   0.176141970
C    7.043798038  -0.778657200  -1.090194030
H    7.287283091  -1.482976305  -1.871437030
C    6.785083962   0.522034140  -1.338280030
H    6.813993814   0.920071815  -2.337796030
C    6.414912804   1.403319766  -0.246380030
O    6.073588004   2.571655559  -0.351837030
N    6.457762121   0.783286382   1.008844970
H    6.155079450   1.342483138   1.797579970
C    6.672978240  -0.551482081   1.292935970
O    6.568456842  -1.034239706   2.404014970
H    7.086536680  -2.285950208   0.331021970
units angstrom
}
""")

S22by5_14_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    2.717766027  -0.578056849   3.494904751
C    2.793508398  -0.571969873   2.415753956
C    2.753054336   0.633650134   1.734349558
H    2.645935858   1.567038531   2.272036098
C    2.855804852   0.624347564   0.333339655
C    2.845637545   1.633662034  -0.673499279
H    2.762013625   2.698030593  -0.533251753
C    2.976224608   0.992808148  -1.884517470
N    3.081930238  -0.360086596  -1.675422891
C    2.997750328  -0.624347564  -0.333339655
C    3.046288127  -1.839842986   0.351754941
H    3.153106953  -2.780217935  -0.172940228
C    2.941516868  -1.796211682   1.733036170
H    2.973148444  -2.718261443   2.297634930
H    3.103876306  -1.056446212  -2.398978775
H    3.012441631   1.398036276  -2.881807744
units angstrom
}
""")

S22by5_14_1p0 = input.process_input("""
molecule dimer {
0 1
C   -0.021074000   1.531861000  -1.363935000
C   -1.274679000   0.974103000  -1.607410000
C   -1.378305000  -0.225698000  -2.308415000
C   -0.228943000  -0.866405000  -2.768794000
C    1.024788000  -0.303517000  -2.531241000
C    1.129000000   0.896679000  -1.829983000
H    0.060074000   2.456563000  -0.809396000
H   -2.165100000   1.465452000  -1.240568000
H   -2.350973000  -0.661612000  -2.492670000
H   -0.310342000  -1.795576000  -3.317270000
H    1.916585000  -0.794084000  -2.899394000
H    2.100035000   1.332676000  -1.640042000
--
0 1
H   -2.941765000   0.895383000   2.223905000
C   -2.022067000   0.425854000   1.901355000
C   -0.814942000   1.074045000   2.106698000
H   -0.785153000   2.044381000   2.585609000
C    0.370429000   0.449285000   1.684746000
C    1.750862000   0.803894000   1.719400000
H    2.187011000   1.699828000   2.127590000
C    2.445136000  -0.231074000   1.135331000
N    1.564646000  -1.213781000   0.755538000
C    0.286121000  -0.826949000   1.061875000
C   -0.928467000  -1.485312000   0.860694000
H   -0.972920000  -2.455485000   0.383401000
C   -2.079285000  -0.841767000   1.287644000
H   -3.038997000  -1.320385000   1.146840000
H    1.807574000  -2.036696000   0.233304000
H    3.502879000  -0.348534000   0.969523000
units angstrom
}
""")

S22by5_14_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    3.693358557  -0.578056849   3.494904751
C    3.769100928  -0.571969873   2.415753956
C    3.728646866   0.633650134   1.734349558
H    3.621528388   1.567038531   2.272036098
C    3.831397382   0.624347564   0.333339655
C    3.821230075   1.633662034  -0.673499279
H    3.737606155   2.698030593  -0.533251753
C    3.951817138   0.992808148  -1.884517470
N    4.057522768  -0.360086596  -1.675422891
C    3.973342858  -0.624347564  -0.333339655
C    4.021880657  -1.839842986   0.351754941
H    4.128699483  -2.780217935  -0.172940228
C    3.917109398  -1.796211682   1.733036170
H    3.948740974  -2.718261443   2.297634930
H    4.079468836  -1.056446212  -2.398978775
H    3.988034161   1.398036276  -2.881807744
units angstrom
}
""")

S22by5_14_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    4.668951087  -0.578056849   3.494904751
C    4.744693458  -0.571969873   2.415753956
C    4.704239396   0.633650134   1.734349558
H    4.597120918   1.567038531   2.272036098
C    4.806989912   0.624347564   0.333339655
C    4.796822605   1.633662034  -0.673499279
H    4.713198685   2.698030593  -0.533251753
C    4.927409668   0.992808148  -1.884517470
N    5.033115298  -0.360086596  -1.675422891
C    4.948935388  -0.624347564  -0.333339655
C    4.997473187  -1.839842986   0.351754941
H    5.104292013  -2.780217935  -0.172940228
C    4.892701928  -1.796211682   1.733036170
H    4.924333504  -2.718261443   2.297634930
H    5.055061366  -1.056446212  -2.398978775
H    4.963626691   1.398036276  -2.881807744
units angstrom
}
""")

S22by5_14_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    6.294938637  -0.578056849   3.494904751
C    6.370681008  -0.571969873   2.415753956
C    6.330226946   0.633650134   1.734349558
H    6.223108468   1.567038531   2.272036098
C    6.432977462   0.624347564   0.333339655
C    6.422810155   1.633662034  -0.673499279
H    6.339186235   2.698030593  -0.533251753
C    6.553397218   0.992808148  -1.884517470
N    6.659102848  -0.360086596  -1.675422891
C    6.574922938  -0.624347564  -0.333339655
C    6.623460737  -1.839842986   0.351754941
H    6.730279563  -2.780217935  -0.172940228
C    6.518689478  -1.796211682   1.733036170
H    6.550321054  -2.718261443   2.297634930
H    6.681048916  -1.056446212  -2.398978775
H    6.589614241   1.398036276  -2.881807744
units angstrom
}
""")

S22by5_15_0p9 = input.process_input("""
molecule dimer {
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    2.995457244   1.318912569   0.115169333
C    3.033773997   0.544134785   1.248235461
H    3.166936649   1.084216460   2.174491246
C    2.913123372  -0.802036026   1.213306349
C    2.965573998  -1.664227788   2.429380731
H    2.009790775  -2.161867438   2.585037720
H    3.726416066  -2.435033978   2.315487569
H    3.189128467  -1.070628980   3.313538183
C    2.718644614  -1.440326451  -0.080379664
O    2.558245305  -2.640081851  -0.255033817
N    2.729839539  -0.560837886  -1.168484485
H    2.554150647  -0.977998743  -2.072617562
C    2.814781928   0.814169728  -1.152798148
O    2.732113465   1.513854058  -2.149163262
H    3.033823338   2.322516737   0.179118562
units angstrom
}
""")

S22by5_15_1p0 = input.process_input("""
molecule dimer {
0 1
N    0.279301000   2.406839000  -0.605752000
C   -1.084857000   2.445746000  -0.551161000
H   -1.659440000   3.023029000  -1.256090000
N   -1.597712000   1.717988000   0.428754000
C   -0.489725000   1.171436000   1.030191000
C   -0.346137000   0.291471000   2.117234000
N   -1.418709000  -0.167777000   2.810144000
H   -1.238875000  -0.959480000   3.404758000
H   -2.291873000  -0.178822000   2.307362000
N    0.885763000  -0.070076000   2.491949000
C    1.935235000   0.407288000   1.796802000
H    2.906033000   0.078841000   2.145818000
N    1.940978000   1.224202000   0.740220000
C    0.695219000   1.577986000   0.406398000
H    0.861007000   2.829804000  -1.310450000
--
0 1
N    1.275461000  -0.647899000  -1.977910000
C    1.413053000  -1.553685000  -0.955067000
H    2.425877000  -1.867078000  -0.746878000
C    0.357598000  -2.023950000  -0.253057000
C    0.482129000  -3.017949000   0.852122000
H    0.175770000  -2.575607000   1.798628000
H   -0.160169000  -3.877041000   0.663950000
H    1.511244000  -3.357277000   0.951366000
C   -0.968471000  -1.529811000  -0.593979000
O   -2.002928000  -1.839696000  -0.019945000
N   -0.995692000  -0.638387000  -1.672042000
H   -1.901406000  -0.250172000  -1.898576000
C    0.068470000  -0.119176000  -2.376376000
O   -0.039788000   0.722701000  -3.253108000
H    2.085329000  -0.276018000  -2.445458000
units angstrom
}
""")

S22by5_15_1p2 = input.process_input("""
molecule dimer {
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    3.951238365   1.318912569   0.115169333
C    3.989555118   0.544134785   1.248235461
H    4.122717770   1.084216460   2.174491246
C    3.868904493  -0.802036026   1.213306349
C    3.921355119  -1.664227788   2.429380731
H    2.965571896  -2.161867438   2.585037720
H    4.682197187  -2.435033978   2.315487569
H    4.144909588  -1.070628980   3.313538183
C    3.674425735  -1.440326451  -0.080379664
O    3.514026426  -2.640081851  -0.255033817
N    3.685620660  -0.560837886  -1.168484485
H    3.509931768  -0.977998743  -2.072617562
C    3.770563049   0.814169728  -1.152798148
O    3.687894586   1.513854058  -2.149163262
H    3.989604459   2.322516737   0.179118562
units angstrom
}
""")

S22by5_15_1p5 = input.process_input("""
molecule dimer {
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    4.907019487   1.318912569   0.115169333
C    4.945336240   0.544134785   1.248235461
H    5.078498892   1.084216460   2.174491246
C    4.824685615  -0.802036026   1.213306349
C    4.877136241  -1.664227788   2.429380731
H    3.921353018  -2.161867438   2.585037720
H    5.637978309  -2.435033978   2.315487569
H    5.100690710  -1.070628980   3.313538183
C    4.630206857  -1.440326451  -0.080379664
O    4.469807548  -2.640081851  -0.255033817
N    4.641401782  -0.560837886  -1.168484485
H    4.465712890  -0.977998743  -2.072617562
C    4.726344171   0.814169728  -1.152798148
O    4.643675708   1.513854058  -2.149163262
H    4.945385581   2.322516737   0.179118562
units angstrom
}
""")

S22by5_15_2p0 = input.process_input("""
molecule dimer {
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    6.499988023   1.318912569   0.115169333
C    6.538304776   0.544134785   1.248235461
H    6.671467428   1.084216460   2.174491246
C    6.417654151  -0.802036026   1.213306349
C    6.470104777  -1.664227788   2.429380731
H    5.514321554  -2.161867438   2.585037720
H    7.230946845  -2.435033978   2.315487569
H    6.693659246  -1.070628980   3.313538183
C    6.223175393  -1.440326451  -0.080379664
O    6.062776084  -2.640081851  -0.255033817
N    6.234370318  -0.560837886  -1.168484485
H    6.058681426  -0.977998743  -2.072617562
C    6.319312707   0.814169728  -1.152798148
O    6.236644244   1.513854058  -2.149163262
H    6.538354117   2.322516737   0.179118562
units angstrom
}
""")

S22by5_16_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    4.749960900   0.000000000   0.000000000
C    3.542697900   0.000000000   0.000000000
H    2.476809900   0.000000000   0.000000000
H    5.813386900   0.000000000   0.000000000
units angstrom
}
""")

S22by5_16_1p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.667578000  -2.124659000
C    0.000000000   0.667578000  -2.124659000
H    0.923621000  -1.232253000  -2.126185000
H   -0.923621000  -1.232253000  -2.126185000
H   -0.923621000   1.232253000  -2.126185000
H    0.923621000   1.232253000  -2.126185000
--
0 1
C    0.000000000   0.000000000   2.900503000
C    0.000000000   0.000000000   1.693240000
H    0.000000000   0.000000000   0.627352000
H    0.000000000   0.000000000   3.963929000
units angstrom
}
""")

S22by5_16_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    5.575564200   0.000000000   0.000000000
C    4.368301200   0.000000000   0.000000000
H    3.302413200   0.000000000   0.000000000
H    6.638990200   0.000000000   0.000000000
units angstrom
}
""")

S22by5_16_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    6.401167500   0.000000000   0.000000000
C    5.193904500   0.000000000   0.000000000
H    4.128016500   0.000000000   0.000000000
H    7.464593500   0.000000000   0.000000000
units angstrom
}
""")

S22by5_16_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    7.777173000   0.000000000   0.000000000
C    6.569910000   0.000000000   0.000000000
H    5.504022000   0.000000000   0.000000000
H    8.840599000   0.000000000   0.000000000
units angstrom
}
""")

S22by5_17_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    3.175061618   0.124369730   0.000000000
H    3.265337861   1.079117991   0.000000000
H    2.221117117   0.000000000   0.000000000
units angstrom
}
""")

S22by5_17_1p0 = input.process_input("""
molecule dimer {
0 1
C    0.780612000  -0.609888000  -1.207543000
C    0.478404000   0.751041000  -1.207904000
C    0.327659000   1.431857000   0.000000000
C    0.478404000   0.751041000   1.207904000
C    0.780612000  -0.609888000   1.207543000
C    0.932151000  -1.289961000   0.000000000
H    0.896669000  -1.137605000  -2.144148000
H    0.357390000   1.278209000  -2.144055000
H    0.091859000   2.487141000   0.000000000
H    0.357390000   1.278209000   2.144055000
H    0.896669000  -1.137605000   2.144148000
H    1.169006000  -2.345167000   0.000000000
--
0 1
O   -2.788527000  -0.274485000   0.000000000
H   -2.622911000  -1.219083000   0.000000000
H   -1.901510000   0.097911000   0.000000000
units angstrom
}
""")

S22by5_17_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    3.915433991   0.124369730   0.000000000
H    4.005710234   1.079117991   0.000000000
H    2.961489490   0.000000000   0.000000000
units angstrom
}
""")

S22by5_17_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    4.655806363   0.124369730   0.000000000
H    4.746082606   1.079117991   0.000000000
H    3.701861862   0.000000000   0.000000000
units angstrom
}
""")

S22by5_17_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    5.889760317   0.124369730   0.000000000
H    5.980036560   1.079117991   0.000000000
H    4.935815816   0.000000000   0.000000000
units angstrom
}
""")

S22by5_18_0p9 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    3.322432676  -0.175158455   0.000000000
H    3.685723470   0.316960994  -0.806073000
H    3.685723470   0.316960994   0.806073000
H    2.324338249   0.000000000   0.000000000
units angstrom
}
""")

S22by5_18_1p0 = input.process_input("""
molecule dimer {
0 1
C   -0.739281000   0.515879000  -1.207108000
C   -1.426144000   0.396545000   0.000000000
C   -0.739281000   0.515879000   1.207108000
C    0.634227000   0.754640000   1.207073000
C    1.321043000   0.873757000   0.000000000
C    0.634227000   0.754640000  -1.207073000
H   -1.271950000   0.420632000  -2.143289000
H   -2.490220000   0.205238000   0.000000000
H   -1.271950000   0.420632000   2.143289000
H    1.166800000   0.847488000   2.143695000
H    2.386359000   1.059631000   0.000000000
H    1.166800000   0.847488000  -2.143695000
--
0 1
N    0.180393000  -2.949123000   0.000000000
H    0.759549000  -3.145948000  -0.806073000
H    0.759549000  -3.145948000   0.806073000
H    0.044417000  -1.944940000   0.000000000
units angstrom
}
""")

S22by5_18_1p2 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    4.097212092  -0.175158455   0.000000000
H    4.460502886   0.316960994  -0.806073000
H    4.460502886   0.316960994   0.806073000
H    3.099117665   0.000000000   0.000000000
units angstrom
}
""")

S22by5_18_1p5 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    4.871991508  -0.175158455   0.000000000
H    5.235282302   0.316960994  -0.806073000
H    5.235282302   0.316960994   0.806073000
H    3.873897081   0.000000000   0.000000000
units angstrom
}
""")

S22by5_18_2p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    6.163290535  -0.175158455   0.000000000
H    6.526581329   0.316960994  -0.806073000
H    6.526581329   0.316960994   0.806073000
H    5.165196108   0.000000000   0.000000000
units angstrom
}
""")

S22by5_19_0p9 = input.process_input("""
molecule dimer {
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    4.308034683   0.304536859   0.000000000
C    3.151543935   0.145763954   0.000000000
H    2.093660645   0.000000000   0.000000000
units angstrom
}
""")

S22by5_19_1p0 = input.process_input("""
molecule dimer {
0 1
C   -0.709774000  -0.990423000   1.207702000
C   -1.406534000  -0.965353000   0.000000000
C   -0.709774000  -0.990423000  -1.207702000
C    0.683965000  -1.040510000  -1.207865000
C    1.380978000  -1.065552000   0.000000000
C    0.683965000  -1.040510000   1.207865000
H   -1.249948000  -0.968628000   2.144051000
H   -2.486920000  -0.923706000   0.000000000
H   -1.249948000  -0.968628000  -2.144051000
H    1.224288000  -1.058075000  -2.144256000
H    2.461589000  -1.102982000   0.000000000
H    1.224288000  -1.058075000   2.144256000
--
0 1
N   -0.003412000   3.535393000   0.000000000
C    0.075196000   2.370704000   0.000000000
H    0.147629000   1.305285000   0.000000000
units angstrom
}
""")

S22by5_19_1p2 = input.process_input("""
molecule dimer {
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    5.005921565   0.304536859   0.000000000
C    3.849430817   0.145763954   0.000000000
H    2.791547527   0.000000000   0.000000000
units angstrom
}
""")

S22by5_19_1p5 = input.process_input("""
molecule dimer {
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    5.703808447   0.304536859   0.000000000
C    4.547317699   0.145763954   0.000000000
H    3.489434409   0.000000000   0.000000000
units angstrom
}
""")

S22by5_19_2p0 = input.process_input("""
molecule dimer {
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    6.866953250   0.304536859   0.000000000
C    5.710462502   0.145763954   0.000000000
H    4.652579212   0.000000000   0.000000000
units angstrom
}
""")

S22by5_20_0p9 = input.process_input("""
molecule dimer {
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    2.189283067   0.000000000  -1.394063000
C    2.189759067   1.207238000  -0.697047000
C    2.189759067   1.207238000   0.697047000
C    2.189283067   0.000000000   1.394063000
C    2.189759067  -1.207238000   0.697047000
C    2.189759067  -1.207238000  -0.697047000
H    2.185453067   0.000000000  -2.475399000
H    2.188807067   2.143565000  -1.238232000
H    2.188807067   2.143565000   1.238232000
H    2.185453067   0.000000000   2.475399000
H    2.188807067  -2.143565000   1.238232000
H    2.188807067  -2.143565000  -1.238232000
units angstrom
}
""")

S22by5_20_1p0 = input.process_input("""
molecule dimer {
0 1
C    0.000000000   0.000000000   1.059035000
C    0.000000000  -1.206008000   1.757674000
C    0.000000000  -1.207177000   3.151591000
C    0.000000000   0.000000000   3.848575000
C    0.000000000   1.207177000   3.151591000
C    0.000000000   1.206008000   1.757674000
H    0.000000000   0.000000000  -0.021580000
H    0.000000000  -2.141639000   1.214422000
H    0.000000000  -2.143566000   3.692995000
H    0.000000000   0.000000000   4.930150000
H    0.000000000   2.143566000   3.692995000
H    0.000000000   2.141639000   1.214422000
--
0 1
C   -1.394063000   0.000000000  -2.454152000
C   -0.697047000   1.207238000  -2.454628000
C    0.697047000   1.207238000  -2.454628000
C    1.394063000   0.000000000  -2.454152000
C    0.697047000  -1.207238000  -2.454628000
C   -0.697047000  -1.207238000  -2.454628000
H   -2.475399000   0.000000000  -2.450322000
H   -1.238232000   2.143565000  -2.453676000
H    1.238232000   2.143565000  -2.453676000
H    2.475399000   0.000000000  -2.450322000
H    1.238232000  -2.143565000  -2.453676000
H   -1.238232000  -2.143565000  -2.453676000
units angstrom
}
""")

S22by5_20_1p2 = input.process_input("""
molecule dimer {
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    2.919149867   0.000000000  -1.394063000
C    2.919625867   1.207238000  -0.697047000
C    2.919625867   1.207238000   0.697047000
C    2.919149867   0.000000000   1.394063000
C    2.919625867  -1.207238000   0.697047000
C    2.919625867  -1.207238000  -0.697047000
H    2.915319867   0.000000000  -2.475399000
H    2.918673867   2.143565000  -1.238232000
H    2.918673867   2.143565000   1.238232000
H    2.915319867   0.000000000   2.475399000
H    2.918673867  -2.143565000   1.238232000
H    2.918673867  -2.143565000  -1.238232000
units angstrom
}
""")

S22by5_20_1p5 = input.process_input("""
molecule dimer {
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    3.649016667   0.000000000  -1.394063000
C    3.649492667   1.207238000  -0.697047000
C    3.649492667   1.207238000   0.697047000
C    3.649016667   0.000000000   1.394063000
C    3.649492667  -1.207238000   0.697047000
C    3.649492667  -1.207238000  -0.697047000
H    3.645186667   0.000000000  -2.475399000
H    3.648540667   2.143565000  -1.238232000
H    3.648540667   2.143565000   1.238232000
H    3.645186667   0.000000000   2.475399000
H    3.648540667  -2.143565000   1.238232000
H    3.648540667  -2.143565000  -1.238232000
units angstrom
}
""")

S22by5_20_2p0 = input.process_input("""
molecule dimer {
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    4.865461333   0.000000000  -1.394063000
C    4.865937333   1.207238000  -0.697047000
C    4.865937333   1.207238000   0.697047000
C    4.865461333   0.000000000   1.394063000
C    4.865937333  -1.207238000   0.697047000
C    4.865937333  -1.207238000  -0.697047000
H    4.861631333   0.000000000  -2.475399000
H    4.864985333   2.143565000  -1.238232000
H    4.864985333   2.143565000   1.238232000
H    4.861631333   0.000000000   2.475399000
H    4.864985333  -2.143565000   1.238232000
H    4.864985333  -2.143565000  -1.238232000
units angstrom
}
""")

S22by5_21_0p9 = input.process_input("""
molecule dimer {
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    2.007797424   0.000000000   0.000000000
N    3.015114828   0.005056388   0.000000000
C    3.796769012   1.132604937   0.000000000
C    5.125653739   0.772354616   0.000000000
C    5.167047225  -0.653193161   0.000000000
C    3.817202589  -1.104920876   0.000000000
C    3.482542920  -2.462094972   0.000000000
C    4.524735226  -3.376178892   0.000000000
C    5.869058665  -2.951641292   0.000000000
C    6.199398544  -1.606705567   0.000000000
H    3.343074787   2.109594763   0.000000000
H    5.961043541   1.451489921   0.000000000
H    2.450153978  -2.785730808   0.000000000
H    4.303017780  -4.434822780   0.000000000
H    6.655123584  -3.694570139   0.000000000
H    7.235724321  -1.294593877   0.000000000
units angstrom
}
""")

S22by5_21_1p0 = input.process_input("""
molecule dimer {
0 1
C    2.511900000   1.625015000   0.000000000
C    2.713009000   0.957854000  -1.208292000
C    3.117782000  -0.376744000  -1.208365000
C    3.321385000  -1.043731000   0.000000000
C    3.117782000  -0.376744000   1.208365000
C    2.713009000   0.957854000   1.208292000
H    2.202404000   2.661136000   0.000000000
H    2.551176000   1.473691000  -2.144590000
H    3.270300000  -0.895141000  -2.144838000
H    3.636814000  -2.078152000   0.000000000
H    3.270300000  -0.895141000   2.144838000
H    2.551176000   1.473691000   2.144590000
--
0 1
H    0.806524000  -0.435887000   0.000000000
N   -0.144241000  -0.768693000   0.000000000
C   -0.516112000  -2.089322000   0.000000000
C   -1.889876000  -2.181449000   0.000000000
C   -2.393232000  -0.847083000   0.000000000
C   -1.264065000   0.019589000   0.000000000
C   -1.389600000   1.411767000   0.000000000
C   -2.672650000   1.936645000   0.000000000
C   -3.805451000   1.097479000   0.000000000
C   -3.679817000  -0.281721000   0.000000000
H    0.231002000  -2.865317000   0.000000000
H   -2.458576000  -3.095605000   0.000000000
H   -0.518873000   2.053952000   0.000000000
H   -2.807757000   3.009786000   0.000000000
H   -4.790599000   1.543937000   0.000000000
H   -4.558019000  -0.914292000   0.000000000
units angstrom
}
""")

S22by5_21_1p2 = input.process_input("""
molecule dimer {
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    2.677063232   0.000000000   0.000000000
N    3.684380636   0.005056388   0.000000000
C    4.466034820   1.132604937   0.000000000
C    5.794919547   0.772354616   0.000000000
C    5.836313033  -0.653193161   0.000000000
C    4.486468397  -1.104920876   0.000000000
C    4.151808728  -2.462094972   0.000000000
C    5.194001034  -3.376178892   0.000000000
C    6.538324473  -2.951641292   0.000000000
C    6.868664352  -1.606705567   0.000000000
H    4.012340595   2.109594763   0.000000000
H    6.630309349   1.451489921   0.000000000
H    3.119419786  -2.785730808   0.000000000
H    4.972283588  -4.434822780   0.000000000
H    7.324389392  -3.694570139   0.000000000
H    7.904990129  -1.294593877   0.000000000
units angstrom
}
""")

S22by5_21_1p5 = input.process_input("""
molecule dimer {
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    3.346329040   0.000000000   0.000000000
N    4.353646444   0.005056388   0.000000000
C    5.135300628   1.132604937   0.000000000
C    6.464185355   0.772354616   0.000000000
C    6.505578841  -0.653193161   0.000000000
C    5.155734205  -1.104920876   0.000000000
C    4.821074536  -2.462094972   0.000000000
C    5.863266842  -3.376178892   0.000000000
C    7.207590281  -2.951641292   0.000000000
C    7.537930160  -1.606705567   0.000000000
H    4.681606403   2.109594763   0.000000000
H    7.299575157   1.451489921   0.000000000
H    3.788685594  -2.785730808   0.000000000
H    5.641549396  -4.434822780   0.000000000
H    7.993655200  -3.694570139   0.000000000
H    8.574255937  -1.294593877   0.000000000
units angstrom
}
""")

S22by5_21_2p0 = input.process_input("""
molecule dimer {
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    4.461772054   0.000000000   0.000000000
N    5.469089458   0.005056388   0.000000000
C    6.250743642   1.132604937   0.000000000
C    7.579628369   0.772354616   0.000000000
C    7.621021855  -0.653193161   0.000000000
C    6.271177219  -1.104920876   0.000000000
C    5.936517550  -2.462094972   0.000000000
C    6.978709856  -3.376178892   0.000000000
C    8.323033295  -2.951641292   0.000000000
C    8.653373174  -1.606705567   0.000000000
H    5.797049417   2.109594763   0.000000000
H    8.415018171   1.451489921   0.000000000
H    4.904128608  -2.785730808   0.000000000
H    6.756992410  -4.434822780   0.000000000
H    9.109098214  -3.694570139   0.000000000
H    9.689698951  -1.294593877   0.000000000
units angstrom
}
""")

S22by5_22_0p9 = input.process_input("""
molecule dimer {
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    1.743489077   0.000000000   0.000000000
C    2.341981491  -1.142898789  -0.483732445
H    2.342838533   0.417604441   0.628041164
C    1.645485086  -1.867622674  -1.447211527
C    2.204739700  -3.035912794  -1.954567993
C    3.449296078  -3.479350313  -1.509647408
C    4.136609561  -2.744696418  -0.547410307
C    3.584309534  -1.574952605  -0.029436748
H    0.681454799  -1.513028491  -1.784467064
H    1.661729182  -3.600082357  -2.699896207
H    3.877956013  -4.387511286  -1.908204233
H    5.102623102  -3.077497147  -0.194005162
H    4.116289930  -1.004251641   0.722333197
units angstrom
}
""")

S22by5_22_1p0 = input.process_input("""
molecule dimer {
0 1
C   -2.007106000   0.763846000  -0.108351000
O   -1.388504000   1.929852000  -0.443121000
H   -0.523812000   1.964652000  -0.006461000
C   -1.463081000  -0.151912000   0.794993000
C   -2.147579000  -1.329509000   1.088368000
C   -3.374321000  -1.603143000   0.489586000
C   -3.914373000  -0.683855000  -0.409103000
C   -3.237050000   0.492961000  -0.709613000
H   -0.510651000   0.056657000   1.264256000
H   -1.715113000  -2.032145000   1.787842000
H   -3.902466000  -2.517387000   0.719795000
H   -4.867073000  -0.882294000  -0.881132000
H   -3.643166000   1.213434000  -1.405759000
--
0 1
O    1.353117000   1.938272000   0.472313000
C    2.036975000   0.786504000   0.149549000
H    1.784285000   2.348749000   1.229711000
C    1.590403000   0.069686000  -0.957415000
C    2.241737000  -1.106977000  -1.312811000
C    3.331567000  -1.566560000  -0.574864000
C    3.769684000  -0.839690000   0.528644000
C    3.122484000   0.338350000   0.896049000
H    0.744551000   0.436798000  -1.521858000
H    1.892146000  -1.664973000  -2.170184000
H    3.833023000  -2.481154000  -0.856667000
H    4.613763000  -1.185010000   1.109263000
H    3.459885000   0.903038000   1.756949000
units angstrom
}
""")

S22by5_22_1p2 = input.process_input("""
molecule dimer {
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    2.324652103   0.000000000   0.000000000
C    2.923144517  -1.142898789  -0.483732445
H    2.924001559   0.417604441   0.628041164
C    2.226648112  -1.867622674  -1.447211527
C    2.785902726  -3.035912794  -1.954567993
C    4.030459104  -3.479350313  -1.509647408
C    4.717772587  -2.744696418  -0.547410307
C    4.165472560  -1.574952605  -0.029436748
H    1.262617825  -1.513028491  -1.784467064
H    2.242892208  -3.600082357  -2.699896207
H    4.459119039  -4.387511286  -1.908204233
H    5.683786128  -3.077497147  -0.194005162
H    4.697452956  -1.004251641   0.722333197
units angstrom
}
""")

S22by5_22_1p5 = input.process_input("""
molecule dimer {
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    2.905815129   0.000000000   0.000000000
C    3.504307543  -1.142898789  -0.483732445
H    3.505164585   0.417604441   0.628041164
C    2.807811138  -1.867622674  -1.447211527
C    3.367065752  -3.035912794  -1.954567993
C    4.611622130  -3.479350313  -1.509647408
C    5.298935613  -2.744696418  -0.547410307
C    4.746635586  -1.574952605  -0.029436748
H    1.843780851  -1.513028491  -1.784467064
H    2.824055234  -3.600082357  -2.699896207
H    5.040282065  -4.387511286  -1.908204233
H    6.264949154  -3.077497147  -0.194005162
H    5.278615982  -1.004251641   0.722333197
units angstrom
}
""")

S22by5_22_2p0 = input.process_input("""
molecule dimer {
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    3.874420172   0.000000000   0.000000000
C    4.472912586  -1.142898789  -0.483732445
H    4.473769628   0.417604441   0.628041164
C    3.776416181  -1.867622674  -1.447211527
C    4.335670795  -3.035912794  -1.954567993
C    5.580227173  -3.479350313  -1.509647408
C    6.267540656  -2.744696418  -0.547410307
C    5.715240629  -1.574952605  -0.029436748
H    2.812385894  -1.513028491  -1.784467064
H    3.792660277  -3.600082357  -2.699896207
H    6.008887108  -4.387511286  -1.908204233
H    7.233554197  -3.077497147  -0.194005162
H    6.247221025  -1.004251641   0.722333197
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

