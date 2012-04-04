"""
| Database of simple molecules, mostly for testing.
| Geometries from nowhere special, and no reference energies defined.

- **cp**  ``'off'``

- **rlxd** ``'off'``

- **subset** [``'h2o'``, ``'nh3'``, ``'ch4'``]

"""
import re
import input

# <<< BASIC Database Module >>>
# Geometries and Reference energies from nowhere special.
# Feel free to add to this database.
dbse = 'BASIC'

# <<< Database Members >>>
HRXN = ['ch4', 'h2o', 'nh3', ]
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, 'ch4'                   )] = ['%s-%s-reagent'      % (dbse, 'ch4')]
RXNM['%s-%s'            % (dbse, 'ch4'                   )] = dict(zip(ACTV['%s-%s' % (dbse, 'ch4')], [+1]))

ACTV['%s-%s'            % (dbse, 'h2o'                   )] = ['%s-%s-reagent'      % (dbse, 'h2o')]
RXNM['%s-%s'            % (dbse, 'h2o'                   )] = dict(zip(ACTV['%s-%s' % (dbse, 'h2o')], [+1]))

ACTV['%s-%s'            % (dbse, 'nh3'                   )] = ['%s-%s-reagent'      % (dbse, 'nh3')]
RXNM['%s-%s'            % (dbse, 'nh3'                   )] = dict(zip(ACTV['%s-%s' % (dbse, 'nh3')], [+1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, 'ch4'                   )] =    0.000
BIND['%s-%s'            % (dbse, 'h2o'                   )] =    0.000
BIND['%s-%s'            % (dbse, 'nh3'                   )] =    0.000

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'ch4'                   )] = 'methane'
TAGL['%s-%s-reagent'    % (dbse, 'ch4'                   )] = 'methane'
TAGL['%s-%s'            % (dbse, 'h2o'                   )] = 'water'
TAGL['%s-%s-reagent'    % (dbse, 'h2o'                   )] = 'water'
TAGL['%s-%s'            % (dbse, 'nh3'                   )] = 'ammonia'
TAGL['%s-%s-reagent'    % (dbse, 'nh3'                   )] = 'ammonia'

# <<< Molecule Specifications >>>
BASIC_ch4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.00014000     1.85916100
H       -0.88855100     0.51306000     1.49468500
H        0.88855100     0.51306000     1.49468500
H        0.00000000    -1.02633900     1.49486800
H        0.00000000     0.00008900     2.94828400
units angstrom
}
""", 0)

BASIC_h2o = input.process_input("""
molecule dimer {
0 1
O       -1.55100700    -0.11452000     0.00000000
H       -1.93425900     0.76250300     0.00000000
H       -0.59967700     0.04071200     0.00000000
units angstrom
}
""", 0)

BASIC_nh3 = input.process_input("""
molecule dimer {
0 1
N       -1.57871800    -0.04661100     0.00000000
H       -2.15862100     0.13639600    -0.80956500
H       -2.15862100     0.13639600     0.80956500
H       -0.84947100     0.65819300     0.00000000
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
    for rgt in ACTV['%s-%s' % (dbse, rxn)]:

        molname = rxnpattern.match(rgt)
        GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))
