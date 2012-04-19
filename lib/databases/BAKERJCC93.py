"""
| Database of molecules that are challenging to optimize.
| Geometries from Baker J. Comput. Chem. 14 1085 (1993), as reported
  in Bakken and Helgaker, J. Chem. Phys. 117, 9160 (2002), with a few
  further corrections.
| No reference energies defined.

- **cp**  ``'off'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``

"""
import re
import input

# <<< BAKERJCC93 Database Module >>>
dbse = 'BAKERJCC93'
isOS = 'true'

# <<< Database Members >>>
HRXN = ['1_3_5_trifluorobenzene', '1_3_5_trisilacyclohexane', '1_3_difluorobenzene', '1_5_difluoronaphthalene', '2_hydroxybicyclopentane', 'ACANIL01', 'acetone', 'acetylene', 'ACHTAR10', 'allene', 'ammonia', 'benzaldehyde', 'benzene', 'benzidine', 'caffeine', 'difuropyrazine', 'dimethylpentane', 'disilyl_ether', 'ethane', 'ethanol', 'furan', 'histidine', 'hydroxysulphane', 'menthone', 'mesityl_oxide', 'methylamine', 'naphthalene', 'neopentane', 'pterin', 'water', ]
HRXN_SM = ['1_3_5_trisilacyclohexane', '2_hydroxybicyclopentane', 'acetone', 'acetylene', 'allene', 'ammonia', 'benzene', 'disilyl_ether', 'ethane', 'ethanol', 'furan', 'hydroxysulphane', 'methylamine', 'neopentane', 'water']
HRXN_LG = ['1_3_difluorobenzene', '1_3_5_trifluorobenzene', '1_5_difluoronaphthalene', 'ACANIL01', 'ACHTAR10', 'benzaldehyde', 'benzidine', 'caffeine', 'difuropyrazine', 'dimethylpentane', 'histidine', 'menthone', 'mesityl_oxide', 'naphthalene', 'pterin']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1_3_5_trifluorobenzene' )] = ['%s-%s-reagent'      % (dbse, '1_3_5_trifluorobenzene')]
RXNM['%s-%s'            % (dbse, '1_3_5_trifluorobenzene' )] = dict(zip(ACTV['%s-%s' % (dbse, '1_3_5_trifluorobenzene')], [+1]))

ACTV['%s-%s'            % (dbse, '1_3_5_trisilacyclohexane' )] = ['%s-%s-reagent'      % (dbse, '1_3_5_trisilacyclohexane')]
RXNM['%s-%s'            % (dbse, '1_3_5_trisilacyclohexane' )] = dict(zip(ACTV['%s-%s' % (dbse, '1_3_5_trisilacyclohexane')], [+1]))

ACTV['%s-%s'            % (dbse, '1_3_difluorobenzene'   )] = ['%s-%s-reagent'      % (dbse, '1_3_difluorobenzene')]
RXNM['%s-%s'            % (dbse, '1_3_difluorobenzene'   )] = dict(zip(ACTV['%s-%s' % (dbse, '1_3_difluorobenzene')], [+1]))

ACTV['%s-%s'            % (dbse, '1_5_difluoronaphthalene' )] = ['%s-%s-reagent'      % (dbse, '1_5_difluoronaphthalene')]
RXNM['%s-%s'            % (dbse, '1_5_difluoronaphthalene' )] = dict(zip(ACTV['%s-%s' % (dbse, '1_5_difluoronaphthalene')], [+1]))

ACTV['%s-%s'            % (dbse, '2_hydroxybicyclopentane' )] = ['%s-%s-reagent'      % (dbse, '2_hydroxybicyclopentane')]
RXNM['%s-%s'            % (dbse, '2_hydroxybicyclopentane' )] = dict(zip(ACTV['%s-%s' % (dbse, '2_hydroxybicyclopentane')], [+1]))

ACTV['%s-%s'            % (dbse, 'ACANIL01'              )] = ['%s-%s-reagent'      % (dbse, 'ACANIL01')]
RXNM['%s-%s'            % (dbse, 'ACANIL01'              )] = dict(zip(ACTV['%s-%s' % (dbse, 'ACANIL01')], [+1]))

ACTV['%s-%s'            % (dbse, 'acetone'               )] = ['%s-%s-reagent'      % (dbse, 'acetone')]
RXNM['%s-%s'            % (dbse, 'acetone'               )] = dict(zip(ACTV['%s-%s' % (dbse, 'acetone')], [+1]))

ACTV['%s-%s'            % (dbse, 'acetylene'             )] = ['%s-%s-reagent'      % (dbse, 'acetylene')]
RXNM['%s-%s'            % (dbse, 'acetylene'             )] = dict(zip(ACTV['%s-%s' % (dbse, 'acetylene')], [+1]))

ACTV['%s-%s'            % (dbse, 'ACHTAR10'              )] = ['%s-%s-reagent'      % (dbse, 'ACHTAR10')]
RXNM['%s-%s'            % (dbse, 'ACHTAR10'              )] = dict(zip(ACTV['%s-%s' % (dbse, 'ACHTAR10')], [+1]))

ACTV['%s-%s'            % (dbse, 'allene'                )] = ['%s-%s-reagent'      % (dbse, 'allene')]
RXNM['%s-%s'            % (dbse, 'allene'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'allene')], [+1]))

ACTV['%s-%s'            % (dbse, 'ammonia'               )] = ['%s-%s-reagent'      % (dbse, 'ammonia')]
RXNM['%s-%s'            % (dbse, 'ammonia'               )] = dict(zip(ACTV['%s-%s' % (dbse, 'ammonia')], [+1]))

ACTV['%s-%s'            % (dbse, 'benzaldehyde'          )] = ['%s-%s-reagent'      % (dbse, 'benzaldehyde')]
RXNM['%s-%s'            % (dbse, 'benzaldehyde'          )] = dict(zip(ACTV['%s-%s' % (dbse, 'benzaldehyde')], [+1]))

ACTV['%s-%s'            % (dbse, 'benzene'               )] = ['%s-%s-reagent'      % (dbse, 'benzene')]
RXNM['%s-%s'            % (dbse, 'benzene'               )] = dict(zip(ACTV['%s-%s' % (dbse, 'benzene')], [+1]))

ACTV['%s-%s'            % (dbse, 'benzidine'             )] = ['%s-%s-reagent'      % (dbse, 'benzidine')]
RXNM['%s-%s'            % (dbse, 'benzidine'             )] = dict(zip(ACTV['%s-%s' % (dbse, 'benzidine')], [+1]))

ACTV['%s-%s'            % (dbse, 'caffeine'              )] = ['%s-%s-reagent'      % (dbse, 'caffeine')]
RXNM['%s-%s'            % (dbse, 'caffeine'              )] = dict(zip(ACTV['%s-%s' % (dbse, 'caffeine')], [+1]))

ACTV['%s-%s'            % (dbse, 'difuropyrazine'        )] = ['%s-%s-reagent'      % (dbse, 'difuropyrazine')]
RXNM['%s-%s'            % (dbse, 'difuropyrazine'        )] = dict(zip(ACTV['%s-%s' % (dbse, 'difuropyrazine')], [+1]))

ACTV['%s-%s'            % (dbse, 'dimethylpentane'       )] = ['%s-%s-reagent'      % (dbse, 'dimethylpentane')]
RXNM['%s-%s'            % (dbse, 'dimethylpentane'       )] = dict(zip(ACTV['%s-%s' % (dbse, 'dimethylpentane')], [+1]))

ACTV['%s-%s'            % (dbse, 'disilyl_ether'         )] = ['%s-%s-reagent'      % (dbse, 'disilyl_ether')]
RXNM['%s-%s'            % (dbse, 'disilyl_ether'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'disilyl_ether')], [+1]))

ACTV['%s-%s'            % (dbse, 'ethane'                )] = ['%s-%s-reagent'      % (dbse, 'ethane')]
RXNM['%s-%s'            % (dbse, 'ethane'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'ethane')], [+1]))

ACTV['%s-%s'            % (dbse, 'ethanol'               )] = ['%s-%s-reagent'      % (dbse, 'ethanol')]
RXNM['%s-%s'            % (dbse, 'ethanol'               )] = dict(zip(ACTV['%s-%s' % (dbse, 'ethanol')], [+1]))

ACTV['%s-%s'            % (dbse, 'furan'                 )] = ['%s-%s-reagent'      % (dbse, 'furan')]
RXNM['%s-%s'            % (dbse, 'furan'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'furan')], [+1]))

ACTV['%s-%s'            % (dbse, 'histidine'             )] = ['%s-%s-reagent'      % (dbse, 'histidine')]
RXNM['%s-%s'            % (dbse, 'histidine'             )] = dict(zip(ACTV['%s-%s' % (dbse, 'histidine')], [+1]))

ACTV['%s-%s'            % (dbse, 'hydroxysulphane'       )] = ['%s-%s-reagent'      % (dbse, 'hydroxysulphane')]
RXNM['%s-%s'            % (dbse, 'hydroxysulphane'       )] = dict(zip(ACTV['%s-%s' % (dbse, 'hydroxysulphane')], [+1]))

ACTV['%s-%s'            % (dbse, 'menthone'              )] = ['%s-%s-reagent'      % (dbse, 'menthone')]
RXNM['%s-%s'            % (dbse, 'menthone'              )] = dict(zip(ACTV['%s-%s' % (dbse, 'menthone')], [+1]))

ACTV['%s-%s'            % (dbse, 'mesityl_oxide'         )] = ['%s-%s-reagent'      % (dbse, 'mesityl_oxide')]
RXNM['%s-%s'            % (dbse, 'mesityl_oxide'         )] = dict(zip(ACTV['%s-%s' % (dbse, 'mesityl_oxide')], [+1]))

ACTV['%s-%s'            % (dbse, 'methylamine'           )] = ['%s-%s-reagent'      % (dbse, 'methylamine')]
RXNM['%s-%s'            % (dbse, 'methylamine'           )] = dict(zip(ACTV['%s-%s' % (dbse, 'methylamine')], [+1]))

ACTV['%s-%s'            % (dbse, 'naphthalene'           )] = ['%s-%s-reagent'      % (dbse, 'naphthalene')]
RXNM['%s-%s'            % (dbse, 'naphthalene'           )] = dict(zip(ACTV['%s-%s' % (dbse, 'naphthalene')], [+1]))

ACTV['%s-%s'            % (dbse, 'neopentane'            )] = ['%s-%s-reagent'      % (dbse, 'neopentane')]
RXNM['%s-%s'            % (dbse, 'neopentane'            )] = dict(zip(ACTV['%s-%s' % (dbse, 'neopentane')], [+1]))

ACTV['%s-%s'            % (dbse, 'pterin'                )] = ['%s-%s-reagent'      % (dbse, 'pterin')]
RXNM['%s-%s'            % (dbse, 'pterin'                )] = dict(zip(ACTV['%s-%s' % (dbse, 'pterin')], [+1]))

ACTV['%s-%s'            % (dbse, 'water'                 )] = ['%s-%s-reagent'      % (dbse, 'water')]
RXNM['%s-%s'            % (dbse, 'water'                 )] = dict(zip(ACTV['%s-%s' % (dbse, 'water')], [+1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1_3_5_trifluorobenzene' )] =    0.000
BIND['%s-%s'            % (dbse, '1_3_5_trisilacyclohexane' )] =    0.000
BIND['%s-%s'            % (dbse, '1_3_difluorobenzene'   )] =    0.000
BIND['%s-%s'            % (dbse, '1_5_difluoronaphthalene' )] =    0.000
BIND['%s-%s'            % (dbse, '2_hydroxybicyclopentane' )] =    0.000
BIND['%s-%s'            % (dbse, 'ACANIL01'              )] =    0.000
BIND['%s-%s'            % (dbse, 'acetone'               )] =    0.000
BIND['%s-%s'            % (dbse, 'acetylene'             )] =    0.000
BIND['%s-%s'            % (dbse, 'ACHTAR10'              )] =    0.000
BIND['%s-%s'            % (dbse, 'allene'                )] =    0.000
BIND['%s-%s'            % (dbse, 'ammonia'               )] =    0.000
BIND['%s-%s'            % (dbse, 'benzaldehyde'          )] =    0.000
BIND['%s-%s'            % (dbse, 'benzene'               )] =    0.000
BIND['%s-%s'            % (dbse, 'benzidine'             )] =    0.000
BIND['%s-%s'            % (dbse, 'caffeine'              )] =    0.000
BIND['%s-%s'            % (dbse, 'difuropyrazine'        )] =    0.000
BIND['%s-%s'            % (dbse, 'dimethylpentane'       )] =    0.000
BIND['%s-%s'            % (dbse, 'disilyl_ether'         )] =    0.000
BIND['%s-%s'            % (dbse, 'ethane'                )] =    0.000
BIND['%s-%s'            % (dbse, 'ethanol'               )] =    0.000
BIND['%s-%s'            % (dbse, 'furan'                 )] =    0.000
BIND['%s-%s'            % (dbse, 'histidine'             )] =    0.000
BIND['%s-%s'            % (dbse, 'hydroxysulphane'       )] =    0.000
BIND['%s-%s'            % (dbse, 'menthone'              )] =    0.000
BIND['%s-%s'            % (dbse, 'mesityl_oxide'         )] =    0.000
BIND['%s-%s'            % (dbse, 'methylamine'           )] =    0.000
BIND['%s-%s'            % (dbse, 'naphthalene'           )] =    0.000
BIND['%s-%s'            % (dbse, 'neopentane'            )] =    0.000
BIND['%s-%s'            % (dbse, 'pterin'                )] =    0.000
BIND['%s-%s'            % (dbse, 'water'                 )] =    0.000

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1_3_5_trifluorobenzene' )] = ''
TAGL['%s-%s-reagent'    % (dbse, '1_3_5_trifluorobenzene' )] = ''
TAGL['%s-%s'            % (dbse, '1_3_5_trisilacyclohexane' )] = ''
TAGL['%s-%s-reagent'    % (dbse, '1_3_5_trisilacyclohexane' )] = ''
TAGL['%s-%s'            % (dbse, '1_3_difluorobenzene'   )] = ''
TAGL['%s-%s-reagent'    % (dbse, '1_3_difluorobenzene'   )] = ''
TAGL['%s-%s'            % (dbse, '1_5_difluoronaphthalene' )] = ''
TAGL['%s-%s-reagent'    % (dbse, '1_5_difluoronaphthalene' )] = ''
TAGL['%s-%s'            % (dbse, '2_hydroxybicyclopentane' )] = ''
TAGL['%s-%s-reagent'    % (dbse, '2_hydroxybicyclopentane' )] = ''
TAGL['%s-%s'            % (dbse, 'ACANIL01'              )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'ACANIL01'              )] = ''
TAGL['%s-%s'            % (dbse, 'acetone'               )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'acetone'               )] = ''
TAGL['%s-%s'            % (dbse, 'acetylene'             )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'acetylene'             )] = ''
TAGL['%s-%s'            % (dbse, 'ACHTAR10'              )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'ACHTAR10'              )] = ''
TAGL['%s-%s'            % (dbse, 'allene'                )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'allene'                )] = ''
TAGL['%s-%s'            % (dbse, 'ammonia'               )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'ammonia'               )] = ''
TAGL['%s-%s'            % (dbse, 'benzaldehyde'          )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'benzaldehyde'          )] = ''
TAGL['%s-%s'            % (dbse, 'benzene'               )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'benzene'               )] = ''
TAGL['%s-%s'            % (dbse, 'benzidine'             )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'benzidine'             )] = ''
TAGL['%s-%s'            % (dbse, 'caffeine'              )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'caffeine'              )] = ''
TAGL['%s-%s'            % (dbse, 'difuropyrazine'        )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'difuropyrazine'        )] = ''
TAGL['%s-%s'            % (dbse, 'dimethylpentane'       )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'dimethylpentane'       )] = ''
TAGL['%s-%s'            % (dbse, 'disilyl_ether'         )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'disilyl_ether'         )] = ''
TAGL['%s-%s'            % (dbse, 'ethane'                )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'ethane'                )] = ''
TAGL['%s-%s'            % (dbse, 'ethanol'               )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'ethanol'               )] = ''
TAGL['%s-%s'            % (dbse, 'furan'                 )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'furan'                 )] = ''
TAGL['%s-%s'            % (dbse, 'histidine'             )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'histidine'             )] = ''
TAGL['%s-%s'            % (dbse, 'hydroxysulphane'       )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'hydroxysulphane'       )] = ''
TAGL['%s-%s'            % (dbse, 'menthone'              )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'menthone'              )] = ''
TAGL['%s-%s'            % (dbse, 'mesityl_oxide'         )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'mesityl_oxide'         )] = ''
TAGL['%s-%s'            % (dbse, 'methylamine'           )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'methylamine'           )] = ''
TAGL['%s-%s'            % (dbse, 'naphthalene'           )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'naphthalene'           )] = ''
TAGL['%s-%s'            % (dbse, 'neopentane'            )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'neopentane'            )] = ''
TAGL['%s-%s'            % (dbse, 'pterin'                )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'pterin'                )] = ''
TAGL['%s-%s'            % (dbse, 'water'                 )] = ''
TAGL['%s-%s-reagent'    % (dbse, 'water'                 )] = ''

# <<< Molecule Specifications >>>
BAKERJCC93_1_3_5_trifluorobenzene = input.process_input("""
molecule {
0 1
F        4.45124771     2.56992907     0.00000000
F       -4.45124771     2.56992907     0.00000000
F        0.00000000    -5.13985813     0.00000000
C        2.27501122     1.31347834     0.00000000
C       -2.27501122     1.31347834     0.00000000
C        0.00000000    -2.62695668     0.00000000
C        2.27446593    -1.31316352     0.00000000
C       -2.27446593    -1.31316352     0.00000000
C        0.00000000     2.62632703     0.00000000
H        4.04176646    -2.33351496     0.00000000
H       -4.04176646    -2.33351496     0.00000000
H        0.00000000     4.66702991     0.00000000
units bohr
}
""", 0)

BAKERJCC93_1_3_5_trisilacyclohexane = input.process_input("""
molecule {
0 1
Si       2.87562701     1.66024403     0.50009833
Si      -2.87562701     1.66024403     0.50009833
Si       0.00000000    -3.32048805     0.50009833
C        0.00000000     3.31617083    -0.65645952
C        2.87188818    -1.65808542    -0.65645952
C       -2.87188818    -1.65808542    -0.65645952
H        0.00000000     5.25402682     0.04550787
H        4.55012070    -2.62701341     0.04550787
H       -4.55012070    -2.62701341     0.04550787
H        0.00000000     3.33620321    -2.71676085
H        2.88923673    -1.66810160    -2.71676085
H       -2.88923673    -1.66810160    -2.71676085
H        5.14953250     2.97308398    -0.46837999
H       -5.14953250     2.97308398    -0.46837999
H        2.91112385     1.68073814     3.29599415
H       -2.91112385     1.68073814     3.29599415
H        0.00000000    -3.36147627     3.29599415
H        0.00000000    -5.94616795    -0.46837999
units bohr
}
""", 0)

BAKERJCC93_1_3_difluorobenzene = input.process_input("""
molecule {
0 1
F        4.45098629     2.53075455     0.00000000
F       -4.45098629     2.53075455     0.00000000
C        2.27459315    -1.35284979     0.00000000
C       -2.27459315    -1.35284979     0.00000000
C        2.27465109     1.27385640     0.00000000
C       -2.27465109     1.27385640     0.00000000
C        0.00000000     2.58727941     0.00000000
C        0.00000000    -2.66641919     0.00000000
H        4.04232694    -2.37256182     0.00000000
H       -4.04232694    -2.37256182     0.00000000
H        0.00000000     4.62804882     0.00000000
H        0.00000000    -4.70730774     0.00000000
units bohr
}
""", 0)

BAKERJCC93_1_5_difluoronaphthalene = input.process_input("""
molecule {
0 1
F        5.77442810     0.00000000     0.00000000
F       -5.77442810     0.00000000     0.00000000
C        0.72785457    -4.70254512     0.00000000
C       -0.72785457     4.70254512     0.00000000
C        3.11062174    -3.60249243     0.00000000
C       -3.11062174     3.60249243     0.00000000
C        3.38479931    -0.98799287     0.00000000
C       -3.38479931     0.98799287     0.00000000
C        1.23776851     0.57124055     0.00000000
C       -1.23776851    -0.57124055     0.00000000
C        1.43014268     3.20907701     0.00000000
C       -1.43014268    -3.20907701     0.00000000
H        0.55204008    -6.73646406     0.00000000
H       -0.55204008     6.73646406     0.00000000
H        4.76445952    -4.80069021     0.00000000
H       -4.76445952     4.80069021     0.00000000
H        3.24999844     4.13948522     0.00000000
H       -3.24999844    -4.13948522     0.00000000
units bohr
}
""", 0)

BAKERJCC93_2_hydroxybicyclopentane = input.process_input("""
molecule {
0 1
O        0.00000000     0.00000000     3.97630549
C        0.61275612     1.71787828    -0.25674160
C       -1.25240609     0.75430367     1.72991074
C       -1.89991796    -1.49181590     0.06218311
C        2.64764592     0.00000000    -1.36190849
C       -0.10732099    -0.53140328    -1.99092676
H       -2.85601026     2.05561017     2.08353099
H        0.13348920     3.49094037    -1.26103342
H        3.57368102    -1.34558615    -0.05933199
H        3.80698111     0.79833379    -2.90613711
H       -1.33579202    -3.34159783     0.85888891
H       -3.90122993    -1.54049270    -0.53915310
H       -0.93405780     0.26002983    -3.74579160
H        1.51218168    -0.82620025     3.41020482
units bohr
}
""", 0)

BAKERJCC93_ACANIL01 = input.process_input("""
molecule {
0 1
O        6.74334167     0.00000000     0.00000000
N        2.75125398    -0.91996681     0.00000000
C       -3.75958919    -3.62046813     0.00000000
C       -1.13660145    -3.38720984     0.00000000
C        0.00427371    -1.00318363     0.00000000
C       -1.53985353     1.15387105     0.00000000
C       -4.16293704     0.91831969     0.00000000
C       -5.26811078    -1.46724271     0.00000000
C        4.57389611     0.80511522     0.00000000
C        4.13207020     3.64054919     0.00000000
H       -4.62306754    -5.47176436     0.00000000
H       -0.00377805    -5.08765397     0.00000000
H       -0.76505606     3.03326534     0.00000000
H       -5.34041651     2.58758240     0.00000000
H       -7.30266577    -1.64802574     0.00000000
H        3.58082506    -2.66479209     0.00000000
H        5.95212032     4.66178191     0.00000000
H        3.08214744     4.23491124     1.70076220
H        3.08214744     4.23491124    -1.70076220
units bohr
}
""", 0)

BAKERJCC93_acetone = input.process_input("""
molecule {
0 1
O        0.00000000     3.46695757     0.00000000
C        0.00000000     1.14032594     0.00000000
C        0.00000000    -0.29542841     2.50138172
C        0.00000000    -0.29542841    -2.50138172
H        0.00000000     1.00440652     4.13754069
H        0.00000000     1.00440652    -4.13754069
H        1.69360304    -1.50630994     2.66984804
H       -1.69360304    -1.50630994     2.66984804
H        1.69360304    -1.50630994    -2.66984804
H       -1.69360304    -1.50630994    -2.66984804
units bohr
}
""", 0)

BAKERJCC93_acetylene = input.process_input("""
molecule {
0 1
C        0.00000000     0.00000000     1.13383600
C        0.00000000     0.00000000    -1.13383600
H        0.00000000     0.00000000     3.02356266
H        0.00000000     0.00000000    -3.02356266
units bohr
}
""", 0)

BAKERJCC93_ACHTAR10 = input.process_input("""
molecule {
0 1
O        0.00000000     0.00000000     3.93735249
O        1.79875939     0.00000000    -0.09531034
N       -4.40589519     1.32037243    -3.31810156
C       -2.43021636    -0.18962157    -2.05696026
C       -0.22185404     1.49597798    -1.20775357
C        1.69726730    -0.59259412     2.46067577
C        3.97685548    -2.11479138     3.27934906
H       -3.68043380     2.27933244    -4.84082518
H       -5.10144333     2.68085421    -2.12147722
H       -3.24985392    -1.18842676    -0.41051393
H       -1.74547418    -1.68142667    -3.35347133
H        0.55351430     2.51912058    -2.85842920
H       -0.88071695     2.99188292     0.10524925
H        5.73529679    -1.04410557     2.94759034
H        4.08562680    -3.90736002     2.21955987
H        3.86856770    -2.56921447     5.31306580
units bohr
}
""", 0)

BAKERJCC93_allene = input.process_input("""
molecule {
0 1
C        0.00000000     0.00000000     0.00000000
C        0.00000000     2.49419295     0.00000000
C        0.00000000    -2.49419295     0.00000000
H        1.76772016    -3.51503166     0.00000000
H       -1.76772016    -3.51503166     0.00000000
H        0.00000000     3.51503166     1.76772016
H        0.00000000     3.51503166    -1.76772016
units bohr
}
""", 0)

BAKERJCC93_ammonia = input.process_input("""
molecule {
0 1
N        0.00000000     0.00000000     0.47690250
H        1.55848945     0.89979432    -0.15896750
H       -1.55848945     0.89979432    -0.15896750
H        0.00000000    -1.79958864    -0.15896750
units bohr
}
""", 0)

BAKERJCC93_benzaldehyde = input.process_input("""
molecule {
0 1
O        6.11695944     0.00000000     0.00000000
C       -0.42811838    -2.25953622     0.00000000
C       -2.92869352    -1.43478712     0.00000000
C       -3.46561640     1.14118082     0.00000000
C       -1.50611491     2.89722764     0.00000000
C        0.99614123     2.07851844     0.00000000
C        1.55290207    -0.51034434     0.00000000
C        4.31002394    -1.46969818     0.00000000
H        4.69277313    -3.52434043     0.00000000
H       -0.04838912    -4.26733408     0.00000000
H       -4.45167820    -2.79426839     0.00000000
H       -5.40516702     1.77808408     0.00000000
H       -1.92653663     4.89495992     0.00000000
H        2.49151439     3.47033786     0.00000000
units bohr
}
""", 0)

BAKERJCC93_benzene = input.process_input("""
molecule {
0 1
C        0.00000000     2.63452745     0.00000000
C        0.00000000    -2.63452745     0.00000000
C        2.28156770     1.31726373     0.00000000
C       -2.28156770     1.31726373     0.00000000
C        2.28156770    -1.31726373     0.00000000
C       -2.28156770    -1.31726373     0.00000000
H        0.00000000     4.67589156     0.00000000
H        0.00000000    -4.67589156     0.00000000
H        4.04944088     2.33794578     0.00000000
H       -4.04944088     2.33794578     0.00000000
H        4.04944088    -2.33794578     0.00000000
H       -4.04944088    -2.33794578     0.00000000
units bohr
}
""", 0)

BAKERJCC93_benzidine = input.process_input("""
molecule {
0 1
N        0.00000000     0.00000000     9.17973038
N        0.00000000     0.00000000    -9.17973038
C       -2.20388942     0.56488223     5.36955702
C        2.20388942    -0.56488223     5.36955702
C       -2.20388942    -0.56488223    -5.36955702
C        2.20388942     0.56488223    -5.36955702
C       -2.20706622     0.56349235     2.73912945
C        2.20706622    -0.56349235     2.73912945
C       -2.20706622    -0.56349235    -2.73912945
C        2.20706622     0.56349235    -2.73912945
C        0.00000000     0.00000000     1.32948630
C        0.00000000     0.00000000    -1.32948630
C        0.00000000     0.00000000     6.67931977
C        0.00000000     0.00000000    -6.67931977
H       -3.93022673     1.02227253     6.36283467
H        3.93022673    -1.02227253     6.36283467
H       -3.93022673    -1.02227253    -6.36283467
H        3.93022673     1.02227253    -6.36283467
H       -3.95573979     1.07384957     1.81596643
H        3.95573979    -1.07384957     1.81596643
H       -3.95573979    -1.07384957    -1.81596643
H        3.95573979     1.07384957    -1.81596643
H        1.67837252    -0.43031314    10.04483176
H       -1.67837252     0.43031314    10.04483176
H        1.67837252     0.43031314   -10.04483176
H       -1.67837252    -0.43031314   -10.04483176
units bohr
}
""", 0)

BAKERJCC93_caffeine = input.process_input("""
molecule {
0 1
O       -1.35796495    -4.55968346     0.00000000
O        6.00359465     0.00000000     0.00000000
N       -4.34699530     0.40790868     0.00000000
N       -2.02147868     4.34704366     0.00000000
N        2.40166495     2.29891253     0.00000000
N        2.38963107    -2.32861610     0.00000000
C       -1.73514100     0.00806819     0.00000000
C       -0.39656652     2.28913440     0.00000000
C       -4.28628286     3.03178701     0.00000000
C       -0.23380597    -2.51268919     0.00000000
C        3.66630626    -0.01011647     0.00000000
C        3.91233427    -4.71649369     0.00000000
C        3.86899427     4.70507045     0.00000000
C       -6.50871497    -1.38799138     0.00000000
H       -6.04146873     4.08346573     0.00000000
H        5.15026261    -4.82916673     1.68145730
H        5.15026261    -4.82916673    -1.68145730
H        2.75182289    -6.45823288     0.00000000
H        5.10374160     4.83609896     1.68332379
H        5.10374160     4.83609896    -1.68332379
H        2.65878836     6.40983420     0.00000000
H       -8.34003564    -0.38023773     0.00000000
H       -6.44634525    -2.62051420     1.68782779
H       -6.44634525    -2.62051420    -1.68782779
units bohr
}
set { guess gwh }
""", 0)

BAKERJCC93_difuropyrazine = input.process_input("""
molecule {
0 1
O        5.24048162     0.00000000     0.00000000
O       -5.24048162     0.00000000     0.00000000
N        1.15705376    -2.55150608     0.00000000
N       -1.15705376     2.55150608     0.00000000
C        1.53596834     2.15160317     0.00000000
C       -1.53596834    -2.15160317     0.00000000
C        2.65703648    -0.27471770     0.00000000
C       -2.65703648     0.27471770     0.00000000
C        5.62186670     2.62201493     0.00000000
C       -5.62186670    -2.62201493     0.00000000
C        3.54881353     4.07756019     0.00000000
C       -3.54881353    -4.07756019     0.00000000
H        7.52697788     3.41239127     0.00000000
H       -7.52697788    -3.41239127     0.00000000
H        3.34537092     6.12657010     0.00000000
H       -3.34537092    -6.12657010     0.00000000
units bohr
}
""", 0)

BAKERJCC93_dimethylpentane = input.process_input("""
molecule {
0 1
C       -1.90302142     1.79989214    -3.12819161
C        0.68098191     1.17008149    -1.92744962
C        0.57347759    -0.44273007     0.55332222
C       -0.57536860     1.07092655     2.79511667
C        0.00000000     0.00000000     5.44078830
C        2.40130119     0.00000000    -3.96848713
C       -0.75740445    -3.03396450     0.24401702
H       -3.17069973     2.76026122    -1.77509550
H       -2.89812692     0.08656618    -3.79174516
H       -1.70835535     3.07391957    -4.77327789
H        1.57127154     2.99847536    -1.42186024
H        2.56484657    -0.84063638     1.06972922
H       -2.64484782     1.25561061     2.54633021
H        0.13997457     3.03676459     2.74695911
H       -0.82887852    -1.89737304     5.71865173
H       -0.77969219     1.22647944     6.94230752
H        2.05566062    -0.15968979     5.78408339
H        1.63584569    -1.79524917    -4.71467976
H        2.65159319     1.28124788    -5.60071709
H        4.31235662    -0.40611038    -3.22785050
H       -0.60282164    -4.20392129     1.96624261
H       -2.79329149    -2.82201336    -0.17320237
H        0.07519864    -4.15853702    -1.30499111
units bohr
}
""", 0)

BAKERJCC93_disilyl_ether = input.process_input("""
molecule {
0 1
Si       0.00000000    -0.06571048     3.03636189
Si       0.00000000    -0.06571048    -3.03636189
O        0.00000000    -0.88817346     0.00000000
H        0.00000000    -2.19565412     4.54756839
H        2.12290049     1.35272566     3.58475023
H       -2.12290049     1.35272566     3.58475023
H        0.00000000    -2.19565412    -4.54756839
H        2.12290049     1.35272566    -3.58475023
H       -2.12290049     1.35272566    -3.58475023
units bohr
}
""", 0)

BAKERJCC93_ethane = input.process_input("""
molecule {
0 1
C        0.00000000     0.00000000     1.45478763
C        0.00000000     0.00000000    -1.45478763
H        1.68084455     0.97043609     2.14455455
H        1.68084455    -0.97043609    -2.14455455
H       -1.68084455     0.97043609     2.14455455
H       -1.68084455    -0.97043609    -2.14455455
H        0.00000000    -1.94087219     2.14455455
H        0.00000000     1.94087219    -2.14455455
units bohr
}
""", 0)

BAKERJCC93_ethanol = input.process_input("""
molecule {
0 1
O        2.94951269     0.00000000     0.00000000
C        0.42864361     0.89070972     0.00000000
C       -1.47274991    -1.22612707     0.00000000
H        4.05795769     1.50458064     0.00000000
H        0.07017562     2.06834349     1.69306899
H        0.07017562     2.06834349    -1.69306899
H       -1.36184741    -2.46035674     1.67199009
H       -1.36184741    -2.46035674    -1.67199009
H       -3.38002050    -0.38513679     0.00000000
units bohr
}
""", 0)

BAKERJCC93_furan = input.process_input("""
molecule {
0 1
O        0.00000000    -2.71155703     0.00000000
C        1.30409645     1.35600277     0.00000000
C       -1.30409645     1.35600277     0.00000000
C        2.07680908    -1.14870311     0.00000000
C       -2.07680908    -1.14870311     0.00000000
H        2.51639050     2.99782755     0.00000000
H       -2.51639050     2.99782755     0.00000000
H        3.99399875    -1.84934869     0.00000000
H       -3.99399875    -1.84934869     0.00000000
units bohr
}
""", 0)

BAKERJCC93_histidine = input.process_input("""
molecule {
0 1
O        3.93683911     0.00000000     5.02858545
O        0.00000000     0.00000000     6.75548572
N       -1.62714005    -0.17063169    -6.38981145
N       -1.55525882     2.75691585    -2.98858739
N       -0.06519044    -3.43699076     1.78152280
C        0.00313112    -0.96382673    -4.49205004
C        0.04394056     0.76526971    -2.47860156
C       -2.44306081     2.04255765    -5.31626413
C        1.61712938     0.55012280    -0.07100455
C        0.25915307    -0.68070217     2.22272340
C        1.63605981    -0.20456172     4.77295177
H        1.09061008    -2.68966721    -4.55294959
H       -1.91653466     4.35407321    -1.95395961
H       -3.77917584     3.21725531    -6.32817718
H        2.24417172     2.47631265     0.44479551
H        3.39934986    -0.46804317    -0.47537143
H       -1.63798632     0.18516799     2.40194434
H       -1.09297021    -4.24259361     3.21643716
H       -1.09919371    -3.74001070     0.16810497
H        0.98612615     0.24935257     8.25422581
units bohr
}
""", 0)

BAKERJCC93_hydroxysulphane = input.process_input("""
molecule {
0 1
S        0.00000000     0.00000000     1.64344454
O        1.55643788     0.00000000    -0.78417924
H        0.70878977    -0.98889634    -2.04698233
H       -2.26522765     0.98889634     1.18771703
units bohr
}
""", 0)

BAKERJCC93_menthone = input.process_input("""
molecule {
0 1
O        0.00000000     0.00000000     4.83502957
C       -5.06597212    -1.27592091     0.49885049
C       -3.60348796    -1.49111229    -2.01995066
C       -1.13779972     0.12182250    -2.02402508
C        0.69335828    -0.53324847     0.24699141
C       -0.81879368    -0.76420189     2.79442766
C       -3.41755812    -2.06413311     2.77868746
C        0.08327139    -0.18247958    -4.66184769
C        3.16977849     1.11788425     0.33780916
C        5.23967937     0.00000000     2.05851212
C        2.74820737     3.91648659     1.03692914
H       -5.73534045     0.69223903     0.75660810
H       -6.80139535    -2.44289264     0.43045930
H       -3.15419510    -3.50140339    -2.39715388
H       -4.86109777    -0.90264082    -3.58603040
H       -1.71463208     2.12647811    -1.82542348
H        1.33530286    -2.48925976    -0.10949068
H       -4.41049264    -1.64891601     4.56938165
H       -3.10227312    -4.12767303     2.77142958
H       -1.27515064     0.19340176    -6.20625544
H        0.83979297    -2.10810531    -4.96157719
H        1.65711962     1.15531285    -4.96195049
H        4.01314574     1.10167735    -1.57473542
H        4.69908810     0.02990650     4.07747056
H        7.03475689     1.05859686     1.90111311
H        5.66887645    -1.98486988     1.56898286
H        4.52277834     5.01677786     0.95132487
H        1.98900684     4.13531008     2.97264568
H        1.40402606     4.85096335    -0.25821233
units bohr
}
""", 0)

BAKERJCC93_mesityl_oxide = input.process_input("""
molecule {
0 1
O        4.30492455     0.00000000     0.00000000
C        0.05024721    -3.82629843     0.00000000
C       -1.35087834    -1.32752917     0.00000000
C       -4.20838872    -1.49335398     0.00000000
C       -0.19920658     0.94023239     0.00000000
C        2.60618767     1.58735088     0.00000000
C        3.31537901     4.37315331     0.00000000
H        1.28461121    -4.00174347     1.67716810
H        1.28461121    -4.00174347    -1.67716810
H       -1.21281465    -5.48980179     0.00000000
H       -5.04695592    -0.57944060     1.68053694
H       -5.04695592    -0.57944060    -1.68053694
H       -4.87911033    -3.47264458     0.00000000
H       -1.42668593     2.59614349     0.00000000
H        2.56758840     5.33389639     1.69384111
H        2.56758840     5.33389639    -1.69384111
H        5.38985873     4.60732325     0.00000000
units bohr
}
set { guess gwh }
""", 0)

BAKERJCC93_methylamine = input.process_input("""
molecule {
0 1
N        1.59169309     0.00000000     0.00000000
C       -1.10781247    -0.03073718     0.00000000
H        2.61432616    -1.63020032     0.00000000
H       -1.81666320     1.93163906     0.00000000
H        2.57804913     1.64911594     0.00000000
H       -1.92979635    -0.95990875     1.69695191
H       -1.92979635    -0.95990875    -1.69695191
units bohr
}
""", 0)

BAKERJCC93_naphthalene = input.process_input("""
molecule {
0 1
C        1.31500993     4.56625993     0.00000000
C       -1.31500993     4.56625993     0.00000000
C        1.31500993    -4.56625993     0.00000000
C       -1.31500993    -4.56625993     0.00000000
C        2.65095410     2.30121210     0.00000000
C       -2.65095410     2.30121210     0.00000000
C        2.65095410    -2.30121210     0.00000000
C       -2.65095410    -2.30121210     0.00000000
C        1.35957848     0.00000000     0.00000000
C       -1.35957848     0.00000000     0.00000000
H        2.32713807     6.33915590     0.00000000
H       -2.32713807     6.33915590     0.00000000
H        2.32713807    -6.33915590     0.00000000
H       -2.32713807    -6.33915590     0.00000000
H        4.69449351     2.36375141     0.00000000
H       -4.69449351     2.36375141     0.00000000
H        4.69449351    -2.36375141     0.00000000
H       -4.69449351    -2.36375141     0.00000000
units bohr
}
""", 0)

BAKERJCC93_neopentane = input.process_input("""
molecule {
0 1
C        0.00000000     0.00000000     0.00000000
C        1.68781269    -1.68781269     1.68781269
C       -1.68781269     1.68781269     1.68781269
C       -1.68781269    -1.68781269    -1.68781269
C        1.68781269     1.68781269    -1.68781269
H        2.93275937    -0.55961452     2.93275937
H       -2.93275937     0.55961452     2.93275937
H        2.93275937     0.55961452    -2.93275937
H       -2.93275937    -0.55961452    -2.93275937
H        0.55961452    -2.93275937     2.93275937
H       -0.55961452     2.93275937     2.93275937
H       -0.55961452    -2.93275937    -2.93275937
H        0.55961452     2.93275937    -2.93275937
H        2.93275937    -2.93275937     0.55961452
H       -2.93275937     2.93275937     0.55961452
H       -2.93275937    -2.93275937    -0.55961452
H        2.93275937     2.93275937    -0.55961452
units bohr
}
""", 0)

BAKERJCC93_pterin = input.process_input("""
molecule {
0 1
O        5.40068710     0.00000000     0.00000000
N        1.67450469    -4.01224809     0.00000000
N       -3.29778810    -2.66298586     0.00000000
N        2.41435003     2.98093954     0.00000000
N       -2.07868639     1.96577816     0.00000000
N       -1.05941931     6.33383022     0.00000000
C        1.04506477    -1.58869528     0.00000000
C       -1.57825490    -0.83595727     0.00000000
C       -0.10078959    -5.76401042     0.00000000
C       -2.66568392    -5.07794600     0.00000000
C        3.13177958     0.52435884     0.00000000
C       -0.33744328     3.67065942     0.00000000
H        0.42296241    -7.73663229     0.00000000
H       -4.11548210    -6.51452882     0.00000000
H        3.70204141     4.43045951     0.00000000
H       -2.96601438     6.68621706     0.00000000
H        0.40817199     7.60076128     0.00000000
units bohr
}
set { guess gwh }
""", 0)

BAKERJCC93_water = input.process_input("""
molecule {
0 1
O        0.00000000    -0.69801390     0.00000000
H        1.48150016     0.34900695     0.00000000
H       -1.48150016     0.34900695     0.00000000
units bohr
}
""", 0)

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
    for rgt in ACTV['%s-%s' % (dbse, rxn)]:

        molname = rxnpattern.match(rgt)
        GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))
