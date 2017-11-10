"""
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.


- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< ALK6 Database Module >>>
dbse = 'ALK6'

# <<< Database Members >>>
HRXN = ['1', '2', '3', '4', '5', '6', ]
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'catLi'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz_catLi') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [1, 1, -1]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'catNa'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz_catNa') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [1, 1, -1]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'catK'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz'),
                                                               '%s-%s-reagent'      % (dbse, 'Bz_catK') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [1, 1, -1]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'Li2'),
                                                               '%s-%s-reagent'      % (dbse, 'Li8') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [4, -1]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'Na2'),
                                                               '%s-%s-reagent'      % (dbse, 'Na8') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [4, -1]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'K2'),
                                                               '%s-%s-reagent'      % (dbse, 'K8') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [4, -1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =    38.4
BIND['%s-%s'            % (dbse, '2'                     )] =    25.0
BIND['%s-%s'            % (dbse, '3'                     )] =    19.2
BIND['%s-%s'            % (dbse, '4'                     )] =    83.2
BIND['%s-%s'            % (dbse, '5'                     )] =    54.6
BIND['%s-%s'            % (dbse, '6'                     )] =    47.1
#Taken from Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H J. Chem. Phys. 2010, 132, 154104.
#(est. CCSD(T)/CBS reference values); all values are in kcal/mol.
# rxn directions left as-is

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Dissociation: Bz ... Li+ --> Bz + Li+"""
TAGL['%s-%s'            % (dbse, '2'                     )] = """Dissociation: Bz ... Na+ --> Bz + Na+"""
TAGL['%s-%s'            % (dbse, '3'                     )] = """Dissociation: Bz ... K+ --> Bz + K+"""
TAGL['%s-%s'            % (dbse, '4'                     )] = """Fragmentation: Li_8 --> 4 Li_2"""
TAGL['%s-%s'            % (dbse, '5'                     )] = """Fragmentation: Na_8 --> 4 Na_2"""
TAGL['%s-%s'            % (dbse, '6'                     )] = """Fragmentation: K_8 --> 4 K_2"""
TAGL['%s-%s-reagent'    % (dbse, 'Bz'                    )] = """benzene"""
TAGL['%s-%s-reagent'    % (dbse, 'Bz_catK'               )] = """Benzene-Potassium Cation Complex"""
TAGL['%s-%s-reagent'    % (dbse, 'Bz_catLi'              )] = """Benzene-Lithium Cation Complex"""
TAGL['%s-%s-reagent'    % (dbse, 'Bz_catNa'              )] = """Benzene-Sodium Cation Complex"""
TAGL['%s-%s-reagent'    % (dbse, 'K2'                    )] = """diatomic potassium"""
TAGL['%s-%s-reagent'    % (dbse, 'K8'                    )] = """potassium cluster"""
TAGL['%s-%s-reagent'    % (dbse, 'Li2'                   )] = """diatomic lithium"""
TAGL['%s-%s-reagent'    % (dbse, 'Li8'                   )] = """lithium cluster"""
TAGL['%s-%s-reagent'    % (dbse, 'Na2'                   )] = """diatomic sodium"""
TAGL['%s-%s-reagent'    % (dbse, 'Na8'                   )] = """sodium cluster"""
TAGL['%s-%s-reagent'    % (dbse, 'catK'                  )] = """potassium cation"""
TAGL['%s-%s-reagent'    % (dbse, 'catLi'                 )] = """lithium cation"""
TAGL['%s-%s-reagent'    % (dbse, 'catNa'                 )] = """sodium cation"""

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'Bz', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.283302421195     1.318265267518     0.000000000000
    C                2.283302421195    -1.318265267518     0.000000000000
    C                0.000000000000    -2.636530535036     0.000000000000
    C               -2.283302421195    -1.318265267518     0.000000000000
    C               -2.283302421195     1.318265267518     0.000000000000
    C               -0.000000000000     2.636530535036     0.000000000000
    H                4.053885505708     2.340511887984     0.000000000000
    H                4.053885505708    -2.340511887984     0.000000000000
    H                0.000000000000    -4.681023775969     0.000000000000
    H               -4.053885505708    -2.340511887984     0.000000000000
    H               -4.053885505708     2.340511887984     0.000000000000
    H               -0.000000000000     4.681023775969     0.000000000000
""")

GEOS['%s-%s-%s' % (dbse, 'Bz_catK', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    C                2.289566419371     1.321881788551    -0.383514001081
    C                2.289566419371    -1.321881788551    -0.383514001081
    C                0.000000000000    -2.643763577102    -0.383514001081
    C               -2.289566419371    -1.321881788551    -0.383514001081
    C               -2.289566419371     1.321881788551    -0.383514001081
    C               -0.000000000000     2.643763577102    -0.383514001081
    H                4.058975575803     2.343450641324    -0.460366236363
    H                4.058975575803    -2.343450641324    -0.460366236363
    H                0.000000000000    -4.686901282648    -0.460366236363
    H               -4.058975575803    -2.343450641324    -0.460366236363
    H               -4.058975575803     2.343450641324    -0.460366236363
    H               -0.000000000000     4.686901282648    -0.460366236363
    K                0.000000000000     0.000000000000     5.063281424668
""")

GEOS['%s-%s-%s' % (dbse, 'Bz_catLi', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    C                2.295473007685     1.325291958905    -0.262038151151
    C                2.295473007685    -1.325291958905    -0.262038151151
    C                0.000000000000    -2.650583917809    -0.262038151151
    C               -2.295473007685    -1.325291958905    -0.262038151151
    C               -2.295473007685     1.325291958905    -0.262038151151
    C               -0.000000000000     2.650583917809    -0.262038151151
    H                4.064863227082     2.346849878375    -0.281984971346
    H                4.064863227082    -2.346849878375    -0.281984971346
    H                0.000000000000    -4.693699756749    -0.281984971346
    H               -4.064863227082    -2.346849878375    -0.281984971346
    H               -4.064863227082     2.346849878375    -0.281984971346
    H               -0.000000000000     4.693699756749    -0.281984971346
    LI               0.000000000000     0.000000000000     3.264138734977
""")

GEOS['%s-%s-%s' % (dbse, 'Bz_catNa', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    C                2.292301278090     1.323460759969    -0.321474274530
    C                2.292301278090    -1.323460759969    -0.321474274530
    C                0.000000000000    -2.646921519938    -0.321474274530
    C               -2.292301278090    -1.323460759969    -0.321474274530
    C               -2.292301278090     1.323460759969    -0.321474274530
    C               -0.000000000000     2.646921519938    -0.321474274530
    H                4.061577204082     2.344952692111    -0.388773066564
    H                4.061577204082    -2.344952692111    -0.388773066564
    H                0.000000000000    -4.689905384222    -0.388773066564
    H               -4.061577204082    -2.344952692111    -0.388773066564
    H               -4.061577204082     2.344952692111    -0.388773066564
    H               -0.000000000000     4.689905384222    -0.388773066564
    NA               0.000000000000     0.000000000000     4.261484046567
""")

GEOS['%s-%s-%s' % (dbse, 'K2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    K                0.000000000000     0.000000000000     3.803514035570
    K                0.000000000000     0.000000000000    -3.803514035570
""")

GEOS['%s-%s-%s' % (dbse, 'K8', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    K                4.958434244000    -4.958434244000     4.958434244000
    K               -4.958434244000     4.958434244000     4.958434244000
    K                4.958434244000     4.958434244000    -4.958434244000
    K               -4.958434244000    -4.958434244000    -4.958434244000
    K                3.056578027000     3.056578027000     3.056578027000
    K               -3.056578027000     3.056578027000    -3.056578027000
    K               -3.056578027000    -3.056578027000     3.056578027000
    K                3.056578027000    -3.056578027000    -3.056578027000
""")

GEOS['%s-%s-%s' % (dbse, 'Li2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    LI               0.000000000000     0.000000000000     2.548418199570
    LI               0.000000000000     0.000000000000    -2.548418199570
""")

GEOS['%s-%s-%s' % (dbse, 'Li8', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    LI               3.240570144106    -3.240570144106     3.240570144106
    LI              -3.240570144106     3.240570144106     3.240570144106
    LI               3.240570144106     3.240570144106    -3.240570144106
    LI              -3.240570144106    -3.240570144106    -3.240570144106
    LI               1.879041556182     1.879041556182     1.879041556182
    LI              -1.879041556182     1.879041556182    -1.879041556182
    LI              -1.879041556182    -1.879041556182     1.879041556182
    LI               1.879041556182    -1.879041556182    -1.879041556182
""")

GEOS['%s-%s-%s' % (dbse, 'Na2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    NA               0.000000000000     0.000000000000     2.849722352481
    NA               0.000000000000     0.000000000000    -2.849722352481
""")

GEOS['%s-%s-%s' % (dbse, 'Na8', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    NA               3.746698430393    -3.746698430393     3.746698430393
    NA              -3.746698430393     3.746698430393     3.746698430393
    NA               3.746698430393     3.746698430393    -3.746698430393
    NA              -3.746698430393    -3.746698430393    -3.746698430393
    NA               2.324757827207     2.324757827207     2.324757827207
    NA              -2.324757827207     2.324757827207    -2.324757827207
    NA              -2.324757827207    -2.324757827207     2.324757827207
    NA               2.324757827207    -2.324757827207    -2.324757827207
""")

GEOS['%s-%s-%s' % (dbse, 'catK', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    K                0.000000000000     0.000000000000     0.000000000000
""")

GEOS['%s-%s-%s' % (dbse, 'catLi', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    LI               0.000000000000     0.000000000000     0.000000000000
""")

GEOS['%s-%s-%s' % (dbse, 'catNa', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    1 1
    NA               0.000000000000     0.000000000000     0.000000000000
""")

