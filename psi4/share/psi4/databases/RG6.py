"""
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.
| Taken from Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H J. Chem. Phys. 2010, 132, 154104.
(experimental reference values); all values are in kcal/mol.




- **cp**  ``'off'`` <erase this comment and after unless on is a valid option> || ``'on'``

- **rlxd** ``'off'`` <erase this comment and after unless on is valid option> || ``'on'``


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

# <<< RG6 Database Module >>>
dbse = 'RG6'

# <<< Database Members >>>
HRXN = ['NeNe', 'ArAr', 'KrKr', 'XeXe', 'RnXe', 'RnRn']
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
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

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, 'NeNe'                  )] =   -0.08  # 10.1039/B509242F  Table 2, exp
BIND['%s-%s'            % (dbse, 'ArAr'                  )] =   -0.28  # "
BIND['%s-%s'            % (dbse, 'KrKr'                  )] =   -0.40  # "
BIND['%s-%s'            % (dbse, 'XeXe'                  )] =   -0.56  # 10.1002/(SICI)1097-461X(1998)66:2<131::AID-QUA4>3.0.CO;2-W  Table VII, exp-corrected
BIND['%s-%s'            % (dbse, 'RnXe'                  )] =   -0.65  # "
BIND['%s-%s'            % (dbse, 'RnRn'                  )] =   -0.79  # "
# Note that Ref has been inverted to D - M - M

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'ArAr'                  )] = """Argon Dimer"""
TAGL['%s-%s-dimer'      % (dbse, 'ArAr'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ArAr'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ArAr'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ArAr'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ArAr'                  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'KrKr'                  )] = """Krypton Dimer"""
TAGL['%s-%s-dimer'      % (dbse, 'KrKr'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'KrKr'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'KrKr'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'KrKr'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'KrKr'                  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'NeNe'                  )] = """Neon Dimer"""
TAGL['%s-%s-dimer'      % (dbse, 'NeNe'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'NeNe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'NeNe'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'NeNe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'NeNe'                  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'RnRn'                  )] = """Radon Dimer"""
TAGL['%s-%s-dimer'      % (dbse, 'RnRn'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'RnRn'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'RnRn'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'RnRn'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'RnRn'                  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'RnXe'                  )] = """Radon-Xenon"""
TAGL['%s-%s-dimer'      % (dbse, 'RnXe'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'RnXe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'RnXe'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'RnXe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'RnXe'                  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'XeXe'                  )] = """Xenon Dimer"""
TAGL['%s-%s-dimer'      % (dbse, 'XeXe'                  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'XeXe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'XeXe'                  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'XeXe'                  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'XeXe'                  )] = """Monomer B from  """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'ArAr', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AR               0.000000000000     0.000000000000     3.552700000000
    --
    0 1
    AR               0.000000000000     0.000000000000    -3.552700000000

""")

GEOS['%s-%s-%s' % (dbse, 'KrKr', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    KR               0.000000000000     0.000000000000     3.788900000000
    --
    0 1
    KR               0.000000000000     0.000000000000    -3.788900000000

""")

GEOS['%s-%s-%s' % (dbse, 'NeNe', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    NE               0.000000000000     0.000000000000     2.919600000000
    --
    0 1
    NE               0.000000000000     0.000000000000    -2.919600000000

""")

GEOS['%s-%s-%s' % (dbse, 'RnRn', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    RN               0.000000000000     0.000000000000    -4.231100000000
    --
    0 1
    RN               0.000000000000     0.000000000000     4.231100000000

""")

GEOS['%s-%s-%s' % (dbse, 'RnXe', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    RN               0.000000000000     0.000000000000    -4.183935000000
    --
    0 1
    XE               0.000000000000     0.000000000000     4.183935000000

""")

GEOS['%s-%s-%s' % (dbse, 'XeXe', 'dimer')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    XE               0.000000000000     0.000000000000    -4.121500000000
    --
    0 1
    XE               0.000000000000     0.000000000000     4.121500000000

""")

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)
