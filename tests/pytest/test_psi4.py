import psi4


def test_psi4_basic():
    """tu1-h2o-energy"""
    #! Sample HF/cc-pVDZ H2O computation

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({'basis': "cc-pVDZ"})
    psi4.energy('scf')

    assert psi4.compare_values(-76.0266327341067125, psi4.get_variable('SCF TOTAL ENERGY'), 6, 'SCF energy')


def test_psi4_cas():
    """casscf-sp"""
    #! CASSCF/6-31G** energy point

    geom = psi4.geometry("""
    O
    H 1 1.00
    H 1 1.00 2 103.1
    """)

    psi4.set_options({
        "basis"           : '6-31G**',
        "reference"       : 'rhf',
        "scf_type"        : 'pk',
        "mcscf_algorithm" : 'ah',
        "qc_module"       : 'detci',
        "nat_orbs"        : True})

    cisd_energy, cisd_wfn = psi4.energy("CISD", return_wfn=True)

    assert psi4.compare_values(-76.2198474477531, cisd_energy, 6, 'CISD Energy')

    psi4.set_options({
        "restricted_docc": [1, 0, 0, 0],
        "active":          [3, 0, 1, 2]})

    casscf_energy = psi4.energy('casscf', ref_wfn=cisd_wfn)

    assert psi4.compare_values(-76.073865006902, casscf_energy, 6, 'CASSCF Energy')


def test_psi4_cc():
    """cc1"""
    #! RHF-CCSD 6-31G** all-electron optimization of the H2O molecule

    h2o = psi4.geometry("""
        O
        H 1 0.97
        H 1 0.97 2 103.0
    """)

    psi4.set_options({"basis": '6-31G**'})

    psi4.optimize('ccsd')

    refnuc   =   9.1654609427539
    refscf   = -76.0229427274435
    refccsd  = -0.20823570806196
    reftotal = -76.2311784355056

    assert psi4.compare_values(refnuc,   h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert psi4.compare_values(refscf,   psi4.get_variable("SCF total energy"), 5, "SCF energy")
    assert psi4.compare_values(refccsd,  psi4.get_variable("CCSD correlation energy"), 4, "CCSD contribution")
    assert psi4.compare_values(reftotal, psi4.get_variable("Current energy"), 7, "Total energy")


def test_psi4_dfmp2():
    """dfmp2-1"""
    #! Density fitted MP2 cc-PVDZ/cc-pVDZ-RI computation of formic acid dimer binding energy
    #! using automatic counterpoise correction.  Monomers are specified using Cartesian coordinates.

    Enuc = 235.946620315069168
    Ecp  = -0.0224119246

    formic_dim = psi4.geometry("""
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
       no_reorient
    """)

    psi4.set_options({
       'basis': 'cc-pvdz',
       'df_basis_scf': 'cc-pvdz-jkfit',
       'df_basis_mp2': 'cc-pvdz-ri',
       # not necessary to specify df_basis* for most basis sets
       'scf_type': 'df',
       'guess': 'sad',
       'd_convergence': 11,
    })

    e_cp = psi4.energy('mp2', bsse_type='cp')

    assert psi4.compare_values(Enuc, formic_dim.nuclear_repulsion_energy(), 7, "Nuclear Repulsion Energy")
    assert psi4.compare_values(Ecp, e_cp, 5, "CP Corrected cc-pVDZ/cc-pVDZ-RI DFMP2")


def test_psi4_sapt():
    """sapt1"""
    #! SAPT0 cc-pVDZ computation of the ethene-ethyne interaction energy, using the cc-pVDZ-JKFIT RI basis for SCF
    #! and cc-pVDZ-RI for SAPT.  Monomer geometries are specified using Cartesian coordinates.

    Eref = [ 85.189064196429101,  -0.00359915058,  0.00362911158,
             -0.00083137117,      -0.00150542374, -0.00230683391 ]

    ethene_ethyne = psi4.geometry("""
         0 1
         C     0.000000    -0.667578    -2.124659
         C     0.000000     0.667578    -2.124659
         H     0.923621    -1.232253    -2.126185
         H    -0.923621    -1.232253    -2.126185
         H    -0.923621     1.232253    -2.126185
         H     0.923621     1.232253    -2.126185
         --
         0 1
         C     0.000000     0.000000     2.900503
         C     0.000000     0.000000     1.693240
         H     0.000000     0.000000     0.627352
         H     0.000000     0.000000     3.963929
         units angstrom
    """)

    # this molecule will crash test if molecule passing broken
    barrier = psi4.geometry("""
     0 1
     He
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "guess": "sad",
        "scf_type": "df",
        "sad_print": 2,
        "d_convergence": 11,
        "puream": True,
        "print": 1})

    psi4.energy('sapt0', molecule=ethene_ethyne)

    Eelst = psi4.get_variable("SAPT ELST ENERGY")
    Eexch = psi4.get_variable("SAPT EXCH ENERGY")
    Eind  = psi4.get_variable("SAPT IND ENERGY")
    Edisp = psi4.get_variable("SAPT DISP ENERGY")
    ET    = psi4.get_variable("SAPT0 TOTAL ENERGY")

    assert psi4.compare_values(Eref[0], ethene_ethyne.nuclear_repulsion_energy(), 9, "Nuclear Repulsion Energy")
    assert psi4.compare_values(Eref[1], Eelst, 6, "SAPT0 Eelst")
    assert psi4.compare_values(Eref[2], Eexch, 6, "SAPT0 Eexch")
    assert psi4.compare_values(Eref[3], Eind, 6, "SAPT0 Eind")
    assert psi4.compare_values(Eref[4], Edisp, 6, "SAPT0 Edisp")
    assert psi4.compare_values(Eref[5], ET, 6, "SAPT0 Etotal")


def test_psi4_scfproperty():
    """scf-property"""
    #! UFH and B3LYP cc-pVQZ properties for the CH2 molecule.

    with open('grid.dat', 'w') as handle:
        handle.write("""\
0.0  0.0  0.0
1.1  1.3  1.4
""")

    ch2 = psi4.geometry("""
        0 3
        c
        h 1 b1
        h 1 b1 2 a1

        b1 = 1.0
        a1 = 125.0
    """)

    # Get a reasonable guess, to save some iterations
    psi4.set_options({
        "scf_type": "pk",
        "basis": "6-31G**",
        "e_convergence": 8,
        "docc": [2, 0, 0, 1],
        "socc": [1, 0, 1, 0],
        "reference": "uhf"})

    ch2.update_geometry()
    assert psi4.compare_values(6.648418918908746, ch2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    props = ['DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES',
             'WIBERG_LOWDIN_INDICES', 'MAYER_INDICES', 'MAYER_INDICES',
             'MO_EXTENTS', 'GRID_FIELD', 'GRID_ESP', 'ESP_AT_NUCLEI',
             'MULTIPOLE(5)', 'NO_OCCUPATIONS']

    psi4.property('scf', properties=props)

    assert psi4.compare_values(psi4.get_variable("CURRENT ENERGY"), -38.91591819679808, 6, "SCF energy")
    assert psi4.compare_values(psi4.get_variable('SCF DIPOLE X'), 0.000000000000, 4, "SCF DIPOLE X")
    assert psi4.compare_values(psi4.get_variable('SCF DIPOLE Y'), 0.000000000000, 4, "SCF DIPOLE Y")
    assert psi4.compare_values(psi4.get_variable('SCF DIPOLE Z'), 0.572697798348, 4, "SCF DIPOLE Z")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE XX'), -7.664066833060, 4, "SCF QUADRUPOLE XX")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE YY'), -6.097755074075, 4, "SCF QUADRUPOLE YY")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE ZZ'), -7.074596012050, 4, "SCF QUADRUPOLE ZZ")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE XY'), 0.000000000000, 4, "SCF QUADRUPOLE XY")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE XZ'), 0.000000000000, 4, "SCF QUADRUPOLE XZ")
    assert psi4.compare_values(psi4.get_variable('SCF QUADRUPOLE YZ'), 0.000000000000, 4, "SCF QUADRUPOLE YZ")

    psi4.property('B3LYP', properties=props)

    assert psi4.compare_values(psi4.get_variable('CURRENT ENERGY'), -39.14134740550916, 6, "B3LYP energy")
    assert psi4.compare_values(psi4.get_variable('B3LYP DIPOLE X'), 0.000000000000, 4, "B3LYP DIPOLE X")
    assert psi4.compare_values(psi4.get_variable('B3LYP DIPOLE Y'), -0.000000000000, 4, "B3LYP DIPOLE Y")
    assert psi4.compare_values(psi4.get_variable('B3LYP DIPOLE Z'), 0.641741521158, 4, "B3LYP DIPOLE Z")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE XX'), -7.616483183211, 4, "B3LYP QUADRUPOLE XX")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE YY'), -6.005896804551, 4, "B3LYP QUADRUPOLE YY")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE ZZ'), -7.021817489904, 4, "B3LYP QUADRUPOLE ZZ")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE XY'), 0.000000000000, 4, "B3LYP QUADRUPOLE XY")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE XZ'), 0.000000000000, 4, "B3LYP QUADRUPOLE XZ")
    assert psi4.compare_values(psi4.get_variable('B3LYP QUADRUPOLE YZ'), -0.000000000000, 4, "B3LYP QUADRUPOLE YZ")
