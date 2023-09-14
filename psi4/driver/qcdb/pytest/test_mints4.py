from utils import *

import qcdb


def test_mints4():
    #! A demonstration of mixed Cartesian/ZMatrix geometry specification, using variables, for
    #! the benzene-hydronium complex.  Atoms can be placed using ZMatrix coordinates, whether they belong
    #! to the same fragment or not.  Note that the Cartesian specification must come before the ZMatrix entries
    #! because the former define absolute positions, while the latter are relative.

    refENuc = 268.6171792624

    refGEOM = \
          [[   0.710500000000,    -0.794637665924,    -1.230622098778],
           [   1.421000000000,    -0.794637665924,     0.000000000000],
           [   0.710500000000,    -0.794637665924,     1.230622098778],
           [  -0.710500000000,    -0.794637665924,     1.230622098778],
           [   1.254500000000,    -0.794637665924,    -2.172857738095],
           [  -1.254500000000,    -0.794637665924,     2.172857738095],
           [  -0.710500000000,    -0.794637665924,    -1.230622098778],
           [  -1.421000000000,    -0.794637665924,     0.000000000000],
           [   2.509000000000,    -0.794637665924,     0.000000000000],
           [   1.254500000000,    -0.794637665924,     2.172857738095],
           [  -1.254500000000,    -0.794637665924,    -2.172857738095],
           [  -2.509000000000,    -0.794637665924,     0.000000000000],
           [   0.000000000000,     3.205362334076,     0.000000000000],
           [   0.494974746831,     3.555362334076,    -0.857321409974],
           [   0.494974746831,     3.555362334076,     0.857321409974],
           [  -0.989949493661,     3.555362334076,     0.000000000000]]

    dimer = qcdb.Molecule("""
    1 1
    # This part is just a normal Cartesian geometry specification for benzene
    C          0.710500000000    -0.794637665924    -1.230622098778
    C          1.421000000000    -0.794637665924     0.000000000000
    C          0.710500000000    -0.794637665924     1.230622098778
    C         -0.710500000000    -0.794637665924     1.230622098778
    H          1.254500000000    -0.794637665924    -2.172857738095
    H         -1.254500000000    -0.794637665924     2.172857738095
    C         -0.710500000000    -0.794637665924    -1.230622098778
    C         -1.421000000000    -0.794637665924     0.000000000000
    H          2.509000000000    -0.794637665924     0.000000000000
    H          1.254500000000    -0.794637665924     2.172857738095
    H         -1.254500000000    -0.794637665924    -2.172857738095
    H         -2.509000000000    -0.794637665924     0.000000000000
    # And the hydronium part is specified using a zmatrix, referencing the benzene coordinates
    X  1  CC  3  30   2  A2
    O  13 R   1  90   2  90
    H  14 OH  13 TDA  1  0
    H  14 OH  15 TDA  13 A1
    H  14 OH  15 TDA  13 -A1

    CC    = 1.421
    CH    = 1.088
    A1    = 120.0
    A2    = 180.0
    OH    = 1.05
    R     = 4.0
    units angstrom
    """)
    dimer.update_geometry()

    assert compare_values(refENuc, dimer.nuclear_repulsion_energy(), 9, "Bz-H3O+: nuclear repulsion energy")

    geom_now = qcdb.mscale(dimer.geometry(), qcdb.constants.bohr2angstroms)
    assert compare_matrices(refGEOM, geom_now, 6, "Bz-H3O+: geometry and orientation")


def test_scf4():
    #! RHF cc-pVDZ energy for water, automatically scanning the symmetric stretch and bending coordinates
    #! using Python's built-in loop mechanisms.  The geometry is apecified using a Z-matrix with variables
    #! that are updated during the potential energy surface scan, and then the same procedure is performed
    #! using polar coordinates, converted to Cartesian coordinates.
    import math

    refENuc = [
        9.78588587740,  9.780670144878629, 8.807297289661147, 8.802603130390768, 8.006633899691952,  8.002366482173423
    ]

    # Define the points on the potential energy surface using standard Python list functions
    Rvals = [0.9, 1.0, 1.1]
    Avals = range(102, 106, 2)

    # Start with a potentital energy scan in Z-matrix coordinates

    h2o = qcdb.Molecule("""
        O
        H 1 R
        H 1 R 2 A
    """)

    print("\n Testing Z-matrix coordinates\n")

    count = 0
    for R in Rvals:
        h2o.set_variable('R', R)  # alternately, h2o.R = R
        for A in Avals:
            h2o.A = A  # alternately, h2o.set_variable('A', A)
            h2o.update_geometry()
            assert compare_values(refENuc[count],
                                  h2o.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy %d" % count)
            count = count + 1

    # And now the same thing, using Python's trigonometry functions, and Cartesian input.  This time
    # we want to reset the Cartesian positions every time the angles and bond lengths change, so we
    # define the geometry inside the loops.  N.B. this requires the basis set to be re-specified after
    # every change of geometry

    print("\n Testing polar coordinates\n")

    count = 0
    for R in Rvals:
        for A in Avals:
            h2o = qcdb.Molecule("""
                O   0.0    0.0    0.0
                H   0.0      R    0.0
                H   0.0  RCosA  RSinA
            """)
            # The non-numeric entries above just define placeholders with names.  They still need
            # to be set, which we do below.
            h2o.R = R
            h2o.set_variable('RCosA', R * math.cos(math.radians(A)))
            h2o.RSinA = R * math.sin(math.radians(A))
            h2o.update_geometry()

            assert compare_values(refENuc[count],
                                  h2o.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy %d" % count)
            count = count + 1
