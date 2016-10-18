import qcdb

geom_h2oa = \
      [[   0.000000000000,     0.000000000000,    -0.079135765807],
       [   0.000000000000,    -0.707106781187,     0.627971015380],
       [   0.000000000000,     0.707106781187,     0.627971015380]]

geom_h2ob = \
      [[   0.000000000000,     0.000000000000,    -0.079135765807],
       [   0.000000000000,    -0.707106781187,     0.627971015380],
       [  -0.000000000000,     0.707106781187,     0.627971015380]]

geom_h2oc = \
      [[   0.000000000000,     0.000000000000,    -0.079135765807],
       [  -0.707106781187,    -0.000000000000,     0.627971015380],
       [   0.707106781187,     0.000000000000,     0.627971015380]]

geom_h2od = \
      [[   0.000000000000,     0.000000000000,     0.079135765807],
       [  -0.707106781187,    -0.000000000000,    -0.627971015380],
       [   0.707106781187,     0.000000000000,    -0.627971015380]]

geom_h2oe = \
      [[   0.000000000000,     0.000000000000,     0.079135765807],
       [   0.707106781187,     0.000000000000,    -0.627971015380],
       [  -0.707106781187,    -0.000000000000,    -0.627971015380]]

geom_h2of = \
      [[   0.000000000000,     0.000000000000,    -0.079135765807],
       [  -0.000000000000,     0.707106781187,     0.627971015380],
       [   0.000000000000,    -0.707106781187,     0.627971015380]]


h2oA = qcdb.Molecule("""
O
H 1 1.0
H 1 1.0 2 90.0
""")

h2oA.update_geometry()
#h2oA.print_out()
geom_now = qcdb.mscale(h2oA.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2oa, geom_now, 6, "H2O A geometry and orientation") #TEST



h2oB = qcdb.Molecule("""
           O          0.000000000000     0.000000000000    -0.079135765807
           H          0.000000000000    -0.707106781187     0.627971015380
           H          0.000000000000     0.707106781187     0.627971015380
""")
h2oB.update_geometry()
#h2oB.print_out()
geom_now = qcdb.mscale(h2oB.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2ob, geom_now, 6, "H2O B geometry and orientation") #TEST


h2oC = qcdb.Molecule("""
           O           0.000000000000    -0.079135765807  0.0
           H          -0.707106781187     0.627971015380  0.0
           H           0.707106781187     0.627971015380  0.0
""")
h2oC.update_geometry()
#h2oC.print_out()
geom_now = qcdb.mscale(h2oC.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2oc, geom_now, 6, "H2O C geometry and orientation") #TEST



h2oD = qcdb.Molecule("""
           O          0.000000000000  0.0    0.079135765807
           H         -0.707106781187  0.0   -0.627971015380
           H          0.707106781187  0.0   -0.627971015380
""")
h2oD.update_geometry()
#h2oD.print_out()
geom_now = qcdb.mscale(h2oD.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2od, geom_now, 6, "H2O D geometry and orientation") #TEST



h2oE = qcdb.Molecule("""
           O          1.0 0.0 0.0
           H          2.0 0.0 0.0
           H          1.0 -1.0 0.0
""")
h2oE.update_geometry()
#h2oE.print_out()
geom_now = qcdb.mscale(h2oE.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2oe, geom_now, 6, "H2O E geometry and orientation") #TEST



h2oF = qcdb.Molecule("""
           O          0.0 0.0 0.0
           H          0.0 0.0 1.0
           H          0.0 -1.0 0.0
""")
h2oF.update_geometry()
#h2oF.print_out()
geom_now = qcdb.mscale(h2oF.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_h2of, geom_now, 6, "H2O F geometry and orientation") #TEST


import qcdb

geom_hofa = \
       [[  0.527716669708,     0.027994251059,     0.000000000000],
        [  0.527716669708,    -0.972005748941,     0.000000000000],
        [ -0.472283330292,     0.027994251059,     0.000000000000]]

geom_hofb = \
       [[  0.353357110938,     0.392946960454,     0.000000000000],
        [  1.060463892125,    -0.314159820733,     0.000000000000],
        [ -0.353749670249,    -0.314159820733,     0.000000000000]]

geom_hofc = \
       [[ -0.353357110938,    -0.392946960454,     0.000000000000],
        [ -1.060463892125,     0.314159820733,     0.000000000000],
        [  0.353749670249,     0.314159820733,     0.000000000000]]

geom_hofd = \
       [[  0.353357110938,     0.392946960454,     0.000000000000],
        [  1.060463892125,    -0.314159820733,     0.000000000000],
        [ -0.353749670249,    -0.314159820733,     0.000000000000]]

geom_hofe = \
       [[ -0.027994251059,     0.527716669708,     0.000000000000],
        [  0.972005748941,     0.527716669708,     0.000000000000],
        [ -0.027994251059,    -0.472283330292,     0.000000000000]]

geom_hoff = \
       [[ -0.527716669708,     0.027994251059,     0.000000000000],
        [ -0.527716669708,    -0.972005748941,     0.000000000000],
        [  0.472283330292,     0.027994251059,     0.000000000000]]


hofA = qcdb.Molecule("""
O
H 1 1.0
F 1 1.0 2 90.0
""")

hofA.update_geometry()
#hofA.print_out()
geom_now = qcdb.mscale(hofA.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hofa, geom_now, 6, "HOF A geometry and orientation") #TEST



hofB = qcdb.Molecule("""
           O          0.000000000000     0.000000000000    -0.079135765807
           H          0.000000000000    -0.707106781187     0.627971015380
           F          0.000000000000     0.707106781187     0.627971015380
""")
hofB.update_geometry()
#hofB.print_out()
geom_now = qcdb.mscale(hofB.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hofb, geom_now, 6, "HOF B geometry and orientation") #TEST


hofC = qcdb.Molecule("""
           O           0.000000000000    -0.079135765807  0.0
           H          -0.707106781187     0.627971015380  0.0
           F           0.707106781187     0.627971015380  0.0
""")
hofC.update_geometry()
#hofC.print_out()
geom_now = qcdb.mscale(hofC.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hofc, geom_now, 6, "HOF C geometry and orientation") #TEST



hofD = qcdb.Molecule("""
           O          0.000000000000  0.0    0.079135765807
           H         -0.707106781187  0.0   -0.627971015380
           F          0.707106781187  0.0   -0.627971015380
""")
hofD.update_geometry()
#hofD.print_out()
geom_now = qcdb.mscale(hofD.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hofd, geom_now, 6, "HOF D geometry and orientation") #TEST



hofE = qcdb.Molecule("""
           O          1.0 0.0 0.0
           H          2.0 0.0 0.0
           F          1.0 -1.0 0.0
""")
hofE.update_geometry()
#hofE.print_out()
geom_now = qcdb.mscale(hofE.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hofe, geom_now, 6, "HOF E geometry and orientation") #TEST



hofF = qcdb.Molecule("""
           O          0.0 0.0 0.0
           H          0.0 0.0 1.0
           F          0.0 -1.0 0.0
""")
hofF.update_geometry()
#hofF.print_out()
geom_now = qcdb.mscale(hofF.geometry(), qcdb.psi_bohr2angstroms)
qcdb.compare_matrices(geom_hoff, geom_now, 6, "HOF F geometry and orientation") #TEST


