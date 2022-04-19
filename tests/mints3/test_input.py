from addons import *

@ctest_labeler("mints;noc1")
def test_mints3():
    ctest_runner(__file__, [
        "Lx-STO-3G.dat",
        "Ly-STO-3G.dat",
        "Lz-STO-3G.dat",
        "Lx-6-311Gss.dat",
        "Ly-6-311Gss.dat",
        "Lz-6-311Gss.dat",
        "Lx-cc-pVTZ.dat",
        "Ly-cc-pVTZ.dat",
        "Lz-cc-pVTZ.dat",
        "Px-STO-3G.dat",
        "Py-STO-3G.dat",
        "Pz-STO-3G.dat",
        "Px-6-311Gss.dat",
        "Py-6-311Gss.dat",
        "Pz-6-311Gss.dat",
        "Px-cc-pVTZ.dat",
        "Py-cc-pVTZ.dat",
        "Pz-cc-pVTZ.dat",
    ])

