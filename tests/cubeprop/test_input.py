from addons import *

@ctest_labeler("quick;cubeprop")
def test_cubeprop():
    ctest_runner(__file__, [
        "Psi_a_1_1-A1.cube.ref",
        "Psi_a_2_2-A1.cube.ref",
        "Psi_a_3_1-B2.cube.ref",
        "Psi_a_4_3-A1.cube.ref",
        "Psi_a_5_1-B1.cube.ref",
        "DUAL.cube.ref",
        "Da.cube.ref",
        "Db.cube.ref",
        "Ds.cube.ref",
        "Dt.cube.ref",
    ])

