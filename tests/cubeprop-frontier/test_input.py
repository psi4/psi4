from addons import *

@ctest_labeler("cubeprop;noc1")
def test_cubeprop_frontier():
    ctest_runner(__file__, extra_infiles=[
        "CH2s_HOMO.cube.ref",
        "CH2s_LUMO.cube.ref",
        "CH2t_3_DOMO.cube.ref",
        "CH2t_4_SOMO.cube.ref",
        "CH2t_5_SOMO.cube.ref",
        "CH2t_6_LVMO.cube.ref",
    ])

