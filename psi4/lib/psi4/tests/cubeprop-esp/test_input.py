from addons import *

@ctest_labeler("quick;cubeprop")
def test_cubeprop_esp():
    ctest_runner(__file__, extra_infiles=[
        "ESP.cube.ref",
        "Dt.cube.ref",
    ])

