from addons import *

@ctest_labeler("quick;smoke;fcidump;noc1")
def test_fcidump():
    ctest_runner(__file__, [
        "Ne.6311G.INTDUMP.ref",
        "Ne.C1.6311G.INTDUMP.ref",
        "Ne.cc-pVDZ.UHF.INTDUMP.ref",
        "Ne.C1.cc-pVDZ.UHF.INTDUMP.ref",
    ])

