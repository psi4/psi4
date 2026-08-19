from addons import *


@ctest_labeler("scf;embpot;gradient")
def test_embpot1_puream():
    ctest_runner(__file__)
