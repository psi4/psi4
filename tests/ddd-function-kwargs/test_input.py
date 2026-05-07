from addons import *

@uusing("qcmanybody")
@ctest_labeler("findif;gradient")
def test_ddd_function_kwargs():
    ctest_runner(__file__)
