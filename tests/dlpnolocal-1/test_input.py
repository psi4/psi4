from addons import *


@ctest_labeler("dlpno;mp2;localization;quick")
def test_dlpnolocal_1():
    ctest_runner(__file__)
