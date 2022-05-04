from addons import *

@ctest_labeler("quick;noc1")
def test_python_databases():
    ctest_runner(__file__)

