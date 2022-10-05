from addons import *

@ctest_labeler("quick;smoke;json")
def test_json_schema_1_gradient():
    ctest_runner(__file__)

