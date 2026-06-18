from addons import *
import pytest

@uusing("ecpint")
@ctest_labeler("properties;df;quick;ecp")
@pytest.mark.xfail(reason="MBIS incompatible with ECP: requires all-electron density. Use all-electron basis or reconstruct density with denspart.")
def test_mbis_ecp():
    ctest_runner(__file__)
