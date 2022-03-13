from pathlib import Path
from addons import *

@ctest_labeler("quick;ci;cas;properties;cart;noc1")
def test_ci_property():
    print(f"{__file__=}")
    print(f"{Path(__file__)=}")
    print(f"{Path(__file__).resolve()=}")
    ctest_runner(__file__, ["grid.dat"])

