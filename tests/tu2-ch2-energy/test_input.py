from pathlib import Path
from addons import *

@ctest_labeler("quick;tutorial")
def test_tu2_ch2_energy():
    print(f"{__file__=}")
    print(f"{Path(__file__)=}")
    print(f"{Path(__file__).resolve()=}")
    ctest_runner(__file__)

