import numpy as np
import pytest

from utils import *

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.mp2]


_ref_h2o_ccpvdz = {
    "df": {
        "HF TOTAL ENERGY": -76.0167614256151865,
        "MP2 SAME-SPIN CORRELATION ENERGY": -0.0527406422061238,
        "MP2 OPPOSITE-SPIN CORRELATION ENERGY": -0.1562926850310142,
        "MP2 CORRELATION ENERGY": -0.2090333272371381,
        "MP2 TOTAL ENERGY": -76.2257947528523232,
        "SCS-MP2 CORRELATION ENERGY": -0.2051314361059251,
        "SCS-MP2 TOTAL ENERGY": -76.2218928617211162,
    },
    "conv": {
        "HF TOTAL ENERGY": -76.01678947133706,
        "MP2 SAME-SPIN CORRELATION ENERGY": -0.05268120425816,
        "MP2 OPPOSITE-SPIN CORRELATION ENERGY": -0.15637564436589,
        "MP2 CORRELATION ENERGY": -0.20905684862405,
        "MP2 TOTAL ENERGY": -76.22584631996111,
        "SCS-MP2 CORRELATION ENERGY": -0.20521117465845,
        "SCS-MP2 TOTAL ENERGY": -76.22200064599551,
    },
}  # yapf: disable
for mp2type in ["df", "conv"]:
    _ref_h2o_ccpvdz[mp2type]["SCF TOTAL ENERGY"] = _ref_h2o_ccpvdz[mp2type]["HF TOTAL ENERGY"]
    _ref_h2o_ccpvdz[mp2type]["CURRENT REFERENCE ENERGY"] = _ref_h2o_ccpvdz[mp2type]["HF TOTAL ENERGY"]
    _ref_h2o_ccpvdz[mp2type]["5050SCS-MP2 CORRELATION ENERGY"] = 0.5 * (
        _ref_h2o_ccpvdz[mp2type]["MP2 SAME-SPIN CORRELATION ENERGY"]
        + _ref_h2o_ccpvdz[mp2type]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
    )
    _ref_h2o_ccpvdz[mp2type]["5050SCS-MP2 TOTAL ENERGY"] = (
        _ref_h2o_ccpvdz[mp2type]["5050SCS-MP2 CORRELATION ENERGY"] + _ref_h2o_ccpvdz[mp2type]["HF TOTAL ENERGY"]
    )


@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"name": "Mp2", "custom": "SCS-MP2", "options": {"mp2_type": "df"}}, id="mp2 (df)"),
        pytest.param({"name": "Mp2", "custom": "MP2", "options": {"mp2_type": "conv"}}, id="mp2 (conv)"),
        pytest.param(
            {
                "name": "Mp2",
                "custom": "SCS-MP2",
                "options": {"mp2_type": "df", "mp2_os_scale": 1.2, "mp2_ss_scale": 0.33333333},
            },
            id="explicit scs mp2 (df)",
        ),
        pytest.param(
            {
                "name": "Mp2",
                "custom": "SCS-MP2",
                "options": {"mp2_type": "conv", "os_scale": 1.2, "ss_scale": 0.33333333},
            },
            id="explicit scs mp2 (conv)",
        ),
        pytest.param(
            {
                "name": "Mp2",
                "custom": "5050SCS-MP2",
                "options": {"mp2_type": "df", "mp2_os_scale": 0.5, "mp2_ss_scale": 0.5},
            },
            id="user-def scs mp2 (df)",
        ),
        pytest.param(
            {"name": "Mp2", "custom": "5050SCS-MP2", "options": {"mp2_type": "conv", "os_scale": 0.5, "ss_scale": 0.5}},
            id="user-def scs mp2 (conv)",
        ),
    ],
)
def test_scsmp2(inp):
    """Formerly known as dfmp2-4"""

    h2o = psi4.geometry(
        """
        O
        H 1 1.0
        H 1 1.0 2 90.0
    """
    )

    psi4.set_options({"basis": "cc-pvdz"})
    psi4.set_options(inp["options"])

    ene, wfn = psi4.energy(inp["name"], return_wfn=True)

    ref_block = _ref_h2o_ccpvdz[inp["options"]["mp2_type"]]
    # ref_corl = ref_block[inp['pv'] + ' CORRELATION ENERGY']
    # ref_tot = ref_block[inp['pv'] + ' TOTAL ENERGY']
    ref_corl = ref_block["MP2 CORRELATION ENERGY"]
    ref_tot = ref_block["MP2 TOTAL ENERGY"]
    ref_custom_corl = ref_block[inp["custom"] + " CORRELATION ENERGY"]
    ref_custom_tot = ref_block[inp["custom"] + " TOTAL ENERGY"]

    for obj in [psi4.core, wfn]:
        for pv in [
            "HF TOTAL ENERGY",
            "SCF TOTAL ENERGY",
            "MP2 SAME-SPIN CORRELATION ENERGY",
            "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
            "MP2 CORRELATION ENERGY",
            "MP2 TOTAL ENERGY",
            "SCS-MP2 CORRELATION ENERGY",
            "SCS-MP2 TOTAL ENERGY",
            "CURRENT REFERENCE ENERGY",
        ]:

            assert compare_values(ref_block[pv], obj.variable(pv), 5, pv)

        if any((x in inp["options"] for x in ["os_scale", "ss_scale", "mp2_os_scale", "mp2_ss_scale"])):
            assert compare_values(
                ref_custom_corl, obj.variable("CUSTOM SCS-MP2 CORRELATION ENERGY"), 5, "custom scsmp2 corl"
            )
            assert compare_values(ref_custom_tot, obj.variable("CUSTOM SCS-MP2 TOTAL ENERGY"), 5, "custom scsmp2 ")

        assert compare_values(ref_corl, obj.variable("CURRENT CORRELATION ENERGY"), 5, "current corl")
        assert compare_values(ref_tot, obj.variable("CURRENT ENERGY"), 5, "current")

    assert compare_values(ref_tot, ene, 5, "return")
    assert compare_values(ref_tot, wfn.energy(), 5, "wfn")
