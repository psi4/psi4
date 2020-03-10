import numpy as np
import pytest

from .utils import *
from .standard_suite_ref import std_suite, answer_hash

import psi4

pytestmark = [pytest.mark.quick, pytest.mark.mp2]


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
)  # yapf: disable
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


@pytest.fixture
def clsd_open_pmols():
    mols = {
        "hf": psi4.core.Molecule.from_string(
            """
                H
                F 1 0.917
              """
        ),
        "bh3p": psi4.core.Molecule.from_string(
            """
                1 2
                B     0.10369114     0.00000000     0.00000000
                H    -1.13269886     0.00000000     0.00000000
                H     3.00000000     0.37149000     0.00000000
                H     3.00000000    -0.37149000     0.00000000
              """
        ),
    }
    return mols


@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene  rhf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2 ene  rhf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene  uhf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  uhf df/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2 ene  uhf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene rohf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene rohf df/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2 ene rohf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 ene  rhf pk/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2 ene  rhf    df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 ene  uhf pk/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  uhf df/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2 ene  uhf    df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 ene rohf pk/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene rohf df/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2 ene rohf    df   ae: * dfmp2",),

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",                   },}, id="mp2 ene  rhf    df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  uhf df/df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",                   },}, id="mp2 ene  uhf    df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene rohf df/df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "true",                   },}, id="mp2 ene rohf    df   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/df   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",                  },}, id="mp2 ene  rhf    df   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  uhf df/df   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",                  },}, id="mp2 ene  uhf    df   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene rohf df/df   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",                  },}, id="mp2 ene rohf    df   ae: * dfocc",),
        # dfocc can't handle pk refs underneath df mp2

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene  rhf pk/conv fc:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/conv fc:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp2 ene  rhf    conv fc:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 ene  rhf pk/conv ae:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/conv ae:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp2 ene  rhf    conv ae:   fnocc",),
 
        # detci must have ae before fc to avoid psio error. cleaning doesn't help
        pytest.param({'driver': 'energy', 'subject': 'hf',   'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'detci', 'freeze_core': 'false', 'scf_type': 'pk',},}, id='mp2 ene  rhf pk/conv ae:   detci'),
        pytest.param({'driver': 'energy', 'subject': 'hf',   'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'detci', 'freeze_core': 'false', 'scf_type': 'df',},}, id='mp2 ene  rhf df/conv ae:   detci'),
        pytest.param({'driver': 'energy', 'subject': 'hf',   'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'detci', 'freeze_core': 'false',                  },}, id='mp2 ene  rhf    conv ae:   detci'),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene  rhf pk/conv fc:   detci",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/conv fc:   detci",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2 ene  rhf    conv fc:   detci",),
      #  pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/conv fc:   detci",),
      #  pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene rohf df/conv fc:   detci",),
##ans        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 ene rohf    conv fc:   detci",),
#long        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 ene rohf    conv ae:   detci",),

 #       pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  rhf    conv fc: * occ",),
     #   pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  rhf    df   fc:   occ",),
#        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  rhf    cd   fc: * occ",),
 #       pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  uhf    conv fc: * occ",),
     #   pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  uhf    df   fc:   occ",),
#        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene  uhf    cd   fc: * occ",),
 #!       pytest.param({'driver': 'energy', 'subject': 'bh3p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ', 'freeze_core': 'true',                     },}, id='mp2 ene rohf    conv fc: * occ',), #
     #   pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene rohf    df   fc:   occ",),
#        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 ene rohf    cd   fc: * occ",),
 #       pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  rhf    conv ae: * occ",),
     #   pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  rhf    df   ae:   occ",),
#        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  rhf    cd   ae: * occ",),
 #       pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  uhf    conv ae: * occ",),
     #   pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  uhf    df   ae:   occ",),
#        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene  uhf    cd   ae: * occ",),
 #!       pytest.param({'driver': 'energy', 'subject': 'bh3p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ', 'freeze_core': 'false',                    },}, id='mp2 ene rohf    conv ae: * occ',),  #
     #   pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene rohf    df   ae:   occ",),
#        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 ene rohf    cd   ae: * occ",),

        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/conv fc: * occ",),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/df   fc:   occ",),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  rhf df/cd   fc: * occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  uhf df/conv fc: * occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  uhf df/df   fc:   occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene  uhf df/cd   fc: * occ",),
        #pytest.param({'driver': 'energy', 'subject': 'bh3p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ', 'freeze_core': 'true',  'scf_type': 'df',},}, id='mp2 ene rohf df/conv fc: * occ',), #
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene rohf df/df   fc:   occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "df",},}, id="mp2 ene rohf df/cd   fc: * occ",),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/conv ae: * occ",),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/df   ae:   occ",),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  rhf df/cd   ae: * occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  uhf df/conv ae: * occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  uhf df/df   ae:   occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene  uhf df/cd   ae: * occ",),
        #pytest.param({'driver': 'energy', 'subject': 'bh3p', 'options': {'reference': 'rohf', 'mp2_type': 'conv', 'qc_module': 'occ', 'freeze_core': 'false', 'scf_type': 'df',},}, id='mp2 ene rohf df/conv ae: * occ',),  #
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene rohf df/df   ae:   occ",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "df",},}, id="mp2 ene rohf df/cd   ae: * occ",),

#        # scf_type set for indexing, not for driver
#       # pytest.param({'driver': 'gradient', 'subject': 'hf',   'options': {'reference': 'rhf',  'mp2_type': 'conv', 'qc_module': 'occ',   'freeze_core': 'true',  'scf_type': 'pk'}}, id='mp2 grad  rhf conv fc: * occ'),
#        pytest.param({"driver": "gradient", "subject": "hf",   "options": { "reference": "rhf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false", "scf_type": "pk", }, }, id="mp2 grad  rhf conv ae: * occ",),
#        pytest.param({"driver": "gradient", "subject": "hf",   "options": { "reference": "rhf", "mp2_type": "df", "qc_module": "occ", "freeze_core": "true", "scf_type": "df", }, }, id="mp2 grad  rhf df   fc:   occ",),
#        pytest.param({"driver": "gradient", "subject": "hf",   "options": { "reference": "rhf", "mp2_type": "df", "qc_module": "dfmp2", "freeze_core": "true", "scf_type": "df", }, }, id="mp2 grad  rhf df   fc: * dfmp2",),
#        pytest.param({"driver": "gradient", "subject": "hf",   "options": { "reference": "rhf", "mp2_type": "df", "qc_module": "occ", "freeze_core": "false", "scf_type": "df", }, }, id="mp2 grad  rhf df   ae:   occ",),
#        pytest.param({"driver": "gradient", "subject": "hf",   "options": { "reference": "rhf", "mp2_type": "df", "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "df", }, }, id="mp2 grad  rhf df   ae: * dfmp2",),
#        pytest.param({"driver": "gradient", "subject": "bh3p", "options": { "reference": "uhf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false", "scf_type": "pk", }, }, id="mp2 grad  uhf conv ae: * occ",),
#        pytest.param({"driver": "gradient", "subject": "bh3p", "options": { "reference": "uhf", "mp2_type": "df", "qc_module": "occ", "freeze_core": "true", "scf_type": "df", }, }, id="mp2 grad  uhf df   fc: * occ",),
#        pytest.param({"driver": "gradient", "subject": "bh3p", "options": { "reference": "uhf", "mp2_type": "df", "qc_module": "occ", "freeze_core": "false", "scf_type": "df", }, }, id="mp2 grad  uhf df   ae: * occ",),
#        # For gradients, this method would be found in the procedures table but would return a ManagedMethodError from proc.py. Normally, this would confuse the analytic-or-findif logic in gradient().
#        #   This scenario is now managed, and this test case makes sure its routing stays managed.
#        pytest.param( { "driver": "gradient", "subject": "bh3p", "options": { "reference": "rohf", "mp2_type": "df", "qc_module": "occ", "freeze_core": "false", "scf_type": "df", "points": 5, }, }, id="mp2 grad rohf df   ae: findif",),
    ],
)  # yapf: disable
def test_mp2_module(inp, clsd_open_pmols, request):
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]
    basis = "cc-pvdz"
    atol = 1.0e-6

    frz = "fc" if inp["options"]["freeze_core"] == "true" else "ae"
    scf_type = inp["options"].get("scf_type", "pk" if inp["options"]["mp2_type"] == "conv" else "df")
    chash = answer_hash(system=inp["subject"], basis=basis, fc=frz, scf_type=scf_type, reference=inp["options"]["reference"], mp2_type=inp["options"]["mp2_type"])
    print("BLOCK", chash)
    ref_block = std_suite[chash]

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 1.0e-3
    chash_conv = answer_hash(
        system=inp["subject"],
        basis=basis,
        fc=frz,
        reference=inp["options"]["reference"],
        mp2_type="conv",
        scf_type="pk",
    )
    ref_block_conv = std_suite[chash_conv]

    psi4.set_options(
        {
            "basis": basis,
            "guess": "sad",
            #'e_convergence': 8,
            #'d_convergence': 7,
            "e_convergence": 10,
            "d_convergence": 9,
        }
    )
    psi4.set_options(inp["options"])

    if inp["driver"] == "energy":
        ene, wfn = psi4.energy("mp2", molecule=subject, return_wfn=True)

    elif inp["driver"] == "gradient":
        grad, wfn = psi4.gradient("mp2", molecule=subject, return_wfn=True)

    # qcvars
    for obj in [psi4.core, wfn]:
        for pv in [
            "HF TOTAL ENERGY",
            "MP2 CORRELATION ENERGY",
            "MP2 TOTAL ENERGY",
#            "MP2 SINGLES ENERGY",
#            "MP2 SAME-SPIN CORRELATION ENERGY",
#            "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
        ]:
            # known contract violators
            if (
                (inp["options"]["qc_module"] == "detci" and pv in ["MP2 SINGLES ENERGY", "MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]) or
                (inp["options"]["qc_module"] == "occ" and inp["options"]["reference"] == "rhf" and inp["options"]["mp2_type"] == "df" and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"])
            ):
                break
            assert compare_values(ref_block[pv], obj.variable(pv), tnm + " " + pv, atol=atol)
            assert compare_values(ref_block_conv[pv], obj.variable(pv), tnm + " " + pv, atol=atol_conv)

        # aliases
        pv = "SCF TOTAL ENERGY"
        assert compare_values(ref_block["HF TOTAL ENERGY"], obj.variable(pv), tnm + " " + pv, atol=atol)

        pv = "CURRENT REFERENCE ENERGY"
        assert compare_values(ref_block["HF TOTAL ENERGY"], obj.variable(pv), tnm + " " + pv, atol=atol)

        pv = "CURRENT CORRELATION ENERGY"
        assert compare_values(ref_block["MP2 CORRELATION ENERGY"], obj.variable(pv), tnm + " " + pv, atol=atol)

    # returns
    assert compare_values(ref_block["MP2 TOTAL ENERGY"], wfn.energy(), tnm + " wfn", atol=atol)

    if inp["driver"] == "energy":
        assert compare_values(ref_block["MP2 TOTAL ENERGY"], ene, tnm + " return")

    elif inp["driver"] == "gradient":
        assert compare_values(ref_block["MP2 TOTAL GRADIENT"], wfn.gradient().np, tnm + " grad wfn", atol=atol)
        assert compare_values(ref_block["MP2 TOTAL GRADIENT"], grad.np, tnm + " grad return", atol=atol)
