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


# yapf: disable
_e1 = (psi4.ValidationError, "please use SCF_TYPE = DF to automatically choose the correct DFJK implementation")
_e2f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'DFMP2'")
_e2a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'UHF' not directable to QC_MODULE 'DFMP2'")
_e3f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'DFMP2'")
_e3a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'DFMP2'")
_e4f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e4a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e5_old = (RuntimeError, "Frozen core/virtual not implemented in OCC module analytic gradients")  # this, if you leave it to occ to raise NYI for RHF/UHF
_e5f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_e5u = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_e6f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e6a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e7f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_e7a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_e8f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_e8a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_e9f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e9a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_e10 = (psi4.ValidationError, "No analytic derivatives for SCF_TYPE CD.")
_e11 = (psi4.ValidationError, "gradients need DF-SCF reference.")
_e12 = (psi4.ValidationError, "gradients need conventional SCF reference.")
# yapf: enable

_nyi1 = "fc conv mp2 gradients NYI"
_nyi2 = "rohf mp2 gradients NYI"
_nyi3 = "cd mp2 gradients NYI"
_nyi4 = "spin components rhf mp2 energies NYI"


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/df   rr dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/df   rr dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/df   rr dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/df   rr dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/df   rr dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/df   rr dfmp2",),

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/df   rr dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/df   rr dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/df   rr dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _e1}, id="mp2  rhf  mem/df   rr dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/df   rr dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/df   rr dfocc",),

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/conv rr occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/conv rr occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/conv rr occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/conv rr occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/conv rr occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/conv rr occ  ",),

        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/conv rr fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/conv rr fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/conv rr fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/conv rr fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/conv rr fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/conv rr fnocc",),
        # below work fine but have to be careful b/c detci can't do fc, then ae w/o psio error, interfere w/detci lines below
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "pk",     },}, id="mp2  rhf   pk/conv rr detci",),
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "direct", },}, id="mp2  rhf drct/conv rr detci",),
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "df",     },}, id="mp2  rhf   df/conv rr detci",),
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "mem_df", },}, id="mp2  rhf  mem/conv rr detci",),
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "disk_df",},}, id="mp2  rhf disk/conv rr detci",),
        # pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "cd",     },}, id="mp2  rhf   cd/conv rr detci",),
        # yapf: enable
    ],
)
def test_mp2_energy_scftype(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]
    _runner_asserter(inp, subject, method, tnm)


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### dfmp2
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2  rhf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2  uhf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2 rohf    df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2  rhf    df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2  uhf    df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2 rohf    df   ae: * dfmp2",),
        ##
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2  rhf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2  uhf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 rohf pk/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2  rhf pk/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2  uhf pk/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 rohf pk/df   ae: * dfmp2",),
        ##
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2  rhf cd/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2  uhf cd/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2 rohf cd/df   fc: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2  rhf cd/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2  uhf cd/df   ae: * dfmp2",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2 rohf cd/df   ae: * dfmp2",),

        ###### fnocc
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   fnocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   fnocc",),

        ###### detci
        # * detci must have ae before fc to avoid psio error. cleaning doesn't help
        # * detci rohf mp2 does not match other programs in the stored reference
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   detci"),
        #pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   detci",),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),
        #pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),

        ###### occ/dfocc
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    df   ae:   dfocc",),
        ##
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  rhf pk/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  uhf pk/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2 rohf pk/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  rhf pk/df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  uhf pk/df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2 rohf pk/df   ae:   dfocc",),
        ##
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2  rhf cd/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2  uhf cd/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2 rohf cd/df   fc:   dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2  rhf cd/df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2  uhf cd/df   ae:   dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2 rohf cd/df   ae:   dfocc",),
        ####
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    cd   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    cd   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    cd   ae: * dfocc",),
        ##
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  rhf pk/cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  uhf pk/cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2 rohf pk/cd   fc: * dfocc",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  rhf pk/cd   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  uhf pk/cd   ae: * dfocc",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2 rohf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2_energy_module(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]
    _runner_asserter(inp, subject, method, tnm)


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### ccenergy
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd  rhf conv fc: * ccenergy",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd  uhf conv fc: * ccenergy",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd rohf conv fc: * ccenergy",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd  rhf conv ae: * ccenergy",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd  uhf conv ae: * ccenergy",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd rohf conv ae: * ccenergy",),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, clsd_open_pmols, request):
    method = "ccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]
    _runner_asserter(inp, subject, method, tnm)


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp2_type) work?

        ###### default qc_module
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2  rhf    conv fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2  uhf    conv fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2 rohf    conv fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2  rhf    conv ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2  uhf    conv ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2 rohf    conv ae: dd     ",),
        ####
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2  rhf    df   fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2  uhf    df   fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2 rohf    df   fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2  rhf    df   ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2  uhf    df   ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2 rohf    df   ae: dd     ",),
        ####
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2  rhf    cd   fc: dd     ", marks=pytest.mark.xfail(reason="no spin for rhf dfocc")),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2  uhf    cd   fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2 rohf    cd   fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2  rhf    cd   ae: dd     ", marks=pytest.mark.xfail(reason="no spin for rhf dfocc")),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2  uhf    cd   ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",                                          "freeze_core": "true",                     },}, id="mp2  rhf         fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",                                          "freeze_core": "true",                     },}, id="mp2  uhf         fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf",                                         "freeze_core": "true",                     },}, id="mp2 rohf         fc: dd     ",),
        pytest.param({"driver": "energy", "subject": "hf",   "options": {"reference": "rhf",                                          "freeze_core": "false",                    },}, id="mp2  rhf         ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "uhf",                                          "freeze_core": "false",                    },}, id="mp2  uhf         ae: dd     ",),
        pytest.param({"driver": "energy", "subject": "bh3p", "options": {"reference": "rohf",                                         "freeze_core": "false",                    },}, id="mp2 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp2_energy_default(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]
    _runner_asserter(inp, subject, method, tnm)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?
        # * test ae and sole fc that differs

        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _e11},       }, id="mp2  rhf   pk/df   rr dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _e11},       }, id="mp2  rhf drct/df   rr dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp2  rhf   df/df   rr dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "mem_df", },                           }, id="mp2  rhf  mem/df   rr dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp2  rhf disk/df   rr dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _e10},       }, id="mp2  rhf   cd/df   rr dfmp2",),

        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _e11},       }, id="mp2  rhf   pk/df   rr dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _e11},       }, id="mp2  rhf drct/df   rr dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp2  rhf   df/df   rr dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _e1, 0: _e1},}, id="mp2  rhf  mem/df   rr dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp2  rhf disk/df   rr dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _e10},       }, id="mp2  rhf   cd/df   rr dfocc",),

        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _e5f},       }, id="mp2 rhf fc pk/conv rr occ  ",),
        ##
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                           }, id="mp2  rhf   pk/conv rr occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                           }, id="mp2  rhf drct/conv rr occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     }, "error": {1: _e12},       }, id="mp2  rhf   df/conv rr occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _e12},       }, id="mp2  rhf  mem/conv rr occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",}, "error": {1: _e12},       }, id="mp2  rhf disk/conv rr occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _e10},       }, id="mp2  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp2_gradient_scftype(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])

    _runner_asserter(inpcopy, subject, method, tnm, dertype=dertype)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests

        ###### dfmp2
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },                    }, id="mp2  rhf    df   fc: * dfmp2",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _e2f},}, id="mp2  uhf    df   fc: * dfmp2",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _e3f},}, id="mp2 rohf    df   fc: * dfmp2",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },                    }, id="mp2  rhf    df   ae: * dfmp2",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _e2a},}, id="mp2  uhf    df   ae: * dfmp2",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _e3a},}, id="mp2 rohf    df   ae: * dfmp2",),

        ###### occ/dfocc
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e5f},}, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e5u},}, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e6f},}, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _e6a},}, id="mp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  rhf    df   fc:   dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  uhf    df   fc:   dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e4f},}, id="mp2 rohf    df   fc:   dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    df   ae:   dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    df   ae:   dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _e4a},}, id="mp2 rohf    df   ae:   dfocc",),
        ####
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e7f},}, id="mp2  rhf    cd   fc: * dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e8f},}, id="mp2  uhf    cd   fc: * dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e9f},}, id="mp2 rohf    cd   fc: * dfocc",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _e7a},}, id="mp2  rhf    cd   ae: * dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _e8a},}, id="mp2  uhf    cd   ae: * dfocc",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _e9a},}, id="mp2 rohf    cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])

    _runner_asserter(inpcopy, subject, method, tnm, dertype=dertype)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp2_type) work? Here we xfail the NYI rather than catch graceful exit.

        ###### default qc_module
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},          }, id="mp2  rhf    conv fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},          }, id="mp2  uhf    conv fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf    conv fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  rhf    conv ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  uhf    conv ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf    conv ae: dd     ",),
        ####
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  rhf    df   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  uhf    df   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf    df   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  rhf    df   ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  uhf    df   ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf    df   ae: dd     ",),
        ####
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi4},}, id="mp2  rhf    cd   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},          }, id="mp2  uhf    cd   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},          }, id="mp2 rohf    cd   fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi4},}, id="mp2  rhf    cd   ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3},          }, id="mp2  uhf    cd   ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3},          }, id="mp2 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",                                          "freeze_core": "true",                     },                               }, id="mp2  rhf         fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",                                          "freeze_core": "true",                     },                               }, id="mp2  uhf         fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf",                                         "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf         fc: dd     ",),
        pytest.param({"driver": "gradient", "subject": "hf",   "options": {"reference": "rhf",                                          "freeze_core": "false",                    },                               }, id="mp2  rhf         ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",                                          "freeze_core": "false",                    },                               }, id="mp2  uhf         ae: dd     ",),
        pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "rohf",                                         "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp2_gradient_default(inp, dertype, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[inp["subject"]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])

    _runner_asserter(inpcopy, subject, method, tnm, dertype=dertype)
