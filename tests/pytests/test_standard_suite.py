import numpy as np
import pytest
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs

import psi4

from .standard_suite_runner import runner_asserter


@pytest.fixture
def clsd_open_pmols():
    return {name: psi4.core.Molecule.from_string(smol, name=name) for name, smol in std_molecules.items()}


# yapf: disable
_p1 = (psi4.ValidationError, "please use SCF_TYPE = DF to automatically choose the correct DFJK implementation")
_p2f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'DFMP2'")
_p2a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'UHF' not directable to QC_MODULE 'DFMP2'")
_p3f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'DFMP2'")
_p3a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'DFMP2'")
_p4f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p4a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'DF', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p5_old = (RuntimeError, "Frozen core/virtual not implemented in OCC module analytic gradients")  # this, if you leave it to occ to raise NYI for RHF/UHF
_p5f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p5u = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_p6f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p6a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CONV', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p7f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p7a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p8f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_p8a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_p9f = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'True', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p9a = (psi4.ManagedMethodError, "Method 'mp2' with MP2_TYPE 'CD', FREEZE_CORE 'False', and REFERENCE 'ROHF' not directable to QC_MODULE 'OCC'")
_p10 = (psi4.ValidationError, "No analytic derivatives for SCF_TYPE CD.")
_p11 = (psi4.ValidationError, "gradients need DF-SCF reference.")
_p12 = (psi4.ValidationError, "gradients need conventional SCF reference.")
_p13 = (psi4.ManagedMethodError, "not directable to QC_MODULE 'FNOCC'")
_p14 = (psi4.ManagedMethodError, "not directable to QC_MODULE 'OCC'")
_p15 = (psi4.ValidationError, "Invalid scf_type for DFCC.")
_p16 = (psi4.ValidationError, "Frozen core is not available for the CC gradients.")
_p17 = (RuntimeError, "Frozen core/virtual not implemented in Orbital-optimized methods")
_p18 = (psi4.ManagedMethodError, "Method 'mp3' with MP_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p19 = (psi4.ManagedMethodError, "Method 'mp3' with MP_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_p20 = (psi4.ManagedMethodError, "Method 'mp2.5' with MP_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p21 = (psi4.ManagedMethodError, "Method 'mp2.5' with MP_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")
_p22 = (psi4.ManagedMethodError, "Method 'lccd' with CC_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'RHF' not directable to QC_MODULE 'OCC'")
_p23 = (psi4.ManagedMethodError, "Method 'lccd' with CC_TYPE 'CONV', FREEZE_CORE 'True', and REFERENCE 'UHF' not directable to QC_MODULE 'OCC'")

# yapf: enable

_nyi1 = pytest.mark.xfail(
    reason="fc conv mp2/mp2.5/mp3/lccd gradients NYI", raises=psi4.ManagedMethodError, strict=True
)
_nyi2 = pytest.mark.xfail(reason="rohf mp2/mp2.5/mp3/lccd gradients NYI", raises=psi4.ManagedMethodError, strict=True)
_nyi3 = pytest.mark.xfail(reason="cd mp2/mp2.5/mp3/lccd gradients NYI", raises=psi4.ManagedMethodError, strict=True)
_nyi4 = pytest.mark.xfail(reason="spin components rhf mp2/mp2.5/mp3/lccd energies NYI", raises=KeyError, strict=True)
_nyi5 = pytest.mark.xfail(reason="df/cd open-shell ccsd energies NYI", raises=psi4.ManagedMethodError, strict=True)
_nyi6 = pytest.mark.xfail(reason="rohf ccsd energies mp2 submethod NYI", raises=KeyError, strict=True)
_nyi7 = pytest.mark.xfail(reason="rohf lccd energies NYI", raises=psi4.ManagedMethodError, strict=True)
_nyi8 = pytest.mark.xfail(reason="rohf lccsd energies NYI", raises=psi4.ValidationError, strict=True)
_nyi9 = pytest.mark.xfail(reason="fc conv orbital-optimized energies NYI", raises=RuntimeError, strict=True)
_nyi10 = pytest.mark.xfail(reason="rohf olccd energies mp2 submethod NYI", raises=KeyError, strict=True)
_nyi11 = pytest.mark.xfail(reason="rohf mp2.5/mp3 energies NYI", raises=psi4.ManagedMethodError, strict=True)
_nyi12 = pytest.mark.xfail(reason="cd scf gradients NYI", raises=psi4.ValidationError, strict=True)

# http://patorjk.com/software/taag/#p=display&c=bash&f=Soft&t=MP3


#
#  ,--.  ,--.,------.    ,------.
#  |  '--'  ||  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .--.  ||  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |  |  ||  |`       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'  `--'`--'        `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                      `---' `---'
#  HF Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### scf_solver
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",                         },}, id="hf  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",                         },}, id="hf  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",                         },}, id="hf rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",                     },}, id="hf  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",                     },}, id="hf  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",                     },}, id="hf rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",                         },}, id="hf  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",                         },}, id="hf  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",                         },}, id="hf rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",                     },}, id="hf  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",                     },}, id="hf  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",                     },}, id="hf rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df",                    },}, id="hf  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df",                    },}, id="hf  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df",                    },}, id="hf rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",                         },}, id="hf  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",                         },}, id="hf  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",                         },}, id="hf rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_hf_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "hf"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    inpcopy["keywords"]["freeze_core"] = "false"

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.  ,--.,------.     ,----.                     ,--.,--.                 ,--.
#  |  '--'  ||  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .--.  ||  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |  |  ||  |`       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'  `--'`--'         `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  HF Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### scf_solver
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",                         },                     }, id="hf  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",                         },                     }, id="hf  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",                         },                     }, id="hf rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",                     },                     }, id="hf  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",                     },                     }, id="hf  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",                     },                     }, id="hf rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",                         },                     }, id="hf  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",                         },                     }, id="hf  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",                         },                     }, id="hf rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",                     },                     }, id="hf  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",                     },                     }, id="hf  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",                     },                     }, id="hf rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df",                    },                     }, id="hf  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df",                    },                     }, id="hf  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df",                    },                     }, id="hf rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",                         }, "marks": {1: _nyi12}}, id="hf  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",                         }, "marks": {1: _nyi12}}, id="hf  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",                         }, "marks": {1: _nyi12}}, id="hf rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_hf_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "hf"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    inpcopy["keywords"]["freeze_core"] = "false"

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------.  ,---.     ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  MP2 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/df   rr dfmp2",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _p1}, id="mp2  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2  rhf   cd/conv rr fnocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "pk",     },             }, id="mp2  rhf   pk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "direct", },             }, id="mp2  rhf drct/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "df",     },             }, id="mp2  rhf   df/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "mem_df", },             }, id="mp2  rhf  mem/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "disk_df",},             }, id="mp2  rhf disk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "cd",     },             }, id="mp2  rhf   cd/conv rr detci",),
        # yapf: enable
    ],
)
def test_mp2_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### dfmp2
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2  rhf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2  uhf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },}, id="mp2 rohf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2  rhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2  uhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },}, id="mp2 rohf    df   ae: * dfmp2",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2  rhf pk/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2  uhf pk/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "pk",},}, id="mp2 rohf pk/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2  rhf pk/df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2  uhf pk/df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "pk",},}, id="mp2 rohf pk/df   ae: * dfmp2",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2  rhf cd/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2  uhf cd/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",  "scf_type": "cd",},}, id="mp2 rohf cd/df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2  rhf cd/df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2  uhf cd/df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false", "scf_type": "cd",},}, id="mp2 rohf cd/df   ae: * dfmp2",),

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   fnocc",),

        ###### detci
        # * detci rohf mp2 does not match other programs in the stored reference
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   detci"),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   detci",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  rhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  uhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2 rohf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  rhf pk/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  uhf pk/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2 rohf pk/df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2  rhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2  uhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2 rohf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2  rhf cd/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2  uhf cd/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2 rohf cd/df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2 rohf    cd   ae: * dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  rhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2  uhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2 rohf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  rhf pk/cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2  uhf pk/cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2 rohf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp2_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2  rhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2  uhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2 rohf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2  rhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2  uhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2 rohf    conv ae: dd     ",            ),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2  rhf    df   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2  uhf    df   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2 rohf    df   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2  rhf    df   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2  uhf    df   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2 rohf    df   ae: dd     ",            ),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2  rhf    cd   fc: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2  uhf    cd   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2 rohf    cd   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2  rhf    cd   ae: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2  uhf    cd   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2 rohf    cd   ae: dd     ",            ),

        ###### default qc_module, mp2_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                     }}, id="mp2  rhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                     }}, id="mp2  uhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }}, id="mp2 rohf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                    }}, id="mp2  rhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                    }}, id="mp2  uhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }}, id="mp2 rohf         ae: dd     ",            ),
        # yapf: enable
    ],
)
def test_mp2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  MP2 Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?
        # * test ae and sole fc that differs

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _p11},       }, id="mp2  rhf   pk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _p11},       }, id="mp2  rhf drct/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp2  rhf   df/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "mem_df", },                           }, id="mp2  rhf  mem/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp2  rhf disk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp2  rhf   cd/df   rr dfmp2",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _p11},       }, id="mp2  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _p11},       }, id="mp2  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp2  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p1, 0: _p1},}, id="mp2  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp2  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp2  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _p5f},       }, id="mp2 rhf fc pk/conv rr occ  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                           }, id="mp2  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                           }, id="mp2  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     }, "error": {1: _p12},       }, id="mp2  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p12},       }, id="mp2  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",}, "error": {1: _p12},       }, id="mp2  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp2  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp2_gradient_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests

        ###### dfmp2
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   },                    }, id="mp2  rhf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _p2f},}, id="mp2  uhf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _p3f},}, id="mp2 rohf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },                    }, id="mp2  rhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _p2a},}, id="mp2  uhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _p3a},}, id="mp2 rohf    df   ae: * dfmp2",),

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p5f},}, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p5u},}, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p6f},}, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p6a},}, id="mp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p4f},}, id="mp2 rohf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p4a},}, id="mp2 rohf    df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p7f},}, id="mp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p8f},}, id="mp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p9f},}, id="mp2 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p7a},}, id="mp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p8a},}, id="mp2  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p9a},}, id="mp2 rohf    cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize(
    "dertype",
    [pytest.param(1, id="grd1", marks=pytest.mark.quick), pytest.param(0, id="grd0", marks=pytest.mark.long),],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp2_type) work? Here we xfail the NYI rather than catch graceful exit.

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},          }, id="mp2  rhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},          }, id="mp2  uhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  rhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  uhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf    conv ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  rhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  uhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  rhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  uhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf    df   ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi4},}, id="mp2  rhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},          }, id="mp2  uhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},          }, id="mp2 rohf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi4},}, id="mp2  rhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3},          }, id="mp2  uhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3},          }, id="mp2 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                     },                               }, id="mp2  rhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                     },                               }, id="mp2  uhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }, "marks": {1: _nyi2},          }, id="mp2 rohf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },                               }, id="mp2  rhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },                               }, id="mp2  uhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }, "marks": {1: _nyi2},          }, id="mp2 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp2_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------.  ,---.     ,-----.    ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    '--. `\    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-..--..--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----''--'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                          `---' `---'
#  MP2.5 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2.5  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2.5  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2.5  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _p1}, id="mp2.5  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2.5  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2.5  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp2.5  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp2.5  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp2.5  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp2.5  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp2.5  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp2.5  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp2p5_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  uhf    df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2.5  rhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2.5  uhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2.5  rhf pk/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2.5  uhf pk/df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2.5  rhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp2.5  uhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2.5  rhf cd/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp2.5  uhf cd/df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp2.5  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp2.5  uhf    cd   ae: * dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2.5  rhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp2.5  uhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2.5  rhf pk/cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp2.5  uhf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2p5_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp2.5  rhf    conv fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp2.5  uhf    conv fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp2.5 rohf    conv fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp2.5  rhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp2.5  uhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp2.5 rohf    conv ae: dd     ", marks=_nyi11),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp2.5  rhf    df   fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp2.5  uhf    df   fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp2.5 rohf    df   fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp2.5  rhf    df   ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp2.5  uhf    df   ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp2.5 rohf    df   ae: dd     ", marks=_nyi11),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp2.5  rhf    cd   fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp2.5  uhf    cd   fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp2.5 rohf    cd   fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp2.5  rhf    cd   ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp2.5  uhf    cd   ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp2.5 rohf    cd   ae: dd     ", marks=_nyi11),

        ###### default qc_module, mp_type
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }}, id="mp2.5  rhf         fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }}, id="mp2.5  uhf         fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }}, id="mp2.5 rohf         fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }}, id="mp2.5  rhf         ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }}, id="mp2.5  uhf         ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }}, id="mp2.5 rohf         ae: dd     ", marks=_nyi11),
        # yapf: enable
    ],
)
def test_mp2p5_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------.  ,---.     ,-----.     ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    |  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    '--. `\    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-..--..--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----''--'`----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  MP2.5 Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?
        # * test ae and sole fc that differs
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _p11},       }, id="mp2.5  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _p11},       }, id="mp2.5  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp2.5  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p1, 0: _p1},}, id="mp2.5  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp2.5  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp2.5  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _p20},       }, id="mp2.5 rhf fc pk/conv rr occ  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                           }, id="mp2.5  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                           }, id="mp2.5  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     }, "error": {1: _p12},       }, id="mp2.5  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p12},       }, id="mp2.5  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",}, "error": {1: _p12},       }, id="mp2.5  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp2.5  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp2p5_gradient_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests
        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p20},}, id="mp2.5  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p21},}, id="mp2.5  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2.5  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2.5  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2.5  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2.5  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2.5  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2.5  uhf    df   ae:   dfocc",),
        # yapf: enable
    ],
)
def test_mp2p5_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize(
    "dertype",
    [pytest.param(1, id="grd1", marks=pytest.mark.quick), pytest.param(0, id="grd0", marks=pytest.mark.long),],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp_type) work? Here we xfail the NYI rather than catch graceful exit.
        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="mp2.5  rhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="mp2.5  uhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                                }, id="mp2.5  rhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                                }, id="mp2.5  uhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf    conv ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp2.5  rhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                                }, id="mp2.5  uhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp2.5  rhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                                }, id="mp2.5  uhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf    df   ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi4}, }, id="mp2.5  rhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},           }, id="mp2.5  uhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi11},}, id="mp2.5 rohf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi4}, }, id="mp2.5  rhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3 },          }, id="mp2.5  uhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi11},}, id="mp2.5 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp2.5  rhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                    },                                }, id="mp2.5  uhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp2.5  rhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                                }, id="mp2.5  uhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp2.5 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp2p5_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp2.5"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------. ,----.     ,------.
#  |   `.'   ||  .--. ''.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  MP3 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp3  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp3  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp3  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _p1}, id="mp3  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp3  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp3  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp3  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="mp3  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="mp3  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp3  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp3  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp3  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },             }, id="mp3  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },             }, id="mp3  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },             }, id="mp3  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="mp3  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="mp3  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },             }, id="mp3  rhf   cd/conv rr fnocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "pk",     },             }, id="mp3  rhf   pk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "direct", },             }, id="mp3  rhf drct/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "df",     },             }, id="mp3  rhf   df/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "mem_df", },             }, id="mp3  rhf  mem/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "disk_df",},             }, id="mp3  rhf disk/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "cd",     },             }, id="mp3  rhf   cd/conv rr detci",),
        # yapf: enable
    ],
)
def test_mp3_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp3  rhf    conv fc:   fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp3  rhf    conv ae:   fnocc",),

        ###### detci
        # * detci rohf mp3 does not match other programs (cfour) in the stored reference
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp3  rhf    conv fc:   detci",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp3 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp3 diff ans", raises=AssertionError)),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp3  rhf    conv ae:   detci"),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp3 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp3 diff ans", raises=AssertionError)),

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  uhf    df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp3  rhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp3  uhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp3  rhf pk/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp3  uhf pk/df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp3  rhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="mp3  uhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp3  rhf cd/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="mp3  uhf cd/df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="mp3  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="mp3  uhf    cd   ae: * dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp3  rhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="mp3  uhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp3  rhf pk/cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="mp3  uhf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp3  rhf    conv fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp3  uhf    conv fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }}, id="mp3 rohf    conv fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp3  rhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp3  uhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }}, id="mp3 rohf    conv ae: dd     ", marks=_nyi11),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp3  rhf    df   fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp3  uhf    df   fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }}, id="mp3 rohf    df   fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp3  rhf    df   ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp3  uhf    df   ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }}, id="mp3 rohf    df   ae: dd     ", marks=_nyi11),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp3  rhf    cd   fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp3  uhf    cd   fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }}, id="mp3 rohf    cd   fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp3  rhf    cd   ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp3  uhf    cd   ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }}, id="mp3 rohf    cd   ae: dd     ", marks=_nyi11),

        ###### default qc_module, mp_type
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }}, id="mp3  rhf         fc: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }}, id="mp3  uhf         fc: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }}, id="mp3 rohf         fc: dd     ", marks=_nyi11),
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }}, id="mp3  rhf         ae: dd     ", marks=_nyi4 ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }}, id="mp3  uhf         ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }}, id="mp3 rohf         ae: dd     ", marks=_nyi11),
        # yapf: enable
    ],
)
def test_mp3_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.   ,--.,------. ,----.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' |  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  MP3 Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?
        # * test ae and sole fc that differs

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _p11},       }, id="mp3  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _p11},       }, id="mp3  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                           }, id="mp3  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p1, 0: _p1},}, id="mp3  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="mp3  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp3  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _p18},       }, id="mp3 rhf fc pk/conv rr occ  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                           }, id="mp3  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                           }, id="mp3  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     }, "error": {1: _p12},       }, id="mp3  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p12},       }, id="mp3  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",}, "error": {1: _p12},       }, id="mp3  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="mp3  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp3_gradient_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p18},}, id="mp3  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p19},}, id="mp3  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp3  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp3  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp3  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp3  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp3  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp3  uhf    df   ae:   dfocc",),
        # yapf: enable
    ],
)
def test_mp3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize(
    "dertype",
    [pytest.param(1, id="grd1", marks=pytest.mark.quick), pytest.param(0, id="grd0", marks=pytest.mark.long),],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp_type) work? Here we xfail the NYI rather than catch graceful exit.

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="mp3  rhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="mp3  uhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                                }, id="mp3  rhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                                }, id="mp3  uhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf    conv ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp3  rhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                                }, id="mp3  uhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp3  rhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                                }, id="mp3  uhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf    df   ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi4}, }, id="mp3  rhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},           }, id="mp3  uhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi11},}, id="mp3 rohf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi4}, }, id="mp3  rhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3 },          }, id="mp3  uhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi11},}, id="mp3 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp3  rhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                    },                                }, id="mp3  uhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }, "marks": {1: _nyi4, 0: _nyi4 },}, id="mp3  rhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                                }, id="mp3  uhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "marks": {1: _nyi2, 0: _nyi11},}, id="mp3 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp3_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "mp3"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.    ,-----. ,-----.,------.      ,------.
#  |  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                     `---' `---'
#  LCCD Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="lccd  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="lccd  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="lccd  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _p1}, id="lccd  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="lccd  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="lccd  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="lccd  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="lccd  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="lccd  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="lccd  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="lccd  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="lccd  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },             }, id="lccd  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },             }, id="lccd  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },             }, id="lccd  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="lccd  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="lccd  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },             }, id="lccd  rhf   cd/conv rr fnocc",),
        # yapf: enable
    ],
)
def test_lccd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="lccd  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="lccd  rhf    conv ae: * fnocc",),

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  rhf    conv fc:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  rhf    conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  uhf    df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="lccd  rhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="lccd  uhf pk/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="lccd  rhf pk/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="lccd  uhf pk/df   ae:   dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="lccd  rhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="lccd  uhf cd/df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="lccd  rhf cd/df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="lccd  uhf cd/df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="lccd  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="lccd  uhf    cd   ae: * dfocc",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="lccd  rhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="lccd  uhf pk/cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="lccd  rhf pk/cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="lccd  uhf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_lccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, cc_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccd  rhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccd  uhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccd rohf    conv fc: dd     ", marks=_nyi7),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccd  rhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccd  uhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccd rohf    conv ae: dd     ", marks=_nyi7),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccd  rhf    df   fc: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccd  uhf    df   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccd rohf    df   fc: dd     ", marks=_nyi7),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccd  rhf    df   ae: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccd  uhf    df   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccd rohf    df   ae: dd     ", marks=_nyi7),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccd  rhf    cd   fc: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccd  uhf    cd   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccd rohf    cd   fc: dd     ", marks=_nyi7),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccd  rhf    cd   ae: dd     ", marks=_nyi4),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccd  uhf    cd   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccd rohf    cd   ae: dd     ", marks=_nyi7),

        ###### default qc_module, cc_type
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }}, id="lccd  rhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }}, id="lccd  uhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }}, id="lccd rohf         fc: dd     ", marks=_nyi7),
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }}, id="lccd  rhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }}, id="lccd  uhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }}, id="lccd rohf         ae: dd     ", marks=_nyi7),
        # yapf: enable
    ],
)
def test_lccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.    ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  |  |   '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |   |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--.'  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `-----' `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  LCCD Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?
        # * test ae and sole fc that differs

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     }, "error": {1: _p11},       }, id="lccd  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", }, "error": {1: _p11},       }, id="lccd  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                           }, id="lccd  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p1, 0: _p1},}, id="lccd  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                           }, id="lccd  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="lccd  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _p22},       }, id="lccd rhf fc pk/conv rr occ  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                           }, id="lccd  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                           }, id="lccd  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     }, "error": {1: _p12},       }, id="lccd  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {1: _p12},       }, id="lccd  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",}, "error": {1: _p12},       }, id="lccd  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     }, "error": {1: _p10},       }, id="lccd  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_lccd_gradient_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p22},}, id="lccd  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p23},}, id="lccd  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="lccd  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="lccd  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="lccd  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="lccd  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="lccd  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="lccd  uhf    df   ae:   dfocc",),
        # yapf: enable
    ],
)
def test_lccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize(
    "dertype",
    [pytest.param(1, id="grd1", marks=pytest.mark.quick), pytest.param(0, id="grd0", marks=pytest.mark.long),],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, cc_type) work? Here we xfail the NYI rather than catch graceful exit.

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="lccd  rhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi1},           }, id="lccd  uhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                                }, id="lccd  rhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                                }, id="lccd  uhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf    conv ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi4, 0: _nyi4 },}, id="lccd  rhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                                }, id="lccd  uhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi4, 0: _nyi4 },}, id="lccd  rhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                                }, id="lccd  uhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf    df   ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi4 },}, id="lccd  rhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3},           }, id="lccd  uhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "marks": {1: _nyi3, 0: _nyi7 },}, id="lccd rohf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi4 },}, id="lccd  rhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3 },          }, id="lccd  uhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "marks": {1: _nyi3, 0: _nyi7 },}, id="lccd rohf    cd   ae: dd     ",),

        ###### default qc_module, cc_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "marks": {1: _nyi1},           }, id="lccd  rhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "marks": {1: _nyi1},           }, id="lccd  uhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                                }, id="lccd  rhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                                }, id="lccd  uhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "marks": {1: _nyi2, 0: _nyi7 },}, id="lccd rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_lccd_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#  ,--.    ,-----. ,-----. ,---.  ,------.      ,------.
#  |  |   '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  LCCSD Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },             }, id="lccsd  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },             }, id="lccsd  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },             }, id="lccsd  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },             }, id="lccsd  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},             }, id="lccsd  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },             }, id="lccsd  rhf   cd/conv rr fnocc",),
        # yapf: enable
    ],
)
def test_lccsd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="lccsd  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="lccsd  rhf    conv ae: * fnocc",),
        # yapf: enable
    ],
)
def test_lccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, cc_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccsd  rhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccsd  uhf    conv fc: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }}, id="lccsd rohf    conv fc: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccsd  rhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccsd  uhf    conv ae: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }}, id="lccsd rohf    conv ae: dd     ", marks=_nyi8),
        ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccsd  rhf    df   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccsd  uhf    df   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }}, id="lccsd rohf    df   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccsd  rhf    df   ae: dd     ",            ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccsd  uhf    df   ae: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }}, id="lccsd rohf    df   ae: dd     ",            ),
        ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccsd  rhf    cd   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccsd  uhf    cd   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }}, id="lccsd rohf    cd   fc: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccsd  rhf    cd   ae: dd     ",            ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccsd  uhf    cd   ae: dd     ",            ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }}, id="lccsd rohf    cd   ae: dd     ",            ),

        ###### default qc_module, cc_type
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }}, id="lccsd  rhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }}, id="lccsd  uhf         fc: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }}, id="lccsd rohf         fc: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }}, id="lccsd  rhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }}, id="lccsd  uhf         ae: dd     ", marks=_nyi8),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }}, id="lccsd rohf         ae: dd     ", marks=_nyi8),
        # yapf: enable
    ],
)
def test_lccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "lccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#   ,-----. ,-----. ,---.  ,------.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                      `---' `---'
#  CCSD Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd  rhf   pk/conv rr ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd  rhf drct/conv rr ccenergy",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd  rhf   df/conv rr ccenergy",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "mem_df", },               }, id="ccsd  rhf  mem/conv rr ccenergy",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd  rhf disk/conv rr ccenergy",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd  rhf   cd/conv rr ccenergy",),

        # TODO reconcile dispute over whether df+conv ccsd directive returns df+conv energy (fnocc) or conv+conv energy (ccenergy)

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd  rhf   pk/conv rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd  rhf drct/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd  rhf   df/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", },               }, id="ccsd  rhf  mem/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd  rhf disk/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd  rhf   cd/conv rr fnocc   ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     }, "error": _p15,}, id="ccsd  rhf   pk/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", }, "error": _p15,}, id="ccsd  rhf drct/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd  rhf   df/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", }, "error": _p1, }, id="ccsd  rhf  mem/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd  rhf disk/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd  rhf   cd/df   rr fnocc   ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd  rhf   pk/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd  rhf drct/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd  rhf   df/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "mem_df", }, "error": _p1, }, id="ccsd  rhf  mem/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd  rhf disk/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd  rhf   cd/df   rr dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "ccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### ccenergy
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd  rhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd  uhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd rohf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd  rhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd  uhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd rohf  conv ae: * ccenergy",),

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd  rhf  conv fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd  rhf  conv ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd  rhf    df fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd  rhf    df ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": _p15,}, id="ccsd  rhf pk/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": _p15,}, id="ccsd  rhf pk/df ae:   fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "cd",},               }, id="ccsd  rhf cd/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",},               }, id="ccsd  rhf cd/df ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd  rhf    cd fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd  rhf    cd ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": _p15,}, id="ccsd  rhf pk/cd fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": _p15,}, id="ccsd  rhf pk/cd ae:   fnocc   ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },               }, id="ccsd  rhf    df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },               }, id="ccsd  rhf    df ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},               }, id="ccsd  rhf pk/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},               }, id="ccsd  rhf pk/df ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},               }, id="ccsd  rhf cd/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},               }, id="ccsd  rhf cd/df ae:   dfocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },               }, id="ccsd  rhf    cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },               }, id="ccsd  rhf    cd ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},               }, id="ccsd  rhf pk/cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},               }, id="ccsd  rhf pk/cd ae:   dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "ccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, cc_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="ccsd  rhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="ccsd  uhf    conv fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }}, id="ccsd rohf    conv fc: dd     ", marks=_nyi6),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="ccsd  rhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="ccsd  uhf    conv ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }}, id="ccsd rohf    conv ae: dd     ", marks=_nyi6),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="ccsd  rhf    df   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="ccsd  uhf    df   fc: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }}, id="ccsd rohf    df   fc: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="ccsd  rhf    df   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="ccsd  uhf    df   ae: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }}, id="ccsd rohf    df   ae: dd     ", marks=_nyi5),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="ccsd  rhf    cd   fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="ccsd  uhf    cd   fc: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }}, id="ccsd rohf    cd   fc: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="ccsd  rhf    cd   ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="ccsd  uhf    cd   ae: dd     ", marks=_nyi5),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }}, id="ccsd rohf    cd   ae: dd     ", marks=_nyi5),

        ###### default qc_module, cc_type
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }}, id="ccsd  rhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }}, id="ccsd  uhf         fc: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }}, id="ccsd rohf         fc: dd     ", marks=_nyi6),
        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }}, id="ccsd  rhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }}, id="ccsd  uhf         ae: dd     ",            ),
        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }}, id="ccsd rohf         ae: dd     ", marks=_nyi6),
        # yapf: enable
    ],
)
def test_ccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "ccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#   ,-----. ,-----. ,---.  ,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  CCSD Gradient


@pytest.mark.parametrize("dertype", [pytest.param(1, id="grd1"), pytest.param(0, id="grd0", marks=pytest.mark.long),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        # * no mixed-type gradients available (like pk+df) so no grad tests

        ###### ccenergy
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true", }, "error": {1: _p16},}, id="ccsd  rhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true", }, "error": {1: _p16},}, id="ccsd  uhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true", }, "error": {1: _p16},}, id="ccsd rohf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",},                    }, id="ccsd  rhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",},                    }, id="ccsd  uhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",},                    }, id="ccsd rohf  conv ae: * ccenergy",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                    }, id="ccsd  rhf    df fc: * dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                    }, id="ccsd  rhf    df ae: * dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "ccsd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'          `-'                                    `---' `---'
#  CCSD(T) Energy


# @pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
# @pytest.mark.parametrize(
#    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
# )
# @pytest.mark.parametrize(
#    "inp",
#    [
#        # yapf: disable
#        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.
#
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd_t_  rhf   pk/conv rr ccenergy",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd_t_  rhf drct/conv rr ccenergy",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd_t_  rhf   df/conv rr ccenergy",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "mem_df", },               }, id="ccsd_t_  rhf  mem/conv rr ccenergy",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd_t_  rhf disk/conv rr ccenergy",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd_t_  rhf   cd/conv rr ccenergy",),
#
#        # TODO reconcile dispute over whether df+conv ccsd_t_ directive returns df+conv energy (fnocc) or conv+conv energy (ccenergy)
#
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd_t_  rhf   pk/conv rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd_t_  rhf drct/conv rr fnocc   ",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd_t_  rhf   df/conv rr fnocc   ",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", },               }, id="ccsd_t_  rhf  mem/conv rr fnocc   ",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd_t_  rhf disk/conv rr fnocc   ",),
#        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd_t_  rhf   cd/conv rr fnocc   ",),
#
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     }, "error": _p15,}, id="ccsd_t_  rhf   pk/df   rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", }, "error": _p15,}, id="ccsd_t_  rhf drct/df   rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd_t_  rhf   df/df   rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", }, "error": _p1, }, id="ccsd_t_  rhf  mem/df   rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd_t_  rhf disk/df   rr fnocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd_t_  rhf   cd/df   rr fnocc   ",),
#
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",     },               }, id="ccsd_t_  rhf   pk/df   rr dfocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "direct", },               }, id="ccsd_t_  rhf drct/df   rr dfocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "df",     },               }, id="ccsd_t_  rhf   df/df   rr dfocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "mem_df", }, "error": _p1, }, id="ccsd_t_  rhf  mem/df   rr dfocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "disk_df",},               }, id="ccsd_t_  rhf disk/df   rr dfocc   ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",     },               }, id="ccsd_t_  rhf   cd/df   rr dfocc   ",),
#        # yapf: enable
#    ],
# )
# def test_ccsd_prt_pr_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
#    method = "ccsd(t)"
#    tnm = request.node.name
#    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]
#
#    inpcopy = {k: v for k, v in inp.items()}
#    inpcopy["driver"] = "energy"
#    inpcopy["call"] = method
#    inpcopy["keywords"]["basis"] = basis
#
#    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### ccenergy
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd_t_  rhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd_t_  uhf  conv fc: * ccenergy",),
        # off from c4 by e-5/e-6 pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },               }, id="ccsd_t_ rohf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd_t_  rhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd_t_  uhf  conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },               }, id="ccsd_t_ rohf  conv ae: * ccenergy",),

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd_t_  rhf  conv fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd_t_  rhf  conv ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd_t_  rhf    df fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd_t_  rhf    df ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": _p15,}, id="ccsd_t_  rhf pk/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": _p15,}, id="ccsd_t_  rhf pk/df ae:   fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "cd",},               }, id="ccsd_t_  rhf cd/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",},               }, id="ccsd_t_  rhf cd/df ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",                   },               }, id="ccsd_t_  rhf    cd fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false",                  },               }, id="ccsd_t_  rhf    cd ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": _p15,}, id="ccsd_t_  rhf pk/cd fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": _p15,}, id="ccsd_t_  rhf pk/cd ae:   fnocc   ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },               }, id="ccsd_t_  rhf    df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },               }, id="ccsd_t_  rhf    df ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},               }, id="ccsd_t_  rhf pk/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},               }, id="ccsd_t_  rhf pk/df ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},               }, id="ccsd_t_  rhf cd/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},               }, id="ccsd_t_  rhf cd/df ae:   dfocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },               }, id="ccsd_t_  rhf    cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },               }, id="ccsd_t_  rhf    cd ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},               }, id="ccsd_t_  rhf pk/cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},               }, id="ccsd_t_  rhf pk/cd ae:   dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "ccsd(t)"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


# @pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
# @pytest.mark.parametrize(
#    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
# )
# @pytest.mark.parametrize(
#    "inp",
#    [
#        # yapf: disable
#        ######## Does the simple interface (default qc_module, scf_type, cc_type) work?
#
#        ###### default qc_module
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                }, id="ccsd_t_  rhf    conv fc: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                }, id="ccsd_t_  uhf    conv fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "marks": _nyi6,}, id="ccsd_t_ rohf    conv fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                }, id="ccsd_t_  rhf    conv ae: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                }, id="ccsd_t_  uhf    conv ae: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "marks": _nyi6,}, id="ccsd_t_ rohf    conv ae: dd     ",),
#        ####
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                }, id="ccsd_t_  rhf    df   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "marks": _nyi5,}, id="ccsd_t_  uhf    df   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "marks": _nyi5,}, id="ccsd_t_ rohf    df   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                }, id="ccsd_t_  rhf    df   ae: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "marks": _nyi5,}, id="ccsd_t_  uhf    df   ae: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "marks": _nyi5,}, id="ccsd_t_ rohf    df   ae: dd     ",),
#        ####
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                }, id="ccsd_t_  rhf    cd   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "marks": _nyi5,}, id="ccsd_t_  uhf    cd   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "marks": _nyi5,}, id="ccsd_t_ rohf    cd   fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                }, id="ccsd_t_  rhf    cd   ae: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "marks": _nyi5,}, id="ccsd_t_  uhf    cd   ae: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "marks": _nyi5,}, id="ccsd_t_ rohf    cd   ae: dd     ",),
#
#        ###### default qc_module, cc_type
#        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                }, id="ccsd_t_  rhf         fc: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "true",                     },                }, id="ccsd_t_  uhf         fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "marks": _nyi6,}, id="ccsd_t_ rohf         fc: dd     ",),
#        pytest.param({"keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                }, id="ccsd_t_  rhf         ae: dd     ",),
#        pytest.param({"keywords": {"reference": "uhf",                                         "freeze_core": "false",                    },                }, id="ccsd_t_  uhf         ae: dd     ",),
#        pytest.param({"keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "marks": _nyi6,}, id="ccsd_t_ rohf         ae: dd     ",),
#        # yapf: enable
#    ],
# )
# def test_ccsd_prt_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
#    method = "ccsd(t)"
#    tnm = request.node.name
#    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]
#
#    inpcopy = {k: v for k, v in inp.items()}
#    if inp.get("marks", False):
#        request.node.add_marker(inp["marks"][dertype])
#    inpcopy["driver"] = "energy"
#    inpcopy["call"] = method
#    inpcopy["keywords"]["basis"] = basis
#
#    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#   ,-----. ,--.    ,-----. ,-----.,------.      ,------.
#  '  .-.  '|  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                              `---' `---'
#  OLCCD Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="olccd  rhf   pk/df   rr dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="olccd  rhf drct/df   rr dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="olccd  rhf   df/df   rr dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": _p1}, id="olccd  rhf  mem/df   rr dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="olccd  rhf disk/df   rr dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="olccd  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },             }, id="olccd  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },             }, id="olccd  rhf drct/conv rr occ  ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },             }, id="olccd  rhf   df/conv rr occ  ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },             }, id="olccd  rhf  mem/conv rr occ  ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},             }, id="olccd  rhf disk/conv rr occ  ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },             }, id="olccd  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_olccd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "olccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": _p17}, id="olccd  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": _p17}, id="olccd  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": _p17}, id="olccd rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },              }, id="olccd  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },              }, id="olccd  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },              }, id="olccd rohf    conv ae: * occ  ",),
        #
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd  rhf  v conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd  uhf  v conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd rohf  v conv ae:   occ  ",),
        # ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd  rhf    df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd  uhf    df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd rohf    df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd  rhf    df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd  uhf    df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd rohf    df   ae:   dfocc",),
        # ##
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd  rhf pk/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd  uhf pk/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd rohf pk/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd  rhf pk/df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd  uhf pk/df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd rohf pk/df   ae:   dfocc",),
        # ##
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="olccd  rhf cd/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="olccd  uhf cd/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "cd",  },}, id="olccd rohf cd/df   fc:   dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="olccd  rhf cd/df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="olccd  uhf cd/df   ae:   dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "false", "scf_type": "cd",  },}, id="olccd rohf cd/df   ae:   dfocc",),
        # ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd  rhf    cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd  uhf    cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="olccd rohf    cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd  rhf    cd   ae: * dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd  uhf    cd   ae: * dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="olccd rohf    cd   ae: * dfocc",),
        # ##
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd  rhf pk/cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd  uhf pk/cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",  "scf_type": "pk",  },}, id="olccd rohf pk/cd   fc: * dfocc",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd  rhf pk/cd   ae: * dfocc",),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd  uhf pk/cd   ae: * dfocc",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false", "scf_type": "pk",  },}, id="olccd rohf pk/cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_olccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "olccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, cc_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="olccd  rhf    conv fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }}, id="olccd  uhf    conv fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }}, id="olccd rohf    conv fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="olccd  rhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }}, id="olccd  uhf    conv ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }}, id="olccd rohf    conv ae: dd     ", marks=_nyi10),
        ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="olccd  rhf    df   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }}, id="olccd  uhf    df   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }}, id="olccd rohf    df   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="olccd  rhf    df   ae: dd     ",             ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }}, id="olccd  uhf    df   ae: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }}, id="olccd rohf    df   ae: dd     ",             ),
        ####
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="olccd  rhf    cd   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }}, id="olccd  uhf    cd   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }}, id="olccd rohf    cd   fc: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="olccd  rhf    cd   ae: dd     ",             ),
        # pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }}, id="olccd  uhf    cd   ae: dd     ",             ),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }}, id="olccd rohf    cd   ae: dd     ",             ),

        ###### default qc_module, cc_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                     }}, id="olccd  rhf         fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                     }}, id="olccd  uhf         fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }}, id="olccd rohf         fc: dd     ", marks=_nyi9 ),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                    }}, id="olccd  rhf         ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                    }}, id="olccd  uhf         ae: dd     ",             ),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }}, id="olccd rohf         ae: dd     ", marks=_nyi10),
        # yapf: enable
    ],
)
def test_olccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    method = "olccd"
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["driver"] = "energy"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis

    runner_asserter(inpcopy, subject, method, basis, tnm)
