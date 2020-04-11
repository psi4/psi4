import numpy as np
import pytest
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs

import psi4

from .standard_suite_runner import runner_asserter

pytestmark = [pytest.mark.quick, pytest.mark.skip]  # skip is temporary until references in place at qcng


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
# yapf: enable

_nyi1 = "fc conv mp2 gradients NYI"
_nyi2 = "rohf mp2 gradients NYI"
_nyi3 = "cd mp2 gradients NYI"
_nyi4 = "spin components rhf mp2 energies NYI"


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py?

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
        # below work fine but have to be careful b/c detci can't do fc, then ae w/o psio error, interfere w/detci lines below
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "pk",     },}, id="mp2  rhf   pk/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "direct", },}, id="mp2  rhf drct/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "df",     },}, id="mp2  rhf   df/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "mem_df", },}, id="mp2  rhf  mem/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "disk_df",},}, id="mp2  rhf disk/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",  "scf_type": "cd",     },}, id="mp2  rhf   cd/conv rr detci",),
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


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
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
        # * detci must have ae before fc to avoid psio error. cleaning doesn't help
        # * detci rohf mp2 does not match other programs in the stored reference
        #pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   detci"),
        #pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   detci",),
        #pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),
        #pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp2 diff ans", raises=AssertionError)),

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


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        ###### ccenergy
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd  rhf conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd  uhf conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },}, id="ccsd rohf conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd  rhf conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd  uhf conv ae: * ccenergy",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },}, id="ccsd rohf conv ae: * ccenergy",),
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


#    +----------------------+----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | name                 | |globals__qc_module| | :py:func:`~psi4.energy()`                                                | :py:func:`~psi4.gradient()`                                  |
#    +                      +                      +------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | _                    | |scf__reference|     | RHF                    | UHF                    | ROHF                   | RHF                | UHF                | ROHF               |
#    +                      +                      +------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | type select [#f1]_   | _                    | CV   | DF       | CD   | CV   | DF       | CD   | CV       | DF   | CD   | CV   | DF   | CD   | CV   | DF   | CD   | CV   | DF   | CD   |
#    +======================+======================+======+==========+======+======+==========+======+==========+======+======+======+======+======+======+======+======+======+======+======+
#    | .. _tlccsd:          | CCENERGY             | D    |          |      | D    |          |      | D        |      |      | D    |      |      | D    |      |      | D    |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | ccsd                 | DETCI                |      |          |      |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | |globals__cc_type|   | DFMP2                |      |          |      |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    |                      | FNOCC                | Y    | D        | D    |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    |                      | OCC                  |      | Y        | Y    |      |          |      |          |      |      |      | D    |      |      |      |      |      |      |      |
#    +----------------------+----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    +----------------------+----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | .. _tlccsdt:         | CCENERGY             | D    |          |      | D    |          |      | D        |      |      | D    |      |      | D    |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | ccsd(t)              | DETCI                |      |          |      |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    | |globals__cc_type|   | DFMP2                |      |          |      |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    |                      | FNOCC                | Y    | D        | D    |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +                      +----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+
#    |                      | OCC                  |      | Y        | Y    |      |          |      |          |      |      |      |      |      |      |      |      |      |      |      |
#    +----------------------+----------------------+------+----------+------+------+----------+------+----------+------+------+------+------+------+------+------+------+------+------+------+


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        # pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        # pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, mp2_type) work?

        ###### default qc_module
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2  rhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2  uhf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     },}, id="mp2 rohf    conv fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2  rhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2  uhf    conv ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    },}, id="mp2 rohf    conv ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2  rhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2  uhf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     },}, id="mp2 rohf    df   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2  rhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2  uhf    df   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    },}, id="mp2 rohf    df   ae: dd     ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2  rhf    cd   fc: dd     ", marks=pytest.mark.xfail(reason="no spin for rhf dfocc")),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2  uhf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     },}, id="mp2 rohf    cd   fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2  rhf    cd   ae: dd     ", marks=pytest.mark.xfail(reason="no spin for rhf dfocc")),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2  uhf    cd   ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    },}, id="mp2 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "true",                     },}, id="mp2  rhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "true",                     },}, id="mp2  uhf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "true",                     },}, id="mp2 rohf         fc: dd     ",),
        pytest.param({"keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },}, id="mp2  rhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },}, id="mp2  uhf         ae: dd     ",),
        pytest.param({"keywords": {"reference": "rohf",                                         "freeze_core": "false",                    },}, id="mp2 rohf         ae: dd     ",),
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


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
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
        pytest.xfail(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
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
        pytest.xfail(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
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
        pytest.xfail(inp["marks"][dertype])

    inpcopy["driver"] = "gradient"
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}

    runner_asserter(inpcopy, subject, method, basis, tnm)
