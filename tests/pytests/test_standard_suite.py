import numpy as np
import pytest
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs

import psi4

from standard_suite_runner import runner_asserter
from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.fixture
def clsd_open_pmols():
    return {name: psi4.core.Molecule.from_string(smol, name=name) for name, smol in std_molecules.items()}

# ADVICE: tests in this file run through methods and check the usual return contracts are fulfilled.
# * tests are arranged by (1) plain methods, then orbital-optimized methods, then DFT methods, (2) method in level-of-theory order, (3) derivative/driver order, (4) "scftype", then "module", then "default"
#   * search sections with `<<<` or scroll for big-letter headers
#   * for new big-letter headers, use:  http://patorjk.com/software/taag/#p=display&c=bash&f=Soft&t=MP3
# * at best, each method+driver has three `def test_<method>_<driver>_<testtype>(...)` test types:
#   * scftype: least important, skip for gradient, write second
#     * purpose: check that if user sets a scf_type, corl method either uses it, computes missing integrals files, or raises an helpful error
#     * scope: only checks rhf ref and dz basis. also, only conv and df, not cd.
#     * a block has single module and type (conv|df). each block cycles through the six common scf_type values
#
#   * module: most important, write first
#     * purpose: check that all routes to computing a method produce the same final result and that all the subcontracts (like storing in wfn and setting reference energy) are fulfilled or excused
#     * scope: full scope, so checks three molecule+basis combinations and all analytic/findif routes
#     * block count: varies by module redundancy and capability
#     * each block fixes: qc_module and corl type and scf_type (>= corl type, so scf+corl combinations among conv+conv, conv+cd, cd+cd, conv+df, cd+df, df+df; mixed are optional)
#     * each block cycles: through up to six ref+aefc combinations
#       * implemented combinations should always be present. ae/fc pairs should always be present
#       * not-yet-implemented combinations may be present to regression test the error messages or may be skipped
#     * checks whether QC values consistent among modules, whether specified ref, types, aefc ran (via QC values), whether specified qc_module ran, whether intended error message thrown
#     * populates the module lines for each method in the capabilities table
#     * findif parameterization important when setting up test to be sure analytic matches 5-point findif. but redundant for day-to-day, so suppressed by default for testing efficiency. use `pytest --runnonroutine` to actually run them
#
#   * default: second-most important, write last
#     * purpose: check that if user doesn't set qc_module or doesn't set corl_type, the calc runs or raises a helpful error
#     * scope: only checks dz basis
#     * block count: four corl type specifications: conv, df, cd, and unspecified
#     * each block fixes: corl type
#     * each block cycles: through six ref+aefc combinations
#     * checks QC values, checks module and type match expectations
#     * always 24 pytest.param lines: 4 "type" blocks namely conv, df, cd, and unspecified (to catch mp2_type, cc_type, etc.). each block with fc r/u/ro and ae r/u/ro
#     * when analytic derivatives missing, new dz entries may be needed at qcng standard_suite_ref.py for findif to check against
#     * populates the summary, non-module line for each method in the capabilities table and also the single- and double-underline default markers
#
# * it's tempting to parametrize and consolidate tests further, but keep in mind that the hard part is managing the differing outcomes (e.g., grd1 vs. grd0)
# * note that doing "good" things like adding a method (particularly to a select_ function) or fixing a missing QCVariable can cause failures here. that's a feature, not a bug.

_MethodError = (psi4.MissingMethodError, psi4.ManagedMethodError)

# yapf: disable
# * tuple is (error_type, string_match_in_error_message, reason_for_human)
# * note that 2nd is regex matched, so raw strings and escape chars may be needed
# depending on DDD routing, _p1/_p10/_p29/_p31 can be `psi4.ValidationError` or `qcengine.exceptions.InputError`
# depending on DDD routing, _p17 can be `RuntimeError` or `qcengine.exeptions.UnknownError`
_p1 = (Exception, "please use SCF_TYPE = DF to automatically choose the correct DFJK implementation", "no mem_df in dfocc")
_p2 = (psi4.MissingMethodError, r"Method=mp2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP2_TYPE=DF, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=DFMP2", "no df open-shell gradients for mp2 by dfmp2")
_p4 = (psi4.MissingMethodError, r"Method=mp2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP2_TYPE=DF, REFERENCE=ROHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=OCC", "no df rohf gradients for mp2 by occ")
_p5_old = (RuntimeError, "Frozen core/virtual not implemented in OCC module analytic gradients")  # this, if you leave it to occ to raise NYI for RHF/UHF
_p5 = (psi4.MissingMethodError, r"Method=mp2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP2_TYPE=CONV, REFERENCE=(R|U|RO)HF, FREEZE_CORE=TRUE, QC_MODULE=(OCC|\(auto\))", "no fc conv gradients for mp2 by occ")
_p6 = (psi4.MissingMethodError, r"Method=mp2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP2_TYPE=(CONV|DF|CD), REFERENCE=ROHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=(OCC|\(auto\))", "no rohf gradients for mp2 by occ")
_p7 = (psi4.MissingMethodError, r"Method=mp2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP2_TYPE=CD, REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=(OCC|\(auto\))", "no cd gradients for mp2 by occ")
_p10 = (Exception, "No analytic derivatives for SCF_TYPE CD.", "no cd scf gradients to underlie post-scf")
_p11 = (psi4.ValidationError, "gradients need DF-SCF reference.", "conv scf needed to underlie df post-scf for mp2/2.5/3/lccd")
_p12 = (psi4.ValidationError, "gradients need conventional SCF reference.", "no df scf for conv gradients for mp2/2.5/3/lccd in occ")
_p15 = (psi4.ValidationError, r"Invalid scf_type='(PK|DIRECT)' for FNOCC energy through `run_fnodfcc`.", "no conventional scf for df/cd cc in fnocc")
_p16 = (psi4.ValidationError, "Frozen core is not available for the CC gradients.", "no fc gradients by ccenergy")
_p17 = (Exception, r"Frozen core\/virtual not implemented in Orbital-optimized methods", "no fc/fv for oo in occ")
_p18 = (psi4.MissingMethodError, r"Method=mp3 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP_TYPE=CONV, REFERENCE=(R|U)HF, FREEZE_CORE=TRUE, QC_MODULE=(OCC|\(auto\))", "no fc conv gradients for mp3 by occ")
_p20 = (psi4.MissingMethodError, r"Method=mp2.5 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions MP_TYPE=CONV, REFERENCE=(R|U)HF, FREEZE_CORE=TRUE, QC_MODULE=(OCC|\(auto\))", "no fc conv gradients for mp2.5 by occ")
_p22 = (psi4.MissingMethodError, r"Method=lccd is not available for requested derivative level \(reqd=1 > avail=0\) under conditions CC_TYPE=CONV, REFERENCE=(R|U)HF, FREEZE_CORE=TRUE, QC_MODULE=(OCC|\(auto\))", "no fc conv gradients for lccd by occ")
_p24 = (psi4.ValidationError, r"Method=mp2 is not available for requested derivative level \(reqd=2 > avail=(1|0)\) under any conditions.", "no hessians for mp2")
_p25 = (psi4.ValidationError, "Only RHF/UHF/RKS/UKS Hessians are currently implemented. SCF_TYPE either CD or OUT_OF_CORE not supported", "no rohf hessians for hf")
_p26 = (psi4.MissingMethodError, r"Method=ccd is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=(R|U)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=OCC", "no conv for ccd by occ")
_p27 = (_MethodError, r"Method=cc2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions CC_TYPE=CONV, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=(CCENERGY|\(auto\))", "no open-shell gradients for cc2 by ccenergy")
_p29 = (Exception, "ROHF reference for DFT is not available.", "no rohf for dft")
_p30 = (psi4.MissingMethodError, r"not available for requested derivative level \(reqd=1 > avail=0\) under conditions SCF_TYPE=CD, REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE)", "no cd gradients")
_p31 = (Exception, r"SCF_TYPE \(CD\) not supported for range-separated functionals", "no cd for lrc in dft")
_p32 = (psi4.MissingMethodError, r"Method=(mp2.5|mp3|remp2|lccd) is not available for any derivative level under conditions (MP|CC)_TYPE=(CONV|DF|CD), REFERENCE=ROHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no rohf mp2.5/mp3/remp2/lccd by occ")
_p33 = (psi4.MissingMethodError, r"Method=(mp2.5|mp3|lccd|ccd|ccsd|ccsd\(t\)|omp2|omp2.5|omp3|olccd|oremp2) is not available for requested derivative level \(reqd=1 > avail=0\) under conditions (MP2|MP|CC)_TYPE=CD, REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no cd gradients by occ")
_p34 = (psi4.MissingMethodError, r"Method=(ccd|ccsd|ccsd\(t\)|a-ccsd\(t\)) is not available for any derivative level under conditions CC_TYPE=(DF|CD), REFERENCE=ROHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no df/cd rohf ccd/ccsd/ccsd(t)/a-ccsd(t) by occ")
_p35 = (psi4.MissingMethodError, r"Method=ccsd\(t\) is not available for requested derivative level \(reqd=1 > avail=0\) under conditions CC_TYPE=CONV, REFERENCE=ROHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=(CCENERGY|\(auto\))", "no conv rohf gradients for ccsd(t) by ccenergy")
_p36 = (psi4.ValidationError, r"Invalid type (DF|CD) for CCENERGY", "no df/cd cc2/cc3/bccd/bccd(t) by psi4")
_p37 = (psi4.MissingMethodError, r"Method=(cc2|cc3) is not available for any derivative level under conditions CC_TYPE=(DF|CD), REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no df/cd cc2/cc3 by psi4")
_p38 = (_MethodError, r"Method=(a-ccsd\(t\)) is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no conv open-shell a-ccsd(t) by ccenergy")
_p39 = (psi4.MissingMethodError, r"Method=(ccd) is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no conv ccd by psi4")
_p40 = (psi4.ValidationError, r"Invalid reference type (U|RO)HF != RHF for FNOCC energy", "no open-shell energies in fnocc")
_p41 = (psi4.MissingMethodError, r"Method=(mp4\(sdq\)|mp4|qcisd|qcisd\(t\)) is not available for any derivative level under conditions (CI|MP)_TYPE=(DF|CD), REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no df/cd mp4(sdq)/mp4/qcisd/qcisd(t) by fnocc")
_p42 = (psi4.MissingMethodError, r"Method=(mp4\(sdq\)|mp4|qcisd|qcisd\(t\)) is not available for any derivative level under conditions (CI|MP)_TYPE=CONV, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no open-shell mp4(sdq)/mp4/qcisd/qcisd(t) by fnocc")
_p43 = (psi4.ValidationError, r"(CD|DF) is not a valid choice", "no df/cd ci by psi4")
_p44 = (psi4.ValidationError, r"Invalid type (DF|CD) for FNOCC energy through `run_(fnocc|cepa)`.", "no df/cd except ccsd/ccsd(t) by fnocc")
_p45 = (psi4.ValidationError, r"Reference UHF for DETCI is not available", "no uhf by detci")
_p46 = (psi4.ValidationError, r"Invalid type (DF|CD) for DETCI energy through `run_detci`", "no df/cd by detci")
_p47 = (psi4.UpgradeHelper, r"Replace method ZAPT with method MP for RHF reference", "retire rhf zapt by detci")
_p48 = (psi4.MissingMethodError, r"Method=cisd is not available for any derivative level under conditions CI_TYPE=CONV, REFERENCE=UHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "no uhf cisd by psi4")

_p98 = (psi4.MissingMethodError, r"Method=ccsd\(t\) is not available for requested derivative level \(reqd=1 > avail=0\) under conditions CC_TYPE=CONV, REFERENCE=(R|U)HF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "temporary: cc(t) disabled w/o qc_module=ccenergy for conv rhf/uhf gradients for ccsd(t) by ccenergy until scaling reworked")
_p99 = (psi4.MissingMethodError, r"Method=(ccsd\(t\)|a-ccsd\(t\)) is not available for any derivative level under conditions CC_TYPE=(DF|CD), REFERENCE=UHF, FREEZE_CORE=(TRUE|FALSE), QC_MODULE=\(auto\)", "temporary: cc(t) disabled w/o qc_module=occ in dfocc until further optimization ")

_w1 = ("CCSD TOTAL GRADIENT",  "wrong semicanonical ae rohf ccsd gradients by ccenergy")
# yapf: enable


#
#  ,--.  ,--.,------.    ,------.
#  |  '--'  ||  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .--.  ||  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |  |  ||  |`       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'  `--'`--'        `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                      `---' `---'
#  <<<  HF Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false",                    },}, id="hf  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false",                    },}, id="hf  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false",                    },}, id="hf rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false",                    },}, id="hf  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false",                    },}, id="hf  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false",                    },}, id="hf rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false",                    },}, id="hf  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false",                    },}, id="hf  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false",                    },}, id="hf rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false",                    },}, id="hf  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false",                    },}, id="hf  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false",                    },}, id="hf rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false",                    },}, id="hf  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false",                    },}, id="hf  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false",                    },}, id="hf rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false",                    },}, id="hf  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false",                    },}, id="hf  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false",                    },}, id="hf rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_hf_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "energy"))


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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",                     }}, id="hf  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",                     }}, id="hf  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",                     }}, id="hf rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",                     }}, id="hf  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",                     }}, id="hf  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",                     }}, id="hf rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",                     }}, id="hf  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",                     }}, id="hf  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",                     }}, id="hf rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",                     }}, id="hf  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",                     }}, id="hf  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                       "freeze_core": "false",                     }}, id="hf rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_hf_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "energy"))


#
#  ,--.  ,--.,------.     ,----.                     ,--.,--.                 ,--.
#  |  '--'  ||  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .--.  ||  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |  |  ||  |`       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'  `--'`--'         `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  HF Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false"                   },                     }, id="hf  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false"                   },                     }, id="hf  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false"                   },                     }, id="hf rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false"                   },                     }, id="hf  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false"                   },                     }, id="hf  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false"                   },                     }, id="hf rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false"                   },                     }, id="hf  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false"                   },                     }, id="hf  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false"                   },                     }, id="hf rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false"                   },                     }, id="hf  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false"                   },                     }, id="hf  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false"                   },                     }, id="hf rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false"                   },                     }, id="hf  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false"                   },                     }, id="hf  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false"                   },                     }, id="hf rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false"                   }, "error": {1: _p10}  }, id="hf  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false"                   }, "error": {1: _p10}  }, id="hf  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false"                   }, "error": {1: _p10}  }, id="hf rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_hf_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false",                    },                               }, id="hf  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false",                    },                               }, id="hf  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false",                    },                               }, id="hf rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false",                    },                               }, id="hf  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false",                    },                               }, id="hf  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false",                    },                               }, id="hf rohf    df   ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p10},           }, id="hf  rhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p10},           }, id="hf  uhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p10},           }, id="hf rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },                               }, id="hf  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },                               }, id="hf  uhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    },                               }, id="hf rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_hf_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "gradient"))


#
#  ,--.  ,--.,------.    ,--.  ,--.                     ,--.
#  |  '--'  ||  .---'    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  .--.  ||  `--,     |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |  |  ||  |`       |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'  `--'`--'        `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  HF Hessian


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2", marks=pytest.mark.d2ints),
        pytest.param(1, id="hes1", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
        pytest.param(0, id="hes0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false",                   },                     }, id="hf  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false",                   },                     }, id="hf  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false",                   }, "error": {2: _p25}  }, id="hf rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false",                   },                     }, id="hf  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false",                   },                     }, id="hf  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false",                   }, "error": {2: _p25}  }, id="hf rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false",                   },                     }, id="hf  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false",                   },                     }, id="hf  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false",                   }, "error": {2: _p25}  }, id="hf rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false",                   },                     }, id="hf  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false",                   },                     }, id="hf  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false",                   }, "error": {2: _p25}  }, id="hf rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false",                   },                     }, id="hf  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false",                   },                     }, id="hf  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false",                   }, "error": {2: _p25}  }, id="hf rohf disk ae:   scf  ",),
        ####
        # only H-by-E available for CD ref Hessians, and loose dflt cholesky_tolerance means they're not close to CONV (at least for dz uhf/rohf, so skipping for now
        # pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false",                   }, "error": {2: _p10, 1: _p10}}, id="hf  rhf   cd ae:   scf  ",),
        # pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false",                   }, "error": {2: _p10, 1: _p10}}, id="hf  uhf   cd ae:   scf  ",),
        # pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false",                   }, "error": {2: _p10, 1: _p10}}, id="hf rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_hf_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "hessian"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2", marks=pytest.mark.d2ints),
        pytest.param(1, id="hes1", marks=pytest.mark.findif),
        pytest.param(0, id="hes0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false",                    },                               }, id="hf  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false",                    },                               }, id="hf  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false",                    }, "error": {2: _p25},           }, id="hf rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false",                    },                               }, id="hf  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false",                    },                               }, id="hf  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false",                    }, "error": {2: _p25},           }, id="hf rohf    df   ae: dd     ",),
        ####
        # only H-by-E available for CD ref Hessians, and loose dflt cholesky_tolerance means they're not close to CONV (at least for dz uhf/rohf, so skipping for now
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {2: _p10, 1: _p10},  }, id="hf  rhf    cd   ae: dd     ",),
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {2: _p10, 1: _p10},  }, id="hf  uhf    cd   ae: dd     ",),
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false",                    }, "error": {2: _p10, 1: _p10},  }, id="hf rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },                               }, id="hf  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },                               }, id="hf  uhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }, "error": {2: _p25},           }, id="hf rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_hf_hessian_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hf", "hessian"))


#
#  ,--.   ,--.,------.  ,---.     ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP2 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp2  rhf   pk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp2  rhf drct/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp2  rhf   df/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="mp2  rhf  mem/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp2  rhf disk/df   rr dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp2  rhf   cd/df   rr dfmp2",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp2  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp2  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp2  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {0: _p1}}, id="mp2  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp2  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp2  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp2  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp2  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp2  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="mp2  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp2  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp2  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp2  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp2  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp2  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="mp2  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp2  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp2  rhf   cd/conv rr fnocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "pk",     },                  }, id="mp2  rhf   pk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "direct", },                  }, id="mp2  rhf drct/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "df",     },                  }, id="mp2  rhf   df/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "mem_df", },                  }, id="mp2  rhf  mem/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "disk_df",},                  }, id="mp2  rhf disk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "cd",     },                  }, id="mp2  rhf   cd/conv rr detci",),
        # yapf: enable
    ],
)
def test_mp2_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "energy"))


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
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp2  rhf    conv ae:   detci"),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp2  rhf    conv fc:   detci",),

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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "energy"))


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
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2  rhf    conv fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2  uhf    conv fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }}, id="mp2 rohf    conv fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2  rhf    conv ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2  uhf    conv ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }}, id="mp2 rohf    conv ae: dd     ",            ),
        ####
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2  rhf    df   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2  uhf    df   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }}, id="mp2 rohf    df   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2  rhf    df   ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2  uhf    df   ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }}, id="mp2 rohf    df   ae: dd     ",            ),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2  rhf    cd   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2  uhf    cd   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }}, id="mp2 rohf    cd   fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2  rhf    cd   ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2  uhf    cd   ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }}, id="mp2 rohf    cd   ae: dd     ",            ),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rhf",                                          "freeze_core": "true",                     }}, id="mp2  rhf         fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "uhf",                                          "freeze_core": "true",                     }}, id="mp2  uhf         fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }}, id="mp2 rohf         fc: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    }}, id="mp2  rhf         ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    }}, id="mp2  uhf         ae: dd     ",            ),
        pytest.param({"xptd": {"qc_module": "dfmp2"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }}, id="mp2 rohf         ae: dd     ",            ),
        # yapf: enable
    ],
)
def test_mp2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "energy"))


#
#  ,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",   "freeze_core": "true",  "scf_type": "pk",      }, "error": {1: _p5},        }, id="mp2 rhf fc pk/conv rr occ  ",),
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _p2}, }, id="mp2  uhf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "true",                   }, "error": {1: _p2}, }, id="mp2 rohf    df   fc: * dfmp2",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  },                    }, id="mp2  rhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _p2}, }, id="mp2  uhf    df   ae: * dfmp2",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "dfmp2", "freeze_core": "false",                  }, "error": {1: _p2}, }, id="mp2 rohf    df   ae: * dfmp2",),

        ###### occ/dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p5}, }, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p5}, }, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p5}, }, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p6}, }, id="mp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="mp2  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p4}, }, id="mp2 rohf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="mp2  uhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p4}, }, id="mp2 rohf    df   ae:   dfocc",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p7}, }, id="mp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p7}, }, id="mp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p7}, }, id="mp2 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p7}, }, id="mp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p7}, }, id="mp2  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    }, "error": {1: _p7}, }, id="mp2 rohf    cd   ae: * dfocc",),
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p5},            }, id="mp2  rhf    conv fc: dd     ",),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p5},            }, id="mp2  uhf    conv fc: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p6},            }, id="mp2 rohf    conv fc: dd     ",),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                               }, id="mp2  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p6},            }, id="mp2 rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "dfmp2"},    "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  rhf    df   fc: dd     ",),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                               }, id="mp2  uhf    df   fc: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p6},            }, id="mp2 rohf    df   fc: dd     ",),
        pytest.param({"xptd": {"qc_module": "dfmp2"},    "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                               }, id="mp2  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p6},            }, id="mp2 rohf    df   ae: dd     ",),
        ####
        pytest.param({                                   "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p7},            }, id="mp2  rhf    cd   fc: dd     ",),
        pytest.param({                                   "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p7},            }, id="mp2  uhf    cd   fc: dd     ",),
        pytest.param({                                   "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p7},            }, id="mp2 rohf    cd   fc: dd     ",),
        pytest.param({                                   "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p7},            }, id="mp2  rhf    cd   ae: dd     ",),
        pytest.param({                                   "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p7},            }, id="mp2  uhf    cd   ae: dd     ",),
        pytest.param({                                   "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p7},            }, id="mp2 rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "dfmp2"},    "keywords": {"reference": "rhf",                                          "freeze_core": "true",                     },                               }, id="mp2  rhf         fc: dd     ",),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "uhf",                                          "freeze_core": "true",                     },                               }, id="mp2  uhf         fc: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }, "error": {1: _p6},            }, id="mp2 rohf         fc: dd     ",),
        pytest.param({"xptd": {"qc_module": "dfmp2"},    "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },                               }, id="mp2  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },                               }, id="mp2  uhf         ae: dd     ",),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }, "error": {1: _p6},            }, id="mp2 rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_mp2_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "gradient"))


#
#  ,--.   ,--.,------.  ,---.     ,--.  ,--.                     ,--.
#  |   `.'   ||  .--. ''.-.  \    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |'.'|  ||  '--' | .-' .'    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |   |  ||  | --' /   '-.    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'   `--'`--'     '-----'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  MP2 Hessian


@pytest.mark.parametrize("dertype", [
    pytest.param(2, id="hes2"),  # no analytic Hessians available
    pytest.param(1, id="hes1", marks=pytest.mark.findif),
    pytest.param(0, id="hes0", marks=[pytest.mark.long, pytest.mark.findif]),
])
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
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {2: _p24, 1: _p5}}, id="mp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {2: _p24, 1: _p5}}, id="mp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {2: _p24, 1: _p5}}, id="mp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {2: _p24}        }, id="mp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {2: _p24}        }, id="mp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    }, "error": {2: _p24, 1: _p6}}, id="mp2 rohf    conv ae: * occ  ",),
        # yapf: enable
    ],
)
def test_mp2_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2", "hessian"))


#
#  ,--.   ,--.,------.  ,---.     ,-----.    ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    '--. `\    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-..--..--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----''--'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                          `---' `---'
#  <<<  MP2.5 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                 }, id="mp2.5  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                 }, id="mp2.5  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                 }, id="mp2.5  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error":{0: _p1}}, id="mp2.5  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                 }, id="mp2.5  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                 }, id="mp2.5  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                 }, id="mp2.5  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                 }, id="mp2.5  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                 }, id="mp2.5  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },                 }, id="mp2.5  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                 }, id="mp2.5  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                 }, id="mp2.5  rhf   cd/conv rr occ  ",),
        # yapf: enable
    ],
)
def test_mp2p5_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "energy"))


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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "energy"))


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
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }},                     id="mp2.5  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }},                     id="mp2.5  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp2.5 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }},                     id="mp2.5  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }},                     id="mp2.5  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp2.5 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }},                     id="mp2.5  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }},                     id="mp2.5  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp2.5 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }},                     id="mp2.5  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }},                     id="mp2.5  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp2.5 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }},                     id="mp2.5  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }},                     id="mp2.5  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp2.5 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }},                     id="mp2.5  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }},                     id="mp2.5  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp2.5 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }},                     id="mp2.5  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }},                     id="mp2.5  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp2.5 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }},                     id="mp2.5  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }},                     id="mp2.5  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp2.5 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp2p5_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "energy"))


#
#  ,--.   ,--.,------.  ,---.     ,-----.     ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    |  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    '--. `\    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-..--..--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----''--'`----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP2.5 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p20},}, id="mp2.5  uhf    conv fc: * occ  ",),
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p20},        }, id="mp2.5  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p20},        }, id="mp2.5  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="mp2.5  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="mp2.5  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="mp2.5  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="mp2.5  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="mp2.5  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="mp2.5  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="mp2.5  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="mp2.5  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="mp2.5  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="mp2.5  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    },                            }, id="mp2.5  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    },                            }, id="mp2.5  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="mp2.5  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="mp2.5  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p32, 0: _p32}}, id="mp2.5 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp2p5_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp2.5", "gradient"))


#
#  ,--.   ,--.,------. ,----.     ,------.
#  |   `.'   ||  .--. ''.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP3 Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp3  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp3  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp3  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {0: _p1}}, id="mp3  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp3  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp3  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp3  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp3  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp3  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="mp3  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp3  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp3  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },                  }, id="mp3  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },                  }, id="mp3  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },                  }, id="mp3  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="mp3  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="mp3  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },                  }, id="mp3  rhf   cd/conv rr fnocc",),

        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "pk",     },                  }, id="mp3  rhf   pk/conv rr detci",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "direct", },                  }, id="mp3  rhf drct/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "df",     },                  }, id="mp3  rhf   df/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "mem_df", },                  }, id="mp3  rhf  mem/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "disk_df",},                  }, id="mp3  rhf disk/conv rr detci",),
        # pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",   "scf_type": "cd",     },                  }, id="mp3  rhf   cd/conv rr detci",),
        # yapf: enable
    ],
)
def test_mp3_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "energy"))


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
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp3 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp3 diff ans", raises=psi4.UpgradeHelper)),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp3  rhf    conv ae:   detci"),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp3 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp3 diff ans", raises=psi4.UpgradeHelper)),

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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "energy"))


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
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }                    }, id="mp3  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }                    }, id="mp3  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp3 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="mp3  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="mp3  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp3 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="mp3  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="mp3  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp3 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="mp3  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="mp3  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp3 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="mp3  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="mp3  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp3 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="mp3  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="mp3  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp3 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="mp3  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }                    }, id="mp3  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p32}}, id="mp3 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="mp3  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }                    }, id="mp3  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p32}}, id="mp3 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp3_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "energy"))


#
#  ,--.   ,--.,------. ,----.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' |  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP3 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p18},}, id="mp3  uhf    conv fc: * occ  ",),
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p18},        }, id="mp3  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p18},        }, id="mp3  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="mp3  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="mp3  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="mp3  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="mp3  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="mp3  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="mp3  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="mp3  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="mp3  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="mp3  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="mp3  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    },                            }, id="mp3  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    },                            }, id="mp3  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="mp3  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="mp3  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p32, 0: _p32}}, id="mp3 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp3_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp3", "gradient"))


#
#  ,--.   ,--.,------.   ,---.  ,-. ,---.  ,------.   ,-----.   ,-.      ,------.
#  |   `.'   ||  .--. ' /    | / .''   .-' |  .-.  \ '  .-.  '  '. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  ||  | `.  `-. |  |  \  :|  | |  |   |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  ||  | .-'    ||  '--'  /'  '-'  '-. |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--' \ '.`-----' `-------'  `-----'--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                               `-'                             `-'                                    `---' `---'
#  <<<  MP4(SDQ) Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
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

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp4_sdq_  rhf    conv fc:   fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp4_sdq_  rhf    conv ae:   fnocc",),
        # yapf: enable
    ],
)
def test_mp4_prsdq_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp4(sdq)", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     },                   }, id="mp4_sdq_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                   }, id="mp4_sdq_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_ rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="mp4_sdq_  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_ rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="mp4_sdq_  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_ rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="mp4_sdq_  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_ rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="mp4_sdq_  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_ rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="mp4_sdq_  rhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_  uhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="mp4_sdq_ rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="mp4_sdq_  rhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_  uhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="mp4_sdq_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp4_prsdq_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp4(sdq)", "energy"))


#
#  ,--.   ,--.,------.   ,---.    ,------.
#  |   `.'   ||  .--. ' /    |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP4 Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="mp4  rhf    conv fc:   fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="mp4  rhf    conv ae:   fnocc",),

        ###### detci
        # * detci rohf mp4 does not match other programs (cfour) in the stored reference
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp4  rhf    conv fc:   detci",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="mp4 rohf    conv fc:   detci", marks=pytest.mark.xfail(reason="detci rohf mp4 diff ans", raises=psi4.UpgradeHelper)),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp4  rhf    conv ae:   detci"),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="mp4 rohf    conv ae:   detci", marks=pytest.mark.xfail(reason="detci rohf mp4 diff ans", raises=psi4.UpgradeHelper)),
        # yapf: enable
    ],
)
def test_mp4_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp4", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     },                   }, id="mp4  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p42}}, id="mp4  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p42}}, id="mp4 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                   }, id="mp4  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p42}}, id="mp4  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p42}}, id="mp4 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4 rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p41}}, id="mp4 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p41}}, id="mp4 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="mp4  rhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p42}}, id="mp4  uhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p42}}, id="mp4 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="mp4  rhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p42}}, id="mp4  uhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p42}}, id="mp4 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_mp4_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "mp4", "energy"))


#
#  ,-------.  ,---.  ,------. ,--------. ,---.     ,------.
#  `--.   /  /  O  \ |  .--. ''--.  .--''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#    /   /  |  .-.  ||  '--' |   |  |    .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#   /   `--.|  | |  ||  | --'    |  |   /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-------'`--' `--'`--'        `--'   '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                `---' `---'
#  <<<  ZAPT2 Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
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
        # RHF DE-DEFINED FOR ZAPT pytest.param({"xptd": {}, "keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "true"}                     }, id="zapt2  rhf   conv fc:   detci"),
        # RHF DE-DEFINED FOR ZAPT pytest.param({"xptd": {}, "keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "detci", "freeze_core": "false"},                   }, id="zapt2  rhf   conv ae:   detci"),
        pytest.param({"xptd": {}, "keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "true"},                    }, id="zapt2 rohf   conv fc:   detci"),
        pytest.param({"xptd": {}, "keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "detci", "freeze_core": "false"},                   }, id="zapt2 rohf   conv ae:   detci"),
        # yapf: enable
    ],
)
def test_zapt2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "zapt2", "energy"))


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
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p47}}, id="zapt2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p45}}, id="zapt2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     },                   }, id="zapt2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p47}}, id="zapt2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p45}}, id="zapt2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    },                   }, id="zapt2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p46}}, id="zapt2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p45}}, id="zapt2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p46}}, id="zapt2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p46}}, id="zapt2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p45}}, id="zapt2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p46}}, id="zapt2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p46}}, id="zapt2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p45}}, id="zapt2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p46}}, id="zapt2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p46}}, id="zapt2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p45}}, id="zapt2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p46}}, id="zapt2 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }, "error": {0: _p47}}, id="zapt2  rhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p45}}, id="zapt2  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     },                   }, id="zapt2 rohf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }, "error": {0: _p47}}, id="zapt2  rhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p45}}, id="zapt2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    },                   }, id="zapt2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_zapt2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "zapt2", "energy"))


#
#   ,-----.,--. ,---.  ,------.      ,------.
#  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                  `---' `---'
#  <<<  CISD Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="cisd  rhf    conv fc:   fnocc"),
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="cisd  rhf    conv ae:   fnocc"),

        ###### detci
        # * detci rohf cisd does not match other programs (cfour) in the stored reference; that is, vcc = tce != guga = detci. detci vals in ref_local so block passes
        pytest.param({                        "keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="cisd  rhf    conv fc:   detci"),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv", "qc_module": "detci", "freeze_core": "true",                   },}, id="cisd rohf    conv fc:   detci"),
        pytest.param({                        "keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="cisd  rhf    conv ae:   detci"),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv", "qc_module": "detci", "freeze_core": "false",                  },}, id="cisd rohf    conv ae:   detci"),
        # yapf: enable
    ],
)
def test_cisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cisd", "energy"))


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, ci_type) work?

        ###### default qc_module
        pytest.param({"xptd": {"qc_module": "fnocc"},               "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "true",                     }                    }, id="cisd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p48}}, id="cisd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "true",                     }                    }, id="cisd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"},               "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "false",                    }                    }, id="cisd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p48}}, id="cisd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "false",                    }                    }, id="cisd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="cisd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="cisd rohf    cd   ae: dd     "),

        ###### default qc_module, ci_type
        pytest.param({"xptd": {"qc_module": "fnocc"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="cisd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p48}}, id="cisd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }                    }, id="cisd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="cisd  rhf         ae: dd     "),
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p48}}, id="cisd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }                    }, id="cisd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cisd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cisd", "energy"))


#
#   ,-----.    ,-----.,--. ,---.  ,------.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  QCISD Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="qcisd  rhf    conv fc:   fnocc"),
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="qcisd  rhf    conv ae:   fnocc"),
        # yapf: enable
    ],
)
def test_qcisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "qcisd", "energy"))


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, ci_type) work?

        ###### default qc_module
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "true",                     },                   }, id="qcisd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "false",                    },                   }, id="qcisd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd rohf    cd   ae: dd     "),

        ###### default qc_module, ci_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="qcisd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="qcisd  rhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_qcisd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "qcisd", "energy"))


#
#   ,-----.    ,-----.,--. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                             `-'          `-'                                    `---' `---'
#  <<<  QCISD(T) Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="qcisd_t_  rhf    conv fc:   fnocc"),
        pytest.param({"keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="qcisd_t_  rhf    conv ae:   fnocc"),
        # yapf: enable
    ],
)
def test_qcisd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "qcisd(t)", "energy"))


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.quick),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, ci_type) work?

        ###### default qc_module
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "true",                     },                   }, id="qcisd_t_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd_t_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd_t_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "false",                    },                   }, id="qcisd_t_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd_t_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd_t_ rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_ rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_ rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="qcisd_t_ rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="qcisd_t_ rohf    cd   ae: dd     "),

        ###### default qc_module, ci_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="qcisd_t_  rhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd_t_  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="qcisd_t_ rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="qcisd_t_  rhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd_t_  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="qcisd_t_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_qcisd_prt_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "qcisd(t)", "energy"))


#
#  ,------. ,-----.,--.    ,------.
#  |  .---''  .--./|  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  `--, |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |`   '  '--'\|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `-----'`--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                        `---' `---'
#  <<<  FCI Energy


# Note: commented lines in fci_energy_module and fci_energy_default run just fine. But the uncommented
#   lines can generate the right table symbols in one run rather than running the expensive ae three
#   times. Former function only runs dz anyways, so nothing lost.

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ROHF FC: gms != detci by 1e-5. detci value stored in psi4 ref_local file
        # pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "detci", "freeze_core": "true",                                   }}, id="fci  rhf    conv fc:   detci"),
        # pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv", "qc_module": "detci", "freeze_core": "true",                                   }}, id="fci rohf    conv fc:   detci"),
        # pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rhf",  "ci_type": "conv", "qc_module": "detci", "freeze_core": "false",                                  }}, id="fci  rhf    conv ae:   detci"),
        # pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv", "qc_module": "detci", "freeze_core": "false",                                  }}, id="fci rohf    conv ae:   detci"),
        # yapf: enable
    ],
)
def test_fci_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "fci", "energy"))


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz", marks=pytest.mark.nonroutine),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Does the simple interface (default qc_module, scf_type, ci_type) work?

        ###### default qc_module
        # pytest.param({"xptd": {"qc_module": "detci"},               "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "true",                     }                    }, id="fci  rhf    conv fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p45}}, id="fci  uhf    conv fc: dd     "),
        # pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "true",                     }                    }, id="fci rohf    conv fc: dd     "),
        # pytest.param({"xptd": {"qc_module": "detci"},               "keywords": {"reference": "rhf",  "ci_type": "conv",                     "freeze_core": "false",                    }                    }, id="fci  rhf    conv ae: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p45}}, id="fci  uhf    conv ae: dd     "),
        # pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf", "ci_type": "conv",                     "freeze_core": "false",                    }                    }, id="fci rohf    conv ae: dd     "),
        # ####
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci  rhf    df   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci  uhf    df   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci rohf    df   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci  rhf    df   ae: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci  uhf    df   ae: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci rohf    df   ae: dd     "),
        # ####
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci  rhf    cd   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci  uhf    cd   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p43}}, id="fci rohf    cd   fc: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci  rhf    cd   ae: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",  "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci  uhf    cd   ae: dd     "),
        # pytest.param({"xptd": {},                                   "keywords": {"reference": "rohf", "ci_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p43}}, id="fci rohf    cd   ae: dd     "),

        ###### default qc_module, ci_type
        pytest.param({"xptd": {"qc_module": "detci"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="fci  rhf         fc: dd     "),  # seconds
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p45}}, id="fci  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }                    }, id="fci rohf         fc: dd     "),  # seconds
        pytest.param({"xptd": {"qc_module": "detci"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="fci  rhf         ae: dd     "),  # 6m
        pytest.param({"xptd": {},                                   "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p45}}, id="fci  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "detci", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }                    }, id="fci rohf         ae: dd     "),  # 30m
        # yapf: enable
    ],
)
def test_fci_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "fci", "energy"))


#
#  ,------. ,------.,--.   ,--.,------.  ,---.     ,------.
#  |  .--. '|  .---'|   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  '--'.'|  `--, |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |\  \ |  `---.|  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--' '--'`------'`--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                `---' `---'
#  <<<  REMP2 Energy

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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  rhf    conv fc:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  rhf    conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  uhf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  uhf    df   ae:   dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="remp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="remp2  uhf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_remp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "remp2", "energy"))


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
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="remp2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="remp2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="remp2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="remp2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="remp2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="remp2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "df",                     "freeze_core": "true",                     }                    }, id="remp2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "df",                     "freeze_core": "true",                     }                    }, id="remp2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "df",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="remp2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "df",                     "freeze_core": "false",                    }                    }, id="remp2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "df",                     "freeze_core": "false",                    }                    }, id="remp2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "df",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="remp2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "cd",                     "freeze_core": "true",                     }                    }, id="remp2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "cd",                     "freeze_core": "true",                     }                    }, id="remp2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "cd",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="remp2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "cd",                     "freeze_core": "false",                    }                    }, id="remp2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "cd",                     "freeze_core": "false",                    }                    }, id="remp2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "cd",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="remp2 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="remp2  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }                    }, id="remp2  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p32}}, id="remp2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="remp2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }                    }, id="remp2  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p32}}, id="remp2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_remp2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "remp2", "energy"))


#
#  ,--.    ,-----. ,-----.,------.      ,------.
#  |  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                     `---' `---'
#  <<<  LCCD Energy


@pytest.mark.parametrize("dertype", [pytest.param(0, id="ene0"),])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are scf_types managed properly by proc.py? Generally skip corl_type=cd, so df & conv only.

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="lccd  rhf   pk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="lccd  rhf drct/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="lccd  rhf   df/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", }, "error": {0: _p1}}, id="lccd  rhf  mem/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="lccd  rhf disk/df   rr dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="lccd  rhf   cd/df   rr dfocc",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "pk",     },                  }, id="lccd  rhf   pk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "direct", },                  }, id="lccd  rhf drct/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "df",     },                  }, id="lccd  rhf   df/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="lccd  rhf  mem/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="lccd  rhf disk/conv rr occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",   "freeze_core": "false",  "scf_type": "cd",     },                  }, id="lccd  rhf   cd/conv rr occ  ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "pk",     },                  }, id="lccd  rhf   pk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "direct", },                  }, id="lccd  rhf drct/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "df",     },                  }, id="lccd  rhf   df/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "mem_df", },                  }, id="lccd  rhf  mem/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "disk_df",},                  }, id="lccd  rhf disk/conv rr fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",  "scf_type": "cd",     },                  }, id="lccd  rhf   cd/conv rr fnocc",),
        # yapf: enable
    ],
)
def test_lccd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "energy"))


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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="lccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="lccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p32}}, id="lccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="lccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="lccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p32}}, id="lccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="lccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="lccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="lccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="lccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="lccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="lccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="lccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="lccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p32}}, id="lccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="lccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="lccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p32}}, id="lccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="lccd  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }                    }, id="lccd  uhf         fc: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p32}}, id="lccd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="lccd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},   "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }                    }, id="lccd  uhf         ae: dd     "),
        pytest.param({"xptd": {},                     "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p32}}, id="lccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_lccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "energy"))


#
#  ,--.    ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  |  |   '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |   |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--.'  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `-----' `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  LCCD Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _p22},}, id="lccd  uhf    conv fc: * occ  ",),
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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p22},        }, id="lccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p22},        }, id="lccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="lccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="lccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="lccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="lccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="lccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="lccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="lccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33},        }, id="lccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="lccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33},        }, id="lccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p32, 0: _p32}}, id="lccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p22},        }, id="lccd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p22},        }, id="lccd  uhf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p32, 0: _p32}}, id="lccd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": {1: "occ"}}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="lccd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="lccd  uhf         ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p32, 0: _p32}}, id="lccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_lccd_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccd", "gradient"))


#
#  ,--.    ,-----. ,-----. ,---.  ,------.      ,------.
#  |  |   '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  LCCSD Energy aka CEPA(0)


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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccsd", "energy"))


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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccsd", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="lccsd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="lccsd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd rohf    conv ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="lccsd  rhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd  uhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd rohf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="lccsd  rhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd  uhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd rohf    df   ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="lccsd  rhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd  uhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd rohf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="lccsd  rhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd  uhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="lccsd  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="lccsd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="lccsd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="lccsd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_lccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "lccsd", "energy"))


#
#   ,-----.,------.,------.   ,---.    ,-. ,--.,-.      ,------.
#  '  .--./|  .---'|  .--. ' /  O  \  / .'/   |'. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  `--, |  '--' ||  .-.  ||  | `|  | |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\|  `---.|  | --' |  | |  ||  |  |  | |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'`------'`--'     `--' `--' \ '. `--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'     `-'                                    `---' `---'
#  <<<  CEPA(1) Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="cepa_1_  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="cepa_1_  rhf    conv ae: * fnocc",),
        # yapf: enable
    ],
)
def test_cepa_pr1_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cepa(1)", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="cepa_1_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="cepa_1_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_ rohf    conv ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="cepa_1_  rhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_  uhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_ rohf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="cepa_1_  rhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_  uhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_ rohf    df   ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="cepa_1_  rhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_  uhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_ rohf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="cepa_1_  rhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_  uhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="cepa_1_  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_1_ rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="cepa_1_  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_1_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cepa_pr1_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cepa(1)", "energy"))


#
#   ,-----.,------.,------.   ,---.    ,-.,----. ,-.      ,------.
#  '  .--./|  .---'|  .--. ' /  O  \  / .''.-.  |'. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  `--, |  '--' ||  .-.  ||  |   .' <  |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\|  `---.|  | --' |  | |  ||  | /'-'  | |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'`------'`--'     `--' `--' \ '.`----' .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'       `-'                                    `---' `---'
#  <<<  CEPA(3) Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="cepa_3_  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="cepa_3_  rhf    conv ae: * fnocc",),
        # yapf: enable
    ],
)
def test_cepa_pr3_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cepa(3)", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="cepa_3_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="cepa_3_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_ rohf    conv ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="cepa_3_  rhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_  uhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_ rohf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="cepa_3_  rhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_  uhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_ rohf    df   ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="cepa_3_  rhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_  uhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_ rohf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="cepa_3_  rhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_  uhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="cepa_3_  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="cepa_3_ rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="cepa_3_  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="cepa_3_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cepa_pr3_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cepa(3)", "energy"))


#
#    ,---.   ,-----.,------. ,------.    ,------.
#   /  O  \ '  .--./|  .--. '|  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  ||  |    |  '--' ||  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  | |  |'  '--'\|  | --' |  |`       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--' `--' `-----'`--'     `--'        `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                      `---' `---'
#  <<<  ACPF Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="acpf  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="acpf  rhf    conv ae: * fnocc",),
        # yapf: enable
    ],
)
def test_acpf_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "acpf", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="acpf  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="acpf  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf rohf    conv ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="acpf  rhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf  uhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf rohf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="acpf  rhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf  uhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf rohf    df   ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="acpf  rhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf  uhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf rohf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="acpf  rhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf  uhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="acpf  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="acpf rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="acpf  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="acpf rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_acpf_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "acpf", "energy"))


#
#    ,---.   ,-----.    ,-----. ,-----.    ,------.
#   /  O  \ '  .-.  '  '  .--./'  .--./    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  ||  | |  |  |  |    |  |        |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  | |  |'  '-'  '-.'  '--'\'  '--'\    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--' `--' `-----'--' `-----' `-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                        `---' `---'
#  <<<  AQCC Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "true",                   },}, id="aqcc  rhf    conv fc: * fnocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc", "freeze_core": "false",                  },}, id="aqcc  rhf    conv ae: * fnocc",),
        # yapf: enable
    ],
)
def test_aqcc_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "aqcc", "energy"))


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
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="aqcc  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="aqcc  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc rohf    conv ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="aqcc  rhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc  uhf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc rohf    df   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="aqcc  rhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc  uhf    df   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc rohf    df   ae: dd     "),
        ####
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p44}}, id="aqcc  rhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc  uhf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc rohf    cd   fc: dd     "),
        pytest.param({                                "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p44}}, id="aqcc  rhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc  uhf    cd   ae: dd     "),
        pytest.param({                                "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="aqcc  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p40}}, id="aqcc rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="aqcc  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p40}}, id="aqcc rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_aqcc_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "aqcc", "energy"))


#
#   ,-----. ,-----.,------.      ,------.
#  '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  CCD Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true", }, "error": {0: _p26},}, id="ccd  rhf    conv fc:   occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="ccd  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },}, id="ccd  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="ccd  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },}, id="ccd  uhf    df   ae:   dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="ccd  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },}, id="ccd  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="ccd  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },}, id="ccd  uhf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_ccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccd", "energy"))


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
        pytest.param({                              "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd  rhf    conv fc: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd  uhf    conv fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd rohf    conv fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd  rhf    conv ae: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd  uhf    conv ae: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                   }, id="ccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                   }, id="ccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                   }, id="ccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                   }, id="ccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                   }, id="ccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                   }, id="ccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                   }, id="ccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                   }, id="ccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({                              "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd  rhf         fc: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd  uhf         fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p39}}, id="ccd rohf         fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd  rhf         ae: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd  uhf         ae: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p39}}, id="ccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccd", "energy"))


#
#   ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCD Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="ccd  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                    }, id="ccd  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="ccd  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                    }, id="ccd  uhf    df   ae:   dfocc",),
        # yapf: enable
    ],
)
def test_ccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccd", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p39, 0: _p39}}, id="ccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p39, 0: _p39}}, id="ccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p39, 0: _p39}}, id="ccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p39, 0: _p39}}, id="ccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p39, 0: _p39}}, id="ccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p39, 0: _p39}}, id="ccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="ccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="ccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="ccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="ccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="ccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="ccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="ccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="ccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p39, 0: _p39}}, id="ccd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p39, 0: _p39}}, id="ccd  uhf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p39, 0: _p39}}, id="ccd rohf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }, "error": {1: _p39, 0: _p39}}, id="ccd  rhf         ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }, "error": {1: _p39, 0: _p39}}, id="ccd  uhf         ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p39, 0: _p39}}, id="ccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccd_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccd", "gradient"))


#
#  ,-----.   ,-----. ,-----.,------.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                       `---' `---'
#  <<<  BCCD Energy


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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd  uhf  conv fc   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd rohf  conv fc sc: * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd  rhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd  uhf  conv ae   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd rohf  conv ae sc: * ccenergy",),
        # yapf: enable
    ],
)
def test_bccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "bccd", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",     }                      }, id="bccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",     }                      }, id="bccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",     }                      }, id="bccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",    }                      }, id="bccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",    }                      }, id="bccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",    }                      }, id="bccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",     }                      }, id="bccd  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "true",     }                      }, id="bccd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf",                                        "freeze_core": "true",     }                      }, id="bccd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",    }                      }, id="bccd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "false",    }                      }, id="bccd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf",                                        "freeze_core": "false",    }                      }, id="bccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_bccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "bccd", "energy"))


#
#   ,-----. ,-----. ,---.     ,------.
#  '  .--./'  .--./'.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |     .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----''-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC2 Energy


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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc2  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc2  uhf  conv fc   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc2 rohf  conv fc sc: * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc2  rhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc2  uhf  conv ae   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc2 rohf  conv ae sc: * ccenergy",),
        # yapf: enable
    ],
)
def test_cc2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc2", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc2 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",     }                    }, id="cc2  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "true",     }                    }, id="cc2  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf",                                        "freeze_core": "true",     }                    }, id="cc2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",    }                    }, id="cc2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "false",    }                    }, id="cc2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf",                                        "freeze_core": "false",    }                    }, id="cc2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cc2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc2", "energy"))


#
#   ,-----. ,-----. ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |     .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\/   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----''-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CC2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                        }, "error": {1: _p16},}, id="cc2  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                        }, "error": {1: _p27},}, id="cc2  uhf  conv fc   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                        }, "error": {1: _p27},}, id="cc2 rohf  conv fc sc: * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                       },                    }, id="cc2  rhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                       }, "error": {1: _p27},}, id="cc2  uhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                       }, "error": {1: _p27},}, id="cc2 rohf  conv ae   : * ccenergy",),
        # yapf: enable
    ],
)
def test_cc2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc2", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable

        ###### default qc_module
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p16},        }, id="cc2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p27},        }, id="cc2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p27},        }, id="cc2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="cc2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p27},        }, id="cc2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p27},        }, id="cc2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p37, 0: _p37}}, id="cc2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p37, 0: _p37}}, id="cc2 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p16},        }, id="cc2  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p27},        }, id="cc2  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p27},        }, id="cc2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="cc2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }, "error": {1: _p27},        }, id="cc2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p27},        }, id="cc2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cc2_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc2", "gradient"))


#
#   ,-----. ,-----. ,---.  ,------.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                      `---' `---'
#  <<<  CCSD Energy


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

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     },                    }, id="ccsd  rhf   pk/conv rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", },                    }, id="ccsd  rhf drct/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },                    }, id="ccsd  rhf   df/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", },                    }, id="ccsd  rhf  mem/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},                    }, id="ccsd  rhf disk/conv rr fnocc   ",),
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },                    }, id="ccsd  rhf   cd/conv rr fnocc   ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",     }, "error": {0: _p15},}, id="ccsd  rhf   pk/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "direct", }, "error": {0: _p15},}, id="ccsd  rhf drct/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "df",     },                    }, id="ccsd  rhf   df/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "mem_df", }, "error": {0: _p1}, }, id="ccsd  rhf  mem/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "disk_df",},                    }, id="ccsd  rhf disk/df   rr fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",     },                    }, id="ccsd  rhf   cd/df   rr fnocc   ",),

        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",     },                    }, id="ccsd  rhf   pk/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "direct", },                    }, id="ccsd  rhf drct/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "df",     },                    }, id="ccsd  rhf   df/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "mem_df", }, "error": {0: _p1}, }, id="ccsd  rhf  mem/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "disk_df",},                    }, id="ccsd  rhf disk/df   rr dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",     },                    }, id="ccsd  rhf   cd/df   rr dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_energy_scftype(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd", "energy"))


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
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="ccsd  rhf  conv fc   : * ccenergy",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="ccsd  uhf  conv fc   : * ccenergy",),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="ccsd rohf  conv fc sd: * ccenergy",),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",  "semicanonical": "true"},                    }, id="ccsd rohf  conv fc sc:   ccenergy",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="ccsd  rhf  conv ae   : * ccenergy",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="ccsd  uhf  conv ae   : * ccenergy",),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="ccsd rohf  conv ae sd: * ccenergy",),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "semicanonical": "true"},                    }, id="ccsd rohf  conv ae sc:   ccenergy",),

        # below run but see #2710
        # pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "df",},                    }, id="ccsd  rhf    df ae   :   ccenergy",),
        # pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "ccenergy", "freeze_core": "false", "scf_type": "df",},                    }, id="ccsd rohf    df ae   :   ccenergy",),

        ###### mrcc
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                         },                    }, id="ccsd  rhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                         },                    }, id="ccsd  uhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                         },                    }, id="ccsd rohf  conv fc sc:   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                        },                    }, id="ccsd  rhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                        },                    }, id="ccsd  uhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                        },                    }, id="ccsd rohf  conv ae sc:   mrcc    ", marks=using("mrcc")),

        ###### fnocc
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "true",                         },                    }, id="ccsd  rhf  conv fc   :   fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false",                        },                    }, id="ccsd  rhf  conv ae   :   fnocc   ",),
        ####
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",                         },                    }, id="ccsd  rhf    df fc   : * fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false",                        },                    }, id="ccsd  rhf    df ae   : * fnocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",      }, "error": {0: _p15},}, id="ccsd  rhf pk/df fc   :   fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",      }, "error": {0: _p15},}, id="ccsd  rhf pk/df ae   :   fnocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "cd",      },                    }, id="ccsd  rhf cd/df fc   :   fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",      },                    }, id="ccsd  rhf cd/df ae   :   fnocc   ",),
        ####
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",                         },                    }, id="ccsd  rhf    cd fc   : * fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false",                        },                    }, id="ccsd  rhf    cd ae   : * fnocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",      }, "error": {0: _p15},}, id="ccsd  rhf pk/cd fc   :   fnocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",      }, "error": {0: _p15},}, id="ccsd  rhf pk/cd ae   :   fnocc   ",),

        ###### dfocc
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                         },                    }, id="ccsd  rhf    df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                         },                    }, id="ccsd  uhf    df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                        },                    }, id="ccsd  rhf    df ae   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                        },                    }, id="ccsd  uhf    df ae   :   dfocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",      },                    }, id="ccsd  rhf pk/df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",      },                    }, id="ccsd  uhf pk/df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",      },                    }, id="ccsd  rhf pk/df ae   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",      },                    }, id="ccsd  uhf pk/df ae   :   dfocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",      },                    }, id="ccsd  rhf cd/df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",      },                    }, id="ccsd  uhf cd/df fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",      },                    }, id="ccsd  rhf cd/df ae   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",      },                    }, id="ccsd  uhf cd/df ae   :   dfocc   ",),
        ####
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                         },                    }, id="ccsd  rhf    cd fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                         },                    }, id="ccsd  uhf    cd fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                        },                    }, id="ccsd  rhf    cd ae   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                        },                    }, id="ccsd  uhf    cd ae   :   dfocc   ",),
        ##
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",      },                    }, id="ccsd  rhf pk/cd fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",      },                    }, id="ccsd  uhf pk/cd fc   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",      },                    }, id="ccsd  rhf pk/cd ae   :   dfocc   ",),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",      },                    }, id="ccsd  uhf pk/cd ae   :   dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="ccsd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="ccsd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }                    }, id="ccsd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="ccsd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="ccsd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="ccsd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "fnocc"},                  "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="ccsd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                   }, id="ccsd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccsd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"},                  "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                   }, id="ccsd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                   }, id="ccsd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccsd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "fnocc"},                  "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                   }, id="ccsd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                   }, id="ccsd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccsd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "fnocc"},                  "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                   }, id="ccsd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                   }, id="ccsd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccsd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }                    }, id="ccsd  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }                    }, id="ccsd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }                    }, id="ccsd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="ccsd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }                    }, id="ccsd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }                    }, id="ccsd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd", "energy"))


#
#   ,-----. ,-----. ,---.  ,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSD Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         }, "error": {1: _p16}}, id="ccsd  rhf  conv fc   : * ccenergy"),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         }, "error": {1: _p16}}, id="ccsd  uhf  conv fc   : * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         }, "error": {1: _p16}}, id="ccsd rohf  conv fc sd: * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",  "semicanonical": "true"}, "error": {1: _p16}}, id="ccsd rohf  conv fc sc: * ccenergy"),  # SEMI-DISABLE
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                   }, id="ccsd  rhf  conv ae   : * ccenergy"),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                   }, id="ccsd  uhf  conv ae   : * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                   }, id="ccsd rohf  conv ae sd: * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false", "semicanonical": "true"}, "wrong": {1: _w1 }}, id="ccsd rohf  conv ae sc: * ccenergy"),  # SEMI-DISABLE

        ###### dfocc
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                                           }, id="ccsd  rhf    df fc   : * dfocc   "),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                                           }, id="ccsd  uhf    df fc   : * dfocc   "),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                                           }, id="ccsd  rhf    df ae   : * dfocc   "),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                                           }, id="ccsd  uhf    df ae   : * dfocc   "),
        # yapf: enable
    ],
)
def test_ccsd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable

        ###### default qc_module
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p16},        }, id="ccsd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p16},        }, id="ccsd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p16},        }, id="ccsd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="ccsd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="ccsd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="ccsd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": {1: "occ", 0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="ccsd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="ccsd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccsd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {1: "occ", 0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="ccsd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="ccsd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccsd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": {0: "fnocc"}},             "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="ccsd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {0: "occ"}},               "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="ccsd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccsd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {0: "fnocc"}},             "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="ccsd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": {0: "occ"}},               "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="ccsd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccsd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p16},        }, id="ccsd  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p16},        }, id="ccsd  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p16},        }, id="ccsd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="ccsd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="ccsd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy", "sdsc": "sd"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   },                            }, id="ccsd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccsd_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd", "gradient"))


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'          `-'                                    `---' `---'
#  <<<  CCSD(T) Energy


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
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf  conv fc   : * ccenergy"),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },                    }, id="ccsd_t_  uhf  conv fc   : * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },                    }, id="ccsd_t_ rohf  conv fc sc: * ccenergy"),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf  conv ae   : * ccenergy"),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },                    }, id="ccsd_t_  uhf  conv ae   : * ccenergy"),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },                    }, id="ccsd_t_ rohf  conv ae sc: * ccenergy"),

        ###### mrcc
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="ccsd_t_  uhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="ccsd_t_ rohf  conv fc sc:   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({                        "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="ccsd_t_  uhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"xptd": {"sdsc": "sc"}, "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="ccsd_t_ rohf  conv ae sc:   mrcc    ", marks=using("mrcc")),

        ###### fnocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf  conv fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "fnocc",    "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf  conv ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf    df fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf    df ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": {0: _p15}}, id="ccsd_t_  rhf pk/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": {0: _p15}}, id="ccsd_t_  rhf pk/df ae:   fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "cd",},                    }, id="ccsd_t_  rhf cd/df fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "cd",},                    }, id="ccsd_t_  rhf cd/df ae:   fnocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf    cd fc: * fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf    cd ae: * fnocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "true",  "scf_type": "pk",}, "error": {0: _p15}}, id="ccsd_t_  rhf pk/cd fc:   fnocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "fnocc",    "freeze_core": "false", "scf_type": "pk",}, "error": {0: _p15}}, id="ccsd_t_  rhf pk/cd ae:   fnocc   ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf    df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="ccsd_t_  uhf    df fc:   dfocc   ",),  # SEMI-DISABLE
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf    df ae:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="ccsd_t_  uhf    df ae:   dfocc   ",),  # SEMI-DISABLE
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="ccsd_t_  rhf pk/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="ccsd_t_  uhf pk/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="ccsd_t_  rhf pk/df ae:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="ccsd_t_  uhf pk/df ae:   dfocc   ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},                    }, id="ccsd_t_  rhf cd/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},                    }, id="ccsd_t_  uhf cd/df fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},                    }, id="ccsd_t_  rhf cd/df ae:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},                    }, id="ccsd_t_  uhf cd/df ae:   dfocc   ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="ccsd_t_  rhf    cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="ccsd_t_  uhf    cd fc:   dfocc   ",),  # SEMI-DISABLE
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="ccsd_t_  rhf    cd ae:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="ccsd_t_  uhf    cd ae:   dfocc   ",),  # SEMI-DISABLE
        ##
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="ccsd_t_  rhf pk/cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="ccsd_t_  uhf pk/cd fc:   dfocc   ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="ccsd_t_  rhf pk/cd ae:   dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="ccsd_t_  uhf pk/cd ae:   dfocc   ",),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd(t)", "energy"))


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
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="ccsd_t_  rhf    conv fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="ccsd_t_  uhf    conv fc: dd     "),
        pytest.param({"sdsc": "sc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="ccsd_t_ rohf    conv fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="ccsd_t_  rhf    conv ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="ccsd_t_  uhf    conv ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="ccsd_t_ rohf    conv ae: dd     "),
        ####
        pytest.param({              "xptd": {"qc_module": "fnocc"},    "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                   }, id="ccsd_t_  rhf    df   fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p99}}, id="ccsd_t_  uhf    df   fc: dd     "),
        pytest.param({              "xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccsd_t_ rohf    df   fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "fnocc"},    "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                   }, id="ccsd_t_  rhf    df   ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p99}}, id="ccsd_t_  uhf    df   ae: dd     "),
        pytest.param({              "xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccsd_t_ rohf    df   ae: dd     "),
        ####
        pytest.param({              "xptd": {"qc_module": "fnocc"},    "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     },                   }, id="ccsd_t_  rhf    cd   fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p99}}, id="ccsd_t_  uhf    cd   fc: dd     "),
        pytest.param({              "xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="ccsd_t_ rohf    cd   fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "fnocc"},    "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    },                   }, id="ccsd_t_  rhf    cd   ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p99}}, id="ccsd_t_  uhf    cd   ae: dd     "),
        pytest.param({              "xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="ccsd_t_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="ccsd_t_  rhf         fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     },                   }, id="ccsd_t_  uhf         fc: dd     "),
        pytest.param({"sdsc": "sc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     },                   }, id="ccsd_t_ rohf         fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="ccsd_t_  rhf         ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    },                   }, id="ccsd_t_  uhf         ae: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    },                   }, id="ccsd_t_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd(t)", "energy"))


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                      `-'          `-'
#  <<<  CCSD(T) Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true", }, "error": {1: _p16},}, id="ccsd_t_  rhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true", }, "error": {1: _p16},}, id="ccsd_t_  uhf  conv fc: * ccenergy",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",},                    }, id="ccsd_t_  rhf  conv ae: * ccenergy",),  # SEMI-DISABLE
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",},                    }, id="ccsd_t_  uhf  conv ae: * ccenergy",),  # SEMI-DISABLE

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                    }, id="ccsd_t_  rhf    df fc: * dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                    }, id="ccsd_t_  uhf    df fc: * dfocc   ",),  # SEMI-DISABLE
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                    }, id="ccsd_t_  rhf    df ae: * dfocc   ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                    }, id="ccsd_t_  uhf    df ae: * dfocc   ",),  # SEMI-DISABLE
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd(t)", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable

        ###### default qc_module
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p98},        }, id="ccsd_t_  rhf    conv fc: dd     "),  # FORMERLY  {1: _p16}
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p98},        }, id="ccsd_t_  uhf    conv fc: dd     "),  # FORMERLY  {1: _p16}
        pytest.param({"sdsc": "sc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p35},        }, id="ccsd_t_ rohf    conv fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p98},        }, id="ccsd_t_  rhf    conv ae: dd     "),  # FORMERLY {}
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p98},        }, id="ccsd_t_  uhf    conv ae: dd     "),  # FORMERLY {}
        pytest.param({"sdsc": "sc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {1: _p35},        }, id="ccsd_t_ rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": {1: "occ", 0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="ccsd_t_  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p99, 0: _p99}}, id="ccsd_t_  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccsd_t_ rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {1: "occ", 0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="ccsd_t_  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},                    "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p99, 0: _p99}}, id="ccsd_t_  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccsd_t_ rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": {          0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="ccsd_t_  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {          0: "occ"}},     "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p99, 0: _p99}}, id="ccsd_t_  uhf    cd   fc: dd     "),  # SOON {1: _p33}
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p34, 0: _p34}}, id="ccsd_t_ rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": {          0: "fnocc"}},   "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="ccsd_t_  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": {          0: "occ"}},     "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p99, 0: _p99}}, id="ccsd_t_  uhf    cd   ae: dd     "),  # SOON {1: _p33}
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p34, 0: _p34}}, id="ccsd_t_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p98},        }, id="ccsd_t_  rhf         fc: dd     "),  # FORMERLY {1: _p16}
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p98},        }, id="ccsd_t_  uhf         fc: dd     "),  # FORMERLY {1: _p16}
        pytest.param({"sdsc": "sc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p35},        }, id="ccsd_t_ rohf         fc: dd     "),
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }, "error": {1: _p98},        }, id="ccsd_t_  rhf         ae: dd     "),  # FORMERLY {}
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }, "error": {1: _p98},        }, id="ccsd_t_  uhf         ae: dd     "),  # FORMERLY {}
        pytest.param({              "xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }, "error": {1: _p35},        }, id="ccsd_t_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "ccsd(t)", "gradient"))


#
#                  ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#   ,--,--.,-----.'  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  ' ,-.  |'-----'|  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  \ '-'  |       '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `--`--'        `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                     `-'          `-'                                    `---' `---'
#  <<<  a-CCSD(T) Energy aka Lambda-CCSD(T) aka CCSD(aT) aka CCSD(T)_L


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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                   },                    }, id="a-ccsd_t_  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                  },                    }, id="a-ccsd_t_  rhf  conv ae   : * ccenergy",),

        ###### mrcc
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="a-ccsd_t_  rhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="a-ccsd_t_  uhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                   },                    }, id="a-ccsd_t_  rohf conv fc sc:   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="a-ccsd_t_  rhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="a-ccsd_t_  uhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                  },                    }, id="a-ccsd_t_  rohf conv ae sc:   mrcc    ", marks=using("mrcc")),

        ###### dfocc
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="a-ccsd_t_  rhf    df fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="a-ccsd_t_  uhf    df fc:   dfocc   ",),  # SEMI-DISABLE
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="a-ccsd_t_  rhf    df ae:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="a-ccsd_t_  uhf    df ae:   dfocc   ",),  # SEMI-DISABLE
        ##
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="a-ccsd_t_  rhf pk/df fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="a-ccsd_t_  uhf pk/df fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="a-ccsd_t_  rhf pk/df ae:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="a-ccsd_t_  uhf pk/df ae:   dfocc   ",),
        ##
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},                    }, id="a-ccsd_t_  rhf cd/df fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "cd",},                    }, id="a-ccsd_t_  uhf cd/df fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},                    }, id="a-ccsd_t_  rhf cd/df ae:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "cd",},                    }, id="a-ccsd_t_  uhf cd/df ae:   dfocc   ",),
        ####
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="a-ccsd_t_  rhf    cd fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",                   },                    }, id="a-ccsd_t_  uhf    cd fc:   dfocc   ",),  # SEMI-DISABLE
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="a-ccsd_t_  rhf    cd ae:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false",                  },                    }, id="a-ccsd_t_  uhf    cd ae:   dfocc   ",),  # SEMI-DISABLE
        ##
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="a-ccsd_t_  rhf pk/cd fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "true",  "scf_type": "pk",},                    }, id="a-ccsd_t_  uhf pk/cd fc:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="a-ccsd_t_  rhf pk/cd ae:   dfocc   ",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ",      "freeze_core": "false", "scf_type": "pk",},                    }, id="a-ccsd_t_  uhf pk/cd ae:   dfocc   ",),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "a-ccsd(t)", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     },                   }, id="a-ccsd_t_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p38}}, id="a-ccsd_t_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p38}}, id="a-ccsd_t_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                   }, id="a-ccsd_t_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p38}}, id="a-ccsd_t_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }, "error": {0: _p38}}, id="a-ccsd_t_ rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="a-ccsd_t_  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p99}}, id="a-ccsd_t_  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="a-ccsd_t_ rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="a-ccsd_t_  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p99}}, id="a-ccsd_t_  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="a-ccsd_t_ rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="a-ccsd_t_  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p99}}, id="a-ccsd_t_  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {0: _p34}}, id="a-ccsd_t_ rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="a-ccsd_t_  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"},      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p99}}, id="a-ccsd_t_  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {0: _p34}}, id="a-ccsd_t_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     },                   }, id="a-ccsd_t_  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p38}}, id="a-ccsd_t_  uhf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p38}}, id="a-ccsd_t_ rohf         fc: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    },                   }, id="a-ccsd_t_  rhf         ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }, "error": {0: _p38}}, id="a-ccsd_t_  uhf         ae: dd     "),
        pytest.param({"xptd": {},                        "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }, "error": {0: _p38}}, id="a-ccsd_t_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "a-ccsd(t)", "energy"))


#
#  ,-----.   ,-----. ,-----.,------.    ,-.,--------.,-.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                       `-'          `-'                                    `---' `---'
#  <<<  BCCD(T) Energy


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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd_t_  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd_t_  uhf  conv fc   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="bccd_t_ rohf  conv fc sc: * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd_t_  rhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd_t_  uhf  conv ae   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="bccd_t_ rohf  conv ae sc: * ccenergy",),
        # yapf: enable
    ],
)
def test_bccd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "bccd(t)", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="bccd_t_  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="bccd_t_  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="bccd_t_ rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="bccd_t_  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="bccd_t_  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="bccd_t_ rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_ rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_ rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p36}}, id="bccd_t_ rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p36}}, id="bccd_t_ rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",     }                    }, id="bccd_t_  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "true",     }                    }, id="bccd_t_  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf",                                        "freeze_core": "true",     }                    }, id="bccd_t_ rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",    }                    }, id="bccd_t_  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "false",    }                    }, id="bccd_t_  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf",                                        "freeze_core": "false",    }                    }, id="bccd_t_ rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_bccd_prt_pr_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "bccd(t)", "energy"))


#
#   ,-----. ,-----.,----.     ,------.
#  '  .--./'  .--./'.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |      .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC3 Energy


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
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc3  rhf  conv fc   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc3  uhf  conv fc   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "true",                         },                    }, id="cc3 rohf  conv fc sc: * ccenergy",),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc3  rhf  conv ae   : * ccenergy",),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc3  uhf  conv ae   : * ccenergy",),
        pytest.param({"sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "ccenergy", "freeze_core": "false",                        },                    }, id="cc3 rohf  conv ae sc: * ccenergy",),

        ###### mrcc
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                        },                    }, id="cc3  rhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "true",                        },                    }, id="cc3  uhf  conv fc   :   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                       },                    }, id="cc3  rhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        pytest.param({              "keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "mrcc",     "freeze_core": "false",                       },                    }, id="cc3  uhf  conv ae   :   mrcc    ", marks=using("mrcc")),
        # yapf: enable
    ],
)
def test_cc3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc3", "energy"))


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
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc3  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc3  uhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",     }                    }, id="cc3 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc3  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc3  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",    }                    }, id="cc3 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3  rhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3 rohf    df   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3  rhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",     }, "error": {0: _p37}}, id="cc3 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                                      "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",    }, "error": {0: _p37}}, id="cc3 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "true",     }                    }, id="cc3  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "true",     }                    }, id="cc3  uhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"}, "sdsc": "sc", "keywords": {"reference": "rohf",                                        "freeze_core": "true",     }                    }, id="cc3 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rhf",                                         "freeze_core": "false",    }                    }, id="cc3  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "uhf",                                         "freeze_core": "false",    }                    }, id="cc3  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "ccenergy"},               "keywords": {"reference": "rohf",                                        "freeze_core": "false",    }                    }, id="cc3 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_cc3_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "cc3", "energy"))



#
#   ,-----. ,--.   ,--.,------.  ,---.     ,------.
#  '  .-.  '|   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                        `---' `---'
#  <<<  OMP2 Energy

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
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2  rhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2  uhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2 rohf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  rhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  uhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2 rohf    df   ae: * dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2 rohf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_omp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2", "energy"))


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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2  rhf    conv ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2  uhf    conv ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2 rohf    conv ae: dd     ",             ),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2  rhf    df   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2  uhf    df   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2 rohf    df   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2  rhf    df   ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2  uhf    df   ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2 rohf    df   ae: dd     ",             ),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2  rhf    cd   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2  uhf    cd   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2 rohf    cd   fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2  rhf    cd   ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2  uhf    cd   ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2 rohf    cd   ae: dd     ",             ),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                     }                    }, id="omp2  rhf         fc: dd     ",             ),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                     }                    }, id="omp2  uhf         fc: dd     ",             ),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                     }                    }, id="omp2 rohf         fc: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    }                    }, id="omp2  rhf         ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    }                    }, id="omp2  uhf         ae: dd     ",             ),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    }                    }, id="omp2 rohf         ae: dd     ",             ),
        # yapf: enable
    ],
)
def test_omp2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2", "energy"))


#
#   ,-----. ,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .-.  '|   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  | |  ||  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '-'  '|  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  OMP2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        ###### occ
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2 rohf    conv ae: * occ  ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2  rhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2  uhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2 rohf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2  rhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2  uhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp2_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2 rohf      df ae: * dfocc",),
        # yapf: enable
    ],
)
def test_omp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp2_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2 rohf    cd   ae: dd     "),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                     },                            }, id="omp2  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                     },                            }, id="omp2  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                     },                            }, id="omp2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                    },                            }, id="omp2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                    },                            }, id="omp2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                    },                            }, id="omp2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_omp2_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2", "gradient"))


#
#   ,-----. ,--.   ,--.,------.  ,---.     ,-----.    ,------.
#  '  .-.  '|   `.'   ||  .--. ''.-.  \    |  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  |'.'|  ||  '--' | .-' .'    '--. `\    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  |   |  ||  | --' /   '-..--..--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `--'   `--'`--'     '-----''--'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                   `---' `---'
#  <<<  OMP2.5 Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5  rhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5  uhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5 rohf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  rhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  uhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5 rohf    df   ae: * dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp2.5 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp2.5 rohf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_omp2p5_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2.5", "energy"))


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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp2.5 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2.5  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2.5  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp2.5 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2.5  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2.5  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp2.5 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2.5  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2.5  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp2.5 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2.5  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2.5  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp2.5 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2.5  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2.5  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp2.5 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="omp2.5  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="omp2.5  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {0: _p17}}, id="omp2.5 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }                    }, id="omp2.5  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }                    }, id="omp2.5  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }                    }, id="omp2.5 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_omp2p5_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2.5", "energy"))


#
#   ,-----. ,--.   ,--.,------.  ,---.     ,-----.     ,----.                     ,--.,--.                 ,--.
#  '  .-.  '|   `.'   ||  .--. ''.-.  \    |  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  | |  ||  |'.'|  ||  '--' | .-' .'    '--. `\    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '-'  '|  |   |  ||  | --' /   '-..--..--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `--'   `--'`--'     '-----''--'`----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  OMP2.5 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        ###### occ
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2.5  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2.5  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp2.5 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5  rhf    conv ae: * occ  ",),  # adz grd0 known fail
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5 rohf    conv ae: * occ  ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2.5  rhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2.5  uhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp2.5 rohf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5  rhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5  uhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp2.5 rohf      df ae: * dfocc",),
        # yapf: enable
    ],
)
def test_omp2p5_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2.5", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2.5  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2.5  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp2.5 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2.5  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2.5  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp2.5 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2.5  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2.5  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp2.5 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2.5  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2.5  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp2.5 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2.5  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2.5  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp2.5 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2.5  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2.5  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp2.5 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp2.5  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp2.5  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp2.5 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="omp2.5  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="omp2.5  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   },                            }, id="omp2.5 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_omp2p5_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp2.5", "gradient"))


#
#   ,-----. ,--.   ,--.,------. ,----.     ,------.
#  '  .-.  '|   `.'   ||  .--. ''.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  |'.'|  ||  '--' |  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  |   |  ||  | --' /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `--'   `--'`--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                        `---' `---'
#  <<<  OMP3 Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  rhf    conv fc:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  uhf    conv fc:   occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3 rohf    conv fc:   occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  rhf    conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  uhf    conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3 rohf    conv ae:   occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3  rhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3  uhf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3 rohf    df   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  rhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  uhf    df   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3 rohf    df   ae: * dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="omp3 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="omp3 rohf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_omp3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp3", "energy"))


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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp3  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp3  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    }                    }, id="omp3 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp3  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp3  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     }                    }, id="omp3 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp3  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp3  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    }                    }, id="omp3 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp3  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp3  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }                    }, id="omp3 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp3  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp3  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }                    }, id="omp3 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                         "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                         "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                        "freeze_core": "true",                     }, "error": {0: _p17}}, id="omp3 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                         "freeze_core": "false",                    }                    }, id="omp3  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                         "freeze_core": "false",                    }                    }, id="omp3  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                        "freeze_core": "false",                    }                    }, id="omp3 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_omp3_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp3", "energy"))


#
#   ,-----. ,--.   ,--.,------. ,----.      ,----.                     ,--.,--.                 ,--.
#  '  .-.  '|   `.'   ||  .--. ''.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  | |  ||  |'.'|  ||  '--' |  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '-'  '|  |   |  ||  | --' /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `--'   `--'`--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  OMP3 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        ###### occ
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp3  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp3  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="omp3 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3 rohf    conv ae: * occ  ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp3  rhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp3  uhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="omp3 rohf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3  rhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3  uhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "mp_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="omp3 rohf      df ae: * dfocc",),
        # yapf: enable
    ],
)
def test_omp3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp3", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp3  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp3  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="omp3 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp3  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp3  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "conv",                     "freeze_core": "false",                    },                            }, id="omp3 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp3  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp3  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "true",                     },                            }, id="omp3 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp3  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp3  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "df",                       "freeze_core": "false",                    },                            }, id="omp3 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp3  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp3  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="omp3 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp3  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp3  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "mp_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="omp3 rohf    cd   ae: dd     "),

        ###### default qc_module, mp_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp3  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp3  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="omp3 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="omp3  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="omp3  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   },                            }, id="omp3 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_omp3_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "omp3", "gradient"))


#
#   ,-----. ,------. ,------.,--.   ,--.,------.  ,---.     ,------.
#  '  .-.  '|  .--. '|  .---'|   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  '--'.'|  `--, |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  |\  \ |  `---.|  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `--' '--'`------'`--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                         `---' `---'
#  <<<  OREMP2 Energy

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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2 rohf    conv ae: * occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2 rohf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  uhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2 rohf    df   ae:   dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="oremp2 rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="oremp2 rohf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_oremp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "oremp2", "energy"))


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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="oremp2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="oremp2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="oremp2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="oremp2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="oremp2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="oremp2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="oremp2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="oremp2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="oremp2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="oremp2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="oremp2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="oremp2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="oremp2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="oremp2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="oremp2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="oremp2 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="oremp2  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="oremp2  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {0: _p17}}, id="oremp2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }                    }, id="oremp2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }                    }, id="oremp2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }                    }, id="oremp2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_oremp2_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "oremp2", "energy"))


#
#   ,-----. ,------. ,------.,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .-.  '|  .--. '|  .---'|   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  | |  ||  '--'.'|  `--, |  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '-'  '|  |\  \ |  `---.|  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `--' '--'`------'`--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  OREMP2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        ###### occ
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="oremp2  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="oremp2  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="oremp2 rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2 rohf    conv ae: * occ  ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="oremp2  rhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="oremp2  uhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="oremp2 rohf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2  rhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2  uhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="oremp2 rohf      df ae: * dfocc",),
        # yapf: enable
    ],
)
def test_oremp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "oremp2", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="oremp2  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="oremp2  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="oremp2 rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="oremp2  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="oremp2  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="oremp2 rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="oremp2  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="oremp2  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="oremp2 rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="oremp2  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="oremp2  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="oremp2 rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="oremp2  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="oremp2  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="oremp2 rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="oremp2  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="oremp2  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="oremp2 rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="oremp2  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="oremp2  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="oremp2 rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="oremp2  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="oremp2  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   },                            }, id="oremp2 rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_oremp2_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "oremp2", "gradient"))


#
#   ,-----. ,--.    ,-----. ,-----.,------.      ,------.
#  '  .-.  '|  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  ||  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '|  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                              `---' `---'
#  <<<  OLCCD Energy


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
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "olccd", "energy"))


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
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd rohf    conv ae: * occ  ",),
        #
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd  rhf  v conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd  uhf  v conv ae:   occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ", "freeze_core": "false",                    "ccl_energy": "true", "tpdm_abcd_type": "compute"},              }, id="olccd rohf  v conv ae:   occ  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd  rhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd  uhf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd rohf    df   fc:   dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  rhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  uhf    df   ae:   dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd rohf    df   ae:   dfocc",),
        ##
        # skipping pk/df for now
        ##
        # skipping cd/df for now
        ####
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd  rhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd  uhf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     },                   }, id="olccd rohf    cd   fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  rhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd  uhf    cd   ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "cd",   "qc_module": "occ", "freeze_core": "false",                    },                   }, id="olccd rohf    cd   ae: * dfocc",),
        ##
        # skipping pk/cd for now
        # yapf: enable
    ],
)
def test_olccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "olccd", "energy"))


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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {0: _p17}}, id="olccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="olccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="olccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    }                    }, id="olccd rohf    conv ae: dd     "),
        ####
        pytest.param({                              "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="olccd  rhf    df   fc: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="olccd  uhf    df   fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     }                    }, id="olccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="olccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="olccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    }                    }, id="olccd rohf    df   ae: dd     "),
        ####
        pytest.param({                              "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="olccd  rhf    cd   fc: dd     "),
        pytest.param({                              "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="olccd  uhf    cd   fc: dd     "),
        pytest.param({                              "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }                    }, id="olccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="olccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="olccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }                    }, id="olccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="olccd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {0: _p17}}, id="olccd  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {0: _p17}}, id="olccd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   }                    }, id="olccd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   }                    }, id="olccd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   }                    }, id="olccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_olccd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "olccd", "energy"))



#
#   ,-----. ,--.    ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  '  .-.  '|  |   '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  | |  ||  |   |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '-'  '|  '--.'  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----' `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  OLCCD Gradient

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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

        ###### occ
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="olccd  rhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="olccd  uhf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ",      "freeze_core": "true", }, "error": {1: _p17, 0: _p17},}, id="olccd rohf    conv fc: * occ  ",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd  rhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd  uhf    conv ae: * occ  ",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "conv", "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd rohf    conv ae: * occ  ",),

        ###### dfocc
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="olccd  rhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="olccd  uhf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ",      "freeze_core": "true", },                             }, id="olccd rohf      df fc: * dfocc",),
        pytest.param({"keywords": {"reference": "rhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd  rhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "uhf",  "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd  uhf      df ae: * dfocc",),
        pytest.param({"keywords": {"reference": "rohf", "cc_type": "df",   "qc_module": "occ",      "freeze_core": "false",},                             }, id="olccd rohf      df ae: * dfocc",),
        # yapf: enable
    ],
)
def test_olccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "olccd", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="olccd  rhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="olccd  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "true",                     }, "error": {1: _p17, 0: _p17}}, id="olccd rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="olccd  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="olccd  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "conv",                     "freeze_core": "false",                    },                            }, id="olccd rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="olccd  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="olccd  uhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "true",                     },                            }, id="olccd rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="olccd  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="olccd  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "df",                       "freeze_core": "false",                    },                            }, id="olccd rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="olccd  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="olccd  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "true",                     }, "error": {1: _p33}         }, id="olccd rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="olccd  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",  "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="olccd  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf", "cc_type": "cd",                       "freeze_core": "false",                    }, "error": {1: _p33}         }, id="olccd rohf    cd   ae: dd     "),

        ###### default qc_module, cc_type
        pytest.param({"xptd": {},                   "keywords": {"reference": "rhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="olccd  rhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "uhf",                                          "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="olccd  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                         "freeze_core": "true",                    }, "error": {1: _p17, 0: _p17}}, id="olccd rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false",                   },                            }, id="olccd  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false",                   },                            }, id="olccd  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "occ"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false",                   },                            }, id="olccd rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_olccd_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "olccd", "gradient"))


#
#   ,---.,--.   ,--.,--.   ,--.,--.  ,--.    ,------.
#  '   .-'\  `.'  / |  |   |  ||  ,'.|  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  `.  `-. \     /  |  |.'.|  ||  |' '  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  .-'    | \   /   |   ,'.   ||  | `   |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----'   `-'    '--'   '--'`--'  `--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                          `---' `---'
#  <<<  SVWN Energy

_psi_grid = {"dft_radial_points": 99, "dft_spherical_points": 590}
_psi_grid_dd = {"dft_radial_points": 99, "dft_spherical_points": 590}


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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="svwn  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="svwn rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_svwn_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "energy"))



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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="svwn  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="svwn  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="svwn rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="svwn  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="svwn  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="svwn rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="svwn  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="svwn  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="svwn rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",},                     }, id="svwn  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",},                     }, id="svwn  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                       "freeze_core": "false",}, "error": {0: _p29}  }, id="svwn rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_svwn_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "energy"))


#
#   ,---.,--.   ,--.,--.   ,--.,--.  ,--.     ,----.                     ,--.,--.                 ,--.
#  '   .-'\  `.'  / |  |   |  ||  ,'.|  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  `.  `-. \     /  |  |.'.|  ||  |' '  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  .-'    | \   /   |   ,'.   ||  | `   |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `-----'   `-'    '--'   '--'`--'  `--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  SVWN Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="svwn rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="svwn rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="svwn rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="svwn rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="svwn  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="svwn rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p10}         }, id="svwn  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p10}         }, id="svwn  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p10, 0: _p29}}, id="svwn rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_svwn_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="svwn rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="svwn rohf    df   ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p10},          }, id="svwn  rhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p10},          }, id="svwn  uhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p10, 0: _p29}  }, id="svwn rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                              }, id="svwn  uhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="svwn rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_svwn_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "gradient"))


#
#   ,---.,--.   ,--.,--.   ,--.,--.  ,--.    ,--.  ,--.                     ,--.
#  '   .-'\  `.'  / |  |   |  ||  ,'.|  |    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  `.  `-. \     /  |  |.'.|  ||  |' '  |    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  .-'    | \   /   |   ,'.   ||  | `   |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `-----'   `-'    '--'   '--'`--'  `--'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  SVWN Hessian


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2", marks=pytest.mark.d2ints),
        pytest.param(1, id="hes1", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
        pytest.param(0, id="hes0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                  }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                  }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                  }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                  }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                  },                                     }, id="svwn  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                  }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf disk ae:   scf  ",),
        ####
        # only H-by-E available for CD ref Hessians, and loose dflt cholesky_tolerance means they're not close to CONV, so skipping for now
        # pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                  }, "error": {2: _p10, 1: _p10}}, id="svwn  rhf   cd ae:   scf  ",),
        # pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                  }, "error": {2: _p10, 1: _p10}}, id="svwn  uhf   cd ae:   scf  ",),
        # pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                  }, "error": {2: _p10, 1: _p10}}, id="svwn rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_svwn_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "hessian"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2", marks=pytest.mark.d2ints),
        pytest.param(1, id="hes1", marks=pytest.mark.findif),
        pytest.param(0, id="hes0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf    df   ae: dd     ",),
        ####
        # only H-by-E available for CD ref Hessians, and loose dflt cholesky_tolerance means they're not close to CONV, so skipping for now
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p10, 1: _p10},  }, id="svwn  rhf    cd   ae: dd     ",),
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p10, 1: _p10},  }, id="svwn  uhf    cd   ae: dd     ",),
        # pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p10, 1: _p10},  }, id="svwn rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                                     }, id="svwn  uhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false", **_psi_grid_dd                   }, "error": {2: _p29, 1: _p29, 0: _p29}}, id="svwn rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_svwn_hessian_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "svwn", "hessian"))


#
#  ,------. ,-----.  ,------.    ,------.
#  |  .--. '|  |) /_ |  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  '--' ||  .-.  \|  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  | --' |  '--' /|  `---.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `------' `------'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  PBE Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="pbe  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="pbe rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_pbe_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "pbe", "energy"))



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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="pbe  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="pbe  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="pbe rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="pbe  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="pbe  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="pbe rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="pbe  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="pbe  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="pbe rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",},                     }, id="pbe  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",},                     }, id="pbe  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                       "freeze_core": "false",}, "error": {0: _p29}  }, id="pbe rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_pbe_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "pbe", "energy"))


#
#  ,------. ,-----.  ,------.     ,----.                     ,--.,--.                 ,--.
#  |  .--. '|  |) /_ |  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  '--' ||  .-.  \|  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  | --' |  '--' /|  `---.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'     `------' `------'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  PBE Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="pbe  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30}         }, id="pbe  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30}         }, id="pbe  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30, 0: _p29}}, id="pbe rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_pbe_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "pbe", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30},        }, id="pbe  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30},        }, id="pbe  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30, 0: _p29}}, id="pbe rohf    cd   ae: dd     "),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                            }, id="pbe  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="pbe rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_pbe_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "pbe", "gradient"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //'-'  ||  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' `----' `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                            `---' `---'
#  <<<  B3LYP Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="b3lyp  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b3lyp rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_b3lyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b3lyp", "energy"))



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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="b3lyp  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="b3lyp  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b3lyp rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="b3lyp  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="b3lyp  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b3lyp rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="b3lyp  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="b3lyp  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b3lyp rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",},                     }, id="b3lyp  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",},                     }, id="b3lyp  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                       "freeze_core": "false",}, "error": {0: _p29}  }, id="b3lyp rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_b3lyp_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b3lyp", "energy"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.      ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //'-'  ||  '--. |  |   |  | --'     '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' `----' `-----' `--'   `--'          `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B3LYP Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="b3lyp  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30}         }, id="b3lyp  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30}         }, id="b3lyp  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30, 0: _p29}}, id="b3lyp rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_b3lyp_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b3lyp", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  uhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  uhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30},        }, id="b3lyp  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30},        }, id="b3lyp  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30, 0: _p29}}, id="b3lyp rohf    cd   ae: dd     "),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                            }, id="b3lyp  uhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}}, id="b3lyp rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_b3lyp_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b3lyp", "gradient"))


#
#             ,-----.   ,---. ,-----.,--.   ,--.    ,------.
#  ,--.   ,--.|  |) /_ | o   \'--,  / \  `.'  /     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |.'.|  ||  .-.  \`..'  | .'  /   .'    \      |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |   .'.   ||  '--' / .'  / /   /   /  .'.  \     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  '--'   '--'`------'  `--'  `--'   '--'   '--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                 `---' `---'
#  <<<  wB97X Energy

# wB97X req
# _psi_grid = {"dft_radial_points": 500, "dft_spherical_points": 1202}

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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="wb97x  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p31}}, id="wb97x  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p31}}, id="wb97x  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="wb97x rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_wb97x_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "wb97x", "energy"))


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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="wb97x  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="wb97x  uhf    conv ae: dd     "),
        pytest.param({"xptd": {                  }, "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="wb97x rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="wb97x  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="wb97x  uhf    df   ae: dd     "),
        pytest.param({"xptd": {                  }, "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="wb97x rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {                  }, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p31}  }, id="wb97x  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {                  }, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p31}  }, id="wb97x  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {                  }, "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="wb97x rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",},                     }, id="wb97x  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",},                     }, id="wb97x  uhf         ae: dd     "),
        pytest.param({"xptd": {                  }, "keywords": {"reference": "rohf",                                       "freeze_core": "false",}, "error": {0: _p29}  }, id="wb97x rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_wb97x_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "wb97x", "energy"))


#
#             ,-----.   ,---. ,-----.,--.   ,--.     ,----.                     ,--.,--.                 ,--.
#  ,--.   ,--.|  |) /_ | o   \'--,  / \  `.'  /     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |.'.|  ||  .-.  \`..'  | .'  /   .'    \      |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |   .'.   ||  '--' / .'  / /   /   /  .'.  \     '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  '--'   '--'`------'  `--'  `--'   '--'   '--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  wB97X Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=[pytest.mark.nonroutine, pytest.mark.findif]),
    ],
)
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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="wb97x rohf   pk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="wb97x rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="wb97x rohf   df ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="wb97x rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                            }, id="wb97x  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {1: _p29, 0: _p29}}, id="wb97x rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30, 0: _p31}}, id="wb97x  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30, 0: _p31}}, id="wb97x  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {1: _p30, 0: _p29}}, id="wb97x rohf   cd ae:   scf  ",),
        # yapf: enable
    ],
)
def test_wb97x_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "wb97x", "gradient"))


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        pytest.param(0, id="grd0", marks=pytest.mark.findif),
    ],
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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  rhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  uhf    conv ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "pk",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="wb97x rohf    conv ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  rhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  uhf    df   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "df",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="wb97x rohf    df   ae: dd     ",),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30, 0: _p31}  }, id="wb97x  rhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30, 0: _p31}  }, id="wb97x  uhf    cd   ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf", "scf_type": "cd",                       "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p30, 0: _p29}  }, id="wb97x rohf    cd   ae: dd     ",),

        ###### default qc_module, mp2_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  rhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                          "freeze_core": "false", **_psi_grid_dd                   },                              }, id="wb97x  uhf         ae: dd     ",),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rohf",                                         "freeze_core": "false", **_psi_grid_dd                   }, "error": {1: _p29, 0: _p29}  }, id="wb97x rohf         ae: dd     ",),
        # yapf: enable
    ],
)
def test_wb97x_gradient_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "wb97x", "gradient"))


#
#  ,-----.   ,---. ,------. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  \|  .--. '|  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \ .-' .'|  '--' ||  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //   '-.|  | --' |  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' '-----'`--'     `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                     `---' `---'
#  <<<  B2PLYP Energy


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
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf   pk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf   pk ae:   scf  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "pk",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  rhf   pk fc:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "pk",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  uhf   pk fc:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "pk",      "freeze_core": "true", **_psi_grid                    }, "error": {0: _p29}}, id="b2plyp rohf   pk fc:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf drct ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "direct",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf drct ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf   df ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf   df ae:   scf  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "df",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  rhf   df fc:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "df",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  uhf   df fc:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "df",      "freeze_core": "true", **_psi_grid                    }, "error": {0: _p29}}, id="b2plyp rohf   df fc:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf  mem ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "mem_df",  "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf  mem ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf disk ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "disk_df", "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf disk ae:   scf  ",),
        ####
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  rhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   },                   }, id="b2plyp  uhf   cd ae:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "false", **_psi_grid                   }, "error": {0: _p29}}, id="b2plyp rohf   cd ae:   scf  ",),
        ##
        pytest.param({"keywords": {"reference": "rhf",  "scf_type": "cd",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  rhf   cd fc:   scf  ",),
        pytest.param({"keywords": {"reference": "uhf",  "scf_type": "cd",      "freeze_core": "true", **_psi_grid                    },                   }, id="b2plyp  uhf   cd fc:   scf  ",),
        pytest.param({"keywords": {"reference": "rohf", "scf_type": "cd",      "freeze_core": "true", **_psi_grid                    }, "error": {0: _p29}}, id="b2plyp rohf   cd fc:   scf  ",),
        # yapf: enable
    ],
)
def test_b2plyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b2plyp", "energy"))



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
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "true", },                     }, id="b2plyp  rhf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "true", },                     }, id="b2plyp  uhf    conv fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "true", }, "error": {0: _p29}  }, id="b2plyp rohf    conv fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="b2plyp  rhf    conv ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "pk",                     "freeze_core": "false",},                     }, id="b2plyp  uhf    conv ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "pk",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b2plyp rohf    conv ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "true", },                     }, id="b2plyp  rhf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "true", },                     }, id="b2plyp  uhf    df   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "true", }, "error": {0: _p29}  }, id="b2plyp rohf    df   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="b2plyp  rhf    df   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "df",                     "freeze_core": "false",},                     }, id="b2plyp  uhf    df   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "df",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b2plyp rohf    df   ae: dd     "),
        ####
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "true", },                     }, id="b2plyp  rhf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "true", },                     }, id="b2plyp  uhf    cd   fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "true", }, "error": {0: _p29}  }, id="b2plyp rohf    cd   fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="b2plyp  rhf    cd   ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",  "scf_type": "cd",                     "freeze_core": "false",},                     }, id="b2plyp  uhf    cd   ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf", "scf_type": "cd",                     "freeze_core": "false",}, "error": {0: _p29}  }, id="b2plyp rohf    cd   ae: dd     "),

        ###### default qc_module, scf_type
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "true", },                     }, id="b2plyp  rhf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "true", },                     }, id="b2plyp  uhf         fc: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                       "freeze_core": "true", }, "error": {0: _p29}  }, id="b2plyp rohf         fc: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "rhf",                                        "freeze_core": "false",},                     }, id="b2plyp  rhf         ae: dd     "),
        pytest.param({"xptd": {"qc_module": "scf"}, "keywords": {"reference": "uhf",                                        "freeze_core": "false",},                     }, id="b2plyp  uhf         ae: dd     "),
        pytest.param({"xptd": {},                   "keywords": {"reference": "rohf",                                       "freeze_core": "false",}, "error": {0: _p29}  }, id="b2plyp rohf         ae: dd     "),
        # yapf: enable
    ],
)
def test_b2plyp_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "b2plyp", "energy"))


#
#  ,-----.   ,---. ,------. ,--.,--.   ,--.,------.      ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  \|  .--. '|  | \  `.'  / |  .--. '    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \ .-' .'|  '--' ||  |  '.    /  |  '--' |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //   '-.|  | --' |  '--. |  |   |  | --'     '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' '-----'`--'     `-----' `--'   `--'          `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B2PLYP  Gradient


#
#   ,-----.
#  '  .--./ ,---. ,--,--,--.,--,--,--. ,---. ,--,--,
#  |  |    | .-. ||        ||        || .-. ||      \
#  '  '--'\' '-' '|  |  |  ||  |  |  |' '-' '|  ||  |
#   `-----' `---' `--`--`--'`--`--`--' `---' `--''--'
#
#  <<<  Common functions


def _processor(inp, dertype, basis, subjects, clsd_open_pmols, request, method, driver):
    """For each pytest.param in this file, takes its input dictionary, inserts additional param
    fields like molecule, basis, dertype, inserts per-test fields like method and driver, sets up
    handling routing like per-basis or per-dertype xfail flags, and sends the result on to
    standard_suite_runner.py

    """
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["keywords"]["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k not in ["error", "wrong", "xptd"]}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("wrong", False) and inp["wrong"].get(dertype, False):
        # clause unused in psi4
        inpcopy["wrong"] = inp["wrong"][dertype]
    if inp.get("error", False) and inp["error"].get(basis, False):
        # clause unused in psi4
        inpcopy["error"] = inp["error"][basis]
    if inp.get("wrong", False) and inp["wrong"].get(basis, False):
        # clause unused in psi4
        inpcopy["wrong"] = inp["wrong"][basis]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    xptd = {}
    for k, v in inp.get("xptd", {}).items():
        if picked_v:= ((v if isinstance(v, (str, bool)) else False) or v.get(dertype, False)):
            xptd[k] = picked_v
    inpcopy["xptd"] = xptd

    inpcopy["driver"] = driver
    inpcopy["call"] = method
    inpcopy["keywords"]["basis"] = basis
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}
    # print("INP", inpcopy)

    return inpcopy, subject, method, basis, tnm
