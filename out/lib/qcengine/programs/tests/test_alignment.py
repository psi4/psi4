import numpy as np
import pytest
import qcelemental as qcel

from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs
from qcengine.testing import using

from .standard_suite_runner import runner_asserter
from .test_standard_suite import _processor


@pytest.fixture
def clsd_open_pmols():
    frame_not_important = {
        name[:-4]: qcel.models.Molecule.from_data(smol, name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }
    frame_part_of_spec = {
        name[:-4] + "-fixed": qcel.models.Molecule.from_data(smol + "\nno_com\nno_reorient\n", name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }

    return {**frame_not_important, **frame_part_of_spec}


#
#  ,--.  ,--.,------.      ,---.  ,--.,--.                                          ,--.
#  |  '--'  ||  .---'     /  O  \ |  |`--' ,---. ,--,--, ,--,--,--. ,---. ,--,--, ,-'  '-.
#  |  .--.  ||  `--,     |  .-.  ||  |,--.| .-. ||      \|        || .-. :|      \'-.  .-'
#  |  |  |  ||  |`       |  | |  ||  ||  |' '-' '|  ||  ||  |  |  |\   --.|  ||  |  |  |
#  `--'  `--'`--'        `--' `--'`--'`--'.`-  / `--''--'`--`--`--' `----'`--''--'  `--'
#                                         `---'
#  <<<  HF Alignment


@pytest.mark.parametrize(
    "scramble",
    # * this parameter alters the input molecule by Cartesian frame and/or atom ordering to imitate real-world inputs.
    # * scramble dictionary are arguments to qcdb.Molecule.scramble() or qcel.models.Molecule.scramble() that computes a shifted/rotated/atom-mapped input molecule from `subjects` below and the transformations to be applied to the reference data.
    # * arguments do_shift=T/F and other boolean values generate random perturbations, so multiple lines or runs makes a fuller test.
    # * specific shifts etc or problematic runs can be reproduced by specifying arrays like the commented example below. note that b/c do_resort is natom-dependent, may need to exclude subjects while debugging.
    # * `id`s below are numbered since `pytest -k` doesn't recognize case and so duplicate entries can be added. the 0, 1, 2, 3, 9 progression is handy for debugging since perturbations are added systematically
    [
        # pytest.param({"do_shift": [ 2.660760432055,  1.477336796939, -2.098045335573], "do_rotate": [[ 0.321861140022,  0.445246880671,  0.835560064745], [-0.447874157747,  0.849136107328, -0.279958229125], [-0.834154749052, -0.28411808546 ,  0.472718487209]], "do_resort": [2, 0, 1]}, id="mTTT3"),
        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": False}, id="srm0"),
        # pytest.param({"do_shift": True, "do_rotate": False, "do_resort": False}, id="Srm1"),
        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": False}, id="SRm2"),
        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": True}, id="srM3"),
        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM8"),
        pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM9"),
    ],
)
@pytest.mark.parametrize(
    "frame",
    # * this parameter alters the input molecule by fix_com and fix_orientation to imitate user signaling frame matters or not.
    [
        pytest.param("fixed"),  # fix_=True (no_com/no_reorient); atres.mol.geom = atin.mol.geom aka scrambled
        pytest.param("free"),  # fix_=False (def)               ; atres.mol.geom = dsl internal orientation
    ],
)
@pytest.mark.parametrize(
    "driver",
    [
        pytest.param("energy", id="ene0"),
        pytest.param("gradient", id="grd1"),
        pytest.param("hessian", id="hes2"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    # this parameter, along with rhf/uhf below covers four different molecules, with nat=2-4.
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"scf_conv": 12},                                                                     }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                   }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {"scf__thresh": 1.e-6},                                                               }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {"scf_type": "pk"},                                                                   }, id="hf  rhf ae: psi4",       marks=using("psi4_derqcsk")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae", "keywords": {"contrl__scftyp": "uhf"},                                                            }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {"scf__uhf": True, "scf__thresh": 1.e-6},                                             }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_type": "pk"},                                               }, id="hf  uhf ae: psi4",       marks=using("psi4_derqcsk")),
        # yapf: enable
    ],
)
def test_hf_alignment(inp, scramble, frame, driver, basis, subjects, clsd_open_pmols, request):
    runner_asserter(
        *_processor(inp, "", basis, subjects, clsd_open_pmols, request, driver, "hf", scramble=scramble, frame=frame)
    )
