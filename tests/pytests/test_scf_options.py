import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"options": {"guess": "core"}}, id="core"),
        pytest.param({"options": {"guess": "gwh"}}, id="gwh"),
        pytest.param({"options": {"guess": "huckel"}}, id="huckel"),
        pytest.param({"options": {"guess": "modhuckel"}}, id="modhuckel"),
        pytest.param({"late_options": {"guess": "read"}}, id="read"),
        pytest.param({"options": {"guess": "sad"}}, id="sad"),
        pytest.param({"options": {"guess": "sadno"}}, id="sadno"),
        pytest.param({"options": {"guess": "sap"}}, id="sap"),
    ],
)
def test_guess_mix_for_broken_symmetry(inp):

    refENuc = 0.17639240356
    refSCF = -0.82648407827446
    refBSSCF = -0.99872135103903

    h2 = psi4.geometry(
        """
        0 1
        H
        H 1 3.0
        symmetry c1
        """
    )

    psi4.set_options({"reference": "uhf", "e_convergence": 12, "basis": "cc-pvdz"})
    psi4.set_options(inp.get("options", {}))

    thisSCF = psi4.energy("scf")
    psi4.set_options(inp.get("late_options", {}))
    psi4.set_options({"guess_mix": True})
    thisBSSCF = psi4.energy("scf")

    assert compare_values(
        refENuc, h2.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy"
    )
    assert compare_values(refSCF, thisSCF, 10, "Reference energy")
    assert compare_values(refBSSCF, thisBSSCF, 10, "Reference broken-symmetry energy")


@pytest.mark.parametrize("ref", ["rhf", "uhf", "rohf"])
@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"options": {"guess": "core"}}, id="core"),
        pytest.param({"options": {"guess": "gwh"}}, id="gwh"),
        pytest.param({"options": {"guess": "huckel"}}, id="huckel"),
        pytest.param({"options": {"guess": "modhuckel"}}, id="modhuckel"),
        pytest.param({"options": {"guess": "sad"}}, id="sad"),
        pytest.param({"options": {"guess": "sadno"}}, id="sadno"),
        pytest.param({"options": {"guess": "sap"}}, id="sap"),
        # pytest.param({'options': {"guess": "sapgau"}}, id="sapgau"),
    ],
)
def test_scf_guess(inp, ref):
    vals = {
        "rhf": {
            "core": -85.33505416384983,
            "gwh": -92.48691177199250,
            "sad": -99.63941801281894,
            "sadno": -99.77439003262768,
            "huckel": -99.82915309755559,
            "modhuckel": -99.80751245771735,
            "sap": -99.97316521623380,
        },
        "uhf": {
            "core": -88.99202896495355,
            "gwh": -94.35585144644561,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641045,
            "huckel": -99.20374047108815,
            "modhuckel": -99.17689049002892,
            "sap": -99.50144270485751,
        },
        "rohf": {
            "core": -88.99202896495353,
            "gwh": -94.35585144644557,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641042,
            "huckel": -99.20374047108817,
            "modhuckel": -99.17689049002892,
            "sap": -99.50144270485754,
        },
    }
    ref_final = {
        # "refnuc": 5.282967161430950,
        "rhf": -100.0584459442408587,
        "uhf": -99.5312257221445549,
        "rohf": -99.5261713512123976,
        "cuhf": -99.5261713512123976,
    }
    options = {
            "reference": ref,
            "basis": "cc-pvtz",
            "scf_type": "pk",
            "df_scf_guess": False,
            "guess": inp["options"]["guess"],
        }
    if ref == "rhf":
        hf = psi4.geometry(
            """
    0 1
    F
    H 1 0.9015
    """
        )
        options["docc"] = [3, 0, 1, 1]
    else:
        hf = psi4.geometry(
            """
    1 2
    F
    H 1 0.9015
    """
        )
        options["docc"] = [3, 0, 0, 1]
        options["socc"] = [0, 0, 1, 0]

    psi4.set_options(options)
    psi4.energy("scf")
    assert compare_values(
        ref_final[ref], psi4.variable("SCF TOTAL ENERGY"), 6, "FINAL SCF ENERGY"
    )
    assert compare_values(
        vals[ref][inp["options"]["guess"]],
        psi4.variable("SCF TOTAL ENERGIES")[0],
        6,
        "INITIAL ITERATION",
    )
