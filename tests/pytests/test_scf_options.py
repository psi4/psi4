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

    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 5e-9})

    thisSCF = psi4.energy("scf")
    psi4.set_options(inp.get("late_options", {}))
    psi4.set_options({"guess_mix": True})
    thisBSSCF = psi4.energy("scf")

    assert compare_values(
        refENuc, h2.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy"
    )
    assert compare_values(refSCF, thisSCF, 10, "Reference energy")
    assert psi4.compare_values(refBSSCF, thisBSSCF, 10, "Reference broken-symmetry energy")


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
            "core": -85.33505416385,
            "gwh": -92.48691177199,
            "sad": -100.28117660005705,
            "sadno": -99.60935769112,
            "huckel": -99.82044374140,
            "modhuckel": -99.79629341034,
            "sap": -99.97316521623,
        },
        "uhf": {
            "core": -88.99202896495,
            "gwh": -94.35585144645,
            "sad": -100.28117660005705,
            "sadno": -98.79741688423,
            "huckel": -99.12283156385,
            "modhuckel": -99.09435147809,
            "sap": -99.50144270486,
        },
        "rohf": {
            "core": -88.99202896495,
            "gwh": -94.35585144645,
            "sad": -100.28117660005705,
            "sadno": -98.79741688423,
            "huckel": -99.12283156385,
            "modhuckel": -99.09435147809,
            "sap": -99.50144270486,
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
    if psi4.core.get_option("scf", "orbital_optimizer_package") == "INTERNAL":
      # until pyooo and SCF TOTAL ENERGIES filled out
      assert compare_values(
        vals[ref][inp["options"]["guess"]],
        psi4.variable("SCF TOTAL ENERGIES")[0],
        6,
        "INITIAL ITERATION",
      )
