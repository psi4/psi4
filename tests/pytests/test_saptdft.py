import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values

hartree_to_kcalmol = constants.conversion_factor("hartree", "kcal/mol")
pytestmark = [pytest.mark.psi, pytest.mark.api]

_sapt_testing_mols = {
    "neutral_water_dimer": """
0 1
8   -0.702196054   -0.056060256   0.009942262
1   -1.022193224   0.846775782   -0.011488714
1   0.257521062   0.042121496   0.005218999
--
0 1
8   2.268880784   0.026340101   0.000508029
1   2.645502399   -0.412039965   0.766632411
1   2.641145101   -0.449872874   -0.744894473
units angstrom
""",
    "hydroxide": """
-1 1
8   -0.702196054   -0.056060256   0.009942262
1   -1.022193224   0.846775782   -0.011488714
--
0 1
8   2.268880784   0.026340101   0.000508029
1   2.645502399   -0.412039965   0.766632411
1   2.641145101   -0.449872874   -0.744894473
units angstrom
""",
}


@pytest.mark.saptdft
@pytest.mark.parametrize(
    "SAPT_DFT_GRAC_COMPUTE, refA, refB, gracA, gracB, geometry, grac_basis",
    [
        (
            "SINGLE",
            0.19807358,
            0.19830016,
            None,
            None,
            'neutral_water_dimer',
            None,
        ),
        (
            "SINGLE",
            0.1307,
            0.19830016,
            0.1307,
            None,
            'neutral_water_dimer',
            None,
        ),
        (
            "ITERATIVE",
            0.2303073898,
            0.19830016,
            None,
            None,
            "hydroxide",
            None,
        ),
        (
            "ITERATIVE",
            0.0924377691,
            0.1306379832,
            None,
            None,
            "hydroxide",
            "aug-cc-pvdz",
        ),
    ],
)
def test_saptdft_auto_grac(SAPT_DFT_GRAC_COMPUTE, refA, refB, gracA, gracB, geometry, grac_basis):
    """
    For SAPT(DFT), one must compute a GRAC shift for each monomer. Ideally,
    this GRAC shift should be close to the experimental Ionization Potential
    (IP) of the monomer. Basis set incompleteness prevents this here.
    e.g., aug-DZ H2O has a shift of 0.1306, compared to 0.1307 experimental IP.
    """
    mol_dimer = psi4.geometry(_sapt_testing_mols[geometry])
    psi4.set_options(
        {
            "basis": "STO-3G",
            "SAPT_DFT_FUNCTIONAL": "pbe0",
            "SAPT_DFT_GRAC_COMPUTE": SAPT_DFT_GRAC_COMPUTE,
        }
    )
    if grac_basis is not None:
        psi4.set_options({"SAPT_DFT_GRAC_BASIS": grac_basis})
    if gracA is not None:
        psi4.set_options({"SAPT_DFT_GRAC_SHIFT_A": gracA})
    if gracB is not None:
        psi4.set_options({"SAPT_DFT_GRAC_SHIFT_B": gracB})
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    compare_values(
        refA,
        psi4.core.variable("SAPT DFT GRAC SHIFT A"),
        8,
        "SAPT DFT GRAC SHIFT A",
    )
    assert compare_values(
        refB,
        psi4.core.variable("SAPT DFT GRAC SHIFT B"),
        8,
        "SAPT DFT GRAC SHIFT B",
    )


if __name__ == "__main__":
    test_saptdft_auto_grac(
            "ITERATIVE",
            0.0924377691,
            0.1306379832,
            None,
            None,
            "hydroxide",
            "aug-cc-pvdz",
        )

