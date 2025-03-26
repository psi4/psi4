import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values

hartree_to_kcalmol = constants.conversion_factor("hartree", "kcal/mol")
pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.mark.saptdft
@pytest.mark.parametrize(
    "SAPT_DFT_GRAC_COMPUTE, refA, refB, gracA, gracB, geometry",
    [
        (
            "SINGLE",
            0.19807358,
            0.19830016,
            0.0,
            0.0,
            """
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
        ),
        (
            "SINGLE",
            0.1307,
            0.19830016,
            0.1307,
            0.0,
            """
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
        ),
        (
            "ITERATIVE",
            0.3258340368,
            0.19830016,
            0.0,
            0.0,
            """
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
        ),
    ],
)
def test_saptdft_auto_grac(SAPT_DFT_GRAC_COMPUTE, refA, refB, gracA, gracB, geometry):
    """
    For SAPT(DFT), one must compute a GRAC shift for each monomer. Ideally,
    this GRAC shift should be close to the experimental Ionization Potential
    (IP) of the monomer. While the present test targets STO-3G, the IP for
    water is 0.1307. Note that using aug-cc-pVDZ, the computed GRAC shift is
    0.13063798506967816, which is close to the experimental value.
    """
    mol_dimer = psi4.geometry(geometry)
    psi4.set_options(
        {
            "basis": "STO-3G",
            "sapt_dft_grac_shift_a": gracA,
            "sapt_dft_grac_shift_b": gracB,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
            "SAPT_DFT_GRAC_COMPUTE": SAPT_DFT_GRAC_COMPUTE,
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    compare_values(
        refA,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_A"),
        8,
        "SAPT_DFT_GRAC_SHIFT_A",
    )
    compare_values(
        refB,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_B"),
        8,
        "SAPT_DFT_GRAC_SHIFT_B",
    )
    return
