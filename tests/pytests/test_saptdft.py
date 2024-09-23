import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values

hartree_to_kcalmol = constants.conversion_factor("hartree", "kcal/mol")
pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.mark.saptdft
def test_saptdft_auto_grac():
    mol_dimer = psi4.geometry(
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
"""
    )
    psi4.set_options(
        {
            "level_shift": 0.060000000000000005,
            "level_shift_cutoff": 0.001,
            "MAXITER": 300,
            "SCF_INITIAL_ACCELERATOR": "ADIIS",
            "basis": "aug-cc-pvdz",
            "sapt_dft_grac_shift_a": -99,
            "sapt_dft_grac_shift_b": -99,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    compare_values(
        0.13053068183319516,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_A"),
        6,
        "SAPT_DFT_GRAC_SHIFT_A",
    )
    compare_values(
        0.13063798506967816,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_B"),
        6,
        "SAPT_DFT_GRAC_SHIFT_B",
    )
    return


@pytest.mark.saptdft
def test_saptdft_monomer_grac():
    mol = psi4.geometry(
        """
0 1
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000

units angstrom
"""
    )
    psi4.set_options(
        {
            "basis": "sto-3g",
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol)
    return


if __name__ == "__main__":
    test_saptdft_auto_grac()
