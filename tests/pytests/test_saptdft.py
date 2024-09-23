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
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
--
0 1
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561

units angstrom
"""
    )
    psi4.set_options(
        {
            # "basis": "sto-3g",
            "basis": "aug-cc-pvtz",
            # "e_convergence": 1e-8,
            # "d_convergence": 1e-8,
            "sapt_dft_grac_shift_a": -99,
            "sapt_dft_grac_shift_b": -99,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)


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


if __name__ == "__main__":
    test_saptdft_auto_grac()
