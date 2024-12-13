import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values
import numpy as np
import qcelemental as qcel
from pprint import pprint as pp

hartree_to_kcalmol = constants.conversion_factor("hartree", "kcal/mol")
pytestmark = [pytest.mark.psi, pytest.mark.api]


psi4.set_memory('60 GB')

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
            "basis": "STO-3G",
            "sapt_dft_grac_shift_a": -99,
            "sapt_dft_grac_shift_b": -99,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    compare_values(
        #  STO-3G target
        0.19807358,
        # aug-cc-pvdz target, 0.1307 (using experimental IP from CCCBDB)
        # 0.13053068183319516,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_A"),
        8,
        "SAPT_DFT_GRAC_SHIFT_A",
    )
    compare_values(
        #  STO-3G target
        0.19830016,
        # aug-cc-pvdz target, 0.1307 (using experimental IP from CCCBDB)
        # 0.13063798506967816,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_B"),
        8,
        "SAPT_DFT_GRAC_SHIFT_B",
    )
    return


@pytest.mark.saptdft
def test_saptdft_auto_grac_iterative():
    mol_dimer = psi4.geometry(
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
"""
    )
    psi4.set_options(
        {
            "SAPT_DFT_GRAC_CONVERGENCE_TIER": "ITERATIVE",
            "basis": "STO-3G",
            "sapt_dft_grac_shift_a": -99,
            "sapt_dft_grac_shift_b": -99,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    compare_values(
        #  STO-3G target
        0.3258340368,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_A"),
        8,
        "SAPT_DFT_GRAC_SHIFT_A",
    )
    compare_values(
        #  STO-3G target
        0.19830016,
        psi4.core.variable("SAPT_DFT_GRAC_SHIFT_B"),
        8,
        "SAPT_DFT_GRAC_SHIFT_B",
    )
    return


@pytest.mark.saptdft
def test_saptdft_external_potential():
    mol_trimer = psi4.geometry(
        """
0 1
O   0.017225   0.031664   0.004802
H  -0.046691  -0.052504   0.962436
H   0.972017   0.055307  -0.185622
--
0 1
O   2.516175   0.894012  -1.014512
H   1.942080   1.572902  -1.410984
H   3.056412   0.561271  -1.739079
symmetry c1
no_reorient
no_com
    """
    )
    fisapt0_external_potential_energies = {
        "SAPT ELST ENERGY": -0.01581004514947182,
        "SAPT ELST10,R ENERGY": -0.01581004514947182,
        "SAPT EXCH ENERGY": 0.012282520736587468,
        "SAPT IND ENERGY": -0.0035613061462424402,
        "SAPT DISP ENERGY": -0.002185724589094623,
        "SAPT TOTAL ENERGY": -0.009274555148221415,
    }

    saptdft_no_external_potential = {
        "SAPT ELST ENERGY": -0.014201712642446296,
        "ELST10,R": -0.014201712642446296,
        "SAPT EXCH ENERGY": 0.014021550175337915,
        "SAPT IND ENERGY": -0.0033383768785273885,
        "SAPT DISP ENERGY": -0.002394920793165888,
        "SAPT TOTAL ENERGY": -0.005913460138801657,
    }
    # External potential containing the third water from the trimer with TIP3P
    # charges
    Chargefield_C = np.array(
        [
            -0.834,
            0.179217,
            2.438389,
            -1.484606,
            0.417,
            -0.194107,
            1.702697,
            -0.966751,
            0.417,
            -0.426657,
            2.563754,
            -2.222683,
        ]
    ).reshape((-1, 4))
    Chargefield_C[:, [1, 2, 3]] /= qcel.constants.bohr2angstroms
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "sapt_dft_grac_shift_a": 0.1307,
            "sapt_dft_grac_shift_b": 0.1307,
        }
    )
    # psi4.energy("fisapt0", external_potentials={"C": Chargefield_C})
    psi4.energy(
        "sapt(dft)",
        external_potentials={"C": Chargefield_C},
    )
    pp(psi4.core.variables())
    for key, value in fisapt0_external_potential_energies.items():
        compare_values(
            value,
            psi4.core.variable(key),
            8,
            key,
        )
    return


def test_fisapt_external_potential():
    mol_trimer = psi4.geometry(
        """
0 1
O   0.017225   0.031664   0.004802
H  -0.046691  -0.052504   0.962436
H   0.972017   0.055307  -0.185622
--
0 1
O   2.516175   0.894012  -1.014512
H   1.942080   1.572902  -1.410984
H   3.056412   0.561271  -1.739079
symmetry c1
no_reorient
no_com
    """
    )
    fisapt0_external_potential_energies = {
        "SAPT ELST ENERGY": -0.01581004514947182,
        "SAPT ELST10,R ENERGY": -0.01581004514947182,
        "SAPT EXCH ENERGY": 0.012282520736587468,
        "SAPT IND ENERGY": -0.0035613061462424402,
        "SAPT DISP ENERGY": -0.002185724589094623,
        "SAPT TOTAL ENERGY": -0.009274555148221415,
    }

    saptdft_no_external_potential = {
        "SAPT ELST ENERGY": -0.014201712642446296,
        "ELST10,R": -0.014201712642446296,
        "SAPT EXCH ENERGY": 0.014021550175337915,
        "SAPT IND ENERGY": -0.0033383768785273885,
        "SAPT DISP ENERGY": -0.002394920793165888,
        "SAPT TOTAL ENERGY": -0.005913460138801657,
    }
    # External potential containing the third water from the trimer with TIP3P
    # charges
    Chargefield_C = np.array(
        [
            -0.834,
            0.179217,
            2.438389,
            -1.484606,
            0.417,
            -0.194107,
            1.702697,
            -0.966751,
            0.417,
            -0.426657,
            2.563754,
            -2.222683,
        ]
    ).reshape((-1, 4))
    Chargefield_C[:, [1, 2, 3]] /= qcel.constants.bohr2angstroms
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
        }
    )
    psi4.energy(
        "fisapt0",
        external_potentials={"C": Chargefield_C},
    )
    pp(psi4.core.variables())
    for key, value in fisapt0_external_potential_energies.items():
        compare_values(
            value,
            psi4.core.variable(key),
            8,
            key,
        )
    return


if __name__ == "__main__":
    # test_saptdft_external_potential()
    test_fisapt_external_potential()
