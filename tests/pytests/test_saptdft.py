import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values
from psi4 import core
import numpy as np
import qcelemental as qcel
from pprint import pprint as pp
# TODO: use logic below to check if einsums is available, then check optoins
# optionally import einsums
# psi4.addons()
# "einsums" in psi4.addons()

# from addons import uusing
import pandas as pd

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
@pytest.mark.dftd4
def test_sapt_dft_compute_ddft_d4():
    """
    Test SAPT(DFT) for correct delta-DFT and -D4 IE terms
    """
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
  units bohr
"""
    )
    dft_functional = "pbe0"
    psi4.set_options(
        {
            "basis": "STO-3G",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "sapt_dft_grac_shift_a": 0.136,
            "sapt_dft_grac_shift_b": 0.136,
            "SAPT_DFT_FUNCTIONAL": dft_functional,
            "SAPT_DFT_DO_DDFT": True,
            "SAPT_DFT_D4_IE": True,
            "SAPT_DFT_D4_TYPE": "SUPERMOLECULAR",
            "SAPT_DFT_DO_DISP": True,
        }
    )
    dft_IE = (
        psi4.energy(dft_functional, bsse_type="CP", molecule=mol_dimer)
        * hartree_to_kcalmol
    )
    # many-body energies
    # dimer    :-150.5070040376297698
    # monomer A: -75.2534797048230359
    # monomer B: -75.2469427300573841

    psi4.energy("SAPT(DFT)")
    # sapt(dft) DFT energies
    # dimer    :-150.5070038696833024
    # monomer A: -75.2469415558837937
    # monomer B: -75.2534791716759486

    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DFT_MONA = psi4.core.variable("DFT MONOMER A ENERGY")
    DFT_MONB = psi4.core.variable("DFT MONOMER B ENERGY")
    # DFT_MONA = psi4.core.variable("DFT MONOMERA")
    # DFT_MONB = psi4.core.variable("DFT MONOMERB")
    DFT_DIMER = psi4.core.variable("DFT DIMER ENERGY")
    # DFT_DIMER = psi4.core.variable("DFT DIMER")

    print(f"{DFT_DIMER=}\n{DFT_MONA=}\n{DFT_MONB=}")
    DFT_IE = (DFT_DIMER - DFT_MONA - DFT_MONB) * hartree_to_kcalmol
    print(f"bsse: {dft_IE=}")
    print(f"SAPT: {DFT_IE=}")
    assert compare_values(dft_IE, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018077857232, DFT_IE, 7, "DFT IE")

    # get dft_IE from SAPT(DFT) Delta DFT term to back-calculate
    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DELTA_HF = psi4.core.variable("SAPT(DFT) DELTA HF")
    DDFT = psi4.core.variable("SAPT(DFT) DELTA DFT")
    DFT_IE_from_dDFT = (DDFT + ELST + EXCH + IND - DELTA_HF) * hartree_to_kcalmol
    assert compare_values(DFT_IE_from_dDFT, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018048599092, DFT_IE, 7, "DFT IE")
    print(f"dDFT: {DFT_IE_from_dDFT=}")


@pytest.mark.saptdft
@pytest.mark.dftd4
@pytest.mark.parametrize(
    "method, d4_type, expected_disp",
    [
        ("SAPT(DFT)-D4(I)", "intermolecular", -0.0041772889),
        ("SAPT(DFT)-D4(S)", "supermolecular", -0.0036056912),
        ("DFT-D4(SAPT)", "gd4_supermolecular", -0.005731715146359108),
    ],
)
def test_saptdft_disp_methods(method, d4_type, expected_disp):
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
units bohr
"""
    )
    psi4.set_options(
        {
            "basis": "STO-3G",
            "sapt_dft_grac_shift_a": 0.136,
            "sapt_dft_grac_shift_b": 0.136,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy(method)
    vars = psi4.core.variables()
    DISP = vars["SAPT DISP ENERGY"]
    assert compare_values(expected_disp, DISP, 8, f"{method} DISP")


@pytest.mark.saptdft
@pytest.mark.dftd4
def test_sapt_dft_compute_ddft_d4_diskdf():
    """
    Test SAPT(DFT) for correct delta-DFT and -D4 IE terms
    """
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
  units bohr
"""
    )
    dft_functional = "pbe0"
    psi4.set_options(
        {
            "basis": "STO-3G",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "scf_type": "disk_df",
            "sapt_dft_grac_shift_a": 0.136,
            "sapt_dft_grac_shift_b": 0.136,
            "SAPT_DFT_FUNCTIONAL": dft_functional,
            "SAPT_DFT_DO_DDFT": True,
            "SAPT_DFT_D4_IE": True,
            "SAPT_DFT_DO_DISP": True,
        }
    )
    dft_IE = (
        psi4.energy(dft_functional, bsse_type="CP", molecule=mol_dimer)
        * hartree_to_kcalmol
    )
    # many-body energies
    # dimer    :-150.5070040376297698
    # monomer A: -75.2534797048230359
    # monomer B: -75.2469427300573841

    psi4.energy("SAPT(DFT)")
    # sapt(dft) DFT energies
    # dimer    :-150.5070038696833024
    # monomer A: -75.2469415558837937
    # monomer B: -75.2534791716759486

    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DFT_MONA = psi4.core.variable("DFT MONOMER A ENERGY")
    DFT_MONB = psi4.core.variable("DFT MONOMER B ENERGY")
    # DFT_MONA = psi4.core.variable("DFT MONOMERA")
    # DFT_MONB = psi4.core.variable("DFT MONOMERB")
    DFT_DIMER = psi4.core.variable("DFT DIMER ENERGY")
    # DFT_DIMER = psi4.core.variable("DFT DIMER")

    print(f"{DFT_DIMER=}\n{DFT_MONA=}\n{DFT_MONB=}")
    DFT_IE = (DFT_DIMER - DFT_MONA - DFT_MONB) * hartree_to_kcalmol
    print(f"bsse: {dft_IE=}")
    print(f"SAPT: {DFT_IE=}")
    assert compare_values(dft_IE, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018077857232, DFT_IE, 7, "DFT IE")

    # get dft_IE from SAPT(DFT) Delta DFT term to back-calculate
    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DELTA_HF = psi4.core.variable("SAPT(DFT) DELTA HF")
    DDFT = psi4.core.variable("SAPT(DFT) DELTA DFT")
    DFT_IE_from_dDFT = (DDFT + ELST + EXCH + IND - DELTA_HF) * hartree_to_kcalmol
    assert compare_values(DFT_IE_from_dDFT, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018048599092, DFT_IE, 7, "DFT IE")
    print(f"dDFT: {DFT_IE_from_dDFT=}")


@pytest.mark.saptdft
@pytest.mark.dftd4
def test_sapt_dft_diskdf():
    """
    Test SAPT(DFT) for correct delta-DFT and -D4 IE terms
    """
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
  units bohr
"""
    )
    psi4.set_memory("280 MB")
    psi4.core.set_num_threads(2)
    dft_functional = "pbe0"
    psi4.set_options(
        {
            "basis": "aug-cc-pvdz",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "scf_type": "disk_df",
            "sapt_dft_grac_shift_a": 0.136,
            "sapt_dft_grac_shift_b": 0.136,
            "SAPT_DFT_FUNCTIONAL": dft_functional,
            "SAPT_DFT_DO_DDFT": True,
            "SAPT_DFT_D4_IE": True,
            "SAPT_DFT_DO_DISP": True,
        }
    )
    dft_IE = (
        psi4.energy(dft_functional, bsse_type="CP", molecule=mol_dimer)
        * hartree_to_kcalmol
    )
    # many-body energies
    # dimer    :-150.5070040376297698
    # monomer A: -75.2534797048230359
    # monomer B: -75.2469427300573841

    psi4.energy("SAPT(DFT)")
    # sapt(dft) DFT energies
    # dimer    :-150.5070038696833024
    # monomer A: -75.2469415558837937
    # monomer B: -75.2534791716759486

    DFT_MONA = psi4.core.variable("DFT MONOMER A ENERGY")
    DFT_MONB = psi4.core.variable("DFT MONOMER B ENERGY")
    DFT_DIMER = psi4.core.variable("DFT DIMER ENERGY")

    print(f"{DFT_DIMER=}\n{DFT_MONA=}\n{DFT_MONB=}")
    DFT_IE = (DFT_DIMER - DFT_MONA - DFT_MONB) * hartree_to_kcalmol
    print(f"bsse: {dft_IE=}")
    print(f"SAPT: {DFT_IE=}")
    assert compare_values(dft_IE, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.91463144301158, DFT_IE, 7, "DFT IE")


@pytest.mark.saptdft
@pytest.mark.dftd4
def test_sapt_dft_compute_ddft_d4_auto_grac():
    """
    Test SAPT(DFT) for correct delta-DFT and -D4 IE terms
    """
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
  units bohr
"""
    )
    dft_functional = "pbe0"
    psi4.set_options(
        {
            "basis": "STO-3G",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "sapt_dft_grac_compute": "SINGLE",
            "SAPT_DFT_FUNCTIONAL": dft_functional,
            "SAPT_DFT_DO_DDFT": True,
            "SAPT_DFT_D4_IE": True,
            "SAPT_DFT_DO_DISP": True,
        }
    )
    dft_IE = (
        psi4.energy(dft_functional, bsse_type="CP", molecule=mol_dimer)
        * hartree_to_kcalmol
    )
    psi4.energy("SAPT(DFT)")
    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DFT_MONA = psi4.core.variable("DFT MONOMER A ENERGY")
    DFT_MONB = psi4.core.variable("DFT MONOMER B ENERGY")
    DFT_DIMER = psi4.core.variable("DFT DIMER ENERGY")

    print(f"{DFT_DIMER=}\n{DFT_MONA=}\n{DFT_MONB=}")
    DFT_IE = (DFT_DIMER - DFT_MONA - DFT_MONB) * hartree_to_kcalmol
    print(f"bsse: {dft_IE=}")
    print(f"SAPT: {DFT_IE=}")
    assert compare_values(dft_IE, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018077857232, DFT_IE, 7, "DFT IE")

    # get dft_IE from SAPT(DFT) Delta DFT term to back-calculate
    ELST = psi4.core.variable("SAPT ELST ENERGY")
    EXCH = psi4.core.variable("SAPT EXCH ENERGY")
    IND = psi4.core.variable("SAPT IND ENERGY")
    DELTA_HF = psi4.core.variable("SAPT(DFT) DELTA HF")
    DDFT = psi4.core.variable("SAPT(DFT) DELTA DFT")
    DFT_IE_from_dDFT = (DDFT + ELST + EXCH + IND - DELTA_HF) * hartree_to_kcalmol
    assert compare_values(DFT_IE_from_dDFT, DFT_IE, 7, "DFT IE")
    assert compare_values(-4.130018048599092, DFT_IE, 7, "DFT IE")
    print(f"dDFT: {DFT_IE_from_dDFT=}")
    compare_values(
        #  STO-3G target
        0.1981702737,
        # aug-cc-pvdz target, 0.1307 (using experimental IP from CCCBDB)
        # 0.13053068183319516,
        psi4.core.variable("SAPT DFT GRAC SHIFT A"),
        8,
        "SAPT DFT GRAC SHIFT A",
    )
    compare_values(
        #  STO-3G target
        0.1983742234,
        # aug-cc-pvdz target, 0.1307 (using experimental IP from CCCBDB)
        # 0.13063798506967816,
        psi4.core.variable("SAPT DFT GRAC SHIFT B"),
        8,
        "SAPT DFT GRAC SHIFT B",
    )
    print(f"{psi4.core.variable('SAPT DFT GRAC SHIFT A')=}")
    print(f"{psi4.core.variable('SAPT DFT GRAC SHIFT B')=}")


@pytest.mark.saptdft
@pytest.mark.dftd4
def test_dftd4():
    """
    Tests SAPT(DFT) module for computing DFT-D4 SAPT-decomposition of energy terms
    """
    mol_dimer = psi4.geometry(
        """
  O -2.930978458   -0.216411437    0.000000000
  H -3.655219777    1.440921844    0.000000000
  H -1.133225297    0.076934530    0.000000000
   --
  O  2.552311356    0.210645882    0.000000000
  H  3.175492012   -0.706268134   -1.433472544
  H  3.175492012   -0.706268134    1.433472544
  units bohr
"""
    )
    psi4.set_options(
        {
            "basis": "STO-3G",
            "sapt_dft_grac_shift_a": 0.136,
            "sapt_dft_grac_shift_b": 0.136,
            "SAPT_DFT_FUNCTIONAL": "pbe0",
            "SAPT_DFT_DO_DDFT": True,
            "SAPT_DFT_D4_IE": True,
            "SAPT_DFT_D4_TYPE": "GD4_SUPERMOLECULAR",
            "SAPT_DFT_DO_DISP": False,
        }
    )
    psi4.energy("SAPT(DFT)", molecule=mol_dimer)
    DISP = psi4.core.variable("SAPT DISP ENERGY")
    assert compare_values(-0.005731715146359108, DISP, 8, "DFT-D4 DISP")


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
            "neutral_water_dimer",
            None,
        ),
        (
            "SINGLE",
            0.1307,
            0.19830016,
            0.1307,
            None,
            "neutral_water_dimer",
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
            0.0924377708,
            0.1306379832,
            None,
            None,
            "hydroxide",
            "aug-cc-pvdz",
        ),
    ],
)
def test_saptdft_auto_grac(
    SAPT_DFT_GRAC_COMPUTE, refA, refB, gracA, gracB, geometry, grac_basis
):
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
    assert compare_values(
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
    return


@pytest.mark.extern
@pytest.mark.saptdft
@pytest.mark.parametrize(
    "test_id, external_pot_keys, energy_method, expected_values, precision",
    [
        (
            "abc",
            ["A", "B", "C"],
            "sapt(dft)",
            {
                "Edisp": -0.002778631330469043,
                "Eelst": -0.04933182514082859,
                "Eexch": 0.01826072035756901,
                "Eind": -0.00783603487368914,
                "Enuc": 37.565065473343004,
                "Etot": -0.04168577098741776,
            },
            8,
        ),
        (
            "abc",
            ["A", "B", "C"],
            "fisapt0",
            {
                "Edisp": -0.002778631330469043,
                "Eelst": -0.04933182514082859,
                "Eexch": 0.01826072035756901,
                "Eind": -0.00783603487368914,
                "Enuc": 37.565065473343004,
                "Etot": -0.04168577098741776,
            },
            8,
        ),
        (
            "ab",
            ["A", "B"],
            "sapt(dft)",
            {
                "Eexch": 0.018716989207357836,
                "Edisp": -0.002806890580921549,
                "Eelst": -0.049844754081429986,
                "Eind": -0.007796394758818957,
                "Enuc": 37.565065473343004,
                "Etot": -0.041731050213812654,
            },
            7,
        ),
        (
            "a",
            ["A"],
            "sapt(dft)",
            {
                "Edisp": -0.002816214841039168,
                "Eelst": -0.029279397864449663,
                "Eexch": 0.01889146242465646,
                "Eind": -0.006356586601442907,
                "Enuc": 37.565065473343004,
                "Etot": -0.019560736882275276,
            },
            8,
        ),
        (
            "b",
            ["B"],
            "sapt(dft)",
            {
                "Edisp": -0.0027835628235023946,
                "Eelst": -0.024996838889304485,
                "Eexch": 0.018658252478543982,
                "Eind": -0.005940905006556488,
                "Enuc": 37.565065473343004,
                "Etot": -0.015063054240819385,
            },
            7,
        ),
        (
            "c",
            ["C"],
            "sapt(dft)",
            {
                "Edisp": -0.0027641749042896808,
                "Eelst": -0.013585309461394246,
                "Eexch": 0.01839802168509688,
                "Eind": -0.004537242936115817,
                "Enuc": 37.565065473343004,
                "Etot": -0.002488705616702865,
            },
            8,
        ),
    ],
)
def test_sapthf_external_potential(
    test_id, external_pot_keys, energy_method, expected_values, precision
):
    """
    Parameterized test for SAPT with external potentials.

    Parameters:
    -----------
    test_id : str
        Identifier for the test case
    external_pot_keys : list
        List of keys from the full external_potentials dict to include in this test
    energy_method : str
        Energy method to use (sapt(dft) or fisapt0)
    expected_values : dict
        Expected energy values to compare against
    precision : int
        Decimal precision for comparison
    """
    # Define the molecule geometry (common to all tests)
    mol = psi4.geometry(
        """
0 1
H 0.0290 -1.1199 -1.5243
O 0.9481 -1.3990 -1.3587
H 1.4371 -0.5588 -1.3099
--
H 1.0088 -1.5240 0.5086
O 1.0209 -1.1732 1.4270
H 1.5864 -0.3901 1.3101
symmetry c1
no_reorient
no_com
        """
    )

    # All potential definitions
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    all_external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) / psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) / psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) / psi_bohr2angstroms],
        ],
        "C": [
            [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
            [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
            [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        ],
    }

    # Select only the potentials needed for this test
    external_potentials = {
        key: all_external_potentials[key] for key in external_pot_keys
    }

    # Set common options
    psi4.set_options(
        {
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )

    # Run the energy calculation
    psi4.energy(
        energy_method,
        external_potentials=external_potentials,
        molecule=mol,
    )

    # Print reference values
    print(f"TEST ID: {test_id}")
    print("REF:")
    import pprint

    pprint.pprint(expected_values)

    # Collect calculated values
    key_labels = [
        ["Eexch", "SAPT EXCH ENERGY"],
        ["Edisp", "SAPT DISP ENERGY"],
        ["Eelst", "SAPT ELST ENERGY"],
        ["Eind", "SAPT IND ENERGY"],
        ["Etot", "SAPT TOTAL ENERGY"],
    ]
    calculated_energies = {k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    calculated_energies["Enuc"] = mol.nuclear_repulsion_energy()

    # Print calculated values
    print("CALC:")
    pprint.pprint(calculated_energies)

    # Compare values
    for k1, k2 in key_labels:
        compare_values(
            expected_values[k1],
            calculated_energies[k1],
            precision,
            k1,
        )

    # Also check nuclear repulsion energy
    compare_values(
        expected_values["Enuc"],
        calculated_energies["Enuc"],
        precision,
        "Enuc",
    )

    return


@pytest.mark.extern
@pytest.mark.saptdft
@pytest.mark.parametrize(
    "test_id, external_pot_keys, energy_method, functional, expected_values, precision",
    [
        (
            "abc",
            ["A", "B", "C"],
            "sapt(dft)",
            "pbe0",
            {
                "Edisp": -0.0030693649240918575,
                "Eelst": -0.04826201906029539,
                "Eexch": 0.021079055937265008,
                "Eind": -0.007991120173136683,
                "Enuc": 37.565065473343004,
                "Etot": -0.03824344798868797,
            },
            8,
        ),
        (
            "abc",
            ["A", "B", "C"],
            "fisapt0",
            "hf",
            {
                "Edisp": -0.0027786313304796356,
                "Eelst": -0.04933182514519132,
                "Eexch": 0.01826072035780734,
                "Eind": -0.007836034873458513,
                "Enuc": 37.565065473343004,
                "Etot": -0.041685770991322126,
            },
            8,
        ),
        (
            "ab",
            ["A", "B"],
            "sapt(dft)",
            "pbe0",
            {
                "Edisp": -0.0031046737482147164,
                "Eelst": -0.048898649876071974,
                "Eexch": 0.021726965401698814,
                "Eind": -0.007951966773147297,
                "Enuc": 37.565065473343004,
                "Etot": -0.03822832499573517,
            },
            8,
        ),
        (
            "a",
            ["A"],
            "sapt(dft)",
            "pbe0",
            {
                "Edisp": -0.0031075464270878045,
                "Eelst": -0.028802617754731585,
                "Eexch": 0.021975641266292034,
                "Eind": -0.006465625562838023,
                "Enuc": 37.565065473343004,
                "Etot": -0.01640014847836538,
            },
            8,
        ),
        (
            "b",
            ["B"],
            "sapt(dft)",
            "pbe0",
            {
                "Edisp": -0.0031243901884552955,
                "Eelst": -0.027251596918624443,
                "Eexch": 0.022466107772087913,
                "Eind": -0.0060653648137171425,
                "Enuc": 37.565065473343004,
                "Etot": -0.013975244148708968,
            },
            8,
        ),
        (
            "c",
            ["C"],
            "sapt(dft)",
            "pbe0",
            {
                "Edisp": -0.0030359862882380568,
                "Eelst": -0.013872053096907422,
                "Eexch": 0.021170116566131776,
                "Eind": -0.004613611052900464,
                "Enuc": 37.565065473343004,
                "Etot": -0.0003515338719141661,
            },
            8,
        ),
    ],
)
def test_saptdft_external_potential(
    test_id, external_pot_keys, energy_method, functional, expected_values, precision
):
    """
    Parameterized test for SAPT with external potentials.

    Parameters:
    -----------
    test_id : str
        Identifier for the test case
    external_pot_keys : list
        List of keys from the full external_potentials dict to include in this test
    energy_method : str
        Energy method to use (sapt(dft) or fisapt0)
    expected_values : dict
        Expected energy values to compare against
    precision : int
        Decimal precision for comparison
    """
    # Define the molecule geometry (common to all tests)
    mol = psi4.geometry(
        """
0 1
H 0.0290 -1.1199 -1.5243
O 0.9481 -1.3990 -1.3587
H 1.4371 -0.5588 -1.3099
--
H 1.0088 -1.5240 0.5086
O 1.0209 -1.1732 1.4270
H 1.5864 -0.3901 1.3101
symmetry c1
no_reorient
no_com
        """
    )

    # All potential definitions
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    all_external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) / psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) / psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) / psi_bohr2angstroms],
        ],
        "C": [
            [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
            [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
            [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        ],
    }

    # Select only the potentials needed for this test
    external_potentials = {
        key: all_external_potentials[key] for key in external_pot_keys
    }

    # Set common options
    psi4.set_options(
        {
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": functional,
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
            "sapt_dft_grac_shift_a": 0.1307,
            "sapt_dft_grac_shift_b": 0.1307,
        }
    )

    # Run the energy calculation
    psi4.energy(
        energy_method,
        external_potentials=external_potentials,
        molecule=mol,
    )

    # Print reference values
    print(f"TEST ID: {test_id}")
    print("REF:")
    import pprint

    pprint.pprint(expected_values)

    # Collect calculated values
    key_labels = [
        ["Eexch", "SAPT EXCH ENERGY"],
        ["Edisp", "SAPT DISP ENERGY"],
        ["Eelst", "SAPT ELST ENERGY"],
        ["Eind", "SAPT IND ENERGY"],
        ["Etot", "SAPT TOTAL ENERGY"],
    ]
    calculated_energies = {k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    calculated_energies["Enuc"] = mol.nuclear_repulsion_energy()

    # Print calculated values
    print("CALC:")
    pprint.pprint(calculated_energies)

    # Compare values
    for k1, k2 in key_labels:
        compare_values(
            expected_values[k1],
            calculated_energies[k1],
            precision,
            k1,
        )

    # Also check nuclear repulsion energy
    compare_values(
        expected_values["Enuc"],
        calculated_energies["Enuc"],
        precision,
        "Enuc",
    )

    return


@pytest.mark.extern
@pytest.mark.saptdft
def test_fisapt0_sapthf_external_potential():
    mol = psi4.geometry(
        """
0 1
H 0.0290 -1.1199 -1.5243
O 0.9481 -1.3990 -1.3587
H 1.4371 -0.5588 -1.3099
--
H 1.0088 -1.5240 0.5086
O 1.0209 -1.1732 1.4270
H 1.5864 -0.3901 1.3101
symmetry c1
no_reorient
no_com
        """
    )
    # All potential definitions
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) / psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) / psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) / psi_bohr2angstroms],
        ],
        "C": [
            [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
            [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
            [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        ],
    }

    # Set common options
    psi4.set_options(
        {
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )

    # Run the FISAPT0 energy calculation
    psi4.energy(
        "fisapt0",
        external_potentials=external_potentials,
        molecule=mol,
    )

    # Collect calculated values
    key_labels = [
        ["Eexch", "SAPT EXCH ENERGY"],
        ["Edisp", "SAPT DISP ENERGY"],
        ["Eelst", "SAPT ELST ENERGY"],
        ["Eind", "SAPT IND ENERGY"],
        ["Etot", "SAPT TOTAL ENERGY"],
    ]
    calculated_fisapt0_energies = {k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    calculated_fisapt0_energies["Enuc"] = mol.nuclear_repulsion_energy()

    # Run the SAPT(HF) energy calculation
    psi4.energy(
        "sapt(dft)",
        external_potentials=external_potentials,
        molecule=mol,
    )

    calculated_sapthf_energies = {k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    calculated_sapthf_energies["Enuc"] = mol.nuclear_repulsion_energy()

    print("FISAPT0:")
    pp(calculated_fisapt0_energies)
    print("SAPT(HF):")
    pp(calculated_sapthf_energies)

    # Compare values
    for k1, k2 in key_labels:
        compare_values(
            calculated_fisapt0_energies[k1],
            calculated_sapthf_energies[k1],
            8,
            k1,
        )

    # Also check nuclear repulsion energy
    compare_values(
        calculated_fisapt0_energies["Enuc"],
        calculated_sapthf_energies["Enuc"],
        8,
        "Enuc",
    )


def test_qcng_embedded_saptdft():
    import qcengine as qcng

    atin = {
        "driver": "energy",
        "extras": {},
        "id": "53500",
        # "keywords": {"SAPT_DFT_GRAC_COMPUTE": "SINGLE"},
        "keywords": {
            # "SAPT_DFT_GRAC_COMPUTE": "SINGLE"
            "SAPT_DFT_GRAC_SHIFT_A": 0.1307,
            "SAPT_DFT_GRAC_SHIFT_B": 0.1307,
        },
        "model": {"basis": "sto-3g", "method": "sapt(dft)"},
        "molecule": {
            "extras": {},
            "fix_com": True,
            "fix_orientation": True,
            "fragment_charges": [0.0, 0.0],
            "fragment_multiplicities": [1, 1],
            "fragments": [[0, 1, 2], [3, 4, 5]],
            "geometry": [
                -1.3269582284372,
                -0.1059385303631,
                0.0187881522475,
                -1.9316652406588,
                1.6001743176504,
                -0.0217105229937,
                0.486644278717,
                0.0795980914346,
                0.009862478759,
                4.287563293074,
                0.0497755770069,
                0.0009600356738,
                4.9992749983517,
                -0.7786426865932,
                1.4487252956894,
                4.9910408984847,
                -0.8501365231326,
                -1.4076465463372,
            ],
            "id": 53204,
            "identifiers": {
                "molecular_formula": "H4O2",
                "molecule_hash": "fa17e4abd32ccf4b5a3db8da6d94507f9ef6941f",
            },
            "molecular_charge": 0.0,
            "molecular_multiplicity": 1,
            "name": "H4O2",
            "provenance": {
                "creator": "QCElemental",
                "routine": "qcelemental.molparse.from_string",
                "version": "0.29.0",
            },
            "schema_name": "qcschema_molecule",
            "schema_version": 2,
            "symbols": ["O", "H", "H", "O", "H", "H"],
            "validated": True,
        },
        "protocols": {"stdout": True},
        "provenance": {
            "creator": "QCElemental",
            "routine": "qcelemental.models.results",
            "version": "0.29.0",
        },
        "schema_name": "qcschema_input",
        "schema_version": 1,
    }

    print("run_qcschema")
    ret_1 = psi4.schema_wrapper.run_qcschema(
        atin,
    )
    print(ret_1)
    assert compare_values(
        -0.00191336,
        ret_1.extras["qcvars"]["SAPT TOTAL ENERGY"],
        4,
        "SAPT(DFT) TOTAL run_qschema",
    )
    print("qcng")
    ret_2 = qcng.compute(
        atin,
        "psi4",
        raise_error=True,
    )
    print(ret_2)
    assert compare_values(
        -0.00191336,
        ret_2.extras["qcvars"]["SAPT TOTAL ENERGY"],
        4,
        "SAPT(DFT) TOTAL qcng",
    )
    return


def test_charge_field_inputs():
    dimer = psi4.geometry("""
    0 1
    8   60.268880784   0.026340101   0.000508029
    1   60.645502399   -0.412039965   0.766632411
    1   60.641145101   -0.449872874   -0.744894473
    --
    0 1
    8   50.268880784   0.026340101   0.000508029
    1   50.645502399   -0.412039965   0.766632411
    1   50.641145101   -0.449872874   -0.744894473
    units angstrom
    symmetry c1
    no_com
    no_reorient
    """)

    Chargefield = np.array(
        [
            0.5972,
            4.802,
            -1.38,
            23.692,
            -0.5679,
            4.723,
            -0.45,
            22.895,
            -0.3662,
            4.229,
            -1.231,
            25.089,
            0.1123,
            3.188,
            -0.936,
            25.072,
            0.009333,
            11.667,
            9.479,
            19.280,
            0.009333,
            10.878,
            8.794,
            19.002,
            0.283950,
            12.371,
            7.698,
            26.868,
            0.283950,
            13.479,
            7.947,
            26.797,
        ]
    ).reshape((-1, 4))
    Chargefield[:, [1, 2, 3]] /= qcel.constants.bohr2angstroms

    psi4.set_options(
        {
            "basis": "jun-cc-pv(D+d)z",
            "freeze_core": "True",
            "scf_type": "df",
            "mp2_type": "df",
            "fisapt_do_fsapt": "false",
            "SAPT_DFT_GRAC_COMPUTE": "single",
        }
    )

    e_a = psi4.energy(
        "sapt(dft)", external_potentials={"A": Chargefield}, molecule=dimer
    )
    e_A = psi4.energy(
        "sapt(dft)", external_potentials={"a": Chargefield}, molecule=dimer
    )
    assert compare_values(e_A, e_a, 7, "e_A==e_a")
    e_b = psi4.energy(
        "sapt(dft)", external_potentials={"b": Chargefield}, molecule=dimer
    )
    e_B = psi4.energy(
        "sapt(dft)", external_potentials={"B": Chargefield}, molecule=dimer
    )
    assert compare_values(e_B, e_b, 7, "e_A==e_a")


def test_einsum_terms():
    """
    built from sapt-dft1 ctest
    """
    Eref_nh = {
        "SAPT ELST ENERGY": -0.22987897,  # mEh
        "SAPT EXCH ENERGY": 0.59560159,  # mEh
        "SAPT IND ENERGY": -0.00010341,  # mEh
        "SAPT DISP ENERGY": 0.00000574,  # mEh
        "CURRENT ENERGY": 0.36562495,  # mEh
    }  # TEST
    mol = psi4.geometry("""
  Ne
  --
  Ne 1 4.5
  units bohr
    """)
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "sapt_dft_grac_shift_a": 0.203293,
            "sapt_dft_grac_shift_b": 0.203293,
            "SAPT_DFT_DO_DHF": False,
            "SAPT_DFT_DO_HYBRID": False,
            "SAPT_DFT_USE_EINSUMS": True,
            "SAPT_DFT_EXCH_DISP_SCALE_SCHEME": "None",
        }
    )
    psi4.energy("sapt(dft)", molecule=mol)
    for k, v in Eref_nh.items():  # TEST
        ref = v
        assert compare_values(
            ref, psi4.variable(k) * 1000, 8, "!hyb, xd=none, !dHF: " + k
        )


def test_saptdft_inf():
    # implement this test
    Eref = {
        "Exch-Ind20,u (A<-B)": 0.00749283,  # TEST
        "Exch-Ind20,u (A->B)": 0.07674835,  # TEST
        "Exch-Ind20,u": 0.08424119,  # TEST
        "Exch-Ind20,u (A<-B) (S^inf)": 0.00749701,  # TEST
        "Exch-Ind20,u (A->B) (S^inf)": 0.07677809,  # TEST
        "Exch-Ind20,u (S^inf)": 0.08427511,  # TEST
        "Exch-Ind20,r (A<-B)": 0.00838528,  # TEST
        "Exch-Ind20,r (A->B)": 0.08663184,  # TEST
        "Exch-Ind20,r": 0.09501712,  # TEST
        "Exch-Ind20,r (A<-B) (S^inf)": 0.00839011,  # TEST
        "Exch-Ind20,r (A->B) (S^inf)": 0.08666563,  # TEST
        "Exch-Ind20,r (S^inf)": 0.09505574,
    }  # TEST
    mol = psi4.geometry("""
  Ne
  --
  Ar 1 6.5
  units bohr
    """)
    psi4.set_options(
        {
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "DO_IND_EXCH_SINF": True,
            "SAPT_DFT_FUNCTIONAL": "HF",
        }
    )
    psi4.energy("sapt(dft)", molecule=mol)
    for k, v in Eref.items():  # TEST
        ref = v / 1000.0
        assert compare_values(
            ref, psi4.variable(k), 6, "sapt(dft) inf exch-ind20: " + k
        )
    psi4.core.clean()
    mol = psi4.geometry("""
0 1
O -0.064997  0.000000  1.202015
H  0.817736  0.000000  1.573527
H  0.070857  0.000000  0.246803
--
0 1
O  0.060546  0.000000 -1.130293
H -0.408968 -0.760178 -1.479405
H -0.408968  0.760178 -1.479405

units angstrom
no_reorient
symmetry c1
""")
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "df_basis_scf": "aug-cc-pvtz-jkfit",
            "df_basis_mp2": "aug-cc-pvtz-ri",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "scf_type": "mem_df",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "DO_DISP_EXCH_SINF": True,
        }
    )
    psi4.energy("sapt(dft)", molecule=mol)
    assert compare_values(
        0.00324766,
        psi4.variable("sapt exch-disp20(s^inf) energy"),
        7,
        "SAPT EXCH-DISP20(S^inf) ENERGY",
    )


if __name__ == "__main__":
    psi4.set_memory("32 GB")
    psi4.set_num_threads(12)
    # test_einsum_terms()

    # pytest this file
    # pytest.main(
    #     [
    #         __file__,
    #         "-v",
    #         "-s",
    #         "-k=test_saptdft_disp_methods",
    #         "--disable-warnings",
    #         # "--maxfail=1",
    #     ]
    # )
    # test_sapt_dft_compute_ddft_d4_diskdf()
    test_sapt_dft_diskdf()
    # test_dftd4()
    # test_saptdft_disp_methods("DFT-D4(SAPT)", "gd4_supermolecular", -0.005731715146359108)
    # test_saptdft_disp_methods(
    #     "SAPT(DFT)-D4(I)", "intermolecular", -0.0041772889
    #     "SAPT(DFT)-D4(S)", "supermolecular", -0.0036056912
    # )
