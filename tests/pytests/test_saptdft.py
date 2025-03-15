import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values
from psi4 import core
import numpy as np
import qcelemental as qcel
from pprint import pprint as pp

hartree_to_kcalmol = constants.conversion_factor("hartree", "kcal/mol")
pytestmark = [pytest.mark.psi, pytest.mark.api]


psi4.set_memory("60 GB")


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
    mol = psi4.geometry(
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
    saptdft_external_potential = {
        "SAPT ELST ENERGY": -0.015497974872921816,
        "ELST10,R": -0.015497974872921816,
        "SAPT EXCH ENERGY": 0.014086924748375487,
        "SAPT IND ENERGY": -0.0036106959685717815,
        "SAPT DISP ENERGY": -0.002398227386784957,
        "SAPT TOTAL ENERGY": -0.007419973479903067,
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
            "SAPT_DFT_FUNCTIONAL": "pbe0",
        }
    )
    psi4.energy(
        "sapt(dft)",
        external_potentials={"C": Chargefield_C},
        molecule=mol,
    )
    pp(psi4.core.variables())
    for key, value in saptdft_external_potential.items():
        compare_values(
            value,
            psi4.core.variable(key),
            8,
            key,
        )
    return


@pytest.mark.saptdft
def test_sapthf_external_potential():
    mol = psi4.geometry(
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
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )
    psi4.energy(
        "sapt(dft)",
        external_potentials={"C": Chargefield_C},
        molecule=mol,
    )
    pp(psi4.core.variables())
    print("REF:")
    pp(fisapt0_external_potential_energies)
    key_labels = [
        ["SAPT ELST10,R ENERGY", "ELST10,R"],
        ["SAPT ELST ENERGY", "SAPT ELST ENERGY"],
        ["SAPT EXCH ENERGY", "SAPT EXCH ENERGY"],
        ["SAPT IND ENERGY", "SAPT IND ENERGY"],
        ["SAPT DISP ENERGY", "SAPT DISP ENERGY"],
        ["SAPT TOTAL ENERGY", "SAPT TOTAL ENERGY"],
    ]
    saptdft_potential_energies = {
        k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    print("CALC:")
    pp(saptdft_potential_energies)
    for k1, k2 in key_labels:
        compare_values(
            fisapt0_external_potential_energies[k1],
            psi4.core.variable(k2),
            8,
            k1,
        )
    return


@pytest.mark.saptdft
def test_sapthf_external_potential2():
    """ """
    mol = psi4.geometry(
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
        "Enuc": 36.965201909136653,
        "SAPT ELST ENERGY": -0.01581004515,
        "SAPT EXCH ENERGY": 0.01228252074,
        "SAPT IND ENERGY": -0.00356130615,
        "SAPT DISP ENERGY": -0.00218572459,
        "SAPT TOTAL ENERGY": -0.00927455515,
        "SAPT ELST10,R ENERGY": -0.01581004514947182,
    }
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
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )
    psi4.energy(
        "sapt(dft)",
        external_potentials={"C": Chargefield_C},
        molecule=mol,
    )
    pp(psi4.core.variables())
    pp(fisapt0_external_potential_energies)
    key_labels = [
        ["SAPT ELST10,R ENERGY", "ELST10,R"],
        ["SAPT ELST ENERGY", "SAPT ELST ENERGY"],
        ["SAPT EXCH ENERGY", "SAPT EXCH ENERGY"],
        ["SAPT IND ENERGY", "SAPT IND ENERGY"],
        ["SAPT DISP ENERGY", "SAPT DISP ENERGY"],
        ["SAPT TOTAL ENERGY", "SAPT TOTAL ENERGY"],
    ]
    for k1, k2 in key_labels:
        compare_values(
            fisapt0_external_potential_energies[k1],
            psi4.core.variable(k2),
            8,
            k1,
        )
    compare_values(
        fisapt0_external_potential_energies["Enuc"],
        mol.nuclear_repulsion_energy(),
        8,
        "Enuc",
    )
    return


def test_fisapt_external_potential():
    mol = psi4.geometry(
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
        molecule=mol,
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


@pytest.mark.saptdft
def test_fisapt0_external_potential_abc():
    """
    Energy test modified from ../fsapt-ext-abc/input.dat by deleting final
    water monomer because SAPT(DFT) does not support trimers at the moment.
    """
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
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) /
             psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) /
             psi_bohr2angstroms],
        ],
        "C": [
            [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
            [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
            [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        ],
    }
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
        external_potentials=external_potentials,
        molecule=mol,
    )
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]  # TEST
    # UPDATED ENERGIES FOR DIMER
    Eref = {
        "Edisp": -0.002778631330469043,
        "Eelst": -0.04933182514082859,
        "Eexch": 0.01826072035756901,
        "Eind": -0.00783603487368914,
        "Enuc": 37.565065473343004,
        "Etot": -0.04168577098741776,
    }

    Epsi = {  # TEST
        "Enuc": mol.nuclear_repulsion_energy(),  # TEST
        "Eelst": core.variable("SAPT ELST ENERGY"),  # TEST
        "Eexch": core.variable("SAPT EXCH ENERGY"),  # TEST
        "Eind": core.variable("SAPT IND ENERGY"),  # TEST
        "Edisp": core.variable("SAPT DISP ENERGY"),  # TESdfs
        "Etot": core.variable("SAPT0 TOTAL ENERGY"),  # TEST
    }  # TEST

    pp(Epsi)
    for key in keys:  # TEST
        compare_values(Eref[key], Epsi[key], 6, key)  # TEST
    return


@pytest.mark.saptdft
def test_fisapt0_external_potential_ab():
    """
    Energy test modified from ../fsapt-ext-abc/input.dat by deleting final
    water monomer because SAPT(DFT) does not support trimers at the moment.
    """
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
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) /
             psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) /
             psi_bohr2angstroms],
        ],
    }
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
        external_potentials=external_potentials,
        molecule=mol,
    )
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]  # TEST
    # UPDATED ENERGIES FOR DIMER
    Eref = {
        'Edisp': -0.002806890580921549,
        'Eelst': -0.049844754081429986,
        'Eexch': 0.018716989207357836,
        'Eind': -0.007796394758818957,
        'Enuc': 37.565065473343004,
        'Etot': -0.041731050213812654
    }

    Epsi = {  # TEST
        "Enuc": mol.nuclear_repulsion_energy(),  # TEST
        "Eelst": core.variable("SAPT ELST ENERGY"),  # TEST
        "Eexch": core.variable("SAPT EXCH ENERGY"),  # TEST
        "Eind": core.variable("SAPT IND ENERGY"),  # TEST
        "Edisp": core.variable("SAPT DISP ENERGY"),  # TESdfs
        "Etot": core.variable("SAPT0 TOTAL ENERGY"),  # TEST
    }  # TEST

    pp(Epsi)
    for key in keys:  # TEST
        compare_values(Eref[key], Epsi[key], 6, key)  # TEST
    return


@pytest.mark.saptdft
def test_sapthf_external_potential_abc():
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
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) /
             psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) /
             psi_bohr2angstroms],
        ],
        "C": [
            [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
            [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
            [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        ],
    }
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )
    psi4.energy(
        "sapt(dft)",
        external_potentials=external_potentials,
        molecule=mol,
    )
    fisapt0_external_potential_energies = {
        "SAPT DISP ENERGY": -0.002185724589094623,
        "SAPT ELST ENERGY": -0.01581004514947182,
        "SAPT ELST10,R ENERGY": -0.01581004514947182,
        "SAPT EXCH ENERGY": 0.012282520736587468,
        "SAPT IND ENERGY": -0.0035613061462424402,
        "SAPT TOTAL ENERGY": -0.009274555148221415,
    }
    # UPDATED ENERGIES FOR DIMER
    fisapt0_external_potential_energies = {
        "Edisp": -0.002778631330469043,
        "Eelst": -0.04933182514082859,
        "Eexch": 0.01826072035756901,
        "Eind": -0.00783603487368914,
        "Enuc": 37.565065473343004,
        "Etot": -0.04168577098741776,
    }
    print("REF:")
    pp(fisapt0_external_potential_energies)
    key_labels = [
        ["Eexch", "SAPT EXCH ENERGY"],
        ["Edisp", "SAPT DISP ENERGY"],
        ["Eelst", "SAPT ELST ENERGY"],
        ["Eind", "SAPT IND ENERGY"],
        ["Etot", "SAPT TOTAL ENERGY"],
    ]
    saptdft_potential_energies = {
        k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    saptdft_potential_energies["Enuc"] = mol.nuclear_repulsion_energy()
    print("CALC:")
    pp(saptdft_potential_energies)
    for k1, k2 in key_labels:
        compare_values(
            fisapt0_external_potential_energies[k1],
            saptdft_potential_energies[k1],
            8,
            k1,
        )
    return


@pytest.mark.saptdft
def test_sapthf_external_potential_ab():
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
    psi_bohr2angstroms = qcel.constants.bohr2angstroms
    external_potentials = {
        "A": [
            [0.417, np.array([-0.5496, -0.6026, 1.5720]) / psi_bohr2angstroms],
            [-0.834, np.array([-1.4545, -0.1932, 1.4677]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.9361, -0.4028, 2.2769]) / psi_bohr2angstroms],
        ],
        "B": [
            [0.417, np.array([-2.5628, -0.8269, -1.6696]) /
             psi_bohr2angstroms],
            [-0.834, np.array([-1.7899, -0.4027, -1.2768]) /
             psi_bohr2angstroms],
            [0.417, np.array([-1.8988, -0.4993, -0.3072]) /
             psi_bohr2angstroms],
        ],
        # "C": [
        #     [0.417, np.array([1.1270, 1.5527, -0.1658]) / psi_bohr2angstroms],
        #     [-0.834, np.array([1.9896, 1.0738, -0.1673]) / psi_bohr2angstroms],
        #     [0.417, np.array([2.6619, 1.7546, -0.2910]) / psi_bohr2angstroms],
        # ],
    }
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": "hf",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
        }
    )
    psi4.energy(
        "sapt(dft)",
        external_potentials=external_potentials,
        molecule=mol,
    )
    fisapt0_external_potential_energies = {
        'Eexch': 0.018716989207357836,
        'Edisp': -0.002806890580921549,
        'Eelst': -0.049844754081429986,
        'Eind': -0.007796394758818957,
        'Enuc': 37.565065473343004,
        'Etot': -0.041731050213812654
    }
    print("REF:")
    pp(fisapt0_external_potential_energies)
    key_labels = [
        ["Eexch", "SAPT EXCH ENERGY"],
        ["Edisp", "SAPT DISP ENERGY"],
        ["Eelst", "SAPT ELST ENERGY"],
        ["Eind", "SAPT IND ENERGY"],
        ["Etot", "SAPT TOTAL ENERGY"],
    ]
    saptdft_potential_energies = {
        k1: psi4.core.variable(k2) for k1, k2 in key_labels}
    saptdft_potential_energies["Enuc"] = mol.nuclear_repulsion_energy()
    print("CALC:")
    pp(saptdft_potential_energies)
    for k1, k2 in key_labels:
        compare_values(
            fisapt0_external_potential_energies[k1],
            saptdft_potential_energies[k1],
            8,
            k1,
        )
    return


if __name__ == "__main__":
    print("test_fisapt0 start")
    test_fisapt0_external_potential_ab()
    print("test_fisapt0 end")
    print("test_saptdft start")
    test_sapthf_external_potential_ab()
    print("test_saptdft end")

    # test_sapthf_external_potential()
    # test_fisapt0_external_potential_abc()
    # test_sapthf_external_potential_abc()
    # test_saptdft_external_potential()
    # test_fisapt_external_potential()
