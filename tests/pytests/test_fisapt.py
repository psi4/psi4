import pytest
import psi4
from psi4 import compare_values, variable
from addons import uusing
import numpy as np

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


# TODO: create test without pandas so it tests
@pytest.mark.fsapt
def test_fsapt_psivars_dict():
    """
    fsapt-psivars: calling fsapt_analysis with psi4 variables after running an
    fisapt0 calcluation requires the user to pass the molecule object
    """
    mol = psi4.geometry(
        """0 1
C 0.00000000 0.00000000 0.00000000
H 1.09000000 0.00000000 0.00000000
H -0.36333333 0.83908239 0.59332085
H -0.36333333 0.09428973 -1.02332709
H -0.36333333 -0.93337212 0.43000624
--
0 1
C 6.44536662 -0.26509169 -0.00000000
H 7.53536662 -0.26509169 -0.00000000
H 6.08203329 0.57399070 0.59332085
H 6.08203329 -0.17080196 -1.02332709
H 6.08203329 -1.19846381 0.43000624
symmetry c1
no_reorient
no_com"""
    )
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
        }
    )
    psi4.core.be_quiet()
    psi4.energy("fisapt0")
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Enuc": 35.07529824960602,
        "Eelst": -3.8035153870907834e-06,
        "Eexch": 1.7912112685446533e-07,
        "Eind": -3.833795151474493e-08,
        "Edisp": -3.288568662589654e-05,
        "Etot": -3.6548418837647605e-05,
    }
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": variable("SAPT ELST ENERGY"),
        "Eexch": variable("SAPT EXCH ENERGY"),
        "Eind": variable("SAPT IND ENERGY"),
        "Edisp": variable("SAPT DISP ENERGY"),
        "Etot": variable("SAPT0 TOTAL ENERGY"),
    }

    for key in keys:
        compare_values(Eref[key], Epsi[key], 6, key)
    fEnergies = psi4.fsapt_analysis(
        molecule=mol,
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "MethylA": [1, 2, 3, 4, 5],
        },
        fragments_b={
            "MethylB": [6, 7, 8, 9, 10],
        },
    )
    fEnergies.pop("Frag1")
    fEnergies.pop("Frag2")
    print(fEnergies)
    fEref = {
        "fEelst": -0.002,
        "fEexch": 0.000,
        "fEindAB": -0.000,
        "fEindBA": -0.000,
        "fEdisp": -0.021,
        "fEedisp": 0.000,
        "fEtot": -0.023,
    }

    # python iterate over zip dictionary keys and values
    for key1, key2 in zip(fEref.keys(), fEnergies.keys()):
        compare_values(fEref[key1], fEnergies[key2][0], 2, key1)


@pytest.mark.fsapt
def test_fsapt_external_potentials():
    """
    fsapt-psivars: calling fsapt_analysis with psi4 variables after running an
    fisapt0 calcluation requires the user to pass the molecule object
    """
    mol = psi4.geometry(
        """
H 0.0290 -1.1199 -1.5243
O 0.9481 -1.3990 -1.3587
H 1.4371 -0.5588 -1.3099
--
H 1.0088 -1.5240 0.5086
O 1.0209 -1.1732 1.4270
H 1.5864 -0.3901 1.3101
--
H -1.0231 1.6243 -0.8743
O -0.5806 2.0297 -0.1111
H -0.9480 1.5096 0.6281
symmetry c1
no_reorient
no_com
"""
    )
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
        }
    )
    psi_bohr2angstroms = psi4.constants.bohr2angstroms
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
    psi4.core.be_quiet()
    psi4.energy("fisapt0", external_potentials=external_potentials)
    print(psi4.core.variables())
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Enuc": 74.2330370461897,
        "Eelst": -0.04919037863747235,
        "Eexch": 0.018239207303845935,
        "Eind": -0.007969545823122322,
        "Edisp": -0.002794948165605119,
        "Etot": -0.04171566532235386,
    }
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": variable("SAPT ELST ENERGY"),
        "Eexch": variable("SAPT EXCH ENERGY"),
        "Eind": variable("SAPT IND ENERGY"),
        "Edisp": variable("SAPT DISP ENERGY"),
        "Etot": variable("SAPT0 TOTAL ENERGY"),
    }

    for key in keys:
        compare_values(Eref[key], Epsi[key], 6, key)
    fEnergies = psi4.fsapt_analysis(
        molecule=mol,
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "w1": [1, 2, 3],
        },
        fragments_b={
            "w3": [4, 5, 6],
        },
    )
    fEnergies.pop("Frag1")
    fEnergies.pop("Frag2")
    fEref = {
        "fEelst": -30.867,
        "fEexch": 11.445,
        "fEindAB": -3.138,
        "fEindBA": -1.863,
        "fEdisp": -1.754,
        "fEedisp": 0.000,
        "fEtot": -26.177,
    }

    for key1, key2 in zip(fEref.keys(), fEnergies.keys()):
        compare_values(fEref[key1], fEnergies[key2][-1], 2, key1)


@pytest.mark.fsapt
@uusing("pandas")
def test_fsapt_psivars():
    """
    fsapt-psivars: calling fsapt_analysis with psi4 variables after running an
    fisapt0 calcluation requires the user to pass the molecule object
    """
    import pandas as pd

    mol = psi4.geometry(
        """0 1
C 0.00000000 0.00000000 0.00000000
H 1.09000000 0.00000000 0.00000000
H -0.36333333 0.83908239 0.59332085
H -0.36333333 0.09428973 -1.02332709
H -0.36333333 -0.93337212 0.43000624
--
0 1
C 6.44536662 -0.26509169 -0.00000000
H 7.53536662 -0.26509169 -0.00000000
H 6.08203329 0.57399070 0.59332085
H 6.08203329 -0.17080196 -1.02332709
H 6.08203329 -1.19846381 0.43000624
symmetry c1
no_reorient
no_com"""
    )
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
        }
    )
    psi4.core.be_quiet()
    psi4.energy("fisapt0")
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Enuc": 35.07529824960602,
        "Eelst": -3.8035153870907834e-06,
        "Eexch": 1.7912112685446533e-07,
        "Eind": -3.833795151474493e-08,
        "Edisp": -3.288568662589654e-05,
        "Etot": -3.6548418837647605e-05,
    }
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": variable("SAPT ELST ENERGY"),
        "Eexch": variable("SAPT EXCH ENERGY"),
        "Eind": variable("SAPT IND ENERGY"),
        "Edisp": variable("SAPT DISP ENERGY"),
        "Etot": variable("SAPT0 TOTAL ENERGY"),
    }

    for key in keys:
        compare_values(Eref[key], Epsi[key], 6, key)
    data = psi4.fsapt_analysis(
        molecule=mol,
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "MethylA": [1, 2, 3, 4, 5],
        },
        fragments_b={
            "MethylB": [6, 7, 8, 9, 10],
        },
    )
    df = pd.DataFrame(data)
    print(df)
    fEnergies = {}
    fkeys = [
        "fEelst",
        "fEexch",
        "fEindAB",
        "fEindBA",
        "fEdisp",
        "fEedisp",
        "fEtot",
    ]

    Energies = df.iloc[0].values[2:]

    for pair in zip(fkeys, Energies):
        fEnergies[pair[0]] = pair[1]

    fEref = {
        "fEelst": -0.002,
        "fEexch": 0.000,
        "fEindAB": -0.000,
        "fEindBA": -0.000,
        "fEdisp": -0.021,
        "fEedisp": 0.000,
        "fEtot": -0.023,
    }

    for key in fkeys:
        compare_values(fEref[key], fEnergies[key], 2, key)


@pytest.mark.fsapt
@uusing("pandas")
def test_fsapt_AtomicOutput():
    """
    fsapt-atomicOutput: calling fsapt_analysis with atomic results
    """
    import pandas as pd

    mol = psi4.geometry(
        """0 1
C 0.00000000 0.00000000 0.00000000
H 1.09000000 0.00000000 0.00000000
H -0.36333333 0.83908239 0.59332085
H -0.36333333 0.09428973 -1.02332709
H -0.36333333 -0.93337212 0.43000624
--
0 1
C 6.44536662 -0.26509169 -0.00000000
H 7.53536662 -0.26509169 -0.00000000
H 6.08203329 0.57399070 0.59332085
H 6.08203329 -0.17080196 -1.02332709
H 6.08203329 -1.19846381 0.43000624
symmetry c1
no_reorient
no_com"""
    )
    psi4.set_options(
        {
            "basis": "jun-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
        }
    )
    psi4.core.be_quiet()
    plan = psi4.energy("fisapt0", return_plan=True)
    atomic_result = psi4.schema_wrapper.run_qcschema(
        plan.plan(wfn_qcvars_only=False),
        clean=True,
        postclean=True,
    )
    qcvars = atomic_result.extras["qcvars"]
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Enuc": 35.07529824960602,
        "Eelst": -3.8035153870907834e-06,
        "Eexch": 1.7912112685446533e-07,
        "Eind": -3.833795151474493e-08,
        "Edisp": -3.288568662589654e-05,
        "Etot": -3.6548418837647605e-05,
    }
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": qcvars["SAPT ELST ENERGY"],
        "Eexch": qcvars["SAPT EXCH ENERGY"],
        "Eind": qcvars["SAPT IND ENERGY"],
        "Edisp": qcvars["SAPT DISP ENERGY"],
        "Etot": qcvars["SAPT0 TOTAL ENERGY"],
    }

    for key in keys:
        compare_values(Eref[key], Epsi[key], 6, key)
    print("Analysis")
    data = psi4.fsapt_analysis(
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "MethylA": [1, 2, 3, 4, 5],
        },
        fragments_b={
            "MethylB": [6, 7, 8, 9, 10],
        },
        atomic_results=atomic_result,
    )
    df = pd.DataFrame(data)
    print(df)
    fEnergies = {}
    fkeys = [
        "fEelst",
        "fEexch",
        "fEindAB",
        "fEindBA",
        "fEdisp",
        "fEedisp",
        "fEtot",
    ]

    Energies = df.iloc[0].values[2:]

    for pair in zip(fkeys, Energies):
        fEnergies[pair[0]] = pair[1]

    fEref = {
        "fEelst": -0.002,
        "fEexch": 0.000,
        "fEindAB": -0.000,
        "fEindBA": -0.000,
        "fEdisp": -0.021,
        "fEedisp": 0.000,
        "fEtot": -0.023,
    }

    for key in fkeys:
        compare_values(fEref[key], fEnergies[key], 2, key)
