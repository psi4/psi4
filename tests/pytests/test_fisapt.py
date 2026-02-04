import pytest
import psi4
from psi4 import compare_values, variable
from addons import uusing
import numpy as np

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.mark.fsapt
def test_fsapt_psivars_dict():
    """
    Test F-SAPT analysis using dictionary output format (no pandas required).

    This test verifies that fsapt_analysis correctly returns F-SAPT energy
    decomposition as a dictionary after running a fisapt0 calculation. The
    molecule object must be passed to fsapt_analysis when using psi4 variables.

    The test validates:
    1. Standard SAPT energy components against reference values
    2. F-SAPT fragment analysis output in dictionary format
    3. Fragment definitions using 1-indexed atom lists
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
    fEnergies = {
        "Elst": fEnergies["Elst"],
        "Exch": fEnergies["Exch"],
        "IndAB": fEnergies["IndAB"],
        "IndBA": fEnergies["IndBA"],
        "Disp": fEnergies["Disp"],
        "EDisp": fEnergies["EDisp"],
        "Total": fEnergies["Total"],
    }
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
    Test F-SAPT analysis with external point charge potentials.

    This test verifies that fisapt0 calculations can incorporate external
    electrostatic potentials (representing environment effects like solvent
    or additional molecules) and that the resulting SAPT energy components
    are correctly computed. The test uses a three-water system where external
    point charges model water molecules.

    The test validates:
    1. Standard SAPT energy components (Eelst, Eexch, Eind, Edisp, Etot)
    2. F-SAPT fragment analysis using dictionary output format
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
    fEnergies = {
        "Elst": fEnergies["Elst"],
        "Exch": fEnergies["Exch"],
        "IndAB": fEnergies["IndAB"],
        "IndBA": fEnergies["IndBA"],
        "Disp": fEnergies["Disp"],
        "EDisp": fEnergies["EDisp"],
        "Total": fEnergies["Total"],
    }
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
    Test F-SAPT analysis with pandas DataFrame output format.

    This test verifies that fsapt_analysis returns data compatible with pandas
    DataFrame construction after running a fisapt0 calculation. Uses the same
    two-methane system as test_fsapt_psivars_dict to validate consistency
    between output formats, and serves as an example of seemless pandas
    integration.

    Requires: pandas

    The test validates:
    1. Standard SAPT energy components against reference values
    2. F-SAPT energies extracted from DataFrame columns
    3. DataFrame structure with proper column names
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

    df_keys = [
        "Elst",
        "Exch",
        "IndAB",
        "IndBA",
        "Disp",
        "EDisp",
        "Total",
    ]

    # Get columns from dataframe that match fkeys
    Energies = df[df_keys].iloc[0].values

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
        print(fEnergies[key], fEref[key])
        compare_values(fEref[key], fEnergies[key], 2, key)


@pytest.mark.fsapt
def test_fsapt_AtomicOutput():
    """
    Test F-SAPT analysis using QCSchema AtomicResult output (no pandas).

    This test verifies that fsapt_analysis works with QCSchema AtomicResult
    objects returned from run_qcschema, using dictionary output format instead
    of pandas. This approach is useful for integration with QCArchive and
    other QCSchema-compatible workflows. Note, QCArchive will flatten
    arrays, so fsapt_ab_size handles reshaping.

    The test validates:
    1. QCSchema plan generation and execution via run_qcschema
    2. F-SAPT analysis from atomic_results parameter
    3. Dictionary output format without pandas dependency
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
    plan = psi4.energy("fisapt0", return_plan=True, molecule=mol)
    atomic_result = psi4.schema_wrapper.run_qcschema(
        plan.plan(wfn_qcvars_only=False),
        clean=True,
        postclean=True,
    )
    fEnergies = psi4.fsapt_analysis(
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "MethylA": [1, 2, 3, 4, 5],
        },
        fragments_b={
            "MethylB": [6, 7, 8, 9, 10],
        },
        atomic_results=atomic_result,
    )
    fEnergies = {
        "Elst": fEnergies["Elst"],
        "Exch": fEnergies["Exch"],
        "IndAB": fEnergies["IndAB"],
        "IndBA": fEnergies["IndBA"],
        "Disp": fEnergies["Disp"],
        "EDisp": fEnergies["EDisp"],
        "Total": fEnergies["Total"],
    }
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
def test_fsapt_output_file():
    """
    Test F-SAPT analysis with file output (writes fsapt.dat).

    This test verifies that fsapt_analysis can write results to a file
    (fsapt.dat) in the specified directory. This is the traditional F-SAPT
    output format for post-processing or visualization tools.

    The test validates:
    1. Standard SAPT energy components against reference values
    2. File output generation in specified directory
    3. Parsing of fsapt.dat file format
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
        }
    )
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
    psi4.fsapt_analysis(
        fragments_a={
            "MethylA": [1, 2, 3, 4, 5],
        },
        fragments_b={
            "MethylB": [6, 7, 8, 9, 10],
        },
        dirname="./fsapt",
    )
    fEnergies = {}
    fkeys = ["fEelst", "fEexch", "fEindAB", "fEindBA", "fEdisp", "fEedisp", "fEtot"]

    with open("./fsapt/fsapt.dat", "r") as fsapt:
        Energies = [float(x) for x in fsapt.readlines()[-2].split()[2:]]

    for pair in zip(fkeys, Energies):
        fEnergies[pair[0]] = pair[1]

    fEref = {
        "fEelst": -0.002,
        "fEexch": 0.000,
        "fEindAB": -0.000,
        "fEindBA": -0.000,
        "fEdisp": 0.000,
        "fEedisp": -0.033,
        "fEtot": -0.036,
    }

    for key in fkeys:
        compare_values(fEref[key], fEnergies[key], 2, key)
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

    for key in fkeys:
        compare_values(fEref[key], fEnergies[key], 2, key)


@pytest.mark.fsapt
def test_fsapt_indices():
    """
    Test F-SAPT fragment index tracking in multi-fragment analysis.

    This test verifies that fsapt_analysis correctly tracks and reports
    which atom indices belong to each fragment pair in the output. Uses
    a more complex system (ethane + N-methylacetamide) with multiple
    user-defined fragments to test index bookkeeping.

    The test validates:
    1. QCSchema plan generation and execution
    2. Correct fragment index assignment in output dictionary
    3. Multi-fragment definitions with links5050 option
    4. DataFrame construction with fragment indices
    5. Molecules of different sizes

    NOTE: This takes a bit longer to run due to size...
    """
    # psi4.set_memory("32 GB")
    # psi4.set_num_threads(12)

    mol = psi4.geometry(
        """
0 1
C   11.54100       27.68600       13.69600
H   12.45900       27.15000       13.44600
C   10.79000       27.96500       12.40600
H   10.55700       27.01400       11.92400
H   9.879000       28.51400       12.64300
H   11.44300       28.56800       11.76200
H   10.90337       27.06487       14.34224
H   11.78789       28.62476       14.21347
--
0 1
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
C   9.372000       27.59000       8.640000
C   11.44600       26.35600       9.091000
C   9.333000       25.25000       9.282000
H   9.874000       26.68900       6.497000
H   9.908000       28.37100       8.093000
H   8.364000       27.46400       8.233000
H   9.317000       27.84600       9.706000
H   9.807000       24.28200       9.160000
H   9.371000       25.57400       10.32900
H   8.328000       25.26700       8.900000
H   11.28800       26.57600       10.14400
H   11.97000       27.14900       8.585000
H   11.93200       25.39300       8.957000
H   10.61998       24.85900       5.366911
units angstrom

symmetry c1
no_reorient
no_com
"""
    )
    psi4.set_options(
        {
            # "basis": "sto-3g",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "guess": "sad",
            # "freeze_core": "true",
        }
    )
    plan = psi4.energy("fisapt0", return_plan=True, molecule=mol)
    atomic_result = psi4.schema_wrapper.run_qcschema(
        plan.plan(wfn_qcvars_only=False),
        clean=True,
        postclean=True,
    )
    data = psi4.fsapt_analysis(
        # NOTE: 1-indexed for fragments_a and fragments_b
        fragments_a={
            "Methyl1_A": [1, 2, 7, 8],
            "Methyl2_A": [3, 4, 5, 6],
        },
        fragments_b={
            "Peptide_B": [9, 10, 11, 16, 26],
            "T-Butyl_B": [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        },
        links5050=True,
        print_output=False,
        atomic_results=atomic_result,
    )
    from pprint import pprint

    pprint(data)
    mol_qcel_dict = mol.to_schema(dtype=2)
    frag1_indices = data["Frag1_indices"]
    frag2_indices = data["Frag2_indices"]
    # Using molecule object for all test to ensure right counts from each
    # fragment are achieved. Note +1 for 1-indexing in fsapt_analysis
    all_A = [i + 1 for i in mol_qcel_dict["fragments"][0]]
    expected_frag1_indices = [
        [1, 2, 7, 8],
        [1, 2, 7, 8],
        [3, 4, 5, 6],
        [3, 4, 5, 6],
        [1, 2, 7, 8],
        [3, 4, 5, 6],
        all_A,
        all_A,
        all_A,
    ]
    all_B = [j + 1 for j in mol_qcel_dict["fragments"][1]]
    expected_frag2_indices = [
        [9, 10, 11, 16, 26],
        [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        [9, 10, 11, 16, 26],
        [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        all_B,
        all_B,
        [9, 10, 11, 16, 26],
        [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        all_B,
    ]
    for i, indices in enumerate(frag1_indices):
        # Assert lists are identical
        e = expected_frag1_indices[i]
        sorted_frag = sorted(indices)
        assert sorted_frag == e, f"Frag1 indices do not match for fragment {
            i
        }: expected {e}, got {sorted_frag}"

    for i, indices in enumerate(frag2_indices):
        e = expected_frag2_indices[i]
        sorted_frag = sorted(indices)
        assert sorted_frag == e, f"Frag2 indices do not match for fragment {
            i
        }: expected {e}, got {sorted_frag}"
    ref_dict = {
        "ClosestContact": [
            12.99840199731447,
            6.708905946098247,
            9.420620786025163,
            3.7293279474020324,
            6.708905946098247,
            3.7293279474020324,
            9.420620786025163,
            3.7293279474020324,
            3.7293279474020324,
        ],
        "Disp": [
            -0.02014331013167378,
            -0.372273895991978,
            -0.0741518363998009,
            -2.432207594127432,
            -0.39241720612365183,
            -2.506359430527233,
            -0.09429514653147468,
            -2.80448149011941,
            -2.898776636650885,
        ],
        "EDisp": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Elst": [
            0.8664961096130384,
            -0.12491437414747963,
            -1.0161642439176575,
            -1.0168685763655532,
            0.7415817354655587,
            -2.0330328202832106,
            -0.14966813430461912,
            -1.1417829505130328,
            -1.291451084817652,
        ],
        "Exch": [
            -0.00016951134321586871,
            0.06302818590927417,
            0.04056059018694508,
            4.214334017378195,
            0.0628586745660583,
            4.25489460756514,
            0.04039107884372921,
            4.277362203287469,
            4.317753282131198,
        ],
        "Frag1": [
            "Methyl1_A",
            "Methyl1_A",
            "Methyl2_A",
            "Methyl2_A",
            "Methyl1_A",
            "Methyl2_A",
            "All",
            "All",
            "All",
        ],
        "Frag2": [
            "Peptide_B",
            "T-Butyl_B",
            "Peptide_B",
            "T-Butyl_B",
            "All",
            "All",
            "Peptide_B",
            "T-Butyl_B",
            "All",
        ],
        "IndAB": [
            -0.027615863638335056,
            -0.015567880312064623,
            -0.09781659817652713,
            -0.26487802957299955,
            -0.04318374395039968,
            -0.36269462774952665,
            -0.12543246181486217,
            -0.2804459098850642,
            -0.40587837169992635,
        ],
        "IndBA": [
            0.0011573940269127988,
            0.03052724256475823,
            -0.002144092785043809,
            -0.13667169442200203,
            0.03168463659167103,
            -0.13881578720704585,
            -0.00098669875813101,
            -0.1061444518572438,
            -0.10713115061537481,
        ],
        "Total": [
            0.8197248185264883,
            -0.4192007219755354,
            -1.1497161810927068,
            0.3637081228885606,
            0.4005240965509529,
            -0.7860080582041462,
            -0.3299913625662185,
            -0.05549259908697479,
            -0.3854839616531933,
        ],
    }
    for key in ['ClosestContact', 'Elst', 'Exch', 'IndAB', 'IndBA', 'Disp', 'EDisp', 'Total']:
        for i, value in enumerate(data[key]):
            f1_f2 = f"{data['Frag1'][i]}-{data['Frag2'][i]}"
            print(
                f1_f2,
                ref_dict[key][i],
                value,
            )
            compare_values(
                ref_dict[key][i],
                value,
                6,
                f"Fragment pair {f1_f2}:{i} for key {key}",
            )
    return


if __name__ == "__main__":
    # test_fsapt_psivars_dict()
    # test_fsapt_external_potentials()
    # test_fsapt_psivars()
    # test_fsapt_psivars_dict()
    # test_fsapt_AtomicOutput()
    # test_fsapt_output_file()
    # test_fsapt_output_file()
    test_fsapt_indices()
    # pytest.main([__file__])
