import pytest
import psi4
from qcelemental import constants
from psi4 import compare_values
from psi4 import core
import numpy as np
from pprint import pprint as pp
from addons import uusing
import os

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


@uusing("pandas")
def test_fsaptdft_timer():
    """Ensure SAPT(DFT) timer CSV output contains expected timing columns."""
    import pandas as pd

    mol = psi4.geometry(_sapt_testing_mols["neutral_water_dimer"])
    np.set_printoptions(precision=10, suppress=True)
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            "FISAPT_FSAPT_FILEPATH": "none",
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
            "SAPT_DFT_DO_FSAPT": "FISAPT",
            "SAPT_DFT_USE_EINSUMS": True,
        }
    )
    psi4.core.clean_timers()
    _, wfn = psi4.energy("sapt(dft)", molecule=mol, return_wfn=True)
    compute_time_saptdft_fi_ein = psi4.core.get_timer_dict()["SAPT(DFT) Energy"]
    psi4.driver.p4util.write_timer_csv("saptdft_fi_useEin_timers.csv")
    df = pd.read_csv("saptdft_fi_useEin_timers.csv")
    os.remove("saptdft_fi_useEin_timers.csv")
    print(f"compute_time_fi_ein: {compute_time_saptdft_fi_ein['wall_time']:.2f}s\n")
    print(df)
    timer_cols = ["timer_name", "wall_time", "user_time", "system_time", "n_calls"]
    for col in timer_cols:
        assert col in df.columns, f"Expected column '{col}' not found in timer CSV"
    return


@pytest.mark.saptdft
@pytest.mark.fsapt
@pytest.mark.saptdft
def test_fsapthf_disp0_fisapt0_psivars():
    """Validate HF SAPT(DFT)+FISAPT energies and fragment terms vs references."""

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
    # Reference energies from FISAPT0/sto-3g
    data = {
        "Disp": [
            -0.003994502630897984,
            -0.0674134405493367,
            -0.013547073146612729,
            -0.41154085645677463,
            -0.07140794318023469,
            -0.42508792960338737,
            -0.017541575777510712,
            -0.47895429700611136,
            -0.49649587278362206,
        ],
        "EDisp": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Elst": [
            0.7173658716748221,
            -0.2055651208882452,
            -0.8178788521844282,
            -0.9342087776241712,
            0.5118007507865769,
            -1.7520876298085994,
            -0.10051298050960611,
            -1.1397738985124164,
            -1.2402868790220225,
        ],
        "Exch": [
            0.000135454393737376,
            0.04719682210839049,
            0.03161592765410519,
            3.8965970121551616,
            0.04733227650212787,
            3.9282129398092667,
            0.031751382047842565,
            3.943793834263552,
            3.9755452163113945,
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
            -0.007097098257764316,
            -0.015628832176699036,
            -0.02607145661295123,
            -0.17470490188100032,
            -0.022725930434463353,
            -0.20077635849395156,
            -0.033168554870715544,
            -0.19033373405769935,
            -0.22350228892841492,
        ],
        "IndBA": [
            0.0003539943187897063,
            0.014741295723273411,
            -0.0017520923518852605,
            -0.0807135744181731,
            0.015095290042063118,
            -0.08246566677005836,
            -0.0013980980330955543,
            -0.06597227869489969,
            -0.06737037672799524,
        ],
        "Total": [
            0.706763719498241,
            -0.22666927578227103,
            -0.8276335466433693,
            2.2954289017738034,
            0.48009444371597,
            1.4677953551304341,
            -0.12086982714512828,
            2.0687596259915324,
            1.9478897988464041,
        ],
    }

    ref_data = data
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Edisp": -0.0007912165332931369,
        "Eelst": -0.0019765265492708295,
        "Eexch": 0.006335438658802855,
        "Eind": -0.0004635353239533062,
        "Enuc": 474.74808217020274,
        "Etot": 0.003104160252285582,
    }
    print("SAPT_DFT_DO_FSAPT = FISAPT0 now testing")
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "FISAPT_FSAPT_FILEPATH": "none",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            "SAPT_DFT_DO_FSAPT": "FISAPT",
            "SAPT_DFT_D4_IE": False,
            "SAPT_DFT_DO_DISP": True,
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
            # Normally on
            "SAPT_DFT_USE_EINSUMS": True,
        }
    )
    _, wfn = psi4.energy("sapt(dft)", molecule=mol, return_wfn=True)
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    pp(Epsi)
    for key in keys:
        compare_values(Eref[key], Epsi[key], 5, key)
    fsapt_data = psi4.fsapt_analysis(
        source=wfn,
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
    )
    for label_key in ["Frag1", "Frag2"]:
        assert list(fsapt_data[label_key]) == ref_data[label_key]
    for col in ["Elst", "Exch", "IndAB", "IndBA", "Disp", "EDisp", "Total"]:
        for i in range(len(ref_data[col])):
            compare_values(
                ref_data[col][i],
                fsapt_data[col][i],
                4,
                f"{ref_data['Frag1'][i]} {ref_data['Frag2'][i]} {col}",
            )


@pytest.mark.saptdft
@pytest.mark.fsapt
@uusing("pandas")
@pytest.mark.saptdft
def test_fsaptdftd4_psivars():
    """Validate SAPT(DFT)-D4(s) scalar variables and FSAPT terms vs references."""
    import pandas as pd

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
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "FISAPT_FSAPT_FILEPATH": "none",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            "SAPT_DFT_DO_FSAPT": "SAPTDFT",
        }
    )
    _, wfn = psi4.energy("sapt(dft)-d4(s)", return_wfn=True)

    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    Eref = {
        "Edisp": -0.004568534767691285,
        "Eelst": -0.0019765266134612602,
        "Eexch": 0.006335438658900877,
        "Eind": -0.0004635353246623952,
        "Enuc": 474.74808217020274,
        "Etot": -0.0006731581723675157,
    }
    Epsi = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    pp(Epsi)
    for key in keys:
        compare_values(Eref[key], Epsi[key], 5, key)
    data = psi4.fsapt_analysis(
        source=wfn,
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
    )
    df = pd.DataFrame(data)
    print("COMPUTED DF")
    print(df)
    pp({k: v.tolist() for k, v in dict(df).items()})
    ref_data = {
        "Disp": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "EDisp": [
            -0.01301450872859415,
            -0.27084948543172516,
            -0.04866416123143305,
            -2.4459445415531187,
            -0.2838639941603193,
            -2.494608702784552,
            -0.0616786699600272,
            -2.716794026984844,
            -2.778472696944871,
        ],
        "Elst": [
            0.7173658713369733,
            -0.20556512078896816,
            -0.8178788520169107,
            -0.9342087766410643,
            0.5118007505480051,
            -1.752087628657975,
            -0.10051298067993741,
            -1.1397738974300324,
            -1.2402868781099698,
        ],
        "Exch": [
            0.00013545439373716334,
            0.047196822108368085,
            0.031615927654092296,
            3.8965970121551456,
            0.047332276502105246,
            3.928212939809238,
            0.03175138204782946,
            3.9437938342635137,
            3.975545216311343,
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
        "Frag1_indices": [
            [1, 2, 7, 8],
            [1, 2, 7, 8],
            [3, 4, 5, 6],
            [3, 4, 5, 6],
            [1, 2, 7, 8],
            [3, 4, 5, 6],
            [1, 2, 7, 8, 3, 4, 5, 6],
            [1, 2, 7, 8, 3, 4, 5, 6],
            [1, 2, 7, 8, 3, 4, 5, 6],
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
        "Frag2_indices": [
            [9, 10, 11, 16, 26],
            [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [9, 10, 11, 16, 26],
            [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [9, 10, 11, 16, 26, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [9, 10, 11, 16, 26, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [9, 10, 11, 16, 26],
            [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [9, 10, 11, 16, 26, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        ],
        "IndAB": [
            -0.007097098275574362,
            -0.01562883221584556,
            -0.026071456678509793,
            -0.17470490231988067,
            -0.02272593049141992,
            -0.20077635899839047,
            -0.033168554954084155,
            -0.19033373453572622,
            -0.22350228948981038,
        ],
        "IndBA": [
            0.0003539943196737901,
            0.014741295760311841,
            -0.001752092356286135,
            -0.08071357462105197,
            0.015095290079985632,
            -0.0824656669773381,
            -0.0013980980366123448,
            -0.06597227886074013,
            -0.06737037689735248,
        ],
        "Total": [
            0.6977436680167427,
            -0.4301052921347477,
            -0.8627506894690783,
            0.2610253266522511,
            0.267638375881995,
            -0.6017253628168272,
            -0.1650070214523356,
            -0.1690799654824966,
            -0.3340869869348322,
        ],
    }

    ref_df = pd.DataFrame(ref_data)
    print("REF")
    print(ref_df)

    for col in ["Elst", "Exch", "IndAB", "IndBA", "Disp", "EDisp", "Total"]:
        for i in range(len(ref_df)):
            compare_values(
                ref_df[col].iloc[i],
                df[col].iloc[i],
                4,
                f"{ref_df['Frag1'].iloc[i]} {ref_df['Frag2'].iloc[i]} {col}",
            )


@pytest.mark.fsapt
@pytest.mark.saptdft
def test_fsaptdftd4_psivars_pbe0_frozen_core():
    """Check PBE0 SAPT(DFT)-D4(i) fragment indices and per-fragment energies."""

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
    print("FSAPT(PBE0)-D4(I)")
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "FISAPT_FSAPT_FILEPATH": "none",
            "SAPT_DFT_FUNCTIONAL": "PBE0",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            "SAPT_DFT_DO_FSAPT": "SAPTDFT",
            "SAPT_DFT_GRAC_SHIFT_A": 0.11652342,
            "SAPT_DFT_GRAC_SHIFT_B": 0.12724880,
        }
    )
    _, wfn = psi4.energy("sapt(dft)-d4(i)", molecule=mol, return_wfn=True)
    fsapt_data = psi4.fsapt_analysis(
        # NOTE: 1-indexed for fragments_a and fragments_b
        source=wfn,
        fragments_a={
            "Methyl1_A": [1, 2, 7, 8],
            "Methyl2_A": range(3, 7),
        },
        fragments_b={
            "Peptide_B": [9, 10, 11, 16, 26],
            "T-Butyl_B": [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        },
        links5050=True,
        print_output=False,
    )
    mol_qcel_dict = mol.to_schema(dtype=2)
    frag1_indices = fsapt_data["Frag1_indices"]
    frag2_indices = fsapt_data["Frag2_indices"]
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
    print(f"{all_A=}")
    print(f"{all_B=}")
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
    data = {
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
        "Elst": [
            0.6885820849831674,
            -0.1340281136412642,
            -0.7952455143967398,
            -1.0514370965090478,
            0.5545539713419032,
            -1.8466826109057877,
            -0.10666342941357243,
            -1.185465210150312,
            -1.2921286395638845,
        ],
        "Exch": [
            0.0001588422784532737,
            0.04729543674199185,
            0.039013303979676645,
            4.045753472042623,
            0.04745427902044513,
            4.0847667760223,
            0.03917214625812992,
            4.093048908784615,
            4.132221055042745,
        ],
        "IndAB": [
            -0.008959884719316583,
            -0.014829489437038195,
            -0.029218977876373193,
            -0.17897054392702005,
            -0.02378937415635478,
            -0.20818952180339326,
            -0.038178862595689776,
            -0.19380003336405824,
            -0.23197889595974802,
        ],
        "IndBA": [
            0.00063667599072615,
            0.019609289365314572,
            -0.0024984807073431003,
            -0.09577726569805996,
            0.020245965356040722,
            -0.09827574640540306,
            -0.0018618047166169503,
            -0.07616797633274539,
            -0.07802978104936234,
        ],
        "Disp": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "EDisp": [
            -0.013635081934517778,
            -0.2986013898239719,
            -0.052012447785209845,
            -3.3086784527471322,
            -0.3122364717584897,
            -3.360690900532342,
            -0.06564752971972762,
            -3.607279842571104,
            -3.6729273722908315,
        ],
        "Total": [
            0.6667826311549935,
            -0.3805542969583513,
            -0.8399621243214389,
            -0.5891099251986289,
            0.2862283341966422,
            -1.4290720495200677,
            -0.17317949316644543,
            -0.9696642221569802,
            -1.1428437153234254,
        ],
    }
    cols = [
        "Frag1",
        "Frag2",
        "Elst",
        "Exch",
        "IndAB",
        "IndBA",
        "Disp",
        "EDisp",
        "Total",
    ]
    ref_data = data

    assert list(fsapt_data["Frag1"]) == ref_data["Frag1"]
    assert list(fsapt_data["Frag2"]) == ref_data["Frag2"]

    for col in cols[2:]:
        for i in range(len(ref_data[col])):
            compare_values(
                ref_data[col][i],
                fsapt_data[col][i],
                4,
                f"{ref_data['Frag1'][i]} {ref_data['Frag2'][i]} {col}",
            )


@pytest.mark.saptdft
@pytest.mark.fsapt
def test_fsaptdft_fisapt0():
    """Confirm fisapt0 and HF SAPT(DFT) produce matching SAPT energy components."""
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

    # Run standard FISAPT0
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
        }
    )
    psi4.energy("fisapt0", molecule=mol)

    # Collect FISAPT0 energies
    fisapt0_energies = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    print("FISAPT0 energies:")
    pp(fisapt0_energies)

    # Clear variables for next calculation
    psi4.core.clean()
    psi4.core.clean_variables()

    # Run SAPT(DFT) with FISAPT option (HF functional to match SAPT0)
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "none",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            # "SAPT_DFT_DO_FSAPT": "FISAPT",
            "SAPT_DFT_DO_FSAPT": "SAPTDFT",
            "SAPT_DFT_D4_IE": False,
            "SAPT_DFT_DO_DISP": True,
            "SAPT_DFT_MP2_DISP_ALG": "FISAPT",
            "SAPT_DFT_USE_EINSUMS": False,
        }
    )
    _, wfn = psi4.energy("sapt(dft)", molecule=mol, return_wfn=True)

    # Collect SAPT(DFT) energies
    saptdft_energies = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    print("SAPT(DFT) with FISAPT energies:")
    pp(saptdft_energies)

    # Compare total energies (5 decimal places = ~0.01 kcal/mol precision)
    keys = ["Enuc", "Eelst", "Eexch", "Eind", "Edisp", "Etot"]
    for key in keys:
        compare_values(
            fisapt0_energies[key],
            saptdft_energies[key],
            5,
            f"Total {key}",
        )


@pytest.mark.saptdft
@pytest.mark.fsapt
def test_fsaptdft_fisapt0_d4():
    """Compare fisapt0-d4 and SAPT(DFT)-D4(i)/FISAPT fragment decompositions."""

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

    # Run standard FISAPT0
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "FISAPT_FSAPT_FILEPATH": "tmp_fisapt",
        }
    )
    _, wfn = psi4.energy("fisapt0-d4", molecule=mol, return_wfn=True)
    with open("tmp_fisapt/fA.dat", "w") as fA:  # TEST
        fA.write("Methyl1_A 1 2 7 8\n")  # TEST
    with open("tmp_fisapt/fB.dat", "w") as fB:  # TEST
        fB.write("Peptide_B  9 10 11 16 26\n")  # TEST
        fB.write("T-Butyl_B  12 13 14 15 17 18 19 20 21 22 23 24 25")  # TEST
    psi4.fsapt_analysis(
        source=wfn,
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
        pdb_dir="tmp_fisapt",
    )
    # remove_fisapt files
    import shutil

    shutil.rmtree("tmp_fisapt")
    # rm pdb
    pdb_files = [
        "Disp.pdb",
        "EDisp.pdb",
        "Elst.pdb",
        "Exch.pdb",
        "IndAB.pdb",
        "IndBA.pdb",
        "Total.pdb",
    ]
    for pdb_file in pdb_files:
        if os.path.exists(pdb_file):
            os.remove(pdb_file)
    # Collect FISAPT0 energies
    fisapt0_energies = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    print("FISAPT0 energies:")
    pp(fisapt0_energies)

    # Clear variables for next calculation
    psi4.core.clean()
    psi4.core.clean_variables()

    # Run SAPT(DFT) with FISAPT option (HF functional to match SAPT0)
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "df",
            "guess": "sad",
            "freeze_core": "true",
            "SAPT_DFT_FUNCTIONAL": "HF",
            "SAPT_DFT_DO_DHF": True,
            "SAPT_DFT_DO_HYBRID": False,
            "SAPT_DFT_DO_FSAPT": "FISAPT",
            "SAPT_DFT_USE_EINSUMS": True,
            "FISAPT_FSAPT_FILEPATH": "tmp",
        }
    )
    _, wfn = psi4.energy("sapt(dft)-d4(i)", molecule=mol, return_wfn=True)

    # Collect SAPT(DFT) energies
    saptdft_energies = {
        "Enuc": mol.nuclear_repulsion_energy(),
        "Eelst": core.variable("SAPT ELST ENERGY"),
        "Eexch": core.variable("SAPT EXCH ENERGY"),
        "Eind": core.variable("SAPT IND ENERGY"),
        "Edisp": core.variable("SAPT DISP ENERGY"),
        "Etot": core.variable("SAPT TOTAL ENERGY"),
    }
    print("SAPT(DFT) with FISAPT energies:")
    pp(saptdft_energies)
    saptdft_fsapt_data = psi4.fsapt_analysis(
        source=wfn,
        fragments_a={
            "Methyl1_A": [1, 2, 7, 8],
        },
        fragments_b={
            "Peptide_B": [9, 10, 11, 16, 26],
            "T-Butyl_B": [12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25],
        },
        links5050=True,
        print_output=False,
        pdb_dir="tmp",
    )
    # remove_fisapt files
    shutil.rmtree("tmp")
    for pdb_file in pdb_files:
        if os.path.exists(pdb_file):
            os.remove(pdb_file)
    ref_data = {
        "Disp": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "EDisp": [
            -0.01301450872859415,
            -0.27084948543172516,
            -0.04866416123143305,
            -2.4459445415531187,
            -0.2838639941603193,
            -2.494608702784552,
            -0.0616786699600272,
            -2.716794026984844,
            -2.778472696944871,
        ],
        "Elst": [
            0.7150991989031183,
            -0.20424514540867733,
            -0.8155054641111548,
            -0.9356354507588378,
            0.510854053494441,
            -1.7511409148699926,
            -0.10040626520803642,
            -1.1398805961675151,
            -1.2402868613755516,
        ],
        "Exch": [
            0.00013680473229194022,
            0.053104025194981225,
            0.030944134785956742,
            3.8913602515981354,
            0.053240829927273164,
            3.9223043863840923,
            0.031080939518248682,
            3.9444642767931164,
            3.9755452163113656,
        ],
        "IndAB": [
            -0.007088768001216014,
            -0.015599384317671102,
            -0.02601507069216381,
            -0.1747990661073405,
            -0.022688152318887114,
            -0.2008141367995043,
            -0.03310383869337982,
            -0.1903984504250116,
            -0.22350228911839143,
        ],
        "IndBA": [
            0.0003529338034038042,
            0.01470702778121176,
            -0.0017519903139542258,
            -0.0806783480564871,
            0.015059961584615564,
            -0.08243033837044134,
            -0.0013990565105504217,
            -0.06597132027527534,
            -0.06737037678582578,
        ],
        "Total": [
            0.6954856149378384,
            -0.4228829336993627,
            -0.8609926109811682,
            0.25430294407805487,
            0.27260268123847575,
            -0.6066896669031133,
            -0.16550699604332975,
            -0.16857998962130782,
            -0.33408698566463757,
        ],
    }
    keys = ["Elst", "Exch", "IndAB", "IndBA", "Disp", "EDisp", "Total"]
    for key in keys:
        for i in range(len(ref_data[key])):
            compare_values(
                ref_data[key][i],
                saptdft_fsapt_data[key][i],
                6,
                f"{saptdft_fsapt_data['Frag1'][i]} {saptdft_fsapt_data['Frag2'][i]} {key}",
            )


if __name__ == "__main__":
    psi4.set_memory("220 GB")
    # psi4.set_num_threads(24)
    psi4.set_num_threads(12)
    test_fsaptdft_fisapt0()
    # pytest.main([
    #     __file__,
    #     "-v",
    #     "-s",
    #     # "-k=test_saptdft_auto_grac",
    #     "--disable-warnings",
    #     # "--maxfail=1",
    # ])
