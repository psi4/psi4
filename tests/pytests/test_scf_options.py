import numpy as np
import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"options": {"guess": "core"}}, id="core"),
        pytest.param({"options": {"guess": "gwh"}}, id="gwh"),
        pytest.param({"options": {"guess": "huckel"}}, id="huckel"),
        pytest.param({"options": {"guess": "modhuckel"}}, id="modhuckel"),
        pytest.param({"late_options": {"guess": "read"}}, id="read"),
        pytest.param({"options": {"guess": "sad"}}, id="sad"),
        pytest.param({"options": {"guess": "sadno"}}, id="sadno"),
        pytest.param({"options": {"guess": "sap"}}, id="sap"),
    ],
)
def test_guess_mix_for_broken_symmetry(inp):

    refENuc = 0.17639240356
    refSCF = -0.82648407827446
    refBSSCF = -0.99872135103903

    h2 = psi4.geometry(
        """
        0 1
        H
        H 1 3.0
        symmetry c1
        """
    )

    psi4.set_options({"reference": "uhf", "e_convergence": 12, "basis": "cc-pvdz"})
    psi4.set_options(inp.get("options", {}))

    thisSCF = psi4.energy("scf")
    psi4.set_options(inp.get("late_options", {}))
    psi4.set_options({"guess_mix": True})
    thisBSSCF = psi4.energy("scf")

    assert compare_values(
        refENuc, h2.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy"
    )
    assert compare_values(refSCF, thisSCF, 10, "Reference energy")
    assert psi4.compare_values(refBSSCF, thisBSSCF, 10, "Reference broken-symmetry energy")


@pytest.mark.parametrize("ref", ["rhf", "uhf", "rohf"])
@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"options": {"guess": "core"}}, id="core"),
        pytest.param({"options": {"guess": "gwh"}}, id="gwh"),
        pytest.param({"options": {"guess": "huckel"}}, id="huckel"),
        pytest.param({"options": {"guess": "modhuckel"}}, id="modhuckel"),
        pytest.param({"options": {"guess": "sad"}}, id="sad"),
        pytest.param({"options": {"guess": "sadno"}}, id="sadno"),
        pytest.param({"options": {"guess": "sap"}}, id="sap"),
        # pytest.param({'options': {"guess": "sapgau"}}, id="sapgau"),
    ],
)
def test_scf_guess(inp, ref):
    vals = {
        "rhf": {
            "core": -85.33505416384983,
            "gwh": -92.48691177199250,
            "sad": -99.63941801281894,
            "sadno": -99.77439003262768,
            "huckel": -99.82915309755559,
            "modhuckel": -99.80751245771735,
            "sap": -99.97316521623380,
        },
        "uhf": {
            "core": -88.99202896495355,
            "gwh": -94.35585144644561,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641045,
            "huckel": -99.20374047108815,
            "modhuckel": -99.17689049002892,
            "sap": -99.50144270485751,
        },
        "rohf": {
           "core": -88.99202896495353,
            "gwh": -94.35585144644557,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641042,
            "huckel": -99.20374047108817,
            "modhuckel": -99.17689049002892,
            "sap": -99.50144270485754,
        },
    }
    ref_final = {
        # "refnuc": 5.282967161430950,
        "rhf": -100.0584459442408587,
        "uhf": -99.5312257221445549,
        "rohf": -99.5261713512123976,
        "cuhf": -99.5261713512123976,
    }
    options = {
            "reference": ref,
            "basis": "cc-pvtz",
            "scf_type": "pk",
            "df_scf_guess": False,
            "guess": inp["options"]["guess"],
        }
    if ref == "rhf":
        hf = psi4.geometry(
            """
    0 1
    F
    H 1 0.9015
    """
        )
        options["docc"] = [3, 0, 1, 1]
    else:
        hf = psi4.geometry(
            """
    1 2
    F
    H 1 0.9015
    """
        )
        options["docc"] = [3, 0, 0, 1]
        options["socc"] = [0, 0, 1, 0]

    psi4.set_options(options)
    psi4.energy("scf")
    assert compare_values(
        ref_final[ref], psi4.variable("SCF TOTAL ENERGY"), 6, "FINAL SCF ENERGY"
    )
    if psi4.core.get_option("scf", "orbital_optimizer_package") == "INTERNAL":  # KP-DIFF-REF
      # until pyooo and SCF TOTAL ENERGIES filled out
      assert compare_values(
        vals[ref][inp["options"]["guess"]],
        psi4.variable("SCF TOTAL ENERGIES")[0],
        6,
        "INITIAL ITERATION",
      )

@pytest.mark.quick
@pytest.mark.parametrize("atom,geom,expected_spin_electrons,expected_atom_alpha", [
    ("hydrogen", "0 2\nH\nsymmetry c1\n", 0.5, [0.5]),
    ("oxygen", "0 3\nO\nsymmetry c1\n", 4.0, [4.0]),
    ("water", "0 1\nO\nH 1 1.0\nH 1 1.0 2 90.0\nsymmetry c1\n", 5.0, [4.0, 0.5, 0.5]),
])
def test_sad_frac_occ_density_occupation(atom, geom, expected_spin_electrons, expected_atom_alpha):
    """SAD fractional occupations should conserve total electrons"""

    atom_mol = psi4.geometry(geom)

    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "reference": "uhf",
            "guess": "sad",
            "sad_functional": "HF",
            "scf_type": "pk",
        }
    )

    base_wfn = psi4.core.Wavefunction.build(atom_mol, psi4.core.get_global_option("BASIS"))
    hf_wfn = psi4.driver.scf_wavefunction_factory("HF", base_wfn, "UHF")
    hf_wfn.initialize()
    hf_wfn.guess()

    S = hf_wfn.S()
    Da = hf_wfn.Da()
    Db = hf_wfn.Db()
    na = Da.vector_dot(S)
    nb = Db.vector_dot(S)

    basis = hf_wfn.basisset()
    atom_from_ao = np.fromiter((basis.function_to_center(ao) for ao in range(basis.nbf())), dtype=int)
    DaS_diag = np.diag(np.asarray(Da) @ np.asarray(S))
    DbS_diag = np.diag(np.asarray(Db) @ np.asarray(S))

    atom_occupations = []
    for atom_idx in range(atom_mol.natom()):
        ao_mask = atom_from_ao == atom_idx
        atom_na = DaS_diag[ao_mask].sum()
        atom_nb = DbS_diag[ao_mask].sum()
        atom_occupations.append((atom_mol.symbol(atom_idx), atom_na, atom_nb))

    #           H H H H H H H       O O O O O O O       H2O H2O H2O H2O H2O H2O H2O
    # PR        na=nb   ntot        na=nb   ntot        O_na=nb O_tot   H_na=nb H_tot
    # pre-3138  0.5     1.0         4.0     8.0         4.0     8.0     0.5     1.0
    # 3138      0.707   1.414       4.464   8.928       4.464   8.928   0.707   1.414
    # 3390      0.5     1.0         4.0     8.0         4.0     8.0     0.5     1.0

    print(f"{na=} {nb=}")
    print("per-atom occupations (alpha beta total):")
    for idx, (symbol, atom_na, atom_nb) in enumerate(atom_occupations):
        print(f"{idx:2d} {symbol:>2s}: {atom_na:12.8f} {atom_nb:12.8f} {atom_na + atom_nb:12.8f}")

    for idx, expected_na in enumerate(expected_atom_alpha):
        symbol, atom_na, atom_nb = atom_occupations[idx]
        assert compare_values(expected_na, atom_na, 10, f"{atom} atom {idx} ({symbol}) SAD alpha occupation")
        assert compare_values(expected_na, atom_nb, 10, f"{atom} atom {idx} ({symbol}) SAD beta occupation")

    # SADGuess returns a spin-restricted density for HF guesses (Db == Da)
    assert compare_values(2 * expected_spin_electrons, na + nb, 10, f"{atom} SAD total density occupation")
    assert compare_values(expected_spin_electrons, na, 10, f"{atom} SAD alpha density occupation")
    assert compare_values(expected_spin_electrons, nb, 10, f"{atom} SAD beta density occupation")


@pytest.mark.quick
@pytest.mark.parametrize("atom,geom,expected_spin_electrons,expected_occ_prefix", [
    ("hydrogen", "0 2\nH\nsymmetry c1\n", 0.5, [0.5]),
    ("oxygen", "0 3\nO\nsymmetry c1\n", 4.0, [1.0, 0.75, 0.75, 0.75, 0.75]),
    ("water", "0 1\nO\nH 1 1.0\nH 1 1.0 2 90.0\nsymmetry c1\n", 5.0, [1.25901924, 0.99495748, 0.84219060, 0.75, 0.75, 0.24318432, 0.16064836]),
])
def test_sad_frac_occ_atomic_natural_occupations(atom, geom, expected_spin_electrons, expected_occ_prefix):
    """SAD fractional occupations should produce the expected atomic natural occupations"""

    # above reference occupations work for pre-3138 and 3390 PRs
    # below are values for 3138
    # H atm [ 0.7071,  0.0, ...] = 0.71
    # O atm [ 1.0000,  0.8660,  0.8660,  0.8660,  0.8660,  0.0, ...] = 4.46
    # water [ 1.6461,  1.0071,  0.9983,  0.8660,  0.8660,  0.2952,  0.1993,  0.0, ...] = 5.88

    atom_mol = psi4.geometry(geom)

    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "reference": "uhf",
            "guess": "sad",
            "sad_functional": "HF",
            "scf_type": "pk",
        }
    )

    base_wfn = psi4.core.Wavefunction.build(atom_mol, psi4.core.get_global_option("BASIS"))
    hf_wfn = psi4.driver.scf_wavefunction_factory("HF", base_wfn, "UHF")
    hf_wfn.initialize()
    hf_wfn.guess()

    Da = np.asarray(hf_wfn.Da())
    S = np.asarray(hf_wfn.S())

    S_eval, S_evec = np.linalg.eigh(S)
    S_half = S_evec @ np.diag(np.sqrt(np.clip(S_eval, a_min=0.0, a_max=None))) @ S_evec.T
    D_orth = S_half @ Da @ S_half
    occ = np.sort(np.linalg.eigvalsh(0.5 * (D_orth + D_orth.T)))[::-1]

    print(f"{atom} alpha natural occupations (top): {occ[:len(expected_occ_prefix) + 3]}")

    assert compare_values(expected_spin_electrons, occ.sum(), 10, f"{atom} SAD alpha natural-occupation sum {occ.sum():.4f}")
    for i, expected_occ in enumerate(expected_occ_prefix):
        assert compare_values(expected_occ, occ[i], 8, f"{atom} SAD alpha natural occupation {i} {occ[i]:.4f}")
    for i in range(len(expected_occ_prefix), len(occ)):
        assert compare_values(0.0, occ[i], 8, f"{atom} SAD alpha virtual natural occupation {i} {occ[i]:.4f}")

# pre-3138
#per-atom occupations (alpha beta total):
# 0  H:   0.50000000   0.50000000   1.00000000
#per-atom occupations (alpha beta total):
# 0  O:   4.00000000   4.00000000   8.00000000
#per-atom occupations (alpha beta total):
# 0  O:   4.00000000   4.00000000   8.00000000
# 1  H:   0.50000000   0.50000000   1.00000000
# 2  H:   0.50000000   0.50000000   1.00000000

# 3138
#per-atom occupations (alpha beta total):
# 0  H:   0.70710678   0.70710678   1.41421356
#per-atom occupations (alpha beta total):
# 0  O:   4.46410162   4.46410162   8.92820323
#per-atom occupations (alpha beta total):
# 0  O:   4.46410162   4.46410162   8.92820323
# 1  H:   0.70710678   0.70710678   1.41421356
# 2  H:   0.70710678   0.70710678   1.41421356


# 3390
# 0  H:   0.50000000   0.50000000   1.00000000
#per-atom occupations (alpha beta total):
# 0  O:   4.00000000   4.00000000   8.00000000
#per-atom occupations (alpha beta total):
# 0  O:   4.00000000   4.00000000   8.00000000
# 1  H:   0.50000000   0.50000000   1.00000000
# 2  H:   0.50000000   0.50000000   1.00000000
