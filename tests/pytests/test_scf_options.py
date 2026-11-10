import numpy as np
import pytest

import psi4

from addons import using, uusing
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
    if (psi4.core.get_option("SCF", "SOSCF")
            and psi4.core.get_option("scf", "second_order_orbital_optimizer_package")
            in ["OTR", "OPENTRUSTREGION"]):
        # This test expects a plain UHF to stop at -0.82648 but guess_mix to reach -0.99872 .
        # The former is a saddle point, and OpenTrustRegion's gradient perturbation is there to
        # leave one, so with it on the plain run descends to -0.99872 too and the test cannot
        # tell the two apart.
        psi4.set_options({"otr_grad_noise": False})

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
            "sap": -99.97314321248359,
        },
        "uhf": {
            "core": -88.99202896495355,
            "gwh": -94.35585144644561,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641045,
            "huckel": -99.20374047108815,
            "modhuckel": -99.17689049002892,
            "sap": -99.50145022872545,
        },
        "rohf": {
           "core": -88.99202896495353,
            "gwh": -94.35585144644557,
            "sad": -99.63941801281894,
            "sadno": -99.11882081641042,
            "huckel": -99.20374047108817,
            "modhuckel": -99.17689049002892,
            "sap": -99.50145022872545,
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
            "sad_frac_occ": True,
            "sad_spin_average": True,
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
            "sad_frac_occ": True,
            "sad_spin_average": True,
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


@pytest.mark.parametrize("soopkg", [
    pytest.param("internal", id="internal"),
    pytest.param("otr", id="otr", marks=using("otr")),
])
def test_soscf_maxiter_budget(soopkg):
    """MAXITER caps the whole SCF, not each stage of it.

    The second-order package takes over partway through, so its own iteration budget has to
    account for what the first-order iterations already spent -- and, for OpenTrustRegion, for
    the macro-iteration that re-evaluates the handoff state psi4 has already counted. Both
    packages should therefore need the same MAXITER to converge the same computation, and
    should report the same iteration count.
    """
    psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
      symmetry c1
    """)
    common = {"basis": "cc-pvdz", "scf_type": "pk", "soscf": True,
              "second_order_orbital_optimizer_package": soopkg, "fail_on_maxiter": True}

    # find the smallest budget that converges
    for maxiter in range(2, 20):
        psi4.set_options({**common, "maxiter": maxiter})
        try:
            psi4.energy("scf")
        except psi4.SCFConvergenceError:
            continue
        break
    else:
        pytest.fail("never converged")

    niter = int(psi4.variable("SCF ITERATIONS"))
    # the budget actually spent is what was reported, not one more or one less
    assert niter == maxiter, f"converged at MAXITER={maxiter} but reported {niter} iterations"

    # and one fewer must fail rather than quietly returning a shorter run
    psi4.set_options({**common, "maxiter": maxiter - 1})
    with pytest.raises(psi4.SCFConvergenceError):
        psi4.energy("scf")


@pytest.mark.scf
@uusing("otr")
def test_otr_precision_floor_polish():
    """OpenTrustRegion stops where the energy stops resolving; DIIS finishes the job.

    Asked for more precision than the objective function can resolve, OpenTrustRegion stops
    short of the thresholds and reports success. OTR_MAX_PRECISION_REACHED_HANDLING = DIIS
    hands the orbitals to first-order iterations, which converge on the orbital gradient
    instead and are not held to the same limit. Nothing in the regular suite reaches that
    floor, so ask for a threshold tight enough to put it there on purpose.
    """
    # methyloxirane, geometry as in cc31 (Pederson et al., CPL)
    psi4.geometry("""
0 1
C     0.210991841737     0.047663662717     0.483960886619
C    -0.960965312103     0.688519555303    -0.115848688934
C     1.558206749694     0.047460114159    -0.180267436891
O    -0.815251456906    -0.734986972268    -0.138532905591
H     0.223550547304    -0.045637654519     1.567373484970
H    -1.759483315918     1.063259112311     0.516330435381
H    -0.862256972568     1.177986074652    -1.079897880020
H     2.088617117792    -0.883541600456     0.026321487274
H     2.169362364888     0.872235726981     0.192209317060
H     1.455347299576     0.149768052462    -1.260354487995
noreorient
""")
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "reference": "rhf",
                      "soscf": True, "second_order_orbital_optimizer_package": "otr",
                      "e_convergence": 1.0e-13, "d_convergence": 1.0e-13,
                      "otr_max_precision_reached_handling": "diis"})

    e, wfn = psi4.energy("scf", return_wfn=True)

    # the floor was reached, so the polish is what finished this SCF
    assert wfn.opentrustregion_stopped_short(), \
        "expected OpenTrustRegion to stop at its precision floor for d_convergence 1e-13"
    assert compare_values(-189.513417237103, e, 9, "SCF energy after the first-order polish")

    # and it took a run of them, not a single step: OpenTrustRegion stops around a 1e-9
    # orbital gradient here and the polish carries it the remaining four orders
    assert int(psi4.variable("SCF ITERATIONS")) >= 10, \
        f'only {int(psi4.variable("SCF ITERATIONS"))} iterations; expected a multi-step polish'

    # FAIL declines to finish it, so the same request becomes non-convergence
    psi4.set_options({"otr_max_precision_reached_handling": "fail"})
    with pytest.raises(psi4.SCFConvergenceError):
        psi4.energy("scf")


@pytest.mark.scf
@uusing("otr")
def test_otr_polish_keeps_its_improvement():
    """A polish that runs out of iterations keeps its orbitals when they are the better ones.

    The polish is allowed to end short, and when it does, whichever of the two solutions has
    the smaller orbital gradient is the one left behind. Both have to be measured the same
    way for that comparison to mean anything, and OpenTrustRegion's cannot be measured after
    the driver diagonalizes: that makes the density a projector onto eigenvectors of the
    current Fock matrix, so the commutator falls to roundoff however unconverged the orbitals
    are, and the comparison then throws the polish's work away every time.
    """
    psi4.geometry("""
0 2
N
O 1 1.158
symmetry c1
""")
    # the SCF of dfmp2-grad2, whose D_CONVERGENCE is tight enough that OpenTrustRegion hits
    # its floor around iteration 7. Its polish converges at iteration 15, so stop two short
    # of that: MAXITER is what ends a polish, and this is a polish ended by it.
    psi4.set_options({"basis": "cc-pvdz", "df_basis_scf": "cc-pvdz-jkfit",
                      "scf_type": "df", "guess": "sad", "reference": "uhf",
                      "soscf": True, "second_order_orbital_optimizer_package": "otr",
                      "e_convergence": 1.0e-10, "d_convergence": 1.0e-10, "maxiter": 13,
                      "otr_max_precision_reached_handling": "diis"})

    with pytest.raises(psi4.SCFConvergenceError) as exc:
        psi4.energy("scf")

    # the SCF is unconverged on its orbital gradient, so the gradient it reports has to be a
    # real measure of these orbitals -- above D_CONVERGENCE, not the roundoff of a density
    # built from the eigenvectors of its own Fock matrix
    assert exc.value.d_conv > 1.0e-10, \
        f"reported Dnorm {exc.value.d_conv:.2e} is below the D_CONVERGENCE it failed to reach"

    # and the polish improved on where OpenTrustRegion stopped, so its orbitals are the ones
    # kept; an equal value would mean OpenTrustRegion's were restored over better ones
    otr_entry_dnorm = exc.value.wfn.otr_outcome_[3]
    assert exc.value.d_conv < otr_entry_dnorm, \
        (f"polish ended at {exc.value.d_conv:.2e} against OpenTrustRegion's "
         f"{otr_entry_dnorm:.2e}; the better orbitals were not the ones kept")
