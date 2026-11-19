import numpy as np
import pytest

from utils import compare_values
from addons import uusing

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

_WATER = """
0 1
O 0.000 0.000 0.117
H 0.000 0.757 -0.469
H 0.000 -0.757 -0.469
symmetry c1
"""


def _charge_like(scale, nbf, seed):
    """Linear potential of the perturbed charge density, V(D) = scale * M Tr(M D) (the form of PCM/PE/DDX)."""
    M = np.random.default_rng(seed).standard_normal((nbf, nbf))
    M = psi4.core.Matrix.from_array(M + M.T)

    def potential(D):
        V = M.clone()
        V.scale(scale * M.vector_dot(D))
        return V

    return potential


def _hx(wfn, xs, singlet=None):
    if singlet is None:
        return [np.asarray(m) for m in wfn.twoel_Hx(xs, True, "MO")]
    return [np.asarray(m) for m in wfn.twoel_Hx_full(xs, True, "MO", singlet)]


def test_external_cpscf_perturbations():
    """All registered perturbations enter; none for RHF triplets; UHF sees the total density change on both spins."""
    mol = psi4.geometry(_WATER)
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "d_convergence": 1e-10, "save_jk": True})
    _, rhf = psi4.energy("scf", molecule=mol, return_wfn=True)
    nbf, nocc = rhf.nso(), rhf.nalpha()
    x = psi4.core.Matrix.from_array(np.random.default_rng(3).standard_normal((nocc, nbf - nocc)))
    v1, v2 = _charge_like(1.0, nbf, 11), _charge_like(0.5, nbf, 12)

    base = _hx(rhf, [x])[0]
    rhf.set_external_cpscf_perturbation("A", v1)
    only_1 = _hx(rhf, [x])[0] - base
    rhf.set_external_cpscf_perturbation("B", v2)
    both = _hx(rhf, [x])[0] - base
    rhf.clear_external_cpscf_perturbations()
    rhf.set_external_cpscf_perturbation("B", v2)
    only_2 = _hx(rhf, [x])[0] - base
    assert np.abs(only_2).max() > 1.0
    assert compare_values(only_1 + only_2, both, 10, "RHF: two perturbations add")

    triplet_with = _hx(rhf, [x], singlet=False)[0]
    rhf.clear_external_cpscf_perturbations()
    triplet_without = _hx(rhf, [x], singlet=False)[0]
    assert compare_values(triplet_without, triplet_with, 10, "RHF triplet: no charge-type response")

    psi4.core.clean()
    psi4.set_options({"reference": "uhf"})
    _, uhf = psi4.energy("scf", molecule=mol, return_wfn=True)
    zero = psi4.core.Matrix.from_array(np.zeros((nocc, nbf - nocc)))
    base_a, base_b = _hx(uhf, [x, zero])
    uhf.set_external_cpscf_perturbation("A", v1)
    with_a, with_b = _hx(uhf, [x, zero])
    assert np.abs(with_a - base_a).max() > 1.0
    assert compare_values(with_a - base_a, with_b - base_b, 10, "UHF: alpha density change acts on both spins")


@uusing("pcmsolver")
def test_pcm_tdscf_uhf_equals_rhf_singlets_and_triplets():
    """For a closed-shell molecule the UHF-TDHF roots are the RHF singlet and triplet roots, also with PCM."""
    mol = psi4.geometry(_WATER)
    pcm = """
    Units = Angstrom
    Medium { SolverType = IEFPCM  Solvent = Water  Nonequilibrium = True }
    Cavity { RadiiSet = UFF  Type = GePol  Area = 0.3  Mode = Implicit }
    """
    nstates = 3

    def roots(reference, triplets="none", n=nstates):
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({"basis": "def2-svp", "scf_type": "pk", "e_convergence": 1e-10, "d_convergence": 1e-9,
                          "reference": reference, "tdscf_states": n, "tdscf_r_convergence": 1e-5,
                          "tdscf_maxiter": 200, "pcm": True})
        if reference == "rhf":
            psi4.set_options({"tdscf_triplets": triplets})
        psi4.pcm_helper(pcm)
        _, wfn = psi4.energy("td-hf", molecule=mol, return_wfn=True)
        return sorted(wfn.variable(f"TD-HF ROOT 0 -> ROOT {i} EXCITATION ENERGY") for i in range(1, n + 1))

    rhf = sorted(roots("rhf", "none") + roots("rhf", "only"))
    uhf = roots("uhf", n=2 * nstates)
    k = 2 * nstates - 2   # the highest RHF roots of each set need not be among the lowest UHF ones
    assert compare_values(np.array(rhf[:k]), np.array(uhf[:k]), 7, "UHF-TDHF roots = RHF singlets + triplets (PCM)")
