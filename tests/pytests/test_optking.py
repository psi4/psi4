import pytest

from utils import *

import psi4

""" Test keyword passing implicitly and explicitly through optking. Test appropriate molecule setting. These tests are also
run in Optking's pytest suite. They are included here as an example of running optking through PsiAPI """

""" frozen, ranged, and external force optimizations """
pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param(
            {
                "name": "hf",
                "options": {"frozen_distance": "1 2 3 4", "scf_type": "pk"},
                "ref_ene": -150.781130356,
            },
            id="frozen_stre",
            marks=pytest.mark.quick
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"frozen_bend": "1 2 3 2 3 4", "scf_type": "pk"},
                "ref_ene": -150.786372411,
            },
            id="frozen_bend",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"frozen_dihedral": "1 2 3 4", "scf_type": "pk"},
                "ref_ene": -150.786766848,
            },
            id="frozen_dihedral",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"ranged_distance": "2 3 1.30 1.35"},
                "ref_ene": -150.7853238,
            },
            id="ranged_stre",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"ranged_bend": "(1 2 3 105.0 110.0) (2 3 4 105.0 110.0)"},
                "ref_ene": -150.7861769,
            },
            id="ranged_bend",
            marks=pytest.mark.quick
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"ranged_dihedral": "(1 2 3 4 100.0 110.0)"},
                "ref_ene": -150.7866419,
            },
            id="ranged_dihedral",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "ext_force_distance": "1 2 '-8.0*(x-0.950)' 3 4 '-8.0*(x-0.950)'"
                },
                "ref_ene": -150.786669,
            },
            id="ext_f_stre",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "ext_force_bend": "1 2 3 '-8.0*(x-105.0)' 2 3 4 '-8.0*(x-105.0)'"
                },
                "ref_ene": -150.786177,
            },
            id="ext_f_bend",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {"ext_force_dihedral": "1 2 3 4 '-8.0*(x-120.0)'"},
                "ref_ene": -150.786647,
            },
            id="ext_f_dihedral",
            marks=pytest.mark.quick
        ),
    ],
)
def test_constraints(inp):
    hooh = psi4.geometry(
        """ 
        H
        O 1 0.90
        O 2 1.40 1 100.0
        H 3 0.90 2 100.0 1 115.0
        no_com
        no_reorient
        """
    )

    psi4_options = {"basis": "cc-PVDZ", "g_convergence": "gau_tight"}
    psi4_options.update(inp["options"])
    psi4.set_options(psi4_options)
    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 3e-8})

    e = psi4.optimize(inp["name"])

    assert psi4.compare_values(e, inp["ref_ene"], 6)  # TEST


""" Dimers """
@pytest.mark.parametrize(
    "inp",
    [
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "frozen_cartesian": " 1 Xyz 4 xYz ",
                    "opt_coordinates": "cartesian",
                },
                "ref_ene": -150.7866491,
            },
            id="frozen_cart1",
            marks=pytest.mark.quick
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "frozen_cartesian": " 2 xyz 3 xyz ",
                    "opt_coordinates": "cartesian",
                },
                "ref_ene": -150.7866390,
            },
            id="frozen_cart2",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "frozen_cartesian": " 1 x 1 y 1 Z 4 x 4 Y 4 z ",
                    "opt_coordinates": "cartesian",
                },
                "ref_ene": -150.7866491,
            },
            id="frozen_cart3",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "frozen_cartesian": " 1 Xyz 4 xYz ",
                    "opt_coordinates": "redundant",
                },
                "ref_ene": -150.7866491,
            },
            id="frozen_cart4",
        ),
        pytest.param(
            {
                "name": "hf",
                "options": {
                    "ext_force_cartesian": "1 x '-2.0*(x-1.0)' 1 y '-2.0*(x-1.0)'"
                },
                "ref_ene": -150.7866742,
            },
            id="ext_f_cartesian",
            marks=pytest.mark.quick
        ),
    ])
def test_cart_constraints(inp):
    hooh = psi4.geometry(
        """
        H  0.90  0.80  0.5
        O  0.00  0.70  0.0
        O  0.00 -0.70  0.0
        H -0.90 -0.80  0.5
        no_com
        no_reorient
        symmetry c1
        """
        )

    psi4_options = {
           "basis": "cc-pvdz",
           "g_convergence": "gau_tight",
           "geom_maxiter": 20,
           "consecutive_backsteps": 1,
       }
    psi4_options.update(inp["options"])
    psi4.set_options(psi4_options)

    if psi4.core.get_global_option("orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 1e-8})

    e = psi4.optimize(inp["name"])

    assert psi4.compare_values(e, inp["ref_ene"], 6)  # TEST

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param(
            {"name": "mp2", "options": {}, "ref_ene": -257.4109749}, id="ne2_dimer"
        ),
        pytest.param(
            {
                "name": "mp2",
                "options": {
                    "frag_ref_atoms": """[[[1]], [[2]]]""",
                    "interfrag_dist_inv": True,
                },
                "ref_ene": -257.4109749,
            },
            id="ne2_dimer2",
        ),
        pytest.param(
            {
                "name": "mp2",
                "options": {
                    "interfrag_coords": str({
                        "Natoms per frag": [1, 1],
                        "A Frag": 1,
                        "A Ref Atoms": [[1]],
                        "A Label": "Ne atom 1",
                        "B Frag": 2,
                        "B Ref Atoms": [[2]],
                        "B Label": "Ne atom 2",
                    })
                },
                "ref_ene": -257.4109749,
            },
            id="ne2_dimer_explicit",
        ),
    ],
)
def test_dimers_ne2_long(inp):
    # Test from long distance start.
    ne2 = psi4.geometry(
        """
      0 1
      Ne  0.0  0.0  0.0
      --
      0 1
      Ne  3.0  0.0  0.0
      nocom
    """
    )

    psi4_options = {
        "basis": "aug-cc-pvdz",
        "geom_maxiter": 30,
        "frag_mode": "MULTI",
        "g_convergence": "gau_verytight",
    }
    psi4_options.update(inp["options"])
    psi4.set_options(psi4_options)

    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        psi4.set_options({"e_convergence": 9, "d_convergence": 3e-8})

    e = psi4.optimize(inp["name"])
    assert psi4.compare_values(e, inp["ref_ene"], 6)
