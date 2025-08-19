import pytest
import pathlib

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
        psi4.set_options({"e_convergence": 9, "d_convergence": 5e-8})

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


@pytest.mark.parametrize("option, test_val", [
    ('ACCEPT_SYMMETRY_BREAKING', False),
    ('ADD_AUXILIARY_BONDS', True),
    ('ALG_GEOM_MAXITER', 30),
    ('AUXILIARY_BOND_FACTOR', 1.2),
    ('BT_DX_CONV', 0.01),
    ('BT_DX_RMS_CHANGE_CONV', 0.01),
    ('BT_MAX_ITER', 10),
    ('BT_PINV_RCOND', 0.01),
    ('CART_HESS_READ', True),
    ('CONJUGATE_GRADIENT_TYPE', "POLAK"),
    ('CONSECUTIVE_BACKSTEPS', 30),
    ('COVALENT_CONNECT', 0.01),
    ('DYNAMIC_LVL', 4),
    ('DYNAMIC_LVL_MAX', 3),
    ('ENSURE_BT_CONVERGENCE', True),
    ('EXT_FORCE_BEND', "2 3 4 `SIN(X)`"),
    ('EXT_FORCE_CARTESIAN', "`1 XY SIN(X)`"),
    ('EXT_FORCE_DIHEDRAL', "`1 2 3 4 SIN(X)`"),
    ('EXT_FORCE_DISTANCE', "`1 2 SIN(X)`"),
    ('EXT_FORCE_OOFP', "`1 2 3 4 SIN(X)`"),
    ('FIX_VAL_NEAR_PI', 0.01),
    ('FLEXIBLE_G_CONVERGENCE', True),
    ('FRAG_MODE', "`MULTI`"),
    # ('FRAG_REF_ATOMS', ),
    ('FREEZE_ALL_DIHEDRALS', True),
    ('FREEZE_INTRAFRAG', True),
    ('FROZEN_BEND', "1 2 3"),
    ('FROZEN_CARTESIAN', "1 XY"),
    ('FROZEN_DIHEDRAL', "1 2 3 4"),
    ('FROZEN_DISTANCE', "1 2"),
    ('FROZEN_OOFP', "1 2 3 4"),
    ('FULL_HESS_EVERY', 6),
    ('GEOM_MAXITER', 7),
    ('G_CONVERGENCE', "CFOUR"),
    ('HESSIAN_FILE', pathlib.Path("random_file.hess")),
    ('HESS_UPDATE', "BOFILL"),
    ('HESS_UPDATE_DEN_TOL', 0.01),
    ('HESS_UPDATE_DQ_TOL', 0.01),
    ('HESS_UPDATE_LIMIT', False),
    ('HESS_UPDATE_LIMIT_MAX', 0.01),
    ('HESS_UPDATE_LIMIT_SCALE', 0.01),
    ('HESS_UPDATE_USE_LAST', 12),
    ('H_BOND_CONNECT', True),
    ('H_GUESS_EVERY', False),
    ('INCLUDE_OOFP', True),
    ('INTERFRAG_COLLINEAR_TOL', 0.01),
    # ('INTERFRAG_COORDS', ), tested elsewhere
    ('INTERFRAG_DIST_INV', True),
    ('INTERFRAG_HESS', "FISCHER_LIKE"),
    ('INTERFRAG_MODE', "PRINCIPAL_AXES"),
    ('INTERFRAG_STEP_LIMIT', 0.01),
    ('INTERFRAG_STEP_LIMIT_MAX', 0.01),
    ('INTERFRAG_STEP_LIMIT_MIN', 0.01),
    ('INTRAFRAG_HESS', "LINDH"),
    ('INTRAFRAG_STEP_LIMIT', 0.01),
    ('INTRAFRAG_STEP_LIMIT_MAX', 0.01),
    ('INTRAFRAG_STEP_LIMIT_MIN', 0.01),
    ('IRC_CONVERGENCE', -0.9),
    ('IRC_DIRECTION', "BACKWARD"),
    ('IRC_MODE', "CONFIRM"),
    ('IRC_POINTS', 30),
    ('IRC_STEP_SIZE', 0.01),
    ('LINEAR_ALGEBRA_TOL', 0.01),
    ('LINEAR_BEND_THRESHOLD', 0.01),
    ('LINESEARCH', True),
    ('LINESEARCH_STEP', 0.01),
    ('MAX_DISP_G_CONVERGENCE', 0.01),
    ('MAX_ENERGY_G_CONVERGENCE', 0.01),
    ('MAX_FORCE_G_CONVERGENCE', 0.01),
    ('OPT_COORDINATES', "CARTESIAN"),
    ('OPT_TYPE', "TS"),
    ('PRINT', 4),
    ('RANGED_BEND', "1 2 3 120 130"),
    ('RANGED_CARTESIAN', "1 XY 1.5 1.7"),
    ('RANGED_DIHEDRAL', "1 2 3 4 120 130"),
    ('RANGED_DISTANCE', "1 2 1.2 1.4"),
    ('RANGED_OOFP', "1 2 3 4 60.0 65.0"),
    ('REDUNDANT_EVAL_TOL', 0.01),
    ('RFO_FOLLOW_ROOT', True),
    ('RFO_NORMALIZATION_MAX', 0.01),
    ('RFO_ROOT', 4),
    ('RMS_DISP_G_CONVERGENCE', 0.01),
    ('RMS_FORCE_G_CONVERGENCE', 0.01),
    ('RSRFO_ALPHA_MAX', 0.01),
    ('SD_HESSIAN', 0.01),
    ('SIMPLE_STEP_SCALING', True),
    ('SMALL_BEND_FIX_THRESHOLD', 0.01),
    ('STEEPEST_DESCENT_TYPE', "BARZILAI_BORWEIN"),
    ('STEP_TYPE', "NR"),
    ('TEST_B', True),
    ('TEST_DERIVATIVE_B', True),
    ('UNFREEZE_DIHEDRALS', "1 2 3 4"),
    ('V3D_TORS_ANGLE_LIM', 0.01),
    ('V3D_TORS_COS_TOL', 0.01),
    ('WRITE_TRAJECTORY', True)]
)
def test_all_keywords(option, test_val):

    import optking
    params = optking.op.OptParams(**{option: test_val})
    processed = params.to_dict(by_alias=True)
    assert processed.get(option) == test_val
