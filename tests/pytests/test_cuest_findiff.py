import numpy as np
import psi4
import pytest

from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.dft, pytest.mark.findif]

# test_cuest.py pins cuEST gradients against stored values from regular Psi4. For
# DFT that yardstick is awkward: Psi4's CPU analytic DFT gradients omit the grid
# weight derivatives, so those stored DFT reference gradients had to be generated
# by finite difference in the first place.
#
# So instead of a reference, each case here computes the gradient twice from the
# same cuEST code path -- analytically (dertype=1) and as a 5-point finite
# difference of the cuEST energy (dertype=0) -- and requires the two to agree.
# That tests the derivative implementation rather than agreement with another
# program, and is insensitive to whichever grid and nuclear-partition conventions
# cuEST uses, since both halves of the comparison use them.
#
# C1 throughout, so finite difference displaces all 3N Cartesians and the two
# gradients compare element by element with no symmetry projection in between.
# Symmetry in the cuEST XC path is test_cuest_symmetry's job.

__geoms = {
    "water": """
    0 1
    O
    H 1 0.9584
    H 1 0.9584 2 104.45
    symmetry c1
    """,
    "water_cation": """
    1 2
    O
    H 1 0.9584
    H 1 0.9584 2 104.45
    symmetry c1
    """,
    "methylamine": """
    units bohr
    N   -1.2443656662    1.5116296128    1.1401094834
    C    1.2072325506    0.2976885350    0.9288948097
    H    2.1975245576    0.3690634856    2.7411987245
    H    2.3661826021    1.2428916800   -0.4976593391
    H    0.9782544525   -1.6841238010    0.3901528276
    H   -2.1240520112    1.4799955998   -0.5728893307
    H   -1.0022350753    3.3701184308    1.5836092757
    symmetry c1
    """,
    "methylamine_cation": """
    units bohr
    1 2
    N   -1.2443656662    1.5116296128    1.1401094834
    C    1.2072325506    0.2976885350    0.9288948097
    H    2.1975245576    0.3690634856    2.7411987245
    H    2.3661826021    1.2428916800   -0.4976593391
    H    0.9782544525   -1.6841238010    0.3901528276
    H   -2.1240520112    1.4799955998   -0.5728893307
    H   -1.0022350753    3.3701184308    1.5836092757
    symmetry c1
    """,
}


@pytest.mark.medlong
@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.parametrize("inp", [

    # LDA, GGA, global hybrid GGA
    pytest.param({"geom": "water",        "methodname": "svwn",   "reference": "rhf"}, id='water_rsvwn'),
    pytest.param({"geom": "water_cation", "methodname": "svwn",   "reference": "uhf"}, id='water_cation_usvwn'),
    pytest.param({"geom": "water",        "methodname": "blyp",   "reference": "rhf"}, id='water_rblyp'),
    pytest.param({"geom": "water_cation", "methodname": "blyp",   "reference": "uhf"}, id='water_cation_ublyp'),
    pytest.param({"geom": "water",        "methodname": "b3lyp",  "reference": "rhf"}, id='water_rb3lyp'),
    pytest.param({"geom": "water_cation", "methodname": "b3lyp",  "reference": "uhf"}, id='water_cation_ub3lyp'),

    # Meta-GGA hybrids: kinetic energy density (tau) terms
    pytest.param({"geom": "water",        "methodname": "m06",    "reference": "rhf"}, id='water_rm06'),
    pytest.param({"geom": "water_cation", "methodname": "m06",    "reference": "uhf"}, id='water_cation_um06'),
    pytest.param({"geom": "water",        "methodname": "pw6b95", "reference": "rhf"}, id='water_rpw6b95'),
    pytest.param({"geom": "water_cation", "methodname": "pw6b95", "reference": "uhf"}, id='water_cation_upw6b95'),

    # Meta-GGA with VV10 nonlocal correlation
    pytest.param({"geom": "water",        "methodname": "b97m-v", "reference": "rhf"}, id='water_rb97m-v'),
    pytest.param({"geom": "water_cation", "methodname": "b97m-v", "reference": "uhf"}, id='water_cation_ub97m-v'),

    # Range-separated (wK). cuEST computes no separate long-range exchange
    # derivative: the DF integral plan carries lrc_exchange_fraction and
    # lr_exchange_omega, so cuEST weights the short- and long-range derivatives
    # internally and folds both into the buffer that becomes
    # gradients_["Coulomb"]. Nothing on the Psi4 side rescales that buffer -- the
    # host passes only a sign and a spin factor -- so an error in the internal
    # weighting surfaces only as a mismatch against finite differences. The four
    # functionals cover four different (alpha, beta) splits, wpbe being the
    # lrc-only corner where alpha is 0 and there is no full-range exchange at all.
    pytest.param({"geom": "water",        "methodname": "wb97x",     "reference": "rhf"}, id='water_rwb97x'),          # 0.158 / 0.842
    pytest.param({"geom": "water_cation", "methodname": "wb97x",     "reference": "uhf"}, id='water_cation_uwb97x'),
    pytest.param({"geom": "water",        "methodname": "wpbe",      "reference": "rhf"}, id='water_rwpbe'),           # 0.0   / 1.0
    pytest.param({"geom": "water_cation", "methodname": "wpbe",      "reference": "uhf"}, id='water_cation_uwpbe'),
    pytest.param({"geom": "water",        "methodname": "wpbe0",     "reference": "rhf"}, id='water_rwpbe0'),          # 0.25  / 0.75
    pytest.param({"geom": "water",        "methodname": "cam-b3lyp", "reference": "rhf"}, id='water_rcam-b3lyp'),      # 0.19  / 0.46

    # Larger system. 7 atoms is 84 displaced energies per case.
    pytest.param({"geom": "methylamine",        "methodname": "b3lyp",  "reference": "rhf"}, id='methylamine_rb3lyp', marks=pytest.mark.long),
    pytest.param({"geom": "methylamine_cation", "methodname": "b3lyp",  "reference": "uhf"}, id='methylamine_cation_ub3lyp', marks=pytest.mark.long),
    pytest.param({"geom": "methylamine",        "methodname": "wb97x",  "reference": "rhf"}, id='methylamine_rwb97x', marks=pytest.mark.long),
    pytest.param({"geom": "methylamine_cation", "methodname": "wb97x",  "reference": "uhf"}, id='methylamine_cation_uwb97x', marks=pytest.mark.long),
])
def test_cuest_dft_findiff(inp, request):
    """cuEST analytic DFT gradient == 5-point finite difference of the cuEST energy."""
    psi4.core.set_num_threads(4)

    molecule = psi4.geometry(__geoms[inp["geom"]])

    psi4.set_options({
        # FD gradient options
        'points': 5,
        'disp_size': 0.005,
        'fd_project': False,
        'scf_type': 'df',
        'dft_nuclear_scheme': 'stratmann',
        'df_basis_scf': 'def2-universal-JKFIT',
        'basis': 'def2-svp',
        # The analytic gradient carries grid weight derivatives while the finite
        # difference lets the grid ride along with the nuclei. Those agree only
        # for a well converged grid.
        'dft_radial_points': 100,
        'dft_spherical_points': 590,
        'maxiter': 300,
        # Differencing energies over a 0.005 bohr step amplifies SCF noise by
        # ~1/h, so the energies must be converged well past the gradient
        # tolerance. Mixed precision would swamp the difference entirely.
        'e_convergence': 10,
        'd_convergence': 9,
        'puream': True,
        'use_cuest': True,
        'cuest_mixed_precision': False,
        'reference': inp['reference'],
    })

    G_analytic = psi4.gradient(inp['methodname'], molecule=molecule, dertype=1)
    psi4.core.clean()
    G_findif = psi4.gradient(inp['methodname'], molecule=molecule, dertype=0)

    assert psi4.compare_values(np.array(G_findif), np.array(G_analytic), 5e-5,
                               f'{request.node.callspec.id} analytic vs 5-point findif gradient')
