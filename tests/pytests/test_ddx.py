import psi4
import pytest
import numpy as np

from psi4 import core

from utils import *
from addons import using, uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]

__geoms = {
    "benzene":
    """
    C -0.886454    0.338645    0.000000
    C  0.508706    0.338645    0.000000
    C  1.206244    1.546396    0.000000
    C  0.508590    2.754905   -0.001199
    C -0.886235    2.754827   -0.001678
    C -1.583836    1.546621   -0.000682
    H -1.436213   -0.613672    0.000450
    H  1.058214   -0.613868    0.001315
    H  2.305924    1.546476    0.000634
    H  1.058790    3.707048   -0.001258
    H -1.436357    3.707108   -0.002631
    H -2.683440    1.546804   -0.000862
    symmetry c1
    units angstrom
    """,
    "fcm":
    """
    C   -0.502049    0.594262    0.000000
    H   -0.145376    1.098660    0.873652
    H   -1.572049    0.594275    0.000000
    Cl   0.084628    1.423927   -1.437034
    F   -0.052065   -0.678535    0.000000
    symmetry c1
    units angstrom
    """,
    "h2o":
    """
    O -0.966135 0.119522 0
    H -0.006135 0.119522 0
    H -1.286590 1.024458 0
    symmetry c1
    units angstrom
    """,
    "nh3":
    """
    N     -0.0000000001    -0.1040380466      0.0000000000
    H     -0.9015844116     0.4818470201     -1.5615900098
    H     -0.9015844116     0.4818470201      1.5615900098
    H      1.8031688251     0.4818470204      0.0000000000
    symmetry c1
    no_reorient
    no_com
    units bohr
    """,
    "h2":
    """
    H 0 0 0.5
    H 0 0 -0.5
    symmetry c1
    units angstrom
    """,
    "h":
    """
    H 0 0 0
    1 1
    symmetry c1
    units angstrom
    """,
}


def _base_test_fock(fock_term, density_matrix, eps=1e-4, tol=1e-6):
    E, V = fock_term(density_matrix)

    perturbation = np.random.random(density_matrix.np.shape)
    perturbation = perturbation + perturbation.T
    delta = np.sum(V.np * perturbation)

    Em, _ = fock_term(core.Matrix.from_array(density_matrix.np - eps * perturbation))
    Ep, _ = fock_term(core.Matrix.from_array(density_matrix.np + eps * perturbation))
    delta_ref = (Ep - Em) / (2 * eps)

    assert abs(delta - delta_ref) < 1e-6


@pytest.mark.quick
@uusing("ddx")
@pytest.mark.parametrize("inp", [
   pytest.param({
       "geom": __geoms["h"],
       "dm": core.Matrix.from_array(np.array([[0.]])),
       "radii": [1.0],  # Angstrom
       "ref": -0.2645886054599999,  # from Gaussian
   }, id='h'),
   pytest.param({
       "geom": __geoms["h2"],
       "dm": core.Matrix.from_array(0.6682326961201372 * np.ones((2, 2))),
       "radii": [1.5873, 1.5873],  # Angstrom
       "ref": -0.0002016948,  # from Gaussian
   }, id='h2'),
])
def test_ddx_fock_build(inp):
    """
    Tests COSMO energy against reference from Gaussian
    and internal consistency of Fock matrix versus the energy.
    """
    from psi4.driver.procrouting.solvent import ddx

    psi4.set_options({
        "df_scf_guess": False,
        #
        "ddx_model": "cosmo",
        "ddx_solvent_epsilon": 1e8,
        "ddx_eta": 0,
        "ddx_lmax": 3,
        "ddx_n_lebedev": 302,
        "ddx_radii": inp["radii"],
        "ddx_solute_spherical_points": 350,
        "ddx_solute_radial_points": 99,
    })

    # build the DDX object to test
    mol = psi4.geometry(inp["geom"])
    basis = psi4.core.BasisSet.build(mol, "BASIS", "sto-3g")
    ddx_options = ddx.get_ddx_options(mol)
    ddx_state = ddx.DdxInterface(mol, ddx_options, basis)

    def get_EV(density_matrix):
        return ddx_state.get_solvation_contributions(density_matrix)

    _base_test_fock(get_EV, inp["dm"])

    E, _ = get_EV(inp["dm"])
    assert compare_values(inp["ref"], E, atol=1e-16, rtol=1e-2)

@pytest.mark.quick
@uusing("ddx")
@pytest.mark.parametrize("inp", [
    pytest.param({
        "geom": __geoms["h2o"],
        "ddx": {"solvent_epsilon": 1e8, "eta": 0, },
        "ref": -75.5946789010,  # from Gaussian
        "solvation": -0.009402,
    }, id='h2o'),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"solvent_epsilon": 1e8, "eta": 0, },
        "ref": -594.993575419,  # from Gaussian
        "solvation": -0.006719,
    }, id='fcm'),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"solvent_epsilon": 2.0, "eta": 0, },
        "ref": -594.990420855,  # from Gaussian
        "solvation": -0.002964,
    }, id='fcmeps'),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"solvent_epsilon": 2.0, "eta": 0.2, },
        "ref": -594.990487330,  # from Gaussian
        "solvation": -0.003041,
    }, id='fcmepseta'),
    #
    pytest.param({
        "geom": __geoms["benzene"],
        "ddx": {"solvent_epsilon": 1e8, "eta": 0, },
        "ref": -229.420391688,  # from Gaussian
        "solvation": -0.005182,
    }, id='benzene'),
])
def test_ddx_rhf_reference(inp):
    mol = psi4.geometry(inp["geom"])
    psi4.set_options({
        "ddx": True,
        "basis": "3-21g",
        "guess": "core",
        "scf_type": "direct",
        #
        "ddx_model": "cosmo",
        "ddx_radii_set": "uff",
        "ddx_lmax": 3,
        "ddx_n_lebedev": 302,
        "ddx_solute_spherical_points": 302,
        "ddx_solute_radial_points": 75,
        "ddx_shift": 0.0,
    })
    for key in inp["ddx"].keys():
        psi4.set_options({"ddx_" + key: inp["ddx"][key]})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    ddx_e = wfn.scalar_variable("dd solvation energy")

    assert compare_values(inp["solvation"], ddx_e, 4, "DDX solvation energy versus Gaussian")
    assert compare_values(inp["ref"], scf_e, 4, "Total SCF energy with DDX versus Gaussian")


@pytest.mark.quick
@uusing("ddx")
@pytest.mark.parametrize("inp", [
    pytest.param({
        "geom": __geoms["h2o"],
        "basis": "cc-pvdz",
        "ddx": {"model": "cosmo", "solvent": "water", "radii_set": "bondi", },
        "ref": -76.0346355428018,
    }, id='h2o-cosmo'),
    pytest.param({
        "geom": __geoms["fcm"],
        "basis": "cc-pvdz",
        "ddx": {"model": "pcm", "solvent": "water", "radii_set": "uff", },
        "ref": -597.9718942424215,
    }, id='fcm-pcm'),
])
def test_ddx_rhf_consistency(inp):
    mol = psi4.geometry(inp["geom"])
    psi4.set_options({
        "ddx": True,
        "basis": inp["basis"],
        #
        "ddx_lmax": 10,
        "ddx_n_lebedev": 302,
        "ddx_solute_spherical_points": 302,
        "ddx_solute_radial_points": 75,
    })
    for key in inp["ddx"].keys():
        psi4.set_options({"ddx_" + key: inp["ddx"][key]})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    assert compare_values(inp["ref"], scf_e, 9, "Total SCF energy with DDX versus reference data")


@uusing("ddx")
@pytest.mark.parametrize("scf_type", ["pk", "out_of_core", "direct", "df", "cd"])
def test_ddx_eri_algorithms(scf_type):
    mol = psi4.geometry(__geoms["nh3"])
    psi4.set_options({
        "scf_type": scf_type,
        "basis": "6-31g",
        "ddx": True,
        "ddx_model": "pcm",
        "ddx_solvent": "water",
        "ddx_radii_set": "uff",
    })
    ref = -56.1715394
    scf_e = psi4.energy('SCF')
    assert compare_values(ref, scf_e, 4, "Total SCF energy with DDX versus reference data")
