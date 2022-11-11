import psi4
import pytest
import numpy as np

from psi4 import core

from utils import compare_values
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
    "methane":
    """
    C  0.268924   -0.099602    0.000000
    H  0.625579   -1.108412    0.000000
    H  0.625597    0.404797    0.873652
    H  0.625597    0.404797   -0.873652
    H -0.801076   -0.099588    0.000000
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
    """
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
   }, id='h', marks=using('ddx')),
   pytest.param({
       "geom": __geoms["h2"],
       "dm": core.Matrix.from_array(0.6682326961201372 * np.ones((2, 2))),
       "radii": [1.5873, 1.5873],  # Angstrom
       "ref": -0.0002016948,  # from Gaussian
   }, id='h2', marks=using('ddx')),
])
def test_ddx_fock_build(inp):
    """
    Tests COSMO energy against reference from Gaussian
    and internal consistency of Fock matrix versus the energy.
    """
    psi4.set_options({
        # Note: DFT grids used for DDX numerical integration
        "dft_spherical_points": 350,
        "dft_nuclear_scheme": "Becke",
        "dft_pruning_scheme": "robust",
        "df_scf_guess": False,
        #
        "ddx__model": "cosmo",
        "ddx__solvent_epsilon": 1e8,
        "ddx__eta": 0,
        "ddx__lmax": 3,
        "ddx__n_lebedev": 302,
        "ddx__radii": inp["radii"],
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
        "ddx": {"model": "cosmo", "solvent_epsilon": 1e8, "eta": 0, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -75.5946789010,  # from Gaussian
        "solvation": -0.009402,
    }, id='h2o', marks=using('ddx')),
    #
    pytest.param({
        "geom": __geoms["methane"],
        "ddx": {"model": "cosmo", "solvent_epsilon": 1e8, "eta": 0, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -39.9764868732,  # from Gaussian
        "solvation": -0.000083,
    }, id='ch4', marks=using('ddx')),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"model": "cosmo", "solvent_epsilon": 1e8, "eta": 0, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -594.993575419,  # from Gaussian
        "solvation": -0.006719,
    }, id='fcm', marks=using('ddx')),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"model": "cosmo", "solvent_epsilon": 2.0, "eta": 0, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -594.990420855,  # from Gaussian
        "solvation": -0.002964,
    }, id='fcmeps', marks=using('ddx')),
    #
    pytest.param({
        "geom": __geoms["fcm"],
        "ddx": {"model": "cosmo", "solvent_epsilon": 2.0, "eta": 0.2, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -594.990487330,  # from Gaussian
        "solvation": -0.003041,
    }, id='fcmepseta', marks=using('ddx')),
    #
    pytest.param({
        "geom": __geoms["benzene"],
        "ddx": {"model": "cosmo", "solvent_epsilon": 1e8, "eta": 0, "lmax": 3,
                "n_lebedev": 302, "radii_set": "uff", },
        "ref": -229.420391688,  # from Gaussian
        "solvation": -0.005182,
    }, id='benzene', marks=using('ddx')),
])
def test_ddx_rhf_reference(inp):
    mol = psi4.geometry(inp["geom"])
    psi4.set_options({
        "ddx": True,
        "basis": "3-21g",
        "guess": "core",
        "scf_type": "direct",
    })
    for key in inp["ddx"].keys():
        psi4.set_options({"ddx__" + key: inp["ddx"][key]})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    ddx_e = wfn.get_variable("ddx energy")

    assert compare_values(inp["solvation"], ddx_e, 4, "DDX solvation energy versus Gaussian")
    assert compare_values(inp["ref"], scf_e, 4, "Total SCF energy with DDX versus Gaussian")


@pytest.mark.quick
@uusing("ddx")
@pytest.mark.parametrize("inp", [
    pytest.param({
        "geom": __geoms["h2o"],
        "basis": "cc-pvdz",
        "ddx": {"model": "cosmo", "solvent": "water", "eta": 0.05,
                "lmax": 10, "n_lebedev": 302, "radii_set": "bondi",
                "radii_scaling": 1.2},
        "solvent": "water",
        "ref": -76.0346373391875,
    }, id='h2o', marks=using('ddx')),
])
def test_ddx_rhf_consistency(inp):
    mol = psi4.geometry(inp["geom"])
    psi4.set_options({
        "ddx": True,
        "basis": inp["basis"],
        "guess": "core",
        "scf_type": "direct",
    })
    for key in inp["ddx"].keys():
        psi4.set_options({"ddx__" + key: inp["ddx"][key]})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    assert compare_values(inp["ref"], scf_e, 9, "Total SCF energy with DDX versus reference data")
