import numpy as np
from numpy import linalg as LA
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

@pytest.fixture
def space1():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5

        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ-f12', molecule=h2o, return_wfn=True)
    space1 = wfn.alpha_orbital_space('p', 'SO', 'ALL')
    return space1

@pytest.fixture
def space2(space1):
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5

        symmetry c1
    """)

    keys = ["BASIS","CABS_BASIS"]
    targets = ["CC-PVDZ-F12","CC-PVDZ-F12-OPTRI"]
    roles = ["ORBITAL","F12"]
    others = ["CC-PVDZ-F12", "CC-PVDZ-F12"]

    combined = psi4.driver.qcdb.libmintsbasisset.BasisSet.pyconstruct_combined(h2o.save_string_xyz(), keys, targets, roles, others)
    combined = psi4.core.BasisSet.construct_from_pydict(h2o, combined, combined["puream"])

    ribs = psi4.core.OrbitalSpace.build_ri_space(combined)
    space2 = psi4.core.OrbitalSpace.build_cabs_space(space1, ribs)
    return space2

def test_OBS_orthonormality(space1, space2):
    obs = space1.basisset()
    cabs = space2.basisset()
    C_obs = np.array(space1.C())
    C_cabs = np.array(space2.C())
    mints = psi4.core.MintsHelper(obs)

    S_obs = np.array(mints.ao_overlap())
    S_mo_obs_obs = LA.multi_dot([C_obs.T, S_obs, C_obs])
    np.testing.assert_allclose(np.dot(S_mo_obs_obs.T, S_mo_obs_obs),\
                               np.eye(48), rtol=1e-05, atol=1e-07)

def test_CABS_orthonormality(space1, space2):
    obs = space1.basisset()
    cabs = space2.basisset()
    C_obs = np.array(space1.C())
    C_cabs = np.array(space2.C())
    mints = psi4.core.MintsHelper(cabs)

    S_cabs = np.array(mints.ao_overlap())
    S_mo_cabs_cabs = LA.multi_dot([C_cabs.T, S_cabs, C_cabs])
    np.testing.assert_allclose(np.dot(S_mo_cabs_cabs.T, S_mo_cabs_cabs),\
                               np.eye(110), rtol=1e-05, atol=1e-07)

def test_OBS_CABS_orthonormality(space1, space2):
    obs = space1.basisset()
    cabs = space2.basisset()
    C_obs = np.array(space1.C())
    C_cabs = np.array(space2.C())
    mints = psi4.core.MintsHelper(cabs)

    S_obs_cabs = np.array(mints.ao_overlap(obs, cabs))
    S_mo_obs_cabs = LA.multi_dot([C_obs.T, S_obs_cabs, C_cabs])
    np.testing.assert_allclose(np.dot(S_mo_obs_cabs.T, S_mo_obs_cabs),\
                               np.zeros((110, 110)), rtol=1e-05, atol=1e-07)

def test_CABS_OBS_orthonormality(space1, space2):
    obs = space1.basisset()
    cabs = space2.basisset()
    C_obs = np.array(space1.C())
    C_cabs = np.array(space2.C())
    mints = psi4.core.MintsHelper(cabs)

    S_cabs_obs = np.array(mints.ao_overlap(cabs, obs))
    S_mo_cabs_obs = LA.multi_dot([C_cabs.T, S_cabs_obs, C_obs])
    np.testing.assert_allclose(np.dot(S_mo_cabs_obs.T, S_mo_cabs_obs),\
                               np.zeros((48, 48)), rtol=1e-05, atol=1e-07)