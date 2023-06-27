import numpy as np
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

@pytest.fixture
def spaces():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5

        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ-f12', molecule=h2o, return_wfn=True)
    obs = wfn.alpha_orbital_space('p', 'SO', 'ALL')

    keys = ["BASIS","CABS_BASIS"]
    targets = ["CC-PVDZ-F12","CC-PVDZ-F12-OPTRI"]
    roles = ["ORBITAL","F12"]
    others = ["CC-PVDZ-F12", "CC-PVDZ-F12"]

    combined = psi4.driver.qcdb.libmintsbasisset.BasisSet.pyconstruct_combined(h2o.save_string_xyz(), keys, targets, roles, others)
    combined = psi4.core.BasisSet.construct_from_pydict(h2o, combined, combined["puream"])

    ribs = psi4.core.OrbitalSpace.build_ri_space(combined)
    cabs = psi4.core.OrbitalSpace.build_cabs_space(obs, ribs)

    return [obs, cabs]

@pytest.mark.parametrize("o1,o2", [(0, 0), (0, 1), (1, 0), (1, 1)])
def test_orthonormality(spaces, o1, o2):
    s1 = spaces[o1]
    s2 = spaces[o2]
    bs1 = s1.basisset()
    bs2 = s2.basisset()
    C1 = np.array(s1.C())
    C2 = np.array(s2.C())
    mints = psi4.core.MintsHelper(bs1)
    S_ao = np.array(mints.ao_overlap(bs1, bs2))
    S_mo = np.linalg.multi_dot([C1.T, S_ao, C2])

    if o1 != o2:
        np.testing.assert_allclose(np.dot(S_mo.T, S_mo), np.zeros((C2.shape[1], C2.shape[1])), rtol=1e-05, atol=1e-07)
    else:
        np.testing.assert_allclose(np.dot(S_mo.T, S_mo), np.eye(C1.shape[1]), rtol=1e-05, atol=1e-07)