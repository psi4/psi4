import numpy as np
import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

@pytest.fixture(params=[(0, 1.0e-6), (2, 1.0e-6), (0, 1.0e-2), (2, 1.0e-2)],
                ids=["all_orbitals", "obs_lindep", "ribs_lindep", "both_lindep"])
def spaces(request):
    """Returns [OBS, CABS, RIBS]. ndrop emulates linear dependencies removed from the
    orbital basis, so that OBS has fewer orbitals than basis functions; ri_tol controls
    how many linear dependencies are removed from the RI basis (issue #3498)."""
    ndrop, ri_tol = request.param

    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5

        symmetry c1
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ-f12', molecule=h2o, return_wfn=True)
    obs = wfn.alpha_orbital_space('p', 'SO', 'ALL')

    if ndrop:
        C = psi4.core.Matrix.from_array(np.array(obs.C())[:, :-ndrop])
        obs = psi4.core.OrbitalSpace(obs.id(), obs.name(), C, obs.basisset(), obs.integral())

    keys = ["BASIS","CABS_BASIS"]
    targets = ["CC-PVDZ-F12","CC-PVDZ-F12-OPTRI"]
    roles = ["ORBITAL","F12"]
    others = ["CC-PVDZ-F12", "CC-PVDZ-F12"]

    combined = psi4.driver.qcdb.libmintsbasisset.BasisSet.pyconstruct_combined(h2o.save_string_xyz(), keys, targets, roles, others)
    combined = psi4.core.BasisSet.construct_from_pydict(h2o, combined, combined["puream"])

    ribs = psi4.core.OrbitalSpace.build_ri_space(combined, ri_tol)
    cabs = psi4.core.OrbitalSpace.build_cabs_space(obs, ribs)

    return [obs, cabs, ribs]

@pytest.mark.parametrize("idx1,idx2", [(0, 0), (0, 1), (1, 0), (1, 1)])
def test_orthonormality(spaces, idx1, idx2):
    s1 = spaces[idx1]
    s2 = spaces[idx2]
    bs1 = s1.basisset()
    bs2 = s2.basisset()
    C1 = np.array(s1.C())
    C2 = np.array(s2.C())
    mints = psi4.core.MintsHelper(bs1)
    S_ao = np.array(mints.ao_overlap(bs1, bs2))
    S_mo = np.linalg.multi_dot([C1.T, S_ao, C2])

    if idx1 != idx2:
        np.testing.assert_allclose(S_mo, np.zeros((C1.shape[1], C2.shape[1])), rtol=1e-05, atol=1e-07)
    else:
        np.testing.assert_allclose(S_mo, np.eye(C1.shape[1]), rtol=1e-05, atol=1e-07)

def test_cabs_dimension(spaces):
    obs, cabs, ribs = spaces

    # CABS is the part of the RI space left over after projecting out the orbitals
    assert cabs.dim()[0] == ribs.dim()[0] - obs.dim()[0]
    assert cabs.C().shape[0] == ribs.basisset().nbf()
