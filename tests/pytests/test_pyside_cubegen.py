"""
Test Python-side generation of cube files
"""
import pytest
import psi4
from utils import compare_cubes

pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_pyside_cubegen():
    mol = psi4.geometry("""
        O 0 0 0
        H 0 0 1.795239827225189
        H 1.693194615993441 0 -0.599043184453037
        symmetry c1
        units au
        """)

    psi4.core.be_quiet()
    psi4.set_options({'basis': "sto-3g",
                      'scf_type': 'pk',
                      'cubeprop_tasks': ['density', 'orbitals']})
    scf_e, wfn = psi4.energy('SCF', return_wfn=True, molecule=mol)
    psi4.cubeprop(wfn)

    cubegen = psi4.core.CubeProperties(wfn)
    Dtot = wfn.Da()
    Dtot.add(wfn.Db())
    cubegen.compute_density(Dtot, "Dtot")

    alpha_orbitals = wfn.Ca_subset("AO", "OCC").np
    # select the three highest occupied orbitals
    occs = alpha_orbitals[:, -3:]
    occs_pm = psi4.core.Matrix.from_array(occs)
    cubegen.compute_orbitals(occs_pm, [0, 2], ["1", "3"], "orbital")

    assert compare_cubes("Dt.cube", "Dtot.cube")
    assert compare_cubes("Psi_a_5_5-A.cube", "orbital_3_3.cube")
    assert compare_cubes("Psi_a_3_3-A.cube", "orbital_1_1.cube")
