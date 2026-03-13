"""
This file tests the ZORA spin-orbit coupling 2p and 3d orbital energy splitting.
"""
import numpy as np
import pytest
import psi4
from addons import uusing

pytestmark = [pytest.mark.api, pytest.mark.quick]

@pytest.fixture
def mo_energies():
    psi4.geometry("""
        0 1
            Kr 0 0 0
        symmetry C1
    """)

    psi4.set_options({
        'reference': 'cghf',
        'scf_type': 'pk',
        'basis': '3-21g',
        'e_convergence': 1e-9,
        'relativistic': 'zora',
        'spin_orbit_coupling': False
    })

    e, wfn = psi4.energy('scf', return_wfn=True)

    psi4.set_options({
        'spin_orbit_coupling': True
    })

    e_so, wfn_so = psi4.energy('scf', return_wfn=True)

    def get_2p3d_energies(wfn):
        """Returns 2p and 3d energies of single atom."""
        # GHF energies are split into alpha/beta due to factory limitations.
        mo_energy_a = wfn.epsilon_a().to_array()
        mo_energy_b = wfn.epsilon_b().to_array()
        mo_energies = np.sort(np.append(mo_energy_a, mo_energy_b))

        mo_e_2p = mo_energies[4:10]
        mo_e_3d = mo_energies[18:28]
        return mo_e_2p, mo_e_3d

    mo_e_2p, mo_e_3d = get_2p3d_energies(wfn)
    mo_e_2p_so, mo_e_3d_so = get_2p3d_energies(wfn_so)

    return (mo_e_2p, mo_e_3d, mo_e_2p_so, mo_e_3d_so)

@uusing("einsums")
def test_2p_splitting(mo_energies):
    """
    Demonstrates the lack of 2p orbital splitting without SOC and shows splitting with SOC.
    """
    mo_e_2p, _, mo_e_2p_so, _ = mo_energies
    assert np.allclose(mo_e_2p, mo_e_2p[0])
    assert not np.allclose(mo_e_2p_so, mo_e_2p_so[0])

@uusing("einsums")
def test_3d_splitting(mo_energies):
    """
    Demonstrates the lack of 3d orbital splitting without SOC and shows splitting with SOC.
    """
    _, mo_e_3d, _, mo_e_3d_so = mo_energies
    assert np.allclose(mo_e_3d, mo_e_3d[0])
    assert not np.allclose(mo_e_3d_so, mo_e_3d_so[0])

@uusing("einsums")
def test_splitting(mo_energies):
    _, _, mo_e_2p_so, mo_e_3d_so = mo_energies
    assert np.allclose(mo_e_2p_so[1], mo_e_2p_so[0])  # 2p 1/2
    assert np.allclose(mo_e_2p_so[2:], mo_e_2p_so[2]) # 2p 3/2

    assert np.allclose(mo_e_3d_so[:4], mo_e_3d_so[0]) # 3d 3/2
    assert np.allclose(mo_e_3d_so[4:], mo_e_3d_so[4]) # 3d 5/2




