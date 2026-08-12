import pytest

import psi4
from os.path import isfile

pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_serialize_wfn():
    """wfn serialization"""

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({'basis': "cc-pVDZ"})
    _, scf_wfn = psi4.energy('scf', return_wfn=True)

    # write the wavefunction to file
    scf_wfn.to_file('pytest_wfn')

    # alternatively store the dict representation of the wavefunction in memory
    wfn_dict = scf_wfn.to_file()

    assert isfile("pytest_wfn.npy")

    # read wavefunction from file
    wfn_new = psi4.core.Wavefunction.from_file('pytest_wfn')
    assert psi4.compare_wavefunctions(scf_wfn, wfn_new, label='Serialization Check(disk)')

    # make a wavefunction from the dict
    wfn_new2 = psi4.core.Wavefunction.from_file(wfn_dict)
    assert psi4.compare_wavefunctions(scf_wfn, wfn_new2, label='Serialization Check(dict)')


def test_restart_scf_serial_wfn():
    """scf restart from wfn file"""

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({'basis': "cc-pVDZ"})
    _, scf_wfn = psi4.energy('scf', return_wfn=True)
    scf_wfn.to_file('pytest_wfn')
    psi4.core.clean()
    psi4.set_options({'maxiter': 1})
    psi4.energy('scf', restart_file='pytest_wfn')
    psi4.core.clean()
    assert psi4.compare_values(-76.0266327341067125, psi4.variable('SCF TOTAL ENERGY'), 6, 'SCF energy')


def test_restart_scf_orbital_file():
    """wfn serialization"""

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({'basis': "cc-pVDZ"})
    _, scf_wfn = psi4.energy('scf', return_wfn=True, write_orbitals='my_mos')
    psi4.core.clean()
    psi4.set_options({'maxiter': 1})
    scf_wfn = psi4.energy('scf', restart_file='my_mos')
    assert psi4.compare_values(-76.0266327341067125, psi4.variable('SCF TOTAL ENERGY'), 6, 'SCF energy')


def test_restart_scf_guess_C_kwarg():
    """SCF restart via the direct guess_C kwarg: plain NumPy arrays (or
    psi4.core.Matrix), handed straight to energy(), reproduce a converged
    calculation in one iteration."""

    h2o = psi4.geometry("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    psi4.set_options({'basis': "cc-pVDZ"})
    e0, wfn0 = psi4.energy('scf', return_wfn=True)

    # Symmetric, restricted reference: a single Matrix suffices (alpha mirrored
    # to beta internally). Passed as a psi4.core.Matrix directly (not .to_array()'d
    # to a bare tuple of per-irrep blocks), since a bare tuple means (Ca, Cb).
    Ca_occ = wfn0.Ca_subset("SO", "OCC")

    psi4.core.clean()
    psi4.set_options({'maxiter': 1})
    e1 = psi4.energy('scf', guess_C=Ca_occ)
    assert psi4.compare_values(e0, e1, 6, 'SCF energy (guess_C kwarg, restricted)')

    # An unrestricted reference needs a (Ca, Cb) tuple.
    psi4.core.clean()
    psi4.set_options({'basis': "cc-pVDZ", 'reference': 'uhf', 'maxiter': 100})
    o2 = psi4.geometry("0 3\nO\nO 1 1.2\n")
    e0u, wfn0u = psi4.energy('scf', molecule=o2, return_wfn=True)
    Ca_occ_u = wfn0u.Ca_subset("SO", "OCC").to_array()
    Cb_occ_u = wfn0u.Cb_subset("SO", "OCC").to_array()

    psi4.core.clean()
    psi4.set_options({'basis': "cc-pVDZ", 'reference': 'uhf', 'maxiter': 2})
    e1u = psi4.energy('scf', molecule=o2, guess_C=(Ca_occ_u, Cb_occ_u))
    assert psi4.compare_values(e0u, e1u, 6, 'SCF energy (guess_C=(Ca, Cb) tuple, UHF)')
