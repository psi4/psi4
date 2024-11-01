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
