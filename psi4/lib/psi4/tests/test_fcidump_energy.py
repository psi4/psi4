import os
import pytest
import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_fcidump_scf_energy():
    """Compare FCIDUMP computed SCF energy against call to energy()"""

    Ne = psi4.geometry("""
      Ne 0 0 0
    """)

    psi4.set_options({'basis': 'cc-pVDZ',
                      'scf_type': 'pk',
                      'reference': 'uhf',
                      'd_convergence': 1e-8,
                      'e_convergence': 1e-8
                     })
    scf_e, scf_wfn = psi4.energy('scf', return_wfn=True)

    psi4.fcidump(scf_wfn, fname='FCIDUMP_SCF', oe_ints=['EIGENVALUES'])
    intdump = psi4.fcidump_from_file('FCIDUMP_SCF')
    e_dict = psi4.energies_from_fcidump(intdump)
    fcidump_e = e_dict['SCF TOTAL ENERGY']

    assert psi4.compare_values(scf_e, fcidump_e, 5, 'SCF energy')


def test_fcidump_mp2_energy():
    """Compare FCIDUMP computed MP2 energy against call to energy()"""

    Ne = psi4.geometry("""
      Ne 0 0 0
    """)

    psi4.set_options({'basis': 'cc-pVDZ',
                      'scf_type': 'pk',
                      'reference': 'uhf',
                      'd_convergence': 1e-8,
                      'e_convergence': 1e-8
                     })
    mp2_e, mp2_wfn = psi4.energy('mp2', return_wfn=True)

    psi4.fcidump(mp2_wfn, fname='FCIDUMP_MP2', oe_ints=['EIGENVALUES'])
    intdump = psi4.fcidump_from_file('FCIDUMP_MP2')
    e_dict = psi4.energies_from_fcidump(intdump)
    fcidump_e = e_dict['SCF TOTAL ENERGY'] + e_dict['MP2 CORRELATION ENERGY']

    assert psi4.compare_values(mp2_e, fcidump_e, 5, 'MP2 energy')
