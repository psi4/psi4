import pytest
import psi4
import itertools
from .utils import compare_integers, compare_values

pytestmark = pytest.mark.quick


def test_no_screening_schwarz():
    """Checks the number of shell quartets screened with Schwarz screening.
    No shell quartets should be screened with a threshold of 0.0"""

    mol = psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)
    psi4.set_options({ "ints_tolerance" : 0.0 ,
                       "screening" : "schwarz" })

    basis = psi4.core.BasisSet.build(mol, target='cc-pVDZ')
    factory = psi4.core.IntegralFactory(basis)
    eri = factory.eri(0)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not eri.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets Schwarz Screened, Cutoff 0')


def test_no_screening_csam():
    """Checks the number of shell quartets screened with CSAM screening.
    No shell quartets should be screened with a threshold of 0.0"""

    mol = psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)
    psi4.set_options({ "ints_tolerance" : 0.0,
                       "screening" : "csam" })

    basis = psi4.core.BasisSet.build(mol, target='cc-pVDZ')
    factory = psi4.core.IntegralFactory(basis)
    eri = factory.eri(0)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not eri.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets CSAM Screened, Cutoff 0')


def test_schwarz_vs_csam_quartets():
    """Checks difference between the number of shell quartets screened with Schwarz and CSAM screening. 
    CSAM is strictly tighter than Schwarz and should screen at least all of the same shell pairs.
    Default threshhold of 1.0E-12 is used"""

    mol = psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)

    psi4.set_options({ "ints_tolerance" : 1e-12})
    basis = psi4.core.BasisSet.build(mol, target='DZ')

    factory = psi4.core.IntegralFactory(basis)
    psi4.set_options({ "ints_tolerance" : 1e-12,
                       "screening" : 'csam' })
    eriCSAM = factory.eri(0)
    psi4.set_options({ "screening" : 'schwarz', 'integral_package': 'libint2' })
    eriSchwarz = factory.eri(0)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count_both = 0
    screen_count_csam = 0
    screen_count_schwarz = 0
    screen_count_none = 0

    for m, n, r, s in quartets:
        screen_schwarz = not eriSchwarz.shell_significant(m, n, r, s)
        screen_csam = not eriCSAM.shell_significant(m, n, r, s)

        if screen_schwarz and screen_csam:
            screen_count_both += 1
        elif screen_csam:
            screen_count_csam += 1
        elif screen_schwarz:
            screen_count_schwarz += 1
        else:
            screen_count_none += 1
    assert compare_integers(75072, screen_count_both, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(2336, screen_count_csam, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(0, screen_count_schwarz, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(27568, screen_count_none, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')


def test_schwarz_vs_csam_energy():
    """Checks difference in Hartree-Fock energy between no screening and CSAM screening, which should be
    insignificant. """

    psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)

    psi4.set_options({'scf_type' : 'direct',
                      'd_convergence' : 1e-12,
                      'screening' : 'schwarz',
                      'ints_tolerance' : 1.0e-12 })
    e_schwarz = psi4.energy('hf/DZ')

    psi4.core.clean()

    psi4.set_options({'scf_type' : 'direct',
                      'd_convergence' : 1e-12,
                      'screening' : 'csam',
                      'ints_tolerance' : 1.0e-12 })
    e_csam = psi4.energy('hf/DZ')

    assert compare_values(e_schwarz, e_csam, 11, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
