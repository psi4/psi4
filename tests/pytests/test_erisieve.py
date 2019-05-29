import pytest
import psi4
import itertools
from .utils import compare_integers

pytestmark = pytest.mark.quick


def test_no_screening_schwarz():
    """Checks the number of shell quartets screened with Schwarz screening.
    No shell quartets should be screened with a threshold of 0.0"""

    psi4.geometry("""
      Ne 0.0 0.0 0.0
      Ne 4.0 0.0 0.0
      Ne 8.0 0.0 0.0
    """)

    _, wfn = psi4.energy('hf/cc-pvdz', return_wfn=True)
    basis = wfn.basisset()
    sieve_schwarz = psi4.core.ERISieve(basis, 0.0, False)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not sieve_schwarz.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets Schwarz Screened, Cutoff 0')


def test_no_screening_csam():
    """Checks the number of shell quartets screened with CSAM screening.
    No shell quartets should be screened with a threshold of 0.0"""

    psi4.geometry("""
      Ne 0.0 0.0 0.0
      Ne 4.0 0.0 0.0
      Ne 8.0 0.0 0.0
    """)

    _, wfn = psi4.energy('hf/cc-pvdz', return_wfn=True)
    basis = wfn.basisset()
    sieve_csam = psi4.core.ERISieve(basis, 0.0, True)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not sieve_csam.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets CSAM Screened, Cutoff 0')


def test_schwarz_vs_csam():
    """Checks difference between the number of shell quartets screened with Schwarz and CSAM screening. 
    CSAM is strictly tighter than Schwarz and should screen at least all of the same shell pairs.
    Default threshhold of 1.0E-12 is used"""

    psi4.geometry("""
      Ne 0.0 0.0 0.0
      Ne 4.0 0.0 0.0
      Ne 8.0 0.0 0.0
    """)

    _, wfn = psi4.energy('hf/cc-pvdz', return_wfn=True)
    basis = wfn.basisset()
    sieve_schwarz = psi4.core.ERISieve(basis, 1.0e-12, False)
    sieve_csam = psi4.core.ERISieve(basis, 1.0e-12, True)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count_both = 0
    screen_count_csam = 0
    screen_count_schwarz = 0
    screen_count_none = 0

    for m, n, r, s in quartets:
        screen_schwarz = not sieve_schwarz.shell_significant(m, n, r, s)
        screen_csam = not sieve_csam.shell_significant(m, n, r, s)

        if screen_schwarz and screen_csam:
            screen_count_both += 1
        elif screen_csam:
            screen_count_csam += 1
        elif screen_schwarz:
            screen_count_schwarz += 1
        else:
            screen_count_none += 1

    assert compare_integers(75792, screen_count_both, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(1344, screen_count_csam, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(0, screen_count_schwarz, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
    assert compare_integers(27840, screen_count_none, 'Schwarz vs CSAM Screening, Cutoff 1.0e-12')
