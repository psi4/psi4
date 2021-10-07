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

def test_no_screening_density():
    """Checks the number of shell quartets screened with Density screening
    No shell quartets should be screened with a threshold of 0.0"""
    mol = psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)
    psi4.set_options({ "ints_tolerance" : 0.0,
                       "scf_type" : "direct" })

    basis = psi4.core.BasisSet.build(mol, target='cc-pVDZ')
    factory = psi4.core.IntegralFactory(basis)
    eri = factory.eri(0)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not eri.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets Density Screened, Cutoff 0')

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

def test_schwarz_vs_density_quartets():
    """Checks difference between the number of shell quartets screened with Schwarz and Density screening. 
    Default threshhold of 1.0E-12 is used"""

    mol = psi4.geometry("""
        0 1
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561
        symmetry c1
        no_reorient
        no_com
    """)

    # Run a quick SCF job to get density matrix
    energy, wfn = psi4.energy('hf/DZ', return_wfn=True)

    psi4.set_options({"ints_tolerance" : 1e-12})
    basis = wfn.basisset()

    factory = psi4.core.IntegralFactory(basis)
    psi4.set_options({ "screening" : 'schwarz', 'integral_package': 'libint2' })
    eriSchwarz = factory.eri(0)

    psi4.set_options({"screening" : 'density' })
    eriDensity = factory.eri(0)
    Da = [wfn.Da()]
    # Update using the Density Matrix
    eriDensity.update_density(Da)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count_schwarz = 0
    screen_count_density = 0

    for m, n, r, s in quartets:
        screen_schwarz = not eriSchwarz.shell_significant(m, n, r, s)
        screen_density = not eriDensity.shell_significant(m, n, r, s)

        if screen_schwarz:
            screen_count_schwarz += 1
        if screen_density:
            screen_count_density += 1
    
    assert compare_integers(screen_count_schwarz, 14504, 'Schwarz Screened Ints Count, Cutoff 1.0e-12')
    assert compare_integers(screen_count_density, 19440, 'Density Screened Ints Count, Cutoff 1.0e-12')

def test_rhf_vs_uhf_screening():
    """Checks difference between the number of shell quartets screened with Density screening in RHF vs UHF. 
    Difference should be 0, mathematically. Default threshhold of 1.0E-12 is used"""

    mol = psi4.geometry("""
        0 1
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561
        symmetry c1
        no_reorient
        no_com
    """)

    psi4.set_options({ "screening" : "density",
                       "ints_tolerance" : 1e-12, 
                       "reference" : "rhf",
                       "integral_package" : 'libint2'})

    # Run an SCF to get RHF Wavefunction
    rhf_energy, rhf_wfn = psi4.energy('hf/DZ', return_wfn=True)
    rhf_basis = rhf_wfn.basisset()
    rhf_factory = psi4.core.IntegralFactory(rhf_basis)
    eriRhf = rhf_factory.eri(0)

    # Run an SCF to get UHF Wavefunction
    psi4.set_options( { "reference" : "uhf" })
    uhf_energy, uhf_wfn = psi4.energy('hf/DZ', return_wfn=True)
    uhf_basis = uhf_wfn.basisset()
    uhf_factory = psi4.core.IntegralFactory(uhf_basis)
    eriUhf = uhf_factory.eri(0)


    Drhf = [rhf_wfn.Da()]
    Duhf = [uhf_wfn.Da(), uhf_wfn.Db()]

    # Update using the Density Matrix
    eriRhf.update_density(Drhf)
    eriUhf.update_density(Duhf)

    shell_inds = range(rhf_basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count_rhf = 0
    screen_count_uhf = 0

    for m, n, r, s in quartets:
        screen_rhf = not eriRhf.shell_significant(m, n, r, s)
        screen_uhf = not eriUhf.shell_significant(m, n, r, s)

        if screen_rhf:
            screen_count_rhf += 1
        if screen_uhf:
            screen_count_uhf += 1
    
    assert compare_integers(screen_count_rhf, 19440, 'RHF Density Screened Ints Count, Cutoff 1.0e-12')
    assert compare_integers(screen_count_uhf, 19440, 'UHF Density Screened Ints Count, Cutoff 1.0e-12')

def test_schwarz_vs_density_energy():
    """Checks difference in Hartree-Fock energy between Schwarz and Density screening (with and without IFB), 
    which should be insignificant.
    """

    mol = psi4.geometry("""
        0 1
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561
        symmetry c1
        no_reorient
        no_com
    """)

    psi4.set_options({'scf_type' : 'direct',
                      'd_convergence' : 1e-6,
                      'df_scf_guess' : False,
                      'screening' : 'schwarz',
                      'ints_tolerance' : 1.0e-12 })
    e_schwarz = psi4.energy('hf/DZ')

    psi4.core.clean()

    psi4.set_options({'scf_type' : 'direct',
                      'd_convergence' : 1e-6,
                      'df_scf_guess' : False,
                      'screening' : 'density',
                      'incfock' : False,
                      'ints_tolerance' : 1.0e-12 })
    e_density = psi4.energy('hf/DZ')

    psi4.core.clean()

    psi4.set_options({'scf_type' : 'direct',
                      'd_convergence' : 1e-6,
                      'df_scf_guess' : False,
                      'screening' : 'density',
                      'incfock' : True,
                      'ints_tolerance' : 1.0e-12 })

    e_incfock = psi4.energy('hf/DZ')

    assert compare_values(e_schwarz, e_density, 10, 'Schwarz vs Density Screening, Cutoff 1.0e-12')
    assert compare_values(e_schwarz, e_incfock, 10, 'Schwarz vs Density Screening, Cutoff 1.0e-12')