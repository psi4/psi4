import pytest
import psi4
import itertools
from utils import compare, compare_integers, compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


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
    """Checks difference between the number of shell quartets computed with Schwarz and Density screening. 
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

    # run schwarz screening calculation
    psi4.set_options({ 
        "scf_type": "direct", 
        "screening" : 'schwarz', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12, 
        "save_jk": True,
        "bench" : 1 

    })
    schwarz_energy, schwarz_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # run density screening calculation
    psi4.set_options({ 
        "scf_type": "direct", 
        "screening" : 'density', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12,
        "save_jk": True,
        "bench" : 1
    })
    density_energy, density_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # prep for comparing results to expected values
    schwarz_computed_shells = schwarz_wfn.jk().computed_shells_per_iter("Quartets")
    density_computed_shells = density_wfn.jk().computed_shells_per_iter("Quartets")

    schwarz_computed_shells_expected = [20290, 20290, 20290, 20290, 20290, 20290, 20290, 20290, 20290]
    density_computed_shells_expected = [13171, 19618, 19665, 19657, 19661, 19661, 19663, 19663, 19663]

    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert(len(schwarz_computed_shells_expected) == schwarz_wfn.iteration_+1)
    assert(len(density_computed_shells_expected) == density_wfn.iteration_+1)

    # actually compare results with expected values
    assert compare(schwarz_computed_shells_expected, schwarz_computed_shells, 'Schwarz Computed Shells Count, Cutoff 1.0e-12')
    assert compare(density_computed_shells_expected, density_computed_shells, 'Density Computed Shells Count, Cutoff 1.0e-12')

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

    # run rhf calculation 
    psi4.set_options({ 
        "scf_type": "direct", 
        "screening" : 'density', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12, 
        "reference" : "rhf",
        "save_jk": True,
        "bench" : 1 

    })
    rhf_energy, rhf_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # run uhf calculation 
    psi4.set_options({ 
        "scf_type": "direct", 
        "screening" : 'density', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12, 
        "reference" : "uhf",
        "save_jk": True,
        "bench" : 1

    })
    uhf_energy, uhf_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # prep for comparing results to expected values
    rhf_computed_shells = rhf_wfn.jk().computed_shells_per_iter("Quartets")
    uhf_computed_shells = uhf_wfn.jk().computed_shells_per_iter("Quartets")
    
    computed_shells_expected = [13171, 19618, 19665, 19657, 19661, 19661, 19663, 19663, 19663]

    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert(len(computed_shells_expected) == rhf_wfn.iteration_+1)
    assert(len(computed_shells_expected) == uhf_wfn.iteration_+1)

    # actually compare results with expected values
    assert compare(computed_shells_expected, rhf_computed_shells, 'Schwarz Computed Shells Count, Cutoff 1.0e-12')
    assert compare(computed_shells_expected, uhf_computed_shells, 'Density Computed Shells Count, Cutoff 1.0e-12')

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
