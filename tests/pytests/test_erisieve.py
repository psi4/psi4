import pytest
import psi4
import itertools
from utils import compare, compare_integers, compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_no_screening_schwarz():
    """Checks the number of shell quartets screened with Schwarz screening.
    no shell quartets should be screened with a threshold of 0.0"""

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
                       "screening": "density",
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

def test_no_screening_none():
    """Checks the number of shell quartets screened with None screening
    No shell quartets should be screened"""
    mol = psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)
    psi4.set_options({ "screening": "none" })

    basis = psi4.core.BasisSet.build(mol, target='cc-pVDZ')
    factory = psi4.core.IntegralFactory(basis)
    eri = factory.eri(0)

    shell_inds = range(basis.nshell())
    quartets = itertools.product(shell_inds, shell_inds, shell_inds, shell_inds)

    screen_count = 0
    for m, n, r, s in quartets:
        if not eri.shell_significant(m, n, r, s):
            screen_count += 1

    assert compare_integers(0, screen_count, 'Quartets None Screened')


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
    """Checks difference in Hartree-Fock energy between Schwarz screening and CSAM screening, which should be
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

def test_schwarz_vs_density_quartets_direct():
    """Checks difference between the number of shell quartets computed with Schwarz and Density screening for DirectJK. 
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

    schwarz_computed_shells_base_expected = [20290, 20290, 20290, 20290, 20290, 20290, 20290, 20290, 20290]
    density_computed_shells_base_expected = [13187, 19683, 19644, 19663, 19661, 19661, 19663, 19663, 19663]

    # we will give a 1% tolerance to account for potential future changes in screening procedures
    # e.g. https://github.com/psi4/psi4/pull/3138
    shell_count_tol = 0.01
    
    schwarz_computed_shells_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in schwarz_computed_shells_base_expected ] 
    density_computed_shells_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in density_computed_shells_base_expected ] 
      
    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert len(schwarz_computed_shells_expected) == schwarz_wfn.iteration_+1
    assert len(density_computed_shells_expected) == density_wfn.iteration_+1

    # actually compare results with expected values
    schwarz_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(schwarz_computed_shells, schwarz_computed_shells_expected) 
    ])
    assert compare(schwarz_pass, True, 'Schwarz Computed Shells Count, Cutoff 1.0e-12')
    
    density_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(density_computed_shells, density_computed_shells_expected) 
    ])
    assert compare(density_pass, True, 'Density Computed Shells Count, Cutoff 1.0e-12')

def test_density_screening_link():
    """Checks difference between the number of shell triplets+quartets computed with Density screening, and a reference value, for DFDirJ+LinK. 
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
        "scf_type": "dfdirj+link", 
        "screening" : 'density', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12, 
        "link_ints_tolerance" : 1e-12, 
        "save_jk": True,
        "bench" : 1 

    })
    density_energy, density_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # prep for comparing results to expected values
    density_computed_triplets = density_wfn.jk().computed_shells_per_iter("Triplets") # shell triplets, from DFDirJ
    density_computed_quartets = density_wfn.jk().computed_shells_per_iter("Quartets") # shell quartets, from LinK

    # reference values, acquired from DFDirJ+LinK from Psi4 v1.8
    density_computed_triplets_base_expected = [17680, 29433, 29488, 29480, 29482, 29482, 29482, 29482, 29482]
    density_computed_quartets_base_expected = [8019, 19341, 19366, 19371, 19371, 19371, 19371, 19371, 19371]    

    # we will give a 1% tolerance to account for potential future changes in screening procedures
    # e.g. https://github.com/psi4/psi4/pull/3138
    shell_count_tol = 0.01
    
    density_computed_triplets_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in density_computed_triplets_base_expected ] 
    density_computed_quartets_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in density_computed_quartets_base_expected ] 
      
    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert len(density_computed_triplets_expected) == density_wfn.iteration_+1
    assert len(density_computed_quartets_expected) == density_wfn.iteration_+1

    # actually compare results with expected values
    triplets_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(density_computed_triplets, density_computed_triplets_expected) 
    ])
    assert compare(triplets_pass, True, 'DFDirJ+LinK Computed Shell Triplets Count, Cutoff 1.0e-12')

    quartets_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(density_computed_quartets, density_computed_quartets_expected) 
    ])
    assert compare(quartets_pass, True, 'DFDirJ+LinK Computed Shell Quartets Count, Cutoff 1.0e-12')

def test_schwarz_screening_cosx():
    """Checks difference between the number of shell triplets+quartets computed with Schwarz screening, and a reference value, for DFDirJ+COSX. 
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
        "scf_type": "dfdirj+cosx", 
        "screening" : 'schwarz', 
        "df_scf_guess" : False,
        "integral_package": 'libint2', 
        "ints_tolerance" : 1e-12, 
        "cosx_ints_tolerance" : 1e-12, 
        "save_jk": True,
        "bench" : 1 

    })
    schwarz_energy, schwarz_wfn = psi4.energy('hf/DZ', return_wfn=True)

    # prep for comparing results to expected values
    schwarz_computed_triplets = schwarz_wfn.jk().computed_shells_per_iter("Triplets") # shell triplets, from DFDirJ
    schwarz_computed_pairs = schwarz_wfn.jk().computed_shells_per_iter("Pairs") # ESP shell pairs, from COSX

    # reference values, acquired from DFDirJ+COSX from Psi4 v1.8
    #schwarz_computed_triplets_expected = [17680, 29433, 29488, 29480, 29482, 29478, 29478, 29478, 29478, 29478]
    schwarz_computed_triplets_base_expected = [17671, 29407, 29469, 29480, 29482, 29478, 29478, 29478, 29478, 29478]
    schwarz_computed_pairs_base_expected = [835082, 864442, 868290, 867307, 867859, 867914, 867930, 867934, 867936, 2543375]

    # we will give a 1% tolerance to account for potential future changes in screening procedures
    # e.g. https://github.com/psi4/psi4/pull/3138
    shell_count_tol = 0.01
    
    schwarz_computed_triplets_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in schwarz_computed_triplets_base_expected ] 
    schwarz_computed_pairs_expected = [ (_ * (1.0 - shell_count_tol), _ * (1.0 + shell_count_tol)) for _ in schwarz_computed_pairs_base_expected ] 
      
    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert len(schwarz_computed_triplets_expected) == schwarz_wfn.iteration_+1
    assert len(schwarz_computed_pairs_expected) == schwarz_wfn.iteration_+1

    # actually compare results with expected values
    triplets_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(schwarz_computed_triplets, schwarz_computed_triplets_expected) 
    ])
    assert compare(triplets_pass, True, 'DFDirJ+COSX Computed Shell Triplets Count, Cutoff 1.0e-12')

    pairs_pass = all([ 
        expected[0] <= actual and actual <= expected[1] 
        for actual, expected in zip(schwarz_computed_pairs, schwarz_computed_pairs_expected) 
    ])
    assert compare(pairs_pass, True, 'DFDirJ+COSX Computed ESP Shell Pairs Count, Cutoff 1.0e-12')

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
    
    computed_shells_expected = [13187, 19683, 19644, 19663, 19661, 19661, 19663, 19663, 19663]

    # compare iteration counts of runs with computed shell quartet array lengths
    # iteration_+1 is used to account for computed_shells arrays including SAD guess results
    assert len(computed_shells_expected) == rhf_wfn.iteration_+1
    assert len(computed_shells_expected) == uhf_wfn.iteration_+1

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

def test_schwarz_vs_none_energy():
    """Checks difference in Hartree-Fock energy between Schwarz screening and no screening, which should be
    insignificant. """

    psi4.geometry("""
        Ne 0.0 0.0 0.0
        Ne 4.0 0.0 0.0
        Ne 8.0 0.0 0.0
    """)

    psi4.set_options({'d_convergence' : 1e-12,
                      'screening' : 'schwarz',
                      'ints_tolerance' : 1.0e-12 })
    e_schwarz = psi4.energy('hf/DZ')

    psi4.core.clean()

    psi4.set_options({'d_convergence' : 1e-12,
                      'screening' : 'none',
                      'ints_tolerance' : 1.0e-12 })
    e_none = psi4.energy('hf/DZ')

    assert compare_values(e_schwarz, e_none, 11, 'Schwarz vs None Screening, Cutoff 1.0e-12')


