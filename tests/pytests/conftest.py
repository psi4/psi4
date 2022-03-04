import pytest


def pytest_configure(config):
    # Register marks to avoid warnings in psi4.test()
    # sync with pyproject.toml
    config.addinivalue_line("markers", "psi: test defined in Psi4 codebase")
    config.addinivalue_line("markers", "cli: test also defined in CTest, usually in Psithon")

    config.addinivalue_line("markers", "check_triplet")
    config.addinivalue_line("markers", "dft")
    config.addinivalue_line("markers", "gga")
    config.addinivalue_line("markers", "hf")
    config.addinivalue_line("markers", "hyb_gga")
    config.addinivalue_line("markers", "hyb_gga_lrc")
    config.addinivalue_line("markers", "lda")
    config.addinivalue_line("markers", "long")
    config.addinivalue_line("markers", "mdi")
    config.addinivalue_line("markers", "mp2")
    config.addinivalue_line("markers", "restricted_singlet")
    config.addinivalue_line("markers", "restricted_triplet")
    config.addinivalue_line("markers", "RPA")
    config.addinivalue_line("markers", "scf")
    config.addinivalue_line("markers", """slow: marks tests as slow (deselect with '-m "not slow"')""")
    config.addinivalue_line("markers", "smoke")
    config.addinivalue_line("markers", "solver")
    config.addinivalue_line("markers", "stress")
    config.addinivalue_line("markers", "TDA")
    config.addinivalue_line("markers", "tdscf")
    config.addinivalue_line("markers", "quick")
    config.addinivalue_line("markers", "unittest")
    config.addinivalue_line("markers", "unrestricted")
    config.addinivalue_line("markers", "cart: test geometries are purely numerical Cartesians and no Z-matrices")
    config.addinivalue_line("markers", "nbody: tests that use the nbody or manybody wrapper")
    config.addinivalue_line("markers", "gradient: tests that run 1st-derivatives wrt energy")
    config.addinivalue_line("markers", "hessian: tests that run 2nd-derivatives wrt energy")
    config.addinivalue_line("markers", "alpha")
    config.addinivalue_line("markers", "cc")
    config.addinivalue_line("markers", "pte")
    config.addinivalue_line("markers", "properties")
    config.addinivalue_line("markers", "dipole")
    config.addinivalue_line("markers", "ghosts")
    config.addinivalue_line("markers", "opt")
    config.addinivalue_line("markers", "df: tests that employ density-fitting")
    #config.addinivalue_line("markers", "")

    config.addinivalue_line("markers", "addon")

    # QCEngine
    config.addinivalue_line("markers", "mp2d: tests using MP2D software; skip if unavailable")
    config.addinivalue_line("markers", "dftd3")
    config.addinivalue_line("markers", "cfour")
    config.addinivalue_line("markers", "gcp")
    config.addinivalue_line("markers", "dftd4: tests using DFTD4 software; skip if unavailable")
    # inherited from QCEngine
    config.addinivalue_line("markers", "dftd3_321")
    config.addinivalue_line("markers", "psi4")

    # run if non-QC software available
    config.addinivalue_line("markers", "memory_profiler")
    config.addinivalue_line("markers", "networkx")

    # run if QC software available
    config.addinivalue_line("markers", "adcc: tests using ADCconnect software; skip if unavailable")
    config.addinivalue_line("markers", "ambit: tests using ambit software; skip if unavailable")
    config.addinivalue_line("markers", "cct3")
    config.addinivalue_line("markers", "chemps2: tests using CheMPS2 software; skip if unavailable")
    config.addinivalue_line("markers", "cppe: tests using cppe software; skip if unavailable")
    config.addinivalue_line("markers", "dkh: tests using dkh software; skip if unavailable")
    config.addinivalue_line("markers", "libefp: tests using LibEFP software; skip if unavailable")
    config.addinivalue_line("markers", "erd: tests using ERD software; skip if unavailable")
    config.addinivalue_line("markers", "fockci: tests using XX software; skip if unavailable")
    config.addinivalue_line("markers", "forte")
    config.addinivalue_line("markers", "gdma: tests using gdma software; skip if unavailable")
    config.addinivalue_line("markers", "gpu_dfcc")
    config.addinivalue_line("markers", "geometric: tests using geomeTRIC software; skip if unavailable")
    config.addinivalue_line("markers", "ipi: tests using i-PI software; skip if unavailable")
    config.addinivalue_line("markers", "mrcc")
    config.addinivalue_line("markers", "openfermionpsi4")
    config.addinivalue_line("markers", "pcmsolver: tests using PCMSolver software; skip if unavailable")
    config.addinivalue_line("markers", "psixas")
    config.addinivalue_line("markers", "resp: tests using resp software; skip if unavailable")
    config.addinivalue_line("markers", "simint: tests using XX software; skip if unavailable")
    config.addinivalue_line("markers", "snsmp2")
    config.addinivalue_line("markers", "v2rdm_casscf")
    config.addinivalue_line("markers", "qcdb")


@pytest.fixture(scope="session", autouse=True)
def set_up_overall(request):
    import psi4

    psi4.set_output_file("pytest_output.dat", False)
    request.addfinalizer(tear_down)


@pytest.fixture(scope="function", autouse=True)
def set_up():
    import psi4

    psi4.core.clean()
    psi4.core.clean_timers()
    psi4.core.clean_options()
    psi4.set_output_file("pytest_output.dat", True)


def tear_down():
    import os
    import glob
    import psi4

    psi4.core.close_outfile()
    patterns = [
        "cavity.*",
        "grid*",
        "pytest_output.*h5",
        "pytest_output.dat",
        "pytest_output.*grad",
        "*pcmsolver.inp",
        "PEDRA.OUT*",
        "timer.dat",
        "FCIDUMP_SCF",
        "FCIDUMP_MP2",
        "*.fchk",
        "*.molden",
    ]
    pytest_scratches = []
    for pat in patterns:
        pytest_scratches.extend(glob.glob(pat))
    for fl in pytest_scratches:
        try:
            os.unlink(fl)
        except (OSError, PermissionError):
            pass
