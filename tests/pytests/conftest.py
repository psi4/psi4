import pytest


def pytest_configure(config):
    # Register marks to avoid warnings in psi4.test()
    # sync with setup.cfg
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
    ]
    pytest_scratches = []
    for pat in patterns:
        pytest_scratches.extend(glob.glob(pat))
    for fl in pytest_scratches:
        try:
            os.unlink(fl)
        except (OSError, PermissionError):
            pass
