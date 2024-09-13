import os
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runnonroutine", action="store_true", default=False, help="run the nonroutine tests in stdsuite"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runnonroutine"):
        # --runnonroutine given in cli: do not skip nonroutine tests
        return
    skip_nonroutine = pytest.mark.skip(reason="need --runnonroutine option to run")
    for item in items:
        if "nonroutine" in item.keywords:
            item.add_marker(skip_nonroutine)


@pytest.fixture(scope="session", autouse=True)
def set_up_overall(request, tmp_path_factory):
    import psi4

    psi4.set_output_file("pytest_output.dat", False)
    os.chdir(tmp_path_factory.getbasetemp())
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
        "pytest_output.log",
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


@pytest.fixture(scope="function")
def snowflake():
    try:
        from qcfractal import FractalSnowflake
        qca_next_branch = False
    except ImportError:
        try:
            from qcfractal.snowflake import FractalSnowflake
            qca_next_branch = True
        except ImportError:
            return None

    if qca_next_branch:
        snowflake = FractalSnowflake()
    else:
        snowflake = FractalSnowflake(logging=True, max_workers=4)
    return snowflake
