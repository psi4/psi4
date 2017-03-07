import pytest

@pytest.fixture(scope="session", autouse=True)
def set_up_overall(request):
    import psi4
    psi4.set_output_file("pytest_output.dat", False)
    request.addfinalizer(tear_down)

@pytest.fixture(scope="function", autouse=True)
def set_up():
    import psi4
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.set_output_file("pytest_output.dat", True)

def tear_down():
    import os
    import glob
    patterns = ['cavity.*', 'grid*', 'pytest_output.*h5',
                'pytest_output.dat',
                '*pcmsolver.inp', 'PEDRA.OUT*', 'timer.dat']
    pytest_scratches = []
    for pat in patterns:
        pytest_scratches.extend(glob.glob(pat))
    for fl in pytest_scratches:
        os.unlink(fl)
