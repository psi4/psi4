import pytest

@pytest.fixture(scope="session", autouse=True)
def set_up_overall(request):
    import psi4
    psi4.set_output_file("output.dat", False)
    request.addfinalizer(tear_down)

@pytest.fixture(scope="function", autouse=True)
def set_up():
    import psi4
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.set_output_file("output.dat", True)

def tear_down():
    pass
