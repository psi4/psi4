import pytest

import psi4
from psi4.driver.driver_util import negotiate_derivative_type

pytestmark = [pytest.mark.psi, pytest.mark.api]

def ordinary():
    pass


def select_fail(mtd, **kwargs):
    raise psi4.ManagedMethodError(['select_fail', "abcdef", 'MP2_TYPE', "CD", "RHF", "labocc"])


def select_pass(mtd, **kwargs):
    ordinary()


mock_proc = {
    'energy': {
        'hf': ordinary,
        'cc': ordinary,
        'mp4': ordinary,
        'EGh_man': select_pass,
        'Egh_man2': select_pass,
        'Eg_man': select_pass,
        'E_man': select_fail,
        'eg_man2': select_fail,
    },
    'gradient': {
        'hf': ordinary,
        'cc': ordinary,
        'EGh_man': select_pass,
        'Egh_man2': select_fail,
        'Eg_man': select_fail,
        'eg_man2': select_fail,
    },
    'hessian': {
        'hf': ordinary,
        'EGh_man': select_fail,
        'Egh_man2': select_fail,
    },
    "properties": {
        "hf": ordinary,
        "cc": ordinary,
    },
}


@pytest.mark.parametrize("inp,out", [
    pytest.param(('hessian', 'hf', None), (2, 2), marks=pytest.mark.d2ints),  # analytic
    pytest.param(('hessian', 'hf', 2), (2, 2), marks=pytest.mark.d2ints),  # analytic
    (('hessian', 'hf', 1), (2, 1)),
    (('hessian', 'hf', 0), (2, 0)),
    (('gradient', 'hf', None), (1, 1)),  # analytic
    (('gradient', 'hf', 1), (1, 1)),  # analytic
    (('gradient', 'hf', 0), (1, 0)),
    (('energy', 'hf', None), (0, 0)),  # analytic
    (('energy', 'hf', 0), (0, 0)),  # analytic
    (('hessian', 'cc', None), (2, 1)),
    (('hessian', 'cc', 1), (2, 1)),
    (('hessian', 'cc', 0), (2, 0)),
    (('gradient', 'cc', None), (1, 1)),  # analytic
    (('gradient', 'cc', 1), (1, 1)),  # analytic
    (('gradient', 'cc', 0), (1, 0)),
    (('energy', 'cc', None), (0, 0)),  # analytic
    (('energy', 'cc', 0), (0, 0)),  # analytic
    (('hessian', 'EGh_man', None), (2, 1)),
    (('hessian', 'EGh_man', 1), (2, 1)),
    (('hessian', 'EGh_man', 0), (2, 0)),
    (('gradient', 'EGh_man', None), (1, 1)),  # analytic
    (('gradient', 'EGh_man', 1), (1, 1)),  # analytic
    (('gradient', 'EGh_man', 0), (1, 0)),
    (('energy', 'EGh_man', None), (0, 0)),  # analytic
    (('energy', 'EGh_man', 0), (0, 0)),  # analytic
    (('hessian', 'Egh_man2', None), (2, 0)),
    (('hessian', 'Egh_man2', 0), (2, 0)),
    (('gradient', 'Egh_man2', None), (1, 0)),
    (('gradient', 'Egh_man2', 0), (1, 0)),
    (('energy', 'Egh_man2', None), (0, 0)),  # analytic
    (('energy', 'Egh_man2', 0), (0, 0)),  # analytic
    (('hessian', 'mp4', None), (2, 0)),
    (('hessian', 'mp4', 0), (2, 0)),
    (('gradient', 'mp4', None), (1, 0)),
    (('gradient', 'mp4', 0), (1, 0)),
    (('energy', 'mp4', None), (0, 0)),  # analytic
    (('energy', 'mp4', 0), (0, 0)),  # analytic
    (('hessian', 'Eg_man', None), (2, 0)),
    (('hessian', 'Eg_man', 0), (2, 0)),
    (('gradient', 'Eg_man', None), (1, 0)),
    (('gradient', 'Eg_man', 0), (1, 0)),
    (('energy', 'Eg_man', None), (0, 0)),  # analytic
    (('energy', 'Eg_man', 0), (0, 0)),  # analytic
    (("properties", "hf", None), ("prop", "prop")),
    (("properties", "cc", 0), ("prop", "prop")),
])
def test_negotiate(inp, out):
    dertype = negotiate_derivative_type(proc=mock_proc, *inp)
    assert dertype == out


@pytest.mark.parametrize("inp", [
    (('hessian', 'cc', 2)),
    (('hessian', 'mp4', 2)),
    (('hessian', 'mp4', 1)),
    (('gradient', 'mp4', 1)),
    (('hessian', 'Eg_man', 2)),
    (('hessian', 'Eg_man', 1)),
    (('gradient', 'Eg_man', 1)),
    (('gradient', 'eg_man2', None)),
    (('energy', 'unavail', None)),
    (('energy', 'unavail', 0)),
    (('energy', 'E_man', None)),
    (('energy', 'E_man', 0)),
    (('hessian', 'Egh_man2', 1)),
    (('gradient', 'Egh_man2', 1)),
    (("properties", "unavail", None)),
    (("properties", "unavail", 0)),
    (("properties", "ccop", 0)),
])
def test_negotiate_missing_error(inp):
    with pytest.raises(psi4.MissingMethodError) as e:
        negotiate_derivative_type(proc=mock_proc, *inp)


@pytest.mark.parametrize("inp", [
    (('gradient', 'hf', 2)),
    (('energy', 'hf', 2)),
    (('energy', 'hf', 1)),
    (('gradient', 'cc', 2)),
    (('energy', 'cc', 2)),
    (('energy', 'cc', 1)),
    (('gradient', 'EGh_man', 2)),
    (('energy', 'EGh_man', 2)),
    (('energy', 'EGh_man', 1)),
    (('gradient', 'Egh_man2', 2)),
    (('energy', 'Egh_man2', 2)),
    (('energy', 'Egh_man2', 1)),
    (('gradient', 'mp4', 2)),
    (('energy', 'mp4', 2)),
    (('energy', 'mp4', 1)),
    (('gradient', 'Eg_man', 2)),
    (('energy', 'Eg_man', 2)),
    (('energy', 'Eg_man', 1)),
    #(("properties", "hf", 2)),  # todo: does not raise b/c prop not in hierarchical sequence. investigate whether should force error.
])
def test_negotiate_excessive_error(inp):
    with pytest.raises(psi4.ValidationError) as e:
        negotiate_derivative_type(proc=mock_proc, *inp)

    assert 'excessive for target calculation' in str(e.value)
