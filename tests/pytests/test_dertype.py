import pytest

import psi4
from psi4.driver.driver_util import negotiate_derivative_type


def ordinary():
    pass


def select_fail(mtd, probe):
    raise psi4.ManagedMethodError('abcdef')


def select_pass(mtd, probe):
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
    },
    'gradient': {
        'hf': ordinary,
        'cc': ordinary,
        'EGh_man': select_pass,
        'Egh_man2': select_fail,
        'Eg_man': select_fail,
    },
    'hessian': {
        'hf': ordinary,
        'EGh_man': select_fail,
        'Egh_man2': select_fail,
    },
}


@pytest.mark.parametrize("inp,out", [
    (('hessian', 'hf', None), "analytic"),
    (('hessian', 'hf', 2), "analytic"),
    (('hessian', 'hf', 1), "2_1"),
    (('hessian', 'hf', 0), "2_0"),
    (('gradient', 'hf', None), "analytic"),
    (('gradient', 'hf', 1), "analytic"),
    (('gradient', 'hf', 0), "1_0"),
    (('energy', 'hf', None), "analytic"),
    (('energy', 'hf', 0), "analytic"),
    (('hessian', 'cc', None), "2_1"),
    (('hessian', 'cc', 1), "2_1"),
    (('hessian', 'cc', 0), "2_0"),
    (('gradient', 'cc', None), "analytic"),
    (('gradient', 'cc', 1), "analytic"),
    (('gradient', 'cc', 0), "1_0"),
    (('energy', 'cc', None), "analytic"),
    (('energy', 'cc', 0), "analytic"),
    (('hessian', 'EGh_man', None), "2_1"),
    (('hessian', 'EGh_man', 1), "2_1"),
    (('hessian', 'EGh_man', 0), "2_0"),
    (('gradient', 'EGh_man', None), "analytic"),
    (('gradient', 'EGh_man', 1), "analytic"),
    (('gradient', 'EGh_man', 0), "1_0"),
    (('energy', 'EGh_man', None), "analytic"),
    (('energy', 'EGh_man', 0), "analytic"),
    (('hessian', 'Egh_man2', None), "2_0"),
    (('hessian', 'Egh_man2', 0), "2_0"),
    (('gradient', 'Egh_man2', None), "1_0"),
    (('gradient', 'Egh_man2', 0), "1_0"),
    (('energy', 'Egh_man2', None), "analytic"),
    (('energy', 'Egh_man2', 0), "analytic"),
    (('hessian', 'mp4', None), "2_0"),
    (('hessian', 'mp4', 0), "2_0"),
    (('gradient', 'mp4', None), "1_0"),
    (('gradient', 'mp4', 0), "1_0"),
    (('energy', 'mp4', None), "analytic"),
    (('energy', 'mp4', 0), "analytic"),
    (('hessian', 'Eg_man', None), "2_0"),
    (('hessian', 'Eg_man', 0), "2_0"),
    (('gradient', 'Eg_man', None), "1_0"),
    (('gradient', 'Eg_man', 0), "1_0"),
    (('energy', 'Eg_man', None), "analytic"),
    (('energy', 'Eg_man', 0), "analytic"),
])
def test_negotiate(inp, out):
    dertype = negotiate_derivative_type(proc=mock_proc, return_strategy=True, *inp)
    assert dertype == out


@pytest.mark.parametrize("inp", [
    (('hessian', 'cc', 2)),
    (('hessian', 'mp4', 2)),
    (('hessian', 'mp4', 1)),
    (('gradient', 'mp4', 1)),
    (('hessian', 'Eg_man', 2)),
    (('hessian', 'Eg_man', 1)),
    (('gradient', 'Eg_man', 1)),
    (('energy', 'unavail', None)),
    (('energy', 'unavail', 0)),
    (('energy', 'E_man', None)),
    (('energy', 'E_man', 0)),
    (('hessian', 'Egh_man2', 1)),
    (('gradient', 'Egh_man2', 1)),
])
def test_negotiate_missing_error(inp):
    with pytest.raises(psi4.MissingMethodError) as e:
        negotiate_derivative_type(proc=mock_proc, return_strategy=True, *inp)


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
])
def test_negotiate_excessive_error(inp):
    with pytest.raises(psi4.ValidationError) as e:
        negotiate_derivative_type(proc=mock_proc, return_strategy=True, *inp)

    assert 'excessive for target calculation' in str(e)
