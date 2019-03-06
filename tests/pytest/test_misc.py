import pytest

import psi4

pytestmark = pytest.mark.quick


def hide_test_xtpl_fn_fn_error():
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy('cbs', scf_basis='cc-pvdz', scf_scheme=psi4.driver_cbs.xtpl_highest_1)

    assert 'Replace extrapolation function with function name' in str(e)


def hide_test_xtpl_cbs_fn_error():
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy(psi4.cbs, scf_basis='cc-pvdz')
        #psi4.energy(psi4.driver.driver_cbs.complete_basis_set, scf_basis='cc-pvdz')

    assert 'Replace cbs or complete_basis_set function with cbs string' in str(e)


@pytest.mark.parametrize("inp,out", [
    ((2, 'C2V'), 2),
    (('A2', 'c2v'), 2),
    (('2', 'C2V'), 2),
])
def test_parse_cotton_irreps(inp, out):
    idx = psi4.driver.driver_util.parse_cotton_irreps(*inp)
    assert idx == out


@pytest.mark.parametrize("inp", [
    ((5, 'cs')),
    (('5', 'cs')),
    ((0, 'cs')),
    (('a2', 'cs')),
])
def test_parse_cotton_irreps_error(inp):
    with pytest.raises(psi4.ValidationError) as e:
        psi4.driver.driver_util.parse_cotton_irreps(*inp)

    assert 'not valid for point group' in str(e)
