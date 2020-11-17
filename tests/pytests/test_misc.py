import pytest
from .addons import using
from .utils import *

import math

import numpy as np
import qcelemental as qcel

import psi4
from psi4.driver import qcdb

pytestmark = pytest.mark.quick


def hide_test_xtpl_fn_fn_error():
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy('cbs', scf_basis='cc-pvdz', scf_scheme=psi4.driver_cbs.xtpl_highest_1)

    assert 'Replace extrapolation function with function name' in str(e.value)


def hide_test_xtpl_cbs_fn_error():
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy(psi4.cbs, scf_basis='cc-pvdz')
        #psi4.energy(psi4.driver.driver_cbs.complete_basis_set, scf_basis='cc-pvdz')

    assert 'Replace cbs or complete_basis_set function with cbs string' in str(e.value)


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

    assert 'not valid for point group' in str(e.value)


# <<<  TODO Deprecated! Delete in Psi4 v1.5  >>>

@using("networkx")
def test_deprecated_qcdb_align_b787():

    soco10 = """
    O  1.0 0.0 0.0
    C  0.0 0.0 0.0
    O -1.0 0.0 0.0
    units ang
    """

    sooc12 = """
    O  1.2 4.0 0.0
    O -1.2 4.0 0.0
    C  0.0 4.0 0.0
    units ang
    """

    ref_rmsd = math.sqrt(2. * 0.2 * 0.2 / 3.)  # RMSD always in Angstroms

    oco10 = qcel.molparse.from_string(soco10)
    oco12 = qcel.molparse.from_string(sooc12)

    oco10_geom_au = oco10['qm']['geom'].reshape((-1, 3)) / qcel.constants.bohr2angstroms
    oco12_geom_au = oco12['qm']['geom'].reshape((-1, 3)) / qcel.constants.bohr2angstroms

    with pytest.warns(FutureWarning) as err:
        rmsd, mill = qcdb.align.B787(
            oco10_geom_au, oco12_geom_au, np.array(['O', 'C', 'O']), np.array(['O', 'O', 'C']), verbose=4, do_plot=False)

    assert compare_values(ref_rmsd, rmsd, 6, 'known rmsd B787')



def test_deprecated_qcdb_align_scramble():
    with pytest.warns(FutureWarning) as err:
        mill = qcdb.align.compute_scramble(4, do_resort=False, do_shift=False, do_rotate=False, deflection=1.0, do_mirror=False)

    assert compare_arrays([0,1,2,3], mill.atommap, 4, 'atommap')

# <<<  TODO Deprecated! Delete when the error messages are removed.  >>>

def test_deprecated_dcft_calls():
    psi4.geometry('He')
    err_substr = "All instances of 'dcft' should be replaced with 'dct'."

    driver_calls = [psi4.energy, psi4.optimize, psi4.gradient, psi4.hessian, psi4.frequencies]

    for call in driver_calls:
        with pytest.raises(psi4.UpgradeHelper) as e:
            call('dcft', basis='cc-pvdz')
        assert err_substr in str(e.value)

    # The errors trapped below are C-side, so they're nameless, Py-side.
    with pytest.raises(Exception) as e:
        psi4.set_module_options('dcft', {'e_convergence': 9})

    assert err_substr in str(e.value)

    with pytest.raises(Exception) as e:
        psi4.set_module_options('dct', {'dcft_functional': 'odc-06'})

    assert err_substr in str(e.value)


def test_deprecated_component_dipole():

    #with pytest.warns(FutureWarning) as e:
    psi4.set_variable("current dipole x", 5)

    with pytest.warns(FutureWarning) as e:
        ans = psi4.variable("current dipole x")

    assert ans == 5
