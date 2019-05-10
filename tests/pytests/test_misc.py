import pytest
from .addons import using_networkx
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


# <<<  TODO Deprecated! Delete in Psi4 v1.5  >>>

@using_networkx
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
