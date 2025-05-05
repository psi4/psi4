import pytest
from addons import uusing
from utils import *

import math

import numpy as np
import qcelemental as qcel

import psi4
from psi4.driver import qcdb

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_variables_deprecated():
    psi4.geometry("he")
    psi4.set_options({"basis": "cc-pvtz", "qc_module": "occ"})
    psi4.energy("mp2")

    qcvars = psi4.core.variables()
    assert "SCS(N)-MP2 CORRELATION ENERGY" in qcvars.keys()
    assert "SCSN-MP2 CORRELATION ENERGY" not in qcvars.keys()

    qcvars = psi4.core.variables(include_deprecated_keys=True)
    assert "SCS(N)-MP2 CORRELATION ENERGY" in qcvars.keys()
    assert "SCSN-MP2 CORRELATION ENERGY" in qcvars.keys()


@pytest.mark.parametrize("call",
    [psi4.energy, psi4.optimize, psi4.gradient, psi4.hessian, psi4.frequencies, psi4.properties])
def test_typo_method_calls(call):
    psi4.geometry('He')
    err_substr = "Did you mean?"

    with pytest.raises(Exception) as e:
        call('ccsdd', basis='cc-pvdz')
    assert err_substr in str(e.value)


def test_xtpl_fn_fn_error():
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy('cbs', scf_basis='cc-pvdz', scf_scheme=psi4.driver_cbs.xtpl_highest_1)

    assert 'Replace extrapolation function with function name' in str(e.value)


@pytest.mark.parametrize("call",
    [psi4.energy, psi4.optimize, psi4.gradient, psi4.hessian, psi4.frequencies, psi4.properties])
def test_xtpl_cbs_fn_error(call):
    psi4.geometry('He')

    with pytest.raises(psi4.UpgradeHelper) as e:
        call(psi4.cbs, scf_basis='cc-pvdz')
        #psi4.energy(psi4.driver.driver_cbs.complete_basis_set, scf_basis='cc-pvdz')

    assert 'Replace cbs or complete_basis_set function with cbs string' in str(e.value)


def test_xtpl_gold_fn_error():
    psi4.geometry('He')
    from psi4.driver.aliases import sherrill_gold_standard

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.energy(sherrill_gold_standard, scf_basis='cc-pvdz')

    assert 'Replace function `energy(sherrill_gold_standard)' in str(e.value)


def test_qmmm_class_error():
    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.QMMM()

    assert 'Replace object with a list of charges and locations in Bohr passed as keyword argument' in str(e.value)


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


def test_slash_in_molecule_name_plus_dfhelper():
    mymol = psi4.core.Molecule.from_arrays(geom=[0, 0, 0, 2, 0, 0], elem=["h", "h"], name="h2/mol")

    # segfaults if any DF (that is, following line commented). runs if DF suppressed (following line active)
    #psi4.set_options({"scf_type": "pk", "df_basis_guess": "false"})

    ene = psi4.energy("b3lyp/cc-pvtz", molecule=mymol)
    assert compare_values(-1.00125358, ene, 5, 'weird mol name ok as file')


# <<<  TODO Deprecated! Delete in Psi4 v1.5  >>>

@uusing("networkx")
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


def test_deprecated_component_dipole():

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.variable("current dipole x")

def test_deprecated_set_module_options():
    err_substr = "instead of `psi4.set_options({<module>__<keys>: <vals>})`"

    with pytest.warns(FutureWarning) as e:
        psi4.set_module_options('scf', {'df_basis_scf': 9})


def test_renamed_qcvars():
    psi4.set_variable("SCS(N)-MP2 TOTAL ENERGY", 3.3)

    with pytest.warns(FutureWarning) as e:
        ans = psi4.variable("SCSN-MP2 TOTAL ENERGY")

    assert ans == 3.3


def test_cancelled_qcvars():
    err_substr = "no direct replacement"

    with pytest.raises(psi4.UpgradeHelper) as e:
        psi4.variable("scsn-mp2 same-spin correlation energy")

    assert err_substr in str(e.value)
