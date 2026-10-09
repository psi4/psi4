"""Cross-check the optional libcint (and simint) integral backends against the default Libint2.

libcint (INTEGRAL_PACKAGE=LIBCINT) is an experimental alternative two-electron
integral engine. Every quantity it produces should match Libint2 to numerical
precision: 4-center ERIs, range-separated erf ERIs, the density-fitting 2-/3-
center integrals, and the resulting SCF/MP2 energies -- for both spherical and
cartesian bases.
"""

import numpy as np
import pytest

from addons import using
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Alternate two-electron engines, each cross-checked against the default Libint2.
ENGINES = [
    pytest.param("libcint", marks=using("libcint")),
    pytest.param("simint", marks=using("simint")),
]


def _water():
    import psi4

    return psi4.geometry(
        """
        O  0.000000  0.000000  0.117790
        H  0.000000  0.755453 -0.471160
        H  0.000000 -0.755453 -0.471160
        units angstrom
        symmetry c1
        """
    )


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("basis,puream", [
    ("cc-pvdz", True),
    ("cc-pvtz", True),
    ("cc-pvdz", False),   # cartesian
    ("6-31G*", False),    # cartesian
])
def test_engine_ao_eri(engine, basis, puream):
    """4-center AO ERIs match Libint2 (spherical and cartesian)."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": basis, "puream": puream})

    def ao_eri(pkg):
        psi4.set_options({"integral_package": pkg})
        wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
        return np.asarray(psi4.core.MintsHelper(wfn.basisset()).ao_eri())

    ref = ao_eri("libint2")
    tst = ao_eri(engine)
    assert compare_values(ref, tst, 11, f"{engine} AO ERI {basis} puream={puream}")


@pytest.mark.parametrize("engine", [p for p in ENGINES if p.values[0] != "simint"])  # simint has no erf
def test_engine_erf_eri(engine):
    """Range-separated erf AO ERIs match Libint2."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": "cc-pvdz", "puream": True})

    def erf(pkg):
        psi4.set_options({"integral_package": pkg})
        wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
        return np.asarray(psi4.core.MintsHelper(wfn.basisset()).ao_erf_eri(0.3))

    assert compare_values(erf("libint2"), erf(engine), 11, f"{engine} AO erf-ERI omega=0.3")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("puream", [True, False])
def test_engine_df_integrals(engine, puream):
    """Density-fitting 2-center (Q|P) and 3-center (Q|mn) match Libint2."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": "cc-pvdz", "puream": puream})
    zero = psi4.core.BasisSet.zero_ao_basis_set()
    prim = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", puream=puream)
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "cc-pvdz-jkfit", puream=puream)

    def tensors(pkg):
        psi4.set_options({"integral_package": pkg})
        m = psi4.core.MintsHelper(prim)
        metric = np.asarray(m.ao_eri(aux, zero, aux, zero))    # (Q|P)
        three = np.asarray(m.ao_eri(aux, zero, prim, prim))    # (Q|mn)
        return metric, three

    m2, t2 = tensors("libint2")
    mc, tc = tensors(engine)
    assert compare_values(m2, mc, 10, f"{engine} DF metric (Q|P) puream={puream}")
    assert compare_values(t2, tc, 10, f"{engine} DF 3-center (Q|mn) puream={puream}")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("puream", [True, False])
def test_engine_scf_conventional(engine, puream):
    """Conventional (PK, 4-center) RHF energy matches Libint2."""
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({
        "basis": "cc-pvdz", "puream": puream, "scf_type": "pk",
        "guess": "core", "df_scf_guess": False,
        "e_convergence": 1e-10, "d_convergence": 1e-10,
    })

    def scf(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy("scf")

    assert compare_values(scf("libint2"), scf(engine), 9, f"{engine} PK-SCF puream={puream}")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("method", ["scf", "mp2"])
def test_engine_df_energies(engine, method):
    """Density-fitted RHF and MP2 energies match Libint2."""
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({
        "basis": "cc-pvdz", "puream": True, "scf_type": "df", "mp2_type": "df",
        "e_convergence": 1e-9, "d_convergence": 1e-9,
    })

    def energy(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy(method)

    assert compare_values(energy("libint2"), energy(engine), 8, f"{engine} DF-{method.upper()}")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("prim_puream,aux_puream", [(False, True), (True, False)])
def test_engine_df_integrals_mixed(engine, prim_puream, aux_puream):
    """DF integrals match Libint2 when primary and auxiliary disagree on cartesian/spherical."""
    import psi4

    psi4.core.clean()
    mol = _water()
    zero = psi4.core.BasisSet.zero_ao_basis_set()
    prim = psi4.core.BasisSet.build(mol, "ORBITAL", "6-31G*", puream=prim_puream)
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_MP2", "cc-pvdz-ri", puream=aux_puream)

    def tensors(pkg):
        psi4.set_options({"integral_package": pkg})
        m = psi4.core.MintsHelper(prim)
        return np.asarray(m.ao_eri(aux, zero, aux, zero)), np.asarray(m.ao_eri(aux, zero, prim, prim))

    m2, t2 = tensors("libint2")
    mc, tc = tensors(engine)
    assert compare_values(m2, mc, 10, f"{engine} DF metric (Q|P) prim/aux puream={prim_puream}/{aux_puream}")
    assert compare_values(t2, tc, 10, f"{engine} DF 3-center (Q|mn) prim/aux puream={prim_puream}/{aux_puream}")


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("basis", ["cc-pvdz", "6-31G*"])
def test_engine_threaded_df_mp2(engine, basis):
    """Threaded DF-MP2 (cloned integral objects; 6-31G* mixes a cartesian primary with a spherical RI aux)."""
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({
        "basis": basis, "scf_type": "df", "mp2_type": "df",
        "e_convergence": 1e-9, "d_convergence": 1e-9,
    })

    def energy(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy("mp2")

    nthreads = psi4.core.get_num_threads()
    psi4.set_num_threads(4)
    try:
        assert compare_values(energy("libint2"), energy(engine), 8, f"{engine} threaded DF-MP2 {basis}")
    finally:
        psi4.set_num_threads(nthreads)
