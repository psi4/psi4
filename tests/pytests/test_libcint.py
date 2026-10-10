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


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("puream", [True, False])
def test_engine_diffuse_high_am(engine, puream):
    """Very diffuse high-l shells (as in auto-generated Cholesky aux bases) aren't screened away.

    libcint's primitive screening underestimates diffuse high-l quartets, e.g., dropping
    (gg|gg) = 0.01 for a g exponent of 2.6e-4 at its default cutoff (scf-auto-cholesky).
    """
    import psi4

    psi4.core.clean()
    mol = psi4.geometry("""
        He 0 0 0
        He 0 0 1.5
        symmetry c1
        """)
    psi4.basis_helper("""
        ****
        He 0
        S 1 1.0
          1.0 1.0
        P 1 1.0
          1.0e-3 1.0
        G 1 1.0
          0.5 1.0
        G 1 1.0
          2.6e-4 1.0
        ****
        """, name="diffuse_g", set_option=True)
    psi4.set_options({"puream": puream})

    def ao_eri(pkg):
        psi4.set_options({"integral_package": pkg})
        wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
        return np.asarray(psi4.core.MintsHelper(wfn.basisset()).ao_eri())

    ref = ao_eri("libint2")
    assert np.abs(ref).max() > 1.e-3
    assert compare_values(ref, ao_eri(engine), 10, f"{engine} AO ERI diffuse g puream={puream}")


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


@pytest.mark.parametrize("engine", [p for p in ENGINES if p.values[0] != "simint"])  # simint segfaults here
def test_engine_screening_none(engine):
    """DF codes initialize the ERI sieve manually under SCREENING NONE, which each engine must support."""
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({
        "basis": "cc-pvdz", "scf_type": "df", "mp2_type": "df", "screening": "none",
        "e_convergence": 1e-9, "d_convergence": 1e-9,
    })

    def energy(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy("mp2")

    assert compare_values(energy("libint2"), energy(engine), 8, f"{engine} DF-MP2 screening=none")


@pytest.mark.parametrize("engine", [pytest.param("libcint", marks=using("libcint"))])
def test_engine_gradient_fallback(engine, tmp_path):
    """LIBCINT computes energies itself and falls back to Libint2 for derivatives; LIBCINT_ONLY refuses."""
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "pk", "e_convergence": 1e-10, "d_convergence": 1e-10})

    def gradient(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return np.asarray(psi4.gradient("scf"))

    ref = gradient("libint2")
    out = tmp_path / "fallback.out"
    psi4.set_output_file(str(out), False)
    try:
        tst = gradient("libcint")
    finally:
        psi4.set_output_file("output.dat", True)
    assert compare_values(ref, tst, 8, "LIBCINT SCF gradient (Libint2 fallback for derivatives)")
    text = out.read_text()
    assert "Two-electron integrals (ERI) from libcint." in text
    assert "Two-electron integrals (ERI 1st deriv) from Libint2 (fallback: libcint doesn't compute them)." in text

    psi4.set_options({"integral_package": "libcint_only"})
    psi4.core.clean()
    msg = "LIBCINT_ONLY can't supply ERI 1st deriv integrals: libcint doesn't compute them"
    with pytest.raises(Exception, match=msg):
        psi4.gradient("scf")


@pytest.mark.parametrize("engine", ENGINES)
def test_engine_sapgau_guess(engine):
    """The SAPGAU guess's 3-center integrals use deliberately non-normalized s shells.

    Stop at the first iteration, whose energy reflects the guess, since a converged
    energy hides a wrong guess.
    """
    import psi4

    psi4.core.clean()
    _water()
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "pk", "guess": "sapgau", "maxiter": 1,
                      "fail_on_maxiter": False})

    def energy(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy("scf")

    ref = energy("libint2")
    assert compare_values(ref, energy(engine), 10, f"{engine} SCF energy at first iteration from SAPGAU guess")


@pytest.mark.parametrize("engine", ENGINES)
def test_engine_diffuse_external_charge(engine):
    """Diffuse external charges use deliberately non-normalized s shells (as do SAPGAU guesses)."""
    import psi4

    psi4.core.clean()
    mol = _water()
    mol.fix_com(True)
    mol.fix_orientation(True)
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "df", "e_convergence": 1e-10, "d_convergence": 1e-10})
    external = [None, [[0.5, [1.5, 0.5, 0.0], 0.8], [-0.5, [-1.0, 1.0, 1.0], 2.0]]]

    def energy(pkg):
        psi4.set_options({"integral_package": pkg})
        psi4.core.clean()
        return psi4.energy("scf", external_potentials=external)

    assert compare_values(energy("libint2"), energy(engine), 9, f"{engine} SCF with diffuse external charges")
