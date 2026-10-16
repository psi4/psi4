"""Cross-check the optional libcint integral backend against the default Libint2.

libcint (INTEGRAL_PACKAGE=LIBCINT) is an experimental alternative two-electron
integral engine. Every quantity it produces should match Libint2 to numerical
precision: 4-center ERIs, range-separated erf ERIs, the density-fitting 2-/3-
center integrals, and the resulting SCF/MP2 energies -- for both spherical and
cartesian bases.
"""

import numpy as np
import pytest

from addons import has_program, uusing
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]


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


@uusing("libcint")
@pytest.mark.parametrize("basis,puream", [
    ("cc-pvdz", True),
    ("cc-pvtz", True),
    ("cc-pvdz", False),   # cartesian
    ("6-31G*", False),    # cartesian
])
def test_libcint_ao_eri(basis, puream):
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
    tst = ao_eri("libcint")
    assert compare_values(ref, tst, 11, f"AO ERI {basis} puream={puream}")


@uusing("libcint")
def test_libcint_erf_eri():
    """Range-separated erf AO ERIs match Libint2."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": "cc-pvdz", "puream": True})

    def erf(pkg):
        psi4.set_options({"integral_package": pkg})
        wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
        return np.asarray(psi4.core.MintsHelper(wfn.basisset()).ao_erf_eri(0.3))

    assert compare_values(erf("libint2"), erf("libcint"), 11, "AO erf-ERI omega=0.3")


@uusing("libcint")
@pytest.mark.parametrize("puream", [True, False])
def test_libcint_df_integrals(puream):
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
    mc, tc = tensors("libcint")
    assert compare_values(m2, mc, 10, f"DF metric (Q|P) puream={puream}")
    assert compare_values(t2, tc, 10, f"DF 3-center (Q|mn) puream={puream}")


@uusing("libcint")
@pytest.mark.parametrize("prim_puream,aux_puream", [(False, True), (True, False)])
def test_libcint_df_integrals_mixed(prim_puream, aux_puream):
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
    mc, tc = tensors("libcint")
    assert compare_values(m2, mc, 10, f"DF metric (Q|P) prim/aux puream={prim_puream}/{aux_puream}")
    assert compare_values(t2, tc, 10, f"DF 3-center (Q|mn) prim/aux puream={prim_puream}/{aux_puream}")


@uusing("libcint")
@pytest.mark.parametrize("basis", ["cc-pvdz", "6-31G*"])
def test_libcint_threaded_df_mp2(basis):
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
        assert compare_values(energy("libint2"), energy("libcint"), 8, f"threaded DF-MP2 {basis}")
    finally:
        psi4.set_num_threads(nthreads)


@uusing("libcint")
def test_libcint_screening_none():
    """DF codes initialize the ERI sieve manually under SCREENING NONE, which the engine must support."""
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

    assert compare_values(energy("libint2"), energy("libcint"), 8, "DF-MP2 screening=none")


@uusing("libcint")
def test_libcint_python_factory_eri():
    """A factory-built ERI object is usable from Python (as in test_erisieve)."""
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": "cc-pvdz", "integral_package": "libcint"})
    bas = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz")
    eri = psi4.core.IntegralFactory(bas).eri(0)
    assert eri.shell_significant(0, 0, 0, 0)
    eri.compute_shell(0, 0, 0, 0)


@uusing("libcint")
def test_libcint_bases_on_different_molecules():
    """Each basis keeps its own molecule's centers: (Q|mn) with Q on a different molecule than m and n."""
    import psi4

    psi4.core.clean()
    mol1 = _water()
    mol2 = psi4.geometry("""
        Ne 0.3 0.2 3.0
        units angstrom
        symmetry c1
        no_com
        no_reorient
        """)
    mol1.update_geometry()
    mol2.update_geometry()
    bs1 = psi4.core.BasisSet.build(mol1, "ORBITAL", "cc-pvdz")
    bs2 = psi4.core.BasisSet.build(mol2, "ORBITAL", "cc-pvdz")
    zero = psi4.core.BasisSet.zero_ao_basis_set()

    def eri(pkg):
        psi4.set_options({"integral_package": pkg})
        return np.asarray(psi4.core.MintsHelper(bs1).ao_eri(bs2, zero, bs1, bs1))

    ref = eri("libint2")
    assert np.abs(ref).max() > 1.e-2
    assert compare_values(ref, eri("libcint"), 10, "(Q|mn) with Q on another molecule")


@uusing("libcint")
def test_libcint_sapgau_guess():
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

    assert compare_values(energy("libint2"), energy("libcint"), 10, "SCF energy at first iteration from SAPGAU guess")


@uusing("libcint")
def test_libcint_diffuse_external_charge():
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

    assert compare_values(energy("libint2"), energy("libcint"), 9, "SCF with diffuse external charges")


@uusing("libcint")
@pytest.mark.parametrize("puream", [True, False])
def test_libcint_scf_conventional(puream):
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

    assert compare_values(scf("libint2"), scf("libcint"), 9, f"PK-SCF puream={puream}")


@uusing("libcint")
@pytest.mark.parametrize("method", ["scf", "mp2"])
def test_libcint_df_energies(method):
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

    assert compare_values(energy("libint2"), energy("libcint"), 8, f"DF-{method.upper()}")


def _diffuse_basis(shell, exp_o, exp_h, puream):
    return f"""
    assign diffuse
    [diffuse]
    {"spherical" if puream else "cartesian"}
    ****
    O 0
    S 1 1.00
      10.0 1.0
    P 1 1.00
      1.0 1.0
    {shell} 1 1.00
      {exp_o} 1.0
    ****
    H 0
    S 1 1.00
      1.0 1.0
    {shell} 1 1.00
      {exp_h} 1.0
    ****
"""


@uusing("libcint")
@pytest.mark.parametrize("puream", [True, False])
@pytest.mark.parametrize("shell,exp_o,exp_h,max_deriv", [
    ("G", "1.0e-4", "3.0e-4", 1),  # Libint2 here has no second derivatives of g functions
    ("F", "1.0e-5", "3.0e-5", 2),
])
def test_libcint_diffuse_high_l(shell, exp_o, exp_h, max_deriv, puream):
    """Diffuse high-l functions match Libint2 in the ERIs and their derivatives.

    libcint's screening of primitive pairs leaves out a factor that is large for
    diffuse primitives of high angular momentum. Without the raised cutoff, the
    one-center integrals of these shells come out exactly zero.
    """
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.basis_helper(_diffuse_basis(shell, exp_o, exp_h, puream), name="diffuse")
    psi4.set_options({"puream": puream})

    # INTEGRAL_PACKAGE is read when the integrals are computed, so compute each
    # one right after building its MintsHelper.
    def both(op):
        return [[np.asarray(m) for m in op(_mints(pkg, mol))] for pkg in ("libint2", "libcint")]

    label = f"diffuse {shell} puream={puream}"
    ref, tst = both(lambda m: [m.ao_eri()])
    assert compare_values(ref[0], tst[0], 10, f"{label} ERI")
    for atom in range(mol.natom()):
        ref, tst = both(lambda m: m.ao_tei_deriv1(atom))
        for xyz, (r, t) in enumerate(zip(ref, tst)):
            assert compare_values(r, t, 10, f"{label} d/d{'xyz'[xyz]}{atom}")
    if max_deriv >= 2:
        ref, tst = both(lambda m: m.ao_tei_deriv2(0, 0))
        for k, (r, t) in enumerate(zip(ref, tst)):
            assert compare_values(r, t, 10, f"{label} d2/d0d0 [{k}]")


def _deriv_setup(basis, puream):
    import psi4

    psi4.core.clean()
    mol = _water()
    psi4.set_options({"basis": basis, "puream": puream})
    return mol


def _mints(pkg, mol):
    import psi4

    psi4.set_options({"integral_package": pkg})
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    return psi4.core.MintsHelper(wfn.basisset())


@uusing("libcint")
@pytest.mark.parametrize("basis,puream,omega", [
    ("cc-pvdz", True, 0.0),
    ("cc-pvdz", False, 0.0),   # cartesian
    ("cc-pvtz", True, 0.0),    # through f
    ("cc-pvdz", True, 0.3),    # range-separated erf
])
def test_libcint_eri_deriv1(basis, puream, omega):
    """First nuclear derivatives of the 4-center ERIs match Libint2, atom by atom."""
    mol = _deriv_setup(basis, puream)
    for atom in range(mol.natom()):
        ref = [np.asarray(m) for m in _mints("libint2", mol).ao_tei_deriv1(atom, omega)]
        tst = [np.asarray(m) for m in _mints("libcint", mol).ao_tei_deriv1(atom, omega)]
        for xyz, (r, t) in enumerate(zip(ref, tst)):
            assert compare_values(r, t, 10, f"d(ab|cd)/d{'xyz'[xyz]}{atom} {basis} puream={puream} omega={omega}")


@uusing("libcint")
@pytest.mark.parametrize("puream", [True, False])
def test_libcint_eri_deriv2(puream):
    """Second nuclear derivatives of the 4-center ERIs match Libint2 for every atom pair."""
    mol = _deriv_setup("6-31G*", puream)
    for a1 in range(mol.natom()):
        for a2 in range(mol.natom()):
            ref = [np.asarray(m) for m in _mints("libint2", mol).ao_tei_deriv2(a1, a2)]
            tst = [np.asarray(m) for m in _mints("libcint", mol).ao_tei_deriv2(a1, a2)]
            for k, (r, t) in enumerate(zip(ref, tst)):
                assert compare_values(r, t, 10, f"d2(ab|cd)/d{a1}d{a2} [{k}] puream={puream}")


@uusing("libcint")
@pytest.mark.parametrize("puream", [True, False])
def test_libcint_df_deriv1(puream):
    """First derivatives of the DF metric (P|Q) and 3-center (P|mn) match Libint2."""
    import psi4

    mol = _deriv_setup("cc-pvdz", puream)
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "cc-pvdz-jkfit", puream=puream)
    for atom in range(mol.natom()):
        for kind in ("metric", "3center"):
            def get(pkg):
                m = _mints(pkg, mol)
                m.set_basisset("DF_BASIS_SCF", aux)
                f = m.ao_metric_deriv1 if kind == "metric" else m.ao_3center_deriv1
                return [np.asarray(x) for x in f(atom, "DF_BASIS_SCF")]

            for xyz, (r, t) in enumerate(zip(get("libint2"), get("libcint"))):
                assert compare_values(r, t, 10, f"DF {kind} d/d{'xyz'[xyz]}{atom} puream={puream}")


def _property(pkg, driver, method, options):
    import psi4

    psi4.core.clean()
    psi4.core.clean_options()
    _water()
    psi4.set_options({"basis": "cc-pvdz", "e_convergence": 1e-10, "d_convergence": 1e-10})
    psi4.set_options(options)
    psi4.set_options({"integral_package": pkg})
    return np.asarray(driver(method))


@uusing("libcint")
@pytest.mark.parametrize("method,options", [
    pytest.param("scf", {"scf_type": "pk"}, id="pk-rhf"),
    pytest.param("scf", {"scf_type": "pk", "puream": False}, id="pk-rhf-cartesian"),
    pytest.param("scf", {"scf_type": "df"}, id="df-rhf"),
    pytest.param("scf", {"scf_type": "pk", "reference": "uhf"}, id="pk-uhf"),
    pytest.param("mp2", {"scf_type": "df", "mp2_type": "df"}, id="df-mp2"),
    # a cartesian primary with spherical fitting bases
    pytest.param("scf", {"basis": "6-31G*", "scf_type": "df"}, id="df-rhf-631gs"),
    pytest.param("mp2", {"basis": "6-31G*", "scf_type": "df", "mp2_type": "df"}, id="df-mp2-631gs"),
    pytest.param("scf", {"scf_type": "df", "screening": "none"}, id="df-rhf-screening-none"),
    pytest.param("wb97x", {"scf_type": "pk", "dft_radial_points": 75, "dft_spherical_points": 302}, id="pk-wb97x"),
    pytest.param("wb97x", {"scf_type": "df", "dft_radial_points": 75, "dft_spherical_points": 302}, id="df-wb97x"),
])
def test_libcint_gradient(method, options):
    """Analytic gradients computed with libcint derivative integrals match Libint2."""
    import psi4

    # gradients carry the SCF convergence error linearly, so converge well past the 1e-9 check
    options = dict(options, e_convergence=1e-12, d_convergence=1e-12)
    ref = _property("libint2", psi4.gradient, method, options)
    tst = _property("libcint", psi4.gradient, method, options)
    assert compare_values(ref, tst, 9, f"{method} gradient {options}")


@uusing("libcint")
@pytest.mark.parametrize("options", [
    pytest.param({"scf_type": "pk"}, id="pk-rhf"),
    pytest.param({"scf_type": "pk", "puream": False}, id="pk-rhf-cartesian"),
    pytest.param({"scf_type": "df"}, id="df-rhf"),
    pytest.param({"basis": "6-31G*", "scf_type": "df"}, id="df-rhf-631gs"),
])
def test_libcint_hessian(options):
    """Analytic RHF Hessians computed with libcint second-derivative integrals match Libint2."""
    import psi4

    ref = _property("libint2", psi4.hessian, "scf", options)
    tst = _property("libcint", psi4.hessian, "scf", options)
    assert compare_values(ref, tst, 8, f"RHF Hessian {options}")


@pytest.mark.skipif(has_program("libcint"), reason="libcint is compiled in")
@pytest.mark.parametrize("method", ["hf", "wb97x"])
def test_libcint_not_compiled_in(method):
    import psi4

    _water()
    psi4.set_options({"basis": "cc-pvdz", "scf_type": "pk", "integral_package": "libcint"})
    with pytest.raises(RuntimeError, match="built without libcint"):
        psi4.energy(method)
