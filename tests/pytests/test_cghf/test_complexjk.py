import pytest
import numpy as np

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]


def _random_complex(shape, rng):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _build_complex_jk():
    """Minimal ComplexDirectJK via ComplexJK.build_JK (SCF_TYPE=DIRECT)."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    return psi4.core.ComplexJK.build_JK(basis, basis)


def _run_compute_D(C_left, C_right=None):
    """Queue one or more C pairs and form D. C_* may be a matrix or a sequence of matrices."""
    if C_right is None:
        C_right = C_left

    lefts = C_left if isinstance(C_left, (list, tuple)) else [C_left]
    rights = C_right if isinstance(C_right, (list, tuple)) else [C_right]
    assert len(lefts) == len(rights)

    jk = _build_complex_jk()
    jk.C_clear()
    for Cl, Cr in zip(lefts, rights):
        jk.C_left_add(Cl)
        jk.C_right_add(Cr)
    jk.compute_D()
    assert len(jk.D()) == len(lefts)
    return jk.D()


def test_complexjk_compute_D_c1_rectangular():
    """C1 non-square C: D = C @ C.conj().T."""
    rng = np.random.default_rng(0)
    C_arr = _random_complex((5, 2), rng)
    C = psi4.core.ComplexMatrix.from_array(C_arr, name="C")

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), C_arr @ C_arr.conj().T)


def test_complexjk_compute_D_c1_asymmetric():
    """C1 non-square C_left != C_right: D = C_left @ C_right.conj().T."""
    rng = np.random.default_rng(2)
    Cl_arr = _random_complex((5, 2), rng)
    Cr_arr = _random_complex((5, 2), rng)
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), Cl_arr @ Cr_arr.conj().T)
    # Sanity: not the same as the symmetric products
    assert not np.allclose(D.to_array(), Cl_arr @ Cl_arr.conj().T)
    assert not np.allclose(D.to_array(), Cr_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_c1_rectangular_D():
    """C_left/C_right with different AO counts yield rectangular D = Cl @ Cr.H."""
    rng = np.random.default_rng(4)
    Cl_arr = _random_complex((5, 2), rng)  # 5 AOs × 2 occ
    Cr_arr = _random_complex((3, 2), rng)  # 3 AOs × 2 occ
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    D_arr = D.to_array()
    assert D_arr.shape == (5, 3)
    np.testing.assert_allclose(D_arr, Cl_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_multi_irrep_rectangular():
    """Multi-irrep non-square C: D_h = C_h @ C_h.conj().T per diagonal tile."""
    rng = np.random.default_rng(1)
    row_sizes = [4, 2, 3]
    col_sizes = [1, 2, 1]
    C = psi4.core.ComplexMatrix("C", row_sizes, col_sizes)

    refs = []
    for h, view in enumerate(C.array_interface()):
        C_h = _random_complex(view.shape, rng)
        view[:] = C_h
        refs.append(C_h @ C_h.conj().T)

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == len(row_sizes)
    D_blocks = D.to_array()
    assert isinstance(D_blocks, list)
    assert len(D_blocks) == len(row_sizes)
    for got, ref, nso in zip(D_blocks, refs, row_sizes):
        assert got.shape == (nso, nso)
        np.testing.assert_allclose(got, ref)


def test_complexjk_compute_D_multiple_C_pairs():
    """Multiple queued C_left/C_right pairs each produce their own D = Cl @ Cr.H."""
    rng = np.random.default_rng(3)
    shapes = [(4, 1), (5, 2), (3, 2)]
    lefts, rights, refs = [], [], []
    for i, shape in enumerate(shapes):
        Cl_arr = _random_complex(shape, rng)
        Cr_arr = _random_complex(shape, rng)
        lefts.append(psi4.core.ComplexMatrix.from_array(Cl_arr, name=f"Cl{i}"))
        rights.append(psi4.core.ComplexMatrix.from_array(Cr_arr, name=f"Cr{i}"))
        refs.append(Cl_arr @ Cr_arr.conj().T)

    Ds = _run_compute_D(lefts, rights)
    assert len(Ds) == len(shapes)
    for D, ref in zip(Ds, refs):
        assert D.num_blocks() == 1
        np.testing.assert_allclose(D.to_array(), ref)


def _h2_sto3g():
    """C1 H2 / STO-3G molecule and orbital basis."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    return mol, basis


def _jk_reference(basis, D):
    """J_pq = sum_rs D_rs (pq|rs), K_pr = sum_qs D_qs (pq|rs) via MintsHelper.ao_eri."""
    mints = psi4.core.MintsHelper(basis)
    I = np.asarray(mints.ao_eri())
    J = np.einsum("pqrs,rs->pq", I, D, optimize=True)
    K = np.einsum("pqrs,qs->pr", I, D, optimize=True)
    return J, K


def test_complexdirectjk_matches_einsum():
    """ComplexDirectJK J/K match explicit einsum contractions of ao_eri with complex D."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(7)
    C_arr = _random_complex((nbf, 1), rng)
    D_ref = C_arr @ C_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    jk.initialize()
    jk.C_clear()
    jk.C_add(psi4.core.ComplexMatrix.from_array(C_arr, name="C"))
    jk.compute()
    jk.finalize()

    np.testing.assert_allclose(jk.D()[0].to_array(), D_ref)
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-10)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-10)


def test_complexdirectjk_asymmetric_matches_einsum():
    """ComplexDirectJK with C_left != C_right still matches einsum (general complex D)."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(11)
    Cl_arr = _random_complex((nbf, 1), rng)
    Cr_arr = _random_complex((nbf, 1), rng)
    D_ref = Cl_arr @ Cr_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    jk.initialize()
    jk.C_clear()
    jk.C_left_add(psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl"))
    jk.C_right_add(psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr"))
    jk.compute()
    jk.finalize()

    np.testing.assert_allclose(jk.D()[0].to_array(), D_ref)
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-10)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-10)


def test_complexdirectjk_real_D_matches_jk():
    """Completely real D: ComplexDirectJK agrees with real Direct JK."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(13)
    C_arr = rng.normal(size=(nbf, 1))

    cjk = psi4.core.ComplexJK.build_JK(basis, basis)
    cjk.initialize()
    cjk.C_clear()
    cjk.C_add(psi4.core.ComplexMatrix.from_array(C_arr.astype(np.complex128), name="C"))
    cjk.compute()
    cjk.finalize()

    rjk = psi4.core.JK.build_JK(basis, basis)
    rjk.initialize()
    rjk.C_clear()
    rjk.C_add(psi4.core.Matrix.from_array(C_arr, name="C"))
    rjk.compute()
    rjk.finalize()

    np.testing.assert_allclose(cjk.D()[0].to_array().real, np.asarray(rjk.D()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.J()[0].to_array().real, np.asarray(rjk.J()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.K()[0].to_array().real, np.asarray(rjk.K()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.J()[0].to_array().imag, 0.0, atol=1e-10)
    np.testing.assert_allclose(cjk.K()[0].to_array().imag, 0.0, atol=1e-10)


def test_complexdirectjk_ghf_spin_blocked_matches_reference():
    """GHF spin-blocked density (2*nbf x 2*nbf) exercises the generalized branch of
    ComplexDirectJK::compute_JK.

    For the ordinary spin-independent Coulomb operator: J only depends on the
    spin-traced total density and is spin-diagonal (J_aa = J_bb = J[D_aa+D_bb],
    J_ab = J_ba = 0); K couples each of the four spin blocks independently
    through the *same* spatial two-electron integrals (K_{sigma,tau} =
    K[D_{sigma,tau}]). This previously silently ignored everything past the
    top-left nbf x nbf (alpha-alpha) block, dropping the beta-beta (and any
    spin-mixing) contribution to the Fock matrix entirely.
    """
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(23)
    # 2 occupied spinors spanning both spin blocks with general complex mixing.
    C_arr = _random_complex((2 * nbf, 2), rng)
    D_ref = C_arr @ C_arr.conj().T

    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    jk.initialize()
    jk.C_clear()
    jk.C_add(psi4.core.ComplexMatrix.from_array(C_arr, name="C"))
    jk.compute()
    jk.finalize()

    D_aa = D_ref[:nbf, :nbf]
    D_bb = D_ref[nbf:, nbf:]
    D_ab = D_ref[:nbf, nbf:]
    D_ba = D_ref[nbf:, :nbf]

    J_tot, _ = _jk_reference(basis, D_aa + D_bb)
    _, K_aa = _jk_reference(basis, D_aa)
    _, K_bb = _jk_reference(basis, D_bb)
    _, K_ab = _jk_reference(basis, D_ab)
    _, K_ba = _jk_reference(basis, D_ba)

    J_ref = np.zeros_like(D_ref)
    K_ref = np.zeros_like(D_ref)
    J_ref[:nbf, :nbf] = J_tot
    J_ref[nbf:, nbf:] = J_tot
    K_ref[:nbf, :nbf] = K_aa
    K_ref[nbf:, nbf:] = K_bb
    K_ref[:nbf, nbf:] = K_ab
    K_ref[nbf:, :nbf] = K_ba

    np.testing.assert_allclose(jk.D()[0].to_array(), D_ref, atol=1e-12)
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-10)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-10)

    # The historical bug left the beta-beta block (and everything else past
    # the top-left nbf x nbf tile) at exactly zero.
    assert not np.allclose(jk.K()[0].to_array()[nbf:, nbf:], 0.0)
    assert not np.allclose(jk.J()[0].to_array()[nbf:, nbf:], 0.0)


def _run_complex_jk(basis, C_arr, screening="NONE", ints_tol=1.0e-12, bench=False):
    """Build/initialize/compute/finalize ComplexDirectJK for a single C."""
    psi4.set_options(
        {
            "SCF_TYPE": "DIRECT",
            "SCREENING": screening,
            "INTS_TOLERANCE": ints_tol,
            "DF_SCF_GUESS": False,
        }
    )
    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    if bench:
        jk.set_bench(1)
    jk.initialize()
    jk.C_clear()
    jk.C_add(psi4.core.ComplexMatrix.from_array(np.asarray(C_arr, dtype=np.complex128), name="C"))
    jk.compute()
    jk.finalize()
    return jk


def _max_eri_am():
    """Highest angular momentum with 4-center ERI support in this libint2 build."""
    amchar = "SPDFGHIKLM"
    max_am = 0
    for am, ch in enumerate(amchar):
        # libint2 uses lowercase shell labels in support keys, e.g. eri_kkkk_d0
        key = f"eri_{ch.lower() * 4}_d0"
        if psi4.core.libint2_supports(key):
            max_am = am
        else:
            break
    return max_am


def _uncontracted_am_basis_string(max_am):
    """Spherical custom basis: one primitive per AM from S through max_am."""
    amchar = "SPDFGHIKLM"
    lines = ["spherical", "****", "H 0"]
    for am in range(max_am + 1):
        lines.append(f"{amchar[am]} 1 1.0")
        lines.append("  1.0 1.0")
    lines.append("****")
    return "\n".join(lines)


def test_complexdirectjk_high_am_single_atom():
    """Single atom with one primitive per shell up through max supported ERI AM."""
    max_am = _max_eri_am()
    assert max_am >= 4  # at least G in a normal Psi4 build
    basis_string = _uncontracted_am_basis_string(max_am)

    mol = psi4.geometry(
        """
        0 1
        He
        symmetry c1
        """
    )

    def basisspec_psi4_yo__anonymous_higham(mol, role):
        mol.set_basis_all_atoms("higham", role=role)
        # One-primitive-per-AM string is written for H; reuse label for He.
        return {"higham": basis_string.replace("H 0", "He 0")}

    psi4.driver.qcdb.libmintsbasisset.basishorde["HIGHAM"] = basisspec_psi4_yo__anonymous_higham
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "HIGHAM")

    assert basis.max_am() == max_am
    assert basis.nshell() == max_am + 1
    nbf = basis.nbf()
    # spherical: sum_{l=0}^{L} (2l+1) = (L+1)^2
    assert nbf == (max_am + 1) ** 2

    rng = np.random.default_rng(17)
    C_arr = _random_complex((nbf, 1), rng)
    D_ref = C_arr @ C_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk = _run_complex_jk(basis, C_arr, screening="NONE")
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-8)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-8)


def test_complexdirectjk_h2_dimer_screening():
    """H2···H2 at ~10 Å: Schwarz screening must match the full einsum reference and skip work."""
    mol = psi4.geometry(
        """
        0 1
        H  0.0  0.0   0.00
        H  0.0  0.0   0.74
        H  0.0  0.0  10.00
        H  0.0  0.0  10.74
        units angstrom
        symmetry c1
        """
    )

    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    assert mol.natom() == 4
    assert basis.nshell() == 4

    nbf = basis.nbf()
    rng = np.random.default_rng(19)
    C_arr = _random_complex((nbf, 2), rng)
    D_ref = C_arr @ C_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk_none = _run_complex_jk(basis, C_arr, screening="NONE", bench=True)
    n_none = jk_none.computed_shells_per_iter()["Quartets"][-1]

    jk_schwarz = _run_complex_jk(basis, C_arr, screening="SCHWARZ", ints_tol=1.0e-10, bench=True)
    n_schwarz = jk_schwarz.computed_shells_per_iter()["Quartets"][-1]

    np.testing.assert_allclose(jk_schwarz.J()[0].to_array(), J_ref, atol=1e-8)
    np.testing.assert_allclose(jk_schwarz.K()[0].to_array(), K_ref, atol=1e-8)
    assert n_schwarz < n_none
    # Without uniqueness, unscreened count is nshell^4
    assert n_none == basis.nshell() ** 4
