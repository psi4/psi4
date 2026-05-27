"""Tests for the psi4.driver.trexio save/load interface."""

import os

import numpy as np
import pytest

import psi4

trexio = pytest.importorskip("trexio")
from psi4.driver import trexio as p4trex  # noqa: E402

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.trexio]


def _pick_back_end():
    """Pick whichever TrexIO backend is compiled in (prefer HDF5)."""
    import tempfile
    for be, suffix in ((trexio.TREXIO_HDF5, ".h5"), (trexio.TREXIO_TEXT, "")):
        d = tempfile.mkdtemp()
        try:
            f = trexio.File(d + "/probe" + suffix, mode="w", back_end=be)
            f.close()
            return be, suffix
        except trexio.Error:
            continue
    pytest.skip("No usable TrexIO back end available")


_BACK_END, _SUFFIX = _pick_back_end()


def _save(wfn, tmp_path, name, **kwargs):
    kwargs.setdefault("back_end", _BACK_END)
    path = str(tmp_path / (name + _SUFFIX))
    p4trex.save(wfn, path, **kwargs)
    return path


def _h2o_rhf(scf_type="pk"):
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    psi4.geometry(
        """
        0 1
        O
        H 1 0.96
        H 1 0.96 2 104.5
        symmetry c1
        """
    )
    psi4.set_options(
        {"basis": "cc-pvdz", "scf_type": scf_type, "e_convergence": 10, "d_convergence": 10}
    )
    return psi4.energy("hf", return_wfn=True)


def _h2o_uhf():
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    psi4.geometry(
        """
        1 2
        O
        H 1 0.96
        H 1 0.96 2 104.5
        symmetry c1
        """
    )
    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "e_convergence": 10,
            "d_convergence": 10,
        }
    )
    return psi4.energy("hf", return_wfn=True)


def test_save_roundtrip_rhf(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_rhf", description="H2O RHF/cc-pVDZ")

    data = p4trex.load(path)

    assert data["nucleus"]["num"] == 3
    assert data["ao"]["num"] == wfn.basisset().nbf()
    assert data["mo"]["num"] == wfn.basisset().nbf()
    assert data["electron"]["up_num"] == wfn.nalpha()
    assert data["electron"]["dn_num"] == wfn.nbeta()
    assert data["mo"]["type"] == "RHF"

    nbf = wfn.basisset().nbf()
    # Canonical AO ordering: convert Psi4 Ca to canonical for comparison.
    from psi4.driver.trexio import _psi4_to_canonical_perm, _canonical_ao_normalization
    perm = _psi4_to_canonical_perm(wfn.basisset())
    ao_norm = _canonical_ao_normalization(wfn.basisset())
    Ca_psi = np.asarray(wfn.Ca())
    Ca_can = np.empty_like(Ca_psi)
    Ca_can[perm] = Ca_psi
    Ca_can = Ca_can / ao_norm[:, None]

    Cmat = np.asarray(data["mo"]["coefficient"]).reshape(nbf, nbf).T
    np.testing.assert_allclose(np.abs(Cmat), np.abs(Ca_can), atol=1e-12)

    eps = np.asarray(data["mo"]["energy"])
    np.testing.assert_allclose(eps, np.asarray(wfn.epsilon_a()), atol=1e-12)

    occ = np.asarray(data["mo"]["occupation"])
    assert np.isclose(occ.sum(), wfn.nalpha() + wfn.nbeta())


def test_save_uhf(tmp_path):
    e, wfn = _h2o_uhf()
    path = _save(wfn, tmp_path, "h2o_uhf")

    data = p4trex.load(path)
    nbf = wfn.basisset().nbf()
    assert data["mo"]["type"] == "UHF"
    assert data["mo"]["num"] == 2 * nbf
    assert data["electron"]["up_num"] == wfn.nalpha()
    assert data["electron"]["dn_num"] == wfn.nbeta()

    spin = np.asarray(data["mo"]["spin"])
    assert (spin[:nbf] == 0).all()
    assert (spin[nbf:] == 1).all()


def test_ao_one_electron_integrals(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_ao1e", save_ao_integrals=True)

    data = p4trex.load(path)
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = np.asarray(data["ao_1e"]["overlap"])
    T = np.asarray(data["ao_1e"]["kinetic"])
    V = np.asarray(data["ao_1e"]["potential_n_e"])
    H = np.asarray(data["ao_1e"]["core_hamiltonian"])

    # Canonical overlap must be unit-diagonal (the whole point of the canonical
    # ao_normalization factor).
    np.testing.assert_allclose(np.diag(S), np.ones(S.shape[0]), atol=1e-12)
    np.testing.assert_allclose(H, T + V, atol=1e-12)

    # Apply the canonical → Psi4 transform and compare to Psi4's own matrices.
    from psi4.driver.trexio import _psi4_to_canonical_perm, _canonical_ao_normalization
    perm = _psi4_to_canonical_perm(wfn.basisset())
    ao_norm = _canonical_ao_normalization(wfn.basisset())

    def _to_psi4(M_can):
        # M_psi4[i,j] = M_can[perm[i], perm[j]] / (ao_norm[perm[i]] * ao_norm[perm[j]])
        M = M_can[np.ix_(perm, perm)]
        scale = ao_norm[perm]
        return M / (scale[:, None] * scale[None, :])

    np.testing.assert_allclose(_to_psi4(S), np.asarray(mints.ao_overlap()), atol=1e-12)
    np.testing.assert_allclose(_to_psi4(T), np.asarray(mints.ao_kinetic()), atol=1e-12)
    np.testing.assert_allclose(_to_psi4(V), np.asarray(mints.ao_potential()), atol=1e-12)


def test_mo_one_electron_integrals(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_mo1e", save_mo_integrals=True)

    data = p4trex.load(path)
    mints = psi4.core.MintsHelper(wfn.basisset())
    Ca = np.asarray(wfn.Ca())
    S_mo_ref = Ca.T @ np.asarray(mints.ao_overlap()) @ Ca
    np.testing.assert_allclose(
        np.asarray(data["mo_1e"]["overlap"]), S_mo_ref, atol=1e-10
    )
    # MO overlap should be identity
    np.testing.assert_allclose(
        np.asarray(data["mo_1e"]["overlap"]), np.eye(Ca.shape[1]), atol=1e-10
    )


def test_ao_eri_sparse(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_eri", save_eri=True)

    data = p4trex.load(path)
    assert "ao_2e" in data
    idx = data["ao_2e"]["indices"]
    val = data["ao_2e"]["values"]
    assert idx.ndim == 2 and idx.shape[1] == 4
    assert val.ndim == 1
    assert idx.shape[0] == val.shape[0]

    mints = psi4.core.MintsHelper(wfn.basisset())
    full = np.asarray(mints.ao_eri())
    nbf = wfn.basisset().nbf()

    # Reconstruct the canonical ERI from the sparse list.
    rebuilt = np.zeros_like(full)
    for (i, j, k, l), v in zip(idx, val):
        for ii, jj in ((i, j), (j, i)):
            for kk, ll in ((k, l), (l, k)):
                rebuilt[ii, jj, kk, ll] = v
                rebuilt[kk, ll, ii, jj] = v

    # Convert canonical back to Psi4 order and compare.
    from psi4.driver.trexio import _psi4_to_canonical_perm, _canonical_ao_normalization
    perm = _psi4_to_canonical_perm(wfn.basisset())
    ao_norm = _canonical_ao_normalization(wfn.basisset())
    inv_scale = 1.0 / ao_norm[perm]
    rebuilt_psi = rebuilt[np.ix_(perm, perm, perm, perm)] * (
        inv_scale[:, None, None, None] * inv_scale[None, :, None, None]
        * inv_scale[None, None, :, None] * inv_scale[None, None, None, :]
    )
    np.testing.assert_allclose(rebuilt_psi, full, atol=1e-12)


def test_mo_eri_sparse(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_mo_eri", save_mo_eri=True)

    data = p4trex.load(path)
    assert "mo_2e" in data
    idx = data["mo_2e"]["indices"]
    val = data["mo_2e"]["values"]

    Ca = np.asarray(wfn.Ca())
    mints = psi4.core.MintsHelper(wfn.basisset())
    ao_eri = np.asarray(mints.ao_eri())
    mo_eri_ref = np.einsum("pqrs,pi,qj,rk,sl->ijkl", ao_eri, Ca, Ca, Ca, Ca, optimize=True)

    mo_eri = np.zeros_like(mo_eri_ref)
    for (i, j, k, l), v in zip(idx, val):
        for ii, jj in ((i, j), (j, i)):
            for kk, ll in ((k, l), (l, k)):
                mo_eri[ii, jj, kk, ll] = v
                mo_eri[kk, ll, ii, jj] = v
    np.testing.assert_allclose(mo_eri, mo_eri_ref, atol=1e-10)


def test_load_wavefunction_reconstructs_basis(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_basis_reco")
    # No basis_name → rebuild basis from the file's basis block.
    rebuilt = p4trex.load_wavefunction(path)
    assert rebuilt.basisset().nbf() == wfn.basisset().nbf()
    np.testing.assert_allclose(
        np.asarray(rebuilt.epsilon_a()), np.asarray(wfn.epsilon_a()), atol=1e-12
    )


def test_overwrite_protection(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_ow")
    with pytest.raises(p4trex.TrexIOError):
        p4trex.save(wfn, path, back_end=_BACK_END)
    p4trex.save(wfn, path, back_end=_BACK_END, overwrite=True)  # must succeed


def test_load_wavefunction(tmp_path):
    e, wfn = _h2o_rhf()
    path = _save(wfn, tmp_path, "h2o_rt")

    rebuilt = p4trex.load_wavefunction(path, reference="rhf", basis_name="cc-pVDZ")
    assert rebuilt.basisset().nbf() == wfn.basisset().nbf()
    assert rebuilt.molecule().natom() == wfn.molecule().natom()
    np.testing.assert_allclose(
        np.asarray(rebuilt.epsilon_a()), np.asarray(wfn.epsilon_a()), atol=1e-12
    )
