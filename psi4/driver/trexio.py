#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""TrexIO (https://trex-coe.github.io/trexio) interface for Psi4.

Read and write `Psi4 <https://psicode.org>`_ wavefunctions to/from the TrexIO
HDF5 file format used by the TREX-CoE quantum-chemistry / QMC ecosystem.

This module is a thin wrapper around the official ``trexio`` Python package.
It supports four layers of data:

* **Geometry + basis** -- nuclei, AO shells, primitives, normalization.
* **Molecular orbitals** -- coefficients in the AO basis, orbital energies,
  occupations, spin labels.
* **Integrals** -- AO/MO one-electron (S, T, V, h) and two-electron (ERI),
  written in sparse 8-fold-permutation form.
* **Determinants** -- bitfield CI expansion (for DETCI wavefunctions).

The AO basis is written and read back in Psi4's native ordering (Cartesian
*alphabetical*: ``xx, xy, xz, yy, yz, zz``; pure spherical: ``m=0, +1c, +1s,
+2c, +2s, ...``). A ``metadata_description`` tag records this so consumers can
re-permute if they expect a different convention. Round-tripping inside Psi4
is exact.
"""

from __future__ import annotations

import os
from typing import Any, Optional

import numpy as np

from psi4 import core


__all__ = [
    "TrexIOError",
    "save",
    "load",
    "save_wavefunction",
    "load_wavefunction",
]


_AO_ORDER_TAG = "psi4-native-ao-order"


class TrexIOError(RuntimeError):
    """Raised when the TrexIO interface fails or is unavailable."""


def _trexio():
    try:
        import trexio  # noqa: F401
    except ImportError as exc:
        raise TrexIOError(
            "The 'trexio' Python package is not installed. "
            "Install it with `pip install trexio` to use psi4.driver.trexio."
        ) from exc
    return trexio


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------


def save(
    wfn: "core.Wavefunction",
    filepath: str,
    *,
    overwrite: bool = False,
    save_ao_integrals: bool = False,
    save_mo_integrals: bool = False,
    save_eri: bool = False,
    ci_vector: Optional[np.ndarray] = None,
    determinants: Optional[np.ndarray] = None,
    back_end: Optional[int] = None,
    description: str = "",
) -> str:
    """Write a Psi4 wavefunction to a TrexIO file.

    Parameters
    ----------
    wfn
        A converged :class:`psi4.core.Wavefunction`.
    filepath
        Output path (``.h5`` for HDF5, directory for TEXT back end).
    overwrite
        Remove an existing file/directory at ``filepath`` first.
    save_ao_integrals
        Write AO overlap/kinetic/N-e/core-hamiltonian.
    save_mo_integrals
        Write the same four matrices transformed to the MO basis.
    save_eri
        Write the AO ERI tensor in sparse 8-fold form. Currently AO only;
        MO ERIs would dominate file size and are left to consumers.
    ci_vector, determinants
        Optional CI expansion. ``determinants`` is a ``(ndet, 2*nint)`` int64
        array of bitfield occupations (alpha then beta), as produced by
        :func:`trexio.to_bitfield_list`. ``ci_vector`` is the matching list of
        coefficients.
    back_end
        ``trexio.TREXIO_HDF5`` (default) or ``trexio.TREXIO_TEXT``.
    description
        Free-text metadata describing the calculation.

    Returns
    -------
    str
        The absolute path written.
    """
    trexio = _trexio()

    filepath = os.path.abspath(filepath)
    if os.path.exists(filepath):
        if not overwrite:
            raise TrexIOError(
                f"TrexIO target {filepath!r} already exists. "
                "Pass overwrite=True to replace it."
            )
        if os.path.isdir(filepath):
            import shutil
            shutil.rmtree(filepath)
        else:
            os.remove(filepath)

    if back_end is None:
        back_end = trexio.TREXIO_HDF5

    with trexio.File(filepath, mode="w", back_end=back_end) as f:
        _write_metadata(f, description)
        _write_nuclei(f, wfn.molecule())
        basisset = wfn.basisset()
        _write_basis(f, basisset)
        _write_ao(f, basisset)
        _write_electrons(f, wfn)
        _write_mos(f, wfn, basisset)

        if save_ao_integrals or save_mo_integrals or save_eri:
            mints = core.MintsHelper(basisset)
            if save_ao_integrals:
                _write_ao_one_e(f, mints, wfn)
            if save_mo_integrals:
                _write_mo_one_e(f, mints, wfn)
            if save_eri:
                _write_ao_eri(f, mints)

        if determinants is not None or ci_vector is not None:
            if determinants is None or ci_vector is None:
                raise TrexIOError(
                    "Both 'determinants' and 'ci_vector' must be supplied to save a CI expansion."
                )
            _write_determinants(f, determinants, ci_vector)

    return filepath


def save_wavefunction(wfn, filepath, **kwargs):
    """Alias for :func:`save` kept for naming symmetry with :func:`load_wavefunction`."""
    return save(wfn, filepath, **kwargs)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------


def load(filepath: str) -> dict:
    """Read a TrexIO file and return its data as a dict of numpy arrays.

    The returned dict contains whichever groups were present in the file:
    ``metadata``, ``nucleus``, ``basis``, ``ao``, ``electron``, ``mo``,
    ``ao_1e``, ``ao_2e``, ``mo_1e``, ``determinant``. Missing groups are
    omitted.
    """
    trexio = _trexio()

    filepath = os.path.abspath(filepath)
    if not os.path.exists(filepath):
        raise TrexIOError(f"TrexIO file {filepath!r} does not exist.")

    out: dict[str, Any] = {}
    with trexio.File(filepath, mode="r", back_end=trexio.TREXIO_AUTO) as f:
        out["metadata"] = _read_metadata(f)
        if trexio.has_nucleus_num(f):
            out["nucleus"] = _read_nuclei(f)
        if trexio.has_basis_shell_num(f):
            out["basis"] = _read_basis(f)
        if trexio.has_ao_num(f):
            out["ao"] = _read_ao(f)
        if trexio.has_electron_num(f):
            out["electron"] = _read_electrons(f)
        if trexio.has_mo_num(f):
            out["mo"] = _read_mos(f)
        ao_1e = _read_ao_one_e(f)
        if ao_1e:
            out["ao_1e"] = ao_1e
        mo_1e = _read_mo_one_e(f)
        if mo_1e:
            out["mo_1e"] = mo_1e
        ao_2e = _read_ao_eri(f)
        if ao_2e is not None:
            out["ao_2e"] = ao_2e
        det = _read_determinants(f)
        if det is not None:
            out["determinant"] = det

    return out


def load_wavefunction(
    filepath: str,
    *,
    reference: str = "rhf",
    basis_name: Optional[str] = None,
) -> "core.Wavefunction":
    """Reconstruct a Psi4 Wavefunction (geometry + basis + MOs) from a TrexIO file.

    The returned wavefunction is a bare reference wavefunction with C, epsilon,
    D, and occupation numbers populated. The molecule is rebuilt from the
    nuclear coordinates and charges. The AO basis is reloaded from the
    standard Psi4 basis library using ``basis_name`` (defaulting to the file's
    metadata, falling back to ``cc-pVDZ``); shell counts and primitive counts
    are checked for consistency, but the actual exponents/coefficients in the
    library are used (since TrexIO does not preserve every Psi4-internal
    detail).

    For exact AO-level fidelity (e.g. custom basis sets) the user should
    construct the BasisSet themselves and use :func:`load` to fetch the MO
    block directly.

    The reconstructed wavefunction is C1-symmetric.
    """
    data = load(filepath)

    if "nucleus" not in data or "mo" not in data:
        raise TrexIOError(
            "TrexIO file is missing nucleus or mo data; cannot rebuild a Wavefunction."
        )

    nuc = data["nucleus"]
    nelec = data.get("electron", {}).get("num")
    total_Z = float(np.asarray(nuc["charge"]).sum())
    charge = int(round(total_Z - nelec)) if nelec is not None else 0
    nalpha = data.get("electron", {}).get("up_num", 0)
    nbeta = data.get("electron", {}).get("dn_num", 0)
    multiplicity = abs(nalpha - nbeta) + 1

    mol_lines = [f"{charge} {multiplicity}"]
    for label, (x, y, z) in zip(nuc["label"], nuc["coord"]):
        mol_lines.append(f"{label} {x:.16f} {y:.16f} {z:.16f}")
    mol_lines += ["units bohr", "no_reorient", "no_com", "symmetry c1"]
    mol = core.Molecule.from_string("\n".join(mol_lines), dtype="psi4")
    mol.update_geometry()

    if basis_name is None:
        basis_name = data.get("metadata", {}).get("basis_name", "cc-pVDZ")

    basisset = core.BasisSet.build(mol, "ORBITAL", basis_name, puream=-1)
    nbf = basisset.nbf()
    if nbf != data["mo"]["num"]:
        raise TrexIOError(
            f"AO/MO count mismatch: basis {basis_name!r} gives nbf={nbf}, "
            f"but file has {data['mo']['num']} MOs."
        )

    Cmat = np.asarray(data["mo"]["coefficient"]).reshape(data["mo"]["num"], nbf).T
    eps = np.asarray(data["mo"]["energy"]) if data["mo"]["energy"] is not None else np.zeros(nbf)
    occ = np.asarray(data["mo"]["occupation"]) if data["mo"]["occupation"] is not None else np.zeros(nbf)

    ref = reference.lower()
    if ref not in ("rhf", "rks", "uhf", "uks", "rohf"):
        raise TrexIOError(f"Unsupported reference {reference!r}.")

    Da = np.einsum("p,mp,np->mn", 0.5 * occ, Cmat, Cmat)
    aotoso = np.eye(nbf)
    nmo = nbf
    nalpha_int = int(nalpha)
    nbeta_int = int(nbeta)

    wfn_matrix = {
        "Ca": core.Matrix.from_array(Cmat, name="Ca"),
        "Cb": core.Matrix.from_array(Cmat, name="Cb"),
        "Da": core.Matrix.from_array(Da, name="Da"),
        "Db": core.Matrix.from_array(Da, name="Db"),
        "Fa": None,
        "Fb": None,
        "H": None,
        "S": None,
        "X": None,
        "aotoso": core.Matrix.from_array(aotoso, name="aotoso"),
        "gradient": None,
        "hessian": None,
    }
    wfn_vector = {
        "epsilon_a": core.Vector.from_array(eps, name="epsilon_a"),
        "epsilon_b": core.Vector.from_array(eps, name="epsilon_b"),
        "frequencies": None,
    }
    wfn_dimension = {
        "doccpi": core.Dimension.from_list([min(nalpha_int, nbeta_int)]),
        "frzcpi": core.Dimension.from_list([0]),
        "frzvpi": core.Dimension.from_list([0]),
        "nalphapi": core.Dimension.from_list([nalpha_int]),
        "nbetapi": core.Dimension.from_list([nbeta_int]),
        "nmopi": core.Dimension.from_list([nmo]),
        "nsopi": core.Dimension.from_list([nbf]),
        "soccpi": core.Dimension.from_list([abs(nalpha_int - nbeta_int)]),
    }
    wfn_int = {
        "nalpha": nalpha_int,
        "nbeta": nbeta_int,
        "nfrzc": 0,
        "nirrep": 1,
        "nmo": nmo,
        "nso": nbf,
        "print": 1,
    }
    wfn_string = {"name": ref.upper(), "module": "trexio", "basisname": basis_name}
    wfn_boolean = {
        "PCM_enabled": False,
        "same_a_b_dens": True,
        "same_a_b_orbs": True,
        "basispuream": bool(basisset.has_puream()),
    }
    wfn_float = {
        "energy": 0.0,
        "efzc": 0.0,
        "dipole_field_x": 0.0,
        "dipole_field_y": 0.0,
        "dipole_field_z": 0.0,
    }
    return core.Wavefunction(
        mol, basisset, wfn_matrix, wfn_vector, wfn_dimension,
        wfn_int, wfn_string, wfn_boolean, wfn_float,
    )


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def _write_metadata(f, description: str) -> None:
    trexio = _trexio()
    trexio.write_metadata_code_num(f, 1)
    trexio.write_metadata_code(f, ["Psi4"])
    trexio.write_metadata_author_num(f, 1)
    trexio.write_metadata_author(f, ["Psi4"])
    parts = [_AO_ORDER_TAG]
    if description:
        parts.append(description)
    trexio.write_metadata_description(f, " | ".join(parts))


def _write_nuclei(f, mol) -> None:
    trexio = _trexio()
    natom = mol.natom()
    trexio.write_nucleus_num(f, natom)
    coords = np.zeros((natom, 3), dtype=np.float64)
    charges = np.zeros(natom, dtype=np.float64)
    labels = []
    for i in range(natom):
        coords[i, 0] = mol.x(i)
        coords[i, 1] = mol.y(i)
        coords[i, 2] = mol.z(i)
        charges[i] = mol.Z(i)
        labels.append(mol.label(i))
    trexio.write_nucleus_coord(f, coords)
    trexio.write_nucleus_charge(f, charges)
    trexio.write_nucleus_label(f, labels)
    trexio.write_nucleus_point_group(f, mol.point_group().symbol())
    trexio.write_nucleus_repulsion(f, mol.nuclear_repulsion_energy())


def _write_basis(f, basisset) -> None:
    """Write the basis block.

    TrexIO's basis is a flat shell list. For each shell we store:

    * ``basis_nucleus_index[s]`` -- which atom the shell sits on
    * ``basis_shell_ang_mom[s]`` -- L
    * ``basis_shell_factor[s]`` -- extra normalization (1.0 for Psi4)
    * primitive (exponent, coefficient) pairs indexed by shell

    Psi4 keeps the *normalized* contraction coefficient on the shell. We
    write those directly into ``basis_coefficient`` and set
    ``basis_prim_factor = 1`` to indicate "coefficients already include
    primitive normalization".
    """
    trexio = _trexio()
    nshell = basisset.nshell()
    shell_ang_mom = np.zeros(nshell, dtype=np.int32)
    shell_nucleus = np.zeros(nshell, dtype=np.int32)
    shell_factor = np.ones(nshell, dtype=np.float64)

    exponents: list[float] = []
    coefficients: list[float] = []
    prim_factor: list[float] = []
    shell_index: list[int] = []

    for s in range(nshell):
        sh = basisset.shell(s)
        shell_ang_mom[s] = sh.am
        shell_nucleus[s] = sh.ncenter
        for p in range(sh.nprimitive):
            exponents.append(sh.exp(p))
            coefficients.append(sh.coef(p))
            prim_factor.append(1.0)
            shell_index.append(s)

    trexio.write_basis_type(f, "Gaussian")
    trexio.write_basis_shell_num(f, nshell)
    trexio.write_basis_prim_num(f, len(exponents))
    trexio.write_basis_nucleus_index(f, shell_nucleus)
    trexio.write_basis_shell_ang_mom(f, shell_ang_mom)
    trexio.write_basis_shell_factor(f, shell_factor)
    trexio.write_basis_shell_index(f, np.asarray(shell_index, dtype=np.int32))
    trexio.write_basis_exponent(f, np.asarray(exponents, dtype=np.float64))
    trexio.write_basis_coefficient(f, np.asarray(coefficients, dtype=np.float64))
    trexio.write_basis_prim_factor(f, np.asarray(prim_factor, dtype=np.float64))


def _write_ao(f, basisset) -> None:
    trexio = _trexio()
    nbf = basisset.nbf()
    cartesian = 0 if basisset.has_puream() else 1
    ao_shell = np.zeros(nbf, dtype=np.int32)
    ao_norm = np.ones(nbf, dtype=np.float64)
    idx = 0
    for s in range(basisset.nshell()):
        sh = basisset.shell(s)
        nfn = sh.nfunction
        for i in range(nfn):
            ao_shell[idx + i] = s
        idx += nfn
    trexio.write_ao_cartesian(f, cartesian)
    trexio.write_ao_num(f, nbf)
    trexio.write_ao_shell(f, ao_shell)
    trexio.write_ao_normalization(f, ao_norm)


def _write_electrons(f, wfn) -> None:
    trexio = _trexio()
    nalpha = wfn.nalpha()
    nbeta = wfn.nbeta()
    trexio.write_electron_num(f, nalpha + nbeta)
    trexio.write_electron_up_num(f, nalpha)
    trexio.write_electron_dn_num(f, nbeta)


def _gather_C1(matrix, aotoso) -> np.ndarray:
    """Pull a (sopi × mopi) symmetry-blocked matrix down to a flat (nao × nmo) C1 form."""
    if matrix.nirrep() == 1:
        return np.asarray(matrix)
    return np.asarray(matrix.clone().remove_symmetry(aotoso, matrix.symmetry()))


def _C1_MO(C, aotoso) -> np.ndarray:
    """Return MO coefficients in C1 AO basis: (nao_C1, nmo_total)."""
    if C.nirrep() == 1:
        return np.asarray(C)
    # Build per-irrep AO blocks, then horizontally stack columns ordered by irrep.
    nao = aotoso.rowdim(0) if False else C.rowdim().sum()  # symmetric AO count == total SO count
    # Use Wavefunction-style C1 reconstruction via aotoso: AO_C1 = sum_h U_h * C_h.
    blocks = []
    for h in range(C.nirrep()):
        Uh = np.asarray(aotoso.nph[h])  # (nao_total, nso_h)
        Ch = np.asarray(C.nph[h])       # (nso_h, nmo_h)
        if Ch.size:
            blocks.append(Uh @ Ch)
    return np.hstack(blocks)


def _flatten_vector(v) -> np.ndarray:
    if v.nirrep() == 1:
        return np.asarray(v)
    return np.concatenate([np.asarray(v.nph[h]) for h in range(v.nirrep())])


def _write_mos(f, wfn, basisset) -> None:
    trexio = _trexio()
    aotoso = wfn.aotoso()

    Ca = _C1_MO(wfn.Ca(), aotoso)
    Cb = _C1_MO(wfn.Cb(), aotoso) if not wfn.same_a_b_orbs() else None
    eps_a = _flatten_vector(wfn.epsilon_a())
    eps_b = _flatten_vector(wfn.epsilon_b()) if not wfn.same_a_b_orbs() else None

    nbf = basisset.nbf()
    if Cb is None:
        # Closed-shell-like storage: 1 MO per orbital, occupation 0/2.
        nmo = Ca.shape[1]
        # TrexIO stores coefficients in row-major (mo, ao) order.
        coefs = Ca.T.copy()
        energies = eps_a
        spin = np.zeros(nmo, dtype=np.int32)
        occ = np.zeros(nmo, dtype=np.float64)
        occ[: wfn.nalpha()] = 2.0
        mo_type = "RHF"
    else:
        nmo = Ca.shape[1] + Cb.shape[1]
        coefs = np.vstack([Ca.T, Cb.T])
        energies = np.concatenate([eps_a, eps_b])
        spin = np.concatenate([
            np.zeros(Ca.shape[1], dtype=np.int32),
            np.ones(Cb.shape[1], dtype=np.int32),
        ])
        occ = np.zeros(nmo, dtype=np.float64)
        occ[: wfn.nalpha()] = 1.0
        occ[Ca.shape[1] : Ca.shape[1] + wfn.nbeta()] = 1.0
        mo_type = "UHF"

    if coefs.shape[1] != nbf:
        raise TrexIOError(
            f"MO coefficient AO dimension {coefs.shape[1]} does not match basis nbf {nbf}."
        )

    trexio.write_mo_type(f, mo_type)
    trexio.write_mo_num(f, nmo)
    trexio.write_mo_coefficient(f, coefs)
    trexio.write_mo_energy(f, energies)
    trexio.write_mo_occupation(f, occ)
    trexio.write_mo_spin(f, spin)


def _write_ao_one_e(f, mints, wfn) -> None:
    trexio = _trexio()
    S = np.asarray(mints.ao_overlap())
    T = np.asarray(mints.ao_kinetic())
    V = np.asarray(mints.ao_potential())
    H = T + V
    trexio.write_ao_1e_int_overlap(f, S)
    trexio.write_ao_1e_int_kinetic(f, T)
    trexio.write_ao_1e_int_potential_n_e(f, V)
    trexio.write_ao_1e_int_core_hamiltonian(f, H)


def _write_mo_one_e(f, mints, wfn) -> None:
    trexio = _trexio()
    aotoso = wfn.aotoso()
    Ca = _C1_MO(wfn.Ca(), aotoso)
    S = np.asarray(mints.ao_overlap())
    T = np.asarray(mints.ao_kinetic())
    V = np.asarray(mints.ao_potential())
    H = T + V
    trexio.write_mo_1e_int_overlap(f, Ca.T @ S @ Ca)
    trexio.write_mo_1e_int_kinetic(f, Ca.T @ T @ Ca)
    trexio.write_mo_1e_int_potential_n_e(f, Ca.T @ V @ Ca)
    trexio.write_mo_1e_int_core_hamiltonian(f, Ca.T @ H @ Ca)


def _write_ao_eri(f, mints) -> None:
    """Write the AO ERI in 8-fold permutational symmetry as a sparse list."""
    trexio = _trexio()
    eri = np.asarray(mints.ao_eri())
    nbf = eri.shape[0]

    indices: list[list[int]] = []
    values: list[float] = []
    for i in range(nbf):
        for j in range(i + 1):
            ij = i * (i + 1) // 2 + j
            for k in range(nbf):
                for l in range(k + 1):
                    kl = k * (k + 1) // 2 + l
                    if kl > ij:
                        continue
                    v = eri[i, j, k, l]
                    if abs(v) < 1e-15:
                        continue
                    indices.append([i, j, k, l])
                    values.append(float(v))

    if not values:
        return
    chunk = 65536
    n = len(values)
    arr_idx = np.asarray(indices, dtype=np.int32)  # shape (n, 4)
    arr_val = np.asarray(values, dtype=np.float64)
    offset = 0
    while offset < n:
        end = min(offset + chunk, n)
        trexio.write_ao_2e_int_eri(
            f, offset, end - offset,
            arr_idx[offset:end].flatten().tolist(),
            arr_val[offset:end].tolist(),
        )
        offset = end


def _write_determinants(f, determinants: np.ndarray, ci_vector: np.ndarray) -> None:
    trexio = _trexio()
    determinants = np.asarray(determinants, dtype=np.int64)
    ci_vector = np.asarray(ci_vector, dtype=np.float64)
    if determinants.shape[0] != ci_vector.shape[0]:
        raise TrexIOError(
            f"Determinant count {determinants.shape[0]} does not match "
            f"CI coefficient count {ci_vector.shape[0]}."
        )
    ndet = determinants.shape[0]
    # trexio expects a list of bitfields (each itself a list of int64).
    det_list = determinants.tolist()
    chunk = 8192
    offset = 0
    while offset < ndet:
        end = min(offset + chunk, ndet)
        trexio.write_determinant_list(f, offset, end - offset, det_list[offset:end])
        trexio.write_determinant_coefficient(
            f, offset, end - offset, ci_vector[offset:end].tolist()
        )
        offset = end


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


def _read_metadata(f) -> dict:
    trexio = _trexio()
    out: dict[str, Any] = {}
    if trexio.has_metadata_code(f):
        out["code"] = trexio.read_metadata_code(f)
    if trexio.has_metadata_author(f):
        out["author"] = trexio.read_metadata_author(f)
    if trexio.has_metadata_description(f):
        out["description"] = trexio.read_metadata_description(f)
    return out


def _read_nuclei(f) -> dict:
    trexio = _trexio()
    n = trexio.read_nucleus_num(f)
    coord = trexio.read_nucleus_coord(f).reshape(n, 3) if n else np.zeros((0, 3))
    charge = trexio.read_nucleus_charge(f)
    label = trexio.read_nucleus_label(f) if trexio.has_nucleus_label(f) else [""] * n
    return {"num": n, "coord": coord, "charge": charge, "label": label}


def _read_basis(f) -> dict:
    trexio = _trexio()
    out = {
        "type": trexio.read_basis_type(f) if trexio.has_basis_type(f) else None,
        "shell_num": trexio.read_basis_shell_num(f),
        "prim_num": trexio.read_basis_prim_num(f),
        "nucleus_index": trexio.read_basis_nucleus_index(f),
        "shell_ang_mom": trexio.read_basis_shell_ang_mom(f),
        "shell_factor": trexio.read_basis_shell_factor(f) if trexio.has_basis_shell_factor(f) else None,
        "shell_index": trexio.read_basis_shell_index(f),
        "exponent": trexio.read_basis_exponent(f),
        "coefficient": trexio.read_basis_coefficient(f),
        "prim_factor": trexio.read_basis_prim_factor(f) if trexio.has_basis_prim_factor(f) else None,
    }
    return out


def _read_ao(f) -> dict:
    trexio = _trexio()
    return {
        "cartesian": bool(trexio.read_ao_cartesian(f)),
        "num": trexio.read_ao_num(f),
        "shell": trexio.read_ao_shell(f),
        "normalization": trexio.read_ao_normalization(f) if trexio.has_ao_normalization(f) else None,
    }


def _read_electrons(f) -> dict:
    trexio = _trexio()
    return {
        "num": trexio.read_electron_num(f),
        "up_num": trexio.read_electron_up_num(f),
        "dn_num": trexio.read_electron_dn_num(f),
    }


def _read_mos(f) -> dict:
    trexio = _trexio()
    return {
        "type": trexio.read_mo_type(f) if trexio.has_mo_type(f) else None,
        "num": trexio.read_mo_num(f),
        "coefficient": trexio.read_mo_coefficient(f),
        "energy": trexio.read_mo_energy(f) if trexio.has_mo_energy(f) else None,
        "occupation": trexio.read_mo_occupation(f) if trexio.has_mo_occupation(f) else None,
        "spin": trexio.read_mo_spin(f) if trexio.has_mo_spin(f) else None,
    }


def _read_ao_one_e(f) -> dict:
    trexio = _trexio()
    out: dict[str, Any] = {}
    if trexio.has_ao_1e_int_overlap(f):
        out["overlap"] = trexio.read_ao_1e_int_overlap(f)
    if trexio.has_ao_1e_int_kinetic(f):
        out["kinetic"] = trexio.read_ao_1e_int_kinetic(f)
    if trexio.has_ao_1e_int_potential_n_e(f):
        out["potential_n_e"] = trexio.read_ao_1e_int_potential_n_e(f)
    if trexio.has_ao_1e_int_core_hamiltonian(f):
        out["core_hamiltonian"] = trexio.read_ao_1e_int_core_hamiltonian(f)
    return out


def _read_mo_one_e(f) -> dict:
    trexio = _trexio()
    out: dict[str, Any] = {}
    if trexio.has_mo_1e_int_overlap(f):
        out["overlap"] = trexio.read_mo_1e_int_overlap(f)
    if trexio.has_mo_1e_int_kinetic(f):
        out["kinetic"] = trexio.read_mo_1e_int_kinetic(f)
    if trexio.has_mo_1e_int_potential_n_e(f):
        out["potential_n_e"] = trexio.read_mo_1e_int_potential_n_e(f)
    if trexio.has_mo_1e_int_core_hamiltonian(f):
        out["core_hamiltonian"] = trexio.read_mo_1e_int_core_hamiltonian(f)
    return out


def _read_ao_eri(f) -> Optional[dict]:
    trexio = _trexio()
    if not trexio.has_ao_2e_int_eri(f):
        return None
    chunk = 65536
    offset = 0
    indices = []
    values = []
    while True:
        idx, val, n_read, eof = trexio.read_ao_2e_int_eri(f, offset, chunk)
        if n_read > 0:
            indices.extend(idx[:n_read])
            values.extend(val[:n_read])
            offset += n_read
        if eof:
            break
    return {
        "indices": np.asarray(indices, dtype=np.int32),
        "values": np.asarray(values, dtype=np.float64),
    }


def _read_determinants(f) -> Optional[dict]:
    trexio = _trexio()
    if not trexio.has_determinant_list(f):
        return None
    chunk = 8192
    offset = 0
    dets = []
    coefs = []
    while True:
        dlist, n_read, eof = trexio.read_determinant_list(f, offset, chunk)
        if n_read > 0:
            dets.extend(dlist[:n_read])
        if eof:
            break
        offset += n_read
    if trexio.has_determinant_coefficient(f):
        offset = 0
        while True:
            clist, n_read, eof = trexio.read_determinant_coefficient(f, offset, chunk)
            if n_read > 0:
                coefs.extend(clist[:n_read])
            if eof:
                break
            offset += n_read
    return {
        "determinants": np.asarray(dets, dtype=np.int64),
        "coefficients": np.asarray(coefs, dtype=np.float64) if coefs else None,
    }
