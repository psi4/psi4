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
file format used by the TREX-CoE quantum-chemistry / QMC ecosystem.

This module is a thin wrapper around the official ``trexio`` Python package.
It supports four layers of data:

* **Geometry + basis** -- nuclei, AO shells, primitives, normalization.
* **Molecular orbitals** -- coefficients in the AO basis, orbital energies,
  occupations, spin labels.
* **Integrals** -- AO/MO one-electron (S, T, V, h) and two-electron (ERI),
  written in sparse 8-fold-permutation form.
* **Determinants** -- bitfield CI expansion (for DETCI wavefunctions).

Files follow the conventions of the TrexIO specification (``trex.org``), so
they interoperate with other TrexIO consumers (TurboRVB, QMC=Chem,
trexio_tools, ...), and are checked against it by trexio-validate, which
recomputes the integrals from the basis stored in the file:

* **Spherical AO ordering**: ``m = 0, +1, -1, +2, -2, ..., +L, -L`` (for
  ``p``: ``z, x, y``). This is Psi4's native order, so no permutation.
* **Solid-harmonic phase**: TrexIO fixes the real solid harmonics to have a
  positive leading coefficient (``S_1^{+1} = +x``, ``S_1^{-1} = +y``).
  Psi4's real solid harmonics already have this phase.
* **Cartesian AO ordering**: alphabetical, i.e. lex descending in
  ``(lx, ly, lz)`` (``xx, xy, xz, yy, yz, zz`` for ``d``). Also Psi4's order.
* **Primitive normalization**: ``basis_coefficient = original_coef``
  (input-file coefficient), ``basis_prim_factor = N_p`` (canonical Gaussian
  primitive normalization). The product reproduces Psi4's normalized
  contraction.
* **Cartesian per-AO factor**: ``ao_normalization[i] = sqrt((2L-1)!! /
  ((2lx-1)!!(2ly-1)!!(2lz-1)!!))`` so every Cartesian AO has unit self-overlap.
  Sphericals use ``ao_normalization = 1``.
* **Two-electron integrals**: physicists' notation ``<ij|kl> = (ik|jl)``.
  Psi4 computes chemists' ``(ij|kl)``; the indices are reordered on write.
* **Effective core potentials** (``ecp`` group): ``nucleus.charge`` is the
  effective charge ``Z - z_core`` (Psi4 already reduces its nuclear charges
  when it builds an ECP basis); ``ecp.power`` is the actual power of ``r``,
  i.e. Psi4's Gaussian-format ``nval - 2`` (libecpint subtracts the 2); and the
  local channel -- Psi4's, and libecpint's, highest angular momentum on the
  centre -- is stored under ``ang_mom = max_ang_mom_plus_1``. The core
  Hamiltonian includes the ECP term, which is also stored on its own as
  ``ao_1e_int.ecp`` / ``mo_1e_int.ecp``. This is the convention of the spec, and
  of pyscf-forge's TrexIO reader and writer.

Round-tripping inside Psi4 is exact.
"""

from __future__ import annotations

import os
from typing import Any, Optional

import numpy as np

from psi4 import core


def _double_factorial(n: int) -> int:
    """``n!!`` with the convention ``(-1)!! = 0!! = 1``."""
    if n <= 0:
        return 1
    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result


def _canonical_ao_normalization(basisset) -> np.ndarray:
    """Per-AO normalization factor ``ao.normalization``, in AO order.

    For sphericals: ``1``.
    For Cartesians: ``sqrt((2L-1)!! / ((2lx-1)!!(2ly-1)!!(2lz-1)!!))`` so that
    every Cartesian AO has unit self-overlap, matching the TrexIO spec.
    """
    nbf = basisset.nbf()
    norm_psi4_order = np.ones(nbf, dtype=np.float64)
    offset = 0
    for s in range(basisset.nshell()):
        sh = basisset.shell(s)
        L = sh.am
        nfn = sh.nfunction
        if not sh.is_pure():
            df_L = _double_factorial(2 * L - 1)
            # Psi4 Cartesian iter: start (a=L,b=0,c=0); next: (a,b-1,c+1) until c==L-a, else (a-1,L-a,0).
            a, b, c = L, 0, 0
            for i in range(nfn):
                df_a = _double_factorial(2 * a - 1)
                df_b = _double_factorial(2 * b - 1)
                df_c = _double_factorial(2 * c - 1)
                norm_psi4_order[offset + i] = np.sqrt(df_L / (df_a * df_b * df_c))
                if c < L - a:
                    b -= 1
                    c += 1
                else:
                    a -= 1
                    c = 0
                    b = L - a
        offset += nfn
    return norm_psi4_order


def _canonicalize_C(C_psi4: np.ndarray, ao_norm: np.ndarray) -> np.ndarray:
    """Rescale (nao, nmo) coefficients for TrexIO's per-AO normalization.

    ``MO_p = sum_i C_psi4[i,p] * AO_psi4[i] = sum_i C_trexio[i,p] * AO_trexio[i]``
    with ``AO_trexio[i] = ao_norm[i] * AO_psi4[i]``, so
    ``C_trexio[i,p] = C_psi4[i,p] / ao_norm[i]``. No permutation: Psi4's AO
    order is already TrexIO's.
    """
    return C_psi4 / ao_norm[:, None]


def _canonicalize_ao_matrix(M_psi4: np.ndarray, ao_norm: np.ndarray) -> np.ndarray:
    """Rescale an AO-basis (nao, nao) operator matrix for TrexIO's per-AO normalization.

    ``<AO_trexio[i]|O|AO_trexio[j]> = ao_norm[i]*ao_norm[j]*<AO_psi4[i]|O|AO_psi4[j]>``.
    """
    return M_psi4 * ao_norm[:, None] * ao_norm[None, :]


__all__ = [
    "TrexIOError",
    "trexio",
    "trexio_from_file",
    "trexio_to_wavefunction",
]


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


def trexio(
    wfn: "core.Wavefunction",
    filepath: str,
    *,
    overwrite: bool = False,
    save_ao_integrals: bool = False,
    save_mo_integrals: bool = False,
    save_eri: bool = False,
    save_mo_eri: bool = False,
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
        Write the AO ERI tensor in sparse 8-fold form.
    save_mo_eri
        Write the MO ERI tensor (transformed from AO) in sparse 8-fold form.
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
        # A molecule is not periodic. The spec does not require saying so, but
        # some readers (pyscf-forge's) read pbc.periodic unconditionally.
        trexio.write_pbc_periodic(f, 0)
        basisset = wfn.basisset()
        _write_basis(f, basisset)
        _write_ecp(f, wfn.molecule(), basisset)
        _write_ao(f, basisset)
        _write_electrons(f, wfn)
        _write_mos(f, wfn, basisset)

        if save_ao_integrals or save_mo_integrals or save_eri or save_mo_eri:
            mints = core.MintsHelper(basisset)
            if save_ao_integrals:
                _write_ao_one_e(f, mints, wfn)
            if save_mo_integrals:
                _write_mo_one_e(f, mints, wfn)
            if save_eri:
                _write_ao_eri(f, mints, basisset)
            if save_mo_eri:
                _write_mo_eri(f, mints, wfn)

        if determinants is not None or ci_vector is not None:
            if determinants is None or ci_vector is None:
                raise TrexIOError(
                    "Both 'determinants' and 'ci_vector' must be supplied to save a CI expansion."
                )
            _write_determinants(f, determinants, ci_vector)

    return filepath


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------


def trexio_from_file(filepath: str) -> dict:
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
        ecp = _read_ecp(f)
        if ecp is not None:
            out["ecp"] = ecp
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
        mo_2e = _read_mo_eri(f)
        if mo_2e is not None:
            out["mo_2e"] = mo_2e
        det = _read_determinants(f)
        if det is not None:
            out["determinant"] = det

    return out


def _basisset_from_trexio(data: dict, mol) -> "core.BasisSet":
    """Reconstruct a Psi4 BasisSet from the basis block in a TrexIO data dict.

    Uses :func:`psi4.core.BasisSet.construct_from_pydict`, which accepts a
    shell-by-shell representation with the *unnormalized* (input-file)
    contraction coefficients. TrexIO canonical stores exactly that in
    ``basis_coefficient``, so no conversion is needed.
    """
    basis = data["basis"]
    ao = data["ao"]
    nuc = data["nucleus"]

    natom = nuc["num"]
    shell_ang_mom = np.asarray(basis["shell_ang_mom"], dtype=np.int64)
    nucleus_index = np.asarray(basis["nucleus_index"], dtype=np.int64)
    shell_index = np.asarray(basis["shell_index"], dtype=np.int64)
    exponents = np.asarray(basis["exponent"], dtype=np.float64)
    coefficients = np.asarray(basis["coefficient"], dtype=np.float64)

    # Group primitives by shell.
    nshell = int(basis["shell_num"])
    shell_prims: list[list[tuple[float, float]]] = [[] for _ in range(nshell)]
    for p in range(int(basis["prim_num"])):
        s = int(shell_index[p])
        shell_prims[s].append((float(exponents[p]), float(coefficients[p])))

    shell_map = []
    for atom in range(natom):
        label = str(nuc["label"][atom]).strip() or "X"
        atom_entry: list = [label, ""]
        for s in range(nshell):
            if int(nucleus_index[s]) != atom:
                continue
            shell_entry: list = [int(shell_ang_mom[s])]
            for exp, coef in shell_prims[s]:
                shell_entry.append([exp, coef])
            atom_entry.append(shell_entry)
        shell_map.append(atom_entry)

    # "BASIS", not "ORBITAL": construct_from_pydict reduces the nuclear charges
    # by the ECP core electrons only for the key "BASIS" (_pybuild_basis maps
    # ORBITAL to BASIS before calling it). With "ORBITAL" an ECP file would
    # rebuild full nuclear charges *and* the ECP, counting the core twice.
    pybs = {
        "key": "BASIS",
        "name": "trexio",
        "blend": "trexio",
        "puream": 0 if ao["cartesian"] else 1,
        "shell_map": shell_map,
        "molecule": mol.to_dict(),
    }
    if "ecp" in data:
        pybs["ecp_shell_map"] = _ecp_shell_map(data)
    return core.BasisSet.construct_from_pydict(mol, pybs, -1)


def _ecp_shell_map(data: dict) -> list:
    """TrexIO's ecp group as construct_from_pydict's ``ecp_shell_map``.

    The inverse of _ecp_items: the channel stored under ang_mom ==
    max_ang_mom_plus_1 is the local one, which Psi4 marks with a negative angular
    momentum (ECPType1); and Psi4 wants the Gaussian-format radial power, the
    TrexIO power plus 2.
    """
    ecp = data["ecp"]
    nuc = data["nucleus"]
    natom = int(nuc["num"])
    channels: dict = {}
    for k in range(ecp["num"]):
        A, l = int(ecp["nucleus_index"][k]), int(ecp["ang_mom"][k])
        channels.setdefault((A, l), []).append(
            [float(ecp["exponent"][k]), float(ecp["coefficient"][k]), int(ecp["power"][k]) + 2])
    ecp_map = []
    for A in range(natom):
        label = str(nuc["label"][A]).strip() or "X"
        entry: list = [label, "", int(ecp["z_core"][A])]
        local = int(ecp["max_ang_mom_plus_1"][A])
        for (B, l), prims in sorted(channels.items(), key=lambda kv: kv[0]):
            if B != A:
                continue
            entry.append([-l if l == local else l] + prims)
        ecp_map.append(entry)
    return ecp_map


def trexio_to_wavefunction(
    filepath: str,
    *,
    reference: str = "rhf",
    basis_name: Optional[str] = None,
) -> "core.Wavefunction":
    """Reconstruct a Psi4 Wavefunction (geometry + basis + MOs) from a TrexIO file.

    The molecule is rebuilt from the nuclear coordinates and charges. The
    AO basis is rebuilt directly from the TrexIO ``basis`` block via
    :func:`psi4.core.BasisSet.construct_from_pydict`, so custom or in-file
    bases round-trip without needing the Psi4 basis library. If
    ``basis_name`` is supplied, the library basis is used instead (handy for
    interop with files produced by other codes).

    The reconstructed wavefunction is C1-symmetric.
    """
    data = trexio_from_file(filepath)

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

    if basis_name is not None:
        basisset = core.BasisSet.build(mol, "ORBITAL", basis_name, puream=-1)
    elif "basis" in data and "ao" in data:
        basisset = _basisset_from_trexio(data, mol)
        basis_name = "trexio"
    else:
        raise TrexIOError(
            "TrexIO file lacks a basis block; pass basis_name= explicitly."
        )
    nbf = basisset.nbf()
    if nbf != data["mo"]["num"]:
        raise TrexIOError(
            f"AO/MO count mismatch: basis {basis_name!r} gives nbf={nbf}, "
            f"but file has {data['mo']['num']} MOs."
        )

    # Coefficients in the file are canonical (nmo, nao_canonical). Convert back
    # to Psi4 AO ordering: AO_psi4[psi4_idx] = AO_can[canon_idx] / ao_norm[canon_idx].
    Cmat_canon = np.asarray(data["mo"]["coefficient"]).reshape(data["mo"]["num"], nbf).T
    ao_norm = _canonical_ao_normalization(basisset)
    Cmat = Cmat_canon * ao_norm[:, None]  # multiply because forward used divide
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
    if description:
        trexio.write_metadata_description(f, description)


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
        labels.append(mol.symbol(i).capitalize())
    trexio.write_nucleus_coord(f, coords)
    trexio.write_nucleus_charge(f, charges)
    trexio.write_nucleus_label(f, labels)
    trexio.write_nucleus_point_group(f, mol.point_group().symbol())
    trexio.write_nucleus_repulsion(f, mol.nuclear_repulsion_energy())


def _write_basis(f, basisset) -> None:
    """Write the basis block in TrexIO canonical convention.

    ``basis_coefficient`` holds the unnormalized (input-file-style)
    contraction coefficients. ``basis_prim_factor`` holds the canonical
    Gaussian primitive normalization ``N_p`` such that ``N_p * orig_coef``
    reproduces Psi4's internally normalized coefficient (i.e. integrals come
    out right when the consumer multiplies ``coefficient * prim_factor``).

    ``basis_shell_factor`` is ``1`` for Gaussian Psi4 bases.
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
            orig = sh.original_coef(p)
            full = sh.coef(p)
            exponents.append(sh.exp(p))
            coefficients.append(orig)
            # f_ks * gamma_ks reproduces Psi4's normalized coefficient.
            # When orig_coef is zero (rare; degenerate contraction), fall back to (1, full).
            prim_factor.append(full / orig if orig != 0.0 else 1.0)
            if orig == 0.0:
                coefficients[-1] = full
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


def _ecp_items(mol, basisset):
    """Psi4's ECP in TrexIO's ecp-group encoding, or None without an ECP.

    TrexIO writes V_ECP = V_local + sum_l dV_l P_l, each term beta r^n exp(-alpha r^2)
    with n the actual power of r (trex.org, "ecp group"), and it stores the local
    channel under ang_mom = max_ang_mom_plus_1. Psi4 keeps the Gaussian-format
    radial power, whose term is r^(nval-2): libecpint's GaussianECP subtracts the 2
    on construction. And the local channel is the shell of highest angular momentum
    on a centre -- the rule libecpint itself applies (L = max l) when Psi4 hands it
    the shells, so it is also what Psi4's integrals mean.

    Z_core is Psi4's own count per atom; nucleus.charge already holds Z_eff, since
    building an ECP orbital basis sets each nuclear charge to Z - ncore. Ghost
    atoms keep charge 0 and have no core removed.

    Returns (max_ang_mom_plus_1, z_core, ang_mom, nucleus_index, exponent,
    coefficient, power).
    """
    if not basisset.has_ECP():
        return None
    natom = mol.natom()
    lmax1 = np.zeros(natom, dtype=np.int32)
    z_core = np.zeros(natom, dtype=np.int32)
    ang, nuc, exps, coefs, powers = [], [], [], [], []
    for A in range(natom):
        if mol.Z(A) > 0:
            z_core[A] = basisset.n_ecp_core(mol.label(A))
        # ecp_shell_on_center returns an index into the ECP shell list, not a shell
        shells = [basisset.ecp_shell(basisset.ecp_shell_on_center(A, i))
                  for i in range(basisset.n_ecp_shell_on_center(A))]
        if not shells:
            continue
        lmax1[A] = max(sh.am for sh in shells)
        for sh in shells:
            for k in range(sh.nprimitive):
                ang.append(sh.am)
                nuc.append(A)
                exps.append(sh.exp(k))
                coefs.append(sh.coef(k))
                powers.append(sh.nval(k) - 2)
    return (lmax1, z_core, np.asarray(ang, dtype=np.int32), np.asarray(nuc, dtype=np.int32),
            np.asarray(exps), np.asarray(coefs), np.asarray(powers, dtype=np.int32))


def _write_ecp(f, mol, basisset) -> None:
    items = _ecp_items(mol, basisset)
    if items is None:
        return
    trexio = _trexio()
    lmax1, z_core, ang, nuc, exps, coefs, powers = items
    trexio.write_ecp_max_ang_mom_plus_1(f, lmax1)
    trexio.write_ecp_z_core(f, z_core)
    trexio.write_ecp_num(f, len(ang))
    trexio.write_ecp_ang_mom(f, ang)
    trexio.write_ecp_nucleus_index(f, nuc)
    trexio.write_ecp_exponent(f, exps)
    trexio.write_ecp_coefficient(f, coefs)
    trexio.write_ecp_power(f, powers)


def _write_ao(f, basisset) -> None:
    """Write AO block in canonical order with the canonical normalization factor."""
    trexio = _trexio()
    nbf = basisset.nbf()
    cartesian = 0 if basisset.has_puream() else 1
    ao_norm_canonical = _canonical_ao_normalization(basisset)

    # ao_shell in Psi4 order, then permute to canonical order.
    ao_shell_psi4 = np.zeros(nbf, dtype=np.int32)
    idx = 0
    for s in range(basisset.nshell()):
        sh = basisset.shell(s)
        for i in range(sh.nfunction):
            ao_shell_psi4[idx + i] = s
        idx += sh.nfunction

    trexio.write_ao_cartesian(f, cartesian)
    trexio.write_ao_num(f, nbf)
    trexio.write_ao_shell(f, ao_shell_psi4)
    trexio.write_ao_normalization(f, ao_norm_canonical)


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
    C, energies, occ, spin, mo_type = _file_mos(wfn)
    nbf = basisset.nbf()
    if C.shape[0] != nbf:
        raise TrexIOError(
            f"MO coefficient AO dimension {C.shape[0]} does not match basis nbf {nbf}."
        )
    coefs = _canonicalize_C(C, _canonical_ao_normalization(basisset)).T.copy()

    trexio.write_mo_type(f, mo_type)
    trexio.write_mo_num(f, coefs.shape[0])
    trexio.write_mo_coefficient(f, coefs)
    trexio.write_mo_energy(f, energies)
    trexio.write_mo_occupation(f, occ)
    trexio.write_mo_spin(f, spin)


def _write_ao_one_e(f, mints, wfn) -> None:
    trexio = _trexio()
    ao_norm = _canonical_ao_normalization(wfn.basisset())

    def _canon(M):
        return _canonicalize_ao_matrix(np.asarray(M), ao_norm)

    S = _canon(mints.ao_overlap())
    T = _canon(mints.ao_kinetic())
    V = _canon(mints.ao_potential())
    H = T + V
    # The core Hamiltonian has to be the one the SCF used, which for an ECP basis
    # is T + V + V_ECP (MintsHelper::so_potential adds so_ecp()).
    if wfn.basisset().has_ECP():
        U = _canon(mints.ao_ecp())
        trexio.write_ao_1e_int_ecp(f, U)
        H = H + U
    trexio.write_ao_1e_int_overlap(f, S)
    trexio.write_ao_1e_int_kinetic(f, T)
    trexio.write_ao_1e_int_potential_n_e(f, V)
    trexio.write_ao_1e_int_core_hamiltonian(f, H)


def _irrep_occupations(occpi, nmopi) -> np.ndarray:
    """0/1 occupations in _C1_MO's irrep-blocked column order.

    Within each irrep Psi4 orders MOs by energy, occupied first, so irrep ``h``
    contributes ``occpi[h]`` ones then ``nmopi[h] - occpi[h]`` zeros.
    """
    return np.concatenate([
        np.concatenate([np.ones(occpi[h]), np.zeros(nmopi[h] - occpi[h])])
        for h in range(nmopi.n())
    ])


def _file_mos(wfn):
    """Every MO-basis quantity in the order the file stores the MOs.

    Returns ``(C, energies, occupations, spin, mo_type)`` with ``C`` of shape
    (nao, mo.num) in Psi4's AO normalization. This is the only place MO order is
    decided: _write_mos and the MO integral writers all go through it, so the
    coefficients, energies, occupations and integrals cannot disagree.

    _C1_MO stacks MOs by irrep, not by energy. Occupations therefore have to be
    built per irrep and carried along when sorting -- marking the first nocc
    columns as occupied is wrong for any molecule with symmetry, and a check on
    the total electron count cannot see it.

    Restricted wavefunctions (RHF, ROHF) have one set, sorted by energy.
    Unrestricted ones store alpha then beta, each sorted by its own energies.
    """
    aotoso = wfn.aotoso()
    nmopi = wfn.nmopi()
    occ_a = _irrep_occupations(wfn.nalphapi(), nmopi)
    occ_b = _irrep_occupations(wfn.nbetapi(), nmopi)

    def _sorted(C, eps, occ):
        order = np.argsort(eps, kind="stable")
        return C[:, order], eps[order], occ[order]

    Ca, eps_a, occ_a = _sorted(_C1_MO(wfn.Ca(), aotoso), _flatten_vector(wfn.epsilon_a()), occ_a)
    if wfn.same_a_b_orbs():
        # One set: nbeta-irrep-occupied orbitals are doubly occupied, the rest of
        # the alpha-occupied ones singly (ROHF). Sorting used alpha energies, so
        # the beta occupations must follow the same permutation.
        order = np.argsort(_flatten_vector(wfn.epsilon_a()), kind="stable")
        occ = occ_a + occ_b[order]
        mo_type = "RHF" if wfn.nalpha() == wfn.nbeta() else "ROHF"
        return Ca, eps_a, occ, np.zeros(Ca.shape[1], dtype=np.int32), mo_type

    Cb, eps_b, occ_b = _sorted(_C1_MO(wfn.Cb(), aotoso), _flatten_vector(wfn.epsilon_b()), occ_b)
    spin = np.concatenate([np.zeros(Ca.shape[1], dtype=np.int32),
                           np.ones(Cb.shape[1], dtype=np.int32)])
    return (np.hstack([Ca, Cb]), np.concatenate([eps_a, eps_b]),
            np.concatenate([occ_a, occ_b]), spin, "UHF")


def _write_mo_one_e(f, mints, wfn) -> None:
    trexio = _trexio()
    # MO 1e ints are invariant under our AO rescaling (C absorbs 1/ao_norm), so
    # they are computed in Psi4's AO basis. They must be mo.num x mo.num: TrexIO
    # reads that many elements, so a smaller array would be read past its end.
    # Spin-forbidden (alpha-beta) elements of a spin-free operator are zero.
    C, _, _, spin, _ = _file_mos(wfn)
    same_spin = spin[:, None] == spin[None, :]
    S = np.asarray(mints.ao_overlap())
    T = np.asarray(mints.ao_kinetic())
    V = np.asarray(mints.ao_potential())

    def _mo(X):
        return np.where(same_spin, C.T @ X @ C, 0.0)

    H = T + V
    if wfn.basisset().has_ECP():
        U = np.asarray(mints.ao_ecp())
        trexio.write_mo_1e_int_ecp(f, _mo(U))
        H = H + U
    trexio.write_mo_1e_int_overlap(f, _mo(S))
    trexio.write_mo_1e_int_kinetic(f, _mo(T))
    trexio.write_mo_1e_int_potential_n_e(f, _mo(V))
    trexio.write_mo_1e_int_core_hamiltonian(f, _mo(H))


def _write_ao_eri(f, mints, basisset) -> None:
    """Write the AO ERI in canonical AO order, 8-fold permutational symmetry, sparse."""
    trexio = _trexio()
    ao_norm = _canonical_ao_normalization(basisset)
    eri = np.asarray(mints.ao_eri())
    nbf = eri.shape[0]
    # Rescale by the ao_norm tensor product (Cartesians only; 1 for sphericals).
    scale = ao_norm
    eri = eri * scale[:, None, None, None] * scale[None, :, None, None] \
              * scale[None, None, :, None] * scale[None, None, None, :]

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
                    # Psi4's (ij|kl) is TrexIO's physicists' <ik|jl>.
                    indices.append([i, k, j, l])
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


def _write_mo_eri(f, mints, wfn) -> None:
    """Write the MO ERI in 8-fold permutational symmetry as a sparse list.

    MO integrals are basis-orientation invariant (the per-AO canonical
    rescaling cancels in the C^T O C transform), so we work directly in
    Psi4-order MOs.
    """
    trexio = _trexio()
    # Over every MO in the file: for unrestricted wavefunctions that is alpha
    # then beta, and all four spin blocks are needed, not just alpha-alpha.
    C, _, _, spin, _ = _file_mos(wfn)
    mo_eri = np.asarray(mints.ao_eri())
    # Transform pqrs in AO → MO via einsum.
    mo_eri = np.einsum("pqrs,pi,qj,rk,sl->ijkl", mo_eri, C, C, C, C, optimize=True)
    nmo = mo_eri.shape[0]

    indices: list[list[int]] = []
    values: list[float] = []
    for i in range(nmo):
        for j in range(i + 1):
            ij = i * (i + 1) // 2 + j
            for k in range(nmo):
                for l in range(k + 1):
                    kl = k * (k + 1) // 2 + l
                    if kl > ij:
                        continue
                    # chemists' (ij|kl) needs spin(i) == spin(j), spin(k) == spin(l)
                    if spin[i] != spin[j] or spin[k] != spin[l]:
                        continue
                    v = mo_eri[i, j, k, l]
                    if abs(v) < 1e-15:
                        continue
                    # Psi4's (ij|kl) is TrexIO's physicists' <ik|jl>.
                    indices.append([i, k, j, l])
                    values.append(float(v))

    if not values:
        return
    chunk = 65536
    n = len(values)
    arr_idx = np.asarray(indices, dtype=np.int32)
    arr_val = np.asarray(values, dtype=np.float64)
    offset = 0
    while offset < n:
        end = min(offset + chunk, n)
        trexio.write_mo_2e_int_eri(
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
    out = {"num": n, "coord": coord, "charge": charge, "label": label}
    if trexio.has_nucleus_repulsion(f):
        out["repulsion"] = float(trexio.read_nucleus_repulsion(f))
    return out


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
    if trexio.has_ao_1e_int_ecp(f):
        out["ecp"] = trexio.read_ao_1e_int_ecp(f)
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
    if trexio.has_mo_1e_int_ecp(f):
        out["ecp"] = trexio.read_mo_1e_int_ecp(f)
    return out


def _read_ecp(f) -> Optional[dict]:
    trexio = _trexio()
    if not trexio.has_ecp_num(f):
        return None
    return {
        "max_ang_mom_plus_1": np.asarray(trexio.read_ecp_max_ang_mom_plus_1(f), dtype=np.int64),
        "z_core": np.asarray(trexio.read_ecp_z_core(f), dtype=np.int64),
        "num": int(trexio.read_ecp_num(f)),
        "ang_mom": np.asarray(trexio.read_ecp_ang_mom(f), dtype=np.int64),
        "nucleus_index": np.asarray(trexio.read_ecp_nucleus_index(f), dtype=np.int64),
        "exponent": np.asarray(trexio.read_ecp_exponent(f), dtype=np.float64),
        "coefficient": np.asarray(trexio.read_ecp_coefficient(f), dtype=np.float64),
        "power": np.asarray(trexio.read_ecp_power(f), dtype=np.int64),
    }


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


def _read_mo_eri(f) -> Optional[dict]:
    trexio = _trexio()
    if not trexio.has_mo_2e_int_eri(f):
        return None
    chunk = 65536
    offset = 0
    indices = []
    values = []
    while True:
        idx, val, n_read, eof = trexio.read_mo_2e_int_eri(f, offset, chunk)
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
