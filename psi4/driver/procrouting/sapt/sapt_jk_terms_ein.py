#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import time
from typing import List, Tuple

import numpy as np

from psi4 import core

from ...p4util import solvers
from .sapt_util import print_sapt_var
import einsums as ein


# Equations come from https://doi.org/10.1063/5.0090688
def localization(
    cache: dict,
    dimer_wfn: core.Wavefunction,
    do_print: bool = True,
) -> None:
    r"""Localize dimer occupied orbitals via Intrinsic Bond Orbitals (IBO).

    Performs IBO localization on the dimer occupied orbitals, separating
    them into frozen-core and active localized subsets. The localized
    orbitals, rotation matrices, and IAO projectors are stored in *cache*.

    Parameters
    ----------
    cache : dict
        SAPT data cache. Must already contain ``'Cocc'``, ``'eps_occ'``.
        Updated in-place with ``'Locc'``, ``'Qocc'``, ``'IAO'``,
        ``'Lfocc'``, and ``'Laocc'``.
    dimer_wfn : core.Wavefunction
        Dimer supermolecular wavefunction (provides basis set and molecule).
    do_print : bool, optional
        Whether to print status output, by default True.
    """
    core.print_out("\n  ==> Localizing Orbitals 1 <== \n\n")
    # Extract monomers to compute frozen core counts
    mol = dimer_wfn.molecule()
    molA = mol.extract_subsets([1], [])
    molB = mol.extract_subsets([2], [])
    nfocc0A = dimer_wfn.basisset().n_frozen_core(
        core.get_option("GLOBALS", "FREEZE_CORE"), molA
    )
    nfocc0B = dimer_wfn.basisset().n_frozen_core(
        core.get_option("GLOBALS", "FREEZE_CORE"), molB
    )
    nfocc_dimer = nfocc0A + nfocc0B

    N_eps_occ = cache["eps_occ"].dimpi()[0]
    Focc = core.Matrix("Focc", N_eps_occ, N_eps_occ)
    for i in range(N_eps_occ):
        Focc.np[i, i] = cache["eps_occ"].np[i]
    ranges = [0, nfocc_dimer, N_eps_occ]  # Separate frozen and active orbitals
    minao = core.BasisSet.build(
        dimer_wfn.molecule(), "BASIS", core.get_global_option("MINAO_BASIS")
    )
    dimer_wfn.set_basisset("MINAO", minao)
    # pybind11 IBO location: ./psi4/src/export_wavefunction.cc
    IBO_loc = core.IBOLocalizer2(
        dimer_wfn.basisset(),
        dimer_wfn.get_basisset("MINAO"),
        cache["Cocc"],
    )
    IBO_loc.print_header()
    ret = IBO_loc.localize(
        cache["Cocc"],
        Focc,
        ranges,
    )
    cache["Locc"] = ret["L"]
    cache["Qocc"] = ret["Q"]
    cache["IAO"] = ret["A"]

    # Extract frozen and active localized orbitals separately
    nn = cache["Cocc"].shape[0]  # number of AO basis functions
    nf = nfocc_dimer
    na = N_eps_occ - nfocc_dimer  # number of active occupied orbitals

    if nf > 0:
        # Store frozen core localized orbitals
        Lfocc = core.Matrix("Lfocc", nn, nf)
        Lfocc.np[:, :] = ret["L"].np[:, :nf]
        cache["Lfocc"] = Lfocc

    # Store active occupied localized orbitals
    Laocc = core.Matrix("Laocc", nn, na)
    Laocc.np[:, :] = ret["L"].np[:, nf:]
    cache["Laocc"] = Laocc
    return


def flocalization(
    cache: dict,
    dimer_wfn: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    do_print: bool = True,
) -> None:
    r"""Localize monomer occupied orbitals separately for F-SAPT partitioning.

    Performs IBO localization independently on monomer A and monomer B
    occupied orbitals. Separates frozen-core and active localized orbitals
    and stores them in *cache*. Handles link-orbital assignments for
    three-body (A-C-B) fragmentation schemes even though I-SAPT is not
    currently implemented for this module.

    Parameters
    ----------
    cache : dict
        SAPT data cache. Must already contain monomer orbital coefficients
        (``'Cocc_A'``, ``'Cocc_B'``) and orbital energies (``'eps_occ_A'``,
        ``'eps_occ_B'``). Updated in-place with ``'Locc_A'``, ``'Locc_B'``,
        ``'Uocc_A'``, ``'Uocc_B'``, ``'Qocc0A'``, ``'Qocc0B'``,
        ``'Lfocc0A'``, ``'Laocc0A'``, ``'Lfocc0B'``, ``'Laocc0B'``,
        ``'Caocc0A'``, and ``'Caocc0B'``.
    dimer_wfn : core.Wavefunction
        Dimer supermolecular wavefunction (provides basis set and molecule).
    do_print : bool, optional
        Whether to print status output, by default True.
    """
    link_assignment = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT").upper()
    core.print_out("  ==> F-SAPT Localization (IBO) <==\n\n")
    core.print_out("  ==> Local orbitals for Monomer A <==\n\n")
    mol = dimer_wfn.molecule()
    molA = mol.extract_subsets([1], [])
    molB = mol.extract_subsets([2], [])
    nfocc0A = dimer_wfn.basisset().n_frozen_core(
        core.get_option("GLOBALS", "FREEZE_CORE"), molA
    )
    nfocc0B = dimer_wfn.basisset().n_frozen_core(
        core.get_option("GLOBALS", "FREEZE_CORE"), molB
    )
    nn = cache["Cocc_A"].shape[0]
    nf = nfocc0A
    nm = cache["Cocc_A"].shape[1]  # total occupied orbitals (frozen + active)
    na = nm - nf  # active occupied orbitals only
    ranges = [0, nf, nm]
    N = cache["eps_occ_A"].shape[0]
    Focc = core.Matrix("Focc", N, N)
    for i in range(N):
        Focc.np[i, i] = cache["eps_occ_A"].np[i]
    IBO_loc = core.IBOLocalizer2(
        dimer_wfn.basisset(),
        dimer_wfn.get_basisset("MINAO"),
        core.Matrix.from_array(cache["Cocc_A"]),
    )
    IBO_loc.print_header()
    ret = IBO_loc.localize(
        core.Matrix.from_array(cache["Cocc_A"]),
        Focc,
        ranges,
    )

    Locc_A = ret["L"]
    Uocc_A = ret["U"]
    Qocc0A = ret["Q"]

    cache["Locc_A"] = Locc_A
    cache["Uocc_A"] = Uocc_A
    cache["Qocc0A"] = Qocc0A

    Lfocc0A = core.Matrix("Lfocc0A", nn, nf)
    Laocc0A = core.Matrix("Laocc0A", nn, na)
    Ufocc0A = core.Matrix("Ufocc0A", nf, nf)
    Uaocc0A = core.Matrix("Uaocc0A", na, na)

    Lfocc0A.np[:, :] = Locc_A.np[:, :nf]
    Laocc0A.np[:, :] = Locc_A.np[:, nf : nf + na]
    Ufocc0A.np[:, :] = Uocc_A.np[:nf, :nf]
    Uaocc0A.np[:, :] = Uocc_A.np[nf : nf + na, nf : nf + na]

    cache["Lfocc0A"] = Lfocc0A
    cache["Laocc0A"] = Laocc0A
    cache["Ufocc0A"] = Ufocc0A
    cache["Uaocc0A"] = Uaocc0A
    # Store active occupied orbitals for dispersion (Caocc0A = Cocc_A[:, nf:])
    Caocc0A = core.Matrix("Caocc0A", nn, na)
    Caocc0A.np[:, :] = cache["Cocc_A"].np[:, nf : nf + na]
    cache["Caocc0A"] = Caocc0A

    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        Locc_A = core.Matrix("Locc_A", nn, nm + 1)
        Locc_A.np[:, :nm] = Locc_A.np[:, :]
        Locc_A.np[:, nm] = cache["thislinkA"].np[:, 0]
        cache["Locc_A"] = Locc_A
    else:
        cache["Locc_A"] = Locc_A

    core.print_out("  ==> Local orbitals for Monomer B <==\n\n")

    nn = cache["Cocc_B"].shape[0]
    nf = nfocc0B
    nm = cache["Cocc_B"].shape[1]  # total occupied orbitals (frozen + active)
    na = nm - nf  # active occupied orbitals only
    ranges = [0, nf, nm]

    N = cache["eps_occ_B"].shape[0]
    Focc = core.Matrix("Focc", N, N)
    for i in range(N):
        Focc.np[i, i] = cache["eps_occ_B"].np[i]

    IBO_loc = core.IBOLocalizer2(
        dimer_wfn.basisset(),
        dimer_wfn.get_basisset("MINAO"),
        core.Matrix.from_array(cache["Cocc_B"]),
    )
    IBO_loc.print_header()
    ret = IBO_loc.localize(
        core.Matrix.from_array(cache["Cocc_B"]),
        Focc,
        ranges,
    )

    Locc_B = ret["L"]
    Uocc_B = ret["U"]
    Qocc0B = ret["Q"]

    cache["Locc_B"] = Locc_B
    cache["Uocc_B"] = Uocc_B
    cache["Qocc0B"] = Qocc0B

    Lfocc0B = core.Matrix("Lfocc0B", nn, nf)
    Laocc0B = core.Matrix("Laocc0B", nn, na)
    Ufocc0B = core.Matrix("Ufocc0B", nf, nf)
    Uaocc0B = core.Matrix("Uaocc0B", na, na)

    Lfocc0B.np[:, :] = Locc_B.np[:, :nf]
    Laocc0B.np[:, :] = Locc_B.np[:, nf : nf + na]
    Ufocc0B.np[:, :] = Uocc_B.np[:nf, :nf]
    Uaocc0B.np[:, :] = Uocc_B.np[nf : nf + na, nf : nf + na]

    cache["Lfocc0B"] = Lfocc0B
    cache["Laocc0B"] = Laocc0B
    cache["Ufocc0B"] = Ufocc0B
    cache["Uaocc0B"] = Uaocc0B
    # Store active occupied orbitals for dispersion (Caocc0B = Cocc_B[:, nf:])
    Caocc0B = core.Matrix("Caocc0B", nn, na)
    Caocc0B.np[:, :] = cache["Cocc_B"].np[:, nf : nf + na]
    cache["Caocc0B"] = Caocc0B

    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        Locc_B = core.Matrix("Locc_B", nn, nm + 1)
        Locc_B.np[:, :nm] = Locc_B.np[:, :]
        Locc_B.np[:, nm] = cache["thislinkB"].np[:, 0]
        cache["Locc_B"] = Locc_B
    else:
        cache["Locc_B"] = Locc_B


def partition(
    cache: dict,
    dimer_wfn: core.Wavefunction,
    do_print: bool = True,
) -> None:
    r"""Partition localized orbitals into monomer A, monomer B, and link fragments.

    Uses IBO charges to assign each localized orbital to a fragment (A, B,
    or linker C). Handles automatic and manual link-orbital selection and
    various link-assignment schemes (SAO, SIAO, etc.).

    Parameters
    ----------
    cache : dict
        SAPT data cache. Must contain ``'Locc'`` and ``'Qocc'`` from a prior
        :func:`localization` call. Updated in-place with fragment-partitioned
        orbitals and assignment vectors.
    dimer_wfn : core.Wavefunction
        Dimer supermolecular wavefunction (provides molecule and basis set).
    do_print : bool, optional
        Whether to print status output, by default True.
    """
    core.print_out("\n  ==> Partitioning <== \n\n")
    # Sizing
    mol = dimer_wfn.molecule()
    natoms = mol.natom()
    n_Locc = cache["Locc"].shape[1]

    # Monomer Atoms
    fragments = mol.get_fragments()
    indA = np.arange(*fragments[0], dtype=int)
    indB = np.arange(*fragments[1], dtype=int)
    indC = (
        np.arange(*fragments[2], dtype=int)
        if len(fragments) == 3
        else np.array([], dtype=int)
    )
    cache["FRAG"] = core.Vector(natoms)
    frag = cache["FRAG"].np
    frag[:] = 0.0
    frag[indA] = 1.0
    frag[indB] = 2.0
    if indC.size:
        frag[indC] = 3.0
    core.print_out("Fragment lookup table:\n")
    cache["FRAG"].print_out()
    core.print_out("   => Atomic Partitioning <= \n\n")
    core.print_out(f"    Monomer A: {len(indA)} atoms\n")
    core.print_out(f"    Monomer B: {len(indB)} atoms\n")
    core.print_out(f"    Monomer C: {len(indC)} atoms\n\n")
    np.set_printoptions(precision=14, suppress=True)

    # Fragment Orbital Charges
    Locc = cache["Locc"].np  # (n_ao x n_occ)
    Qocc = cache["Qocc"].np  # (n_atom x n_occ) orbital populations per atom

    n_ao, n_occ = Locc.shape

    QF = core.Matrix(3, n_Locc).np
    QF.fill(0.0)
    QF[0, :] = Qocc[indA, :].sum(axis=0)
    QF[1, :] = Qocc[indB, :].sum(axis=0)
    if indC.size:
        QF[2, :] = Qocc[indC, :].sum(axis=0)

    # --- link identification ---
    link_orbs: List[int] = []
    link_atoms: List[Tuple[int, int]] = []
    link_types: List[str] = []

    def top_two_atoms_for_orb(a: int) -> tuple[int, int]:
        A_sorted = np.argsort(Qocc[:, a])[::-1]
        A1, A2 = int(A_sorted[0]), int(A_sorted[1])
        return (A1, A2) if A1 < A2 else (A2, A1)

    link_sel = core.get_option("FISAPT", "FISAPT_LINK_SELECTION").upper()
    if link_sel == "AUTOMATIC":
        delta = float(core.get_option("FISAPT", "FISAPT_CHARGE_COMPLETENESS"))
        for a in range(n_occ):
            if np.any(QF[:, a] > delta):
                continue
            if QF[0, a] + QF[2, a] > delta:
                link_orbs.append(a)
                link_types.append("AC")
            elif QF[1, a] + QF[2, a] > delta:
                link_orbs.append(a)
                link_types.append("BC")
            elif QF[0, a] + QF[1, a] > delta:
                link_orbs.append(a)
                link_types.append("AB")
            else:
                raise ValueError(
                    "FISAPT: 3c-2e style bond encountered (no single/pair exceeds delta)."
                )
        for a in link_orbs:
            link_atoms.append(top_two_atoms_for_orb(a))
    elif link_sel == "MANUAL":
        if not core.get_option("FISAPT", "FISAPT_MANUAL_LINKS"):
            raise ValueError(
                "FISAPT: MANUAL selection requires manual_links (0-based atom pairs)."
            )
        S = set(indA.tolist())
        T = set(indB.tolist())
        U = set(indC.tolist())
        for A1, A2 in core.get_option("FISAPT", "FISAPT_MANUAL_LINKS"):
            prod = Qocc[A1, :] * Qocc[A2, :]
            a = int(np.argmax(prod))
            link_orbs.append(a)
            A1_, A2_ = (A1, A2) if A1 < A2 else (A2, A1)
            link_atoms.append((A1_, A2_))
            if (A1_ in S) and (A2_ in U):
                link_types.append("AC")
            elif (A1_ in T) and (A2_ in U):
                link_types.append("BC")
            elif (A1_ in S) and (A2_ in T):
                link_types.append("AB")
            else:
                raise ValueError("FISAPT: manual pair is not AB, AC, or BC.")
    else:
        raise ValueError("FISAPT: Unrecognized FISAPT_LINK_SELECTION.")
    link_orbs = np.array(link_orbs, dtype=int)

    # --- Z per fragment originals ---
    ZA = core.Vector(natoms)
    ZB = core.Vector(natoms)
    ZC = core.Vector(natoms)
    ZA.np[:] = 0.0
    ZB.np[:] = 0.0
    ZC.np[:] = 0.0

    Z_all = np.array([mol.Z(i) for i in range(natoms)], dtype=float)
    ZA.np[indA] = Z_all[indA]
    ZB.np[indB] = Z_all[indB]
    if indC.size:
        ZC.np[indC] = Z_all[indC]

    cache["ZA"] = ZA
    cache["ZB"] = ZB
    cache["ZC"] = ZC
    cache["ZA_orig"] = core.Vector.from_array(ZA.np.copy())
    cache["ZB_orig"] = core.Vector.from_array(ZB.np.copy())
    cache["ZC_orig"] = core.Vector.from_array(ZC.np.copy())

    # --- link assignment (C vs AB vs SAO*/SIAO*) ---
    orbsA: List[int] = []
    orbsB: List[int] = []
    orbsC: List[int] = []
    orbsL: List[int] = []
    typesL: List[str] = []

    la = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT").upper()
    valid = {"AB", "C", "SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"}
    if la not in valid:
        raise ValueError("FISAPT: FISAPT_LINK_ASSIGNMENT not recognized.")

    # --- link assignment (C vs AB vs SAO*/SIAO*) ---
    orbsA: List[int] = []
    orbsB: List[int] = []
    orbsC: List[int] = []
    orbsL: List[int] = []
    typesL: List[str] = []

    la = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT").upper()
    valid = {"AB", "C", "SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"}
    if la not in valid:
        raise ValueError("FISAPT: FISAPT_LINK_ASSIGNMENT not recognized.")

    if la in {"C", "SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"}:
        for a, (A1, A2), t in zip(link_orbs, link_atoms, link_types):
            typesL.append(t)
            if t == "AC":
                ZA.np[A1] -= 1.0
                ZC.np[A1] += 1.0
                orbsC.append(a)
                orbsL.append(a)
            elif t == "BC":
                ZB.np[A1] -= 1.0
                ZC.np[A1] += 1.0
                orbsC.append(a)
                orbsL.append(a)
            elif t == "AB":
                ZA.np[A1] -= 1.0
                ZC.np[A1] += 1.0
                ZB.np[A2] -= 1.0
                ZC.np[A2] += 1.0
                orbsC.append(a)
                orbsL.append(a)
    elif la == "AB":
        for a, (A1, A2), t in zip(link_orbs, link_atoms, link_types):
            if t == "AC":
                ZA.np[A1] += 1.0
                ZC.np[A1] -= 1.0
                orbsA.append(a)
            elif t == "BC":
                ZB.np[A1] += 1.0
                ZC.np[A1] -= 1.0
                orbsB.append(a)
            elif t == "AB":
                raise ValueError(
                    "FISAPT: AB link requires LINK_ASSIGNMENT C in this scheme."
                )

    # --- electron counts per fragment; enforce closed-shell ---
    fragment_charges = mol.get_fragment_charges()
    qA, qB = int(fragment_charges[0]), int(fragment_charges[1])
    qC = int(fragment_charges[2]) if len(fragment_charges) == 3 else 0

    def i_round(x: float) -> int:
        # protect against bankers rounding
        return int(np.floor(x + 0.5))

    ZA2 = i_round(float(ZA.np.sum()))
    ZB2 = i_round(float(ZB.np.sum()))
    ZC2 = i_round(float(ZC.np.sum()))
    EA2, EB2, EC2 = ZA2 - qA, ZB2 - qB, ZC2 - qC

    if EA2 % 2 or EB2 % 2 or EC2 % 2:
        raise ValueError(
            "FISAPT: fragment charge incompatible with singlet (odd electron count)."
        )

    NA2, NB2, NC2 = EA2 // 2, EB2 // 2, EC2 // 2
    if (NA2 + NB2 + NC2) != n_occ:
        raise ValueError(
            "FISAPT: sum of fragment electrons incompatible with total electrons."
        )

    RA2 = NA2 - len(orbsA)
    RB2 = NB2 - len(orbsB)
    RC2 = NC2 - len(orbsC)

    # --- greedy fill using QF weights (C then A then B), excluding taken orbs ---
    taken = set(orbsA) | set(orbsB) | set(orbsC)

    def take_top(weights: np.ndarray, k: int, taken_set: set) -> list[int]:
        if k <= 0:
            return []
        order = np.argsort(weights)[::-1]
        picked: List[int] = []
        for a in order:
            if int(a) in taken_set:
                continue
            picked.append(int(a))
            taken_set.add(int(a))
            if len(picked) == k:
                break
        return picked

    orbsC += take_top(QF[2, :], RC2, taken)
    orbsA += take_top(QF[0, :], RA2, taken)
    orbsB += take_top(QF[1, :], RB2, taken)

    # --- sort & link ordering swap like C++ ---
    orbsA = np.array(sorted(set(orbsA)), dtype=int)
    orbsB = np.array(sorted(set(orbsB)), dtype=int)
    orbsC = np.array(sorted(set(orbsC)), dtype=int)
    orbsL = np.array(orbsL, dtype=int)
    if orbsL.size > 1 and orbsL[0] > orbsL[1]:
        orbsL[[0, 1]] = orbsL[[1, 0]]
        typesL[0], typesL[1] = typesL[1], typesL[0]

    # --- build LoccA/B/C/L as psi4 Matrices (column extracts) ---
    def cols(M_np: np.ndarray, idx: np.ndarray) -> core.Matrix:
        if idx.size == 0:
            return np.zeros(M_np.shape[0])
        return core.Matrix.from_array(M_np[:, idx])

    def extract_columns(cols, A: core.Matrix) -> core.Matrix:
        cols = np.asarray(cols, dtype=int)
        if cols.size == 0:
            return None
        A2 = A[:, cols]
        return core.Matrix.from_array(A2)

    cache["LoccA"] = extract_columns(orbsA, Locc)
    cache["LoccB"] = extract_columns(orbsB, Locc)
    cache["LoccC"] = extract_columns(orbsC, Locc)
    cache["LoccL"] = extract_columns(orbsL, Locc)

    cache["QF"] = QF
    # --- summary numbers  ---
    ZA_int, ZB_int, ZC_int = (
        i_round(ZA.np.sum()),
        i_round(ZB.np.sum()),
        i_round(ZC.np.sum()),
    )
    YA, YB, YC = int(2 * orbsA.size), int(2 * orbsB.size), int(2 * orbsC.size)

    core.print_out("   => Partition Summary <= \n\n")
    core.print_out(
        f"    Monomer A: {ZA_int - YA:2d} charge, {ZA_int:3d} protons, {
            YA:3d} electrons, {len(orbsA):3d} docc\n"
    )
    core.print_out(
        f"    Monomer B: {ZB_int - YB:2d} charge, {ZB_int:3d} protons, {
            YB:3d} electrons, {len(orbsB):3d} docc\n"
    )
    core.print_out(
        f"    Monomer C: {ZC_int - YC:2d} charge, {ZC_int:3d} protons, {
            YC:3d} electrons, {len(orbsC):3d} docc\n"
    )
    return cache


def build_sapt_jk_cache(
    wfn_dimer: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    jk: core.JK,
    do_print: bool = True,
    external_potentials: dict = None,
) -> dict:
    r"""Construct the dimer-centered basis set (DCBS) cache of integrals and
    matrices for SAPT(DFT).

    Builds all one- and two-electron integrals needed for the electrostatics,
    exchange, and induction components of SAPT(DFT). Density matrices (Eq. 5,
    7), Coulomb/exchange matrices, nuclear potentials, and overlap integrals
    are computed and stored in the returned dictionary.

    .. math::

        \mathbf{P}^X = \mathbf{C}^{X,\text{occ}}(\mathbf{C}^{X,\text{occ}})^\dagger
        \quad (\text{Eq. 5})

    .. math::

        \mathbf{P}^{X,\text{vir}} = \mathbf{C}^{X,\text{vir}}(\mathbf{C}^{X,\text{vir}})^\dagger
        \quad (\text{Eq. 7})

    Parameters
    ----------
    wfn_dimer : core.Wavefunction
        Dimer supermolecular wavefunction.
    wfn_A : core.Wavefunction
        Monomer A wavefunction in the dimer-centered basis set (DCBS).
    wfn_B : core.Wavefunction
        Monomer B wavefunction in the dimer-centered basis set (DCBS).
    jk : core.JK
        Psi4 JK integral engine for computing Coulomb and exchange matrices.
    do_print : bool, optional
        Whether to print status output, by default True.
    external_potentials : dict, optional
        Dictionary of external potential objects keyed by ``'A'`` and ``'B'``.

    Returns
    -------
    dict
        Cache dictionary containing orbital coefficients, density matrices,
        Coulomb (``'J_A'``, ``'J_B'``, ``'J_O'``), exchange (``'K_A'``, ``'K_B'``,
        ``'K_O'``), nuclear potential (``'V_A'``, ``'V_B'``), overlap (``'S'``),
        orbital energies, and nuclear repulsion interaction energy.
    """
    core.print_out("\n  ==> Preparing SAPT Data Cache <== \n\n")
    jk.print_header()

    cache = {}
    cache["wfn_A"] = wfn_A
    cache["wfn_B"] = wfn_B

    # NOTE: scf_A from FISAPT0 and SAPT(DFT) wfn_A have slightly different
    # coefficients, so numerical exactness is not achieved everywhere, but the
    # pytests for final values are still quite robust

    # First grab the orbitals as psi4.core.Matrix objects
    cache["Cocc_A"] = wfn_A.Ca_subset("AO", "OCC")
    cache["Cocc_A"].name = "Cocc_A"
    cache["Cvir_A"] = wfn_A.Ca_subset("AO", "VIR")
    cache["Cvir_A"].name = "Cvir_A"

    cache["Cocc_B"] = wfn_B.Ca_subset("AO", "OCC")
    cache["Cocc_B"].name = "Cocc_B"
    cache["Cvir_B"] = wfn_B.Ca_subset("AO", "VIR")
    cache["Cvir_B"].name = "Cvir_B"

    cache["eps_occ_A"] = wfn_A.epsilon_a_subset("AO", "OCC")
    cache["eps_vir_A"] = wfn_A.epsilon_a_subset("AO", "VIR")
    cache["eps_occ_B"] = wfn_B.epsilon_a_subset("AO", "OCC")
    cache["eps_vir_B"] = wfn_B.epsilon_a_subset("AO", "VIR")

    # localization
    if core.get_option("SAPT", "SAPT_DFT_DO_FSAPT"):
        cache["Cfocc"] = wfn_dimer.Ca_subset("AO", "FROZEN_OCC")
        cache["eps_all"] = wfn_dimer.epsilon_a_subset("AO", "ALL")

        cache["Call"] = wfn_dimer.Ca_subset("AO", "ALL")
        cache["Cocc"] = wfn_dimer.Ca_subset("AO", "OCC")
        cache["Cvir"] = wfn_dimer.Ca_subset("AO", "VIR")

        cache["eps_occ"] = wfn_dimer.epsilon_a_subset("AO", "OCC")
        cache["eps_vir"] = wfn_dimer.epsilon_a_subset("AO", "VIR")

        cache["Caocc"] = wfn_dimer.Ca_subset("AO", "ACTIVE_OCC")
        cache["Cavir"] = wfn_dimer.Ca_subset("AO", "ACTIVE_VIR")
        cache["Cfvir"] = wfn_dimer.Ca_subset("AO", "FROZEN_VIR")

        cache["eps_focc"] = wfn_dimer.epsilon_a_subset("AO", "FROZEN_OCC")
        cache["eps_aocc"] = wfn_dimer.epsilon_a_subset("AO", "ACTIVE_OCC")
        cache["eps_avir"] = wfn_dimer.epsilon_a_subset("AO", "ACTIVE_VIR")
        cache["eps_fvir"] = wfn_dimer.epsilon_a_subset("AO", "FROZEN_VIR")

    # Build the densities as HF takes an extra "step", Eq. 5
    cache["D_A"] = chain_gemm_einsums([cache["Cocc_A"], cache["Cocc_A"]], ["N", "T"])
    cache["D_B"] = chain_gemm_einsums([cache["Cocc_B"], cache["Cocc_B"]], ["N", "T"])
    # Eq. 7
    cache["P_A"] = chain_gemm_einsums([cache["Cvir_A"], cache["Cvir_A"]], ["N", "T"])
    cache["P_B"] = chain_gemm_einsums([cache["Cvir_B"], cache["Cvir_B"]], ["N", "T"])

    # Potential ints - store as psi4.core.Matrix
    mints = core.MintsHelper(wfn_A.basisset())
    cache["V_A"] = mints.ao_potential()
    mints = core.MintsHelper(wfn_B.basisset())
    cache["V_B"] = mints.ao_potential()

    # External Potentials need to add to V_A and V_B
    if external_potentials:
        if external_potentials.get("A") is not None:
            ext_A = wfn_A.external_pot().computePotentialMatrix(wfn_A.basisset())
            cache["V_A"].add(ext_A)
        if external_potentials.get("B") is not None:
            ext_B = wfn_B.external_pot().computePotentialMatrix(wfn_B.basisset())
            cache["V_B"].add(ext_B)

    # Anything else we might need
    # S corresponds to the overlap matrix, S^{AO}
    cache["S"] = wfn_A.S().clone()
    cache["S"].name = "S"

    # J and K matrices
    jk.C_clear()

    # Normal J/K for Monomer A
    jk.C_left_add(wfn_A.Ca_subset("SO", "OCC"))
    jk.C_right_add(wfn_A.Ca_subset("SO", "OCC"))

    # Normal J/K for Monomer B
    jk.C_left_add(wfn_B.Ca_subset("SO", "OCC"))
    jk.C_right_add(wfn_B.Ca_subset("SO", "OCC"))

    DB_S_CA = chain_gemm_einsums([cache["D_B"], cache["S"], cache["Cocc_A"]])
    jk.C_left_add(DB_S_CA)
    jk.C_right_add(cache["Cocc_A"])

    jk.compute()

    # Clone them as the JK object will overwrite. Store as psi4.core.Matrix
    cache["J_A"] = jk.J()[0].clone()
    cache["K_A"] = jk.K()[0].clone()
    cache["J_B"] = jk.J()[1].clone()
    cache["K_B"] = jk.K()[1].clone()
    cache["J_O"] = jk.J()[2].clone()
    # K_O needs transpose
    K_O = jk.K()[2].clone().transpose()
    cache["K_O"] = core.Matrix.from_array(K_O.np)
    cache["K_O"].name = "K_O"

    monA_nr = wfn_A.molecule().nuclear_repulsion_energy()
    monB_nr = wfn_B.molecule().nuclear_repulsion_energy()
    dimer_nr = wfn_A.molecule().extract_subsets([1, 2]).nuclear_repulsion_energy()

    cache["extern_extern_IE"] = 0.0
    if external_potentials:
        dimer_nr += wfn_dimer.external_pot().computeNuclearEnergy(wfn_dimer.molecule())
        if external_potentials.get("A") is not None:
            monA_nr += wfn_A.external_pot().computeNuclearEnergy(wfn_A.molecule())
        if external_potentials.get("B") is not None:
            monB_nr += wfn_B.external_pot().computeNuclearEnergy(wfn_B.molecule())
        if (
            external_potentials.get("A") is not None
            and external_potentials.get("B") is not None
        ):
            cache["extern_extern_IE"] = (
                wfn_A.external_pot().computeExternExternInteraction(
                    wfn_B.external_pot()
                )
            )

    cache["nuclear_repulsion_energy"] = dimer_nr - monA_nr - monB_nr
    return cache


def electrostatics(cache: dict, do_print: bool = True) -> tuple[dict, float]:
    r"""Compute the first-order electrostatic interaction energy :math:`E^{(1)}_{\text{elst}}`.

    Evaluates the Coulombic interaction between unperturbed monomer charge
    distributions (Eq. 4 of Xie et al.):

    .. math::

        E^{(1)}_{\text{elst}} = 2\mathbf{P}^A \cdot \mathbf{V}^B
        + 2\mathbf{P}^B \cdot \mathbf{V}^A
        + 4\mathbf{P}^B \cdot \mathbf{J}^A + V_{\text{nuc}}

    Parameters
    ----------
    cache : dict
        SAPT data cache from :func:`build_sapt_jk_cache`.
    do_print : bool, optional
        Whether to print the result, by default True.

    Returns
    -------
    tuple[dict, float]
        A dictionary ``{'Elst10,r': float}`` and the extern-extern
        interaction energy (zero if no external potentials).
    """
    if do_print:
        core.print_out("\n  ==> E10 Electrostatics <== \n\n")

    # Eq. 4
    Elst10 = 2.0 * ein.core.dot(cache["D_A"].np, cache["V_B"].np)
    Elst10 += 2.0 * ein.core.dot(cache["D_B"].np, cache["V_A"].np)
    Elst10 += 4.0 * ein.core.dot(cache["D_B"].np, cache["J_A"].np)
    Elst10 += cache["nuclear_repulsion_energy"]

    if do_print:
        core.print_out(print_sapt_var("Elst10,r ", Elst10, short=True))
        core.print_out("\n")

    # External Potentials interacting with each other (V_A_ext, V_B_ext)
    extern_extern_ie = 0
    if cache.get("extern_extern_IE"):
        extern_extern_ie = cache["extern_extern_IE"]
        core.print_out(print_sapt_var("Extern-Extern ", extern_extern_ie, short=True))
        core.print_out("\n")

    return {"Elst10,r": Elst10}, extern_extern_ie


def felst(
    cache: dict,
    sapt_elst: dict,
    dimer_wfn: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    jk: core.JK,
    do_print: bool = True,
) -> dict:
    r"""Compute F-SAPT partitioned electrostatic interaction energy.

    Decomposes the total :math:`E^{(1)}_{\text{elst}}` (Eq. 4) into
    atom-pair and orbital-pair contributions for functional-group
    analysis (F-SAPT). Nuclear-nuclear, nuclear-electron, and
    electron-electron terms are accumulated into the ``Elst_AB``
    breakdown matrix stored in *cache*.

    Parameters
    ----------
    cache : dict
        SAPT data cache with localized orbitals from :func:`flocalization`.
    sapt_elst : dict
        Total SAPT electrostatic energies from :func:`electrostatics`.
    dimer_wfn : core.Wavefunction
        Dimer supermolecular wavefunction.
    wfn_A : core.Wavefunction
        Monomer A wavefunction.
    wfn_B : core.Wavefunction
        Monomer B wavefunction.
    jk : core.JK
        JK integral engine.
    do_print : bool, optional
        Whether to print output, by default True.

    Returns
    -------
    dict
        Updated *cache* dictionary with ``'Elst_AB'`` breakdown matrix.
    """
    core.timer_on("F-SAPT Elst Setup")
    if do_print:
        core.print_out("  ==> F-SAPT Electrostatics <==\n\n")

    link_assignment = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT").upper()
    mol = dimer_wfn.molecule()  # dimer molecule
    dimer_basis = dimer_wfn.basisset()
    nA_atoms = mol.natom()
    nB_atoms = mol.natom()

    # Sizing
    L0A = (
        cache["Locc_A"]
        if link_assignment not in {"SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"}
        else cache["Locc_A"]
    )
    L0B = (
        cache["Locc_B"]
        if link_assignment not in {"SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"}
        else cache["Locc_B"]
    )
    na = L0A.np.shape[1]
    nb = L0B.np.shape[1]

    # Initialize breakdown matrix (nA_atoms + na + 1, nB_atoms + nb + 1)
    Elst_AB = np.zeros((nA_atoms + na + 1, nB_atoms + nb + 1))

    # Terms for total
    Elst1_terms = np.zeros(4)  # [0]: a-B, [1]: A-b, [2]: a-b, [3]: nuc

    # Nuclear-nuclear interactions (A <-> B)
    ZA = cache["ZA"]
    ZB = cache["ZB"]

    # Vectorized nuclear-nuclear interactions
    # Build distance matrix
    for A in range(nA_atoms):
        for B in range(nB_atoms):
            if A == B:
                continue
            R = mol.xyz(A).distance(mol.xyz(B))
            if R == 0:
                continue
            E = ZA.np[A] * ZB.np[B] / R
            Elst_AB[A, B] = E
            Elst1_terms[3] += E

    # External A - atom B interactions
    if "A" in cache.get("external_potentials", {}):
        ext_pot_A = cache["external_potentials"]["A"]
        for B in range(nB_atoms):
            atom_mol = core.Molecule([core.Atom(ZB.np[B])])
            atom_mol.set_geometry([mol.xyz(B)])
            interaction = ext_pot_A.computeNuclearEnergy(atom_mol)
            Elst_AB[nA_atoms + na, B] = interaction
            Elst1_terms[3] += interaction

    # External B - atom A interactions
    if "B" in cache.get("external_potentials", {}):
        ext_pot_B = cache["external_potentials"]["B"]
        for A in range(nA_atoms):
            atom_mol = core.Molecule([core.Atom(ZA.np[A])])
            atom_mol.set_geometry([mol.xyz(A)])
            interaction = ext_pot_B.computeNuclearEnergy(atom_mol)
            Elst_AB[A, nB_atoms + nb] = interaction
            Elst1_terms[3] += interaction

    core.timer_off("F-SAPT Elst Setup")
    # => a <-> b (electron-electron interactions via DFHelper) <= //

    # Get auxiliary basis for density fitting
    aux_basis = dimer_wfn.get_basisset("DF_BASIS_SCF")

    # Create DFHelper object
    dfh = core.DFHelper(dimer_basis, aux_basis)

    # Set memory following fisapt.cc logic
    # Note: In C++, doubles_ is the total memory budget in doubles
    # Here we use a reasonable default or get from options if available
    memory_doubles = core.get_memory() // 8
    dfh.set_memory(memory_doubles)
    dfh.set_method("DIRECT_iaQ")
    dfh.set_nthreads(core.get_num_threads())
    dfh.initialize()
    dfh.print_header()

    # Create Matrix objects from numpy arrays for L0A and L0B
    L0A = core.Matrix.from_array(L0A.np)
    L0B = core.Matrix.from_array(L0B.np)

    # Add orbital spaces
    dfh.add_space("a", L0A)
    dfh.add_space("b", L0B)

    # Add transformations for diagonal blocks (a,a|Q) and (b,b|Q)
    dfh.add_transformation("Aaa", "a", "a")
    dfh.add_transformation("Abb", "b", "b")

    # Perform the transformation
    dfh.transform()

    # Extract diagonal 3-index integrals (vectorized)
    nQ = aux_basis.nbf()
    QaC = np.zeros((na, nQ))
    QbC = np.zeros((nb, nQ))

    # Process in batches for better memory efficiency
    batch_size = min(100, max(na, nb))

    # Extract Aaa diagonal elements
    for start_a in range(0, na, batch_size):
        end_a = min(start_a + batch_size, na)
        for a in range(start_a, end_a):
            tensor = dfh.get_tensor("Aaa", [a, a + 1], [a, a + 1], [0, nQ])
            QaC[a, :] = tensor.np.flatten()

    # Extract Abb diagonal elements
    for start_b in range(0, nb, batch_size):
        end_b = min(start_b + batch_size, nb)
        for b in range(start_b, end_b):
            tensor = dfh.get_tensor("Abb", [b, b + 1], [b, b + 1], [0, nQ])
            QbC[b, :] = tensor.np.flatten()

    # Compute electrostatic interaction: Elst10_3 = 4.0 * QaC @ QbC.T
    Elst10_3 = 4.0 * np.dot(QaC, QbC.T)

    # Store in breakdown matrix and accumulate total
    Elst1_terms[2] += np.sum(Elst10_3)
    Elst_AB[nA_atoms : nA_atoms + na, nB_atoms : nB_atoms + nb] += Elst10_3

    # Store QaC and QbC in cache for reuse in f-induction
    cache["Vlocc0A"] = core.Matrix.from_array(QaC)
    cache["Vlocc0B"] = core.Matrix.from_array(QbC)

    # Clear DFHelper spaces for next use
    dfh.clear_spaces()

    core.timer_on("F-SAPT Elst Final")
    # => A <-> b (nuclei A interacting with orbitals b) <= //
    L0B_mat = core.Matrix.from_array(L0B.np)
    L0B_mat.name = "L0B_mat"

    L0A_mat = core.Matrix.from_array(L0A.np)
    ext_pot = core.ExternalPotential()
    for A in range(nA_atoms):
        if ZA.np[A] == 0.0:
            continue

        ext_pot.clear()
        atom_pos = mol.xyz(A)
        ext_pot.addCharge(ZA.np[A], atom_pos[0], atom_pos[1], atom_pos[2])

        Vtemp = ext_pot.computePotentialMatrix(dimer_basis)
        Vtemp_mat = Vtemp.clone()
        Vtemp_mat.name = "Vtemp_mat"

        Vbb = chain_gemm_einsums([L0B_mat, Vtemp_mat, L0B_mat], ["T", "N", "N"])
        Vbb.name = "Vbb"

        # Vectorized diagonal extraction
        diag_Vbb = np.diag(Vbb.np)
        E_vec = 2.0 * diag_Vbb
        Elst1_terms[1] += np.sum(E_vec)
        Elst_AB[A, nB_atoms : nB_atoms + nb] += E_vec

    # Add external-A <-> orbital b interaction
    if "A" in cache.get("external_potentials", {}):
        ext_pot_A = cache["external_potentials"]["A"]
        Vtemp = ext_pot_A.computePotentialMatrix(dimer_basis)

        Vtemp_mat = Vtemp.clone()
        Vbb = chain_gemm_einsums([L0B_mat, Vtemp_mat, L0B_mat], ["T", "N", "N"])

        # Vectorized diagonal extraction
        diag_Vbb = np.diag(Vbb.np)
        E_vec = 2.0 * diag_Vbb
        Elst1_terms[1] += np.sum(E_vec)
        Elst_AB[nA_atoms + na, nB_atoms : nB_atoms + nb] += E_vec

    # => a <-> B (orbitals a interacting with nuclei B) <= //

    for B in range(nB_atoms):
        if ZB.np[B] == 0.0:
            continue

        ext_pot.clear()
        atom_pos = mol.xyz(B)
        ext_pot.addCharge(ZB.np[B], atom_pos[0], atom_pos[1], atom_pos[2])

        Vtemp = ext_pot.computePotentialMatrix(dimer_basis)

        Vtemp_mat = Vtemp.clone()
        Vaa = chain_gemm_einsums([L0A_mat, Vtemp_mat, L0A_mat], ["T", "N", "N"])

        # Vectorized diagonal extraction
        diag_Vaa = np.diag(Vaa.np)
        E_vec = 2.0 * diag_Vaa
        Elst1_terms[0] += np.sum(E_vec)
        Elst_AB[nA_atoms : nA_atoms + na, B] += E_vec

    # Add orbital a <-> external-B interaction
    if "B" in cache.get("external_potentials", {}):
        ext_pot_B = cache["external_potentials"]["B"]
        Vtemp = ext_pot_B.computePotentialMatrix(dimer_basis)

        Vtemp_mat = Vtemp.clone()
        Vaa = chain_gemm_einsums([L0A_mat, Vtemp_mat, L0A_mat], ["T", "N", "N"])

        # Vectorized diagonal extraction
        diag_Vaa = np.diag(Vaa.np)
        E_vec = 2.0 * diag_Vaa
        Elst1_terms[0] += np.sum(E_vec)
        Elst_AB[nA_atoms : nA_atoms + na, nB_atoms + nb] += E_vec

    # Clear DFHelper for next use
    dfh.clear_spaces()
    cache["dfh"] = dfh  # Store DFHelper in cache for potential reuse
    Elst10 = np.sum(Elst1_terms)
    core.print_out(f"    Elst10,r            = {Elst10 * 1000:.8f} [mEh]\n")
    # Ensure that partition matches SAPT elst energy. Should be equal to
    # numerical precision and effectively free to check assertion here.
    assert abs(Elst10 - sapt_elst) < 1e-8, (
        f"FELST: Localized Elst10,r does not match SAPT Elst10,r!\n{Elst10 =}, {sapt_elst}"
    )

    # Add extern-extern contribution if both external potentials exist
    if "A" in cache.get("external_potentials", {}) and "B" in cache.get(
        "external_potentials", {}
    ):
        ext_pot_A = cache["external_potentials"]["A"]
        ext_pot_B = cache["external_potentials"]["B"]
        ext_ext = ext_pot_A.computeExternExternInteraction(ext_pot_B) * 2.0
        Elst_AB[nA_atoms + na, nB_atoms + nb] += ext_ext

    # Store breakdown matrix in cache
    cache["Elst_AB"] = core.Matrix.from_array(Elst_AB)
    core.timer_off("F-SAPT Elst Final")
    return cache


def fexch(
    cache: dict,
    sapt_exch10_s2: float,
    sapt_exch10: float,
    dimer_wfn: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    jk: core.JK,
    do_print: bool = True,
) -> dict:
    r"""Compute the F-SAPT first-order exchange partitioning.

    Uses the Exch10(S^2) approximation with orbital partitioning and
    follows Eq. 6 of Xie et al. (2022). Note S = S^{AO}

    .. math::

        E^{(1)}_{\text{exch}}(S^2) = -2(\mathbf{P}^A \mathbf{S} \mathbf{P}^B \mathbf{S} \mathbf{P}^{A,\text{vir}}) \cdot \boldsymbol{\omega}^B
          - 2(\mathbf{P}^B \mathbf{S} \mathbf{P}^A \mathbf{S} \mathbf{P}^{B,\text{vir}}) \cdot \boldsymbol{\omega}^A
          - 2(\mathbf{P}^{A,\text{vir}} \mathbf{S} \mathbf{P}^B) \cdot \mathbf{K}[\mathbf{P}^A \mathbf{S} \mathbf{P}^{B,\text{vir}}]

    Parameters
    ----------
    cache : dict
        SAPT/F-SAPT cache containing localized orbitals and intermediates.
    sapt_exch10_s2 : float
        Total SAPT first-order exchange energy in the :math:`S^2` approximation.
    sapt_exch10 : float
        Total SAPT first-order exchange energy used for optional scaling.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction used for basis and molecular metadata.
    wfn_A : core.Wavefunction
        Monomer A wavefunction.
    wfn_B : core.Wavefunction
        Monomer B wavefunction.
    jk : core.JK
        JK object for Coulomb/exchange intermediates.
    do_print : bool, optional
        Whether to print exchange diagnostics, by default True.

    Returns
    -------
    dict
        Updated cache with ``Exch_AB`` matrix.
    """
    if do_print:
        core.print_out("  ==> F-SAPT Exchange <==\n\n")

    mol = dimer_wfn.molecule()
    nA_atoms = nB_atoms = mol.natom()
    na = cache["Locc_A"].shape[1]
    nb = cache["Locc_B"].shape[1]
    nr = cache["Cvir_A"].shape[1]
    ns = cache["Cvir_B"].shape[1]

    link_assignment = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT")
    na1 = na
    nb1 = nb
    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        na1 = na + 1
        nb1 = nb + 1

    Exch10_2 = 0.0
    Exch10_2_terms = [0.0, 0.0, 0.0]

    Exch_AB = np.zeros((nA_atoms + na1 + 1, nB_atoms + nb1 + 1))

    S = cache["S"]
    V_A = cache["V_A"]
    J_A = cache["J_A"]
    V_B = cache["V_B"]
    J_B = cache["J_B"]

    LoccA = cache["Locc_A"].clone()
    LoccA.name = "LoccA"
    LoccB = cache["Locc_B"].clone()
    LoccB.name = "LoccB"
    CvirA = cache["Cvir_A"]
    CvirB = cache["Cvir_B"]
    CvirA.name = "CvirA"
    CvirB.name = "CvirB"

    dfh = cache["dfh"]

    dfh.add_space("a", LoccA)
    dfh.add_space("r", CvirA)
    dfh.add_space("b", LoccB)
    dfh.add_space("s", CvirB)

    dfh.add_transformation("Aar", "a", "r")
    dfh.add_transformation("Abs", "b", "s")

    dfh.transform()

    W_A = V_A.clone()
    ein.core.axpy(2.0, J_A.np, W_A.np)
    W_A.name = "W_A"
    W_B = V_B.clone()
    ein.core.axpy(2.0, J_B.np, W_B.np)
    W_B.name = "W_B"

    WAbs = chain_gemm_einsums([LoccB, W_A, CvirB], ["T", "N", "N"])
    WBar = chain_gemm_einsums([LoccA, W_B, CvirA], ["T", "N", "N"])
    WAbs.name = "WAbs"
    WBar.name = "WBar"

    Sab = chain_gemm_einsums([LoccA, S, LoccB], ["T", "N", "N"])
    Sba = chain_gemm_einsums([LoccB, S, LoccA], ["T", "N", "N"])
    Sas = chain_gemm_einsums([LoccA, S, CvirB], ["T", "N", "N"])
    Sas.name = "Sas"
    Sab.name = "Sab"

    LoccB.name = "LoccB"
    CvirA.name = "CvirA"
    Sbr = chain_gemm_einsums([LoccB, S, CvirA], ["T", "N", "N"])

    Sab.name = "Sab"
    Sba.name = "Sba"
    Sas.name = "Sas"
    Sbr.name = "Sbr"

    WBab = chain_gemm_einsums([WBar, Sbr], ["N", "T"])
    WAba = chain_gemm_einsums([WAbs, Sas], ["N", "T"])
    WBab.name = "WBab"
    WAba.name = "WAba"

    E_exch1 = np.zeros((na, nb))
    E_exch2 = np.zeros((na, nb))

    for a in range(na):
        for b in range(nb):
            E_exch1[a, b] = -2.0 * Sab.np[a, b] * WBab.np[a, b]
            E_exch2[a, b] = -2.0 * Sba.np[b, a] * WAba.np[b, a]

    nQ = dimer_wfn.get_basisset("DF_BASIS_SCF").nbf()
    TrQ = core.Matrix("TrQ", nr, nQ)
    TsQ = core.Matrix("TsQ", ns, nQ)
    TbQ = core.Matrix("TbQ", nb, nQ)
    TaQ = core.Matrix("TaQ", na, nQ)

    dfh.add_disk_tensor("Bab", (na, nb, nQ))

    for a in range(na):
        TrQ.np[:, :] = dfh.get_tensor("Aar", [a, a + 1], [0, nr], [0, nQ]).np.reshape(
            nr, nQ
        )
        TbQ.np[:, :] = np.dot(Sbr.np, TrQ.np)
        dfh.write_disk_tensor("Bab", TbQ, [a, a + 1])

    dfh.add_disk_tensor("Bba", (nb, na, nQ))

    for b in range(nb):
        TsQ.np[:, :] = dfh.get_tensor("Abs", [b, b + 1], [0, ns], [0, nQ]).np.reshape(
            ns, nQ
        )
        TaQ.np[:, :] = np.dot(Sas.np, TsQ.np)
        dfh.write_disk_tensor("Bba", TaQ, [b, b + 1])

    E_exch3 = np.zeros((na, nb))

    for a in range(na):
        TbQ.np[:, :] = dfh.get_tensor("Bab", [a, a + 1], [0, nb], [0, nQ]).np.reshape(
            nb, nQ
        )
        for b in range(nb):
            TaQ_slice = dfh.get_tensor(
                "Bba", [b, b + 1], [a, a + 1], [0, nQ]
            ).np.reshape(nQ)
            E_exch3[a, b] = -2.0 * np.dot(TbQ.np[b, :], TaQ_slice)

    for a in range(na):
        for b in range(nb):
            Exch_AB[a + nA_atoms, b + nB_atoms] = (
                E_exch1[a, b] + E_exch2[a, b] + E_exch3[a, b]
            )
            Exch10_2_terms[0] += E_exch1[a, b]
            Exch10_2_terms[1] += E_exch2[a, b]
            Exch10_2_terms[2] += E_exch3[a, b]

    Exch10_2 = sum(Exch10_2_terms)

    if do_print:
        core.print_out(f"    Exch10(S^2)         = {Exch10_2 * 1000:18.10f} [mEh]\n")
        core.print_out(
            f"    Exch10(S^2)-true    = {sapt_exch10_s2 * 1000:18.10f} [mEh]\n"
        )
        core.print_out(f"    Exch10-true         = {sapt_exch10 * 1000:18.10f} [mEh]\n")
        core.print_out("\n")

    if core.get_option("FISAPT", "FISAPT_FSAPT_EXCH_SCALE"):
        scale = sapt_exch10 / Exch10_2
        Exch_AB *= scale
        if do_print:
            core.print_out(
                f"    Scaling F-SAPT Exch10(S^2) by {scale:11.3E} to match Exch10\n\n"
            )

    cache["Exch_AB"] = core.Matrix.from_array(Exch_AB)
    dfh.clear_spaces()
    return cache


def build_ind_pot(vars: dict) -> core.Matrix:
    r"""Build the induction potential in the MO basis for one monomer due to the other.

    Constructs :math:`\tilde{\boldsymbol{\omega}}^X` (Eq. 16 of Xie et al. 2022).
    By swapping A/B labels in ``vars``, the potential for either monomer can
    be computed.

    .. math::

        \tilde{\boldsymbol{\omega}}^X = (\mathbf{C}^{Y,\text{occ}})^\dagger \boldsymbol{\omega}^X \mathbf{C}^{Y,\text{vir}}

    Parameters
    ----------
    vars : dict
        Dictionary containing the matrices required for the induction
        potential build, including ``V_B``, ``J_B``, ``Cocc_A``, and ``Cvir_A``.

    Returns
    -------
    core.Matrix
        Induction potential in the occupied-virtual MO block.
    """
    w_B = vars["V_B"].clone()
    ein.core.axpy(2.0, vars["J_B"].np, w_B.np)
    return chain_gemm_einsums(
        [vars["Cocc_A"], w_B, vars["Cvir_A"]],
        ["T", "N", "N"],
    )


def build_exch_ind_pot_AB(vars: dict) -> core.Matrix:
    r"""Build the exchange-induction potential for monomer A due to monomer B.

    Constructs the exchange-induction operator in the MO basis following
    Eq. 17 of Xie et al. (2022), involving overlap-weighted density matrix
    products and Coulomb/exchange contractions.

    Parameters
    ----------
    vars : dict
        Dictionary of AO and MO intermediates required for the A<-B
        exchange-induction construction.

    Returns
    -------
    core.Matrix
        Exchange-induction potential for monomer A in the occupied-virtual
        MO block.
    """

    K_B = vars["K_B"]
    J_O = vars["J_O"]
    K_O = vars["K_O"]
    J_P_B = vars["J_P_B"]
    J_A = vars["J_A"]
    K_A = vars["K_A"]
    J_B = vars["J_B"]
    D_A = vars["D_A"]
    D_B = vars["D_B"]
    S = vars["S"]
    V_B = vars["V_B"]
    V_A = vars["V_A"]

    # Exch-Ind Potential A
    EX_A = K_B.clone()
    EX_A.scale(-1.0)
    ein.core.axpy(-2.0, J_O.np, EX_A.np)
    ein.core.axpy(1.0, K_O.np, EX_A.np)
    ein.core.axpy(2.0, J_P_B.np, EX_A.np)

    # Apply all the axpy operations to EX_A
    S_DB, S_DB_VA, S_DB_VA_DB_S = chain_gemm_einsums(
        [S, D_B, V_A, D_B, S], return_tensors=[True, True, False, True]
    )
    S_DB_JA, S_DB_JA_DB_S = chain_gemm_einsums(
        [S_DB, J_A, D_B, S], return_tensors=[True, False, True]
    )
    S_DB_S_DA, S_DB_S_DA_VB = chain_gemm_einsums(
        [S_DB, S, D_A, V_B],
        return_tensors=[False, True, True],
    )
    ein.core.axpy(-1.0, S_DB_VA.np, EX_A.np)
    ein.core.axpy(-2.0, S_DB_JA.np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([S_DB, K_A]).np, EX_A.np)
    ein.core.axpy(1.0, S_DB_S_DA_VB.np, EX_A.np)
    ein.core.axpy(2.0, chain_gemm_einsums([S_DB_S_DA, J_B]).np, EX_A.np)
    ein.core.axpy(1.0, S_DB_VA_DB_S.np, EX_A.np)
    ein.core.axpy(2.0, S_DB_JA_DB_S.np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([S_DB, K_O], ["N", "T"]).np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([V_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(-2.0, chain_gemm_einsums([J_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([K_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([V_B, D_A, S, D_B, S]).np, EX_A.np)
    ein.core.axpy(2.0, chain_gemm_einsums([J_B, D_A, S, D_B, S]).np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([K_O, D_B, S]).np, EX_A.np)

    EX_A_MO = chain_gemm_einsums(
        [vars["Cocc_A"], EX_A, vars["Cvir_A"]],
        ["T", "N", "N"],
    )
    return EX_A_MO


def build_exch_ind_pot_BA(vars: dict) -> core.Matrix:
    r"""Build the exchange-induction potential for monomer B due to monomer A.

    Analogous to :func:`build_exch_ind_pot_AB` with A/B roles swapped,
    following Eq. 17 of Xie et al. (2022).

    Parameters
    ----------
    vars : dict
        Dictionary of AO and MO intermediates required for the B<-A
        exchange-induction construction.

    Returns
    -------
    core.Matrix
        Exchange-induction potential for monomer B in the occupied-virtual
        MO block.
    """

    K_B = vars["K_B"]
    J_O = vars["J_O"]
    K_O = vars["K_O"]
    J_P_A = vars["J_P_A"]
    J_A = vars["J_A"]
    K_A = vars["K_A"]
    J_B = vars["J_B"]
    D_A = vars["D_A"]
    D_B = vars["D_B"]
    S = vars["S"]
    V_B = vars["V_B"]
    V_A = vars["V_A"]

    EX_B = K_A.clone()
    EX_B.scale(-1.0)
    ein.core.axpy(-2.0, J_O.np, EX_B.np)
    ein.core.axpy(1.0, K_O.np, EX_B.np.T)
    ein.core.axpy(2.0, J_P_A.np, EX_B.np)

    S_DA, S_DA_VB, S_DA_VB_DA_S = chain_gemm_einsums(
        [S, D_A, V_B, D_A, S], return_tensors=[True, True, False, True]
    )
    S_DA_JB, S_DA_JB_DA_S = chain_gemm_einsums(
        [S_DA, J_B, D_A, S], return_tensors=[True, False, True]
    )
    S_DA_S_DB, S_DA_S_DB_VA = chain_gemm_einsums(
        [S_DA, S, D_B, V_A],
        return_tensors=[False, True, True],
    )

    # Apply all the axpy operations to EX_B
    ein.core.axpy(-1.0, S_DA_VB.np, EX_B.np)
    ein.core.axpy(-2.0, S_DA_JB.np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([S_DA, K_B]).np, EX_B.np)
    ein.core.axpy(1.0, S_DA_S_DB_VA.np, EX_B.np)
    ein.core.axpy(2.0, chain_gemm_einsums([S_DA_S_DB, J_A]).np, EX_B.np)
    ein.core.axpy(1.0, S_DA_VB_DA_S.np, EX_B.np)
    ein.core.axpy(2.0, S_DA_JB_DA_S.np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([S_DA, K_O]).np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([V_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(-2.0, chain_gemm_einsums([J_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([K_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([V_A, D_B, S, D_A, S]).np, EX_B.np)
    ein.core.axpy(2.0, chain_gemm_einsums([J_A, D_B, S, D_A, S]).np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([K_O, D_A, S], ["T", "N", "N"]).np, EX_B.np)

    EX_B_MO = chain_gemm_einsums(
        [vars["Cocc_B"], EX_B, vars["Cvir_B"]],
        ["T", "N", "N"],
    )
    return EX_B_MO


def build_exch_ind_pot_avg(vars: dict) -> core.Matrix:
    r"""Build the averaged exchange-induction potential for link-orbital SAO/SIAO methods.

    Uses the older :func:`core.triplet` API for matrix triple products.
    This variant handles the SAO/SIAO link-orbital partitioning where
    average exchange-induction potentials are needed.

    Parameters
    ----------
    vars : dict
        Dictionary of AO and MO intermediates required to construct the
        averaged exchange-induction potential.

    Returns
    -------
    core.Matrix
        Averaged exchange-induction potential in the occupied-virtual
        MO block.
    """
    Ca = vars["Cocc_A"]
    Cr = vars["Cvir_A"]

    S = vars["S"]

    D_A = vars["D_A"]
    J_A = vars["J_A"]
    K_A = vars["K_A"]
    V_A = vars["V_A"]
    D_B = vars["D_B"]
    J_B = vars["J_B"]
    K_B = vars["K_B"]
    V_B = vars["V_B"]
    D_X = vars["D_X"]
    J_X = vars["J_X"]
    K_X = vars["K_X"]
    D_Y = vars["D_Y"]
    J_Y = vars["J_Y"]
    K_Y = vars["K_Y"]

    J_O = vars["J_O"]
    K_O = vars["K_O"]
    K_AOY = vars["K_AOY"]

    J_P = vars["J_P"]
    J_PYAY = vars["J_PYAY"]

    W = core.Matrix.from_array(-K_B.np)

    T = core.triplet(S, D_B, J_A, False, False, False)
    W.np[:] += -2.0 * T.np

    W.np[:] += K_O.np

    W.np[:] += -2.0 * J_O.np

    T = core.triplet(S, D_B, K_A, False, False, False)
    W.np[:] += T.np

    T = core.triplet(J_B, D_B, S, False, False, False)
    W.np[:] += -2.0 * T.np

    T = core.triplet(K_B, D_B, S, False, False, False)
    W.np[:] += T.np
    T = core.triplet(K_Y, D_Y, S, False, False, False)
    W.np[:] += T.np

    T1 = core.triplet(S, D_B, J_A, False, False, False)
    T = core.triplet(T1, D_B, S, False, False, False)
    W.np[:] += 2.0 * T.np
    T1 = core.triplet(S, D_Y, J_A, False, False, False)
    T = core.triplet(T1, D_Y, S, False, False, False)
    W.np[:] += 2.0 * T.np

    T1 = core.triplet(J_B, D_A, S, False, False, False)
    T = core.triplet(T1, D_B, S, False, False, False)
    W.np[:] += 2.0 * T.np

    T = core.triplet(K_O, D_B, S, False, False, False)
    W.np[:] += -1.0 * T.np
    T = core.triplet(K_AOY, D_Y, S, False, False, False)
    W.np[:] += -1.0 * T.np

    W.np[:] += 2.0 * J_P.np
    W.np[:] += 2.0 * J_PYAY.np

    T1 = core.triplet(S, D_B, S, False, False, False)
    T = core.triplet(T1, D_A, J_B, False, False, False)
    W.np[:] += 2.0 * T.np

    T = core.triplet(S, D_B, K_O, False, False, True)
    W.np[:] += -1.0 * T.np
    T = core.triplet(S, D_Y, K_AOY, False, False, True)
    W.np[:] += -1.0 * T.np

    T = core.triplet(S, D_B, V_A, False, False, False)
    W.np[:] += -1.0 * T.np

    T = core.triplet(V_B, D_B, S, False, False, False)
    W.np[:] += -1.0 * T.np

    T1 = core.triplet(S, D_B, V_A, False, False, False)
    T = core.triplet(T1, D_B, S, False, False, False)
    W.np[:] += T.np
    T1 = core.triplet(S, D_Y, V_A, False, False, False)
    T = core.triplet(T1, D_Y, S, False, False, False)
    W.np[:] += T.np

    T1 = core.triplet(V_B, D_A, S, False, False, False)
    T = core.triplet(T1, D_B, S, False, False, False)
    W.np[:] += T.np

    T1 = core.triplet(S, D_B, S, False, False, False)
    T = core.triplet(T1, D_A, V_B, False, False, False)
    W.np[:] += T.np

    return core.triplet(Ca, W, Cr, True, False, False)


def find(
    cache: dict,
    scalars: dict,
    dimer_wfn: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    jk: core.JK,
    do_print: bool = True,
) -> dict:
    r"""Compute the F-SAPT induction partitioning.

    Partitions the second-order induction and exchange-induction energies into
    atomic pair contributions for F-SAPT analysis. Computes both uncoupled and
    (optionally) coupled induction using the CPSCF solver.

    Parameters
    ----------
    cache : dict
        SAPT data cache containing orbital coefficients, density matrices, and integrals.
    scalars : dict
        Reference scalar energies for validation.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction.
    wfn_A : core.Wavefunction
        Monomer A wavefunction.
    wfn_B : core.Wavefunction
        Monomer B wavefunction.
    jk : core.JK
        JK integral engine.
    do_print : bool, optional
        Whether to print results, by default True.

    Returns
    -------
    dict
        Updated cache with ``Ind_AB`` and ``IndAB_AB``/``IndBA_AB`` matrices.
    """
    if do_print:
        core.print_out("  ==> F-SAPT Induction <==\n\n")

    ind_scale = core.get_option("FISAPT", "FISAPT_FSAPT_IND_SCALE")
    link_assignment = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT")

    mol = dimer_wfn.molecule()
    nA = mol.natom()
    nB = mol.natom()
    na = cache["Locc_A"].shape[1]
    nb = cache["Locc_B"].shape[1]
    nr = cache["Cvir_A"].shape[1]
    ns = cache["Cvir_B"].shape[1]

    na1 = na
    nb1 = nb
    # for the SAOn/SIAOn variants, we sometimes need na1 = na+1 (with link
    # orbital) and sometimes na (without) - be careful with this!
    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        na1 = na + 1
        nb1 = nb + 1

    Locc_A = cache["Locc_A"].clone()
    Locc_A.name = "LoccA"
    Locc_B = cache["Locc_B"].clone()
    Locc_B.name = "LoccB"

    Uocc_A = cache["Uocc_A"]
    Uocc_B = cache["Uocc_B"]

    Cocc_A = cache["Cocc_A"]
    Cocc_B = cache["Cocc_B"]
    Cvir_A = cache["Cvir_A"]
    Cvir_B = cache["Cvir_B"]

    # Cvir_A.set_name("Cvir_A")
    # print(Cvir_A)

    eps_occ_A = cache["eps_occ_A"]
    eps_occ_B = cache["eps_occ_B"]
    eps_vir_A = cache["eps_vir_A"]
    eps_vir_B = cache["eps_vir_B"]

    # Collect relevant variables
    S = cache["S"]
    D_A = cache["D_A"]
    V_A = cache["V_A"]
    J_A = cache["J_A"]
    K_A = cache["K_A"]
    D_B = cache["D_B"]
    V_B = cache["V_B"]
    J_B = cache["J_B"]
    K_B = cache["K_B"]
    J_O = cache["J_O"]
    K_O = cache["K_O"]
    J_P_A = cache["J_P_A"]
    J_P_B = cache["J_P_B"]

    aux_basis = dimer_wfn.get_basisset("DF_BASIS_SCF")
    nQ = aux_basis.nbf()

    dfh = cache["dfh"]

    # ESPs - external potential entries
    dfh.add_disk_tensor("WBar", (nB + nb1 + 1, na, nr))
    dfh.add_disk_tensor("WAbs", (nA + na1 + 1, nb, ns))

    # Nuclear Contribution to ESPs
    ext_pot = core.ExternalPotential()
    ZA = cache["ZA"].np
    for A in range(nA):
        ext_pot.clear()
        atom_pos = mol.xyz(A)
        ext_pot.addCharge(ZA[A], atom_pos[0], atom_pos[1], atom_pos[2])
        Vtemp = ext_pot.computePotentialMatrix(dimer_wfn.basisset())
        Vbs = core.Matrix.from_array(
            chain_gemm_einsums([Cocc_B, Vtemp, Cvir_B], ["T", "N", "N"])
        )
        # Vbs_A doesn't agree... Cocc_B and Cvir_B
        dfh.write_disk_tensor("WAbs", Vbs, (A, A + 1))

    ZB = cache["ZB"].np
    for B in range(nB):
        ext_pot.clear()
        atom_pos = mol.xyz(B)
        ext_pot.addCharge(ZB[B], atom_pos[0], atom_pos[1], atom_pos[2])
        Vtemp = ext_pot.computePotentialMatrix(dimer_wfn.basisset())
        Var = core.Matrix.from_array(
            chain_gemm_einsums([Cocc_A, Vtemp, Cvir_A], ["T", "N", "N"])
        )
        dfh.write_disk_tensor("WBar", Var, (B, B + 1))

    dfh.add_space("a", core.Matrix.from_array(Cocc_A))
    dfh.add_space("r", core.Matrix.from_array(Cvir_A))
    dfh.add_space("b", core.Matrix.from_array(Cocc_B))
    dfh.add_space("s", core.Matrix.from_array(Cvir_B))

    dfh.add_transformation("Aar", "a", "r")
    dfh.add_transformation("Abs", "b", "s")

    dfh.transform()

    RaC = cache["Vlocc0A"]  # na x nQ
    RbD = cache["Vlocc0B"]  # nb x nQ

    TsQ = core.Matrix("TsQ", ns, nQ)
    T1As = core.Matrix("T1As", na1, ns)
    # print(f"{na1 = }, {nb1 = }, {ns = }, {nr = }")
    # print(f"{na1 = }, {nb1 = }, {ns = }, {nQ = }")
    for B in range(nb):
        # print(f"{TsQ.np.shape =}")
        # TODO: CONTINUE HERE
        # fill_tensor is not working properly with 2D slices yet...
        # dfh.fill_tensor("Abs", TsQ, [B, B + 1])
        dfh.fill_tensor("Abs", TsQ, [B, B + 1], [0, ns], [0, nQ])
        TsQ = core.Matrix.from_array(TsQ.np[0, :, :])
        T1As.gemm(False, True, 2.0, RaC, TsQ, 0.0)
        for A in range(na1):
            row_view = core.Matrix.from_array(T1As.np[A : A + 1, :])
            dfh.write_disk_tensor("WAbs", row_view, (nA + A, nA + A + 1), (B, B + 1))

    TrQ = core.Matrix("TrQ", nr, nQ)
    T1Br = core.Matrix("T1Br", nb1, nr)
    for A in range(na):
        # dfh.fill_tensor("Abs", TsQ, [B, B + 1], [0, ns], [0, nQ])
        # TsQ = core.Matrix.from_array(TsQ.np[0, :, :])
        dfh.fill_tensor("Aar", TrQ, [A, A + 1], [0, nr], [0, nQ])
        TrQ = core.Matrix.from_array(TrQ.np[0, :, :])
        T1Br.gemm(False, True, 2.0, RbD, TrQ, 0.0)
        for B in range(nb1):
            row_view = core.Matrix.from_array(T1Br.np[B : B + 1, :])
            dfh.write_disk_tensor("WBar", row_view, (nB + B, nB + B + 1), (A, A + 1))

    xA = core.Matrix("xA", na, nr)
    xB = core.Matrix("xB", nb, ns)
    wB = core.Matrix("wB", na, nr)
    wA = core.Matrix("wA", nb, ns)

    uAT = core.Matrix("uAT", nb, ns)
    wAT = core.Matrix("wAT", nb, ns)
    uBT = core.Matrix("uBT", na, nr)
    wBT = core.Matrix("wBT", na, nr)

    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        D_X = core.doublet(cache["thislinkA"], cache["thislinkA"], False, True)
        D_Y = core.doublet(cache["thislinkB"], cache["thislinkB"], False, True)
        J_X = cache["JLA"]
        K_X = cache["KLA"]
        J_Y = cache["JLB"]
        K_Y = cache["KLB"]

        K_AOY = cache["K_AOY"]
        K_XOB = core.Matrix.from_array(cache["K_XOB"].np.T)
        J_P_YAY = cache["J_P_YAY"]
        J_P_XBX = cache["J_P_XBX"]

        mapA = {
            "Cocc_A": Locc_A,
            "Cvir_A": Cvir_A,
            "S": S,
            "D_A": D_A,
            "V_A": V_A,
            "J_A": J_A,
            "K_A": K_A,
            "D_B": D_B,
            "V_B": V_B,
            "J_B": J_B,
            "K_B": K_B,
            "D_X": D_X,
            "J_X": J_X,
            "K_X": K_X,
            "D_Y": D_Y,
            "J_Y": J_Y,
            "K_Y": K_Y,
            "J_O": J_O,
            "K_O": K_O,
            "K_AOY": K_AOY,
            "J_P": J_P_A,
            "J_PYAY": J_P_YAY,
        }

        raise NotImplementedError("find() not ready yet for link orbitals")
        wBT = build_ind_pot(mapA)
        uBT = build_exch_ind_pot_avg(mapA)

        K_O.np = K_O.np
        K_O.np[:] = K_O.np.T

        mapB = {
            "Cocc_A": Locc_B,
            "Cvir_A": Cvir_B,
            "S": S,
            "D_A": D_B,
            "V_A": V_B,
            "J_A": J_B,
            "K_A": K_B,
            "D_B": D_A,
            "V_B": V_A,
            "J_B": J_A,
            "K_B": K_A,
            "D_X": D_Y,
            "J_X": J_Y,
            "K_X": K_Y,
            "D_Y": D_X,
            "J_Y": J_X,
            "K_Y": K_X,
            "J_O": J_O,
            "K_O": K_O,
            "K_AOY": K_XOB,
            "J_P": J_P_B,
            "J_PYAY": J_P_XBX,
        }

        wAT = build_ind_pot(mapB)
        uAT = build_exch_ind_pot_avg(mapB)

        K_O.np[:] = K_O.np.T

    else:
        mapA = {
            "S": S,
            "J_O": J_O,
            "K_O": K_O,
            "Cocc_A": Locc_A,
            "Cvir_A": Cvir_A,
            "D_A": D_A,
            "V_A": V_A,
            "J_A": J_A,
            "K_A": K_A,
            "J_P_A": J_P_A,
            "Cocc_B": Locc_B,
            "Cvir_B": Cvir_B,
            "D_B": D_B,
            "V_B": V_B,
            "J_B": J_B,
            "K_B": K_B,
            "J_P_B": J_P_B,
        }

        # V_B and J_B are equivalent, but Cocc_A and Cvir_A differ from FISAPT0... is there a critical disagreement that would only surface here? Seems unlikely.
        # Locc_A magnitudes are about the same with different signs, but that is okay.
        # Cvir_A does not have the same magnitude for most terms... This is an issue.
        wBT = build_ind_pot(
            {
                "V_B": V_B,
                "J_B": J_B,
                "Cocc_A": Locc_A,
                "Cvir_A": Cvir_A,
            }
        )
        wAT = build_ind_pot(
            {
                "V_B": V_A,
                "J_B": J_A,
                "Cocc_A": Locc_B,
                "Cvir_A": Cvir_B,
            }
        )
        uBT = build_exch_ind_pot_AB(mapA)
        uAT = build_exch_ind_pot_BA(mapA)

    wBT.name = "wBT"
    uBT.name = "uBT"
    wAT.name = "wAT"
    uAT.name = "uAT"
    V_B.name = "V_B"
    J_B.name = "J_B"

    Ind20u_AB_terms = core.Matrix("Ind20 [A<-B] (a x B)", na, nB + nb1 + 1)
    Ind20u_BA_terms = core.Matrix("Ind20 [B<-A] (A x b)", nA + na1 + 1, nb)
    Ind20u_AB_termsp = Ind20u_AB_terms.np
    Ind20u_BA_termsp = Ind20u_BA_terms.np

    Ind20u_AB = 0.0
    Ind20u_BA = 0.0

    ExchInd20u_AB_terms = core.Matrix("ExchInd20 [A<-B] (a x B)", na, nB + nb1 + 1)
    ExchInd20u_BA_terms = core.Matrix("ExchInd20 [B<-A] (A x b)", nA + na1 + 1, nb)
    ExchInd20u_AB_termsp = ExchInd20u_AB_terms.np
    ExchInd20u_BA_termsp = ExchInd20u_BA_terms.np

    ExchInd20u_AB = 0.0
    ExchInd20u_BA = 0.0

    # sna = snB = snb = snA = 0
    # sExchInd20u_AB_terms = core.Matrix("sExchInd20 [A<-B] (a x B)", sna, snB + snb + 1)
    # sExchInd20u_BA_terms = core.Matrix("sExchInd20 [B<-A] (A x b)", snA + sna + 1, snb)
    # sExchInd20u_AB_termsp = sExchInd20u_AB_terms.np
    # sExchInd20u_BA_termsp = sExchInd20u_BA_terms.np
    #
    # sExchInd20u_AB = 0.0
    # sExchInd20u_BA = 0.0

    Indu_AB_terms = core.Matrix("Ind [A<-B] (a x B)", na, nB + nb1 + 1)
    Indu_BA_terms = core.Matrix("Ind [B<-A] (A x b)", nA + na1 + 1, nb)
    Indu_AB = 0.0
    Indu_BA = 0.0

    # Commented out terms are for sSAPT0 scaling... do we really want this?
    # sIndu_AB_terms = core.Matrix("sInd [A<-B] (a x B)", sna, snB + snb + 1)
    # sIndu_BA_terms = core.Matrix("sInd [B<-A] (A x b)", snA + sna + 1, snb)
    # sIndu_AB_termsp = sIndu_AB_terms.np
    # sIndu_BA_termsp = sIndu_BA_terms.np
    # sIndu_AB = 0.0
    # sIndu_BA = 0.0

    # ==> A <- B Uncoupled <==
    if dimer_wfn.has_potential_variable("B"):
        Var = core.triplet(Cocc_A, cache["VB_extern"], Cvir_A, True, False, False)
        dfh.write_disk_tensor("WBar", Var, (nB + nb1, nB + nb1 + 1))
    else:
        Var = core.Matrix("zero", na, nr)
        Var.zero()
        dfh.write_disk_tensor("WBar", Var, (nB + nb1, nB + nb1 + 1))

    for B in range(nB + nb1 + 1):  # add one for external potential
        # ESP
        dfh.fill_tensor("WBar", wB, [B, B + 1])
        # Uncoupled
        for a in range(na):
            for r in range(nr):
                # fill_tensor wB as (1, na, nr), so we take first index only
                xA.np[a, r] = wB.np[0, a, r] / (eps_occ_A.np[a] - eps_vir_A.np[r])

        x2A = core.doublet(Uocc_A, xA, True, False)
        x2Ap = x2A.np

        for a in range(na):
            Jval = 2.0 * np.dot(x2Ap[a, :], wBT.np[a, :])
            Kval = 2.0 * np.dot(x2Ap[a, :], uBT.np[a, :])
            Ind20u_AB += Jval
            ExchInd20u_AB_termsp[a, B] = Kval
            ExchInd20u_AB += Kval
            Ind20u_AB_termsp[a, B] = Jval
            # if core.get_option("SAPT", "SSAPT0_SCALE"):
            #     sExchInd20u_AB_termsp[a, B] = Kval
            #     sExchInd20u_AB += Kval
            #     sIndu_AB_termsp[a, B] = Jval + Kval
            #     sIndu_AB += Jval + Kval

            Indu_AB_terms.np[a, B] = Jval + Kval
            Indu_AB += Jval + Kval

    # ==> B <- A Uncoupled <==
    if dimer_wfn.has_potential_variable("A"):
        Vbs = core.triplet(Cocc_B, cache["VA_extern"], Cvir_B, True, False, False)
        dfh.write_disk_tensor("WAbs", Vbs, (nA + na1, nA + na1 + 1))
    else:
        Vbs = core.Matrix("zero", nb, ns)
        Vbs.zero()
        dfh.write_disk_tensor("WAbs", Vbs, (nA + na1, nA + na1 + 1))

    for A in range(nA + na1 + 1):
        dfh.fill_tensor("WAbs", wA, [A, A + 1])
        for b in range(nb):
            for s in range(ns):
                xB.np[b, s] = wA.np[0, b, s] / (eps_occ_B.np[b] - eps_vir_B.np[s])

        x2B = core.doublet(Uocc_B, xB, True, False)
        x2Bp = x2B.np

        for b in range(nb):
            Jval = 2.0 * np.dot(x2Bp[b, :], wAT.np[b, :])
            Kval = 2.0 * np.dot(x2Bp[b, :], uAT.np[b, :])
            Ind20u_BA_termsp[A, b] = Jval
            Ind20u_BA += Jval
            ExchInd20u_BA_termsp[A, b] = Kval
            ExchInd20u_BA += Kval
            # if core.get_option("SAPT", "SSAPT0_SCALE"):
            #     sExchInd20u_BA_termsp[A, b] = Kval
            #     sExchInd20u_BA += Kval
            #     sIndu_BA_termsp[A, b] = Jval + Kval
            #     sIndu_BA += Jval + Kval

            Indu_BA_terms.np[A, b] = Jval + Kval
            Indu_BA += Jval + Kval

    # Currently Ind20 and Exch-Ind are qualitatively coming out with wrong sign even...
    if do_print:
        core.print_out(
            f"    Ind20,u (A<-B)          = {Ind20u_AB * 1000:18.8f} [mEh]\n"
        )
        core.print_out(
            f"    Ind20,u (B<-A)          = {Ind20u_BA * 1000:18.8f} [mEh]\n"
        )
        assert (
            abs(scalars["Ind20,u (A<-B)"] - Ind20u_AB) < 1e-8
        ), f"Ind20u_AB mismatch: {1000 * scalars['Ind20,u (A<-B)']:.8f} vs {
            1000 * Ind20u_AB:.8f}"
        assert (
            abs(scalars["Ind20,u (A->B)"] - Ind20u_BA) < 1e-8
        ), f"Ind20u_BA mismatch: {1000 * scalars['Ind20,u (A->B)']:.8f} vs {
            1000 * Ind20u_BA:.8f}"
        core.print_out(
            f"    Ind20,u                 = {
                Ind20u_AB + Ind20u_BA * 1000:18.8f} [mEh]\n"
        )
        core.print_out(
            f"    Exch-Ind20,u (A<-B)     = {ExchInd20u_AB * 1000:18.8f} [mEh]\n"
        )
        core.print_out(
            f"    Exch-Ind20,u (B<-A)     = {ExchInd20u_BA * 1000:18.8f} [mEh]\n"
        )
        assert (
            abs(scalars["Exch-Ind20,u (A<-B)"] - ExchInd20u_AB) < 1e-8
        ), f"ExchInd20u_AB mismatch: {1000 * scalars['Exch-Ind20,u (A<-B)']:.8f} vs {
            1000 * ExchInd20u_AB:.8f}"
        assert (
            abs(scalars["Exch-Ind20,u (A->B)"] - ExchInd20u_BA) < 1e-8
        ), f"ExchInd20u_BA mismatch: {1000 * scalars['Exch-Ind20,u (A->B)']:.8f} vs {
            1000 * ExchInd20u_BA:.8f}"
        core.print_out(
            f"    Exch-Ind20,u            = {ExchInd20u_AB + ExchInd20u_BA * 1000:18.8f} [mEh]\n\n"
        )

    # Induction scaling
    if ind_scale:
        dHF = scalars.get("Delta HF Correction", 0.0)
        IndHF = scalars["Ind20,r"] + scalars["Exch-Ind20,r"] + dHF
        IndSAPT0 = scalars["Ind20,r"] + scalars["Exch-Ind20,r"]

        Sdelta = IndHF / IndSAPT0

        # NOTE: if doing ind_resp, logic below needs adjusted
        SrAB = (scalars["Ind20,r (A<-B)"] + scalars["Exch-Ind20,r (A<-B)"]) / (
            scalars["Ind20,u (A<-B)"] + scalars["Exch-Ind20,u (A<-B)"]
        )
        SrBA = (scalars["Ind20,r (A->B)"] + scalars["Exch-Ind20,r (A->B)"]) / (
            scalars["Ind20,u (A->B)"] + scalars["Exch-Ind20,u (A->B)"]
        )

        if do_print:
            core.print_out(f"    Scaling for delta HF        = {Sdelta:11.3E}\n")
            core.print_out(f"    Scaling for response (A<-B) = {SrAB:11.3E}\n")
            core.print_out(f"    Scaling for response (A->B) = {SrBA:11.3E}\n")
            core.print_out(f"    Scaling for total (A<-B)    = {Sdelta * SrAB:11.3E}\n")
            core.print_out(f"    Scaling for total (A->B)    = {Sdelta * SrBA:11.3E}\n")
            core.print_out("\n")

        # Apply scaling to all terms
        Indu_AB_terms.scale(Sdelta * SrAB)
        Indu_BA_terms.scale(Sdelta * SrBA)
        Ind20u_AB_terms.scale(Sdelta * SrAB)
        ExchInd20u_AB_terms.scale(Sdelta * SrAB)
        Ind20u_BA_terms.scale(Sdelta * SrBA)
        ExchInd20u_BA_terms.scale(Sdelta * SrBA)

        # Apply SSAPT0 scaling if enabled
        # if "sExch-Ind20,r" in scalars:
        #     sIndu_AB_terms.scale(sSdelta * sSrAB)
        #     sIndu_BA_terms.scale(sSdelta * sSrBA)

    IndAB_AB = core.Matrix("IndAB_AB", nA + na1 + 1, nB + nb1 + 1)
    IndBA_AB = core.Matrix("IndBA_AB", nA + na1 + 1, nB + nb1 + 1)

    # Final assembly might be wrong... backtrace time
    for a in range(na):
        for B in range(nB + nb1 + 1):
            IndAB_AB.np[a + nA, B] = Ind20u_AB_termsp[a, B] + ExchInd20u_AB_termsp[a, B]
    for A in range(nA + na1 + 1):
        for b in range(nb):
            IndBA_AB.np[A, b + nB] = Ind20u_BA_termsp[A, b] + ExchInd20u_BA_termsp[A, b]

    cache["IndAB_AB"] = IndAB_AB
    cache["IndBA_AB"] = IndBA_AB

    # if core.get_option("SAPT", "SSAPT0_SCALE"):
    #     cache["sExchInd20u_AB"] = sExchInd20u_AB
    #     cache["sExchInd20u_BA"] = sExchInd20u_BA
    #     cache["sIndu_AB"] = sIndu_AB
    #     cache["sIndu_BA"] = sIndu_BA

    """
    Ind20,u (A<-B)      =    -0.000005862 [mEh]
    Ind20,u (B<-A)      =    -0.000003086 [mEh]
    Ind20,u             =    -0.000008949 [mEh]
    Exch-Ind20,u (A<-B) =     0.000000887 [mEh]
    Exch-Ind20,u (B<-A) =     0.000000291 [mEh]
    Exch-Ind20,u        =     0.000001178 [mEh]
    """

    # NOT IMPLEMENTED YET
    # if (ind_resp) {
    #     outfile->Printf("  COUPLED INDUCTION (You asked for it!):\n\n");
    dfh.clear_all()
    return cache


def fdisp0(
    cache: dict,
    scalars: dict,
    dimer_wfn: core.Wavefunction,
    wfn_A: core.Wavefunction,
    wfn_B: core.Wavefunction,
    jk: core.JK,
    do_print: bool = True,
) -> dict:
    r"""Compute the F-SAPT0 dispersion partitioning.

    Partitions the second-order dispersion and exchange-dispersion energies
    into atomic pair contributions for F-SAPT analysis. Note, this does not use
    DFT energies and is a Hartree-Fock (SAPT0) dispersion energy only.

    Parameters
    ----------
    cache : dict
        SAPT data cache containing orbital coefficients, density matrices, and integrals.
    scalars : dict
        Reference scalar energies for validation.
    dimer_wfn : core.Wavefunction
        Dimer wavefunction.
    wfn_A : core.Wavefunction
        Monomer A wavefunction.
    wfn_B : core.Wavefunction
        Monomer B wavefunction.
    jk : core.JK
        JK integral engine.
    do_print : bool, optional
        Whether to print results, by default True.

    Returns
    -------
    dict
        Updated cache with ``Disp_AB`` matrix and ``Disp20,u``/``Exch-Disp20,u`` scalar energies.
    """
    if do_print:
        core.print_out("  ==> F-SAPT0 Dispersion <==\n\n")

    core.timer_on("F-SAPT Disp Setup")
    # ind_scale = core.get_option("FISAPT", "FISAPT_FSAPT_IND_SCALE")
    link_assignment = core.get_option("FISAPT", "FISAPT_LINK_ASSIGNMENT")

    mol = dimer_wfn.molecule()
    dimer_basis = dimer_wfn.basisset()
    nA = mol.natom()
    nB = mol.natom()
    nfa = cache["Lfocc0A"].shape[1]
    nfb = cache["Lfocc0B"].shape[1]
    # Use active occupied dimensions (excluding frozen core) to match C++ FISAPT fdisp
    na = cache["Caocc0A"].shape[1]
    nb = cache["Caocc0B"].shape[1]
    nr = cache["Cvir_A"].shape[1]
    ns = cache["Cvir_B"].shape[1]
    nn = cache["Cocc_A"].shape[0]  # number of AO basis functions

    na1 = na
    nb1 = nb

    # Disp_AB = core.Matrix("Disp_AB", nA + nfa + na1 + 1, nB + nfb + nb1 + 1)
    snA = 0
    snfa = 0
    sna = 0
    snB = 0
    snfb = 0
    snb = 0
    # if options_.get_bool("FISAPT", "FISAPT_SSAPT0_SCALE"):
    #     snA = nA
    #     snfa = nfa
    #     sna = na
    #     snB = nB
    #     snfb = nfb
    #     snb = nb

    if link_assignment in ["SAO0", "SAO1", "SAO2", "SIAO0", "SIAO1", "SIAO2"]:
        na1 = na + 1
        nb1 = nb + 1

    Locc_A = cache["Locc_A"].clone()
    Locc_A.name = "LoccA"
    Locc_B = cache["Locc_B"].clone()
    Locc_B.name = "LoccB"

    # Use active occupied orbitals (excluding frozen core) to match C++ FISAPT fdisp
    Cocc_A = cache["Caocc0A"]
    Cocc_B = cache["Caocc0B"]
    Cvir_A = cache["Cvir_A"]
    Cvir_B = cache["Cvir_B"]

    # Use only active occupied orbital energies (skip frozen core)
    # nfa and nfb are already defined at the start of fdisp0
    eps_occ_A = core.Vector.from_array(cache["eps_occ_A"].np[nfa:])
    eps_occ_B = core.Vector.from_array(cache["eps_occ_B"].np[nfb:])
    eps_vir_A = cache["eps_vir_A"]
    eps_vir_B = cache["eps_vir_B"]

    # Collect relevant variables
    S = cache["S"]
    D_A = cache["D_A"]
    P_A = cache["P_A"]
    V_A = cache["V_A"]
    J_A = cache["J_A"]
    K_A = cache["K_A"]
    D_B = cache["D_B"]
    P_B = cache["P_B"]
    V_B = cache["V_B"]
    J_B = cache["J_B"]
    K_B = cache["K_B"]
    K_O = cache["K_O"]

    aux_basis = dimer_wfn.get_basisset("DF_BASIS_SCF")
    nQ = aux_basis.nbf()

    # => Auxiliary C matrices <= //
    # Cr1 = (I - D_B * S) * Cvir_A
    Cr1 = chain_gemm_einsums([D_B, S, Cvir_A])
    ein.core.axpy(-1.0, Cvir_A.np, Cr1.np)

    # Cs1 = (I - D_A * S) * Cvir_B
    Cs1 = chain_gemm_einsums([D_A, S, Cvir_B])
    ein.core.axpy(-1.0, Cvir_B.np, Cs1.np)

    # Ca2 = D_B * S * Cocc_A
    Ca2 = chain_gemm_einsums([D_B, S, Cocc_A])

    # Cb2 = D_A * S * Cocc_B
    Cb2 = chain_gemm_einsums([D_A, S, Cocc_B])

    # Cr3 = 2 * (D_B * S * Cvir_A - D_A * S * D_B * S * Cvir_A)
    Cr3 = chain_gemm_einsums([D_B, S, Cvir_A])
    CrX = chain_gemm_einsums([D_A, S, D_B, S, Cvir_A])
    Cr3.subtract(CrX)
    Cr3.scale(2.0)

    # Cs3 = 2 * (D_A * S * Cvir_B - D_B * S * D_A * S * Cvir_B)
    Cs3 = chain_gemm_einsums([D_A, S, Cvir_B])
    CsX = chain_gemm_einsums([D_B, S, D_A, S, Cvir_B])
    Cs3.subtract(CsX)
    Cs3.scale(2.0)

    # Ca4 = -2 * D_A * S * D_B * S * Cocc_A
    Ca4 = chain_gemm_einsums([D_A, S, D_B, S, Cocc_A])
    Ca4.scale(-2.0)

    # Cb4 = -2 * D_B * S * D_A * S * Cocc_B
    Cb4 = chain_gemm_einsums([D_B, S, D_A, S, Cocc_B])
    Cb4.scale(-2.0)

    # => Auxiliary V matrices <= #

    # Jbr = 2.0 * Cocc_B.T @ J_A @ Cvir_A
    Jbr = chain_gemm_einsums([Cocc_B, J_A, Cvir_A], ["T", "N", "N"])
    Jbr.scale(2.0)

    # Kbr = -1.0 * Cocc_B.T @ K_A @ Cvir_A
    Kbr = chain_gemm_einsums([Cocc_B, K_A, Cvir_A], ["T", "N", "N"])
    Kbr.scale(-1.0)

    # Jas = 2.0 * Cocc_A.T @ J_B @ Cvir_B
    Jas = chain_gemm_einsums([Cocc_A, J_B, Cvir_B], ["T", "N", "N"])
    Jas.scale(2.0)

    # Kas = -1.0 * Cocc_A.T @ K_B @ Cvir_B
    Kas = chain_gemm_einsums([Cocc_A, K_B, Cvir_B], ["T", "N", "N"])
    Kas.scale(-1.0)

    # KOas = 1.0 * Cocc_A.T @ K_O @ Cvir_B
    KOas = chain_gemm_einsums([Cocc_A, K_O, Cvir_B], ["T", "N", "N"])

    # KObr = 1.0 * Cocc_B.T @ K_O.T @ Cvir_A
    # Note: K_O is transposed (second 'T' in the transpose list)
    KObr = chain_gemm_einsums([Cocc_B, K_O, Cvir_A], ["T", "T", "N"])

    # JBas = -2.0 * (Cocc_A.T @ S @ D_B) @ J_A @ Cvir_B
    temp_JBas = chain_gemm_einsums([Cocc_A, S, D_B], ["T", "N", "N"])
    JBas = chain_gemm_einsums([temp_JBas, J_A, Cvir_B], ["N", "N", "N"])
    JBas.scale(-2.0)

    # JAbr = -2.0 * (Cocc_B.T @ S @ D_A) @ J_B @ Cvir_A
    temp_JAbr = chain_gemm_einsums([Cocc_B, S, D_A], ["T", "N", "N"])
    JAbr = chain_gemm_einsums([temp_JAbr, J_B, Cvir_A], ["N", "N", "N"])
    JAbr.scale(-2.0)

    # Jbs = 4.0 * Cocc_B.T @ J_A @ Cvir_B
    Jbs = chain_gemm_einsums([Cocc_B, J_A, Cvir_B], ["T", "N", "N"])
    Jbs.scale(4.0)

    # Jar = 4.0 * Cocc_A.T @ J_B @ Cvir_A
    Jar = chain_gemm_einsums([Cocc_A, J_B, Cvir_A], ["T", "N", "N"])
    Jar.scale(4.0)

    # JAas = -2.0 * (Cocc_A.T @ J_B @ D_A) @ S @ Cvir_B
    temp_JAas = chain_gemm_einsums([Cocc_A, J_B, D_A], ["T", "N", "N"])
    JAas = chain_gemm_einsums([temp_JAas, S, Cvir_B], ["N", "N", "N"])
    JAas.scale(-2.0)

    # JBbr = -2.0 * (Cocc_B.T @ J_A @ D_B) @ S @ Cvir_A
    temp_JBbr = chain_gemm_einsums([Cocc_B, J_A, D_B], ["T", "N", "N"])
    JBbr = chain_gemm_einsums([temp_JBbr, S, Cvir_A], ["N", "N", "N"])
    JBbr.scale(-2.0)

    # Get your signs right Hesselmann!
    # Vbs = 2.0 * Cocc_B.T @ V_A @ Cvir_B
    Vbs = chain_gemm_einsums([Cocc_B, V_A, Cvir_B], ["T", "N", "N"])
    Vbs.scale(2.0)

    # Var = 2.0 * Cocc_A.T @ V_B @ Cvir_A
    Var = chain_gemm_einsums([Cocc_A, V_B, Cvir_A], ["T", "N", "N"])
    Var.scale(2.0)

    # VBas = -1.0 * (Cocc_A.T @ S @ D_B) @ V_A @ Cvir_B
    temp_VBas = chain_gemm_einsums([Cocc_A, S, D_B], ["T", "N", "N"])
    VBas = chain_gemm_einsums([temp_VBas, V_A, Cvir_B], ["N", "N", "N"])
    VBas.scale(-1.0)

    # VAbr = -1.0 * (Cocc_B.T @ S @ D_A) @ V_B @ Cvir_A
    temp_VAbr = chain_gemm_einsums([Cocc_B, S, D_A], ["T", "N", "N"])
    VAbr = chain_gemm_einsums([temp_VAbr, V_B, Cvir_A], ["N", "N", "N"])
    VAbr.scale(-1.0)

    # VRas = 1.0 * (Cocc_A.T @ V_B @ P_A) @ S @ Cvir_B
    temp_VRas = chain_gemm_einsums([Cocc_A, V_B, P_A], ["T", "N", "N"])
    VRas = chain_gemm_einsums([temp_VRas, S, Cvir_B], ["N", "N", "N"])

    # VSbr = 1.0 * (Cocc_B.T @ V_A @ P_B) @ S @ Cvir_A
    temp_VSbr = chain_gemm_einsums([Cocc_B, V_A, P_B], ["T", "N", "N"])
    VSbr = chain_gemm_einsums([temp_VSbr, S, Cvir_A], ["N", "N", "N"])

    # Sas = Cocc_A.T @ S @ Cvir_B
    Sas = chain_gemm_einsums([Cocc_A, S, Cvir_B], ["T", "N", "N"])

    # Sbr = Cocc_B.T @ S @ Cvir_A
    Sbr = chain_gemm_einsums([Cocc_B, S, Cvir_A], ["T", "N", "N"])

    # Qbr = Jbr + Kbr + KObr + JAbr + JBbr + VAbr + VSbr
    Qbr = Jbr.clone()
    Qbr.add(Kbr)
    Qbr.add(KObr)
    Qbr.add(JAbr)
    Qbr.add(JBbr)
    Qbr.add(VAbr)
    Qbr.add(VSbr)

    # Qas = Jas + Kas + KOas + JAas + JBas + VBas + VRas
    Qas = Jas.clone()
    Qas.add(Kas)
    Qas.add(KOas)
    Qas.add(JAas)
    Qas.add(JBas)
    Qas.add(VBas)
    Qas.add(VRas)

    # SBar = Cocc_A.T @ S @ D_B @ S @ Cvir_A
    SBar = chain_gemm_einsums([Cocc_A, S, D_B, S, Cvir_A], ["T", "N", "N", "N", "N"])

    # SAbs = Cocc_B.T @ S @ D_A @ S @ Cvir_B
    SAbs = chain_gemm_einsums([Cocc_B, S, D_A, S, Cvir_B], ["T", "N", "N", "N", "N"])

    # Qar = Jar + Var
    Qar = Jar.clone()
    Qar.add(Var)

    # Qbs = Jbs + Vbs
    Qbs = Jbs.clone()
    Qbs.add(Vbs)

    # => Integrals from DFHelper <= #

    # Build list of orbital space matrices for DF transformations
    # Order: Cocc_A, Cvir_A, Cocc_B, Cvir_B, Cr1, Cs1, Ca2, Cb2, Cr3, Cs3, Ca4, Cb4
    # Convert einsums RuntimeTensorD objects to core.Matrix objects for DFHelper
    # RuntimeTensorD supports buffer protocol, so np.asarray() can convert to numpy
    orbital_spaces = [
        core.Matrix.from_array(Cocc_A),  # 0: 'a'
        core.Matrix.from_array(Cvir_A),  # 1: 'r'
        core.Matrix.from_array(Cocc_B),  # 2: 'b'
        core.Matrix.from_array(Cvir_B),  # 3: 's'
        core.Matrix.from_array(Cr1),  # 4: 'r1'
        core.Matrix.from_array(Cs1),  # 5: 's1'
        core.Matrix.from_array(Ca2),  # 6: 'a2'
        core.Matrix.from_array(Cb2),  # 7: 'b2'
        core.Matrix.from_array(Cr3),  # 8: 'r3'
        core.Matrix.from_array(Cs3),  # 9: 's3'
        core.Matrix.from_array(Ca4),  # 10: 'a4'
        core.Matrix.from_array(Cb4),  # 11: 'b4'
    ]

    # Calculate total columns for memory allocation
    ncol = sum(mat.shape[1] for mat in orbital_spaces)
    # All should have same number of rows (AO basis)
    nrows = orbital_spaces[0].shape[0]

    # Initialize DFHelper
    aux_basis = dimer_wfn.get_basisset("DF_BASIS_SCF")
    dfh = core.DFHelper(dimer_basis, aux_basis)

    # Set memory: total available minus space needed for orbital matrices
    # Note: In C++, doubles_ is the total memory budget in doubles
    # Here we use a reasonable default or get from options if available
    memory_doubles = core.get_memory() // 8
    orbital_memory = nrows * ncol
    dfh.set_memory(memory_doubles - orbital_memory)
    # print set memory in GB
    core.print_out(
        f"    Setting DFHelper memory to {
            (memory_doubles - orbital_memory) * 8 / 1e9:.3f} GB\n"
    )

    dfh.set_method("DIRECT_iaQ")
    dfh.set_nthreads(core.get_num_threads())
    dfh.initialize()
    dfh.print_header()

    # Add orbital spaces
    dfh.add_space("a", orbital_spaces[0])  # Cocc_A
    dfh.add_space("r", orbital_spaces[1])  # Cvir_A
    dfh.add_space("b", orbital_spaces[2])  # Cocc_B
    dfh.add_space("s", orbital_spaces[3])  # Cvir_B
    dfh.add_space("r1", orbital_spaces[4])  # Cr1
    dfh.add_space("s1", orbital_spaces[5])  # Cs1
    dfh.add_space("a2", orbital_spaces[6])  # Ca2
    dfh.add_space("b2", orbital_spaces[7])  # Cb2
    dfh.add_space("r3", orbital_spaces[8])  # Cr3
    dfh.add_space("s3", orbital_spaces[9])  # Cs3
    dfh.add_space("a4", orbital_spaces[10])  # Ca4
    dfh.add_space("b4", orbital_spaces[11])  # Cb4

    # Add DF transformations
    # Format: (name, left_space, right_space) -> computes (left|right) integrals
    dfh.add_transformation("Aar", "r", "a")  # (r|a) virtuals_A x occupied_A
    dfh.add_transformation("Abs", "s", "b")  # (s|b) virtuals_B x occupied_B
    dfh.add_transformation("Bas", "s1", "a")  # (s1|a) Cs1 x occupied_A
    dfh.add_transformation("Bbr", "r1", "b")  # (r1|b) Cr1 x occupied_B
    dfh.add_transformation("Cas", "s", "a2")  # (s|a2) virtuals_B x Ca2
    dfh.add_transformation("Cbr", "r", "b2")  # (r|b2) virtuals_A x Cb2
    dfh.add_transformation("Dar", "r3", "a")  # (r3|a) Cr3 x occupied_A
    dfh.add_transformation("Dbs", "s3", "b")  # (s3|b) Cs3 x occupied_B
    dfh.add_transformation("Ear", "r", "a4")  # (r|a4) virtuals_A x Ca4
    dfh.add_transformation("Ebs", "s", "b4")  # (s|b4) virtuals_B x Cb4

    # TODO: Handle link orbital spaces for parallel/perpendicular coupling (lines 6950-7018)
    # For now, skip this and proceed with standard dispersion calculation

    # Perform DF transformations
    dfh.transform()

    # Clear spaces now that transformations are done
    dfh.clear_spaces()

    # => Memory blocking setup

    # Number of threads (single-threaded in Python)
    nT = 1

    # Calculate overhead for work arrays
    overhead = 0
    overhead += 5 * nT * na * nb  # Tab, Vab, T2ab, V2ab, Iab work arrays
    # For link orbitals with parperp, we'd need more, but we're skipping that
    overhead += (
        2 * na * ns + 2 * nb * nr + 2 * na * nr + 2 * nb * ns
    )  # S and Q matrices
    # E_disp20 and E_exch_disp20 thread work and final
    overhead += 2 * na * nb * (nT + 1)
    # sE_exch_disp20 thread work and final
    overhead += 1 * sna * snb * (nT + 1)
    overhead += 1 * (nA + nfa + na) * (nB + nfb + nb)  # Disp_AB
    overhead += 1 * (snA + snfa + sna) * (snB + snfb + snb)  # sDisp_AB
    overhead += 12 * nn * nn  # D, V, J, K, P, C matrices for A and B

    # Available memory for dispersion calculation
    total_memory = core.get_memory() // 8  # Convert bytes to doubles
    rem = total_memory - overhead

    core.print_out(
        f"    {total_memory} doubles - {overhead} overhead leaves {rem} for dispersion\n"
    )

    if rem < 0:
        raise Exception("Too little static memory for fdisp0")

    # Calculate cost per r or s virtual orbital
    # Each r needs: Aar, Bbr, Cbr, Dar (each is na x nQ or nb x nQ)
    cost_r = 2 * na * nQ + 2 * nb * nQ
    # Factor of 2 because we hold both r and s slices
    max_r_l = rem // (2 * cost_r)
    max_s_l = max_r_l
    max_r = min(max_r_l, nr)
    max_s = min(max_s_l, ns)

    if max_r < 1 or max_s < 1:
        raise Exception("Too little dynamic memory for fdisp0")

    nrblocks = (nr + max_r - 1) // max_r  # Ceiling division
    nsblocks = (ns + max_s - 1) // max_s

    core.print_out(
        f"    Processing a single (r,s) pair requires {cost_r * 2} doubles\n"
    )
    core.print_out(f"    {nr} values of r processed in {nrblocks} blocks of {max_r}\n")
    core.print_out(
        f"    {ns} values of s processed in {nsblocks} blocks of {max_s}\n\n"
    )

    # => Compute Far = Dar + Ear and Fbs = Dbs + Ebs
    # These represent combined D and E DF integrals that will be reused in the main loop

    # Add disk tensor for Far
    dfh.add_disk_tensor("Far", (nr, na, nQ))

    # Loop over r blocks to compute Far = Dar + Ear
    for rstart in range(0, nr, max_r):
        nrblock = min(max_r, nr - rstart)

        # Allocate matrices to hold the tensor slices
        Dar = core.Matrix("Dar block", nrblock * na, nQ)
        Ear = core.Matrix("Ear block", nrblock * na, nQ)

        # Fill Dar and Ear from disk tensors
        dfh.fill_tensor("Dar", Dar, [rstart, rstart + nrblock], [0, na], [0, nQ])
        dfh.fill_tensor("Ear", Ear, [rstart, rstart + nrblock], [0, na], [0, nQ])

        # Compute Far = Dar + Ear (element-wise addition)
        Dar.np[:, :] += Ear.np[:, :]

        # Write Far back to disk (Dar now contains Dar + Ear)
        dfh.write_disk_tensor("Far", Dar, (rstart, rstart + nrblock))

    # Add disk tensor for Fbs
    dfh.add_disk_tensor("Fbs", (ns, nb, nQ))

    # Loop over s blocks to compute Fbs = Dbs + Ebs
    for sstart in range(0, ns, max_s):
        nsblock = min(max_s, ns - sstart)

        # Allocate matrices to hold the tensor slices
        Dbs = core.Matrix("Dbs block", nsblock * nb, nQ)
        Ebs = core.Matrix("Ebs block", nsblock * nb, nQ)

        # Fill Dbs and Ebs from disk tensors
        dfh.fill_tensor("Dbs", Dbs, [sstart, sstart + nsblock], [0, nb], [0, nQ])
        dfh.fill_tensor("Ebs", Ebs, [sstart, sstart + nsblock], [0, nb], [0, nQ])

        # Compute Fbs = Dbs + Ebs (element-wise addition)
        Dbs.np[:, :] += Ebs.np[:, :]

        # Write Fbs back to disk (Dbs now contains Dbs + Ebs)
        dfh.write_disk_tensor("Fbs", Dbs, (sstart, sstart + nsblock))

    E_disp20_comp = core.Matrix("E_disp20", na, nb)
    E_exch_disp20_comp = core.Matrix("E_exch_disp20", na, nb)

    # => MO to LO Transformation
    Uaocc_A = cache["Uaocc0A"]
    Uaocc_B = cache["Uaocc0B"]
    UAp = Uaocc_A.np
    UBp = Uaocc_B.np

    # Orbital energies (already numpy arrays)
    # In the dispersion formula: indices a,b are occupied and r,s are virtual
    eap = eps_occ_A  # occupied energies for monomer A (index a)
    ebp = eps_occ_B  # occupied energies for monomer B (index b)
    erp = eps_vir_A  # virtual energies for monomer A (index r)
    esp = eps_vir_B  # virtual energies for monomer B (index s)

    # => Work arrays for inner loop
    Tab = core.Matrix("Tab", na, nb)
    Vab = core.Matrix("Vab", na, nb)
    T2ab = core.Matrix("T2ab", na, nb)
    V2ab = core.Matrix("V2ab", na, nb)
    Iab = core.Matrix("Iab", na, nb)

    # => Main r,s loop <= //
    # Allocate and fill r-block tensors
    Aar = core.Matrix("Aar block", nrblock * na, nQ)
    Far = core.Matrix("Far block", nrblock * na, nQ)
    Bbr = core.Matrix("Bbr block", nrblock * nb, nQ)
    Cbr = core.Matrix("Cbr block", nrblock * nb, nQ)

    # Allocate and fill s-block tensors
    Abs = core.Matrix("Abs block", nsblock * nb, nQ)
    Fbs = core.Matrix("Fbs block", nsblock * nb, nQ)
    Bas = core.Matrix("Bas block", nsblock * na, nQ)
    Cas = core.Matrix("Cas block", nsblock * na, nQ)
    core.timer_off("F-SAPT Disp Setup")

    core.timer_on("F-SAPT Disp Compute")
    for rstart in range(0, nr, max_r):
        nrblock = min(max_r, nr - rstart)

        dfh.fill_tensor("Aar", Aar, [rstart, rstart + nrblock], [0, na], [0, nQ])
        dfh.fill_tensor("Far", Far, [rstart, rstart + nrblock], [0, na], [0, nQ])
        dfh.fill_tensor("Bbr", Bbr, [rstart, rstart + nrblock], [0, nb], [0, nQ])
        dfh.fill_tensor("Cbr", Cbr, [rstart, rstart + nrblock], [0, nb], [0, nQ])

        # Get numpy pointers for r-block tensors and reshape to 3D
        # Tensors are stored as 2D with shape (nrblock * nX, nQ) and need to be (nrblock, nX, nQ)
        Aarp = Aar.np.reshape(nrblock, na, nQ)
        Farp = Far.np.reshape(nrblock, na, nQ)
        Bbrp = Bbr.np.reshape(nrblock, nb, nQ)
        Cbrp = Cbr.np.reshape(nrblock, nb, nQ)

        for sstart in range(0, ns, max_s):
            nsblock = min(max_s, ns - sstart)

            dfh.fill_tensor("Abs", Abs, [sstart, sstart + nsblock], [0, nb], [0, nQ])
            dfh.fill_tensor("Fbs", Fbs, [sstart, sstart + nsblock], [0, nb], [0, nQ])
            dfh.fill_tensor("Bas", Bas, [sstart, sstart + nsblock], [0, na], [0, nQ])
            dfh.fill_tensor("Cas", Cas, [sstart, sstart + nsblock], [0, na], [0, nQ])

            # Get numpy pointers for s-block tensors and reshape to 3D
            # Tensors are stored as 2D with shape (nsblock * nX, nQ) and need to be (nsblock, nX, nQ)
            Absp = Abs.np.reshape(nsblock, nb, nQ)
            Fbsp = Fbs.np.reshape(nsblock, nb, nQ)
            Basp = Bas.np.reshape(nsblock, na, nQ)
            Casp = Cas.np.reshape(nsblock, na, nQ)

            nrs = nrblock * nsblock

            # => RS inner loop <= //
            for rs in range(nrs):
                r = rs // nsblock
                s = rs % nsblock

                # Get pointers to work arrays and energy matrices
                Tabp = Tab.np
                Vabp = Vab.np
                T2abp = T2ab.np
                V2abp = V2ab.np
                Iabp = Iab.np
                E_disp20Tp = E_disp20_comp.np
                E_exch_disp20Tp = E_exch_disp20_comp.np

                # => Amplitudes, Disp20 <= //

                # Vab = Aar[r] @ Abs[s].T
                # Extract slices for r-th and s-th orbitals
                # Store these as we need them for Exch-Disp20 too
                Aar_r = Aarp[r, :, :]
                Abs_s = Absp[s, :, :]
                # Use einsum to match C++ DGEMM('N', 'T', ...) more closely
                np.einsum("aQ,bQ->ab", Aar_r, Abs_s, out=Vabp, optimize=True)

                # Compute amplitudes Tab[a,b] = Vab[a,b] / (ea + eb - er - es)
                for a in range(na):
                    for b in range(nb):
                        Tabp[a, b] = Vabp[a, b] / (
                            eap.np[a]
                            + ebp.np[b]
                            - erp.np[r + rstart]
                            - esp.np[s + sstart]
                        )

                # Transform to localized orbital basis
                # T2ab = UA.T @ Tab @ UB
                Iabp[:, :] = Tabp @ UBp
                T2abp[:, :] = UAp.T @ Iabp

                # V2ab = UA.T @ Vab @ UB
                Iabp[:, :] = Vabp @ UBp
                V2abp[:, :] = UAp.T @ Iabp

                # Accumulate Disp20
                for a in range(na):
                    for b in range(nb):
                        E_disp20Tp[a, b] += 4.0 * T2abp[a, b] * V2abp[a, b]

                # => Exch-Disp20 <= //

                # > Q1-Q3 < //
                # Vab = Bas[s] @ Bbr[r].T + Cas[s] @ Cbr[r].T + Aar[r] @ Fbs[s].T + Far[r] @ Abs[s].T
                # Extract slices for r-th and s-th orbitals
                Bas_s = Basp[s, :, :]
                Bbr_r = Bbrp[r, :, :]
                Cas_s = Casp[s, :, :]
                Cbr_r = Cbrp[r, :, :]
                Far_r = Farp[r, :, :]
                Fbs_s = Fbsp[s, :, :]

                Vabp[:, :] = Bas_s @ Bbr_r.T
                Vabp[:, :] += Cas_s @ Cbr_r.T
                Vabp[:, :] += Aar_r @ Fbs_s.T
                Vabp[:, :] += Far_r @ Abs_s.T

                # > V,J,K < //
                # Add outer product contributions using DGER equivalent
                # C_DGER(na, nb, 1.0, &Sasp[0][s + sstart], ns, &Qbrp[0][r + rstart], nr, Vabp[0], nb);
                Vabp[:, :] += np.outer(Sas.np[:, s + sstart], Qbr.np[:, r + rstart])

                # C_DGER(na, nb, 1.0, &Qasp[0][s + sstart], ns, &Sbrp[0][r + rstart], nr, Vabp[0], nb);
                Vabp[:, :] += np.outer(Qas.np[:, s + sstart], Sbr.np[:, r + rstart])

                # C_DGER(na, nb, 1.0, &Qarp[0][r + rstart], nr, &SAbsp[0][s + sstart], ns, Vabp[0], nb);
                Vabp[:, :] += np.outer(Qar.np[:, r + rstart], SAbs.np[:, s + sstart])

                # C_DGER(na, nb, 1.0, &SBarp[0][r + rstart], nr, &Qbsp[0][s + sstart], ns, Vabp[0], nb);
                Vabp[:, :] += np.outer(SBar.np[:, r + rstart], Qbs.np[:, s + sstart])

                # Transform to localized orbital basis
                Iabp[:, :] = Vabp @ UBp
                V2abp[:, :] = UAp.T @ Iabp

                # Accumulate ExchDisp20
                for a in range(na):
                    for b in range(nb):
                        E_exch_disp20Tp[a, b] -= 2.0 * T2abp[a, b] * V2abp[a, b]

    core.timer_off("F-SAPT Disp Compute")
    # => Accumulate thread results <= //
    E_disp20 = core.Matrix("E_disp20", nA + nfa + na1 + 1, nB + nfb + nb1 + 1)
    E_exch_disp20 = core.Matrix("E_exch_disp20", nA + nfa + na1 + 1, nB + nfb + nb1 + 1)

    # Single-threaded, so just use the first (and only) thread result
    # E_disp20.copy(E_disp20_comp)
    # E_exch_disp20.copy(E_exch_disp20_comp)
    for a in range(na):
        for b in range(nb):
            E_disp20.np[a + nfa + nA, b + nfb + nB] = E_disp20_comp.np[a, b]
            E_exch_disp20.np[a + nfa + nA, b + nfb + nB] = E_exch_disp20_comp.np[a, b]

    # => Populate cache['E'] matrix <= //
    # Store energy matrices and scalars
    # cache['E_DISP20'] = E_disp20
    # cache['E_EXCH_DISP20'] = E_exch_disp20
    # add E_disp20 and E_exch_disp20
    # Disp_AB = core.Matrix("DISP_AB", na, nb)
    Disp_AB = core.Matrix("Disp_AB", nA + nfa + na1 + 1, nB + nfb + nb1 + 1)
    Disp_AB.np[:, :] = E_disp20.np + E_exch_disp20.np
    cache["Disp_AB"] = Disp_AB
    # => Output printing <= //
    # {Elst10*1000:.8f} [mEh]
    Disp20 = np.sum(E_disp20.np)
    ExchDisp20 = np.sum(E_exch_disp20.np)

    cache["Exch-Disp20,u"] = ExchDisp20
    cache["Disp20,u"] = Disp20
    # if do_print:
    #     core.print_out(f"    Disp20              = {Disp20 * 1000:.8f} [mEh]\n")
    #     core.print_out(f"    Exch-Disp20         = {ExchDisp20 * 1000:.8f} [mEh]\n")
    #     core.print_out("\n")
    #     assert abs(scalars['Disp20,u'] - Disp20) < 1e-6, f"Disp20 scalar mismatch! {scalars['Disp20,u'] = } {Disp20 = }"
    #     assert abs(scalars['Exch-Disp20,u'] - ExchDisp20) < 1e-6, f"ExchDisp20 scalar mismatch!\nRef: {scalars['Exch-Disp20,u']:.4e}\nAct: {ExchDisp20:.4e}"
    return cache


def chain_gemm_einsums(
    tensors: list[core.Matrix],
    transposes: list[str] = None,
    prefactors_C: list[float] = None,
    prefactors_AB: list[float] = None,
    return_tensors: list[bool] = None,
) -> core.Matrix | list[core.Matrix]:
    """
    Computes a chain of einsum matrix multiplications

    Parameters
    ----------
    tensors : list[core.Matrix]
        List of tensors to be contracted.
    transposes : list[str], optional
        List of transpose operations for each tensor, where "N" means no transpose and "T" means transpose.
    prefactors_C : list[float], optional
        List of prefactors for the resulting tensors in the chain.
    prefactors_AB : list[float], optional
        List of prefactors for the tensors being multiplied in the chain.
    return_tensors : list[bool], optional
        List indicating which intermediate tensors should be returned. If None,
        only the final tensor is returned. Note that these are only
        intermediate tensors and final tensor; hence, the length of this list
        should be one less than the number of tensors.
    """
    # initialization "computed_tensors" with the first tensor of the chain
    computed_tensors = [tensors[0]]
    N = len(tensors)
    if transposes is None:
        transposes = ["N"] * N
    if prefactors_C is None:
        prefactors_C = [0.0] * (N - 1)
    if prefactors_AB is None:
        prefactors_AB = [1.0] * (N - 1)
    try:
        for i in range(len(tensors) - 1):
            A = computed_tensors[-1]
            B = tensors[i + 1]

            # For intermediate results (i > 0), always use 'N' for T1 since A is a computed intermediate
            T1 = transposes[i] if i == 0 else "N"
            T2 = transposes[i + 1]
            A_size = A.shape[0]
            if T1 == "T":
                A_size = A.shape[1]
            B_size = B.shape[1]
            if T2 == "T":
                B_size = B.shape[0]

            # Initialize output as psi4.core.Matrix with zeros
            C = core.Matrix(A_size, B_size)
            C.zero()
            # Use ein.core.gemm to write to C.np
            ein.core.gemm(T1, T2, prefactors_AB[i], A.np, B.np, prefactors_C[i], C.np)
            computed_tensors.append(C)
    except Exception as e:
        raise ValueError(
            f"Error in einsum_chain_gemm: {e}\n{i=}\n{A=}\n{B=}\n{T1=}\n{T2=}"
        )
    if return_tensors is None:
        return computed_tensors[-1]
    returned_tensors = []
    for i, r in enumerate(return_tensors):
        if r:
            returned_tensors.append(computed_tensors[i + 1])
    return returned_tensors


def exchange(cache: dict, jk: core.JK, do_print: bool = True) -> dict:
    r"""Compute the first-order exchange energy :math:`E^{(1)}_{\text{exch}}`.

    Evaluates both the :math:`S^2` approximation (Eq. 6) and the
    :math:`S^\infty` (Eq. 9) first-order exchange energies from
    Xie et al. (2022).

    The :math:`S^2` approximation is:

    .. math::

        E^{(1)}_{\text{exch}}(S^2) = -2(\mathbf{P}^A \mathbf{S} \mathbf{P}^B \mathbf{S} \mathbf{P}^{A,\text{vir}}) \cdot \boldsymbol{\omega}^B
          - 2(\mathbf{P}^B \mathbf{S} \mathbf{P}^A \mathbf{S} \mathbf{P}^{B,\text{vir}}) \cdot \boldsymbol{\omega}^A
          - 2(\mathbf{P}^{A,\text{vir}} \mathbf{S} \mathbf{P}^B) \cdot \mathbf{K}[\mathbf{P}^A \mathbf{S} \mathbf{P}^{B,\text{vir}}]

    The :math:`S^\infty` exchange energy uses the full inverse overlap metric (Eq. 9):

    .. math::

        E^{(1)}_{\text{exch}} = -2\mathbf{P}^A \cdot \mathbf{K}^B
          + 2\mathbf{T}^{AB} \cdot (\mathbf{h}^A + \mathbf{h}^B)
          + 2\mathbf{T}^{AA} \cdot \mathbf{h}^B
          + 2\mathbf{T}^{BB} \cdot \mathbf{h}^A
          + 2\mathbf{T}^{BB} \cdot \mathbf{W}^{AB}
          + 2\mathbf{T}^{AA} \cdot \mathbf{W}^{AB}
          + 2\mathbf{T}^{BB} \cdot \mathbf{W}^{AA}
          + 2\mathbf{T}^{AB} \cdot \mathbf{W}^{AB}

    where the intermediates are (Eq. 8, 10):

    .. math::

        \boldsymbol{\omega}^X = 2\mathbf{J}^X + \mathbf{V}^X

        \mathbf{h}^X = \mathbf{V}^X + 2\mathbf{J}^X - \mathbf{K}^X

    Parameters
    ----------
    cache : dict
        SAPT data cache from :func:`build_sapt_jk_cache`.
    jk : core.JK
        JK integral engine for computing Coulomb and exchange matrices.
    do_print : bool, optional
        Whether to print the result, by default True.

    Returns
    -------
    dict
        Dictionary with keys ``'Exch10(S^2)'`` and ``'Exch10'`` mapping to
        the :math:`S^2` and :math:`S^\infty` exchange energies, respectively.
    """

    if do_print:
        core.print_out("\n  ==> E10 Exchange Einsums <== \n\n")

    # Eq. 10: h^A = V^A + 2*J^A - K^A
    h_A = cache["V_A"].clone()
    ein.core.axpy(2.0, cache["J_A"].np, h_A.np)
    ein.core.axpy(-1.0, cache["K_A"].np, h_A.np)

    # Eq. 10: h^B = V^B + 2*J^B - K^B
    h_B = cache["V_B"].clone()
    ein.core.axpy(2.0, cache["J_B"].np, h_B.np)
    ein.core.axpy(-1.0, cache["K_B"].np, h_B.np)

    # Eq. 8: omega^A = V^A + 2*J^A
    w_A = cache["V_A"].clone()
    ein.core.axpy(2.0, cache["J_A"].np, w_A.np)

    # Eq. 8: omega^B = V^B + 2*J^B
    w_B = cache["V_B"].clone()
    ein.core.axpy(2.0, cache["J_B"].np, w_B.np)

    # Build inverse exchange metric
    nocc_A = cache["Cocc_A"].shape[1]
    nocc_B = cache["Cocc_B"].shape[1]
    SAB = chain_gemm_einsums(
        [cache["Cocc_A"], cache["S"], cache["Cocc_B"]],
        ["T", "N", "N"],
    )

    num_occ = nocc_A + nocc_B

    Sab = core.Matrix(num_occ, num_occ)
    Sab.np[:nocc_A, nocc_A:] = SAB.np
    Sab.np[nocc_A:, :nocc_A] = SAB.np.T
    Sab.np[np.diag_indices_from(Sab.np)] += 1
    Sab.power(-1.0, 1.0e-14)
    Sab.np[np.diag_indices_from(Sab.np)] -= 1.0

    Tmo_AA = core.Matrix.from_array(Sab.np[:nocc_A, :nocc_A])
    Tmo_BB = core.Matrix.from_array(Sab.np[nocc_A:, nocc_A:])
    Tmo_AB = core.Matrix.from_array(Sab.np[:nocc_A, nocc_A:])

    T_AA = chain_gemm_einsums(
        [cache["Cocc_A"], Tmo_AA, cache["Cocc_A"]], ["N", "N", "T"]
    )
    T_BB = chain_gemm_einsums(
        [cache["Cocc_B"], Tmo_BB, cache["Cocc_B"]], ["N", "N", "T"]
    )
    T_AB = chain_gemm_einsums(
        [cache["Cocc_A"], Tmo_AB, cache["Cocc_B"]], ["N", "N", "T"]
    )

    S = cache["S"]
    D_A = cache["D_A"]
    P_A = cache["P_A"]
    D_B = cache["D_B"]
    P_B = cache["P_B"]

    # Compute the J and K matrices
    jk.C_clear()

    jk.C_left_add(core.Matrix.from_array(cache["Cocc_A"]))
    jk.C_right_add(chain_gemm_einsums([cache["Cocc_A"], Tmo_AA]))

    jk.C_left_add(core.Matrix.from_array(cache["Cocc_B"]))
    jk.C_right_add(chain_gemm_einsums([cache["Cocc_A"], Tmo_AB]))

    jk.C_left_add(core.Matrix.from_array(cache["Cocc_A"]))
    jk.C_right_add(chain_gemm_einsums([P_B, S, cache["Cocc_A"]]))
    # This also works... you can choose to form the density-like matrix either
    # way..., just remember that the C_right_add has an adjoint (transpose, and switch matmul order)
    # jk.C_left_add(core.Matrix.from_array(einsum_chain_gemm([D_A, S, cache['Cvir_B']])))
    # jk.C_right_add(core.Matrix.from_array(cache['Cvir_B']))
    jk.compute()

    JT_A, JT_AB, Jij = jk.J()
    KT_A, KT_AB, Kij = jk.K()

    # Eq. 6: E^(1)_exch(S^2) — three-term S^2 exchange
    Exch_s2 = 0.0

    # Save some intermediate tensors to avoid recomputation in the next steps
    DA_S_DB_S_PA = chain_gemm_einsums([D_A, S, D_B, S, P_A])
    Exch_s2 -= 2.0 * ein.core.dot(w_B.np, DA_S_DB_S_PA.np)

    DB_S_DA_S_PB = chain_gemm_einsums([D_B, S, D_A, S, P_B])
    Exch_s2 -= 2.0 * ein.core.dot(w_A.np, DB_S_DA_S_PB.np)
    Exch_s2 -= 2.0 * ein.core.dot(Kij.np, chain_gemm_einsums([P_A, S, D_B]).np)

    if do_print:
        core.print_out(print_sapt_var("Exch10(S^2) ", Exch_s2, short=True))
        core.print_out("\n")

    # Eq. 9: E^(1)_exch(S^inf) — full inverse-overlap exchange
    Exch10 = 0.0
    Exch10 -= 2.0 * ein.core.dot(D_A.np, cache["K_B"].np)
    Exch10 += 2.0 * ein.core.dot(T_AA.np, h_B.np)
    Exch10 += 2.0 * ein.core.dot(T_BB.np, h_A.np)
    Exch10 += 2.0 * ein.core.dot(T_AB.np, h_A.np + h_B.np)
    Exch10 += 4.0 * ein.core.dot(T_BB.np, JT_AB.np - 0.5 * KT_AB.np)
    Exch10 += 4.0 * ein.core.dot(T_AA.np, JT_AB.np - 0.5 * KT_AB.np.T)
    Exch10 += 4.0 * ein.core.dot(T_BB.np, JT_A.np - 0.5 * KT_A.np)
    Exch10 += 4.0 * ein.core.dot(T_AB.np, JT_AB.np - 0.5 * KT_AB.np.T)

    if do_print:
        core.set_variable("Exch10", Exch10)
        core.print_out(print_sapt_var("Exch10", Exch10, short=True))
        core.print_out("\n")

    return {"Exch10(S^2)": Exch_s2, "Exch10": Exch10}


def induction(
    cache: dict,
    jk: core.JK,
    do_print: bool = True,
    maxiter: int = 12,
    conv: float = 1.0e-8,
    do_response: bool = True,
    Sinf: bool = False,
    sapt_jk_B: core.JK | None = None,
) -> dict:
    r"""Compute second-order induction and exchange-induction energies.

    Evaluates the uncoupled and (optionally) coupled second-order induction
    energy :math:`E^{(2)}_{\text{ind}}` and exchange-induction energy
    :math:`E^{(2)}_{\text{exch-ind}}` from Xie et al. (2022).

    The induction energy for monomer A polarized by B is (Eq. 14):

    .. math::

        E^{(2)}_{\text{ind}}(A \leftarrow B) = 2\mathbf{x}^A \cdot \tilde{\boldsymbol{\omega}}^B

    where the induction potential :math:`\tilde{\omega}^B` is given by (Eq. 16):

    .. math::

        \tilde{\boldsymbol{\omega}}^A = (\mathbf{C}^{B,\text{occ}})^\dagger \boldsymbol{\omega}^A \mathbf{C}^{B,\text{vir}}

    and the uncoupled response amplitudes are (Eq. 20):

    .. math::

        (x^A)^a_r = -(\tilde{\omega}^B)^a_r / (\epsilon_r - \epsilon_a)

    For monomer B polarized by A, the formulas are analogous with A and B swapped.

    Parameters
    ----------
    cache : dict
        SAPT data cache from :func:`build_sapt_jk_cache`.
    jk : core.JK
        JK integral engine for Coulomb and exchange matrices.
    do_print : bool, optional
        Whether to print results, by default True.
    maxiter : int, optional
        Maximum CPSCF iterations for coupled induction, by default 12.
    conv : float, optional
        Convergence threshold for CPSCF solver, by default 1.0e-8.
    do_response : bool, optional
        Whether to compute coupled (CPSCF) induction, by default True.
    Sinf : bool, optional
        Whether to include :math:`S^\infty` exchange-induction, by default False.
    sapt_jk_B : core.JK or None, optional
        Separate JK object for monomer B, by default None (uses same as A).

    Returns
    -------
    dict
        Dictionary containing induction energies with keys such as
        ``'Ind20,u (A<-B)'``, ``'Ind20,u (A->B)'``, ``'Ind20,u'``,
        ``'Exch-Ind20,u (A<-B)'``, ``'Exch-Ind20,u (A->B)'``,
        ``'Exch-Ind20,u'``, and coupled variants (``'Ind20,r'``, etc.)
        when ``do_response=True``.
    """

    if do_print:
        core.print_out("\n  ==> E20 Induction Einsums <== \n\n")

    # Build Induction and Exchange-Induction potentials
    S = cache["S"]

    D_A = cache["D_A"]
    V_A = cache["V_A"]

    J_A = cache["J_A"]
    K_A = cache["K_A"]

    D_B = cache["D_B"]
    V_B = cache["V_B"]
    J_B = cache["J_B"]
    K_B = cache["K_B"]

    K_O = cache["K_O"]
    J_O = cache["J_O"]

    # Set up matrix multiplication plans
    plan_matmul_tt = ein.core.compile_plan("ij", "ik", "kj")

    # Prepare JK calculations
    jk.C_clear()

    DB_S, DB_S_CA = chain_gemm_einsums(
        [D_B, S, cache["Cocc_A"]], return_tensors=[True, True]
    )
    jk.C_left_add(core.Matrix.from_array(DB_S_CA))
    jk.C_right_add(core.Matrix.from_array(cache["Cocc_A"]))

    jk.C_left_add(
        core.Matrix.from_array(chain_gemm_einsums([DB_S, D_A, S, cache["Cocc_B"]]))
    )
    jk.C_right_add(core.Matrix.from_array(cache["Cocc_B"]))

    DA_S, DA_S_DB_S_CA = chain_gemm_einsums(
        [D_A, S, D_B, S, cache["Cocc_A"]],
        return_tensors=[True, False, False, True],
    )
    jk.C_left_add(core.Matrix.from_array(DA_S_DB_S_CA))
    jk.C_right_add(core.Matrix.from_array(cache["Cocc_A"]))

    jk.compute()

    J_Ot, J_P_B, J_P_A = jk.J()
    K_Ot, K_P_B, K_P_A = jk.K()

    # Save for later usage in find()
    cache["J_P_A"] = J_P_A
    cache["J_P_B"] = J_P_B

    # Eq. 17: exchange-induction potential for A due to B
    EX_A = K_B.clone()
    EX_A.scale(-1.0)
    ein.core.axpy(-2.0, J_O.np, EX_A.np)
    ein.core.axpy(1.0, K_O.np, EX_A.np)
    ein.core.axpy(2.0, J_P_B.np, EX_A.np)

    # Apply all the axpy operations to EX_A
    S_DB, S_DB_VA, S_DB_VA_DB_S = chain_gemm_einsums(
        [S, D_B, V_A, D_B, S], return_tensors=[True, True, False, True]
    )
    S_DB_JA, S_DB_JA_DB_S = chain_gemm_einsums(
        [S_DB, J_A, D_B, S], return_tensors=[True, False, True]
    )
    S_DB_S_DA, S_DB_S_DA_VB = chain_gemm_einsums(
        [S_DB, S, D_A, V_B],
        return_tensors=[False, True, True],
    )
    ein.core.axpy(-1.0, S_DB_VA.np, EX_A.np)
    ein.core.axpy(-2.0, S_DB_JA.np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([S_DB, K_A]).np, EX_A.np)
    ein.core.axpy(1.0, S_DB_S_DA_VB.np, EX_A.np)
    ein.core.axpy(2.0, chain_gemm_einsums([S_DB_S_DA, J_B]).np, EX_A.np)
    ein.core.axpy(1.0, S_DB_VA_DB_S.np, EX_A.np)
    ein.core.axpy(2.0, S_DB_JA_DB_S.np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([S_DB, K_O], ["N", "T"]).np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([V_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(-2.0, chain_gemm_einsums([J_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([K_B, D_B, S]).np, EX_A.np)
    ein.core.axpy(1.0, chain_gemm_einsums([V_B, D_A, S, D_B, S]).np, EX_A.np)
    ein.core.axpy(2.0, chain_gemm_einsums([J_B, D_A, S, D_B, S]).np, EX_A.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([K_O, D_B, S]).np, EX_A.np)

    EX_A_MO_1 = chain_gemm_einsums(
        [cache["Cocc_A"], EX_A, cache["Cvir_A"]],
        ["T", "N", "N"],
    )
    mapA = {
        "S": S,
        "J_O": J_O,
        "K_O": K_O,
        "Cocc_A": cache["Cocc_A"],
        "Cvir_A": cache["Cvir_A"],
        "D_A": D_A,
        "V_A": V_A,
        "J_A": J_A,
        "K_A": K_A,
        "J_P_A": J_P_A,
        "Cocc_B": cache["Cocc_B"],
        "Cvir_B": cache["Cvir_B"],
        "D_B": D_B,
        "V_B": V_B,
        "J_B": J_B,
        "K_B": K_B,
        "J_P_B": J_P_B,
    }
    EX_A_MO = build_exch_ind_pot_AB(mapA)
    assert np.allclose(EX_A_MO, EX_A_MO_1), "EX_A_MO and EX_A_MO_1 do not match!"

    # Eq. 17: exchange-induction potential for B due to A
    EX_B = K_A.clone()
    EX_B.scale(-1.0)
    ein.core.axpy(-2.0, J_O.np, EX_B.np)
    ein.core.axpy(1.0, K_O.np, EX_B.np.T)
    ein.core.axpy(2.0, J_P_A.np, EX_B.np)
    cache["J_P_A"] = J_P_A
    cache["J_P_B"] = J_P_B

    S_DA, S_DA_VB, S_DA_VB_DA_S = chain_gemm_einsums(
        [S, D_A, V_B, D_A, S], return_tensors=[True, True, False, True]
    )
    S_DA_JB, S_DA_JB_DA_S = chain_gemm_einsums(
        [S_DA, J_B, D_A, S], return_tensors=[True, False, True]
    )
    S_DA_S_DB, S_DA_S_DB_VA = chain_gemm_einsums(
        [S_DA, S, D_B, V_A],
        return_tensors=[False, True, True],
    )

    # Apply all the axpy operations to EX_B
    ein.core.axpy(-1.0, S_DA_VB.np, EX_B.np)
    ein.core.axpy(-2.0, S_DA_JB.np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([S_DA, K_B]).np, EX_B.np)
    ein.core.axpy(1.0, S_DA_S_DB_VA.np, EX_B.np)
    ein.core.axpy(2.0, chain_gemm_einsums([S_DA_S_DB, J_A]).np, EX_B.np)
    ein.core.axpy(1.0, S_DA_VB_DA_S.np, EX_B.np)
    ein.core.axpy(2.0, S_DA_JB_DA_S.np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([S_DA, K_O]).np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([V_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(-2.0, chain_gemm_einsums([J_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([K_A, D_A, S]).np, EX_B.np)
    ein.core.axpy(1.0, chain_gemm_einsums([V_A, D_B, S, D_A, S]).np, EX_B.np)
    ein.core.axpy(2.0, chain_gemm_einsums([J_A, D_B, S, D_A, S]).np, EX_B.np)
    ein.core.axpy(-1.0, chain_gemm_einsums([K_O, D_A, S], ["T", "N", "N"]).np, EX_B.np)

    EX_B_MO_1 = chain_gemm_einsums(
        [cache["Cocc_B"], EX_B, cache["Cvir_B"]],
        ["T", "N", "N"],
    )
    EX_B_MO = build_exch_ind_pot_BA(mapA)
    assert np.allclose(EX_B_MO, EX_B_MO_1), "EX_B_MO and EX_B_MO_1 do not match!"

    # Eq. 8: omega^A = V^A + 2*J^A
    w_A = V_A.clone()
    w_A.name = "w_A"
    ein.core.axpy(2.0, J_A.np, w_A.np)

    # Eq. 8: omega^B = V^B + 2*J^B
    w_B = V_B.clone()
    w_B.name = "w_B"
    ein.core.axpy(2.0, J_B.np, w_B.np)

    w_B_MOA_1 = chain_gemm_einsums(
        [cache["Cocc_A"], w_B, cache["Cvir_A"]],
        ["T", "N", "N"],
    )
    w_A_MOB_1 = chain_gemm_einsums(
        [cache["Cocc_B"], w_A, cache["Cvir_B"]],
        ["T", "N", "N"],
    )

    # Eq. 16: induction potential omega_B in MO basis of A
    w_B_MOA = build_ind_pot(
        {
            "V_B": V_B,
            "J_B": J_B,
            "Cocc_A": cache["Cocc_A"],
            "Cvir_A": cache["Cvir_A"],
        }
    )
    w_B_MOA.name = "w_B_MOA"
    # Eq. 16: induction potential omega_A in MO basis of B
    w_A_MOB = build_ind_pot(
        {
            "V_B": V_A,
            "J_B": J_A,
            "Cocc_A": cache["Cocc_B"],
            "Cvir_A": cache["Cvir_B"],
        }
    )
    w_A_MOB.name = "w_A_MOB"
    assert np.allclose(w_B_MOA, w_B_MOA_1), "w_B_MOA and w_B_MOA_1 do not match!"
    assert np.allclose(w_A_MOB, w_A_MOB_1), "w_A_MOB and w_A_MOB_1 do not match!"

    # Do uncoupled induction calculations
    core.print_out("   => Uncoupled Induction <= \n\n")

    # Create uncoupled response vectors by element-wise division
    unc_x_B_MOA = w_B_MOA.clone()
    unc_x_A_MOB = w_A_MOB.clone()

    eps_occ_A = cache["eps_occ_A"]
    eps_vir_A = cache["eps_vir_A"]
    eps_occ_B = cache["eps_occ_B"]
    eps_vir_B = cache["eps_vir_B"]

    # Eq. 20
    for r in range(unc_x_B_MOA.shape[0]):
        for a in range(unc_x_B_MOA.shape[1]):
            unc_x_B_MOA.np[r, a] /= eps_occ_A.np[r] - eps_vir_A.np[a]

    # Eq. 20
    for r in range(unc_x_A_MOB.shape[0]):
        for a in range(unc_x_A_MOB.shape[1]):
            unc_x_A_MOB.np[r, a] /= eps_occ_B.np[r] - eps_vir_B.np[a]

    # Eq. 14: E^(2)_ind(A<-B) = 2 * x^A . omega_tilde^B
    unc_ind_ab = 2.0 * ein.core.dot(unc_x_B_MOA.np, w_B_MOA.np)
    unc_ind_ba = 2.0 * ein.core.dot(unc_x_A_MOB.np, w_A_MOB.np)
    unc_indexch_ab = 2.0 * ein.core.dot(unc_x_B_MOA.np, EX_A_MO.np)
    unc_indexch_ba = 2.0 * ein.core.dot(unc_x_A_MOB.np, EX_B_MO.np)

    ret = {}
    ret["Ind20,u (A<-B)"] = unc_ind_ab
    ret["Ind20,u (A->B)"] = unc_ind_ba
    ret["Ind20,u"] = unc_ind_ab + unc_ind_ba
    ret["Exch-Ind20,u (A<-B)"] = unc_indexch_ab
    ret["Exch-Ind20,u (A->B)"] = unc_indexch_ba
    ret["Exch-Ind20,u"] = unc_indexch_ba + unc_indexch_ab

    plist = [
        "Ind20,u (A<-B)",
        "Ind20,u (A->B)",
        "Ind20,u",
        "Exch-Ind20,u (A<-B)",
        "Exch-Ind20,u (A->B)",
        "Exch-Ind20,u",
    ]

    if do_print:
        for name in plist:
            core.print_out(print_sapt_var(name, ret[name], short=True))
            core.print_out("\n")

    # Exch-Ind without S^2 (Sinf calculations)
    if Sinf:
        nocc_A = cache["Cocc_A"].shape[1]
        nocc_B = cache["Cocc_B"].shape[1]
        SAB = core.triplet(
            cache["Cocc_A"], cache["S"], cache["Cocc_B"], True, False, False
        )
        num_occ = nocc_A + nocc_B

        Sab = core.Matrix(num_occ, num_occ)
        Sab.np[:nocc_A, nocc_A:] = SAB.np
        Sab.np[nocc_A:, :nocc_A] = SAB.np.T
        Sab.np[np.diag_indices_from(Sab.np)] += 1
        Sab.power(-1.0, 1.0e-14)

        Tmo_AA = core.Matrix.from_array(Sab.np[:nocc_A, :nocc_A])
        Tmo_BB = core.Matrix.from_array(Sab.np[nocc_A:, nocc_A:])
        Tmo_AB = core.Matrix.from_array(Sab.np[:nocc_A, nocc_A:])

        T_A = core.triplet(cache["Cocc_A"], Tmo_AA, cache["Cocc_A"], False, False, True)
        T_B = core.triplet(cache["Cocc_B"], Tmo_BB, cache["Cocc_B"], False, False, True)
        T_AB = core.triplet(
            cache["Cocc_A"], Tmo_AB, cache["Cocc_B"], False, False, True
        )

        sT_A = core.Matrix.chain_dot(
            cache["Cvir_A"],
            unc_x_B_MOA,
            Tmo_AA,
            cache["Cocc_A"],
            trans=[False, True, False, True],
        )
        sT_B = core.Matrix.chain_dot(
            cache["Cvir_B"],
            unc_x_A_MOB,
            Tmo_BB,
            cache["Cocc_B"],
            trans=[False, True, False, True],
        )
        sT_AB = core.Matrix.chain_dot(
            cache["Cvir_A"],
            unc_x_B_MOA,
            Tmo_AB,
            cache["Cocc_B"],
            trans=[False, True, False, True],
        )
        sT_BA = core.Matrix.chain_dot(
            cache["Cvir_B"],
            unc_x_A_MOB,
            Tmo_AB,
            cache["Cocc_A"],
            trans=[False, True, True, True],
        )

        jk.C_clear()

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_A"], Tmo_AA))
        jk.C_right_add(cache["Cocc_A"])

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_B"], Tmo_BB))
        jk.C_right_add(cache["Cocc_B"])

        jk.C_left_add(core.Matrix.chain_dot(cache["Cocc_A"], Tmo_AB))
        jk.C_right_add(cache["Cocc_B"])

        jk.compute()

        J_AA_inf, J_BB_inf, J_AB_inf = jk.J()
        K_AA_inf, K_BB_inf, K_AB_inf = jk.K()

        # A <- B
        EX_AA_inf = V_B.clone()
        EX_AA_inf.axpy(
            -1.00, core.Matrix.chain_dot(S, T_AB, V_B, trans=[False, True, False])
        )
        EX_AA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_B, V_B))
        EX_AA_inf.axpy(2.00, J_AB_inf)
        EX_AA_inf.axpy(
            -2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf, trans=[False, True, False])
        )
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AB_inf))
        EX_AA_inf.axpy(2.00, J_BB_inf)
        EX_AA_inf.axpy(
            -2.00, core.Matrix.chain_dot(S, T_AB, J_BB_inf, trans=[False, True, False])
        )
        EX_AA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_BB_inf))
        EX_AA_inf.axpy(-1.00, K_AB_inf.transpose())
        EX_AA_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, True, True])
        )
        EX_AA_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_B, K_AB_inf, trans=[False, False, True])
        )
        EX_AA_inf.axpy(-1.00, K_BB_inf)
        EX_AA_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_AB, K_BB_inf, trans=[False, True, False])
        )
        EX_AA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_BB_inf))

        EX_AB_inf = V_A.clone()
        EX_AB_inf.axpy(
            -1.00, core.Matrix.chain_dot(S, T_AB, V_A, trans=[False, True, False])
        )
        EX_AB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_B, V_A))
        EX_AB_inf.axpy(2.00, J_AA_inf)
        EX_AB_inf.axpy(
            -2.00, core.Matrix.chain_dot(S, T_AB, J_AA_inf, trans=[False, True, False])
        )
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AA_inf))
        EX_AB_inf.axpy(2.00, J_AB_inf)
        EX_AB_inf.axpy(
            -2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf, trans=[False, True, False])
        )
        EX_AB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_B, J_AB_inf))
        EX_AB_inf.axpy(-1.00, K_AA_inf)
        EX_AB_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_AB, K_AA_inf, trans=[False, True, False])
        )
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_AA_inf))
        EX_AB_inf.axpy(-1.00, K_AB_inf)
        EX_AB_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, True, False])
        )
        EX_AB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_B, K_AB_inf))

        # B <- A
        EX_BB_inf = V_A.clone()
        EX_BB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_A))
        EX_BB_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_A, V_A))
        EX_BB_inf.axpy(2.00, J_AB_inf)
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf))
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AB_inf))
        EX_BB_inf.axpy(2.00, J_AA_inf)
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AA_inf))
        EX_BB_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AA_inf))
        EX_BB_inf.axpy(-1.00, K_AB_inf)
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf))
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_AB_inf))
        EX_BB_inf.axpy(-1.00, K_AA_inf)
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_AA_inf))
        EX_BB_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_AA_inf))

        EX_BA_inf = V_B.clone()
        EX_BA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_AB, V_B))
        EX_BA_inf.axpy(-1.00, core.Matrix.chain_dot(S, T_A, V_B))
        EX_BA_inf.axpy(2.00, J_BB_inf)
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_BB_inf))
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_BB_inf))
        EX_BA_inf.axpy(2.00, J_AB_inf)
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_AB, J_AB_inf))
        EX_BA_inf.axpy(-2.00, core.Matrix.chain_dot(S, T_A, J_AB_inf))
        EX_BA_inf.axpy(-1.00, K_BB_inf)
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_AB, K_BB_inf))
        EX_BA_inf.axpy(1.00, core.Matrix.chain_dot(S, T_A, K_BB_inf))
        EX_BA_inf.axpy(-1.00, K_AB_inf.transpose())
        EX_BA_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_AB, K_AB_inf, trans=[False, False, True])
        )
        EX_BA_inf.axpy(
            1.00, core.Matrix.chain_dot(S, T_A, K_AB_inf, trans=[False, False, True])
        )

        unc_ind_ab_total = 2.0 * (
            sT_A.vector_dot(EX_AA_inf) + sT_AB.vector_dot(EX_AB_inf)
        )
        unc_ind_ba_total = 2.0 * (
            sT_B.vector_dot(EX_BB_inf) + sT_BA.vector_dot(EX_BA_inf)
        )
        unc_indexch_ab_inf = unc_ind_ab_total - unc_ind_ab
        unc_indexch_ba_inf = unc_ind_ba_total - unc_ind_ba

        ret["Exch-Ind20,u (A<-B) (S^inf)"] = unc_indexch_ab_inf
        ret["Exch-Ind20,u (A->B) (S^inf)"] = unc_indexch_ba_inf
        ret["Exch-Ind20,u (S^inf)"] = unc_indexch_ba_inf + unc_indexch_ab_inf

        if do_print:
            for name in plist[3:]:
                name = name + " (S^inf)"

                core.print_out(print_sapt_var(name, ret[name], short=True))
                core.print_out("\n")

    # Do coupled induction calculations
    if do_response:
        core.print_out("\n   => Coupled Induction <= \n\n")

        cphf_r_convergence = core.get_option("SAPT", "CPHF_R_CONVERGENCE")
        x_B_MOA, x_A_MOB = _sapt_cpscf_solve(
            cache,
            jk,
            w_B_MOA.np,
            w_A_MOB.np,
            maxiter,
            cphf_r_convergence,
            sapt_jk_B=sapt_jk_B,
        )
        x_B_MOA = core.Matrix.from_array(x_B_MOA)
        x_A_MOB = core.Matrix.from_array(x_A_MOB)

        ind_ab = 2.0 * ein.core.dot(x_B_MOA.np, w_B_MOA.np)
        ind_ba = 2.0 * ein.core.dot(x_A_MOB.np, w_A_MOB.np)
        indexch_ab = 2.0 * ein.core.dot(x_B_MOA.np, EX_A_MO.np)
        indexch_ba = 2.0 * ein.core.dot(x_A_MOB.np, EX_B_MO.np)

        ret["Ind20,r (A<-B)"] = ind_ab
        ret["Ind20,r (A->B)"] = ind_ba
        ret["Ind20,r"] = ind_ab + ind_ba
        ret["Exch-Ind20,r (A<-B)"] = indexch_ab
        ret["Exch-Ind20,r (A->B)"] = indexch_ba
        ret["Exch-Ind20,r"] = indexch_ba + indexch_ab

        if do_print:
            core.print_out("\n")
            for name in plist:
                name = name.replace(",u", ",r")
                core.print_out(print_sapt_var(name, ret[name], short=True))
                core.print_out("\n")

        # Exch-Ind without S^2
        if Sinf:
            cT_A = core.Matrix.chain_dot(
                cache["Cvir_A"],
                x_B_MOA,
                Tmo_AA,
                cache["Cocc_A"],
                trans=[False, True, False, True],
            )
            cT_B = core.Matrix.chain_dot(
                cache["Cvir_B"],
                x_A_MOB,
                Tmo_BB,
                cache["Cocc_B"],
                trans=[False, True, False, True],
            )
            cT_AB = core.Matrix.chain_dot(
                cache["Cvir_A"],
                x_B_MOA,
                Tmo_AB,
                cache["Cocc_B"],
                trans=[False, True, False, True],
            )
            cT_BA = core.Matrix.chain_dot(
                cache["Cvir_B"],
                x_A_MOB,
                Tmo_AB,
                cache["Cocc_A"],
                trans=[False, True, True, True],
            )

            ind_ab_total = 2.0 * (
                cT_A.vector_dot(EX_AA_inf) + cT_AB.vector_dot(EX_AB_inf)
            )
            ind_ba_total = 2.0 * (
                cT_B.vector_dot(EX_BB_inf) + cT_BA.vector_dot(EX_BA_inf)
            )
            indexch_ab_inf = ind_ab_total - ind_ab
            indexch_ba_inf = ind_ba_total - ind_ba

            ret["Exch-Ind20,r (A<-B) (S^inf)"] = indexch_ab_inf
            ret["Exch-Ind20,r (A->B) (S^inf)"] = indexch_ba_inf
            ret["Exch-Ind20,r (S^inf)"] = indexch_ba_inf + indexch_ab_inf

            if do_print:
                for name in plist[3:]:
                    name = name.replace(",u", ",r") + " (S^inf)"

                    core.print_out(print_sapt_var(name, ret[name], short=True))
                    core.print_out("\n")

    return ret


def _sapt_cpscf_solve(
    cache: dict,
    jk: core.JK,
    rhsA: np.ndarray,
    rhsB: np.ndarray,
    maxiter: int,
    conv: float,
    sapt_jk_B: core.JK | None = None,
) -> list:
    r"""Solve the coupled-perturbed SCF (CPSCF) equations for SAPT induction.

    Implements the coupled-perturbed Kohn-Sham (CPKS) or Hartree-Fock (CPHF)
    equations using a conjugate-gradient solver. The CPSCF response
    vectors :math:`\mathbf{x}^A` and :math:`\mathbf{x}^B` satisfy (Eq. 21-26
    of Xie et al. 2022):

    .. math::

        \mathbf{H}^{(1)} \mathbf{x}^A = \mathbf{\omega}^{B}

    and
    .. math::

        \mathbf{H}^{(1)} \mathbf{x}^B = \mathbf{\omega}^{A}

    where :math:`\mathbf{H}^{(1)}` includes exchange-correlation
    and exact exchange contributions.

    Parameters
    ----------
    cache : dict
        SAPT data cache containing wavefunctions and orbital energies.
    jk : core.JK
        JK integral engine for monomer A.
    rhsA : np.ndarray
        Right-hand side vector for monomer A response (:math:`-\tilde{\omega}^B`).
    rhsB : np.ndarray
        Right-hand side vector for monomer B response (:math:`-\tilde{\omega}^A`).
    maxiter : int
        Maximum number of CPSCF iterations.
    conv : float
        Convergence threshold (relative residual norm).
    sapt_jk_B : core.JK or None, optional
        Separate JK object for monomer B, by default None (uses same as A).

    Returns
    -------
    list
        Converged response vectors ``[x_A, x_B]`` as numpy arrays.
    """

    cache["wfn_A"].set_jk(jk)
    if sapt_jk_B:
        cache["wfn_B"].set_jk(sapt_jk_B)
    else:
        cache["wfn_B"].set_jk(jk)

    def setup_P_X(eps_occ, eps_vir, name="P_X"):
        P_X = ein.utils.tensor_factory(
            name, [eps_occ.shape[0], eps_vir.shape[0]], np.float64, "einsums"
        )

        ones_occ = ein.utils.tensor_factory(
            "ones_occ", [eps_occ.shape[0]], np.float64, "einsums"
        )
        ones_vir = ein.utils.tensor_factory(
            "ones_vir", [eps_vir.shape[0]], np.float64, "einsums"
        )
        ones_occ.set_all(1.0)
        ones_vir.set_all(1.0)
        plan_outer = ein.core.compile_plan("ia", "i", "a")
        plan_outer.execute(0.0, P_X, 1.0, eps_occ.np, ones_vir)
        eps_vir_2D = ein.utils.tensor_factory(
            "eps_vir_2D", [eps_occ.shape[0], eps_vir.shape[0]], np.float64, "einsums"
        )
        plan_outer.execute(0.0, eps_vir_2D, 1.0, ones_occ, eps_vir.np)
        ein.core.axpy(-1.0, eps_vir_2D, P_X)
        return P_X

    # Make a preconditioner function
    P_A = setup_P_X(cache["eps_occ_A"], cache["eps_vir_A"])
    P_B = setup_P_X(cache["eps_occ_B"], cache["eps_vir_B"])

    # Preconditioner function
    def apply_precon(x_vec, act_mask):
        if act_mask[0]:
            pA = x_vec[0].copy()
            pA /= P_A
        else:
            pA = False

        if act_mask[1]:
            pB = x_vec[1].copy()
            pB /= P_B
        else:
            pB = False
        return [pA, pB]

    # Hx function
    def hessian_vec(x_vec, act_mask):
        # TODO: to convert to einsums fully here, would need to re-write
        # cphf_HX, onel_Hx, and twoel_Hx functions in libscf_solver/uhf.cc
        if act_mask[0]:
            xA = cache["wfn_A"].cphf_Hx([core.Matrix.from_array(x_vec[0])])[0].np
        else:
            xA = False

        if act_mask[1]:
            xB = cache["wfn_B"].cphf_Hx([core.Matrix.from_array(x_vec[1])])[0].np
        else:
            xB = False

        return [xA, xB]

    # Manipulate the printing
    sep_size = 51
    core.print_out("   " + ("-" * sep_size) + "\n")
    core.print_out("   " + "SAPT Coupled Induction Solver".center(sep_size) + "\n")
    core.print_out("   " + ("-" * sep_size) + "\n")
    core.print_out("    Maxiter             = %11d\n" % maxiter)
    core.print_out("    Convergence         = %11.3E\n" % conv)
    core.print_out("   " + ("-" * sep_size) + "\n")

    tstart = time.time()
    core.print_out(
        "     %4s %12s     %12s     %9s\n" % ("Iter", "(A<-B)", "(B->A)", "Time [s]")
    )
    core.print_out("   " + ("-" * sep_size) + "\n")

    start_resid = [ein.core.dot(rhsA, rhsA), ein.core.dot(rhsB, rhsB)]

    def pfunc(niter, x_vec, r_vec):
        if niter == 0:
            niter = "Guess"
        else:
            niter = "%5d" % niter
        # Compute IndAB
        valA = (ein.core.dot(r_vec[0], r_vec[0]) / start_resid[0]) ** 0.5
        if valA < conv:
            cA = "*"
        else:
            cA = " "

        # Compute IndBA
        valB = (ein.core.dot(r_vec[1], r_vec[1]) / start_resid[1]) ** 0.5
        if valB < conv:
            cB = "*"
        else:
            cB = " "

        core.print_out(
            "    %5s %15.6e%1s %15.6e%1s %9d\n"
            % (niter, valA, cA, valB, cB, time.time() - tstart)
        )
        return [valA, valB]

    # Compute the solver
    vecs, resid = solvers.cg_solver_ein(
        [rhsA, rhsB],
        hessian_vec,
        apply_precon,
        maxiter=maxiter,
        rcond=conv,
        printlvl=0,
        printer=pfunc,
    )
    core.print_out("   " + ("-" * sep_size) + "\n")

    return vecs
