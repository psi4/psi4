"""
SAD (Superposition of Atomic Densities) initial guess.

The atomic SCF is a single-atom restricted (spin-averaged) calculation,
done one atom at a time. The algorithm:

1. Build S, T, V (and ECP if present) integrals via Psi4's MintsHelper.
2. X = S^(-1/2) with small-eigenvalue truncation.
3. SCF iteration with degeneracy-aware aufbau + optimal damping:
     - Diagonalize F in the orthonormal basis: F' = X^T F X.
     - Cluster eigenvalues that are degenerate to within a tolerance and
       distribute electrons equally within each cluster.
     - Build the alpha-density Da from the occupations.
     - Use Psi4's JK builder to get J(Da) and K(Da); add Vxc if a DFT
       functional was requested.
     - F_alpha = H + 2 J - alpha_x K + Vxc.
     - ODA line search: pick lambda in [0, 1] minimizing
       E((1-lambda) D + lambda D_tent). Two JK builds per iteration but
       guarantees monotone descent.

This works for both spherical and Cartesian basis sets; we cluster on
orbital energies, not on angular-momentum labels.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from psi4 import core


# Tolerance for grouping orbital energies into degenerate clusters
# (Hartree). Loose enough to accept accidental Cartesian splittings.
_DEGEN_TOL = 1.0e-4


# =====================================================================
#  Configuration access
# =====================================================================


def _opt_str(key: str) -> str:
    return core.get_option("SCF", key)


def _opt_int(key: str) -> int:
    return core.get_option("SCF", key)


def _opt_double(key: str) -> float:
    return core.get_option("SCF", key)


def _print_(level: int, msg: str) -> None:
    if _opt_int("SAD_PRINT") >= level:
        core.print_out(msg + "\n")


# =====================================================================
#  Linear-algebra helpers
# =====================================================================


def _orthogonalizer(S: np.ndarray, tol: float = 1.0e-10) -> np.ndarray:
    """Build X = S^(-1/2) with small-eigenvalue truncation."""
    sval, svec = np.linalg.eigh(S)
    keep = sval > tol
    return svec[:, keep] @ np.diag(1.0 / np.sqrt(sval[keep]))


def _diag_in_ortho_basis(F: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Solve F C = S C E by diagonalizing X^T F X (X = S^(-1/2)).
    Returns (C, E) with C in the AO basis and E sorted ascending."""
    Fp = X.T @ F @ X
    E, Cp = np.linalg.eigh(Fp)
    C = X @ Cp
    return C, E


def _degeneracy_clusters(eps: np.ndarray, tol: float) -> List[List[int]]:
    """Cluster sorted eigenvalues into degeneracy groups (within tol)."""
    clusters: List[List[int]] = []
    for i, e in enumerate(eps):
        if clusters and (e - eps[clusters[-1][0]]) < tol:
            clusters[-1].append(i)
        else:
            clusters.append([i])
    return clusters


def _aufbau_with_degeneracy(
    eps: np.ndarray, n_per_spin: float, tol: float = _DEGEN_TOL
) -> Dict[int, float]:
    """Fill aufbau-style across orbitals grouped by degeneracy.

    Each cluster of size k can hold k electrons per spin (one per spatial
    orbital per spin, closed-shell limit). The partial-filling cluster
    gets the remaining electrons split equally across its members.

    Returns a dict orbital_index -> occ_per_spin in [0, 1].
    """
    occ: Dict[int, float] = {}
    remaining = float(n_per_spin)
    eps_tol = 1.0e-12
    for cluster in _degeneracy_clusters(eps, tol):
        capacity = float(len(cluster))
        if remaining >= capacity - eps_tol:
            for i in cluster:
                occ[i] = 1.0
            remaining -= capacity
        elif remaining > eps_tol:
            share = remaining / capacity
            for i in cluster:
                occ[i] = share
            remaining = 0.0
        else:
            break
    return occ


def _build_density(C: np.ndarray, occ: Dict[int, float]) -> np.ndarray:
    """Alpha-density Da = sum_i occ_per_spin[i] * |C_i><C_i| (per spin)."""
    nbf = C.shape[0]
    Da = np.zeros((nbf, nbf))
    for i, o in occ.items():
        if o > 0:
            Da += o * np.outer(C[:, i], C[:, i])
    return Da


def _extract_occ(
    C: np.ndarray, E: np.ndarray, occ: Dict[int, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """Pick out the occupied AO orbitals and their energies in order."""
    idx = sorted(i for i, o in occ.items() if o > 0)
    if not idx:
        return np.zeros((C.shape[0], 0)), np.zeros(0)
    return C[:, idx], E[np.array(idx)]


# =====================================================================
#  DFT helpers
# =====================================================================


def _maybe_build_functional(name: str):
    if name.upper() == "HF":
        return None, 1.0
    from psi4.driver.procrouting.dft import build_superfunctional

    sup, _ = build_superfunctional(name, restricted=True)
    return sup, sup.x_alpha()


def _maybe_build_vbase(sup, basis: core.BasisSet):
    if sup is None:
        return None
    vb = core.VBase.build(basis, sup, "RV")
    vb.initialize()
    return vb


def _build_jk(basis: core.BasisSet, fit_basis):
    use_df = _opt_str("SAD_SCF_TYPE") in ("DF", "MEM_DF", "DISK_DF")
    aux = fit_basis if (use_df and fit_basis is not None) else core.BasisSet.zero_ao_basis_set()
    # Use the full ctor with explicit do_wK and memory so the JK factory picks the
    # right backend (MemDFJK / DirectJK / ...) instead of falling back to DiskDFJK
    # with the "simple constructor" warning.
    doubles = int(core.get_memory() / 2 / 8)
    jk = core.JK.build_JK(basis, aux, False, doubles)
    jk.initialize()
    return jk


# =====================================================================
#  Per-atom SCF
# =====================================================================


def _hermite_cubic_min(f0: float, g0: float, f1: float, g1: float) -> Tuple[float, float]:
    """Minimize the cubic Hermite interpolant on [0, 1] given function value
    and derivative at the two endpoints. Returns (lambda*, f(lambda*))."""
    d = g0 + g1 + 2.0 * (f0 - f1)
    c = 3.0 * (f1 - f0) - 2.0 * g0 - g1
    b = g0
    candidates = [0.0, 1.0]
    if abs(d) > 1e-14:
        disc = c * c - 3.0 * b * d
        if disc >= 0:
            sd = float(np.sqrt(disc))
            for sign in (-1.0, 1.0):
                root = (-c + sign * sd) / (3.0 * d)
                if 0.0 <= root <= 1.0:
                    candidates.append(root)
    elif abs(c) > 1e-14:
        root = -b / (2.0 * c)
        if 0.0 <= root <= 1.0:
            candidates.append(root)

    def feval(x: float) -> float:
        return f0 + b * x + c * x * x + d * x * x * x

    best = min(candidates, key=feval)
    return best, feval(best)


def _atomic_scf(
    atomic_basis: core.BasisSet, fit_basis, functional_name: str
) -> Dict[str, object]:
    """Single-atom restricted spin-averaged SCF with degeneracy-aware
    aufbau and ODA. Works with both spherical and Cartesian shells.

    Convergence uses a cubic Hermite line search between the current
    density Da and the aufbau-tentative density Dat, using function
    values and gradients (tr[(Dat - Da) F]) at both endpoints. The
    cubic is appropriate because the energy is genuinely quadratic in
    D for HF and weakly higher-order under DFT.

    The Fock matrix is built density-form (JK is driven by a Cholesky
    factor of Da, so the mixed ODA density is handled exactly).
    """
    nbf = atomic_basis.nbf()
    mol = atomic_basis.molecule()
    Z = int(round(mol.Z(0)))
    n_per_spin = 0.5 * Z

    mints = core.MintsHelper(atomic_basis)
    S = np.asarray(mints.ao_overlap())
    H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())
    if atomic_basis.has_ECP():
        H = H + np.asarray(mints.ao_ecp())

    X = _orthogonalizer(S)

    sup, alpha_x = _maybe_build_functional(functional_name)
    vbase = _maybe_build_vbase(sup, atomic_basis)
    jk = _build_jk(atomic_basis, fit_basis)

    def fock_and_energy(Da: np.ndarray) -> Tuple[np.ndarray, float]:
        """Build F_alpha and the total energy at the given alpha density.

        JK is driven by a Cholesky factor of Da, which lets us evaluate
        F and E at the ODA-mixed density even though it doesn't
        factorize as C @ diag(occ) @ C.T."""
        Da_mat = core.Matrix.from_array(Da)
        L = Da_mat.partial_cholesky_factorize(1.0e-12)
        jk.C_clear()
        jk.C_left_add(L)
        jk.compute()
        J = np.asarray(jk.J()[0])
        K = np.asarray(jk.K()[0])

        F = H + 2.0 * J - alpha_x * K
        E_xc = 0.0
        if vbase is not None:
            vbase.set_D([Da_mat])
            V_xc = core.Matrix("V_xc", nbf, nbf)
            vbase.compute_V([V_xc])
            F = F + np.asarray(V_xc)
            E_xc = vbase.quadrature_values().get("FUNCTIONAL", 0.0)

        D_total = 2.0 * Da
        # E = D.H + D.J(Da) - alpha_x Da.K(Da) + Exc  (restricted closed-shell)
        E = (
            float(np.einsum("ij,ij->", D_total, H))
            + float(np.einsum("ij,ij->", D_total, J))
            - alpha_x * float(np.einsum("ij,ij->", Da, K))
            + E_xc
        )
        return F, E

    # Initial guess: aufbau over the core-Hamiltonian eigenstates.
    C, eps = _diag_in_ortho_basis(H, X)
    occ = _aufbau_with_degeneracy(eps, n_per_spin)
    Da = _build_density(C, occ)
    F, E = fock_and_energy(Da)

    E_tol = _opt_double("SAD_E_CONVERGENCE")
    D_tol = _opt_double("SAD_D_CONVERGENCE")
    max_iter = _opt_int("SAD_MAXITER")
    converged = False

    for it in range(1, max_iter + 1):
        # Aufbau tentative density from the current Fock matrix.
        C_t, eps_t = _diag_in_ortho_basis(F, X)
        occ_t = _aufbau_with_degeneracy(eps_t, n_per_spin)
        Da_tent = _build_density(C_t, occ_t)
        F_tent, E_tent = fock_and_energy(Da_tent)

        # Cubic Hermite line search on E((1-lam)*Da + lam*Da_tent).
        dDa = Da_tent - Da
        # Total-density gradients: d E / d (Da_total scaling). For
        # closed-shell, total D = 2 Da, and dE/d Da_total = F_alpha, so
        # tr[ΔD F] = 2 tr[ΔDa F_alpha].
        g0 = 2.0 * float(np.einsum("ij,ij->", dDa, F))
        g1 = 2.0 * float(np.einsum("ij,ij->", dDa, F_tent))
        lam, _ = _hermite_cubic_min(E, g0, E_tent, g1)

        Da_new = (1.0 - lam) * Da + lam * Da_tent
        F_new, E_new = fock_and_energy(Da_new)

        dE = abs(E_new - E)
        dD = float(np.linalg.norm(Da_new - Da))
        _print_(
            2,
            f"  @Atomic SCF iter {it:3d}  E = {E_new:20.14f}"
            f"   dE = {dE:14.6e}   |dD| = {dD:14.6e}   lambda = {lam:6.4f}",
        )
        if it > 1 and dE < E_tol and dD < D_tol:
            converged = True

        Da = Da_new
        F = F_new
        E = E_new

        if converged:
            break
    else:
        core.print_out(
            f"\n  WARNING: SAD atomic SCF did not converge in {max_iter} iterations "
            f"for Z={Z}.\n"
        )

    # Final orbital expansion: diagonalize the converged F, take the
    # occupied orbitals (with their fractional occupations) for the
    # Huckel collection.
    C_final, eps_final = _diag_in_ortho_basis(F, X)
    occ_final = _aufbau_with_degeneracy(eps_final, n_per_spin)
    C_occ, E_occ = _extract_occ(C_final, eps_final, occ_final)

    if vbase is not None:
        vbase.finalize()
    jk.finalize()

    D_total = 2.0 * Da  # spin-averaged closed-shell
    return {
        "D": D_total,
        "C_occ": C_occ,
        "E_occ": E_occ,
        "n_electrons": float(np.einsum("ij,ij->", D_total, S)),
        "converged": converged,
    }


# =====================================================================
#  Unique-atom dedup
# =====================================================================


def _identify_unique_atoms(
    molecule: core.Molecule, atomic_bases: Sequence[core.BasisSet]
) -> Tuple[List[int], List[int], int]:
    """Match symmetry-equivalent atoms (same Z, identical atomic basis)."""
    n = molecule.natom()
    unique_of = list(range(n))
    for l in range(n - 1):
        for m in range(l + 1, n):
            if unique_of[m] != m:
                continue
            if molecule.Z(l) != molecule.Z(m):
                continue
            b_l, b_m = atomic_bases[l], atomic_bases[m]
            if (
                b_l.nbf() != b_m.nbf()
                or b_l.nshell() != b_m.nshell()
                or b_l.nprimitive() != b_m.nprimitive()
                or b_l.max_am() != b_m.max_am()
                or b_l.max_nprimitive() != b_m.max_nprimitive()
                or b_l.has_puream() != b_m.has_puream()
                or b_l.n_ecp_core() != b_m.n_ecp_core()
            ):
                continue
            unique_of[m] = l

    rank_of = [0] * n
    n_unique = 0
    for A in range(n):
        if unique_of[A] == A:
            rank_of[A] = n_unique
            n_unique += 1
        else:
            rank_of[A] = rank_of[unique_of[A]]
    return unique_of, rank_of, n_unique


# =====================================================================
#  Top-level entry points
# =====================================================================


def compute_sad_guess(
    basis: core.BasisSet,
    atomic_bases: Sequence[core.BasisSet],
    atomic_fit_bases=None,
    AO2SO: "core.Matrix | None" = None,
) -> Dict[str, object]:
    """Build the molecular SAD guess.

    Returns a dict with Da, Db (core.Matrix in the SO basis), Ca, Cb
    (core.Matrix from the partial Cholesky of Da), HuckelC, HuckelE
    (in the AO basis), and n_electrons (float)."""
    molecule = basis.molecule()
    functional_name = _opt_str("SAD_FUNCTIONAL")

    n_atoms = molecule.natom()
    z_per_atom = [int(round(molecule.Z(A))) for A in range(n_atoms)]
    unique_of, rank_of, _n_unique = _identify_unique_atoms(molecule, atomic_bases)

    _print_(1, f"  Performing Atomic SCF Computations  (functional = {functional_name})")
    per_atom: Dict[int, Dict[str, object]] = {}
    for A in range(n_atoms):
        if unique_of[A] != A:
            continue
        Z = z_per_atom[A]
        if Z == 0:
            continue
        nbf_a = atomic_bases[A].nbf()
        if Z > 2 * nbf_a:
            raise RuntimeError(
                f"SAD: atom {molecule.symbol(A)} has more electrons than basis "
                f"functions ({Z} > 2*{nbf_a})."
            )
        _print_(1, f"\n  Atom {A}  Z = {Z}  nbf = {nbf_a}")
        fit = atomic_fit_bases[A] if atomic_fit_bases is not None else None
        per_atom[rank_of[A]] = _atomic_scf(atomic_bases[A], fit, functional_name)

    # Assemble molecular DAO (block-diagonal by atom) and Huckel orbital collection.
    nbf_mol = basis.nbf()
    DAO = np.zeros((nbf_mol, nbf_mol))
    huckel_cols: List[np.ndarray] = []
    huckel_es: List[float] = []
    offset = 0
    for A in range(n_atoms):
        nbf_a = atomic_bases[A].nbf()
        if z_per_atom[A] > 0:
            res = per_atom[rank_of[A]]
            D_atom = res["D"]
            # Half-density goes to Da_ = Db_; total D = 2 * DAO.
            DAO[offset:offset + nbf_a, offset:offset + nbf_a] += 0.5 * D_atom

            C_atom = res["C_occ"]
            E_atom = res["E_occ"]
            for c in range(C_atom.shape[1]):
                col = np.zeros(nbf_mol)
                col[offset:offset + nbf_a] = C_atom[:, c]
                huckel_cols.append(col)
                huckel_es.append(float(E_atom[c]))
        offset += nbf_a

    HuckelC_np = (
        np.column_stack(huckel_cols) if huckel_cols else np.zeros((nbf_mol, 0))
    )
    HuckelE_np = np.array(huckel_es) if huckel_es else np.zeros(0)

    if AO2SO is None:
        # Fall back to a C1-only identity transform when the caller didn't
        # pass a symmetry transform.
        nbf_dim = core.Dimension([basis.nbf()])
        AO2SO = core.Matrix("AO2SO (C1)", nbf_dim, nbf_dim)
        for i in range(basis.nbf()):
            AO2SO.set(0, i, i, 1.0)
    Da_mat = core.Matrix("Da SAD", AO2SO.coldim(), AO2SO.coldim())
    DAO_mat = core.Matrix.from_array(DAO)
    Da_mat.apply_symmetry(DAO_mat, AO2SO)
    Db_mat = Da_mat  # closed-shell SAD

    tol = _opt_double("SAD_CHOL_TOLERANCE")
    Ca_mat = Da_mat.partial_cholesky_factorize(tol)
    Ca_mat.name = "Ca SAD"
    Cb_mat = Ca_mat

    n_elec_total = sum(
        float(per_atom[rank_of[A]]["n_electrons"])
        for A in range(n_atoms)
        if z_per_atom[A] > 0
    )
    _print_(1, f"\n  Number of electrons in SAD guess:  {n_elec_total:12.6f}\n")

    return {
        "Da": Da_mat,
        "Db": Db_mat,
        "Ca": Ca_mat,
        "Cb": Cb_mat,
        "HuckelC": core.Matrix.from_array(HuckelC_np, name="C_Huckel (MINAO)"),
        "HuckelE": core.Vector.from_array(HuckelE_np, name="E_Huckel (MINAO)"),
        "n_electrons": n_elec_total,
    }


def populate_sad_guess(
    wfn,
    atomic_bases: Sequence[core.BasisSet],
    atomic_fit_bases,
    natorb: bool,
) -> None:
    """Populate an HF wavefunction's Da/Db/Ca/Cb/occupations from a SAD guess.

    Called from C++ HF::guess(). When ``natorb`` is True, the SAD density is
    diagonalized to give natural orbitals (doi:10.1021/acs.jctc.8b01089);
    otherwise the SAD density's partial Cholesky factor is taken as the
    orbital guess and the alpha/beta occupations are set to its irrep
    dimensions. For ROHF and CUHF the total density ``Dt_`` (and ROHF's
    ``Ct_``) are also updated.
    """
    sad = compute_sad_guess(wfn.basisset(), atomic_bases, atomic_fit_bases, AO2SO=wfn.aotoso())

    if natorb:
        Dhelp = sad["Da"].clone()
        Dhelp.scale(-1.0)  # negative so most-occupied eigenvalues come first
        Dhelp.transform(wfn.S())
        X = wfn.S().clone()
        X.power(-0.5, 1.0e-10)
        Dhelp.transform(X)

        Cno = core.Matrix("SAD NO temp", Dhelp.rowdim(), Dhelp.coldim())
        Dhelp.diagonalize(Cno, wfn.epsilon_a())
        wfn.epsilon_b().copy(wfn.epsilon_a())
        wfn.Ca().gemm(False, False, 1.0, X, Cno, 0.0)
        wfn.Cb().copy(wfn.Ca())

        # ROHF needs Ct_ = X^T S Ca.
        if hasattr(wfn, "Ct"):
            Ct = wfn.Ct()
            # X already holds S^(-1/2); we want X.T @ S @ Ca but X.T == X for a symmetric
            # S^(-1/2). Use linalg.triplet if available; here Matrix.gemm chains suffice.
            tmp = core.Matrix("tmp", wfn.S().rowdim(), wfn.Ca().coldim())
            tmp.gemm(False, False, 1.0, wfn.S(), wfn.Ca(), 0.0)
            Ct.gemm(True, False, 1.0, X, tmp, 0.0)
    else:
        wfn.Da().copy(sad["Da"])
        wfn.Db().copy(sad["Db"])

        Ca_sad = sad["Ca"]
        Cb_sad = sad["Cb"]
        Ca = wfn.Ca()
        Cb = wfn.Cb()
        # Allowed MO count per irrep is the orbital-basis MO dimension (Ca's column count).
        nmopi = Ca.coldim()
        sad_dim = core.Dimension(Ca.nirrep(), "SAD Dimensions")

        for h in range(Ca.nirrep()):
            nso = Ca_sad.rowdim()[h]
            nmo_sad = Ca_sad.coldim()[h]
            nmo = min(nmo_sad, nmopi[h])
            sad_dim[h] = nmo
            if nso == 0 or nmo == 0:
                continue
            Ca.nph[h][:, :nmo] = Ca_sad.nph[h][:, :nmo]
            Cb.nph[h][:, :nmo] = Cb_sad.nph[h][:, :nmo]

        wfn.set_nalphapi(sad_dim)
        wfn.set_nbetapi(sad_dim)

        # ROHF / CUHF: refresh the total density Dt_ = Da_ + Db_.
        if hasattr(wfn, "Dt"):
            Dt = wfn.Dt()
            Dt.copy(wfn.Da())
            Dt.add(wfn.Db())


def populate_huckel_guess(
    wfn,
    atomic_bases: Sequence[core.BasisSet],
    atomic_fit_bases,
    updated_rule: bool,
) -> None:
    """Populate an HF wavefunction's Fa/Fb from the GWH (Huckel) guess."""
    Fhuckel = compute_huckel_guess(
        wfn.basisset(), atomic_bases, atomic_fit_bases, updated_rule, AO2SO=wfn.aotoso()
    )
    wfn.Fa().copy(Fhuckel)
    wfn.Fb().copy(Fhuckel)


def _gwh_k(Ei: float, Ej: float, updated_rule: bool) -> float:
    k = 1.75
    if not updated_rule:
        return k
    delta = (Ei - Ej) / (Ei + Ej)
    d2 = delta * delta
    return k + d2 + d2 * d2 * (1 - k)


def compute_huckel_guess(
    basis: core.BasisSet,
    atomic_bases: Sequence[core.BasisSet],
    atomic_fit_bases,
    updated_rule: bool,
    AO2SO: "core.Matrix | None" = None,
) -> core.Matrix:
    """Construct the GWH Fock matrix from atomic occupied MOs."""
    sad = compute_sad_guess(basis, atomic_bases, atomic_fit_bases, AO2SO=AO2SO)
    Chu = np.asarray(sad["HuckelC"])
    Ehu = np.asarray(sad["HuckelE"])
    nbf = basis.nbf()
    nhu = Chu.shape[1]

    mints = core.MintsHelper(basis)
    S = np.asarray(mints.ao_overlap())
    SChu = S @ Chu
    ChuSChu = Chu.T @ SChu

    huckelmo = np.zeros((nhu, nhu))
    for i in range(nhu):
        huckelmo[i, i] = Ehu[i]
        for j in range(nhu):
            huckelmo[i, j] = (
                0.5
                * _gwh_k(Ehu[i], Ehu[j], updated_rule)
                * ChuSChu[i, j]
                * (Ehu[i] + Ehu[j])
            )

    huckelao = SChu @ huckelmo @ SChu.T
    if AO2SO is None:
        nbf_dim = core.Dimension([basis.nbf()], "AO")
        AO2SO = core.Matrix("AO2SO (C1)", nbf_dim, nbf_dim)
        for i in range(basis.nbf()):
            AO2SO.set(0, i, i, 1.0)
    huckel = core.Matrix("Huckel SO matrix", AO2SO.coldim(), AO2SO.coldim())
    huckel.apply_symmetry(core.Matrix.from_array(huckelao), AO2SO)
    return huckel
