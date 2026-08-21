/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "local_qia.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

namespace psi {

namespace {

/// Solve A X = B in place on B, with A symmetric and both stored row-major.
///
/// A is symmetric here (an overlap block), so it needs no transposition to be
/// read as Fortran-ordered; B does, which is what the two copies are for.
void solve_in_place(double* A, double* B, int n, int nrhs) {
    if (n == 0 || nrhs == 0) return;

    std::vector<double> B_fortran(static_cast<size_t>(n) * nrhs);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < nrhs; k++) {
            B_fortran[static_cast<size_t>(k) * n + i] = B[static_cast<size_t>(i) * nrhs + k];
        }
    }

    std::vector<int> ipiv(n);
    int errcode = C_DGESV(n, nrhs, A, n, ipiv.data(), B_fortran.data(), n);
    if (errcode != 0) {
        throw PSIEXCEPTION(
            "LocalQiaBuilder: the Boughton-Pulay refit of the LMO coefficients hit a singular "
            "overlap block. Either the primary basis is linearly dependent on this domain or the "
            "coefficient tolerance is admitting a degenerate set; set_bp_refit(false) transforms "
            "with the raw coefficient slice instead.");
    }

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < nrhs; k++) {
            B[static_cast<size_t>(i) * nrhs + k] = B_fortran[static_cast<size_t>(k) * n + i];
        }
    }
}

}  // namespace

LocalQiaBuilder::LocalQiaBuilder(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux)
    : primary_(primary),
      aux_(aux),
      ints_tolerance_(1.0e-10),
      coef_tol_lmo_(1.0e-4),
      coef_tol_pao_(1.0e-4),
      bp_refit_(true),
      domains_set_(false),
      shell_pairs_computed_(0),
      shell_pairs_offered_(0) {
    if (!primary_ || !aux_) {
        throw PSIEXCEPTION("LocalQiaBuilder: both a primary and an auxiliary basis set are required.");
    }
    if (primary_->molecule()->natom() != aux_->molecule()->natom()) {
        throw PSIEXCEPTION(
            "LocalQiaBuilder: the primary and auxiliary basis sets are on different molecules. The "
            "domains are expressed per atom, so the two have to agree on what the atoms are.");
    }

    natom_ = aux_->molecule()->natom();
    nbf_ = primary_->nbf();
    naux_ = aux_->nbf();

    nthread_ = Process::environment.get_n_threads();

    // Which primary functions and shells sit on each atom. Fixed by the basis
    // sets, so this is the one map that never depends on the demand.
    atom_to_bf_.assign(natom_, {});
    for (int mu = 0; mu < nbf_; mu++) {
        atom_to_bf_[primary_->function_to_center(mu)].push_back(mu);
    }
    atom_to_shell_.assign(natom_, {});
    for (int M = 0; M < primary_->nshell(); M++) {
        atom_to_shell_[primary_->shell_to_center(M)].push_back(M);
    }

    // The auxiliary function range of each atom. Shells are ordered by center,
    // so an atom's auxiliary functions are contiguous, which is what lets a
    // caller scatter a whole block with one slice; check it rather than
    // assuming it, because a caller that trusted a wrong range would get a
    // silently misplaced block rather than an error.
    atom_aux_range_.assign(natom_, {0, 0});
    std::vector<int> aux_count(natom_, 0);
    std::vector<int> aux_first(natom_, naux_);
    std::vector<int> aux_last(natom_, 0);
    for (int Q = 0; Q < aux_->nshell(); Q++) {
        int center = aux_->shell_to_center(Q);
        int start = aux_->shell(Q).function_index();
        int nq = aux_->shell(Q).nfunction();
        aux_count[center] += nq;
        aux_first[center] = std::min(aux_first[center], start);
        aux_last[center] = std::max(aux_last[center], start + nq);
    }
    for (int A = 0; A < natom_; A++) {
        if (aux_count[A] == 0) continue;
        if (aux_last[A] - aux_first[A] != aux_count[A]) {
            throw PSIEXCEPTION(
                "LocalQiaBuilder: the auxiliary basis functions of an atom are not contiguous, "
                "which the per-atom blocking assumes.");
        }
        atom_aux_range_[A] = {aux_first[A], aux_last[A]};
    }
}

void LocalQiaBuilder::set_nthreads(int nthread) {
    if (nthread < 1) {
        throw PSIEXCEPTION("LocalQiaBuilder: set_nthreads needs a positive thread count.");
    }
    nthread_ = nthread;
}

void LocalQiaBuilder::set_ints_tolerance(double tol) {
    ints_tolerance_ = tol;
    metric_shell_diag_.clear();
}

void LocalQiaBuilder::set_coef_tolerance_lmo(double tol) { coef_tol_lmo_ = tol; }

void LocalQiaBuilder::set_coef_tolerance_pao(double tol) { coef_tol_pao_ = tol; }

void LocalQiaBuilder::set_bp_refit(bool on) { bp_refit_ = on; }

void LocalQiaBuilder::set_domains(const std::vector<std::vector<int>>& atom_to_lmos,
                                  const std::vector<std::vector<int>>& atom_to_paos) {
    if (static_cast<int>(atom_to_lmos.size()) != natom_ || static_cast<int>(atom_to_paos.size()) != natom_) {
        throw PSIEXCEPTION(
            "LocalQiaBuilder: set_domains needs one LMO list and one PAO list per atom of the "
            "auxiliary basis's molecule.");
    }
    atom_to_lmos_ = atom_to_lmos;
    atom_to_paos_ = atom_to_paos;
    domains_set_ = true;
}

void LocalQiaBuilder::build_metric_shell_diag() {
    // Only the diagonal of the auxiliary metric is wanted, so only the diagonal
    // shell blocks are computed: nshell integral blocks rather than the
    // nshell^2 a full metric would take. The screening estimate reads it as a
    // per-shell maximum anyway.
    metric_shell_diag_.assign(aux_->nshell(), 0.0);

    auto zero = BasisSet::zero_ao_basis_set();
    auto factory = std::make_shared<IntegralFactory>(aux_, zero, aux_, zero);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread_);
    eris[0] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    for (int thread = 1; thread < nthread_; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(eris.front()->clone());
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread_)
    for (int Q = 0; Q < aux_->nshell(); Q++) {
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        int nq = aux_->shell(Q).nfunction();
        eris[thread]->compute_shell(Q, 0, Q, 0);
        const double* buffer = eris[thread]->buffer();
        double diag = 0.0;
        for (int q = 0; q < nq; q++) {
            diag = std::max(diag, buffer[q * nq + q]);
        }
        metric_shell_diag_[Q] = diag;
    }
}

void LocalQiaBuilder::build_bf_domains(SharedMatrix C_lmo, SharedMatrix C_pao) {
    atom_in_domain1_.assign(natom_, {});
    atom_in_domain2_.assign(natom_, {});
    domain_shells1_.assign(natom_, {});
    domain_shells2_.assign(natom_, {});
    domain_bfs1_.assign(natom_, {});
    domain_bfs2_.assign(natom_, {});

    double** Clmo = C_lmo->pointer();
    double** Cpao = C_pao->pointer();

    // Which atoms an orbital lives on, by the coefficient tolerance. Atom
    // granularity rather than basis-function granularity, matching DLPNO: a
    // shell has to be computed whole in any case, and keeping whole atoms is
    // what makes the refit's overlap block well conditioned.
    auto orbital_atoms = [&](double** C, int col, double tol, std::vector<char>& into) {
        for (int A = 0; A < natom_; A++) {
            if (into[A]) continue;
            for (int mu : atom_to_bf_[A]) {
                if (std::fabs(C[mu][col]) > tol) {
                    into[A] = 1;
                    break;
                }
            }
        }
    };

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread_)
    for (int A = 0; A < natom_; A++) {
        if (atom_to_lmos_[A].empty() || atom_to_paos_[A].empty()) continue;

        std::vector<char> in1(natom_, 0);
        std::vector<char> in2(natom_, 0);
        for (int i : atom_to_lmos_[A]) orbital_atoms(Clmo, i, coef_tol_lmo_, in1);
        for (int u : atom_to_paos_[A]) orbital_atoms(Cpao, u, coef_tol_pao_, in2);

        for (int B = 0; B < natom_; B++) {
            if (in1[B]) {
                domain_shells1_[A].insert(domain_shells1_[A].end(), atom_to_shell_[B].begin(),
                                          atom_to_shell_[B].end());
                domain_bfs1_[A].insert(domain_bfs1_[A].end(), atom_to_bf_[B].begin(), atom_to_bf_[B].end());
            }
            if (in2[B]) {
                domain_shells2_[A].insert(domain_shells2_[A].end(), atom_to_shell_[B].begin(),
                                          atom_to_shell_[B].end());
                domain_bfs2_[A].insert(domain_bfs2_[A].end(), atom_to_bf_[B].begin(), atom_to_bf_[B].end());
            }
        }
        atom_in_domain1_[A] = std::move(in1);
        atom_in_domain2_[A] = std::move(in2);
    }
}

std::vector<SharedMatrix> LocalQiaBuilder::compute(SharedMatrix C_lmo, SharedMatrix C_pao, SharedMatrix S) {
    if (!domains_set_) {
        throw PSIEXCEPTION("LocalQiaBuilder: call set_domains() before compute().");
    }
    if (!C_lmo || !C_pao) {
        throw PSIEXCEPTION("LocalQiaBuilder: compute() needs both coefficient matrices.");
    }
    if (C_lmo->nrow() != nbf_ || C_pao->nrow() != nbf_) {
        throw PSIEXCEPTION(
            "LocalQiaBuilder: the coefficient matrices must have one row per primary basis "
            "function.");
    }
    if (bp_refit_ && (!S || S->nrow() != nbf_ || S->ncol() != nbf_)) {
        throw PSIEXCEPTION(
            "LocalQiaBuilder: the Boughton-Pulay refit needs the (nbf, nbf) AO overlap. Pass it, or "
            "set_bp_refit(false) to transform with the raw coefficient slice.");
    }

    const int nlmo_total = C_lmo->ncol();
    const int npao_total = C_pao->ncol();
    for (int A = 0; A < natom_; A++) {
        for (int i : atom_to_lmos_[A]) {
            if (i < 0 || i >= nlmo_total) {
                throw PSIEXCEPTION("LocalQiaBuilder: an LMO index in the domains is out of range.");
            }
        }
        for (int u : atom_to_paos_[A]) {
            if (u < 0 || u >= npao_total) {
                throw PSIEXCEPTION("LocalQiaBuilder: a PAO index in the domains is out of range.");
            }
        }
    }

    build_bf_domains(C_lmo, C_pao);

    const bool screening = ints_tolerance_ > 0.0;
    if (screening && metric_shell_diag_.empty()) build_metric_shell_diag();

    // S C_lmo once for the whole run: the refit's right-hand side is a slice of
    // it, and slicing is cheaper than re-forming it per atom.
    SharedMatrix SC_lmo;
    if (bp_refit_) SC_lmo = linalg::doublet(S, C_lmo, false, false);

    // Per-atom transform matrices. Both depend only on the atom, so they are
    // built once here rather than once per auxiliary shell as DLPNO does - an
    // atom carries several auxiliary shells and the refit is a linear solve.
    std::vector<SharedMatrix> X(natom_);          // (|bfs1|, |lmos|)
    std::vector<SharedMatrix> Cpao_dom(natom_);   // (|bfs2|, |paos|)
    std::vector<SharedMatrix> out(natom_);

    double** Clmo = C_lmo->pointer();
    double** Cpao = C_pao->pointer();

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread_)
    for (int A = 0; A < natom_; A++) {
        const auto& lmos = atom_to_lmos_[A];
        const auto& paos = atom_to_paos_[A];
        if (lmos.empty() || paos.empty()) continue;
        if (atom_aux_range_[A].second == atom_aux_range_[A].first) continue;

        const auto& bfs1 = domain_bfs1_[A];
        const auto& bfs2 = domain_bfs2_[A];
        const int nbf1 = bfs1.size();
        const int nbf2 = bfs2.size();
        const int nlmo = lmos.size();
        const int npao = paos.size();

        out[A] = std::make_shared<Matrix>("(Q|iu) domain", atom_aux_range_[A].second - atom_aux_range_[A].first,
                                          nlmo * npao);

        // A domain can come out empty if the coefficient tolerance admits
        // nothing for any requested orbital. The block is then exactly zero,
        // which the freshly allocated matrix already is, and there is no
        // transform to set up.
        if (nbf1 == 0 || nbf2 == 0) continue;

        Cpao_dom[A] = std::make_shared<Matrix>("C_pao domain", nbf2, npao);
        double** Cp = Cpao_dom[A]->pointer();
        for (int n = 0; n < nbf2; n++) {
            for (int u = 0; u < npao; u++) Cp[n][u] = Cpao[bfs2[n]][paos[u]];
        }

        X[A] = std::make_shared<Matrix>("C_lmo domain", nbf1, nlmo);
        double** Xp = X[A]->pointer();
        // The refit is the identity when the domain is the whole basis - it
        // solves S X = S C - so skip it there. That is not only faster: running
        // it anyway would push every element through the full overlap's
        // condition number for a result already known exactly, which is what an
        // exactness check would then be measuring.
        const bool refit = bp_refit_ && nbf1 < nbf_;
        if (refit) {
            double** SC = SC_lmo->pointer();
            for (int m = 0; m < nbf1; m++) {
                for (int i = 0; i < nlmo; i++) Xp[m][i] = SC[bfs1[m]][lmos[i]];
            }
            auto S_dom = std::make_shared<Matrix>("S domain", nbf1, nbf1);
            double** Sd = S_dom->pointer();
            double** Sp = S->pointer();
            for (int m = 0; m < nbf1; m++) {
                for (int n = 0; n < nbf1; n++) Sd[m][n] = Sp[bfs1[m]][bfs1[n]];
            }
            solve_in_place(Sd[0], Xp[0], nbf1, nlmo);
        } else {
            for (int m = 0; m < nbf1; m++) {
                for (int i = 0; i < nlmo; i++) Xp[m][i] = Clmo[bfs1[m]][lmos[i]];
            }
        }
    }

    auto factory = std::make_shared<IntegralFactory>(aux_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread_);
    eris[0] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    for (int thread = 1; thread < nthread_; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(eris.front()->clone());
    }

    size_t computed = 0;
    size_t offered = 0;

    // One auxiliary shell per task. Shells of the same atom write disjoint rows
    // of that atom's block, so nothing is shared across tasks but the integral
    // engines, which are per thread.
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthread_) reduction(+ : computed, offered)
    for (int Q = 0; Q < aux_->nshell(); Q++) {
        const int centerQ = aux_->shell_to_center(Q);
        if (!out[centerQ] || !X[centerQ]) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        const int nq = aux_->shell(Q).nfunction();
        const int qstart = aux_->shell(Q).function_index();

        const auto& bfs1 = domain_bfs1_[centerQ];
        const auto& bfs2 = domain_bfs2_[centerQ];
        const int nbf1 = bfs1.size();
        const int nbf2 = bfs2.size();
        const int nlmo = atom_to_lmos_[centerQ].size();
        const int npao = atom_to_paos_[centerQ].size();

        // Global basis function to position within this atom's domain, -1 when
        // outside it. Dense because the scatter below is in the inner loop.
        std::vector<int> pos1(nbf_, -1);
        std::vector<int> pos2(nbf_, -1);
        for (int m = 0; m < nbf1; m++) pos1[bfs1[m]] = m;
        for (int n = 0; n < nbf2; n++) pos2[bfs2[n]] = n;

        // (mn|Q) for this shell's auxiliary functions, zeroed because every
        // screened-out or out-of-domain shell pair contributes exactly zero.
        std::vector<double> block(static_cast<size_t>(nq) * nbf1 * nbf2, 0.0);

        for (int M : domain_shells1_[centerQ]) {
            const int nm = primary_->shell(M).nfunction();
            const int mstart = primary_->shell(M).function_index();
            const int centerM = primary_->shell_to_center(M);

            for (int N : domain_shells2_[centerQ]) {
                const int nn = primary_->shell(N).nfunction();
                const int nstart = primary_->shell(N).function_index();
                const int centerN = primary_->shell_to_center(N);

                // (MN|Q) and (NM|Q) are the same integrals. When both orderings
                // are in domain, compute one and scatter it twice.
                const bool MN_symmetry =
                    atom_in_domain1_[centerQ][centerN] && atom_in_domain2_[centerQ][centerM];
                if (N < M && MN_symmetry) continue;
                offered++;

                if (screening &&
                    metric_shell_diag_[Q] * eris[thread]->shell_pair_value(M, N) <
                        ints_tolerance_ * ints_tolerance_) {
                    continue;
                }
                computed++;

                eris[thread]->compute_shell(Q, 0, M, N);
                const double* buffer = eris[thread]->buffer();

                for (int q = 0, index = 0; q < nq; q++) {
                    double* dest = block.data() + static_cast<size_t>(q) * nbf1 * nbf2;
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            dest[static_cast<size_t>(pos1[mstart + m]) * nbf2 + pos2[nstart + n]] = buffer[index];
                        }
                    }
                }

                if (N > M && MN_symmetry) {
                    for (int q = 0, index = 0; q < nq; q++) {
                        double* dest = block.data() + static_cast<size_t>(q) * nbf1 * nbf2;
                        for (int m = 0; m < nm; m++) {
                            for (int n = 0; n < nn; n++, index++) {
                                dest[static_cast<size_t>(pos1[nstart + n]) * nbf2 + pos2[mstart + m]] =
                                    buffer[index];
                            }
                        }
                    }
                }
            }
        }

        // X^T (mn|Q) C_pao, straight into the atom's output rows: the second
        // axis of a row is already (lmo, pao) row-major, which is what the
        // second GEMM writes.
        std::vector<double> half(static_cast<size_t>(nlmo) * nbf2);
        double** Xp = X[centerQ]->pointer();
        double** Cp = Cpao_dom[centerQ]->pointer();
        double** dest = out[centerQ]->pointer();
        const int row0 = atom_aux_range_[centerQ].first;

        for (int q = 0; q < nq; q++) {
            const double* src = block.data() + static_cast<size_t>(q) * nbf1 * nbf2;
            C_DGEMM('T', 'N', nlmo, nbf2, nbf1, 1.0, Xp[0], nlmo, const_cast<double*>(src), nbf2, 0.0,
                    half.data(), nbf2);
            C_DGEMM('N', 'N', nlmo, npao, nbf2, 1.0, half.data(), nbf2, Cp[0], npao, 0.0,
                    dest[qstart + q - row0], npao);
        }
    }

    shell_pairs_computed_ = computed;
    shell_pairs_offered_ = offered;

    return out;
}

}  // namespace psi
