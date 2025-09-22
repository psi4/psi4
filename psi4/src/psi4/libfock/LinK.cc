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

#include "jk.h"
#include "SplitJK.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

LinK::LinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {
    timer_on("LinK: Setup");

    // => General Setup <= //

    // thread count
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif

    // set default lr_symmetric_ value
    lr_symmetric_ = true;

    // set up LinK integral tolerance
    if (options["LINK_INTS_TOLERANCE"].has_changed() && options.get_str("SCREENING") != "NONE") {
        linK_ints_cutoff_ = options.get_double("LINK_INTS_TOLERANCE");
    } else {
        linK_ints_cutoff_ = cutoff_;
    }

    timer_off("LinK: Setup");
}

LinK::~LinK() {}

size_t LinK::num_computed_shells() {
    return num_computed_shells_;
}

void LinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> LinK: Linear Exchange K <==\n\n");

        outfile->Printf("    K Screening Cutoff:%11.0E\n", linK_ints_cutoff_);
    }
}

// build the K matrix using Ochsenfelds's Linear Exchange (LinK) algorithm
// To follow this code, compare with figure 1 of DOI: 10.1063/1.476741
void LinK::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    // LinK does not support non-symmetric matrices
    if (!lr_symmetric_) {
        throw PSIEXCEPTION("Non-symmetric K matrix builds are currently not supported in the LinK algorithm.");
    }

    // ==> Prep Auxiliary Quantities <== //

    // => Sizing <= //
    int nshell = primary_->nshell();
    int nbf = primary_->nbf();
    int nthread = nthreads_;

    // => Atom Blocking <= //
    std::vector<int> shell_endpoints_for_atom;
    std::vector<int> basis_endpoints_for_shell;

    int atomic_ind = -1;
    for (int P = 0; P < nshell; P++) {
        if (primary_->shell(P).ncenter() > atomic_ind) {
            shell_endpoints_for_atom.push_back(P);
            atomic_ind++;
        }
        basis_endpoints_for_shell.push_back(primary_->shell_to_basis_function(P));
    }
    shell_endpoints_for_atom.push_back(nshell);
    basis_endpoints_for_shell.push_back(nbf);

    size_t natom = shell_endpoints_for_atom.size() - 1;

    size_t max_functions_per_atom = 0L;
    for (size_t atom = 0; atom < natom; atom++) {
        size_t size = 0L;
        for (int P = shell_endpoints_for_atom[atom]; P < shell_endpoints_for_atom[atom + 1]; P++) {
            size += primary_->shell(P).nfunction();
        }
        max_functions_per_atom = std::max(max_functions_per_atom, size);
    }

    if (debug_) {
        outfile->Printf("  ==> LinK: Atom Blocking <==\n\n");
        for (size_t atom = 0; atom < natom; atom++) {
            outfile->Printf("  Atom: %3d, Atom Start: %4d, Atom End: %4d\n", atom, shell_endpoints_for_atom[atom],
                            shell_endpoints_for_atom[atom + 1]);
            for (int P = shell_endpoints_for_atom[atom]; P < shell_endpoints_for_atom[atom + 1]; P++) {
                int size = primary_->shell(P).nfunction();
                int off = primary_->shell(P).function_index();
                int off2 = basis_endpoints_for_shell[P];
                outfile->Printf("    Shell: %4d, Size: %4d, Offset: %4d, Offset2: %4d\n", P, size, off,
                                off2);
            }
        }
        outfile->Printf("\n");
    }

    // ==> Prep Atom Pairs <== //
    // Atom-pair blocking inherited from DirectJK code
    // TODO: Test shell-pair blocking

    std::vector<std::pair<int, int>> atom_pairs;
    for (size_t Patom = 0; Patom < natom; Patom++) {
        for (size_t Qatom = 0; Qatom <= Patom; Qatom++) {
            bool found = false;
            for (int P = shell_endpoints_for_atom[Patom]; P < shell_endpoints_for_atom[Patom + 1]; P++) {
                for (int Q = shell_endpoints_for_atom[Qatom]; Q < shell_endpoints_for_atom[Qatom + 1]; Q++) {
                    if (eri_computers[0]->shell_pair_significant(P, Q)) {
                        found = true;
                        atom_pairs.emplace_back(Patom, Qatom);
                        break;
                    }
                }
                if (found) break;
            }
        }
    }

    // ==> Start "Pre-ordering and pre-selection to find significant elements in P_uv" in Fig. 1 of paper <== //

    // ==> Prep Bra-Bra Shell Pairs <== //

    // A comparator used for sorting integral screening values
    auto screen_compare = [](const std::pair<int, double> &a,
                                    const std::pair<int, double> &b) { return a.second > b.second; };

    std::vector<std::vector<int>> significant_bras(nshell);
    double max_integral = eri_computers[0]->max_integral();

#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PQ_shell_values;
        for (size_t Q = 0; Q < nshell; Q++) {
            double pq_pq = std::sqrt(eri_computers[0]->shell_ceiling2(P, Q, P, Q));
            double schwarz_value = std::sqrt(pq_pq * max_integral);
            if (schwarz_value >= cutoff_) {
                PQ_shell_values.emplace_back(Q, schwarz_value);
            }
        }
        std::sort(PQ_shell_values.begin(), PQ_shell_values.end(), screen_compare);

        for (const auto& value : PQ_shell_values) {
            significant_bras[P].push_back(value.first);
        }
    }

    // ==> Prep Bra-Ket Shell Pairs <== //

    // => Calculate Shell Ceilings <= //
    std::vector<double> shell_ceilings(nshell, 0.0);

    // sqrt(Umax|Umax) in Ochsenfeld Eq. 3
#pragma omp parallel for
    for (int P = 0; P < nshell; P++) {
        for (int Q = 0; Q <= P; Q++) {
            double val = std::sqrt(eri_computers[0]->shell_ceiling2(P, Q, P, Q));
            shell_ceilings[P] = std::max(shell_ceilings[P], val);
#pragma omp critical
            shell_ceilings[Q] = std::max(shell_ceilings[Q], val);
        }
    }

    std::vector<std::vector<int>> significant_kets(nshell);

    // => Use shell ceilings to compute significant ket-shells for each bra-shell <= //
    // => Implementation of Eq. 4 in paper <= //
#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PR_shell_values;
        for (size_t R = 0; R < nshell; R++) {
            double screen_val = shell_ceilings[P] * shell_ceilings[R] * eri_computers[0]->shell_pair_max_density(P, R);
            if (screen_val >= linK_ints_cutoff_) {
                PR_shell_values.emplace_back(R, screen_val);
            }
        }
        std::sort(PR_shell_values.begin(), PR_shell_values.end(), screen_compare);

        for (const auto& value : PR_shell_values) {
            significant_kets[P].push_back(value.first);
        }
    }

    size_t natom_pair = atom_pairs.size();

    // ==> End "Pre-ordering and pre-selection to find significant elements in P_uv" in Fig. 1 of paper <== //

    // ==> Intermediate Buffers <== //

    // Temporary buffers used during the K contraction process to
    // Take full advantage of permutational symmetry of ERIs
    std::vector<std::vector<SharedMatrix>> KT;

    // To prevent race conditions, give every thread a buffer
    for (int thread = 0; thread < nthread; thread++) {
        std::vector<SharedMatrix> K2;
        for (size_t ind = 0; ind < D.size(); ind++) {
            // (pq|rs) can be contracted into Kpr, Kps, Kqr, Kqs (hence the 4)
            K2.push_back(std::make_shared<Matrix>("KT (linK)", 4 * max_functions_per_atom, nbf));
        }
        KT.push_back(K2);
    }

    // ==> Start "Loop over significant 'bra'-shell pairs uh" in Fig. 1 of paper <== //
    // Number of computed shell quartets is tracked for benchmarking purposes
    num_computed_shells_ = 0L;
    size_t computed_shells = 0L;

    // ==> Integral Formation Loop <== //

#pragma omp parallel for num_threads(nthread) schedule(dynamic) reduction(+ : computed_shells)
    for (size_t ipair = 0L; ipair < natom_pair; ipair++) { // O(N) shell-pairs in asymptotic limit
        int Patom = atom_pairs[ipair].first;
        int Qatom = atom_pairs[ipair].second;

        // Number of shells per atom
        int nPshell = shell_endpoints_for_atom[Patom + 1] - shell_endpoints_for_atom[Patom];
        int nQshell = shell_endpoints_for_atom[Qatom + 1] - shell_endpoints_for_atom[Qatom];

        // First shell per atom
        int Pstart = shell_endpoints_for_atom[Patom];
        int Qstart = shell_endpoints_for_atom[Qatom];

        // Number of basis functions per atom
        int nPbasis = basis_endpoints_for_shell[Pstart + nPshell] - basis_endpoints_for_shell[Pstart];
        int nQbasis = basis_endpoints_for_shell[Qstart + nQshell] - basis_endpoints_for_shell[Qstart];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Keep track of contraction indices for stripeout (Towards end of this function)
        std::vector<std::unordered_set<int>> P_stripeout_list(nPshell);
        std::vector<std::unordered_set<int>> Q_stripeout_list(nQshell);

        bool touched = false;
        for (int P = Pstart; P < Pstart + nPshell; P++) {
            for (int Q = Qstart; Q < Qstart + nQshell; Q++) {

                if (Q > P) continue;
                if (!eri_computers[0]->shell_pair_significant(P, Q)) continue;

                int dP = P - Pstart;
                int dQ = Q - Qstart;

                // => Start "Formation of Significant Shell Pair List ML" in Fig. 1 of Paper <= //

                // Significant ket shell pairs RS for bra shell pair PQ
                // represents the merge of ML_P and ML_Q (mini-lists) as defined in Oschenfeld
                // Unordered set structure allows for automatic merging as new elements are added
                std::unordered_set<int> ML_PQ;

                // Form ML_P as part of ML_PQ
                for (const int R : significant_kets[P]) {
                    bool is_significant = false;
                    for (const int S : significant_bras[R]) {
                        double screen_val = eri_computers[0]->shell_pair_max_density(P, R) * std::sqrt(eri_computers[0]->shell_ceiling2(P, Q, R, S));

                        if (screen_val >= linK_ints_cutoff_) {
                            if (!is_significant) is_significant = true;
                            int RS = (R >= S) ? (R * nshell + S) : (S * nshell + R);
                            if (RS > P * nshell + Q) continue;
                            ML_PQ.emplace(RS);
                            Q_stripeout_list[dQ].emplace(S);
                        }
                        else break;
                    }
                    if (!is_significant) break;
                }

                // Form ML_Q as part of ML_PQ
                for (const int R : significant_kets[Q]) {
                    bool is_significant = false;
                    for (const int S : significant_bras[R]) {
                        double screen_val = eri_computers[0]->shell_pair_max_density(Q, R) * std::sqrt(eri_computers[0]->shell_ceiling2(P, Q, R, S));

                        if (screen_val >= linK_ints_cutoff_) {
                            if (!is_significant) is_significant = true;
                            int RS = (R >= S) ? (R * nshell + S) : (S * nshell + R);
                            if (RS > P * nshell + Q) continue;
                            ML_PQ.emplace(RS);
                            P_stripeout_list[dP].emplace(S);
                        }
                        else break;
                    }
                    if (!is_significant) break;
                }

                // Loop over significant RS pairs
                for (const int RS : ML_PQ) {

                    int R = RS / nshell;
                    int S = RS % nshell;

                    if (!eri_computers[0]->shell_pair_significant(R, S)) continue;
                    if (!eri_computers[0]->shell_significant(P, Q, R, S)) continue;

                    if (eri_computers[thread]->compute_shell(P, Q, R, S) == 0)
                        continue;
                    computed_shells++;

                    const double* buffer = eri_computers[thread]->buffer();

                    // Number of basis functions in shells P, Q, R, S
                    int shell_P_nfunc = primary_->shell(P).nfunction();
                    int shell_Q_nfunc = primary_->shell(Q).nfunction();
                    int shell_R_nfunc = primary_->shell(R).nfunction();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    // Basis Function Starting index for shell
                    int shell_P_start = primary_->shell(P).function_index();
                    int shell_Q_start = primary_->shell(Q).function_index();
                    int shell_R_start = primary_->shell(R).function_index();
                    int shell_S_start = primary_->shell(S).function_index();

                    // Basis Function offset from first basis function in the atom
                    int shell_P_offset = basis_endpoints_for_shell[P] - basis_endpoints_for_shell[Pstart];
                    int shell_Q_offset = basis_endpoints_for_shell[Q] - basis_endpoints_for_shell[Qstart];

                    for (size_t ind = 0; ind < D.size(); ind++) {
                        double** Kp = K[ind]->pointer();
                        double** Dp = D[ind]->pointer();
                        double** KTp = KT[thread][ind]->pointer();
                        const double* buffer2 = buffer;

                        if (!touched) {
                            ::memset((void*)KTp[0L * max_functions_per_atom], '\0', nPbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[1L * max_functions_per_atom], '\0', nPbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[2L * max_functions_per_atom], '\0', nQbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[3L * max_functions_per_atom], '\0', nQbasis * nbf * sizeof(double));
                        }

                        // Four pointers needed for PR, PS, QR, QS
                        double* K1p = KTp[0L * max_functions_per_atom];
                        double* K2p = KTp[1L * max_functions_per_atom];
                        double* K3p = KTp[2L * max_functions_per_atom];
                        double* K4p = KTp[3L * max_functions_per_atom];

                        double prefactor = 1.0;
                        if (P == Q) prefactor *= 0.5;
                        if (R == S) prefactor *= 0.5;
                        if (P == R && Q == S) prefactor *= 0.5;

                        // => Computing integral contractions to K buffers <= //
                        for (int p = 0; p < shell_P_nfunc; p++) {
                            for (int q = 0; q < shell_Q_nfunc; q++) {
                                for (int r = 0; r < shell_R_nfunc; r++) {
                                    for (int s = 0; s < shell_S_nfunc; s++) {

                                        K1p[(p + shell_P_offset) * nbf + r + shell_R_start] +=
                                            prefactor * (Dp[q + shell_Q_start][s + shell_S_start]) * (*buffer2);
                                        K2p[(p + shell_P_offset) * nbf + s + shell_S_start] +=
                                            prefactor * (Dp[q + shell_Q_start][r + shell_R_start]) * (*buffer2);
                                        K3p[(q + shell_Q_offset) * nbf + r + shell_R_start] +=
                                            prefactor * (Dp[p + shell_P_start][s + shell_S_start]) * (*buffer2);
                                        K4p[(q + shell_Q_offset) * nbf + s + shell_S_start] +=
                                            prefactor * (Dp[p + shell_P_start][r + shell_R_start]) * (*buffer2);

                                        buffer2++;
                                    }
                                }
                            }
                        }
                    }
                    touched = true;
                }
            }
        }

        // => Master shell quartet loops <= //

        if (!touched) continue;

        // => Stripe out (Writing to K matrix) <= //
        for (auto& KTmat : KT[thread]) {
            KTmat->scale(2.0);
        }

        for (size_t ind = 0; ind < D.size(); ind++) {
            double** KTp = KT[thread][ind]->pointer();
            double** Kp = K[ind]->pointer();

            double* K1p = KTp[0L * max_functions_per_atom];
            double* K2p = KTp[1L * max_functions_per_atom];
            double* K3p = KTp[2L * max_functions_per_atom];
            double* K4p = KTp[3L * max_functions_per_atom];

            // K_PR and K_PS
            for (int P = Pstart; P < Pstart + nPshell; P++) {
                int dP = P - Pstart;
                int shell_P_start = primary_->shell(P).function_index();
                int shell_P_nfunc = primary_->shell(P).nfunction();
                int shell_P_offset = basis_endpoints_for_shell[P] - basis_endpoints_for_shell[Pstart];
                for (const int S : P_stripeout_list[dP]) {
                    int shell_S_start = primary_->shell(S).function_index();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    for (int p = 0; p < shell_P_nfunc; p++) {
                        for (int s = 0; s < shell_S_nfunc; s++) {
#pragma omp atomic
                            Kp[shell_P_start + p][shell_S_start + s] += K1p[(p + shell_P_offset) * nbf + s + shell_S_start];
#pragma omp atomic
                            Kp[shell_P_start + p][shell_S_start + s] += K2p[(p + shell_P_offset) * nbf + s + shell_S_start];
                        }
                    }

                }
            }

            // K_QR and K_QS
            for (int Q = Qstart; Q < Qstart + nQshell; Q++) {
                int dQ = Q - Qstart;
                int shell_Q_start = primary_->shell(Q).function_index();
                int shell_Q_nfunc = primary_->shell(Q).nfunction();
                int shell_Q_offset = basis_endpoints_for_shell[Q] - basis_endpoints_for_shell[Qstart];
                for (const int S : Q_stripeout_list[dQ]) {
                    int shell_S_start = primary_->shell(S).function_index();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    for (int q = 0; q < shell_Q_nfunc; q++) {
                        for (int s = 0; s < shell_S_nfunc; s++) {
#pragma omp atomic
                            Kp[shell_Q_start + q][shell_S_start + s] += K3p[(q + shell_Q_offset) * nbf + s + shell_S_start];
#pragma omp atomic
                            Kp[shell_Q_start + q][shell_S_start + s] += K4p[(q + shell_Q_offset) * nbf + s + shell_S_start];
                        }
                    }

                }
            }

        }  // End stripe out

    }  // End master task list

    for (auto& Kmat : K) {
        Kmat->hermitivitize();
    }

    num_computed_shells_ = computed_shells;
}

}  // namespace psi
