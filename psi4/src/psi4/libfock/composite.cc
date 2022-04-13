#include "composite.h"

#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <vector>
#include <algorithm>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

SplitJKBase::SplitJKBase(std::shared_ptr<BasisSet> primary, Options& options) : primary_(primary), options_(options) {

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

}

CompositeJK::CompositeJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, 
                         std::string jtype, std::string ktype, Options& options)
                        : JK(primary), auxiliary_(auxiliary), jtype_(jtype), ktype_(ktype), options_(options) { common_init(); }

void CompositeJK::common_init() {

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    if (jtype_ == "DIRECT_DF") {
        jalgo_ = std::make_shared<DirectDFJ>(primary_, auxiliary_, options_);
    } else {
        throw PSIEXCEPTION("J BUILD TYPE " + jtype_ + " IS NOT SUPPORTED IN COMPOSITE JK!");
    }

    if (ktype_ == "LINK") {
        kalgo_ = std::make_shared<LinK>(primary_, options_);
    } else {
        throw PSIEXCEPTION("K BUILD TYPE " + ktype_ + " IS NOT SUPPORTED IN COMPOSITE JK!");
    }
}

size_t CompositeJK::memory_estimate() {
    return 0;   // only Direct-based integral algorithms are currently implemented in Composite JK
}

void CompositeJK::preiterations() { return; }
void CompositeJK::postiterations() { return; }

void CompositeJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> CompositeJK: Mix/Match JK Builds <==\n\n");
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        if (do_J_) outfile->Printf("    J Algorithm:       %11s\n", jtype_.c_str());
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        if (do_K_) outfile->Printf("    K Algorithm:       %11s\n", ktype_.c_str());
        outfile->Printf("    Integrals threads: %11d\n", nthread_);
        outfile->Printf("    Incremental Fock:  %11s\n", incfock_ ? "Yes" : "No");
        outfile->Printf("\n");
    }
    jalgo_->print_header();
    kalgo_->print_header();
}

void CompositeJK::compute_JK() {

    if (!lr_symmetric_) {
        throw PSIEXCEPTION("Non-symmetric K matrix builds are currently not supported in Composite JK.");
    }

    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Preprocessing");
        incfock_setup();
        int reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        double dconv = options_.get_double("D_CONVERGENCE");
        double Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= dconv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);
        
        if (!initial_iteration_ && (Dnorm >= dconv)) incfock_count_ += 1;
        timer_off("CompositeJK: INCFOCK Preprocessing");
    }

    // Matrices to use/build depending on whether or not incremental Fock build is performed in the iteration
    std::vector<SharedMatrix>& D_ref = (do_incfock_iter_ ? delta_D_ao_ : D_ao_);
    std::vector<SharedMatrix>& J_ref = (do_incfock_iter_ ? delta_J_ao_ : J_ao_);
    std::vector<SharedMatrix>& K_ref = (do_incfock_iter_ ? delta_K_ao_ : K_ao_);

    jalgo_->build_G_component(D_ref, J_ref);
    kalgo_->build_G_component(D_ref, K_ref);

    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Postprocessing");
        incfock_postiter();
        timer_off("CompositeJK: INCFOCK Postprocessing");
    }

    if (initial_iteration_) initial_iteration_ = false;
}

DirectDFJ::DirectDFJ(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
                        : SplitJKBase(primary, options), auxiliary_(auxiliary) {
    cutoff_ = options_.get_double("INTS_TOLERANCE");
    build_metric();
    build_ints();
}

void DirectDFJ::print_header() {
    if (print_) {
        outfile->Printf("  ==> Direct Density-Fitted J <==\n\n");
        outfile->Printf("    Primary Basis: %11s\n", primary_->name().c_str());
        outfile->Printf("    Auxiliary Basis: %11s\n", auxiliary_->name().c_str());
        outfile->Printf("    ERI Screening Cutoff: %11.0E\n", cutoff_);
        outfile->Printf("\n");
    }
}

void DirectDFJ::build_metric() {
    timer_on("DirectDFJ: Build Metric");

    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_fitting_metric();
    Jmet_ = metric->get_metric();

    timer_off("DirectDFJ: Build Metric");
}

void DirectDFJ::build_ints() {
    timer_on("DirectDFJ: Build Ints");

    auto zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);
    ints_.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
    for (int thread = 1; thread < nthread_; thread++) {
        ints_.push_back(std::shared_ptr<TwoBodyAOInt>(ints_[0]->clone()));
    }

    timer_off("DirectDFJ: Build Ints");
}

void DirectDFJ::build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {

    timer_on("DirectDFJ: J");

    // => Zeroing <= //

    for (auto& Jmat : J) {
        Jmat->zero();
    }

    // => Sizing <= //
    int pri_nshell = primary_->nshell();
    int aux_nshell = auxiliary_->nshell();
    int nmat = D.size();

    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Get significant primary shells <=
    const auto& shell_pairs = ints_[0]->shell_pairs();

    // Weigand 2002 doi: 10.1039/b204199p (Figure 1)
    // The gamma intermediates defined in Figure 1
    // gammaP = (P|uv) * Duv
    // (P|Q) * gammaQ = gammaP
    SharedMatrix gamma = std::make_shared<Matrix>(nmat, naux);
    gamma->zero();
    double* gamp = gamma->pointer()[0];

    int max_nbf_per_shell = 0;
    for (int P = 0; P < pri_nshell; P++) {
        max_nbf_per_shell = std::max(max_nbf_per_shell, primary_->shell(P).nfunction());
    }

    // Temporary buffers for J to minimize race tonditions
    std::vector<std::vector<SharedMatrix>> JT;

    for (int thread = 0; thread < nthread_; thread++) {
        std::vector<SharedVector> gp2;
        std::vector<SharedMatrix> J2;
        for (size_t ind = 0; ind < nmat; ind++) {
            J2.push_back(std::make_shared<Matrix>(max_nbf_per_shell, max_nbf_per_shell));
            J2[ind]->zero();
        }
        JT.push_back(J2);
    }

    // Solve for gammaP = (P|rs)*Drs
#pragma omp parallel for num_threads(nthread_) schedule(dynamic)
    for (int P = 0; P < aux_nshell; P++) {

        int p_start = auxiliary_->shell_to_basis_function(P);
        int num_p = auxiliary_->shell(P).nfunction();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (const auto& UV : shell_pairs) {

            int U = UV.first;
            int V = UV.second;

            int u_start = primary_->shell_to_basis_function(U);
            int num_u = primary_->shell(U).nfunction();

            int v_start = primary_->shell_to_basis_function(V);
            int num_v = primary_->shell(V).nfunction();

            ints_[thread]->compute_shell(P, 0, U, V);

            const double* buffer = ints_[thread]->buffer();

            double prefactor = 2.0;
            if (U == V) prefactor *= 0.5;

            for (int i = 0; i < D.size(); i++) {
                double** Dp = D[i]->pointer();
                double *Puv = const_cast<double *>(buffer);

                /*
                for (int p = p_start; p < p_start + num_p; p++) {
                    int dp = p - p_start;
                    for (int u = u_start; u < u_start + num_u; u++) {
                        int du = u - u_start;
                        for (int v = v_start; v < v_start + num_v; v++) {
                            int dv = v - v_start;
                            gamp[i * naux + p] += prefactor * (*Puv) * Dp[u][v];
                            Puv++;
                        }
                    }
                }
                */

                std::vector<double> Dbuff(num_u * num_v, 0.0);
                for (int u = u_start; u < u_start + num_u; u++) {
                    int du = u - u_start;
                    for (int v = v_start; v < v_start + num_v; v++) {
                        int dv = v - v_start;
                        Dbuff[du * num_v + dv] = Dp[u][v];
                    }
                }
                C_DGEMV('N', num_p, num_u * num_v, prefactor, Puv, num_u * num_v, Dbuff.data(), 1, 1.0, &(gamp[i * naux + p_start]), 1);
            }
        }
    }

    // Solve for gammaQ, (P|Q) * gammaQ = gammaP
    SharedMatrix Jmet_copy = Jmet_->clone();

    std::vector<int> ipiv(naux);
    C_DGESV(naux, nmat, Jmet_copy->pointer()[0], naux, ipiv.data(), gamp, naux);

#pragma omp parallel for num_threads(nthread_) schedule(dynamic)
    for (int UVidx = 0; UVidx < shell_pairs.size(); UVidx++) {

        int U = shell_pairs[UVidx].first;
        int V = shell_pairs[UVidx].second;

        int u_start = primary_->shell_to_basis_function(U);
        int num_u = primary_->shell(U).nfunction();

        int v_start = primary_->shell_to_basis_function(V);
        int num_v = primary_->shell(V).nfunction();

        double prefactor = 2.0;
        if (U == V) prefactor *= 0.5;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int Q = 0; Q < aux_nshell; Q++) {

            int q_start = auxiliary_->shell_to_basis_function(Q);
            int num_q = auxiliary_->shell(Q).nfunction();

            ints_[thread]->compute_shell(Q, 0, U, V);

            const double* buffer = ints_[thread]->buffer();

            for (int i = 0; i < D.size(); i++) {
                double* JTp = JT[thread][i]->pointer()[0];
                double* Quv = const_cast<double *>(buffer);

                /*
                for (int q = q_start; q < q_start + num_q; q++) {
                    int dq = q - q_start;
                    for (int u = u_start; u < u_start + num_u; u++) {
                        int du = u - u_start;
                        for (int v = v_start; v < v_start + num_v; v++) {
                            int dv = v - v_start;
                            JTp[du * num_v + dv] += prefactor * (*Quv) * gamp[i * naux + q];
                            Quv++;
                        }
                    }
                }
                */
                C_DGEMV('T', num_q, num_u * num_v, prefactor, Quv, num_u * num_v, &(gamp[i * naux + q_start]), 1, 1.0, JTp, 1);
            }
        }

        // => Stripeout <= //

        for (int i = 0; i < D.size(); i++) {
            double* JTp = JT[thread][i]->pointer()[0];
            double** Jp = J[i]->pointer();
            for (int u = u_start; u < u_start + num_u; u++) {
                int du = u - u_start;
                for (int v = v_start; v < v_start + num_v; v++) {
                    int dv = v - v_start;
                    Jp[u][v] += JTp[du * num_v + dv];
                }
            }
            JT[thread][i]->zero();
        }
    }

    for (auto& Jmat : J) {
        Jmat->hermitivitize();
    }

    timer_off("DirectDFJ: J");
}

LinK::LinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJKBase(primary, options) {
    cutoff_ = options.get_double("INTS_TOLERANCE");
    if (options_["LINK_INTS_TOLERANCE"].has_changed()) {
        linK_ints_cutoff_ = options_.get_double("LINK_INTS_TOLERANCE");
    } else {
        linK_ints_cutoff_ = options_.get_double("INTS_TOLERANCE");
    }

    build_ints();
}

void LinK::print_header() {
    if (print_) {
        outfile->Printf("  ==> Linear Exchange (LinK) <==\n\n");
        outfile->Printf("    Primary Basis: %11s\n", primary_->name().c_str());
        outfile->Printf("    ERI Screening Cutoff: %11.0E\n", cutoff_);
        outfile->Printf("    Density Screening Cutoff: %11.0E\n", linK_ints_cutoff_);
        outfile->Printf("\n");
    }
}

void LinK::build_ints() {
    timer_on("LinK: Build Ints");

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    ints_.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    for (int thread = 1; thread < nthread_; thread++) {
        ints_.push_back(std::shared_ptr<TwoBodyAOInt>(ints_[0]->clone()));
    }

    timer_off("LinK: Build Ints");
}

// To follow this code, compare with figure 1 of DOI: 10.1063/1.476741
void LinK::build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& K) {

    timer_on("LinK: K");

    for (auto& integral : ints_) {
        integral->update_density(D);
    }
    // ==> Prep Auxiliary Quantities <== //

    // => Zeroing <= //
    for (auto& Kmat : K) {
        Kmat->zero();
    }

    // => Sizing <= //
    int nshell = primary_->nshell();
    int nbf = primary_->nbf();
    int nthread = nthread_;

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
                    if (ints_[0]->shell_pair_significant(P, Q)) {
                        found = true;
                        atom_pairs.emplace_back(Patom, Qatom);
                        break;
                    }
                }
                if (found) break;
            }
        }
    }

    // ==> Prep Bra-Bra Shell Pairs <== //

    // A comparator used for sorting integral screening values
    auto screen_compare = [](const std::pair<int, double> &a, 
                                    const std::pair<int, double> &b) { return a.second > b.second; };

    std::vector<std::vector<int>> significant_bras(nshell);
    double max_integral = ints_[0]->max_integral();

#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PQ_shell_values;
        for (size_t Q = 0; Q < nshell; Q++) {
            double pq_pq = std::sqrt(ints_[0]->shell_ceiling2(P, Q, P, Q));
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
            double val = std::sqrt(ints_[0]->shell_ceiling2(P, Q, P, Q));
            shell_ceilings[P] = std::max(shell_ceilings[P], val);
#pragma omp critical
            shell_ceilings[Q] = std::max(shell_ceilings[Q], val);
        }
    }

    std::vector<std::vector<int>> significant_kets(nshell);

    // => Use shell ceilings to compute significant ket-shells for each bra-shell <= //
#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PR_shell_values;
        for (size_t R = 0; R < nshell; R++) {
            double screen_val = shell_ceilings[P] * shell_ceilings[R] * ints_[0]->shell_pair_max_density(P, R);
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

    // Number of computed shell quartets is tracked for benchmarking purposes
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
                if (!ints_[0]->shell_pair_significant(P, Q)) continue;

                int dP = P - Pstart;
                int dQ = Q - Qstart;

                // => Formation of Significant Shell Pair List ML <= //

                // Significant ket shell pairs RS for bra shell pair PQ
                // represents the merge of ML_P and ML_Q (mini-lists) as defined in Oschenfeld
                // Unordered set structure allows for automatic merging as new elements are added
                std::unordered_set<int> ML_PQ;

                // Form ML_P as part of ML_PQ
                for (const int R : significant_kets[P]) {
                    bool is_significant = false;
                    for (const int S : significant_bras[R]) {
                        double screen_val = ints_[0]->shell_pair_max_density(P, R) * std::sqrt(ints_[0]->shell_ceiling2(P, Q, R, S));

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
                        double screen_val = ints_[0]->shell_pair_max_density(Q, R) * std::sqrt(ints_[0]->shell_ceiling2(P, Q, R, S));

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

                    if (!ints_[0]->shell_pair_significant(R, S)) continue;
                    if (!ints_[0]->shell_significant(P, Q, R, S)) continue;

                    if (ints_[thread]->compute_shell(P, Q, R, S) == 0)
                        continue;
                    computed_shells++;

                    const double* buffer = ints_[thread]->buffer();

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
        Kmat->scale(2.0);
        Kmat->hermitivitize();
    }

    if (bench_) {
        auto mode = std::ostream::app;
        auto printer = PsiOutStream("bench.dat", mode);
        size_t ntri = nshell * (nshell + 1L) / 2L;
        size_t possible_shells = ntri * (ntri + 1L) / 2L;
        printer.Printf("(LinK) Computed %20zu Shell Quartets out of %20zu, (%11.3E ratio)\n", computed_shells,
                        possible_shells, computed_shells / (double)possible_shells);
    }

    timer_off("LinK: K");
}

}