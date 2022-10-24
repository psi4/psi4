/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

DFJLinK::DFJLinK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options) : JK(primary), auxiliary_(auxiliary), options_(options) { 
    timer_on("DFJLinK::Setup");
    common_init(); 
    timer_off("DFJLinK::Setup");
}

DFJLinK::~DFJLinK() {}

void DFJLinK::common_init() {

    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif

    incfock_ = options_.get_bool("INCFOCK");
    incfock_count_ = 0;
    do_incfock_iter_ = false;
    if (options_.get_int("INCFOCK_FULL_FOCK_EVERY") <= 0) {
        throw PSIEXCEPTION("Invalid input for option INCFOCK_FULL_FOCK_EVERY (<= 0)");
    }
    density_screening_ = options_.get_str("SCREENING") == "DENSITY";
   
    set_cutoff(options_.get_double("INTS_TOLERANCE"));

    // => Direct Density-Fitted Coulomb Setup <= //

    // pre-compute coulomb fitting metric
    timer_on("DFJLinK: Coulomb Metric");
    FittingMetric J_metric_obj(auxiliary_, true);
    J_metric_obj.form_fitting_metric();
    J_metric_ = J_metric_obj.get_metric();
    timer_off("DFJLinK: Coulomb Metric");

    // pre-construct per-thread TwoBodyAOInt objects for computing 3- and 4-index ERIs
    timer_on("DFJLinK: ERI Computers");
    eri_computers_["4-Center"].emplace({}); 
    eri_computers_["3-Center"].emplace({});
    
    eri_computers_["4-Center"].resize(nthreads_);
    eri_computers_["3-Center"].resize(nthreads_);
    
    auto zero = BasisSet::zero_ao_basis_set();
    IntegralFactory rifactory(auxiliary_, zero, primary_, primary_);
    IntegralFactory factory(primary_, primary_, primary_, primary_);
    eri_computers_["4-Center"][0] = std::shared_ptr<TwoBodyAOInt>(factory.eri());
    eri_computers_["3-Center"][0] = std::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    for(int rank = 1; rank < nthreads_; rank++) {
        eri_computers_["4-Center"][rank] = std::shared_ptr<TwoBodyAOInt>(eri_computers_["4-Center"].front()->clone());
        eri_computers_["3-Center"][rank] = std::shared_ptr<TwoBodyAOInt>(eri_computers_["3-Center"].front()->clone());
    }
    timer_off("DFJLinK: ERI Computers");

    // => Linear Exchange Setup <= //
    
    // set up LinK integral tolerance
    if (options_["LINK_INTS_TOLERANCE"].has_changed()) {
        linK_ints_cutoff_ = options_.get_double("LINK_INTS_TOLERANCE");
    } else {
        linK_ints_cutoff_ = options_.get_double("INTS_TOLERANCE");
    }
}
size_t DFJLinK::num_computed_shells() { 
    //no bench data returned - to come in a future update
    return JK::num_computed_shells(); 
}

size_t DFJLinK::memory_estimate() {
    return 0;  // Memory is O(N^2), which psi4 counts as effectively 0
}

void DFJLinK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> DFJLinK: Density-Fitted J and Linear Exchange K <==\n\n");

        outfile->Printf("    J tasked:           %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:           %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:          %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:              %11.3E\n", omega_);
        outfile->Printf("    Integrals threads:  %11d\n", nthreads_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Incremental Fock :  %11s\n", (incfock_ ? "Yes" : "No"));
        outfile->Printf("    Screening Type:   %11s\n", screen_type.c_str());
        outfile->Printf("    Screening Cutoff: %11.0E\n", cutoff_);
    }
}

void DFJLinK::preiterations() {}

void DFJLinK::incfock_setup() {

    // The prev_D_ao_ condition is used to handle stability analysis case
    if (initial_iteration_ || prev_D_ao_.size() != D_ao_.size()) {
        initial_iteration_ = true;

        prev_D_ao_.clear();
        delta_D_ao_.clear();

        if (do_wK_) {
            prev_wK_ao_.clear();
            delta_wK_ao_.clear();
        }

        if (do_J_) {
            prev_J_ao_.clear();
            delta_J_ao_.clear();
        }

        if (do_K_) {
            prev_K_ao_.clear();
            delta_K_ao_.clear();
        }
    
        for (size_t N = 0; N < D_ao_.size(); N++) {
            prev_D_ao_.push_back(std::make_shared<Matrix>("D Prev", D_ao_[N]->nrow(), D_ao_[N]->ncol()));
            delta_D_ao_.push_back(std::make_shared<Matrix>("Delta D", D_ao_[N]->nrow(), D_ao_[N]->ncol()));

            if (do_wK_) {
                prev_wK_ao_.push_back(std::make_shared<Matrix>("wK Prev", wK_ao_[N]->nrow(), wK_ao_[N]->ncol()));
                delta_wK_ao_.push_back(std::make_shared<Matrix>("Delta wK", wK_ao_[N]->nrow(), wK_ao_[N]->ncol()));
            }
                
            if (do_J_) {
                prev_J_ao_.push_back(std::make_shared<Matrix>("J Prev", J_ao_[N]->nrow(), J_ao_[N]->ncol()));
                delta_J_ao_.push_back(std::make_shared<Matrix>("Delta J", J_ao_[N]->nrow(), J_ao_[N]->ncol()));
            }
        
            if (do_K_) {
                prev_K_ao_.push_back(std::make_shared<Matrix>("K Prev", K_ao_[N]->nrow(), K_ao_[N]->ncol()));
                delta_K_ao_.push_back(std::make_shared<Matrix>("Delta K", K_ao_[N]->nrow(), K_ao_[N]->ncol()));
            }
        }
    } else {
        for (size_t N = 0; N < D_ao_.size(); N++) {
            delta_D_ao_[N]->copy(D_ao_[N]);
            delta_D_ao_[N]->subtract(prev_D_ao_[N]);
        }
    }
}
void DFJLinK::incfock_postiter() {
    if (do_incfock_iter_) {
        for (size_t N = 0; N < D_ao_.size(); N++) {

            if (do_wK_) {
                prev_wK_ao_[N]->add(delta_wK_ao_[N]);
                wK_ao_[N]->copy(prev_wK_ao_[N]);
            }

            if (do_J_) {
                prev_J_ao_[N]->add(delta_J_ao_[N]);
                J_ao_[N]->copy(prev_J_ao_[N]);
            }

            if (do_K_) {
                prev_K_ao_[N]->add(delta_K_ao_[N]);
                K_ao_[N]->copy(prev_K_ao_[N]);
            }

            prev_D_ao_[N]->copy(D_ao_[N]);
        }
    } else {
        for (size_t N = 0; N < D_ao_.size(); N++) {
            if (do_wK_) prev_wK_ao_[N]->copy(wK_ao_[N]);
            if (do_J_) prev_J_ao_[N]->copy(J_ao_[N]);
            if (do_K_) prev_K_ao_[N]->copy(K_ao_[N]);
            prev_D_ao_[N]->copy(D_ao_[N]);
        }
    }
}

void DFJLinK::compute_JK() {
 
    int njk = D_ao_.size();

    if (incfock_) {
        timer_on("DFJLinK: INCFOCK Preprocessing");
        incfock_setup();
        int reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        double incfock_conv = options_.get_double("INCFOCK_CONVERGENCE");
        double Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= incfock_conv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);
        
        if (!initial_iteration_ && (Dnorm >= incfock_conv)) incfock_count_ += 1;
        timer_off("DFJLinK: INCFOCK Preprocessing");
    }
    
    std::vector<SharedMatrix>& D_ref = (do_incfock_iter_ ? delta_D_ao_ : D_ao_);
    std::vector<SharedMatrix>& J_ref = (do_incfock_iter_ ? delta_J_ao_ : J_ao_);
    std::vector<SharedMatrix>& K_ref = (do_incfock_iter_ ? delta_K_ao_ : K_ao_);
    std::vector<SharedMatrix>& wK_ref = (do_incfock_iter_ ? delta_wK_ao_ : wK_ao_);

    if (density_screening_) {
        for (auto eri_computer : eri_computers_["4-Center"]) {
            eri_computer->update_density(D_ref);
	}
    }


    //if (density_screening_) {
    //    eri_computers_["4-Center"][0]->update_density(D_ref);
    //    for (int thread = 1; thread < nthreads_; thread++) {
//	    eri_computers_["4-Center"][thread] = std::shared_ptr<TwoBodyAOInt>(eri_computers_["4-Center"][0]->clone());
	//}
	//eri_computers_["3-Center"][0]->update_density(D_ref);
    //}

    if (do_wK_) throw PSIEXCEPTION("DFJLinK does not support wK integrals yet!");

    // D_eff, the effective pseudo-density matrix is either:
    //   (1) the regular density: D_eff == D_lr = C_lo x C*ro
    //   (2) the difference density: D_eff == dD_lr = (C_lo x C_ro)_{iter} - (C_lo x C_ro)_{iter - 1}
    //
    std::vector<SharedMatrix> D_eff(njk);

    if (do_J_) {
        timer_on("DFJLinK: J");
        //build_J(D_eff, J_ao_);
        for (auto& Jmat : J_ref) Jmat->zero();
	build_J(D_ref, J_ref);
        timer_off("DFJLinK: J");
    }
    
    if (do_K_) {
        timer_on("DFJLinK: K");
        //build_K(D_eff, K_ao_);
        for (auto& Kmat : K_ref) Kmat->zero();
        build_K(D_ref, K_ref);
        timer_off("DFJLinK: K");
    }
    
    if (incfock_) {
        timer_on("DFJLinK: INCFOCK Postprocessing");
        incfock_postiter();
        timer_off("DFJLinK: INCFOCK Postprocessing");
    }

    if (initial_iteration_) initial_iteration_ = false;
}

void DFJLinK::postiterations() {}

void DFJLinK::build_J(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& J) {
    
    timer_on("Setup");

    // => Sizing <= //
    int njk = D.size();
    int nbf = primary_->nbf();
    int nshell = primary_->nshell();
    int nbf_aux = auxiliary_->nbf();
    int nshell_aux = auxiliary_->nshell();

    // benchmarking 
    size_t nshellpair = eri_computers_["3-Center"][0]->shell_pairs().size();
    size_t nshelltriplet = nshell_aux * nshellpair;
    size_t computed_triplets1 = 0, computed_triplets2 = 0;

    // screening threshold
    double thresh2 = options_.get_double("INTS_TOLERANCE") * options_.get_double("INTS_TOLERANCE");

    // per-thread G Vector buffers (for accumulating thread contributions to G)
    // G is the contraction of the density matrix with the 3-index ERIs
    std::vector<std::vector<SharedVector>> GT(njk, std::vector<SharedVector>(nthreads_));

    // H is the contraction of G with the inverse coulomb metric
    std::vector<SharedVector> H(njk);

    // per-thread J Matrix buffers (for accumulating thread contributions to J)
    std::vector<std::vector<SharedMatrix>> JT(njk, std::vector<SharedMatrix>(nthreads_));

    // initialize per-thread objects
    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthreads_; thread++) {
            JT[jki][thread] = std::make_shared<Matrix>(nbf, nbf);
            GT[jki][thread] = std::make_shared<Vector>(nbf_aux);
        }
        H[jki] = std::make_shared<Vector>(nbf_aux);
    }

    // diagonal shell maxima of J_metric_ for screening
    std::vector<double> J_metric_shell_diag(nshell_aux, 0.0);
    for (size_t s = 0; s < nshell_aux; s++) {
        int bf_start = auxiliary_->shell(s).function_index();
        int bf_end = bf_start + auxiliary_->shell(s).nfunction();
        for (size_t bf = bf_start; bf < bf_end; bf++) {
            J_metric_shell_diag[s] = std::max(J_metric_shell_diag[s], J_metric_->get(bf, bf));
        }
    }

    // shell maxima of D for screening
    Matrix Dshell(nshell, nshell);
    auto Dshellp = Dshell.pointer();

    for(size_t M = 0; M < nshell; M++) {
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        for(size_t N = 0; N < nshell; N++) {
            int nn = primary_->shell(N).nfunction();
            int nstart = primary_->shell(N).function_index();
            for(size_t jki = 0; jki < njk; jki++) {
                auto Dp = D[jki]->pointer();
                for(size_t m = mstart; m < mstart + nm; m++) {
                    for(size_t n = nstart; n < nstart + nn; n++) {
                        Dshellp[M][N] = std::max(Dshellp[M][N], std::abs(Dp[m][n]));
                    }
                }
            }
        }
    }

    timer_off("Setup");

    //  => First Contraction <= //

    // contract D with three-index DF ERIs to get G:
    // G_{p} = D_{mn} * (mn|p)

    timer_on("ERI1");
#pragma omp parallel for schedule(guided) num_threads(nthreads_) reduction(+ : computed_triplets1)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = eri_computers_["3-Center"][rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(Dshellp[M][N] * Dshellp[M][N] * J_metric_shell_diag[P] * eri_computers_["3-Center"][rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets1++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();
        eri_computers_["3-Center"][rank]->compute_shell(P, 0, M, N);
        const auto & buffer = eri_computers_["3-Center"][rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto GTp = GT[jki][rank]->pointer();
            auto Dp = D[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        GTp[p] += buffer[index] * Dp[m][n];
                        if (N != M) GTp[p] += buffer[index] * Dp[n][m];
                    }
                }
            }

        }

    }

    timer_off("ERI1");

    //  => Second Contraction <= //

    //  linear solve for H:
    //  G_{p} = H_{q} (q|p)

    timer_on("Metric");

    std::vector<int> ipiv(nbf_aux);

    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthreads_; thread++) {
            H[jki]->add(*GT[jki][thread]);
        }
        C_DGESV(nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux);
    }


    // I believe C_DSYSV should be faster than C_GESV, but I've found the opposite to be true.
    // This performance issue should be investigated, but is not consequential here.
    // The cost of either linear solve is dwarfed by the actual integral computation.
    //
    //std::vector<double> work(3 * nbf_aux);
    //int errcode = C_DSYSV('U', nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux, work.data(), 3 * nbf_aux);

    // shell maxima of H for screening
    Vector H_shell_max(nshell_aux);
    auto H_shell_maxp = H_shell_max.pointer();

    for(size_t jki = 0; jki < njk; jki++) {

        auto Hp = H[jki]->pointer();

        for (int P = 0; P < nshell_aux; P++) {
            int np = auxiliary_->shell(P).nfunction();
            int pstart = auxiliary_->shell(P).function_index();
            for (int p = pstart; p < pstart + np; p++) {
                H_shell_maxp[P] = std::max(H_shell_maxp[P], std::abs(Hp[p]));
            }
        }

    }

    timer_off("Metric");

    //  => Third Contraction <= //

    // contract H with three-index DF ERIs to get J
    // J_{mn} = H_{p} (mn|p)

    timer_on("ERI2");

#pragma omp parallel for schedule(guided) num_threads(nthreads_) reduction(+ : computed_triplets2)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = eri_computers_["3-Center"][rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(H_shell_maxp[P] * H_shell_maxp[P] * J_metric_shell_diag[P] * eri_computers_["3-Center"][rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets2++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();

        eri_computers_["3-Center"][rank]->compute_shell(P, 0, M, N);
        const auto & buffer = eri_computers_["3-Center"][rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto JTp = JT[jki][rank]->pointer();
            auto Hp = H[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        JTp[m][n] += buffer[index] * Hp[p];
                        if (N != M) JTp[n][m] += buffer[index] * Hp[p];
                    }
                }
            }

        }

    }

    timer_off("ERI2");

    if (bench_) {
        auto mode = std::ostream::app;
        PsiOutStream printer("bench.dat", mode);
        printer.Printf("DFJ ERI Shells: %zu,%zu,%zu\n", computed_triplets1, computed_triplets2, nshelltriplet);
    }

    for(size_t jki = 0; jki < njk; jki++) {
        for (size_t thread = 0; thread < nthreads_; thread++) {
            J[jki]->add(JT[jki][thread]);
        }
        J[jki]->hermitivitize();
    }

}

// To follow this code, compare with figure 1 of DOI: 10.1063/1.476741
void DFJLinK::build_K(std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& K) {

    if (!lr_symmetric_) {
        throw PSIEXCEPTION("Non-symmetric K matrix builds are currently not supported in the LinK algorithm.");
    }

    timer_on("build_linK()");

    // ==> Prep Auxiliary Quantities <== //

    // => Zeroing <= //
    for (auto& Kmat : K) {
        Kmat->zero();
    }

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
                    if (eri_computers_["4-Center"][0]->shell_pair_significant(P, Q)) {
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
    double max_integral = eri_computers_["4-Center"][0]->max_integral();

#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PQ_shell_values;
        for (size_t Q = 0; Q < nshell; Q++) {
            double pq_pq = std::sqrt(eri_computers_["4-Center"][0]->shell_ceiling2(P, Q, P, Q));
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
            double val = std::sqrt(eri_computers_["4-Center"][0]->shell_ceiling2(P, Q, P, Q));
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
            double screen_val = shell_ceilings[P] * shell_ceilings[R] * eri_computers_["4-Center"][0]->shell_pair_max_density(P, R);
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
                if (!eri_computers_["4-Center"][0]->shell_pair_significant(P, Q)) continue;

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
                        double screen_val = eri_computers_["4-Center"][0]->shell_pair_max_density(P, R) * std::sqrt(eri_computers_["4-Center"][0]->shell_ceiling2(P, Q, R, S));

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
                        double screen_val = eri_computers_["4-Center"][0]->shell_pair_max_density(Q, R) * std::sqrt(eri_computers_["4-Center"][0]->shell_ceiling2(P, Q, R, S));

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

                    if (!eri_computers_["4-Center"][0]->shell_pair_significant(R, S)) continue;
                    if (!eri_computers_["4-Center"][0]->shell_significant(P, Q, R, S)) continue;

                    if (eri_computers_["4-Center"][thread]->compute_shell(P, Q, R, S) == 0)
                        continue;
                    computed_shells++;

                    const double* buffer = eri_computers_["4-Center"][thread]->buffer();

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
    timer_off("build_linK()");
}

}  // namespace psi
