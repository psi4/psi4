
/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "jk_grad.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace psi;

namespace psi {
namespace scfgrad {

DFJCOSKGrad::DFJCOSKGrad(int deriv, std::shared_ptr<MintsHelper> mints, Options& options) :
    JKGrad(deriv,mints->get_basisset("ORBITAL")), auxiliary_(mints->get_basisset("DF_BASIS_SCF")), mints_(mints), options_(options)
{
    common_init(); 
}
DFJCOSKGrad::~DFJCOSKGrad() {}
void DFJCOSKGrad::common_init() {
    ints_num_threads_ = 1;
#ifdef _OPENMP
    ints_num_threads_ = Process::environment.get_n_threads();
#endif

}
void DFJCOSKGrad::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> DFJCOSKGrad: Chain-of-Spheres SCF Gradients <==\n\n");

        outfile->Printf("    Gradient:          %11d\n", deriv_);
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", ints_num_threads_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    J Screening Type:   %11s\n", screen_type.c_str());
        outfile->Printf("    J Screening Cutoff: %11.0E\n", cutoff_);
        outfile->Printf("    K Screening Cutoff: %11.0E\n", options_.get_double("COSX_INTS_TOLERANCE"));
        outfile->Printf("    K Density Cutoff:   %11.0E\n", options_.get_double("COSX_DENSITY_TOLERANCE"));
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", options_.get_double("COSX_BASIS_TOLERANCE"));
        outfile->Printf("    K Overlap Fitting:  %11s\n", (options_.get_bool("COSX_OVERLAP_FITTING") ? "Yes" : "No"));
        outfile->Printf("\n");
    }
}
void DFJCOSKGrad::compute_gradient() {
    if (!do_J_ && !do_K_ && !do_wK_) return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");
    
    // => Set up gradients <= //
    int natom = primary_->molecule()->natom();
    gradients_.clear();
    if (do_J_) {
        gradients_["Coulomb"] = std::make_shared<Matrix>("Coulomb Gradient", natom, 3);
    }
    if (do_K_) {
        gradients_["Exchange"] = std::make_shared<Matrix>("Exchange Gradient", natom, 3);
    }
    if (do_wK_) {
        gradients_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Gradient", natom, 3);
    }

    if (do_J_) {
        build_JGrad();
        outfile->Printf("coulomb gradients successfully calculated.\n");
    }
    if (do_K_) {
        gradients_["Exchange"] = Process::environment.arrays["COSK_EXCHANGE_GRADIENT"];
        if (Ca_ == Cb_){
            gradients_["Exchange"]->scale(2.0);
        }

        outfile->Printf("exchange gradients successfully calculated.\n");
    }
    if (do_wK_) {
        throw PSIEXCEPTION("No wK for DFJCOSK.");
    }

    // Printing
    gradients_["Coulomb"]->print();
    if (do_K_) {
        gradients_["Exchange"]->print();
     }
}


void DFJCOSKGrad::build_JGrad() {
    if (!do_J_) return;
    timer_on("J Grad");
    int nthread_df = ints_num_threads_;

    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(ints_num_threads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)); //eri(1) ??
    for (int t = 1; t < ints_num_threads_; t++) {
        eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    const std::vector<std::pair<int, int>>& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();
    int nshell_aux = auxiliary_->nshell();
    int ntriplets = nshell_aux * npairs;

    // contract c_A
    // c_A = (B|pq) Dt_pq
    timer_on("contract c_A");

    SharedVector c;
    double* cp;

    c = std::make_shared<Vector>("c", naux);
    cp = c->pointer();
    std::vector<SharedVector> cT(nthread_df);
    for(size_t thread = 0; thread < nthread_df; thread++) {
        cT[thread] = std::make_shared<Vector>(naux);
    }


    double** Dtp = Dt_->pointer();


#pragma omp parallel for schedule(guided) num_threads(nthread_df)
    for (long int PMN = 0L; PMN < ntriplets; PMN++) {
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        int P = PMN / npairs;
        int MN = PMN % npairs;
        int M = shell_pairs[MN].first;
        int N = shell_pairs[MN].second;

        eri[thread]->compute_shell(P, 0, M, N);

        const auto & buffer = eri[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int Pstart = auxiliary_->shell(P).function_index();

        int nM = primary_->shell(M).nfunction();
        int Mstart = primary_->shell(M).function_index();

        int nN = primary_->shell(N).nfunction();
        int Nstart = primary_->shell(N).function_index();

        auto cTp = cT[thread]->pointer();

        for (int p = Pstart, index=0; p < Pstart + nP; p++) {
            for (int m = Mstart; m < Mstart + nM; m++) {
                for (int n = Nstart; n < Nstart + nN; n++, index++) {
                    cTp[p] += buffer[index] * Dtp[m][n];
                    if (N != M) cTp[p] += buffer[index] * Dtp[n][m];
                    
                }
            }
        }
    }

    timer_off("contract c_A");
    timer_on("contract d_A");
    // contract d_A
    // d_A = (A|B)^-1 c_A

    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_full_eig_inverse(condition_);
    SharedMatrix J = metric->get_metric();
    double** Jp = J->pointer();

    // add up per-thread contributions
    for(size_t thread = 0; thread < nthread_df; thread++) {
        c->add(*cT[thread]);
    }

    auto d = std::make_shared<Vector>("d", naux);
    double* dp = d->pointer();

    C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);


    timer_off("contract d_A");
    timer_on("metric grad");
    // metrix grad term
    // J = 0.5 * (A|B)^x d_B d_A

    std::map<std::string, SharedMatrix> densities;

    auto D = std::make_shared<Matrix>("D", naux, naux);
    auto Dp = D->pointer();
    C_DGER(naux, naux, 1, dp, 1, dp, 1, Dp[0], naux);
    densities["Coulomb"] = D;

    auto results = mints_->metric_grad(densities, "DF_BASIS_SCF");

    for (const auto& kv: results) {
        gradients_[kv.first] = kv.second;
    }
 

    timer_off("metric grad");
    timer_on("eri grad");
    // eri grad term
    // J += (A|pq)^x d_A Dt_pq

    int natom = primary_->molecule()->natom();
    std::vector<SharedMatrix> Jtemps;
    for (int t = 0; t < ints_num_threads_; t++) {
        Jtemps.push_back(std::make_shared<Matrix>("Jtemp", natom, 3));
    }

    // > Integrals < //
#pragma omp parallel for schedule(guided) num_threads(nthread_df)
    for (long int PMN = 0L; PMN < ntriplets; PMN++) {
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        int P = PMN / npairs;
        int MN = PMN % npairs;
        int M = shell_pairs[MN].first;
        int N = shell_pairs[MN].second;

        eri[thread]->compute_shell_deriv1(P, 0, M, N);

        const double* buffer = eri[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int cP = auxiliary_->shell(P).ncartesian();
        int aP = auxiliary_->shell(P).ncenter();
        int oP = auxiliary_->shell(P).function_index();

        int nM = primary_->shell(M).nfunction();
        int cM = primary_->shell(M).ncartesian();
        int aM = primary_->shell(M).ncenter();
        int oM = primary_->shell(M).function_index();

        int nN = primary_->shell(N).nfunction();
        int cN = primary_->shell(N).ncartesian();
        int aN = primary_->shell(N).ncenter();
        int oN = primary_->shell(N).function_index();

        const auto& buffers = eri[thread]->buffers();
        const double* Px = buffers[0];
        const double* Py = buffers[1];
        const double* Pz = buffers[2];
        const double* Mx = buffers[3];
        const double* My = buffers[4];
        const double* Mz = buffers[5];
        const double* Nx = buffers[6];
        const double* Ny = buffers[7];
        const double* Nz = buffers[8];

        double perm = (M == N ? 1.0 : 2.0);

        double** grad_Jp = Jtemps[thread]->pointer();

        for (int p = oP, index = 0; p < nP + oP; p++) {
            for (int m = oM; m < nM + oM; m++) {
                for (int n = oN; n < nN + oN; n++, index++) {
                    //  J^x = (A|pq)^x d_A Dt_pq
                    double Ival = 1.0 * perm * dp[p] * Dtp[m][n];
                    grad_Jp[aP][0] += Ival * Px[index];
                    grad_Jp[aP][1] += Ival * Py[index];
                    grad_Jp[aP][2] += Ival * Pz[index];
                    grad_Jp[aM][0] += Ival * Mx[index];
                    grad_Jp[aM][1] += Ival * My[index];
                    grad_Jp[aM][2] += Ival * Mz[index];
                    grad_Jp[aN][0] += Ival * Nx[index];
                    grad_Jp[aN][1] += Ival * Ny[index];
                    grad_Jp[aN][2] += Ival * Nz[index];
                }
            }
        }
    }

    timer_off("eri grad");
    // => Temporary Gradient Reduction <= //

    for (int t = 0; t < ints_num_threads_; t++) {
        gradients_["Coulomb"]->add(Jtemps[t]);
    }
    timer_off("J Grad");
}

void DFJCOSKGrad::compute_hessian(){
    throw PSIEXCEPTION("DFJCOSK analytic Hessians not implemented.");
}

}
}