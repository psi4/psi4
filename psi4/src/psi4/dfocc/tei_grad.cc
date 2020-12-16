/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

/** Standard library includes */
#include <fstream>
#include "psi4/libqt/qt.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libpsi4util/process.h"
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::tei_grad(std::string aux_type, std::map<std::string, SharedMatrix>& gradients) {
    //===========================================================================================
    //============================== Two-electron Gradient ======================================
    //===========================================================================================
    // Two-electron gradients
    int df_ints_num_threads_ = 1;
#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif

    std::string aux_name, intermed_name, intermed_short;
    if (aux_type == "JK") {
        aux_name = "SCF";
        intermed_name = "RefSep";
        intermed_short = "RefSep";
    } else if (aux_type == "RI") {
        aux_name = "CC";
        intermed_name = "Correlation";
        intermed_short = "Corr";
    } else {
        throw PSIEXCEPTION("Contact a developer. An unrecognized aux_type was passed to DFOCC::tei_grad.");
    }

    // Read in the basis set informations
    auto primary_ = get_basisset("ORBITAL");
    auto auxiliary_ = get_basisset("DF_BASIS_" + aux_name);
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    auto naux = auxiliary_->nbf();

    //===========================================================================================
    //============================== 3-Index Gradient ===========================================
    //===========================================================================================
    // Read gQso
    auto gQso = std::make_shared<Tensor2d>("" + intermed_name + " 3-Index TPDM (Q|nn)", naux, nso_, nso_);
    gQso->read(psio_, PSIF_DFOCC_DENS, true, true);

    // (Q | mu nu)^X
    auto idx3_short = "3-Index:" + intermed_short;
    timer_on("Grad: " + idx3_short);
    gradients[idx3_short] = matrix_factory()->create_shared_matrix(idx3_short + " Gradient", natom, 3);
    gradients[idx3_short]->zero();

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory2(
        new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri(1)));
    }

    const auto &shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //
    int max_rows;
    max_rows = auxiliary_->nshell();

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary_->nshell());

    // => Temporary Gradients <= //
    std::vector<SharedMatrix> Jtemps2;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jtemps2.push_back(std::make_shared<Matrix>("Jtemp2", natom, 3));
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P, 0, M, N);

            const double *buffer = eri[thread]->buffer();
            const auto& buffers = eri[thread]->buffers();

            int nP = auxiliary_->shell(P).nfunction();
            int cP = auxiliary_->shell(P).ncartesian();
            int aP = auxiliary_->shell(P).ncenter();
            // int oP = auxiliary_->shell(P).function_index() - pstart;
            int oP = auxiliary_->shell(P).function_index();

            int nM = primary_->shell(M).nfunction();
            int cM = primary_->shell(M).ncartesian();
            int aM = primary_->shell(M).ncenter();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int cN = primary_->shell(N).ncartesian();
            int aN = primary_->shell(N).ncenter();
            int oN = primary_->shell(N).function_index();

            int ncart = cP * cM * cN;
            const double *Px = buffers[0];
            const double *Py = buffers[1];
            const double *Pz = buffers[2];
            const double *Mx = buffers[3];
            const double *My = buffers[4];
            const double *Mz = buffers[5];
            const double *Nx = buffers[6];
            const double *Ny = buffers[7];
            const double *Nz = buffers[8];

            double perm = (M == N ? 1.0 : 2.0);

            double **grad_Jp;
            grad_Jp = Jtemps2[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        double Ival = 1.0 * perm * gQso->get(p + oP, (m + oM) * nso_ + (n + oN));
                        grad_Jp[aP][0] += Ival * (*Px);
                        grad_Jp[aP][1] += Ival * (*Py);
                        grad_Jp[aP][2] += Ival * (*Pz);
                        grad_Jp[aM][0] += Ival * (*Mx);
                        grad_Jp[aM][1] += Ival * (*My);
                        grad_Jp[aM][2] += Ival * (*Mz);
                        grad_Jp[aN][0] += Ival * (*Nx);
                        grad_Jp[aN][1] += Ival * (*Ny);
                        grad_Jp[aN][2] += Ival * (*Nz);

                        Px++;
                        Py++;
                        Pz++;
                        Mx++;
                        My++;
                        Mz++;
                        Nx++;
                        Ny++;
                        Nz++;
                    }
                }
            }
        }
    }

    // => Temporary Gradient Reduction <= //
    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients[idx3_short]->add(Jtemps2[t]);
    }

    timer_off("Grad: " + idx3_short);
}  // end

}  // namespace dfoccwave
}  // namespace psi

