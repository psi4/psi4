/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sieve.h"
#include "dfocc.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/mintshelper.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ref_grad() {
    // outfile->Printf("\tref_grad is starting... \n");
    /********************************************************************************************/
    /************************** Build Intermediates and TPDM ************************************/
    /********************************************************************************************/
    timer_on("DF-SCF integrals");
    df_ref();
    fock_so();
    timer_off("DF-SCF integrals");

    // Build SO-Basis OPDM
    Dso = SharedTensor2d(new Tensor2d("SO-basis Density Matrix", nso_, nso_));
    Dso->gemm(false, true, CoccA, CoccA, 2.0, 0.0);

    // Build TPDM OO Blok in the MO Basis
    // G_ij^Q = 4*\delta_ij \sum_{k} c_kk^Q - 2*c_ij^Q
    gQoo_ref = SharedTensor2d(new Tensor2d("DF_BASIS_SCF 3-Index TPDM <O|O>", nQ_ref, naoccA * naoccA));
    gQoo_ref->copy(cQooA);
    gQoo_ref->scale(-2.0);
    for (int Q = 0; Q < nQ_ref; Q++) {
        for (int i = 0; i < noccA; ++i) {
            double summ = 0.0;
            int ii = oo_pair_idxAA(i, i);
            for (int j = 0; j < noccA; ++j) {
                int jj = oo_pair_idxAA(j, j);
                summ += cQooA->get(Q, jj);
            }
            gQoo_ref->add(Q, ii, 4.0 * summ);
        }
    }

    // Backtransform MO basis 3-Index TPDM to SO-basis 3-index TPDM
    gQso_ref = SharedTensor2d(new Tensor2d("3-Index TPDM", nQ_ref, nso2_));
    gQon_ref = SharedTensor2d(new Tensor2d("DF_BASIS_SCF G_imu^Q", nQ_ref, nso_ * noccA));
    // G_im^Q = \sum_{j} G_ij^Q * Cnj
    gQon_ref->contract(false, true, nQ_ref * noccA, nso_, noccA, gQoo_ref, CoccA, 1.0, 0.0);
    // G_mn^Q = \sum_{i} Cmi * G_in^Q
    gQso_ref->contract233(false, false, nso_, nso_, CoccA, gQon_ref, 1.0, 0.0);
    // gQso_ref->print();

    // Build G(P,Q) : 2-Index TPDM
    Gaux_ref = SharedTensor2d(new Tensor2d("2-Index TPDM", nQ_ref, nQ_ref));
    // Gaux_ref->gemm(false, true, cQso, gQso_ref, 0.5, 0.0); // SO basis
    Gaux_ref->gemm(false, true, cQooA, gQoo_ref, 0.5, 0.0);  // MO basis
    // Gaux_ref->print();

    // Build Wmn = 2*\sum_{i} e_i * Cmi Cni
    for (int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
    for (int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
    Wso = SharedTensor2d(new Tensor2d("SO-basis GFM", nso_, nso_));
    for (int mu = 0; mu < nso_; mu++) {
        for (int nu = 0; nu < nso_; nu++) {
            double summ = 0.0;
            for (int i = 0; i < noccA; i++) {
                summ += 2.0 * FockA->get(i, i) * CmoA->get(mu, i) * CmoA->get(nu, i);
            }
            Wso->set(mu, nu, summ);
        }
    }

    /********************************************************************************************/
    /************************** Gradient ********************************************************/
    /********************************************************************************************/
    std::map<std::string, SharedMatrix> gradients;
    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Core");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("3index");
    gradient_terms.push_back("Metric");
    gradient_terms.push_back("Total");

    /********************************************************************************************/
    /************************** Nuclear Gradient ************************************************/
    /********************************************************************************************/
    // => Nuclear Gradient <= //
    gradients["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1(dipole_field_strength_).clone());
    gradients["Nuclear"]->set_name("Nuclear Gradient");
    gradients["Nuclear"]->print_atom_vector();

    /********************************************************************************************/
    /************************** One-Electron Gradient *******************************************/
    /********************************************************************************************/
    auto mints = std::make_shared<MintsHelper>(basisset_, options_);

    timer_on("Grad: T V Perturb");
    auto D = std::make_shared<Matrix>("AO-basis OPDM", nmo_, nmo_);
    Dso->to_shared_matrix(D);
    gradients["Core"] = mints->core_hamiltonian_grad(D);
    gradients["Core"]->print_atom_vector();
    timer_off("Grad: T V Perturb");

    /********************************************************************************************/
    /************************** Overlap Gradient ************************************************/
    /********************************************************************************************/
    timer_on("Grad: S");
    auto W = std::make_shared<Matrix>("AO-basis Energy-Weighted OPDM", nmo_, nmo_);
    Wso->to_shared_matrix(W);
    gradients["Overlap"] = mints->overlap_grad(W);
    gradients["Overlap"]->scale(-1.0);
    gradients["Overlap"]->print_atom_vector();
    timer_off("Grad: S");

    /********************************************************************************************/
    /************************** Two-electron Gradient *******************************************/
    /********************************************************************************************/
    // Two-electron gradients
    int df_ints_num_threads_ = 1;
#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif

    // Read in the basis set informations
    std::shared_ptr<BasisSet> primary_ = get_basisset("ORBITAL");
    std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_SCF");
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    // auxiliary_->print();
    int nbasis = primary_->nbf();

    /********************************************************************************************/
    /************************** Metric Gradient *************************************************/
    /********************************************************************************************/
    // JPQ_X
    timer_on("Grad: Metric");
    gradients["Metric"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["Metric"]->set_name("Metric Gradient");
    gradients["Metric"]->zero();

    // => Sizing <= //
    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory(
        new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //
    std::vector<SharedMatrix> Jtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jtemps.push_back(SharedMatrix(new Matrix("Jtemp", natom, 3)));
    }

    std::vector<std::pair<int, int> > PQ_pairs;
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

    int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        Jint[thread]->compute_shell_deriv1(P, 0, Q, 0);
        const double* buffer = Jint[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int cP = auxiliary_->shell(P).ncartesian();
        int aP = auxiliary_->shell(P).ncenter();
        int oP = auxiliary_->shell(P).function_index();

        int nQ = auxiliary_->shell(Q).nfunction();
        int cQ = auxiliary_->shell(Q).ncartesian();
        int aQ = auxiliary_->shell(Q).ncenter();
        int oQ = auxiliary_->shell(Q).function_index();

        int ncart = cP * cQ;
        const double* Px = buffer + 0 * ncart;
        const double* Py = buffer + 1 * ncart;
        const double* Pz = buffer + 2 * ncart;
        const double* Qx = buffer + 3 * ncart;
        const double* Qy = buffer + 4 * ncart;
        const double* Qz = buffer + 5 * ncart;

        double perm = (P == Q ? 1.0 : 2.0);

        double** grad_Jp;
        grad_Jp = Jtemps[thread]->pointer();

        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) {
                double Uval = perm * Gaux_ref->get(p + oP, q + oQ);
                grad_Jp[aP][0] -= Uval * (*Px);
                grad_Jp[aP][1] -= Uval * (*Py);
                grad_Jp[aP][2] -= Uval * (*Pz);
                grad_Jp[aQ][0] -= Uval * (*Qx);
                grad_Jp[aQ][1] -= Uval * (*Qy);
                grad_Jp[aQ][2] -= Uval * (*Qz);

                Px++;
                Py++;
                Pz++;
                Qx++;
                Qy++;
                Qz++;
            }
        }
    }

    // => Temporary Gradient Reduction <= //
    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients["Metric"]->add(Jtemps[t]);
    }

    gradients["Metric"]->print_atom_vector();
    timer_off("Grad: Metric");

    /********************************************************************************************/
    /************************** 3-Index Gradient ************************************************/
    /********************************************************************************************/
    // (Q | mu nu)^X
    timer_on("Grad: 3index");

    // int naux = nQ_ref;
    gradients["3index"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["3index"]->set_name("3index Gradient");
    gradients["3index"]->zero();

    // => Sizing <= //
    // int natom = primary_->molecule()->natom();
    // int nso = primary_->nbf();
    // int naux = auxiliary_->nbf();

    std::shared_ptr<ERISieve> sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, 0.0));
    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
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

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory2(
        new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri(1)));
    }

    // => Temporary Gradients <= //
    std::vector<SharedMatrix> Jtemps2;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jtemps2.push_back(SharedMatrix(new Matrix("Jtemp2", natom, 3)));
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

            const double* buffer = eri[thread]->buffer();

            int nP = auxiliary_->shell(P).nfunction();
            int cP = auxiliary_->shell(P).ncartesian();
            int aP = auxiliary_->shell(P).ncenter();
            int oP = auxiliary_->shell(P).function_index() - pstart;

            int nM = primary_->shell(M).nfunction();
            int cM = primary_->shell(M).ncartesian();
            int aM = primary_->shell(M).ncenter();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int cN = primary_->shell(N).ncartesian();
            int aN = primary_->shell(N).ncenter();
            int oN = primary_->shell(N).function_index();

            int ncart = cP * cM * cN;
            const double* Px = buffer + 0 * ncart;
            const double* Py = buffer + 1 * ncart;
            const double* Pz = buffer + 2 * ncart;
            const double* Mx = buffer + 3 * ncart;
            const double* My = buffer + 4 * ncart;
            const double* Mz = buffer + 5 * ncart;
            const double* Nx = buffer + 6 * ncart;
            const double* Ny = buffer + 7 * ncart;
            const double* Nz = buffer + 8 * ncart;

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Jp;
            grad_Jp = Jtemps2[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        double Ival = 1.0 * perm * gQso_ref->get(p + oP, (m + oM) * nso + (n + oN));
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
        gradients["3index"]->add(Jtemps2[t]);
    }

    gradients["3index"]->print_atom_vector();
    timer_off("Grad: 3index");

    /********************************************************************************************/
    /************************** Total Gradient **************************************************/
    /********************************************************************************************/
    // => Total Gradient <= //
    SharedMatrix total = SharedMatrix(gradients["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients.count(gradient_terms[i])) {
            total->add(gradients[gradient_terms[i]]);
        }
    }

    gradients["Total"] = total;
    gradients["Total"]->set_name("Total Gradient");

    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients.count(gradient_terms[i])) {
                gradients[gradient_terms[i]]->print_atom_vector();
            }
        }
    } else {
        gradients["Total"]->print_atom_vector();
    }

    // outfile->Printf("\tref_grad is done. \n");
}  // end

}  // namespace dfoccwave
}  // namespace psi
