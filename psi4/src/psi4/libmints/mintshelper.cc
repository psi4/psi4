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

#include "x2cint.h"

#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"

#include "psi4/libmints/matrix.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_dkh
#include <DKH/DKH_MANGLE.h>
#define F_DKH DKH_MANGLE_MODULE(dkh_main, dkh, DKH_MAIN, DKH)

extern "C" {
void F_DKH(double *S, double *V, double *T, double *pVp, int *nbf, int *dkh_order);
}
#endif

namespace psi {

/**
 * IWLWriter functor for use with SO TEIs
 **/
class PSI_API IWLWriter {
    IWL &writeto_;
    size_t count_;
    int &current_buffer_count_;

    Label *plabel_;
    Value *pvalue_;

   public:
    IWLWriter(IWL &writeto) : writeto_(writeto), count_(0), current_buffer_count_(writeto_.index()) {
        plabel_ = writeto_.labels();
        pvalue_ = writeto_.values();
    }

    void operator()(int i, int j, int k, int l, int, int, int, int, int, int, int, int, double value) {
        int current_label_position = 4 * current_buffer_count_;

        // Save the labels
        plabel_[current_label_position++] = i;
        plabel_[current_label_position++] = j;
        plabel_[current_label_position++] = k;
        plabel_[current_label_position] = l;

        // Save the value
        pvalue_[current_buffer_count_++] = value;

        // Increment overall counter
        count_++;

        // If our IWL buffer is full dump to disk.
        if (current_buffer_count_ == writeto_.ints_per_buffer()) {
            writeto_.last_buffer() = 0;
            writeto_.buffer_count() = current_buffer_count_;
            writeto_.put();
            current_buffer_count_ = 0;
        }
    }

    size_t count() const { return count_; }
};

MintsHelper::MintsHelper(std::shared_ptr<BasisSet> basis, Options &options, int print)
    : options_(options), print_(print) {
    init_helper(basis);
}

MintsHelper::MintsHelper(std::shared_ptr<Wavefunction> wavefunction)
    : options_(wavefunction->options()), print_(wavefunction->get_print()) {
    init_helper(wavefunction);
}

MintsHelper::~MintsHelper() {}

void MintsHelper::init_helper(std::shared_ptr<Wavefunction> wavefunction) {
    if (wavefunction->basisset().get() == 0) {
        outfile->Printf("  Wavefunction does not have a basisset!");
        throw PSIEXCEPTION("Wavefunction does not have a basisset, what did you do?!");
    }

    psio_ = wavefunction->psio();
    basisset_ = wavefunction->basisset();
    molecule_ = basisset_->molecule();

    // Make sure molecule is valid.
    molecule_->update_geometry();

    common_init();
}

void MintsHelper::init_helper(std::shared_ptr<BasisSet> basis) {
    basisset_ = basis;
    molecule_ = basis->molecule();
    psio_ = _default_psio_lib_;

    // Make sure molecule is valid.
    molecule_->update_geometry();

    common_init();
}

void MintsHelper::common_init() {
    // Print the molecule.
    if (print_) molecule_->print();

    // Print the basis set
    if (print_) basisset_->print_detail();

    // How many threads?
    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    // Create integral factory
    integral_ = std::make_shared<IntegralFactory>(basisset_);

    // Get the SO basis object.
    sobasis_ = std::make_shared<SOBasisSet>(basisset_, integral_);

    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis_->dimension();

    // Create a matrix factory and initialize it
    factory_ = std::make_shared<MatrixFactory>();
    factory_->init_with(dimension, dimension);

    // Integral cutoff
    cutoff_ = Process::environment.options.get_double("INTS_TOLERANCE");
}

std::shared_ptr<PetiteList> MintsHelper::petite_list() const {
    auto pt = std::make_shared<PetiteList>(basisset_, integral_);
    return pt;
}

std::shared_ptr<PetiteList> MintsHelper::petite_list(bool val) const {
    auto pt = std::make_shared<PetiteList>(basisset_, integral_, val);
    return pt;
}

std::shared_ptr<BasisSet> MintsHelper::basisset() const { return basisset_; }

std::shared_ptr<SOBasisSet> MintsHelper::sobasisset() const { return sobasis_; }

std::shared_ptr<MatrixFactory> MintsHelper::factory() const { return factory_; }

std::shared_ptr<IntegralFactory> MintsHelper::integral() const { return integral_; }

int MintsHelper::nbf() const { return basisset_->nbf(); }

void MintsHelper::integrals() {
    if (print_) {
        outfile->Printf(" MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");
    }

    // Get ERI object
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i) tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->eri()));
    auto eri = std::make_shared<TwoBodySOInt>(tb, integral_);

    //// Print out some useful information
    if (print_) {
        outfile->Printf("   Calculation information:\n");
        outfile->Printf("      Number of threads:              %4d\n", nthread_);
        outfile->Printf("      Number of atoms:                %4d\n", molecule_->natom());
        outfile->Printf("      Number of AO shells:            %4d\n", basisset_->nshell());
        outfile->Printf("      Number of SO shells:            %4d\n", sobasis_->nshell());
        outfile->Printf("      Number of primitives:           %4d\n", basisset_->nprimitive());
        outfile->Printf("      Number of atomic orbitals:      %4d\n", basisset_->nao());
        outfile->Printf("      Number of basis functions:      %4d\n\n", basisset_->nbf());
        outfile->Printf("      Number of irreps:               %4d\n", sobasis_->nirrep());
        outfile->Printf("      Integral cutoff                 %4.2e\n", cutoff_);
        outfile->Printf("      Number of functions per irrep: [");
        for (int i = 0; i < sobasis_->nirrep(); ++i) {
            outfile->Printf("%4d ", sobasis_->nfunction_in_irrep(i));
        }
        outfile->Printf("]\n\n");
    }

    // Compute one-electron integrals.
    one_electron_integrals();

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio_.get(), PSIF_SO_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Let the user know what we're doing.
    if (print_) {
        outfile->Printf("      Computing two-electron integrals...");
    }

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        eri->compute_shell(shellIter, writer);
    }

    // Flush out buffers.
    ERIOUT.flush(1);

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    if (print_) {
        outfile->Printf("done\n");
        outfile->Printf(
            "      Computed %lu non-zero two-electron integrals.\n"
            "        Stored in file %d.\n\n",
            writer.count(), PSIF_SO_TEI);
    }
}

void MintsHelper::integrals_erf(double w) {
    double omega = (w == -1.0 ? options_.get_double("OMEGA_ERF") : w);

    IWL ERIOUT(psio_.get(), PSIF_SO_ERF_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Get ERI object
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i) tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->erf_eri(omega)));
    auto erf = std::make_shared<TwoBodySOInt>(tb, integral_);

    // Let the user know what we're doing.
    outfile->Printf("      Computing non-zero ERF integrals (omega = %.3f)...", omega);

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) erf->compute_shell(shellIter, writer);

    // Flush the buffers
    ERIOUT.flush(1);

    // Keep the integrals around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    outfile->Printf("done\n");
    outfile->Printf(
        "      Computed %lu non-zero ERF integrals.\n"
        "        Stored in file %d.\n\n",
        writer.count(), PSIF_SO_ERF_TEI);
}

void MintsHelper::integrals_erfc(double w) {
    double omega = (w == -1.0 ? options_.get_double("OMEGA_ERF") : w);

    IWL ERIOUT(psio_.get(), PSIF_SO_ERFC_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Get ERI object
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i)
        tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->erf_complement_eri(omega)));
    auto erf = std::make_shared<TwoBodySOInt>(tb, integral_);

    // Let the user know what we're doing.
    outfile->Printf("      Computing non-zero ERFComplement integrals...");

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) erf->compute_shell(shellIter, writer);

    // Flush the buffers
    ERIOUT.flush(1);

    // Keep the integrals around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    outfile->Printf("done\n");
    outfile->Printf(
        "      Computed %lu non-zero ERFComplement integrals.\n"
        "        Stored in file %d.\n\n",
        writer.count(), PSIF_SO_ERFC_TEI);
}

void MintsHelper::one_electron_integrals() {
    //    outfile->Printf( " OEINTS: Wrapper to libmints.\n   by Justin Turney\n\n");
    //
    //    // Print out some useful information
    //    outfile->Printf( "   Calculation information:\n");
    //    outfile->Printf( "      Number of atoms:                %4d\n", molecule_->natom());
    //    outfile->Printf( "      Number of AO shells:            %4d\n", basisset_->nshell());
    //    outfile->Printf( "      Number of SO shells:            %4d\n", sobasis_->nshell());
    //    outfile->Printf( "      Number of primitives:           %4d\n", basisset_->nprimitive());
    //    outfile->Printf( "      Number of atomic orbitals:      %4d\n", basisset_->nao());
    //    outfile->Printf( "      Number of basis functions:      %4d\n\n", basisset_->nbf());
    //    outfile->Printf( "      Number of irreps:               %4d\n", sobasis_->nirrep());
    //    outfile->Printf( "      Number of functions per irrep: [");
    //    for (int i=0; i<sobasis_->nirrep(); ++i) {
    //        outfile->Printf( "%4d ", sobasis_->nfunction_in_irrep(i));
    //    }
    //    outfile->Printf( "]\n\n");

    // Compute and dump one-electron SO integrals.

    if (options_.get_str("RELATIVISTIC") == "NO" || options_.get_str("RELATIVISTIC") == "DKH") {
        // Overlap
        so_overlap()->save(psio_, PSIF_OEI);

        // Kinetic
        so_kinetic()->save(psio_, PSIF_OEI);

        // Potential -- DKH perturbation added to potential integrals if needed.
        so_potential()->save(psio_, PSIF_OEI);
    } else if (options_.get_str("RELATIVISTIC") == "X2C") {
        outfile->Printf(" OEINTS: Using relativistic (X2C) overlap, kinetic, and potential integrals.\n");

        if (!rel_basisset_) {
            throw PSIEXCEPTION("OEINTS: X2C requested, but relativistic basis was not set.");
        }
        X2CInt x2cint;
        SharedMatrix so_overlap_x2c = so_overlap();
        SharedMatrix so_kinetic_x2c = so_kinetic();
        SharedMatrix so_potential_x2c = so_potential();
        x2cint.compute(basisset_, rel_basisset_, so_overlap_x2c, so_kinetic_x2c, so_potential_x2c);

        // Overlap
        so_overlap_x2c->save(psio_, PSIF_OEI);

        // Kinetic
        so_kinetic_x2c->save(psio_, PSIF_OEI);

        // Potential
        so_potential_x2c->save(psio_, PSIF_OEI);
    }

    // Dipoles
    std::vector<SharedMatrix> dipole_mats = so_dipole();
    for (SharedMatrix m : dipole_mats) {
        m->save(psio_, PSIF_OEI);
    }

    // Quadrupoles
    std::vector<SharedMatrix> quadrupole_mats = so_quadrupole();
    for (SharedMatrix m : quadrupole_mats) {
        m->save(psio_, PSIF_OEI);
    }

    if (print_) {
        outfile->Printf(
            " OEINTS: Overlap, kinetic, potential, dipole, and quadrupole integrals\n"
            "         stored in file %d.\n\n",
            PSIF_OEI);
    }
}

void MintsHelper::integral_gradients() {
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_derivatives", __FILE__, __LINE__);
}

void MintsHelper::integral_hessians() {
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_hessians", __FILE__, __LINE__);
}

void MintsHelper::one_body_ao_computer(std::vector<std::shared_ptr<OneBodyAOInt>> ints, SharedMatrix out, bool symm) {
    // Grab basis info
    std::shared_ptr<BasisSet> bs1 = ints[0]->basis1();
    std::shared_ptr<BasisSet> bs2 = ints[0]->basis2();

    // Limit to the number of incoming onbody ints
    size_t nthread = nthread_;
    if (nthread > ints.size()) {
        nthread = ints.size();
    }

    // Grab the buffers
    std::vector<const double *> ints_buff(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        ints_buff[thread] = ints[thread]->buffer();
    }

    double **outp = out->pointer();

// Loop it
#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t MU = 0; MU < bs1->nshell(); ++MU) {
        const size_t num_mu = bs1->shell(MU).nfunction();
        const size_t index_mu = bs1->shell(MU).function_index();

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        if (symm) {
            // Triangular
            for (size_t NU = 0; NU <= MU; ++NU) {
                const size_t num_nu = bs2->shell(NU).nfunction();
                const size_t index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        outp[nu][mu] = outp[mu][nu] = ints_buff[rank][index++];
                    }
                }
            }  // End NU
        }      // End Symm
        else {
            // Rectangular
            for (size_t NU = 0; NU < bs2->nshell(); ++NU) {
                const size_t num_nu = bs2->shell(NU).nfunction();
                const size_t index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        // printf("%zu %zu | %zu %zu | %lf\n", MU, NU, mu, nu, ints_buff[rank][index]);
                        outp[mu][nu] = ints_buff[rank][index++];
                    }
                }
            }  // End NU
        }      // End Rectangular
    }          // End Mu
}
void MintsHelper::grad_two_center_computer(std::vector<std::shared_ptr<OneBodyAOInt>> ints, SharedMatrix D,
                                           SharedMatrix out) {
    // Grab basis info
    std::shared_ptr<BasisSet> bs1 = ints[0]->basis1();
    std::shared_ptr<BasisSet> bs2 = ints[0]->basis2();
    if (bs1 != bs2) {
        throw PSIEXCEPTION("BasisSets must be the same for deriv1");
    }

    if (D->nirrep() > 1) {
        throw PSIEXCEPTION("Density must be of C1 symmetry");
    }

    // Limit to the number of incoming onbody ints
    size_t nthread = nthread_;
    if (nthread > ints.size()) {
        nthread = ints.size();
    }

    // Grab the buffers
    std::vector<const double *> ints_buff(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        ints_buff[thread] = ints[thread]->buffer();
    }

    double **outp = out->pointer();
    double **Dp = D->pointer();

#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t P = 0; P < basisset_->nshell(); P++) {
        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        for (size_t Q = 0; Q <= P; Q++) {
            ints[rank]->compute_shell_deriv1(P, Q);

            size_t nP = basisset_->shell(P).nfunction();
            size_t oP = basisset_->shell(P).function_index();
            size_t aP = basisset_->shell(P).ncenter();

            size_t nQ = basisset_->shell(Q).nfunction();
            size_t oQ = basisset_->shell(Q).function_index();
            size_t aQ = basisset_->shell(Q).ncenter();

            size_t offset = nP * nQ;
            const double *ref = ints_buff[rank];
            double perm = (P == Q ? 1.0 : 2.0);

            // Px
            double Px = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Px += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aP][0] += Px;

            // Py
            double Py = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Py += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aP][1] += Py;

            // Pz
            double Pz = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Pz += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aP][2] += Pz;

            // Qx
            double Qx = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Qx += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aQ][0] += Qx;

            // Qy
            double Qy = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Qy += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aQ][1] += Qy;

            // Qz
            double Qz = 0.0;
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    Qz += perm * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
#pragma omp atomic update
            outp[aQ][2] += Qz;
        }
    }
}

SharedMatrix MintsHelper::ao_overlap() {
    // Overlap
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_overlap()));
    }
    auto overlap_mat = std::make_shared<Matrix>(PSIF_AO_S, basisset_->nbf(), basisset_->nbf());
    one_body_ao_computer(ints_vec, overlap_mat, true);

    // Is this needed?
    //    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

// JWM 4/3/2015
SharedMatrix MintsHelper::ao_overlap(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2) {
    // Overlap
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_overlap()));
    }
    auto overlap_mat = std::make_shared<Matrix>(PSIF_AO_S, bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, overlap_mat, false);
    return overlap_mat;
}

SharedMatrix MintsHelper::ao_kinetic() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_kinetic()));
    }
    auto kinetic_mat = std::make_shared<Matrix>("AO-basis Kinetic Ints", basisset_->nbf(), basisset_->nbf());
    one_body_ao_computer(ints_vec, kinetic_mat, true);
    return kinetic_mat;
}

SharedMatrix MintsHelper::ao_kinetic(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2) {
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_kinetic()));
    }
    auto kinetic_mat = std::make_shared<Matrix>("AO-basis Kinetic Ints", bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, kinetic_mat, false);
    return kinetic_mat;
}

SharedMatrix MintsHelper::ao_potential() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential()));
    }
    SharedMatrix potential_mat =
        std::make_shared<Matrix>("AO-basis Potential Ints", basisset_->nbf(), basisset_->nbf());
    one_body_ao_computer(ints_vec, potential_mat, true);
    return potential_mat;
}

SharedMatrix MintsHelper::ao_potential(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2) {
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_potential()));
    }
    auto potential_mat = std::make_shared<Matrix>("AO-basis Potential Ints", bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, potential_mat, false);
    return potential_mat;
}

SharedMatrix MintsHelper::ao_ecp() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_ecp()));
    }
    auto ecp_mat = std::make_shared<Matrix>("AO-basis ECP Ints", basisset_->nbf(), basisset_->nbf());
    one_body_ao_computer(ints_vec, ecp_mat, true);
    return ecp_mat;
}

SharedMatrix MintsHelper::ao_ecp(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2) {
    IntegralFactory factory(bs1, bs2, bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_ecp()));
    }
    auto ecp_mat = std::make_shared<Matrix>("AO-basis ECP Ints", bs1->nbf(), bs2->nbf());
    one_body_ao_computer(ints_vec, ecp_mat, false);
    return ecp_mat;
}

SharedMatrix MintsHelper::ao_pvp() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_rel_potential()));
    }
    auto pVp_mat = std::make_shared<Matrix>("AO-basis pVp Ints", basisset_->nbf(), basisset_->nbf());
    one_body_ao_computer(ints_vec, pVp_mat, true);
    return pVp_mat;
}

SharedMatrix MintsHelper::ao_dkh(int dkh_order) {
#ifdef USING_dkh
    MintsHelper decon(rel_basisset_);
    SharedMatrix S = decon.ao_overlap();
    SharedMatrix T = decon.ao_kinetic();
    SharedMatrix Torig = T->clone();
    SharedMatrix V = decon.ao_potential();
    SharedMatrix Vorig = V->clone();
    SharedMatrix pVp = decon.ao_pvp();
    SharedMatrix H_dk = T->clone();
    H_dk->zero();

    double *Sp = S->pointer()[0];
    double *Tp = T->pointer()[0];
    double *Vp = V->pointer()[0];
    double *pVpp = pVp->pointer()[0];

    if (dkh_order < 1) dkh_order = 2;
    if (dkh_order > 4) dkh_order = 4;

    outfile->Printf("    Computing %d-order Douglas-Kroll-Hess integrals.\n", dkh_order);

    int nbf = rel_basisset_->nbf();
    //    rel_basisset_->print_detail();

    // Call DKH code from Markus Reiher
    F_DKH(Sp, Vp, Tp, pVpp, &nbf, &dkh_order);

    H_dk->add(V);
    H_dk->add(T);
    H_dk->subtract(Vorig);
    H_dk->subtract(Torig);

    H_dk->set_name("AO-basis DKH Ints");
    //    H_dk->print();

    // DKH is solved in the decontracted basis ... transform to the contracted.

    SharedMatrix S_inv = S->clone();
    S_inv->general_invert();

    SharedMatrix S_cd = ao_overlap(basisset_, rel_basisset_);
    //    S_cd->print();

    auto D = std::make_shared<Matrix>("D", nbf, basisset_->nbf());
    //  Form D = S_uu^{-1} S_uc.  Notice that we transpose S_cu
    D->gemm(false, true, 1.0, S_inv, S_cd, 0.0);

    auto H_dk_cc = std::make_shared<Matrix>("AO-basis DKH Ints", basisset_->nbf(), basisset_->nbf());
    H_dk_cc->transform(H_dk, D);
    //    H_dk_cc->print();

    return H_dk_cc;
#else
    outfile->Printf("    Douglas-Kroll-Hess integrals of order %d requested but are not available.\n", dkh_order);
    throw PSIEXCEPTION("Douglas-Kroll-Hess integrals requested but were not compiled in.");
#endif
}

SharedMatrix MintsHelper::so_dkh(int dkh_order) {
    SharedMatrix dkh = factory_->create_shared_matrix("SO Douglas-Kroll-Hess Integrals");
    dkh->apply_symmetry(ao_dkh(dkh_order), petite_list()->aotoso());
    return dkh;
}

SharedMatrix MintsHelper::ao_helper(const std::string &label, std::shared_ptr<TwoBodyAOInt> ints) {
    std::shared_ptr<BasisSet> bs1 = ints->basis1();
    std::shared_ptr<BasisSet> bs2 = ints->basis2();
    std::shared_ptr<BasisSet> bs3 = ints->basis3();
    std::shared_ptr<BasisSet> bs4 = ints->basis4();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();
    int nbf4 = bs4->nbf();

    auto I = std::make_shared<Matrix>(label, nbf1 * nbf2, nbf3 * nbf4);
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();

    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {
                for (int Q = 0; Q < bs4->nshell(); Q++) {
                    ints->compute_shell(M, N, P, Q);

                    for (int m = 0, index = 0; m < bs1->shell(M).nfunction(); m++) {
                        for (int n = 0; n < bs2->shell(N).nfunction(); n++) {
                            for (int p = 0; p < bs3->shell(P).nfunction(); p++) {
                                for (int q = 0; q < bs4->shell(Q).nfunction(); q++, index++) {
                                    Ip[(bs1->shell(M).function_index() + m) * nbf2 + bs2->shell(N).function_index() + n]
                                      [(bs3->shell(P).function_index() + p) * nbf4 + bs4->shell(Q).function_index() +
                                       q] = buffer[index];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3, nbf4};
    I->set_numpy_shape(nshape);

    return I;
}

SharedMatrix MintsHelper::ao_shell_getter(const std::string &label, std::shared_ptr<TwoBodyAOInt> ints, int M, int N,
                                          int P, int Q) {
    int mfxn = basisset_->shell(M).nfunction();
    int nfxn = basisset_->shell(N).nfunction();
    int pfxn = basisset_->shell(P).nfunction();
    int qfxn = basisset_->shell(Q).nfunction();
    auto I = std::make_shared<Matrix>(label, mfxn * nfxn, pfxn * qfxn);
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();

    ints->compute_shell(M, N, P, Q);

    for (int m = 0, index = 0; m < mfxn; m++) {
        for (int n = 0; n < nfxn; n++) {
            for (int p = 0; p < pfxn; p++) {
                for (int q = 0; q < qfxn; q++, index++) {
                    Ip[m * nfxn + n][p * qfxn + q] = buffer[index];
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{mfxn, nfxn, pfxn, qfxn};
    I->set_numpy_shape(nshape);

    return I;
}

SharedMatrix MintsHelper::ao_erf_eri(double omega, std::shared_ptr<IntegralFactory> input_factory) {
    std::shared_ptr<IntegralFactory> factory;
    if (input_factory) {
        factory = input_factory;
    } else {
        factory = integral_;
    }
    return ao_helper("AO ERF ERI Integrals", std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega)));
}

SharedMatrix MintsHelper::ao_eri(std::shared_ptr<IntegralFactory> input_factory) {
    std::shared_ptr<IntegralFactory> factory;
    if (input_factory) {
        factory = input_factory;
    } else {
        factory = integral_;
    }

    return ao_helper("AO ERI Tensor", std::shared_ptr<TwoBodyAOInt>(factory->eri()));
}

SharedMatrix MintsHelper::ao_eri(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.eri());
    return ao_helper("AO ERI Tensor", ints);
}

SharedMatrix MintsHelper::ao_eri_shell(int M, int N, int P, int Q) {
    if (eriInts_ == 0) {
        eriInts_ = std::shared_ptr<TwoBodyAOInt>(integral_->eri());
    }
    return ao_shell_getter("AO ERI Tensor", eriInts_, M, N, P, Q);
}

SharedMatrix MintsHelper::ao_erfc_eri(double omega) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->erf_complement_eri(omega));
    return ao_helper("AO ERFC ERI Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12(std::shared_ptr<CorrelationFactor> corr) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12(corr));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12(std::shared_ptr<CorrelationFactor> corr, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12(corr));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_scaled(std::shared_ptr<CorrelationFactor> corr) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12_scaled(corr));
    return ao_helper("AO F12 Scaled Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_scaled(std::shared_ptr<CorrelationFactor> corr, std::shared_ptr<BasisSet> bs1,
                                        std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                                        std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12_scaled(corr));
    return ao_helper("AO F12 Scaled Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::shared_ptr<CorrelationFactor> corr) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12_squared(corr));
    return ao_helper("AO F12 Squared Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::shared_ptr<CorrelationFactor> corr, std::shared_ptr<BasisSet> bs1,
                                         std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                                         std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12_squared(corr));
    return ao_helper("AO F12 Squared Tensor", ints);
}

SharedMatrix MintsHelper::ao_3coverlap_helper(const std::string &label, std::shared_ptr<ThreeCenterOverlapInt> ints) {
    std::shared_ptr<BasisSet> bs1 = ints->basis1();
    std::shared_ptr<BasisSet> bs2 = ints->basis2();
    std::shared_ptr<BasisSet> bs3 = ints->basis3();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();

    auto I = std::make_shared<Matrix>(label, nbf1 * nbf2, nbf3);
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();

    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {
                ints->compute_shell(M, N, P);
                int Mfi = bs1->shell(M).function_index();
                int Nfi = bs2->shell(N).function_index();
                int Pfi = bs3->shell(P).function_index();

                for (int m = Mfi, index = 0; m < (Mfi + bs1->shell(M).nfunction()); m++) {
                    for (int n = Nfi; n < (Nfi + bs2->shell(N).nfunction()); n++) {
                        for (int p = Pfi; p < (Pfi + bs3->shell(P).nfunction()); p++) {
                            Ip[m * nbf2 + n][p] = buffer[index++];
                        }
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3};
    I->set_numpy_shape(nshape);

    return I;
}
SharedMatrix MintsHelper::ao_3coverlap() {
    std::vector<SphericalTransform> trans;
    for (int i = 0; i <= basisset_->max_am(); i++) {
        trans.push_back(SphericalTransform(i));
    }
    std::shared_ptr<ThreeCenterOverlapInt> ints =
        std::make_shared<ThreeCenterOverlapInt>(trans, basisset_, basisset_, basisset_);
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}

SharedMatrix MintsHelper::ao_3coverlap(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                       std::shared_ptr<BasisSet> bs3) {
    int max_am = std::max(std::max(bs1->max_am(), bs2->max_am()), bs3->max_am());
    std::vector<SphericalTransform> trans;
    for (int i = 0; i <= max_am; i++) {
        trans.push_back(SphericalTransform(i));
    }
    auto ints = std::make_shared<ThreeCenterOverlapInt>(trans, bs1, bs2, bs3);
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12g12(std::shared_ptr<CorrelationFactor> corr) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12g12(corr));
    return ao_helper("AO F12G12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_double_commutator(std::shared_ptr<CorrelationFactor> corr) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12_double_commutator(corr));
    return ao_helper("AO F12 Double Commutator Tensor", ints);
}

SharedMatrix MintsHelper::mo_erf_eri(double omega, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_erf_eri(omega), C1, C2, C3, C4);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_erfc_eri(double omega, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                      SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_erfc_eri(omega), C1, C2, C3, C4);
    mo_ints->set_name("MO ERFC ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2,
                                 SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_squared(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2,
                                         SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_squared(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Squared Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12g12(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2,
                                    SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12g12(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12G12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_double_commutator(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1,
                                                   SharedMatrix C2, SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_double_commutator(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Double Commutator Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri(SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), C1, C2, C3, C4);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_erf_eri(double omega, SharedMatrix Co, SharedMatrix Cv) {
    SharedMatrix mo_ints = mo_eri_helper(ao_erf_eri(omega), Co, Cv);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri(SharedMatrix Co, SharedMatrix Cv) {
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), Co, Cv);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                        SharedMatrix C4) {
    int nso = basisset_->nbf();
    int n1 = C1->colspi()[0];
    int n2 = C2->colspi()[0];
    int n3 = C3->colspi()[0];
    int n4 = C4->colspi()[0];

    double **C1p = C1->pointer();
    double **C2p = C2->pointer();
    double **C3p = C3->pointer();
    double **C4p = C4->pointer();

    double **Isop = Iso->pointer();
    auto I2 = std::make_shared<Matrix>("MO ERI Tensor", n1 * nso, nso * nso);
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', n1, nso * (size_t)nso * nso, nso, 1.0, C1p[0], n1, Isop[0], nso * (size_t)nso * nso, 0.0, I2p[0],
            nso * (size_t)nso * nso);

    Iso.reset();
    auto I3 = std::make_shared<Matrix>("MO ERI Tensor", n1 * nso, nso * n3);
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', n1 * (size_t)nso * nso, n3, nso, 1.0, I2p[0], nso, C3p[0], n3, 0.0, I3p[0], n3);

    I2.reset();
    auto I4 = std::make_shared<Matrix>("MO ERI Tensor", nso * n1, n3 * nso);
    double **I4p = I4->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * n1 + i][j * nso + n] = I3p[i * nso + m][n * n3 + j];
                }
            }
        }
    }

    I3.reset();
    auto I5 = std::make_shared<Matrix>("MO ERI Tensor", n2 * n1, n3 * nso);
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', n2, n1 * (size_t)n3 * nso, nso, 1.0, C2p[0], n2, I4p[0], n1 * (size_t)n3 * nso, 0.0, I5p[0],
            n1 * (size_t)n3 * nso);

    I4.reset();
    auto I6 = std::make_shared<Matrix>("MO ERI Tensor", n2 * n1, n3 * n4);
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', n2 * (size_t)n1 * n3, n4, nso, 1.0, I5p[0], nso, C4p[0], n4, 0.0, I6p[0], n4);

    I5.reset();
    auto Imo = std::make_shared<Matrix>("MO ERI Tensor", n1 * n2, n3 * n4);
    double **Imop = Imo->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int a = 0; a < n2; a++) {
                for (int b = 0; b < n4; b++) {
                    Imop[i * n2 + a][j * n4 + b] = I6p[a * n1 + i][j * n4 + b];
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{n1, n2, n3, n4};
    Imo->set_numpy_shape(nshape);

    return Imo;
}

SharedMatrix MintsHelper::mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv) {
    int nso = basisset_->nbf();
    int nocc = Co->colspi()[0];
    int nvir = Cv->colspi()[0];

    double **Cop = Co->pointer();
    double **Cvp = Cv->pointer();

    double **Isop = Iso->pointer();
    auto I2 = std::make_shared<Matrix>("MO ERI Tensor", nocc * nso, nso * nso);
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', nocc, nso * (size_t)nso * nso, nso, 1.0, Cop[0], nocc, Isop[0], nso * (size_t)nso * nso, 0.0,
            I2p[0], nso * (size_t)nso * nso);

    Iso.reset();
    auto I3 = std::make_shared<Matrix>("MO ERI Tensor", nocc * nso, nso * nocc);
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', nocc * (size_t)nso * nso, nocc, nso, 1.0, I2p[0], nso, Cop[0], nocc, 0.0, I3p[0], nocc);

    I2.reset();
    auto I4 = std::make_shared<Matrix>("MO ERI Tensor", nso * nocc, nocc * nso);
    double **I4p = I4->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * nocc + i][j * nso + n] = I3p[i * nso + m][n * nocc + j];
                }
            }
        }
    }

    I3.reset();
    auto I5 = std::make_shared<Matrix>("MO ERI Tensor", nvir * nocc, nocc * nso);
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', nvir, nocc * (size_t)nocc * nso, nso, 1.0, Cvp[0], nvir, I4p[0], nocc * (size_t)nocc * nso, 0.0,
            I5p[0], nocc * (size_t)nocc * nso);

    I4.reset();
    auto I6 = std::make_shared<Matrix>("MO ERI Tensor", nvir * nocc, nocc * nvir);
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', nvir * (size_t)nocc * nocc, nvir, nso, 1.0, I5p[0], nso, Cvp[0], nvir, 0.0, I6p[0], nvir);

    I5.reset();
    auto Imo = std::make_shared<Matrix>("MO ERI Tensor", nocc * nvir, nocc * nvir);
    double **Imop = Imo->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = 0; a < nvir; a++) {
                for (int b = 0; b < nvir; b++) {
                    Imop[i * nvir + a][j * nvir + b] = I6p[a * nocc + i][j * nvir + b];
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{nocc, nvir, nocc, nvir};
    Imo->set_numpy_shape(nshape);

    return Imo;
}

SharedMatrix MintsHelper::mo_spin_eri(SharedMatrix Co, SharedMatrix Cv) {
    int n1 = Co->colspi()[0];
    int n2 = Cv->colspi()[0];
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), Co, Cv);
    SharedMatrix mo_spin_ints = mo_spin_eri_helper(mo_ints, n1, n2);
    mo_ints.reset();
    mo_spin_ints->set_name("MO Spin ERI Tensor");
    return mo_spin_ints;
}

SharedMatrix MintsHelper::mo_spin_eri_helper(SharedMatrix Iso, int n1, int n2) {
    int n12 = n1 * 2;
    int n22 = n2 * 2;

    double **Isop = Iso->pointer();
    auto Ispin = std::make_shared<Matrix>("MO ERI Tensor", 4 * n1 * n1, 4 * n2 * n2);
    double **Ispinp = Ispin->pointer();

    double first, second;
    int mask1, mask2;
    for (int i = 0; i < n12; i++) {
        for (int j = 0; j < n12; j++) {
            for (int k = 0; k < n22; k++) {
                for (int l = 0; l < n22; l++) {
                    mask1 = (i % 2 == k % 2) * (j % 2 == l % 2);
                    mask2 = (i % 2 == l % 2) * (j % 2 == k % 2);

                    first = Isop[i / 2 * n2 + k / 2][j / 2 * n2 + l / 2];
                    second = Isop[i / 2 * n2 + l / 2][j / 2 * n2 + k / 2];
                    Ispinp[i * n12 + j][k * n22 + l] = first * mask1 - second * mask2;
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{n12, n12, n22, n22};
    Ispin->set_numpy_shape(nshape);

    return Ispin;
}

SharedMatrix MintsHelper::so_overlap() {
    if (factory_->nirrep() == 1) {
        SharedMatrix ret = ao_overlap();
        ret->set_name(PSIF_SO_S);
        return ret;
    } else {
        SharedMatrix overlap_mat(factory_->create_matrix(PSIF_SO_S));
        overlap_mat->apply_symmetry(ao_overlap(), petite_list()->aotoso());
        return overlap_mat;
    }
}

SharedMatrix MintsHelper::so_kinetic() {
    if (factory_->nirrep() == 1) {
        SharedMatrix ret = ao_kinetic();
        ret->set_name(PSIF_SO_T);
        return ret;
    } else {
        SharedMatrix kinetic_mat(factory_->create_matrix(PSIF_SO_T));
        kinetic_mat->apply_symmetry(ao_kinetic(), petite_list()->aotoso());
        return kinetic_mat;
    }
}

SharedMatrix MintsHelper::so_ecp() {
    if (!basisset_->has_ECP()) {
        SharedMatrix ecp_mat = factory_->create_shared_matrix("SO Basis ECP");
        ecp_mat->zero();
        outfile->Printf("\n\tWarning! ECP integrals requested, but no ECP basis detected.  Returning zeros.\n");
        return ecp_mat;
    }

    if (factory_->nirrep() == 1) {
        SharedMatrix ecp_mat(ao_ecp());
        ecp_mat->set_name("AO Basis ECP");
        return ecp_mat;
    } else {
        SharedMatrix ecp_mat = factory_->create_shared_matrix("SO Basis ECP");
        ecp_mat->apply_symmetry(ao_ecp(), petite_list()->aotoso());
        return ecp_mat;
    }
}

SharedMatrix MintsHelper::so_potential(bool include_perturbations) {
    // No symmetry
    SharedMatrix potential_mat;
    if (factory_->nirrep() == 1) {
        potential_mat = ao_potential();
        potential_mat->set_name(PSIF_SO_V);
    } else {
        potential_mat = factory_->create_shared_matrix(PSIF_SO_V);
        potential_mat->apply_symmetry(ao_potential(), petite_list()->aotoso());
    }

    // Add ECPs, if needed
    if (basisset_->has_ECP()) {
        potential_mat->add(so_ecp());
    }

    // Handle addition of any perturbations here and not in SCF code.
    if (include_perturbations) {
        if (options_.get_bool("PERTURB_H")) {
            std::string perturb_with = options_.get_str("PERTURB_WITH");
            Vector3 lambda(0.0, 0.0, 0.0);

            if (perturb_with == "DIPOLE_X")
                lambda[0] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE_Y")
                lambda[1] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE_Z")
                lambda[2] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE") {
                if (options_["PERTURB_DIPOLE"].size() != 3)
                    throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
                for (int n = 0; n < 3; ++n) lambda[n] = options_["PERTURB_DIPOLE"][n].to_double();
            } else {
                outfile->Printf("  MintsHelper doesn't understand the requested perturbation, might be done in SCF.");
            }

            OperatorSymmetry msymm(1, molecule_, integral_, factory_);
            std::vector<SharedMatrix> dipoles = msymm.create_matrices("Dipole");
            OneBodySOInt *so_dipole = integral_->so_dipole();
            so_dipole->compute(dipoles);

            if (lambda[0] != 0.0) {
                if (msymm.component_symmetry(0) != 0) {
                    outfile->Printf("  WARNING: Requested mu(x) perturbation, but mu(x) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(x).\n", lambda[0]);
                    dipoles[0]->scale(lambda[0]);
                    potential_mat->add(dipoles[0]);
                }
            }
            if (lambda[1] != 0.0) {
                if (msymm.component_symmetry(1) != 0) {
                    outfile->Printf("  WARNING: Requested mu(y) perturbation, but mu(y) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(y).\n", lambda[1]);
                    dipoles[1]->scale(lambda[1]);
                    potential_mat->add(dipoles[1]);
                }
            }
            if (lambda[2] != 0.0) {
                if (msymm.component_symmetry(2) != 0) {
                    outfile->Printf("  WARNING: Requested mu(z) perturbation, but mu(z) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(z).\n", lambda[2]);
                    dipoles[2]->scale(lambda[2]);
                    potential_mat->add(dipoles[2]);
                }
            }
        }

        if (options_.get_str("RELATIVISTIC") == "DKH") {
            int dkh_order = options_.get_int("DKH_ORDER");
            SharedMatrix dkh = so_dkh(dkh_order);

            outfile->Printf("    Adding Douglas-Kroll-Hess corrections to the potential integrals.\n");

            potential_mat->add(dkh);
        }
    }

    return potential_mat;
}

std::vector<SharedMatrix> MintsHelper::so_dipole() {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    return dipole;
}

std::vector<SharedMatrix> MintsHelper::so_quadrupole() {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Quadrupole");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector<SharedMatrix> MintsHelper::so_traceless_quadrupole() {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Traceless Quadrupole");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_traceless_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector<SharedMatrix> MintsHelper::so_nabla() {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::P, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla = msymm.create_matrices("SO Nabla");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_nabla());
    ints->compute(nabla);

    return nabla;
}

std::vector<SharedMatrix> MintsHelper::so_angular_momentum() {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::L, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> am = msymm.create_matrices("SO Angular Momentum");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_angular_momentum());
    ints->compute(am);

    return am;
}

std::vector<SharedMatrix> MintsHelper::ao_angular_momentum() {
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> angmom;

    angmom.push_back(std::make_shared<Matrix>("AO Lx", basisset_->nbf(), basisset_->nbf()));
    angmom.push_back(std::make_shared<Matrix>("AO Ly", basisset_->nbf(), basisset_->nbf()));
    angmom.push_back(std::make_shared<Matrix>("AO Lz", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_angular_momentum());
    ints->compute(angmom);

    return angmom;
}

std::vector<SharedMatrix> MintsHelper::ao_dipole() {
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole;

    dipole.push_back(std::make_shared<Matrix>("AO Mux", basisset_->nbf(), basisset_->nbf()));
    dipole.push_back(std::make_shared<Matrix>("AO Muy", basisset_->nbf(), basisset_->nbf()));
    dipole.push_back(std::make_shared<Matrix>("AO Muz", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_dipole());
    ints->compute(dipole);

    return dipole;
}

std::vector<SharedMatrix> MintsHelper::ao_quadrupole() {
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole;

    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole XX", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole XY", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole XZ", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole YY", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole YZ", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Quadrupole ZZ", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector<SharedMatrix> MintsHelper::ao_traceless_quadrupole() {
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole;

    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole XX", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole XY", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole XZ", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole YY", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole YZ", basisset_->nbf(), basisset_->nbf()));
    quadrupole.push_back(std::make_shared<Matrix>("AO Traceless Quadrupole ZZ", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_traceless_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector<SharedMatrix> MintsHelper::ao_efp_multipole_potential(const std::vector<double> &origin, int deriv) {
    if (origin.size() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    Vector3 v3origin(origin[0], origin[1], origin[2]);

    std::vector<SharedMatrix> mult;
    mult.push_back(std::make_shared<Matrix>("AO EFP Charge 0", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole X", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole Y", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole Z", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XX", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole YY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole ZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole YZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXX", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YYY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole ZZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XYY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YYZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XYZ", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_efp_multipole_potential(deriv));
    ints->set_origin(v3origin);
    ints->compute(mult);

    return mult;
}

std::vector<SharedMatrix> MintsHelper::electric_field(const std::vector<double> &origin, int deriv) {
    if (origin.size() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    Vector3 v3origin(origin[0], origin[1], origin[2]);

    std::vector<SharedMatrix> field;
    field.push_back(std::make_shared<Matrix>("Ex integrals", basisset_->nbf(), basisset_->nbf()));
    field.push_back(std::make_shared<Matrix>("Ey integrals", basisset_->nbf(), basisset_->nbf()));
    field.push_back(std::make_shared<Matrix>("Ez integrals", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->electric_field(deriv));
    ints->set_origin(v3origin);
    ints->compute(field);

    return field;
}

std::vector<SharedMatrix> MintsHelper::ao_nabla() {
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla;

    nabla.push_back(std::make_shared<Matrix>("AO Px", basisset_->nbf(), basisset_->nbf()));
    nabla.push_back(std::make_shared<Matrix>("AO Py", basisset_->nbf(), basisset_->nbf()));
    nabla.push_back(std::make_shared<Matrix>("AO Pz", basisset_->nbf(), basisset_->nbf()));

    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_nabla());
    ints->compute(nabla);

    return nabla;
}

std::shared_ptr<CdSalcList> MintsHelper::cdsalcs(int needed_irreps, bool project_out_translations,
                                                 bool project_out_rotations) {
    return std::make_shared<CdSalcList>(molecule_, needed_irreps, project_out_translations, project_out_rotations);
}

SharedMatrix MintsHelper::mo_transform(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                       SharedMatrix C4) {
    // Attempts to transform integrals in the most efficient manner. Will transpose left, right
    // and left_right where left and right are (left|right) indices. Does not consider the fimal
    // perturbation eg (12|34) -> (13|24), therefore integrals of type (oo|vv) will not be computed
    // in the optimal order. However, the first transformed index is guaranteed to be the smallest.

    if ((C1->nirrep() + C2->nirrep() + C3->nirrep() + C4->nirrep()) > 4) {
        throw PSIEXCEPTION("MO Transform: Incoming orbitals must be C1 symmetry.");
    }
    int nso = C1->rowspi()[0];

    // Check C dimensions
    int dim_check = 0;
    dim_check += (C2->rowspi()[0] != nso);
    dim_check += (C3->rowspi()[0] != nso);
    dim_check += (C4->rowspi()[0] != nso);

    if (dim_check) {
        throw PSIEXCEPTION("MO Transform: Eigenvector lengths of the C matrices are not identical.");
    }

    if ((Iso->nirrep()) > 1) {
        throw PSIEXCEPTION("MO Transform: The ERI has more than one irrep.");
    }

    // Make sure I is square and of correction dimesion
    int Irows = Iso->rowspi()[0];
    int Icols = Iso->colspi()[0];

    if (((nso * nso) != Irows) || ((nso * nso) != Icols)) {
        throw PSIEXCEPTION("MO Transform: ERI shape does match that of the C matrices.");
    }

    int n1 = C1->colspi()[0];
    int n2 = C2->colspi()[0];
    int n3 = C3->colspi()[0];
    int n4 = C4->colspi()[0];
    std::vector<int> nshape{n1, n2, n3, n4};

    double **C1p = C1->pointer();
    double **C2p = C2->pointer();
    double **C3p = C3->pointer();
    double **C4p = C4->pointer();

    int shape_left = n1 * n2;
    int shape_right = n3 * n4;

    // Transform smallest indices first
    bool transpose_left = (n1 > n2);
    bool transpose_right = (n3 > n4);
    bool transpose_left_right = ((transpose_left ? n2 : n1) > (transpose_right ? n4 : n3));

    int tmp_n;
    double **tmp_p;
    if (transpose_left) {
        tmp_n = n1;
        n1 = n2;
        n2 = tmp_n;
        tmp_p = C1p;
        C1p = C2p;
        C2p = tmp_p;
    }
    if (transpose_right) {
        tmp_n = n3;
        n3 = n4;
        n4 = tmp_n;
        tmp_p = C3p;
        C3p = C4p;
        C4p = tmp_p;
    }
    if (transpose_left_right) {
        tmp_n = n1;
        n1 = n3;
        n3 = tmp_n;
        tmp_n = n2;
        n2 = n4;
        n4 = tmp_n;
        tmp_p = C1p;
        C1p = C3p;
        C3p = tmp_p;
        tmp_p = C2p;
        C2p = C4p;
        C4p = tmp_p;
    }

    double **Isop = Iso->pointer();
    auto I2 = std::make_shared<Matrix>("MO ERI Tensor", n1 * nso, nso * nso);
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', n1, nso * (size_t)nso * nso, nso, 1.0, C1p[0], n1, Isop[0], nso * (size_t)nso * nso, 0.0, I2p[0],
            nso * (size_t)nso * nso);

    Iso.reset();
    auto I3 = std::make_shared<Matrix>("MO ERI Tensor", n1 * nso, nso * n3);
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', n1 * (size_t)nso * nso, n3, nso, 1.0, I2p[0], nso, C3p[0], n3, 0.0, I3p[0], n3);

    I2.reset();
    auto I4 = std::make_shared<Matrix>("MO ERI Tensor", nso * n1, n3 * nso);
    double **I4p = I4->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * n1 + i][j * nso + n] = I3p[i * nso + m][n * n3 + j];
                }
            }
        }
    }

    I3.reset();
    auto I5 = std::make_shared<Matrix>("MO ERI Tensor", n2 * n1, n3 * nso);
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', n2, n1 * (size_t)n3 * nso, nso, 1.0, C2p[0], n2, I4p[0], n1 * (size_t)n3 * nso, 0.0, I5p[0],
            n1 * (size_t)n3 * nso);

    I4.reset();
    auto I6 = std::make_shared<Matrix>("MO ERI Tensor", n2 * n1, n3 * n4);
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', n2 * (size_t)n1 * n3, n4, nso, 1.0, I5p[0], nso, C4p[0], n4, 0.0, I6p[0], n4);

    I5.reset();
    auto Imo = std::make_shared<Matrix>("MO ERI Tensor", shape_left, shape_right);
    double **Imop = Imo->pointer();

    // Currently 2143, need to transform back
    int left, right;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int a = 0; a < n2; a++) {
                // Tranpose left
                if (transpose_left) {
                    left = a * n1 + i;
                } else {
                    left = i * n2 + a;
                }

                for (int b = 0; b < n4; b++) {
                    right = j * n4 + b;

                    // Transpose right
                    if (transpose_right) {
                        right = b * n3 + j;
                    } else {
                        right = j * n4 + b;
                    }

                    // Transpose left_right
                    if (transpose_left_right) {
                        Imop[right][left] = I6p[a * n1 + i][j * n4 + b];
                    } else {
                        Imop[left][right] = I6p[a * n1 + i][j * n4 + b];
                    }
                }
            }
        }
    }
    // Build numpy and final matrix shape
    Imo->set_numpy_shape(nshape);

    return Imo;
}
SharedMatrix MintsHelper::potential_grad(SharedMatrix D) {
    // Potential derivs
    int natom = basisset_->molecule()->natom();
    auto V = std::make_shared<Matrix>("Potential Gradient", natom, 3);

    // Build temps
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    std::vector<SharedMatrix> Vtemps;
    for (size_t i = 0; i < nthread_; i++) {
        Vtemps.push_back(V->clone());
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
    }

    // Lower Triangle
    std::vector<std::pair<int, int>> PQ_pairs;
    for (int P = 0; P < basisset_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

    double **Dp = D->pointer();

#pragma omp parallel for schedule(dynamic) num_threads(nthread_)
    for (size_t PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
        size_t P = PQ_pairs[PQ].first;
        size_t Q = PQ_pairs[PQ].second;

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        ints_vec[rank]->compute_shell_deriv1(P, Q);
        const double *buffer = ints_vec[rank]->buffer();

        size_t nP = basisset_->shell(P).nfunction();
        size_t oP = basisset_->shell(P).function_index();
        size_t aP = basisset_->shell(P).ncenter();

        size_t nQ = basisset_->shell(Q).nfunction();
        size_t oQ = basisset_->shell(Q).function_index();
        size_t aQ = basisset_->shell(Q).ncenter();

        double perm = (P == Q ? 1.0 : 2.0);

        double **Vp = Vtemps[rank]->pointer();

        for (size_t A = 0; A < natom; A++) {
            const double *ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
            const double *ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
            const double *ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    double Vval = perm * Dp[p + oP][q + oQ];
                    Vp[A][0] += Vval * (*ref0++);
                    Vp[A][1] += Vval * (*ref1++);
                    Vp[A][2] += Vval * (*ref2++);
                }
            }
        }
    }

    // Sum it up
    for (size_t t = 0; t < nthread_; t++) {
        V->axpy(1.0, Vtemps[t]);
    }
    return V;
}

SharedMatrix MintsHelper::kinetic_grad(SharedMatrix D) {
    // Overlap
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_kinetic(1)));
    }
    SharedMatrix kinetic_mat(new Matrix("Kinetic Gradient", basisset_->molecule()->natom(), 3));
    grad_two_center_computer(ints_vec, D, kinetic_mat);
    return kinetic_mat;
}
SharedMatrix MintsHelper::overlap_grad(SharedMatrix D) {
    // Overlap
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_overlap(1)));
    }
    SharedMatrix overlap_mat(new Matrix("Overlap Gradient", basisset_->molecule()->natom(), 3));
    grad_two_center_computer(ints_vec, D, overlap_mat);
    return overlap_mat;
}
SharedMatrix MintsHelper::perturb_grad(SharedMatrix D) {
    double xlambda = 0.0;
    double ylambda = 0.0;
    double zlambda = 0.0;

    std::string perturb_with = options_.get_str("PERTURB_WITH");
    if (perturb_with == "DIPOLE_X")
        xlambda = options_.get_double("PERTURB_MAGNITUDE");
    else if (perturb_with == "DIPOLE_Y")
        ylambda = options_.get_double("PERTURB_MAGNITUDE");
    else if (perturb_with == "DIPOLE_Z")
        zlambda = options_.get_double("PERTURB_MAGNITUDE");
    else if (perturb_with == "DIPOLE") {
        if (options_["PERTURB_DIPOLE"].size() != 3)
            throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
        xlambda = options_["PERTURB_DIPOLE"][0].to_double();
        ylambda = options_["PERTURB_DIPOLE"][1].to_double();
        zlambda = options_["PERTURB_DIPOLE"][2].to_double();
    } else {
        std::string msg("Gradients for a ");
        msg += perturb_with;
        msg += " perturbation are not available yet.\n";
        throw PSIEXCEPTION(msg);
    }

    int natoms = basisset_->molecule()->natom();
    auto perturbation_gradient = std::make_shared<Matrix>("Perturbation Gradient", natoms, 3);
    auto dipole_gradients = dipole_grad(D);
    double lambdas[3] = {xlambda, ylambda, zlambda};
    C_DGEMM('n', 't', 3 * natoms, 1, 3, 1.0, dipole_gradients->pointer()[0], 3, &lambdas[0], 3, 0.0,
            perturbation_gradient->pointer()[0], 1);
    return perturbation_gradient;
}

SharedMatrix MintsHelper::dipole_grad(SharedMatrix D) {
    // Computes skeleton (Hellman-Feynman like) dipole derivatives for each perturbation
    double **Dp = D->pointer();

    int natom = molecule_->natom();
    auto ret = std::make_shared<Matrix>("Dipole dervatives (pert*component, i.e. 3Nx3)", 3 * natom, 3);
    double **Pp = ret->pointer();

    std::shared_ptr<OneBodyAOInt> Dint(integral_->ao_dipole(1));
    const double *buffer = Dint->buffer();

    for (int P = 0; P < basisset_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            Dint->compute_shell_deriv1(P, Q);

            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();

            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            const double *ref = buffer;
            double prefac = (P == Q ? 1.0 : 2.0);

            /*
             * Mu X derivatives
             */
            // Px
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 0][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Py
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 1][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Pz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 2][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 0][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 1][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 2][0] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            /*
             * Mu Y derivatives
             */
            // Px
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 0][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Py
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 1][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Pz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 2][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 0][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 1][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 2][1] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            /*
             * Mu Z derivatives
             */
            // Px
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 0][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Py
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 1][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Pz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aP + 2][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 0][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 1][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Pp[3 * aQ + 2][2] += prefac * Dp[p + oP][q + oQ] * (*ref++);
                }
            }
        }
    }
    return ret;
}

SharedMatrix MintsHelper::core_hamiltonian_grad(SharedMatrix D) {
    auto ret = kinetic_grad(D);
    ret->set_name("Core Hamiltonian Gradient");
    ret->add(potential_grad(D));

    if (options_.get_bool("PERTURB_H")) {
        ret->add(perturb_grad(D));
    }
    return ret;
}

void MintsHelper::play() {}

/* 1st and 2nd derivatives of OEI in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_overlap_kinetic_deriv1_helper(const std::string &type, int atom) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::shared_ptr<OneBodyAOInt> GInt;

    if (type == "OVERLAP") {
        std::shared_ptr<OneBodyAOInt> Int(integral_->ao_overlap(1));
        GInt = Int;
    } else {
        std::shared_ptr<OneBodyAOInt> Int(integral_->ao_kinetic(1));
        GInt = Int;
    }

    std::shared_ptr<BasisSet> bs1 = GInt->basis1();
    std::shared_ptr<BasisSet> bs2 = GInt->basis2();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_" << type << "_deriv1_" << atom << cartcomp[p];
        grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1, nbf2)));
    }

    const double *buffer = GInt->buffer();

    for (int P = 0; P < bs1->nshell(); P++)
        for (int Q = 0; Q < bs2->nshell(); Q++) {
            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();

            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            if (aP != atom && aQ != atom) continue;

            GInt->compute_shell_deriv1(P, Q);
            int offset = 0;

            if (aP == atom) {
                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[0]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[1]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[2]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;
            } else {
                offset += 3 * nP * nQ;
            }

            if (aQ == atom) {
                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[0]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[1]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        grad[2]->add(p + oP, q + oQ, buffer[p * nQ + q + offset]);
                    }
                }
                offset += nP * nQ;
            }

            else {
                offset += 3 * nP * nQ;
            }
        }

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_potential_deriv1_helper(int atom) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(1));

    std::shared_ptr<BasisSet> bs1 = Vint->basis1();
    std::shared_ptr<BasisSet> bs2 = Vint->basis2();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();

    int natom = basisset_->molecule()->natom();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_potential_deriv1_" << atom << cartcomp[p];
        grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1, nbf2)));
    }

    const double *buffer = Vint->buffer();

    for (int P = 0; P < bs1->nshell(); P++)
        for (int Q = 0; Q < bs2->nshell(); Q++) {
            int nP = bs1->shell(P).nfunction();
            int oP = bs1->shell(P).function_index();
            int aP = bs1->shell(P).ncenter();

            int nQ = bs2->shell(Q).nfunction();
            int oQ = bs2->shell(Q).function_index();
            int aQ = bs2->shell(Q).ncenter();

            Vint->compute_shell_deriv1(P, Q);

            const double *ref0 = &buffer[3 * atom * nP * nQ + 0 * nP * nQ];
            const double *ref1 = &buffer[3 * atom * nP * nQ + 1 * nP * nQ];
            const double *ref2 = &buffer[3 * atom * nP * nQ + 2 * nP * nQ];
            for (int p = 0; p < nP; p++)
                for (int q = 0; q < nQ; q++) {
                    grad[0]->set(p + oP, q + oQ, (*ref0++));
                    grad[1]->set(p + oP, q + oQ, (*ref1++));
                    grad[2]->set(p + oP, q + oQ, (*ref2++));
                }
        }

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_potential_deriv2_helper(int atom1, int atom2) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("x");
    cartcomp.push_back("y");
    cartcomp.push_back("z");

    std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(2));

    std::shared_ptr<BasisSet> bs1 = Vint->basis1();
    std::shared_ptr<BasisSet> bs2 = Vint->basis2();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int natom = molecule_->natom();

    std::vector<SharedMatrix> grad;
    for (int a = 0, ab = 0; a < 3; a++)
        for (int b = 0; b < 3; b++, ab++) {
            std::stringstream sstream;
            sstream << "ao_potential_deriv2_" << atom1 << atom2 << cartcomp[a] << cartcomp[b];
            grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1, nbf2)));
            grad[ab]->zero();
        }

    const double *buffer = Vint->buffer();

    for (int P = 0; P < bs1->nshell(); P++) {
        for (int Q = 0; Q < bs2->nshell(); Q++) {
            int nP = bs1->shell(P).nfunction();
            int oP = bs1->shell(P).function_index();
            int aP = bs1->shell(P).ncenter();

            int nQ = bs2->shell(Q).nfunction();
            int oQ = bs2->shell(Q).function_index();
            int aQ = bs2->shell(Q).ncenter();

            size_t offset = nP * nQ;

            std::map<std::string, double> grad_map;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    std::string grad_key = cartcomp[a] + cartcomp[b];
                    grad_map[grad_key] = 0;
                }

            for (int atom = 0; atom < natom; atom++) {
                double Z = molecule_->Z(atom);
                Vint->set_origin(molecule_->xyz(atom));

                std::pair<int, int> target_atoms(atom1, atom2);
                std::vector<std::pair<int, int>> vec_pairs;
                std::vector<int> vec_pos;
                std::map<int, std::string> pos_map;
                std::string strings;
                std::string key;

                vec_pairs.push_back(std::make_pair(aP, aP));
                vec_pairs.push_back(std::make_pair(aQ, aQ));
                vec_pairs.push_back(std::make_pair(atom, atom));
                vec_pairs.push_back(std::make_pair(atom, aP));
                vec_pairs.push_back(std::make_pair(aP, aQ));
                vec_pairs.push_back(std::make_pair(atom, aQ));

                for (std::vector<std::pair<int, int>>::iterator it = vec_pairs.begin(); it != vec_pairs.end(); ++it) {
                    if ((*it).first == target_atoms.first && (*it).second == target_atoms.second)
                        vec_pos.push_back(it - vec_pairs.begin());
                }

                pos_map[0] = "AA";
                pos_map[1] = "BB";
                pos_map[2] = "CC";
                pos_map[3] = "CA";
                pos_map[4] = "AB";
                pos_map[5] = "CB";

                if (vec_pos.empty()) continue;

                Vint->compute_shell_deriv2(P, Q);

                const double *CxAx = buffer + 0 * offset;
                const double *CxAy = buffer + 1 * offset;
                const double *CxAz = buffer + 2 * offset;
                const double *CyAx = buffer + 3 * offset;
                const double *CyAy = buffer + 4 * offset;
                const double *CyAz = buffer + 5 * offset;
                const double *CzAx = buffer + 6 * offset;
                const double *CzAy = buffer + 7 * offset;
                const double *CzAz = buffer + 8 * offset;
                const double *AxAx = buffer + 9 * offset;
                const double *AxAy = buffer + 10 * offset;
                const double *AxAz = buffer + 11 * offset;
                const double *AyAy = buffer + 12 * offset;
                const double *AyAz = buffer + 13 * offset;
                const double *AzAz = buffer + 14 * offset;
                const double *BxBx = buffer + 15 * offset;
                const double *BxBy = buffer + 16 * offset;
                const double *BxBz = buffer + 17 * offset;
                const double *ByBy = buffer + 18 * offset;
                const double *ByBz = buffer + 19 * offset;
                const double *BzBz = buffer + 20 * offset;
                const double *CxCx = buffer + 21 * offset;
                const double *CxCy = buffer + 22 * offset;
                const double *CxCz = buffer + 23 * offset;
                const double *CyCy = buffer + 24 * offset;
                const double *CyCz = buffer + 25 * offset;
                const double *CzCz = buffer + 26 * offset;

                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        std::map<std::string, double> hess_map;

                        hess_map["CxAx"] = Z * (*CxAx);
                        hess_map["CxAy"] = Z * (*CxAy);
                        hess_map["CxAz"] = Z * (*CxAz);
                        hess_map["CyAx"] = Z * (*CyAx);
                        hess_map["CyAy"] = Z * (*CyAy);
                        hess_map["CyAz"] = Z * (*CyAz);
                        hess_map["CzAx"] = Z * (*CzAx);
                        hess_map["CzAy"] = Z * (*CzAy);
                        hess_map["CzAz"] = Z * (*CzAz);
                        hess_map["AxAx"] = Z * (*AxAx);
                        hess_map["AxAy"] = Z * (*AxAy);
                        hess_map["AxAz"] = Z * (*AxAz);
                        hess_map["AyAy"] = Z * (*AyAy);
                        hess_map["AyAz"] = Z * (*AyAz);
                        hess_map["AzAz"] = Z * (*AzAz);
                        hess_map["BxBx"] = Z * (*BxBx);
                        hess_map["BxBy"] = Z * (*BxBy);
                        hess_map["BxBz"] = Z * (*BxBz);
                        hess_map["ByBy"] = Z * (*ByBy);
                        hess_map["ByBz"] = Z * (*ByBz);
                        hess_map["BzBz"] = Z * (*BzBz);
                        hess_map["CxCx"] = Z * (*CxCx);
                        hess_map["CxCy"] = Z * (*CxCy);
                        hess_map["CxCz"] = Z * (*CxCz);
                        hess_map["CyCy"] = Z * (*CyCy);
                        hess_map["CyCz"] = Z * (*CyCz);
                        hess_map["CzCz"] = Z * (*CzCz);

                        hess_map["AyAx"] = hess_map["AxAy"];
                        hess_map["AzAx"] = hess_map["AxAz"];
                        hess_map["AzAy"] = hess_map["AyAz"];
                        hess_map["ByBx"] = hess_map["BxBy"];
                        hess_map["BzBx"] = hess_map["BxBz"];
                        hess_map["BzBy"] = hess_map["ByBz"];
                        hess_map["CyCx"] = hess_map["CxCy"];
                        hess_map["CzCx"] = hess_map["CxCz"];
                        hess_map["CzCy"] = hess_map["CyCz"];

                        for (auto pos : vec_pos) {
                            strings = pos_map[pos];
                            for (int a = 0; a < 3; a++)
                                for (int b = 0; b < 3; b++) {
                                    key = strings[0] + cartcomp[a] + strings[1] + cartcomp[b];
                                    std::string grad_key = cartcomp[a] + cartcomp[b];
                                    if (pos < 3 && a <= b) {
                                        grad_map[grad_key] += hess_map[key];
                                    }
                                    if (pos == 3) {
                                        if (aP == atom && a == b)
                                            grad_map[grad_key] += 2.0 * hess_map[key];
                                        else
                                            grad_map[grad_key] += hess_map[key];
                                    }
                                    if (pos == 4) {
                                        std::string key1 = pos_map[2][0] + cartcomp[a] + pos_map[2][1] + cartcomp[b];
                                        std::string key2 = pos_map[3][0] + cartcomp[a] + pos_map[3][1] + cartcomp[b];
                                        std::string key3 = pos_map[1][0] + cartcomp[a] + pos_map[1][1] + cartcomp[b];
                                        if (aP == aQ && a == b)
                                            grad_map[grad_key] +=
                                                2.0 * (hess_map[key1] + hess_map[key2] - hess_map[key3]);
                                        else
                                            grad_map[grad_key] += hess_map[key1] + hess_map[key2] - hess_map[key3];
                                    }
                                    if (pos == 5) {
                                        std::string key1 = pos_map[3][0] + cartcomp[b] + pos_map[3][1] + cartcomp[a];
                                        std::string key2 = pos_map[0][0] + cartcomp[a] + pos_map[0][1] + cartcomp[b];
                                        std::string key3 = pos_map[1][0] + cartcomp[a] + pos_map[1][1] + cartcomp[b];
                                        if (aQ == atom && a == b)
                                            grad_map[grad_key] +=
                                                2.0 * (hess_map[key1] + hess_map[key2] - hess_map[key3]);
                                        else
                                            grad_map[grad_key] += hess_map[key1] + hess_map[key2] - hess_map[key3];
                                    }
                                }
                        }

                        for (int a = 0, ab = 0; a < 3; a++)
                            for (int b = 0; b < 3; b++, ab++) {
                                std::string grad_key = cartcomp[a] + cartcomp[b];
                                grad[ab]->add(p + oP, q + oQ, grad_map[grad_key]);
                                grad_map[grad_key] = 0;
                            }

                        ++CxAx;
                        ++CxAy;
                        ++CxAz;
                        ++CyAx;
                        ++CyAy;
                        ++CyAz;
                        ++CzAx;
                        ++CzAy;
                        ++CzAz;
                        ++AxAx;
                        ++AxAy;
                        ++AxAz;
                        ++AyAy;
                        ++AyAz;
                        ++AzAz;
                        ++BxBx;
                        ++BxBy;
                        ++BxBz;
                        ++ByBy;
                        ++ByBz;
                        ++BzBz;
                        ++CxCx;
                        ++CxCy;
                        ++CxCz;
                        ++CyCy;
                        ++CyCz;
                        ++CzCz;
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2};
    for (int p = 0; p < 9; p++) grad[p]->set_numpy_shape(nshape);

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_overlap_kinetic_deriv2_helper(const std::string &type, int atom1, int atom2) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::shared_ptr<OneBodyAOInt> GInt;

    if (type == "OVERLAP") {
        std::shared_ptr<OneBodyAOInt> Int(integral_->ao_overlap(2));
        GInt = Int;
    } else {
        std::shared_ptr<OneBodyAOInt> Int(integral_->ao_kinetic(2));
        GInt = Int;
    }

    std::shared_ptr<BasisSet> bs1 = GInt->basis1();
    std::shared_ptr<BasisSet> bs2 = GInt->basis2();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
            std::stringstream sstream;
            sstream << "ao_" << type << "_deriv2_" << atom1 << atom2 << cartcomp[p] << cartcomp[q];
            grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1, nbf2)));
        }

    const double *buffer = GInt->buffer();

    for (int P = 0; P < bs1->nshell(); P++) {
        for (int Q = 0; Q < bs2->nshell(); Q++) {
            int nP = bs1->shell(P).nfunction();
            int oP = bs1->shell(P).function_index();
            int aP = bs1->shell(P).ncenter();

            int nQ = bs2->shell(Q).nfunction();
            int oQ = bs2->shell(Q).function_index();
            int aQ = bs2->shell(Q).ncenter();

            if (aP != atom1 && aQ != atom1 && aP != atom2 && aQ != atom2) continue;

            GInt->compute_shell_deriv2(P, Q);

            size_t offset = nP * nQ;

            const double *pxx = buffer + 0 * offset;
            const double *pxy = buffer + 1 * offset;
            const double *pxz = buffer + 2 * offset;
            const double *pyy = buffer + 3 * offset;
            const double *pyz = buffer + 4 * offset;
            const double *pzz = buffer + 5 * offset;

            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    double tmpxx = (*pxx);
                    double tmpxy = (*pxy);
                    double tmpxz = (*pxz);
                    double tmpyy = (*pyy);
                    double tmpyz = (*pyz);
                    double tmpzz = (*pzz);

                    if (atom1 == atom2) {
                        if (aP == aQ) {
                            grad[0]->set(p + oP, q + oQ, 0);
                            grad[1]->set(p + oP, q + oQ, tmpxy);
                            grad[2]->set(p + oP, q + oQ, tmpxz);
                            grad[3]->set(p + oP, q + oQ, -tmpxy);
                            grad[4]->set(p + oP, q + oQ, 0);
                            grad[5]->set(p + oP, q + oQ, tmpyz);
                            grad[6]->set(p + oP, q + oQ, -tmpxz);
                            grad[7]->set(p + oP, q + oQ, -tmpyz);
                            grad[8]->set(p + oP, q + oQ, 0);
                        } else {
                            grad[0]->set(p + oP, q + oQ, tmpxx);
                            grad[1]->set(p + oP, q + oQ, tmpxy);
                            grad[2]->set(p + oP, q + oQ, tmpxz);
                            grad[4]->set(p + oP, q + oQ, tmpyy);
                            grad[5]->set(p + oP, q + oQ, tmpyz);
                            grad[8]->set(p + oP, q + oQ, tmpzz);
                        }
                    } else {
                        if (aP == atom1 && aQ == atom2) {
                            grad[0]->set(p + oP, q + oQ, -1.0 * tmpxx);
                            grad[1]->set(p + oP, q + oQ, -1.0 * tmpxy);
                            grad[2]->set(p + oP, q + oQ, -1.0 * tmpxz);
                            grad[3]->set(p + oP, q + oQ, -1.0 * tmpxy);
                            grad[4]->set(p + oP, q + oQ, -1.0 * tmpyy);
                            grad[5]->set(p + oP, q + oQ, -1.0 * tmpyz);
                            grad[6]->set(p + oP, q + oQ, -1.0 * tmpxz);
                            grad[7]->set(p + oP, q + oQ, -1.0 * tmpyz);
                            grad[8]->set(p + oP, q + oQ, -1.0 * tmpzz);
                        }
                    }

                    ++pxx;
                    ++pxy;
                    ++pxz;
                    ++pyy;
                    ++pyz;
                    ++pzz;
                }
            }
        }
    }

    return grad;
}

/* 1st and 2nd derivatives of TEI in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_tei_deriv1(int atom, double omega,
                                                     std::shared_ptr<IntegralFactory> input_factory) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::shared_ptr<IntegralFactory> factory;
    if (input_factory) {
        factory = input_factory;
    } else {
        factory = integral_;
    }

    std::shared_ptr<TwoBodyAOInt> ints;
    if (omega == 0.0) {
        ints = std::shared_ptr<TwoBodyAOInt>(factory->eri(1));
    } else {
        ints = std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega, 1));
    }

    std::shared_ptr<BasisSet> bs1 = ints->basis1();
    std::shared_ptr<BasisSet> bs2 = ints->basis2();
    std::shared_ptr<BasisSet> bs3 = ints->basis3();
    std::shared_ptr<BasisSet> bs4 = ints->basis4();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();
    int nbf4 = bs4->nbf();

    int natom = basisset_->molecule()->natom();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_tei_deriv1_" << atom << cartcomp[p];
        grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1 * nbf2, nbf3 * nbf4)));
    }

    const double *buffer = ints->buffer();

    for (int P = 0; P < bs1->nshell(); P++) {
        for (int Q = 0; Q < bs2->nshell(); Q++) {
            for (int R = 0; R < bs3->nshell(); R++) {
                for (int S = 0; S < bs4->nshell(); S++) {
                    int Psize = bs1->shell(P).nfunction();
                    int Qsize = bs2->shell(Q).nfunction();
                    int Rsize = bs3->shell(R).nfunction();
                    int Ssize = bs4->shell(S).nfunction();

                    int Pncart = bs1->shell(P).ncartesian();
                    int Qncart = bs2->shell(Q).ncartesian();
                    int Rncart = bs3->shell(R).ncartesian();
                    int Sncart = bs4->shell(S).ncartesian();

                    int Poff = bs1->shell(P).function_index();
                    int Qoff = bs2->shell(Q).function_index();
                    int Roff = bs3->shell(R).function_index();
                    int Soff = bs4->shell(S).function_index();

                    int Pcenter = bs1->shell(P).ncenter();
                    int Qcenter = bs2->shell(Q).ncenter();
                    int Rcenter = bs3->shell(R).ncenter();
                    int Scenter = bs4->shell(S).ncenter();

                    size_t stride = Pncart * Qncart * Rncart * Sncart;
                    size_t delta;

                    delta = 0L;

                    if (Pcenter != atom && Qcenter != atom && Rcenter != atom && Scenter != atom) continue;

                    if (Pcenter == atom && Qcenter == atom && Rcenter == atom && Scenter == atom) continue;

                    ints->compute_shell_deriv1(P, Q, R, S);

                    double Ax, Ay, Az;
                    double Bx, By, Bz;
                    double Cx, Cy, Cz;
                    double Dx, Dy, Dz;
                    double X = 0, Y = 0, Z = 0;

                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    int i = (Poff + p) * nbf2 + Qoff + q;
                                    int j = (Roff + r) * nbf4 + Soff + s;

                                    Ax = buffer[0 * stride + delta];
                                    Ay = buffer[1 * stride + delta];
                                    Az = buffer[2 * stride + delta];
                                    Cx = buffer[3 * stride + delta];
                                    Cy = buffer[4 * stride + delta];
                                    Cz = buffer[5 * stride + delta];
                                    Dx = buffer[6 * stride + delta];
                                    Dy = buffer[7 * stride + delta];
                                    Dz = buffer[8 * stride + delta];

                                    Bx = -(Ax + Cx + Dx);
                                    By = -(Ay + Cy + Dy);
                                    Bz = -(Az + Cz + Dz);

                                    if (Pcenter == atom) {
                                        X += Ax;
                                        Y += Ay;
                                        Z += Az;
                                    }

                                    if (Qcenter == atom) {
                                        X += Bx;
                                        Y += By;
                                        Z += Bz;
                                    }

                                    if (Rcenter == atom) {
                                        X += Cx;
                                        Y += Cy;
                                        Z += Cz;
                                    }

                                    if (Scenter == atom) {
                                        X += Dx;
                                        Y += Dy;
                                        Z += Dz;
                                    }

                                    grad[0]->set(i, j, X);
                                    grad[1]->set(i, j, Y);
                                    grad[2]->set(i, j, Z);

                                    X = 0, Y = 0, Z = 0;
                                    delta++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3, nbf4};
    for (int p = 0; p < 3; p++) grad[p]->set_numpy_shape(nshape);

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_tei_deriv2(int atom1, int atom2) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("x");
    cartcomp.push_back("y");
    cartcomp.push_back("z");

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    for (int thread = 0; thread < nthreads; thread++) ints.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->eri(2)));

    std::shared_ptr<BasisSet> bs1 = ints[0]->basis1();
    std::shared_ptr<BasisSet> bs2 = ints[0]->basis2();
    std::shared_ptr<BasisSet> bs3 = ints[0]->basis3();
    std::shared_ptr<BasisSet> bs4 = ints[0]->basis4();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();
    int nbf4 = bs4->nbf();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
            std::stringstream sstream;
            sstream << "ao_tei_deriv2_" << atom1 << atom2 << cartcomp[p] << cartcomp[q];
            grad.push_back(SharedMatrix(new Matrix(sstream.str(), nbf1 * nbf2, nbf3 * nbf4)));
        }

    std::vector<std::vector<int>> shell_quartets;

    for (int P = 0; P < bs1->nshell(); P++)
        for (int Q = 0; Q < bs2->nshell(); Q++)
            for (int R = 0; R < bs3->nshell(); R++)
                for (int S = 0; S < bs4->nshell(); S++) {
                    std::vector<int> tmp;
                    tmp.push_back(P);
                    tmp.push_back(Q);
                    tmp.push_back(R);
                    tmp.push_back(S);
                    shell_quartets.push_back(tmp);
                }

#pragma omp parallel for num_threads(nthreads) schedule(dynamic)

    for (int i = 0; i < shell_quartets.size(); i++) {
        int P = shell_quartets[i][0];
        int Q = shell_quartets[i][1];
        int R = shell_quartets[i][2];
        int S = shell_quartets[i][3];

        int Psize = bs1->shell(P).nfunction();
        int Qsize = bs2->shell(Q).nfunction();
        int Rsize = bs3->shell(R).nfunction();
        int Ssize = bs4->shell(S).nfunction();

        int Pncart = bs1->shell(P).ncartesian();
        int Qncart = bs2->shell(Q).ncartesian();
        int Rncart = bs3->shell(R).ncartesian();
        int Sncart = bs4->shell(S).ncartesian();

        int Poff = bs1->shell(P).function_index();
        int Qoff = bs2->shell(Q).function_index();
        int Roff = bs3->shell(R).function_index();
        int Soff = bs4->shell(S).function_index();

        int Pcenter = bs1->shell(P).ncenter();
        int Qcenter = bs2->shell(Q).ncenter();
        int Rcenter = bs3->shell(R).ncenter();
        int Scenter = bs4->shell(S).ncenter();

        size_t stride = Pncart * Qncart * Rncart * Sncart;
        size_t delta;

        delta = 0L;

        std::pair<int, int> atoms(atom1, atom2);
        std::vector<std::pair<int, int>> vec_pairs;
        std::vector<int> vec_pos;
        std::map<int, std::string> pos_map;
        std::string strings;
        std::string key;
        std::map<std::string, double> grad_map;

        vec_pairs.push_back(std::make_pair(Pcenter, Pcenter));
        vec_pairs.push_back(std::make_pair(Qcenter, Qcenter));
        vec_pairs.push_back(std::make_pair(Rcenter, Rcenter));
        vec_pairs.push_back(std::make_pair(Scenter, Scenter));
        vec_pairs.push_back(std::make_pair(Pcenter, Qcenter));
        vec_pairs.push_back(std::make_pair(Pcenter, Rcenter));
        vec_pairs.push_back(std::make_pair(Pcenter, Scenter));
        vec_pairs.push_back(std::make_pair(Qcenter, Rcenter));
        vec_pairs.push_back(std::make_pair(Qcenter, Scenter));
        vec_pairs.push_back(std::make_pair(Rcenter, Scenter));

        for (std::vector<std::pair<int, int>>::iterator it = vec_pairs.begin(); it != vec_pairs.end(); ++it) {
            if ((*it).first == atoms.first && (*it).second == atoms.second) vec_pos.push_back(it - vec_pairs.begin());
        }

        pos_map[0] = "AA";
        pos_map[1] = "BB";
        pos_map[2] = "CC";
        pos_map[3] = "DD";
        pos_map[4] = "AB";
        pos_map[5] = "AC";
        pos_map[6] = "AD";
        pos_map[7] = "BC";
        pos_map[8] = "BD";
        pos_map[9] = "CD";

        for (int p = 0, pq = 0; p < 3; p++)
            for (int q = 0; q < 3; q++, pq++) {
                std::string grad_key = cartcomp[p] + cartcomp[q];
                grad_map[grad_key] = 0;
            }

        if (vec_pos.empty()) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        ints[thread]->compute_shell_deriv2(P, Q, R, S);

        const double *buffer = ints[thread]->buffer();
        std::unordered_map<std::string, double> hess_map;

        for (int p = 0; p < Psize; p++) {
            for (int q = 0; q < Qsize; q++) {
                for (int r = 0; r < Rsize; r++) {
                    for (int s = 0; s < Ssize; s++) {
                        int i = (Poff + p) * nbf2 + Qoff + q;
                        int j = (Roff + r) * nbf4 + Soff + s;

                        hess_map["AxAx"] = buffer[9 * stride + delta];
                        hess_map["AxAy"] = buffer[10 * stride + delta];
                        hess_map["AxAz"] = buffer[11 * stride + delta];
                        hess_map["AxCx"] = buffer[12 * stride + delta];
                        hess_map["AxCy"] = buffer[13 * stride + delta];
                        hess_map["AxCz"] = buffer[14 * stride + delta];
                        hess_map["AxDx"] = buffer[15 * stride + delta];
                        hess_map["AxDy"] = buffer[16 * stride + delta];
                        hess_map["AxDz"] = buffer[17 * stride + delta];
                        hess_map["AyAy"] = buffer[18 * stride + delta];
                        hess_map["AyAz"] = buffer[19 * stride + delta];
                        hess_map["AyCx"] = buffer[20 * stride + delta];
                        hess_map["AyCy"] = buffer[21 * stride + delta];
                        hess_map["AyCz"] = buffer[22 * stride + delta];
                        hess_map["AyDx"] = buffer[23 * stride + delta];
                        hess_map["AyDy"] = buffer[24 * stride + delta];
                        hess_map["AyDz"] = buffer[25 * stride + delta];
                        hess_map["AzAz"] = buffer[26 * stride + delta];
                        hess_map["AzCx"] = buffer[27 * stride + delta];
                        hess_map["AzCy"] = buffer[28 * stride + delta];
                        hess_map["AzCz"] = buffer[29 * stride + delta];
                        hess_map["AzDx"] = buffer[30 * stride + delta];
                        hess_map["AzDy"] = buffer[31 * stride + delta];
                        hess_map["AzDz"] = buffer[32 * stride + delta];
                        hess_map["CxCx"] = buffer[33 * stride + delta];
                        hess_map["CxCy"] = buffer[34 * stride + delta];
                        hess_map["CxCz"] = buffer[35 * stride + delta];
                        hess_map["CxDx"] = buffer[36 * stride + delta];
                        hess_map["CxDy"] = buffer[37 * stride + delta];
                        hess_map["CxDz"] = buffer[38 * stride + delta];
                        hess_map["CyCy"] = buffer[39 * stride + delta];
                        hess_map["CyCz"] = buffer[40 * stride + delta];
                        hess_map["CyDx"] = buffer[41 * stride + delta];
                        hess_map["CyDy"] = buffer[42 * stride + delta];
                        hess_map["CyDz"] = buffer[43 * stride + delta];
                        hess_map["CzCz"] = buffer[44 * stride + delta];
                        hess_map["CzDx"] = buffer[45 * stride + delta];
                        hess_map["CzDy"] = buffer[46 * stride + delta];
                        hess_map["CzDz"] = buffer[47 * stride + delta];
                        hess_map["DxDx"] = buffer[48 * stride + delta];
                        hess_map["DxDy"] = buffer[49 * stride + delta];
                        hess_map["DxDz"] = buffer[50 * stride + delta];
                        hess_map["DyDy"] = buffer[51 * stride + delta];
                        hess_map["DyDz"] = buffer[52 * stride + delta];
                        hess_map["DzDz"] = buffer[53 * stride + delta];

                        // Translational invariance relationships

                        hess_map["AxBx"] = -(hess_map["AxAx"] + hess_map["AxCx"] + hess_map["AxDx"]);
                        hess_map["AxBy"] = -(hess_map["AxAy"] + hess_map["AxCy"] + hess_map["AxDy"]);
                        hess_map["AxBz"] = -(hess_map["AxAz"] + hess_map["AxCz"] + hess_map["AxDz"]);
                        hess_map["AyBx"] = -(hess_map["AxAy"] + hess_map["AyCx"] + hess_map["AyDx"]);
                        hess_map["AyBy"] = -(hess_map["AyAy"] + hess_map["AyCy"] + hess_map["AyDy"]);
                        hess_map["AyBz"] = -(hess_map["AyAz"] + hess_map["AyCz"] + hess_map["AyDz"]);
                        hess_map["AzBx"] = -(hess_map["AxAz"] + hess_map["AzCx"] + hess_map["AzDx"]);
                        hess_map["AzBy"] = -(hess_map["AyAz"] + hess_map["AzCy"] + hess_map["AzDy"]);
                        hess_map["AzBz"] = -(hess_map["AzAz"] + hess_map["AzCz"] + hess_map["AzDz"]);
                        hess_map["BxCx"] = -(hess_map["AxCx"] + hess_map["CxCx"] + hess_map["CxDx"]);
                        hess_map["BxCy"] = -(hess_map["AxCy"] + hess_map["CxCy"] + hess_map["CyDx"]);
                        hess_map["BxCz"] = -(hess_map["AxCz"] + hess_map["CxCz"] + hess_map["CzDx"]);
                        hess_map["ByCx"] = -(hess_map["AyCx"] + hess_map["CxCy"] + hess_map["CxDy"]);
                        hess_map["ByCy"] = -(hess_map["AyCy"] + hess_map["CyCy"] + hess_map["CyDy"]);
                        hess_map["ByCz"] = -(hess_map["AyCz"] + hess_map["CyCz"] + hess_map["CzDy"]);
                        hess_map["BzCx"] = -(hess_map["AzCx"] + hess_map["CxCz"] + hess_map["CxDz"]);
                        hess_map["BzCy"] = -(hess_map["AzCy"] + hess_map["CyCz"] + hess_map["CyDz"]);
                        hess_map["BzCz"] = -(hess_map["AzCz"] + hess_map["CzCz"] + hess_map["CzDz"]);
                        hess_map["BxDx"] = -(hess_map["AxDx"] + hess_map["CxDx"] + hess_map["DxDx"]);
                        hess_map["BxDy"] = -(hess_map["AxDy"] + hess_map["CxDy"] + hess_map["DxDy"]);
                        hess_map["BxDz"] = -(hess_map["AxDz"] + hess_map["CxDz"] + hess_map["DxDz"]);
                        hess_map["ByDx"] = -(hess_map["AyDx"] + hess_map["CyDx"] + hess_map["DxDy"]);
                        hess_map["ByDy"] = -(hess_map["AyDy"] + hess_map["CyDy"] + hess_map["DyDy"]);
                        hess_map["ByDz"] = -(hess_map["AyDz"] + hess_map["CyDz"] + hess_map["DyDz"]);
                        hess_map["BzDx"] = -(hess_map["AzDx"] + hess_map["CzDx"] + hess_map["DxDz"]);
                        hess_map["BzDy"] = -(hess_map["AzDy"] + hess_map["CzDy"] + hess_map["DyDz"]);
                        hess_map["BzDz"] = -(hess_map["AzDz"] + hess_map["CzDz"] + hess_map["DzDz"]);

                        hess_map["BxBx"] = hess_map["AxAx"] + hess_map["AxCx"] + hess_map["AxDx"] + hess_map["AxCx"] +
                                           hess_map["CxCx"] + hess_map["CxDx"] + hess_map["AxDx"] + hess_map["CxDx"] +
                                           hess_map["DxDx"];

                        hess_map["ByBy"] = hess_map["AyAy"] + hess_map["AyCy"] + hess_map["AyDy"] + hess_map["AyCy"] +
                                           hess_map["CyCy"] + hess_map["CyDy"] + hess_map["AyDy"] + hess_map["CyDy"] +
                                           hess_map["DyDy"];

                        hess_map["BzBz"] = hess_map["AzAz"] + hess_map["AzCz"] + hess_map["AzDz"] + hess_map["AzCz"] +
                                           hess_map["CzCz"] + hess_map["CzDz"] + hess_map["AzDz"] + hess_map["CzDz"] +
                                           hess_map["DzDz"];

                        hess_map["BxBy"] = hess_map["AxAy"] + hess_map["AxCy"] + hess_map["AxDy"] + hess_map["AyCx"] +
                                           hess_map["CxCy"] + hess_map["CxDy"] + hess_map["AyDx"] + hess_map["CyDx"] +
                                           hess_map["DxDy"];

                        hess_map["BxBz"] = hess_map["AxAz"] + hess_map["AxCz"] + hess_map["AxDz"] + hess_map["AzCx"] +
                                           hess_map["CxCz"] + hess_map["CxDz"] + hess_map["AzDx"] + hess_map["CzDx"] +
                                           hess_map["DxDz"];

                        hess_map["ByBz"] = hess_map["AyAz"] + hess_map["AyCz"] + hess_map["AyDz"] + hess_map["AzCy"] +
                                           hess_map["CyCz"] + hess_map["CyDz"] + hess_map["AzDy"] + hess_map["CzDy"] +
                                           hess_map["DyDz"];

                        hess_map["AyAx"] = hess_map["AxAy"];
                        hess_map["AzAx"] = hess_map["AxAz"];
                        hess_map["AzAy"] = hess_map["AyAz"];

                        hess_map["ByBx"] = hess_map["BxBy"];
                        hess_map["BzBx"] = hess_map["BxBz"];
                        hess_map["BzBy"] = hess_map["ByBz"];

                        hess_map["CyCx"] = hess_map["CxCy"];
                        hess_map["CzCx"] = hess_map["CxCz"];
                        hess_map["CzCy"] = hess_map["CyCz"];

                        hess_map["DyDx"] = hess_map["DxDy"];
                        hess_map["DzDx"] = hess_map["DxDz"];
                        hess_map["DzDy"] = hess_map["DyDz"];

                        for (std::vector<int>::iterator it = vec_pos.begin(); it != vec_pos.end(); ++it) {
                            strings = pos_map[*it];
                            for (int p = 0; p < 3; p++)
                                for (int q = 0; q < 3; q++) {
                                    key = strings[0] + cartcomp[p] + strings[1] + cartcomp[q];
                                    std::string grad_key = cartcomp[p] + cartcomp[q];
                                    if (atom1 == atom2 && (*it) <= 3)
                                        grad_map[grad_key] += hess_map[key];
                                    else
                                        grad_map[grad_key] += 2.0 * hess_map[key];
                                }
                        }

                        for (int p = 0, pq = 0; p < 3; p++)
                            for (int q = 0; q < 3; q++, pq++) {
                                std::string grad_key = cartcomp[p] + cartcomp[q];
                                grad[pq]->set(i, j, grad_map[grad_key]);
                                grad_map[grad_key] = 0;
                            }
                        delta++;
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3, nbf4};
    for (int p = 0; p < 9; p++) grad[p]->set_numpy_shape(nshape);

    return grad;
}

/* OEI derivatives in both ao and mo basis */

std::vector<SharedMatrix> MintsHelper::ao_oei_deriv1(const std::string &oei_type, int atom) {
    std::vector<SharedMatrix> ao_grad;

    if (oei_type == "OVERLAP")
        ao_grad = ao_overlap_kinetic_deriv1_helper("OVERLAP", atom);
    else if (oei_type == "KINETIC")
        ao_grad = ao_overlap_kinetic_deriv1_helper("KINETIC", atom);
    else if (oei_type == "POTENTIAL")
        ao_grad = ao_potential_deriv1_helper(atom);
    else
        throw PSIEXCEPTION("Not a valid choice of OEI");

    return ao_grad;
}

std::vector<SharedMatrix> MintsHelper::ao_oei_deriv2(const std::string &oei_type, int atom1, int atom2) {
    std::vector<SharedMatrix> ao_grad_12;
    std::vector<SharedMatrix> ao_grad_21;

    if (oei_type == "OVERLAP") {
        ao_grad_12 = ao_overlap_kinetic_deriv2_helper("OVERLAP", atom1, atom2);
        if (atom1 != atom2) ao_grad_21 = ao_overlap_kinetic_deriv2_helper("OVERLAP", atom2, atom1);
    } else if (oei_type == "KINETIC") {
        ao_grad_12 = ao_overlap_kinetic_deriv2_helper("KINETIC", atom1, atom2);
        if (atom1 != atom2) ao_grad_21 = ao_overlap_kinetic_deriv2_helper("KINETIC", atom2, atom1);
    } else if (oei_type == "POTENTIAL") {
        ao_grad_12 = ao_potential_deriv2_helper(atom1, atom2);
        if (atom1 != atom2) ao_grad_21 = ao_potential_deriv2_helper(atom2, atom1);
    } else
        throw PSIEXCEPTION("Not a valid choice of OEI");

    for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
            int pq = p * 3 + q;
            int qp = q * 3 + p;

            if (atom1 == atom2) {
                if (q < p) {
                    ao_grad_12[pq]->add(ao_grad_12[qp]);
                    ao_grad_12[qp] = ao_grad_12[pq];
                }
            } else
                ao_grad_12[pq]->add(ao_grad_21[qp]);
        }

    return ao_grad_12;
}

std::vector<SharedMatrix> MintsHelper::mo_oei_deriv1(const std::string &oei_type, int atom, SharedMatrix C1,
                                                     SharedMatrix C2) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::vector<SharedMatrix> ao_grad;
    ao_grad = ao_oei_deriv1(oei_type, atom);

    // Assuming C1 symmetry
    int nbf1 = ao_grad[0]->rowdim();
    int nbf2 = ao_grad[0]->coldim();

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "mo_" << oei_type << "_deriv1_" << atom << cartcomp[p];
        SharedMatrix temp(new Matrix(sstream.str(), nbf1, nbf2));
        temp->transform(C1, ao_grad[p], C2);
        mo_grad.push_back(temp);
    }
    return mo_grad;
}

std::vector<SharedMatrix> MintsHelper::mo_oei_deriv2(const std::string &oei_type, int atom1, int atom2, SharedMatrix C1,
                                                     SharedMatrix C2) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::vector<SharedMatrix> ao_grad;
    ao_grad = ao_oei_deriv2(oei_type, atom1, atom2);

    // Assuming C1 symmetry
    int nbf1 = ao_grad[0]->rowdim();
    int nbf2 = ao_grad[0]->coldim();

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0, pq = 0; p < 3; p++)
        for (int q = 0; q < 3; q++, pq++) {
            std::stringstream sstream;
            sstream << "mo_" << oei_type << "_deriv2_" << atom1 << atom2 << cartcomp[p] << cartcomp[q];
            SharedMatrix temp(new Matrix(sstream.str(), nbf1, nbf2));
            temp->transform(C1, ao_grad[pq], C2);
            mo_grad.push_back(temp);
        }
    return mo_grad;
}

/*  TEI derivatives in  MO basis */

std::vector<SharedMatrix> MintsHelper::mo_tei_deriv1(int atom, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                                     SharedMatrix C4) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::vector<SharedMatrix> ao_grad = ao_tei_deriv1(atom);

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "mo_tei_deriv1_" << atom << cartcomp[p];
        SharedMatrix temp = mo_eri_helper(ao_grad[p], C1, C2, C3, C4);
        temp->set_name(sstream.str());
        mo_grad.push_back(temp);
    }
    return mo_grad;
}

std::vector<SharedMatrix> MintsHelper::mo_tei_deriv2(int atom1, int atom2, SharedMatrix C1, SharedMatrix C2,
                                                     SharedMatrix C3, SharedMatrix C4) {
    std::vector<std::string> cartcomp;
    cartcomp.push_back("X");
    cartcomp.push_back("Y");
    cartcomp.push_back("Z");

    std::vector<SharedMatrix> ao_grad = ao_tei_deriv2(atom1, atom2);
    std::vector<SharedMatrix> mo_grad;
    for (int p = 0, pq = 0; p < 3; p++)
        for (int q = 0; q < 3; q++, pq++) {
            std::stringstream sstream;
            sstream << "mo_tei_deriv2_" << atom1 << atom2 << cartcomp[p] << cartcomp[q];
            SharedMatrix temp = mo_eri_helper(ao_grad[pq], C1, C2, C3, C4);
            temp->set_name(sstream.str());
            mo_grad.push_back(temp);
        }
    return mo_grad;
}

}  // namespace psi
