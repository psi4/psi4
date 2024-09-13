/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#ifdef USING_ecpint
#include "psi4/libmints/ecpint.h"
#endif
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "electricfield.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
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

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_scf.h>
#include <brian_geom_opt.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern brianInt brianRestrictionType;

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

MintsHelper::MintsHelper(std::shared_ptr<BasisSet> basis,
                         std::map<std::string, std::shared_ptr<psi::BasisSet>> basissets, Options &options, int print)
    : options_(options), print_(print) {
    init_helper(basis, basissets);
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
    basissets_ = wavefunction->basissets();
    molecule_ = basisset_->molecule();

    // Make sure molecule is valid.
    molecule_->update_geometry();

    common_init();
}

void MintsHelper::init_helper(std::shared_ptr<BasisSet> basis,
                              std::map<std::string, std::shared_ptr<psi::BasisSet>> basissets) {
    basisset_ = basis;
    basissets_ = basissets;
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

std::shared_ptr<BasisSet> MintsHelper::get_basisset(std::string label) {
    // This may be slightly confusing, but better than changing this in 800 other places
    if (label == "ORBITAL") {
        return basisset_;
    } else if (not basisset_exists(label)) {
        outfile->Printf("Could not find requested basisset (%s).", label.c_str());
        throw PSIEXCEPTION("MintsHelper::get_basisset: Requested basis set (" + label + ") was not set!\n");
    } else {
        return basissets_[label];
    }
}
void MintsHelper::set_basisset(std::string label, std::shared_ptr<BasisSet> basis) {
    if (label == "ORBITAL") {
        throw PSIEXCEPTION("Cannot set the ORBITAL basis after the Wavefunction is built!");
    } else {
        basissets_[label] = basis;
    }
}

bool MintsHelper::basisset_exists(std::string label) {
    if (basissets_.count(label) == 0)
        return false;
    else
        return true;
}

std::shared_ptr<MatrixFactory> MintsHelper::factory() const { return factory_; }

std::shared_ptr<IntegralFactory> MintsHelper::integral() const { return integral_; }

int MintsHelper::nbf() const { return basisset_->nbf(); }

void MintsHelper::integrals() {
    if (print_) {
        outfile->Printf(" MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");
    }

    // Get ERI object
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb(nthread_);
    tb[0] = std::shared_ptr<TwoBodyAOInt>(integral_->eri());
    for (int i = 1; i < nthread_; ++i) tb[i] = std::shared_ptr<TwoBodyAOInt>(tb.front()->clone());
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
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb(nthread_);
    tb[0] = std::shared_ptr<TwoBodyAOInt>(integral_->erf_eri(omega));
    for (int i = 1; i < nthread_; ++i) tb[i] = std::shared_ptr<TwoBodyAOInt>(tb.front()->clone());
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
    std::vector<std::shared_ptr<TwoBodyAOInt>> tb(nthread_);
    tb[0] = std::shared_ptr<TwoBodyAOInt>(integral_->erf_complement_eri(omega));
    for (int i = 1; i < nthread_; ++i)
        tb[i] = std::shared_ptr<TwoBodyAOInt>(tb.front()->clone());
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

    // Limit to the number of incoming onebody ints
    size_t nthread = nthread_;
    if (nthread > ints.size()) {
        nthread = ints.size();
    }

    double **outp = out->pointer();

    const auto &shell_pairs = ints[0]->shellpairs();
    size_t n_pairs = shell_pairs.size();

    // Loop it
#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t p = 0; p < n_pairs; ++p) {
        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto mu = shell_pairs[p].first;
        auto nu = shell_pairs[p].second;
        const size_t num_mu = bs1->shell(mu).nfunction();
        const size_t index_mu = bs1->shell(mu).function_index();
        const size_t num_nu = bs2->shell(nu).nfunction();
        const size_t index_nu = bs2->shell(nu).function_index();

        ints[rank]->compute_shell(mu, nu);
        const auto *ints_buff = ints[rank]->buffers()[0];

        size_t index = 0;
        if (symm) {
            for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                    outp[nu][mu] = outp[mu][nu] = ints_buff[index++];
                }
            }
        } else {
            for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                    outp[mu][nu] = ints_buff[index++];
                }
            }
        }
    }
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

    double **outp = out->pointer();
    double **Dp = D->pointer();

    const auto &shell_pairs = ints[0]->shellpairs();
    size_t n_pairs = shell_pairs.size();

    // Loop it
#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t p = 0; p < n_pairs; ++p) {
        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        ints[rank]->compute_shell_deriv1(P, Q);
        const auto &ints_buff = ints[rank]->buffers();

        size_t nP = basisset_->shell(P).nfunction();
        size_t oP = basisset_->shell(P).function_index();
        size_t aP = basisset_->shell(P).ncenter();

        size_t nQ = basisset_->shell(Q).nfunction();
        size_t oQ = basisset_->shell(Q).function_index();
        size_t aQ = basisset_->shell(Q).ncenter();

        size_t offset = nP * nQ;
        double perm = (P == Q ? 1.0 : 2.0);

        // Px
        double Px = 0.0;
        const double *ref = ints_buff[0];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Px += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aP][0] += Px;

        // Py
        double Py = 0.0;
        ref = ints_buff[1];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Py += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aP][1] += Py;

        // Pz
        double Pz = 0.0;
        ref = ints_buff[2];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Pz += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aP][2] += Pz;

        // Qx
        double Qx = 0.0;
        ref = ints_buff[3];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Qx += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aQ][0] += Qx;

        // Qy
        double Qy = 0.0;
        ref = ints_buff[4];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Qy += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aQ][1] += Qy;

        // Qz
        double Qz = 0.0;
        ref = ints_buff[5];
        for (size_t p = 0; p < nP; p++) {
            for (size_t q = 0; q < nQ; q++) {
                Qz += perm * Dp[p + oP][q + oQ] * (*ref++);
            }
        }
#pragma omp atomic update
        outp[aQ][2] += Qz;
    }
}

SharedMatrix MintsHelper::ao_overlap() {
    // Overlap
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_overlap()));
    }
    auto overlap_mat = std::make_shared<Matrix>(PSIF_AO_S, basisset_->nbf(), basisset_->nbf());

#ifdef USING_BrianQC
    if (brianEnable) {
        brianInt integralType = BRIAN_INTEGRAL_TYPE_OVERLAP;
        brianSCFBuild1e(&brianCookie, &integralType, overlap_mat->get_pointer());
        checkBrian();

        return overlap_mat;
    }
#endif

    one_body_ao_computer(ints_vec, overlap_mat, true);

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
    one_body_ao_computer(ints_vec, overlap_mat, bs1 == bs2);
    return overlap_mat;
}

SharedMatrix MintsHelper::ao_kinetic() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_kinetic()));
    }
    auto kinetic_mat = std::make_shared<Matrix>("AO-basis Kinetic Ints", basisset_->nbf(), basisset_->nbf());

#ifdef USING_BrianQC
    if (brianEnable) {
        brianInt integralType = BRIAN_INTEGRAL_TYPE_KINETIC;
        brianSCFBuild1e(&brianCookie, &integralType, kinetic_mat->get_pointer());
        checkBrian();

        return kinetic_mat;
    }
#endif

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
    one_body_ao_computer(ints_vec, kinetic_mat, bs1 == bs2);
    return kinetic_mat;
}

SharedMatrix MintsHelper::ao_potential() {
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential()));
    }
    SharedMatrix potential_mat =
        std::make_shared<Matrix>("AO-basis Potential Ints", basisset_->nbf(), basisset_->nbf());

#ifdef USING_BrianQC
    if (brianEnable) {
        brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
        brianSCFBuild1e(&brianCookie, &integralType, potential_mat->get_pointer());
        checkBrian();

        return potential_mat;
    }
#endif

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
    one_body_ao_computer(ints_vec, potential_mat, bs1 == bs2);
    return potential_mat;
}

#ifdef USING_ecpint
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
    one_body_ao_computer(ints_vec, ecp_mat, bs1 == bs2);
    return ecp_mat;
}
#endif

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
    MintsHelper decon(get_basisset("BASIS_RELATIVISTIC"));
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

    int nbf = get_basisset("BASIS_RELATIVISTIC")->nbf();
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

    SharedMatrix S_cd = ao_overlap(basisset_, get_basisset("BASIS_RELATIVISTIC"));
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

    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {
                for (int Q = 0; Q < bs4->nshell(); Q++) {
                    ints->compute_shell(M, N, P, Q);
                    const double *buffer = ints->buffer();

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

    ints->compute_shell(M, N, P, Q);
    const double *buffer = ints->buffer();

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

SharedMatrix MintsHelper::ao_f12(std::vector<std::pair<double, double>> exp_coeff) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12(exp_coeff));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12(exp_coeff));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::vector<std::pair<double, double>> exp_coeff) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12_squared(exp_coeff));
    return ao_helper("AO F12 Squared Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::vector<std::pair<double, double>> exp_coeff,
                                         std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                         std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12_squared(exp_coeff));
    return ao_helper("AO F12 Squared Tensor", ints);
}

std::vector<std::pair<double, double>> MintsHelper::f12_cgtg(double exponent) {
    // The fitting coefficients and the exponents
    std::vector<std::pair<double, double>> exp_coeff = {};
    std::vector<double> coeffs = {-0.31442480597241274, -0.30369575353387201, -0.16806968430232927,
                                  -0.098115812152857612, -0.060246640234342785, -0.037263541968504843};
    std::vector<double> exps = {0.22085085450735284, 1.0040191632019282, 3.6212173098378728,
                                12.162483236221904, 45.855332448029337, 254.23460688554644};

    for (int i = 0; i < exps.size(); i++){
        auto exp_scaled = (exponent * exponent) * exps[i];
        exp_coeff.push_back(std::make_pair(exp_scaled, coeffs[i]));
    }
    
    return exp_coeff;
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

    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {
                ints->compute_shell(M, N, P);
                const double *buffer = ints->buffers()[0];
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
    std::shared_ptr<ThreeCenterOverlapInt> ints =
        std::make_shared<ThreeCenterOverlapInt>(basisset_, basisset_, basisset_);
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}

SharedMatrix MintsHelper::ao_3coverlap(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                       std::shared_ptr<BasisSet> bs3) {
    auto ints = std::make_shared<ThreeCenterOverlapInt>(bs1, bs2, bs3);
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12g12(std::vector<std::pair<double, double>> exp_coeff) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12g12(exp_coeff));
    return ao_helper("AO F12G12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12g12(std::vector<std::pair<double, double>> exp_coeff,
                                    std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                    std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12g12(exp_coeff));
    return ao_helper("AO F12G12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff) {
    std::shared_ptr<TwoBodyAOInt> ints(integral_->f12_double_commutator(exp_coeff));
    return ao_helper("AO F12 Double Commutator Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff,
                                         std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                         std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr<TwoBodyAOInt> ints(intf.f12_double_commutator(exp_coeff));
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

SharedMatrix MintsHelper::mo_f12(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1, SharedMatrix C2,
                                 SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12(exp_coeff), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_squared(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1,
                                         SharedMatrix C2, SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_squared(exp_coeff), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Squared Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12g12(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1, SharedMatrix C2,
                                    SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12g12(exp_coeff), C1, C2, C3, C4);
    mo_ints->set_name("MO F12G12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1,
                                                   SharedMatrix C2, SharedMatrix C3, SharedMatrix C4) {
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_double_commutator(exp_coeff), C1, C2, C3, C4);
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

bool MintsHelper::are_ints_cached(const std::string &label, bool include_perturbation) {
    auto it = cached_oe_ints_.find(std::make_pair(label, include_perturbation));
    return it != cached_oe_ints_.end();
}

void MintsHelper::cache_ao_to_so_ints(SharedMatrix ao_ints, const std::string &label, bool include_perturbation) {
    auto p = std::make_pair(label, include_perturbation);
    if (factory_->nirrep() == 1) {
        cached_oe_ints_[p] = ao_ints;
    } else {
        SharedMatrix so_ints_sym(factory_->create_matrix(label));
        so_ints_sym->apply_symmetry(ao_ints, petite_list()->aotoso());
        cached_oe_ints_[p] = so_ints_sym;
    }
}

SharedMatrix MintsHelper::so_overlap_nr() {
    std::string label(PSIF_SO_S);
    SharedMatrix overlap_mat;
    if (factory_->nirrep() == 1) {
        overlap_mat = ao_overlap();
        overlap_mat->set_name(label);
    } else {
        overlap_mat = factory_->create_shared_matrix(label);
        overlap_mat->apply_symmetry(ao_overlap(), petite_list()->aotoso());
    }
    return overlap_mat;
}

SharedMatrix MintsHelper::so_kinetic_nr() {
    std::string label(PSIF_SO_T);
    SharedMatrix kinetic_mat;
    if (factory_->nirrep() == 1) {
        kinetic_mat = ao_kinetic();
        kinetic_mat->set_name(label);
    } else {
        kinetic_mat = factory_->create_shared_matrix(label);
        kinetic_mat->apply_symmetry(ao_kinetic(), petite_list()->aotoso());
    }
    return kinetic_mat;
}

SharedMatrix MintsHelper::so_potential_nr(bool include_perturbations) {
    std::string label(PSIF_SO_V);
    // No symmetry
    SharedMatrix potential_mat;
    if (factory_->nirrep() == 1) {
        potential_mat = ao_potential();
        potential_mat->set_name(label);
    } else {
        potential_mat = factory_->create_shared_matrix(label);
        potential_mat->apply_symmetry(ao_potential(), petite_list()->aotoso());
    }

    // Add ECPs, if needed
    if (basisset_->has_ECP()) {
#ifdef USING_ecpint
        potential_mat->add(so_ecp());
#endif
    }

    // Handle addition of any perturbations here and not in SCF code.
    if (include_perturbations) {
        if (options_.get_bool("PERTURB_H")) {
            add_dipole_perturbation(potential_mat);
        }
    }

    return potential_mat;
}

SharedMatrix MintsHelper::so_overlap(bool include_perturbations) {
    std::string label(PSIF_SO_S);
    auto p = std::make_pair(label, include_perturbations);
    if (!are_ints_cached(label, include_perturbations)) {
        if (options_.get_str("RELATIVISTIC") == "X2C") {
            // generate so_overlap, so_kinetic, so_potential and cache them
            compute_so_x2c_ints(include_perturbations);
        } else {
            cached_oe_ints_[p] = so_overlap_nr();
        }
    }
    return cached_oe_ints_[p];
}

SharedMatrix MintsHelper::so_kinetic(bool include_perturbations) {
    std::string label(PSIF_SO_T);
    auto p = std::make_pair(label, include_perturbations);
    if (!are_ints_cached(label, include_perturbations)) {
        if (options_.get_str("RELATIVISTIC") == "X2C") {
            // generate so_overlap, so_kinetic, so_potential and cache them
            compute_so_x2c_ints(include_perturbations);
        } else {
            cached_oe_ints_[p] = so_kinetic_nr();
        }
    }
    return cached_oe_ints_[p];
}

#ifdef USING_ecpint
SharedMatrix MintsHelper::so_ecp() {
    std::string label(PSIF_SO_ECP);
    if (!basisset_->has_ECP()) {
        SharedMatrix ecp_mat = factory_->create_shared_matrix(label);
        ecp_mat->zero();
        outfile->Printf("\n\tWarning! ECP integrals requested, but no ECP basis detected.  Returning zeros.\n");
        return ecp_mat;
    }
    if (!are_ints_cached(label, false)) {
        cache_ao_to_so_ints(ao_ecp(), label, false);
    }
    auto p = std::make_pair(label, false);
    return cached_oe_ints_[p];
}
#endif

SharedMatrix MintsHelper::so_potential(bool include_perturbations) {
    std::string label(PSIF_SO_V);
    auto p = std::make_pair(label, include_perturbations);
    if (!are_ints_cached(label, include_perturbations)) {
        if (options_.get_str("RELATIVISTIC") == "X2C") {
            // generate so_overlap, so_kinetic, so_potential and cache them
            compute_so_x2c_ints(include_perturbations);
        } else {
            cached_oe_ints_[p] = so_potential_nr(include_perturbations);
            // Add DKH correction if requested
            if (options_.get_str("RELATIVISTIC") == "DKH") {
                if (include_perturbations) {
                    int dkh_order = options_.get_int("DKH_ORDER");
                    SharedMatrix dkh = so_dkh(dkh_order);
                    outfile->Printf("    Adding Douglas-Kroll-Hess corrections to the potential integrals.\n");
                    cached_oe_ints_[p]->add(dkh);
                }
            }
        }
    }
    return cached_oe_ints_[p];
}

void MintsHelper::add_dipole_perturbation(SharedMatrix potential_mat) {
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
    std::unique_ptr<OneBodySOInt> so_dipole = integral_->so_dipole();
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

void MintsHelper::compute_so_x2c_ints(bool include_perturbations) {
    outfile->Printf(" OEINTS: Using relativistic (X2C) overlap, kinetic, and potential integrals.\n");

    if (!basisset_exists("BASIS_RELATIVISTIC")) {
        throw PSIEXCEPTION("OEINTS: X2C requested, but relativistic basis (BASIS_RELATIVISTIC) was not set.");
    }
    SharedMatrix so_overlap_x2c = so_overlap_nr();
    SharedMatrix so_kinetic_x2c = so_kinetic_nr();
    SharedMatrix so_potential_x2c = so_potential_nr(include_perturbations);

    std::vector<double> lambda(3, 0.0);

    if (include_perturbations) {
        if (options_.get_bool("PERTURB_H")) {
            std::string perturb_with = options_.get_str("PERTURB_WITH");
            outfile->Printf("\n  perturb_with = %s", perturb_with.c_str());
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
        }
    }

    X2CInt x2cint;
    x2cint.compute(molecule_, basisset_, get_basisset("BASIS_RELATIVISTIC"), so_overlap_x2c, so_kinetic_x2c,
                   so_potential_x2c, lambda);

    // Overwrite cached integrals
    cached_oe_ints_[std::make_pair(PSIF_SO_S, include_perturbations)] = so_overlap_x2c;
    cached_oe_ints_[std::make_pair(PSIF_SO_T, include_perturbations)] = so_kinetic_x2c;
    cached_oe_ints_[std::make_pair(PSIF_SO_V, include_perturbations)] = so_potential_x2c;
}

std::vector<SharedMatrix> MintsHelper::so_dipole() const {
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

std::vector<SharedMatrix> MintsHelper::so_nabla() const {
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::P, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla = msymm.create_matrices("SO Nabla");

    std::shared_ptr<OneBodySOInt> ints(integral_->so_nabla());
    ints->compute(nabla);

    return nabla;
}

std::vector<SharedMatrix> MintsHelper::so_angular_momentum() const {
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

std::vector<SharedMatrix> MintsHelper::ao_multipoles(int order, const std::vector<double> &origin) {
    if (origin.size() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    Vector3 v3origin(origin[0], origin[1], origin[2]);
    std::vector<SharedMatrix> ret;
    for (int l = 1; l <= order; ++l) {
        for (int ii = 0; ii <= l; ii++) {
            int lx = l - ii;
            for (int lz = 0; lz <= ii; lz++) {
                int ly = ii - lz;
                std::stringstream sstream;
                if (l == 1) {
                    sstream << "Dipole ";
                } else if (l == 2) {
                    sstream << "Quadrupole ";
                } else if (l == 3) {
                    sstream << "Octupole ";
                } else if (l == 4) {
                    sstream << "Hexadecapole ";
                } else {
                    int n = (1 << l);
                    sstream << n << "-pole ";
                }
                std::string name = sstream.str();
                for (int xval = 0; xval < lx; ++xval) name += "X";
                for (int yval = 0; yval < ly; ++yval) name += "Y";
                for (int zval = 0; zval < lz; ++zval) name += "Z";
                auto mat = std::make_shared<Matrix>(name, factory_->norb(), factory_->norb());
                ret.push_back(mat);
            }
        }
    }
    std::shared_ptr<OneBodyAOInt> multipole_int(integral_->ao_multipoles(order));
    multipole_int->set_origin(v3origin);
    multipole_int->compute(ret);
    return ret;
}

std::vector<SharedMatrix> MintsHelper::ao_efp_multipole_potential(const std::vector<double> &origin, int deriv) {
    std::vector<SharedMatrix> ret = ao_multipole_potential(3, origin, deriv);
    // EFP expects the following order of Cartesian components
    //       | // Charge
    //  0    |      0
    //       | // Dipole
    //  1    |      X
    //  2    |      Y
    //  3    |      Z
    //       | // Quadrupole
    //  4    |      XX
    //  5    |      YY
    //  6    |      ZZ
    //  7    |      XY
    //  8    |      XZ
    //  9    |      YZ
    //       | // Octupole
    // 10    |      XXX
    // 11    |      YYY
    // 12    |      ZZZ
    // 13    |      XXY
    // 14    |      XXZ
    // 15    |      XYY
    // 16    |      YYZ
    // 17    |      XZZ
    // 18    |      YZZ
    // 19    |      XYZ
    // Using this mapping, one can convert alphabetical ordering
    // to EFP ordering of the components
    std::vector<int> map_components{
        0,                                      // charge
        1,  2,  3,                              // dipole
        4,  7,  9,  5,  6,  8,                  // quadrupole
        10, 16, 19, 11, 12, 13, 17, 15, 18, 14  // octupole
    };
    std::vector<SharedMatrix> ret_reordered;
    for (size_t i = 0; i < 20; ++i) {
        ret_reordered.push_back(std::move(ret[map_components[i]]));
    }
    return ret_reordered;
}

std::vector<SharedMatrix> MintsHelper::ao_multipole_potential(int order, const std::vector<double> &origin, int deriv) {
    if (origin.size() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    Vector3 v3origin(origin[0], origin[1], origin[2]);

    std::vector<SharedMatrix> ret;
    for (int l = 0; l <= order; ++l) {
        for (int ii = 0; ii <= l; ii++) {
            int lx = l - ii;
            for (int lz = 0; lz <= ii; lz++) {
                int ly = ii - lz;
                std::string name = "AO Multipole Potential ";
                for (int xval = 0; xval < lx; ++xval) name += "X";
                for (int yval = 0; yval < ly; ++yval) name += "Y";
                for (int zval = 0; zval < lz; ++zval) name += "Z";
                if (lx == 0 && ly == 0 && lz == 0) name += "0";
                auto mat = std::make_shared<Matrix>(name, basisset_->nbf(), basisset_->nbf());
                ret.push_back(mat);
            }
        }
    }
    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_multipole_potential(order, deriv));
    ints->set_origin(v3origin);
    ints->compute(ret);
    return ret;
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

SharedMatrix MintsHelper::induction_operator(SharedMatrix coords, SharedMatrix moments) {
    SharedMatrix mat = std::make_shared<Matrix>("Induction operator", basisset_->nbf(), basisset_->nbf());
    ContractOverDipolesFunctor dipfun(moments, mat);
    auto field_integrals_ = std::unique_ptr<ElectricFieldInt>(static_cast<ElectricFieldInt *>(integral_->electric_field().release()));
    field_integrals_->compute_with_functor(dipfun, coords);
    mat->scale(-1.0);

    return mat;
}

SharedMatrix MintsHelper::electric_field_value(SharedMatrix coords, SharedMatrix D) {
    auto field_integrals_ = std::unique_ptr<ElectricFieldInt>(static_cast<ElectricFieldInt *>(integral_->electric_field().release()));

    SharedMatrix efields = std::make_shared<Matrix>("efields", coords->nrow(), 3);
    auto fieldfun = ContractOverDensityFieldFunctor(efields, D);
    field_integrals_->compute_with_functor(fieldfun, coords);

    return efields;
}

SharedMatrix MintsHelper::ao_potential_erf(const std::vector<double> &origin, double omega, int deriv) {
    SharedMatrix int_erf = std::make_shared<Matrix>("AO Potential Erf", basisset_->nbf(), basisset_->nbf());
    Vector3 v3origin(origin[0], origin[1], origin[2]);
    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_potential_erf(omega, deriv));
    ints->set_origin(v3origin);
    ints->compute(int_erf);
    return int_erf;
}

SharedMatrix MintsHelper::ao_potential_erf_complement(const std::vector<double> &origin, double omega, int deriv) {
    SharedMatrix int_erfc = std::make_shared<Matrix>("AO Potential Erf Complement", basisset_->nbf(), basisset_->nbf());
    Vector3 v3origin(origin[0], origin[1], origin[2]);
    std::shared_ptr<OneBodyAOInt> ints(integral_->ao_potential_erf_complement(omega, deriv));
    ints->set_origin(v3origin);
    ints->compute(int_erfc);
    return int_erfc;
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

SharedVector MintsHelper::electrostatic_potential_value(SharedVector charges, SharedMatrix coords, SharedMatrix D) {
    if (coords->ncol() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    if (coords->nrow() != charges->dim()) {
        throw PSIEXCEPTION("Dimension mismatch charges and coordinates.");
    }

    auto potential_integrals_ = std::unique_ptr<PCMPotentialInt>(static_cast<PCMPotentialInt *>(integral_->pcm_potentialint().release()));
    std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
    for (size_t i = 0; i < coords->nrow(); ++i) {
        Zxyz.push_back({charges->pointer()[i], {coords->pointer()[i][0],
                                                coords->pointer()[i][1],
                                                coords->pointer()[i][2]}});
    }
    potential_integrals_->set_charge_field(Zxyz);

    SharedVector potvalues = std::make_shared<Vector>("potential values", coords->nrow());
    ContractOverDensityFunctor contract_density_functor(potvalues->dim(), potvalues->pointer(), D);
    potential_integrals_->compute(contract_density_functor);
    return potvalues;
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
#ifdef USING_BrianQC
    if (brianEnable) {
        int densityCount = (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) ? 1 : 2;

        SharedMatrix dummyInput, dummyOutput;
        if (densityCount > 1) {
            dummyInput = std::make_shared<Matrix>("dummy", basisset_->nbf(), basisset_->nbf());
            dummyOutput = std::make_shared<Matrix>("dummy", basisset_->molecule()->natom(), 3);
        }

        brianInt integralType = BRIAN_INTEGRAL_TYPE_NUCLEAR;
        brianOPTBuildGradient1eDeriv(&brianCookie, &integralType, D->get_pointer(),
                                     (densityCount > 1 ? dummyInput->get_pointer() : nullptr), nullptr, nullptr,
                                     V->get_pointer(), (densityCount > 1 ? dummyOutput->get_pointer() : nullptr));

        return V;
    }
#endif

    // Build temps
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    std::vector<SharedMatrix> Vtemps;
    for (size_t i = 0; i < nthread_; i++) {
        Vtemps.push_back(V->clone());
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
    }

    const auto &shell_pairs = ints_vec[0]->shellpairs();

    double **Dp = D->pointer();

#pragma omp parallel for schedule(dynamic) num_threads(nthread_)
    for (size_t PQ = 0L; PQ < shell_pairs.size(); PQ++) {
        size_t P = shell_pairs[PQ].first;
        size_t Q = shell_pairs[PQ].second;

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        ints_vec[rank]->compute_shell_deriv1(P, Q);
        const auto &buffers = ints_vec[rank]->buffers();

        size_t nP = basisset_->shell(P).nfunction();
        size_t oP = basisset_->shell(P).function_index();
        size_t aP = basisset_->shell(P).ncenter();

        size_t nQ = basisset_->shell(Q).nfunction();
        size_t oQ = basisset_->shell(Q).function_index();
        size_t aQ = basisset_->shell(Q).ncenter();

        double perm = (P == Q ? 1.0 : 2.0);

        double **Vp = Vtemps[rank]->pointer();

        // The three buffers have Px, Py, Pz, then the next three are Qx, Qy, Qz, with the next
        // 3*natom containing derivatives with respect to each atom's nuclear charge position.
        for (size_t X = 0; X < 2 + natom; X++) {
            auto A = (X == 0) ? aP : ((X == 1) ? aQ : X - 2);

            const double *ref0 = buffers[3 * X + 0];
            const double *ref1 = buffers[3 * X + 1];
            const double *ref2 = buffers[3 * X + 2];
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
    // Kinetic
    SharedMatrix kinetic_mat(new Matrix("Kinetic Gradient", basisset_->molecule()->natom(), 3));
#ifdef USING_BrianQC
    if (brianEnable) {
        int densityCount = (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) ? 1 : 2;

        SharedMatrix dummyInput, dummyOutput;
        if (densityCount > 1) {
            dummyInput = std::make_shared<Matrix>("dummy", basisset_->nbf(), basisset_->nbf());
            dummyOutput = std::make_shared<Matrix>("dummy", basisset_->molecule()->natom(), 3);
        }

        brianInt integralType = BRIAN_INTEGRAL_TYPE_KINETIC;
        brianOPTBuildGradient1eDeriv(
            &brianCookie, &integralType, D->get_pointer(), (densityCount > 1 ? dummyInput->get_pointer() : nullptr),
            nullptr, nullptr, kinetic_mat->get_pointer(), (densityCount > 1 ? dummyOutput->get_pointer() : nullptr));

        return kinetic_mat;
    }
#endif
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_kinetic(1)));
    }
    grad_two_center_computer(ints_vec, D, kinetic_mat);
    return kinetic_mat;
}
SharedMatrix MintsHelper::overlap_grad(SharedMatrix D) {
    // Overlap
    SharedMatrix overlap_mat(new Matrix("Overlap Gradient", basisset_->molecule()->natom(), 3));
#ifdef USING_BrianQC
    if (brianEnable) {
        int densityCount = (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) ? 1 : 2;

        SharedMatrix dummyInput, dummyOutput;
        if (densityCount > 1) {
            dummyInput = std::make_shared<Matrix>("dummy", basisset_->nbf(), basisset_->nbf());
            dummyOutput = std::make_shared<Matrix>("dummy", basisset_->molecule()->natom(), 3);
        }

        brianInt integralType = BRIAN_INTEGRAL_TYPE_OVERLAP;
        brianOPTBuildGradient1eDeriv(&brianCookie, &integralType, nullptr, nullptr, D->get_pointer(),
                                     (densityCount > 1 ? dummyInput->get_pointer() : nullptr),
                                     overlap_mat->get_pointer(),
                                     (densityCount > 1 ? dummyOutput->get_pointer() : nullptr));

        return overlap_mat;
    }
#endif
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++) {
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_overlap(1)));
    }
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

#ifdef USING_ecpint
SharedMatrix MintsHelper::effective_core_potential_grad(SharedMatrix D) {
    int natom = basisset_->molecule()->natom();
    auto grad = std::make_shared<Matrix>("Effective Core Potential Gradient", natom, 3);

    // Build temps
    std::vector<std::shared_ptr<OneBodyAOInt>> ecp_ints_vec;
    std::vector<SharedMatrix> gradtemps;
    for (size_t i = 0; i < nthread_; i++) {
        gradtemps.push_back(grad->clone());
        ecp_ints_vec.push_back(std::shared_ptr<ECPInt>(dynamic_cast<ECPInt*>(integral_->ao_ecp(1).release())));
    }

    // Lower Triangle
    std::vector<std::pair<int, int>> PQ_pairs;
    for (int P = 0; P < basisset_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

    // Make a list of all ECP centers
    std::set<int> ecp_centers;
    for (int ecp_shell = 0; ecp_shell < basisset_->n_ecp_shell(); ++ecp_shell){
        const GaussianShell &ecp = basisset_->ecp_shell(ecp_shell);
        ecp_centers.insert(ecp.ncenter());
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

        ecp_ints_vec[rank]->compute_shell_deriv1(P, Q);
        const auto &buffers = ecp_ints_vec[rank]->buffers();

        size_t nP = basisset_->shell(P).nfunction();
        size_t oP = basisset_->shell(P).function_index();
        size_t aP = basisset_->shell(P).ncenter();

        size_t nQ = basisset_->shell(Q).nfunction();
        size_t oQ = basisset_->shell(Q).function_index();
        size_t aQ = basisset_->shell(Q).ncenter();

        // Make a list of all ECP centers and the current basis function pair
        std::set<int> all_centers(ecp_centers.begin(), ecp_centers.end());
        all_centers.insert(aP);
        all_centers.insert(aQ);

        double perm = (P == Q ? 1.0 : 2.0);

        double **Vp = gradtemps[rank]->pointer();

        size_t size = nP * nQ;
        for (const int center : all_centers) {
            const double *ref0 = buffers[3 * center + 0];
            const double *ref1 = buffers[3 * center + 1];
            const double *ref2 = buffers[3 * center + 2];
            for (size_t p = 0; p < nP; p++) {
                for (size_t q = 0; q < nQ; q++) {
                    double Vval = perm * Dp[p + oP][q + oQ];
                    Vp[center][0] += Vval * (*ref0++);
                    Vp[center][1] += Vval * (*ref1++);
                    Vp[center][2] += Vval * (*ref2++);
                }
            }
        }
    }

    // Sum it up
    for (size_t t = 0; t < nthread_; t++) {
        grad->axpy(1.0, gradtemps[t]);
    }
    return grad;
}
#endif

SharedMatrix MintsHelper::multipole_grad(SharedMatrix D, int order, const std::vector<double> &origin) {
    if (origin.size() != 3) throw PSIEXCEPTION("Origin argument must have length 3.");
    // Computes skeleton (Hellman-Feynman like) multipole derivatives for each perturbation
    double **Dp = D->pointer();

    int natom = molecule_->natom();
    int nmult = (order + 1) * (order + 2) * (order + 3) / 6 - 1;
    auto ret = std::make_shared<Matrix>("Multipole dervatives (pert*component, i.e. 3NxN_mult)", 3 * natom, nmult);
    double **Pp = ret->pointer();

    std::shared_ptr<OneBodyAOInt> Mint(integral_->ao_multipoles(order, 1));
    Vector3 v3origin(origin[0], origin[1], origin[2]);
    Mint->set_origin(v3origin);

    const auto &shell_pairs = Mint->shellpairs();
    size_t n_pairs = shell_pairs.size();

    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;

        Mint->compute_shell_deriv1(P, Q);
        const auto &buffers = Mint->buffers();

        const auto &shellP = basisset_->shell(P);
        const auto &shellQ = basisset_->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        double prefac = (P == Q ? 1.0 : 2.0);

        // loop over the multipole components (x, y, z, xx, xy, xz, ...)
        for (int chunk = 0; chunk < nmult; ++chunk) {
            // loop over the derivative components (x, y, z)
            for (int comp = 0; comp < 3; ++comp) {
                // ordering in buffers is (xPx, xPy, xPz, xQx, xQy, xQz, yPx, yPy, ..., xxPx, xxPy, ...)
                // bra derivatives on atom aP
                const double *ref_bra = buffers[6 * chunk + comp];
                // ket derivatives on atom aQ (3 elements offset)
                const double *ref_ket = buffers[6 * chunk + comp + 3];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Pp[3 * aP + comp][chunk] += prefac * Dp[p + oP][q + oQ] * (*ref_bra++);
                        Pp[3 * aQ + comp][chunk] += prefac * Dp[p + oP][q + oQ] * (*ref_ket++);
                    }
                }
            }
        }
    }
    return ret;
}

SharedMatrix MintsHelper::dipole_grad(SharedMatrix D) {
    // Computes skeleton (Hellman-Feynman like) dipole derivatives for each perturbation
    // call the more general routine for arbitrary order multipoles
    return multipole_grad(D, 1, {0.0, 0.0, 0.0});
}

SharedMatrix MintsHelper::core_hamiltonian_grad(SharedMatrix D) {
    auto ret = kinetic_grad(D);
    ret->set_name("Core Hamiltonian Gradient");
    ret->add(potential_grad(D));

    if (options_.get_bool("PERTURB_H")) {
        ret->add(perturb_grad(D));
    }
    if (basisset_->n_ecp_shell()) {
#ifdef USING_ecpint
        ret->add(effective_core_potential_grad(D));
#endif
    }
    return ret;
}

std::map<std::string, SharedMatrix> MintsHelper::metric_grad(std::map<std::string, SharedMatrix> &D,
                                                             const std::string &aux_name) {
    // Construct integral factory.
    auto auxiliary = get_basisset(aux_name);
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary, BasisSet::zero_ao_basis_set(), auxiliary,
                                                       BasisSet::zero_ao_basis_set());
    std::vector<std::shared_ptr<TwoBodyAOInt>> Jint(nthread_);
    Jint[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1));
    for (int t = 1; t < nthread_; t++) {
        Jint[t] = std::shared_ptr<TwoBodyAOInt>(Jint.front()->clone());
    }

    // Construct temporary matrices for each thread
    int natom = basisset_->molecule()->natom();
    std::map<std::string, std::vector<SharedMatrix>> temps;
    for (auto kv : D) {
        temps[kv.first] = std::vector<SharedMatrix>();
        auto &temp = temps[kv.first];
        for (int j = 0; j < nthread_; j++) {
            temp.push_back(std::make_shared<Matrix>("temp", natom, 3));
        }
    }

    // Construct pairs of aux AO shells
    std::vector<std::pair<int, int>> PQ_pairs;
    for (int P = 0; P < auxiliary->nshell(); P++) {
        for (auto Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

    // Perform threading contraction of "densities" against metric derivative integrals.
#pragma omp parallel for schedule(dynamic) num_threads(nthread_)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
        auto P = PQ_pairs[PQ].first;
        auto Q = PQ_pairs[PQ].second;
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        Jint[thread]->compute_shell_deriv1(P, 0, Q, 0);
        const auto buffers = Jint[thread]->buffers();

        int nP = auxiliary->shell(P).nfunction();
        int cP = auxiliary->shell(P).ncartesian();
        int aP = auxiliary->shell(P).ncenter();
        int oP = auxiliary->shell(P).function_index();

        int nQ = auxiliary->shell(Q).nfunction();
        int cQ = auxiliary->shell(Q).ncartesian();
        int aQ = auxiliary->shell(Q).ncenter();
        int oQ = auxiliary->shell(Q).function_index();

        int ncart = cP * cQ;
        const double *Px = buffers[0];
        const double *Py = buffers[1];
        const double *Pz = buffers[2];
        const double *Qx = buffers[3];
        const double *Qy = buffers[4];
        const double *Qz = buffers[5];

        double perm = (P == Q ? 1.0 : 2.0);

        std::vector<std::pair<double **, double **>> read_write_pairs;
        for (auto kv : D) {
            auto temp1 = D[kv.first]->pointer();
            auto temp2 = temps[kv.first][thread]->pointer();
            auto pair = std::make_pair(temp1, temp2);
            read_write_pairs.push_back(pair);
        }

        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) {
                for (auto &pair : read_write_pairs) {
                    double val = 0.5 * perm * pair.first[p + oP][q + oQ];
                    auto grad_mat = pair.second;
                    grad_mat[aP][0] -= val * (*Px);
                    grad_mat[aP][1] -= val * (*Py);
                    grad_mat[aP][2] -= val * (*Pz);
                    grad_mat[aQ][0] -= val * (*Qx);
                    grad_mat[aQ][1] -= val * (*Qy);
                    grad_mat[aQ][2] -= val * (*Qz);
                }

                Px++;
                Py++;
                Pz++;
                Qx++;
                Qy++;
                Qz++;
            }
        }
    }

    // Sum results across the various threads
    std::map<std::string, SharedMatrix> gradient_contributions;
    for (auto kv : temps) {
        auto gradient_contribution = std::make_shared<Matrix>(kv.first + " Gradient", natom, 3);
        for (auto thread_matrix : kv.second) {
            gradient_contribution->add(thread_matrix);
        }
        gradient_contributions[kv.first] = gradient_contribution;
    }
    return gradient_contributions;
}

// TODO: DFMP2 might be able to use this, if we resolve the following:
//  1. Should we move density back transform into this loop?
//  2. Should we explicitly hermitivitize during contraction against derivative integral?
//  3. Can we force a relation between intermed_name and gradient_name, to simplify the argument list?
SharedMatrix MintsHelper::three_idx_grad(const std::string &aux_name, const std::string &intermed_name,
                                         const std::string &gradient_name) {
    // Construct integral factory.
    auto primary = get_basisset("ORBITAL");
    auto auxiliary = get_basisset(aux_name);
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary, BasisSet::zero_ao_basis_set(), primary, primary);
    std::vector<std::shared_ptr<TwoBodyAOInt>> Jint(nthread_);
    Jint[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1));
    for (int t = 1; t < nthread_; t++) {
        Jint[t] = std::shared_ptr<TwoBodyAOInt>(Jint.front()->clone());
    }

    const auto &shell_pairs = Jint[0]->shell_pairs();
    int npairs = shell_pairs.size();  // Number of pairs of primary orbital shells

    // => Memory Constraints <= //
    const auto nprim = primary->nbf();
    const auto naux = auxiliary->nbf();
    const auto ntri = (nprim * (nprim + 1)) / 2;
    auto row_cost = sizeof(double) * (nprim * nprim + ntri);
    // Assume we can devote 80% of Psi's memory to this. 80% was pulled from a hat.
    auto max_rows = 0.8 * Process::environment.get_memory() / row_cost;

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary->nshell(); P++) {
        int nP = auxiliary->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary->nshell());

    // Construct temporary matrices for each thread
    int natom = basisset_->molecule()->natom();
    std::vector<SharedMatrix> temps;
    for (int j = 0; j < nthread_; j++) {
        temps.push_back(std::make_shared<Matrix>("temp", natom, 3));
    }

    psio_address next_Pmn = PSIO_ZERO;
    // Individual block reads may very well not use all of this memory.
    auto temp = std::vector<double>(static_cast<size_t>(max_rows) * ntri);
    auto data = temp.data();

    // Perform threaded contraction of "densities" against 3-index derivative integrals.
    // Loop over blocks. Each block is all (P|mn) belonging to certain aux. orbital shells.
    // Large aux. basis sets may require multiple blocks.
    for (int block = 0; block < Pstarts.size() - 1; block++) {
        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary->nshell() ? naux : auxiliary->shell(Pstop).function_index());
        int np = pstop - pstart;

        // Read values from disk. We assume that only the "lower triangle" of (P|mn) is stored.
        psio_->read(PSIF_AO_TPDM, intermed_name.c_str(), (char *)temp.data(), sizeof(double) * np * ntri, next_Pmn,
                    &next_Pmn);

        // Now get them into a matrix, not just the lower triangle.
        // This is the price of only storing the lower triangle on disk.
        auto idx3_matrix = std::make_shared<Matrix>(np, nprim * nprim);
        auto idx3p = idx3_matrix->pointer();
#pragma omp parallel for
        for (int aux = 0; aux < np; aux++) {
            auto elt = &data[ntri * aux];
            for (int p = 0; p < nprim; p++) {
                for (int q = 0; q <= p; q++) {
                    idx3p[aux][p * nprim + q] = *elt;
                    idx3p[aux][q * nprim + p] = *elt;
                    elt++;
                }
            }
        }

        // For each block, loop over aux. shell, then primary shell pairs
#pragma omp parallel for schedule(dynamic) num_threads(nthread_)
        for (long int PMN = 0L; PMN < static_cast<long int>(NP) * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            Jint[thread]->compute_shell_deriv1(P, 0, M, N);

            const auto &buffers = Jint[thread]->buffers();

            int nP = auxiliary->shell(P).nfunction();
            int cP = auxiliary->shell(P).ncartesian();
            int aP = auxiliary->shell(P).ncenter();
            int oP = auxiliary->shell(P).function_index() - pstart;

            int nM = primary->shell(M).nfunction();
            int cM = primary->shell(M).ncartesian();
            int aM = primary->shell(M).ncenter();
            int oM = primary->shell(M).function_index();

            int nN = primary->shell(N).nfunction();
            int cN = primary->shell(N).ncartesian();
            int aN = primary->shell(N).ncenter();
            int oN = primary->shell(N).function_index();

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

            auto grad_Jp = temps[thread]->pointer();

            // Within each of those, then loop over each function in the shell.
            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        double Ival = 1.0 * perm * idx3p[p + oP][(m + oM) * nprim + (n + oN)];
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

    // Sum results across the various threads
    auto idx3_grad = std::make_shared<Matrix>(intermed_name + " Gradient", natom, 3);
    for (const auto &thread_contribution : temps) {
        idx3_grad->add(thread_contribution);
    }

    return idx3_grad;
}

void MintsHelper::play() {}

/* 1st and 2nd derivatives of OEI in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_overlap_kinetic_deriv1_helper(const std::string &type, int atom) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
        grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
    }

    const auto &shell_pairs = GInt->shellpairs();
    size_t n_pairs = shell_pairs.size();

    // Loop it
    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        const auto &shellP = basisset_->shell(P);
        const auto &shellQ = basisset_->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        if (aP != atom && aQ != atom) continue;

        GInt->compute_shell_deriv1(P, Q);
        const auto &buffers = GInt->buffers();
        double scale = P == Q ? 0.5 : 1.0;

        if (aP == atom) {
            // Px
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, scale * buffers[0][p * nQ + q]);
                    grad[0]->add(q + oQ, p + oP, scale * buffers[0][p * nQ + q]);
                }
            }

            // Py
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(p + oP, q + oQ, scale * buffers[1][p * nQ + q]);
                    grad[1]->add(q + oQ, p + oP, scale * buffers[1][p * nQ + q]);
                }
            }

            // Pz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(p + oP, q + oQ, scale * buffers[2][p * nQ + q]);
                    grad[2]->add(q + oQ, p + oP, scale * buffers[2][p * nQ + q]);
                }
            }
        }

        if (aQ == atom) {
            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, scale * buffers[3][p * nQ + q]);
                    grad[0]->add(q + oQ, p + oP, scale * buffers[3][p * nQ + q]);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(p + oP, q + oQ, scale * buffers[4][p * nQ + q]);
                    grad[1]->add(q + oQ, p + oP, scale * buffers[4][p * nQ + q]);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(p + oP, q + oQ, scale * buffers[5][p * nQ + q]);
                    grad[2]->add(q + oQ, p + oP, scale * buffers[5][p * nQ + q]);
                }
            }
        }
    }

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_potential_deriv1_helper(int atom) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
        grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
    }

    const auto &shell_pairs = Vint->shellpairs();
    size_t n_pairs = shell_pairs.size();

    // Loop it
    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        const auto &shellP = bs1->shell(P);
        const auto &shellQ = bs2->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        Vint->compute_shell_deriv1(P, Q);
        const auto &buffers = Vint->buffers();

        double scale = P == Q ? 0.5 : 1.0;

        const double *ref0 = buffers[3 * atom + 6];
        const double *ref1 = buffers[3 * atom + 7];
        const double *ref2 = buffers[3 * atom + 8];
        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) {
                grad[0]->add(p + oP, q + oQ, scale * (*ref0));
                grad[1]->add(p + oP, q + oQ, scale * (*ref1));
                grad[2]->add(p + oP, q + oQ, scale * (*ref2));
                grad[0]->add(q + oQ, p + oP, scale * (*ref0++));
                grad[1]->add(q + oQ, p + oP, scale * (*ref1++));
                grad[2]->add(q + oQ, p + oP, scale * (*ref2++));
            }
        }

        if (aP == atom) {
            ref0 = buffers[0];
            ref1 = buffers[1];
            ref2 = buffers[2];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, scale * (*ref0));
                    grad[1]->add(p + oP, q + oQ, scale * (*ref1));
                    grad[2]->add(p + oP, q + oQ, scale * (*ref2));
                    grad[0]->add(q + oQ, p + oP, scale * (*ref0++));
                    grad[1]->add(q + oQ, p + oP, scale * (*ref1++));
                    grad[2]->add(q + oQ, p + oP, scale * (*ref2++));
                }
            }
        }

        if (aQ == atom) {
            ref0 = buffers[3];
            ref1 = buffers[4];
            ref2 = buffers[5];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, scale * (*ref0));
                    grad[1]->add(p + oP, q + oQ, scale * (*ref1));
                    grad[2]->add(p + oP, q + oQ, scale * (*ref2));
                    grad[0]->add(q + oQ, p + oP, scale * (*ref0++));
                    grad[1]->add(q + oQ, p + oP, scale * (*ref1++));
                    grad[2]->add(q + oQ, p + oP, scale * (*ref2++));
                }
            }
        }
    }

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_overlap_half_deriv1_helper(const std::string &half_der_side, int atom) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    std::shared_ptr<OneBodyAOInt> GInt(integral_->ao_overlap(1));

    std::shared_ptr<BasisSet> bs1 = GInt->basis1();
    std::shared_ptr<BasisSet> bs2 = GInt->basis2();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_overlap_half_deriv1_" << atom << cartcomp[p];
        grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
    }

    const auto &shell_pairs = GInt->shellpairs();
    size_t n_pairs = shell_pairs.size();

    // Loop it
    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        const auto &shellP = basisset_->shell(P);
        const auto &shellQ = basisset_->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        if (aP != atom && aQ != atom) continue;

        GInt->compute_shell_deriv1(P, Q);
        const auto &buffers = GInt->buffers();
        int offset = 0;
        if (aP == atom && half_der_side == "LEFT") {
            // Px
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, buffers[0][p * nQ + q]);
                }
            }

            // Py
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(p + oP, q + oQ, buffers[1][p * nQ + q]);
                }
            }

            // Pz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(p + oP, q + oQ, buffers[2][p * nQ + q]);
                }
            }
        }

        if (aQ == atom && half_der_side == "LEFT") {
            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(q + oQ, p + oP, buffers[3][p * nQ + q]);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(q + oQ, p + oP, buffers[4][p * nQ + q]);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(q + oQ, p + oP, buffers[5][p * nQ + q]);
                }
            }
        }

        if (aP == atom && half_der_side == "RIGHT") {
            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(q + oQ, p + oP, buffers[0][p * nQ + q]);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(q + oQ, p + oP, buffers[1][p * nQ + q]);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(q + oQ, p + oP, buffers[2][p * nQ + q]);
                }
            }
        }

        if (aQ == atom && half_der_side == "RIGHT") {
            // Qx
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[0]->add(p + oP, q + oQ, buffers[3][p * nQ + q]);
                }
            }

            // Qy
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[1]->add(p + oP, q + oQ, buffers[4][p * nQ + q]);
                }
            }

            // Qz
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    grad[2]->add(p + oP, q + oQ, buffers[5][p * nQ + q]);
                }
            }
        }
    }

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_potential_deriv2_helper(int atom1, int atom2) {
    /* NOTE: the x, y, and z in this vector must remain lowercase for this function */
    std::array<std::string, 3> cartcomp{{"x", "y", "z"}};

    std::shared_ptr<OneBodyAOInt> Vint(integral_->ao_potential(2));

    std::shared_ptr<BasisSet> bs1 = Vint->basis1();
    std::shared_ptr<BasisSet> bs2 = Vint->basis2();

    // Sets up the field of partial charges
    std::vector<std::pair<double, std::array<double, 3>>> full_params;
    const auto &mol = bs1->molecule();
    for (int A = 0; A < mol->natom(); A++) {
        full_params.push_back({(double)mol->Z(A), {mol->x(A), mol->y(A), mol->z(A)}});
    }

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int natom = molecule_->natom();

    std::vector<SharedMatrix> hess;
    for (int a = 0, ab = 0; a < 3; a++)
        for (int b = 0; b < 3; b++, ab++) {
            std::stringstream sstream;
            sstream << "ao_potential_deriv2_" << atom1 << atom2 << cartcomp[a] << cartcomp[b];
            hess.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
        }

    const double scale = (atom1 == atom2 ? 2.0 : 1.0);

    auto upper_triangle_index = [](int matrix_dim, long i, long j) {
        return std::min(i, j) * (2 * matrix_dim - std::min(i, j) - 1) / 2 + std::max(i, j);
    };

    const auto &shell_pairs = Vint->shellpairs();
    size_t n_pairs = shell_pairs.size();

    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        const auto &shellP = bs1->shell(P);
        const auto &shellQ = bs2->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        bool do_full_field;
        int num_buffers;
        if (atom1 != aP && atom2 != aQ) {
            // Neither bra nor ket are atoms of interest - only compute external charge contributions for atom1 and
            // atom2
            do_full_field = false;
            num_buffers = 12;
            std::dynamic_pointer_cast<PotentialInt>(Vint)->set_charge_field(
                {{(double)mol->Z(atom1), {mol->x(atom1), mol->y(atom1), mol->z(atom1)}},
                 {(double)mol->Z(atom2), {mol->x(atom2), mol->y(atom2), mol->z(atom2)}}});
        } else {
            // One of the atoms of interest is in the bra or ket - do a full computation
            do_full_field = true;
            num_buffers = 3 * (natom + 2);
            std::dynamic_pointer_cast<PotentialInt>(Vint)->set_charge_field(full_params);
        }

        Vint->compute_shell_deriv2(P, Q);
        const auto &buffers = Vint->buffers();

        double perm = P == Q ? 0.5 : 1.0;
        // clang-format off
        // The derivatives emerge in upper triangular order 
        // AxAx, AxAy, AxAz, AxBx, AxBy, AxBz, AxC1x, AxC2x, ... AxCNz
        //       AyAy, AyAz, AyBx, AyBy, AyBz, AyC1x, AyC2x, ... AyCNz
        //             AzAz, AzBx, AzBy, AyBz, AzC1x, AzC2x, ... AzCNz
        //                                                         .
        //                                                         .
        //                                                         .
        //                                               CNyCNz CNzCNz
        //                                                      CNzCNz
        // clang-format on
        // For each shell pair, we get 3 buffers in the bra, 3 in the ket and 3 for each atom

        if (aP == atom1 && aP == atom2) {
            const double *AxAx = buffers[upper_triangle_index(num_buffers, 0, 0)];
            const double *AxAy = buffers[upper_triangle_index(num_buffers, 0, 1)];
            const double *AxAz = buffers[upper_triangle_index(num_buffers, 0, 2)];
            const double *AyAy = buffers[upper_triangle_index(num_buffers, 1, 1)];
            const double *AyAz = buffers[upper_triangle_index(num_buffers, 1, 2)];
            const double *AzAz = buffers[upper_triangle_index(num_buffers, 2, 2)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * (*AxAx));
                    hess[1]->add(p + oP, q + oQ, perm * (*AxAy));
                    hess[2]->add(p + oP, q + oQ, perm * (*AxAz));
                    hess[4]->add(p + oP, q + oQ, perm * (*AyAy));
                    hess[5]->add(p + oP, q + oQ, perm * (*AyAz));
                    hess[8]->add(p + oP, q + oQ, perm * (*AzAz));
                    hess[0]->add(q + oQ, p + oP, perm * (*AxAx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*AxAy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*AxAz++));
                    hess[4]->add(q + oQ, p + oP, perm * (*AyAy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*AyAz++));
                    hess[8]->add(q + oQ, p + oP, perm * (*AzAz++));
                }
            }
        }

        if (aQ == atom1 && aQ == atom2) {
            const double *BxBx = buffers[upper_triangle_index(num_buffers, 3, 3)];
            const double *BxBy = buffers[upper_triangle_index(num_buffers, 3, 4)];
            const double *BxBz = buffers[upper_triangle_index(num_buffers, 3, 5)];
            const double *ByBy = buffers[upper_triangle_index(num_buffers, 4, 4)];
            const double *ByBz = buffers[upper_triangle_index(num_buffers, 4, 5)];
            const double *BzBz = buffers[upper_triangle_index(num_buffers, 5, 5)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * (*BxBx));
                    hess[1]->add(p + oP, q + oQ, perm * (*BxBy));
                    hess[2]->add(p + oP, q + oQ, perm * (*BxBz));
                    hess[4]->add(p + oP, q + oQ, perm * (*ByBy));
                    hess[5]->add(p + oP, q + oQ, perm * (*ByBz));
                    hess[8]->add(p + oP, q + oQ, perm * (*BzBz));
                    hess[0]->add(q + oQ, p + oP, perm * (*BxBx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*BxBy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*BxBz++));
                    hess[4]->add(q + oQ, p + oP, perm * (*ByBy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*ByBz++));
                    hess[8]->add(q + oQ, p + oP, perm * (*BzBz++));
                }
            }
        }

        if (aP == atom1 && aQ == atom2) {
            const double *AxBx = buffers[upper_triangle_index(num_buffers, 0, 3)];
            const double *AyBx = buffers[upper_triangle_index(num_buffers, 1, 3)];
            const double *AzBx = buffers[upper_triangle_index(num_buffers, 2, 3)];
            const double *AxBy = buffers[upper_triangle_index(num_buffers, 0, 4)];
            const double *AyBy = buffers[upper_triangle_index(num_buffers, 1, 4)];
            const double *AzBy = buffers[upper_triangle_index(num_buffers, 2, 4)];
            const double *AxBz = buffers[upper_triangle_index(num_buffers, 0, 5)];
            const double *AyBz = buffers[upper_triangle_index(num_buffers, 1, 5)];
            const double *AzBz = buffers[upper_triangle_index(num_buffers, 2, 5)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * scale * (*AxBx));
                    hess[1]->add(p + oP, q + oQ, perm * (*AxBy));
                    hess[2]->add(p + oP, q + oQ, perm * (*AxBz));
                    hess[3]->add(p + oP, q + oQ, perm * (*AyBx));
                    hess[4]->add(p + oP, q + oQ, perm * scale * (*AyBy));
                    hess[5]->add(p + oP, q + oQ, perm * (*AyBz));
                    hess[6]->add(p + oP, q + oQ, perm * (*AzBx));
                    hess[7]->add(p + oP, q + oQ, perm * (*AzBy));
                    hess[8]->add(p + oP, q + oQ, perm * scale * (*AzBz));
                    hess[0]->add(q + oQ, p + oP, perm * scale * (*AxBx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*AxBy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*AxBz++));
                    hess[3]->add(q + oQ, p + oP, perm * (*AyBx++));
                    hess[4]->add(q + oQ, p + oP, perm * scale * (*AyBy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*AyBz++));
                    hess[6]->add(q + oQ, p + oP, perm * (*AzBx++));
                    hess[7]->add(q + oQ, p + oP, perm * (*AzBy++));
                    hess[8]->add(q + oQ, p + oP, perm * scale * (*AzBz++));
                }
            }
        }

        if (atom1 == atom2) {
            int C = do_full_field ? 3 * atom1 + 6 : 6;
            const double *CxCx = buffers[upper_triangle_index(num_buffers, C + 0, C + 0)];
            const double *CxCy = buffers[upper_triangle_index(num_buffers, C + 0, C + 1)];
            const double *CxCz = buffers[upper_triangle_index(num_buffers, C + 0, C + 2)];
            const double *CyCy = buffers[upper_triangle_index(num_buffers, C + 1, C + 1)];
            const double *CyCz = buffers[upper_triangle_index(num_buffers, C + 1, C + 2)];
            const double *CzCz = buffers[upper_triangle_index(num_buffers, C + 2, C + 2)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * (*CxCx));
                    hess[1]->add(p + oP, q + oQ, perm * (*CxCy));
                    hess[2]->add(p + oP, q + oQ, perm * (*CxCz));
                    hess[4]->add(p + oP, q + oQ, perm * (*CyCy));
                    hess[5]->add(p + oP, q + oQ, perm * (*CyCz));
                    hess[8]->add(p + oP, q + oQ, perm * (*CzCz));
                    hess[0]->add(q + oQ, p + oP, perm * (*CxCx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*CxCy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*CxCz++));
                    hess[4]->add(q + oQ, p + oP, perm * (*CyCy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*CyCz++));
                    hess[8]->add(q + oQ, p + oP, perm * (*CzCz++));
                }
            }
        }

        if (aP == atom1) {
            int C = do_full_field ? 3 * atom2 + 6 : 9;
            const double *AxCx = buffers[upper_triangle_index(num_buffers, 0, C + 0)];
            const double *AyCx = buffers[upper_triangle_index(num_buffers, 1, C + 0)];
            const double *AzCx = buffers[upper_triangle_index(num_buffers, 2, C + 0)];
            const double *AxCy = buffers[upper_triangle_index(num_buffers, 0, C + 1)];
            const double *AyCy = buffers[upper_triangle_index(num_buffers, 1, C + 1)];
            const double *AzCy = buffers[upper_triangle_index(num_buffers, 2, C + 1)];
            const double *AxCz = buffers[upper_triangle_index(num_buffers, 0, C + 2)];
            const double *AyCz = buffers[upper_triangle_index(num_buffers, 1, C + 2)];
            const double *AzCz = buffers[upper_triangle_index(num_buffers, 2, C + 2)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * scale * (*AxCx));
                    hess[1]->add(p + oP, q + oQ, perm * (*AxCy));
                    hess[2]->add(p + oP, q + oQ, perm * (*AxCz));
                    hess[3]->add(p + oP, q + oQ, perm * (*AyCx));
                    hess[4]->add(p + oP, q + oQ, perm * scale * (*AyCy));
                    hess[5]->add(p + oP, q + oQ, perm * (*AyCz));
                    hess[6]->add(p + oP, q + oQ, perm * (*AzCx));
                    hess[7]->add(p + oP, q + oQ, perm * (*AzCy));
                    hess[8]->add(p + oP, q + oQ, perm * scale * (*AzCz));
                    hess[0]->add(q + oQ, p + oP, perm * scale * (*AxCx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*AxCy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*AxCz++));
                    hess[3]->add(q + oQ, p + oP, perm * (*AyCx++));
                    hess[4]->add(q + oQ, p + oP, perm * scale * (*AyCy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*AyCz++));
                    hess[6]->add(q + oQ, p + oP, perm * (*AzCx++));
                    hess[7]->add(q + oQ, p + oP, perm * (*AzCy++));
                    hess[8]->add(q + oQ, p + oP, perm * scale * (*AzCz++));
                }
            }
        }

        if (aQ == atom2) {
            int C = do_full_field ? 3 * atom1 + 6 : 6;
            const double *CxBx = buffers[upper_triangle_index(num_buffers, C + 0, 3)];
            const double *CxBy = buffers[upper_triangle_index(num_buffers, C + 0, 4)];
            const double *CxBz = buffers[upper_triangle_index(num_buffers, C + 0, 5)];
            const double *CyBx = buffers[upper_triangle_index(num_buffers, C + 1, 3)];
            const double *CyBy = buffers[upper_triangle_index(num_buffers, C + 1, 4)];
            const double *CyBz = buffers[upper_triangle_index(num_buffers, C + 1, 5)];
            const double *CzBx = buffers[upper_triangle_index(num_buffers, C + 2, 3)];
            const double *CzBy = buffers[upper_triangle_index(num_buffers, C + 2, 4)];
            const double *CzBz = buffers[upper_triangle_index(num_buffers, C + 2, 5)];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    hess[0]->add(p + oP, q + oQ, perm * scale * (*CxBx));
                    hess[1]->add(p + oP, q + oQ, perm * (*CxBy));
                    hess[2]->add(p + oP, q + oQ, perm * (*CxBz));
                    hess[3]->add(p + oP, q + oQ, perm * (*CyBx));
                    hess[4]->add(p + oP, q + oQ, perm * scale * (*CyBy));
                    hess[5]->add(p + oP, q + oQ, perm * (*CyBz));
                    hess[6]->add(p + oP, q + oQ, perm * (*CzBx));
                    hess[7]->add(p + oP, q + oQ, perm * (*CzBy));
                    hess[8]->add(p + oP, q + oQ, perm * scale * (*CzBz));
                    hess[0]->add(q + oQ, p + oP, perm * scale * (*CxBx++));
                    hess[1]->add(q + oQ, p + oP, perm * (*CxBy++));
                    hess[2]->add(q + oQ, p + oP, perm * (*CxBz++));
                    hess[3]->add(q + oQ, p + oP, perm * (*CyBx++));
                    hess[4]->add(q + oQ, p + oP, perm * scale * (*CyBy++));
                    hess[5]->add(q + oQ, p + oP, perm * (*CyBz++));
                    hess[6]->add(q + oQ, p + oP, perm * (*CzBx++));
                    hess[7]->add(q + oQ, p + oP, perm * (*CzBy++));
                    hess[8]->add(q + oQ, p + oP, perm * scale * (*CzBz++));
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2};
    for (int p = 0; p < 9; p++) hess[p]->set_numpy_shape(nshape);

    return hess;
}

std::vector<SharedMatrix> MintsHelper::ao_overlap_kinetic_deriv2_helper(const std::string &type, int atom1, int atom2) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
    for (int p = 0; p < 3; p++) {
        for (int q = 0; q < 3; q++) {
            std::stringstream sstream;
            sstream << "ao_" << type << "_deriv2_" << atom1 << atom2 << cartcomp[p] << cartcomp[q];
            grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
        }
    }

    const auto &shell_pairs = GInt->shellpairs();
    size_t n_pairs = shell_pairs.size();

    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        const auto &shellP = bs1->shell(P);
        const auto &shellQ = bs2->shell(Q);

        int nP = shellP.nfunction();
        int oP = shellP.function_index();
        int aP = shellP.ncenter();

        int nQ = shellQ.nfunction();
        int oQ = shellQ.function_index();
        int aQ = shellQ.ncenter();

        if (aP != atom1 && aQ != atom1 && aP != atom2 && aQ != atom2) continue;

        GInt->compute_shell_deriv2(P, Q);
        const auto &buffers = GInt->buffers();

        double perm = P == Q ? 0.5 : 1.0;

        // This code makes use of the translational invariance relations
        //     ∂^2 S   ∂^2 S       ∂^2 S
        //     ----- = -----  =  - -----
        //     ∂A ∂A   ∂B ∂B       ∂A ∂B
        // and writes everything in terms of derivs w.r.t. center A only
        const double *pxx = buffers[0];
        const double *pxy = buffers[1];
        const double *pxz = buffers[2];
        const double *pyy = buffers[6];
        const double *pyz = buffers[7];
        const double *pzz = buffers[11];

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
                        grad[0]->set(p + oP, q + oQ, perm * 0);
                        grad[1]->set(p + oP, q + oQ, perm * tmpxy);
                        grad[2]->set(p + oP, q + oQ, perm * tmpxz);
                        grad[3]->set(p + oP, q + oQ, perm * -tmpxy);
                        grad[4]->set(p + oP, q + oQ, perm * 0);
                        grad[5]->set(p + oP, q + oQ, perm * tmpyz);
                        grad[6]->set(p + oP, q + oQ, perm * -tmpxz);
                        grad[7]->set(p + oP, q + oQ, perm * -tmpyz);
                        grad[8]->set(p + oP, q + oQ, perm * 0);
                        grad[0]->set(q + oQ, p + oP, perm * 0);
                        grad[1]->set(q + oQ, p + oP, perm * tmpxy);
                        grad[2]->set(q + oQ, p + oP, perm * tmpxz);
                        grad[3]->set(q + oQ, p + oP, perm * -tmpxy);
                        grad[4]->set(q + oQ, p + oP, perm * 0);
                        grad[5]->set(q + oQ, p + oP, perm * tmpyz);
                        grad[6]->set(q + oQ, p + oP, perm * -tmpxz);
                        grad[7]->set(q + oQ, p + oP, perm * -tmpyz);
                        grad[8]->set(q + oQ, p + oP, perm * 0);
                    } else {
                        grad[0]->set(p + oP, q + oQ, perm * tmpxx);
                        grad[1]->set(p + oP, q + oQ, perm * tmpxy);
                        grad[2]->set(p + oP, q + oQ, perm * tmpxz);
                        grad[4]->set(p + oP, q + oQ, perm * tmpyy);
                        grad[5]->set(p + oP, q + oQ, perm * tmpyz);
                        grad[8]->set(p + oP, q + oQ, perm * tmpzz);
                        grad[0]->set(q + oQ, p + oP, perm * tmpxx);
                        grad[1]->set(q + oQ, p + oP, perm * tmpxy);
                        grad[2]->set(q + oQ, p + oP, perm * tmpxz);
                        grad[4]->set(q + oQ, p + oP, perm * tmpyy);
                        grad[5]->set(q + oQ, p + oP, perm * tmpyz);
                        grad[8]->set(q + oQ, p + oP, perm * tmpzz);
                    }
                } else {
                    if (aP == atom1 && aQ == atom2) {
                        grad[0]->set(p + oP, q + oQ, perm * -tmpxx);
                        grad[1]->set(p + oP, q + oQ, perm * -tmpxy);
                        grad[2]->set(p + oP, q + oQ, perm * -tmpxz);
                        grad[3]->set(p + oP, q + oQ, perm * -tmpxy);
                        grad[4]->set(p + oP, q + oQ, perm * -tmpyy);
                        grad[5]->set(p + oP, q + oQ, perm * -tmpyz);
                        grad[6]->set(p + oP, q + oQ, perm * -tmpxz);
                        grad[7]->set(p + oP, q + oQ, perm * -tmpyz);
                        grad[8]->set(p + oP, q + oQ, perm * -tmpzz);
                        grad[0]->set(q + oQ, p + oP, perm * -tmpxx);
                        grad[1]->set(q + oQ, p + oP, perm * -tmpxy);
                        grad[2]->set(q + oQ, p + oP, perm * -tmpxz);
                        grad[3]->set(q + oQ, p + oP, perm * -tmpxy);
                        grad[4]->set(q + oQ, p + oP, perm * -tmpyy);
                        grad[5]->set(q + oQ, p + oP, perm * -tmpyz);
                        grad[6]->set(q + oQ, p + oP, perm * -tmpxz);
                        grad[7]->set(q + oQ, p + oP, perm * -tmpyz);
                        grad[8]->set(q + oQ, p + oP, perm * -tmpzz);
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

    return grad;
}

/* 1st derivatives of electric dipole integrals in the AO basis */
std::vector<SharedMatrix> MintsHelper::ao_elec_dip_deriv1_helper(int atom) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    std::shared_ptr<OneBodyAOInt> Dint(integral_->ao_multipoles(1, 1));

    std::shared_ptr<BasisSet> bs1 = Dint->basis1();
    std::shared_ptr<BasisSet> bs2 = Dint->basis2();
    const auto bs1_equiv_bs2 = (bs1 == bs2);

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_mu" << cartcomp[p] << "_deriv1_";
        for (int q = 0; q < 3; q++) {
            sstream << atom << cartcomp[q];
            grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1, nbf2));
        }
    }

    const auto &shell_pairs = Dint->shellpairs();
    size_t n_pairs = shell_pairs.size();

    for (size_t p = 0; p < n_pairs; ++p) {
        auto P = shell_pairs[p].first;
        auto Q = shell_pairs[p].second;
        int nP = basisset_->shell(P).nfunction();
        int oP = basisset_->shell(P).function_index();
        int aP = basisset_->shell(P).ncenter();

        int nQ = basisset_->shell(Q).nfunction();
        int oQ = basisset_->shell(Q).function_index();
        int aQ = basisset_->shell(Q).ncenter();

        if (aP != atom && aQ != atom) continue;

        Dint->compute_shell_deriv1(P, Q);
        const auto &buffers = Dint->buffers();

        for (int mu_cart = 0; mu_cart < 3; mu_cart++) {
            for (int atom_cart = 0; atom_cart < 3; atom_cart++) {
                const double *bufferP = buffers[6 * mu_cart + atom_cart];
                const double *bufferQ = buffers[6 * mu_cart + atom_cart + 3];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        if (atom == aP) {
                            grad[3 * mu_cart + atom_cart]->add(p + oP, q + oQ, *bufferP);
                            if (bs1_equiv_bs2 && P != Q) {
                                grad[3 * mu_cart + atom_cart]->add(q + oQ, p + oP, *bufferP);
                            }
                            bufferP++;
                        }
                        if (atom == aQ) {
                            grad[3 * mu_cart + atom_cart]->add(p + oP, q + oQ, *bufferQ);
                            if (bs1_equiv_bs2 && P != Q) {
                                grad[3 * mu_cart + atom_cart]->add(q + oQ, p + oP, *bufferQ);
                            }
                            bufferQ++;
                        }
                    }
                }
            }
        }
    }

    return grad;
}

/* 1st and 2nd derivatives of TEI in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_tei_deriv1(int atom, double omega,
                                                     std::shared_ptr<IntegralFactory> input_factory) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
        grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1 * nbf2, nbf3 * nbf4));
    }

    const auto &buffers = ints->buffers();
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

                                    Ax = buffers[0][delta];
                                    Ay = buffers[1][delta];
                                    Az = buffers[2][delta];
                                    Bx = buffers[3][delta];
                                    By = buffers[4][delta];
                                    Bz = buffers[5][delta];
                                    Cx = buffers[6][delta];
                                    Cy = buffers[7][delta];
                                    Cz = buffers[8][delta];
                                    Dx = buffers[9][delta];
                                    Dy = buffers[10][delta];
                                    Dz = buffers[11][delta];

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

/* 1st and 2nd derivatives of metric in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_metric_deriv1(int atom, const std::string &aux_name) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    auto aux = get_basisset(aux_name);
    auto factory =
        std::make_shared<IntegralFactory>(aux, BasisSet::zero_ao_basis_set(), aux, BasisSet::zero_ao_basis_set());

    auto ints = std::shared_ptr<TwoBodyAOInt>(factory->eri(1));
    auto naux = aux->nbf();
    int natom = basisset_->molecule()->natom();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_metric_deriv1_" << atom << cartcomp[p];
        grad.push_back(std::make_shared<Matrix>(sstream.str(), naux, naux));
    }

    const auto &buffers = ints->buffers();
    for (int P = 0; P < aux->nshell(); P++) {
        for (int Q = 0; Q < aux->nshell(); Q++) {
            int Psize = aux->shell(P).nfunction();
            int Qsize = aux->shell(Q).nfunction();

            int Poff = aux->shell(P).function_index();
            int Qoff = aux->shell(Q).function_index();

            int Pcenter = aux->shell(P).ncenter();
            int Qcenter = aux->shell(Q).ncenter();

            size_t delta;

            delta = 0L;

            if (Pcenter != atom && Qcenter != atom) continue;

            if (Pcenter == atom && Qcenter == atom) continue;

            ints->compute_shell_deriv1(P, 0, Q, 0);

            double Ax, Ay, Az;
            double Bx, By, Bz;
            double X = 0, Y = 0, Z = 0;

            for (int p = 0; p < Psize; p++) {
                for (int q = 0; q < Qsize; q++) {
                    int i = p + Poff;
                    int j = q + Qoff;

                    Ax = buffers[0][delta];
                    Ay = buffers[1][delta];
                    Az = buffers[2][delta];
                    Bx = buffers[3][delta];
                    By = buffers[4][delta];
                    Bz = buffers[5][delta];

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

                    grad[0]->set(i, j, X);
                    grad[1]->set(i, j, Y);
                    grad[2]->set(i, j, Z);

                    X = 0, Y = 0, Z = 0;
                    delta++;
                }
            }
        }
    }

    return grad;
}

/* 1st and 2nd derivatives of TEI in AO basis  */

std::vector<SharedMatrix> MintsHelper::ao_3center_deriv1(int atom, const std::string &aux_name) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    auto aux = get_basisset(aux_name);
    auto factory = std::make_shared<IntegralFactory>(aux, BasisSet::zero_ao_basis_set(), basisset_, basisset_);

    auto ints = std::shared_ptr<TwoBodyAOInt>(factory->eri(1));
    auto naux = aux->nbf();
    int natom = basisset_->molecule()->natom();

    std::vector<SharedMatrix> grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "ao_3center_deriv1_" << atom << cartcomp[p];
        grad.push_back(std::make_shared<Matrix>(sstream.str(), naux, nbf() * nbf()));
    }

    const auto &buffers = ints->buffers();
    for (int P = 0; P < aux->nshell(); P++) {
        for (int Q = 0; Q < basisset_->nshell(); Q++) {
            for (int R = 0; R < basisset_->nshell(); R++) {
                int Psize = aux->shell(P).nfunction();
                int Qsize = basisset_->shell(Q).nfunction();
                int Rsize = basisset_->shell(R).nfunction();

                int Pncart = aux->shell(P).ncartesian();
                int Qncart = basisset_->shell(Q).ncartesian();
                int Rncart = basisset_->shell(R).ncartesian();

                int Poff = aux->shell(P).function_index();
                int Qoff = basisset_->shell(Q).function_index();
                int Roff = basisset_->shell(R).function_index();

                int Pcenter = aux->shell(P).ncenter();
                int Qcenter = basisset_->shell(Q).ncenter();
                int Rcenter = basisset_->shell(R).ncenter();

                size_t delta;

                delta = 0L;

                if (Pcenter != atom && Qcenter != atom && Rcenter != atom) continue;

                if (Pcenter == atom && Qcenter == atom && Rcenter == atom) continue;

                ints->compute_shell_deriv1(P, 0, Q, R);

                double Ax, Ay, Az;
                double Bx, By, Bz;
                double Cx, Cy, Cz;
                double X = 0, Y = 0, Z = 0;

                for (int p = 0; p < Psize; p++) {
                    for (int q = 0; q < Qsize; q++) {
                        for (int r = 0; r < Rsize; r++) {
                            int i = Poff + p;
                            int j = (Qoff + q) * nbf() + Roff + r;

                            Ax = buffers[0][delta];
                            Ay = buffers[1][delta];
                            Az = buffers[2][delta];
                            Bx = buffers[3][delta];
                            By = buffers[4][delta];
                            Bz = buffers[5][delta];
                            Cx = buffers[6][delta];
                            Cy = buffers[7][delta];
                            Cz = buffers[8][delta];

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

    // Build numpy and final matrix shape
    std::vector<int> nshape{naux, nbf(), nbf()};
    for (int p = 0; p < 3; p++) grad[p]->set_numpy_shape(nshape);

    return grad;
}

std::vector<SharedMatrix> MintsHelper::ao_tei_deriv2(int atom1, int atom2) {
    /* NOTE: the x, y, and z in this vector must remain lowercase for this function */
    std::array<std::string, 3> cartcomp{{"x", "y", "z"}};

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints(nthreads);
    ints[0] = std::shared_ptr<TwoBodyAOInt>(integral_->eri(2));
    for (int thread = 1; thread < nthreads; thread++) ints[thread] = std::shared_ptr<TwoBodyAOInt>(ints.front()->clone());

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
            grad.push_back(std::make_shared<Matrix>(sstream.str(), nbf1 * nbf2, nbf3 * nbf4));
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
        const auto &buffers = ints[thread]->buffers();

        std::unordered_map<std::string, double> hess_map;

        for (int p = 0; p < Psize; p++) {
            for (int q = 0; q < Qsize; q++) {
                for (int r = 0; r < Rsize; r++) {
                    for (int s = 0; s < Ssize; s++) {
                        int i = (Poff + p) * nbf2 + Qoff + q;
                        int j = (Roff + r) * nbf4 + Soff + s;
                        hess_map["AxAx"] = buffers[0][delta];
                        hess_map["AxAy"] = buffers[1][delta];
                        hess_map["AxAz"] = buffers[2][delta];
                        hess_map["AxBx"] = buffers[3][delta];
                        hess_map["AxBy"] = buffers[4][delta];
                        hess_map["AxBz"] = buffers[5][delta];
                        hess_map["AxCx"] = buffers[6][delta];
                        hess_map["AxCy"] = buffers[7][delta];
                        hess_map["AxCz"] = buffers[8][delta];
                        hess_map["AxDx"] = buffers[9][delta];
                        hess_map["AxDy"] = buffers[10][delta];
                        hess_map["AxDz"] = buffers[11][delta];
                        hess_map["AyAy"] = buffers[12][delta];
                        hess_map["AyAz"] = buffers[13][delta];
                        hess_map["AyBx"] = buffers[14][delta];
                        hess_map["AyBy"] = buffers[15][delta];
                        hess_map["AyBz"] = buffers[16][delta];
                        hess_map["AyCx"] = buffers[17][delta];
                        hess_map["AyCy"] = buffers[18][delta];
                        hess_map["AyCz"] = buffers[19][delta];
                        hess_map["AyDx"] = buffers[20][delta];
                        hess_map["AyDy"] = buffers[21][delta];
                        hess_map["AyDz"] = buffers[22][delta];
                        hess_map["AzAz"] = buffers[23][delta];
                        hess_map["AzBx"] = buffers[24][delta];
                        hess_map["AzBy"] = buffers[25][delta];
                        hess_map["AzBz"] = buffers[26][delta];
                        hess_map["AzCx"] = buffers[27][delta];
                        hess_map["AzCy"] = buffers[28][delta];
                        hess_map["AzCz"] = buffers[29][delta];
                        hess_map["AzDx"] = buffers[30][delta];
                        hess_map["AzDy"] = buffers[31][delta];
                        hess_map["AzDz"] = buffers[32][delta];
                        hess_map["BxBx"] = buffers[33][delta];
                        hess_map["BxBy"] = buffers[34][delta];
                        hess_map["BxBz"] = buffers[35][delta];
                        hess_map["BxCx"] = buffers[36][delta];
                        hess_map["BxCy"] = buffers[37][delta];
                        hess_map["BxCz"] = buffers[38][delta];
                        hess_map["BxDx"] = buffers[39][delta];
                        hess_map["BxDy"] = buffers[40][delta];
                        hess_map["BxDz"] = buffers[41][delta];
                        hess_map["ByBy"] = buffers[42][delta];
                        hess_map["ByBz"] = buffers[43][delta];
                        hess_map["ByCx"] = buffers[44][delta];
                        hess_map["ByCy"] = buffers[45][delta];
                        hess_map["ByCz"] = buffers[46][delta];
                        hess_map["ByDx"] = buffers[47][delta];
                        hess_map["ByDy"] = buffers[48][delta];
                        hess_map["ByDz"] = buffers[49][delta];
                        hess_map["BzBz"] = buffers[50][delta];
                        hess_map["BzCx"] = buffers[51][delta];
                        hess_map["BzCy"] = buffers[52][delta];
                        hess_map["BzCz"] = buffers[53][delta];
                        hess_map["BzDx"] = buffers[54][delta];
                        hess_map["BzDy"] = buffers[55][delta];
                        hess_map["BzDz"] = buffers[56][delta];
                        hess_map["CxCx"] = buffers[57][delta];
                        hess_map["CxCy"] = buffers[58][delta];
                        hess_map["CxCz"] = buffers[59][delta];
                        hess_map["CxDx"] = buffers[60][delta];
                        hess_map["CxDy"] = buffers[61][delta];
                        hess_map["CxDz"] = buffers[62][delta];
                        hess_map["CyCy"] = buffers[63][delta];
                        hess_map["CyCz"] = buffers[64][delta];
                        hess_map["CyDx"] = buffers[65][delta];
                        hess_map["CyDy"] = buffers[66][delta];
                        hess_map["CyDz"] = buffers[67][delta];
                        hess_map["CzCz"] = buffers[68][delta];
                        hess_map["CzDx"] = buffers[69][delta];
                        hess_map["CzDy"] = buffers[70][delta];
                        hess_map["CzDz"] = buffers[71][delta];
                        hess_map["DxDx"] = buffers[72][delta];
                        hess_map["DxDy"] = buffers[73][delta];
                        hess_map["DxDz"] = buffers[74][delta];
                        hess_map["DyDy"] = buffers[75][delta];
                        hess_map["DyDz"] = buffers[76][delta];
                        hess_map["DzDz"] = buffers[77][delta];

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

std::vector<SharedMatrix> MintsHelper::ao_overlap_half_deriv1(const std::string &half_der_side, int atom) {
    std::vector<SharedMatrix> ao_grad;

    if (half_der_side == "LEFT" || half_der_side == "RIGHT")
        ao_grad = ao_overlap_half_deriv1_helper(half_der_side, atom);
    else
        throw PSIEXCEPTION("Not a valid choice of half derivative side: must be LEFT or RIGHT");

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
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    std::vector<SharedMatrix> ao_grad;
    ao_grad = ao_oei_deriv1(oei_type, atom);

    // Assuming C1 symmetry
    int nbf1 = ao_grad[0]->rowdim();
    int nbf2 = ao_grad[0]->coldim();

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "mo_" << oei_type << "_deriv1_" << atom << cartcomp[p];
        auto temp = std::make_shared<Matrix>(sstream.str(), nbf1, nbf2);
        temp->transform(C1, ao_grad[p], C2);
        mo_grad.push_back(temp);
    }
    return mo_grad;
}

std::vector<SharedMatrix> MintsHelper::mo_oei_deriv2(const std::string &oei_type, int atom1, int atom2, SharedMatrix C1,
                                                     SharedMatrix C2) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
            auto temp = std::make_shared<Matrix>(sstream.str(), nbf1, nbf2);
            temp->transform(C1, ao_grad[pq], C2);
            mo_grad.push_back(temp);
        }
    return mo_grad;
}

std::vector<SharedMatrix> MintsHelper::mo_overlap_half_deriv1(const std::string &half_der_side, int atom,
                                                              SharedMatrix C1, SharedMatrix C2) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    std::vector<SharedMatrix> ao_grad;
    ao_grad = ao_overlap_half_deriv1(half_der_side, atom);

    // Assuming C1 symmetry
    int nbf1 = ao_grad[0]->rowdim();
    int nbf2 = ao_grad[0]->coldim();

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0; p < 3; p++) {
        std::stringstream sstream;
        sstream << "mo_overlap_half_deriv1_" << atom << cartcomp[p];
        auto temp = std::make_shared<Matrix>(sstream.str(), nbf1, nbf2);
        temp->transform(C1, ao_grad[p], C2);
        mo_grad.push_back(temp);
    }
    return mo_grad;
}

/* Electric dipole moment derivatives in both AO and MO basis */

std::vector<SharedMatrix> MintsHelper::ao_elec_dip_deriv1(int atom) {
    std::vector<SharedMatrix> ao_grad = ao_elec_dip_deriv1_helper(atom);

    return ao_grad;
}

std::vector<SharedMatrix> MintsHelper::mo_elec_dip_deriv1(int atom, SharedMatrix C1, SharedMatrix C2) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

    std::vector<SharedMatrix> ao_grad = ao_elec_dip_deriv1(atom);

    // Assuming C1 symmetry
    int nbf1 = ao_grad[0]->rowdim();
    int nbf2 = ao_grad[0]->coldim();

    std::vector<SharedMatrix> mo_grad;
    for (int p = 0; p < 9; p++) {
        std::stringstream sstream;
        sstream << "mo_elec_dip_deriv1_" << atom << cartcomp[p];
        auto temp = std::make_shared<Matrix>(sstream.str(), nbf1, nbf2);
        temp->transform(C1, ao_grad[p], C2);
        mo_grad.push_back(temp);
    }

    return mo_grad;
}

/*  TEI derivatives in  MO basis */

std::vector<SharedMatrix> MintsHelper::mo_tei_deriv1(int atom, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                                     SharedMatrix C4) {
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
    std::array<std::string, 3> cartcomp{{"X", "Y", "Z"}};

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
