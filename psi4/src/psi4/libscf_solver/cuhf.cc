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

#include "cuhf.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"

#include "psi4/libfock/jk.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

using namespace psi;

namespace psi {
namespace scf {

CUHF::CUHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

CUHF::CUHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

CUHF::~CUHF() {}

void CUHF::common_init() {
    name_ = "CUHF";

    Fa_ = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_ = SharedMatrix(factory_->create_matrix("F beta"));
    Fp_ = SharedMatrix(factory_->create_matrix("F charge"));
    Fm_ = SharedMatrix(factory_->create_matrix("F spin"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Dp_ = SharedMatrix(factory_->create_matrix("D charge"));
    Dt_ = SharedMatrix(factory_->create_matrix("D total"));
    Da_old_ = SharedMatrix(factory_->create_matrix("D alpha old"));
    Db_old_ = SharedMatrix(factory_->create_matrix("D beta  old"));
    Dt_old_ = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_ = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_ = SharedMatrix(factory_->create_matrix("C beta"));
    Cno_ = SharedMatrix(factory_->create_matrix("C NOs"));
    Cno_temp_ = SharedMatrix(factory_->create_matrix("C NO temp"));
    J_ = SharedMatrix(factory_->create_matrix("J total"));
    Ka_ = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_ = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());
    No_ = SharedVector(factory_->create_vector());
    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;

    subclass_init();
}

void CUHF::damping_update(double damping_percentage) {
    Da_->scale(1.0 - damping_percentage);
    Da_->axpy(damping_percentage, Da_old_);
    Db_->scale(1.0 - damping_percentage);
    Db_->axpy(damping_percentage, Db_old_);
    Dt_->copy(Da_);
    Dt_->add(Db_);
}

void CUHF::finalize() {
    // Form lagrangian
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < Lagrangian_->rowdim(h); ++m) {
            for (int n = 0; n < Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i = 0; i < nalphapi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);
                }
                for (int i = 0; i < nbetapi_[h]; ++i)
                    sum += epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);

                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dt_.reset();
    Da_old_.reset();
    Db_old_.reset();
    Dt_old_.reset();
    Dp_.reset();
    Fp_.reset();
    Fm_.reset();
    Cno_.reset();
    Cno_temp_.reset();
    No_.reset();

    HF::finalize();
}

void CUHF::save_density_and_energy() {
    Da_old_->copy(Dt_);
    Db_old_->copy(Dt_);
    Dt_old_->copy(Dt_);
}

void CUHF::form_G() {
    // Push the C matrix on
    std::vector<SharedMatrix>& C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();
    J_->copy(J[0]);
    J_->add(J[1]);
    Ka_ = K[0];
    Kb_ = K[1];
}

void CUHF::compute_spin_contamination() {
    double dN = 0.0;

    for (int h = 0; h < S_->nirrep(); h++) {
        int nbf = S_->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0) continue;

        auto Ht = std::make_shared<Matrix>("H Temp", nbf, nb);
        auto Ft = std::make_shared<Matrix>("F Temp", na, nb);

        double** Sp = S_->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N', 'N', nbf, nb, nbf, 1.0, Sp[0], nbf, Cbp[0], nmo, 0.0, Htp[0], nb);
        C_DGEMM('T', 'N', na, nb, nbf, 1.0, Cap[0], nmo, Htp[0], nb, 0.0, Ftp[0], nb);

        for (long int ab = 0; ab < (long int)na * nb; ab++) dN += Ftp[0][ab] * Ftp[0][ab];
    }

    double dS = (double)nbeta_ - (double)dN;

    double nm = (nalpha_ - nbeta_) / 2.0;
    double S2 = nm * (nm + 1.0);

    outfile->Printf("\n  @Spin Contamination Metric: %8.5F\n", dS);
    outfile->Printf("  @S^2 Expected:              %8.5F\n", S2);
    outfile->Printf("  @S^2 Observed:              %8.5F\n", S2 + dS);
}

void CUHF::form_initial_F() {
    // Form the initial Fock matrix to get initial orbitals
    Fp_->copy(J_);
    Fp_->scale(2.0);
    Fp_->subtract(Ka_);
    Fp_->subtract(Kb_);
    Fp_->scale(0.5);

    Fa_->copy(H_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }
    Fa_->add(Fp_);

    // Just reuse alpha for beta
    Fb_->copy(Fa_);

    if (debug_) {
        outfile->Printf("Initial Fock alpha matrix:\n");
        Fa_->print("outfile");
        outfile->Printf("Initial Fock beta matrix:\n");
        Fb_->print("outfile");
    }
}

void CUHF::form_F() {
    // Form (rho_a + rho_b) / 2
    Dp_->copy(Dt_);
    Dp_->scale(-0.5);  // This is a hack to get the eigenvectors in the
                       // order that I want
    if (debug_) {
        outfile->Printf("Charge Density Matrix (SO Basis):\n");
        Dp_->print();
    }

    // Transfrom to the orthonormal basis
    Dp_->transform(S_);
    Dp_->transform(X_);
    if (debug_) {
        outfile->Printf("Charge Density Matrix (Orthonormal Basis):\n");
        Dp_->print();
    }

    // Diagonalize the charge density and form the natural orbitals
    Dp_->diagonalize(Cno_temp_, No_);
    if (debug_) {
        outfile->Printf("CUHF Natural Orbital Occupations:\n");
        No_->print();
    }
    Cno_->gemm(false, false, 1.0, X_, Cno_temp_, 0.0);

    // Now we form the contributions to the Fock matrix from
    // the charge and spin densities
    Fp_->copy(J_);
    Fp_->scale(2.0);
    Fp_->subtract(Ka_);
    Fp_->subtract(Kb_);
    Fp_->scale(0.5);

    Fm_->copy(Ka_);
    Fm_->subtract(Kb_);
    Fm_->scale(-0.5);

    // Transform the spin density contributions to the NO basis
    Fm_->transform(Cno_);

    // Zero the core-virtual contributions
    //
    //            [ Fm_cc Fm_co   0   ]
    // Fm_tilde = [ Fm_oc Fm_oo Fm_ov ]
    //            [   0   Fm_vo Fm_vv ]
    //
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < nbetapi_[h]; ++i) {
            for (int j = nalphapi_[h]; j < nmopi_[h]; ++j) {
                Fm_->set(h, i, j, 0.0);
                Fm_->set(h, j, i, 0.0);
            }
        }
    }

    // Return to the SO basis
    Fm_->back_transform(Cno_);
    Fm_->transform(S_);

    // Build the modified alpha and beta Fock matrices
    Fa_->copy(H_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }
    Fa_->add(Fp_);
    Fa_->add(Fm_);

    Fb_->copy(H_);
    for (const auto& Vext : external_potentials_) {
        Fb_->add(Vext);
    }
    Fb_->add(Fp_);
    Fb_->subtract(Fm_);

    if (debug_) {
        Fa_->print("outfile");
        Fb_->print("outfile");
    }
}

void CUHF::form_C(double shift) {
    if (shift == 0.0) {
        diagonalize_F(Fa_, Ca_, epsilon_a_);
        diagonalize_F(Fb_, Cb_, epsilon_b_);
    } else {
        auto shifted_F = SharedMatrix(factory_->create_matrix("F"));

        auto Cvir = Ca_subset("SO", "VIR");
        auto SCvir = std::make_shared<Matrix>(nirrep_, S_->rowspi(), Cvir->colspi());
        SCvir->gemm(false, false, 1.0, S_, Cvir, 0.0);
        shifted_F->gemm(false, true, shift, SCvir, SCvir, 0.0);
        shifted_F->add(Fa_);
        diagonalize_F(shifted_F, Ca_, epsilon_a_);

        Cvir = Cb_subset("SO", "VIR");
        SCvir = std::make_shared<Matrix>(nirrep_, S_->rowspi(), Cvir->colspi());
        SCvir->gemm(false, false, 1.0, S_, Cvir, 0.0);
        shifted_F->gemm(false, true, shift, SCvir, SCvir, 0.0);
        shifted_F->add(Fb_);
        diagonalize_F(shifted_F, Cb_, epsilon_b_);
    }
    find_occupation();
    if (debug_) {
        Ca_->print("outfile");
        Cb_->print("outfile");
    }
}

void CUHF::form_D() {
    Da_->zero();
    Db_->zero();

    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** Cb = Cb_->pointer(h);
        double** Da = Da_->pointer(h);
        double** Db = Db_->pointer(h);

        C_DGEMM('N', 'T', nso, nso, na, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, Da[0], nso);
        C_DGEMM('N', 'T', nso, nso, nb, 1.0, Cb[0], nmo, Cb[0], nmo, 0.0, Db[0], nso);
    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        outfile->Printf("in CUHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

double CUHF::compute_initial_E() { return nuclearrep_ + Dt_->vector_dot(H_); }

double CUHF::compute_E() {
    double DH = Dt_->vector_dot(H_);
    double DT = Dt_->vector_dot(T_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);

    double one_electron_E = DH;
    double two_electron_E = 0.5 * (DFa + DFb - one_electron_E);
    double Eelec = 0.5 * (DH + DFa + DFb);

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = DT;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["VV10_E"] = 0.0;
    energies_["-D"] = 0.0;

    // outfile->Printf( "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

bool CUHF::stability_analysis() {
    throw PSIEXCEPTION("CUHF stability analysis has not been implemented yet.  Sorry :(");
    return false;
}

std::shared_ptr<CUHF> CUHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    std::shared_ptr<Wavefunction> wfn = Wavefunction::c1_deep_copy(basis);
    auto hf_wfn = std::make_shared<CUHF>(wfn, functional_, wfn->options(), wfn->psio());

    // now just have to copy the matrices that UHF initializes
    // include only those that are not temporary (some deleted in finalize())
    if (Ca_) hf_wfn->Ca_ = Ca_subset("AO", "ALL");
    if (Cb_) hf_wfn->Cb_ = Cb_subset("AO", "ALL");
    if (Da_) hf_wfn->Da_ = Da_subset("AO");
    if (Db_) hf_wfn->Db_ = Db_subset("AO");
    if (Fa_) hf_wfn->Fa_ = Fa_subset("AO");
    if (Fb_) hf_wfn->Fb_ = Fb_subset("AO");
    if (epsilon_a_) hf_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nalphapi_, "AO", "ALL");
    if (epsilon_b_) hf_wfn->epsilon_b_ = epsilon_subset_helper(epsilon_b_, nbetapi_, "AO", "ALL");
    // H_ ans X_ reset in the HF constructor, copy them over here
    SharedMatrix SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}

void CUHF::compute_SAD_guess(bool natorb) {
    // Form the SAD guess
    HF::compute_SAD_guess(natorb);
    if (!natorb) {
        // Form the total density used in energy evaluation
        Dt_->copy(Da_);
        Dt_->add(Db_);
    }
}

void CUHF::setup_potential() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("CUHF: Cannot compute XC components!");
    }
}
}  // namespace scf
}  // namespace psi
