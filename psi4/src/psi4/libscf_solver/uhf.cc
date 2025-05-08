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

#ifdef USING_OpenOrbitalOptimizer
#include <openorbitaloptimizer/scfsolver.hpp>
#endif

#include "uhf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "psi4/physconst.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_types.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianEnableDFT;

#endif

namespace psi {
namespace scf {

UHF::UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

UHF::UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
         std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

UHF::~UHF() {}

void UHF::common_init() {
    name_ = "UHF";

    mix_performed_ = false;

    // TODO: Move that to the base object
    step_scale_ = options_.get_double("FOLLOW_STEP_SCALE");
    step_increment_ = options_.get_double("FOLLOW_STEP_INCREMENT");

    Fa_ = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_ = SharedMatrix(factory_->create_matrix("F beta"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Da_old_ = SharedMatrix(factory_->create_matrix("Old alpha SCF density"));
    Db_old_ = SharedMatrix(factory_->create_matrix("Old beta SCF density"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_ = SharedMatrix(factory_->create_matrix("alpha MO coefficients (C)"));
    Cb_ = SharedMatrix(factory_->create_matrix("beta MO coefficients (C)"));
    Ga_ = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_ = SharedMatrix(factory_->create_matrix("G beta"));
    Va_ = SharedMatrix(factory_->create_matrix("V alpha"));
    Vb_ = SharedMatrix(factory_->create_matrix("V beta"));
    J_ = SharedMatrix(factory_->create_matrix("J total"));
    Ka_ = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_ = SharedMatrix(factory_->create_matrix("K beta"));
    wKa_ = SharedMatrix(factory_->create_matrix("wK alpha"));
    wKb_ = SharedMatrix(factory_->create_matrix("wK beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("alpha orbital energies");
    epsilon_b_ = SharedVector(factory_->create_vector());
    epsilon_b_->set_name("beta orbital energies");

    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;

    subclass_init();
}

void UHF::finalize() {
    // Form lagrangian
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < Lagrangian_->rowdim(h); ++m) {
            for (int n = 0; n < Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i = 0; i < nalphapi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);
                }
                for (int i = 0; i < nbetapi_[h]; ++i) {
                    sum += epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);
                }
                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Da_old_.reset();
    Db_old_.reset();
    Ga_.reset();
    Gb_.reset();

    compute_nos();

    HF::finalize();
}

void UHF::save_density_and_energy() {
    Da_old_->copy(Da_);
    Db_old_->copy(Db_);
}
void UHF::form_V() {
    // // Push the C matrix on
    // std::vector<SharedMatrix> & C = potential_->C();
    // C.clear();
    // C.push_back(Ca_subset("SO", "OCC"));
    // C.push_back(Cb_subset("SO", "OCC"));

    // // Run the potential object
    // potential_->compute();

    // // Pull the V matrices off
    // const std::vector<SharedMatrix> & V = potential_->V();
    // Va_ = V[0];
    // Vb_ = V[1];
    potential_->set_D({Da_, Db_});
    potential_->compute_V({Va_, Vb_});
    // Vb_ = Va_;
}
void UHF::form_G() {
    if (functional_->needs_xc()) {
        form_V();
        Ga_->copy(Va_);
        Gb_->copy(Vb_);
    } else {
        Ga_->zero();
        Gb_->zero();
    }

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
    const std::vector<SharedMatrix>& wK = jk_->wK();
    J_->copy(J[0]);
    J_->add(J[1]);
    if (functional_->is_x_hybrid()) {
        Ka_ = K[0];
        Kb_ = K[1];
    }
    if (functional_->is_x_lrc()) {
        wKa_ = wK[0];
        wKb_ = wK[1];
    }
    Ga_->add(J_);
    Gb_->add(J_);

    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    if (functional_->is_x_hybrid() && !(functional_->is_x_lrc() && jk_->get_wcombine())) {
        Ga_->axpy(-alpha, Ka_);
        Gb_->axpy(-alpha, Kb_);
    } else {
        Ka_->zero();
        Kb_->zero();
    }

    if (functional_->is_x_lrc()) {
        if (jk_->get_wcombine()) {
            Ga_->axpy(-1.0, wKa_);
            Gb_->axpy(-1.0, wKb_);
        } else {
            Ga_->axpy(-beta, wKa_);
            Gb_->axpy(-beta, wKb_);
        }
    } else {
        wKa_->zero();
        wKb_->zero();
    }
}

void UHF::form_F() {
    Fa_->copy(H_);
    Fa_->add(Ga_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }

    Fb_->copy(H_);
    Fb_->add(Gb_);
    for (const auto& Vext : external_potentials_) {
        Fb_->add(Vext);
    }

    if (debug_) {
        Fa_->print("outfile");
        Fb_->print("outfile");
    }
}

void UHF::form_C(double shift) {
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
    if (options_.get_bool("GUESS_MIX") && !mix_performed_) {
        if (Ca_->nirrep() != 1) {
            outfile->Printf("  Mixing alpha HOMO/LUMO orbitals (%d,%d)\n", nalpha_, nalpha_ + 1);
            throw InputException("Warning: cannot mix alpha HOMO/LUMO orbitals. Run in C1 symmetry.",
                                 "to 'symmetry c1'", __FILE__, __LINE__);
        }

        // SAD doesn't have orbitals in iteration 0, other guesses do
        bool have_orbitals = !sad_ || (sad_ && iteration_ > 0);
        if (have_orbitals) {
            Ca_->rotate_columns(0, nalpha_ - 1, nalpha_, pc_pi * 0.25);
            if (nbeta_ > 0) {
                outfile->Printf("  Mixing beta HOMO/LUMO orbitals (%d,%d)\n", nbeta_, nbeta_ + 1);
                Cb_->rotate_columns(0, nbeta_ - 1, nbeta_, -pc_pi * 0.25);
            }
            mix_performed_ = true;

            // Since we've changed the orbitals, delete the DIIS history
            // so that we don't fall back to spin-restricted orbitals
            if (initialized_diis_manager_) {
                diis_manager_.attr("delete_diis_file")();
                diis_manager_ = py::none();
                initialized_diis_manager_ = false;
            }
        }
    }
    find_occupation();
    if (debug_) {
        Ca_->print("outfile");
        Cb_->print("outfile");
    }
}

void UHF::form_D() {
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

    if (debug_) {
        outfile->Printf("in UHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}
void UHF::damping_update(double damping_percentage) {
    Da_->scale(1.0 - damping_percentage);
    Da_->axpy(damping_percentage, Da_old_);
    Db_->scale(1.0 - damping_percentage);
    Db_->axpy(damping_percentage, Db_old_);
}

double UHF::compute_initial_E() {
    auto Dt = Da_->clone();
    Dt->add(Db_);
    return nuclearrep_ + 0.5 * (Dt->vector_dot(H_));
}

double UHF::compute_E() {
    // E_DFT = 2.0 D*H + D*J - \alpha D*K + E_xc
    double kinetic_E = Da_->vector_dot(T_);
    kinetic_E += Db_->vector_dot(T_);
    double one_electron_E = Da_->vector_dot(H_);
    one_electron_E += Db_->vector_dot(H_);
    double coulomb_E = Da_->vector_dot(J_);
    coulomb_E += Db_->vector_dot(J_);

    double XC_E = 0.0;
    double VV10_E = 0.0;
    if (functional_->needs_xc()) {
        XC_E = potential_->quadrature_values()["FUNCTIONAL"];
    }
    if (functional_->needs_vv10()) {
        VV10_E = potential_->quadrature_values()["VV10"];
    }

    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    double exchange_E = 0.0;
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha * Da_->vector_dot(Ka_);
        exchange_E -= alpha * Db_->vector_dot(Kb_);
    }
    if (functional_->is_x_lrc()) {
        if (jk_->get_do_wK() && jk_->get_wcombine()) {
            exchange_E -= Da_->vector_dot(wKa_);
            exchange_E -= Db_->vector_dot(wKb_);
        } else {
            exchange_E -= beta * Da_->vector_dot(wKa_);
            exchange_E -= beta * Db_->vector_dot(wKb_);
        }
    }

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = kinetic_E;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * (coulomb_E + exchange_E);
    energies_["XC"] = XC_E;
    energies_["VV10"] = VV10_E;
    energies_["-D"] = scalar_variable("-D Energy");
    double dashD_E = energies_["-D"];

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * coulomb_E;
    Etotal += 0.5 * exchange_E;
    Etotal += XC_E;
    Etotal += VV10_E;
    Etotal += dashD_E;

    return Etotal;
}
std::vector<SharedMatrix> UHF::onel_Hx(std::vector<SharedMatrix> x_vec) {
    if ((x_vec.size() % 2) != 0) {
        throw PSIEXCEPTION("UHF::onel_Hx expect incoming vector to alternate A/B");
    }

    // This is a bypass for C1 input
    std::vector<bool> c1_input_;

    bool needs_ao = false;
    bool needs_so = false;
    for (size_t i = 0; i < x_vec.size(); i++) {
        if ((x_vec[i]->nirrep() == 1) && (nirrep_ != 1)) {
            c1_input_.push_back(true);
            needs_ao = true;
        } else {
            c1_input_.push_back(false);
            needs_so = true;
        }
    }

    SharedMatrix Caocc_ao, Cavir_ao, Fa_ao, Caocc_so, Cavir_so;
    SharedMatrix Cbocc_ao, Cbvir_ao, Fb_ao, Cbocc_so, Cbvir_so;

    if (needs_ao) {
        Caocc_ao = Ca_subset("AO", "OCC");
        Cbocc_ao = Cb_subset("AO", "OCC");
        Cavir_ao = Ca_subset("AO", "VIR");
        Cbvir_ao = Cb_subset("AO", "VIR");
        Fa_ao = matrix_subset_helper(Fa_, Ca_, "AO", "Fock");
        Fb_ao = matrix_subset_helper(Fb_, Cb_, "AO", "Fock");
    }

    if (needs_so) {
        Caocc_so = Ca_subset("SO", "OCC");
        Cavir_so = Ca_subset("SO", "VIR");
        Cbocc_so = Cb_subset("SO", "OCC");
        Cbvir_so = Cb_subset("SO", "VIR");
    }

    // Compute Fij x_ia - Fab x_ia
    SharedMatrix Fa, Fb, Cao, Cbo, Cav, Cbv;
    std::vector<SharedMatrix> ret;
    for (size_t i = 0; i < x_vec.size() / 2; i++) {
        if (c1_input_[i]) {
            if ((x_vec[2 * i]->rowspi() != Caocc_ao->colspi()) || (x_vec[2 * i]->colspi() != Cavir_ao->colspi())) {
                throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Fa = Fa_ao;
            Fb = Fb_ao;

            Cao = Caocc_ao;
            Cbo = Cbocc_ao;
            Cav = Cavir_ao;
            Cbv = Cbvir_ao;

        } else {
            if ((x_vec[2 * i + 1]->rowspi() != Cbocc_so->colspi()) ||
                (x_vec[2 * i + 1]->colspi() != Cbvir_so->colspi())) {
                throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Fa = Fa_;
            Fb = Fb_;

            Cao = Caocc_so;
            Cbo = Cbocc_so;
            Cav = Cavir_so;
            Cbv = Cbvir_so;
        }
        // Alpha
        SharedMatrix tmp1 = linalg::triplet(Cao, Fa, Cao, true, false, false);
        SharedMatrix result = linalg::doublet(tmp1, x_vec[2 * i], false, false);

        SharedMatrix tmp2 = linalg::triplet(x_vec[2 * i], Cav, Fa, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Cav, 1.0);

        ret.push_back(result);

        // Beta
        tmp1 = linalg::triplet(Cbo, Fb, Cbo, true, false, false);
        result = linalg::doublet(tmp1, x_vec[2 * i + 1], false, false);

        tmp2 = linalg::triplet(x_vec[2 * i + 1], Cbv, Fb, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Cbv, 1.0);

        ret.push_back(result);
    }

    return ret;
}
std::vector<SharedMatrix> UHF::twoel_Hx(std::vector<SharedMatrix> x_vec, bool combine, std::string return_basis) {
    if ((x_vec.size() % 2) != 0) {
        throw PSIEXCEPTION("UHF::twoel_Hx expect incoming vector to alternate A/B");
    }
    // This is a bypass for C1 input
    std::vector<bool> c1_input_;

    bool needs_ao = false;
    bool needs_so = false;
    for (size_t i = 0; i < x_vec.size(); i++) {
        if ((x_vec[i]->nirrep() == 1) && (nirrep_ != 1)) {
            c1_input_.push_back(true);
            needs_ao = true;
        } else {
            c1_input_.push_back(false);
            needs_so = true;
        }
    }

    SharedMatrix Caocc_ao, Cavir_ao, Caocc_so, Cavir_so;
    SharedMatrix Cbocc_ao, Cbvir_ao, Cbocc_so, Cbvir_so;

    if (needs_ao) {
        Caocc_ao = Ca_subset("AO", "OCC");
        Cavir_ao = Ca_subset("AO", "VIR");
        Cbocc_ao = Cb_subset("AO", "OCC");
        Cbvir_ao = Cb_subset("AO", "VIR");
    } else {
        Caocc_so = Ca_subset("SO", "OCC");
        Cavir_so = Ca_subset("SO", "VIR");
        Cbocc_so = Cb_subset("SO", "OCC");
        Cbvir_so = Cb_subset("SO", "VIR");
    }

    // Setup jk
    auto& Cl = jk_->C_left();
    auto& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    int nvecs = x_vec.size() / 2;

    // We actually want to compute all alpha then all beta, its smart enough to figure out the half transform
    SharedMatrix Cao, Cbo, Cav, Cbv;
    for (size_t i = 0; i < nvecs; i++) {
        if (c1_input_[i]) {
            if ((x_vec[2 * i]->rowspi() != Caocc_ao->colspi()) || (x_vec[2 * i]->colspi() != Cavir_ao->colspi())) {
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Cao = Caocc_ao;
            Cav = Cavir_ao;
        } else {
            if ((x_vec[2 * i]->rowspi() != Caocc_so->colspi()) || (x_vec[2 * i]->colspi() != Cavir_so->colspi())) {
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Cao = Caocc_so;
            Cav = Cavir_so;
        }

        Cl.push_back(Cao);

        SharedMatrix R = linalg::doublet(Cav, x_vec[2 * i], false, true);
        R->scale(-1.0);
        Cr.push_back(R);
    }
    for (size_t i = 0; i < nvecs; i++) {
        if (c1_input_[i]) {
            if ((x_vec[2 * i + 1]->rowspi() != Cbocc_ao->colspi()) ||
                (x_vec[2 * i + 1]->colspi() != Cbvir_ao->colspi())) {
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Cbo = Cbocc_ao;
            Cbv = Cbvir_ao;
        } else {
            if ((x_vec[2 * i + 1]->rowspi() != Cbocc_so->colspi()) ||
                (x_vec[2 * i + 1]->colspi() != Cbvir_so->colspi())) {
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Cbo = Cbocc_so;
            Cbv = Cbvir_so;
        }
        Cl.push_back(Cbo);

        SharedMatrix R = linalg::doublet(Cbv, x_vec[2 * i + 1], false, true);
        R->scale(-1.0);
        Cr.push_back(R);
    }

    // Compute JK
    // J_mn = I_mnpq Ca_pi Ca_qa Xa_ai, I_mnpq Cb_pi Cb_qa Xb_ai pairs
    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();
    const std::vector<SharedMatrix>& wK = jk_->wK();

    std::vector<SharedMatrix> Vx;
    if (functional_->needs_xc()) {
        std::vector<SharedMatrix> Dx;
        // Gotta reorder the wizardry
        for (size_t i = 0; i < nvecs; i++) {
            auto Dx_a = linalg::doublet(Cl[i], Cr[i], false, true);
            auto Dx_b = linalg::doublet(Cl[nvecs + i], Cr[nvecs + i], false, true);
            Vx.push_back(std::make_shared<Matrix>("Vax temp", Dx_a->rowspi(), Dx_a->colspi(), Dx_a->symmetry()));
            Vx.push_back(std::make_shared<Matrix>("Vbx temp", Dx_b->rowspi(), Dx_b->colspi(), Dx_b->symmetry()));
            Dx.push_back(Dx_a);
            Dx.push_back(Dx_b);
        }
        potential_->compute_Vx(Dx, Vx);
    }

    std::vector<SharedMatrix> V_ext_pert;
    for (const auto& pert : external_cpscf_perturbations_) {
        if (print_ > 1) outfile->Printf("Adding external CPSCF contribution %s.\n", pert.first.c_str());
        for (size_t i = 0; i < nvecs; i++) {
            auto Dx_a = linalg::doublet(Cl[i], Cr[i], false, true);
            auto Dx_b = linalg::doublet(Cl[nvecs + i], Cr[nvecs + i], false, true);
            V_ext_pert.push_back(pert.second(Dx_a));
            V_ext_pert.push_back(pert.second(Dx_b));
        }
    }

    Cl.clear();
    Cr.clear();

    // Build return vector
    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    std::vector<SharedMatrix> ret;
    if (combine) {
        for (size_t i = 0; i < nvecs; i++) {
            // Scale the Xa -> alpha Coulomb terms.
            J[i]->scale(2.0);
            // Add the Xb -> alpha Coulomb terms.
            J[i]->axpy(2.0, J[nvecs + i]);
            // X -> alpha terms in the AO/SO basis are the same as the X -> beta.
            J[nvecs + i]->copy(J[i]);
            if (functional_->needs_xc()) {
                J[i]->axpy(2.0, Vx[2 * i]);
                J[nvecs + i]->axpy(2.0, Vx[2 * i + 1]);
            }
            if (functional_->is_x_hybrid()) {
                J[i]->axpy(-alpha, K[i]);
                J[i]->axpy(-alpha, K[i]->transpose());

                J[nvecs + i]->axpy(-alpha, K[nvecs + i]);
                J[nvecs + i]->axpy(-alpha, K[nvecs + i]->transpose());
            }
            if (functional_->is_x_lrc()) {
                J[i]->axpy(-beta, wK[i]);
                J[i]->axpy(-beta, wK[i]->transpose());

                J[nvecs + i]->axpy(-beta, wK[nvecs + i]);
                J[nvecs + i]->axpy(-beta, wK[nvecs + i]->transpose());
            }
            if (V_ext_pert.size()) {
                J[i]->axpy(2.0, V_ext_pert[2 * i]);
                J[nvecs + i]->axpy(2.0, V_ext_pert[2 * i + 1]);
            }
            ret.push_back(J[i]);
            ret.push_back(J[nvecs + i]);
        }
    } else {
        if (jk_->get_wcombine()) {
            throw PSIEXCEPTION(
                "UHF::twoel_Hx user asked for wcombine but combine==false in SCF::twoel_Hx. Please set wcombine false "
                "in your input.");
        }
        for (size_t i = 0; i < nvecs; i++) {
            // Add opposite-spin terms to alpha.
            J[i]->add(J[nvecs + i]);
            // Before MO transformation, beta and alpha are the same.
            J[nvecs + i]->copy(J[i]);
            // Augment same-spin Coulomb terms with exchange-correlation.
            if (functional_->needs_xc()) {
                J[i]->add(Vx[2 * i]);
                J[nvecs + i]->add(Vx[2 * i + 1]);
            }
            if (V_ext_pert.size()) {
                J[i]->add(V_ext_pert[2 * i]);
                J[nvecs + i]->add(V_ext_pert[2 * i + 1]);
            }
            ret.push_back(J[i]);
            ret.push_back(J[nvecs + i]);
            if (functional_->is_x_hybrid()) {
                K[i]->scale(alpha);
                K[nvecs + i]->scale(alpha);
                if (functional_->is_x_lrc()) {
                    K[i]->axpy(beta, wK[i]);
                    K[nvecs + i]->axpy(beta, wK[nvecs + i]);
                }
                ret.push_back(K[i]);
                ret.push_back(K[nvecs + i]);
            } else if (functional_->is_x_lrc()) {
                wK[i]->scale(beta);
                wK[nvecs + i]->scale(beta);
                ret.push_back(wK[i]);
                ret.push_back(wK[nvecs + i]);
            }
        }
    }

    // Transform if needed
    if (return_basis == "SO") {
        /* pass */
    } else if (return_basis == "MO") {
        for (size_t i = 0; i < nvecs; i++) {
            if (c1_input_[i]) {
                ret[2 * i] = linalg::triplet(Caocc_ao, ret[2 * i], Cavir_ao, true, false, false);
                ret[2 * i + 1] = linalg::triplet(Cbocc_ao, ret[2 * i + 1], Cbvir_ao, true, false, false);
            } else {
                ret[2 * i] = linalg::triplet(Caocc_so, ret[2 * i], Cavir_so, true, false, false);
                ret[2 * i + 1] = linalg::triplet(Cbocc_so, ret[2 * i + 1], Cbvir_so, true, false, false);
            }
        }
    } else {
        throw PSIEXCEPTION("SCF::twoel_Hx: return_basis option not understood.");
    }

    return ret;
}
std::vector<SharedMatrix> UHF::cphf_Hx(std::vector<SharedMatrix> x_vec) {
    // Compute quantities
    auto onel = onel_Hx(x_vec);
    auto twoel = twoel_Hx(x_vec, true, "MO");

    for (size_t i = 0; i < onel.size(); i++) {
        onel[i]->add(twoel[i]);
    }

    return onel;
}
std::vector<SharedMatrix> UHF::cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol, int max_iter,
                                          int print_lvl) {
    if ((x_vec.size() % 2) != 0) {
        throw PSIEXCEPTION("UHF::cphf_solve expect incoming vector to alternate A/B");
    }

    std::time_t start, stop;
    start = std::time(nullptr);
    cphf_converged_ = false;
    cphf_nfock_builds_ = 0;

    // => Figure out the type of perturbation tensor <= //
    // This helps get around difficulties with non-totally symmetric perturbations
    std::vector<bool> c1_input_;

    bool needs_ao = false;
    bool needs_so = false;
    for (size_t i = 0; i < x_vec.size(); i++) {
        if ((x_vec[i]->nirrep() == 1) && (nirrep_ != 1)) {
            c1_input_.push_back(true);
            needs_ao = true;
        } else {
            c1_input_.push_back(false);
            needs_so = true;
        }
    }

    // => Build preconditioner <= //
    SharedMatrix Precon_ao_a, Precon_ao_b, Precon_so_a, Precon_so_b;

    if (needs_ao) {
        // MO (C1) Fock Matrix (Inactive Fock in Helgaker's language)
        auto Caocc_ao = Ca_subset("AO", "ALL");
        auto Cbocc_ao = Cb_subset("AO", "ALL");
        auto Fa_ao = matrix_subset_helper(Fa_, Ca_, "AO", "Fock");
        auto Fb_ao = matrix_subset_helper(Fb_, Cb_, "AO", "Fock");
        auto IFock_ao_a = linalg::triplet(Caocc_ao, Fa_ao, Caocc_ao, true, false, false);
        auto IFock_ao_b = linalg::triplet(Cbocc_ao, Fb_ao, Cbocc_ao, true, false, false);
        Precon_ao_a = std::make_shared<Matrix>("Precon", nalpha_, nmo_ - nalpha_);
        Precon_ao_b = std::make_shared<Matrix>("Precon", nbeta_, nmo_ - nbeta_);

        auto denom_ap = Precon_ao_a->pointer()[0];
        auto f_ap = IFock_ao_a->pointer();
        for (size_t i = 0, target = 0; i < nalpha_; i++) {
            for (size_t a = nalpha_; a < nmo_; a++) {
                denom_ap[target++] = -f_ap[i][i] + f_ap[a][a];
            }
        }

        auto denom_bp = Precon_ao_b->pointer()[0];
        auto f_bp = IFock_ao_b->pointer();
        for (size_t i = 0, target = 0; i < nbeta_; i++) {
            for (size_t a = nbeta_; a < nmo_; a++) {
                denom_bp[target++] = -f_bp[i][i] + f_bp[a][a];
            }
        }
    }

    if (needs_so) {
        // Grab occ and vir orbitals
        Dimension virpi_a = nmopi_ - nalphapi_;
        Dimension virpi_b = nmopi_ - nbetapi_;

        // MO Fock Matrix (Inactive Fock in Helgaker's language)
        SharedMatrix IFock_a = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
        SharedMatrix IFock_b = linalg::triplet(Cb_, Fb_, Cb_, true, false, false);
        Precon_so_a = std::make_shared<Matrix>("Alpha Precon", nirrep_, nalphapi_, virpi_a);
        Precon_so_b = std::make_shared<Matrix>("Beta Precon", nirrep_, nbetapi_, virpi_b);

        for (size_t h = 0; h < nirrep_; h++) {
            if (virpi_a[h] && nalphapi_[h]) {
                double* denom_ap = Precon_so_a->pointer(h)[0];
                double** f_ap = IFock_a->pointer(h);
                for (size_t i = 0, target = 0, max_i = nalphapi_[h], max_a = nmopi_[h]; i < max_i; i++) {
                    for (size_t a = max_i; a < max_a; a++) {
                        denom_ap[target++] = -f_ap[i][i] + f_ap[a][a];
                    }
                }
            }

            if (virpi_b[h] && nbetapi_[h]) {
                double* denom_bp = Precon_so_b->pointer(h)[0];
                double** f_bp = IFock_b->pointer(h);
                for (size_t i = 0, target = 0, max_i = nbetapi_[h], max_a = nmopi_[h]; i < max_i; i++) {
                    for (size_t a = max_i; a < max_a; a++) {
                        denom_bp[target++] = -f_bp[i][i] + f_bp[a][a];
                    }
                }
            }
        }
    }

    // => Header <= //
    if (print_lvl) {
        outfile->Printf("\n");
        outfile->Printf("   ==> Coupled-Perturbed %s Solver <==\n\n", options_.get_str("REFERENCE").c_str());
        outfile->Printf("    Maxiter             = %11d\n", max_iter);
        outfile->Printf("    Convergence         = %11.3E\n", conv_tol);
        outfile->Printf("    Number of equations = %11ld\n", x_vec.size());
        outfile->Printf("   -----------------------------------------------------\n");
        outfile->Printf("     %4s %14s %12s  %6s  %6s\n", "Iter", "Residual RMS", "Max RMS", "Remain", "Time [s]");
        outfile->Printf("   -----------------------------------------------------\n");
    }

    // => Initial state <= //

    // What vectors do we need?
    int nvecs = x_vec.size() / 2;
    int nvecs_ab = x_vec.size();
    std::vector<SharedMatrix> ret_vec, r_vec, z_vec, p_vec;
    std::vector<double> resid(nvecs), resid_denom(nvecs), rms(nvecs), rzpre(nvecs);
    std::vector<bool> active(nvecs);
    std::vector<int> active_map(nvecs);

    // => Initial CG guess <= //
    for (size_t i = 0; i < nvecs; i++) {
        ret_vec.push_back(x_vec[2 * i]->clone());
        ret_vec.push_back(x_vec[2 * i + 1]->clone());

        if (c1_input_[i]) {
            ret_vec[2 * i]->apply_denominator(Precon_ao_a);
            ret_vec[2 * i + 1]->apply_denominator(Precon_ao_b);
        } else {
            ret_vec[2 * i]->apply_denominator(Precon_so_a);
            ret_vec[2 * i + 1]->apply_denominator(Precon_so_b);
        }

        r_vec.push_back(x_vec[2 * i]->clone());
        r_vec.push_back(x_vec[2 * i + 1]->clone());
    }

    // Calc hessian vector product, find residual and conditioned residual
    std::vector<SharedMatrix> Ax_vec = cphf_Hx(ret_vec);

    double max_rms = 0.0;
    double mean_rms = 0.0;
    int nremain = 0;
    for (size_t i = 0; i < nvecs; i++) {
        r_vec[2 * i]->subtract(Ax_vec[2 * i]);
        r_vec[2 * i + 1]->subtract(Ax_vec[2 * i + 1]);

        resid[i] = r_vec[2 * i]->sum_of_squares();
        resid[i] += r_vec[2 * i + 1]->sum_of_squares();

        resid_denom[i] = x_vec[2 * i]->sum_of_squares();
        resid_denom[i] += x_vec[2 * i + 1]->sum_of_squares();

        // Compute residuals
        if (resid_denom[i] < 1.e-14) {
            resid_denom[i] = 1.e-14;  // Prevent rel denom from being too small
        }
        rms[i] = std::sqrt(resid[i] / resid_denom[i]);
        mean_rms += rms[i];
        if (rms[i] > max_rms) {
            max_rms = rms[i];
        }
        active[i] = true;
        active_map[i] = nremain;

        // p and z vectors
        z_vec.push_back(r_vec[2 * i]->clone());
        z_vec.push_back(r_vec[2 * i + 1]->clone());

        if (c1_input_[i]) {
            z_vec[2 * i]->apply_denominator(Precon_ao_a);
            z_vec[2 * i + 1]->apply_denominator(Precon_ao_b);
        } else {
            z_vec[2 * i]->apply_denominator(Precon_so_a);
            z_vec[2 * i + 1]->apply_denominator(Precon_so_b);
        }
        p_vec.push_back(z_vec[2 * i]->clone());
        p_vec.push_back(z_vec[2 * i + 1]->clone());
        nremain++;
    }
    mean_rms /= (double)nvecs;
    cphf_nfock_builds_ += nremain;

    stop = std::time(nullptr);
    if (print_lvl > 1) {
        outfile->Printf("    %5s %14.3e %12.3e %7d %9ld\n", "Guess", mean_rms, max_rms, nremain, stop - start);
    }

    // => CG iterations <= //
    for (int cg_iter = 1; cg_iter < max_iter; cg_iter++) {
        // Build an active vector
        std::vector<SharedMatrix> active_p_vec;
        for (size_t i = 0; i < nremain; i++) {
            // outfile->Printf("Giving vec %d", active_map[i]);
            active_p_vec.push_back(p_vec[2 * active_map[i]]);
            active_p_vec.push_back(p_vec[2 * active_map[i] + 1]);
        }

        // Calc hessian vector product
        std::vector<SharedMatrix> Ap_vec = cphf_Hx(active_p_vec);

        max_rms = 0.0;
        mean_rms = 0.0;
        nremain = 0;
        int new_remain = 0;
        // Find factors and scale
        for (size_t i = 0; i < nvecs; i++) {
            if (!active[i]) {
                mean_rms += rms[i];
                continue;
            }

            // Compute update
            rzpre[i] = r_vec[2 * i]->vector_dot(z_vec[2 * i]);
            rzpre[i] += r_vec[2 * i + 1]->vector_dot(z_vec[2 * i + 1]);

            double tmp_denom = p_vec[2 * i]->vector_dot(Ap_vec[2 * nremain]);
            tmp_denom += p_vec[2 * i + 1]->vector_dot(Ap_vec[2 * nremain + 1]);
            double alpha = rzpre[i] / tmp_denom;

            if (std::isnan(alpha)) {
                outfile->Printf("RHF::CPHF Warning CG alpha is zero/nan for vec %lu. Stopping vec.\n", i);
                active[i] = false;
                alpha = 0.0;
            }

            // Update vectors
            ret_vec[2 * i]->axpy(alpha, p_vec[2 * i]);
            r_vec[2 * i]->axpy(-alpha, Ap_vec[2 * nremain]);

            ret_vec[2 * i + 1]->axpy(alpha, p_vec[2 * i + 1]);
            r_vec[2 * i + 1]->axpy(-alpha, Ap_vec[2 * nremain + 1]);

            // Get residual
            resid[i] = r_vec[2 * i]->sum_of_squares();
            resid[i] += r_vec[2 * i + 1]->sum_of_squares();

            rms[i] = std::sqrt(resid[i] / resid_denom[i]);
            if (rms[i] > max_rms) {
                max_rms = rms[i];
            }
            mean_rms += rms[i];

            // Figure out what we need to do.
            if (rms[i] < conv_tol) {
                active[i] = false;
            } else {
                active_map[new_remain] = i;
                new_remain++;
            }
            nremain++;
        }
        cphf_nfock_builds_ += nremain;
        nremain = new_remain;
        mean_rms /= (double)nvecs;

        stop = std::time(nullptr);
        if (print_lvl) {
            outfile->Printf("    %5d %14.3e %12.3e %7d %9ld\n", cg_iter, mean_rms, max_rms, nremain, stop - start);
        }

        // Check convergence
        if ((max_rms < conv_tol) && (!nremain)) {
            break;
        }

        // Update p and z
        for (size_t i = 0; i < nvecs; i++) {
            if (!active[i]) continue;
            z_vec[2 * i]->copy(r_vec[2 * i]);
            z_vec[2 * i + 1]->copy(r_vec[2 * i + 1]);

            if (c1_input_[i]) {
                z_vec[2 * i]->apply_denominator(Precon_ao_a);
                z_vec[2 * i + 1]->apply_denominator(Precon_ao_b);
            } else {
                z_vec[2 * i]->apply_denominator(Precon_so_a);
                z_vec[2 * i + 1]->apply_denominator(Precon_so_b);
            }

            double tmp_numer = r_vec[2 * i]->vector_dot(z_vec[2 * i]);
            tmp_numer += r_vec[2 * i + 1]->vector_dot(z_vec[2 * i + 1]);
            double beta = tmp_numer / rzpre[i];

            p_vec[2 * i]->scale(beta);
            p_vec[2 * i]->add(z_vec[2 * i]);

            p_vec[2 * i + 1]->scale(beta);
            p_vec[2 * i + 1]->add(z_vec[2 * i + 1]);
        }
    }

    // Convergence
    if (!nremain) {
        cphf_converged_ = true;
    }

    // Print out tail
    if (print_lvl > 1) {
        outfile->Printf("   -----------------------------------------------------\n");
        outfile->Printf("\n");
        if (nremain) {
            outfile->Printf("    Warning! %d equations did not converge!\n\n", nremain);
        } else {
            outfile->Printf("    Solver has converged.\n\n");
        }
    }

    return ret_vec;
}
int UHF::soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) {
    std::time_t start, stop;
    start = std::time(nullptr);

    // => Build gradient and preconditioner <= //

    // Grab occ and vir orbitals
    auto Cocc_a = Ca_subset("SO", "OCC");
    auto Cvir_a = Ca_subset("SO", "VIR");
    auto Gradient_a = linalg::triplet(Cocc_a, Fa_, Cvir_a, true, false, false);

    auto Cocc_b = Cb_subset("SO", "OCC");
    auto Cvir_b = Cb_subset("SO", "VIR");
    auto Gradient_b = linalg::triplet(Cocc_b, Fb_, Cvir_b, true, false, false);

    // Make sure the MO gradient is reasonably small
    if ((Gradient_a->absmax() > 0.3) || (Gradient_b->absmax() > 0.3)) {
        if (print_ > 1) {
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }
    auto ret_x =
        cphf_solve({Gradient_a, Gradient_b}, soscf_conv, soscf_max_iter, soscf_print ? 2 : 0);

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, ret_x[0]);
    rotate_orbitals(Cb_, ret_x[1]);

    return cphf_nfock_builds_;
}

void UHF::compute_nos() {
    // Compute UHF NOs and NOONs [J. Chem. Phys. 88, 4926 (1988)] -- TDC, 8/15

    // Build S^1/2
    SharedMatrix SHalf = S_->clone();
    SHalf->power(0.5);

    // Diagonalize S^1/2 Dt S^1/2
    SharedMatrix SDS = factory_->create_shared_matrix("S^1/2 Dt S^1/2");
    SDS->copy(Da_);
    SDS->add(Db_);
    SDS->transform(SHalf);

    SharedMatrix UHF_NOs = factory_->create_shared_matrix("UHF NOs");
    SharedVector UHF_NOONs(factory_->create_vector());
    SDS->diagonalize(UHF_NOs, UHF_NOONs, descending);

    // Print the NOONs -- code ripped off from OEProp::compute_no_occupations()
    int max_num;
    if (options_.get_str("PRINT_NOONS") == "ALL")
        max_num = nmo_;
    else
        max_num = to_integer(options_.get_str("PRINT_NOONS"));

    std::vector<std::tuple<double, int, int> > metric;
    for (int h = 0; h < UHF_NOONs->nirrep(); h++)
        for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++)
            metric.push_back(std::tuple<double, int, int>(UHF_NOONs->get(h, i), h, i));

    std::sort(metric.begin(), metric.end(), std::greater<std::tuple<double, int, int> >());
    int offset = nalpha_;
    int start_occ = offset - max_num;
    start_occ = (start_occ < 0 ? 0 : start_occ);
    int stop_vir = offset + max_num + 1;
    stop_vir = (int)((size_t)stop_vir >= metric.size() ? metric.size() : stop_vir);
    std::vector<std::string> labels = basisset_->molecule()->irrep_labels();
    outfile->Printf("\n  UHF NO Occupations:\n");
    for (int index = start_occ; index < stop_vir; index++) {
        if (index < offset) {
            outfile->Printf("  HONO-%-2d: %4d%3s %9.7f\n", offset - index - 1, std::get<2>(metric[index]) + 1,
                            labels[std::get<1>(metric[index])].c_str(), std::get<0>(metric[index]));
        } else {
            outfile->Printf("  LUNO+%-2d: %4d%3s %9.7f\n", index - offset, std::get<2>(metric[index]) + 1,
                            labels[std::get<1>(metric[index])].c_str(), std::get<0>(metric[index]));
        }
    }
    outfile->Printf("\n");

    if (options_.get_bool("SAVE_UHF_NOS")) {
        // Save the NOs to Ca and Cb. The resulting orbitals will be restricted.

        outfile->Printf("  Saving the UHF Natural Orbitals.\n");

        SharedMatrix SHalf_inv = S_->clone();
        SHalf_inv->power(-0.5);
        Ca_->gemm(false, false, 1.0, SHalf_inv, UHF_NOs, 0.0);

        double actv_threshold = 0.02;

        // Transform the average Fock matrix to the NO basis
        SharedMatrix F_UHF_NOs = factory_->create_shared_matrix("Fock Matrix");
        F_UHF_NOs->copy(Fa_);
        F_UHF_NOs->add(Fb_);
        F_UHF_NOs->transform(Ca_);

        // Sort orbitals according to type (core,active,virtual) and energy
        std::vector<std::tuple<int, double, int, int> > sorted_nos;
        for (int h = 0; h < UHF_NOONs->nirrep(); h++) {
            for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++) {
                double noon = UHF_NOONs->get(h, i);
                int type = 0;  // core      NO >= 1.98
                if (noon < actv_threshold) {
                    type = 2;  // virtual   NO < 0.02
                } else if (noon < 2.0 - actv_threshold) {
                    type = 1;  // active    0.02 <= NO < 1.98
                }
                double epsilon = F_UHF_NOs->get(h, i, i);
                sorted_nos.push_back(std::tuple<int, double, int, int>(type, epsilon, h, i));
            }
        }
        std::sort(sorted_nos.begin(), sorted_nos.end());

        // Build the final set of UHF NOs
        std::vector<int> irrep_count(nirrep_, 0);

        for (size_t i = 0; i < sorted_nos.size(); i++) {
            int h = std::get<2>(sorted_nos[i]);
            int Ca_p = std::get<3>(sorted_nos[i]);
            int Cb_p = irrep_count[h];
            for (int mu = 0; mu < Ca_->colspi(h); mu++) {
                double value = Ca_->get(h, mu, Ca_p);
                Cb_->set(h, mu, Cb_p, value);
            }
            irrep_count[h] += 1;
        }

        // Copy sorted orbitals to Ca
        Ca_->copy(Cb_);

        // Suggest an active space
        Dimension corepi(nirrep_);
        Dimension actvpi(nirrep_);
        for (size_t i = 0; i < sorted_nos.size(); i++) {
            int type = std::get<0>(sorted_nos[i]);
            int h = std::get<2>(sorted_nos[i]);
            if (type == 0) corepi[h] += 1;
            if (type == 1) actvpi[h] += 1;
        }
        outfile->Printf("\n  Active Space from UHF-NOs (NO threshold = %.4f):\n\n", actv_threshold);

        outfile->Printf("    restricted_docc = [");
        for (int h = 0; h < nirrep_; h++) {
            outfile->Printf("%s%d", h ? "," : "", corepi[h]);
        }
        outfile->Printf("]\n");

        outfile->Printf("    active = [");
        for (int h = 0; h < nirrep_; h++) {
            outfile->Printf("%s%d", h ? "," : "", actvpi[h]);
        }
        outfile->Printf("]\n");
    }
}

std::shared_ptr<UHF> UHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
    auto hf_wfn = std::make_shared<UHF>(wfn, functional_, wfn->options(), wfn->psio());

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
    auto SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}

void UHF::setup_potential() {
    if (functional_->needs_xc()) {
        potential_ = std::make_shared<UV>(functional_, basisset_, options_);
        potential_->initialize();
    } else {
        potential_ = nullptr;
    }
}

void UHF::openorbital_scf() {
#ifndef USING_OpenOrbitalOptimizer
  throw PSIEXCEPTION("OpenOrbitalOptimizer support has not been enabled in this Psi4 build!\n");
#else
  std::function<OpenOrbitalOptimizer::FockBuilderReturn<double, double>(const OpenOrbitalOptimizer::DensityMatrix<double, double> &)> fock_builder = [&](const OpenOrbitalOptimizer::DensityMatrix<double, double> & dm) {
    // Grab the orbitals and occupations
    std::vector<arma::mat> orbitals = dm.first;
    std::vector<arma::vec> occupations = dm.second;
    assert(orbitals.size() == 2*nirrep_);
    assert(occupations.size() == 2*nirrep_);

    // Throw away zero occupations and calculate size of the matrix
    Dimension namopi(nirrep_), nbmopi(nirrep_);
    for(int h=0; h<nirrep_; h++) {
      if(nsopi_[h]==0) {
        namopi[h]=0;
        nbmopi[h]=0;
        continue;
      }

      // spin-up
      arma::uvec idx(arma::find(occupations[h]!=0.0));
      orbitals[h] = orbitals[h].cols(idx);
      occupations[h] = occupations[h](idx);
      namopi[h] = idx.n_elem;
      // This interface can't handle negative occupations
      if(idx.n_elem>0 and arma::min(occupations[h])<0.0) {
          throw PSIEXCEPTION("Negative orbital occupations not supported in Psi4 interface!\n");
      }

      // spin-down
      int hd=h+nirrep_;
      idx=(arma::find(occupations[hd]!=0.0));
      orbitals[hd] = orbitals[hd].cols(idx);
      occupations[hd] = occupations[hd](idx);
      nbmopi[h] = idx.n_elem;
      // This interface can't handle negative occupations
      if(idx.n_elem>0 and arma::min(occupations[hd])<0.0) {
          throw PSIEXCEPTION("Negative orbital occupations not supported in Psi4 interface!\n");
      }
    }

    // Form the dummy orbitals for the jk object
    auto Cadummy = std::make_shared<Matrix>("Dummy alpha orbitals", nsopi_, namopi);
    for(int h=0;h<nirrep_;h++) {
      if(namopi[h]==0)
        // Skip case of nothing to do
        continue;
      // Get the block of X
      const arma::mat Xblock(X_->to_armadillo_matrix(h));
      arma::mat Cablock = Xblock*orbitals[h]*arma::diagmat(arma::sqrt(occupations[h]));
      Cadummy->from_armadillo_matrix(Cablock,h);
    }
    auto Cbdummy = std::make_shared<Matrix>("Dummy beta orbitals", nsopi_, nbmopi);
    for(int h=0;h<nirrep_;h++) {
      if(nbmopi[h]==0)
        // Skip case of nothing to do
        continue;
      // Get the block of X
      const arma::mat Xblock(X_->to_armadillo_matrix(h));
      int hd=h+nirrep_;
      arma::mat Cbblock = Xblock*orbitals[hd]*arma::diagmat(arma::sqrt(occupations[hd]));
      Cbdummy->from_armadillo_matrix(Cbblock,h);
    }

    std::vector<SharedMatrix>& jkC = jk_->C_left();
    jkC.clear();
    jkC.push_back(Cadummy);
    jkC.push_back(Cbdummy);
    jk_->compute();
    const std::vector<SharedMatrix>& Jvec = jk_->J();
    const std::vector<SharedMatrix>& Kvec = jk_->K();
    const std::vector<SharedMatrix>& wKvec = jk_->wK();

    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
      // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
      alpha = 1.0;
      beta = 1.0;
    }
#endif

    std::vector<arma::mat> Vxca(nirrep_), Vxcb(nirrep_);
    if (functional_->needs_xc()) {
      auto Padummy = std::make_shared<Matrix>("Dummy alpha density", nsopi_, nsopi_);
      for(int h=0;h<nirrep_;h++) {
        if(namopi[h]==0)
          // Skip case of nothing to do
          continue;
        // Get the block of X
        arma::mat Cablock = Cadummy->to_armadillo_matrix(h);
        Padummy->from_armadillo_matrix(Cablock*Cablock.t(),h);
      }

      auto Pbdummy = std::make_shared<Matrix>("Dummy beta density", nsopi_, nsopi_);
      for(int h=0;h<nirrep_;h++) {
        if(nbmopi[h]==0)
          // Skip case of nothing to do
          continue;
        // Get the block of X
        const arma::mat Xblock(X_->to_armadillo_matrix(h));
        int hd=h+nirrep_;
        arma::mat Cbblock = Cbdummy->to_armadillo_matrix(h);
        Pbdummy->from_armadillo_matrix(Cbblock*Cbblock.t(),h);
      }

      potential_->set_D({Padummy, Pbdummy});
      potential_->compute_V({Va_, Vb_});
      for(int h=0;h<nirrep_;h++) {
        if(nsopi_[h]==0)
          // Skip case of nothing to do
          continue;
        Vxca[h] = Va_->to_armadillo_matrix(h);
        Vxcb[h] = Vb_->to_armadillo_matrix(h);
      }
    }

    // Build the Fock matrix and components of the total energy in each block
    std::vector<arma::mat> fock(2*nirrep_);
    double Ecore=0.0, Ecoul=0.0, Eexch=0.0;
    for(int h=0;h<nirrep_;h++) {
      if(nsopi_[h]==0)
        // Skip case of nothing to do
        continue;

      const arma::mat Xblock(X_->to_armadillo_matrix(h));
      arma::mat J_AO(Jvec[0]->to_armadillo_matrix(h));
      J_AO += Jvec[1]->to_armadillo_matrix(h);
      const arma::mat coreH(H_->to_armadillo_matrix(h));
      const arma::mat Cablock(Cadummy->to_armadillo_matrix(h));
      const arma::mat Cbblock(Cbdummy->to_armadillo_matrix(h));
      arma::mat Ka_AO;
      arma::mat Kb_AO;

      if (functional_->is_x_hybrid() && !(functional_->is_x_lrc() && jk_->get_wcombine())) {
        Ka_AO = -alpha*(Kvec[0]->to_armadillo_matrix(h));
        Kb_AO = -alpha*(Kvec[1]->to_armadillo_matrix(h));
      } else {
        Ka_AO.zeros(coreH.n_rows, coreH.n_cols);
        Kb_AO.zeros(coreH.n_rows, coreH.n_cols);
      }

      if (functional_->is_x_lrc()) {
        const arma::mat wKa_AO(wKvec[0]->to_armadillo_matrix(h));
        const arma::mat wKb_AO(wKvec[1]->to_armadillo_matrix(h));
        if (jk_->get_wcombine()) {
          Ka_AO -= wKa_AO;
          Kb_AO -= wKb_AO;
        } else {
          Ka_AO -= beta*wKa_AO;
          Kb_AO -= beta*wKb_AO;
        }
      }

      // Minus sign in K has already been taken into account above
      int hb=h+nirrep_;
      if (functional_->needs_xc()) {
        fock[h] = Xblock.t()*(coreH+J_AO+Ka_AO+Vxca[h])*Xblock;
        fock[hb] = Xblock.t()*(coreH+J_AO+Kb_AO+Vxcb[h])*Xblock;
      } else {
        fock[h] = Xblock.t()*(coreH+J_AO+Ka_AO)*Xblock;
        fock[hb] = Xblock.t()*(coreH+J_AO+Kb_AO)*Xblock;
      }

      arma::mat Pa_AO(Cablock*Cablock.t());
      arma::mat Pb_AO(Cbblock*Cbblock.t());
      arma::mat P_AO(Pa_AO+Pb_AO);
      Ecore += arma::trace(P_AO*coreH);
      Ecoul += 0.5*arma::trace(J_AO*P_AO);
      Eexch += 0.5*(arma::trace(Ka_AO*Pa_AO) + arma::trace(Kb_AO*Pb_AO));
    }
    double XC_E = 0.0;
    double VV10_E = 0.0;
    if (functional_->needs_xc()) {
        XC_E = potential_->quadrature_values()["FUNCTIONAL"];
    }
    if (functional_->needs_vv10()) {
        VV10_E = potential_->quadrature_values()["VV10"];
    }
    double Etot = Ecore+Ecoul+Eexch+nuclearrep_+XC_E+VV10_E;

    return std::make_pair(Etot,fock);
  };

  arma::uword nirrep(nirrep_);
  // Two types of particles, each of which has nirrep blocks
  arma::uvec number_of_blocks_per_particle_type({nirrep, nirrep});
  // Orbitals in each nirrep block have maximal occupation 1
  arma::vec maximum_occupation(2*nirrep,arma::fill::value(1.0));
  // Number of electrons is na+nb
  arma::vec number_of_particles({(double) nalpha_, (double) nbeta_});

  // Descriptions of the blocks: character table
  std::vector<std::string> block_descriptions(2*nirrep);
  CharacterTable ct = molecule_->point_group()->char_table();
  for(int h=0; h<nirrep_; h++) {
    int hd = h+nirrep_;
    block_descriptions[h] = std::string("alpha ") + ct.gamma(h).symbol();
    block_descriptions[hd] = std::string("beta ")  + ct.gamma(h).symbol();
  }
  double E_tol = options_.get_double("E_CONVERGENCE");
  double start_diis = options_.get_double("SCF_INITIAL_START_DIIS_TRANSITION");
  double finish_diis = options_.get_double("SCF_INITIAL_FINISH_DIIS_TRANSITION");
  int maxvecs = options_.get_double("DIIS_MAX_VECS");
  int maxiter = options_.get_int("MAXITER");
  bool fail_on_maxiter = options_.get_bool("FAIL_ON_MAXITER");

  // Get the orbital guess
  std::vector<arma::mat> orbitals(2*nirrep);
  std::vector<arma::vec> occupations(2*nirrep);
  for(int h=0;h<nirrep_;h++) {
    if(nsopi_[h]==0)
      continue;

    auto Xblock(X_->to_armadillo_matrix(h));
    auto Sblock(S_->to_armadillo_matrix(h));
    auto Cablock(Ca_->to_armadillo_matrix(h));
    auto Cbblock(Cb_->to_armadillo_matrix(h));

    if(Cablock.n_cols) {
      int hd=h+nirrep_;
      orbitals[h] = Xblock.t()*Sblock*Cablock;
      occupations[h].zeros(Cablock.n_cols);
      if(nalphapi_[h]>0)
        occupations[h].subvec(0,nalphapi_[h]-1).ones();

      orbitals[hd] = Xblock.t()*Sblock*Cbblock;
      occupations[hd].zeros(Cbblock.n_cols);
      if(nbetapi_[h]>0)
        occupations[hd].subvec(0,nbetapi_[h]-1).ones();
    }
  };

  std::function<void(const std::map<std::string,std::any> &)> callback_function = [&](const std::map<std::string,std::any> & data) {
    std::string reference = options_.get_str("REFERENCE");
    if(options_.get_str("SCF_TYPE").ends_with("DF"))
      reference = "DF-" + reference;

    int iiter = std::any_cast<size_t>(data.at("iter"));
    double E = std::any_cast<double>(data.at("E"));
    double dE = std::any_cast<double>(data.at("dE"));
    double Dnorm = std::any_cast<double>(data.at("diis_error"));
    std::string step = std::any_cast<std::string>(data.at("step"));

    outfile->Printf("   @%s iter %3i: %20.14f   %12.5e   %-11.5e %s\n", reference.c_str(), iiter, E, dE, Dnorm, step.c_str());
  };

  OpenOrbitalOptimizer::SCFSolver<double, double> scfsolver(number_of_blocks_per_particle_type, maximum_occupation, number_of_particles, fock_builder, block_descriptions);
  scfsolver.maximum_iterations(maxiter); // mod
  scfsolver.verbosity(5);  // mod
  scfsolver.convergence_threshold(E_tol);  // mod
  scfsolver.maximum_history_length(maxvecs); // mod
  scfsolver.callback_function(callback_function); // mod
  scfsolver.diis_epsilon(start_diis); // mod
  scfsolver.diis_threshold(finish_diis); // mod
  scfsolver.initialize_with_orbitals(orbitals, occupations);
  scfsolver.run();
  if(fail_on_maxiter and not scfsolver.converged())
    throw PSIEXCEPTION("SCF did not converge and FAIL_ON_MAXITER is set to true.\n");

  // Update the orbitals with OOO's solution
  auto solution = scfsolver.get_solution();
  orbitals = solution.first;
  occupations = solution.second;
  for(int h=0;h<nirrep_;h++) {
    if(nsopi_[h]==0)
      continue;
    int hd=h+nirrep_;
    const arma::mat Xblock(X_->to_armadillo_matrix(h));
    const arma::mat Sblock(S_->to_armadillo_matrix(h));
    arma::mat Cablock(Xblock*orbitals[h]);
    Ca_->from_armadillo_matrix(Cablock,h);
    arma::mat Cbblock(Xblock*orbitals[hd]);
    Cb_->from_armadillo_matrix(Cbblock,h);

    // We hope that the occupations are integer...
    arma::uvec nonzero(arma::find(occupations[h]>0.0));
    double diff(arma::norm(occupations[h](nonzero)-arma::ones<arma::vec>(nonzero.n_elem),2));
    if(diff!=0.0) {
      fprintf(stderr,"Non-integer alpha occupations in symmetry block %i\n",h);
    }
    nalphapi_[h] = (int) nonzero.n_elem;

    nonzero=arma::find(occupations[hd]>0.0);
    diff=arma::norm(occupations[hd](nonzero)-arma::ones<arma::vec>(nonzero.n_elem),2);
    if(diff!=0.0) {
      fprintf(stderr,"Non-integer beta occupations in symmetry block %i\n",h);
    }
    nbetapi_[h] = (int) nonzero.n_elem;
  }

  // Form the density matrix
  form_D();
  // Form the two-electron part
  form_G();
  // Compute the energy
  compute_E();
#endif
}

}  // namespace scf
}  // namespace psi
