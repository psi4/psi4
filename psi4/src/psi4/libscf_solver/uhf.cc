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
#include "psi4/libdiis/diisentry.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#include "stability.h"

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

    // TODO: Move that to the base object
    step_scale_ = options_.get_double("FOLLOW_STEP_SCALE");
    step_increment_ = options_.get_double("FOLLOW_STEP_INCREMENT");

    Fa_ = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_ = SharedMatrix(factory_->create_matrix("F beta"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Dt_ = SharedMatrix(factory_->create_matrix("D total"));
    Da_old_ = SharedMatrix(factory_->create_matrix("Old alpha SCF density"));
    Db_old_ = SharedMatrix(factory_->create_matrix("Old beta SCF density"));
    Dt_old_ = SharedMatrix(factory_->create_matrix("D total old"));
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
}

void UHF::finalize() {
    // Form lagrangian
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < Lagrangian_->rowdim(h); ++m) {
            for (int n = 0; n < Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i = 0; i < doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i) +
                           epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);
                }
                for (int i = doccpi_[h]; i < doccpi_[h] + soccpi_[h]; ++i)
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);

                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dt_.reset();
    Da_old_.reset();
    Db_old_.reset();
    Dt_old_.reset();
    Ga_.reset();
    Gb_.reset();

    compute_nos();

    HF::finalize();
}

void UHF::save_density_and_energy() {
    Da_old_->copy(Da_);
    Db_old_->copy(Db_);
    Dt_old_->copy(Dt_);
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

    if (alpha != 0.0) {
        Ga_->axpy(-alpha, Ka_);
        Gb_->axpy(-alpha, Kb_);
    } else {
        Ka_->zero();
        Kb_->zero();
    }

    if (functional_->is_x_lrc()) {
        Ga_->axpy(-beta, wKa_);
        Gb_->axpy(-beta, wKb_);
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

void UHF::form_C() {
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    if (options_.get_bool("GUESS_MIX") && (iteration_ == 0)) {
        if (Ca_->nirrep() == 1) {
            outfile->Printf("  Mixing alpha HOMO/LUMO orbitals (%d,%d)\n\n", nalpha_, nalpha_ + 1);
            Ca_->rotate_columns(0, nalpha_ - 1, nalpha_, pc_pi * 0.25);
            Cb_->rotate_columns(0, nbeta_ - 1, nbeta_, -pc_pi * 0.25);
        } else {
            throw InputException("Warning: cannot mix alpha HOMO/LUMO orbitals. Run in C1 symmetry.",
                                 "to 'symmetry c1'", __FILE__, __LINE__);
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

    Dt_->copy(Da_);
    Dt_->add(Db_);

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
    Dt_->copy(Da_);
    Dt_->add(Db_);
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E() {
    Dt_->copy(Da_);
    Dt_->add(Db_);
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double UHF::compute_E() {
    // E_DFT = 2.0 D*H + D*J - \alpha D*K + E_xc
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

    double exchange_E = 0.0;
    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha * Da_->vector_dot(Ka_);
        exchange_E -= alpha * Db_->vector_dot(Kb_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -= beta * Da_->vector_dot(wKa_);
        exchange_E -= beta * Db_->vector_dot(wKb_);
    }

    energies_["Nuclear"] = nuclearrep_;
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

    // Compute Fij x_ia - Fab x_ia

    SharedMatrix Ca_occ = Ca_subset("SO", "OCC");
    SharedMatrix Ca_vir = Ca_subset("SO", "VIR");
    SharedMatrix Cb_occ = Cb_subset("SO", "OCC");
    SharedMatrix Cb_vir = Cb_subset("SO", "VIR");

    std::vector<SharedMatrix> ret;
    for (size_t i = 0; i < x_vec.size() / 2; i++) {
        if ((x_vec[2 * i]->rowspi() != Ca_occ->colspi()) || (x_vec[2 * i]->colspi() != Ca_vir->colspi())) {
            throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
        }

        if ((x_vec[2 * i + 1]->rowspi() != Cb_occ->colspi()) || (x_vec[2 * i + 1]->colspi() != Cb_vir->colspi())) {
            throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
        }

        // Alpha
        SharedMatrix tmp1 = linalg::triplet(Ca_occ, Fa_, Ca_occ, true, false, false);
        SharedMatrix result = linalg::doublet(tmp1, x_vec[2 * i], false, false);

        SharedMatrix tmp2 = linalg::triplet(x_vec[2 * i], Ca_vir, Fa_, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Ca_vir, 1.0);

        ret.push_back(result);

        // Beta
        tmp1 = linalg::triplet(Cb_occ, Fb_, Cb_occ, true, false, false);
        result = linalg::doublet(tmp1, x_vec[2 * i + 1], false, false);

        tmp2 = linalg::triplet(x_vec[2 * i + 1], Cb_vir, Fb_, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Cb_vir, 1.0);

        ret.push_back(result);
    }

    return ret;
}
std::vector<SharedMatrix> UHF::twoel_Hx(std::vector<SharedMatrix> x_vec, bool combine, std::string return_basis) {
    if ((x_vec.size() % 2) != 0) {
        throw PSIEXCEPTION("UHF::onel_Hx expect incoming vector to alternate A/B");
    }

    // Compute Fij x_ia - Fab x_ia

    SharedMatrix Ca_occ = Ca_subset("SO", "OCC");
    SharedMatrix Ca_vir = Ca_subset("SO", "VIR");
    SharedMatrix Cb_occ = Cb_subset("SO", "OCC");
    SharedMatrix Cb_vir = Cb_subset("SO", "VIR");

    // Setup jk
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    int nvecs = x_vec.size() / 2;

    // We actually want to compute all alpha then all beta, its smart enough to figure out the half transform
    for (size_t i = 0; i < nvecs; i++) {
        if ((x_vec[2 * i]->rowspi() != Ca_occ->colspi()) || (x_vec[2 * i]->colspi() != Ca_vir->colspi())) {
            throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
        }

        Cl.push_back(Ca_occ);

        SharedMatrix R = linalg::doublet(Ca_vir, x_vec[2 * i], false, true);
        R->scale(-1.0);
        Cr.push_back(R);
    }
    for (size_t i = 0; i < nvecs; i++) {
        if ((x_vec[2 * i + 1]->rowspi() != Cb_occ->colspi()) || (x_vec[2 * i + 1]->colspi() != Cb_vir->colspi())) {
            throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
        }
        Cl.push_back(Cb_occ);

        SharedMatrix R = linalg::doublet(Cb_vir, x_vec[2 * i + 1], false, true);
        R->scale(-1.0);
        Cr.push_back(R);
    }

    // Compute JK
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

    Cl.clear();
    Cr.clear();

    // Build return vector
    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();
    std::vector<SharedMatrix> ret;
    if (combine) {
        for (size_t i = 0; i < nvecs; i++) {
            J[i]->scale(2.0);
            J[i]->axpy(2.0, J[nvecs + i]);
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
            ret.push_back(J[i]);
            ret.push_back(J[nvecs + i]);
        }
    } else {
        for (size_t i = 0; i < nvecs; i++) {
            J[i]->add(J[nvecs + i]);
            J[nvecs + i]->copy(J[i]);
            if (functional_->needs_xc()) {
                J[i]->add(Vx[2 * i]);
                J[nvecs + i]->add(Vx[2 * i + 1]);
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
            ret[2 * i] = linalg::triplet(Ca_occ, ret[2 * i], Ca_vir, true, false, false);
            ret[2 * i + 1] = linalg::triplet(Cb_occ, ret[2 * i + 1], Cb_vir, true, false, false);
        }
    } else {
        throw PSIEXCEPTION("SCF::twoel_Hx: return_basis option not understood.");
    }

    return ret;
}
std::vector<SharedMatrix> UHF::cphf_Hx(std::vector<SharedMatrix> x_vec) {
    // Compute quantities
    std::vector<SharedMatrix> onel = onel_Hx(x_vec);
    std::vector<SharedMatrix> twoel = twoel_Hx(x_vec, true, "MO");

    for (size_t i = 0; i < onel.size(); i++) {
        onel[i]->add(twoel[i]);
    }

    return onel;
}
std::vector<SharedMatrix> UHF::cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol, int max_iter,
                                          int print_lvl) {
    if ((x_vec.size() % 2) != 0) {
        throw PSIEXCEPTION("UHF::onel_Hx expect incoming vector to alternate A/B");
    }

    std::time_t start, stop;
    start = std::time(nullptr);
    cphf_converged_ = false;
    cphf_nfock_builds_ = 0;

    // => Build preconditioner <= //

    // Grab occ and vir orbitals
    Dimension virpi_a = nmopi_ - nalphapi_;
    Dimension virpi_b = nmopi_ - nbetapi_;

    // MO Fock Matrix (Inactive Fock in Helgaker's language)
    SharedMatrix IFock_a = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
    SharedMatrix IFock_b = linalg::triplet(Cb_, Fb_, Cb_, true, false, false);

    auto Precon_a = std::make_shared<Matrix>("Alpha Precon", nirrep_, nalphapi_, virpi_a);
    auto Precon_b = std::make_shared<Matrix>("Beta Precon", nirrep_, nbetapi_, virpi_b);

    for (size_t h = 0; h < nirrep_; h++) {
        if (virpi_a[h] && nalphapi_[h]) {
            double* denom_ap = Precon_a->pointer(h)[0];
            double** f_ap = IFock_a->pointer(h);
            for (size_t i = 0, target = 0; i < nalphapi_[h]; i++) {
                for (size_t a = nalphapi_[h]; a < nmopi_[h]; a++) {
                    denom_ap[target++] = -f_ap[i][i] + f_ap[a][a];
                }
            }
        }

        if (virpi_b[h] && nbetapi_[h]) {
            double* denom_bp = Precon_b->pointer(h)[0];
            double** f_bp = IFock_b->pointer(h);
            for (size_t i = 0, target = 0; i < nbetapi_[h]; i++) {
                for (size_t a = nbetapi_[h]; a < nmopi_[h]; a++) {
                    denom_bp[target++] = -f_bp[i][i] + f_bp[a][a];
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

        ret_vec[2 * i]->apply_denominator(Precon_a);
        ret_vec[2 * i + 1]->apply_denominator(Precon_b);

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

        z_vec[2 * i]->apply_denominator(Precon_a);
        z_vec[2 * i + 1]->apply_denominator(Precon_b);
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
            z_vec[2 * i]->apply_denominator(Precon_a);

            z_vec[2 * i + 1]->copy(r_vec[2 * i + 1]);
            z_vec[2 * i + 1]->apply_denominator(Precon_b);

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
    SharedMatrix Cocc_a = Ca_subset("SO", "OCC");
    SharedMatrix Cvir_a = Ca_subset("SO", "VIR");
    SharedMatrix Gradient_a = linalg::triplet(Cocc_a, Fa_, Cvir_a, true, false, false);

    SharedMatrix Cocc_b = Cb_subset("SO", "OCC");
    SharedMatrix Cvir_b = Cb_subset("SO", "VIR");
    SharedMatrix Gradient_b = linalg::triplet(Cocc_b, Fb_, Cvir_b, true, false, false);

    // Make sure the MO gradient is reasonably small
    if ((Gradient_a->absmax() > 0.3) || (Gradient_b->absmax() > 0.3)) {
        if (print_ > 1) {
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }
    std::vector<SharedMatrix> ret_x =
        cphf_solve({Gradient_a, Gradient_b}, soscf_conv, soscf_max_iter, soscf_print ? 2 : 0);

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, ret_x[0]);
    rotate_orbitals(Cb_, ret_x[1]);

    return cphf_nfock_builds_;
}

double UHF::compute_orbital_gradient(bool save_fock, int max_diis_vectors) {
    SharedMatrix gradient_a = form_FDSmSDF(Fa_, Da_);
    SharedMatrix gradient_b = form_FDSmSDF(Fb_, Db_);

    if (save_fock) {
        if (initialized_diis_manager_ == false) {
            diis_manager_ = std::make_shared<DIISManager>(max_diis_vectors, "HF DIIS vector", DIISManager::LargestError,
                                                          DIISManager::OnDisk);
            diis_manager_->set_error_vector_size(2, DIISEntry::Matrix, gradient_a.get(), DIISEntry::Matrix,
                                                 gradient_b.get());
            diis_manager_->set_vector_size(2, DIISEntry::Matrix, Fa_.get(), DIISEntry::Matrix, Fb_.get());
            initialized_diis_manager_ = true;
        }

        diis_manager_->add_entry(4, gradient_a.get(), gradient_b.get(), Fa_.get(), Fb_.get());
    }

    if (options_.get_bool("DIIS_RMS_ERROR")) {
        return std::sqrt(0.5 * (std::pow(gradient_a->rms(), 2) + std::pow(gradient_b->rms(), 2)));
    } else {
        return std::max(gradient_a->absmax(), gradient_b->absmax());
    }
}

bool UHF::diis() { return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get()); }

bool UHF::stability_analysis() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("Stability analysis not yet supported for XC functionals.");
    }
    auto stab = std::make_shared<UStab>(shared_from_this(), options_);
    stab->compute_energy();
    SharedMatrix eval_sym = stab->analyze();
    outfile->Printf("    Lowest UHF->UHF stability eigenvalues: \n");
    std::vector<std::pair<double, int> > eval_print;
    for (int h = 0; h < eval_sym->nirrep(); ++h) {
        for (int i = 0; i < eval_sym->rowdim(h); ++i) {
            eval_print.push_back(std::make_pair(eval_sym->get(h, i, 0), h));
        }
    }
    print_stability_analysis(eval_print);

    // And now, export the eigenvalues to a PSI4 array, mainly for testing purposes

    Process::environment.arrays["SCF STABILITY EIGENVALUES"] = eval_sym;
    if (stab->is_unstable() && options_.get_str("STABILITY_ANALYSIS") == "FOLLOW") {
        if (attempt_number_ == 1) {
            stab_val = stab->get_eigval();
        } else if (stab_val - stab->get_eigval() < 1e-4) {
            // We probably fell on the same minimum, increase step_scale_
            outfile->Printf("    Negative eigenvalue similar to previous one, wavefunction\n");
            outfile->Printf("    likely to be in the same minimum.\n");
            step_scale_ += step_increment_;
            outfile->Printf("    Modifying FOLLOW_STEP_SCALE to %f.\n", step_scale_);
        } else {
            stab_val = stab->get_eigval();
        }
        //     outfile->Printf( "OLD ORBS");
        //     Ca_->print();
        stab->rotate_orbs(step_scale_);
        //     outfile->Printf( "NEW ORBS");
        //     Ca_->print();

        // Ask politely SCF control for a new set of iterations
        return true;
    } else {
        outfile->Printf("    Stability analysis over.\n");
        // We are done, no more iterations
        return false;
    }
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
    std::shared_ptr<Wavefunction> wfn = Wavefunction::c1_deep_copy(basis);
    auto hf_wfn = std::make_shared<UHF>(wfn, functional_, wfn->options(), wfn->psio());

    // now just have to copy the matrices that UHF initializes
    // include only those that are not temporary (some deleted in finalize())
    if (Ca_) hf_wfn->Ca_ = Ca_subset("AO", "ALL");
    if (Cb_) hf_wfn->Cb_ = Cb_subset("AO", "ALL");
    if (Da_) hf_wfn->Da_ = Da_subset("AO");
    if (Db_) hf_wfn->Db_ = Db_subset("AO");
    if (Fa_) hf_wfn->Fa_ = Fa_subset("AO");
    if (Fb_) hf_wfn->Fb_ = Fb_subset("AO");
    if (epsilon_a_) hf_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nsopi_, "AO", "ALL");
    if (epsilon_b_) hf_wfn->epsilon_b_ = epsilon_subset_helper(epsilon_b_, nsopi_, "AO", "ALL");
    // H_ ans X_ reset in the HF constructor, copy them over here
    SharedMatrix SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}
}  // namespace scf
}  // namespace psi
