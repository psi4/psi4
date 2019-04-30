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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psifiles.h"
#include "psi4/physconst.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdiis/diisentry.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#include "rhf.h"

namespace psi {
namespace scf {

RHF::RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

RHF::RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
         std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

RHF::~RHF() {}

void RHF::common_init() {
    name_ = "RHF";

    if (multiplicity_ != 1) throw PSIEXCEPTION("RHF: RHF reference is only for singlets.");

    // Allocate matrix memory
    Fa_ = SharedMatrix(factory_->create_matrix("F"));
    Fb_ = Fa_;
    Ca_ = SharedMatrix(factory_->create_matrix("MO coefficients (C)"));
    Cb_ = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("orbital energies");
    epsilon_b_ = epsilon_a_;
    Da_ = SharedMatrix(factory_->create_matrix("SCF density"));
    Db_ = Da_;
    Lagrangian_ = SharedMatrix(factory_->create_matrix("X"));
    D_ = Da_;
    Dold_ = SharedMatrix(factory_->create_matrix("D old"));
    Va_ = SharedMatrix(factory_->create_matrix("V"));
    Vb_ = Va_;
    G_ = SharedMatrix(factory_->create_matrix("G"));
    J_ = SharedMatrix(factory_->create_matrix("J"));
    K_ = SharedMatrix(factory_->create_matrix("K"));
    wK_ = SharedMatrix(factory_->create_matrix("wK"));

    same_a_b_dens_ = true;
    same_a_b_orbs_ = true;
}

void RHF::finalize() {
    // Form lagrangian
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < Lagrangian_->rowdim(h); ++m) {
            for (int n = 0; n < Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i = 0; i < doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);
                }
                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dold_.reset();
    G_.reset();
    J_.reset();
    K_.reset();
    wK_.reset();

    HF::finalize();
}

SharedMatrix RHF::Da() const { return D_; }

void RHF::save_density_and_energy() {
    Dold_->copy(D_);  // Save previous density
}

void forPermutation(int depth, std::vector<int>& array, std::vector<int>& indices, int curDepth,
                    std::vector<std::vector<int> >& finalindex) {
    int length = array.size();
    if (curDepth == 0) {
        finalindex.push_back(indices);
        return;
    }
    for (int i = 0; i < length; i++) {
        bool isgood = true;
        for (int j = length - 1; j >= curDepth && isgood; j--) {
            if (indices[j] == array[i]) isgood = false;
        }
        if (isgood) {
            indices[curDepth - 1] = array[i];
            forPermutation(depth, array, indices, curDepth - 1, finalindex);
        }
    }
}
void RHF::form_V() {
    // Push the C matrix on
    // std::vector<SharedMatrix> & C = potential_->C();
    // C.clear();
    // C.push_back(Ca_subset("SO", "OCC"));

    // // Run the potential object
    // potential_->compute();

    // // Pull the V matrices off
    // const std::vector<SharedMatrix> & V = potential_->V();
    // Va_ = V[0];
    potential_->set_D({D_});
    potential_->compute_V({Va_});
    Vb_ = Va_;
}
void RHF::form_G() {
    if (functional_->needs_xc()) {
        form_V();
        G_->copy(Va_);
    } else {
        G_->zero();
    }

    /// Push the C matrix on
    std::vector<SharedMatrix>& C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();
    const std::vector<SharedMatrix>& wK = jk_->wK();
    J_ = J[0];
    if (functional_->is_x_hybrid()) {
        K_ = K[0];
    }
    if (functional_->is_x_lrc()) {
        wK_ = wK[0];
    }

    G_->axpy(2.0, J_);

    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();

    if (alpha != 0.0) {
        G_->axpy(-alpha, K_);
    } else {
        K_->zero();
    }

    if (functional_->is_x_lrc()) {
        G_->axpy(-beta, wK_);
    } else {
        wK_->zero();
    }
}

double RHF::compute_orbital_gradient(bool save_fock, int max_diis_vectors) {
    // Conventional DIIS (X'[FDS - SDF]X, where X levels things out)
    SharedMatrix gradient = form_FDSmSDF(Fa_, Da_);

    if (save_fock) {
        if (initialized_diis_manager_ == false) {
            if (scf_type_ == "DIRECT") {
                diis_manager_ = std::make_shared<DIISManager>(max_diis_vectors, "HF DIIS vector",
                                                              DIISManager::LargestError, DIISManager::InCore);
            } else {
                diis_manager_ = std::make_shared<DIISManager>(max_diis_vectors, "HF DIIS vector",
                                                              DIISManager::LargestError, DIISManager::OnDisk);
            }
            diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, gradient.get());
            diis_manager_->set_vector_size(1, DIISEntry::Matrix, Fa_.get());
            initialized_diis_manager_ = true;
        }
        diis_manager_->add_entry(2, gradient.get(), Fa_.get());
    }

    if (options_.get_bool("DIIS_RMS_ERROR")) {
        return gradient->rms();
    } else {
        return gradient->absmax();
    }
}

bool RHF::diis() { return diis_manager_->extrapolate(1, Fa_.get()); }

void RHF::form_F() {
    Fa_->copy(H_);
    Fa_->add(G_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }

    if (debug_) {
        Fa_->print();
        J_->print();
        K_->print();
        if (functional_->needs_xc()) {
            Va_->print();
        }
        G_->print();
    }
}

void RHF::form_C() {
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    find_occupation();
}

void RHF::form_D() {
    D_->zero();

    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = doccpi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** D = D_->pointer(h);

        C_DGEMM('N', 'T', nso, nso, na, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, D[0], nso);
    }

    if (debug_) {
        outfile->Printf("in RHF::form_D:\n");
        D_->print();
    }
}

void RHF::damping_update(double damping_percentage) {
    D_->scale(1.0 - damping_percentage);
    D_->axpy(damping_percentage, Dold_);
}

double RHF::compute_initial_E() {
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E() {
    double one_electron_E = 2.0 * D_->vector_dot(H_);
    double coulomb_E = 2.0 * D_->vector_dot(J_);

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
        exchange_E -= alpha * Da_->vector_dot(K_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -= beta * Da_->vector_dot(wK_);
    }

    double two_electron_E = D_->vector_dot(Fa_) - 0.5 * one_electron_E;

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = coulomb_E + exchange_E;
    energies_["XC"] = XC_E;
    energies_["VV10"] = VV10_E;
    energies_["-D"] = scalar_variable("-D Energy");
    double dashD_E = energies_["-D"];

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += coulomb_E;
    Etotal += exchange_E;
    Etotal += XC_E;
    Etotal += VV10_E;
    Etotal += dashD_E;

    return Etotal;
}
std::vector<SharedMatrix> RHF::onel_Hx(std::vector<SharedMatrix> x_vec) {
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

    SharedMatrix Cocc_ao, Cvir_ao, F_ao, Cocc_so, Cvir_so;
    if (needs_ao) {
        Cocc_ao = Ca_subset("AO", "OCC");
        Cvir_ao = Ca_subset("AO", "VIR");
        F_ao = matrix_subset_helper(Fa_, Ca_, "AO", "Fock");
    }
    if (needs_so) {
        Cocc_so = Ca_subset("SO", "OCC");
        Cvir_so = Ca_subset("SO", "VIR");
    }

    // Compute Fij x_ia - Fab x_ia
    std::vector<SharedMatrix> ret;
    SharedMatrix F, Co, Cv;
    for (size_t i = 0; i < x_vec.size(); i++) {
        if (c1_input_[i]) {
            if ((x_vec[i]->rowspi()[0] != nalpha_) || (x_vec[i]->colspi()[0] != (nmo_ - nalpha_))) {
                throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            F = F_ao;
            Co = Cocc_ao;
            Cv = Cvir_ao;

        } else {
            if ((x_vec[i]->rowspi() != Cocc_so->colspi()) || (x_vec[i]->colspi() != Cvir_so->colspi())) {
                throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            F = Fa_;
            Co = Cocc_so;
            Cv = Cvir_so;
        }

        SharedMatrix tmp1 = linalg::triplet(Co, F, Co, true, false, false);
        SharedMatrix result = linalg::doublet(tmp1, x_vec[i], false, false);

        SharedMatrix tmp2 = linalg::triplet(x_vec[i], Cv, F, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Cv, 1.0);

        ret.push_back(result);
    }

    return ret;
}
std::vector<SharedMatrix> RHF::twoel_Hx(std::vector<SharedMatrix> x_vec, bool combine, std::string return_basis) {
    // Make sure we have a JK object
    if (!jk_) {
        throw PSIEXCEPTION("RHF::twoel_Hx: JK object is not initialized, please set option SAVE_JK to True.");
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

    SharedMatrix Cocc_ao, Cvir_ao, Cocc_so, Cvir_so;
    if (needs_ao) {
        Cocc_ao = Ca_subset("AO", "OCC");
        Cvir_ao = Ca_subset("AO", "VIR");
    }
    if (needs_so) {
        Cocc_so = Ca_subset("SO", "OCC");
        Cvir_so = Ca_subset("SO", "VIR");
    }

    // Compute Ipqrs,K_rs -> J and IprqsK_rs -> K
    // K = Co x Cv.T
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    SharedMatrix Co, Cv;
    for (size_t i = 0; i < x_vec.size(); i++) {
        if (c1_input_[i]) {
            if ((x_vec[i]->rowspi()[0] != nalpha_) || (x_vec[i]->colspi()[0] != (nmo_ - nalpha_))) {
                throw PSIEXCEPTION("SCF::onel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Co = Cocc_ao;
            Cv = Cvir_ao;
        } else {
            if ((x_vec[i]->rowspi() != Cocc_so->colspi()) || (x_vec[i]->colspi() != Cvir_so->colspi())) {
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
            }
            Co = Cocc_so;
            Cv = Cvir_so;
        }

        Cl.push_back(Co);
        SharedMatrix R = linalg::doublet(Cv, x_vec[i], false, true);
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
        for (size_t i = 0; i < x_vec.size(); i++) {
            Dx.push_back(linalg::doublet(Cl[i], Cr[i], false, true));
            Vx.push_back(std::make_shared<Matrix>("Vx Temp", Dx[i]->rowspi(), Dx[i]->colspi()));
        }
        potential_->compute_Vx(Dx, Vx);
    }

    Cl.clear();
    Cr.clear();

    // Build return vector, ohyea thats fun
    double alpha = functional_->x_alpha();
    double beta = functional_->x_beta();
    std::vector<SharedMatrix> ret;
    if (combine) {
        // Cocc_ni (4 * J[D]_nm - K[D]_nm - K[D]_mn) C_vir_ma
        for (size_t i = 0; i < x_vec.size(); i++) {
            J[i]->scale(4.0);
            if (functional_->is_x_hybrid()) {
                J[i]->axpy(-alpha, K[i]);
                J[i]->axpy(-alpha, K[i]->transpose());
            }
            if (functional_->is_x_lrc()) {
                J[i]->axpy(-beta, wK[i]);
                J[i]->axpy(-beta, wK[i]->transpose());
            }
            if (functional_->needs_xc()) {
                J[i]->axpy(4.0, Vx[i]);
            }
            ret.push_back(J[i]);
        }
    } else {
        for (size_t i = 0; i < x_vec.size(); i++) {
            // always have a J-like piece (optionally include Xc)
            if (functional_->needs_xc()) {
                J[i]->add(Vx[i]);
            }
            ret.push_back(J[i]);
            // may have K^HF
            if (functional_->is_x_hybrid()) {
                K[i]->scale(alpha);
                // may have K^HF + wK
                if (functional_->is_x_lrc()) {
                    K[i]->axpy(beta, wK[i]);
                }
                ret.push_back(K[i]);
                // could also have just wK
            } else if (functional_->is_x_lrc()) {
                wK[i]->scale(beta);
                ret.push_back(wK[i]);
            }  // also may have not k-like piece
        }
    }

    // Transform if needed
    if (return_basis == "SO") {
        /* pass */
    } else if (return_basis == "MO") {
        for (size_t i = 0; i < ret.size(); i++) {
            if (c1_input_[i]) {
                ret[i] = linalg::triplet(Cocc_ao, ret[i], Cvir_ao, true, false, false);
            } else {
                ret[i] = linalg::triplet(Cocc_so, ret[i], Cvir_so, true, false, false);
            }
        }
    } else {
        throw PSIEXCEPTION("SCF::twoel_Hx: return_basis option not understood.");
    }

    return ret;
}
std::vector<SharedMatrix> RHF::cphf_Hx(std::vector<SharedMatrix> x_vec) {
    // Compute quantities
    std::vector<SharedMatrix> onel = onel_Hx(x_vec);
    std::vector<SharedMatrix> twoel = twoel_Hx(x_vec, true, "MO");

    for (size_t i = 0; i < onel.size(); i++) {
        onel[i]->add(twoel[i]);
    }

    return onel;
}
std::vector<SharedMatrix> RHF::cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol, int max_iter,
                                          int print_lvl) {
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
    SharedMatrix Precon_ao, Precon_so;

    if (needs_ao) {
        // MO (C1) Fock Matrix (Inactive Fock in Helgaker's language)
        SharedMatrix Cocc_ao = Ca_subset("AO", "ALL");
        SharedMatrix F_ao = matrix_subset_helper(Fa_, Ca_, "AO", "Fock");
        SharedMatrix IFock_ao = linalg::triplet(Cocc_ao, F_ao, Cocc_ao, true, false, false);
        Precon_ao = std::make_shared<Matrix>("Precon", nalpha_, nmo_ - nalpha_);

        double* denomp = Precon_ao->pointer()[0];
        double** fp = IFock_ao->pointer();

        for (size_t i = 0, target = 0; i < nalpha_; i++) {
            for (size_t a = nalpha_; a < nmo_; a++) {
                denomp[target++] = -fp[i][i] + fp[a][a];
            }
        }
    }

    if (needs_so) {
        // MO Fock Matrix (Inactive Fock in Helgaker's language)
        Dimension virpi = nmopi_ - nalphapi_;
        SharedMatrix IFock_so = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
        Precon_so = std::make_shared<Matrix>("Precon", nirrep_, doccpi_, virpi);

        for (size_t h = 0; h < nirrep_; h++) {
            if (!doccpi_[h] || !virpi[h]) continue;
            double* denomp = Precon_so->pointer(h)[0];
            double** fp = IFock_so->pointer(h);

            for (size_t i = 0, target = 0; i < doccpi_[h]; i++) {
                for (size_t a = doccpi_[h]; a < nmopi_[h]; a++) {
                    denomp[target++] = -fp[i][i] + fp[a][a];
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
    int nvecs = x_vec.size();
    std::vector<SharedMatrix> ret_vec, r_vec, z_vec, p_vec;
    std::vector<double> resid(nvecs), resid_denom(nvecs), rms(nvecs), rzpre(nvecs);
    std::vector<bool> active(nvecs);
    std::vector<int> active_map(nvecs);

    // => Initial CG guess <= //
    for (size_t i = 0; i < nvecs; i++) {
        ret_vec.push_back(x_vec[i]->clone());
        if (c1_input_[i]) {
            ret_vec[i]->apply_denominator(Precon_ao);
        } else {
            ret_vec[i]->apply_denominator(Precon_so);
        }
        r_vec.push_back(x_vec[i]->clone());
    }

    // Calc hessian vector product, find residual and conditioned residual
    std::vector<SharedMatrix> Ax_vec = cphf_Hx(ret_vec);

    double max_rms = 0.0;
    double mean_rms = 0.0;
    int nremain = 0;
    for (size_t i = 0; i < nvecs; i++) {
        r_vec[i]->subtract(Ax_vec[i]);
        resid[i] = r_vec[i]->sum_of_squares();
        resid_denom[i] = x_vec[i]->sum_of_squares();

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
        z_vec.push_back(r_vec[i]->clone());
        if (c1_input_[i]) {
            z_vec[i]->apply_denominator(Precon_ao);
        } else {
            z_vec[i]->apply_denominator(Precon_so);
        }
        p_vec.push_back(z_vec[i]->clone());
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
            active_p_vec.push_back(p_vec[active_map[i]]);
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
            rzpre[i] = r_vec[i]->vector_dot(z_vec[i]);
            double alpha = rzpre[i] / p_vec[i]->vector_dot(Ap_vec[nremain]);

            if (std::isnan(alpha)) {
                outfile->Printf("RHF::CPHF Warning CG alpha is zero/nan for vec %lu. Stopping vec.\n", i);
                active[i] = false;
                alpha = 0.0;
            }

            // Update vectors
            ret_vec[i]->axpy(alpha, p_vec[i]);
            r_vec[i]->axpy(-alpha, Ap_vec[nremain]);

            // Get residual
            resid[i] = r_vec[i]->sum_of_squares();
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
            z_vec[i]->copy(r_vec[i]);
            // z_vec[i]->apply_denominator(Precon);
            if (c1_input_[i]) {
                z_vec[i]->apply_denominator(Precon_ao);
            } else {
                z_vec[i]->apply_denominator(Precon_so);
            }

            double beta = r_vec[i]->vector_dot(z_vec[i]) / rzpre[i];

            p_vec[i]->scale(beta);
            p_vec[i]->add(z_vec[i]);
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

    // => Cleanup <= //
    Precon_ao.reset();
    Precon_so.reset();

    return ret_vec;
}

int RHF::soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print) {
    int fock_builds;
    std::time_t start, stop;
    start = std::time(nullptr);

    // => Build gradient and preconditioner <= //

    // Grab occ and vir orbitals
    SharedMatrix Cocc = Ca_subset("SO", "OCC");
    SharedMatrix Cvir = Ca_subset("SO", "VIR");

    // Gradient RHS
    SharedMatrix Gradient = linalg::triplet(Cocc, Fa_, Cvir, true, false, false);

    // Make sure the MO gradient is reasonably small
    if (Gradient->absmax() > 0.3) {
        if (print_ > 1) {
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }

    std::vector<SharedMatrix> ret_x = cphf_solve({Gradient}, soscf_conv, soscf_max_iter, soscf_print ? 2 : 0);

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, ret_x[0]);

    return cphf_nfock_builds_;
}

bool RHF::stability_analysis() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("Stability analysis not yet supported for XC functionals.");
    }
    if (scf_type_ == "DF" || scf_type_ == "CD") {
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    } else {
#define ID(x) ints.DPD_ID(x)
        // Build the Fock Matrix
        auto moF = std::make_shared<Matrix>("MO basis fock matrix", nmopi_, nmopi_);
        moF->transform(Fa_, Ca_);

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        IntegralTransform ints(shared_from_this(), spaces, IntegralTransform::TransformationType::Restricted,
                               IntegralTransform::OutputType::DPDOnly, IntegralTransform::MOOrdering::QTOrder,
                               IntegralTransform::FrozenOrbitals::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Asing, Atrip, I;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "MO Ints (OV|OV)");
        // Singlet A_ia_jb = 4 (ia|jb)
        global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "RHF Singlet Hessian (IA|JB)", 4.0);
        // Triplet A_ia_jb = -(ib|ja)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"),
                                    "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0,
                               "MO Ints (OO|VV)");
        // Triplet A_ia_jb -= (ij|ab)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"),
                                    "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "RHF Triplet Hessian (IA|JB)");
        for (int h = 0; h < Atrip.params->nirreps; ++h) {
            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            for (int ia = 0; ia < Atrip.params->rowtot[h]; ++ia) {
                int iabs = Atrip.params->roworb[h][ia][0];
                int aabs = Atrip.params->roworb[h][ia][1];
                int isym = Atrip.params->psym[iabs];
                int asym = Atrip.params->qsym[aabs];
                int irel = iabs - Atrip.params->poff[isym];
                int arel = aabs - Atrip.params->qoff[asym] + doccpi_[asym];
                for (int jb = 0; jb < Atrip.params->coltot[h]; ++jb) {
                    int jabs = Atrip.params->colorb[h][jb][0];
                    int babs = Atrip.params->colorb[h][jb][1];
                    int jsym = Atrip.params->rsym[jabs];
                    int bsym = Atrip.params->ssym[babs];
                    int jrel = jabs - Atrip.params->roff[jsym];
                    int brel = babs - Atrip.params->soff[bsym] + doccpi_[bsym];
                    // Triplet A_ia_jb += delta_ij F_ab - delta_ab F_ij
                    if ((iabs == jabs) && (asym == bsym)) Atrip.matrix[h][ia][jb] += moF->get(asym, arel, brel);
                    if ((aabs == babs) && (isym == jsym)) Atrip.matrix[h][ia][jb] -= moF->get(isym, irel, jrel);
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Atrip, h);
        }
        // Singlet A += Triplet A
        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_axpy(&Atrip, &Asing, 1.0);
        global_dpd_->buf4_close(&Atrip);
        global_dpd_->buf4_close(&Asing);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> > singlet_eval_sym;
        std::vector<std::pair<double, int> > triplet_eval_sym;

        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "RHF Triplet Hessian (IA|JB)");
        for (int h = 0; h < Asing.params->nirreps; ++h) {
            int dim = Asing.params->rowtot[h];
            if (dim == 0) continue;
            double* evals = init_array(dim);
            double** evecs = block_matrix(dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Asing, h);
            global_dpd_->buf4_mat_irrep_rd(&Asing, h);
            sq_rsp(dim, dim, Asing.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Asing, h);

            int mindim = dim < 5 ? dim : 5;
            for (int i = 0; i < mindim; i++) singlet_eval_sym.push_back(std::make_pair(evals[i], h));

            zero_arr(evals, dim);
            zero_mat(evecs, dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            sq_rsp(dim, dim, Atrip.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Atrip, h);

            for (int i = 0; i < mindim; i++) triplet_eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(evecs);
            free(evals);
        }

        outfile->Printf("    Lowest singlet (RHF->RHF) stability eigenvalues:-\n");
        print_stability_analysis(singlet_eval_sym);
        outfile->Printf("    Lowest triplet (RHF->UHF) stability eigenvalues:-\n");
        print_stability_analysis(triplet_eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }

    // FOLLOW is not implemented for RHF
    return false;
}

std::shared_ptr<RHF> RHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    std::shared_ptr<Wavefunction> wfn = Wavefunction::c1_deep_copy(basis);
    auto hf_wfn = std::make_shared<RHF>(wfn, functional_, wfn->options(), wfn->psio());
    // now just have to copy the matrices that RHF initializes
    // include only those that are not temporary (some deleted in finalize())
    if (Ca_) {
        hf_wfn->Ca_ = Ca_subset("AO", "ALL");
        hf_wfn->Cb_ = hf_wfn->Ca_;
    }
    if (Da_) {
        hf_wfn->Da_ = Da_subset("AO");
        hf_wfn->Db_ = hf_wfn->Da_;
        hf_wfn->D_ = hf_wfn->Da_;
    }
    if (Fa_) {
        hf_wfn->Fa_ = Fa_subset("AO");
        hf_wfn->Fb_ = hf_wfn->Fa_;
    }
    if (epsilon_a_) {
        hf_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nsopi_, "AO", "ALL");
        hf_wfn->epsilon_b_ = hf_wfn->epsilon_a_;
    }
    // H_ ans X_ reset in the HF constructor, copy them over here
    SharedMatrix SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}
}  // namespace scf
}  // namespace psi
