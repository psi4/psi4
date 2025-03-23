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
#include "psi4/libdpd/dpd.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#ifdef USING_OpenOrbitalOptimizer
#include <openorbitaloptimizer/scfsolver.hpp>
#endif

#include "rhf.h"

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
    Dold_ = SharedMatrix(factory_->create_matrix("D old"));
    Va_ = SharedMatrix(factory_->create_matrix("V"));
    Vb_ = Va_;
    G_ = SharedMatrix(factory_->create_matrix("G"));
    J_ = SharedMatrix(factory_->create_matrix("J"));
    K_ = SharedMatrix(factory_->create_matrix("K"));
    wK_ = SharedMatrix(factory_->create_matrix("wK"));

    same_a_b_dens_ = true;
    same_a_b_orbs_ = true;

    subclass_init();
}

void RHF::finalize() {
    // Form lagrangian
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < Lagrangian_->rowdim(h); ++m) {
            for (int n = 0; n < Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i = 0; i < nalphapi_[h]; ++i) {
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

void RHF::save_density_and_energy() {
    Dold_->copy(Da_);  // Save previous density
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
    potential_->set_D({Da_});
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

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    if (functional_->is_x_hybrid() && !(functional_->is_x_lrc() && jk_->get_wcombine())) {
        G_->axpy(-alpha, K_);
    } else {
        K_->zero();
    }

    if (functional_->is_x_lrc()) {
        if (jk_->get_wcombine()) {
            G_->axpy(-1.0, wK_);
        } else {
            G_->axpy(-beta, wK_);
        }
    } else {
        wK_->zero();
    }
}

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

void RHF::form_C(double shift) {
    if (shift == 0.0) {
        diagonalize_F(Fa_, Ca_, epsilon_a_);
    } else {
        auto shifted_F = SharedMatrix(factory_->create_matrix("F"));
        auto Cvir = Ca_subset("SO", "VIR");

        auto SCvir = std::make_shared<Matrix>(nirrep_, S_->rowspi(), Cvir->colspi());
        SCvir->gemm(false, false, 1.0, S_, Cvir, 0.0);
        shifted_F->gemm(false, true, shift, SCvir, SCvir, 0.0);
        shifted_F->add(Fa_);
        diagonalize_F(shifted_F, Ca_, epsilon_a_);
    }
    find_occupation();
}

void RHF::form_D() {
    Da_->zero();

    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];

        if (nso == 0 || nmo == 0) continue;

        auto Ca = Ca_->pointer(h);
        auto D = Da_->pointer(h);

        C_DGEMM('N', 'T', nso, nso, na, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, D[0], nso);
    }

    if (debug_) {
        outfile->Printf("in RHF::form_D:\n");
        Da_->print();
    }
}

void RHF::damping_update(double damping_percentage) {
    Da_->scale(1.0 - damping_percentage);
    Da_->axpy(damping_percentage, Dold_);
}

double RHF::compute_initial_E() {
    double Etotal = nuclearrep_ + Da_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E() {
    double one_electron_E = 2.0 * Da_->vector_dot(H_);
    double kinetic_E = 2.0 * Da_->vector_dot(T_);
    double coulomb_E = 2.0 * Da_->vector_dot(J_);

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

#ifdef USING_BrianQC
    if (brianEnable and brianEnableDFT) {
        // BrianQC multiplies with the exact exchange factors inside the Fock building, so we must not do it here
        alpha = 1.0;
        beta = 1.0;
    }
#endif

    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha * Da_->vector_dot(K_);
    }
    if (functional_->is_x_lrc()) {
        if (jk_->get_do_wK() && jk_->get_wcombine()) {
            exchange_E -= Da_->vector_dot(wK_);
        } else {
            exchange_E -= beta * Da_->vector_dot(wK_);
        }
    }

    double two_electron_E = Da_->vector_dot(Fa_) - 0.5 * one_electron_E;

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = kinetic_E;
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

        auto tmp1 = linalg::triplet(Co, F, Co, true, false, false);
        auto result = linalg::doublet(tmp1, x_vec[i], false, false);

        auto tmp2 = linalg::triplet(x_vec[i], Cv, F, false, true, false);
        result->gemm(false, false, -1.0, tmp2, Cv, 1.0);

        ret.push_back(result);
    }

    return ret;
}
std::vector<SharedMatrix> RHF::twoel_Hx_full(std::vector<SharedMatrix> x_vec, bool combine, std::string return_basis, bool singlet) {
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
                throw PSIEXCEPTION("SCF::twoel_Hx incoming rotation matrices must have shape (occ x vir).");
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
            Vx.push_back(std::make_shared<Matrix>("Vx Temp", Dx[i]->rowspi(), Dx[i]->colspi(), Dx[i]->symmetry()));
        }
        potential_->compute_Vx_full(Dx, Vx, singlet);
    }

    std::vector<SharedMatrix> V_ext_pert;
    for (const auto& pert : external_cpscf_perturbations_) {
        if (print_ > 1) outfile->Printf("Adding external CPSCF contribution %s.\n", pert.first.c_str());
        for (size_t i = 0; i < x_vec.size(); i++) {
            auto Dx = linalg::doublet(Cl[i], Cr[i], false, true);
            V_ext_pert.push_back(pert.second(Dx));
        }
    }

    Cl.clear();
    Cr.clear();

    // Build return vector, ohyea thats fun
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
    if (!singlet) {
        for (size_t i = 0; i < x_vec.size(); i++) {
            J[i]->zero();
        }
    }
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
            if (V_ext_pert.size()) {
                J[i]->axpy(4.0, V_ext_pert[i]);
            }
            ret.push_back(J[i]);
        }
    } else {
        if (jk_->get_wcombine()) {
            throw PSIEXCEPTION(
                "RHF::twoel_Hx user asked for wcombine but combine==false in SCF::twoel_Hx. Please set wcombine false "
                "in your input.");
        }
        for (size_t i = 0; i < x_vec.size(); i++) {
            // always have a J-like piece (optionally include Xc)
            if (functional_->needs_xc()) {
                J[i]->add(Vx[i]);
            }
            if (V_ext_pert.size()) {
                J[i]->add(V_ext_pert[i]);
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
        auto Cocc_ao = Ca_subset("AO", "ALL");
        auto F_ao = matrix_subset_helper(Fa_, Ca_, "AO", "Fock");
        auto IFock_ao = linalg::triplet(Cocc_ao, F_ao, Cocc_ao, true, false, false);
        Precon_ao = std::make_shared<Matrix>("Precon", nalpha_, nmo_ - nalpha_);

        auto denomp = Precon_ao->pointer()[0];
        auto fp = IFock_ao->pointer();

        for (size_t i = 0, target = 0; i < nalpha_; i++) {
            for (size_t a = nalpha_; a < nmo_; a++) {
                denomp[target++] = -fp[i][i] + fp[a][a];
            }
        }
    }

    if (needs_so) {
        // MO Fock Matrix (Inactive Fock in Helgaker's language)
        auto virpi = nmopi_ - nalphapi_;
        auto IFock_so = linalg::triplet(Ca_, Fa_, Ca_, true, false, false);
        Precon_so = std::make_shared<Matrix>("Precon", nirrep_, nalphapi_, virpi);

        for (size_t h = 0; h < nirrep_; h++) {
            if (!nalphapi_[h] || !virpi[h]) continue;
            auto denomp = Precon_so->pointer(h)[0];
            auto fp = IFock_so->pointer(h);

            for (size_t i = 0, target = 0; i < nalphapi_[h]; i++) {
                for (size_t a = nalphapi_[h]; a < nmopi_[h]; a++) {
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
                int arel = aabs - Atrip.params->qoff[asym] + nalphapi_[asym];
                for (int jb = 0; jb < Atrip.params->coltot[h]; ++jb) {
                    int jabs = Atrip.params->colorb[h][jb][0];
                    int babs = Atrip.params->colorb[h][jb][1];
                    int jsym = Atrip.params->rsym[jabs];
                    int bsym = Atrip.params->ssym[babs];
                    int jrel = jabs - Atrip.params->roff[jsym];
                    int brel = babs - Atrip.params->soff[bsym] + nalphapi_[bsym];
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
            if (DSYEV_ascending(dim, Asing.matrix[h], evals, evecs) != 0){
                throw PSIEXCEPTION("DSYEV diagonalizer failed in RHF stability check!");
            }
            global_dpd_->buf4_mat_irrep_close(&Asing, h);

            int mindim = dim < 5 ? dim : 5;
            for (int i = 0; i < mindim; i++) singlet_eval_sym.push_back(std::make_pair(evals[i], h));

            zero_arr(evals, dim);
            zero_mat(evecs, dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            if (DSYEV_ascending(dim, Atrip.matrix[h], evals, evecs) != 0){
                throw PSIEXCEPTION("DSYEV diagonalizer failed in RHF stability check!");
            }
            global_dpd_->buf4_mat_irrep_close(&Atrip, h);

            for (int i = 0; i < mindim; i++) triplet_eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(evecs);
            free(evals);
        }

        outfile->Printf("    Lowest singlet (RHF->RHF) stability eigenvalues:\n");
        print_stability_analysis(singlet_eval_sym);
        outfile->Printf("    Lowest triplet (RHF->UHF) stability eigenvalues:\n");
        print_stability_analysis(triplet_eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }

    // FOLLOW is not implemented for RHF
    return false;
}

std::shared_ptr<RHF> RHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
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
    }
    if (Fa_) {
        hf_wfn->Fa_ = Fa_subset("AO");
        hf_wfn->Fb_ = hf_wfn->Fa_;
    }
    if (epsilon_a_) {
        hf_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nalphapi_, "AO", "ALL");
        hf_wfn->epsilon_b_ = hf_wfn->epsilon_a_;
    }
    // H_ ans X_ reset in the HF constructor, copy them over here
    auto SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}

void RHF::setup_potential() {
    if (functional_->needs_xc()) {
        potential_ = std::make_shared<RV>(functional_, basisset_, options_);
        potential_->initialize();
    } else {
        potential_ = nullptr;
    }
}

void RHF::openorbital_scf() {
#ifndef USING_OpenOrbitalOptimizer
  throw PSIEXCEPTION("OpenOrbitalOptimizer support has not been enabled in this Psi4 build!\n");
#else
  std::function<OpenOrbitalOptimizer::FockBuilderReturn<double, double>(const OpenOrbitalOptimizer::DensityMatrix<double, double> &)> fock_builder = [&](const OpenOrbitalOptimizer::DensityMatrix<double, double> & dm) {
    // Grab the orbitals and occupations
    std::vector<arma::mat> orbitals = dm.first;
    std::vector<arma::vec> occupations = dm.second;
    assert(orbitals.size() == nirrep_);
    assert(occupations.size() == nirrep_);

    // Throw away zero occupations and calculate size of the matrix
    Dimension nmopi(nirrep_);
    for(int h=0; h<nirrep_; h++) {
      if(nsopi_[h]==0) {
        nmopi[h]=0;
        continue;
      }
      arma::uvec idx(arma::find(occupations[h]!=0.0));
      orbitals[h] = orbitals[h].cols(idx);
      occupations[h] = occupations[h](idx);
      nmopi[h] = idx.n_elem;
      // This interface can't handle negative occupations
      if(idx.n_elem>0 and arma::min(occupations[h])<0.0) {
          throw PSIEXCEPTION("Negative orbital occupations not supported in Psi4 interface!\n");
      }
    }

    // Form the dummy orbitals for the jk object
    auto Cdummy = std::make_shared<Matrix>("Dummy orbitals", nsopi_, nmopi);
    for(int h=0;h<nirrep_;h++) {
      if(nmopi[h]==0)
        // Skip case of nothing to do
        continue;
      // Get the block of X
      const arma::mat Xblock(X_->to_armadillo_matrix(h));
      printf("h=%i: Xblock is %i x %i, orbitals is %i x %i, occupations is %i\n",h,Xblock.n_rows,Xblock.n_cols,orbitals[h].n_rows,orbitals[h].n_cols,occupations[h].n_elem); fflush(stdout);
      arma::mat Cblock = Xblock*orbitals[h]*arma::diagmat(arma::sqrt(occupations[h]));
      Cdummy->from_armadillo_matrix(Cblock,h);
    }

    std::vector<SharedMatrix>& jkC = jk_->C_left();
    jkC.clear();
    jkC.push_back(Cdummy);
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

    // Build the Fock matrix and components of the total energy in each block
    std::vector<arma::mat> fock(nirrep_);
    double Ecore=0.0, Ecoul=0.0, Eexch=0.0;
    for(int h=0;h<nirrep_;h++) {
      if(nsopi_[h]==0)
        // Skip case of nothing to do
        continue;
      const arma::mat Xblock(X_->to_armadillo_matrix(h));
      const arma::mat J_AO(Jvec[0]->to_armadillo_matrix(h));
      const arma::mat coreH(H_->to_armadillo_matrix(h));
      const arma::mat Cblock(Cdummy->to_armadillo_matrix(h));
      arma::mat K_AO(Kvec[0]->to_armadillo_matrix(h));

      if (functional_->is_x_hybrid() && !(functional_->is_x_lrc() && jk_->get_wcombine())) {
        K_AO *= -alpha;
      } else {
        K_AO.zeros();
      }

      if (functional_->is_x_lrc()) {
        const arma::mat wK_AO(wKvec[0]->to_armadillo_matrix(h));
        if (jk_->get_wcombine()) {
          K_AO -= wK_AO;
        } else {
          K_AO -= beta*wK_AO;
        }
      }

      // Minus sign in K has already been taken into account above
      fock[h] = Xblock.t()*(coreH+J_AO+0.5*K_AO)*Xblock;

      arma::mat P_AO(Cblock*Cblock.t());
      Ecore += arma::trace(P_AO*coreH);
      Ecoul += 0.5*arma::trace(J_AO*P_AO);
      Eexch += 0.25*arma::trace(K_AO*P_AO);
    }
    double Etot = Ecore+Ecoul+Eexch+nuclearrep_;
    printf("Enucrep = %.6f\nEcore = %.6f\nEcoul = %.6f\nEexch = %.6f\n",nuclearrep_,Ecore,Ecoul,Eexch);

    return std::make_pair(Etot,fock);
  };

  arma::uword nirrep(nirrep_);
  // Only one particle type, which has nirrep blocks
  arma::uvec number_of_blocks_per_particle_type({nirrep});
  // Orbitals in each nirrep block have maximal occupation 2
  arma::vec maximum_occupation(nirrep,arma::fill::value(2.0));
  // Number of electrons is na+nb
  arma::vec number_of_particles({(double) (nalpha_+nbeta_)});

  // Descriptions of the blocks: character table
  std::vector<std::string> block_descriptions(nirrep);
  CharacterTable ct = molecule_->point_group()->char_table();
  for(int h=0; h<nirrep_; h++)
    block_descriptions[h] = ct.gamma(h).symbol();

  double E_tol = options_.get_double("E_CONVERGENCE");

  // Get the orbital guess
  std::vector<arma::mat> orbitals(nirrep);
  std::vector<arma::vec> occupations(nirrep);
  for(int h=0;h<nirrep_;h++) {
    printf("irrep %i nalphapi=%i nsopi=%i\n",h,nalphapi_[h],nsopi_[h]);
    if(nsopi_[h]==0)
      continue;

    auto Xblock(X_->to_armadillo_matrix(h));
    auto Sblock(S_->to_armadillo_matrix(h));
    auto Cblock(Ca_->to_armadillo_matrix(h));

    printf("irrep %i\n",h);
    arma::mat Smo(Cblock.t()*Sblock*Cblock);
    Smo.print("Smo");

    if(Cblock.n_cols) {
      orbitals[h] = Xblock.t()*Sblock*Cblock;
      occupations[h].zeros(Cblock.n_cols);
      if(nalphapi_[h]>0)
        occupations[h].subvec(0,nalphapi_[h]-1).ones();
      occupations[h] *= 2;
    }
    Xblock.print("X");
    Sblock.print("S");
    Cblock.print("C");
    orbitals[h].print("orbitals");
    occupations[h].t().print("occupations");
  };

  OpenOrbitalOptimizer::SCFSolver<double, double> scfsolver(number_of_blocks_per_particle_type, maximum_occupation, number_of_particles, fock_builder, block_descriptions);
  scfsolver.verbosity(5);  // mod
  scfsolver.convergence_threshold(E_tol);  // mod
  scfsolver.initialize_with_orbitals(orbitals, occupations);
  scfsolver.run();

  // Update the orbitals with OOO's solution
  auto solution = scfsolver.get_solution();
  orbitals = solution.first;
  for(int h=0;h<nirrep_;h++) {
    if(nsopi_[h]==0)
      continue;
    const arma::mat Xblock(X_->to_armadillo_matrix(h));
    const arma::mat Sblock(S_->to_armadillo_matrix(h));
    arma::mat Cblock(Xblock*orbitals[h]);
    Ca_->from_armadillo_matrix(Cblock,h));
  };

#endif
}

}  // namespace scf
}  // namespace psi
