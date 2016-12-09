/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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


#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

#include "cubature.h"
#include "points.h"
#include "v.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/integral.h"

#include <cstdlib>
#include <sstream>
#include <string>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

using ULI = unsigned long int;

namespace psi {

VBase::VBase(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary,
             Options& options)
    : options_(options), primary_(primary), functional_(functional) {
    common_init();
}
VBase::~VBase() {}
void VBase::common_init() {
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    num_threads_ = 1;
#ifdef _OPENMP
    num_threads_ = omp_get_max_threads();
#endif
}
std::shared_ptr<VBase> VBase::build_V(std::shared_ptr<BasisSet> primary,
                                      std::shared_ptr<SuperFunctional> functional, Options& options,
                                      const std::string& type) {

    std::shared_ptr<VBase> v;
    if (type == "RV") {
        if (!functional->is_unpolarized()) {
            throw PSIEXCEPTION("Passed in functional was polarized for RV reference.");
        }
        v = std::shared_ptr<VBase>(new RV(functional, primary, options));
    } else if (type == "UV") {
        if (functional->is_unpolarized()) {
            throw PSIEXCEPTION("Passed in functional was unpolarized for UV reference.");
        }
        v = std::shared_ptr<VBase>(new UV(functional, primary, options));
    } else {
        throw PSIEXCEPTION("V: V type is not recognized");
    }

    return v;
}
void VBase::set_D(std::vector<SharedMatrix> Dvec) {
    if (Dvec.size() > 2) {
        throw PSIEXCEPTION("VBase::set_D: Can only set up to two D vectors.");
    }

    // Build AO2USO matrix, if needed
    if (!AO2USO_ && (Dvec[0]->nirrep() != 1)) {
        std::shared_ptr<IntegralFactory> integral(
            new IntegralFactory(primary_, primary_, primary_, primary_));
        std::shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
        AO2USO_ = SharedMatrix(pet->aotoso());
        USO2AO_ = AO2USO_->transpose();
    }

    if (AO2USO_){
        nbf_ = AO2USO_->rowspi()[0];
    } else {
        nbf_ = Dvec[0]->rowspi()[0];
    }

    // Allocate the densities
    if (D_AO_.size() != Dvec.size()) {
        D_AO_.clear();
        for (size_t i = 0; i < Dvec.size(); i++) {
            D_AO_.push_back(SharedMatrix(new Matrix("D AO temp", nbf_, nbf_)));
        }
    }

    // Copy over the AO
    for (size_t i = 0; i < Dvec.size(); i++) {
        if (Dvec[i]->nirrep() != 1) {
            D_AO_[i]->remove_symmetry(Dvec[i], USO2AO_);
        } else {
            D_AO_[i]->copy(Dvec[i]);
        }
    }
}
void VBase::initialize() {
    timer_on("V: Grid");
    grid_ = std::shared_ptr<DFTGrid>(new DFTGrid(primary_->molecule(), primary_, options_));
    timer_off("V: Grid");
}
SharedMatrix VBase::compute_gradient() {
    throw PSIEXCEPTION("VBase: gradient not implemented for this V instance.");
}
SharedMatrix VBase::compute_hessian() {
    throw PSIEXCEPTION("VBase: hessian not implemented for this V instance.");
}
void VBase::compute_V(std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("VBase: deriv not implemented for this V instance.");
}
void VBase::compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("VBase: deriv not implemented for this Vx instance.");
}
void VBase::finalize() { grid_.reset(); }
void VBase::print_header() const {
    outfile->Printf("  ==> DFT Potential <==\n\n");
    functional_->print("outfile", print_);
    grid_->print("outfile", print_);
}
std::shared_ptr<BlockOPoints> VBase::get_block(int block) { return grid_->blocks()[block]; }
size_t VBase::nblocks() { return grid_->blocks().size(); }

RV::RV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary,
       Options& options)
    : VBase(functional, primary, options) {}
RV::~RV() {}
void RV::initialize() {
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread
        std::shared_ptr<PointFunctions> point_tmp =
            std::shared_ptr<PointFunctions>(new RKSFunctions(primary_, max_points, max_functions));
        point_tmp->set_ansatz(functional_->ansatz());
        point_workers_.push_back(point_tmp);

        // Need a functional worker per thread
        functional_workers_.push_back(functional_->build_worker());
    }
}
void RV::finalize() { VBase::finalize(); }
void RV::print_header() const { VBase::print_header(); }
void RV::compute_V(std::vector<SharedMatrix> ret) {
    if ((D_AO_.size() != 1) || (ret.size() != 1)) {
        throw PSIEXCEPTION("V: RKS should have only one D/V Matrix");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0]);
    }

    // Per thread temporaries
    std::vector<SharedMatrix> V_local;
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++) {
        V_local.push_back(SharedMatrix(new Matrix("V Temp", max_functions, max_functions)));
        Q_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature Temp", max_points)));
    }

    SharedMatrix V_AO(new Matrix("V AO Temp", nbf_, nbf_));
    double** Vp = V_AO->pointer();

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);


    // Traverse the blocks of points
    #pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** V2p = V_local[rank]->pointer();
        double* QTp = Q_temp[rank]->pointer();

        // Scratch
        double** Tp = pworker->scratch()[0]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double * x = block->x();
        double * y = block->y();
        double * z = block->z();
        double * w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute Rho, Phi, etc
        // timer_on("V: Properties");
        pworker->compute_points(block);
        // timer_off("V: Properties");

        // Compute functional values
        // timer_on("V: Functional");
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        // timer_off("V: Functional");

        if (debug_ > 4) {
            block->print("outfile", debug_);
            pworker->print("outfile", debug_);
        }

        // timer_on("V: V_XC");
        double** phi = pworker->basis_value("PHI")->pointer();
        double * rho_a = pworker->point_value("RHO_A")->pointer();
        double * zk = vals["V"]->pointer();
        double * v_rho_a = vals["V_RHO_A"]->pointer();


        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints, w, 1, zk, 1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq[rank] += C_DDOT(npoints, w, 1, rho_a, 1);
        rhoaxq[rank] += C_DDOT(npoints, QTp, 1, x, 1);
        rhoayq[rank] += C_DDOT(npoints, QTp, 1, y, 1);
        rhoazq[rank] += C_DDOT(npoints, QTp, 1, z, 1);

        // => LSDA contribution (symmetrized) <= //
        // timer_on("V: LSDA");
        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(Tp[P]), '\0', nlocal * sizeof(double));
            C_DAXPY(nlocal, 0.5 * v_rho_a[P] * w[P], phi[P], 1, Tp[P], 1);
        }
        // timer_off("V: LSDA");

        // => GGA contribution (symmetrized) <= //
        if (ansatz >= 1) {
            // timer_on("V: GGA");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double * rho_ax = pworker->point_value("RHO_AX")->pointer();
            double * rho_ay = pworker->point_value("RHO_AY")->pointer();
            double * rho_az = pworker->point_value("RHO_AZ")->pointer();
            double * v_sigma_aa = vals["V_GAMMA_AA"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P]), phix[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P]), phiy[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_az[P]), phiz[P], 1, Tp[P], 1);
            }
            // timer_off("V: GGA");
        }

        // Single GEMM slams GGA+LSDA together (man but GEM's hot!)
        // timer_on("V: LSDA");
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], max_functions, Tp[0], max_functions,
                0.0, V2p[0], max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m];
            }
        }
        // timer_off("V: LSDA");

        // => Meta contribution <= //
        if (ansatz >= 2) {
            // timer_on("V: Meta");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double * v_tau_a = vals["V_TAU_A"]->pointer();

            double** phi_w[3];
            phi_w[0] = phix;
            phi_w[1] = phiy;
            phi_w[2] = phiz;

            for (int i = 0; i < 3; i++) {
                double** phiw = phi_w[i];
                for (int P = 0; P < npoints; P++) {
                    ::memset(static_cast<void*>(Tp[P]), '\0', nlocal * sizeof(double));
                    C_DAXPY(nlocal, v_tau_a[P] * w[P], phiw[P], 1, Tp[P], 1);
                }
                C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phiw[0], max_functions, Tp[0],
                        max_functions, 1.0, V2p[0], max_functions);
            }
            // timer_off("V: Meta");
        }

        // => Unpacking <= //
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                #pragma omp atomic update
                Vp[mg][ng] += V2p[ml][nl];
                #pragma omp atomic update
                Vp[ng][mg] += V2p[ml][nl];
            }
            #pragma omp atomic update
            Vp[mg][mg] += V2p[ml][ml];
        }
        // timer_off("V: V_XC");
    }

    // Set the result
    if (AO2USO_){
        ret[0]->apply_symmetry(V_AO, AO2USO_);
    } else {
        ret[0]->copy(V_AO);
    }

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"]      = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"]     = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"]     = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"]     = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"]      = quad_values_["RHO_A"];
    quad_values_["RHO_BX"]     = quad_values_["RHO_AX"];
    quad_values_["RHO_BY"]     = quad_values_["RHO_AY"];
    quad_values_["RHO_BZ"]     = quad_values_["RHO_AZ"];

    if (debug_) {
        outfile->Printf("   => Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",
                        quad_values_["RHO_AX"], quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",
                        quad_values_["RHO_BX"], quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }
}
void RV::compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    if (D_AO_.size() != 1) {
        throw PSIEXCEPTION("Vx: RKS should have only one D Matrix");
    }
    if ((Dx.size() != ret.size()) || (Dx.size() == 0)) {
        throw PSIEXCEPTION("Vx: RKS input and output sizes should be the same.");
    }

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();
    if (ansatz >= 2){
        throw PSIEXCEPTION("Vx: RKS does not support rotated V builds for MGGA's");
    }

    int old_point_deriv = point_workers_[0]->deriv();
    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set pointers to SCF density
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0]);
    }

    std::vector<SharedMatrix> Dx_vec;
    for (size_t i = 0; i < Dx.size(); i++){
        if (Dx[i]->nirrep() != 1) {
            SharedMatrix Dx_mat = SharedMatrix(new Matrix("D AO temp", nbf_, nbf_));
            Dx_mat->remove_symmetry(Dx[i], USO2AO_);
            Dx_vec.push_back(Dx_mat);
        } else {
            Dx_vec.push_back(Dx[i]);
        }
    }


    // Per [R]ank quantities
    std::vector<SharedMatrix> R_Vx_local, R_Dx_local;
    std::vector<std::shared_ptr<Vector>> R_rho_k, R_rho_k_x, R_rho_k_y, R_rho_k_z, R_gamma_k;
    for (size_t i = 0; i < num_threads_; i++) {
        R_Vx_local.push_back(SharedMatrix(new Matrix("Vx Temp", max_functions, max_functions)));
        R_Dx_local.push_back(SharedMatrix(new Matrix("Dk Temp", max_functions, max_functions)));

        R_rho_k.push_back(std::shared_ptr<Vector>(new Vector("Rho K Temp", max_points)));

        if (ansatz >= 1){
            R_rho_k_x.push_back(std::shared_ptr<Vector>(new Vector("RHO K X Temp", max_points)));
            R_rho_k_y.push_back(std::shared_ptr<Vector>(new Vector("RHO K Y Temp", max_points)));
            R_rho_k_z.push_back(std::shared_ptr<Vector>(new Vector("Rho K Z Temp", max_points)));
            R_gamma_k.push_back(std::shared_ptr<Vector>(new Vector("Gamma K Temp", max_points)));
        }

        functional_workers_[i]->set_deriv(2);
        functional_workers_[i]->allocate();
    }

    // Output quantities
    std::vector<SharedMatrix> Vx_AO;
    for (size_t i = 0; i < Dx.size(); i++){
        Vx_AO.push_back(SharedMatrix(new Matrix("Vx AO Temp", nbf_, nbf_)));
    }

    // Traverse the blocks of points
    #pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        // => Setup <= //
        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** Vx_localp = R_Vx_local[rank]->pointer();
        double** Dx_localp = R_Dx_local[rank]->pointer();

        // => Compute blocks <= //
        double** Tp = pworker->scratch()[0]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double * w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute Rho, Phi, etc
        pworker->compute_points(block);

        // Compute functional values
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);

        // => Grab quantities <= //
        // LDA
        double** phi = pworker->basis_value("PHI")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* v2_rho2 = vals["V_RHO_A_RHO_A"]->pointer();
        double* rho_k = R_rho_k[rank]->pointer();

        // GGA
        double* rho_k_x;
        double* rho_k_y;
        double* rho_k_z;
        double* gamma_k;
        double** phi_x;
        double** phi_y;
        double** phi_z;
        double * rho_x;
        double * rho_y;
        double * rho_z;
        if (ansatz >= 1) {
            rho_k_x = R_rho_k_x[rank]->pointer();
            rho_k_y = R_rho_k_y[rank]->pointer();
            rho_k_z = R_rho_k_z[rank]->pointer();
            gamma_k = R_gamma_k[rank]->pointer();
            phi_x = pworker->basis_value("PHI_X")->pointer();
            phi_y = pworker->basis_value("PHI_Y")->pointer();
            phi_z = pworker->basis_value("PHI_Z")->pointer();
            rho_x = pworker->point_value("RHO_AX")->pointer();
            rho_y = pworker->point_value("RHO_AY")->pointer();
            rho_z = pworker->point_value("RHO_AZ")->pointer();
        }

        // Meta
        // Forget that!

        // Loop over perturbation tensors
        for (size_t dindex = 0; dindex < Dx_vec.size(); dindex++) {
            double** Dxp = Dx_vec[dindex]->pointer();

            // => Build Rotated Densities <= //
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < nlocal; nl++) {
                    int ng = function_map[nl];
                    Dx_localp[ml][nl] = Dxp[mg][ng];
                }
            }

            // Rho_a = D^k_xy phi_xa phi_ya
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], max_functions, Dx_localp[0],
                    max_functions, 0.0, Tp[0], max_functions);
            C_DGEMM('N', 'T', npoints, nlocal, nlocal, 1.0, phi[0], max_functions, Dx_localp[0],
                    max_functions, 1.0, Tp[0], max_functions);
            for (int P = 0; P < npoints; P++) {
                rho_k[P] = 0.5 * C_DDOT(nlocal, phi[P], 1, Tp[P], 1);
            }
            if (ansatz >= 1) {
                for (int P = 0; P < npoints; P++) {
                    rho_k_x[P] = C_DDOT(nlocal, phi_x[P], 1, Tp[P], 1);
                    rho_k_y[P] = C_DDOT(nlocal, phi_y[P], 1, Tp[P], 1);
                    rho_k_z[P] = C_DDOT(nlocal, phi_z[P], 1, Tp[P], 1);
                    gamma_k[P] = rho_k_x[P] * rho_x[P];
                    gamma_k[P] += rho_k_y[P] * rho_y[P];
                    gamma_k[P] += rho_k_z[P] * rho_z[P];
                    gamma_k[P] *= 2;
                }
            }
            // => LSDA contribution (symmetrized) <= //
            // timer_on("V: LSDA");
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tp[P]), '\0', nlocal * sizeof(double));
                if (rho_a[P] < 1.e-6) continue;
                C_DAXPY(nlocal, 0.5 * v2_rho2[P] * w[P] * rho_k[P], phi[P], 1, Tp[P], 1);
            }

            // => GGA contribution <= //
            if (ansatz >= 1) {
                double* v_gamma = vals["V_GAMMA_AA"]->pointer();
                double* v2_gamma_gamma = vals["V_GAMMA_AA_GAMMA_AA"]->pointer();
                double* v2_rho_gamma = vals["V_RHO_A_GAMMA_AA"]->pointer();
                double tmp_val = 0.0, v2_val = 0.0;

                for (int P = 0; P < npoints; P++) {
                    if (rho_a[P] < 1.e-6) continue;

                    // V contributions
                    C_DAXPY(nlocal, (0.5 * w[P] * v2_rho_gamma[P] * gamma_k[P]), phi[P], 1, Tp[P],
                            1);

                    // W contributions
                    v2_val = (v2_rho_gamma[P] * rho_k[P] + v2_gamma_gamma[P] * gamma_k[P]);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_x[P] + v2_val * rho_x[P]);
                    C_DAXPY(nlocal, tmp_val, phi_x[P], 1, Tp[P], 1);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_y[P] + v2_val * rho_y[P]);
                    C_DAXPY(nlocal, tmp_val, phi_y[P], 1, Tp[P], 1);

                    tmp_val = 2.0 * w[P] * (v_gamma[P] * rho_k_z[P] + v2_val * rho_z[P]);
                    C_DAXPY(nlocal, tmp_val, phi_z[P], 1, Tp[P], 1);
                }
            }

            // Put it all together
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], max_functions, Tp[0],
                    max_functions, 0.0, Vx_localp[0], max_functions);

            // Symmetrization (V is Hermitian)
            for (int m = 0; m < nlocal; m++) {
                for (int n = 0; n <= m; n++) {
                    Vx_localp[m][n] = Vx_localp[n][m] = Vx_localp[m][n] + Vx_localp[n][m];
                }
            }

            // => Unpacking <= //
            double** Vxp = Vx_AO[dindex]->pointer();
            for (int ml = 0; ml < nlocal; ml++) {
                int mg = function_map[ml];
                for (int nl = 0; nl < ml; nl++) {
                    int ng = function_map[nl];
                    #pragma omp atomic update
                    Vxp[mg][ng] += Vx_localp[ml][nl];
                    #pragma omp atomic update
                    Vxp[ng][mg] += Vx_localp[ml][nl];
                }
                #pragma omp atomic update
                Vxp[mg][mg] += Vx_localp[ml][ml];
            }
            // timer_off("V: V_XC");
        }
    }

    // Set the result
    for (size_t i = 0; i < Dx.size(); i++){
        if (Dx[i]->nirrep() != 1) {
            ret[i]->apply_symmetry(Vx_AO[i], AO2USO_);
        } else {
            ret[i]->copy(Vx_AO[i]);
        }
    }

    // Reset the workers
    for (size_t i = 0; i < num_threads_; i++) {
         functional_workers_[i]->set_deriv(old_func_deriv);
         functional_workers_[i]->allocate();
    }

}
SharedMatrix RV::compute_gradient()
{
    if ((D_AO_.size() != 1))
        throw PSIEXCEPTION("V: RKS should have only one D Matrix");

    // How many atoms?
    int natom = primary_->molecule()->natom();

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0]);
        point_workers_[i]->set_deriv((functional_->is_gga() || functional_->is_meta() ? 2 : 1));
    }

    // Per thread temporaries
    std::vector<SharedMatrix> V_local, G_local;
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++){
        G_local.push_back(SharedMatrix(new Matrix("G Temp", natom, 3)));
        V_local.push_back(SharedMatrix(new Matrix("V Temp", max_functions, max_functions)));
        Q_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature Tempt", max_points)));
    }

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);

    // Traverse the blocks of points
    #pragma omp parallel for private (rank) schedule(dynamic) num_threads(num_threads_)
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** V2p = V_local[rank]->pointer();
        double** Gp = G_local[rank]->pointer();
        double* QTp = Q_temp[rank]->pointer();
        double** Dp = pworker->D_scratch()[0]->pointer();

        // Scratch
        double** Tp = pworker->scratch()[0]->pointer();
        SharedMatrix U_local(pworker->scratch()[0]->clone());
        double** Up = U_local->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // timer_on("V: Properties");
        pworker->compute_points(block);
        // timer_off("V: Properties");

        // timer_on("V: Functional");
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        // timer_off("V: Functional");

        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* zk = vals["V"]->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();

        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints,w,1,zk,1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq[rank]  += C_DDOT(npoints,w,1,rho_a,1);
        rhoaxq[rank] += C_DDOT(npoints,QTp,1,x,1);
        rhoayq[rank] += C_DDOT(npoints,QTp,1,y,1);
        rhoazq[rank] += C_DDOT(npoints,QTp,1,z,1);

        // => LSDA Contribution <= //
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tp[P], 1);
        }

        // => GGA Contribution (Term 1) <= //
        if (fworker->is_gga()) {
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P]), phi_x[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P]), phi_y[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P]), phi_z[P], 1, Tp[P], 1);
            }

        }

        // => Synthesis <= //
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,Tp[0],max_functions,Dp[0],max_functions,0.0,Up[0],max_functions);

        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_x[0][ml],max_functions);
            Gp[A][1] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_y[0][ml],max_functions);
            Gp[A][2] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_z[0][ml],max_functions);
        }

        // => GGA Contribution (Term 2) <= //
        if (fworker->is_gga()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();

            C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi[0],max_functions,Dp[0],max_functions,0.0,Up[0],max_functions);

            // x
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xx[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
            }

            // y
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
            }

            // z
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_zz[0][ml],max_functions);
            }

        }

        // => Meta Contribution <= //
        if (fworker->is_meta()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* v_tau_a = vals["V_TAU_A"]->pointer();

            double** phi_i[3];
            phi_i[0] = phi_x;
            phi_i[1] = phi_y;
            phi_i[2] = phi_z;

            double** phi_ij[3][3];
            phi_ij[0][0] = phi_xx;
            phi_ij[0][1] = phi_xy;
            phi_ij[0][2] = phi_xz;
            phi_ij[1][0] = phi_xy;
            phi_ij[1][1] = phi_yy;
            phi_ij[1][2] = phi_yz;
            phi_ij[2][0] = phi_xz;
            phi_ij[2][1] = phi_yz;
            phi_ij[2][2] = phi_zz;

            for (int i = 0; i < 3; i++) {
                double*** phi_j = phi_ij[i];
                C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi_i[i][0],max_functions,Dp[0],max_functions,0.0,Up[0],max_functions);
                for (int P = 0; P < npoints; P++) {
                    ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                    C_DAXPY(nlocal, -2.0 * w[P] * (v_tau_a[P]), Up[P], 1, Tp[P], 1);
                }
                for (int ml = 0; ml < nlocal; ml++) {
                    int A = primary_->function_to_center(function_map[ml]);
                    Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_j[0][0][ml],max_functions);
                    Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_j[1][0][ml],max_functions);
                    Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_j[2][0][ml],max_functions);
                }
            }
        }
        U_local.reset();
    }

    // Sum up the matrix
    SharedMatrix G(new Matrix("XC Gradient", natom, 3));
    for (auto const &val: G_local){
        G->add(val);
    }

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"]      = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"]     = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"]     = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"]     = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"]      = quad_values_["RHO_A"];
    quad_values_["RHO_BX"]     = quad_values_["RHO_AX"];
    quad_values_["RHO_BY"]     = quad_values_["RHO_AY"];
    quad_values_["RHO_BZ"]     = quad_values_["RHO_AZ"];

    if (debug_) {
        outfile->Printf("   => XC Gradient: Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",
                        quad_values_["RHO_AX"], quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",
                        quad_values_["RHO_BX"], quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_deriv(old_deriv);
    }

    // RKS
    G->scale(2.0);

    return G;
}


SharedMatrix RV::compute_hessian()
{
    if(functional_->is_gga() || functional_->is_meta())
        throw PSIEXCEPTION("Hessians for GGA and meta GGA functionals are not yet implemented.");

    if ((D_AO_.size() != 1))
        throw PSIEXCEPTION("V: RKS should have only one D Matrix");

    // Build the target Hessian Matrix
    int natom = primary_->molecule()->natom();
    SharedMatrix H(new Matrix("XC Hessian", 3*natom,3*natom));
    double** Hp = H->pointer();

    // Thread info
    int rank = 0;

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();
    int old_func_deriv = functional_->deriv();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    int derivlev = (functional_->is_gga() || functional_->is_meta()) ? 3 : 2;
    functional_->set_deriv(derivlev);

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0]);
        point_workers_[i]->set_deriv(derivlev);
        functional_workers_[i]->set_deriv(derivlev);
        functional_workers_[i]->allocate();
    }

    // Per thread temporaries
    std::vector<SharedMatrix> V_local;
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++){
        V_local.push_back(SharedMatrix(new Matrix("V Temp", max_functions, max_functions)));
        Q_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature Tempt", max_points)));
    }

    std::shared_ptr<Vector> QT(new Vector("Quadrature Temp", max_points));
    double* QTp = QT->pointer();
    const std::vector<std::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();

    for (size_t Q = 0; Q < blocks.size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** V2p = V_local[rank]->pointer();
        double* QTp = Q_temp[rank]->pointer();
        double** Dp = pworker->D_scratch()[0]->pointer();

        // Scratch
        double** Tp = pworker->scratch()[0]->pointer();
        SharedMatrix U_local(pworker->scratch()[0]->clone());
        double** Up = U_local->pointer();

        // ACS TODO: these need to be threaded eventually, to fit in with the new infrastructure
        SharedVector Tmpx(new Vector("Tx", max_functions));
        SharedVector Tmpy(new Vector("Ty", max_functions));
        SharedVector Tmpz(new Vector("Tz", max_functions));
        double *pTx = Tmpx->pointer();
        double *pTy = Tmpy->pointer();
        double *pTz = Tmpz->pointer();
        SharedMatrix Tx(U_local->clone());
        SharedMatrix Ty(U_local->clone());
        SharedMatrix Tz(U_local->clone());
        double **pTx2 = Tx->pointer();
        double **pTy2 = Ty->pointer();
        double **pTz2 = Tz->pointer();

        std::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        pworker->compute_points(block);
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);

        double** phi    = pworker->basis_value("PHI")->pointer();
        double** phi_x  = pworker->basis_value("PHI_X")->pointer();
        double** phi_y  = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z  = pworker->basis_value("PHI_Z")->pointer();
        double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
        double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
        double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
        double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
        double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
        double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
        double* v_rho_a  = vals["V_RHO_A"]->pointer();
        double* v_rho_aa = vals["V_RHO_A_RHO_A"]->pointer();

        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 4.0 * w[P] * v_rho_aa[P], Tp[P], 1, Up[P], 1);
        }


        // => LSDA Contribution <= //


        /*
         *                        m             n  ∂^2 F
         *  H_mn <- 4 D_ab ɸ_a ɸ_b  D_cd ɸ_c ɸ_d   ------
         *                                         ∂ ρ^2
         */

        // T = ɸ D
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi[0],max_functions,Dp[0],max_functions,0.0,Tp[0],max_functions);


        // Compute rho, to filter out small values
        for (int P = 0; P < npoints; P++) {
            double rho = C_DDOT(nlocal,phi[P],1,Tp[P],1);
            if(fabs(rho) < 1E-8){
                v_rho_a[P] = 0.0;
                v_rho_aa[P] = 0.0;
            }
//            outfile->Printf("%f %f %f\n", w[P], v_rho_a[P], v_rho_aa[P]);
        }

        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 4.0 * w[P] * v_rho_aa[P], Tp[P], 1, Up[P], 1);
        }

        for (int ml = 0; ml < nlocal; ml++) {
            pTx[ml] = C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_x[0][ml],max_functions);
            pTy[ml] = C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_y[0][ml],max_functions);
            pTz[ml] = C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_z[0][ml],max_functions);
        }

        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            double mx = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_x[0][ml],max_functions);
            double my = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_y[0][ml],max_functions);
            double mz = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_z[0][ml],max_functions);
            for (int nl = 0; nl < nlocal; nl++) {
                int B = primary_->function_to_center(function_map[nl]);
                Hp[3*A+0][3*B+0] += mx*pTx[nl];
                Hp[3*A+0][3*B+1] += mx*pTy[nl];
                Hp[3*A+0][3*B+2] += mx*pTz[nl];
                Hp[3*A+1][3*B+0] += my*pTx[nl];
                Hp[3*A+1][3*B+1] += my*pTy[nl];
                Hp[3*A+1][3*B+2] += my*pTz[nl];
                Hp[3*A+2][3*B+0] += mz*pTx[nl];
                Hp[3*A+2][3*B+1] += mz*pTy[nl];
                Hp[3*A+2][3*B+2] += mz*pTz[nl];
            }
        }

        /*
         *                        mn  ∂ F
         *  H_mn <- 2 D_ab ɸ_a ɸ_b    ---
         *                            ∂ ρ
         */
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], Tp[P], 1, Up[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            double Txx = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_xx[0][ml],max_functions);
            double Txy = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_xy[0][ml],max_functions);
            double Txz = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_xz[0][ml],max_functions);
            double Tyy = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_yy[0][ml],max_functions);
            double Tyz = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_yz[0][ml],max_functions);
            double Tzz = C_DDOT(npoints,&Up[0][ml],max_functions,&phi_zz[0][ml],max_functions);
            Hp[3*A+0][3*A+0] += Txx;
            Hp[3*A+0][3*A+1] += Txy;
            Hp[3*A+0][3*A+2] += Txz;
            Hp[3*A+1][3*A+0] += Txy;
            Hp[3*A+1][3*A+1] += Tyy;
            Hp[3*A+1][3*A+2] += Tyz;
            Hp[3*A+2][3*A+0] += Txz;
            Hp[3*A+2][3*A+1] += Tyz;
            Hp[3*A+2][3*A+2] += Tzz;
        }

        /*
         *                    m    n  ∂ F
         *  H_mn <- 2 D_ab ɸ_a  ɸ_b   ---
         *                            ∂ ρ
         */
        // T = ɸ_x D
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi_x[0],max_functions,Dp[0],max_functions,0.0,pTx2[0],max_functions);
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi_y[0],max_functions,Dp[0],max_functions,0.0,pTy2[0],max_functions);
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi_z[0],max_functions,Dp[0],max_functions,0.0,pTz2[0],max_functions);
        // x derivatives
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], pTx2[P], 1, Up[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Hp[3*A+0][3*A+0] += C_DDOT(npoints,&pTx2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+0][3*A+1] += C_DDOT(npoints,&pTy2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+0][3*A+2] += C_DDOT(npoints,&pTz2[0][ml],max_functions,&Up[0][ml],max_functions);;
        }
        // y derivatives
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], pTy2[P], 1, Up[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Hp[3*A+1][3*A+0] += C_DDOT(npoints,&pTx2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+1][3*A+1] += C_DDOT(npoints,&pTy2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+1][3*A+2] += C_DDOT(npoints,&pTz2[0][ml],max_functions,&Up[0][ml],max_functions);;
        }
        // x derivatives
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Up[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, 2.0 * w[P] * v_rho_a[P], pTz2[P], 1, Up[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Hp[3*A+2][3*A+0] += C_DDOT(npoints,&pTx2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+2][3*A+1] += C_DDOT(npoints,&pTy2[0][ml],max_functions,&Up[0][ml],max_functions);;
            Hp[3*A+2][3*A+2] += C_DDOT(npoints,&pTz2[0][ml],max_functions,&Up[0][ml],max_functions);;
        }

    }

    if (debug_) {
        outfile->Printf( "   => XC Hessian: Numerical Integrals <=\n\n");
        outfile->Printf( "    Functional Value:  %24.16E\n",quad_values_["FUNCTIONAL"]);
        outfile->Printf( "    <\\rho_a>        :  %24.16E\n",quad_values_["RHO_A"]);
        outfile->Printf( "    <\\rho_b>        :  %24.16E\n",quad_values_["RHO_B"]);
        outfile->Printf( "    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",quad_values_["RHO_AX"],quad_values_["RHO_AY"],quad_values_["RHO_AZ"]);
        outfile->Printf( "    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",quad_values_["RHO_BX"],quad_values_["RHO_BY"],quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_deriv(old_deriv);
    }
    functional_->set_deriv(old_func_deriv);

    // RKS
    H->scale(2.0);
    H->hermitivitize();

    return H;
}

UV::UV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary,
       Options& options)
    : VBase(functional, primary, options) {}
UV::~UV() {}
void UV::initialize() {
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    for (size_t i = 0; i < num_threads_; i++) {
        // Need a points worker per thread
        std::shared_ptr<PointFunctions> point_tmp =
            std::shared_ptr<PointFunctions>(new UKSFunctions(primary_, max_points, max_functions));
        point_tmp->set_ansatz(functional_->ansatz());
        point_workers_.push_back(point_tmp);

        // Need a functional worker per thread
        functional_workers_.push_back(functional_->build_worker());
    }
}
void UV::finalize()
{
    VBase::finalize();
}
void UV::print_header() const
{
    VBase::print_header();
}
void UV::compute_V(std::vector<SharedMatrix> ret)
{
    if ((D_AO_.size() != 2) || (ret.size() != 2))
        throw PSIEXCEPTION("V: UKS should have two D/V Matrices");

    // Thread info
    int rank = 0;

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0], D_AO_[1]);
    }

    // Per thread temporaries
    std::vector<SharedMatrix> Va_local, Vb_local;
    std::vector<std::shared_ptr<Vector>> Qa_temp, Qb_temp;
    for (size_t i = 0; i < num_threads_; i++){
        Va_local.push_back(SharedMatrix(new Matrix("Va Temp", max_functions, max_functions)));
        Vb_local.push_back(SharedMatrix(new Matrix("Vb Temp", max_functions, max_functions)));
        Qa_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature A Temp", max_points)));
        Qb_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature B Temp", max_points)));
    }

    SharedMatrix Va_AO(new Matrix("Va Temp", nbf_, nbf_));
    SharedMatrix Vb_AO(new Matrix("Vb Temp", nbf_, nbf_));
    double** Vap = Va_AO->pointer();
    double** Vbp = Vb_AO->pointer();

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);
    std::vector<double> rhobq(num_threads_);
    std::vector<double> rhobxq(num_threads_);
    std::vector<double> rhobyq(num_threads_);
    std::vector<double> rhobzq(num_threads_);

    // Loop over grid
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double** Va2p = Va_local[rank]->pointer();
        double** Vb2p = Vb_local[rank]->pointer();
        double* QTap = Qa_temp[rank]->pointer();
        double* QTbp = Qb_temp[rank]->pointer();

        // Scratch
        double** Tap = pworker->scratch()[0]->pointer();
        double** Tbp = pworker->scratch()[1]->pointer();

        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // timer_on("V: Properties");
        pworker->compute_points(block);
        // timer_off("V: Properties");

        // timer_on("V: Functional");
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        // timer_off("V: Functional");

        if (debug_ > 3) {
            block->print("outfile", debug_);
            pworker->print("outfile", debug_);
        }

        // timer_on("V: V_XC");
        double** phi = pworker->basis_value("PHI")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* rho_b = pworker->point_value("RHO_B")->pointer();
        double* zk = vals["V"]->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_b = vals["V_RHO_B"]->pointer();

        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints, w, 1, zk, 1);
        for (int P = 0; P < npoints; P++) {
            QTap[P] = w[P] * rho_a[P];
            QTbp[P] = w[P] * rho_b[P];
        }
        rhoaq[rank] += C_DDOT(npoints, w, 1, rho_a, 1);
        rhoaxq[rank] += C_DDOT(npoints, QTap, 1, x, 1);
        rhoayq[rank] += C_DDOT(npoints, QTap, 1, y, 1);
        rhoazq[rank] += C_DDOT(npoints, QTap, 1, z, 1);
        rhobq[rank] += C_DDOT(npoints, w, 1, rho_b, 1);
        rhobxq[rank] += C_DDOT(npoints, QTbp, 1, x, 1);
        rhobyq[rank] += C_DDOT(npoints, QTbp, 1, y, 1);
        rhobzq[rank] += C_DDOT(npoints, QTbp, 1, z, 1);

        // => LSDA contribution (symmetrized) <= //
        // timer_on("V: LSDA");
        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(Tap[P]), '\0', nlocal * sizeof(double));
            ::memset(static_cast<void*>(Tbp[P]), '\0', nlocal * sizeof(double));
            C_DAXPY(nlocal, 0.5 * v_rho_a[P] * w[P], phi[P], 1, Tap[P], 1);
            C_DAXPY(nlocal, 0.5 * v_rho_b[P] * w[P], phi[P], 1, Tbp[P], 1);
        }
        // timer_off("V: LSDA");

        // => GGA contribution (symmetrized) <= //
        if (ansatz >= 1) {
            // timer_on("V: GGA");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double * rho_ax = pworker->point_value("RHO_AX")->pointer();
            double * rho_ay = pworker->point_value("RHO_AY")->pointer();
            double * rho_az = pworker->point_value("RHO_AZ")->pointer();
            double * rho_bx = pworker->point_value("RHO_BX")->pointer();
            double * rho_by = pworker->point_value("RHO_BY")->pointer();
            double * rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double * v_sigma_aa = vals["V_GAMMA_AA"]->pointer();
            double * v_sigma_ab = vals["V_GAMMA_AB"]->pointer();
            double * v_sigma_bb = vals["V_GAMMA_BB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P] + v_sigma_ab[P] * rho_bx[P]), phix[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P] + v_sigma_ab[P] * rho_by[P]), phiy[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_az[P] + v_sigma_ab[P] * rho_bz[P]), phiz[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_bx[P] + v_sigma_ab[P] * rho_ax[P]), phix[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_by[P] + v_sigma_ab[P] * rho_ay[P]), phiy[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_bz[P] + v_sigma_ab[P] * rho_az[P]), phiz[P], 1, Tbp[P], 1);
            }
            // timer_off("V: GGA");
        }

        // timer_on("V: LSDA");
        // Single GEMM slams GGA+LSDA together (man but GEM's hot!)
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], max_functions, Tap[0],
                max_functions, 0.0, Va2p[0], max_functions);
        C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], max_functions, Tbp[0],
                max_functions, 0.0, Vb2p[0], max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                Va2p[m][n] = Va2p[n][m] = Va2p[m][n] + Va2p[n][m];
                Vb2p[m][n] = Vb2p[n][m] = Vb2p[m][n] + Vb2p[n][m];
            }
        }
        // timer_off("V: LSDA");

        // => Meta contribution <= //
        if (ansatz >= 2) {
            // timer_on("V: Meta");
            double** phix = pworker->basis_value("PHI_X")->pointer();
            double** phiy = pworker->basis_value("PHI_Y")->pointer();
            double** phiz = pworker->basis_value("PHI_Z")->pointer();
            double * v_tau_a = vals["V_TAU_A"]->pointer();
            double * v_tau_b = vals["V_TAU_B"]->pointer();

            double** phi[3];
            phi[0] = phix;
            phi[1] = phiy;
            phi[2] = phiz;

            double* v_tau[2];
            v_tau[0] = v_tau_a;
            v_tau[1] = v_tau_b;

            double** V_val[2];
            V_val[0] = Va2p;
            V_val[1] = Vb2p;

            for (int s = 0; s < 2; s++) {
                double** V2p = V_val[s];
                double* v_taup = v_tau[s];
                for (int i = 0; i < 3; i++) {
                    double** phiw = phi[i];
                    for (int P = 0; P < npoints; P++) {
                        ::memset(static_cast<void*>(Tap[P]), '\0', nlocal * sizeof(double));
                        C_DAXPY(nlocal, v_taup[P] * w[P], phiw[P], 1, Tap[P], 1);
                    }
                    C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phiw[0], max_functions, Tap[0],
                            max_functions, 1.0, V2p[0], max_functions);
                }
            }

            // timer_off("V: Meta");
        }

        // => Unpacking <= //
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                # pragma omp atomic update
                Vap[mg][ng] += Va2p[ml][nl];
                # pragma omp atomic update
                Vap[ng][mg] += Va2p[ml][nl];
                # pragma omp atomic update
                Vbp[mg][ng] += Vb2p[ml][nl];
                # pragma omp atomic update
                Vbp[ng][mg] += Vb2p[ml][nl];
            }
            # pragma omp atomic update
            Vap[mg][mg] += Va2p[ml][ml];
            # pragma omp atomic update
            Vbp[mg][mg] += Vb2p[ml][ml];
        }
        // timer_off("V: V_XC");
    }

    // Set the result
    if (AO2USO_){
        ret[0]->apply_symmetry(Va_AO, AO2USO_);
        ret[1]->apply_symmetry(Vb_AO, AO2USO_);
    } else {
        ret[0]->copy(Va_AO);
        ret[1]->copy(Vb_AO);
    }

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"]      = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"]     = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"]     = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"]     = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"]      = std::accumulate(rhobq.begin(), rhobq.end(), 0.0);
    quad_values_["RHO_BX"]     = std::accumulate(rhobxq.begin(), rhobxq.end(), 0.0);
    quad_values_["RHO_BY"]     = std::accumulate(rhobyq.begin(), rhobyq.end(), 0.0);
    quad_values_["RHO_BZ"]     = std::accumulate(rhobzq.begin(), rhobzq.end(), 0.0);

    if (debug_) {
        outfile->Printf("   => Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",
                        quad_values_["RHO_AX"], quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",
                        quad_values_["RHO_BX"], quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }
}
SharedMatrix UV::compute_gradient()
{
    if ((D_AO_.size() != 2))
        throw PSIEXCEPTION("V: UKS should have two D Matrices");

    int rank = 0;

    // Build the target gradient Matrix
    int natom = primary_->molecule()->natom();
    SharedMatrix G(new Matrix("XC Gradient", natom,3));
    double** Gp = G->pointer();

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Set Hessian derivative level in properties
    int old_deriv = point_workers_[0]->deriv();

    // Setup the pointers
    for (size_t i = 0; i < num_threads_; i++){
        point_workers_[i]->set_pointers(D_AO_[0], D_AO_[1]);
        point_workers_[i]->set_deriv((functional_->is_gga() || functional_->is_meta() ? 2 : 1));
    }

    // Thread scratch
    std::vector<std::shared_ptr<Vector>> Q_temp;
    for (size_t i = 0; i < num_threads_; i++){
        Q_temp.push_back(std::shared_ptr<Vector>(new Vector("Quadrature Temp", max_points)));
    }

    std::vector<double> functionalq(num_threads_);
    std::vector<double> rhoaq(num_threads_);
    std::vector<double> rhoaxq(num_threads_);
    std::vector<double> rhoayq(num_threads_);
    std::vector<double> rhoazq(num_threads_);
    std::vector<double> rhobq(num_threads_);
    std::vector<double> rhobxq(num_threads_);
    std::vector<double> rhobyq(num_threads_);
    std::vector<double> rhobzq(num_threads_);


    // timer_off("V: V_XC");
    for (size_t Q = 0; Q < grid_->blocks().size(); Q++) {

        // Get thread info
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = point_workers_[rank];
        double* QTp = Q_temp[rank]->pointer();

        double** Tap = pworker->scratch()[0]->pointer();
        double** Tbp = pworker->scratch()[1]->pointer();
        double** Dap = pworker->D_scratch()[0]->pointer();
        double** Dbp = pworker->D_scratch()[1]->pointer();

        SharedMatrix Ua_local(pworker->scratch()[0]->clone());
        double** Uap = Ua_local->pointer();
        SharedMatrix Ub_local(pworker->scratch()[1]->clone());
        double** Ubp = Ub_local->pointer();

        // Grid info
        std::shared_ptr<BlockOPoints> block = grid_->blocks()[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Compute grid and functional
        // timer_on("V: Properties");
        pworker->compute_points(block);
        // timer_off("V: Properties");

        // timer_on("V: Functional");
        std::map<std::string, SharedVector>& vals = fworker->compute_functional(pworker->point_values(), npoints);
        // timer_off("V: Functional");

        // More pointers
        double** phi = pworker->basis_value("PHI")->pointer();
        double** phi_x = pworker->basis_value("PHI_X")->pointer();
        double** phi_y = pworker->basis_value("PHI_Y")->pointer();
        double** phi_z = pworker->basis_value("PHI_Z")->pointer();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* rho_b = pworker->point_value("RHO_B")->pointer();
        double* zk = vals["V"]->pointer();
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_b = vals["V_RHO_B"]->pointer();

        // => Quadrature values <= //
        functionalq[rank] += C_DDOT(npoints, w, 1, zk, 1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq[rank] += C_DDOT(npoints, w, 1, rho_a, 1);
        rhoaxq[rank] += C_DDOT(npoints, QTp, 1, x, 1);
        rhoayq[rank] += C_DDOT(npoints, QTp, 1, y, 1);
        rhoazq[rank] += C_DDOT(npoints, QTp, 1, z, 1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_b[P];
        }
        rhobq[rank] += C_DDOT(npoints, w, 1, rho_b, 1);
        rhobxq[rank] += C_DDOT(npoints, QTp, 1, x, 1);
        rhobyq[rank] += C_DDOT(npoints, QTp, 1, y, 1);
        rhobzq[rank] += C_DDOT(npoints, QTp, 1, z, 1);

        // => LSDA Contribution <= //
        for (int P = 0; P < npoints; P++) {
            ::memset((void*)Tap[P], '\0', sizeof(double) * nlocal);
            ::memset((void*)Tbp[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tap[P], 1);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_b[P], phi[P], 1, Tbp[P], 1);
        }

        // => GGA Contribution (Term 1) <= //
        if (fworker->is_gga()) {
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* rho_bx = pworker->point_value("RHO_BX")->pointer();
            double* rho_by = pworker->point_value("RHO_BY")->pointer();
            double* rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]),
                        phi_x[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]),
                        phi_y[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]),
                        phi_z[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]),
                        phi_x[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]),
                        phi_y[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]),
                        phi_z[P], 1, Tbp[P], 1);
            }
        }

        // => Synthesis <= //
        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, Tap[0], max_functions, Dap[0],
                max_functions, 0.0, Uap[0], max_functions);
        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, Tbp[0], max_functions, Dbp[0],
                max_functions, 0.0, Ubp[0], max_functions);

        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_x[0][ml], max_functions);
            Gp[A][1] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_y[0][ml], max_functions);
            Gp[A][2] += C_DDOT(npoints, &Uap[0][ml], max_functions, &phi_z[0][ml], max_functions);
            Gp[A][0] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_x[0][ml], max_functions);
            Gp[A][1] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_y[0][ml], max_functions);
            Gp[A][2] += C_DDOT(npoints, &Ubp[0][ml], max_functions, &phi_z[0][ml], max_functions);
        }

        // => GGA Contribution (Term 2) <= //
        if (fworker->is_gga()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* rho_ax = pworker->point_value("RHO_AX")->pointer();
            double* rho_ay = pworker->point_value("RHO_AY")->pointer();
            double* rho_az = pworker->point_value("RHO_AZ")->pointer();
            double* rho_bx = pworker->point_value("RHO_BX")->pointer();
            double* rho_by = pworker->point_value("RHO_BY")->pointer();
            double* rho_bz = pworker->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], max_functions, Dap[0],
                    max_functions, 0.0, Uap[0], max_functions);
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], max_functions, Dbp[0],
                    max_functions, 0.0, Ubp[0], max_functions);

            // x
            for (int P = 0; P < npoints; P++) {
                ::memset((void*)Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*)Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]),
                        Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal,
                        -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]),
                        Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] +=
                    C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xx[0][ml], max_functions);
                Gp[A][1] +=
                    C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xy[0][ml], max_functions);
                Gp[A][2] +=
                    C_DDOT(npoints, &Tap[0][ml], max_functions, &phi_xz[0][ml], max_functions);
                Gp[A][0] +=
                    C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xx[0][ml], max_functions);
                Gp[A][1] +=
                    C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xy[0][ml], max_functions);
                Gp[A][2] +=
                    C_DDOT(npoints, &Tbp[0][ml], max_functions, &phi_xz[0][ml], max_functions);
            }

            // y
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]), Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]), Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][0] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
            }

            // z
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]), Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]), Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_zz[0][ml],max_functions);
                Gp[A][0] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_zz[0][ml],max_functions);
            }

        }

        // => Meta Contribution <= //
        if (fworker->is_meta()) {
            double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
            double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
            double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
            double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
            double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
            double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
            double* v_tau_a = vals["V_TAU_A"]->pointer();
            double* v_tau_b = vals["V_TAU_B"]->pointer();

            double** phi_i[3];
            phi_i[0] = phi_x;
            phi_i[1] = phi_y;
            phi_i[2] = phi_z;

            double** phi_ij[3][3];
            phi_ij[0][0] = phi_xx;
            phi_ij[0][1] = phi_xy;
            phi_ij[0][2] = phi_xz;
            phi_ij[1][0] = phi_xy;
            phi_ij[1][1] = phi_yy;
            phi_ij[1][2] = phi_yz;
            phi_ij[2][0] = phi_xz;
            phi_ij[2][1] = phi_yz;
            phi_ij[2][2] = phi_zz;

            double** Ds[2];
            Ds[0] = Dap;
            Ds[1] = Dap;

            double* v_tau_s[2];
            v_tau_s[0] = v_tau_a;
            v_tau_s[1] = v_tau_b;

            for (int s = 0; s < 2; s++) {
//                double** Dp = Ds[s];
                double* v_tau = v_tau_s[s];
                for (int i = 0; i < 3; i++) {
                    double*** phi_j = phi_ij[i];
                    C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi_i[i][0],max_functions,Dap[0],max_functions,0.0,Uap[0],max_functions);
                    for (int P = 0; P < npoints; P++) {
                        ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                        C_DAXPY(nlocal, -2.0 * w[P] * (v_tau[P]), Uap[P], 1, Tap[P], 1);
                    }
                    for (int ml = 0; ml < nlocal; ml++) {
                        int A = primary_->function_to_center(function_map[ml]);
                        Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_j[0][0][ml],max_functions);
                        Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_j[1][0][ml],max_functions);
                        Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_j[2][0][ml],max_functions);
                    }
                }
            }
        }
        Ua_local.reset();
        Ub_local.reset();

    }
    // timer_off("V: V_XC");

    quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    quad_values_["RHO_A"]      = std::accumulate(rhoaq.begin(), rhoaq.end(), 0.0);
    quad_values_["RHO_AX"]     = std::accumulate(rhoaxq.begin(), rhoaxq.end(), 0.0);
    quad_values_["RHO_AY"]     = std::accumulate(rhoayq.begin(), rhoayq.end(), 0.0);
    quad_values_["RHO_AZ"]     = std::accumulate(rhoazq.begin(), rhoazq.end(), 0.0);
    quad_values_["RHO_B"]      = std::accumulate(rhobq.begin(), rhobq.end(), 0.0);
    quad_values_["RHO_BX"]     = std::accumulate(rhobxq.begin(), rhobxq.end(), 0.0);
    quad_values_["RHO_BY"]     = std::accumulate(rhobyq.begin(), rhobyq.end(), 0.0);
    quad_values_["RHO_BZ"]     = std::accumulate(rhobzq.begin(), rhobzq.end(), 0.0);

    if (debug_) {
        outfile->Printf("   => XC Gradient: Numerical Integrals <=\n\n");
        outfile->Printf("    Functional Value:  %24.16E\n", quad_values_["FUNCTIONAL"]);
        outfile->Printf("    <\\rho_a>        :  %24.16E\n", quad_values_["RHO_A"]);
        outfile->Printf("    <\\rho_b>        :  %24.16E\n", quad_values_["RHO_B"]);
        outfile->Printf("    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",
                        quad_values_["RHO_AX"], quad_values_["RHO_AY"], quad_values_["RHO_AZ"]);
        outfile->Printf("    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",
                        quad_values_["RHO_BX"], quad_values_["RHO_BY"], quad_values_["RHO_BZ"]);
    }

    for (size_t i = 0; i < num_threads_; i++) {
        point_workers_[i]->set_deriv(old_deriv);
    }

    return G;
}

}
