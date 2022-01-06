/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "points.h"
#include "cubature.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include "gau2grid/gau2grid.h"

#include <cmath>

namespace psi {

SAPFunctions::SAPFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions)
    : PointFunctions(primary, max_points, max_functions) {
    current_basis_map_ = &basis_values_;
}
SAPFunctions::~SAPFunctions() {}
std::vector<SharedMatrix> SAPFunctions::scratch() {
    std::vector<SharedMatrix> vec;
    vec.push_back(temp_);
    return vec;
}
void SAPFunctions::build_temps() { temp_ = std::make_shared<Matrix>("Temp", max_points_, max_functions_); }
void SAPFunctions::allocate() {
    BasisFunctions::allocate();
    point_values_.clear();
    build_temps();
}
void SAPFunctions::compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    // => Build basis function values <= //
    block_index_ = block->index();
    if (!force_compute && cache_map_ && (cache_map_->find(block->index()) != cache_map_->end())) {
        current_basis_map_ = &(*cache_map_)[block->index()];
    } else {
        current_basis_map_ = &basis_values_;
        BasisFunctions::compute_functions(block);
    }
}
void SAPFunctions::print(std::string out, int print) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => SAPFunctions <=\n\n");
    printer->Printf("    Point Values:\n");
    for (std::map<std::string, std::shared_ptr<Vector> >::const_iterator it = point_values_.begin();
         it != point_values_.end(); it++) {
        printer->Printf("    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");
    BasisFunctions::print(out, print);
}
std::vector<SharedMatrix> SAPFunctions::D_scratch() {
    throw PSIEXCEPTION("SAPFunctions::density matrices are not appropriate. Read the source.");
}
void SAPFunctions::compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    throw PSIEXCEPTION("SAPFunctions::orbitals are not appropriate. Read the source.");
}
void SAPFunctions::set_Cs(SharedMatrix /*Ca_AO*/, SharedMatrix /*Cb_AO*/) {
    throw PSIEXCEPTION("SAPFunctions::orbitals are not appropriate. Read the source.");
}
void SAPFunctions::set_Cs(SharedMatrix /*Ca_AO*/) {
    throw PSIEXCEPTION("SAPFunctions::orbitals are not appropriate. Read the source.");
}
void SAPFunctions::set_pointers(SharedMatrix /*Ca_AO*/, SharedMatrix /*Cb_AO*/) {
    throw PSIEXCEPTION("SAPFunctions::orbitals are not appropriate. Read the source.");
}
void SAPFunctions::set_pointers(SharedMatrix /*Ca_AO*/) {
    throw PSIEXCEPTION("SAPFunctions::orbitals are not appropriate. Read the source.");
}

RKSFunctions::RKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions)
    : PointFunctions(primary, max_points, max_functions) {
    set_ansatz(0);
    current_basis_map_ = &basis_values_;
}
RKSFunctions::~RKSFunctions() {}
std::vector<SharedMatrix> RKSFunctions::scratch() {
    std::vector<SharedMatrix> vec;
    vec.push_back(temp_);
    return vec;
}
std::vector<SharedMatrix> RKSFunctions::D_scratch() {
    std::vector<SharedMatrix> vec;
    vec.push_back(D_local_);
    return vec;
}
void RKSFunctions::build_temps() {
    temp_ = std::make_shared<Matrix>("Temp", max_points_, max_functions_);
    D_local_ = std::make_shared<Matrix>("Dlocal", max_functions_, max_functions_);
}
void RKSFunctions::allocate() {
    BasisFunctions::allocate();

    point_values_.clear();

    if (ansatz_ >= 0) {
        point_values_["RHO_A"] = std::make_shared<Vector>("RHO_A", max_points_);
    }

    if (ansatz_ >= 1) {
        point_values_["RHO_AX"] = std::make_shared<Vector>("RHO_AX", max_points_);
        point_values_["RHO_AY"] = std::make_shared<Vector>("RHO_AY", max_points_);
        point_values_["RHO_AZ"] = std::make_shared<Vector>("RHO_AZ", max_points_);
        point_values_["GAMMA_AA"] = std::make_shared<Vector>("GAMMA_AA", max_points_);
    }

    if (ansatz_ >= 2) {
        point_values_["RHO_XX"] = std::make_shared<Vector>("RHO_XX", max_points_);
        point_values_["RHO_YY"] = std::make_shared<Vector>("RHO_YY", max_points_);
        point_values_["RHO_ZZ"] = std::make_shared<Vector>("RHO_ZZ", max_points_);
        point_values_["TAU_A"] = std::make_shared<Vector>("TAU_A", max_points_);
    }
    build_temps();
}
void RKSFunctions::set_pointers(SharedMatrix D_AO) { D_AO_ = D_AO; }
void RKSFunctions::set_pointers(SharedMatrix /*Da_AO*/, SharedMatrix /*Db_AO*/) {
    throw PSIEXCEPTION("RKSFunctions::unrestricted pointers are not appropriate. Read the source.");
}
void RKSFunctions::compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    if (!D_AO_) throw PSIEXCEPTION("RKSFunctions: call set_pointers.");

    // => Build basis function values <= //
    block_index_ = block->index();
    if (!force_compute && cache_map_ && (cache_map_->find(block->index()) != cache_map_->end())) {
        current_basis_map_ = &(*cache_map_)[block->index()];
    } else {
        current_basis_map_ = &basis_values_;
        BasisFunctions::compute_functions(block);
    }

    // => Global information <= //
    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal = function_map.size();

    double** Tp = temp_->pointer();

    // => Build local D matrix <= //
    double** Dp = D_AO_->pointer();
    double** D2p = D_local_->pointer();

    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        for (int nl = 0; nl <= ml; nl++) {
            int ng = function_map[nl];

            double Dval = Dp[mg][ng];

            D2p[ml][nl] = Dval;
            D2p[nl][ml] = Dval;
        }
    }

    // => Build LSDA quantities <= //
    double** phip = basis_value("PHI")->pointer();
    double* rhoap = point_value("RHO_A")->pointer();
    size_t coll_funcs = basis_value("PHI")->ncol();

    // Rho_a = 2.0 * D_xy phi_xa phi_ya
    C_DGEMM('N', 'N', npoints, nlocal, nlocal, 2.0, phip[0], coll_funcs, D2p[0], nglobal, 0.0, Tp[0], nglobal);
    for (int P = 0; P < npoints; P++) {
        rhoap[P] = C_DDOT(nlocal, phip[P], 1, Tp[P], 1);
    }

    // => Build GGA quantities <= //
    // Rho^l_a = D_xy phi_xa phi^l_ya
    if (ansatz_ >= 1) {
        double** phixp = basis_value("PHI_X")->pointer();
        double** phiyp = basis_value("PHI_Y")->pointer();
        double** phizp = basis_value("PHI_Z")->pointer();
        double* rhoaxp = point_value("RHO_AX")->pointer();
        double* rhoayp = point_value("RHO_AY")->pointer();
        double* rhoazp = point_value("RHO_AZ")->pointer();
        double* gammaaap = point_value("GAMMA_AA")->pointer();

        for (int P = 0; P < npoints; P++) {
            // 2.0 for Px D P + P D Px
            double rho_x = 2.0 * C_DDOT(nlocal, phixp[P], 1, Tp[P], 1);
            double rho_y = 2.0 * C_DDOT(nlocal, phiyp[P], 1, Tp[P], 1);
            double rho_z = 2.0 * C_DDOT(nlocal, phizp[P], 1, Tp[P], 1);
            rhoaxp[P] = rho_x;
            rhoayp[P] = rho_y;
            rhoazp[P] = rho_z;
            gammaaap[P] = rho_x * rho_x + rho_y * rho_y + rho_z * rho_z;
        }
    }

    // => Build Meta quantities <= //
    if (ansatz_ >= 2) {
        double** phixp = basis_value("PHI_X")->pointer();
        double** phiyp = basis_value("PHI_Y")->pointer();
        double** phizp = basis_value("PHI_Z")->pointer();
        double* taup = point_value("TAU_A")->pointer();

        std::fill(taup, taup + npoints, 0.0);

        double** phi[3];
        phi[0] = phixp;
        phi[1] = phiyp;
        phi[2] = phizp;

        for (int x = 0; x < 3; x++) {
            double** phic = phi[x];
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phic[0], coll_funcs, D2p[0], nglobal, 0.0, Tp[0], nglobal);
            for (int P = 0; P < npoints; P++) {
                taup[P] += C_DDOT(nlocal, phic[P], 1, Tp[P], 1);
            }
        }

        // Kinetic terms
        // double** phixxp = basis_value("PHI_XX")->pointer();
        // double** phixyp = basis_value("PHI_XY")->pointer();
        // double** phixzp = basis_value("PHI_XZ")->pointer();
        // double** phiyyp = basis_value("PHI_YY")->pointer();
        // double** phiyzp = basis_value("PHI_YZ")->pointer();
        // double** phizzp = basis_value("PHI_ZZ")->pointer();

        // double* rhoxxp = point_values_["RHO_XX"]->pointer();
        // double* rhoyyp = point_values_["RHO_YY"]->pointer();
        // double* rhozzp = point_values_["RHO_ZZ"]->pointer();

        // double* laplp = point_values_["LAPL_RHO_A"]->pointer();

        // // Diagonal terms phi^xx_a D_ab phi_b
        // for (int P = 0; P < npoints; P++) {
        //      rhoxxp[P] = 2.0 * C_DDOT(nlocal,phixxp[P],1,Tp[P],1);
        //      rhoyyp[P] = 2.0 * C_DDOT(nlocal,phiyyp[P],1,Tp[P],1);
        //      rhozzp[P] = 2.0 * C_DDOT(nlocal,phizzp[P],1,Tp[P],1);
        // }

        // // Cross terms phi^x_a D_ab phi^x_b
        // C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phixp[0],coll_funcs,D2p[0],nglobal,0.0,Tp[0],nglobal);
        // for (int P = 0; P < npoints; P++) {
        //      rhoxxp[P] += 2.0 * C_DDOT(nlocal,phixp[P],1,Tp[P],1);
        // }

        // C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phiyp[0],coll_funcs,D2p[0],nglobal,0.0,Tp[0],nglobal);
        // for (int P = 0; P < npoints; P++) {
        //      rhoyyp[P] += 2.0 * C_DDOT(nlocal,phiyp[P],1,Tp[P],1);
        // }

        // C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phizp[0],coll_funcs,D2p[0],nglobal,0.0,Tp[0],nglobal);
        // for (int P = 0; P < npoints; P++) {
        //      rhozzp[P] += 2.0 * C_DDOT(nlocal,phizp[P],1,Tp[P],1);
        // }

        // // Put it together
        // for (int P = 0; P < npoints; P++) {
        //     laplp[P]  = rhoxxp[P];
        //     laplp[P] += rhoyyp[P];
        //     laplp[P] += rhozzp[P];
        // }
    }
}

void RKSFunctions::set_Cs(SharedMatrix C_AO) {
    C_AO_ = C_AO;
    C_local_ = std::make_shared<Matrix>("C local", max_functions_, C_AO_->colspi()[0]);
    orbital_values_["PSI_A"] = std::make_shared<Matrix>("PSI_A", C_AO_->colspi()[0], max_points_);
    orbital_values_["PSI_B"] = orbital_value("PSI_A");
}
void RKSFunctions::set_Cs(SharedMatrix /*Ca_AO*/, SharedMatrix /*Cb_AO*/) {
    throw PSIEXCEPTION("RKSFunctions::unrestricted pointers are not appropriate. Read the source.");
}
void RKSFunctions::compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    // => Build basis function values <= //
    block_index_ = block->index();
    if (!force_compute && cache_map_ && (cache_map_->find(block->index()) != cache_map_->end())) {
        current_basis_map_ = &(*cache_map_)[block->index()];
    } else {
        current_basis_map_ = &basis_values_;
        BasisFunctions::compute_functions(block);
    }
    // timer_off("Functions: Points");

    // => Global information <= //

    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal = function_map.size();

    // => Build local C matrix <= //

    int na = C_AO_->colspi()[0];
    double** Cap = C_AO_->pointer();
    double** Ca2p = C_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        C_DCOPY(na, Cap[mg], 1, Ca2p[ml], 1);
    }

    // => Build orbitals <= //

    double** phip = basis_value("PHI")->pointer();
    double** psiap = orbital_value("PSI_A")->pointer();
    size_t coll_funcs = basis_value("PHI")->ncol();

    C_DGEMM('T', 'T', na, npoints, nlocal, 1.0, Ca2p[0], na, phip[0], coll_funcs, 0.0, psiap[0], max_points_);
}

void RKSFunctions::print(std::string out, int print) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    std::string ans;
    if (ansatz_ == 0) {
        ans = "LSDA";
    } else if (ansatz_ == 1) {
        ans = "GGA";
    } else if (ansatz_ == 2) {
        ans = "Meta-GGA";
    }

    printer->Printf("   => RKSFunctions: %s Ansatz <=\n\n", ans.c_str());

    printer->Printf("    Point Values:\n");
    for (std::map<std::string, std::shared_ptr<Vector> >::const_iterator it = point_values_.begin();
         it != point_values_.end(); it++) {
        printer->Printf("    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");

    BasisFunctions::print(out, print);
}

UKSFunctions::UKSFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions)
    : PointFunctions(primary, max_points, max_functions) {
    set_ansatz(0);
    current_basis_map_ = &basis_values_;
}
UKSFunctions::~UKSFunctions() {}
std::vector<SharedMatrix> UKSFunctions::scratch() {
    std::vector<SharedMatrix> vec;
    vec.push_back(tempa_);
    vec.push_back(tempb_);
    return vec;
}
std::vector<SharedMatrix> UKSFunctions::D_scratch() {
    std::vector<SharedMatrix> vec;
    vec.push_back(Da_local_);
    vec.push_back(Db_local_);
    return vec;
}
void UKSFunctions::build_temps() {
    tempa_ = std::make_shared<Matrix>("Temp", max_points_, max_functions_);
    Da_local_ = std::make_shared<Matrix>("Dlocal", max_functions_, max_functions_);
    tempb_ = std::make_shared<Matrix>("Temp", max_points_, max_functions_);
    Db_local_ = std::make_shared<Matrix>("Dlocal", max_functions_, max_functions_);
}
void UKSFunctions::allocate() {
    BasisFunctions::allocate();

    point_values_.clear();

    if (ansatz_ >= 0) {
        point_values_["RHO_A"] = std::make_shared<Vector>("RHO_A", max_points_);
        point_values_["RHO_B"] = std::make_shared<Vector>("RHO_B", max_points_);
    }

    if (ansatz_ >= 1) {
        point_values_["RHO_AX"] = std::make_shared<Vector>("RHO_AX", max_points_);
        point_values_["RHO_AY"] = std::make_shared<Vector>("RHO_AY", max_points_);
        point_values_["RHO_AZ"] = std::make_shared<Vector>("RHO_AZ", max_points_);
        point_values_["RHO_BX"] = std::make_shared<Vector>("RHO_BX", max_points_);
        point_values_["RHO_BY"] = std::make_shared<Vector>("RHO_BY", max_points_);
        point_values_["RHO_BZ"] = std::make_shared<Vector>("RHO_BZ", max_points_);
        point_values_["GAMMA_AA"] = std::make_shared<Vector>("GAMMA_AA", max_points_);
        point_values_["GAMMA_AB"] = std::make_shared<Vector>("GAMMA_AB", max_points_);
        point_values_["GAMMA_BB"] = std::make_shared<Vector>("GAMMA_BB", max_points_);
    }

    if (ansatz_ >= 2) {
        point_values_["TAU_A"] = std::make_shared<Vector>("TAU_A", max_points_);
        point_values_["TAU_B"] = std::make_shared<Vector>("TAU_A", max_points_);
    }
    build_temps();
}
void UKSFunctions::set_pointers(SharedMatrix /*Da_AO*/) {
    throw PSIEXCEPTION("UKSFunctions::restricted pointers are not appropriate. Read the source.");
}
void UKSFunctions::set_pointers(SharedMatrix Da_AO, SharedMatrix Db_AO) {
    Da_AO_ = Da_AO;
    Db_AO_ = Db_AO;
}
void UKSFunctions::compute_points(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    if (!Da_AO_) throw PSIEXCEPTION("UKSFunctions: call set_pointers.");

    // => Build basis function values <= //
    block_index_ = block->index();
    if (!force_compute && cache_map_ && (cache_map_->find(block->index()) != cache_map_->end())) {
        current_basis_map_ = &(*cache_map_)[block->index()];
    } else {
        current_basis_map_ = &basis_values_;
        BasisFunctions::compute_functions(block);
    }

    // => Global information <= //
    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal = function_map.size();

    double** Tap = tempa_->pointer();
    double** Tbp = tempb_->pointer();

    // => Build local D matrix <= //
    double** Dap = Da_AO_->pointer();
    double** Da2p = Da_local_->pointer();
    double** Dbp = Db_AO_->pointer();
    double** Db2p = Db_local_->pointer();

    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        for (int nl = 0; nl <= ml; nl++) {
            int ng = function_map[nl];

            double Daval = Dap[mg][ng];
            double Dbval = Dbp[mg][ng];

            Da2p[ml][nl] = Daval;
            Da2p[nl][ml] = Daval;
            Db2p[ml][nl] = Dbval;
            Db2p[nl][ml] = Dbval;
        }
    }

    // => Build LSDA quantities <= //
    double** phip = basis_value("PHI")->pointer();
    double* rhoap = point_value("RHO_A")->pointer();
    double* rhobp = point_value("RHO_B")->pointer();
    size_t coll_funcs = basis_value("PHI")->ncol();

    C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phip[0], coll_funcs, Da2p[0], nglobal, 0.0, Tap[0], nglobal);
    for (int P = 0; P < npoints; P++) {
        rhoap[P] = C_DDOT(nlocal, phip[P], 1, Tap[P], 1);
    }

    C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phip[0], coll_funcs, Db2p[0], nglobal, 0.0, Tbp[0], nglobal);
    for (int P = 0; P < npoints; P++) {
        rhobp[P] = C_DDOT(nlocal, phip[P], 1, Tbp[P], 1);
    }

    // => Build GGA quantities <= //
    if (ansatz_ >= 1) {
        double** phixp = basis_value("PHI_X")->pointer();
        double** phiyp = basis_value("PHI_Y")->pointer();
        double** phizp = basis_value("PHI_Z")->pointer();
        double* rhoaxp = point_value("RHO_AX")->pointer();
        double* rhoayp = point_value("RHO_AY")->pointer();
        double* rhoazp = point_value("RHO_AZ")->pointer();
        double* rhobxp = point_value("RHO_BX")->pointer();
        double* rhobyp = point_value("RHO_BY")->pointer();
        double* rhobzp = point_value("RHO_BZ")->pointer();
        double* gammaaap = point_value("GAMMA_AA")->pointer();
        double* gammaabp = point_value("GAMMA_AB")->pointer();
        double* gammabbp = point_value("GAMMA_BB")->pointer();

        for (int P = 0; P < npoints; P++) {
            // 2.0 for Px D P + P D Px
            double rhoa_x = 2.0 * C_DDOT(nlocal, phixp[P], 1, Tap[P], 1);
            double rhoa_y = 2.0 * C_DDOT(nlocal, phiyp[P], 1, Tap[P], 1);
            double rhoa_z = 2.0 * C_DDOT(nlocal, phizp[P], 1, Tap[P], 1);
            double rhob_x = 2.0 * C_DDOT(nlocal, phixp[P], 1, Tbp[P], 1);
            double rhob_y = 2.0 * C_DDOT(nlocal, phiyp[P], 1, Tbp[P], 1);
            double rhob_z = 2.0 * C_DDOT(nlocal, phizp[P], 1, Tbp[P], 1);
            rhoaxp[P] = rhoa_x;
            rhoayp[P] = rhoa_y;
            rhoazp[P] = rhoa_z;
            rhobxp[P] = rhob_x;
            rhobyp[P] = rhob_y;
            rhobzp[P] = rhob_z;
            gammaaap[P] = rhoa_x * rhoa_x + rhoa_y * rhoa_y + rhoa_z * rhoa_z;
            gammaabp[P] = rhoa_x * rhob_x + rhoa_y * rhob_y + rhoa_z * rhob_z;
            gammabbp[P] = rhob_x * rhob_x + rhob_y * rhob_y + rhob_z * rhob_z;
        }
    }

    // => Build Meta quantities <= //
    if (ansatz_ >= 2) {
        double** phixp = basis_value("PHI_X")->pointer();
        double** phiyp = basis_value("PHI_Y")->pointer();
        double** phizp = basis_value("PHI_Z")->pointer();
        double* tauap = point_value("TAU_A")->pointer();
        double* taubp = point_value("TAU_B")->pointer();

        std::fill(tauap, tauap + npoints, 0.0);
        std::fill(taubp, taubp + npoints, 0.0);

        double** phi[3];
        phi[0] = phixp;
        phi[1] = phiyp;
        phi[2] = phizp;

        double* tau[2];
        tau[0] = tauap;
        tau[1] = taubp;

        double** D[2];
        D[0] = Da2p;
        D[1] = Db2p;

        double** T[2];
        T[0] = Tap;
        T[1] = Tbp;

        for (int x = 0; x < 3; x++) {
            for (int t = 0; t < 2; t++) {
                double** phic = phi[x];
                double** Dc = D[t];
                double** Tc = T[t];
                double* tauc = tau[t];
                C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phic[0], coll_funcs, Dc[0], nglobal, 0.0, Tc[0],
                        nglobal);
                for (int P = 0; P < npoints; P++) {
                    tauc[P] += 0.5 * C_DDOT(nlocal, phic[P], 1, Tc[P], 1);
                }
            }
        }
    }
}

void UKSFunctions::set_Cs(SharedMatrix /*Ca_AO*/) {
    throw PSIEXCEPTION("UKSFunctions::restricted pointers are not appropriate. Read the source.");
}
void UKSFunctions::set_Cs(SharedMatrix Ca_AO, SharedMatrix Cb_AO) {
    Ca_AO_ = Ca_AO;
    Cb_AO_ = Cb_AO;
    Ca_local_ = std::make_shared<Matrix>("Ca local", max_functions_, Ca_AO_->colspi()[0]);
    Cb_local_ = std::make_shared<Matrix>("Cb local", max_functions_, Cb_AO_->colspi()[0]);
    orbital_values_["PSI_A"] = std::make_shared<Matrix>("PSI_A", Ca_AO_->colspi()[0], max_points_);
    orbital_values_["PSI_B"] = std::make_shared<Matrix>("PSI_B", Cb_AO_->colspi()[0], max_points_);
}
void UKSFunctions::compute_orbitals(std::shared_ptr<BlockOPoints> block, bool force_compute) {
    // => Build basis function values <= //
    block_index_ = block->index();
    if (!force_compute && cache_map_ && (cache_map_->find(block->index()) != cache_map_->end())) {
        current_basis_map_ = &(*cache_map_)[block->index()];
    } else {
        current_basis_map_ = &basis_values_;
        BasisFunctions::compute_functions(block);
    }

    // => Global information <= //

    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal = function_map.size();

    // => Build local C matrix <= //

    int na = Ca_AO_->colspi()[0];
    double** Cap = Ca_AO_->pointer();
    double** Ca2p = Ca_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        C_DCOPY(na, Cap[mg], 1, Ca2p[ml], 1);
    }

    int nb = Cb_AO_->colspi()[0];
    double** Cbp = Cb_AO_->pointer();
    double** Cb2p = Cb_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        C_DCOPY(na, Cbp[mg], 1, Cb2p[ml], 1);
    }

    // => Build orbitals <= //

    double** phip = basis_value("PHI")->pointer();
    double** psiap = orbital_value("PSI_A")->pointer();
    double** psibp = orbital_value("PSI_B")->pointer();
    size_t coll_funcs = basis_value("PHI")->ncol();

    C_DGEMM('T', 'T', na, npoints, nlocal, 1.0, Ca2p[0], na, phip[0], coll_funcs, 0.0, psiap[0], max_points_);
    C_DGEMM('T', 'T', nb, npoints, nlocal, 1.0, Cb2p[0], nb, phip[0], coll_funcs, 0.0, psibp[0], max_points_);
}

void UKSFunctions::print(std::string out, int print) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    std::string ans;
    if (ansatz_ == 0) {
        ans = "LSDA";
    } else if (ansatz_ == 1) {
        ans = "GGA";
    } else if (ansatz_ == 2) {
        ans = "Meta-GGA";
    }

    printer->Printf("   => UKSFunctions: %s Ansatz <=\n\n", ans.c_str());

    printer->Printf("    Point Values:\n");
    for (std::map<std::string, std::shared_ptr<Vector> >::const_iterator it = point_values_.begin();
         it != point_values_.end(); it++) {
        printer->Printf("    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");

    BasisFunctions::print(out, print);
}

PointFunctions::PointFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions)
    : BasisFunctions(primary, max_points, max_functions) {
    set_ansatz(0);
}
PointFunctions::~PointFunctions() {}
SharedVector PointFunctions::point_value(const std::string& key) { return point_values_[key]; }

SharedMatrix PointFunctions::orbital_value(const std::string& key) { return orbital_values_[key]; }

BasisFunctions::BasisFunctions(std::shared_ptr<BasisSet> primary, int max_points, int max_functions)
    : primary_(primary), max_points_(max_points), max_functions_(max_functions) {
    if (!primary_->has_puream()) {
        puream_ = false;
        return;
    }

    puream_ = true;
    set_deriv(0);
}
BasisFunctions::~BasisFunctions() {}
void BasisFunctions::allocate() {
    basis_values_.clear();
    basis_temps_.clear();

    int max_am = primary_->max_am();
    int max_cart = (max_am + 1) * (max_am + 2) / 2;

    if (deriv_ >= 0) {
        basis_values_["PHI"] = std::make_shared<Matrix>("PHI", max_points_, max_functions_);
        basis_temps_["PHI"] = std::make_shared<Matrix>("PHI", max_points_, max_functions_);
    }

    if (deriv_ >= 1) {
        basis_values_["PHI_X"] = std::make_shared<Matrix>("PHI_X", max_points_, max_functions_);
        basis_values_["PHI_Y"] = std::make_shared<Matrix>("PHI_Y", max_points_, max_functions_);
        basis_values_["PHI_Z"] = std::make_shared<Matrix>("PHI_Z", max_points_, max_functions_);
        basis_temps_["PHI_X"] = std::make_shared<Matrix>("PHI_X", max_points_, max_functions_);
        basis_temps_["PHI_Y"] = std::make_shared<Matrix>("PHI_Y", max_points_, max_functions_);
        basis_temps_["PHI_Z"] = std::make_shared<Matrix>("PHI_Z", max_points_, max_functions_);
    }

    if (deriv_ >= 2) {
        basis_values_["PHI_XX"] = std::make_shared<Matrix>("PHI_XX", max_points_, max_functions_);
        basis_values_["PHI_XY"] = std::make_shared<Matrix>("PHI_XY", max_points_, max_functions_);
        basis_values_["PHI_XZ"] = std::make_shared<Matrix>("PHI_XZ", max_points_, max_functions_);
        basis_values_["PHI_YY"] = std::make_shared<Matrix>("PHI_YY", max_points_, max_functions_);
        basis_values_["PHI_YZ"] = std::make_shared<Matrix>("PHI_YZ", max_points_, max_functions_);
        basis_values_["PHI_ZZ"] = std::make_shared<Matrix>("PHI_ZZ", max_points_, max_functions_);
        basis_temps_["PHI_XX"] = std::make_shared<Matrix>("PHI_XX", max_points_, max_functions_);
        basis_temps_["PHI_XY"] = std::make_shared<Matrix>("PHI_XY", max_points_, max_functions_);
        basis_temps_["PHI_XZ"] = std::make_shared<Matrix>("PHI_XZ", max_points_, max_functions_);
        basis_temps_["PHI_YY"] = std::make_shared<Matrix>("PHI_YY", max_points_, max_functions_);
        basis_temps_["PHI_YZ"] = std::make_shared<Matrix>("PHI_YZ", max_points_, max_functions_);
        basis_temps_["PHI_ZZ"] = std::make_shared<Matrix>("PHI_ZZ", max_points_, max_functions_);
    }

    if (deriv_ >= 3) throw PSIEXCEPTION("BasisFunctions: Only up to Hessians are currently supported");
}
void BasisFunctions::compute_functions(std::shared_ptr<BlockOPoints> block) {
    // Pull out data
    int npoints = block->npoints();
    double* x = block->x();
    double* y = block->y();
    double* z = block->z();
    std::vector<double> xyz(npoints * 3);
    ::memcpy(xyz.data(), x, sizeof(double) * npoints);
    ::memcpy(xyz.data() + npoints, y, sizeof(double) * npoints);
    ::memcpy(xyz.data() + 2 * npoints, z, sizeof(double) * npoints);

    const std::vector<int>& shells = block->shells_local_to_global();

    // Declare tmps
    std::vector<double> center(3, 0.0);

    // Declare pointers
    double *tmpp, *tmp_xp, *tmp_yp, *tmp_zp;
    double *tmp_xxp, *tmp_xyp, *tmp_xzp, *tmp_yyp, *tmp_yzp, *tmp_zzp;
    double *valuesp, *values_xp, *values_yp, *values_zp;
    double *values_xxp, *values_xyp, *values_xzp, *values_yyp, *values_yzp, *values_zzp;

    if (deriv_ >= 0) {
        tmpp = basis_temps_["PHI"]->pointer()[0];
        valuesp = basis_values_["PHI"]->pointer()[0];
    }
    if (deriv_ >= 1) {
        tmp_xp = basis_temps_["PHI_X"]->pointer()[0];
        tmp_yp = basis_temps_["PHI_Y"]->pointer()[0];
        tmp_zp = basis_temps_["PHI_Z"]->pointer()[0];
        values_xp = basis_values_["PHI_X"]->pointer()[0];
        values_yp = basis_values_["PHI_Y"]->pointer()[0];
        values_zp = basis_values_["PHI_Z"]->pointer()[0];
    }
    if (deriv_ >= 2) {
        tmp_xxp = basis_temps_["PHI_XX"]->pointer()[0];
        tmp_xyp = basis_temps_["PHI_XY"]->pointer()[0];
        tmp_xzp = basis_temps_["PHI_XZ"]->pointer()[0];
        tmp_yyp = basis_temps_["PHI_YY"]->pointer()[0];
        tmp_yzp = basis_temps_["PHI_YZ"]->pointer()[0];
        tmp_zzp = basis_temps_["PHI_ZZ"]->pointer()[0];
        values_xxp = basis_values_["PHI_XX"]->pointer()[0];
        values_xyp = basis_values_["PHI_XY"]->pointer()[0];
        values_xzp = basis_values_["PHI_XZ"]->pointer()[0];
        values_yyp = basis_values_["PHI_YY"]->pointer()[0];
        values_yzp = basis_values_["PHI_YZ"]->pointer()[0];
        values_zzp = basis_values_["PHI_ZZ"]->pointer()[0];
    }

    int nvals = 0;
    for (size_t Qlocal = 0; Qlocal < shells.size(); Qlocal++) {
        int Qglobal = shells[Qlocal];
        const GaussianShell& Qshell = primary_->shell(Qglobal);
        Vector3 v = Qshell.center();
        int L = Qshell.am();
        int nQ = Qshell.nfunction();
        int nprim = Qshell.nprimitive();
        const double* alpha = Qshell.exps();
        const double* norm = Qshell.coefs();

        // Copy over centerp to a double*
        center[0] = v[0];
        center[1] = v[1];
        center[2] = v[2];

        // Make new pointers, gg computes along rows so we need to skip down `nval` rows.
        size_t row_shift = static_cast<size_t>(nvals) * npoints;
        double* phi_start = tmpp + row_shift;
        const int order = (int)puream_ ? GG_SPHERICAL_GAUSSIAN : GG_CARTESIAN_CCA;

        // Copmute collocation
        if (deriv_ == 0) {
            gg_collocation(L, npoints, xyz.data(), 1, nprim, norm, alpha, center.data(), order, phi_start);
        } else if (deriv_ == 1) {
            double* phi_x_start = tmp_xp + row_shift;
            double* phi_y_start = tmp_yp + row_shift;
            double* phi_z_start = tmp_zp + row_shift;
            gg_collocation_deriv1(L, npoints, xyz.data(), 1, nprim, norm, alpha, center.data(), order, phi_start,
                                  phi_x_start, phi_y_start, phi_z_start);

        } else if (deriv_ == 2) {
            double* phi_x_start = tmp_xp + row_shift;
            double* phi_y_start = tmp_yp + row_shift;
            double* phi_z_start = tmp_zp + row_shift;
            double* phi_xx_start = tmp_xxp + row_shift;
            double* phi_xy_start = tmp_xyp + row_shift;
            double* phi_xz_start = tmp_xzp + row_shift;
            double* phi_yy_start = tmp_yyp + row_shift;
            double* phi_yz_start = tmp_yzp + row_shift;
            double* phi_zz_start = tmp_zzp + row_shift;
            gg_collocation_deriv2(L, npoints, xyz.data(), 1, nprim, norm, alpha, center.data(), order, phi_start,
                                  phi_x_start, phi_y_start, phi_z_start, phi_xx_start, phi_xy_start, phi_xz_start,
                                  phi_yy_start, phi_yz_start, phi_zz_start);
        }

        if (puream_) {
            nvals += 2 * L + 1;
        } else {
            nvals += nQ;  // Cartesian is already computed
        }
    }

    // GG spits it out tranpose of what we need
    int nso = max_functions_;
    gg_fast_transpose(nso, npoints, tmpp, valuesp);
    if (deriv_ >= 1) {
        gg_fast_transpose(nso, npoints, tmp_xp, values_xp);
        gg_fast_transpose(nso, npoints, tmp_yp, values_yp);
        gg_fast_transpose(nso, npoints, tmp_zp, values_zp);
    }
    if (deriv_ >= 2) {
        gg_fast_transpose(nso, npoints, tmp_xxp, values_xxp);
        gg_fast_transpose(nso, npoints, tmp_xyp, values_xyp);
        gg_fast_transpose(nso, npoints, tmp_xzp, values_xzp);
        gg_fast_transpose(nso, npoints, tmp_yyp, values_yyp);
        gg_fast_transpose(nso, npoints, tmp_yzp, values_yzp);
        gg_fast_transpose(nso, npoints, tmp_zzp, values_zzp);
    }
}
void BasisFunctions::print(std::string out, int print) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => BasisFunctions: Derivative = %d, Max Points = %d <=\n\n", deriv_, max_points_);

    printer->Printf("    Basis Values:\n");
    for (std::map<std::string, SharedMatrix>::const_iterator it = basis_values_.begin(); it != basis_values_.end();
         it++) {
        printer->Printf("    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");
}

}  // namespace psi
