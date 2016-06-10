/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <cmath>
#include "points.h"
#include "cubature.h"
#include "psiconfig.h"
#include "libparallel/ParallelPrinter.h"
namespace psi {

RKSFunctions::RKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions) :
    PointFunctions(primary,max_points,max_functions)
{
    set_ansatz(0);
}
RKSFunctions::~RKSFunctions()
{
}
std::vector<SharedMatrix> RKSFunctions::scratch()
{
    std::vector<SharedMatrix> vec;
    vec.push_back(temp_);
    return vec;
}
std::vector<SharedMatrix> RKSFunctions::D_scratch()
{
    std::vector<SharedMatrix> vec;
    vec.push_back(D_local_);
    return vec;
}
void RKSFunctions::build_temps()
{
    temp_ = SharedMatrix(new Matrix("Temp",max_points_,max_functions_));
    D_local_ = SharedMatrix(new Matrix("Dlocal",max_functions_,max_functions_));
}
void RKSFunctions::allocate()
{
    BasisFunctions::allocate();

    point_values_.clear();

    if (ansatz_ >= 0) {
        point_values_["RHO_A"] = boost::shared_ptr<Vector>(new Vector("RHO_A", max_points_));
        point_values_["RHO_B"] = point_values_["RHO_A"];
    }

    if (ansatz_ >= 1) {
        point_values_["RHO_AX"] = boost::shared_ptr<Vector>(new Vector("RHO_AX", max_points_));
        point_values_["RHO_AY"] = boost::shared_ptr<Vector>(new Vector("RHO_AY", max_points_));
        point_values_["RHO_AZ"] = boost::shared_ptr<Vector>(new Vector("RHO_AZ", max_points_));
        point_values_["RHO_BX"] = point_values_["RHO_AX"];
        point_values_["RHO_BY"] = point_values_["RHO_AY"];
        point_values_["RHO_BZ"] = point_values_["RHO_AZ"];
        point_values_["GAMMA_AA"] = boost::shared_ptr<Vector>(new Vector("GAMMA_AA", max_points_));
        point_values_["GAMMA_AB"] = point_values_["GAMMA_AA"];
        point_values_["GAMMA_BB"] = point_values_["GAMMA_AA"];
    }

    if (ansatz_ >= 2) {
        point_values_["TAU_A"] = boost::shared_ptr<Vector>(new Vector("TAU_A", max_points_));
        point_values_["TAU_B"] = point_values_["TAU_A"];
    }
}
void RKSFunctions::set_pointers(SharedMatrix D_AO)
{
    D_AO_ = D_AO;
    build_temps();
}
void RKSFunctions::set_pointers(SharedMatrix /*Da_AO*/, SharedMatrix /*Db_AO*/)
{
    throw PSIEXCEPTION("RKSFunctions::unrestricted pointers are not appropriate. Read the source.");
}
void RKSFunctions::compute_points(boost::shared_ptr<BlockOPoints> block)
{
    if (!D_AO_) 
        throw PSIEXCEPTION("RKSFunctions: call set_pointers.");

    // => Build basis function values <= //
    timer_on("Points");
    BasisFunctions::compute_functions(block);
    timer_off("Points");

    // => Global information <= //
    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal  = function_map.size();

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
    double** phip = basis_values_["PHI"]->pointer();
    double* rhoap = point_values_["RHO_A"]->pointer();

    C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phip[0],nglobal,D2p[0],nglobal,0.0,Tp[0],nglobal);
    for (int P = 0; P < npoints; P++) {
        rhoap[P] = C_DDOT(nlocal,phip[P],1,Tp[P],1);
    }

    // => Build GGA quantities <= //
    if (ansatz_ >= 1) {

        double** phixp = basis_values_["PHI_X"]->pointer();
        double** phiyp = basis_values_["PHI_Y"]->pointer();
        double** phizp = basis_values_["PHI_Z"]->pointer();
        double* rhoaxp = point_values_["RHO_AX"]->pointer();
        double* rhoayp = point_values_["RHO_AY"]->pointer();
        double* rhoazp = point_values_["RHO_AZ"]->pointer();
        double* gammaaap = point_values_["GAMMA_AA"]->pointer();

        for (int P = 0; P < npoints; P++) {
            double rho_x = 2.0 * C_DDOT(nlocal,phixp[P],1,Tp[P],1);
            double rho_y = 2.0 * C_DDOT(nlocal,phiyp[P],1,Tp[P],1);
            double rho_z = 2.0 * C_DDOT(nlocal,phizp[P],1,Tp[P],1);
            rhoaxp[P] = rho_x;
            rhoayp[P] = rho_y;
            rhoazp[P] = rho_z;
            gammaaap[P] = rho_x * rho_x + rho_y * rho_y + rho_z * rho_z;
        }
    }

    // => Build Meta quantities <= //
    if (ansatz_ >= 2) {
        double** phixp = basis_values_["PHI_X"]->pointer();
        double** phiyp = basis_values_["PHI_Y"]->pointer();
        double** phizp = basis_values_["PHI_Z"]->pointer();
        double* taup = point_values_["TAU_A"]->pointer();

        ::memset((void*) taup, '\0', sizeof(double) * npoints);

        double** phi[3];
        phi[0] = phixp;
        phi[1] = phiyp;
        phi[2] = phizp;

        for (int x = 0; x < 3; x++) {
            double** phic = phi[x];
            C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phic[0],nglobal,D2p[0],nglobal,0.0,Tp[0],nglobal);
            for (int P = 0; P < npoints; P++) {
                taup[P] += C_DDOT(nlocal, phic[P], 1, Tp[P], 1);
            }
        }
    }
}

void RKSFunctions::set_Cs(SharedMatrix C_AO)
{
    C_AO_ = C_AO;
    C_local_ = boost::shared_ptr<Matrix>(new Matrix("C local", max_functions_, C_AO_->colspi()[0]));
    orbital_values_["PSI_A"] = boost::shared_ptr<Matrix>(new Matrix("PSI_A", C_AO_->colspi()[0], max_points_));
    orbital_values_["PSI_B"] = orbital_values_["PSI_A"];
}
void RKSFunctions::set_Cs(SharedMatrix /*Ca_AO*/, SharedMatrix /*Cb_AO*/)
{
    throw PSIEXCEPTION("RKSFunctions::unrestricted pointers are not appropriate. Read the source.");
}
void RKSFunctions::compute_orbitals(boost::shared_ptr<BlockOPoints> block)
{
    // => Build basis function values <= //

    timer_on("Points");
    BasisFunctions::compute_functions(block);
    timer_off("Points");

    // => Global information <= //

    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal  = function_map.size();

    // => Build local C matrix <= //

    int na = C_AO_->colspi()[0];
    double** Cap  = C_AO_->pointer();
    double** Ca2p = C_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        ::memcpy(Ca2p[ml],Cap[mg],na*sizeof(double));
    }

    // => Build orbitals <= //

    double** phip = basis_values_["PHI"]->pointer();
    double** psiap = orbital_values_["PSI_A"]->pointer();

    C_DGEMM('T','T',na,npoints,nlocal,1.0,Ca2p[0],na,phip[0],nglobal,0.0,psiap[0],max_points_);
}


void RKSFunctions::print(std::string out, int print) const
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));
   std::string ans;
    if (ansatz_ == 0) {
        ans = "LSDA";
    } else if (ansatz_ == 1) {
        ans = "GGA";
    } else if (ansatz_ == 2) {
        ans = "Meta-GGA";
    }

    printer->Printf( "   => RKSFunctions: %s Ansatz <=\n\n", ans.c_str());

    printer->Printf( "    Point Values:\n");
    for (std::map<std::string, boost::shared_ptr<Vector> >::const_iterator it = point_values_.begin();
        it != point_values_.end(); it++) {
        printer->Printf( "    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");

    BasisFunctions::print(out,print);
}

UKSFunctions::UKSFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions) :
    PointFunctions(primary,max_points, max_functions)
{
    set_ansatz(0);
}
UKSFunctions::~UKSFunctions()
{
}
std::vector<SharedMatrix> UKSFunctions::scratch()
{
    std::vector<SharedMatrix> vec;
    vec.push_back(tempa_);
    vec.push_back(tempb_);
    return vec;
}
std::vector<SharedMatrix> UKSFunctions::D_scratch()
{
    std::vector<SharedMatrix> vec;
    vec.push_back(Da_local_);
    vec.push_back(Db_local_);
    return vec;
}
void UKSFunctions::build_temps()
{
    tempa_ = SharedMatrix(new Matrix("Temp",max_points_,max_functions_));
    Da_local_ = SharedMatrix(new Matrix("Dlocal",max_functions_,max_functions_));
    tempb_ = SharedMatrix(new Matrix("Temp",max_points_,max_functions_));
    Db_local_ = SharedMatrix(new Matrix("Dlocal",max_functions_,max_functions_));
}
void UKSFunctions::allocate()
{
    BasisFunctions::allocate();

    point_values_.clear();

    if (ansatz_ >= 0) {
        point_values_["RHO_A"] = boost::shared_ptr<Vector>(new Vector("RHO_A", max_points_));
        point_values_["RHO_B"] = boost::shared_ptr<Vector>(new Vector("RHO_B", max_points_));
    }

    if (ansatz_ >= 1) {
        point_values_["RHO_AX"] = boost::shared_ptr<Vector>(new Vector("RHO_AX", max_points_));
        point_values_["RHO_AY"] = boost::shared_ptr<Vector>(new Vector("RHO_AY", max_points_));
        point_values_["RHO_AZ"] = boost::shared_ptr<Vector>(new Vector("RHO_AZ", max_points_));
        point_values_["RHO_BX"] = boost::shared_ptr<Vector>(new Vector("RHO_BX", max_points_));
        point_values_["RHO_BY"] = boost::shared_ptr<Vector>(new Vector("RHO_BY", max_points_));
        point_values_["RHO_BZ"] = boost::shared_ptr<Vector>(new Vector("RHO_BZ", max_points_));
        point_values_["GAMMA_AA"] = boost::shared_ptr<Vector>(new Vector("GAMMA_AA", max_points_));
        point_values_["GAMMA_AB"] = boost::shared_ptr<Vector>(new Vector("GAMMA_AB", max_points_));
        point_values_["GAMMA_BB"] = boost::shared_ptr<Vector>(new Vector("GAMMA_BB", max_points_));
    }

    if (ansatz_ >= 2) {
        point_values_["TAU_A"] = boost::shared_ptr<Vector>(new Vector("TAU_A", max_points_));
        point_values_["TAU_B"] = boost::shared_ptr<Vector>(new Vector("TAU_A", max_points_));
    }
}
void UKSFunctions::set_pointers(SharedMatrix /*Da_AO*/)
{
    throw PSIEXCEPTION("UKSFunctions::restricted pointers are not appropriate. Read the source.");
}
void UKSFunctions::set_pointers(SharedMatrix Da_AO, SharedMatrix Db_AO)
{
    Da_AO_ = Da_AO;
    Db_AO_ = Db_AO;
    build_temps();
}
void UKSFunctions::compute_points(boost::shared_ptr<BlockOPoints> block)
{
    if (!Da_AO_) 
        throw PSIEXCEPTION("UKSFunctions: call set_pointers.");

    // => Build basis function values <= //
    timer_on("Points");
    BasisFunctions::compute_functions(block);
    timer_off("Points");

    // => Global information <= //
    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal  = function_map.size();

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
    double** phip = basis_values_["PHI"]->pointer();
    double* rhoap = point_values_["RHO_A"]->pointer();
    double* rhobp = point_values_["RHO_B"]->pointer();

    C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phip[0],nglobal,Da2p[0],nglobal,0.0,Tap[0],nglobal);
    for (int P = 0; P < npoints; P++) {
        rhoap[P] = C_DDOT(nlocal,phip[P],1,Tap[P],1);
    }

    C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phip[0],nglobal,Db2p[0],nglobal,0.0,Tbp[0],nglobal);
    for (int P = 0; P < npoints; P++) {
        rhobp[P] = C_DDOT(nlocal,phip[P],1,Tbp[P],1);
    }

    // => Build GGA quantities <= //
    if (ansatz_ >= 1) {

        double** phixp = basis_values_["PHI_X"]->pointer();
        double** phiyp = basis_values_["PHI_Y"]->pointer();
        double** phizp = basis_values_["PHI_Z"]->pointer();
        double* rhoaxp = point_values_["RHO_AX"]->pointer();
        double* rhoayp = point_values_["RHO_AY"]->pointer();
        double* rhoazp = point_values_["RHO_AZ"]->pointer();
        double* rhobxp = point_values_["RHO_BX"]->pointer();
        double* rhobyp = point_values_["RHO_BY"]->pointer();
        double* rhobzp = point_values_["RHO_BZ"]->pointer();
        double* gammaaap = point_values_["GAMMA_AA"]->pointer();
        double* gammaabp = point_values_["GAMMA_AB"]->pointer();
        double* gammabbp = point_values_["GAMMA_BB"]->pointer();

        for (int P = 0; P < npoints; P++) {
            double rhoa_x = 2.0 * C_DDOT(nlocal,phixp[P],1,Tap[P],1);
            double rhoa_y = 2.0 * C_DDOT(nlocal,phiyp[P],1,Tap[P],1);
            double rhoa_z = 2.0 * C_DDOT(nlocal,phizp[P],1,Tap[P],1);
            double rhob_x = 2.0 * C_DDOT(nlocal,phixp[P],1,Tbp[P],1);
            double rhob_y = 2.0 * C_DDOT(nlocal,phiyp[P],1,Tbp[P],1);
            double rhob_z = 2.0 * C_DDOT(nlocal,phizp[P],1,Tbp[P],1);
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
        double** phixp = basis_values_["PHI_X"]->pointer();
        double** phiyp = basis_values_["PHI_Y"]->pointer();
        double** phizp = basis_values_["PHI_Z"]->pointer();
        double* tauap = point_values_["TAU_A"]->pointer();
        double* taubp = point_values_["TAU_B"]->pointer();

        ::memset((void*) tauap, '\0', sizeof(double) * npoints);
        ::memset((void*) taubp, '\0', sizeof(double) * npoints);

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
                double*  tauc = tau[t];
                C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phic[0],nglobal,Dc[0],nglobal,0.0,Tc[0],nglobal);
                for (int P = 0; P < npoints; P++) {
                    tauc[P] += C_DDOT(nlocal, phic[P], 1, Tc[P], 1);
                }
            }
        }
    }
}

void UKSFunctions::set_Cs(SharedMatrix /*Ca_AO*/)
{
    throw PSIEXCEPTION("UKSFunctions::restricted pointers are not appropriate. Read the source.");
}
void UKSFunctions::set_Cs(SharedMatrix Ca_AO, SharedMatrix Cb_AO)
{
    Ca_AO_ = Ca_AO;
    Cb_AO_ = Cb_AO;
    Ca_local_ = boost::shared_ptr<Matrix>(new Matrix("Ca local", max_functions_, Ca_AO_->colspi()[0]));
    Cb_local_ = boost::shared_ptr<Matrix>(new Matrix("Cb local", max_functions_, Cb_AO_->colspi()[0]));
    orbital_values_["PSI_A"] = boost::shared_ptr<Matrix>(new Matrix("PSI_A", Ca_AO_->colspi()[0], max_points_));
    orbital_values_["PSI_B"] = boost::shared_ptr<Matrix>(new Matrix("PSI_B", Cb_AO_->colspi()[0], max_points_));
}
void UKSFunctions::compute_orbitals(boost::shared_ptr<BlockOPoints> block)
{
    // => Build basis function values <= //

    timer_on("Points");
    BasisFunctions::compute_functions(block);
    timer_off("Points");

    // => Global information <= //

    int npoints = block->npoints();
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nglobal = max_functions_;
    int nlocal  = function_map.size();

    // => Build local C matrix <= //

    int na = Ca_AO_->colspi()[0];
    double** Cap  = Ca_AO_->pointer();
    double** Ca2p = Ca_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        ::memcpy(Ca2p[ml],Cap[mg],na*sizeof(double));
    }

    int nb = Cb_AO_->colspi()[0];
    double** Cbp  = Cb_AO_->pointer();
    double** Cb2p = Cb_local_->pointer();
    for (int ml = 0; ml < nlocal; ml++) {
        int mg = function_map[ml];
        ::memcpy(Cb2p[ml],Cbp[mg],nb*sizeof(double));
    }

    // => Build orbitals <= //

    double** phip = basis_values_["PHI"]->pointer();
    double** psiap = orbital_values_["PSI_A"]->pointer();
    double** psibp = orbital_values_["PSI_B"]->pointer();

    C_DGEMM('T','T',na,npoints,nlocal,1.0,Ca2p[0],na,phip[0],nglobal,0.0,psiap[0],max_points_);
    C_DGEMM('T','T',nb,npoints,nlocal,1.0,Cb2p[0],nb,phip[0],nglobal,0.0,psibp[0],max_points_);
}


void UKSFunctions::print(std::string out, int print) const
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));
   std::string ans;
    if (ansatz_ == 0) {
        ans = "LSDA";
    } else if (ansatz_ == 1) {
        ans = "GGA";
    } else if (ansatz_ == 2) {
        ans = "Meta-GGA";
    }

    printer->Printf( "   => UKSFunctions: %s Ansatz <=\n\n", ans.c_str());

    printer->Printf( "    Point Values:\n");
    for (std::map<std::string, boost::shared_ptr<Vector> >::const_iterator it = point_values_.begin();
        it != point_values_.end(); it++) {
        printer->Printf( "    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");

    BasisFunctions::print(out,print);
}

PointFunctions::PointFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions) :
    BasisFunctions(primary,max_points, max_functions)
{
    set_ansatz(0);
}
PointFunctions::~PointFunctions()
{
}
SharedVector PointFunctions::point_value(const std::string& key)
{
    return point_values_[key];
}

SharedMatrix PointFunctions::orbital_value(const std::string& key)
{
    return orbital_values_[key];
}

BasisFunctions::BasisFunctions(boost::shared_ptr<BasisSet> primary, int max_points, int max_functions) :
    primary_(primary), max_points_(max_points), max_functions_(max_functions)
{
    build_spherical();
    set_deriv(0);
}
BasisFunctions::~BasisFunctions()
{
}
void BasisFunctions::build_spherical()
{
    if (!primary_->has_puream()) {
        puream_ = false;
        return;
    }

    puream_ = true;

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));

    for (int L = 0; L <= primary_->max_am(); L++) {
        std::vector<boost::tuple<int,int,double> > comp;
        boost::shared_ptr<SphericalTransformIter> trans(fact->spherical_transform_iter(L));
        for (trans->first(); !trans->is_done();trans->next()) {
            comp.push_back(boost::tuple<int,int,double>(
                trans->pureindex(),
                trans->cartindex(),
                trans->coef()));
        }
        spherical_transforms_.push_back(comp);
    }
}
void BasisFunctions::allocate()
{
    basis_values_.clear();
    basis_temps_.clear();

    int max_am = primary_->max_am();
    int max_cart = (max_am + 1) * (max_am + 2) / 2;

    if (deriv_ >= 0) {
        basis_values_["PHI"] = SharedMatrix (new Matrix("PHI", max_points_, max_functions_));
        basis_temps_["PHI"] = SharedMatrix (new Matrix("PHI", max_points_, max_cart));
    }

    if (deriv_ >= 1) {
        basis_values_["PHI_X"] = SharedMatrix (new Matrix("PHI_X", max_points_, max_functions_));
        basis_values_["PHI_Y"] = SharedMatrix (new Matrix("PHI_Y", max_points_, max_functions_));
        basis_values_["PHI_Z"] = SharedMatrix (new Matrix("PHI_Z", max_points_, max_functions_));
        basis_temps_["PHI_X"] = SharedMatrix (new Matrix("PHI_X", max_points_, max_cart));
        basis_temps_["PHI_Y"] = SharedMatrix (new Matrix("PHI_Y", max_points_, max_cart));
        basis_temps_["PHI_Z"] = SharedMatrix (new Matrix("PHI_Z", max_points_, max_cart));
    }

    if (deriv_ >= 2) {
        basis_values_["PHI_XX"] = SharedMatrix (new Matrix("PHI_XX", max_points_, max_functions_));
        basis_values_["PHI_XY"] = SharedMatrix (new Matrix("PHI_XY", max_points_, max_functions_));
        basis_values_["PHI_XZ"] = SharedMatrix (new Matrix("PHI_XZ", max_points_, max_functions_));
        basis_values_["PHI_YY"] = SharedMatrix (new Matrix("PHI_YY", max_points_, max_functions_));
        basis_values_["PHI_YZ"] = SharedMatrix (new Matrix("PHI_YZ", max_points_, max_functions_));
        basis_values_["PHI_ZZ"] = SharedMatrix (new Matrix("PHI_ZZ", max_points_, max_functions_));
        basis_temps_["PHI_XX"] = SharedMatrix (new Matrix("PHI_XX", max_points_, max_cart));
        basis_temps_["PHI_XY"] = SharedMatrix (new Matrix("PHI_XY", max_points_, max_cart));
        basis_temps_["PHI_XZ"] = SharedMatrix (new Matrix("PHI_XZ", max_points_, max_cart));
        basis_temps_["PHI_YY"] = SharedMatrix (new Matrix("PHI_YY", max_points_, max_cart));
        basis_temps_["PHI_YZ"] = SharedMatrix (new Matrix("PHI_YZ", max_points_, max_cart));
        basis_temps_["PHI_ZZ"] = SharedMatrix (new Matrix("PHI_ZZ", max_points_, max_cart));
    }

    if (deriv_ >= 3)
        throw PSIEXCEPTION("BasisFunctions: Only up to Hessians are currently supported");
}
SharedMatrix BasisFunctions::basis_value(const std::string& key)
{
    return basis_values_[key];
}
void BasisFunctions::compute_functions(boost::shared_ptr<BlockOPoints> block)
{
    int max_am = primary_->max_am();
    int max_cart = (max_am + 1) * (max_am + 2) / 2;

    int nso = max_functions_;

    int npoints = block->npoints();
    double *restrict x = block->x();
    double *restrict y = block->y();
    double *restrict z = block->z();

    int maxL = primary_->max_am();

    double *restrict xc_pow = new double[maxL + 3];
    double *restrict yc_pow = new double[maxL + 3];
    double *restrict zc_pow = new double[maxL + 3];

    const std::vector<int>& shells = block->shells_local_to_global();

    int nsig_functions = block->functions_local_to_global().size();

    if (deriv_ == 0) {
        double** cartp = basis_temps_["PHI"]->pointer();
        double** purep = basis_values_["PHI"]->pointer();

        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(purep[P]),'\0',nsig_functions*sizeof(double));
        }

        int function_offset = 0;
        for (size_t Qlocal = 0; Qlocal < shells.size(); Qlocal++) {
            int Qglobal = shells[Qlocal];
            const GaussianShell& Qshell = primary_->shell(Qglobal);
            Vector3 v     = Qshell.center();
            int L         = Qshell.am();
            int nQ        = Qshell.nfunction();
            int nprim     = Qshell.nprimitive();
            const double *alpha = Qshell.exps();
            const double *norm  = Qshell.coefs();

            const std::vector<boost::tuple<int,int,double> >& transform = spherical_transforms_[L];

            xc_pow[0] = 1.0;
            yc_pow[0] = 1.0;
            zc_pow[0] = 1.0;

            // Computation of points
            for (int P = 0; P < npoints; P++) {

                double xc = x[P] - v[0];
                double yc = y[P] - v[1];
                double zc = z[P] - v[2];

                for (int LL = 1; LL < L + 1; LL++) {
                    xc_pow[LL] = xc_pow[LL - 1] * xc;
                    yc_pow[LL] = yc_pow[LL - 1] * yc;
                    zc_pow[LL] = zc_pow[LL - 1] * zc;
                }

                double R2 = xc * xc + yc * yc + zc * zc;
                double S0 = 0.0;
                for (int K = 0; K < nprim; K++) {
                    S0 += norm[K] * exp(-alpha[K] * R2);
                }

                for (int i=0, index = 0; i<=L; ++i) {
                    int l = L-i;
                    for (int j=0; j<=i; ++j, ++index) {
                        int m = i-j;
                        int n = j;

                        cartp[P][index] = S0 * xc_pow[l] * yc_pow[m] * zc_pow[n];
                    }
                }
            }

            // Spherical transform
            if (puream_) {
                for (size_t index = 0; index < transform.size(); index++) {
                    int pureindex = boost::get<0>(transform[index]);
                    int cartindex = boost::get<1>(transform[index]);
                    double coef   = boost::get<2>(transform[index]);

                    C_DAXPY(npoints,coef,&cartp[0][cartindex],max_cart,&purep[0][pureindex + function_offset],nso);
                }
            } else {
                for (int q = 0; q < nQ; q++) {
                    C_DCOPY(npoints,&cartp[0][q],max_cart,&purep[0][q + function_offset],nso);
                }
            }

            function_offset += nQ;
        }
    } else if (deriv_ == 1) {
        double** cartp = basis_temps_["PHI"]->pointer();
        double** cartxp = basis_temps_["PHI_X"]->pointer();
        double** cartyp = basis_temps_["PHI_Y"]->pointer();
        double** cartzp = basis_temps_["PHI_Z"]->pointer();
        double** purep = basis_values_["PHI"]->pointer();
        double** purexp = basis_values_["PHI_X"]->pointer();
        double** pureyp = basis_values_["PHI_Y"]->pointer();
        double** purezp = basis_values_["PHI_Z"]->pointer();

        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(purep[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purexp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(pureyp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purezp[P]),'\0',nsig_functions*sizeof(double));
        }

        int function_offset = 0;
        for (size_t Qlocal = 0; Qlocal < shells.size(); Qlocal++) {
            int Qglobal = shells[Qlocal];
            const GaussianShell& Qshell = primary_->shell(Qglobal);
            Vector3 v     = Qshell.center();
            int L         = Qshell.am();
            int nQ        = Qshell.nfunction();
            int nprim     = Qshell.nprimitive();
            const double *alpha = Qshell.exps();
            const double *norm  = Qshell.coefs();

            const std::vector<boost::tuple<int,int,double> >& transform = spherical_transforms_[L];

            xc_pow[0] = 0.0;
            yc_pow[0] = 0.0;
            zc_pow[0] = 0.0;
            xc_pow[1] = 1.0;
            yc_pow[1] = 1.0;
            zc_pow[1] = 1.0;

            for (int P = 0; P < npoints; P++) {

                double xc = x[P] - v[0];
                double yc = y[P] - v[1];
                double zc = z[P] - v[2];

                for (int LL = 2; LL < L + 2; LL++) {
                    xc_pow[LL] = xc_pow[LL - 1] * xc;
                    yc_pow[LL] = yc_pow[LL - 1] * yc;
                    zc_pow[LL] = zc_pow[LL - 1] * zc;
                }

                double R2 = xc * xc + yc * yc + zc * zc;
                double V1 = 0.0;
                double V2 = 0.0;
//                double V3 = 0.0;
                double T1,T2;
                for (int K = 0; K < nprim; K++) {
                    T1 =  norm[K] * exp(-alpha[K] * R2);
                    T2 =  -2.0 * alpha[K] * T1;
                    V1 += T1;
                    V2 += T2;
                }
                double S0 = V1;
                double SX = V2 * xc;
                double SY = V2 * yc;
                double SZ = V2 * zc;

                for (int i=0, index = 0; i<=L; ++i) {
                    int l = L-i+1;
                    for (int j=0; j<=i; ++j, ++index) {
                        int m = i-j+1;
                        int n = j+1;

                        int lp = l - 1;
                        int mp = m - 1;
                        int np = n - 1;

                        double xyz = xc_pow[l] * yc_pow[m] * zc_pow[n];
                        cartp[P][index] = S0 * xyz;
                        cartxp[P][index] = S0 * lp * xc_pow[l-1] * yc_pow[m] * zc_pow[n] + SX * xyz;
                        cartyp[P][index] = S0 * mp * xc_pow[l] * yc_pow[m-1] * zc_pow[n] + SY * xyz;
                        cartzp[P][index] = S0 * np * xc_pow[l] * yc_pow[m] * zc_pow[n-1] + SZ * xyz;
                    }
                }
            }

            // Spherical transform
            if (puream_) {
                for (size_t index = 0; index < transform.size(); index++) {
                    int pureindex = boost::get<0>(transform[index]);
                    int cartindex = boost::get<1>(transform[index]);
                    double coef   = boost::get<2>(transform[index]);

                    C_DAXPY(npoints,coef,&cartp[0][cartindex],max_cart,&purep[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartxp[0][cartindex],max_cart,&purexp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartyp[0][cartindex],max_cart,&pureyp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartzp[0][cartindex],max_cart,&purezp[0][pureindex + function_offset],nso);
                }
            } else {
                for (int q = 0; q < nQ; q++) {
                    C_DCOPY(npoints,&cartp[0][q],max_cart,&purep[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartxp[0][q],max_cart,&purexp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartyp[0][q],max_cart,&pureyp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartzp[0][q],max_cart,&purezp[0][q + function_offset],nso);
                }
            }

            function_offset += nQ;
        }
    } else if (deriv_ == 2) {
        double** cartp = basis_temps_["PHI"]->pointer();
        double** cartxp = basis_temps_["PHI_X"]->pointer();
        double** cartyp = basis_temps_["PHI_Y"]->pointer();
        double** cartzp = basis_temps_["PHI_Z"]->pointer();
        double** cartxxp = basis_temps_["PHI_XX"]->pointer();
        double** cartxyp = basis_temps_["PHI_XY"]->pointer();
        double** cartxzp = basis_temps_["PHI_XZ"]->pointer();
        double** cartyyp = basis_temps_["PHI_YY"]->pointer();
        double** cartyzp = basis_temps_["PHI_YZ"]->pointer();
        double** cartzzp = basis_temps_["PHI_ZZ"]->pointer();
        double** purep = basis_values_["PHI"]->pointer();
        double** purexp = basis_values_["PHI_X"]->pointer();
        double** pureyp = basis_values_["PHI_Y"]->pointer();
        double** purezp = basis_values_["PHI_Z"]->pointer();
        double** purexxp = basis_values_["PHI_XX"]->pointer();
        double** purexyp = basis_values_["PHI_XY"]->pointer();
        double** purexzp = basis_values_["PHI_XZ"]->pointer();
        double** pureyyp = basis_values_["PHI_YY"]->pointer();
        double** pureyzp = basis_values_["PHI_YZ"]->pointer();
        double** purezzp = basis_values_["PHI_ZZ"]->pointer();

        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(purep[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purexp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(pureyp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purezp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purexxp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purexyp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purexzp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(pureyyp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(pureyzp[P]),'\0',nsig_functions*sizeof(double));
            ::memset(static_cast<void*>(purezzp[P]),'\0',nsig_functions*sizeof(double));
        }

        int function_offset = 0;
        for (size_t Qlocal = 0; Qlocal < shells.size(); Qlocal++) {
            int Qglobal = shells[Qlocal];
            const GaussianShell& Qshell = primary_->shell(Qglobal);
            Vector3 v     = Qshell.center();
            int L         = Qshell.am();
            int nQ        = Qshell.nfunction();
            int nprim     = Qshell.nprimitive();
            const double *alpha = Qshell.exps();
            const double *norm  = Qshell.coefs();

            const std::vector<boost::tuple<int,int,double> >& transform = spherical_transforms_[L];

            xc_pow[0] = 0.0;
            yc_pow[0] = 0.0;
            zc_pow[0] = 0.0;
            xc_pow[1] = 0.0;
            yc_pow[1] = 0.0;
            zc_pow[1] = 0.0;
            xc_pow[2] = 1.0;
            yc_pow[2] = 1.0;
            zc_pow[2] = 1.0;

            for (int P = 0; P < npoints; P++) {

                double xc = x[P] - v[0];
                double yc = y[P] - v[1];
                double zc = z[P] - v[2];

                for (int LL = 3; LL < L + 3; LL++) {
                    xc_pow[LL] = xc_pow[LL - 1] * xc;
                    yc_pow[LL] = yc_pow[LL - 1] * yc;
                    zc_pow[LL] = zc_pow[LL - 1] * zc;
                }

                double R2 = xc * xc + yc * yc + zc * zc;
                double V1 = 0.0;
                double V2 = 0.0;
                double V3 = 0.0;
                double T1,T2,T3;
                for (int K = 0; K < nprim; K++) {
                    T1 =  norm[K] * exp(-alpha[K] * R2);
                    T2 =  -2.0 * alpha[K] * T1;
                    T3 =  -2.0 * alpha[K] * T2;
                    V1 += T1;
                    V2 += T2;
                    V3 += T3;
                }
                double S = V1;
                double SX = V2 * xc;
                double SY = V2 * yc;
                double SZ = V2 * zc;
                double SXY = V3 * xc * yc;
                double SXZ = V3 * xc * zc;
                double SYZ = V3 * yc * zc;
                double SXX = V3 * xc * xc + V2;
                double SYY = V3 * yc * yc + V2;
                double SZZ = V3 * zc * zc + V2;

                for (int i=0, index = 0; i<=L; ++i) {
                    int l = L-i+2;
                    for (int j=0; j<=i; ++j, ++index) {
                        int m = i-j+2;
                        int n = j+2;

                        int lp = l - 2;
                        int mp = m - 2;
                        int np = n - 2;

                        double A = xc_pow[l] * yc_pow[m] * zc_pow[n];
                        double AX = lp * xc_pow[l-1] * yc_pow[m] * zc_pow[n];
                        double AY = mp * xc_pow[l] * yc_pow[m-1] * zc_pow[n];
                        double AZ = np * xc_pow[l] * yc_pow[m] * zc_pow[n-1];
                        double AXY = lp * mp * xc_pow[l-1] * yc_pow[m-1] * zc_pow[n];
                        double AXZ = lp * np * xc_pow[l-1] * yc_pow[m] * zc_pow[n-1];
                        double AYZ = mp * np * xc_pow[l] * yc_pow[m-1] * zc_pow[n-1];
                        double AXX = lp * (lp - 1) * xc_pow[l-2] * yc_pow[m] * zc_pow[n];
                        double AYY = mp * (mp - 1) * xc_pow[l] * yc_pow[m-2] * zc_pow[n];
                        double AZZ = np * (np - 1) * xc_pow[l] * yc_pow[m] * zc_pow[n-2];

                        cartp[P][index] = S * A;
                        cartxp[P][index] = S * AX + SX * A; 
                        cartyp[P][index] = S * AY + SY * A; 
                        cartzp[P][index] = S * AZ + SZ * A; 
                        cartxxp[P][index] = SXX * A + SX * AX + SX * AX + S * AXX;
                        cartyyp[P][index] = SYY * A + SY * AY + SY * AY + S * AYY;
                        cartzzp[P][index] = SZZ * A + SZ * AZ + SZ * AZ + S * AZZ;
                        cartxyp[P][index] = SXY * A + SX * AY + SY * AX + S * AXY;
                        cartxzp[P][index] = SXZ * A + SX * AZ + SZ * AX + S * AXZ;
                        cartyzp[P][index] = SYZ * A + SY * AZ + SZ * AY + S * AYZ;
                    }
                }
            }

            // Spherical transform
            if (puream_) {
                for (size_t index = 0; index < transform.size(); index++) {
                    int pureindex = boost::get<0>(transform[index]);
                    int cartindex = boost::get<1>(transform[index]);
                    double coef   = boost::get<2>(transform[index]);

                    C_DAXPY(npoints,coef,&cartp[0][cartindex],max_cart,&purep[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartxp[0][cartindex],max_cart,&purexp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartyp[0][cartindex],max_cart,&pureyp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartzp[0][cartindex],max_cart,&purezp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartxxp[0][cartindex],max_cart,&purexxp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartxyp[0][cartindex],max_cart,&purexyp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartxzp[0][cartindex],max_cart,&purexzp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartyyp[0][cartindex],max_cart,&pureyyp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartyzp[0][cartindex],max_cart,&pureyzp[0][pureindex + function_offset],nso);
                    C_DAXPY(npoints,coef,&cartzzp[0][cartindex],max_cart,&purezzp[0][pureindex + function_offset],nso);
                }
            } else {
                for (int q = 0; q < nQ; q++) {
                    C_DCOPY(npoints,&cartp[0][q],max_cart,&purep[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartxp[0][q],max_cart,&purexp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartyp[0][q],max_cart,&pureyp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartzp[0][q],max_cart,&purezp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartxxp[0][q],max_cart,&purexxp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartxyp[0][q],max_cart,&purexyp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartxzp[0][q],max_cart,&purexzp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartyyp[0][q],max_cart,&pureyyp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartyzp[0][q],max_cart,&pureyzp[0][q + function_offset],nso);
                    C_DCOPY(npoints,&cartzzp[0][q],max_cart,&purezzp[0][q + function_offset],nso);
                }
            }

            function_offset += nQ;
        }
    }

    delete[] xc_pow;
    delete[] yc_pow;
    delete[] zc_pow;
}
void BasisFunctions::print(std::string out, int print) const
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));
   printer->Printf( "   => BasisFunctions: Derivative = %d, Max Points = %d <=\n\n", deriv_, max_points_);

    printer->Printf( "    Basis Values:\n");
    for (std::map<std::string, SharedMatrix >::const_iterator it = basis_values_.begin();
        it != basis_values_.end(); it++) {
        printer->Printf( "    %s\n", (*it).first.c_str());
        if (print > 3) {
            (*it).second->print();
        }
    }
    printer->Printf("\n\n");
}

} // Namespace psi