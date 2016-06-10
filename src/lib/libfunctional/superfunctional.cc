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

#include <libmints/vector.h>
#include <libdisp/dispersion.h>
#include "superfunctional.h"
#include "functional.h"
#include "libparallel/ParallelPrinter.h"
using namespace psi;

namespace psi {

SuperFunctional::SuperFunctional()
{
    common_init();
}
SuperFunctional::~SuperFunctional()
{
}
void SuperFunctional::common_init()
{
    max_points_ = 0;
    deriv_ = 0;
    name_ = "";
    description_ = "";
    citation_ = "";
    x_omega_ = 0.0;
    c_omega_ = 0.0;
    x_alpha_ = 0.0;
    c_alpha_ = 0.0;
    c_ss_alpha_ = 0.0;
    c_os_alpha_ = 0.0;
}
boost::shared_ptr<SuperFunctional> SuperFunctional::blank()
{
    return boost::shared_ptr<SuperFunctional>(new SuperFunctional());
}
void SuperFunctional::print(std::string out, int level) const
{
    if (level < 1) return;
    boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
             boost::shared_ptr<OutFile>(new OutFile(out)));
    printer->Printf( "   => %s Composite Functional <= \n\n", name_.c_str());

    printer->Printf( "%s", description_.c_str());
    printer->Printf( "\n");
    
    printer->Printf( "%s", citation_.c_str());
    printer->Printf( "\n");
    
    printer->Printf( "    Points   = %14d\n", max_points_);
    printer->Printf( "    Deriv    = %14d\n", deriv_);
    printer->Printf( "    GGA      = %14s\n", (is_gga() ? "TRUE" : "FALSE"));
    printer->Printf( "    Meta     = %14s\n", (is_meta() ? "TRUE" : "FALSE"));
    printer->Printf( "\n");

    printer->Printf( "    X_LRC        = %14s\n", (is_x_lrc() ? "TRUE" : "FALSE"));
    printer->Printf( "    X_Hybrid     = %14s\n", (is_x_hybrid() ? "TRUE" : "FALSE"));
    printer->Printf( "    X_Alpha      = %14.6E\n", x_alpha_);
    printer->Printf( "    X_Omega      = %14.6E\n", x_omega_);
    printer->Printf( "    C_LRC        = %14s\n", (is_c_lrc() ? "TRUE" : "FALSE"));
    printer->Printf( "    C_Hybrid     = %14s\n", (is_c_hybrid() ? "TRUE" : "FALSE"));
    printer->Printf( "    C_Alpha      = %14.6E\n", c_alpha_);
    printer->Printf( "    C_Omega      = %14.6E\n", c_omega_);
    //printer->Printf( "    C_SCS_Hybrid = %14s\n", (is_c_scs_hybrid() ? "TRUE" : "FALSE"));
    //printer->Printf( "    C_SS_Alpha   = %14.6E\n", c_ss_alpha_);
    //printer->Printf( "    C_OS_Alpha   = %14.6E\n", c_os_alpha_);
    printer->Printf( "\n");
    
    printer->Printf( "   => Exchange Functionals <=\n\n");
    for (int i = 0; i < x_functionals_.size(); i++) {
        printer->Printf( "    %6.4f %7s", (1.0 - x_alpha_) * x_functionals_[i]->alpha(),
            x_functionals_[i]->name().c_str());
        if (x_functionals_[i]->omega()) {
            printer->Printf( " [omega = %6.4f]", x_functionals_[i]->omega());
        }
        printer->Printf("\n");
    }    
    if (x_omega_) {
        printer->Printf( "    %6.4f %7s [omega = %6.4f]\n", (1.0 - x_alpha_), "HF,LR", x_omega_);
    }
    if (x_alpha_) {
        printer->Printf( "    %6.4f %7s \n", x_alpha_, "HF");
    }
    printer->Printf( "\n");
     
    printer->Printf( "   => Correlation Functionals <=\n\n");
    for (int i = 0; i < c_functionals_.size(); i++) {
        printer->Printf( "    %6.4f %7s", (1.0 - c_alpha_) * c_functionals_[i]->alpha(),
            c_functionals_[i]->name().c_str());
        if (c_functionals_[i]->omega()) {
            printer->Printf( " [omega = %6.4f]", c_functionals_[i]->omega());
        }
        printer->Printf("\n");
    }
    
     // Not currently defined   
    if (c_omega_) {
        printer->Printf( "    %6.4f %7s [omega = %6.4f]\n", (1.0 - c_alpha_), "MP2,LR", c_omega_);
    }
    if (c_alpha_) {
        printer->Printf( "    %6.4f %7s \n", c_alpha_, "MP2");
    }
    if (c_ss_alpha_) {
        printer->Printf( "    %6.4f %s \n", c_ss_alpha_, "Same-Spin SCS-DF-MP2");
    } 
    if (c_os_alpha_) {
        printer->Printf( "    %6.4f %s \n", c_os_alpha_, "Opposite-Spin SCS-DF-MP2");
    } 
    printer->Printf( "\n");

    if (level > 1) {
        for (int i = 0; i < x_functionals_.size(); i++) {
            x_functionals_[i]->print(out,level);
        }
        for (int i = 0; i < c_functionals_.size(); i++) {
            c_functionals_[i]->print(out,level);
        }
    }

    if (dispersion_) {
        dispersion_->print();
    }
}
void SuperFunctional::add_x_functional(boost::shared_ptr<Functional> fun) 
{
    x_functionals_.push_back(fun);
}
void SuperFunctional::add_c_functional(boost::shared_ptr<Functional> fun) 
{
    c_functionals_.push_back(fun);
}
boost::shared_ptr<Functional> SuperFunctional::c_functional(const std::string& name)
{
    for (int Q = 0; Q < c_functionals_.size(); Q++) {
        if (name == c_functionals_[Q]->name())
            return c_functionals_[Q];
    }
    throw PSIEXCEPTION("Functional not found within SuperFunctional");
}
boost::shared_ptr<Functional> SuperFunctional::x_functional(const std::string& name)
{
    for (int Q = 0; Q < x_functionals_.size(); Q++) {
        if (name == x_functionals_[Q]->name())
            return x_functionals_[Q];
    }
    throw PSIEXCEPTION("Functional not found within SuperFunctional");
}
bool SuperFunctional::is_gga() const 
{
    for (int i = 0; i < x_functionals_.size(); i++) {
        if (x_functionals_[i]->is_gga()) 
            return true;
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        if (c_functionals_[i]->is_gga()) 
            return true;
    }
    return false;
}
bool SuperFunctional::is_meta() const 
{
    for (int i = 0; i < x_functionals_.size(); i++) {
        if (x_functionals_[i]->is_meta()) 
            return true;
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        if (c_functionals_[i]->is_meta()) 
            return true;
    }
    return false;
}
void SuperFunctional::partition_gks()
{
    // DFAs need to know about omega
    for (int i = 0; i < x_functionals_.size(); i++) {
        x_functionals_[i]->set_omega(x_omega_);
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        c_functionals_[i]->set_omega(c_omega_);
    }
    // Do nothing for alpha
}
void SuperFunctional::allocate()
{
    values_.clear();

    std::vector<std::string> list;

    // LSDA
    if (deriv_ >= 0) {
        list.push_back("V");
    }   
    if (deriv_ >= 1) {
        list.push_back("V_RHO_A");
        list.push_back("V_RHO_B");
    }   
    if (deriv_ >= 2) {
        list.push_back("V_RHO_A_RHO_A");
        list.push_back("V_RHO_A_RHO_B");
        list.push_back("V_RHO_B_RHO_B");
    }   
 
    // GGA 
    if (is_gga()) {
        if (deriv_ >= 1) {
            list.push_back("V_GAMMA_AA");
            list.push_back("V_GAMMA_AB");
            list.push_back("V_GAMMA_BB");
        }   
        if (deriv_ >= 2) {
            list.push_back("V_GAMMA_AA_GAMMA_AA");
            list.push_back("V_GAMMA_AA_GAMMA_AB");
            list.push_back("V_GAMMA_AA_GAMMA_BB");
            list.push_back("V_GAMMA_AB_GAMMA_AB");
            list.push_back("V_GAMMA_AB_GAMMA_BB");
            list.push_back("V_GAMMA_BB_GAMMA_BB");
        }   
    }
 
    // Meta 
    if (is_meta()) {
        if (deriv_ >= 1) {
            list.push_back("V_TAU_A");
            list.push_back("V_TAU_B");
        }   
        if (deriv_ >= 2) {
            list.push_back("V_TAU_A_TAU_A");
            list.push_back("V_TAU_A_TAU_B");
            list.push_back("V_TAU_B_TAU_B");
        }   
    }

    // LSDA-GGA cross
    if (is_gga()) {
        if (deriv_ >= 2) {
            list.push_back("V_RHO_A_GAMMA_AA");
            list.push_back("V_RHO_A_GAMMA_AB");
            list.push_back("V_RHO_A_GAMMA_BB");
            list.push_back("V_RHO_B_GAMMA_AA");
            list.push_back("V_RHO_B_GAMMA_AB");
            list.push_back("V_RHO_B_GAMMA_BB");
        }   
    }

    // LSDA-Meta cross
    if (is_meta()) {
        if (deriv_ >= 2) {
            list.push_back("V_RHO_A_TAU_A");
            list.push_back("V_RHO_A_TAU_B");
            list.push_back("V_RHO_B_TAU_A");
            list.push_back("V_RHO_B_TAU_B");
        }   
    }

    // GGA-Meta cross
    if (is_gga() && is_meta()) {
        if (deriv_ >= 2) {
            list.push_back("V_GAMMA_AA_TAU_A");
            list.push_back("V_GAMMA_AA_TAU_B");
            list.push_back("V_GAMMA_AB_TAU_A");
            list.push_back("V_GAMMA_AB_TAU_B");
            list.push_back("V_GAMMA_BB_TAU_A");
            list.push_back("V_GAMMA_BB_TAU_B");
        }   
    }

    for (int i = 0; i < list.size(); i++) {
        values_[list[i]] = SharedVector(new Vector(list[i],max_points_));
    }
}
std::map<std::string, SharedVector>& SuperFunctional::compute_functional(const std::map<std::string, SharedVector>& vals, int npoints)
{
    npoints = (npoints == -1 ? vals.find("RHO_A")->second->dimpi()[0] : npoints);
    
    for (std::map<std::string, SharedVector>::const_iterator it = values_.begin();
        it != values_.end(); ++it) {
        ::memset((void*)((*it).second->pointer()),'\0',sizeof(double) * npoints);
    }

    for (int i = 0; i < x_functionals_.size(); i++) {
        x_functionals_[i]->compute_functional(vals, values_, npoints, deriv_, (1.0 - x_alpha_));
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        c_functionals_[i]->compute_functional(vals, values_, npoints, deriv_, (1.0 - c_alpha_));
    }
    
    return values_;
}
void SuperFunctional::test_functional(SharedVector rho_a, 
                                      SharedVector rho_b,
                                      SharedVector gamma_aa,
                                      SharedVector gamma_ab,
                                      SharedVector gamma_bb,
                                      SharedVector tau_a,
                                      SharedVector tau_b)
{
    std::map<std::string, SharedVector> props;
    props["RHO_A"] = rho_a;
    props["RHO_B"] = rho_b;
    props["GAMMA_AA"] = gamma_aa;
    props["GAMMA_AB"] = gamma_ab;
    props["GAMMA_BB"] = gamma_bb;
    props["TAU_A"] = tau_a;
    props["TAU_B"] = tau_b;
    compute_functional(props);
}
SharedVector SuperFunctional::value(const std::string& key) 
{
    return values_[key];
}
int SuperFunctional::ansatz() const 
{
    if (is_meta()) return 2;
    if (is_gga())  return 1;
    return 0;
}

}