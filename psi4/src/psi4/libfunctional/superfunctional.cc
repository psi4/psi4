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

#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "superfunctional.h"
#include "functional.h"
#include "LibXCfunctional.h"
#include <cmath>
using namespace psi;

namespace psi {

SuperFunctional::SuperFunctional() { common_init(); }
SuperFunctional::~SuperFunctional() {}
void SuperFunctional::common_init() {
    max_points_ = 0;
    deriv_ = 0;
    name_ = "";
    description_ = "";
    citation_ = "";

    x_omega_ = 0.0;
    c_omega_ = 0.0;
    x_alpha_ = 0.0;
    x_beta_ = 0.0;
    c_alpha_ = 0.0;

    needs_grac_ = false;
    grac_shift_ = 0.0;
    grac_alpha_ = 0.5;
    grac_beta_ = 40.0;

    needs_vv10_ = false;
    vv10_b_ = 0.0;
    vv10_c_ = 0.0;

    libxc_xc_func_ = false;
    locked_ = false;
}
std::shared_ptr<SuperFunctional> SuperFunctional::blank() {
    return std::shared_ptr<SuperFunctional>(new SuperFunctional());
}
std::shared_ptr<SuperFunctional> SuperFunctional::XC_build(std::string name, bool unpolarized) {
    // Only allow build from full XC kernals
    if (name.find("_XC_") == std::string::npos) {
        throw PSIEXCEPTION("XC_build requires full _XC_ functional names");
    }

    // Build the superfuncitonal
    std::shared_ptr<SuperFunctional> sup = std::shared_ptr<SuperFunctional>(new SuperFunctional());

    // Build LibXC functional
    LibXCFunctional* xc_func = new LibXCFunctional(name, unpolarized);

    // Copy params
    sup->set_name(xc_func->name());
    sup->set_description(xc_func->description());
    sup->set_citation(xc_func->citation());
    sup->set_x_omega(xc_func->omega());
    sup->set_x_alpha(xc_func->global_exchange());
    sup->set_x_beta(xc_func->lr_exchange());
    if (xc_func->needs_vv10()){
        sup->set_vv10_b(xc_func->vv10_b());
        sup->set_vv10_c(xc_func->vv10_c());
    }
    sup->add_c_functional(static_cast<std::shared_ptr<Functional>>(xc_func));
    sup->libxc_xc_func_ = true;

    return sup;
}
std::shared_ptr<SuperFunctional> SuperFunctional::build_worker() {
    // Build the superfuncitonal
    std::shared_ptr<SuperFunctional> sup = std::shared_ptr<SuperFunctional>(new SuperFunctional());

    // Clone over parts
    for (int i = 0; i < x_functionals_.size(); i++) {
        sup->add_x_functional(x_functionals_[i]->build_worker());
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        sup->add_c_functional(c_functionals_[i]->build_worker());
    }

    sup->deriv_ = deriv_;
    sup->max_points_ = max_points_;
    sup->libxc_xc_func_ = libxc_xc_func_;
    if (needs_vv10_) {
        sup->vv10_b_ = vv10_b_;
        sup->vv10_c_ = vv10_c_;
        sup->needs_vv10_ = true;
    }
    if (needs_grac_) {
        sup->needs_grac_ = true;
        sup->grac_shift_ = grac_shift_;
        sup->set_grac_functional(grac_functional_->build_worker());
    }
    sup->allocate();

    return sup;
}
void SuperFunctional::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer =
        (out == "outfile" ? outfile : std::shared_ptr<OutFile>(new OutFile(out)));
    printer->Printf("   => Composite Functional: %s <= \n\n", name_.c_str());

    if (description_ != "") {
        printer->Printf("%s", description_.c_str());
        printer->Printf("\n");
    }

    printer->Printf("%s", citation_.c_str());
    printer->Printf("\n\n");

    printer->Printf("    Deriv            = %14d\n", deriv_);
    printer->Printf("    GGA              = %14s\n", (is_gga() ? "TRUE" : "FALSE"));
    printer->Printf("    Meta             = %14s\n", (is_meta() ? "TRUE" : "FALSE"));
    printer->Printf("\n");

    printer->Printf("    Exchange Hybrid  = %14s\n", (is_x_hybrid() ? "TRUE" : "FALSE"));
    printer->Printf("    Exchange Alpha   = %14.6f\n", x_alpha_);
    printer->Printf("\n");

    printer->Printf("    Exchange LRC     = %14s\n", (is_x_lrc() ? "TRUE" : "FALSE"));
    printer->Printf("    Exchange Beta    = %14.6f\n", x_beta_);
    printer->Printf("    Exchange Omega   = %14.6f\n", x_omega_);
    if (is_c_lrc() || is_c_hybrid()) {
        printer->Printf("\n");
        printer->Printf("    MP2 Hybrid       = %14s\n", (is_c_hybrid() ? "TRUE" : "FALSE"));
        printer->Printf("    MP2 Alpha        = %14.6f\n", c_alpha_);
        printer->Printf("\n");
        printer->Printf("    MP2 LRC          = %14s\n", (is_c_lrc() ? "TRUE" : "FALSE"));
        printer->Printf("    MP2 Omega        = %14.6f\n", c_omega_);
    }
    printer->Printf("\n");

    if (libxc_xc_func_){
        // Well thats nasty
        std::vector<std::tuple<std::string, int, double>> mix_data =
            dynamic_cast<LibXCFunctional*>(c_functionals_[0].get())
                ->get_mix_data();

        int nxc = 0;
        int nexch = 0;
        int ncorr = 0;
        for (int i = 0; i < mix_data.size(); i++) {
            int val = std::get<1>(mix_data[i]);
            if (val == 0) {
                nexch++;
            } else if (val == 1) {
                ncorr++;
            } else if (val == 2) {
                nxc++;
            } else {
                throw PSIEXCEPTION("Functional type not understood");
            }
        }

        if (nxc) {
            printer->Printf("   => Exchange-Correlation Functionals <=\n\n");
            for (int i = 0; i < mix_data.size(); i++) {
                if (std::get<1>(mix_data[i]) != 2) continue;

                printer->Printf("    %6.4f   %7s", std::get<2>(mix_data[i]),
                                std::get<0>(mix_data[i]).c_str());
                printer->Printf("\n");
            }
            printer->Printf("\n");
        }

        if (nexch) {
            printer->Printf("   => Exchange Functionals <=\n\n");
            for (int i = 0; i < mix_data.size(); i++) {
                if (std::get<1>(mix_data[i]) != 0) continue;

                printer->Printf("    %6.4f   %7s", std::get<2>(mix_data[i]),
                                std::get<0>(mix_data[i]).c_str());
                if (c_functionals_[0]->omega()) {
                    printer->Printf(" [omega = %6.4f]", c_functionals_[0]->omega());
                }
                printer->Printf("\n");
            }
            printer->Printf("\n");
        }

        if ((x_omega_ + x_alpha_) > 0.0) {
            printer->Printf("   => Exact (HF) Exchange <=\n\n");
            if (x_omega_) {
                printer->Printf("    %6.4f   %7s [omega = %6.4f]\n", (x_beta_), "HF,LR",
                                x_omega_);
            }
            if (x_alpha_) {
                printer->Printf("    %6.4f   %7s \n", x_alpha_, "HF");
            }
            printer->Printf("\n");
        }

        if (ncorr) {
            printer->Printf("   => Correlation Functionals <=\n\n");
            for (int i = 0; i < mix_data.size(); i++) {
                if (std::get<1>(mix_data[i]) != 1) continue;

                printer->Printf("    %6.4f   %7s", std::get<2>(mix_data[i]),
                                std::get<0>(mix_data[i]).c_str());
                printer->Printf("\n");
            }
            printer->Printf("\n");
        }

    } else {
        printer->Printf("   => Exchange Functionals <=\n\n");
        for (int i = 0; i < x_functionals_.size(); i++) {
            printer->Printf("    %6.4f   %7s", x_functionals_[i]->alpha(),
                            x_functionals_[i]->name().c_str());
            if (x_functionals_[i]->omega()) {
                printer->Printf(" [omega = %6.4f]", x_functionals_[i]->omega());
            }
            printer->Printf("\n");
        }
        printer->Printf("\n");

        printer->Printf("   => Correlation Functionals <=\n\n");
        for (int i = 0; i < c_functionals_.size(); i++) {
            printer->Printf("    %6.4f   %7s", (1.0 - c_alpha_) * c_functionals_[i]->alpha(),
                            c_functionals_[i]->name().c_str());
            if (c_functionals_[i]->omega()) {
                printer->Printf(" [omega = %6.4f]", c_functionals_[i]->omega());
            }
            printer->Printf("\n");
        }
    }

    if (c_omega_) {
        printer->Printf("    %6.4f   %7s [omega = %6.4f]\n", (1.0 - c_alpha_), "MP2,LR", c_omega_);
    }
    if (c_alpha_) {
        printer->Printf("    %6.4f   %7s \n", c_alpha_, "MP2");
    }
    printer->Printf("\n");

    if (needs_grac_) {
        printer->Printf("   => Asymptotic Correction <=\n\n");
        printer->Printf("    Functional       = %14s\n", grac_functional_->name().c_str());
        printer->Printf("    Bulk Shift       = %14.6f\n", grac_shift_);
        printer->Printf("    GRAC Alpha       = %14.6f\n", grac_alpha_);
        printer->Printf("    GRAC Beta        = %14.6f\n", grac_beta_);
        printer->Printf("\n");
    }

    if (needs_vv10_){
        printer->Printf("   => VV10 Non-Local Parameters <=\n\n");
        printer->Printf("    VV10 beta        = %16.4E\n", vv10_b_);
        printer->Printf("    VV10 C           = %16.4E\n", vv10_c_);
        printer->Printf("\n");
    }

    if (level > 1) {
        for (int i = 0; i < x_functionals_.size(); i++) {
            x_functionals_[i]->print(out, level);
        }
        for (int i = 0; i < c_functionals_.size(); i++) {
            c_functionals_[i]->print(out, level);
        }
        if (grac_functional_){
            grac_functional_->print(out, level);
        }
    }
}
void SuperFunctional::can_edit() {
    if (libxc_xc_func_) {
        throw PSIEXCEPTION("Cannot set parameter on full LibXC XC builds\n");
    }
    if (locked_) {
        throw PSIEXCEPTION("The SuperFunctional is locked and cannot be edited.\n");
    }
}
void SuperFunctional::set_x_omega(double omega) {
    can_edit();
    x_omega_ = omega;
}
void SuperFunctional::set_c_omega(double omega) {
    can_edit();
    c_omega_ = omega;
}
void SuperFunctional::set_x_alpha(double alpha) {
    can_edit();
    x_alpha_ = alpha;
}
void SuperFunctional::set_x_beta(double beta) {
    can_edit();
    x_beta_ = beta;
}
void SuperFunctional::set_c_alpha(double alpha) {
    can_edit();
    c_alpha_ = alpha;
}
void SuperFunctional::set_vv10_b(double vv10_b) {
    can_edit();
    needs_vv10_ = true;
    vv10_b_ = vv10_b;
}
void SuperFunctional::set_vv10_c(double vv10_c) {
    can_edit();
    needs_vv10_ = true;
    vv10_c_ = vv10_c;
}
void SuperFunctional::set_grac_alpha(double grac_alpha) {
    can_edit();
    grac_alpha_ = grac_alpha;
}
void SuperFunctional::set_grac_beta(double grac_beta) {
    can_edit();
    grac_beta_ = grac_beta;
}
void SuperFunctional::set_grac_shift(double grac_shift) {
    can_edit();
    if (!grac_functional_){
        throw PSIEXCEPTION("Set the GRAC functional before setting the shift.");
    }
    needs_grac_ = true;
    grac_shift_ = grac_shift;
}
void SuperFunctional::add_x_functional(std::shared_ptr<Functional> fun) {
    can_edit();
    x_functionals_.push_back(fun);
}
void SuperFunctional::add_c_functional(std::shared_ptr<Functional> fun) {
    can_edit();
    c_functionals_.push_back(fun);
}
std::shared_ptr<Functional> SuperFunctional::c_functional(const std::string& name) {
    for (int Q = 0; Q < c_functionals_.size(); Q++) {
        if (name == c_functionals_[Q]->name()) return c_functionals_[Q];
    }
    throw PSIEXCEPTION("Functional not found within SuperFunctional");
}
std::shared_ptr<Functional> SuperFunctional::x_functional(const std::string& name) {
    for (int Q = 0; Q < x_functionals_.size(); Q++) {
        if (name == x_functionals_[Q]->name()) return x_functionals_[Q];
    }
    throw PSIEXCEPTION("Functional not found within SuperFunctional");
}
bool SuperFunctional::is_gga() const {
    for (int i = 0; i < x_functionals_.size(); i++) {
        if (x_functionals_[i]->is_gga())
            return true;
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        if (c_functionals_[i]->is_gga())
            return true;
    }
    if (needs_grac_){
        return true;
    }
    return false;
}
bool SuperFunctional::is_meta() const {
    for (int i = 0; i < x_functionals_.size(); i++) {
        if (x_functionals_[i]->is_meta()) return true;
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        if (c_functionals_[i]->is_meta()) return true;
    }
    return false;
}
bool SuperFunctional::is_unpolarized() const {
    // Need to make sure they are all the same
    std::vector<bool> bool_arr;

    for (int i = 0; i < x_functionals_.size(); i++) {
        bool_arr.push_back(x_functionals_[i]->is_unpolarized());
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        bool_arr.push_back(c_functionals_[i]->is_unpolarized());
    }

    size_t num_true = 0;
    for (int i = 0; i < bool_arr.size(); i++) {
        num_true += bool_arr[i];
    }

    if (num_true == 0) {
        return false;
    } else if (num_true == bool_arr.size()) {
        return true;
    } else {
        outfile->Printf("Mix of polarized and unpolarized functionals detected.\n");
        throw PSIEXCEPTION("All functionals must either be polarized or unpolarized.");
    }
}
void SuperFunctional::allocate() {
    // Make sure were either polarized or not
    bool is_polar = !is_unpolarized();
    values_.clear();

    // Set the lock, after allocation no more changes are allowed
    set_lock(true);

    std::vector<std::string> list;

    // LSDA
    if (deriv_ >= 0) {
        list.push_back("V");
    }
    if (deriv_ >= 1) {
        list.push_back("V_RHO_A");
        if (is_polar){
            list.push_back("V_RHO_B");
        }
    }
    if (deriv_ >= 2) {
        list.push_back("V_RHO_A_RHO_A");

        if (is_polar){
            list.push_back("V_RHO_A_RHO_B");
            list.push_back("V_RHO_B_RHO_B");
        }
    }

    // GGA
    if (is_gga()) {
        if (deriv_ >= 1) {
            list.push_back("V_GAMMA_AA");

            if (is_polar){
                list.push_back("V_GAMMA_AB");
                list.push_back("V_GAMMA_BB");
            }
        }
        if (deriv_ >= 2) {
            list.push_back("V_GAMMA_AA_GAMMA_AA");

            if (is_polar){
                list.push_back("V_GAMMA_AA_GAMMA_AB");
                list.push_back("V_GAMMA_AA_GAMMA_BB");
                list.push_back("V_GAMMA_AB_GAMMA_AB");
                list.push_back("V_GAMMA_AB_GAMMA_BB");
                list.push_back("V_GAMMA_BB_GAMMA_BB");
            }
        }
    }

    // Meta
    if (is_meta()) {
        if (deriv_ >= 1) {
            list.push_back("V_TAU_A");
            // list.push_back("V_LAPL_A");

            if (is_polar){
                list.push_back("V_TAU_B");
                // list.push_back("V_LAPL_B");
            }
        }
        if (deriv_ >= 2) {
            list.push_back("V_TAU_A_TAU_A");

            if (is_polar){
                list.push_back("V_TAU_A_TAU_B");
                list.push_back("V_TAU_B_TAU_B");
            }
        }
    }

    // LSDA-GGA cross
    if (is_gga()) {
        if (deriv_ >= 2) {
            list.push_back("V_RHO_A_GAMMA_AA");

            if (is_polar){
                list.push_back("V_RHO_A_GAMMA_AB");
                list.push_back("V_RHO_A_GAMMA_BB");
                list.push_back("V_RHO_B_GAMMA_AA");
                list.push_back("V_RHO_B_GAMMA_AB");
                list.push_back("V_RHO_B_GAMMA_BB");
            }
        }
    }

    // LSDA-Meta cross
    if (is_meta()) {
        if (deriv_ >= 2) {
            list.push_back("V_RHO_A_TAU_A");

            if (is_polar){
                list.push_back("V_RHO_A_TAU_B");
                list.push_back("V_RHO_B_TAU_A");
                list.push_back("V_RHO_B_TAU_B");
            }
        }
    }

    // GGA-Meta cross
    if (is_gga() && is_meta()) {
        if (deriv_ >= 2) {
            list.push_back("V_GAMMA_AA_TAU_A");
            if (is_polar){
                list.push_back("V_GAMMA_AA_TAU_B");
                list.push_back("V_GAMMA_AB_TAU_A");
                list.push_back("V_GAMMA_AB_TAU_B");
                list.push_back("V_GAMMA_BB_TAU_A");
                list.push_back("V_GAMMA_BB_TAU_B");
            }
        }
    }

    for (int i = 0; i < list.size(); i++) {
        values_[list[i]] = SharedVector(new Vector(list[i], max_points_));
    }

    if (needs_grac_) {
        ac_values_["V_RHO_A"] = SharedVector(new Vector("V_RHO_A", max_points_));
        ac_values_["V_GAMMA_AA"] = SharedVector(new Vector("V_GAMMA_AA", max_points_));
        if (is_polar) {
            throw PSIEXCEPTION("GRAC is not implemented for UKS functionals.");
            ac_values_["V_RHO_B"] = SharedVector(new Vector("V_RHO_B", max_points_));
            ac_values_["V_GAMMA_AB"] = SharedVector(new Vector("V_GAMMA_AB", max_points_));
            ac_values_["V_GAMMA_BB"] = SharedVector(new Vector("V_GAMMA_BB", max_points_));
        }
    }
}
std::map<std::string, SharedVector>& SuperFunctional::compute_functional(
    const std::map<std::string, SharedVector>& vals, int npoints) {
    npoints = (npoints == -1 ? vals.find("RHO_A")->second->dimpi()[0] : npoints);

    for (std::map<std::string, SharedVector>::const_iterator it = values_.begin();
         it != values_.end(); ++it) {
        ::memset((void*)((*it).second->pointer()), '\0', sizeof(double) * npoints);
    }

    for (int i = 0; i < x_functionals_.size(); i++) {
        x_functionals_[i]->compute_functional(vals, values_, npoints, deriv_);
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        c_functionals_[i]->compute_functional(vals, values_, npoints, deriv_);
    }

    // Apply the grac shift, only valid for gradient computations
    if (needs_grac_ && (deriv_ == 1)) {
        for (std::map<std::string, SharedVector>::const_iterator it = ac_values_.begin();
             it != ac_values_.end(); ++it) {
            ::memset((void*)((*it).second->pointer()), '\0', sizeof(double) * npoints);
        }
        grac_functional_->compute_functional(vals, ac_values_, npoints, 1);

        if (is_unpolarized()){

            double* rho = vals.find("RHO_A")->second->pointer();
            double* rho_x = vals.find("RHO_AX")->second->pointer();
            double* rho_y = vals.find("RHO_AY")->second->pointer();
            double* rho_z = vals.find("RHO_AZ")->second->pointer();

            double* v_rho = values_["V_RHO_A"]->pointer();
            double* v_gamma = values_["V_GAMMA_AA"]->pointer();

            double* grac_v_rho = ac_values_["V_RHO_A"]->pointer();
            double* grac_v_gamma = ac_values_["V_GAMMA_AA"]->pointer();

            const double galpha = -1.0 * grac_alpha_;
            const double gbeta = grac_beta_;
            const double gshift = grac_shift_;
            const double pow43 = 4.0 / 3.0;
            double denx;

            # pragma omp simd
            for (size_t i = 0; i < npoints; i++){

                // Gotta be careful of this, specially for CP results
                if (rho[i] < 1.e-14){
                    denx = 1.e8; // Will force grac_fx to 1
                } else {
                    denx = std::fabs(rho_x[i] + rho_y[i] + rho_z[i]) / std::pow(rho[i], pow43);
                }

                double grac_fx = 1.0 / (1.0 + std::exp(galpha * (denx - gbeta)));

                v_rho[i] = ((1.0 - grac_fx) * (v_rho[i] - gshift)) + (grac_fx * grac_v_rho[i]);
                v_gamma[i] = (1.0 - grac_fx) * v_gamma[i];

                // We neglect the gradient of denx with v_gamma, as it requires the laplacian and
                // there is virtually no difference. See DOI: 10.1002/cphc.200700504
            }
        }


        // This is turned off by allocate for now, this doesnt appear to be quite correct.
        else{

            throw PSIEXCEPTION("GRAC is not implemented for UKS functionals.");
            double* rho_a = vals.find("RHO_A")->second->pointer();
            double* rho_a_x = vals.find("RHO_AX")->second->pointer();
            double* rho_a_y = vals.find("RHO_AY")->second->pointer();
            double* rho_a_z = vals.find("RHO_AZ")->second->pointer();

            double* rho_b = vals.find("RHO_B")->second->pointer();
            double* rho_b_x = vals.find("RHO_BX")->second->pointer();
            double* rho_b_y = vals.find("RHO_BY")->second->pointer();
            double* rho_b_z = vals.find("RHO_BZ")->second->pointer();

            double* v_rho_a = values_["V_RHO_A"]->pointer();
            double* v_rho_b = values_["V_RHO_B"]->pointer();
            double* v_gamma_aa = values_["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = values_["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = values_["V_GAMMA_BB"]->pointer();

            double* grac_v_rho_a = ac_values_["V_RHO_A"]->pointer();
            double* grac_v_rho_b = ac_values_["V_RHO_B"]->pointer();
            double* grac_v_gamma_aa = ac_values_["V_GAMMA_AA"]->pointer();
            double* grac_v_gamma_ab = ac_values_["V_GAMMA_AB"]->pointer();
            double* grac_v_gamma_bb = ac_values_["V_GAMMA_BB"]->pointer();

            const double galpha = -1.0 * grac_alpha_;
            const double gbeta = grac_beta_;
            const double gshift = grac_shift_;

            # pragma omp simd
            for (size_t i = 0; i < npoints; i++){
                double denx = std::fabs(rho_a_x[i] + rho_a_y[i] + rho_a_z[i] + rho_b_x[i] +
                                          rho_b_y[i] + rho_b_z[i]) /
                                std::pow((rho_a[i] + rho_b[i]), 4.0 / 3.0);
                double grac_fx = 1.0 / (1.0 + std::exp(galpha * (denx - gbeta)));

                v_rho_a[i] = (1.0 - grac_fx) * (v_rho_a[i] - gshift) + (grac_fx * grac_v_rho_a[i]);
                v_rho_b[i] = (1.0 - grac_fx) * (v_rho_b[i] - gshift) + (grac_fx * grac_v_rho_b[i]);

                v_gamma_aa[i] = ((1.0 - grac_fx) * v_gamma_aa[i]) + (grac_fx * grac_v_gamma_aa[i]);
                v_gamma_ab[i] = ((1.0 - grac_fx) * v_gamma_ab[i]) + (grac_fx * grac_v_gamma_ab[i]);
                v_gamma_bb[i] = ((1.0 - grac_fx) * v_gamma_bb[i]) + (grac_fx * grac_v_gamma_bb[i]);
            }
        }
    }

    return values_;
}
void SuperFunctional::test_functional(SharedVector rho_a, SharedVector rho_b, SharedVector gamma_aa,
                                      SharedVector gamma_ab, SharedVector gamma_bb,
                                      SharedVector tau_a, SharedVector tau_b) {
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
