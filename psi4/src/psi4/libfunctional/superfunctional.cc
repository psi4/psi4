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

#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "superfunctional.h"
#include "functional.h"
#include "LibXCfunctional.h"

#include <cmath>
#include <cstdlib>

// using namespace psi;

namespace psi {

SuperFunctional::SuperFunctional() { common_init(); }
SuperFunctional::~SuperFunctional() {}
void SuperFunctional::common_init() {
    max_points_ = 0;
    deriv_ = 0;
    name_ = "";
    description_ = "";
    citation_ = "";
    xclib_description_ = "";

    x_omega_ = 0.0;
    c_omega_ = 0.0;
    x_alpha_ = 0.0;
    x_beta_ = 1.0;
    c_alpha_ = 0.0;
    c_os_alpha_ = 0.0;
    c_ss_alpha_ = 0.0;

    needs_grac_ = false;
    grac_shift_ = 0.0;
    grac_alpha_ = 0.5;
    grac_beta_ = 40.0;

    needs_vv10_ = false;
    vv10_b_ = 0.0;
    vv10_c_ = 0.0;
    vv10_beta_ = 0.0;

    libxc_xc_func_ = false;
    locked_ = false;
    density_tolerance_ = 0.0;
}
std::shared_ptr<SuperFunctional> SuperFunctional::blank() { return std::make_shared<SuperFunctional>(); }
std::shared_ptr<SuperFunctional> SuperFunctional::XC_build(std::string name, bool unpolarized, const std::optional<std::map<std::string, double>>& tweakers_) {
    // Only allow build from full XC kernels
    if (name.find("XC_") == std::string::npos) {
        throw PSIEXCEPTION("XC_build requires full XC_ functional names");
    }

    // Build the superfuncitonal
    auto sup = std::make_shared<SuperFunctional>();

    // Build LibXC functional
    LibXCFunctional* xc_func = new LibXCFunctional(name, unpolarized);

    // Copy params
    sup->set_name(xc_func->name());
    sup->set_description(xc_func->description());
    sup->set_citation(xc_func->citation());
    sup->set_xclib_description(xc_func->xclib_description());
    sup->set_x_omega(xc_func->omega());
    sup->set_x_alpha(xc_func->global_exchange());
    sup->set_x_beta(xc_func->lr_exchange());

    // User tweakers
    if (tweakers_.value().size()) {
        xc_func->set_tweak(tweakers_.value(), true);
    }

    if (xc_func->needs_vv10()) {
        sup->set_vv10_b(xc_func->vv10_b());
        sup->set_vv10_c(xc_func->vv10_c());
    }
    sup->add_c_functional(static_cast<std::shared_ptr<Functional>>(xc_func));
    sup->libxc_xc_func_ = true;

    return sup;
}
std::shared_ptr<SuperFunctional> SuperFunctional::build_polarized() {
    // Build the superfunctional
    auto sup = std::make_shared<SuperFunctional>();

    // Clone over parts
    for (int i = 0; i < x_functionals_.size(); i++) {
        sup->add_x_functional(x_functionals_[i]->build_polarized());
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        sup->add_c_functional(c_functionals_[i]->build_polarized());
    }

    sup->deriv_ = deriv_;
    sup->max_points_ = max_points_;
    sup->libxc_xc_func_ = libxc_xc_func_;
    if (needs_vv10_) {
        sup->needs_vv10_ = true;
        sup->vv10_b_ = vv10_b_;
        sup->vv10_c_ = vv10_c_;
        sup->vv10_beta_ = vv10_beta_;
    }
    if (needs_grac_) {
        sup->needs_grac_ = true;
        sup->grac_shift_ = grac_shift_;
        sup->grac_alpha_ = grac_alpha_;
        sup->grac_beta_ = grac_beta_;
        sup->set_grac_x_functional(grac_x_functional_->build_polarized());
        sup->set_grac_c_functional(grac_c_functional_->build_polarized());
    }
    sup->name_ = name_;
    sup->description_ = description_;
    sup->citation_ = citation_;
    sup->xclib_description_ = xclib_description_;
    sup->x_omega_ = x_omega_;
    sup->c_omega_ = c_omega_;
    sup->x_alpha_ = x_alpha_;
    sup->x_beta_ = x_beta_;
    sup->c_alpha_ = c_alpha_;
    sup->c_ss_alpha_ = c_ss_alpha_;
    sup->c_os_alpha_ = c_os_alpha_;
    sup->density_tolerance_ = density_tolerance_;
    sup->allocate();

    return sup;
}
std::shared_ptr<SuperFunctional> SuperFunctional::build_worker() {
    // Build the superfunctional
    auto sup = std::make_shared<SuperFunctional>();

    // Clone over parts
    for (int i = 0; i < x_functionals_.size(); i++) {
        sup->add_x_functional(x_functionals_[i]->build_worker());
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        sup->add_c_functional(c_functionals_[i]->build_worker());
    }

    // Workers dont need omega or alpha
    sup->deriv_ = deriv_;
    sup->max_points_ = max_points_;
    sup->libxc_xc_func_ = libxc_xc_func_;
    if (needs_vv10_) {
        sup->needs_vv10_ = true;
        sup->vv10_b_ = vv10_b_;
        sup->vv10_c_ = vv10_c_;
        sup->vv10_beta_ = vv10_beta_;
    }
    if (needs_grac_) {
        sup->needs_grac_ = true;
        sup->grac_shift_ = grac_shift_;
        sup->grac_alpha_ = grac_alpha_;
        sup->grac_beta_ = grac_beta_;
        sup->set_grac_x_functional(grac_x_functional_->build_worker());
        sup->set_grac_c_functional(grac_c_functional_->build_worker());
    }
    sup->allocate();

    return sup;
}
void SuperFunctional::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));

    if (xclib_description_ != "") {
        printer->Printf("%s\n\n", xclib_description_.c_str());
    }

    printer->Printf("   => Composite Functional: %s <= \n\n", name_.c_str());

    if (description_ != "") {
        printer->Printf("%s", description_.c_str());
        printer->Printf("\n");
    }

    printer->Printf("%s", citation_.c_str());
    printer->Printf("\n\n");

    printer->Printf("    Deriv               = %14d\n", deriv_);
    printer->Printf("    GGA                 = %14s\n", (is_gga() ? "TRUE" : "FALSE"));
    printer->Printf("    Meta                = %14s\n", (is_meta() ? "TRUE" : "FALSE"));
    printer->Printf("\n");

    printer->Printf("    Exchange Hybrid     = %14s\n", (is_x_hybrid() ? "TRUE" : "FALSE"));
    printer->Printf("    MP2 Hybrid          = %14s\n",
                    ((is_c_lrc() || is_c_hybrid() || is_c_scs_hybrid()) ? "TRUE" : "FALSE"));
    printer->Printf("\n");

    if (libxc_xc_func_) {
        // Well thats nasty
        std::vector<std::tuple<std::string, int, double>> mix_data =
            dynamic_cast<LibXCFunctional*>(c_functionals_[0].get())->get_mix_data();

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

                printer->Printf("    %6.4f   %14s", std::get<2>(mix_data[i]), std::get<0>(mix_data[i]).c_str());
                printer->Printf("\n");
            }
            printer->Printf("\n");
        }

        if (nexch) {
            printer->Printf("   => Exchange Functionals <=\n\n");
            for (int i = 0; i < mix_data.size(); i++) {
                if (std::get<1>(mix_data[i]) != 0) continue;

                printer->Printf("    %6.4f   %14s", std::get<2>(mix_data[i]), std::get<0>(mix_data[i]).c_str());
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
                printer->Printf("    %6.4f   %14s [omega = %6.4f]\n", (x_beta_), "HF,LR", x_omega_);
            }
            if (x_alpha_) {
                printer->Printf("    %6.4f   %14s \n", x_alpha_, "HF");
            }
            printer->Printf("\n");
        }

        if (ncorr) {
            printer->Printf("   => Correlation Functionals <=\n\n");
            for (int i = 0; i < mix_data.size(); i++) {
                if (std::get<1>(mix_data[i]) != 1) continue;

                printer->Printf("    %6.4f   %14s", std::get<2>(mix_data[i]), std::get<0>(mix_data[i]).c_str());
                printer->Printf("\n");
            }
            printer->Printf("\n");
        }

    } else {
        printer->Printf("   => Exchange Functionals <=\n\n");
        for (int i = 0; i < x_functionals_.size(); i++) {
            printer->Printf("    %6.4f   %14s", x_functionals_[i]->alpha(), x_functionals_[i]->name().c_str());
            if (x_functionals_[i]->omega()) {
                printer->Printf(" [omega = %6.4f]", x_functionals_[i]->omega());
            }
            printer->Printf("\n");
        }
        printer->Printf("\n");

        if ((x_omega_ + x_alpha_) > 0.0) {
            printer->Printf("   => Exact (HF) Exchange <=\n\n");
            if (x_omega_) {
                printer->Printf("    %6.4f   %14s [omega = %6.4f]\n", (x_beta_), "HF,LR", x_omega_);
            }
            if (x_alpha_) {
                printer->Printf("    %6.4f   %14s \n", x_alpha_, "HF");
            }
            printer->Printf("\n");
        }

        printer->Printf("   => Correlation Functionals <=\n\n");
        for (int i = 0; i < c_functionals_.size(); i++) {
            printer->Printf("    %6.4f   %14s", c_functionals_[i]->alpha(), c_functionals_[i]->name().c_str());
            if (c_functionals_[i]->omega()) {
                printer->Printf(" [omega = %6.4f]", c_functionals_[i]->omega());
            }
            printer->Printf("\n");
        }
        printer->Printf("\n");
    }

    if (is_c_lrc() || is_c_hybrid() || is_c_scs_hybrid()) {
        printer->Printf("   => MP2 Correlation <=\n\n");
        if (c_omega_) {
            printer->Printf("    %6.4f   %7s [omega = %6.4f]\n", (1.0 - c_alpha_), "MP2,LR", c_omega_);
        }
        if (c_alpha_) {
            if (is_c_scs_hybrid()) {
                printer->Printf("    MP2 SCS Hybrid      = %14s\n", "TRUE");
                printer->Printf("    MP2 OS Alpha        = %14.6f\n", c_os_alpha_);
                printer->Printf("    MP2 SS Alpha        = %14.6f\n", c_ss_alpha_);
            } else {
                printer->Printf("    %6.4f   %7s \n", c_alpha_, "MP2");
            }
            printer->Printf("\n");
        }
    }

    print_density_threshold("outfile", 1);

    if (needs_grac_) {
        printer->Printf("   => Asymptotic Correction <=\n\n");
        printer->Printf("    X Functional        = %14s\n", grac_x_functional_->name().c_str());
        printer->Printf("    C Functional        = %14s\n", grac_c_functional_->name().c_str());
        printer->Printf("    Bulk Shift          = %14.6f\n", grac_shift_);
        printer->Printf("    GRAC Alpha          = %14.6f\n", grac_alpha_);
        printer->Printf("    GRAC Beta           = %14.6f\n", grac_beta_);
        printer->Printf("\n");
    }

    if (needs_vv10_) {
        printer->Printf("   => VV10 Non-Local Parameters <=\n\n");
        printer->Printf("    VV10 B              = %14.4E\n", vv10_b_);
        printer->Printf("    VV10 C              = %14.4E\n", vv10_c_);
        printer->Printf("\n");
    }

    if (level > 1) {
        for (int i = 0; i < x_functionals_.size(); i++) {
            x_functionals_[i]->print(out, level);
        }
        for (int i = 0; i < c_functionals_.size(); i++) {
            c_functionals_[i]->print(out, level);
        }
        if (needs_grac_) {
            if (grac_x_functional_) {
                grac_x_functional_->print(out, level);
            }
            if (grac_c_functional_) {
                grac_c_functional_->print(out, level);
            }
        }
    }
}
void SuperFunctional::can_edit() {
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
    vv10_beta_ = 1.0 / 32.0 * std::pow((3.0 / (vv10_b_ * vv10_b_)), (3.0 / 4.0));
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
    if (!grac_x_functional_) {
        throw PSIEXCEPTION("Set the GRAC functional before setting the shift.");
    }
    needs_grac_ = true;
    grac_shift_ = grac_shift;
}
void SuperFunctional::set_c_ss_alpha(double alpha) {
    can_edit();
    c_ss_alpha_ = alpha;
}
void SuperFunctional::set_c_os_alpha(double alpha) {
    can_edit();
    c_os_alpha_ = alpha;
}
void SuperFunctional::add_x_functional(std::shared_ptr<Functional> fun) {
    can_edit();
    x_functionals_.push_back(fun);
}
void SuperFunctional::add_c_functional(std::shared_ptr<Functional> fun) {
    can_edit();
    c_functionals_.push_back(fun);
}
void SuperFunctional::set_density_tolerance(double cut) {
    density_tolerance_ = cut;
    for (int Q = 0; Q < c_functionals_.size(); Q++) {
        c_functionals_[Q]->set_density_cutoff(density_tolerance_);
    };
    for (int Q = 0; Q < x_functionals_.size(); Q++) {
        x_functionals_[Q]->set_density_cutoff(density_tolerance_);
    };
}
void SuperFunctional::print_density_threshold(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => LibXC Density Thresholds  <==\n\n");
    double val = 0.0;
    for (int Q = 0; Q < c_functionals_.size(); Q++) {
        val = c_functionals_[Q]->query_density_cutoff();
        printer->Printf("    %s:  %6.2E \n", c_functionals_[Q]->name().c_str(), val);
    };
    for (int Q = 0; Q < x_functionals_.size(); Q++) {
        val = x_functionals_[Q]->query_density_cutoff();
        printer->Printf("    %s:  %6.2E \n", x_functionals_[Q]->name().c_str(), val);
    };
    printer->Printf("\n");
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
        if (x_functionals_[i]->is_gga()) return true;
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        if (c_functionals_[i]->is_gga()) return true;
    }
    if (needs_grac_ || needs_vv10_) {
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

    // Temporaries
    list.push_back("Q_TMP");

    // LSDA
    if (deriv_ >= 0) {
        list.push_back("V");
    }
    if (deriv_ >= 1) {
        list.push_back("V_RHO_A");
        if (is_polar) {
            list.push_back("V_RHO_B");
        }
    }
    if (deriv_ >= 2) {
        list.push_back("V_RHO_A_RHO_A");

        if (is_polar) {
            list.push_back("V_RHO_A_RHO_B");
            list.push_back("V_RHO_B_RHO_B");
        }
    }

    // GGA
    if (is_gga()) {
        if (deriv_ >= 1) {
            list.push_back("V_GAMMA_AA");

            if (is_polar) {
                list.push_back("V_GAMMA_AB");
                list.push_back("V_GAMMA_BB");
            }
        }
        if (deriv_ >= 2) {
            list.push_back("V_GAMMA_AA_GAMMA_AA");

            if (is_polar) {
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

            if (is_polar) {
                list.push_back("V_TAU_B");
                // list.push_back("V_LAPL_B");
            }
        }
        if (deriv_ >= 2) {
            list.push_back("V_TAU_A_TAU_A");

            if (is_polar) {
                list.push_back("V_TAU_A_TAU_B");
                list.push_back("V_TAU_B_TAU_B");
            }
        }
    }

    // LSDA-GGA cross
    if (is_gga()) {
        if (deriv_ >= 2) {
            list.push_back("V_RHO_A_GAMMA_AA");

            if (is_polar) {
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

            if (is_polar) {
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
            if (is_polar) {
                list.push_back("V_GAMMA_AA_TAU_B");
                list.push_back("V_GAMMA_AB_TAU_A");
                list.push_back("V_GAMMA_AB_TAU_B");
                list.push_back("V_GAMMA_BB_TAU_A");
                list.push_back("V_GAMMA_BB_TAU_B");
            }
        }
    }

    for (int i = 0; i < list.size(); i++) {
        values_[list[i]] = std::make_shared<Vector>(list[i], max_points_);
    }

    if (needs_grac_) {
        ac_values_["V"] = std::make_shared<Vector>("V", max_points_);  // Not actually used
        ac_values_["V_RHO_A"] = std::make_shared<Vector>("V_RHO_A", max_points_);
        ac_values_["V_GAMMA_AA"] = std::make_shared<Vector>("V_GAMMA_AA", max_points_);
        if (is_polar) {
            throw PSIEXCEPTION("GRAC is not implemented for UKS functionals.");
            ac_values_["V_RHO_B"] = std::make_shared<Vector>("V_RHO_B", max_points_);
            ac_values_["V_GAMMA_AB"] = std::make_shared<Vector>("V_GAMMA_AB", max_points_);
            ac_values_["V_GAMMA_BB"] = std::make_shared<Vector>("V_GAMMA_BB", max_points_);
        }
    }

    if (needs_vv10_) {
        vv_values_["W0"] = std::make_shared<Vector>("W0", max_points_);
        vv_values_["KAPPA"] = std::make_shared<Vector>("KAPPA", max_points_);
        vv_values_["GRID_WX"] = std::make_shared<Vector>("W_X_GRID", max_points_);
        vv_values_["GRID_WY"] = std::make_shared<Vector>("W_Y_GRID", max_points_);
        vv_values_["GRID_WZ"] = std::make_shared<Vector>("W_Z_GRID", max_points_);
    }
}
std::map<std::string, SharedVector>& SuperFunctional::compute_functional(
    const std::map<std::string, SharedVector>& vals, int npoints, bool singlet) {
    if (!singlet && is_unpolarized() and deriv_ == 2) {
        // We want triplet spin integration for the second derivatives.

        // First, make a UKS version of this functional.
        auto UKS = build_polarized();

        // Now compute the UKS version of this functional. Naturally, the input values will need to be UKS-ified first.
        std::map<std::string, SharedVector> UKS_vals;
        if (true) {
            UKS_vals["RHO_A"] = std::make_shared<Vector>(std::move((vals.at("RHO_A")->clone())));
            UKS_vals["RHO_A"]->scale(0.5); // Un-spinsum
            UKS_vals["RHO_B"] = UKS_vals["RHO_A"];
        }
        if (vals.count("RHO_AX")) {
            UKS_vals["RHO_AX"] = std::make_shared<Vector>(std::move(vals.at("RHO_AX")->clone()));
            UKS_vals["RHO_AX"]->scale(0.5); // Un-spinsum
            UKS_vals["RHO_BX"] = UKS_vals["RHO_AX"];
            UKS_vals["RHO_AY"] = std::make_shared<Vector>(std::move(vals.at("RHO_AY")->clone()));
            UKS_vals["RHO_AY"]->scale(0.5); // Un-spinsum
            UKS_vals["RHO_BY"] = UKS_vals["RHO_AY"];
            UKS_vals["RHO_AZ"] = std::make_shared<Vector>(std::move(vals.at("RHO_AZ")->clone()));
            UKS_vals["RHO_AZ"]->scale(0.5); // Un-spinsum
            UKS_vals["RHO_BZ"] = UKS_vals["RHO_AZ"];
            UKS_vals["GAMMA_AA"] = std::make_shared<Vector>(std::move(vals.at("GAMMA_AA")->clone()));
            UKS_vals["GAMMA_AA"]->scale(0.25); // Un-spinsum
            UKS_vals["GAMMA_AB"] = UKS_vals["GAMMA_AA"];
            UKS_vals["GAMMA_BB"] = UKS_vals["GAMMA_AA"];
        }

        auto temp = UKS->compute_functional(UKS_vals, npoints);
        values_ = std::move(UKS->compute_functional(UKS_vals, npoints));

        // Now we take the magic triplet combinations.
        if (true) {
            values_.erase("V_RHO_B_RHO_B");
            values_.at("V_RHO_A_RHO_A")->axpby(-0.5, 0.5, *values_.at("V_RHO_A_RHO_B"));
            values_.erase("V_RHO_A_RHO_B");
        }
        if (vals.count("RHO_AX")) {
            values_.erase("V_GAMMA_BB");
            values_.at("V_GAMMA_AA")->axpby(-0.25, 0.5, *values_.at("V_GAMMA_AB"));
            values_.erase("V_GAMMA_AB");

            values_.erase("V_RHO_A_GAMMA_AB");
            values_.erase("V_RHO_B_GAMMA_AA");
            values_.erase("V_RHO_B_GAMMA_AB");
            values_.erase("V_RHO_B_GAMMA_BB");
            values_.at("V_RHO_A_GAMMA_AA")->axpby(-0.25, 0.25, *values_.at("V_RHO_A_GAMMA_BB"));
            values_.erase("V_RHO_A_GAMMA_BB");

            values_.erase("V_GAMMA_AA_GAMMA_AB");
            values_.erase("V_GAMMA_AB_GAMMA_AB");
            values_.erase("V_GAMMA_AB_GAMMA_BB");
            values_.erase("V_GAMMA_BB_GAMMA_BB");
            values_.at("V_GAMMA_AA_GAMMA_AA")->axpby(-0.125, 0.125, *values_.at("V_GAMMA_AA_GAMMA_BB"));
            values_.erase("V_GAMMA_AA_GAMMA_BB");
        }
        return values_;
    }
    npoints = (npoints == -1 ? vals.find("RHO_A")->second->dimpi()[0] : npoints);

    // Zero out values
    for (auto kv : values_) {
        kv.second->zero();
    }

    for (int i = 0; i < x_functionals_.size(); i++) {
        x_functionals_[i]->compute_functional(vals, values_, npoints, deriv_);
    }
    for (int i = 0; i < c_functionals_.size(); i++) {
        c_functionals_[i]->compute_functional(vals, values_, npoints, deriv_);
    }

    // Apply the grac shift, only valid for gradient computations
    if (needs_grac_ && (deriv_ == 1)) {
        for (std::map<std::string, SharedVector>::const_iterator it = ac_values_.begin(); it != ac_values_.end();
             ++it) {
            ::memset((void*)((*it).second->pointer()), '\0', sizeof(double) * npoints);
        }
        if (grac_x_functional_) {
            grac_x_functional_->compute_functional(vals, ac_values_, npoints, 1);
        }
        if (grac_c_functional_) {
            grac_c_functional_->compute_functional(vals, ac_values_, npoints, 1);
        }

        if (is_unpolarized()) {
            double* rho = vals.find("RHO_A")->second->pointer();
            double* sigma = vals.find("GAMMA_AA")->second->pointer();

            double* v_rho = values_["V_RHO_A"]->pointer();
            double* v_gamma = values_["V_GAMMA_AA"]->pointer();

            double* grac_v_rho = ac_values_["V_RHO_A"]->pointer();
            double* grac_v_gamma = ac_values_["V_GAMMA_AA"]->pointer();

            const double galpha = -1.0 * grac_alpha_;
            const double gbeta = grac_beta_;
            const double gshift = grac_shift_;
            const double pow43 = 4.0 / 3.0;
            double denx;

#pragma omp simd
            for (size_t i = 0; i < npoints; i++) {
                if (rho[i] < 1.e-16) {
                    denx = 1.e2;  // Will force grac_fx to 1
                } else {
                    denx = std::pow(sigma[i], 0.5) / std::pow(rho[i], pow43);
                }

                double grac_fx = 1.0 / (1.0 + std::exp(galpha * (denx - gbeta)));

                double sr_grac_fx = (1.0 - grac_fx);
                v_rho[i] = sr_grac_fx * (v_rho[i] - gshift) + (grac_fx * grac_v_rho[i]);
                v_gamma[i] = sr_grac_fx * v_gamma[i];

                // We neglect the gradient of denx with v_gamma, as it requires the laplacian and
                // there is virtually no difference. See DOI: 10.1002/cphc.200700504
            }
        }

        // This is turned off by allocate for now, this doesnt appear to be quite correct.
        else {
            throw PSIEXCEPTION("GRAC is not implemented for UKS functionals.");
            // double* rho_a = vals.find("RHO_A")->second->pointer();
            // double* rho_a_x = vals.find("RHO_AX")->second->pointer();
            // double* rho_a_y = vals.find("RHO_AY")->second->pointer();
            // double* rho_a_z = vals.find("RHO_AZ")->second->pointer();

            // double* rho_b = vals.find("RHO_B")->second->pointer();
            // double* rho_b_x = vals.find("RHO_BX")->second->pointer();
            // double* rho_b_y = vals.find("RHO_BY")->second->pointer();
            // double* rho_b_z = vals.find("RHO_BZ")->second->pointer();

            // double* v_rho_a = values_["V_RHO_A"]->pointer();
            // double* v_rho_b = values_["V_RHO_B"]->pointer();
            // double* v_gamma_aa = values_["V_GAMMA_AA"]->pointer();
            // double* v_gamma_ab = values_["V_GAMMA_AB"]->pointer();
            // double* v_gamma_bb = values_["V_GAMMA_BB"]->pointer();

            // double* grac_v_rho_a = ac_values_["V_RHO_A"]->pointer();
            // double* grac_v_rho_b = ac_values_["V_RHO_B"]->pointer();
            // double* grac_v_gamma_aa = ac_values_["V_GAMMA_AA"]->pointer();
            // double* grac_v_gamma_ab = ac_values_["V_GAMMA_AB"]->pointer();
            // double* grac_v_gamma_bb = ac_values_["V_GAMMA_BB"]->pointer();

            // const double galpha = -1.0 * grac_alpha_;
            // const double gbeta = grac_beta_;
            // const double gshift = grac_shift_;

            // # pragma omp simd
            // for (size_t i = 0; i < npoints; i++){
            //     double denx = std::fabs(rho_a_x[i] + rho_a_y[i] + rho_a_z[i] + rho_b_x[i] +
            //                               rho_b_y[i] + rho_b_z[i]) /
            //                     std::pow((rho_a[i] + rho_b[i]), 4.0 / 3.0);
            //     double grac_fx = 1.0 / (1.0 + std::exp(galpha * (denx - gbeta)));

            //     v_rho_a[i] = (1.0 - grac_fx) * (v_rho_a[i] - gshift) + (grac_fx * grac_v_rho_a[i]);
            //     v_rho_b[i] = (1.0 - grac_fx) * (v_rho_b[i] - gshift) + (grac_fx * grac_v_rho_b[i]);

            //     v_gamma_aa[i] = ((1.0 - grac_fx) * v_gamma_aa[i]) + (grac_fx * grac_v_gamma_aa[i]);
            //     v_gamma_ab[i] = ((1.0 - grac_fx) * v_gamma_ab[i]) + (grac_fx * grac_v_gamma_ab[i]);
            //     v_gamma_bb[i] = ((1.0 - grac_fx) * v_gamma_bb[i]) + (grac_fx * grac_v_gamma_bb[i]);
            // }
        }
    }

    return values_;
}
std::map<std::string, SharedVector> SuperFunctional::compute_vv10_cache(const std::map<std::string, SharedVector>& vals,
                                                                        std::shared_ptr<BlockOPoints> block,
                                                                        double rho_thresh, int npoints, bool internal) {
    npoints = (npoints == -1 ? vals.find("RHO_A")->second->dimpi()[0] : npoints);

    // Precompute prefactors
    const double Wp_pref = (4.0 / 3.0) * M_PI;
    const double Wg_pref = vv10_c_;
    const double kappa_pref = (1.5 * vv10_b_ * M_PI) / std::pow(9.0 * M_PI, 1.0 / 6.0);
    const double kappa_pow = 1.0 / 6.0;

    size_t nact = 0;
    // rho_thresh = 0.0;

    // Zero out the vectors
    vv_values_["W0"]->zero();
    vv_values_["KAPPA"]->zero();

    double* w0p = vv_values_["W0"]->pointer();
    double* kappap = vv_values_["KAPPA"]->pointer();
    double* rhop = vals.find("RHO_A")->second->pointer();
    double* gammap = vals.find("GAMMA_AA")->second->pointer();

    // Eh, worth a shot
    // clang 10 on Mac objects at runtime: #pragma omp simd
    for (size_t i = 0; i < npoints; i++) {
        if (rhop[i] < rho_thresh) continue;

        double Wp = Wp_pref * rhop[i];

        double Wg_tmp = gammap[i] / (rhop[i] * rhop[i]);
        double Wg = Wg_pref * Wg_tmp * Wg_tmp;

        w0p[i] = std::sqrt(Wp + Wg);
        kappap[i] = kappa_pref * std::pow(rhop[i], kappa_pow);
        nact++;
    }

    if (internal) {
        return vv_values_;
    }

    // printf("Nact/Ntot %zu / %d\n", nact, npoints);
    // Sieve the results
    auto w_vec = std::make_shared<Vector>("W Grid points", nact);
    auto x_vec = std::make_shared<Vector>("X Grid points", nact);
    auto y_vec = std::make_shared<Vector>("Y Grid points", nact);
    auto z_vec = std::make_shared<Vector>("Z Grid points", nact);
    auto rho_vec = std::make_shared<Vector>("RHO points", nact);
    auto w0_vec = std::make_shared<Vector>("W0 points", nact);
    auto kappa_vec = std::make_shared<Vector>("KAPPA points", nact);

    double* w_vecp = w_vec->pointer();
    double* x_vecp = x_vec->pointer();
    double* y_vecp = y_vec->pointer();
    double* z_vecp = z_vec->pointer();
    double* rho_vecp = rho_vec->pointer();
    double* w0_vecp = w0_vec->pointer();
    double* kappa_vecp = kappa_vec->pointer();

    size_t sieve_pos = 0;
    for (size_t i = 0; i < npoints; i++) {
        if (rhop[i] < rho_thresh) continue;

        w_vecp[sieve_pos] = block->w()[i];
        x_vecp[sieve_pos] = block->x()[i];
        y_vecp[sieve_pos] = block->y()[i];
        z_vecp[sieve_pos] = block->z()[i];

        rho_vecp[sieve_pos] = rhop[i];
        w0_vecp[sieve_pos] = w0p[i];
        kappa_vecp[sieve_pos] = kappap[i];

        sieve_pos++;
    }

    std::map<std::string, SharedVector> ret;

    ret["W"] = w_vec;
    ret["X"] = x_vec;
    ret["Y"] = y_vec;
    ret["Z"] = z_vec;
    ret["RHO"] = rho_vec;
    ret["W0"] = w0_vec;
    ret["KAPPA"] = kappa_vec;

    return ret;
}
double SuperFunctional::compute_vv10_kernel(const std::map<std::string, SharedVector>& vals,
                                            const std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                                            std::shared_ptr<BlockOPoints> block, int npoints, bool do_grad) {
    // Kernel between left (*this) and right (vv10_cache) grids

    // Compute the vv10 cache in place
    const double l_thresh = 1.e-12;
    const size_t l_npoints = block->npoints();
    compute_vv10_cache(vals, block, l_thresh, block->npoints(), true);

    // Grab values to update
    double vv10_e = 0.0;
    double* v_rho = values_["V_RHO_A"]->pointer();
    double* v_gamma = values_["V_GAMMA_AA"]->pointer();
    double* x_grid = vv_values_["GRID_WX"]->pointer();
    double* y_grid = vv_values_["GRID_WY"]->pointer();
    double* z_grid = vv_values_["GRID_WZ"]->pointer();
    std::fill(x_grid, x_grid + l_npoints, 0.0);
    std::fill(y_grid, y_grid + l_npoints, 0.0);
    std::fill(z_grid, z_grid + l_npoints, 0.0);

    // Constants
    const double vv10_beta = vv10_beta_;

    // Get left points
    const double* l_x = block->x();
    const double* l_y = block->y();
    const double* l_z = block->z();
    const double* l_w = block->w();
    const double* l_rho = vals.find("RHO_A")->second->pointer();
    const double* l_gamma = vals.find("GAMMA_AA")->second->pointer();
    const double* l_W0 = vv_values_["W0"]->pointer();
    const double* l_kappa = vv_values_["KAPPA"]->pointer();

    for (size_t i = 0; i < l_npoints; i++) {
        // Add Phi agnostic quantities
        vv10_e += l_w[i] * l_rho[i] * vv10_beta;
        v_rho[i] = vv10_beta;
        v_gamma[i] = 0.0;

        if (l_rho[i] < l_thresh) continue;

        // Compute interior kernel
        double phi = 0.0;
        double U = 0.0;
        double W = 0.0;
        double xc = 0.0;
        double yc = 0.0;
        double zc = 0.0;
        for (auto r_block : vv10_cache) {
            // Get right points
            const double* r_x = r_block["X"]->pointer();
            const double* r_y = r_block["Y"]->pointer();
            const double* r_z = r_block["Z"]->pointer();
            const double* r_w = r_block["W"]->pointer();
            const double* r_rho = r_block["RHO"]->pointer();
            const double* r_W0 = r_block["W0"]->pointer();
            const double* r_kappa = r_block["KAPPA"]->pointer();

            const size_t r_npoints = r_block["KAPPA"]->dimpi()[0];

            // Interior Kernel
            if (do_grad) {
#pragma omp simd reduction(+ : phi, U, W, xc, yc, zc)
                for (size_t j = 0; j < r_npoints; j++) {
                    // Distance between grid points
                    const double d_x = l_x[i] - r_x[j];
                    const double d_y = l_y[i] - r_y[j];
                    const double d_z = l_z[i] - r_z[j];
                    const double R2 = d_x * d_x + d_y * d_y + d_z * d_z;

                    // g/gp values
                    const double g = l_W0[i] * R2 + l_kappa[i];
                    const double gp = r_W0[j] * R2 + r_kappa[j];
                    const double gs = g + gp;

                    // Sum the kernel
                    const double phi_kernel = (-1.5 * r_w[j] * r_rho[j]) / (g * gp * gs);

                    // Dumb question, does FMA do subtraction?
                    phi += phi_kernel;
                    const double tmp_U = -1.0 * phi_kernel * ((1.0 / g) + (1.0 / gs));
                    U += tmp_U;
                    W += tmp_U * R2;

                    // Grid contribution
                    const double Q = -2.0 * phi_kernel * (l_W0[i] / g + r_W0[j] / gp + (l_W0[i] + r_W0[j]) / gs);
                    xc += Q * d_x;
                    yc += Q * d_y;
                    zc += Q * d_z;
                }

            } else {
#pragma omp simd reduction(+ : phi, U, W)
                for (size_t j = 0; j < r_npoints; j++) {
                    // Distance between grid points
                    const double d_x = l_x[i] - r_x[j];
                    const double d_y = l_y[i] - r_y[j];
                    const double d_z = l_z[i] - r_z[j];
                    const double R2 = d_x * d_x + d_y * d_y + d_z * d_z;

                    // g/gp values
                    const double g = l_W0[i] * R2 + l_kappa[i];
                    const double gp = r_W0[j] * R2 + r_kappa[j];
                    const double gs = g + gp;

                    // Sum the kernel
                    const double phi_kernel = (-1.5 * r_w[j] * r_rho[j]) / (g * gp * gs);

                    // Dumb question, does FMA do subtraction?
                    phi += phi_kernel;
                    const double tmp_U = -1.0 * phi_kernel * ((1.0 / g) + (1.0 / gs));
                    U += tmp_U;
                    W += tmp_U * R2;
                }
            }

        }  // End r blocks
        // Mathematica for the win
        const double kappa_dn = l_kappa[i] / (6.0 * l_rho[i]);
        const double w0_dgamma = vv10_c_ * l_gamma[i] / (l_W0[i] * std::pow(l_rho[i], 4.0));
        const double w0_drho =
            2.0 / l_W0[i] * (M_PI / 3.0 - vv10_c_ * (l_gamma[i] * l_gamma[i]) / std::pow(l_rho[i], 5.0));

        // Sum it all together
        vv10_e += 0.5 * l_w[i] * l_rho[i] * phi;
        v_rho[i] += phi + l_rho[i] * (kappa_dn * U + w0_drho * W);
        v_gamma[i] += l_rho[i] * w0_dgamma * W;
        if (do_grad) {
            x_grid[i] += l_rho[i] * l_w[i] * xc;
            y_grid[i] += l_rho[i] * l_w[i] * yc;
            z_grid[i] += l_rho[i] * l_w[i] * zc;
        }
    }

    // printf("Nact/Ntot Ext %zu / %zu\n", nact, l_npoints);
    return vv10_e;
}
void SuperFunctional::test_functional(SharedVector rho_a, SharedVector rho_b, SharedVector gamma_aa,
                                      SharedVector gamma_ab, SharedVector gamma_bb, SharedVector tau_a,
                                      SharedVector tau_b) {
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
int SuperFunctional::ansatz() const {
    if (is_meta()) return 2;
    if (is_gga()) return 1;
    return 0;
}
}  // namespace psi
