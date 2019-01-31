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

#include "functional.h"
#include "LibXCfunctional.h"

#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <cmath>
#include <string>
#include <algorithm>

// LibXC helper utility for setter functions, not really supposed to do this
#include <xc.h>

using namespace psi;

namespace psi {

LibXCFunctional::LibXCFunctional(std::string xc_name, bool unpolarized) {
    xc_func_name_ = xc_name;
    func_id_ = xc_functional_get_number(xc_name.c_str());
    unpolarized_ = unpolarized;
    lr_exch_ = 0.0;
    global_exch_ = 0.0;

    xc_functional_ = std::unique_ptr<xc_func_type>(new xc_func_type);

    // Build the functional
    int polar_value;
    if (unpolarized_) {
        polar_value = 1;
    } else {
        polar_value = 2;
    }

    // Replace in C++14
    xc_functional_ = std::unique_ptr<xc_func_type>(new xc_func_type);

    if (xc_func_init(xc_functional_.get(), func_id_, polar_value) != 0) {
        outfile->Printf("Functional '%d' not found\n", xc_name.c_str());
        throw PSIEXCEPTION("Could not find required LibXC functional");
    }

    // Extract citation information
    name_ = xc_name;
    description_ = "    " + std::string(xc_functional_->info->name);
    for (size_t i = 0; i < 5; i++) {
        if (xc_functional_->info->refs[i]) {
            if (i != 0) {
                citation_ += "\n";
            }
            citation_ += "    ";
            citation_ += xc_functional_->info->refs[i]->ref;
        }
    }

    // Extract variables
    if (xc_functional_->info->family == XC_FAMILY_HYB_GGA || xc_functional_->info->family == XC_FAMILY_HYB_MGGA) {
        /* Range separation? */
        lrc_ = false;
        if (xc_functional_->info->flags & XC_FLAGS_HYB_CAMY) {
            outfile->Printf("Functional '%d' is a HYB_CAMY functional which is not supported in Psi4\n",
                            xc_name.c_str());
            throw PSIEXCEPTION("HYB_CAMY functionals not supported in Psi4 at present");
        }
        if (xc_functional_->info->flags & XC_FLAGS_HYB_LC) {
            outfile->Printf("Functional '%d' is a HYB_LC functional which is not supported in Psi4\n", xc_name.c_str());
            throw PSIEXCEPTION("HYB_LC functionals not supported in Psi4 at present");
        }
        if (xc_functional_->info->flags & XC_FLAGS_HYB_LCY) {
            outfile->Printf("Functional '%d' is a HYB_LCY functional which is not supported in Psi4\n",
                            xc_name.c_str());
            throw PSIEXCEPTION("HYB_LCY functionals not supported in Psi4 at present");
        }
        if (xc_functional_->info->flags & XC_FLAGS_HYB_CAM) {
            lrc_ = true;
            double alpha, beta;
            xc_hyb_cam_coef(xc_functional_.get(), &omega_, &alpha, &beta);

            /*
              The values alpha and beta have a different meaning in
              psi4 and libxc.

              In libxc, alpha is the contribution from full exact
              exchange (at all ranges), and beta is the contribution
              from short-range only exchange, yielding alpha exact
              exchange at the long range and alpha+beta in the short
              range in total.

              In Psi4, alpha is the amount of exchange at all ranges,
              while beta is the difference between the amount of
              exchange in the long range and in the short range,
              meaning alpha+beta at the long range, and alpha only at
              the short range.

              These differences amount to the transform

              SR      = LibXC_ALPHA + LibXC_BETA = Psi4_ALPHA
              LR      = LibXC_ALPHA              = Psi4_ALPHA + Psi4_BETA
              LR - SR =             - LibXC_BETA =              Psi4_BETA
            */

            global_exch_ = alpha + beta;
            lr_exch_ = -1.0 * beta;
        }

        if (!lrc_) {
            // Global hybrid
            global_exch_ = xc_hyb_exx_coef(xc_functional_.get());
        }
    }

    // Figure out the family
    int family = xc_functional_->info->family;

    std::vector<int> gga_vec = {XC_FAMILY_GGA, XC_FAMILY_HYB_GGA};
    if (std::find(gga_vec.begin(), gga_vec.end(), family) != gga_vec.end()) {
        gga_ = true;
    }

    std::vector<int> meta_vec = {XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA};
    if (std::find(meta_vec.begin(), meta_vec.end(), family) != meta_vec.end()) {
        gga_ = true;
        meta_ = true;
    }

    // Set any other parameters
    user_omega_ = false;
    exc_ = xc_functional_->info->flags & XC_FLAGS_HAVE_EXC;
    vxc_ = xc_functional_->info->flags & XC_FLAGS_HAVE_VXC;
    fxc_ = xc_functional_->info->flags & XC_FLAGS_HAVE_FXC;

    // VV10 info, ONLY for passing up the chain
    needs_vv10_ = false;
    vv10_c_ = 0.0;
    vv10_b_ = 0.0;
    if (xc_functional_->info->flags & XC_FLAGS_VV10) {
        xc_nlc_coef(xc_functional_.get(), &vv10_b_, &vv10_c_);
        needs_vv10_ = true;
    }
}
LibXCFunctional::~LibXCFunctional() { xc_func_end(xc_functional_.get()); }
std::shared_ptr<Functional> LibXCFunctional::build_worker() {
    // Build functional
    auto func = std::make_shared<LibXCFunctional>(xc_func_name_, unpolarized_);

    // User omega
    if (user_omega_) {
        func->set_omega(omega_);
    }

    // User tweakers
    if (user_tweakers_.size()) {
        func->set_tweak(user_tweakers_);
    }

    func->set_alpha(alpha_);
    func->set_gga(gga_);
    func->set_meta(meta_);
    func->set_lsda_cutoff(lsda_cutoff_);
    func->set_meta_cutoff(meta_cutoff_);
    func->exc_ = exc_;
    func->vxc_ = vxc_;
    func->fxc_ = fxc_;

    return static_cast<std::shared_ptr<Functional>>(func);
}
void LibXCFunctional::set_omega(double omega) {
    omega_ = omega;
    user_omega_ = true;
    if (xc_func_name_ == "XC_GGA_X_WPBEH") {
        xc_gga_x_wpbeh_set_params(xc_functional_.get(), omega);
    } else if (xc_func_name_ == "XC_GGA_X_HJS_PBE") {
        xc_gga_x_hjs_set_params(xc_functional_.get(), omega);
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_LRC_WPBEH") {
        xc_gga_x_wpbeh_set_params(xc_functional_->func_aux[0], omega);
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97X") {
        xc_functional_->cam_omega = omega;
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97") {
        xc_functional_->cam_omega = omega;
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97X_V") {
        xc_functional_->cam_omega = omega;
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97X_D") {
        xc_functional_->cam_omega = omega;
    } else if (xc_func_name_ == "XC_HYB_MGGA_X_M11") {
        xc_functional_->cam_omega = omega;
    } else {
        outfile->Printf("LibXCfunctional: set_omega is not defined for functional %s\n.", xc_func_name_.c_str());
        throw PSIEXCEPTION("LibXCfunctional: set_omega not defined for input functional");
    }
}
std::map<std::string, double> LibXCFunctional::query_libxc(const std::string& functional) {
    std::map<std::string, double> params;

    if (functional == "XC_HYB_CAM_COEF") {
        double omega, alpha, beta;
        xc_hyb_cam_coef(xc_functional_.get(), &omega, &alpha, &beta);
        // LibXC and Psi4 conventions for alpha and beta differ, see
        // above for full description
        params["OMEGA"] = omega;
        params["ALPHA"] = alpha + beta;
        params["BETA"] = -beta;
    } else if (functional == "XC_NLC_COEF") {
        double nlc_b, nlc_c;
        xc_nlc_coef(xc_functional_.get(), &nlc_b, &nlc_c);
        params["NLC_B"] = nlc_b;
        params["NLC_C"] = nlc_c;
    } else if (functional == "XC_HYB_EXX_COEF") {
        double mixing = xc_hyb_exx_coef(xc_functional_.get());
        params["MIXING"] = mixing;
    } else {
        outfile->Printf("LibXCFunctional: query_libxc unknown function to query parameters for: %s\n.",
                        functional.c_str());
        throw PSIEXCEPTION("LibXCFunctional: query_libxc unknown functional.");
    }

    return params;
}
void LibXCFunctional::set_tweak(std::vector<double> values) {
    bool failed = true;
    size_t vsize = values.size();
    if (xc_func_name_ == "XC_GGA_X_B86") {
        if (vsize == 3) {
            // (XC(func_type) *p, FLOAT beta, FLOAT gamma, FLOAT omega);
            xc_gga_x_b86_set_params(xc_functional_.get(), values[0], values[1], values[2]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_X_B88") {
        if (vsize == 2) {
            // (XC(func_type) *p, FLOAT beta, FLOAT gamma);
            xc_gga_x_b88_set_params(xc_functional_.get(), values[0], values[1]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_X_PBE") {
        if (vsize == 2) {
            //  (XC(func_type) *p, FLOAT kappa, FLOAT mu);
            xc_gga_x_pbe_set_params(xc_functional_.get(), values[0], values[1]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_C_PBE") {
        if (vsize == 1) {
            // (XC(func_type) *p, FLOAT beta);
            xc_gga_c_pbe_set_params(xc_functional_.get(), values[0]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_X_PW91") {
        if (vsize == 7) {
            // (XC(func_type) *p, FLOAT a, FLOAT b, FLOAT c, FLOAT d, FLOAT f, FLOAT alpha, FLOAT expo);
            xc_gga_x_pw91_set_params(xc_functional_.get(), values[0], values[1], values[2], values[3], values[4],
                                     values[5], values[6]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_X_RPBE") {
        if (vsize == 2) {
            // (XC(func_type) *p, FLOAT kappa, FLOAT mu);
            xc_gga_x_rpbe_set_params(xc_functional_.get(), values[0], values[1]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_X_OPTX") {
        if (vsize == 3) {
            // (XC(func_type) *p, FLOAT a, FLOAT b, FLOAT gamma);
            xc_gga_x_optx_set_params(xc_functional_.get(), values[0], values[1], values[2]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_GGA_C_LYP") {
        if (vsize == 4) {
            // (XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d);
            xc_gga_c_lyp_set_params(xc_functional_.get(), values[0], values[1], values[2], values[3]);
            failed = false;
        }
    } else if ((xc_func_name_ == "XC_HYB_GGA_XC_HSE03") || (xc_func_name_ == "XC_HYB_GGA_XC_HSE06")) {
        if (vsize == 3) {
            // "Mixing parameter beta", "Screening parameter omega_HF", "Screening parameter omega_PBE"
            xc_func_set_ext_params(xc_functional_.get(), values.data());
            failed = false;
        }
    } else if (xc_func_name_ == "XC_MGGA_X_TPSS") {
        if (vsize == 5) {
            // (xc_func_type *p, double b, double c, double e, double kappa, double mu, double BLOC_a, double BLOC_bu);
            xc_mgga_x_tpss_set_params(xc_functional_.get(), values[0], values[1], values[2], values[3], values[4], 2.0,
                                      0.0);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_MGGA_C_TPSS") {
        if (vsize == 6) {
            // (xc_func_type *p, double beta, double d, double C0_0, double C0_1, double C0_2, double C0_3);
            xc_mgga_c_tpss_set_params(xc_functional_.get(), values[0], values[1], values[2], values[3], values[4],
                                      values[5]);
            failed = false;
        }
    } else if (xc_func_name_ == "XC_MGGA_C_BC95") {
        if (vsize == 2) {
            // (XC(func_type) *p, FLOAT css, FLOAT copp);
            xc_mgga_c_bc95_set_params(xc_functional_.get(), values[0], values[1]);
            failed = false;
        }
        // } else if (xc_func_name_ == "XC_MGGA_C_PKZB") {
        //     if (vsize == 6) {
        //         // ((XC(func_type) *p, FLOAT beta, FLOAT d, FLOAT C0_0, FLOAT C0_1, FLOAT C0_2, FLOAT
        //         // C0_3);
        //         xc_mgga_c_pkzb_set_params(xc_functional_.get(), values[0], values[1], values[2], values[3],
        //                                   values[4], values[5]);
        //         failed = false;
        //     }
    } else {
        throw PSIEXCEPTION(
            "LibXCfunctional: set_tweak: There are no known tweaks for this functional, please double check "
            "the functional form and add them if required.");
    }

    // Did we match fully?
    if (failed) {
        throw PSIEXCEPTION(
            "LibXCfunctional: set_tweak: Mismatch in size of tweaker vector and expected number of "
            "input parameters.");
    }

    user_tweakers_ = values;
}
std::vector<std::tuple<std::string, int, double>> LibXCFunctional::get_mix_data() {
    std::vector<std::tuple<std::string, int, double>> ret;

    if (xc_functional_->mix_coef == nullptr) {
        std::string name = std::string(xc_functional_->info->name);
        int kind = xc_functional_->info->kind;
        double coef = 1.0;
        ret.push_back(std::make_tuple(name, kind, coef));

    } else {
        for (size_t i = 0; i < xc_functional_->n_func_aux; i++) {
            std::string name = std::string(xc_functional_->func_aux[i]->info->name);
            int kind = xc_functional_->func_aux[i]->info->kind;
            double coef = xc_functional_->mix_coef[i];

            ret.push_back(std::make_tuple(name, kind, coef));
        }
    }
    return ret;
}
void LibXCFunctional::compute_functional(const std::map<std::string, SharedVector>& in,
                                         const std::map<std::string, SharedVector>& out, int npoints, int deriv) {
    // => Input variables <= //

    if ((deriv >= 1) & (!vxc_)) {
        std::string error = "LibXCfunctional: No derivative implemented for " + name_ + ".";
        throw PSIEXCEPTION(error.c_str());
    }
    if ((deriv >= 2) & (!fxc_)) {
        std::string error = "LibXCfunctional: No second derivative implemented for " + name_ + ".";
        throw PSIEXCEPTION(error.c_str());
    }
    if (deriv >= 3) {
        throw PSIEXCEPTION("LibXCfunctional: Third derivatives are not implemented!");
    }

    // => Input variables <= //

    double* rho_ap = nullptr;
    double* rho_bp = nullptr;
    double* gamma_aap = nullptr;
    double* gamma_abp = nullptr;
    double* gamma_bbp = nullptr;
    double* tau_ap = nullptr;
    double* tau_bp = nullptr;
    double* lapl_ap = nullptr;
    double* lapl_bp = nullptr;

    if (true) {
        rho_ap = in.find("RHO_A")->second->pointer();

        if (!unpolarized_) {
            rho_bp = in.find("RHO_B")->second->pointer();
        }
    }
    if (gga_) {
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        if (!unpolarized_) {
            gamma_abp = in.find("GAMMA_AB")->second->pointer();
            gamma_bbp = in.find("GAMMA_BB")->second->pointer();
        }
    }
    if (meta_) {
        tau_ap = in.find("TAU_A")->second->pointer();
        // lapl_ap = in.find("LAPL_RHO_A")->second->pointer();
        if (!unpolarized_) {
            tau_bp = in.find("TAU_B")->second->pointer();
            // lapl_bp = in.find("LAPL_RHO_B")->second->pointer();
        }
    }

    // => Outut variables <= //

    double* v = nullptr;

    double* v_rho_a = nullptr;
    double* v_rho_b = nullptr;
    double* v_gamma_aa = nullptr;
    double* v_gamma_ab = nullptr;
    double* v_gamma_bb = nullptr;
    double* v_tau_a = nullptr;
    double* v_tau_b = nullptr;
    double* v_lapl_a = nullptr;
    double* v_lapl_b = nullptr;

    double* v_rho_a_rho_a = nullptr;
    double* v_rho_a_rho_b = nullptr;
    double* v_rho_b_rho_b = nullptr;
    double* v_gamma_aa_gamma_aa = nullptr;
    double* v_gamma_aa_gamma_ab = nullptr;
    double* v_gamma_aa_gamma_bb = nullptr;
    double* v_gamma_ab_gamma_ab = nullptr;
    double* v_gamma_ab_gamma_bb = nullptr;
    double* v_gamma_bb_gamma_bb = nullptr;
    double* v_tau_a_tau_a = nullptr;
    double* v_tau_a_tau_b = nullptr;
    double* v_tau_b_tau_b = nullptr;
    double* v_rho_a_gamma_aa = nullptr;
    double* v_rho_a_gamma_ab = nullptr;
    double* v_rho_a_gamma_bb = nullptr;
    double* v_rho_b_gamma_aa = nullptr;
    double* v_rho_b_gamma_ab = nullptr;
    double* v_rho_b_gamma_bb = nullptr;
    double* v_rho_a_tau_a = nullptr;
    double* v_rho_a_tau_b = nullptr;
    double* v_rho_b_tau_a = nullptr;
    double* v_rho_b_tau_b = nullptr;
    double* v_gamma_aa_tau_a = nullptr;
    double* v_gamma_aa_tau_b = nullptr;
    double* v_gamma_ab_tau_a = nullptr;
    double* v_gamma_ab_tau_b = nullptr;
    double* v_gamma_bb_tau_a = nullptr;
    double* v_gamma_bb_tau_b = nullptr;

    if (deriv >= 0) {
        // Energy doesnt make sense for all functionals
        if (exc_) {
            v = out.find("V")->second->pointer();
        }
    }
    if (deriv >= 1) {
        if (true) {
            v_rho_a = out.find("V_RHO_A")->second->pointer();
            if (!unpolarized_) {
                v_rho_b = out.find("V_RHO_B")->second->pointer();
            }
        }
        if (gga_) {
            v_gamma_aa = out.find("V_GAMMA_AA")->second->pointer();
            if (!unpolarized_) {
                v_gamma_ab = out.find("V_GAMMA_AB")->second->pointer();
                v_gamma_bb = out.find("V_GAMMA_BB")->second->pointer();
            }
        }
        if (meta_) {
            v_tau_a = out.find("V_TAU_A")->second->pointer();
            if (!unpolarized_) {
                v_tau_b = out.find("V_TAU_B")->second->pointer();
            }
            // v_lapl_a = out.find("V_LAPL_A")->second->pointer();
            // v_lapl_b = out.find("V_LAPL_B")->second->pointer();
        }
    }
    if (deriv >= 2) {
        if (true) {
            v_rho_a_rho_a = out.find("V_RHO_A_RHO_A")->second->pointer();

            if (!unpolarized_) {
                v_rho_a_rho_b = out.find("V_RHO_A_RHO_B")->second->pointer();
                v_rho_b_rho_b = out.find("V_RHO_B_RHO_B")->second->pointer();
            }
        }
        if (gga_) {
            v_gamma_aa_gamma_aa = out.find("V_GAMMA_AA_GAMMA_AA")->second->pointer();

            if (!unpolarized_) {
                v_gamma_aa_gamma_ab = out.find("V_GAMMA_AA_GAMMA_AB")->second->pointer();
                v_gamma_aa_gamma_bb = out.find("V_GAMMA_AA_GAMMA_BB")->second->pointer();
                v_gamma_ab_gamma_ab = out.find("V_GAMMA_AB_GAMMA_AB")->second->pointer();
                v_gamma_ab_gamma_bb = out.find("V_GAMMA_AB_GAMMA_BB")->second->pointer();
                v_gamma_bb_gamma_bb = out.find("V_GAMMA_BB_GAMMA_BB")->second->pointer();
            }
        }
        if (meta_) {
            v_tau_a_tau_a = out.find("V_TAU_A_TAU_A")->second->pointer();

            if (!unpolarized_) {
                v_tau_a_tau_b = out.find("V_TAU_A_TAU_B")->second->pointer();
                v_tau_b_tau_b = out.find("V_TAU_B_TAU_B")->second->pointer();
            }
        }
        if (gga_) {
            v_rho_a_gamma_aa = out.find("V_RHO_A_GAMMA_AA")->second->pointer();

            if (!unpolarized_) {
                v_rho_a_gamma_ab = out.find("V_RHO_A_GAMMA_AB")->second->pointer();
                v_rho_a_gamma_bb = out.find("V_RHO_A_GAMMA_BB")->second->pointer();
                v_rho_b_gamma_aa = out.find("V_RHO_B_GAMMA_AA")->second->pointer();
                v_rho_b_gamma_ab = out.find("V_RHO_B_GAMMA_AB")->second->pointer();
                v_rho_b_gamma_bb = out.find("V_RHO_B_GAMMA_BB")->second->pointer();
            }
        }
        if (meta_) {
            v_rho_a_tau_a = out.find("V_RHO_A_TAU_A")->second->pointer();

            if (!unpolarized_) {
                v_rho_a_tau_b = out.find("V_RHO_A_TAU_B")->second->pointer();
                v_rho_b_tau_a = out.find("V_RHO_B_TAU_A")->second->pointer();
                v_rho_b_tau_b = out.find("V_RHO_B_TAU_B")->second->pointer();
            }
        }
        if (gga_ && meta_) {
            v_gamma_aa_tau_a = out.find("V_GAMMA_AA_TAU_A")->second->pointer();
            if (!unpolarized_) {
                v_gamma_aa_tau_b = out.find("V_GAMMA_AA_TAU_B")->second->pointer();
                v_gamma_ab_tau_a = out.find("V_GAMMA_AB_TAU_A")->second->pointer();
                v_gamma_ab_tau_b = out.find("V_GAMMA_AB_TAU_B")->second->pointer();
                v_gamma_bb_tau_a = out.find("V_GAMMA_BB_TAU_A")->second->pointer();
                v_gamma_bb_tau_b = out.find("V_GAMMA_BB_TAU_B")->second->pointer();
            }
        }
    }

    if (unpolarized_) {
        if (deriv == 0) {
            throw PSIEXCEPTION("LibXCFunction deriv=0 is not implemented, call deriv >=1");
        }

        // Compute deriv
        if (deriv >= 1) {
            // Allocate
            std::vector<double> fv(npoints);
            std::vector<double> fv_rho(npoints);

            // GGA
            std::vector<double> fgamma, fv_gamma;
            if (gga_) {
                fgamma.resize(npoints);
                fv_gamma.resize(npoints);
            }

            // Meta
            std::vector<double> flapl, fv_lapl, fv_tau;
            if (meta_) {
                flapl.resize(npoints);
                fv_lapl.resize(npoints);
                fv_tau.resize(npoints);
            }

            double* fvp = nullptr;
            if (exc_) {
                fvp = fv.data();
            }

            // Compute
            if (meta_) {
                xc_mgga_exc_vxc(xc_functional_.get(), npoints, rho_ap, gamma_aap, flapl.data(), tau_ap, fvp,
                                fv_rho.data(), fv_gamma.data(), fv_lapl.data(), fv_tau.data());
            } else if (gga_) {
                xc_gga_exc_vxc(xc_functional_.get(), npoints, rho_ap, gamma_aap, fvp, fv_rho.data(), fv_gamma.data());

            } else {
                xc_lda_exc_vxc(xc_functional_.get(), npoints, rho_ap, fvp, fv_rho.data());
            }
            // printf("%s | %lf %lf\n", xc_func_name_.c_str(), fv_rho[0], fv_gamma[0]);

            // Re-apply
            if (exc_) {
                for (size_t i = 0; i < npoints; i++) {
                    v[i] += alpha_ * fv[i] * rho_ap[i];
                }
            }

            C_DAXPY(npoints, alpha_, fv_rho.data(), 1, v_rho_a, 1);

            if (gga_) {
                C_DAXPY(npoints, alpha_, fv_gamma.data(), 1, v_gamma_aa, 1);
            }

            if (meta_) {
                C_DAXPY(npoints, 0.5 * alpha_, fv_tau.data(), 1, v_tau_a, 1);
            }
        }

        // Compute second derivative
        if (deriv >= 2) {
            if (meta_) {
                throw PSIEXCEPTION(
                    "Second derivative for meta functionals is not yet "
                    "available");

            } else if (gga_) {
                std::vector<double> fv2_rho2(npoints);
                std::vector<double> fv2_rho_gamma(npoints);
                std::vector<double> fv2_gamma2(npoints);

                xc_gga_fxc(xc_functional_.get(), npoints, rho_ap, gamma_aap, fv2_rho2.data(), fv2_rho_gamma.data(),
                           fv2_gamma2.data());

                C_DAXPY(npoints, alpha_, fv2_rho2.data(), 1, v_rho_a_rho_a, 1);
                C_DAXPY(npoints, alpha_, fv2_gamma2.data(), 1, v_gamma_aa_gamma_aa, 1);
                C_DAXPY(npoints, alpha_, fv2_rho_gamma.data(), 1, v_rho_a_gamma_aa, 1);

            } else {
                std::vector<double> fv2_rho2(npoints);

                xc_lda_fxc(xc_functional_.get(), npoints, rho_ap, fv2_rho2.data());

                C_DAXPY(npoints, alpha_, fv2_rho2.data(), 1, v_rho_a_rho_a, 1);
            }
        }

    } else {  // End unpolarized

        // Allocate input data
        std::vector<double> frho(npoints * 2);
        std::vector<double> fv(npoints);
        std::vector<double> fv_rho(npoints * 2);

        C_DCOPY(npoints, rho_ap, 1, frho.data(), 2);
        C_DCOPY(npoints, rho_bp, 1, (frho.data() + 1), 2);

        std::vector<double> fgamma, fv_gamma;
        if (gga_) {
            fgamma.resize(npoints * 3);
            fv_gamma.resize(npoints * 3);

            C_DCOPY(npoints, gamma_aap, 1, fgamma.data(), 3);
            C_DCOPY(npoints, gamma_abp, 1, (fgamma.data() + 1), 3);
            C_DCOPY(npoints, gamma_bbp, 1, (fgamma.data() + 2), 3);
        }

        std::vector<double> ftau, flapl, fv_lapl, fv_tau;
        if (meta_) {
            ftau.resize(npoints * 2);
            flapl.resize(npoints * 2);
            fv_lapl.resize(npoints * 2);
            fv_tau.resize(npoints * 2);

            C_DCOPY(npoints, tau_ap, 1, ftau.data(), 2);
            C_DCOPY(npoints, tau_bp, 1, (ftau.data() + 1), 2);
        }

        // Compute first deriv
        if (deriv == 0) {
            throw PSIEXCEPTION("LibXCFunction deriv=0 is not implemented, call deriv >=1");
        }
        if (deriv >= 1) {
            // Special cases
            double* fvp;
            if (exc_) {
                fvp = fv.data();
            } else {
                fvp = nullptr;
            }

            if (meta_) {
                xc_mgga_exc_vxc(xc_functional_.get(), npoints, frho.data(), fgamma.data(), flapl.data(), ftau.data(),
                                fvp, fv_rho.data(), fv_gamma.data(), fv_lapl.data(), fv_tau.data());

            } else if (gga_) {
                xc_gga_exc_vxc(xc_functional_.get(), npoints, frho.data(), fgamma.data(), fvp, fv_rho.data(),
                               fv_gamma.data());

            } else {
                xc_lda_exc_vxc(xc_functional_.get(), npoints, frho.data(), fvp, fv_rho.data());
            }

            // Re-apply
            if (exc_) {
                for (size_t i = 0; i < npoints; i++) {
                    v[i] += alpha_ * fv[i] * (rho_ap[i] + rho_bp[i]);
                }
            }

            C_DAXPY(npoints, alpha_, fv_rho.data(), 2, v_rho_a, 1);
            C_DAXPY(npoints, alpha_, (fv_rho.data() + 1), 2, v_rho_b, 1);

            if (gga_) {
                C_DAXPY(npoints, alpha_, fv_gamma.data(), 3, v_gamma_aa, 1);
                C_DAXPY(npoints, alpha_, (fv_gamma.data() + 1), 3, v_gamma_ab, 1);
                C_DAXPY(npoints, alpha_, (fv_gamma.data() + 2), 3, v_gamma_bb, 1);
            }

            if (meta_) {
                C_DAXPY(npoints, 0.5 * alpha_, fv_tau.data(), 2, v_tau_a, 1);
                C_DAXPY(npoints, 0.5 * alpha_, (fv_tau.data() + 1), 2, v_tau_b, 1);
            }
        }

        // Compute second deriv
        if (deriv >= 2) {
            if (meta_) {
                throw PSIEXCEPTION("Second derivative for meta functionals is not yet available");

            } else if (gga_) {
                std::vector<double> fv2_rho2(npoints * 3);
                std::vector<double> fv2_rhogamma(npoints * 6);
                std::vector<double> fv2_gamma2(npoints * 6);

                xc_gga_fxc(xc_functional_.get(), npoints, frho.data(), fgamma.data(), fv2_rho2.data(),
                           fv2_rhogamma.data(), fv2_gamma2.data());

                for (size_t i = 0; i < npoints; i++) {
                    // v2rho2(3)       = (u_u, u_d, d_d)
                    v_rho_a_rho_a[i] += alpha_ * fv2_rho2[3 * i];
                    v_rho_a_rho_b[i] += alpha_ * fv2_rho2[3 * i + 1];
                    v_rho_b_rho_b[i] += alpha_ * fv2_rho2[3 * i + 2];

                    // v2gamma2(6)     = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)
                    v_gamma_aa_gamma_aa[i] += alpha_ * fv2_gamma2[6 * i];
                    v_gamma_aa_gamma_ab[i] += alpha_ * fv2_gamma2[6 * i + 1];
                    v_gamma_aa_gamma_bb[i] += alpha_ * fv2_gamma2[6 * i + 2];
                    v_gamma_ab_gamma_ab[i] += alpha_ * fv2_gamma2[6 * i + 3];
                    v_gamma_ab_gamma_bb[i] += alpha_ * fv2_gamma2[6 * i + 4];
                    v_gamma_bb_gamma_bb[i] += alpha_ * fv2_gamma2[6 * i + 5];

                    // v2rhogamma(6)   = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
                    v_rho_a_gamma_aa[i] += alpha_ * fv2_rhogamma[6 * i];
                    v_rho_a_gamma_ab[i] += alpha_ * fv2_rhogamma[6 * i + 1];
                    v_rho_a_gamma_bb[i] += alpha_ * fv2_rhogamma[6 * i + 2];
                    v_rho_b_gamma_aa[i] += alpha_ * fv2_rhogamma[6 * i + 3];
                    v_rho_b_gamma_ab[i] += alpha_ * fv2_rhogamma[6 * i + 4];
                    v_rho_b_gamma_bb[i] += alpha_ * fv2_rhogamma[6 * i + 5];
                }

            } else {
                std::vector<double> fv2_rho2(npoints * 3);

                xc_lda_fxc(xc_functional_.get(), npoints, frho.data(), fv2_rho2.data());

                for (size_t i = 0; i < npoints; i++) {
                    // v2rho2(3)       = (u_u, u_d, d_d)
                    v_rho_a_rho_a[i] += alpha_ * fv2_rho2[3 * i];
                    v_rho_a_rho_b[i] += alpha_ * fv2_rho2[3 * i + 1];
                    v_rho_b_rho_b[i] += alpha_ * fv2_rho2[3 * i + 2];
                }
            }
        }
        if (deriv > 2) {
            throw PSIEXCEPTION("TRYING TO COPMUTE DERIV > 3 ");
        }
    }  // End polarized
}

}  // End namespace
