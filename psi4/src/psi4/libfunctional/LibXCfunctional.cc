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

#include <cmath>
#include <string>
#include "LibXCfunctional.h"
#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"

using namespace psi;

namespace psi {

LibXCFunctional::LibXCFunctional(std::string xc_name, bool unpolarized) {
    xc_func_name_ = xc_name;
    func_id_ = xc_functional_get_number(xc_name.c_str());
    unpolarized_ = unpolarized;
    lr_exch_ = 0.0;
    global_exch_ = 0.0;

    // Build the functional
    int polar_value;
    if (unpolarized_) {
        polar_value = 1;
    } else {
        polar_value = 2;
    }
    if (xc_func_init(&xc_functional_, func_id_, polar_value) != 0) {
        outfile->Printf("Functional '%d' not found\n", xc_name.c_str());
        throw PSIEXCEPTION("Could not find required LIBXC functional");
    }

    // Extract citation information
    name_ = xc_name;
    description_ = "    " + std::string(xc_functional_.info->name);
    for (size_t i = 0; i < 5; i++) {
        if (xc_functional_.info->refs[i]) {
            if (i != 0) {
                citation_ += "\n";
            }
            citation_ += "    ";
            citation_ += xc_functional_.info->refs[i]->ref;
        }
    }

    // Extract variables
    if (xc_functional_.info->family == XC_FAMILY_HYB_GGA ||
        xc_functional_.info->family == XC_FAMILY_HYB_MGGA) {
        /* Range separation? */
        int rangesep = 0;
        if (xc_functional_.info->flags & XC_FLAGS_HYB_CAM) rangesep++;
        if (xc_functional_.info->flags & XC_FLAGS_HYB_CAMY) rangesep++;
        if (xc_functional_.info->flags & XC_FLAGS_HYB_LC) rangesep++;
        if (xc_functional_.info->flags & XC_FLAGS_HYB_LCY) rangesep++;
        if (rangesep) {
            lrc_ = true;

            double alpha, beta;
            xc_hyb_cam_coef(&xc_functional_, &omega_, &alpha, &beta);

            lr_exch_ = alpha;
            global_exch_ = alpha + beta;
            if (std::abs(1.0 - lr_exch_) > 1.e14) {
                throw PSIEXCEPTION(
                    "PSI Currently cannot computation functionals with less than 100%% long range "
                    "exact exchange.\n");
            }

        } else {
            global_exch_ = xc_hyb_exx_coef(&xc_functional_);
        }
    }

    // Figure out the family
    int family = xc_functional_.info->family;

    std::vector<int> gga_vec = {XC_FAMILY_GGA, XC_FAMILY_HYB_GGA};
    if (std::find(gga_vec.begin(), gga_vec.end(), family) != gga_vec.end()) {
        gga_ = true;
    }

    std::vector<int> meta_vec = {XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA};
    if (std::find(meta_vec.begin(), meta_vec.end(), family) != meta_vec.end()) {
        gga_ = true;
        meta_ = true;
    }
}
LibXCFunctional::~LibXCFunctional() { xc_func_end(&xc_functional_); }
std::shared_ptr<Functional> LibXCFunctional::build_worker() {
    // Build functional
    std::shared_ptr<LibXCFunctional> func(new LibXCFunctional(xc_func_name_, unpolarized_));

    // Tweak
    if (omega_ != 0.0) {
        func->set_omega(omega_);
    }

    func->set_alpha(alpha_);
    func->set_gga(gga_);
    func->set_meta(meta_);
    func->set_lsda_cutoff(lsda_cutoff_);
    func->set_meta_cutoff(meta_cutoff_);

    return static_cast<std::shared_ptr<Functional>>(func);
}
void LibXCFunctional::set_omega(double omega) {
    omega_ = omega;
    if (xc_func_name_ == "XC_GGA_X_WPBEH") {
        xc_gga_x_wpbeh_set_params(&xc_functional_, omega);
    } else if (xc_func_name_ == "XC_GGA_X_HJS_PBE") {
        xc_gga_x_hjs_set_params(&xc_functional_, omega);
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97X") {
        xc_functional_.cam_omega = omega;
    } else if (xc_func_name_ == "XC_HYB_GGA_XC_WB97") {
        xc_functional_.cam_omega = omega;
    } else {
        outfile->Printf("set_omega is not defined for functional %s\n.", xc_func_name_.c_str());
        throw PSIEXCEPTION("set_omega not defined for input functional");
    }
}
std::vector<std::tuple<std::string, int, double>> LibXCFunctional::get_mix_data() {
    std::vector<std::tuple<std::string, int, double>> ret;

    if (xc_functional_.mix_coef == nullptr) {
        std::string name = std::string(xc_functional_.info->name);
        int kind = xc_functional_.info->kind;
        double coef = 1.0;
        ret.push_back(std::make_tuple(name, kind, coef));

    } else {
        for (size_t i = 0; i < xc_functional_.n_func_aux; i++) {
            std::string name = std::string(xc_functional_.func_aux[i]->info->name);
            int kind = xc_functional_.func_aux[i]->info->kind;
            double coef = xc_functional_.mix_coef[i];

            ret.push_back(std::make_tuple(name, kind, coef));
        }
    }
    return ret;
}
void LibXCFunctional::compute_functional(const std::map<std::string, SharedVector>& in,
                                         const std::map<std::string, SharedVector>& out,
                                         int npoints, int deriv) {
    // => Input variables <= //

    double* rho_ap = NULL;
    double* rho_bp = NULL;
    double* gamma_aap = NULL;
    double* gamma_abp = NULL;
    double* gamma_bbp = NULL;
    double* tau_ap = NULL;
    double* tau_bp = NULL;
    double* lapl_ap = NULL;
    double* lapl_bp = NULL;

    // gga_ = true;

    if (true) {
        rho_ap = in.find("RHO_A")->second->pointer();
        rho_bp = in.find("RHO_B")->second->pointer();
    }
    if (gga_) {
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        gamma_abp = in.find("GAMMA_AB")->second->pointer();
        gamma_bbp = in.find("GAMMA_BB")->second->pointer();
    }
    if (meta_) {
        tau_ap = in.find("TAU_A")->second->pointer();
        tau_bp = in.find("TAU_B")->second->pointer();
        // lapl_ap = in.find("LAPL_RHO_A")->second->pointer();
        // lapl_bp = in.find("LAPL_RHO_B")->second->pointer();
    }

    // => Outut variables <= //

    double* v = NULL;

    double* v_rho_a = NULL;
    double* v_rho_b = NULL;
    double* v_gamma_aa = NULL;
    double* v_gamma_ab = NULL;
    double* v_gamma_bb = NULL;
    double* v_tau_a = NULL;
    double* v_tau_b = NULL;
    double* v_lapl_a = NULL;
    double* v_lapl_b = NULL;

    double* v_rho_a_rho_a = NULL;
    double* v_rho_a_rho_b = NULL;
    double* v_rho_b_rho_b = NULL;
    double* v_gamma_aa_gamma_aa = NULL;
    double* v_gamma_aa_gamma_ab = NULL;
    double* v_gamma_aa_gamma_bb = NULL;
    double* v_gamma_ab_gamma_ab = NULL;
    double* v_gamma_ab_gamma_bb = NULL;
    double* v_gamma_bb_gamma_bb = NULL;
    double* v_tau_a_tau_a = NULL;
    double* v_tau_a_tau_b = NULL;
    double* v_tau_b_tau_b = NULL;
    double* v_rho_a_gamma_aa = NULL;
    double* v_rho_a_gamma_ab = NULL;
    double* v_rho_a_gamma_bb = NULL;
    double* v_rho_b_gamma_aa = NULL;
    double* v_rho_b_gamma_ab = NULL;
    double* v_rho_b_gamma_bb = NULL;
    double* v_rho_a_tau_a = NULL;
    double* v_rho_a_tau_b = NULL;
    double* v_rho_b_tau_a = NULL;
    double* v_rho_b_tau_b = NULL;
    double* v_gamma_aa_tau_a = NULL;
    double* v_gamma_aa_tau_b = NULL;
    double* v_gamma_ab_tau_a = NULL;
    double* v_gamma_ab_tau_b = NULL;
    double* v_gamma_bb_tau_a = NULL;
    double* v_gamma_bb_tau_b = NULL;

    if (deriv >= 0) {
        v = out.find("V")->second->pointer();
    }
    if (deriv >= 1) {
        if (true) {
            v_rho_a = out.find("V_RHO_A")->second->pointer();
            v_rho_b = out.find("V_RHO_B")->second->pointer();
        }
        if (gga_) {
            v_gamma_aa = out.find("V_GAMMA_AA")->second->pointer();
            v_gamma_ab = out.find("V_GAMMA_AB")->second->pointer();
            v_gamma_bb = out.find("V_GAMMA_BB")->second->pointer();
        }
        if (meta_) {
            v_tau_a = out.find("V_TAU_A")->second->pointer();
            v_tau_b = out.find("V_TAU_B")->second->pointer();
            // v_lapl_a = out.find("V_LAPL_A")->second->pointer();
            // v_lapl_b = out.find("V_LAPL_B")->second->pointer();
        }
    }
    if (deriv >= 2) {
        if (true) {
            v_rho_a_rho_a = out.find("V_RHO_A_RHO_A")->second->pointer();
            v_rho_a_rho_b = out.find("V_RHO_A_RHO_B")->second->pointer();
            v_rho_b_rho_b = out.find("V_RHO_B_RHO_B")->second->pointer();
        }
        if (gga_) {
            v_gamma_aa_gamma_aa = out.find("V_GAMMA_AA_GAMMA_AA")->second->pointer();
            v_gamma_aa_gamma_ab = out.find("V_GAMMA_AA_GAMMA_AB")->second->pointer();
            v_gamma_aa_gamma_bb = out.find("V_GAMMA_AA_GAMMA_BB")->second->pointer();
            v_gamma_ab_gamma_ab = out.find("V_GAMMA_AB_GAMMA_AB")->second->pointer();
            v_gamma_ab_gamma_bb = out.find("V_GAMMA_AB_GAMMA_BB")->second->pointer();
            v_gamma_bb_gamma_bb = out.find("V_GAMMA_BB_GAMMA_BB")->second->pointer();
        }
        if (meta_) {
            v_tau_a_tau_a = out.find("V_TAU_A_TAU_A")->second->pointer();
            v_tau_a_tau_b = out.find("V_TAU_A_TAU_B")->second->pointer();
            v_tau_b_tau_b = out.find("V_TAU_B_TAU_B")->second->pointer();
        }
        if (gga_) {
            v_rho_a_gamma_aa = out.find("V_RHO_A_GAMMA_AA")->second->pointer();
            v_rho_a_gamma_ab = out.find("V_RHO_A_GAMMA_AB")->second->pointer();
            v_rho_a_gamma_bb = out.find("V_RHO_A_GAMMA_BB")->second->pointer();
            v_rho_b_gamma_aa = out.find("V_RHO_B_GAMMA_AA")->second->pointer();
            v_rho_b_gamma_ab = out.find("V_RHO_B_GAMMA_AB")->second->pointer();
            v_rho_b_gamma_bb = out.find("V_RHO_B_GAMMA_BB")->second->pointer();
        }
        if (meta_) {
            v_rho_a_tau_a = out.find("V_RHO_A_TAU_A")->second->pointer();
            v_rho_a_tau_b = out.find("V_RHO_A_TAU_B")->second->pointer();
            v_rho_b_tau_a = out.find("V_RHO_B_TAU_A")->second->pointer();
            v_rho_b_tau_b = out.find("V_RHO_B_TAU_B")->second->pointer();
        }
        if (gga_ && meta_) {
            v_gamma_aa_tau_a = out.find("V_GAMMA_AA_TAU_A")->second->pointer();
            v_gamma_aa_tau_b = out.find("V_GAMMA_AA_TAU_B")->second->pointer();
            v_gamma_ab_tau_a = out.find("V_GAMMA_AB_TAU_A")->second->pointer();
            v_gamma_ab_tau_b = out.find("V_GAMMA_AB_TAU_B")->second->pointer();
            v_gamma_bb_tau_a = out.find("V_GAMMA_BB_TAU_A")->second->pointer();
            v_gamma_bb_tau_b = out.find("V_GAMMA_BB_TAU_B")->second->pointer();
        }
    }

    if (unpolarized_) {
        if (deriv == 0) {
            throw PSIEXCEPTION("LibXCFunction deriv=0 is not implemented, call deriv >=1");
        }

        // Allocate
        std::vector<double> frho(npoints);
        std::vector<double> fv(npoints);
        std::vector<double> fv_rho(npoints);

        C_DAXPY(npoints, 2.0, rho_ap, 1, frho.data(), 1);

        // GGA
        std::vector<double> fsigma, fv_sigma;
        if (gga_) {
            fsigma.resize(npoints);
            fv_sigma.resize(npoints);
            C_DAXPY(npoints, 4.0, gamma_aap, 1, fsigma.data(), 1);
        }

        // Meta
        std::vector<double> flapl, fv_lapl, fv_tau;
        if (meta_) {
            flapl.resize(npoints);
            fv_lapl.resize(npoints);
            fv_tau.resize(npoints);
        }

        // Compute deriv
        if (deriv >= 1) {
            if (meta_) {
                xc_mgga_exc_vxc(&xc_functional_, npoints, frho.data(), fsigma.data(), flapl.data(),
                                tau_ap, fv.data(), fv_rho.data(), fv_sigma.data(), fv_lapl.data(),
                                fv_tau.data());
            } else if (gga_) {
                xc_gga_exc_vxc(&xc_functional_, npoints, frho.data(), fsigma.data(), fv.data(),
                               fv_rho.data(), fv_sigma.data());

            } else {
                xc_lda_exc_vxc(&xc_functional_, npoints, frho.data(), fv.data(), fv_rho.data());
            }

            // Re-apply
            for (size_t i = 0; i < npoints; i++) {
                v[i] += alpha_ * fv[i] * (2.0 * rho_ap[i]);
            }
            C_DAXPY(npoints, alpha_, fv_rho.data(), 1, v_rho_a, 1);

            if (gga_) {
                C_DAXPY(npoints, 2.0 * alpha_, fv_sigma.data(), 1, v_gamma_aa, 1);
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
                std::vector<double> fv2_rhosigma(npoints);
                std::vector<double> fv2_sigma2(npoints);

                xc_gga_fxc(&xc_functional_, npoints, frho.data(), fsigma.data(), fv2_rho2.data(),
                           fv2_rhosigma.data(), fv2_sigma2.data());

                C_DAXPY(npoints, alpha_, fv2_rho2.data(), 1, v_rho_a_rho_a, 1);
                C_DAXPY(npoints, alpha_, fv2_sigma2.data(), 1, v_gamma_aa_gamma_aa, 1);
                C_DAXPY(npoints, alpha_, fv2_rhosigma.data(), 1, v_rho_a_gamma_aa, 1);

            } else {
                std::vector<double> fv2_rho2(npoints);

                xc_lda_fxc(&xc_functional_, npoints, frho.data(), fv2_rho2.data());

                C_DAXPY(npoints, alpha_, fv2_rho2.data(), 1, v_rho_a_rho_a, 1);
            }
        }

    } else {  // End unpolarized
        if (deriv == 0) {
            throw PSIEXCEPTION("LibXCFunction deriv=0 is not implemented, call deriv >=1");
        }

        // Allocate input data
        std::vector<double> frho(npoints * 2);
        std::vector<double> fv(npoints);
        std::vector<double> fv_rho(npoints * 2);

        C_DAXPY(npoints, 1.0, rho_ap, 1, frho.data(), 2);
        C_DAXPY(npoints, 1.0, rho_bp, 1, (frho.data() + 1), 2);

        std::vector<double> fsigma, fv_sigma;
        if (gga_) {
            fsigma.resize(npoints * 3);
            fv_sigma.resize(npoints * 3);

            C_DAXPY(npoints, 1.0, gamma_aap, 1, fsigma.data(), 3);
            C_DAXPY(npoints, 1.0, gamma_abp, 1, (fsigma.data() + 1), 3);
            C_DAXPY(npoints, 1.0, gamma_bbp, 1, (fsigma.data() + 2), 3);
        }

        std::vector<double> ftau, flapl, fv_lapl, fv_tau;
        if (meta_) {
            ftau.resize(npoints * 2);
            flapl.resize(npoints * 2);
            fv_lapl.resize(npoints * 2);
            fv_tau.resize(npoints * 2);

            C_DAXPY(npoints, 0.5, tau_ap, 1, ftau.data(), 2);
            C_DAXPY(npoints, 0.5, tau_bp, 1, (ftau.data() + 1), 2);
        }

        // Compute first deriv
        if (deriv >= 1) {
            if (meta_) {
                xc_mgga_exc_vxc(&xc_functional_, npoints, frho.data(), fsigma.data(), flapl.data(),
                                ftau.data(), fv.data(), fv_rho.data(), fv_sigma.data(),
                                fv_lapl.data(), fv_tau.data());

            } else if (gga_) {
                xc_gga_exc_vxc(&xc_functional_, npoints, frho.data(), fsigma.data(), fv.data(),
                               fv_rho.data(), fv_sigma.data());

            } else {
                xc_lda_exc_vxc(&xc_functional_, npoints, frho.data(), fv.data(), fv_rho.data());
            }

            // Re-apply
            for (size_t i = 0; i < npoints; i++) {
                v[i] += alpha_ * fv[i] * (rho_ap[i] + rho_bp[i]);
            }

            C_DAXPY(npoints, alpha_, fv_rho.data(), 2, v_rho_a, 1);
            C_DAXPY(npoints, alpha_, (fv_rho.data() + 1), 2, v_rho_b, 1);

            if (gga_) {
                C_DAXPY(npoints, alpha_, fv_sigma.data(), 3, v_gamma_aa, 1);
                C_DAXPY(npoints, alpha_, (fv_sigma.data() + 1), 3, v_gamma_ab, 1);
                C_DAXPY(npoints, alpha_, (fv_sigma.data() + 2), 3, v_gamma_bb, 1);
            }

            if (meta_) {
                C_DAXPY(npoints, 0.5 * alpha_, fv_tau.data(), 2, v_tau_a, 1);
                C_DAXPY(npoints, 0.5 * alpha_, (fv_tau.data() + 1), 2, v_tau_b, 1);
            }
        }

        // Compute first deriv
        if (deriv >= 2) {
            if (meta_) {
                throw PSIEXCEPTION("TRYING TO COMPUTE MGGA FUNCTIONAL");

            } else if (gga_) {
                // spin polarized
                std::vector<double> fv(npoints);
                std::vector<double> frho(npoints * 2);
                std::vector<double> fsigma(npoints * 3);

                for (size_t i = 0; i < npoints; i++) {
                    frho[2 * i] = rho_ap[i];
                    frho[2 * i + 1] = rho_bp[i];

                    fsigma[3 * i] = gamma_aap[i];
                    fsigma[3 * i + 1] = gamma_abp[i];
                    fsigma[3 * i + 2] = gamma_bbp[i];
                }

                std::vector<double> fv2_rho2(npoints * 3);
                std::vector<double> fv2_rhosigma(npoints * 6);
                std::vector<double> fv2_sigma2(npoints * 6);

                xc_gga_fxc(&xc_functional_, npoints, frho.data(), fsigma.data(), fv2_rho2.data(),
                           fv2_rhosigma.data(), fv2_sigma2.data());

                for (size_t i = 0; i < npoints; i++) {
                    // v2rho2(3)       = (u_u, u_d, d_d)
                    v_rho_a_rho_a[i] += alpha_ * fv2_rho2[3 * i];
                    v_rho_a_rho_b[i] += alpha_ * fv2_rho2[3 * i + 1];
                    v_rho_b_rho_b[i] += alpha_ * fv2_rho2[3 * i + 2];

                    // v2sigma2(6)     = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd,
                    // dd_dd)
                    v_gamma_aa_gamma_aa[i] += alpha_ * fv2_sigma2[6 * i];
                    v_gamma_aa_gamma_ab[i] += alpha_ * fv2_sigma2[6 * i + 1];
                    v_gamma_aa_gamma_bb[i] += alpha_ * fv2_sigma2[6 * i + 2];
                    v_gamma_ab_gamma_ab[i] += alpha_ * fv2_sigma2[6 * i + 3];
                    v_gamma_ab_gamma_bb[i] += alpha_ * fv2_sigma2[6 * i + 4];
                    v_gamma_bb_gamma_bb[i] += alpha_ * fv2_sigma2[6 * i + 5];

                    // v2rhosigma(6)   = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
                    v_rho_a_gamma_aa[i] += alpha_ * fv2_rhosigma[6 * i];
                    v_rho_a_gamma_ab[i] += alpha_ * fv2_rhosigma[6 * i + 1];
                    v_rho_a_gamma_bb[i] += alpha_ * fv2_rhosigma[6 * i + 2];
                    v_rho_b_gamma_ab[i] += alpha_ * fv2_rhosigma[6 * i + 3];
                    v_rho_b_gamma_bb[i] += alpha_ * fv2_rhosigma[6 * i + 4];
                    v_rho_b_gamma_bb[i] += alpha_ * fv2_rhosigma[6 * i + 5];
                }

            } else {
                // spin polarized
                std::vector<double> frho(npoints * 2);

                for (size_t i = 0; i < npoints; i++) {
                    frho[2 * i] = rho_ap[i];
                    frho[2 * i + 1] = rho_bp[i];
                }

                std::vector<double> fv2_rho2(npoints * 3);

                xc_lda_fxc(&xc_functional_, npoints, frho.data(), fv2_rho2.data());

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
