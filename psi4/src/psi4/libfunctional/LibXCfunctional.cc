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
#include "psi4/libmints/vector.h"
#include "LibXCfunctional.h"
#include "psi4/psi4-dec.h"

using namespace psi;

namespace psi {

LibXCFunctional::LibXCFunctional(std::string xc_name)
{

    func_id_ = xc_functional_get_number(xc_name.c_str());
    if(xc_func_init(&xc_functional_, func_id_, XC_POLARIZED) != 0){
        outfile->Printf("Functional '%d' not found\n", func_id_);
        throw PSIEXCEPTION("Could not find required LIBXC functional");
    }

    name_ = std::string(xc_functional_.info->name);
    // outfile->Printf("I found the LIBXC functional, the name is %s", name_.c_str());
    // description_ = str(xc_functional_->name);
    // citation_ = std::string(xc_functional_.info->refs[0]);

    // int kind = xc_functional_->kind;
    // int number = xc_functional_->number;
    int family = xc_functional_.info->family;
    outfile->Printf("The family of %s is %d\n", name_.c_str(), family);

    gga_ = (family >= 2);
    meta_ = (family >= 4);

}
LibXCFunctional::~LibXCFunctional()
{
    xc_func_end(&xc_functional_);
}
void LibXCFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{

    // Overall scale factor
    double scale = alpha_ * alpha;

    // => Input variables <= //

    double* rho_ap = NULL;
    double* rho_bp = NULL;
    double* gamma_aap = NULL;
    double* gamma_abp = NULL;
    double* gamma_bbp = NULL;
    double* tau_ap = NULL;
    double* tau_bp = NULL;

    // gga_ = true;

    if (true) {
        rho_ap = in.find("RHO_A")->second->pointer();
        rho_bp = in.find("RHO_B")->second->pointer();
    }
    if (gga_) {
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        gamma_abp = in.find("GAMMA_AB")->second->pointer();
        gamma_bbp = in.find("GAMMA_BB")->second->pointer();
        // outfile->Printf("% 9.7f %9.7f %9.7f\n", gamma_aap[0], gamma_abp[0], gamma_bbp[0]);
    }
    if (meta_)  {
        tau_ap = in.find("TAU_A")->second->pointer();
        tau_bp = in.find("TAU_B")->second->pointer();
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
        }
    }

    if (deriv == 0){
        if (meta_){
            // outfile->Printf("Executing MGGA");
            // xc_mgga_exc(&xc_functional_, npoints, rho_ap, gamma_aap, v);
            throw PSIEXCEPTION("TRYING TO COMPUTE MGGA FUNCTIONAL");

        } else if (gga_) {
            // outfile->Printf("Executing GGA");
            throw PSIEXCEPTION("NYI");
            xc_gga_exc(&xc_functional_, npoints, rho_ap, gamma_aap, v);

        } else{
            // outfile->Printf("Executing LDA");
            throw PSIEXCEPTION("NYI");
            xc_lda_exc(&xc_functional_, npoints, rho_ap, v);
        }
    } else if (deriv == 1) {
        if (meta_){
            outfile->Printf("Executing MGGA");
            // xc_mgga_exc(&xc_functional_, npoints, rho_ap, gamma_aap, v);
            throw PSIEXCEPTION("TRYING TO COMPUTE MGGA FUNCTIONAL");

        } else if (gga_) {
            // xc_gga_exc_vxc(&xc_functional_, npoints, rho_ap, gamma_aap, v, v_rho_a, v_gamma_aa);
            outfile->Printf("Executing GGA\n");

            // spin polarized
            std::vector<double> fv(npoints);
            std::vector<double> frho(npoints*2);
            std::vector<double> fsigma(npoints*3);

            for (size_t i=0; i < npoints; i++){
                frho[2 * i] = rho_ap[i];
                frho[2 * i + 1] = rho_bp[i];

                fsigma[3 * i] = gamma_aap[i];
                fsigma[3 * i + 1] = gamma_abp[i];
                fsigma[3 * i + 2] = gamma_bbp[i];
            }


            std::vector<double> fv_rho(npoints*2);
            std::vector<double> fv_sigma(npoints*3);

            xc_gga_exc_vxc(&xc_functional_, npoints, frho.data(), fsigma.data(), fv.data(), fv_rho.data(), fv_sigma.data());

            for (size_t i=0; i < npoints; i++){
                v[i] += fv[i] * (rho_ap[i] + rho_bp[i]);
                v_rho_a[i] += fv_rho[2 * i];
                v_rho_b[i] += fv_rho[2 * i + 1];
                v_gamma_aa[i] += fv_sigma[3 * i];
                v_gamma_ab[i] += fv_sigma[3 * i + 1];
                v_gamma_bb[i] += fv_sigma[3 * i + 2];

            }

        } else{
            // xc_lda_exc_vxc(&xc_functional_, npoints, rho_ap, v, v_rho_a);
            outfile->Printf("Executing LDA\n");

            // spin polarized
            std::vector<double> fv(npoints);
            std::vector<double> frho(npoints*2);
            std::vector<double> fv_rho(npoints*2);

            for (size_t i=0; i < npoints; i++){
                frho[2 * i] = rho_ap[i];
                frho[2 * i + 1] = rho_bp[i];
            }

            xc_lda_exc_vxc(&xc_functional_, npoints, frho.data(), fv.data(), fv_rho.data());


            for (size_t i=0; i < npoints; i++){
                v[i] += fv[i] * (rho_ap[i] + rho_bp[i]);
                v_rho_a[i] += fv_rho[2 * i];
                v_rho_b[i] += fv_rho[2 * i + 1];
            }
        }
    } else {
        throw PSIEXCEPTION("TRYING TO COPMUTE DERIV > 1 ");
    }



}

} // End namespace