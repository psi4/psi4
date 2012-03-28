#include <libmints/vector.h>
#include "wpbex_functional.h"
#include "utility.h"
#include <psi4-dec.h>
#include <cmath>

using namespace psi;

namespace psi {

wPBEXFunctional::wPBEXFunctional()
{
    common_init();
}
wPBEXFunctional::~wPBEXFunctional()
{
}
void wPBEXFunctional::common_init()
{
    name_ = "   wPBE_X";
    description_ = "    Short-Range PBE Exchange Hole (Exact Model)";
    citation_ = "   NULL"; 
    
    alpha_ = 1.0;
    omega_ = 0.3;
    
    gga_ = true;
    lrc_ = true;
    meta_ = false; 

    _K0_ = 3.0 * pow(3.0 / (4.0 * M_PI), 1.0/3.0);
    _k0_ = pow(6.0 * M_PI * M_PI, 1.0/3.0);
    
}
void wPBEXFunctional::set_parameter(const std::string& key, double val) 
{
    throw PSIEXCEPTION("wPBEX functional does not contain any adjustable parameters");
}
void wPBEXFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    compute_sigma_functional(in,out,npoints,deriv,alpha,true);
    compute_sigma_functional(in,out,npoints,deriv,alpha,false);
}
void wPBEXFunctional::compute_sigma_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("wPBEXFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

    // => Input variables (spin-polarized) <= //

    double* rho_s = NULL;
    double* gamma_s = NULL;
    double* tau_s = NULL;
    rho_s = in.find(spin ? "RHO_A" : "RHO_B")->second->pointer();
    gamma_s = in.find(spin ? "GAMMA_AA" : "GAMMA_BB")->second->pointer();

    // => Output variables <= //

    double* v = NULL;
    double* v_rho = NULL;
    double* v_gamma = NULL;
    
    v = out.find(spin ? "V" : "V")->second->pointer();
    if (deriv >=1) {
        v_rho = out.find(spin ? "V_RHO_A" : "V_RHO_B")->second->pointer();
        v_gamma = out.find(spin ? "V_GAMMA_AA" : "V_GAMMA_BB")->second->pointer();
    }
     
    // => Main Loop over points <= //
    for (int Q = 0; Q < npoints; Q++) {
    
        // => Primitive variables <= //
        double rho;
        double gamma;

        rho = rho_s[Q];
        if (rho < lsda_cutoff_) {
            continue;
        } 
        gamma = gamma_s[Q];

        // Powers of rho
        double rho13 = pow(rho,1.0/3.0);
        double rho43 = rho * rho13;
        double rho73 = rho * rho * rho13;

        //=>  Factors <= //

        // > LSDA < //
        double E, E_rho;
        E = - 0.5 * _K0_ * rho43;
        E_rho = -4.0/6.0 * rho13;
    
        // > Mixed LSDA/GGA spin-enhancement factor < //

        // s
        double s, s_rho, s_gamma;
        s = sqrt(gamma) / rho43;
        s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
        s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) / rho43; 

        // Spin unpolarized intermediates
        double P = 2.0 * rho;
        double S = 1.0 / (2.0 * _k0_) * s;

        // Omega is all weird in this case
        double omega = omega_ * pow(2.0, -1.0/6.0);

        // F
        double F, F_P, F_S;

        // Call Dr. Scuseria's nice function
        wpbe_F(rho, s, omega, &F, &F_P, &F_S);

        double F_rho = 2.0 * F_P;
        double F_s = 1.0 / (2.0 * _k0_) * F_S;

        // => Assembly <= //
        v[Q] += A * E * F;
        if (deriv >= 1) {
            v_rho[Q] += A * (F * E_rho + E * F_rho + E * F_s * s_rho);
            v_gamma[Q] += A * (E * F_s * s_gamma);
        }
    }
}

}
