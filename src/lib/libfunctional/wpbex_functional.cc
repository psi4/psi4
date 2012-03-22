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
void wPBEXFunctional::computeRKSFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    computeFunctional(in,out,npoints,deriv,alpha,true, true);
    computeFunctional(in,out,npoints,deriv,alpha,true, false);
}
void wPBEXFunctional::computeUKSFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    computeFunctional(in,out,npoints,deriv,alpha,true, true);
    computeFunctional(in,out,npoints,deriv,alpha,false, false);
}
void wPBEXFunctional::computeFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool in_alpha, bool out_alpha)
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
    rho_s = in.find(in_alpha ? "RHO_A" : "RHO_B")->second->pointer();
    gamma_s = in.find(in_alpha ? "GAMMA_AA" : "GAMMA_BB")->second->pointer();

    // => Output variables <= //

    double* v = NULL;
    double* v_rho = NULL;
    double* v_gamma = NULL;
    
    v = out.find(out_alpha ? "V" : "V")->second->pointer();
    if (deriv >=1) {
        v_rho = out.find(out_alpha ? "V_RHO_A" : "V_RHO_B")->second->pointer();
        v_gamma = out.find(out_alpha ? "V_GAMMA_AA" : "V_GAMMA_BB")->second->pointer();
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
        E = - _K0_ * rho43;
        E_rho = -4.0/3.0 * rho13;
    
        // > Mixed LSDA/GGA spin-enhancement factor < //

        // s
        double s, s_rho, s_gamma;
        s = sqrt(gamma) / rho43;
        s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
        s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) * rho43; 

        // F
        double F, F_rho, F_s;
        F = 1.0;
        F_rho = 0.0;
        F_s = 0.0;
        // TODO

        // => Assembly <= //
        v[Q] += E * F;
        if (deriv >= 1) {
            v_rho[Q] += A * (F * E_rho + E * F_rho + E * F_s * s_rho);
            v_gamma[Q] += A * (E_rho * F_s * s_gamma);
        }
    }
}

}
