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
    description_ = "    Short-Range PBE Exchange Hole (HJS Model)\n";
    citation_ = "    Henderson, Janesko, and Scuseria, J. Chem. Phys., 128, 194105 (2008)\nWeintrab, Henderson, and Scuseria, J. Chem. Theory. Comput., 5, 754 (2009)\n"; 
    
    alpha_ = 1.0;
    omega_ = 0.3;
    
    gga_ = true;
    lrc_ = true;
    meta_ = false; 

    _K0_ = 3.0 * pow(3.0 / (4.0 * M_PI), 1.0/3.0);
    _k0_ = pow(6.0 * M_PI * M_PI, 1.0/3.0);
    _s0_ = 2.0; 
    _pi12_ = sqrt(M_PI);
    _s_min_tol_ = 1.0E-10;
    _nu_min_tol_ = 1.0E-10;

    B88_ = false;
}
void wPBEXFunctional::set_parameter(const std::string& key, double val) 
{
    parameters_[key] = val;
    if (key == "A") {
        _A_ = val;
    } else if (key == "B") {
        _B_ = val;
    } else if (key == "C") {
        _C_ = val;
    } else if (key == "D") {
        _D_ = val;
    } else if (key == "E") {
        _E_ = val;
    } else if (key.substr(0,2) == "Ha") {
        int index = atoi(key.substr(2).c_str());
        if (_Ha_.size() < index + 1) {
            _Ha_.resize(index + 1);
            _Ha_[index] = val;
        }
    } else if (key.substr(0,2) == "Hb") {
        int index = atoi(key.substr(2).c_str());
        if (_Hb_.size() < index + 1) {
            _Hb_.resize(index + 1);
            _Hb_[index] = val;
        }
    } else {
        throw PSIEXCEPTION("Error, unknown HJS exchange functional parameter");    
    }
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

        // =>  Factors <= //

        // > LSDA < //
        double E, E_rho;
        E = - 0.5 * _K0_ * rho43;
        E_rho = - 4.0/6.0 * _K0_ * rho13;
    
        // > Mixed LSDA/GGA spin-enhancement factor < //

        // Targets
        double F;
        double F_rho;
        double F_s;
    
        // Intermediates for S

        // s
        double s, s_rho, s_gamma;
        s = sqrt(gamma) / rho43;
        s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
        s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) / rho43; 

        // Spin unpolarized intermediate
        double S = 1.0 / (2.0 * _k0_) * s;
        double S_s = 1.0 / (2.0 * _k0_);

        // Rescaled reduced density gradient (B88)
        double Sigma = S;
        double Sigma_S = 1.0;
        if (B88_) {
            if (S < 1.0E2) {
                double xi = 1.0 / (exp(20.0) - 1.0);
                Sigma = -log((exp(-S) + xi)/(1.0 + xi));
                Sigma_S = 1.0 / (xi * exp(S) + 1);    
            } else {
                Sigma = 20.0;
                Sigma_S = 0.0;
            }
        } 

        // Intermediates for RHO 

        double nu = omega_ / (_k0_ * rho13);
        double nu_rho = -1.0/3.0 * omega_ / (_k0_ * rho43); 

        double F_nu;   
        double F_Sigma;
        hjs_F(Sigma, nu, &F, &F_Sigma, &F_nu);

        double F_S = F_Sigma * Sigma_S;
        F_s = F_S * S_s;
        F_rho = F_nu * nu_rho;       

        // => Assembly <= //
        v[Q] += A * E * F;
        if (deriv >= 1) {
            v_rho[Q] += A * (F * E_rho + E * F_rho + E * F_s * s_rho);
            v_gamma[Q] += A * (E * F_s * s_gamma);
        }
    }
}

void wPBEXFunctional::hjs_F(double s, double nu, double* Fhjs, double *Fhjs_s, double* Fhjs_nu) 
{
    // ==> Small s (LSDA) limit <== //
    if (s <= _s_min_tol_) {

        *Fhjs_s = 0.0;

        // ==> Small nu limit <== //
        if (nu <= _nu_min_tol_) {
            *Fhjs= 1.0;
            *Fhjs_nu = -2.3632718012073548E+00; // Magic limit
            return;
        }

        double a = 1.0/2.0 * nu;
        double a_nu = 1.0/2.0; 

        // Easy part
        double F1k = 0.0;
        double F1k_a = 0.0;
        // Hard part
        double F2k = 0.0;
        double F2k_a = 0.0;
        
        // Common intermediates
        double x = -1.0/(4.0 * a * a);
        double expv = exp(x);
        double erfv = erf(1.0/(2.0*a));

        // Easy part (seems very stable)
        F1k = 1.0 - 8.0 / 3.0 * a * _pi12_ * erfv;
        F1k_a = -8.0/3.0 * (_pi12_ * erfv - expv / a);

        // Hard part (quite unstable)
        if (a < 10.0) {
            // Conventional
            double b = expv - 1.0;
            double c = 2.0* a * a * b + 0.5; 
            F2k = -8.0 / 3.0 * a * (2.0 * a * (b - c));
            F2k_a = -4.0 / (3.0 * x * x) * (b - c) +
                4.0 / (3.0 * x) * (expv - 1.0 / (2.0 * x * x) * (expv - 1.0) + 1.0 / (2.0 * x) * expv);

            // Chain rule
            F2k_a *= 0.5 / (a * a * a); 

        } else {
            // Taylor Series (sixth terms guarantees 10^-15 relative accuracy for a > 10.0)
            int N = 7;
            double xp = 1.0;
            double xpm1 = 0.0;
            int factp1 = 1;
            int factp2 = 2;

            for (int i = 0; i < N; i++) {
                F2k += 4.0 / 3.0 * (xp / factp1 + 0.5 * xp / factp2);
                F2k_a += 4.0 / 3.0 * (i * xpm1 / factp1 + 0.5 * i * xpm1 / factp2);
                factp1 *= i + 2;
                factp2 *= i + 3;
                xpm1 = xp;
                xp *= x;
            }
        
            // Chain rule
            F2k_a *= 0.5 / (a * a * a); 
        }
        
        *Fhjs = F1k + F2k;
        *Fhjs_nu = (F1k_a + F2k_a) * a_nu;
        return;
    }

    // ==> HJS Implementation of wPBE <== //

    //  > H  and H_s < //
    
    double H; 
    double H_s; 
    {
        double HN = 0.0;
        double HD = 0.0;
        double HN_s = 0.0;
        double HD_s = 0.0;
        
        double temp1, temp0;
        
        temp1 = 1.0;
        temp0 = 0.0;
        for (int i = 0; i < _Ha_.size(); i++) {
            HN += _Ha_[i] * temp1;
            HN_s += _Ha_[i] * i * temp0;
            temp0 = temp1;
            temp1 *= s;
        } 
        temp1 = 1.0;
        temp0 = 0.0;
        for (int i = 0; i < _Hb_.size(); i++) {
            HD += _Hb_[i] * temp1;
            HD_s += _Hb_[i] * i * temp0;
            temp0 = temp1;
            temp1 *= s;
        }
    
        H = HN / HD;
        H_s = (HN_s * HD - HN * HD_s) / (HD * HD);
    }
   
    //  > xi < //
    
    double xi = H*(s*s);
    
    //  > eta < //
    
    double eta = _A_+xi;
    
    //  > lambda < //
    
    double lambda = _D_+xi;
    
    //  > R < //
    
    double R = (_A_*pow(lambda,7.0/2.0)*(1.2E1/5.0))/(sqrt(eta)+sqrt(xi));
    
    //  > F < //
    
    double F = (xi*(-1.0/2.0))/_C_-((s*s)*(1.0/2.7E1))/(_C_*((s*s)*1.0/(_s0_*_s0_)+1.0))+1.0;
    
    //  > G < //
    
    double G = -(-R+_A_*(lambda*lambda*lambda)*(6.0/5.0)+_B_*(lambda*lambda)*(4.0/1.5E1)+pow(lambda,7.0/2.0)*_pi12_*(4.0/5.0)+_C_*F*lambda*(2.0/5.0))/_E_;
    
    //  > chi < //
    
    double chi = nu*1.0/sqrt(lambda+nu*nu);
    
    //  > L < //
    
    double L = (lambda+nu*nu)*(chi+1.0);
    
    //  > U < //
    
    double U = (_B_*(-4.0/9.0))/L-_C_*F*1.0/(L*L)*(chi*(1.0/2.0)+1.0)*(4.0/9.0)-_E_*G*1.0/(L*L*L)*(chi*(9.0/8.0)+(chi*chi)*(3.0/8.0)+1.0)*(8.0/9.0);
    
    //  > Dxe < //
    
    double Dxe = (sqrt(xi+nu*nu)+sqrt(eta+nu*nu))*(nu+sqrt(xi+nu*nu));
    
    //  > Dex < //
    
    double Dex = (sqrt(xi+nu*nu)+sqrt(eta+nu*nu))*(nu+sqrt(eta+nu*nu));
    
    //  > Q < //
    
    double Q = _A_*(eta/Dex+xi/Dxe);
    
    //  > Dlx < //
    
    double Dlx = (sqrt(xi+nu*nu)+sqrt(lambda+nu*nu))*(nu+sqrt(lambda+nu*nu));
    
    //  > Dle < //
    
    double Dle = (sqrt(eta+nu*nu)+sqrt(lambda+nu*nu))*(nu+sqrt(lambda+nu*nu));
    
    //  > Xx < //
    
    double Xx = (lambda-xi)/Dlx;
    
    //  > Xe < //
    
    double Xe = -(eta-lambda)/Dle;
    
    //  > Lx < //
    
    double Lx = log(-Xx+1.0);
    if (fabs(Xx) < 0.001) {
        double temp = Xx;
        Lx = 0.0;
        for (int k = 1; k <= 5; k++) {
            Lx -= temp / (double) k;
            temp *= Xx;
        }
    }   
 
    //  > Le < //
    
    double Le = log(-Xe+1.0);
    if (fabs(Xe) < 0.001) {
        double temp = Xe;
        Le = 0.0;
        for (int k = 1; k <= 5; k++) {
            Le -= temp / (double) k;
            temp *= Xe;
        }
    }   
    
    //  > V < //
    
    double V = Le*eta*-2.0+Lx*xi*2.0;
    
    //  > Z < //
    
    double Z = Q+U+V;
    
    //  > Fhjs < //
    
    *Fhjs = Z;
    
    //  > xi_s < //
    
    double xi_s = H*s*2.0;
    
    //  > xi_H < //
    
    double xi_H = s*s;
    
    //  > eta_xi < //
    
    double eta_xi = 1.0;
    
    //  > lambda_xi < //
    
    double lambda_xi = 1.0;
    
    //  > R_xi < //
    
    double R_xi = _A_*pow(lambda,7.0/2.0)*1.0/sqrt(xi)*1.0/pow(sqrt(eta)+sqrt(xi),2.0)*(-6.0/5.0);
    
    //  > R_eta < //
    
    double R_eta = _A_*1.0/sqrt(eta)*pow(lambda,7.0/2.0)*1.0/pow(sqrt(eta)+sqrt(xi),2.0)*(-6.0/5.0);
    
    //  > R_lambda < //
    
    double R_lambda = (_A_*pow(lambda,5.0/2.0)*(4.2E1/5.0))/(sqrt(eta)+sqrt(xi));
    
    //  > F_s < //
    
    double F_s = (s*(-2.0/2.7E1))/(_C_*((s*s)*1.0/(_s0_*_s0_)+1.0))+((s*s*s)*1.0/(_s0_*_s0_)*1.0/pow((s*s)*1.0/(_s0_*_s0_)+1.0,2.0)*(2.0/2.7E1))/_C_;
    
    //  > F_xi < //
    
    double F_xi = (-1.0/2.0)/_C_;
    
    //  > G_lambda < //
    
    double G_lambda = -(_B_*lambda*(8.0/1.5E1)+_A_*(lambda*lambda)*(1.8E1/5.0)+pow(lambda,5.0/2.0)*_pi12_*(1.4E1/5.0)+_C_*F*(2.0/5.0))/_E_;
    
    //  > G_R < //
    
    double G_R = 1.0/_E_;
    
    //  > G_F < //
    
    double G_F = (_C_*lambda*(-2.0/5.0))/_E_;
    
    //  > chi_nu < //
    
    double chi_nu = -(nu*nu)*1.0/pow(lambda+nu*nu,3.0/2.0)+1.0/sqrt(lambda+nu*nu);
    
    //  > chi_lambda < //
    
    double chi_lambda = nu*1.0/pow(lambda+nu*nu,3.0/2.0)*(-1.0/2.0);
    
    //  > L_nu < //
    
    double L_nu = nu*(chi+1.0)*2.0;
    
    //  > L_lambda < //
    
    double L_lambda = chi+1.0;
    
    //  > L_chi < //
    
    double L_chi = lambda+nu*nu;
    
    //  > U_F < //
    
    double U_F = _C_*1.0/(L*L)*(chi*(1.0/2.0)+1.0)*(-4.0/9.0);
    
    //  > U_G < //
    
    double U_G = _E_*1.0/(L*L*L)*(chi*(9.0/8.0)+(chi*chi)*(3.0/8.0)+1.0)*(-8.0/9.0);
    
    //  > U_chi < //
    
    double U_chi = _C_*F*1.0/(L*L)*(-2.0/9.0)-_E_*G*1.0/(L*L*L)*(chi*(3.0/4.0)+9.0/8.0)*(8.0/9.0);
    
    //  > U_L < //
    
    double U_L = _B_*1.0/(L*L)*(4.0/9.0)+_C_*F*1.0/(L*L*L)*(chi*(1.0/2.0)+1.0)*(8.0/9.0)+_E_*G*1.0/(L*L*L*L)*(chi*(9.0/8.0)+(chi*chi)*(3.0/8.0)+1.0)*(8.0/3.0);
    
    //  > Dxe_nu < //
    
    double Dxe_nu = (nu+sqrt(xi+nu*nu))*(nu*1.0/sqrt(eta+nu*nu)+nu*1.0/sqrt(xi+nu*nu))+(nu*1.0/sqrt(xi+nu*nu)+1.0)*(sqrt(xi+nu*nu)+sqrt(eta+nu*nu));
    
    //  > Dxe_xi < //
    
    double Dxe_xi = 1.0/sqrt(xi+nu*nu)*(nu+sqrt(xi+nu*nu))*(1.0/2.0)+(sqrt(xi+nu*nu)+sqrt(eta+nu*nu))*1.0/sqrt(xi+nu*nu)*(1.0/2.0);
    
    //  > Dxe_eta < //
    
    double Dxe_eta = 1.0/sqrt(eta+nu*nu)*(nu+sqrt(xi+nu*nu))*(1.0/2.0);
    
    //  > Dex_nu < //
    
    double Dex_nu = (nu+sqrt(eta+nu*nu))*(nu*1.0/sqrt(eta+nu*nu)+nu*1.0/sqrt(xi+nu*nu))+(nu*1.0/sqrt(eta+nu*nu)+1.0)*(sqrt(xi+nu*nu)+sqrt(eta+nu*nu));
    
    //  > Dex_xi < //
    
    double Dex_xi = 1.0/sqrt(xi+nu*nu)*(nu+sqrt(eta+nu*nu))*(1.0/2.0);
    
    //  > Dex_eta < //
    
    double Dex_eta = 1.0/sqrt(eta+nu*nu)*(nu+sqrt(eta+nu*nu))*(1.0/2.0)+(sqrt(xi+nu*nu)+sqrt(eta+nu*nu))*1.0/sqrt(eta+nu*nu)*(1.0/2.0);
    
    //  > Q_xi < //
    
    double Q_xi = _A_/Dxe;
    
    //  > Q_eta < //
    
    double Q_eta = _A_/Dex;
    
    //  > Q_Dxe < //
    
    double Q_Dxe = -_A_*1.0/(Dxe*Dxe)*xi;
    
    //  > Q_Dex < //
    
    double Q_Dex = -_A_*1.0/(Dex*Dex)*eta;
    
    //  > Dlx_nu < //
    
    double Dlx_nu = (nu+sqrt(lambda+nu*nu))*(nu*1.0/sqrt(lambda+nu*nu)+nu*1.0/sqrt(xi+nu*nu))+(nu*1.0/sqrt(lambda+nu*nu)+1.0)*(sqrt(xi+nu*nu)+sqrt(lambda+nu*nu));
    
    //  > Dlx_xi < //
    
    double Dlx_xi = 1.0/sqrt(xi+nu*nu)*(nu+sqrt(lambda+nu*nu))*(1.0/2.0);
    
    //  > Dlx_lambda < //
    
    double Dlx_lambda = 1.0/sqrt(lambda+nu*nu)*(nu+sqrt(lambda+nu*nu))*(1.0/2.0)+(sqrt(xi+nu*nu)+sqrt(lambda+nu*nu))*1.0/sqrt(lambda+nu*nu)*(1.0/2.0);
    
    //  > Dle_nu < //
    
    double Dle_nu = (nu+sqrt(lambda+nu*nu))*(nu*1.0/sqrt(eta+nu*nu)+nu*1.0/sqrt(lambda+nu*nu))+(nu*1.0/sqrt(lambda+nu*nu)+1.0)*(sqrt(eta+nu*nu)+sqrt(lambda+nu*nu));
    
    //  > Dle_eta < //
    
    double Dle_eta = 1.0/sqrt(eta+nu*nu)*(nu+sqrt(lambda+nu*nu))*(1.0/2.0);
    
    //  > Dle_lambda < //
    
    double Dle_lambda = 1.0/sqrt(lambda+nu*nu)*(nu+sqrt(lambda+nu*nu))*(1.0/2.0)+(sqrt(eta+nu*nu)+sqrt(lambda+nu*nu))*1.0/sqrt(lambda+nu*nu)*(1.0/2.0);
    
    //  > Xx_xi < //
    
    double Xx_xi = -1.0/Dlx;
    
    //  > Xx_lambda < //
    
    double Xx_lambda = 1.0/Dlx;
    
    //  > Xx_Dlx < //
    
    double Xx_Dlx = -1.0/(Dlx*Dlx)*(lambda-xi);
    
    //  > Xe_eta < //
    
    double Xe_eta = -1.0/Dle;
    
    //  > Xe_lambda < //
    
    double Xe_lambda = 1.0/Dle;
    
    //  > Xe_Dle < //
    
    double Xe_Dle = 1.0/(Dle*Dle)*(eta-lambda);
    
    //  > Lx_Xx < //
    
    double Lx_Xx = 1.0/(Xx-1.0);
    
    //  > Le_Xe < //
    
    double Le_Xe = 1.0/(Xe-1.0);
    
    //  > V_xi < //
    
    double V_xi = Lx*2.0;
    
    //  > V_eta < //
    
    double V_eta = Le*-2.0;
    
    //  > V_Lx < //
    
    double V_Lx = xi*2.0;
    
    //  > V_Le < //
    
    double V_Le = eta*-2.0;
    
    //  > Z_U < //
    
    double Z_U = 1.0;
    
    //  > Z_Q < //
    
    double Z_Q = 1.0;
    
    //  > Z_V < //
    
    double Z_V = 1.0;
    
    //  > Fhjs_Z < //
    
    double Fhjs_Z = 1.0;
    
    //  > Fhjs_V < //
    
    double Fhjs_V = Fhjs_Z*Z_V;
    
    //  > Fhjs_Le < //
    
    double Fhjs_Le = Fhjs_V*V_Le;
    
    //  > Fhjs_Lx < //
    
    double Fhjs_Lx = Fhjs_V*V_Lx;
    
    //  > Fhjs_Xe < //
    
    double Fhjs_Xe = Fhjs_Le*Le_Xe;
    
    //  > Fhjs_Xx < //
    
    double Fhjs_Xx = Fhjs_Lx*Lx_Xx;
    
    //  > Fhjs_Dle < //
    
    double Fhjs_Dle = Fhjs_Xe*Xe_Dle;
    
    //  > Fhjs_Dlx < //
    
    double Fhjs_Dlx = Fhjs_Xx*Xx_Dlx;
    
    //  > Fhjs_Q < //
    
    double Fhjs_Q = Fhjs_Z*Z_Q;
    
    //  > Fhjs_Dex < //
    
    double Fhjs_Dex = Fhjs_Q*Q_Dex;
    
    //  > Fhjs_Dxe < //
    
    double Fhjs_Dxe = Fhjs_Q*Q_Dxe;
    
    //  > Fhjs_U < //
    
    double Fhjs_U = Fhjs_Z*Z_U;
    
    //  > Fhjs_L < //
    
    double Fhjs_L = Fhjs_U*U_L;
    
    //  > Fhjs_chi < //
    
    double Fhjs_chi = Fhjs_L*L_chi+Fhjs_U*U_chi;
    
    //  > Fhjs_G < //
    
    double Fhjs_G = Fhjs_U*U_G;
    
    //  > Fhjs_F < //
    
    double Fhjs_F = Fhjs_G*G_F+Fhjs_U*U_F;
    
    //  > Fhjs_R < //
    
    double Fhjs_R = Fhjs_G*G_R;
    
    //  > Fhjs_lambda < //
    
    double Fhjs_lambda = Fhjs_chi*chi_lambda+Dle_lambda*Fhjs_Dle+Dlx_lambda*Fhjs_Dlx+Fhjs_G*G_lambda+Fhjs_L*L_lambda+Fhjs_R*R_lambda+Fhjs_Xe*Xe_lambda+Fhjs_Xx*Xx_lambda;
    
    //  > Fhjs_eta < //
    
    double Fhjs_eta = Dex_eta*Fhjs_Dex+Dle_eta*Fhjs_Dle+Dxe_eta*Fhjs_Dxe+Fhjs_Q*Q_eta+Fhjs_R*R_eta+Fhjs_V*V_eta+Fhjs_Xe*Xe_eta;
    
    //  > Fhjs_xi < //
    
    double Fhjs_xi = Fhjs_eta*eta_xi+Fhjs_lambda*lambda_xi+Dex_xi*Fhjs_Dex+Dlx_xi*Fhjs_Dlx+Dxe_xi*Fhjs_Dxe+F_xi*Fhjs_F+Fhjs_Q*Q_xi+Fhjs_R*R_xi+Fhjs_V*V_xi+Fhjs_Xx*Xx_xi;
    
    //  > Fhjs_H < //
    
    double Fhjs_H = Fhjs_xi*xi_H;
    
    //  > Fhjs_s < //
    
    *Fhjs_s = Fhjs_xi*xi_s+F_s*Fhjs_F+Fhjs_H*H_s;
    if ((s > 1.0E5 && nu > 1.0E5) || s > 1.0E10 || nu > 1.0E8) {
        *Fhjs_s = 0.0;
    }
    
    //  > Fhjs_nu < //
    
    *Fhjs_nu = Fhjs_chi*chi_nu+Dex_nu*Fhjs_Dex+Dle_nu*Fhjs_Dle+Dlx_nu*Fhjs_Dlx+Dxe_nu*Fhjs_Dxe+Fhjs_L*L_nu;
     
}


}
