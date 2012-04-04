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
    citation_ = "    Henderson, Janesko, and Scuseria, J. Chem. Phys., 128, 194105 (2008)\n"; 
    
    alpha_ = 1.0;
    omega_ = 0.3;
    
    gga_ = true;
    lrc_ = true;
    meta_ = false; 

    hjs_ = false;

    _K0_ = 3.0 * pow(3.0 / (4.0 * M_PI), 1.0/3.0);
    _k0_ = pow(6.0 * M_PI * M_PI, 1.0/3.0);
    _s0_ = 2.0; 
    _pi12_ = sqrt(M_PI);
    _s_min_tol_ = 1.0E-10;
    _nu_min_tol_ = 1.0E-10;
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

        //=>  Factors <= //

        // > LSDA < //
        double E, E_rho;
        E = - 0.5 * _K0_ * rho43;
        E_rho = - 4.0/6.0 * _K0_ * rho13;
    
        // > Mixed LSDA/GGA spin-enhancement factor < //

        // s
        double s, s_rho, s_gamma;
        s = sqrt(gamma) / rho43;
        s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
        s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) / rho43; 

        // Spin unpolarized intermediates
        double P = 2.0 * rho;
        double P_rho = 2.0;
        double S = 1.0 / (2.0 * _k0_) * s;
        double S_s = 1.0 / (2.0 * _k0_);

        double F;
        double F_rho;
        double F_s;

        // HJS Kernel
        if (hjs_) {
            // Nu intermediate
            double nu = omega_ / (_k0_ * rho13);
            double nu_rho = -1.0/3.0 * omega_ / (_k0_ * rho43); 
            double F_nu;   
            double F_S;
            hjs_F(S, nu, &F, &F_S, &F_nu);
            F_s = F_S * S_s;
            F_rho = F_nu * nu_rho;       
        } else {
            double F_P;
            double F_S;
            hse_F(P, S, omega_, &F, &F_P, &F_S);
            F_s = F_S * S_s;
            F_rho = F_P * P_rho;
        } 

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
    // ==> Small s limit <== //
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
    
    //  > alpha < //
    
    double alpha = -lambda/(lambda+nu*sqrt(lambda+nu*nu)+nu*nu);
    
    //  > Q < //
    
    double Q = -_A_*((nu*2.0)/(sqrt(xi+nu*nu)+sqrt(eta+nu*nu))-1.0);
    
    //  > V < //
    
    double V = eta*log((nu+sqrt(eta+nu*nu))/(nu+sqrt(lambda+nu*nu)))*-2.0+xi*log((nu+sqrt(xi+nu*nu))/(nu+sqrt(lambda+nu*nu)))*2.0;
    
    //  > U < //
    
    double U = (_B_*alpha*(4.0/9.0))/lambda-_C_*F*1.0/(lambda*lambda)*((alpha*alpha)*(3.0/2.0)+(alpha*alpha*alpha)*(1.0/2.0))*(4.0/9.0)+_E_*G*1.0/(lambda*lambda*lambda)*((alpha*alpha*alpha)*(5.0/2.0)+(alpha*alpha*alpha*alpha)*(1.5E1/8.0)+(alpha*alpha*alpha*alpha*alpha)*(3.0/8.0))*(8.0/9.0);
    
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
    
    //  > alpha_nu < //
    
    double alpha_nu = lambda*(nu*2.0+(nu*nu)*1.0/sqrt(lambda+nu*nu)+sqrt(lambda+nu*nu))*1.0/pow(lambda+nu*sqrt(lambda+nu*nu)+nu*nu,2.0);
    
    //  > alpha_lambda < //
    
    double alpha_lambda = -1.0/(lambda+nu*sqrt(lambda+nu*nu)+nu*nu)+lambda*(nu*1.0/sqrt(lambda+nu*nu)*(1.0/2.0)+1.0)*1.0/pow(lambda+nu*sqrt(lambda+nu*nu)+nu*nu,2.0);
    
    //  > Q_nu < //
    
    double Q_nu = -_A_*(2.0/(sqrt(xi+nu*nu)+sqrt(eta+nu*nu))-nu*1.0/pow(sqrt(xi+nu*nu)+sqrt(eta+nu*nu),2.0)*(nu*1.0/sqrt(eta+nu*nu)+nu*1.0/sqrt(xi+nu*nu))*2.0);
    
    //  > Q_xi < //
    
    double Q_xi = _A_*nu*1.0/pow(sqrt(xi+nu*nu)+sqrt(eta+nu*nu),2.0)*1.0/sqrt(xi+nu*nu);
    
    //  > Q_eta < //
    
    double Q_eta = _A_*nu*1.0/pow(sqrt(xi+nu*nu)+sqrt(eta+nu*nu),2.0)*1.0/sqrt(eta+nu*nu);
    
    //  > V_nu < //
    
    double V_nu = (eta*((nu*1.0/sqrt(eta+nu*nu)+1.0)/(nu+sqrt(lambda+nu*nu))-(nu*1.0/sqrt(lambda+nu*nu)+1.0)*(nu+sqrt(eta+nu*nu))*1.0/pow(nu+sqrt(lambda+nu*nu),2.0))*(nu+sqrt(lambda+nu*nu))*-2.0)/(nu+sqrt(eta+nu*nu))+(xi*((nu*1.0/sqrt(xi+nu*nu)+1.0)/(nu+sqrt(lambda+nu*nu))-(nu*1.0/sqrt(lambda+nu*nu)+1.0)*1.0/pow(nu+sqrt(lambda+nu*nu),2.0)*(nu+sqrt(xi+nu*nu)))*(nu+sqrt(lambda+nu*nu))*2.0)/(nu+sqrt(xi+nu*nu));
    
    //  > V_xi < //
    
    double V_xi = log((nu+sqrt(xi+nu*nu))/(nu+sqrt(lambda+nu*nu)))*2.0+(xi*1.0/sqrt(xi+nu*nu))/(nu+sqrt(xi+nu*nu));
    
    //  > V_eta < //
    
    double V_eta = log((nu+sqrt(eta+nu*nu))/(nu+sqrt(lambda+nu*nu)))*-2.0-(eta*1.0/sqrt(eta+nu*nu))/(nu+sqrt(eta+nu*nu));
    
    //  > V_lambda < //
    
    double V_lambda = (eta*1.0/sqrt(lambda+nu*nu))/(nu+sqrt(lambda+nu*nu))-(xi*1.0/sqrt(lambda+nu*nu))/(nu+sqrt(lambda+nu*nu));
    
    //  > U_lambda < //
    
    double U_lambda = _B_*alpha*1.0/(lambda*lambda)*(-4.0/9.0)+_C_*F*1.0/(lambda*lambda*lambda)*((alpha*alpha)*(3.0/2.0)+(alpha*alpha*alpha)*(1.0/2.0))*(8.0/9.0)-_E_*G*1.0/(lambda*lambda*lambda*lambda)*((alpha*alpha*alpha)*(5.0/2.0)+(alpha*alpha*alpha*alpha)*(1.5E1/8.0)+(alpha*alpha*alpha*alpha*alpha)*(3.0/8.0))*(8.0/3.0);
    
    //  > U_F < //
    
    double U_F = _C_*1.0/(lambda*lambda)*((alpha*alpha)*(3.0/2.0)+(alpha*alpha*alpha)*(1.0/2.0))*(-4.0/9.0);
    
    //  > U_G < //
    
    double U_G = _E_*1.0/(lambda*lambda*lambda)*((alpha*alpha*alpha)*(5.0/2.0)+(alpha*alpha*alpha*alpha)*(1.5E1/8.0)+(alpha*alpha*alpha*alpha*alpha)*(3.0/8.0))*(8.0/9.0);
    
    //  > U_alpha < //
    
    double U_alpha = (_B_*(4.0/9.0))/lambda+_E_*G*1.0/(lambda*lambda*lambda)*((alpha*alpha)*(1.5E1/2.0)+(alpha*alpha*alpha)*(1.5E1/2.0)+(alpha*alpha*alpha*alpha)*(1.5E1/8.0))*(8.0/9.0)-_C_*F*1.0/(lambda*lambda)*(alpha*3.0+(alpha*alpha)*(3.0/2.0))*(4.0/9.0);
    
    //  > Z_Q < //
    
    double Z_Q = 1.0;
    
    //  > Z_V < //
    
    double Z_V = 1.0;
    
    //  > Z_U < //
    
    double Z_U = 1.0;
    
    //  > Fhjs_Z < //
    
    double Fhjs_Z = 1.0;
    
    //  > Fhjs_U < //
    
    double Fhjs_U = Fhjs_Z*Z_U;
    
    //  > Fhjs_V < //
    
    double Fhjs_V = Fhjs_Z*Z_V;
    
    //  > Fhjs_Q < //
    
    double Fhjs_Q = Fhjs_Z*Z_Q;
    
    //  > Fhjs_alpha < //
    
    double Fhjs_alpha = Fhjs_U*U_alpha;
    
    //  > Fhjs_G < //
    
    double Fhjs_G = Fhjs_U*U_G;
    
    //  > Fhjs_F < //
    
    double Fhjs_F = Fhjs_G*G_F+Fhjs_U*U_F;
    
    //  > Fhjs_R < //
    
    double Fhjs_R = Fhjs_G*G_R;
    
    //  > Fhjs_lambda < //
    
    double Fhjs_lambda = Fhjs_alpha*alpha_lambda+Fhjs_G*G_lambda+Fhjs_R*R_lambda+Fhjs_U*U_lambda+Fhjs_V*V_lambda;
    
    //  > Fhjs_eta < //
    
    double Fhjs_eta = Fhjs_Q*Q_eta+Fhjs_R*R_eta+Fhjs_V*V_eta;
    
    //  > Fhjs_xi < //
    
    double Fhjs_xi = Fhjs_eta*eta_xi+Fhjs_lambda*lambda_xi+F_xi*Fhjs_F+Fhjs_Q*Q_xi+Fhjs_R*R_xi+Fhjs_V*V_xi;
    
    //  > Fhjs_H < //
    
    double Fhjs_H = Fhjs_xi*xi_H;
    
    //  > Fhjs_s < //
    
    *Fhjs_s = Fhjs_xi*xi_s+F_s*Fhjs_F+Fhjs_H*H_s;
    
    //  > Fhjs_nu < //
    
    *Fhjs_nu = Fhjs_alpha*alpha_nu+Fhjs_Q*Q_nu+Fhjs_V*V_nu;
    
}
void wPBEXFunctional::hse_F(double rho, double s, double omega, double *fx, double *d1rfx, double *d1sfx)
{
    /* Initialized data */
    static double zero = 0.;
    static double nine = 9.;
    static double ten = 10.;
    static double fifteen = 15.;
    static double sixteen = 16.;
    static double r36 = 36.;
    static double r64 = 64.;
    static double r81 = 81.;
    static double r256 = 256.;
    static double r384 = 384.;
    static double one = 1.;
    static double r27 = 27.;
    static double r128 = 128.;
    static double r144 = 144.;
    static double r288 = 288.;
    static double r324 = 324.;
    static double two = 2.;
    static double r729 = 729.;
    static double r20 = 20.;
    static double r32 = 32.;
    static double r243 = 243.;
    static double r2187 = 2187.;
    static double r6561 = 6561.;
    static double r40 = 40.;
    static double r12 = 12.;
    static double r25 = 25.;
    static double r30 = 30.;
    static double three = 3.;
    static double r54 = 54.;
    static double r75 = 75.;
    static double r105 = 105.;
    static double r135 = 135.;
    static double r1215 = 1215.;
    static double r15309 = 15309.;
    static double four = 4.;
    static double five = 5.;
    static double six = 6.;
    static double seven = 7.;
    static double eight = 8.;

    /* System generated locals */
    double d__1;

    /* Local variables */
    static double d1rexpei, d1sexpei, piexperf, a, b, c__, d__, e, f, h__,
	     w, x, a2, a3, a4, f2, f3, f4, f5, f6, f7, f8, f9, s2, s3, s4, s5,
	     w2, w3, w4, w5, w6, w7, w8, s6, t1, a12, a32, a52, f12, f13, f14,
	     f23, f43, f32, f18, f72, f34, a72, f94, eg, t10, piexperfd1, f98;
    static double pi, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, eb1, ha1, 
	    ha2, ha3, ha4, ha5, fc1, fc2, pi2, np1, np2, expfcutoff, g_a__, 
	    g_b__, f2d1, f3d1, f4d1, f5d1, f6d1, f1516, f8d1, f9d1, dhs, 
	    d1rpiexperf, xkf, d1spiexperf, t2t9, ega1, ega2, ega3, t10d1, 
	    d1sf, d1sh, dhs2, dhs3, dhs4, d1rw, pi_23__, hden, d1sf2, d1sf3, 
	    d1sf4, d1sf5, d1sf6, d1sf7, d1sf8, d1sf9, dhs72, d1rf2, dhs92, 
	    d1rf3, d1rf4, d1rf5, d1rf6, d1rf7, d1rf8, d1rf9, d1st1, d1rt1, 
	    hsbw;
    static double dhsw, hnum, srpi, d1seg, d1rt10, d1st10, hsbw2, hsbw3, 
	    hsbw4, dhsw2, term1, term2, term3, term4, term5;
    static double d1rnp1, hsbw12, dhsbw, hsbw32, d1rnp2, hsbw52, expei, 
	    hsbw72, dhsw52, dhsw72, d1sg_a__, d1sg_b__, d1sdhs, dhsbw2, 
	    dhsbw3, dhsbw4, dhsbw5, expei1, expei2, expei3, expei4, d1rt2t9, 
	    d1st2t9, dhsbw12, dhsbw32, term1d1, dhsbw52, dhsbw72, hsbwa94, 
	    dhsbw92, egscut, xkfrho, d1shden, expeid1, hsbwa942, hsbwa943, 
	    hsbwa945, d1rhsbw, d1shsbw, d1rdhsw, d1shnum, three_13__, 
	    hsbwa9412, d1rterm1, d1sterm1, d1sterm2, d1sterm3, d1rterm3, 
	    d1sterm4, d1rterm4, d1sterm5, d1rterm5, wcutoff;

    /* General constants */
    f12 = .5;
    f13 = one / three;
    f14 = .25;
    f18 = .125;
    f23 = two * f13;
    f43 = two * f23;
    f32 = 1.5;
    f72 = 3.5;
    f34 = .75;
    f94 = 2.25;
    f98 = 1.125;
    f1516 = fifteen / sixteen;
    pi = acos(-one);
    pi2 = pi * pi;
    pi_23__ = pow(pi2, f13);
    srpi = sqrt(pi);
    three_13__ = pow(three, f13);
    /* Constants from fit */
    ea1 = -1.128223946706117;
    ea2 = 1.452736265762971;
    ea3 = -1.243162299390327;
    ea4 = .971824836115601;
    ea5 = -.568861079687373;
    ea6 = .246880514820192;
    ea7 = -.065032363850763;
    ea8 = .008401793031216;
    eb1 = 1.455915450052607;
    /* Constants for PBE hole */
    a = 1.0161144;
    b = -.37170836;
    c__ = -.077215461;
    d__ = .57786348;
    e = -.051955731;
    x = -eight / nine;
    /* Constants for fit of H(s) (PBE) */
    ha1 = .00979681;
    ha2 = .0410834;
    ha3 = .18744;
    ha4 = .00120824;
    ha5 = .0347188;
    /* Constants for F(H) (PBE) */
    fc1 = 6.4753871;
    fc2 = .4796583;
    /* Constants for polynomial expansion for EG for small s */
    ega1 = -.0262841788;
    ega2 = -.07117647788;
    ega3 = .08534541323;
    /* Constants for large x expansion of exp(x)*ei(-x) */
    expei1 = 4.0364;
    expei2 = 1.15198;
    expei3 = 5.03627;
    expei4 = 4.1916;
    /* Cutoff criterion below which to use polynomial expansion */
    egscut = .08;
    wcutoff = 14.;
    expfcutoff = 700.;
    /* Calculate prelim variables */
    d__1 = three * pi2 * rho;
    xkf = pow(d__1, f13);
    xkfrho = xkf * rho;
    a2 = a * a;
    a3 = a2 * a;
    a4 = a3 * a;
    a12 = sqrt(a);
    a32 = a12 * a;
    a52 = a32 * a;
    a72 = a52 * a;
    w = omega / xkf;
    w2 = w * w;
    w3 = w2 * w;
    w4 = w2 * w2;
    w5 = w3 * w2;
    w6 = w5 * w;
    w7 = w6 * w;
    w8 = w7 * w;
    d1rw = -(one / (three * rho)) * w;
    x = -eight / nine;
    s2 = s * s;
    s3 = s2 * s;
    s4 = s2 * s2;
    s5 = s4 * s;
    s6 = s5 * s;
    /* Calculate wPBE enhancement factor */
    hnum = ha1 * s2 + ha2 * s4;
    hden = one + ha3 * s4 + ha4 * s5 + ha5 * s6;
    h__ = hnum / hden;
    d1shnum = two * ha1 * s + four * ha2 * s3;
    d1shden = four * ha3 * s3 + five * ha4 * s4 + six * ha5 * s5;
    d1sh = (hden * d1shnum - hnum * d1shden) / (hden * hden);
    f = fc1 * h__ + fc2;
    d1sf = fc1 * d1sh;
    /* Change exponent of Gaussian if we're using the simple approx. */
    if (w > wcutoff) {
	eb1 = 2.;
    }
    /* Calculate helper variables (should be moved later on...) */
    hsbw = s2 * h__ + eb1 * w2;
    hsbw2 = hsbw * hsbw;
    hsbw3 = hsbw2 * hsbw;
    hsbw4 = hsbw3 * hsbw;
    hsbw12 = sqrt(hsbw);
    hsbw32 = hsbw12 * hsbw;
    hsbw52 = hsbw32 * hsbw;
    hsbw72 = hsbw52 * hsbw;
    d1shsbw = d1sh * s2 + two * s * h__;
    d1rhsbw = two * eb1 * d1rw * w;
    dhsbw = d__ + s2 * h__ + eb1 * w2;
    dhsbw2 = dhsbw * dhsbw;
    dhsbw3 = dhsbw2 * dhsbw;
    dhsbw4 = dhsbw3 * dhsbw;
    dhsbw5 = dhsbw4 * dhsbw;
    dhsbw12 = sqrt(dhsbw);
    dhsbw32 = dhsbw12 * dhsbw;
    dhsbw52 = dhsbw32 * dhsbw;
    dhsbw72 = dhsbw52 * dhsbw;
    dhsbw92 = dhsbw72 * dhsbw;
    hsbwa94 = f94 * hsbw / a;
    hsbwa942 = hsbwa94 * hsbwa94;
    hsbwa943 = hsbwa942 * hsbwa94;
    hsbwa945 = hsbwa943 * hsbwa942;
    hsbwa9412 = sqrt(hsbwa94);
    dhs = d__ + s2 * h__;
    dhs2 = dhs * dhs;
    dhs3 = dhs2 * dhs;
    dhs4 = dhs3 * dhs;
    dhs72 = dhs3 * sqrt(dhs);
    dhs92 = dhs72 * dhs;
    d1sdhs = two * s * h__ + s2 * d1sh;
    dhsw = dhs + w2;
    dhsw2 = dhsw * dhsw;
    dhsw52 = sqrt(dhsw) * dhsw2;
    dhsw72 = dhsw52 * dhsw;
    d1rdhsw = two * d1rw * w;
    if (s > egscut) {
	d__1 = f32 * s * sqrt(h__ / a);
	g_a__ = srpi * (fifteen * e + six * c__ * (one + f * s2) * dhs + four 
		* b * dhs2 + eight * a * dhs3) * (one / (sixteen * dhs72)) - 
		f34 * pi * sqrt(a) * exp(f94 * h__ * s2 / a) * (one - erf(
		d__1));
	d__1 = f32 * sqrt(h__ / a) * s;
	d1sg_a__ = one / r32 * srpi * (r36 * (two * h__ + d1sh * s) / (a12 * 
		sqrt(h__ / a)) + one / dhs92 * (-eight * a * d1sdhs * dhs3 - 
		r105 * d1sdhs * e - r30 * c__ * d1sdhs * dhs * (one + s2 * f) 
		+ r12 * dhs2 * (-b * d1sdhs + c__ * s * (d1sf * s + two * f)
		)) - r54 * exp(f94 * h__ * s2 / a) * srpi * s * (two * h__ + 
		d1sh * s) * erf(d__1) / a12);
	g_b__ = f1516 * srpi * s2 / dhs72;
	d1sg_b__ = fifteen * srpi * s * (four * dhs - seven * d1sdhs * s) / 
		(r32 * dhs92);
	eg = -(f34 * pi + g_a__) / g_b__;
	d1seg = (-four * d1sg_a__ * g_b__ + d1sg_b__ * (four * g_a__ + three *
		 pi)) / (four * g_b__ * g_b__);
    } else {
	eg = ega1 + ega2 * s2 + ega3 * s4;
	d1seg = two * ega2 * s + four * ega3 * s3;
    }
    /* Calculate the terms needed in any case */
    term2 = (dhs2 * b + dhs * c__ + two * e + dhs * s2 * c__ * f + two * s2 * 
	    eg) / (two * dhs3);
    d1sterm2 = (-six * d1sdhs * (eg * s2 + e) + dhs2 * (-d1sdhs * b + s * 
	    c__ * (d1sf * s + two * f)) + two * dhs * (two * eg * s - 
	    d1sdhs * c__ + s2 * (d1seg - d1sdhs * c__ * f))) / (two * dhs4);
    term3 = -w * (four * dhsw2 * b + six * dhsw * c__ + fifteen * e + six * 
	    dhsw * s2 * c__ * f + fifteen * s2 * eg) / (eight * dhs * dhsw52);
    d1sterm3 = w * (two * d1sdhs * dhsw * (four * dhsw2 * b + six * dhsw * 
	    c__ + fifteen * e + three * s2 * (five * eg + two * dhsw * c__ * 
	    f)) + dhs * (r75 * d1sdhs * (eg * s2 + e) + four * dhsw2 * (
	    d1sdhs * b - three * s * c__ * (d1sf * s + two * f)) - six * 
	    dhsw * (-three * d1sdhs * c__ + s * (ten * eg + five * d1seg * 
	    s - three * d1sdhs * s * c__ * f)))) / (sixteen * dhs2 * dhsw72);
    d1rterm3 = (-two * d1rw * dhsw * (four * dhsw2 * b + six * dhsw * c__ + 
	    fifteen * e + three * s2 * (five * eg + two * dhsw * c__ * f)) + 
	    w * d1rdhsw * (r75 * (eg * s2 + e) + two * dhsw * (two * dhsw * b 
	    + nine * c__ + nine * s2 * c__ * f))) / (sixteen * dhs * dhsw72);
    term4 = -w3 * (dhsw * c__ + five * e + dhsw * s2 * c__ * f + five * s2 * 
	    eg) / (two * dhs2 * dhsw52);
    d1sterm4 = w3 * (four * d1sdhs * dhsw * (dhsw * c__ + five * e + s2 * (
	    five * eg + dhsw * c__ * f)) + dhs * (r25 * d1sdhs * (eg * s2 + e)
	     - two * dhsw2 * s * c__ * (d1sf * s + two * f) + dhsw * (three 
	    * d1sdhs * c__ + s * (-r20 * eg - ten * d1seg * s + three * 
	    d1sdhs * s * c__ * f)))) / (four * dhs3 * dhsw72);
    d1rterm4 = w2 * (-six * d1rw * dhsw * (dhsw * c__ + five * e + s2 * (five 
	    * eg + dhsw * c__ * f)) + w * d1rdhsw * (r25 * (eg * s2 + e) + 
	    three * dhsw * c__ * (one + s2 * f))) / (four * dhs2 * dhsw72);
    term5 = -w5 * (e + s2 * eg) / (dhs3 * dhsw52);
    d1sterm5 = w5 * (six * d1sdhs * dhsw * (eg * s2 + e) + dhs * (-two * dhsw 
	    * s * (two * eg + d1seg * s) + five * d1sdhs * (eg * s2 + e))) /
	     (two * dhs4 * dhsw72);
    d1rterm5 = w4 * five * (eg * s2 + e) * (-two * d1rw * dhsw + d1rdhsw * w) 
	    / (two * dhs3 * dhsw72);
    if (s > 0. || w > 0.) {
	t10 = f12 * a * log(hsbw / dhsbw);
	t10d1 = f12 * a * (one / hsbw - one / dhsbw);
	d1st10 = d1shsbw * t10d1;
	d1rt10 = d1rhsbw * t10d1;
    }
    /* Calculate exp(x)*f(x) depending on size of x */
    if (hsbwa94 < expfcutoff) {
	piexperf = pi * exp(hsbwa94) * erfc(hsbwa9412);
	d__1 = -hsbwa94;
	expei = exp(hsbwa94) * Ei(d__1);
    } else {
	piexperf = pi * (one / (srpi * hsbwa9412) - one / (two * sqrt(pi * 
		hsbwa943)) + three / (four * sqrt(pi * hsbwa945)));
	expei = -(one / hsbwa94) * (hsbwa942 + expei1 * hsbwa94 + expei2) / (
		hsbwa942 + expei3 * hsbwa94 + expei4);
    }

    /* Calculate the derivatives (based on the orig. expression) */
    /* --> Is this ok? ==> seems to be ok... */
    piexperfd1 = -(three * srpi * sqrt(hsbw / a)) / (two * hsbw) + nine * 
	    piexperf / (four * a);
    d1spiexperf = d1shsbw * piexperfd1;
    d1rpiexperf = d1rhsbw * piexperfd1;
    expeid1 = f14 * (four / hsbw + nine * expei / a);
    d1sexpei = d1shsbw * expeid1;
    d1rexpei = d1rhsbw * expeid1;
    if (w == zero) {
        /* Fall back to original expression for the PBE hole */
	t1 = -f12 * a * expei;
	d1st1 = -f12 * a * d1sexpei;
	d1rt1 = -f12 * a * d1rexpei;
	if (s > 0.) {
	    term1 = t1 + t10;
	    d1sterm1 = d1st1 + d1st10;
	    d1rterm1 = d1rt1 + d1rt10;
	    *fx = x * (term1 + term2);
	    *d1sfx = x * (d1sterm1 + d1sterm2);
	    *d1rfx = x * d1rterm1;
	} else {
	    *fx = 1.;
            /* TODO    This is checked to be true for term1 */
            /*         How about the other terms??? */
	    *d1sfx = 0.;
	    *d1rfx = 0.;
	}
    } else if (w > wcutoff) {
        /* Use simple Gaussian approximation for large w */
        /* print *,rho,s," LARGE w" */
	term1 = -f12 * a * (expei + log(dhsbw) - log(hsbw));
	term1d1 = -a / (two * dhsbw) - f98 * expei;
	d1sterm1 = d1shsbw * term1d1;
	d1rterm1 = d1rhsbw * term1d1;
	*fx = x * (term1 + term2 + term3 + term4 + term5);
	*d1sfx = x * (d1sterm1 + d1sterm2 + d1sterm3 + d1sterm4 + d1sterm5);
	*d1rfx = x * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5);
    } else {
        /* For everything else, use the full blown expression */
        /* First, we calculate the polynomials for the first term */
	np1 = -f32 * ea1 * a12 * w + r27 * ea3 * w3 / (eight * a12) - r243 * 
		ea5 * w5 / (r32 * a32) + r2187 * ea7 * w7 / (r128 * a52);
	d1rnp1 = -f32 * ea1 * d1rw * a12 + r81 * ea3 * d1rw * w2 / (eight * 
		a12) - r1215 * ea5 * d1rw * w4 / (r32 * a32) + r15309 * ea7 * 
		d1rw * w6 / (r128 * a52);
	np2 = -a + f94 * ea2 * w2 - r81 * ea4 * w4 / (sixteen * a) + r729 * 
		ea6 * w6 / (r64 * a2) - r6561 * ea8 * w8 / (r256 * a3);
	d1rnp2 = f12 * (nine * ea2 * d1rw * w) - r81 * ea4 * d1rw * w3 / (
		four * a) + r2187 * ea6 * d1rw * w5 / (r32 * a2) - r6561 * 
		ea8 * d1rw * w7 / (r32 * a3);
        /* The first term is */
	t1 = f12 * (np1 * piexperf + np2 * expei);
	d1st1 = f12 * (d1spiexperf * np1 + d1sexpei * np2);
	d1rt1 = f12 * (d1rnp2 * expei + d1rpiexperf * np1 + d1rexpei * np2 + 
		d1rnp1 * piexperf);
        /* The factors for the main polynomial in w and their derivatives */
	f2 = f12 * ea1 * srpi * a / dhsbw12;
	f2d1 = -ea1 * srpi * a / (four * dhsbw32);
	d1sf2 = d1shsbw * f2d1;
	d1rf2 = d1rhsbw * f2d1;
	f3 = f12 * ea2 * a / dhsbw;
	f3d1 = -ea2 * a / (two * dhsbw2);
	d1sf3 = d1shsbw * f3d1;
	d1rf3 = d1rhsbw * f3d1;
	f4 = ea3 * srpi * (-f98 / hsbw12 + f14 * a / dhsbw32);
	f4d1 = ea3 * srpi * (nine / (sixteen * hsbw32) - three * a / (eight * 
		dhsbw52));
	d1sf4 = d1shsbw * f4d1;
	d1rf4 = d1rhsbw * f4d1;
	f5 = ea4 * (one / r128) * (-r144 * (one / hsbw) + r64 * (one / dhsbw2)
		 * a);
	f5d1 = ea4 * (f98 / hsbw2 - a / dhsbw3);
	d1sf5 = d1shsbw * f5d1;
	d1rf5 = d1rhsbw * f5d1;
	f6 = ea5 * (three * srpi * (three * dhsbw52 * (nine * hsbw - two * a) 
		+ four * hsbw32 * a2)) / (r32 * dhsbw52 * hsbw32 * a);
	f6d1 = ea5 * srpi * (r27 / (r32 * hsbw52) - r81 / (r64 * hsbw32 * a) 
		- fifteen * a / (sixteen * dhsbw72));
	d1sf6 = d1shsbw * f6d1;
	d1rf6 = d1rhsbw * f6d1;
	f7 = ea6 * (r32 * a / dhsbw3 + (-r36 + r81 * s2 * h__ / a) / hsbw2) / 
		r32;
	d1sf7 = ea6 * (three * (r27 * d1sh * dhsbw4 * hsbw * s2 + eight * 
		d1shsbw * a * (three * dhsbw4 - four * hsbw3 * a) + r54 * 
		dhsbw4 * s * (hsbw - d1shsbw * s) * h__)) / (r32 * dhsbw4 * 
		hsbw3 * a);
	d1rf7 = ea6 * d1rhsbw * (f94 / hsbw3 - three * a / dhsbw4 - r81 * s2 *
		 h__ / (sixteen * hsbw3 * a));
	f8 = ea7 * (-three * srpi * (-r40 * hsbw52 * a3 + nine * dhsbw72 * (
		r27 * hsbw2 - six * hsbw * a + four * a2))) / (r128 * dhsbw72 
		* hsbw52 * a2);
	f8d1 = ea7 * srpi * (r135 / (r64 * hsbw72) + r729 / (r256 * hsbw32 * 
		a2) - r243 / (r128 * hsbw52 * a) - r105 * a / (r32 * dhsbw92))
		;
	d1sf8 = d1shsbw * f8d1;
	d1rf8 = d1rhsbw * f8d1;
	f9 = (r324 * ea6 * eb1 * dhsbw4 * hsbw * a + ea8 * (r384 * hsbw3 * a3 
		+ dhsbw4 * (-r729 * hsbw2 + r324 * hsbw * a - r288 * a2))) / (
		r128 * dhsbw4 * hsbw3 * a2);
	f9d1 = -(r81 * ea6 * eb1 / (sixteen * hsbw3 * a)) + ea8 * (r27 / (
		four * hsbw4) + r729 / (r128 * hsbw2 * a2) - r81 / (sixteen * 
		hsbw3 * a) - r12 * a / dhsbw5);
	d1sf9 = d1shsbw * f9d1;
	d1rf9 = d1rhsbw * f9d1;
	t2t9 = f2 * w + f3 * w2 + f4 * w3 + f5 * w4 + f6 * w5 + f7 * w6 + f8 *
		 w7 + f9 * w8;
	d1st2t9 = d1sf2 * w + d1sf3 * w2 + d1sf4 * w3 + d1sf5 * w4 + d1sf6 * 
		w5 + d1sf7 * w6 + d1sf8 * w7 + d1sf9 * w8;
	d1rt2t9 = d1rw * f2 + d1rf2 * w + two * d1rw * f3 * w + d1rf3 * w2 + 
		three * d1rw * f4 * w2 + d1rf4 * w3 + four * d1rw * f5 * w3 + 
		d1rf5 * w4 + five * d1rw * f6 * w4 + d1rf6 * w5 + six * d1rw *
		 f7 * w5 + d1rf7 * w6 + seven * d1rw * f8 * w6 + d1rf8 * w7 + 
		eight * d1rw * f9 * w7 + d1rf9 * w8;
        /*       The final value of term1 for 0 < omega < wcutoff is: */
	term1 = t1 + t2t9 + t10;
	d1sterm1 = d1st1 + d1st2t9 + d1st10;
	d1rterm1 = d1rt1 + d1rt2t9 + d1rt10;
        /*       The final value for the enhancement factor and its */
        /*       derivatives is: */
	*fx = x * (term1 + term2 + term3 + term4 + term5);
	*d1sfx = x * (d1sterm1 + d1sterm2 + d1sterm3 + d1sterm4 + d1sterm5);
	*d1rfx = x * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5);
    }
}


}
