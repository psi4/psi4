#include <libmints/vector.h>
#include "cfunctional.h"
#include "utility.h"
#include <psi4-dec.h>
#include <cmath>

using namespace psi;

namespace psi {

CFunctional::CFunctional()
{
    common_init();
}
CFunctional::~CFunctional()
{
}
void CFunctional::common_init()
{
    gga_type_ = GGA_None;
    meta_type_ = Meta_None;
    lsda_type_ = LSDA_None;

    // PW92
    _c0_ = 1.0/4.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(M_PI,-1.0/3.0);
    _two13_ = pow(2.0, 1.0/3.0);
    _d2fz0_ = 1.709921;

    _c0a_ = 0.016887;
    _a1a_ = 0.11125;
    _b1a_ = 10.357;
    _b2a_ = 3.6231;
    _b3a_ = 0.88026;
    _b4a_ = 0.49671;

    _c0p_ = 0.031091;
    _a1p_ = 0.21370;
    _b1p_ = 7.5957;
    _b2p_ = 3.5876;
    _b3p_ = 1.6382;
    _b4p_ = 0.49294;
   
    _c0f_ = 0.015545;
    _a1f_ = 0.20548; 
    _b1f_ = 14.1189; 
    _b2f_ = 6.1977; 
    _b3f_ = 3.3662;
    _b4f_ = 0.62517;
 
    _B97_ss_gamma_ = 0.0380;
    _B97_os_gamma_ = 0.0031;
}
void CFunctional::set_parameter(const std::string& key, double val) 
{
    parameters_[key] = val;

    if (key == "B97_ss_gamma") {
        _B97_ss_gamma_ = val;
    } else if (key.substr(0,8) == "B97_ss_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(8).c_str());
        if (_B97_ss_a_.size() < index + 1) {
            _B97_ss_a_.resize(index + 1);
            _B97_ss_a_[index] = val;
        }
    } else if (key == "B97_os_gamma") {
        _B97_os_gamma_ = val;
    } else if (key.substr(0,8) == "B97_os_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(8).c_str());
        if (_B97_os_a_.size() < index + 1) {
            _B97_os_a_.resize(index + 1);
            _B97_os_a_[index] = val;
        }
    } else {
        throw PSIEXCEPTION("Error, unknown generalized correlation functional parameter");    
    }
}
void CFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    compute_ss_functional(in,out,npoints,deriv,alpha,true);
    compute_ss_functional(in,out,npoints,deriv,alpha,false);
    compute_os_functional(in,out,npoints,deriv,alpha);
}
void CFunctional::compute_ss_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("CFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

    // => Input variables (spin-polarized) <= //

    double* rho_s = NULL;
    double* gamma_s = NULL;
    double* tau_s = NULL;
    rho_s = in.find(spin ? "RHO_A" : "RHO_B")->second->pointer();
    if (gga_) {
        gamma_s = in.find(spin ? "GAMMA_AA" : "GAMMA_BB")->second->pointer();
    }
    if (meta_) {
        tau_s = in.find(spin ? "TAU_A" : "TAU_B")->second->pointer();
    }

    // => Output variables <= //

    double* v = NULL;
    double* v_rho = NULL;
    double* v_gamma = NULL;
    double* v_tau = NULL;
    
    v = out.find(spin ? "V" : "V")->second->pointer();
    if (deriv >= 1) {
        v_rho = out.find(spin ? "V_RHO_A" : "V_RHO_B")->second->pointer();
        if (gga_) {
            v_gamma = out.find(spin ? "V_GAMMA_AA" : "V_GAMMA_BB")->second->pointer();
        }
        if (meta_) {
            v_tau = out.find(spin ? "V_TAU_A" : "V_TAU_B")->second->pointer();
        }
    }
     
    // => Main Loop over points <= //
    for (int Q = 0; Q < npoints; Q++) {
    
        // => Primitive variables <= //
        double rho;
        double gamma;
        double tau;

        rho = rho_s[Q];
        if (rho < lsda_cutoff_) {
            continue;
        } 
        if (gga_) {
            gamma = gamma_s[Q];
        }
        if (meta_) {
            tau = tau_s[Q];
        }

        // => Powers of rho <= //
        double rho13 = pow(rho, 1.0/3.0);
        double rho43 = rho * rho13;
        double rho73 = rho * rho43;

        // => LDA Part <= //
        
        double E;
        double E_rho;
        if (lsda_type_ == LSDA_None) {
            throw PSIEXCEPTION("WTF");
        } else if (lsda_type_ == PW92) {
            
            double r = _c0_ * pow(rho, -1.0/3.0);
            double r_rho = -1.0/3.0 * _c0_ * pow(rho, -4.0/3.0);
            
            double EcP;
            double EcP_r;
            G2(r, &EcP, &EcP_r);
           
            E = rho * EcP;
            E_rho = EcP_r * r_rho; 
        }

        // => GGA Part <= //
        double s, s_rho, s_gamma;
        if (gga_) {
            s = sqrt(gamma) / rho43;
            s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
            s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) / rho43; 
        } else {
            s = 0.0;
            s_rho = 0.0;
            s_gamma = 0.0;
        }

        double Fs;
        double Fs_s;
        if (gga_type_ == GGA_None) {
            Fs = 1.0;
            Fs_s = 0.0;
        } else if (gga_type_ == B97) {
            Fs = 0.0;
            Fs_s = 0.0;
            double s2 = s * s;
            double gs2 = 1.0 + _B97_ss_gamma_ * s2; 
            double g = _B97_ss_gamma_ * s2 / gs2; 
            double g_s = 2.0 * _B97_ss_gamma_ * s / (gs2 * gs2);

            double buf = 1.0; 
            double buf2 = 0.0; 

            int size = _B97_ss_a_.size();
            for (int A = 0; A < size; A++) {
                double a = _B97_ss_a_[A];
                Fs += a * buf;
                Fs_s += A * a * buf2 * g_s; 
                buf2 = buf;
                buf *= g;                  
            }
        }

        // => Meta <= //
        double D, D_rho, D_gamma, D_tau;
        if (meta_type_ == Meta_None) {
            D = 1.0;
            D_gamma = 0.0;
            D_rho = 0.0;
            D_tau = 0.0;
        } else if (meta_type_ == B95) {
            D = 1.0 - gamma / (8.0 * tau * rho);
            D_gamma = - 1.0 / (8.0 * tau * rho * rho);
            D_rho = gamma / (8.0 * tau * rho * rho);
            D_tau = gamma / (8.0 * tau * tau * rho);
        }
 
        // => Assembly <= //
        v[Q] += A * E * Fs * D; 
        if (deriv >= 1) {
            v_rho[Q] += A * (Fs * D  * (E_rho) +
                             E  * D  * (Fs_s * s_rho) +
                             E  * Fs * (D_rho));
            if (gga_) {
                v_gamma[Q] += A * (E  * D  * (Fs_s * s_gamma) +
                                   E  * Fs * (D_gamma));
            }
            if (meta_) {
                v_tau[Q] += A * (E * Fs * D_tau);
            }
        }
    }
}
void CFunctional::compute_os_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("CFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

    // => Input variables (ab-type, no meta, no cross terms) <= //

    double* rho_ap = NULL;
    double* rho_bp = NULL;
    double* gamma_aap = NULL;
    double* gamma_bbp = NULL;

    rho_ap = in.find("RHO_A")->second->pointer();
    rho_bp = in.find("RHO_B")->second->pointer();
    if (gga_) {
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        gamma_bbp = in.find("GAMMA_BB")->second->pointer();
    }

    // => Output variables <= //

    double* v = NULL;
    double* v_rho_a = NULL;
    double* v_rho_b = NULL;
    double* v_gamma_aa = NULL;
    double* v_gamma_bb = NULL;
    
    v = out.find("V")->second->pointer();
    if (deriv >= 1) {
        v_rho_a = out.find("V_RHO_A")->second->pointer();
        v_rho_b = out.find("V_RHO_B")->second->pointer();
        if (gga_) {
            v_gamma_aa = out.find("V_GAMMA_AA")->second->pointer();
            v_gamma_bb = out.find("V_GAMMA_BB")->second->pointer();
        }
    }
     
    // => Main Loop over points <= //
    for (int Q = 0; Q < npoints; Q++) {
    
        // => Primitive variables <= //
        double rho_a;
        double rho_b;
        double gamma_aa;
        double gamma_bb;

        rho_a = rho_ap[Q];
        rho_b = rho_bp[Q];
        if (rho_a < lsda_cutoff_ || rho_b < lsda_cutoff_) {
            continue;
        } 
        if (gga_) {
            gamma_aa = gamma_aap[Q];
            gamma_bb = gamma_bbp[Q];
        }

        // => Powers of rho <= //
        double rho_a13 = pow(rho_a, 1.0/3.0);
        double rho_a43 = rho_a * rho_a13;
        double rho_a73 = rho_a * rho_a43;
        double rho_a83 = rho_a73 * rho_a13;
        double rho_b13 = pow(rho_b, 1.0/3.0);
        double rho_b43 = rho_b * rho_b13;
        double rho_b73 = rho_b * rho_b43;
        double rho_b83 = rho_b73 * rho_a13;

        // => LDA Part <= //
        
        double E;
        double E_rho_a;
        double E_rho_b;
        if (lsda_type_ == LSDA_None) {
            throw PSIEXCEPTION("WTF");
        } else if (lsda_type_ == PW92) {

            double z = (rho_a - rho_b) / (rho_a + rho_b);

            double f = (pow(1.0 + z, 4.0/3.0) + (1.0 - z, 4.0/3.0) - 2.0) / (2.0 * (_two13_ - 1.0));    

            double r = _c0_ * pow(rho_a + rho_b, -1.0/3.0);
            
            double EcP;
            double EcP_r;

            double EcF;
            double EcF_r;

            double Ac;
            double Ac_r;

            G1(r, &Ac , &Ac_r );
            G2(r, &EcP, &EcP_r);
            G3(r, &EcF, &EcF_r);

            double Q = EcP + Ac * f * (1.0 - z * z * z * z) / _d2fz0_ + (EcF - EcP) * f * z * z * z * z;

            double ra = _c0_ * pow(rho_a, -1.0/3.0);
            double rb = _c0_ * pow(rho_b, -1.0/3.0);
            
            double EcPa;
            double EcPa_ra;

            double EcPb;
            double EcPb_rb;

            G2(ra, &EcPa, &EcPa_ra);
            G2(rb, &EcPb, &EcPb_rb);

            double Ea = rho_a * EcPa;
            double Eb = rho_b * EcPb;

            E = (rho_a + rho_b) * Q - Ea - Eb;

            double E_Ea = -1.0;
            double E_Eb = -1.0;

            double Ea_EcPa = rho_a;
            double Eb_EcPb = rho_b;

            double E_EcPa = Ea_EcPa * E_Ea; 
            double E_EcPb = Eb_EcPb * E_Eb; 

            double E_ra = EcPa_ra * E_EcPa; 
            double E_rb = EcPb_rb * E_EcPb; 

            double E_Q = (rho_a + rho_b);

            double Q_EcF = f * z * z * z * z;

            double E_EcF = Q_EcF * E_Q;

            double Q_EcP = 1.0 - f * z * z * z * z;

            double E_EcP = Q_EcP * E_Q;

            double Q_Ac = f * (1.0 - z * z * z * z) / _d2fz0_;

            double E_Ac = Q_Ac * E_Q;

            double E_r = EcP_r * E_EcP + EcF_r * E_EcF + Ac_r * E_Ac;

            double Q_f = Ac * (1.0 - z * z * z * z) / _d2fz0_ + (EcF - EcP) * z * z * z * z;

            double E_f = Q_f * E_Q;

            double Q_z = -4.0 * Ac * f * z * z * z / _d2fz0_ + 4.0 * (EcF - EcP) * f * z * z * z; 
            double f_z = (4.0/3.0 * pow(1.0 + z, 1.0/3.0) + 4.0/3.0 * pow(1.0 - z, 1.0/3.0)) / (2.0 * (_two13_ - 1.0));

            double E_z = Q_z * E_Q + f_z * Q_f;

            double r_rho_a = -1.0/3.0 * _c0_ * pow(rho_a + rho_b, -4.0/3.0);
            double r_rho_b = r_rho_a;

            double ra_rho_a = -1.0/3.0 * _c0_ * pow(rho_a, -4.0/3.0);
            double rb_rho_b = -1.0/3.0 * _c0_ * pow(rho_b, -4.0/3.0);

            double z_rho_a =  2.0 * rho_b / ((rho_a + rho_b) * (rho_a + rho_b));
            double z_rho_b = -2.0 * rho_a / ((rho_a + rho_b) * (rho_a + rho_b));

            E_rho_a = r_rho_a * E_r + z_rho_a * E_z + ra_rho_a * E_ra;
            E_rho_b = r_rho_b * E_r + z_rho_b * E_z + rb_rho_b * E_rb;
        }

        // => GGA Part <= //
        double s2, s2_rho_a, s2_gamma_aa, s2_rho_b, s2_gamma_bb;
        if (gga_) {
            s2 = gamma_aa / rho_a83 + gamma_bb / rho_b83;
            s2_rho_a = -8.0/3.0 * gamma_aa / (rho_a83 * rho_a); 
            s2_rho_b = -8.0/3.0 * gamma_bb / (rho_b83 * rho_b); 
            s2_gamma_aa = 1.0 / rho_a83; 
            s2_gamma_bb = 1.0 / rho_b83; 
        } else {
            s2 = 0.0;
            s2_rho_a = 0.0;
            s2_rho_b = 0.0;
            s2_gamma_aa = 0.0;
            s2_gamma_bb = 0.0;
        }

        double Fs2;
        double Fs2_s2;
        if (gga_type_ == GGA_None) {
            Fs2 = 1.0;
            Fs2_s2 = 0.0;
        } else if (gga_type_ == B97) {
            Fs2 = 0.0;
            Fs2_s2 = 0.0;
            double gs2 = 1.0 + _B97_os_gamma_ * s2; 
            double g = _B97_os_gamma_ * s2 / gs2; 
            double g_s2 = _B97_os_gamma_ / (gs2 * gs2);

            double buf = 1.0; 
            double buf2 = 0.0; 

            int size = _B97_os_a_.size();
            for (int A = 0; A < size; A++) {
                double a = _B97_os_a_[A];
                Fs2 += a * buf;
                Fs2_s2 += A * a * buf2 * g_s2; 
                buf2 = buf;
                buf *= g;                  
            }
        }

        // => Assembly <= //
        v[Q] += A * E * Fs2; 
        if (deriv >= 1) {
            v_rho_a[Q] += A * (Fs2 * E_rho_a +
                               E  * Fs2_s2 * s2_rho_a);
            v_rho_b[Q] += A * (Fs2 * E_rho_b +
                               E  * Fs2_s2 * s2_rho_b);
            if (gga_) {
                v_gamma_aa[Q] += A * (E  * Fs2_s2 * s2_gamma_aa);
                v_gamma_bb[Q] += A * (E  * Fs2_s2 * s2_gamma_bb);
            }
        }
    }

}
void CFunctional::G1(double r, double *G, double *G_r)
{
    double r12 = sqrt(r);
    double r32 = r * r12;
    double r2  = r * r;
    
    double A1 = _c0a_;    
    double a1 = _a1a_;
    double b1 = _b1a_;
    double b2 = _b2a_;
    double b3 = _b3a_;
    double b4 = _b4a_;

    double D = 2.0 * A1 * (b1 * r12 + b2 * r + b3 * r32 + b4 * r2);
    double D_r = 2.0 * A1 * (1.0/2.0 * b1 / r12 + b2 + 3.0/2.0 * b3 * r12 + 2.0 * b4 * r);
    
    double Q = 1.0 + 1.0 / D;
    double R = log(Q); 
    double L = -2.0 * A1 * (1.0 + a1 * r);

    // Minus: GOTCHA
    *G = - L * R; 

    double L_r = -2.0 * A1 * a1; 
    double R_r = - 1.0 / (Q * D * D) * D_r;

    // Minus: GOTCHA
    *G_r = - (L_r * R + L * R_r);
}
void CFunctional::G2(double r, double *G, double *G_r)
{
    double r12 = sqrt(r);
    double r32 = r * r12;
    double r2  = r * r;
    
    double A1 = _c0p_;    
    double a1 = _a1p_;
    double b1 = _b1p_;
    double b2 = _b2p_;
    double b3 = _b3p_;
    double b4 = _b4p_;

    double D = 2.0 * A1 * (b1 * r12 + b2 * r + b3 * r32 + b4 * r2);
    double D_r = 2.0 * A1 * (1.0/2.0 * b1 / r12 + b2 + 3.0/2.0 * b3 * r12 + 2.0 * b4 * r);
    
    double Q = 1.0 + 1.0 / D;
    double R = log(Q); 
    double L = -2.0 * A1 * (1.0 + a1 * r);

    *G = L * R; 

    double L_r = -2.0 * A1 * a1; 
    double R_r = - 1.0 / (Q * D * D) * D_r;

    *G_r = (L_r * R + L * R_r);
}
void CFunctional::G3(double r, double *G, double *G_r)
{
    double r12 = sqrt(r);
    double r32 = r * r12;
    double r2  = r * r;
    
    double A1 = _c0f_;    
    double a1 = _a1f_;
    double b1 = _b1f_;
    double b2 = _b2f_;
    double b3 = _b3f_;
    double b4 = _b4f_;

    double D = 2.0 * A1 * (b1 * r12 + b2 * r + b3 * r32 + b4 * r2);
    double D_r = 2.0 * A1 * (1.0/2.0 * b1 / r12 + b2 + 3.0/2.0 * b3 * r12 + 2.0 * b4 * r);
    
    double Q = 1.0 + 1.0 / D;
    double R = log(Q); 
    double L = -2.0 * A1 * (1.0 + a1 * r);

    *G = L * R; 

    double L_r = -2.0 * A1 * a1; 
    double R_r = - 1.0 / (Q * D * D) * D_r;

    *G_r = (L_r * R + L * R_r);
}

}
