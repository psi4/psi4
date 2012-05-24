#include <libmints/vector.h>
#include "xfunctional.h"
#include "utility.h"
#include <psi4-dec.h>
#include <cmath>

using namespace psi;

namespace psi {

XFunctional::XFunctional()
{
    common_init();
}
XFunctional::~XFunctional()
{
}
void XFunctional::common_init()
{
    gga_type_ = GGA_None;
    meta_type_ = Meta_None;
    sr_type_ = SR_None;
 
    _K0_ = 3.0 * pow(3.0 / (4.0 * M_PI), 1.0/3.0);
    _C0_ = 3.0 / 5.0 * pow(6.0 * M_PI * M_PI, 2.0/3.0);   
    _k0_ = pow(6.0 * M_PI * M_PI, 1.0/3.0);
    _pi12_ = sqrt(M_PI);
    
    _B88_d_ = 0.0042;
    _B88_a_ = 1.0;

    _PBE_kp_ = 0.804;
    _PBE_mu_ = 0.2195149727645171;

    _PW91_a1_ = 0.19645 / (2.0 * _k0_);
    _PW91_a2_ = 7.7956 / (2.0 * _k0_);
    _PW91_a3_ = 0.2743 / (4.0 * _k0_ * _k0_);
    _PW91_a4_ = 0.1508 / (4.0 * _k0_ * _k0_);
    _PW91_a5_ = 100.0 / (4.0 * _k0_ * _k0_);
    _PW91_a6_ = 0.004 / (16.0 * _k0_ * _k0_ * _k0_ * _k0_);

    _B97_gamma_ = 0.004;
}
void XFunctional::set_parameter(const std::string& key, double val) 
{
    parameters_[key] = val;

    if (key == "B88_d") {
        _B88_d_ = val;
    } else if (key == "B88_a") {
        _B88_a_ = val;
    } else if (key == "PW91_a1") {
        _PW91_a1_ = val;
    } else if (key == "PW91_a2") {
        _PW91_a2_ = val;
    } else if (key == "PW91_a3") {
        _PW91_a3_ = val;
    } else if (key == "PW91_a4") {
        _PW91_a4_ = val;
    } else if (key == "PW91_a5") {
        _PW91_a5_ = val;
    } else if (key == "PW91_a6") {
        _PW91_a6_ = val;
    } else if (key == "B97_gamma") {
        _B97_gamma_ = val;
    } else if (key == "PBE_kp") {
        _PBE_kp_ = val;
    } else if (key == "PBE_mu") {
        _PBE_mu_ = val;
    } else if (key.substr(0,5) == "B97_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(5).c_str());
        if (_B97_a_.size() < index + 1) {
            _B97_a_.resize(index + 1);
            _B97_a_[index] = val;
        }
    } else if (key.substr(0,6) == "Meta_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(6).c_str());
        if (_Meta_a_.size() < index + 1) {
            _Meta_a_.resize(index + 1);
            _Meta_a_[index] = val;
        }
    } else {
        throw PSIEXCEPTION("Error, unknown generalized exchange functional parameter");    
    }
}
void XFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    compute_sigma_functional(in,out,npoints,deriv,alpha,true);
    compute_sigma_functional(in,out,npoints,deriv,alpha,false);
}
void XFunctional::compute_sigma_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("XFunctional: 2nd and higher partials not implemented yet.");
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

        // Powers of rho
        double rho13 = pow(rho,1.0/3.0);
        double rho23 = rho13 * rho13;
        double rho43 = rho * rho13;
        double rho53 = rho * rho23;
        double rho73 = rho * rho * rho13;

        //=>  Factors <= //

        // > LSDA < //
        double E, E_rho;
        E = - 0.5 * _K0_ * rho43;
        E_rho = -4.0/6.0 * _K0_ * rho13;
    
        // > GGA < //

        // s
        double s, s_rho, s_gamma, s_tau;
        if (gga_) {
            s = sqrt(gamma) / rho43;
            s_rho = - 4.0 / 3.0 * sqrt(gamma) / rho73;
            s_gamma = 1.0 / 2.0 * pow(gamma,-1.0/2.0) / rho43; 
            s_tau = 0.0;
        } else {
            s = 0.0;
            s_rho = 0.0;
            s_gamma = 0.0;
            s_tau = 0.0;
        }

        // Fs
        double Fs, Fs_s;
        switch (gga_type_) {
            case GGA_None: {
                Fs = 1.0;
                Fs_s = 0.0;
                break;
            }
            case B88: {
                double s2p1 = s * s + 1.0;
                double s2p1_12 = sqrt(s2p1);
                double asinhs = log(s + s2p1_12);

                double N = 2.0 / _K0_ * _B88_a_ * _B88_d_ * s * s;
                double D = 1.0 + 6.0 * _B88_d_ * s * asinhs;

                double N_s = 4.0 / _K0_ * _B88_a_ * _B88_d_ * s;
                double D_s = 6.0 * _B88_d_ * asinhs + 6.0 * _B88_d_ * s / s2p1_12;
                
                Fs = 1.0 + N / D;
                Fs_s = (N_s * D - D_s * N) / (D * D);
                break;
            }
            case PBE: {
                double mus2 = 1.0 + _PBE_mu_ * s * s / (4.0 * _k0_ * _k0_ * _PBE_kp_);
                Fs = 1.0 + _PBE_kp_ - _PBE_kp_ / mus2; 
                Fs_s = 2.0 / (mus2 * mus2) * _PBE_mu_ * s / (4.0 * _k0_ * _k0_); 
                break;
            }
            case RPBE: {
                double mus2 = 1.0 + _PBE_mu_ * s * s / (4.0 * _k0_ * _k0_ * _PBE_kp_);
                Fs = 1.0 + _PBE_kp_ * (1.0 - exp(- _PBE_mu_ * s * s / _PBE_kp_); 
                Fs_s = 2.0 * _PBE_mu_ * s * exp(- _PBE_mu_ * s * s / _PBE_kp_); 
                break;
            }
            case PW91: {
                double bs = _PW91_a2_ * s;
                double bs2 = bs * bs;
                double bs2p1 = 1.0 + bs2;        
                double bs2p1_12 = sqrt(bs2p1);
                double aasinhv = _PW91_a1_ * log(bs + bs2p1_12);   
                double dexpv = _PW91_a4_ * exp(-_PW91_a5_ * s * s);

                double N = 1.0 + aasinhv * s + (_PW91_a3_ - dexpv) * s * s;
                double D = 1.0 + aasinhv * s  + _PW91_a6_ * s * s * s * s;
                
                double N_s =  aasinhv + _PW91_a1_ * _PW91_a2_ * s / bs2p1_12 + 2.0 * s * (_PW91_a3_ - dexpv) + 2.0 * s * s * s * _PW91_a5_ * dexpv;
                double D_s =  aasinhv + _PW91_a1_ * _PW91_a2_ * s / bs2p1_12 + 4.0 * _PW91_a6_ * s * s * s; 
            
                Fs = N / D;
                Fs_s = (N_s * D - D_s * N) / (D * D);
                break;
            }
            case B97: {
                int size = _B97_a_.size();
                Fs = 0.0;
                Fs_s = 0.0;
                double s2 = s * s;
                double gs2 = 1.0 + _B97_gamma_ * s2; 
                double g = _B97_gamma_ * s2 / gs2; 
                double g_s = 2.0 * _B97_gamma_ * s / (gs2 * gs2);
                
                double buf = 1.0; 
                double buf2 = 0.0; 

                for (int A = 0; A < size; A++) {
                    double a = _B97_a_[A];
                    Fs += a * buf;
                    Fs_s += A * a * buf2 * g_s; 
                    buf2 = buf;
                    buf *= g;                  
                }
                break;
            }
        }
        
        // > Meta < //

        // w
        double w, w_rho, w_gamma, w_tau;
        if (meta_) {
            double tauLSDA = _C0_ * rho53;
            double denom = tauLSDA + tau;
            w = (tauLSDA - tau) / (denom);
            w_rho = 10.0/3.0 * _C0_ * rho23 * tau / (denom * denom);
            w_gamma = 0.0;
            w_tau = -2.0 * tauLSDA / (denom * denom);        
        } else {
            w = 0.0;
            w_rho = 0.0;
            w_gamma = 0.0;
            w_tau = 0.0;
        }

        // Fw
        double Fw, Fw_w;
        switch (meta_type_) {
            case Meta_None: {
                Fw = 1.0;   
                Fw_w = 0.0;
                break;
            }
            default: {
                int size = _Meta_a_.size();
                Fw = 0.0;
                Fw_w = 0.0;
            
                double buf = 1.0; 
                double buf2 = 0.0; 

                for (int A = 0; A < size; A++) {
                    double a = _Meta_a_[A];
                    Fw += a * buf;
                    Fw_w += A * a * buf2; 
                    buf2 = buf;
                    buf *= w;                  
                }
                break;
            }
        }

        // > SR < //

        // k
        double k, k_rho, k_gamma, k_tau;
        switch (sr_type_) {
            case SR_None: {
                k = 0.0;
                k_rho = 0.0;
                k_gamma = 0.0;
                k_tau = 0.0;
                break;
            }
            case LSDA: {
                k = _k0_ * rho13;
                k_rho = 1.0/3.0 * _k0_ / rho23;
                k_gamma = 0.0; 
                k_tau = 0.0;
                break;
            }
            case GGA: {
                double Fs12 = sqrt(Fs);
                k = _k0_ * rho13 / Fs12;
                k_rho = 1.0/3.0 * _k0_ / (Fs12 * rho23) 
                    - 1.0 / 2.0 * k / (Fs * Fs12) * Fs_s * s_rho;
                k_gamma = - 1.0 / 2.0 * k / (Fs * Fs12) * Fs_s * s_gamma;
                k_tau = 0.0;
                break;
            }
            case Meta: { 
                double Fs12 = sqrt(Fs);
                double Fw12 = sqrt(Fw);
                k = _k0_ * rho13 / (Fs12 * Fw12);
                k_rho = 1.0/3.0 * _k0_ / (Fs12 * Fw12 *rho23) 
                    - 1.0 / 2.0 * k / (Fs * Fs12 * Fw12) * Fs_s * s_rho;
                    - 1.0 / 2.0 * k / (Fw * Fs12 * Fw12) * Fw_w * w_rho;
                k_gamma = - 1.0 / 2.0 * k / (Fs * Fs12 * Fw12) * Fs_s * s_gamma;
                k_tau   = - 1.0 / 2.0 * k / (Fw * Fs12 * Fw12) * Fw_w * w_tau;
                break;
            }
        }
 
        // Fk
        double Fk, Fk_k;
        switch (sr_type_) {
            case SR_None: {
                Fk = 1.0;
                Fk_k = 0.0;
                break;
            } 
            default: {    
                double a = omega_ / (2.0 * k);
                double a_k = - omega_ / (2.0 * k * k); 

                if (a < 1.0E-10) {
                    Fk = 1.0;
                    Fk_k = 0.0;
                    break;
                }

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
                
                Fk = F1k + F2k;
                Fk_k = (F1k_a + F2k_a) * a_k;
                break;
            }
        }
    
        // => Assembly <= //
        v[Q] += A * E * Fs * Fw * Fk; 
        if (deriv >= 1) {
            v_rho[Q] += A * (Fs * Fw * Fk * (E_rho) +
                             E  * Fw * Fk * (Fs_s * s_rho) +
                             E  * Fs * Fk * (Fw_w * w_rho) + 
                             E  * Fs * Fw * (Fk_k * k_rho));
            if (gga_) {
                v_gamma[Q] += A * (E  * Fw * Fk * (Fs_s * s_gamma) +
                                   E  * Fs * Fk * (Fw_w * w_gamma) + 
                                   E  * Fs * Fw * (Fk_k * k_gamma));
            }
            if (meta_) {
                v_tau[Q] += A * (E  * Fw * Fk * (Fs_s * s_tau) +
                                 E  * Fs * Fk * (Fw_w * w_tau) + 
                                 E  * Fs * Fw * (Fk_k * k_tau));
            }
        }
    }
}

}
