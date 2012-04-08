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

    Truhlar_version_ = false;
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
    if (!Truhlar_version_) {
        compute_ss_functional(in,out,npoints,deriv,alpha,true);
        compute_ss_functional(in,out,npoints,deriv,alpha,false);
    }
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
        double rho83 = rho73 * rho13;

        // => LDA Part <= //
        
        double E;
        double E_rho;
        if (lsda_type_ == LSDA_None) {
            throw PSIEXCEPTION("WTF");
        } else if (lsda_type_ == PW92) {
            double E_z;
            PW92_C(rho, 1.0, &E, &E_rho, &E_z);    
        }

        // => GGA Part <= //
        double s2, s2_rho, s2_gamma;
        if (gga_) {
            s2 = gamma / rho83;
            s2_rho = - 8.0 / 3.0 * gamma / (rho83 * rho);
            s2_gamma = 1.0 / rho83; 
        } else {
            s2 = 0.0;
            s2_rho = 0.0;
            s2_gamma = 0.0;
        }

        double Fs2;
        double Fs2_s2;
        if (gga_type_ == GGA_None) {
            Fs2 = 1.0;
            Fs2_s2 = 0.0;
        } else if (gga_type_ == B97) {
            Fs2 = 0.0;
            Fs2_s2 = 0.0;
            double gs2 = 1.0 + _B97_ss_gamma_ * s2; 
            double g = _B97_ss_gamma_ * s2 / gs2; 
            double g_s2 = _B97_ss_gamma_ / (gs2 * gs2);

            double buf = 1.0; 
            double buf2 = 0.0; 

            int size = _B97_ss_a_.size();
            for (int A = 0; A < size; A++) {
                double a = _B97_ss_a_[A];
                Fs2 += a * buf;
                Fs2_s2 += A * a * buf2 * g_s2; 
                buf2 = buf;
                buf *= g;                  
            }
        }

        // => Meta <= //
        double D, D_rho, D_gamma, D_tau;
        if (meta_type_ == Meta_None) {
            D = 1.0;
            D_rho = 0.0;
            D_gamma = 0.0;
            D_tau = 0.0;
        } else if (meta_type_ == B95) {
            D = 1.0 - gamma / (4.0 * tau * rho);
            D_rho = gamma / (4.0 * tau * rho * rho);
            D_gamma = - 1.0 / (4.0 * tau * rho);
            D_tau = gamma / (4.0 * tau * tau * rho);
        }
 
        // => Assembly <= //
        v[Q] += A * E * Fs2 * D; 
        if (deriv >= 1) {
            v_rho[Q] += A * (Fs2 * D   * (E_rho) +
                             E   * D   * (Fs2_s2 * s2_rho) +
                             E   * Fs2 * (D_rho));
            if (gga_) {
                v_gamma[Q] += A * (E  * D   * (Fs2_s2 * s2_gamma) +
                                   E  * Fs2 * (D_gamma));
            }
            if (meta_) {
                v_tau[Q] += A * (E * Fs2 * D_tau);
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
    double* tau_ap = NULL;
    double* tau_bp = NULL;

    rho_ap = in.find("RHO_A")->second->pointer();
    rho_bp = in.find("RHO_B")->second->pointer();
    if (gga_) {
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        gamma_bbp = in.find("GAMMA_BB")->second->pointer();
    }
    if (meta_) {
        tau_ap = in.find("TAU_A")->second->pointer();
        tau_bp = in.find("TAU_B")->second->pointer();
    }

    // => Output variables <= //

    double* v = NULL;
    double* v_rho_a = NULL;
    double* v_rho_b = NULL;
    double* v_gamma_aa = NULL;
    double* v_gamma_bb = NULL;
    double* v_tau_a = NULL;
    double* v_tau_b = NULL;
    
    v = out.find("V")->second->pointer();
    if (deriv >= 1) {
        v_rho_a = out.find("V_RHO_A")->second->pointer();
        v_rho_b = out.find("V_RHO_B")->second->pointer();
        if (gga_) {
            v_gamma_aa = out.find("V_GAMMA_AA")->second->pointer();
            v_gamma_bb = out.find("V_GAMMA_BB")->second->pointer();
        }
        if (meta_) {
            v_tau_a = out.find("V_TAU_A")->second->pointer();
            v_tau_b = out.find("V_TAU_B")->second->pointer();
        }
    }
     
    // => Main Loop over points <= //
    for (int Q = 0; Q < npoints; Q++) {
    
        // => Primitive variables <= //
        double rho_a;
        double rho_b;
        double gamma_aa;
        double gamma_bb;
        double tau_a;
        double tau_b;

        rho_a = rho_ap[Q];
        rho_b = rho_bp[Q];
        if (rho_a < lsda_cutoff_ || rho_b < lsda_cutoff_) {
            continue;
        } 
        if (gga_) {
            gamma_aa = gamma_aap[Q];
            gamma_bb = gamma_bbp[Q];
        }
        if (meta_) {
            tau_a = tau_ap[Q];
            tau_b = tau_bp[Q];
        }

        // => Powers of rho <= //
        double rho_a13 = pow(rho_a, 1.0/3.0);
        double rho_a43 = rho_a * rho_a13;
        double rho_a73 = rho_a * rho_a43;
        double rho_a83 = rho_a73 * rho_a13;
        double rho_b13 = pow(rho_b, 1.0/3.0);
        double rho_b43 = rho_b * rho_b13;
        double rho_b73 = rho_b * rho_b43;
        double rho_b83 = rho_b73 * rho_b13;

        // => LDA Part <= //
        
        double E;
        double E_rho_a;
        double E_rho_b;
        if (lsda_type_ == LSDA_None) {
            throw PSIEXCEPTION("WTF");
        } else if (lsda_type_ == PW92) {

            double rho = rho_a + rho_b;
            double rho_rho_a = 1.0;
            double rho_rho_b = 1.0;

            double z = (rho_a - rho_b) / (rho_a + rho_b);
            double z_rho_a =  2.0 * rho_b / ((rho_a + rho_b) * (rho_a + rho_b));
            double z_rho_b = -2.0 * rho_a / ((rho_a + rho_b) * (rho_a + rho_b));

            double Eab;
            double Eab_rho; 
            double Eab_z; 

            PW92_C(rho, z, &Eab, &Eab_rho, &Eab_z);

            double Ea;
            double Ea_rho_a; 
            double Ea_z; 

            PW92_C(rho_a, 1.0, &Ea, &Ea_rho_a, &Ea_z);

            double Eb;
            double Eb_rho_b; 
            double Eb_z; 

            PW92_C(rho_b, 1.0, &Eb, &Eb_rho_b, &Eb_z);

            E = Eab - Ea - Eb;
            E_rho_a = Eab_rho * rho_rho_a + Eab_z * z_rho_a - Ea_rho_a;
            E_rho_b = Eab_rho * rho_rho_b + Eab_z * z_rho_b - Eb_rho_b;
        }

        // => GGA Part <= //
        double s2, s2_rho_a, s2_gamma_aa, s2_rho_b, s2_gamma_bb;
        if (gga_) {
            s2 = 1.0/2.0 * gamma_aa / rho_a83 + 1.0/2.0 * gamma_bb / rho_b83;
            s2_rho_a = -4.0/3.0 * gamma_aa / (rho_a83 * rho_a); 
            s2_rho_b = -4.0/3.0 * gamma_bb / (rho_b83 * rho_b); 
            s2_gamma_aa = 1.0/2.0 / rho_a83; 
            s2_gamma_bb = 1.0/2.0 / rho_b83; 
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
        if(!Truhlar_version_) {
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
        } else {
            double M;
            double M_rho_a;
            double M_rho_b;
            double M_gamma_a;
            double M_gamma_b;
            double M_tau_a;
            double M_tau_b;

            M05c(rho_a,rho_b,gamma_aa,gamma_bb,tau_a,tau_b, &M, &M_rho_a, &M_rho_b, &M_gamma_a, &M_gamma_b, &M_tau_a, &M_tau_b);
            
            v[Q] += M;
            v_rho_a[Q] += M_rho_a;
            v_rho_b[Q] += M_rho_b;
            v_gamma_aa[Q] += M_gamma_a;
            v_gamma_bb[Q] += M_gamma_b;
            v_tau_a[Q] += M_tau_a;
            v_tau_b[Q] += M_tau_b;
        }
    }

}
void CFunctional::PW92_C(double rho, double z, double* PW92, double* PW92_rho, double* PW92_z)
{
    //  > r < //
    
    double r = _c0_*1.0/pow(rho,1.0/3.0);
    
    //  > Ac < //
    
    double Ac = _c0a_*log((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0)*(_a1a_*r+1.0)*-2.0;
    
    //  > EcP < //
    
    double EcP = _c0p_*log((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0)*(_a1p_*r+1.0)*-2.0;
    
    //  > EcF < //
    
    double EcF = _c0f_*log((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0)*(_a1f_*r+1.0)*-2.0;
    
    //  > f < //
    
    double f = (pow(z+1.0,4.0/3.0)+pow(-z+1.0,4.0/3.0)-2.0)/(_two13_*2.0-2.0);
    
    //  > E < //
    
    double E = rho*(EcP+f*(z*z*z*z)*(EcF-EcP)+(Ac*f*(z*z*z*z-1.0))/_d2fz0_);
    
    //  > PW92 < //
    
    *PW92 = E;
    
    //  > r_rho < //
    
    double r_rho = _c0_*1.0/pow(rho,4.0/3.0)*(-1.0/3.0);
    
    //  > Ac_r < //
    
    double Ac_r = _a1a_*_c0a_*log((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0)*-2.0+((_a1a_*r+1.0)*(_b2a_+_b4a_*r*2.0+_b1a_*1.0/sqrt(r)*(1.0/2.0)+_b3a_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r),2.0))/((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0);
    
    //  > EcP_r < //
    
    double EcP_r = _a1p_*_c0p_*log((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0)*-2.0+((_a1p_*r+1.0)*(_b2p_+_b4p_*r*2.0+_b1p_*1.0/sqrt(r)*(1.0/2.0)+_b3p_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r),2.0))/((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0);
    
    //  > EcF_r < //
    
    double EcF_r = _a1f_*_c0f_*log((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0)*-2.0+((_a1f_*r+1.0)*(_b2f_+_b4f_*r*2.0+_b1f_*1.0/sqrt(r)*(1.0/2.0)+_b3f_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r),2.0))/((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0);
    
    //  > f_z < //
    
    double f_z = (pow(z+1.0,1.0/3.0)*(4.0/3.0)-pow(-z+1.0,1.0/3.0)*(4.0/3.0))/(_two13_*2.0-2.0);
    
    //  > E_rho < //
    
    double E_rho = EcP+f*(z*z*z*z)*(EcF-EcP)+(Ac*f*(z*z*z*z-1.0))/_d2fz0_;
    
    //  > E_z < //
    
    double E_z = rho*(f*(z*z*z)*(EcF-EcP)*4.0+(Ac*f*(z*z*z)*4.0)/_d2fz0_);
    
    //  > E_Ac < //
    
    double E_Ac = (f*rho*(z*z*z*z-1.0))/_d2fz0_;
    
    //  > E_EcP < //
    
    double E_EcP = -rho*(f*(z*z*z*z)-1.0);
    
    //  > E_EcF < //
    
    double E_EcF = f*rho*(z*z*z*z);
    
    //  > E_f < //
    
    double E_f = rho*((z*z*z*z)*(EcF-EcP)+(Ac*(z*z*z*z-1.0))/_d2fz0_);
    
    //  > PW92_E < //
    
    double PW92_E = 1.0;
    
    //  > PW92_f < //
    
    double PW92_f = E_f*PW92_E;
    
    //  > PW92_EcF < //
    
    double PW92_EcF = E_EcF*PW92_E;
    
    //  > PW92_EcP < //
    
    double PW92_EcP = E_EcP*PW92_E;
    
    //  > PW92_Ac < //
    
    double PW92_Ac = E_Ac*PW92_E;
    
    //  > PW92_r < //
    
    double PW92_r = Ac_r*PW92_Ac+EcF_r*PW92_EcF+EcP_r*PW92_EcP;
    
    //  > PW92_rho < //
    
    *PW92_rho = PW92_r*r_rho+E_rho*PW92_E;
    
    //  > PW92_z < //
    
    *PW92_z = PW92_f*f_z+E_z*PW92_E;
}
void CFunctional::PW92_V(double r, double z, double* PW92, double* PW92_r, double* PW92_z)
{
    //  > Ac < //
    
    double Ac = _c0a_*log((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0)*(_a1a_*r+1.0)*-2.0;
    
    //  > EcP < //
    
    double EcP = _c0p_*log((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0)*(_a1p_*r+1.0)*-2.0;
    
    //  > EcF < //
    
    double EcF = _c0f_*log((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0)*(_a1f_*r+1.0)*-2.0;
    
    //  > f < //
    
    double f = (pow(z+1.0,4.0/3.0)+pow(-z+1.0,4.0/3.0)-2.0)/(_two13_*2.0-2.0);
    
    //  > E < //
    
    double E = (EcP+f*(z*z*z*z)*(EcF-EcP)+(Ac*f*(z*z*z*z-1.0))/_d2fz0_);
    
    //  > PW92 < //
    
    *PW92 = E;
    
    //  > Ac_r < //
    
    double Ac_r = _a1a_*_c0a_*log((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0)*-2.0+((_a1a_*r+1.0)*(_b2a_+_b4a_*r*2.0+_b1a_*1.0/sqrt(r)*(1.0/2.0)+_b3a_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r),2.0))/((1.0/2.0)/(_c0a_*(_b2a_*r+_b1a_*sqrt(r)+_b3a_*pow(r,3.0/2.0)+_b4a_*(r*r)))+1.0);
    
    //  > EcP_r < //
    
    double EcP_r = _a1p_*_c0p_*log((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0)*-2.0+((_a1p_*r+1.0)*(_b2p_+_b4p_*r*2.0+_b1p_*1.0/sqrt(r)*(1.0/2.0)+_b3p_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r),2.0))/((1.0/2.0)/(_c0p_*(_b2p_*r+_b1p_*sqrt(r)+_b3p_*pow(r,3.0/2.0)+_b4p_*(r*r)))+1.0);
    
    //  > EcF_r < //
    
    double EcF_r = _a1f_*_c0f_*log((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0)*-2.0+((_a1f_*r+1.0)*(_b2f_+_b4f_*r*2.0+_b1f_*1.0/sqrt(r)*(1.0/2.0)+_b3f_*sqrt(r)*(3.0/2.0))*1.0/pow(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r),2.0))/((1.0/2.0)/(_c0f_*(_b2f_*r+_b1f_*sqrt(r)+_b3f_*pow(r,3.0/2.0)+_b4f_*(r*r)))+1.0);
    
    //  > f_z < //
    
    double f_z = (pow(z+1.0,1.0/3.0)*(4.0/3.0)-pow(-z+1.0,1.0/3.0)*(4.0/3.0))/(_two13_*2.0-2.0);
    
    
    //  > E_z < //
    
    double E_z = (f*(z*z*z)*(EcF-EcP)*4.0+(Ac*f*(z*z*z)*4.0)/_d2fz0_);
    
    //  > E_Ac < //
    
    double E_Ac = (f*(z*z*z*z-1.0))/_d2fz0_;
    
    //  > E_EcP < //
    
    double E_EcP = -(f*(z*z*z*z)-1.0);
    
    //  > E_EcF < //
    
    double E_EcF = f*(z*z*z*z);
    
    //  > E_f < //
    
    double E_f = ((z*z*z*z)*(EcF-EcP)+(Ac*(z*z*z*z-1.0))/_d2fz0_);
    
    //  > PW92_E < //
    
    double PW92_E = 1.0;
    
    //  > PW92_f < //
    
    double PW92_f = E_f*PW92_E;
    
    //  > PW92_EcF < //
    
    double PW92_EcF = E_EcF*PW92_E;
    
    //  > PW92_EcP < //
    
    double PW92_EcP = E_EcP*PW92_E;
    
    //  > PW92_Ac < //
    
    double PW92_Ac = E_Ac*PW92_E;
    
    //  > PW92_r < //
    
    *PW92_r = Ac_r*PW92_Ac+EcF_r*PW92_EcF+EcP_r*PW92_EcP;
    
    //  > PW92_z < //
    
    *PW92_z = PW92_f*f_z+E_z*PW92_E;
}

void CFunctional::M05c(double RA, double RB, double GA, double GB, double TA, double TB, double* F, double* F_rho_a, double* F_rho_b, double* F_gamma_aa, double* F_gamma_bb, double* F_tau_a, double* F_tau_b)
{
      double Pi, F6, F43, Pi34, F13, RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ, P, EUEG,Denom, DenPA, DenPB, DenGA, DenGB;
      double PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA,ChiAP,ChiAG;
      double PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB,ChiBP,ChiBG;
      double sopp0, sopp1,sopp2, sopp3, sopp4, COpp;
      double U, W, dUdChiA,dUdChiB,dUdPA,dUdPB,dUdGA,dUdGB, dWdU,dWdPA,dWdPB, dWdGA,dWdGB,EUEGPA,EUEGPB;
      double DTol,F1, F2, F3, F4;

      F1 = 1.0;
      F2 = 2.0;
      F3 = 3.0;
      F4 = 4.0;

      COpp = 0.0031E0;
      sopp0= 1.00000E0;
      sopp1= 3.78569E0;
      sopp2= -14.15261E0;
      sopp3= -7.46589E0;
      sopp4= 17.94491E0;

      DTol =  1.0E-8;
      Pi = F4*atan(F1);
      F6=6.0E0;
      F43 = F4 / F3;
      Pi34 = F3 / (F4*Pi);
      F13 = F1 / F3;

      *F = 0.0;
      *F_rho_a = 0.0;
      *F_gamma_aa = 0.0;
      *F_tau_a = 0.0;
      *F_rho_b = 0.0;
      *F_gamma_bb = 0.0;
      *F_tau_b = 0.0;

      PA = RA;
      TauA = TA;
      GAA = GA; 
      if (RA > DTol && TA > DTol) {
            m05css(DTol,PA,GAA,TauA,&FA,&FPA,&FGA,&FTA,&EUA,&ChiA,&EUPA,&ChiAP,&ChiAG);
            *F += FA;
            *F_rho_a += FPA;
            *F_gamma_aa += FGA;
            *F_tau_a += FTA;
      }
      PB = RB;
      TauB = TB;
      GBB = GB; 
      if (RB > DTol && TB > DTol) {
            m05css(DTol,PB,GBB,TauB,&FB,&FPB,&FGB,&FTB,&EUB,&ChiB,&EUPB,&ChiBP,&ChiBG);
            *F += FB;
            *F_rho_b += FPB;
            *F_gamma_bb += FGB;
            *F_tau_b += FTB;
      }
      if (RA > DTol && RB > DTol) {
          P = PA + PB;
          RS = pow(Pi34/P , F13); 
          RSP = -RS/(F3*P);
          Zeta = (PA-PB)/P;
          dZdA = (F1-Zeta)/P;
          dZdB = (-F1-Zeta)/P;
          PW92_V(RS,Zeta,&PotLC,&dLdS,&dLdZ);
          EUEG = P*PotLC - EUA - EUB;
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB));
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)));
          *F += EUEG*W;
          dUdChiA =COpp/pow(F1 + COpp*(ChiA+ChiB),2.0);
          dUdChiB =COpp/pow(F1 + COpp*(ChiA+ChiB),2.0);
          dUdPA= dUdChiA*ChiAP;
          dUdPB= dUdChiB*ChiBP;
          dUdGA= dUdChiA*ChiAG;
          dUdGB= dUdChiB*ChiBG;
          dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4));
          dWdPA= dWdU*dUdPA;
          dWdPB= dWdU*dUdPB;
          dWdGA= dWdU*dUdGA;
          dWdGB= dWdU*dUdGB;
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA;
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB;
          *F_rho_a += EUEGPA*W + EUEG*dWdPA; 
          *F_gamma_aa += EUEG*dWdGA; 
          *F_rho_b += EUEGPB*W + EUEG*dWdPB;
          *F_gamma_bb += EUEG*dWdGB;
      }
} 


void CFunctional::m05css(double DTol,double PX,double GX,double TX,double* F,double* FP,double* FG,double* FT,double* EUEG,double* Chi,double* EUEGP, double* ChiP,double* ChiG)
{
      double Pt25, F1, F2, F3, F4, F5, F6, F8, F11;
      double sss0,sss1, sss2, sss3, sss4, Css;
      double Pi, Pi34, F13, F23, F43, F53, F83, F113;
      double RS, FDUEG, D, Fscc, RSP, dFsccP, dFsccG;
      double E, W, U, dFsccT, dUdChi, dWdU, dWdP, dWdG; 
      double PotLC,dLdS,dLdZ;
 
      Pt25 = 0.25E0;
      F1 = 1.0E0;
      F2 = 2.0E0;
      F3 = 3.0E0;
      F4 = 4.0E0;
      F5 = 5.0E0;
      F6 = 6.0E0;
      F8 = 8.0E0;
      F11 =11.0E0;

      Css = 0.06E0;
      sss0=  1.00000E0;
      sss1=  3.77344E0;
      sss2=  -26.04463E0;
      sss3=  30.69913E0;
      sss4=  -9.22695E0;
      
      if ((PX <= DTol)) {
        *EUEG = 0.0;
        *Chi = 0.0;
        *EUEGP = 0.0;
        *ChiP = 0.0;
        *ChiG = 0.0;
        *F  = 0.0;
        *FP = 0.0;
        *FG = 0.0;
        *FT = 0.0;
        return;
      }      
        Pi = F4*atan(F1);
        Pi34 = F3 / (F4*Pi);
        F13 = F1 / F3;
        F23 = F2 / F3;
        F43 = F2 * F23;
        F53 = F5 / F3;
        F83 = F8 / F3;
        F113 = F11 / F3;
        FDUEG = (F3/F5)*pow(F6*Pi*Pi,F23);
        RS = pow(Pi34/PX, F13);
        PW92_V(RS,F1,&PotLC,&dLdS,&dLdZ);
        *EUEG = PX*PotLC;
        D = TX - Pt25*GX/PX;
        *Chi = GX/pow(PX,F83);
        U = Css* *Chi/(F1 + Css* *Chi);
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)));
        Fscc=D/TX;
        E = Fscc*W* *EUEG;
        *F = E;
        RSP = -RS/(F3*PX);
        *ChiG = F1/pow(PX,F83);
        *ChiP = -F83* *Chi/PX;
        dFsccP=Pt25*GX/(TX*PX*PX);
        dFsccG=-Pt25/(TX*PX);
        dFsccT=Pt25*GX/(PX*TX*TX);
        dUdChi=Css/pow((F1+Css* *Chi),2.0);
        dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4));
        dWdP=dWdU*dUdChi* *ChiP;
        dWdG=dWdU*dUdChi* *ChiG;
        *EUEGP = PotLC + PX*dLdS*RSP;
        *FP = (dFsccP*W* *EUEG  + Fscc*dWdP* *EUEG + Fscc*W* *EUEGP);
        *FG = (dFsccG*W* *EUEG + Fscc*dWdG* *EUEG);
        *FT = (dFsccT*W* *EUEG);
}

}
