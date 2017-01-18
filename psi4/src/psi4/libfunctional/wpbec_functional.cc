/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/libmints/vector.h"
#include "wpbec_functional.h"
#include "utility.h"
#include "psi4/psi4-dec.h"
#include <cmath>

using namespace psi;

namespace psi {

wPBECFunctional::wPBECFunctional()
{
}
wPBECFunctional::~wPBECFunctional()
{
}
void wPBECFunctional::common_init()
{
    alpha_ = 1.0;

    meta_ = false;

    //lsda_cutoff_ = 1.0E-18;

    if (type_ == pw92c_type) {
        name_ = "   PW92C-New";
        description_ = "    New Implementation of PW92C in wPBEc-sr.\n";
        citation_ = "   Perdew and Yang, PRB, 45, 13244 (1992).\n";
        gga_ = true;
        lrc_ = false;
        omega_ = 0.0;
    } else if (type_ == pbec_type) {
        name_ = "   PBEC-New";
        description_ = "    New Implementation of PBEC in wPBEc-sr.\n";
        citation_ = "   Perdew, Burke, and Ernzerhof, PRL, 77, 3865 (1996).\n";
        gga_ = true;
        lrc_ = false;
        omega_ = 0.0;
    } else if (type_ == pw92c_sr_type) {
        name_ = "   PW92C-SR-New";
        description_ = "    New Implementation of PW92C-SR in wPBEc-sr.\n";
        citation_ = "   Paziani, Moroni, Gori-Giorgi, and Bachelet, PRB, 73, 155111 (2006).\n";
        gga_ = true;
        lrc_ = true;
        omega_ = 0.3;
    } else if (type_ == pbec_sr_type) {
        name_ = "   PBEC-SR-New";
        description_ = "    New Implementation of PBEC-SR in wPBEc-sr.\n";
        citation_ = "   Goll, Werner, Stoll, Leninger, Gori-Giorgi, and Savin, Chem. Phys., 329, 276 (2006).\n";
        gga_ = true;
        lrc_ = true;
        omega_ = 0.3;
    } else {
        throw PSIEXCEPTION("Bad wPBEC_Type.");
    }
}
void wPBECFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("wPBECFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

    // => Input variables (spin-polarized) <= //

    double* rho_ap = in.find("RHO_A")->second->pointer();
    double* rho_bp = in.find("RHO_B")->second->pointer();
    double* gamma_aap = in.find("GAMMA_AA")->second->pointer();
    double* gamma_abp = in.find("GAMMA_AB")->second->pointer();
    double* gamma_bbp = in.find("GAMMA_BB")->second->pointer();

    // => Output variables <= //

    double* v = NULL;
    double* v_rho_a = NULL;
    double* v_rho_b = NULL;
    double* v_gamma_aa = NULL;
    double* v_gamma_ab = NULL;
    double* v_gamma_bb = NULL;

    v = out.find("V")->second->pointer();
    if (deriv >=1) {
        v_rho_a = out.find("V_RHO_A")->second->pointer();
        v_rho_b = out.find("V_RHO_B")->second->pointer();
        v_gamma_aa = out.find("V_GAMMA_AA")->second->pointer();
        v_gamma_ab = out.find("V_GAMMA_AB")->second->pointer();
        v_gamma_bb = out.find("V_GAMMA_BB")->second->pointer();
    }

    // => Main Loop over points <= //
    for (int Q = 0; Q < npoints; Q++) {

        double rho_a = rho_ap[Q];
        double rho_b = rho_bp[Q];
        double gamma_aa = gamma_aap[Q];
        double gamma_ab = gamma_abp[Q];
        double gamma_bb = gamma_bbp[Q];

        double rho = rho_a + rho_b;
        double z = (rho_a - rho_b) / (rho_a + rho_b);
        double s = gamma_aa + 2.0 * gamma_ab + gamma_bb;

        if (rho < lsda_cutoff_) {
            continue;
        }

        // Condition the z variable
        const double z_eps = 5.0E-16;
        if (rho_a > rho_b && fabs(  1.0 - z) < z_eps) {
            z =  (1.0 - z_eps);
        }
        if (rho_b > rho_a && fabs(- 1.0 + z) < z_eps) {
            z = -(1.0 - z_eps);
        }

        double F = 0.0;
        double F_rho = 0.0;
        double F_z = 0.0;
        double F_s = 0.0;

        if (type_ == pw92c_type) {
            pw92c_f(rho,z,&F,&F_rho,&F_z);
        } else if (type_ == pbec_type) {
            pbec_f(rho,z,s,&F,&F_rho,&F_z,&F_s);
        } else if (type_ == pw92c_sr_type) {
            pw92c_sr_f(omega_,rho,z,&F,&F_rho,&F_z);
        } else if (type_ == pbec_sr_type) {
            pbec_sr_f(omega_,rho,z,s,&F,&F_rho,&F_z,&F_s);
        } else {
            throw PSIEXCEPTION("Bad wPBEC_Type.");
        }

        v[Q] += A * F;
        if (deriv >= 1) {
            double rho_rho_a = 1.0;
            double rho_rho_b = 1.0;
            double z_rho_a =   2.0 * rho_b / (rho * rho);
            double z_rho_b = - 2.0 * rho_a / (rho * rho);
            double s_gamma_aa = 1.0;
            double s_gamma_ab = 2.0;
            double s_gamma_bb = 1.0;

            v_rho_a[Q]    += A * F_rho * rho_rho_a + F_z * z_rho_a;
            v_rho_b[Q]    += A * F_rho * rho_rho_b + F_z * z_rho_b;
            v_gamma_aa[Q] += A * F_s * s_gamma_aa;
            v_gamma_ab[Q] += A * F_s * s_gamma_ab;
            v_gamma_bb[Q] += A * F_s * s_gamma_bb;
        }

        //outfile->Printf("rho = %11.3E, gamma = %11.3E, %11.3E %11.3E %11.3E\n", rho_a, gamma_aa, v[Q], v_rho_a[Q], v_gamma_aa[Q]);
    }
}
void wPBECFunctional::pw92c_eps(
    double rho,
    double z,
    double* eps,
    double* eps_rho,
    double* eps_z)
{
    // => Parameters <= //

//#define PW92_OLD
#ifndef PW92_OLD
    const double c0       =   6.2035049089939998E-01;
    const double two13    =   1.2599210498948732E+00;
    const double d2       =   1.7099209341613668E+00;
    const double c0a      =   1.6886900000000000E-02;
    const double a1a      =   1.1125000000000000E-01;
    const double b1a      =   1.0356999999999999E+01;
    const double b2a      =   3.6231000000000000E+00;
    const double b3a      =   8.8026000000000004E-01;
    const double b4a      =   4.9670999999999998E-01;
    const double c0p      =   3.1090699999999999E-02;
    const double a1p      =   2.1370000000000000E-01;
    const double b1p      =   7.5956999999999999E+00;
    const double b2p      =   3.5876000000000001E+00;
    const double b3p      =   1.6382000000000001E+00;
    const double b4p      =   4.9293999999999999E-01;
    const double c0f      =   1.5545349999999999E-02;
    const double a1f      =   2.0548000000000000E-01;
    const double b1f      =   1.4118900000000000E+01;
    const double b2f      =   6.1977000000000002E+00;
    const double b3f      =   3.3662000000000001E+00;
    const double b4f      =   6.2517000000000000E-01;
#else
    const double c0       =   6.2035049089939998E-01;
    const double two13    =   1.2599210498948732E+00;
    const double d2       =   1.7099210000000000E+00;
    const double c0a      =   1.6886999999999999E-02;
    const double a1a      =   1.1125000000000000E-01;
    const double b1a      =   1.0356999999999999E+01;
    const double b2a      =   3.6231000000000000E+00;
    const double b3a      =   8.8026000000000004E-01;
    const double b4a      =   4.9670999999999998E-01;
    const double c0p      =   3.1091000000000001E-02;
    const double a1p      =   2.1370000000000000E-01;
    const double b1p      =   7.5956999999999999E+00;
    const double b2p      =   3.5876000000000001E+00;
    const double b3p      =   1.6382000000000001E+00;
    const double b4p      =   4.9293999999999999E-01;
    const double c0f      =   1.5545000000000000E-02;
    const double a1f      =   2.0548000000000000E-01;
    const double b1f      =   1.4118900000000000E+01;
    const double b2f      =   6.1977000000000002E+00;
    const double b3f      =   3.3662000000000001E+00;
    const double b4f      =   6.2517000000000000E-01;
#endif

    // ==> Code <== //

    //  > r_s < //

    double r_s = c0*1.0/pow(rho,1.0/3.0);

    //  > Ac < //

    double Ac = c0a*log((1.0/2.0)/(c0a*(b2a*r_s+b1a*sqrt(r_s)+b3a*pow(r_s,3.0/2.0)+b4a*(r_s*r_s)))+1.0)*(a1a*r_s+1.0)*-2.0;

    //  > EcP < //

    double EcP = c0p*log((1.0/2.0)/(c0p*(b2p*r_s+b1p*sqrt(r_s)+b3p*pow(r_s,3.0/2.0)+b4p*(r_s*r_s)))+1.0)*(a1p*r_s+1.0)*-2.0;

    //  > EcF < //

    double EcF = c0f*log((1.0/2.0)/(c0f*(b2f*r_s+b1f*sqrt(r_s)+b3f*pow(r_s,3.0/2.0)+b4f*(r_s*r_s)))+1.0)*(a1f*r_s+1.0)*-2.0;

    //  > f < //

    double f = (pow(z+1.0,4.0/3.0)+pow(-z+1.0,4.0/3.0)-2.0)/(two13*2.0-2.0);

    //  > Z < //

    double Z = EcP+f*(z*z*z*z)*(EcF-EcP)+(Ac*f*(z*z*z*z-1.0))/d2;

    //  > eps < //

    *eps = Z;

    //  > r_s_rho < //

    double r_s_rho = c0*1.0/pow(rho,4.0/3.0)*(-1.0/3.0);

    //  > Ac_r_s < //

    double Ac_r_s = a1a*c0a*log((1.0/2.0)/(c0a*(b2a*r_s+b1a*sqrt(r_s)+b3a*pow(r_s,3.0/2.0)+b4a*(r_s*r_s)))+1.0)*-2.0+((a1a*r_s+1.0)*(b2a+b4a*r_s*2.0+b1a*1.0/sqrt(r_s)*(1.0/2.0)+b3a*sqrt(r_s)*(3.0/2.0))*1.0/pow(b2a*r_s+b1a*sqrt(r_s)+b3a*pow(r_s,3.0/2.0)+b4a*(r_s*r_s),2.0))/((1.0/2.0)/(c0a*(b2a*r_s+b1a*sqrt(r_s)+b3a*pow(r_s,3.0/2.0)+b4a*(r_s*r_s)))+1.0);

    //  > EcP_r_s < //

    double EcP_r_s = a1p*c0p*log((1.0/2.0)/(c0p*(b2p*r_s+b1p*sqrt(r_s)+b3p*pow(r_s,3.0/2.0)+b4p*(r_s*r_s)))+1.0)*-2.0+((a1p*r_s+1.0)*(b2p+b4p*r_s*2.0+b1p*1.0/sqrt(r_s)*(1.0/2.0)+b3p*sqrt(r_s)*(3.0/2.0))*1.0/pow(b2p*r_s+b1p*sqrt(r_s)+b3p*pow(r_s,3.0/2.0)+b4p*(r_s*r_s),2.0))/((1.0/2.0)/(c0p*(b2p*r_s+b1p*sqrt(r_s)+b3p*pow(r_s,3.0/2.0)+b4p*(r_s*r_s)))+1.0);

    //  > EcF_r_s < //

    double EcF_r_s = a1f*c0f*log((1.0/2.0)/(c0f*(b2f*r_s+b1f*sqrt(r_s)+b3f*pow(r_s,3.0/2.0)+b4f*(r_s*r_s)))+1.0)*-2.0+((a1f*r_s+1.0)*(b2f+b4f*r_s*2.0+b1f*1.0/sqrt(r_s)*(1.0/2.0)+b3f*sqrt(r_s)*(3.0/2.0))*1.0/pow(b2f*r_s+b1f*sqrt(r_s)+b3f*pow(r_s,3.0/2.0)+b4f*(r_s*r_s),2.0))/((1.0/2.0)/(c0f*(b2f*r_s+b1f*sqrt(r_s)+b3f*pow(r_s,3.0/2.0)+b4f*(r_s*r_s)))+1.0);

    //  > f_z < //

    double f_z = (pow(z+1.0,1.0/3.0)*(4.0/3.0)-pow(-z+1.0,1.0/3.0)*(4.0/3.0))/(two13*2.0-2.0);

    //  > Z_z < //

    double Z_z = f*(z*z*z)*(EcF-EcP)*4.0+(Ac*f*(z*z*z)*4.0)/d2;

    //  > Z_Ac < //

    double Z_Ac = (f*(z*z*z*z-1.0))/d2;

    //  > Z_EcP < //

    double Z_EcP = -f*(z*z*z*z)+1.0;

    //  > Z_EcF < //

    double Z_EcF = f*(z*z*z*z);

    //  > Z_f < //

    double Z_f = (z*z*z*z)*(EcF-EcP)+(Ac*(z*z*z*z-1.0))/d2;

    //  > eps_Z < //

    double eps_Z = 1.0;

    //  > eps_f < //

    double eps_f = Z_f*eps_Z;

    //  > eps_EcF < //

    double eps_EcF = Z_EcF*eps_Z;

    //  > eps_EcP < //

    double eps_EcP = Z_EcP*eps_Z;

    //  > eps_Ac < //

    double eps_Ac = Z_Ac*eps_Z;

    //  > eps_r_s < //

    double eps_r_s = Ac_r_s*eps_Ac+EcF_r_s*eps_EcF+EcP_r_s*eps_EcP;

    //  > eps_rho < //

    *eps_rho = eps_r_s*r_s_rho;

    //  > eps_z < //

    *eps_z = Z_z*eps_Z+eps_f*f_z;
}
void wPBECFunctional::pw92c_sr_eps(
    double omega,
    double rho,
    double z,
    double* eps,
    double* eps_rho,
    double* eps_z,
    double* eps_sr,
    double* eps_sr_rho,
    double* eps_sr_z)
{
    // => Parameters <= //

    const double Q0       =  -6.2181381739309802E-02;
    const double Qa       =   5.8460500000000000E+00;
    const double Qb       =   7.4495253826340555E+00;
    const double Qc       =   3.9174400000000000E+00;
    const double Qd       =   3.4485100000000002E+00;
    const double p1       =  -2.2655000000000002E-02;
    const double p2       =   4.3190000000000001E-01;
    const double p3       =   4.0000000000000001E-02;
    const double Bg       =   2.0711000000000035E-02;
    const double Cg       =   8.1930600000000006E-02;
    const double Dg       =  -1.2771299999999999E-02;
    const double Eg       =   1.8589800000000001E-03;
    const double Fg       =   7.5241100000000005E-01;
    const double Ad       =   7.8494900000000001E-01;
    const double q1a      =  -3.8800000000000001E-01;
    const double q2a      =   6.7600000000000005E-01;
    const double q3a      =   5.4700000000000004E-01;
    const double t1a      =  -4.9500000000000002E+00;
    const double t2a      =   1.0000000000000000E+00;
    const double t3a      =   3.1000000000000000E-01;
    const double c0       =   6.2035049089939998E-01;
    const double sqrt_2pi =   2.5066282746310002E+00;
    const double alpha2   =   2.7150535898260325E-01;
    const double two53    =   3.1748021039363992E+00;

    //  > ec < //

    double ec;
    double ec_rho;
    double ec_z;
    pw92c_eps(rho,z,&ec,&ec_rho,&ec_z);
    *eps     = ec;
    *eps_rho = ec_rho;
    *eps_z   = ec_z;

    // ==> Code <== //

    //  > P2 < //

    double P2 = pow(z+1.0,2.0/3.0)*(1.0/2.0)+pow(-z+1.0,2.0/3.0)*(1.0/2.0);

    //  > P8 < //

    double P8 = pow(z+1.0,8.0/3.0)*(1.0/2.0)+pow(-z+1.0,8.0/3.0)*(1.0/2.0);

    //  > r_s < //

    double r_s = c0*1.0/pow(rho,1.0/3.0);

    //  > g0 < //

    double g0 = exp(-Fg*r_s)*(Bg*r_s+Cg*(r_s*r_s)+Dg*(r_s*r_s*r_s)+Eg*(r_s*r_s*r_s*r_s)+1.0)*(1.0/2.0);

    //  > b0 < //

    double b0 = Ad*r_s;

    //  > D2 < //

    double D2 = 1.0/(r_s*r_s)*exp(-q3a*r_s)*(q1a*r_s+q2a*(r_s*r_s));

    //  > D3 < //

    double D3 = 1.0/(r_s*r_s*r_s)*exp(-r_s*t3a)*(r_s*t1a+(r_s*r_s)*t2a);

    //  > r_sp < //

    double r_sp = r_s*pow(2.0/(z+1.0),1.0/3.0);

    //  > r_sm < //

    double r_sm = r_s*pow(-2.0/(z-1.0),1.0/3.0);

    //  > g2p < //

    double g2p = (1.0/(r_sp*r_sp)*two53*(p1*r_sp+1.0)*(1.0/5.0))/(alpha2*(p2*r_sp+p3*(r_sp*r_sp)+1.0));

    //  > g2m < //

    double g2m = (1.0/(r_sm*r_sm)*two53*(p1*r_sm+1.0)*(1.0/5.0))/(alpha2*(p2*r_sm+p3*(r_sm*r_sm)+1.0));

    //  > C2 < //

    double C2 = 1.0/(r_s*r_s*r_s)*((z*z)*3.0-3.0)*(g0-1.0/2.0)*(1.0/8.0);

    //  > C3 < //

    double C3 = (g0*1.0/(r_s*r_s*r_s)*(z*z-1.0))/sqrt_2pi;

    //  > C4 < //

    double C4 = 1.0/(r_s*r_s*r_s)*(g2m*pow(z*(1.0/2.0)-1.0/2.0,2.0)+g2p*pow(z*(1.0/2.0)+1.0/2.0,2.0)-D2*(z*z-1.0)-(P8*1.0/(r_s*r_s)*(1.0/5.0))/alpha2)*(-9.0/6.4E1);

    //  > C5 < //

    double C5 = (1.0/(r_s*r_s*r_s)*(g2m*pow(z*(1.0/2.0)-1.0/2.0,2.0)+g2p*pow(z*(1.0/2.0)+1.0/2.0,2.0)-D3*(z*z-1.0))*(-9.0/4.0E1))/sqrt_2pi;

    //  > a1 < //

    double a1 = C3*(b0*b0*b0*b0*b0*b0)*4.0+C5*(b0*b0*b0*b0*b0*b0*b0*b0);

    //  > a2 < //

    double a2 = C2*(b0*b0*b0*b0*b0*b0)*4.0+C4*(b0*b0*b0*b0*b0*b0*b0*b0)+(b0*b0*b0*b0)*ec*6.0;

    //  > a3 < //

    double a3 = C3*(b0*b0*b0*b0*b0*b0*b0*b0);

    //  > a4 < //

    double a4 = C2*(b0*b0*b0*b0*b0*b0*b0*b0)+(b0*b0*b0*b0*b0*b0)*ec*4.0;

    //  > a5 < //

    double a5 = (b0*b0*b0*b0*b0*b0*b0*b0)*ec;

    //  > eta < //

    double eta = (omega*sqrt(r_s))/P2;

    //  > Q < //

    double Q = Q0*log((Qa*eta+Qb*(eta*eta)+Qc*(eta*eta*eta)+1.0)/(Qa*eta+Qd*(eta*eta)+1.0));

    //  > eps_lr < //

    double eps_lr = 1.0/pow((b0*b0)*(omega*omega)+1.0,4.0)*((P2*P2*P2)*Q+a1*(omega*omega*omega)+a2*(omega*omega*omega*omega)+a3*(omega*omega*omega*omega*omega)+a4*(omega*omega*omega*omega*omega*omega)+a5*(omega*omega*omega*omega*omega*omega*omega*omega));

    //  > Z < //

    double Z = ec-eps_lr;

    //  > eps_sr < //

    *eps_sr = Z;

    //  > P2_z < //

    double P2_z = 1.0/pow(z+1.0,1.0/3.0)*(1.0/3.0)-1.0/pow(-z+1.0,1.0/3.0)*(1.0/3.0);

    //  > P8_z < //

    double P8_z = pow(z+1.0,5.0/3.0)*(4.0/3.0)-pow(-z+1.0,5.0/3.0)*(4.0/3.0);

    //  > r_s_rho < //

    double r_s_rho = c0*1.0/pow(rho,4.0/3.0)*(-1.0/3.0);

    //  > g0_r_s < //

    double g0_r_s = exp(-Fg*r_s)*(Bg+Cg*r_s*2.0+Dg*(r_s*r_s)*3.0+Eg*(r_s*r_s*r_s)*4.0)*(1.0/2.0)-Fg*exp(-Fg*r_s)*(Bg*r_s+Cg*(r_s*r_s)+Dg*(r_s*r_s*r_s)+Eg*(r_s*r_s*r_s*r_s)+1.0)*(1.0/2.0);

    //  > b0_r_s < //

    double b0_r_s = Ad;

    //  > D2_r_s < //

    double D2_r_s = 1.0/(r_s*r_s*r_s)*exp(-q3a*r_s)*(q1a*r_s+q2a*(r_s*r_s))*-2.0+1.0/(r_s*r_s)*exp(-q3a*r_s)*(q1a+q2a*r_s*2.0)-q3a*1.0/(r_s*r_s)*exp(-q3a*r_s)*(q1a*r_s+q2a*(r_s*r_s));

    //  > D3_r_s < //

    double D3_r_s = 1.0/(r_s*r_s*r_s*r_s)*exp(-r_s*t3a)*(r_s*t1a+(r_s*r_s)*t2a)*-3.0+1.0/(r_s*r_s*r_s)*exp(-r_s*t3a)*(t1a+r_s*t2a*2.0)-1.0/(r_s*r_s*r_s)*t3a*exp(-r_s*t3a)*(r_s*t1a+(r_s*r_s)*t2a);

    //  > r_sp_z < //

    double r_sp_z = r_s*1.0/pow(2.0/(z+1.0),2.0/3.0)*1.0/pow(z+1.0,2.0)*(-2.0/3.0);

    //  > r_sp_r_s < //

    double r_sp_r_s = pow(2.0/(z+1.0),1.0/3.0);

    //  > r_sm_z < //

    double r_sm_z = r_s*1.0/pow(-2.0/(z-1.0),2.0/3.0)*1.0/pow(z-1.0,2.0)*(2.0/3.0);

    //  > r_sm_r_s < //

    double r_sm_r_s = pow(-2.0/(z-1.0),1.0/3.0);

    //  > g2p_r_sp < //

    double g2p_r_sp = (p1*1.0/(r_sp*r_sp)*two53*(1.0/5.0))/(alpha2*(p2*r_sp+p3*(r_sp*r_sp)+1.0))-(1.0/(r_sp*r_sp*r_sp)*two53*(p1*r_sp+1.0)*(2.0/5.0))/(alpha2*(p2*r_sp+p3*(r_sp*r_sp)+1.0))-(1.0/(r_sp*r_sp)*two53*(p2+p3*r_sp*2.0)*(p1*r_sp+1.0)*1.0/pow(p2*r_sp+p3*(r_sp*r_sp)+1.0,2.0)*(1.0/5.0))/alpha2;

    //  > g2m_r_sm < //

    double g2m_r_sm = (p1*1.0/(r_sm*r_sm)*two53*(1.0/5.0))/(alpha2*(p2*r_sm+p3*(r_sm*r_sm)+1.0))-(1.0/(r_sm*r_sm*r_sm)*two53*(p1*r_sm+1.0)*(2.0/5.0))/(alpha2*(p2*r_sm+p3*(r_sm*r_sm)+1.0))-(1.0/(r_sm*r_sm)*two53*(p2+p3*r_sm*2.0)*(p1*r_sm+1.0)*1.0/pow(p2*r_sm+p3*(r_sm*r_sm)+1.0,2.0)*(1.0/5.0))/alpha2;

    //  > C2_z < //

    double C2_z = 1.0/(r_s*r_s*r_s)*z*(g0-1.0/2.0)*(3.0/4.0);

    //  > C2_r_s < //

    double C2_r_s = 1.0/(r_s*r_s*r_s*r_s)*((z*z)*3.0-3.0)*(g0-1.0/2.0)*(-3.0/8.0);

    //  > C2_g0 < //

    double C2_g0 = 1.0/(r_s*r_s*r_s)*((z*z)*3.0-3.0)*(1.0/8.0);

    //  > C3_z < //

    double C3_z = (g0*1.0/(r_s*r_s*r_s)*z*2.0)/sqrt_2pi;

    //  > C3_r_s < //

    double C3_r_s = (g0*1.0/(r_s*r_s*r_s*r_s)*(z*z-1.0)*-3.0)/sqrt_2pi;

    //  > C3_g0 < //

    double C3_g0 = (1.0/(r_s*r_s*r_s)*(z*z-1.0))/sqrt_2pi;

    //  > C4_z < //

    double C4_z = 1.0/(r_s*r_s*r_s)*(D2*z*-2.0+g2m*(z*(1.0/2.0)-1.0/2.0)+g2p*(z*(1.0/2.0)+1.0/2.0))*(-9.0/6.4E1);

    //  > C4_P8 < //

    double C4_P8 = (1.0/(r_s*r_s*r_s*r_s*r_s)*(9.0/3.2E2))/alpha2;

    //  > C4_r_s < //

    double C4_r_s = 1.0/(r_s*r_s*r_s*r_s)*(g2m*pow(z*(1.0/2.0)-1.0/2.0,2.0)+g2p*pow(z*(1.0/2.0)+1.0/2.0,2.0)-D2*(z*z-1.0)-(P8*1.0/(r_s*r_s)*(1.0/5.0))/alpha2)*(2.7E1/6.4E1)-(P8*1.0/(r_s*r_s*r_s*r_s*r_s*r_s)*(9.0/1.6E2))/alpha2;

    //  > C4_D2 < //

    double C4_D2 = 1.0/(r_s*r_s*r_s)*(z*z-1.0)*(9.0/6.4E1);

    //  > C4_g2p < //

    double C4_g2p = 1.0/(r_s*r_s*r_s)*pow(z*(1.0/2.0)+1.0/2.0,2.0)*(-9.0/6.4E1);

    //  > C4_g2m < //

    double C4_g2m = 1.0/(r_s*r_s*r_s)*pow(z*(1.0/2.0)-1.0/2.0,2.0)*(-9.0/6.4E1);

    //  > C5_z < //

    double C5_z = (1.0/(r_s*r_s*r_s)*(D3*z*-2.0+g2m*(z*(1.0/2.0)-1.0/2.0)+g2p*(z*(1.0/2.0)+1.0/2.0))*(-9.0/4.0E1))/sqrt_2pi;

    //  > C5_r_s < //

    double C5_r_s = (1.0/(r_s*r_s*r_s*r_s)*(g2m*pow(z*(1.0/2.0)-1.0/2.0,2.0)+g2p*pow(z*(1.0/2.0)+1.0/2.0,2.0)-D3*(z*z-1.0))*(2.7E1/4.0E1))/sqrt_2pi;

    //  > C5_D3 < //

    double C5_D3 = (1.0/(r_s*r_s*r_s)*(z*z-1.0)*(9.0/4.0E1))/sqrt_2pi;

    //  > C5_g2p < //

    double C5_g2p = (1.0/(r_s*r_s*r_s)*pow(z*(1.0/2.0)+1.0/2.0,2.0)*(-9.0/4.0E1))/sqrt_2pi;

    //  > C5_g2m < //

    double C5_g2m = (1.0/(r_s*r_s*r_s)*pow(z*(1.0/2.0)-1.0/2.0,2.0)*(-9.0/4.0E1))/sqrt_2pi;

    //  > a1_b0 < //

    double a1_b0 = C3*(b0*b0*b0*b0*b0)*2.4E1+C5*(b0*b0*b0*b0*b0*b0*b0)*8.0;

    //  > a1_C3 < //

    double a1_C3 = (b0*b0*b0*b0*b0*b0)*4.0;

    //  > a1_C5 < //

    double a1_C5 = b0*b0*b0*b0*b0*b0*b0*b0;

    //  > a2_ec < //

    double a2_ec = (b0*b0*b0*b0)*6.0;

    //  > a2_b0 < //

    double a2_b0 = C2*(b0*b0*b0*b0*b0)*2.4E1+C4*(b0*b0*b0*b0*b0*b0*b0)*8.0+(b0*b0*b0)*ec*2.4E1;

    //  > a2_C2 < //

    double a2_C2 = (b0*b0*b0*b0*b0*b0)*4.0;

    //  > a2_C4 < //

    double a2_C4 = b0*b0*b0*b0*b0*b0*b0*b0;

    //  > a3_b0 < //

    double a3_b0 = C3*(b0*b0*b0*b0*b0*b0*b0)*8.0;

    //  > a3_C3 < //

    double a3_C3 = b0*b0*b0*b0*b0*b0*b0*b0;

    //  > a4_ec < //

    double a4_ec = (b0*b0*b0*b0*b0*b0)*4.0;

    //  > a4_b0 < //

    double a4_b0 = C2*(b0*b0*b0*b0*b0*b0*b0)*8.0+(b0*b0*b0*b0*b0)*ec*2.4E1;

    //  > a4_C2 < //

    double a4_C2 = b0*b0*b0*b0*b0*b0*b0*b0;

    //  > a5_ec < //

    double a5_ec = b0*b0*b0*b0*b0*b0*b0*b0;

    //  > a5_b0 < //

    double a5_b0 = (b0*b0*b0*b0*b0*b0*b0)*ec*8.0;

    //  > eta_P2 < //

    double eta_P2 = -1.0/(P2*P2)*omega*sqrt(r_s);

    //  > eta_r_s < //

    double eta_r_s = (omega*1.0/sqrt(r_s)*(1.0/2.0))/P2;

    //  > Q_eta < //

    double Q_eta = (Q0*((Qa+Qb*eta*2.0+Qc*(eta*eta)*3.0)/(Qa*eta+Qd*(eta*eta)+1.0)-(Qa+Qd*eta*2.0)*1.0/pow(Qa*eta+Qd*(eta*eta)+1.0,2.0)*(Qa*eta+Qb*(eta*eta)+Qc*(eta*eta*eta)+1.0))*(Qa*eta+Qd*(eta*eta)+1.0))/(Qa*eta+Qb*(eta*eta)+Qc*(eta*eta*eta)+1.0);

    //  > eps_lr_P2 < //

    double eps_lr_P2 = (P2*P2)*Q*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0)*3.0;

    //  > eps_lr_b0 < //

    double eps_lr_b0 = b0*(omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,5.0)*((P2*P2*P2)*Q+a1*(omega*omega*omega)+a2*(omega*omega*omega*omega)+a3*(omega*omega*omega*omega*omega)+a4*(omega*omega*omega*omega*omega*omega)+a5*(omega*omega*omega*omega*omega*omega*omega*omega))*-8.0;

    //  > eps_lr_a1 < //

    double eps_lr_a1 = (omega*omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > eps_lr_a2 < //

    double eps_lr_a2 = (omega*omega*omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > eps_lr_a3 < //

    double eps_lr_a3 = (omega*omega*omega*omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > eps_lr_a4 < //

    double eps_lr_a4 = (omega*omega*omega*omega*omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > eps_lr_a5 < //

    double eps_lr_a5 = (omega*omega*omega*omega*omega*omega*omega*omega)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > eps_lr_Q < //

    double eps_lr_Q = (P2*P2*P2)*1.0/pow((b0*b0)*(omega*omega)+1.0,4.0);

    //  > Z_ec < //

    double Z_ec = 1.0;

    //  > Z_eps_lr < //

    double Z_eps_lr = -1.0;

    //  > eps_sr_Z < //

    double eps_sr_Z = 1.0;

    //  > eps_sr_eps_lr < //

    double eps_sr_eps_lr = Z_eps_lr*eps_sr_Z;

    //  > eps_sr_Q < //

    double eps_sr_Q = eps_lr_Q*eps_sr_eps_lr;

    //  > eps_sr_eta < //

    double eps_sr_eta = Q_eta*eps_sr_Q;

    //  > eps_sr_a5 < //

    double eps_sr_a5 = eps_lr_a5*eps_sr_eps_lr;

    //  > eps_sr_a4 < //

    double eps_sr_a4 = eps_lr_a4*eps_sr_eps_lr;

    //  > eps_sr_a3 < //

    double eps_sr_a3 = eps_lr_a3*eps_sr_eps_lr;

    //  > eps_sr_a2 < //

    double eps_sr_a2 = eps_lr_a2*eps_sr_eps_lr;

    //  > eps_sr_a1 < //

    double eps_sr_a1 = eps_lr_a1*eps_sr_eps_lr;

    //  > eps_sr_C5 < //

    double eps_sr_C5 = a1_C5*eps_sr_a1;

    //  > eps_sr_C4 < //

    double eps_sr_C4 = a2_C4*eps_sr_a2;

    //  > eps_sr_C3 < //

    double eps_sr_C3 = a1_C3*eps_sr_a1+a3_C3*eps_sr_a3;

    //  > eps_sr_C2 < //

    double eps_sr_C2 = a2_C2*eps_sr_a2+a4_C2*eps_sr_a4;

    //  > eps_sr_g2m < //

    double eps_sr_g2m = C4_g2m*eps_sr_C4+C5_g2m*eps_sr_C5;

    //  > eps_sr_g2p < //

    double eps_sr_g2p = C4_g2p*eps_sr_C4+C5_g2p*eps_sr_C5;

    //  > eps_sr_r_sm < //

    double eps_sr_r_sm = eps_sr_g2m*g2m_r_sm;

    //  > eps_sr_r_sp < //

    double eps_sr_r_sp = eps_sr_g2p*g2p_r_sp;

    //  > eps_sr_D3 < //

    double eps_sr_D3 = C5_D3*eps_sr_C5;

    //  > eps_sr_D2 < //

    double eps_sr_D2 = C4_D2*eps_sr_C4;

    //  > eps_sr_b0 < //

    double eps_sr_b0 = a1_b0*eps_sr_a1+a2_b0*eps_sr_a2+a3_b0*eps_sr_a3+a4_b0*eps_sr_a4+a5_b0*eps_sr_a5+eps_lr_b0*eps_sr_eps_lr;

    //  > eps_sr_g0 < //

    double eps_sr_g0 = C2_g0*eps_sr_C2+C3_g0*eps_sr_C3;

    //  > eps_sr_r_s < //

    double eps_sr_r_s = C2_r_s*eps_sr_C2+C3_r_s*eps_sr_C3+C4_r_s*eps_sr_C4+C5_r_s*eps_sr_C5+D2_r_s*eps_sr_D2+D3_r_s*eps_sr_D3+b0_r_s*eps_sr_b0+eps_sr_eta*eta_r_s+eps_sr_g0*g0_r_s+eps_sr_r_sm*r_sm_r_s+eps_sr_r_sp*r_sp_r_s;

    //  > eps_sr_P8 < //

    double eps_sr_P8 = C4_P8*eps_sr_C4;

    //  > eps_sr_P2 < //

    double eps_sr_P2 = eps_lr_P2*eps_sr_eps_lr+eps_sr_eta*eta_P2;

    //  > eps_sr_ec < //

    double eps_sr_ec = Z_ec*eps_sr_Z+a2_ec*eps_sr_a2+a4_ec*eps_sr_a4+a5_ec*eps_sr_a5;

    //  > eps_sr_rho < //

    *eps_sr_rho = ec_rho*eps_sr_ec+eps_sr_r_s*r_s_rho;

    //  > eps_sr_z < //

    *eps_sr_z = C2_z*eps_sr_C2+C3_z*eps_sr_C3+C4_z*eps_sr_C4+C5_z*eps_sr_C5+P2_z*eps_sr_P2+P8_z*eps_sr_P8+ec_z*eps_sr_ec+eps_sr_r_sm*r_sm_z+eps_sr_r_sp*r_sp_z;
}
void wPBECFunctional::pw92c_f(
    double rho,
    double z,
    double* f,
    double* f_rho,
    double* f_z)
{
    double e;
    double e_rho;
    double e_z;

    pw92c_eps(rho,z,&e,&e_rho,&e_z);

    *f = rho * e;
    *f_rho = e + rho * e_rho;
    *f_z = rho * e_z;
}
void wPBECFunctional::pbec_f(
    double rho,
    double z,
    double s,
    double* f,
    double* f_rho,
    double* f_z,
    double* f_s)
{
    // => Parameters <= //

    const double B_PBE    =   6.6724550603149219E-02;
    const double G        =   3.1090690869654901E-02;
    const double cf       =   1.9846863952198559E+00;

    //  > pw92c_eps < //

    double ec;
    double ec_rho;
    double ec_z;
    pw92c_eps(rho,z,&ec,&ec_rho,&ec_z);

    // ==> Code <== //

    //  > ks < //

    double ks = cf*pow(rho,1.0/6.0);

    //  > P < //

    double P = pow(z+1.0,2.0/3.0)*(1.0/2.0)+pow(-z+1.0,2.0/3.0)*(1.0/2.0);

    //  > t2 < //

    double t2 = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho)*s*(1.0/4.0);

    //  > A < //

    double A = B_PBE/(G*(exp(-(1.0/(P*P*P)*ec)/G)-1.0));

    //  > H < //

    double H = G*(P*P*P)*log((B_PBE*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > Z < //

    double Z = rho*(H+ec);

    //  > f < //

    *f = Z;

    //  > ks_rho < //

    double ks_rho = cf*1.0/pow(rho,5.0/6.0)*(1.0/6.0);

    //  > P_z < //

    double P_z = 1.0/pow(z+1.0,1.0/3.0)*(1.0/3.0)-1.0/pow(-z+1.0,1.0/3.0)*(1.0/3.0);

    //  > t2_rho < //

    double t2_rho = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho*rho)*s*(-1.0/2.0);

    //  > t2_s < //

    double t2_s = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho)*(1.0/4.0);

    //  > t2_ks < //

    double t2_ks = 1.0/(P*P)*1.0/(ks*ks*ks)*1.0/(rho*rho)*s*(-1.0/2.0);

    //  > t2_P < //

    double t2_P = 1.0/(P*P*P)*1.0/(ks*ks)*1.0/(rho*rho)*s*(-1.0/2.0);

    //  > A_ec < //

    double A_ec = B_PBE*1.0/(G*G)*1.0/(P*P*P)*exp(-(1.0/(P*P*P)*ec)/G)*1.0/pow(exp(-(1.0/(P*P*P)*ec)/G)-1.0,2.0);

    //  > A_P < //

    double A_P = B_PBE*1.0/(G*G)*1.0/(P*P*P*P)*ec*exp(-(1.0/(P*P*P)*ec)/G)*1.0/pow(exp(-(1.0/(P*P*P)*ec)/G)-1.0,2.0)*-3.0;

    //  > H_P < //

    double H_P = G*(P*P)*log((B_PBE*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0)*3.0;

    //  > H_t2 < //

    double H_t2 = (G*(P*P*P)*((B_PBE*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+(A*B_PBE*t2)/(G*(A*t2+(A*A)*(t2*t2)+1.0))-(B_PBE*t2*(A+(A*A)*t2*2.0)*(A*t2+1.0)*1.0/pow(A*t2+(A*A)*(t2*t2)+1.0,2.0))/G))/((B_PBE*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > H_A < //

    double H_A = (G*(P*P*P)*((B_PBE*(t2*t2))/(G*(A*t2+(A*A)*(t2*t2)+1.0))-(B_PBE*t2*(t2+A*(t2*t2)*2.0)*(A*t2+1.0)*1.0/pow(A*t2+(A*A)*(t2*t2)+1.0,2.0))/G))/((B_PBE*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > Z_rho < //

    double Z_rho = H+ec;

    //  > Z_ec < //

    double Z_ec = rho;

    //  > Z_H < //

    double Z_H = rho;

    //  > f_Z < //

    double f_Z = 1.0;

    //  > f_H < //

    double f_H = Z_H*f_Z;

    //  > f_A < //

    double f_A = H_A*f_H;

    //  > f_t2 < //

    double f_t2 = H_t2*f_H;

    //  > f_P < //

    double f_P = A_P*f_A+H_P*f_H+f_t2*t2_P;

    //  > f_ks < //

    double f_ks = f_t2*t2_ks;

    //  > f_ec < //

    double f_ec = A_ec*f_A+Z_ec*f_Z;

    //  > f_rho < //

    *f_rho = Z_rho*f_Z+ec_rho*f_ec+f_ks*ks_rho+f_t2*t2_rho;

    //  > f_z < //

    *f_z = P_z*f_P+ec_z*f_ec;

    //  > f_s < //

    *f_s = f_t2*t2_s;
}
void wPBECFunctional::pw92c_sr_f(
    double omega,
    double rho,
    double z,
    double* f,
    double* f_rho,
    double* f_z)
{
    double e;
    double e_rho;
    double e_z;

    double e_sr;
    double e_sr_rho;
    double e_sr_z;

    pw92c_sr_eps(omega,rho,z,&e,&e_rho,&e_z,&e_sr,&e_sr_rho,&e_sr_z);

    *f = rho * e_sr;
    *f_rho = e_sr + rho * e_sr_rho;
    *f_z = rho * e_sr_z;
}
void wPBECFunctional::pbec_sr_f(
    double omega,
    double rho,
    double z,
    double s,
    double* f,
    double* f_rho,
    double* f_z,
    double* f_s)
{
    // => Parameters <= //

    const double a_c      =   2.7799999999999998E+00;
    const double B_PBE    =   6.6724550603149219E-02;
    const double G        =   3.1090690869654901E-02;
    const double cf       =   1.9846863952198559E+00;

    //  > pw92c_eps < //

    double ec;
    double ec_rho;
    double ec_z;
    double ec_sr;
    double ec_sr_rho;
    double ec_sr_z;

    pw92c_sr_eps(omega,rho,z,&ec,&ec_rho,&ec_z,&ec_sr,&ec_sr_rho,&ec_sr_z);

    if (ec_sr >= -1.0E-20) {
        *f = 0.0;
        *f_rho = 0.0;
        *f_z = 0.0;
        *f_s = 0.0;
        return;
    }

    // ==> Code <== //

    //  > B < //

    double B = B_PBE*pow(ec_sr/ec,a_c);

    //  > ks < //

    double ks = cf*pow(rho,1.0/6.0);

    //  > P < //

    double P = pow(z+1.0,2.0/3.0)*(1.0/2.0)+pow(-z+1.0,2.0/3.0)*(1.0/2.0);

    //  > t2 < //

    double t2 = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho)*s*(1.0/4.0);

    //  > V < //

    double V = -(1.0/(P*P*P)*ec_sr)/G;

    //  > U < //

    double U = exp(V)-1.0;
    if (fabs(V) < 1.0E-5) {
        U = V + V*V/2.0 + V*V*V/6.0 + V*V*V*V/24.0;
    }

    //  > A < //

    double A = B/(G*U);

    //  > H < //

    double H = G*(P*P*P)*log((B*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > Z < //

    double Z = rho*(H+ec_sr);

    //  > f < //

    *f = Z;

    //  > B_ec < //

    double B_ec = -B_PBE*a_c*1.0/(ec*ec)*ec_sr*pow(ec_sr/ec,a_c-1.0);

    //  > B_ec_sr < //

    double B_ec_sr = (B_PBE*a_c*pow(ec_sr/ec,a_c-1.0))/ec;

    //  > ks_rho < //

    double ks_rho = cf*1.0/pow(rho,5.0/6.0)*(1.0/6.0);

    //  > P_z < //

    double P_z = 1.0/pow(z+1.0,1.0/3.0)*(1.0/3.0)-1.0/pow(-z+1.0,1.0/3.0)*(1.0/3.0);

    //  > t2_rho < //

    double t2_rho = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho*rho)*s*(-1.0/2.0);

    //  > t2_s < //

    double t2_s = 1.0/(P*P)*1.0/(ks*ks)*1.0/(rho*rho)*(1.0/4.0);

    //  > t2_ks < //

    double t2_ks = 1.0/(P*P)*1.0/(ks*ks*ks)*1.0/(rho*rho)*s*(-1.0/2.0);

    //  > t2_P < //

    double t2_P = 1.0/(P*P*P)*1.0/(ks*ks)*1.0/(rho*rho)*s*(-1.0/2.0);

    //  > V_ec_sr < //

    double V_ec_sr = -1.0/(P*P*P)/G;

    //  > V_P < //

    double V_P = (1.0/(P*P*P*P)*ec_sr*3.0)/G;

    //  > U_V < //

    double U_V = exp(V);

    //  > A_B < //

    double A_B = 1.0/(G*U);

    //  > A_U < //

    double A_U = -(B*1.0/(U*U))/G;

    //  > H_B < //

    double H_B = ((P*P*P)*t2*(A*t2+1.0))/(((B*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0)*(A*t2+(A*A)*(t2*t2)+1.0));

    //  > H_P < //

    double H_P = G*(P*P)*log((B*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0)*3.0;

    //  > H_t2 < //

    double H_t2 = (G*(P*P*P)*((B*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+(A*B*t2)/(G*(A*t2+(A*A)*(t2*t2)+1.0))-(B*t2*(A+(A*A)*t2*2.0)*(A*t2+1.0)*1.0/pow(A*t2+(A*A)*(t2*t2)+1.0,2.0))/G))/((B*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > H_A < //

    double H_A = (G*(P*P*P)*((B*(t2*t2))/(G*(A*t2+(A*A)*(t2*t2)+1.0))-(B*t2*(t2+A*(t2*t2)*2.0)*(A*t2+1.0)*1.0/pow(A*t2+(A*A)*(t2*t2)+1.0,2.0))/G))/((B*t2*(A*t2+1.0))/(G*(A*t2+(A*A)*(t2*t2)+1.0))+1.0);

    //  > Z_rho < //

    double Z_rho = H+ec_sr;

    //  > Z_ec_sr < //

    double Z_ec_sr = rho;

    //  > Z_H < //

    double Z_H = rho;

    //  > f_Z < //

    double f_Z = 1.0;

    //  > f_H < //

    double f_H = Z_H*f_Z;

    //  > f_A < //

    double f_A = H_A*f_H;

    //  > f_U < //

    double f_U = A_U*f_A;

    //  > f_V < //

    double f_V = U_V*f_U;

    //  > f_t2 < //

    double f_t2 = H_t2*f_H;

    //  > f_P < //

    double f_P = H_P*f_H+V_P*f_V+f_t2*t2_P;

    //  > f_ks < //

    double f_ks = f_t2*t2_ks;

    //  > f_B < //

    double f_B = A_B*f_A+H_B*f_H;

    //  > f_ec_sr < //

    double f_ec_sr = B_ec_sr*f_B+V_ec_sr*f_V+Z_ec_sr*f_Z;

    //  > f_ec < //

    double f_ec = B_ec*f_B;

    //  > f_rho < //

    *f_rho = Z_rho*f_Z+ec_rho*f_ec+ec_sr_rho*f_ec_sr+f_ks*ks_rho+f_t2*t2_rho;

    //  > f_z < //

    *f_z = P_z*f_P+ec_z*f_ec+ec_sr_z*f_ec_sr;

    //  > f_s < //

    *f_s = f_t2*t2_s;
}

}
