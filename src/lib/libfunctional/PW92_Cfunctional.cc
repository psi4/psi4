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

#include <libmints/vector.h>
#include "PW92_Cfunctional.h"
#include "utility.h"
#include <cmath>

using namespace psi;

namespace psi {

PW92_CFunctional::PW92_CFunctional()
{
    name_ = "PW92_C";
    description_ = "    PW92 LSDA Correlation\n";
    citation_ = "    J.P. Perdew and Y. Wang, Phys. Rev. B., 45(23), 13244, 1992\n";
    alpha_ = 1.0;
    omega_ = 0.0;
    lrc_ = false;
    gga_ = false;
    meta_ = false;
    parameters_["c"] =   6.2035049089939986E-01;
    parameters_["two_13"] =   1.2599210498948732E+00;
    parameters_["d2fz0"] =   1.7099210000000000E+00;
    parameters_["Aa"] =   1.6886999999999999E-02;
    parameters_["a1a"] =   1.1125000000000000E-01;
    parameters_["b1a"] =   1.0356999999999999E+01;
    parameters_["b2a"] =   3.6231000000000000E+00;
    parameters_["b3a"] =   8.8026000000000004E-01;
    parameters_["b4a"] =   4.9670999999999998E-01;
    parameters_["c0p"] =   3.1091000000000001E-02;
    parameters_["a1p"] =   2.1370000000000000E-01;
    parameters_["b1p"] =   7.5956999999999999E+00;
    parameters_["b2p"] =   3.5876000000000001E+00;
    parameters_["b3p"] =   1.6382000000000001E+00;
    parameters_["b4p"] =   4.9293999999999999E-01;
    parameters_["c0f"] =   1.5545000000000000E-02;
    parameters_["a1f"] =   2.0548000000000000E-01;
    parameters_["b1f"] =   1.4118900000000000E+01;
    parameters_["b2f"] =   6.1977000000000002E+00;
    parameters_["b3f"] =   3.3662000000000001E+00;
    parameters_["b4f"] =   6.2517000000000000E-01;
}
PW92_CFunctional::~PW92_CFunctional()
{
}
void PW92_CFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    double c = parameters_["c"];
    double two_13 = parameters_["two_13"];
    double d2fz0 = parameters_["d2fz0"];
    double Aa = parameters_["Aa"];
    double a1a = parameters_["a1a"];
    double b1a = parameters_["b1a"];
    double b2a = parameters_["b2a"];
    double b3a = parameters_["b3a"];
    double b4a = parameters_["b4a"];
    double c0p = parameters_["c0p"];
    double a1p = parameters_["a1p"];
    double b1p = parameters_["b1p"];
    double b2p = parameters_["b2p"];
    double b3p = parameters_["b3p"];
    double b4p = parameters_["b4p"];
    double c0f = parameters_["c0f"];
    double a1f = parameters_["a1f"];
    double b1f = parameters_["b1f"];
    double b2f = parameters_["b2f"];
    double b3f = parameters_["b3f"];
    double b4f = parameters_["b4f"];

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

    if (true) {
        rho_ap = in.find("RHO_A")->second->pointer();
        rho_bp = in.find("RHO_B")->second->pointer();
    }
    if (gga_) {  
        gamma_aap = in.find("GAMMA_AA")->second->pointer();
        gamma_abp = in.find("GAMMA_AB")->second->pointer();
        gamma_bbp = in.find("GAMMA_BB")->second->pointer();
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

    // => Loop over points <= //

    for (int Q = 0; Q < npoints; Q++) {

        // Input variables 
        double rho_a;
        double rho_b;
        double gamma_aa;
        double gamma_ab;
        double gamma_bb;
        double tau_a;
        double tau_b;

        if (true) {
            rho_a = rho_ap[Q];
            rho_b = rho_bp[Q];
        }        
        if (gga_) {
            gamma_aa = gamma_aap[Q];
            gamma_ab = gamma_abp[Q];
            gamma_bb = gamma_bbp[Q];
        }        
        if (meta_) {
            tau_a = tau_ap[Q];
            tau_b = tau_bp[Q];
        }        

        // Definitions (asymptotics to prevent numerical problems)
        if (rho_a < lsda_cutoff_ && rho_b < lsda_cutoff_) {
            continue;
        } else if (rho_a < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t14893 = rho_a+rho_b;
                double t14894 = 1.0/pow(t14893,1.0/3.0);
                double t14895 = c*t14894;
                double t14896 = sqrt(t14895);
                double t14897 = pow(t14895,3.0/2.0);
                double t14898 = c*c;
                double t14899 = 1.0/pow(t14893,2.0/3.0);
                double t14900 = 1.0/c0p;
                double t14901 = b1p*t14896;
                double t14902 = b3p*t14897;
                double t14903 = b4p*t14898*t14899;
                double t14904 = b2p*c*t14894;
                double t14905 = t14901+t14902+t14903+t14904;
                double t14906 = 1.0/t14905;
                double t14907 = t14900*t14906*(1.0/2.0);
                double t14908 = t14907+1.0;
                double t14909 = log(t14908);
                double t14910 = a1p*c*t14894;
                double t14911 = t14910+1.0;
                double t14912 = c0p*t14911*t14909*2.0;
                v[Q] += scale * (-t14893*(t14912-((pow(2.0,1.0/3.0)*2.0-2.0)*(t14912-c0f*log((1.0/2.0)/(c0f*(b1f*t14896+b3f*t14897+b2f*c*t14894+b4f*t14898*t14899))+1.0)*(a1f*c*t14894+1.0)*2.0))/(two_13*2.0-2.0)));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t14914 = rho_a+rho_b;
                double t14915 = 1.0/pow(t14914,1.0/3.0);
                double t14916 = c*t14915;
                double t14917 = sqrt(t14916);
                double t14918 = b1f*t14917;
                double t14919 = pow(t14916,3.0/2.0);
                double t14920 = b3f*t14919;
                double t14921 = c*c;
                double t14922 = 1.0/pow(t14914,2.0/3.0);
                double t14923 = b4f*t14921*t14922;
                double t14924 = b2f*c*t14915;
                double t14925 = t14920+t14923+t14924+t14918;
                double t14926 = 1.0/pow(t14914,4.0/3.0);
                double t14927 = b1p*t14917;
                double t14928 = b3p*t14919;
                double t14929 = b4p*t14921*t14922;
                double t14930 = b2p*c*t14915;
                double t14931 = t14930+t14927+t14928+t14929;
                double t14932 = 1.0/pow(t14914,5.0/3.0);
                double t14933 = 1.0/sqrt(t14916);
                double t14934 = 1.0/c0f;
                double t14935 = 1.0/t14925;
                double t14936 = t14934*t14935*(1.0/2.0);
                double t14937 = t14936+1.0;
                double t14938 = 1.0/c0p;
                double t14939 = 1.0/t14931;
                double t14940 = t14938*t14939*(1.0/2.0);
                double t14941 = t14940+1.0;
                double t14942 = a1p*c*t14915;
                double t14943 = t14942+1.0;
                double t14944 = 1.0/t14941;
                double t14945 = 1.0/(t14931*t14931);
                double t14946 = b4p*t14921*t14932*(2.0/3.0);
                double t14947 = b2p*c*t14926*(1.0/3.0);
                double t14948 = b1p*c*t14933*t14926*(1.0/6.0);
                double t14949 = b3p*c*t14917*t14926*(1.0/2.0);
                double t14950 = t14946+t14947+t14948+t14949;
                double t14951 = t14950*t14943*t14944*t14945;
                double t14952 = log(t14941);
                double t14953 = pow(2.0,1.0/3.0);
                double t14954 = t14953*2.0;
                double t14955 = t14954-2.0;
                double t14956 = two_13*2.0;
                double t14957 = t14956-2.0;
                double t14958 = 1.0/t14957;
                double t14959 = log(t14937);
                double t14960 = a1f*c*t14915;
                double t14961 = t14960+1.0;
                v_rho_a[Q] += scale * (t14914*(-t14951+t14955*t14958*(t14951-(1.0/(t14925*t14925)*t14961*(b2f*c*t14926*(1.0/3.0)+b4f*t14921*t14932*(2.0/3.0)+b1f*c*t14933*t14926*(1.0/6.0)+b3f*c*t14917*t14926*(1.0/2.0)))/t14937+a1f*c*c0f*t14926*t14959*(2.0/3.0)-a1p*c*c0p*t14952*t14926*(2.0/3.0))+a1p*c*c0p*t14952*t14926*(2.0/3.0))-c0p*t14943*t14952*2.0-t14955*t14958*(c0f*t14961*t14959*2.0-c0p*t14943*t14952*2.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t14963 = rho_a+rho_b;
                double t14964 = 1.0/pow(t14963,1.0/3.0);
                double t14965 = c*t14964;
                double t14966 = sqrt(t14965);
                double t14967 = b1f*t14966;
                double t14968 = pow(t14965,3.0/2.0);
                double t14969 = b3f*t14968;
                double t14970 = c*c;
                double t14971 = 1.0/pow(t14963,2.0/3.0);
                double t14972 = b4f*t14970*t14971;
                double t14973 = b2f*c*t14964;
                double t14974 = t14972+t14973+t14967+t14969;
                double t14975 = 1.0/pow(t14963,4.0/3.0);
                double t14976 = b1p*t14966;
                double t14977 = b3p*t14968;
                double t14978 = b4p*t14970*t14971;
                double t14979 = b2p*c*t14964;
                double t14980 = t14976+t14977+t14978+t14979;
                double t14981 = 1.0/pow(t14963,5.0/3.0);
                double t14982 = 1.0/sqrt(t14965);
                double t14983 = 1.0/c0f;
                double t14984 = 1.0/t14974;
                double t14985 = t14983*t14984*(1.0/2.0);
                double t14986 = t14985+1.0;
                double t14987 = 1.0/c0p;
                double t14988 = 1.0/t14980;
                double t14989 = t14987*t14988*(1.0/2.0);
                double t14990 = t14989+1.0;
                double t14991 = a1p*c*t14964;
                double t14992 = t14991+1.0;
                double t14993 = 1.0/t14990;
                double t14994 = 1.0/(t14980*t14980);
                double t14995 = b4p*t14970*t14981*(2.0/3.0);
                double t14996 = b2p*c*t14975*(1.0/3.0);
                double t14997 = b1p*c*t14982*t14975*(1.0/6.0);
                double t14998 = b3p*c*t14966*t14975*(1.0/2.0);
                double t14999 = t14995+t14996+t14997+t14998;
                double t15000 = t14992*t14993*t14994*t14999;
                double t15001 = log(t14990);
                double t15002 = pow(2.0,1.0/3.0);
                double t15003 = t15002*2.0;
                double t15004 = t15003-2.0;
                double t15005 = two_13*2.0;
                double t15006 = t15005-2.0;
                double t15007 = 1.0/t15006;
                double t15008 = log(t14986);
                double t15009 = a1f*c*t14964;
                double t15010 = t15009+1.0;
                v_rho_b[Q] += scale * (t14963*(-t15000+t15004*t15007*(t15000-(1.0/(t14974*t14974)*t15010*(b2f*c*t14975*(1.0/3.0)+b4f*t14970*t14981*(2.0/3.0)+b1f*c*t14982*t14975*(1.0/6.0)+b3f*c*t14966*t14975*(1.0/2.0)))/t14986+a1f*c*c0f*t14975*t15008*(2.0/3.0)-a1p*c*c0p*t14975*t15001*(2.0/3.0))+a1p*c*c0p*t14975*t15001*(2.0/3.0))-c0p*t14992*t15001*2.0-t15004*t15007*(c0f*t15010*t15008*2.0-c0p*t14992*t15001*2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t15017 = rho_a+rho_b;
                double t15018 = 1.0/pow(t15017,1.0/3.0);
                double t15019 = c*t15018;
                double t15020 = sqrt(t15019);
                double t15021 = b1f*t15020;
                double t15022 = pow(t15019,3.0/2.0);
                double t15023 = b3f*t15022;
                double t15024 = c*c;
                double t15025 = 1.0/pow(t15017,2.0/3.0);
                double t15026 = b4f*t15024*t15025;
                double t15027 = b2f*c*t15018;
                double t15028 = t15021+t15023+t15026+t15027;
                double t15029 = 1.0/pow(t15017,7.0/3.0);
                double t15030 = 1.0/pow(t15017,8.0/3.0);
                double t15031 = 1.0/sqrt(t15019);
                double t15032 = b1p*t15020;
                double t15033 = b3p*t15022;
                double t15034 = b4p*t15024*t15025;
                double t15035 = b2p*c*t15018;
                double t15036 = t15032+t15033+t15034+t15035;
                double t15037 = 1.0/pow(t15019,3.0/2.0);
                double t15038 = a1f*c*t15018;
                double t15039 = t15038+1.0;
                double t15040 = 1.0/c0f;
                double t15041 = 1.0/t15028;
                double t15042 = t15040*t15041*(1.0/2.0);
                double t15043 = t15042+1.0;
                double t15044 = 1.0/t15043;
                double t15045 = 1.0/pow(t15017,4.0/3.0);
                double t15054 = 1.0/pow(t15017,5.0/3.0);
                double t15056 = b4f*t15024*t15054*(2.0/3.0);
                double t15057 = b2f*c*t15045*(1.0/3.0);
                double t15058 = b1f*c*t15031*t15045*(1.0/6.0);
                double t15059 = b3f*c*t15020*t15045*(1.0/2.0);
                double t15046 = t15056+t15057+t15058+t15059;
                double t15047 = a1p*c*t15018;
                double t15048 = t15047+1.0;
                double t15049 = 1.0/c0p;
                double t15050 = 1.0/t15036;
                double t15051 = t15050*t15049*(1.0/2.0);
                double t15052 = t15051+1.0;
                double t15053 = 1.0/t15052;
                double t15061 = b4p*t15024*t15054*(2.0/3.0);
                double t15062 = b2p*c*t15045*(1.0/3.0);
                double t15063 = b1p*c*t15031*t15045*(1.0/6.0);
                double t15064 = b3p*c*t15020*t15045*(1.0/2.0);
                double t15055 = t15061+t15062+t15063+t15064;
                double t15060 = t15046*t15046;
                double t15065 = t15055*t15055;
                double t15066 = 1.0/(t15028*t15028);
                double t15067 = 1.0/(t15036*t15036);
                double t15068 = b4p*t15030*t15024*(1.0E1/9.0);
                double t15069 = b2p*c*t15029*(4.0/9.0);
                double t15070 = b1p*c*t15031*t15029*(2.0/9.0);
                double t15071 = b3p*c*t15020*t15029*(2.0/3.0);
                double t15072 = b3p*t15030*t15031*t15024*(1.0/1.2E1);
                double t15073 = t15070+t15071+t15072+t15068+t15069-b1p*t15030*t15024*t15037*(1.0/3.6E1);
                double t15074 = t15053*t15073*t15048*t15067;
                double t15075 = 1.0/(t15036*t15036*t15036);
                double t15076 = 1.0/(t15052*t15052);
                double t15077 = 1.0/(t15036*t15036*t15036*t15036);
                double t15078 = t15065*t15048*t15049*t15076*t15077*(1.0/2.0);
                double t15079 = log(t15052);
                double t15080 = a1p*c*t15053*t15045*t15055*t15067*(2.0/3.0);
                double t15081 = pow(2.0,1.0/3.0);
                double t15082 = t15081*2.0;
                double t15083 = t15082-2.0;
                double t15084 = two_13*2.0;
                double t15085 = t15084-2.0;
                double t15086 = 1.0/t15085;
                double t15087 = log(t15043);
                v_rho_a_rho_a[Q] += scale * (t15017*(t15080+t15074+t15078-t15083*t15086*(t15080+t15074+t15078-t15053*t15065*t15048*t15075*2.0+t15060*t15044*1.0/(t15028*t15028*t15028)*t15039*2.0-t15044*t15039*t15066*(b2f*c*t15029*(4.0/9.0)+b4f*t15030*t15024*(1.0E1/9.0)+b1f*c*t15031*t15029*(2.0/9.0)+b3f*c*t15020*t15029*(2.0/3.0)-b1f*t15030*t15024*t15037*(1.0/3.6E1)+b3f*t15030*t15031*t15024*(1.0/1.2E1))-t15040*t15060*1.0/(t15043*t15043)*1.0/(t15028*t15028*t15028*t15028)*t15039*(1.0/2.0)+a1f*c*c0f*t15029*t15087*(8.0/9.0)-a1p*c*c0p*t15029*t15079*(8.0/9.0)-a1f*c*t15044*t15045*t15046*t15066*(2.0/3.0))-t15053*t15065*t15048*t15075*2.0-a1p*c*c0p*t15029*t15079*(8.0/9.0))-t15083*t15086*(t15044*t15046*t15039*t15066-t15053*t15055*t15048*t15067-a1f*c*c0f*t15045*t15087*(2.0/3.0)+a1p*c*c0p*t15045*t15079*(2.0/3.0))*2.0-t15053*t15055*t15048*t15067*2.0+a1p*c*c0p*t15045*t15079*(4.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t15089 = rho_a+rho_b;
                double t15090 = 1.0/pow(t15089,1.0/3.0);
                double t15091 = c*t15090;
                double t15092 = sqrt(t15091);
                double t15093 = b1f*t15092;
                double t15094 = pow(t15091,3.0/2.0);
                double t15095 = b3f*t15094;
                double t15096 = c*c;
                double t15097 = 1.0/pow(t15089,2.0/3.0);
                double t15098 = b4f*t15096*t15097;
                double t15099 = b2f*c*t15090;
                double t15100 = t15093+t15095+t15098+t15099;
                double t15101 = 1.0/pow(t15089,7.0/3.0);
                double t15102 = 1.0/pow(t15089,8.0/3.0);
                double t15103 = 1.0/sqrt(t15091);
                double t15104 = b1p*t15092;
                double t15105 = b3p*t15094;
                double t15106 = b4p*t15096*t15097;
                double t15107 = b2p*c*t15090;
                double t15108 = t15104+t15105+t15106+t15107;
                double t15109 = 1.0/pow(t15091,3.0/2.0);
                double t15110 = a1f*c*t15090;
                double t15111 = t15110+1.0;
                double t15112 = 1.0/c0f;
                double t15113 = 1.0/t15100;
                double t15114 = t15112*t15113*(1.0/2.0);
                double t15115 = t15114+1.0;
                double t15116 = 1.0/t15115;
                double t15117 = 1.0/pow(t15089,4.0/3.0);
                double t15126 = 1.0/pow(t15089,5.0/3.0);
                double t15128 = b4f*t15126*t15096*(2.0/3.0);
                double t15129 = b2f*c*t15117*(1.0/3.0);
                double t15130 = b1f*c*t15103*t15117*(1.0/6.0);
                double t15131 = b3f*c*t15117*t15092*(1.0/2.0);
                double t15118 = t15130+t15131+t15128+t15129;
                double t15119 = a1p*c*t15090;
                double t15120 = t15119+1.0;
                double t15121 = 1.0/c0p;
                double t15122 = 1.0/t15108;
                double t15123 = t15121*t15122*(1.0/2.0);
                double t15124 = t15123+1.0;
                double t15125 = 1.0/t15124;
                double t15133 = b4p*t15126*t15096*(2.0/3.0);
                double t15134 = b2p*c*t15117*(1.0/3.0);
                double t15135 = b1p*c*t15103*t15117*(1.0/6.0);
                double t15136 = b3p*c*t15117*t15092*(1.0/2.0);
                double t15127 = t15133+t15134+t15135+t15136;
                double t15132 = t15118*t15118;
                double t15137 = t15127*t15127;
                double t15138 = 1.0/(t15100*t15100);
                double t15139 = 1.0/(t15108*t15108);
                double t15140 = b4p*t15102*t15096*(1.0E1/9.0);
                double t15141 = b2p*c*t15101*(4.0/9.0);
                double t15142 = b1p*c*t15101*t15103*(2.0/9.0);
                double t15143 = b3p*c*t15101*t15092*(2.0/3.0);
                double t15144 = b3p*t15102*t15103*t15096*(1.0/1.2E1);
                double t15145 = t15140+t15141+t15142+t15143+t15144-b1p*t15102*t15109*t15096*(1.0/3.6E1);
                double t15146 = t15120*t15125*t15145*t15139;
                double t15147 = 1.0/(t15108*t15108*t15108);
                double t15148 = 1.0/(t15124*t15124);
                double t15149 = 1.0/(t15108*t15108*t15108*t15108);
                double t15150 = t15120*t15121*t15137*t15148*t15149*(1.0/2.0);
                double t15151 = log(t15124);
                double t15152 = a1p*c*t15125*t15117*t15127*t15139*(2.0/3.0);
                double t15153 = pow(2.0,1.0/3.0);
                double t15154 = t15153*2.0;
                double t15155 = t15154-2.0;
                double t15156 = two_13*2.0;
                double t15157 = t15156-2.0;
                double t15158 = 1.0/t15157;
                double t15159 = log(t15115);
                v_rho_a_rho_b[Q] += scale * (t15089*(t15150+t15152+t15146-t15155*t15158*(t15150+t15152+t15146-t15120*t15125*t15137*t15147*2.0+1.0/(t15100*t15100*t15100)*t15111*t15132*t15116*2.0-t15111*t15116*t15138*(b2f*c*t15101*(4.0/9.0)+b4f*t15102*t15096*(1.0E1/9.0)+b1f*c*t15101*t15103*(2.0/9.0)+b3f*c*t15101*t15092*(2.0/3.0)-b1f*t15102*t15109*t15096*(1.0/3.6E1)+b3f*t15102*t15103*t15096*(1.0/1.2E1))-1.0/(t15100*t15100*t15100*t15100)*t15111*t15112*t15132*1.0/(t15115*t15115)*(1.0/2.0)+a1f*c*c0f*t15101*t15159*(8.0/9.0)-a1p*c*c0p*t15101*t15151*(8.0/9.0)-a1f*c*t15116*t15117*t15118*t15138*(2.0/3.0))-t15120*t15125*t15137*t15147*2.0-a1p*c*c0p*t15101*t15151*(8.0/9.0))-t15155*t15158*(t15111*t15116*t15118*t15138-t15120*t15125*t15127*t15139-a1f*c*c0f*t15117*t15159*(2.0/3.0)+a1p*c*c0p*t15151*t15117*(2.0/3.0))*2.0-t15120*t15125*t15127*t15139*2.0+a1p*c*c0p*t15151*t15117*(4.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t15161 = rho_a+rho_b;
                double t15162 = 1.0/pow(t15161,1.0/3.0);
                double t15163 = c*t15162;
                double t15164 = sqrt(t15163);
                double t15165 = b1f*t15164;
                double t15166 = pow(t15163,3.0/2.0);
                double t15167 = b3f*t15166;
                double t15168 = c*c;
                double t15169 = 1.0/pow(t15161,2.0/3.0);
                double t15170 = b4f*t15168*t15169;
                double t15171 = b2f*c*t15162;
                double t15172 = t15170+t15171+t15165+t15167;
                double t15173 = 1.0/pow(t15161,7.0/3.0);
                double t15174 = 1.0/pow(t15161,8.0/3.0);
                double t15175 = 1.0/sqrt(t15163);
                double t15176 = b1p*t15164;
                double t15177 = b3p*t15166;
                double t15178 = b4p*t15168*t15169;
                double t15179 = b2p*c*t15162;
                double t15180 = t15176+t15177+t15178+t15179;
                double t15181 = 1.0/pow(t15163,3.0/2.0);
                double t15182 = a1f*c*t15162;
                double t15183 = t15182+1.0;
                double t15184 = 1.0/c0f;
                double t15185 = 1.0/t15172;
                double t15186 = t15184*t15185*(1.0/2.0);
                double t15187 = t15186+1.0;
                double t15188 = 1.0/t15187;
                double t15189 = 1.0/pow(t15161,4.0/3.0);
                double t15198 = 1.0/pow(t15161,5.0/3.0);
                double t15200 = b4f*t15168*t15198*(2.0/3.0);
                double t15201 = b2f*c*t15189*(1.0/3.0);
                double t15202 = b1f*c*t15175*t15189*(1.0/6.0);
                double t15203 = b3f*c*t15164*t15189*(1.0/2.0);
                double t15190 = t15200+t15201+t15202+t15203;
                double t15191 = a1p*c*t15162;
                double t15192 = t15191+1.0;
                double t15193 = 1.0/c0p;
                double t15194 = 1.0/t15180;
                double t15195 = t15193*t15194*(1.0/2.0);
                double t15196 = t15195+1.0;
                double t15197 = 1.0/t15196;
                double t15205 = b4p*t15168*t15198*(2.0/3.0);
                double t15206 = b2p*c*t15189*(1.0/3.0);
                double t15207 = b1p*c*t15175*t15189*(1.0/6.0);
                double t15208 = b3p*c*t15164*t15189*(1.0/2.0);
                double t15199 = t15205+t15206+t15207+t15208;
                double t15204 = t15190*t15190;
                double t15209 = t15199*t15199;
                double t15210 = 1.0/(t15172*t15172);
                double t15211 = 1.0/(t15180*t15180);
                double t15212 = b4p*t15174*t15168*(1.0E1/9.0);
                double t15213 = b2p*c*t15173*(4.0/9.0);
                double t15214 = b1p*c*t15173*t15175*(2.0/9.0);
                double t15215 = b3p*c*t15164*t15173*(2.0/3.0);
                double t15216 = b3p*t15174*t15175*t15168*(1.0/1.2E1);
                double t15217 = t15212+t15213+t15214+t15215+t15216-b1p*t15181*t15174*t15168*(1.0/3.6E1);
                double t15218 = t15211*t15217*t15192*t15197;
                double t15219 = 1.0/(t15180*t15180*t15180);
                double t15220 = 1.0/(t15196*t15196);
                double t15221 = 1.0/(t15180*t15180*t15180*t15180);
                double t15222 = t15220*t15221*t15209*t15192*t15193*(1.0/2.0);
                double t15223 = log(t15196);
                double t15224 = a1p*c*t15211*t15197*t15189*t15199*(2.0/3.0);
                double t15225 = pow(2.0,1.0/3.0);
                double t15226 = t15225*2.0;
                double t15227 = t15226-2.0;
                double t15228 = two_13*2.0;
                double t15229 = t15228-2.0;
                double t15230 = 1.0/t15229;
                double t15231 = log(t15187);
                v_rho_b_rho_b[Q] += scale * (t15161*(t15222+t15224+t15218-t15230*t15227*(t15222+t15224+t15218-t15209*t15192*t15219*t15197*2.0+t15204*1.0/(t15172*t15172*t15172)*t15183*t15188*2.0-t15210*t15183*t15188*(b2f*c*t15173*(4.0/9.0)+b4f*t15174*t15168*(1.0E1/9.0)+b1f*c*t15173*t15175*(2.0/9.0)+b3f*c*t15164*t15173*(2.0/3.0)-b1f*t15181*t15174*t15168*(1.0/3.6E1)+b3f*t15174*t15175*t15168*(1.0/1.2E1))-t15204*1.0/(t15172*t15172*t15172*t15172)*t15183*t15184*1.0/(t15187*t15187)*(1.0/2.0)+a1f*c*c0f*t15231*t15173*(8.0/9.0)-a1p*c*c0p*t15223*t15173*(8.0/9.0)-a1f*c*t15210*t15190*t15188*t15189*(2.0/3.0))-t15209*t15192*t15219*t15197*2.0-a1p*c*c0p*t15223*t15173*(8.0/9.0))-t15230*t15227*(t15210*t15190*t15183*t15188-t15211*t15192*t15197*t15199-a1f*c*c0f*t15231*t15189*(2.0/3.0)+a1p*c*c0p*t15223*t15189*(2.0/3.0))*2.0-t15211*t15192*t15197*t15199*2.0+a1p*c*c0p*t15223*t15189*(4.0/3.0));
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t15258 = rho_a+rho_b;
                double t15259 = 1.0/pow(t15258,1.0/3.0);
                double t15260 = c*t15259;
                double t15261 = sqrt(t15260);
                double t15262 = pow(t15260,3.0/2.0);
                double t15263 = c*c;
                double t15264 = 1.0/pow(t15258,2.0/3.0);
                double t15265 = 1.0/c0p;
                double t15266 = b1p*t15261;
                double t15267 = b3p*t15262;
                double t15268 = b4p*t15263*t15264;
                double t15269 = b2p*c*t15259;
                double t15270 = t15266+t15267+t15268+t15269;
                double t15271 = 1.0/t15270;
                double t15272 = t15271*t15265*(1.0/2.0);
                double t15273 = t15272+1.0;
                double t15274 = log(t15273);
                double t15275 = a1p*c*t15259;
                double t15276 = t15275+1.0;
                double t15277 = c0p*t15274*t15276*2.0;
                v[Q] += scale * (-t15258*(t15277-((pow(2.0,1.0/3.0)*2.0-2.0)*(t15277-c0f*log((1.0/2.0)/(c0f*(b1f*t15261+b3f*t15262+b2f*c*t15259+b4f*t15263*t15264))+1.0)*(a1f*c*t15259+1.0)*2.0))/(two_13*2.0-2.0)));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t15279 = rho_a+rho_b;
                double t15280 = 1.0/pow(t15279,1.0/3.0);
                double t15281 = c*t15280;
                double t15282 = sqrt(t15281);
                double t15283 = b1f*t15282;
                double t15284 = pow(t15281,3.0/2.0);
                double t15285 = b3f*t15284;
                double t15286 = c*c;
                double t15287 = 1.0/pow(t15279,2.0/3.0);
                double t15288 = b4f*t15286*t15287;
                double t15289 = b2f*c*t15280;
                double t15290 = t15283+t15285+t15288+t15289;
                double t15291 = 1.0/pow(t15279,4.0/3.0);
                double t15292 = b1p*t15282;
                double t15293 = b3p*t15284;
                double t15294 = b4p*t15286*t15287;
                double t15295 = b2p*c*t15280;
                double t15296 = t15292+t15293+t15294+t15295;
                double t15297 = 1.0/pow(t15279,5.0/3.0);
                double t15298 = 1.0/sqrt(t15281);
                double t15299 = 1.0/c0f;
                double t15300 = 1.0/t15290;
                double t15301 = t15300*t15299*(1.0/2.0);
                double t15302 = t15301+1.0;
                double t15303 = 1.0/c0p;
                double t15304 = 1.0/t15296;
                double t15305 = t15303*t15304*(1.0/2.0);
                double t15306 = t15305+1.0;
                double t15307 = a1p*c*t15280;
                double t15308 = t15307+1.0;
                double t15309 = 1.0/t15306;
                double t15310 = 1.0/(t15296*t15296);
                double t15311 = b4p*t15286*t15297*(2.0/3.0);
                double t15312 = b2p*c*t15291*(1.0/3.0);
                double t15313 = b1p*c*t15291*t15298*(1.0/6.0);
                double t15314 = b3p*c*t15282*t15291*(1.0/2.0);
                double t15315 = t15311+t15312+t15313+t15314;
                double t15316 = t15310*t15315*t15308*t15309;
                double t15317 = log(t15306);
                double t15318 = pow(2.0,1.0/3.0);
                double t15319 = t15318*2.0;
                double t15320 = t15319-2.0;
                double t15321 = two_13*2.0;
                double t15322 = t15321-2.0;
                double t15323 = 1.0/t15322;
                double t15324 = log(t15302);
                double t15325 = a1f*c*t15280;
                double t15326 = t15325+1.0;
                v_rho_a[Q] += scale * (t15279*(-t15316+t15320*t15323*(t15316-(1.0/(t15290*t15290)*t15326*(b2f*c*t15291*(1.0/3.0)+b4f*t15286*t15297*(2.0/3.0)+b1f*c*t15291*t15298*(1.0/6.0)+b3f*c*t15282*t15291*(1.0/2.0)))/t15302+a1f*c*c0f*t15324*t15291*(2.0/3.0)-a1p*c*c0p*t15317*t15291*(2.0/3.0))+a1p*c*c0p*t15317*t15291*(2.0/3.0))-c0p*t15308*t15317*2.0-t15320*t15323*(c0f*t15324*t15326*2.0-c0p*t15308*t15317*2.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t15328 = rho_a+rho_b;
                double t15329 = 1.0/pow(t15328,1.0/3.0);
                double t15330 = c*t15329;
                double t15331 = sqrt(t15330);
                double t15332 = b1f*t15331;
                double t15333 = pow(t15330,3.0/2.0);
                double t15334 = b3f*t15333;
                double t15335 = c*c;
                double t15336 = 1.0/pow(t15328,2.0/3.0);
                double t15337 = b4f*t15335*t15336;
                double t15338 = b2f*c*t15329;
                double t15339 = t15332+t15334+t15337+t15338;
                double t15340 = 1.0/pow(t15328,4.0/3.0);
                double t15341 = b1p*t15331;
                double t15342 = b3p*t15333;
                double t15343 = b4p*t15335*t15336;
                double t15344 = b2p*c*t15329;
                double t15345 = t15341+t15342+t15343+t15344;
                double t15346 = 1.0/pow(t15328,5.0/3.0);
                double t15347 = 1.0/sqrt(t15330);
                double t15348 = 1.0/c0f;
                double t15349 = 1.0/t15339;
                double t15350 = t15348*t15349*(1.0/2.0);
                double t15351 = t15350+1.0;
                double t15352 = 1.0/c0p;
                double t15353 = 1.0/t15345;
                double t15354 = t15352*t15353*(1.0/2.0);
                double t15355 = t15354+1.0;
                double t15356 = a1p*c*t15329;
                double t15357 = t15356+1.0;
                double t15358 = 1.0/t15355;
                double t15359 = 1.0/(t15345*t15345);
                double t15360 = b4p*t15335*t15346*(2.0/3.0);
                double t15361 = b2p*c*t15340*(1.0/3.0);
                double t15362 = b1p*c*t15340*t15347*(1.0/6.0);
                double t15363 = b3p*c*t15331*t15340*(1.0/2.0);
                double t15364 = t15360+t15361+t15362+t15363;
                double t15365 = t15364*t15357*t15358*t15359;
                double t15366 = log(t15355);
                double t15367 = pow(2.0,1.0/3.0);
                double t15368 = t15367*2.0;
                double t15369 = t15368-2.0;
                double t15370 = two_13*2.0;
                double t15371 = t15370-2.0;
                double t15372 = 1.0/t15371;
                double t15373 = log(t15351);
                double t15374 = a1f*c*t15329;
                double t15375 = t15374+1.0;
                v_rho_b[Q] += scale * (t15328*(-t15365+t15372*t15369*(t15365-(1.0/(t15339*t15339)*t15375*(b2f*c*t15340*(1.0/3.0)+b4f*t15335*t15346*(2.0/3.0)+b1f*c*t15340*t15347*(1.0/6.0)+b3f*c*t15331*t15340*(1.0/2.0)))/t15351+a1f*c*c0f*t15340*t15373*(2.0/3.0)-a1p*c*c0p*t15340*t15366*(2.0/3.0))+a1p*c*c0p*t15340*t15366*(2.0/3.0))-c0p*t15357*t15366*2.0-t15372*t15369*(c0f*t15373*t15375*2.0-c0p*t15357*t15366*2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t15382 = rho_a+rho_b;
                double t15383 = 1.0/pow(t15382,1.0/3.0);
                double t15384 = c*t15383;
                double t15385 = sqrt(t15384);
                double t15386 = b1f*t15385;
                double t15387 = pow(t15384,3.0/2.0);
                double t15388 = b3f*t15387;
                double t15389 = c*c;
                double t15390 = 1.0/pow(t15382,2.0/3.0);
                double t15391 = b4f*t15390*t15389;
                double t15392 = b2f*c*t15383;
                double t15393 = t15391+t15392+t15386+t15388;
                double t15394 = 1.0/pow(t15382,7.0/3.0);
                double t15395 = 1.0/pow(t15382,8.0/3.0);
                double t15396 = 1.0/sqrt(t15384);
                double t15397 = b1p*t15385;
                double t15398 = b3p*t15387;
                double t15399 = b4p*t15390*t15389;
                double t15400 = b2p*c*t15383;
                double t15401 = t15400+t15397+t15398+t15399;
                double t15402 = 1.0/pow(t15384,3.0/2.0);
                double t15403 = a1f*c*t15383;
                double t15404 = t15403+1.0;
                double t15405 = 1.0/c0f;
                double t15406 = 1.0/t15393;
                double t15407 = t15405*t15406*(1.0/2.0);
                double t15408 = t15407+1.0;
                double t15409 = 1.0/t15408;
                double t15410 = 1.0/pow(t15382,4.0/3.0);
                double t15419 = 1.0/pow(t15382,5.0/3.0);
                double t15421 = b4f*t15419*t15389*(2.0/3.0);
                double t15422 = b2f*c*t15410*(1.0/3.0);
                double t15423 = b1f*c*t15410*t15396*(1.0/6.0);
                double t15424 = b3f*c*t15410*t15385*(1.0/2.0);
                double t15411 = t15421+t15422+t15423+t15424;
                double t15412 = a1p*c*t15383;
                double t15413 = t15412+1.0;
                double t15414 = 1.0/c0p;
                double t15415 = 1.0/t15401;
                double t15416 = t15414*t15415*(1.0/2.0);
                double t15417 = t15416+1.0;
                double t15418 = 1.0/t15417;
                double t15426 = b4p*t15419*t15389*(2.0/3.0);
                double t15427 = b2p*c*t15410*(1.0/3.0);
                double t15428 = b1p*c*t15410*t15396*(1.0/6.0);
                double t15429 = b3p*c*t15410*t15385*(1.0/2.0);
                double t15420 = t15426+t15427+t15428+t15429;
                double t15425 = t15411*t15411;
                double t15430 = t15420*t15420;
                double t15431 = 1.0/(t15393*t15393);
                double t15432 = 1.0/(t15401*t15401);
                double t15433 = b4p*t15395*t15389*(1.0E1/9.0);
                double t15434 = b2p*c*t15394*(4.0/9.0);
                double t15435 = b1p*c*t15394*t15396*(2.0/9.0);
                double t15436 = b3p*c*t15385*t15394*(2.0/3.0);
                double t15437 = b3p*t15395*t15396*t15389*(1.0/1.2E1);
                double t15438 = t15433+t15434+t15435+t15436+t15437-b1p*t15402*t15395*t15389*(1.0/3.6E1);
                double t15439 = t15413*t15432*t15418*t15438;
                double t15440 = 1.0/(t15401*t15401*t15401);
                double t15441 = 1.0/(t15417*t15417);
                double t15442 = 1.0/(t15401*t15401*t15401*t15401);
                double t15443 = t15430*t15413*t15414*t15441*t15442*(1.0/2.0);
                double t15444 = log(t15417);
                double t15445 = a1p*c*t15410*t15420*t15432*t15418*(2.0/3.0);
                double t15446 = pow(2.0,1.0/3.0);
                double t15447 = t15446*2.0;
                double t15448 = t15447-2.0;
                double t15449 = two_13*2.0;
                double t15450 = t15449-2.0;
                double t15451 = 1.0/t15450;
                double t15452 = log(t15408);
                v_rho_a_rho_a[Q] += scale * (t15382*(t15443+t15445+t15439-t15451*t15448*(t15443+t15445+t15439-t15430*t15413*t15440*t15418*2.0+t15404*t15425*t15409*1.0/(t15393*t15393*t15393)*2.0-t15404*t15431*t15409*(b2f*c*t15394*(4.0/9.0)+b4f*t15395*t15389*(1.0E1/9.0)+b1f*c*t15394*t15396*(2.0/9.0)+b3f*c*t15385*t15394*(2.0/3.0)-b1f*t15402*t15395*t15389*(1.0/3.6E1)+b3f*t15395*t15396*t15389*(1.0/1.2E1))-t15404*t15405*t15425*1.0/(t15408*t15408)*1.0/(t15393*t15393*t15393*t15393)*(1.0/2.0)+a1f*c*c0f*t15452*t15394*(8.0/9.0)-a1p*c*c0p*t15444*t15394*(8.0/9.0)-a1f*c*t15410*t15411*t15431*t15409*(2.0/3.0))-t15430*t15413*t15440*t15418*2.0-a1p*c*c0p*t15444*t15394*(8.0/9.0))-t15451*t15448*(t15411*t15404*t15431*t15409-t15420*t15413*t15432*t15418-a1f*c*c0f*t15410*t15452*(2.0/3.0)+a1p*c*c0p*t15410*t15444*(2.0/3.0))*2.0-t15420*t15413*t15432*t15418*2.0+a1p*c*c0p*t15410*t15444*(4.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t15454 = rho_a+rho_b;
                double t15455 = 1.0/pow(t15454,1.0/3.0);
                double t15456 = c*t15455;
                double t15457 = sqrt(t15456);
                double t15458 = b1f*t15457;
                double t15459 = pow(t15456,3.0/2.0);
                double t15460 = b3f*t15459;
                double t15461 = c*c;
                double t15462 = 1.0/pow(t15454,2.0/3.0);
                double t15463 = b4f*t15461*t15462;
                double t15464 = b2f*c*t15455;
                double t15465 = t15460+t15463+t15464+t15458;
                double t15466 = 1.0/pow(t15454,7.0/3.0);
                double t15467 = 1.0/pow(t15454,8.0/3.0);
                double t15468 = 1.0/sqrt(t15456);
                double t15469 = b1p*t15457;
                double t15470 = b3p*t15459;
                double t15471 = b4p*t15461*t15462;
                double t15472 = b2p*c*t15455;
                double t15473 = t15470+t15471+t15472+t15469;
                double t15474 = 1.0/pow(t15456,3.0/2.0);
                double t15475 = a1f*c*t15455;
                double t15476 = t15475+1.0;
                double t15477 = 1.0/c0f;
                double t15478 = 1.0/t15465;
                double t15479 = t15477*t15478*(1.0/2.0);
                double t15480 = t15479+1.0;
                double t15481 = 1.0/t15480;
                double t15482 = 1.0/pow(t15454,4.0/3.0);
                double t15491 = 1.0/pow(t15454,5.0/3.0);
                double t15493 = b4f*t15461*t15491*(2.0/3.0);
                double t15494 = b2f*c*t15482*(1.0/3.0);
                double t15495 = b1f*c*t15482*t15468*(1.0/6.0);
                double t15496 = b3f*c*t15482*t15457*(1.0/2.0);
                double t15483 = t15493+t15494+t15495+t15496;
                double t15484 = a1p*c*t15455;
                double t15485 = t15484+1.0;
                double t15486 = 1.0/c0p;
                double t15487 = 1.0/t15473;
                double t15488 = t15486*t15487*(1.0/2.0);
                double t15489 = t15488+1.0;
                double t15490 = 1.0/t15489;
                double t15498 = b4p*t15461*t15491*(2.0/3.0);
                double t15499 = b2p*c*t15482*(1.0/3.0);
                double t15500 = b1p*c*t15482*t15468*(1.0/6.0);
                double t15501 = b3p*c*t15482*t15457*(1.0/2.0);
                double t15492 = t15500+t15501+t15498+t15499;
                double t15497 = t15483*t15483;
                double t15502 = t15492*t15492;
                double t15503 = 1.0/(t15465*t15465);
                double t15504 = 1.0/(t15473*t15473);
                double t15505 = b4p*t15461*t15467*(1.0E1/9.0);
                double t15506 = b2p*c*t15466*(4.0/9.0);
                double t15507 = b1p*c*t15466*t15468*(2.0/9.0);
                double t15508 = b3p*c*t15457*t15466*(2.0/3.0);
                double t15509 = b3p*t15461*t15467*t15468*(1.0/1.2E1);
                double t15510 = t15505+t15506+t15507+t15508+t15509-b1p*t15461*t15474*t15467*(1.0/3.6E1);
                double t15511 = t15510*t15504*t15490*t15485;
                double t15512 = 1.0/(t15473*t15473*t15473);
                double t15513 = 1.0/(t15489*t15489);
                double t15514 = 1.0/(t15473*t15473*t15473*t15473);
                double t15515 = t15502*t15513*t15514*t15485*t15486*(1.0/2.0);
                double t15516 = log(t15489);
                double t15517 = a1p*c*t15504*t15490*t15482*t15492*(2.0/3.0);
                double t15518 = pow(2.0,1.0/3.0);
                double t15519 = t15518*2.0;
                double t15520 = t15519-2.0;
                double t15521 = two_13*2.0;
                double t15522 = t15521-2.0;
                double t15523 = 1.0/t15522;
                double t15524 = log(t15480);
                v_rho_a_rho_b[Q] += scale * (t15454*(t15511+t15515+t15517-t15520*t15523*(t15511+t15515+t15517-t15502*t15512*t15490*t15485*2.0+t15481*1.0/(t15465*t15465*t15465)*t15476*t15497*2.0-t15503*t15481*t15476*(b2f*c*t15466*(4.0/9.0)+b4f*t15461*t15467*(1.0E1/9.0)+b1f*c*t15466*t15468*(2.0/9.0)+b3f*c*t15457*t15466*(2.0/3.0)-b1f*t15461*t15474*t15467*(1.0/3.6E1)+b3f*t15461*t15467*t15468*(1.0/1.2E1))-1.0/(t15480*t15480)*1.0/(t15465*t15465*t15465*t15465)*t15476*t15477*t15497*(1.0/2.0)+a1f*c*c0f*t15524*t15466*(8.0/9.0)-a1p*c*c0p*t15516*t15466*(8.0/9.0)-a1f*c*t15503*t15481*t15482*t15483*(2.0/3.0))-t15502*t15512*t15490*t15485*2.0-a1p*c*c0p*t15516*t15466*(8.0/9.0))-t15520*t15523*(t15503*t15481*t15483*t15476-t15504*t15490*t15492*t15485-a1f*c*c0f*t15524*t15482*(2.0/3.0)+a1p*c*c0p*t15516*t15482*(2.0/3.0))*2.0-t15504*t15490*t15492*t15485*2.0+a1p*c*c0p*t15516*t15482*(4.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t15526 = rho_a+rho_b;
                double t15527 = 1.0/pow(t15526,1.0/3.0);
                double t15528 = c*t15527;
                double t15529 = sqrt(t15528);
                double t15530 = b1f*t15529;
                double t15531 = pow(t15528,3.0/2.0);
                double t15532 = b3f*t15531;
                double t15533 = c*c;
                double t15534 = 1.0/pow(t15526,2.0/3.0);
                double t15535 = b4f*t15533*t15534;
                double t15536 = b2f*c*t15527;
                double t15537 = t15530+t15532+t15535+t15536;
                double t15538 = 1.0/pow(t15526,7.0/3.0);
                double t15539 = 1.0/pow(t15526,8.0/3.0);
                double t15540 = 1.0/sqrt(t15528);
                double t15541 = b1p*t15529;
                double t15542 = b3p*t15531;
                double t15543 = b4p*t15533*t15534;
                double t15544 = b2p*c*t15527;
                double t15545 = t15541+t15542+t15543+t15544;
                double t15546 = 1.0/pow(t15528,3.0/2.0);
                double t15547 = a1f*c*t15527;
                double t15548 = t15547+1.0;
                double t15549 = 1.0/c0f;
                double t15550 = 1.0/t15537;
                double t15551 = t15550*t15549*(1.0/2.0);
                double t15552 = t15551+1.0;
                double t15553 = 1.0/t15552;
                double t15554 = 1.0/pow(t15526,4.0/3.0);
                double t15563 = 1.0/pow(t15526,5.0/3.0);
                double t15565 = b4f*t15533*t15563*(2.0/3.0);
                double t15566 = b2f*c*t15554*(1.0/3.0);
                double t15567 = b1f*c*t15540*t15554*(1.0/6.0);
                double t15568 = b3f*c*t15554*t15529*(1.0/2.0);
                double t15555 = t15565+t15566+t15567+t15568;
                double t15556 = a1p*c*t15527;
                double t15557 = t15556+1.0;
                double t15558 = 1.0/c0p;
                double t15559 = 1.0/t15545;
                double t15560 = t15558*t15559*(1.0/2.0);
                double t15561 = t15560+1.0;
                double t15562 = 1.0/t15561;
                double t15570 = b4p*t15533*t15563*(2.0/3.0);
                double t15571 = b2p*c*t15554*(1.0/3.0);
                double t15572 = b1p*c*t15540*t15554*(1.0/6.0);
                double t15573 = b3p*c*t15554*t15529*(1.0/2.0);
                double t15564 = t15570+t15571+t15572+t15573;
                double t15569 = t15555*t15555;
                double t15574 = t15564*t15564;
                double t15575 = 1.0/(t15537*t15537);
                double t15576 = 1.0/(t15545*t15545);
                double t15577 = b4p*t15533*t15539*(1.0E1/9.0);
                double t15578 = b2p*c*t15538*(4.0/9.0);
                double t15579 = b1p*c*t15540*t15538*(2.0/9.0);
                double t15580 = b3p*c*t15529*t15538*(2.0/3.0);
                double t15581 = b3p*t15540*t15533*t15539*(1.0/1.2E1);
                double t15582 = t15580+t15581+t15577+t15578+t15579-b1p*t15533*t15546*t15539*(1.0/3.6E1);
                double t15583 = t15562*t15582*t15557*t15576;
                double t15584 = 1.0/(t15545*t15545*t15545);
                double t15585 = 1.0/(t15561*t15561);
                double t15586 = 1.0/(t15545*t15545*t15545*t15545);
                double t15587 = t15574*t15557*t15558*t15585*t15586*(1.0/2.0);
                double t15588 = log(t15561);
                double t15589 = a1p*c*t15562*t15554*t15564*t15576*(2.0/3.0);
                double t15590 = pow(2.0,1.0/3.0);
                double t15591 = t15590*2.0;
                double t15592 = t15591-2.0;
                double t15593 = two_13*2.0;
                double t15594 = t15593-2.0;
                double t15595 = 1.0/t15594;
                double t15596 = log(t15552);
                v_rho_b_rho_b[Q] += scale * (t15526*(t15583+t15587+t15589-t15592*t15595*(t15583+t15587+t15589-t15562*t15574*t15557*t15584*2.0+t15553*1.0/(t15537*t15537*t15537)*t15548*t15569*2.0-t15553*t15548*t15575*(b2f*c*t15538*(4.0/9.0)+b4f*t15533*t15539*(1.0E1/9.0)+b1f*c*t15540*t15538*(2.0/9.0)+b3f*c*t15529*t15538*(2.0/3.0)-b1f*t15533*t15546*t15539*(1.0/3.6E1)+b3f*t15540*t15533*t15539*(1.0/1.2E1))-1.0/(t15552*t15552)*1.0/(t15537*t15537*t15537*t15537)*t15548*t15549*t15569*(1.0/2.0)+a1f*c*c0f*t15538*t15596*(8.0/9.0)-a1p*c*c0p*t15538*t15588*(8.0/9.0)-a1f*c*t15553*t15554*t15555*t15575*(2.0/3.0))-t15562*t15574*t15557*t15584*2.0-a1p*c*c0p*t15538*t15588*(8.0/9.0))-t15592*t15595*(t15553*t15555*t15548*t15575-t15562*t15564*t15557*t15576-a1f*c*c0f*t15554*t15596*(2.0/3.0)+a1p*c*c0p*t15554*t15588*(2.0/3.0))*2.0-t15562*t15564*t15557*t15576*2.0+a1p*c*c0p*t15554*t15588*(4.0/3.0));
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t14259 = rho_a+rho_b;
                double t14260 = 1.0/pow(t14259,1.0/3.0);
                double t14261 = c*t14260;
                double t14262 = sqrt(t14261);
                double t14263 = pow(t14261,3.0/2.0);
                double t14264 = c*c;
                double t14265 = 1.0/pow(t14259,2.0/3.0);
                double t14266 = 1.0/c0p;
                double t14267 = b1p*t14262;
                double t14268 = b3p*t14263;
                double t14269 = b4p*t14264*t14265;
                double t14270 = b2p*c*t14260;
                double t14271 = t14270+t14267+t14268+t14269;
                double t14272 = 1.0/t14271;
                double t14273 = t14272*t14266*(1.0/2.0);
                double t14274 = t14273+1.0;
                double t14275 = log(t14274);
                double t14276 = a1p*c*t14260;
                double t14277 = t14276+1.0;
                double t14278 = c0p*t14275*t14277*2.0;
                double t14279 = rho_a-rho_b;
                double t14280 = t14279*t14279;
                double t14281 = 1.0/t14259;
                double t14282 = t14281*t14279;
                double t14283 = two_13*2.0;
                double t14284 = t14283-2.0;
                double t14285 = 1.0/t14284;
                double t14286 = 1.0/(t14259*t14259*t14259*t14259);
                double t14287 = t14280*t14280;
                double t14288 = t14282+1.0;
                double t14289 = pow(t14288,4.0/3.0);
                double t14290 = -t14282+1.0;
                double t14291 = pow(t14290,4.0/3.0);
                double t14292 = t14291+t14289-2.0;
                v[Q] += scale * (-t14259*(t14278-t14292*t14285*t14286*t14287*(t14278-c0f*log((1.0/2.0)/(c0f*(b1f*t14262+b3f*t14263+b2f*c*t14260+b4f*t14264*t14265))+1.0)*(a1f*c*t14260+1.0)*2.0)+(Aa*t14292*t14285*log((1.0/2.0)/(Aa*(b1a*t14262+b3a*t14263+b2a*c*t14260+b4a*t14264*t14265))+1.0)*(t14286*t14287-1.0)*(a1a*c*t14260+1.0)*2.0)/d2fz0));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t14294 = rho_a+rho_b;
                double t14295 = 1.0/pow(t14294,1.0/3.0);
                double t14296 = c*t14295;
                double t14297 = sqrt(t14296);
                double t14298 = b1p*t14297;
                double t14299 = pow(t14296,3.0/2.0);
                double t14300 = b3p*t14299;
                double t14301 = c*c;
                double t14302 = 1.0/pow(t14294,2.0/3.0);
                double t14303 = b4p*t14301*t14302;
                double t14304 = b2p*c*t14295;
                double t14305 = t14300+t14303+t14304+t14298;
                double t14306 = 1.0/pow(t14294,4.0/3.0);
                double t14307 = 1.0/c0p;
                double t14308 = 1.0/t14305;
                double t14309 = t14307*t14308*(1.0/2.0);
                double t14310 = t14309+1.0;
                double t14311 = a1p*c*t14295;
                double t14312 = t14311+1.0;
                double t14313 = rho_a-rho_b;
                double t14314 = t14313*t14313;
                double t14315 = 1.0/t14294;
                double t14316 = t14313*t14315;
                double t14317 = two_13*2.0;
                double t14318 = t14317-2.0;
                double t14319 = 1.0/t14318;
                double t14320 = 1.0/c0f;
                double t14321 = b1f*t14297;
                double t14322 = b3f*t14299;
                double t14323 = b4f*t14301*t14302;
                double t14324 = b2f*c*t14295;
                double t14325 = t14321+t14322+t14323+t14324;
                double t14326 = 1.0/t14325;
                double t14327 = t14320*t14326*(1.0/2.0);
                double t14328 = t14327+1.0;
                double t14329 = log(t14328);
                double t14330 = a1f*c*t14295;
                double t14331 = t14330+1.0;
                double t14332 = log(t14310);
                double t14342 = c0f*t14331*t14329*2.0;
                double t14343 = c0p*t14312*t14332*2.0;
                double t14333 = t14342-t14343;
                double t14334 = t14316+1.0;
                double t14335 = pow(t14334,4.0/3.0);
                double t14336 = -t14316+1.0;
                double t14337 = pow(t14336,4.0/3.0);
                double t14338 = t14335+t14337-2.0;
                double t14339 = 1.0/(t14294*t14294);
                double t14370 = t14313*t14339;
                double t14340 = t14315-t14370;
                double t14341 = 1.0/(t14294*t14294*t14294*t14294);
                double t14344 = t14314*t14314;
                double t14345 = 1.0/pow(t14294,5.0/3.0);
                double t14346 = 1.0/sqrt(t14296);
                double t14347 = 1.0/t14310;
                double t14348 = 1.0/(t14305*t14305);
                double t14349 = b4p*t14301*t14345*(2.0/3.0);
                double t14350 = b2p*c*t14306*(1.0/3.0);
                double t14351 = b1p*c*t14306*t14346*(1.0/6.0);
                double t14352 = b3p*c*t14306*t14297*(1.0/2.0);
                double t14353 = t14350+t14351+t14352+t14349;
                double t14354 = t14312*t14353*t14347*t14348;
                double t14355 = 1.0/(t14294*t14294*t14294*t14294*t14294);
                double t14356 = 1.0/d2fz0;
                double t14357 = 1.0/Aa;
                double t14358 = b1a*t14297;
                double t14359 = b3a*t14299;
                double t14360 = b4a*t14301*t14302;
                double t14361 = b2a*c*t14295;
                double t14362 = t14360+t14361+t14358+t14359;
                double t14363 = 1.0/t14362;
                double t14364 = t14363*t14357*(1.0/2.0);
                double t14365 = t14364+1.0;
                double t14366 = log(t14365);
                double t14367 = a1a*c*t14295;
                double t14368 = t14367+1.0;
                double t14369 = pow(t14334,1.0/3.0);
                double t14371 = t14340*t14369*(4.0/3.0);
                double t14372 = pow(t14336,1.0/3.0);
                double t14373 = t14371-t14340*t14372*(4.0/3.0);
                double t14374 = t14341*t14344;
                double t14375 = t14374-1.0;
                v_rho_a[Q] += scale * (-t14343-t14294*(t14354-t14341*t14344*t14319*t14338*(t14354-(t14331*1.0/(t14325*t14325)*(b2f*c*t14306*(1.0/3.0)+b4f*t14301*t14345*(2.0/3.0)+b1f*c*t14306*t14346*(1.0/6.0)+b3f*c*t14306*t14297*(1.0/2.0)))/t14328+a1f*c*c0f*t14306*t14329*(2.0/3.0)-a1p*c*c0p*t14332*t14306*(2.0/3.0))-a1p*c*c0p*t14332*t14306*(2.0/3.0)+t14341*t14333*t14344*t14319*t14373-t14333*t14344*t14319*t14355*t14338*4.0+t14313*t14314*t14341*t14333*t14319*t14338*4.0-Aa*t14319*t14338*t14356*t14366*t14368*(t14344*t14355*4.0-t14313*t14314*t14341*4.0)*2.0+Aa*t14319*t14373*t14356*t14366*t14375*t14368*2.0+(1.0/(t14362*t14362)*t14319*t14338*t14356*t14375*t14368*(b2a*c*t14306*(1.0/3.0)+b4a*t14301*t14345*(2.0/3.0)+b1a*c*t14306*t14346*(1.0/6.0)+b3a*c*t14306*t14297*(1.0/2.0)))/t14365-Aa*a1a*c*t14306*t14319*t14338*t14356*t14366*t14375*(2.0/3.0))-t14341*t14333*t14344*t14319*t14338-Aa*t14319*t14338*t14356*t14366*t14375*t14368*2.0);
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t14377 = rho_a+rho_b;
                double t14378 = 1.0/pow(t14377,1.0/3.0);
                double t14379 = c*t14378;
                double t14380 = sqrt(t14379);
                double t14381 = b1p*t14380;
                double t14382 = pow(t14379,3.0/2.0);
                double t14383 = b3p*t14382;
                double t14384 = c*c;
                double t14385 = 1.0/pow(t14377,2.0/3.0);
                double t14386 = b4p*t14384*t14385;
                double t14387 = b2p*c*t14378;
                double t14388 = t14381+t14383+t14386+t14387;
                double t14389 = 1.0/pow(t14377,4.0/3.0);
                double t14390 = 1.0/c0p;
                double t14391 = 1.0/t14388;
                double t14392 = t14390*t14391*(1.0/2.0);
                double t14393 = t14392+1.0;
                double t14394 = a1p*c*t14378;
                double t14395 = t14394+1.0;
                double t14396 = rho_a-rho_b;
                double t14397 = t14396*t14396;
                double t14398 = 1.0/t14377;
                double t14399 = t14396*t14398;
                double t14400 = two_13*2.0;
                double t14401 = t14400-2.0;
                double t14402 = 1.0/t14401;
                double t14403 = 1.0/c0f;
                double t14404 = b1f*t14380;
                double t14405 = b3f*t14382;
                double t14406 = b4f*t14384*t14385;
                double t14407 = b2f*c*t14378;
                double t14408 = t14404+t14405+t14406+t14407;
                double t14409 = 1.0/t14408;
                double t14410 = t14403*t14409*(1.0/2.0);
                double t14411 = t14410+1.0;
                double t14412 = log(t14411);
                double t14413 = a1f*c*t14378;
                double t14414 = t14413+1.0;
                double t14415 = log(t14393);
                double t14426 = c0f*t14412*t14414*2.0;
                double t14427 = c0p*t14415*t14395*2.0;
                double t14416 = t14426-t14427;
                double t14417 = t14399+1.0;
                double t14418 = pow(t14417,4.0/3.0);
                double t14419 = -t14399+1.0;
                double t14420 = pow(t14419,4.0/3.0);
                double t14421 = t14420+t14418-2.0;
                double t14422 = 1.0/(t14377*t14377);
                double t14423 = t14422*t14396;
                double t14424 = t14423+t14398;
                double t14425 = 1.0/(t14377*t14377*t14377*t14377);
                double t14428 = t14397*t14397;
                double t14429 = 1.0/pow(t14377,5.0/3.0);
                double t14430 = 1.0/sqrt(t14379);
                double t14431 = 1.0/t14393;
                double t14432 = 1.0/(t14388*t14388);
                double t14433 = b4p*t14384*t14429*(2.0/3.0);
                double t14434 = b2p*c*t14389*(1.0/3.0);
                double t14435 = b1p*c*t14430*t14389*(1.0/6.0);
                double t14436 = b3p*c*t14380*t14389*(1.0/2.0);
                double t14437 = t14433+t14434+t14435+t14436;
                double t14438 = 1.0/(t14377*t14377*t14377*t14377*t14377);
                double t14439 = 1.0/d2fz0;
                double t14440 = 1.0/Aa;
                double t14441 = b1a*t14380;
                double t14442 = b3a*t14382;
                double t14443 = b4a*t14384*t14385;
                double t14444 = b2a*c*t14378;
                double t14445 = t14441+t14442+t14443+t14444;
                double t14446 = 1.0/t14445;
                double t14447 = t14440*t14446*(1.0/2.0);
                double t14448 = t14447+1.0;
                double t14449 = log(t14448);
                double t14450 = a1a*c*t14378;
                double t14451 = t14450+1.0;
                double t14452 = pow(t14417,1.0/3.0);
                double t14453 = t14424*t14452*(4.0/3.0);
                double t14454 = pow(t14419,1.0/3.0);
                double t14455 = t14453-t14424*t14454*(4.0/3.0);
                double t14456 = t14425*t14428;
                double t14457 = t14456-1.0;
                v_rho_b[Q] += scale * (-t14427+t14377*(-t14431*t14432*t14437*t14395+a1p*c*c0p*t14415*t14389*(2.0/3.0)+t14402*t14421*t14416*t14428*t14438*4.0+t14402*t14416*t14425*t14428*t14455+t14402*t14421*t14425*t14428*(t14431*t14432*t14437*t14395-(t14414*1.0/(t14408*t14408)*(b2f*c*t14389*(1.0/3.0)+b4f*t14384*t14429*(2.0/3.0)+b1f*c*t14430*t14389*(1.0/6.0)+b3f*c*t14380*t14389*(1.0/2.0)))/t14411+a1f*c*c0f*t14412*t14389*(2.0/3.0)-a1p*c*c0p*t14415*t14389*(2.0/3.0))+t14402*t14421*t14416*t14425*t14396*t14397*4.0+Aa*t14402*t14421*t14451*t14439*t14449*(t14428*t14438*4.0+t14425*t14396*t14397*4.0)*2.0+Aa*t14402*t14451*t14455*t14439*t14457*t14449*2.0-(t14402*t14421*t14451*1.0/(t14445*t14445)*t14439*t14457*(b2a*c*t14389*(1.0/3.0)+b4a*t14384*t14429*(2.0/3.0)+b1a*c*t14430*t14389*(1.0/6.0)+b3a*c*t14380*t14389*(1.0/2.0)))/t14448+Aa*a1a*c*t14402*t14421*t14439*t14457*t14449*t14389*(2.0/3.0))-t14402*t14421*t14416*t14425*t14428-Aa*t14402*t14421*t14451*t14439*t14457*t14449*2.0);
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t14464 = rho_a+rho_b;
                double t14465 = 1.0/pow(t14464,1.0/3.0);
                double t14466 = c*t14465;
                double t14467 = sqrt(t14466);
                double t14468 = b1p*t14467;
                double t14469 = pow(t14466,3.0/2.0);
                double t14470 = b3p*t14469;
                double t14471 = c*c;
                double t14472 = 1.0/pow(t14464,2.0/3.0);
                double t14473 = b4p*t14471*t14472;
                double t14474 = b2p*c*t14465;
                double t14475 = t14470+t14473+t14474+t14468;
                double t14476 = 1.0/pow(t14464,7.0/3.0);
                double t14477 = 1.0/pow(t14464,8.0/3.0);
                double t14478 = 1.0/sqrt(t14466);
                double t14479 = a1p*c*t14465;
                double t14480 = t14479+1.0;
                double t14481 = 1.0/c0p;
                double t14482 = 1.0/t14475;
                double t14483 = t14481*t14482*(1.0/2.0);
                double t14484 = t14483+1.0;
                double t14485 = 1.0/t14484;
                double t14486 = 1.0/pow(t14464,4.0/3.0);
                double t14526 = 1.0/pow(t14464,5.0/3.0);
                double t14528 = b4p*t14471*t14526*(2.0/3.0);
                double t14529 = b2p*c*t14486*(1.0/3.0);
                double t14530 = b1p*c*t14486*t14478*(1.0/6.0);
                double t14531 = b3p*c*t14467*t14486*(1.0/2.0);
                double t14487 = t14530+t14531+t14528+t14529;
                double t14488 = rho_a-rho_b;
                double t14489 = 1.0/t14464;
                double t14490 = t14488*t14489;
                double t14491 = two_13*2.0;
                double t14492 = t14491-2.0;
                double t14493 = 1.0/t14492;
                double t14494 = 1.0/c0f;
                double t14495 = b1f*t14467;
                double t14496 = b3f*t14469;
                double t14497 = b4f*t14471*t14472;
                double t14498 = b2f*c*t14465;
                double t14499 = t14495+t14496+t14497+t14498;
                double t14500 = 1.0/t14499;
                double t14501 = t14500*t14494*(1.0/2.0);
                double t14502 = t14501+1.0;
                double t14503 = log(t14502);
                double t14504 = a1f*c*t14465;
                double t14505 = t14504+1.0;
                double t14506 = log(t14484);
                double t14514 = c0f*t14503*t14505*2.0;
                double t14515 = c0p*t14506*t14480*2.0;
                double t14507 = t14514-t14515;
                double t14508 = t14488*t14488;
                double t14509 = t14490+1.0;
                double t14510 = pow(t14509,4.0/3.0);
                double t14511 = -t14490+1.0;
                double t14512 = pow(t14511,4.0/3.0);
                double t14513 = t14510+t14512-2.0;
                double t14516 = 1.0/(t14464*t14464);
                double t14520 = t14516*t14488;
                double t14517 = -t14520+t14489;
                double t14518 = 1.0/(t14464*t14464*t14464*t14464);
                double t14519 = pow(t14509,1.0/3.0);
                double t14521 = t14517*t14519*(4.0/3.0);
                double t14522 = pow(t14511,1.0/3.0);
                double t14555 = t14522*t14517*(4.0/3.0);
                double t14523 = t14521-t14555;
                double t14524 = 1.0/(t14464*t14464*t14464*t14464*t14464);
                double t14525 = t14508*t14508;
                double t14527 = 1.0/(t14475*t14475);
                double t14532 = 1.0/t14502;
                double t14533 = 1.0/(t14499*t14499);
                double t14534 = b4f*t14471*t14526*(2.0/3.0);
                double t14535 = b2f*c*t14486*(1.0/3.0);
                double t14536 = b1f*c*t14486*t14478*(1.0/6.0);
                double t14537 = b3f*c*t14467*t14486*(1.0/2.0);
                double t14538 = t14534+t14535+t14536+t14537;
                double t14539 = t14505*t14532*t14533*t14538;
                double t14540 = a1p*c*c0p*t14506*t14486*(2.0/3.0);
                double t14556 = t14480*t14527*t14485*t14487;
                double t14557 = a1f*c*c0f*t14503*t14486*(2.0/3.0);
                double t14541 = t14540-t14556+t14539-t14557;
                double t14542 = t14487*t14487;
                double t14543 = 1.0/pow(t14466,3.0/2.0);
                double t14544 = b4p*t14471*t14477*(1.0E1/9.0);
                double t14545 = b2p*c*t14476*(4.0/9.0);
                double t14546 = b1p*c*t14476*t14478*(2.0/9.0);
                double t14547 = b3p*c*t14467*t14476*(2.0/3.0);
                double t14548 = b3p*t14471*t14477*t14478*(1.0/1.2E1);
                double t14549 = t14544+t14545+t14546+t14547+t14548-b1p*t14471*t14543*t14477*(1.0/3.6E1);
                double t14550 = 1.0/(t14475*t14475*t14475);
                double t14551 = t14550*t14542*t14480*t14485*2.0;
                double t14552 = t14538*t14538;
                double t14553 = 1.0/(t14484*t14484);
                double t14554 = 1.0/(t14475*t14475*t14475*t14475);
                double t14558 = t14516*2.0;
                double t14559 = 1.0/(t14464*t14464*t14464);
                double t14577 = t14559*t14488*2.0;
                double t14560 = t14558-t14577;
                double t14561 = t14517*t14517;
                double t14562 = a1p*c*c0p*t14506*t14476*(8.0/9.0);
                double t14563 = 1.0/(t14464*t14464*t14464*t14464*t14464*t14464);
                double t14564 = 1.0/d2fz0;
                double t14565 = 1.0/Aa;
                double t14566 = b1a*t14467;
                double t14567 = b3a*t14469;
                double t14568 = b4a*t14471*t14472;
                double t14569 = b2a*c*t14465;
                double t14570 = t14566+t14567+t14568+t14569;
                double t14571 = 1.0/t14570;
                double t14572 = t14571*t14565*(1.0/2.0);
                double t14573 = t14572+1.0;
                double t14574 = log(t14573);
                double t14575 = a1a*c*t14465;
                double t14576 = t14575+1.0;
                double t14578 = t14522*t14560*(4.0/3.0);
                double t14579 = 1.0/pow(t14509,2.0/3.0);
                double t14580 = t14561*t14579*(4.0/9.0);
                double t14581 = 1.0/pow(t14511,2.0/3.0);
                double t14582 = t14561*t14581*(4.0/9.0);
                double t14583 = t14580+t14582+t14578-t14560*t14519*(4.0/3.0);
                double t14584 = t14524*t14525*4.0;
                double t14596 = t14508*t14518*t14488*4.0;
                double t14585 = t14584-t14596;
                double t14586 = 1.0/t14573;
                double t14587 = t14525*t14518;
                double t14588 = t14587-1.0;
                double t14589 = 1.0/(t14570*t14570);
                double t14590 = b4a*t14471*t14526*(2.0/3.0);
                double t14591 = b2a*c*t14486*(1.0/3.0);
                double t14592 = b1a*c*t14486*t14478*(1.0/6.0);
                double t14593 = b3a*c*t14467*t14486*(1.0/2.0);
                double t14594 = t14590+t14591+t14592+t14593;
                double t14595 = t14594*t14594;
                v_rho_a_rho_a[Q] += scale * (-t14464*(t14551+t14562-t14480*t14527*t14485*t14549-t14513*t14525*t14518*t14493*(t14551+t14562-t14480*t14527*t14485*t14549-t14505*t14532*t14552*1.0/(t14499*t14499*t14499)*2.0+t14505*t14532*t14533*(b2f*c*t14476*(4.0/9.0)+b4f*t14471*t14477*(1.0E1/9.0)+b1f*c*t14476*t14478*(2.0/9.0)+b3f*c*t14467*t14476*(2.0/3.0)-b1f*t14471*t14543*t14477*(1.0/3.6E1)+b3f*t14471*t14477*t14478*(1.0/1.2E1))+1.0/(t14502*t14502)*t14505*t14552*t14494*1.0/(t14499*t14499*t14499*t14499)*(1.0/2.0)-a1f*c*c0f*t14503*t14476*(8.0/9.0)-t14542*t14480*t14481*t14553*t14554*(1.0/2.0)+a1f*c*t14532*t14533*t14538*t14486*(2.0/3.0)-a1p*c*t14527*t14485*t14486*t14487*(2.0/3.0))-t14513*t14541*t14524*t14525*t14493*8.0-t14523*t14524*t14507*t14525*t14493*8.0+t14523*t14541*t14525*t14518*t14493*2.0+t14513*t14507*t14525*t14563*t14493*2.0E1-t14542*t14480*t14481*t14553*t14554*(1.0/2.0)+t14513*t14507*t14508*t14518*t14493*1.2E1+t14507*t14525*t14518*t14493*t14583-t14513*t14524*t14507*t14508*t14493*t14488*3.2E1+t14513*t14541*t14508*t14518*t14493*t14488*8.0+t14523*t14507*t14508*t14518*t14493*t14488*8.0-a1p*c*t14527*t14485*t14486*t14487*(2.0/3.0)-Aa*t14523*t14564*t14493*t14574*t14576*t14585*4.0+Aa*t14564*t14493*t14574*t14583*t14576*t14588*2.0+Aa*t14513*t14564*t14493*t14574*t14576*(t14525*t14563*2.0E1+t14508*t14518*1.2E1-t14524*t14508*t14488*3.2E1)*2.0-t14513*t14564*t14493*t14576*t14585*t14594*t14586*t14589*2.0+t14523*t14564*t14493*t14576*t14594*t14586*t14588*t14589*2.0+t14513*1.0/(t14570*t14570*t14570)*t14564*t14493*t14576*t14586*t14595*t14588*2.0-t14513*t14564*t14493*t14576*t14586*t14588*t14589*(b2a*c*t14476*(4.0/9.0)+b4a*t14471*t14477*(1.0E1/9.0)+b1a*c*t14476*t14478*(2.0/9.0)+b3a*c*t14467*t14476*(2.0/3.0)-b1a*t14471*t14543*t14477*(1.0/3.6E1)+b3a*t14471*t14477*t14478*(1.0/1.2E1))+Aa*a1a*c*t14513*t14564*t14493*t14574*t14486*t14585*(4.0/3.0)+Aa*a1a*c*t14513*t14564*t14493*t14574*t14476*t14588*(8.0/9.0)-Aa*a1a*c*t14523*t14564*t14493*t14574*t14486*t14588*(4.0/3.0)-t14513*1.0/(t14570*t14570*t14570*t14570)*t14564*1.0/(t14573*t14573)*t14493*t14565*t14576*t14595*t14588*(1.0/2.0)-a1a*c*t14513*t14564*t14493*t14486*t14594*t14586*t14588*t14589*(2.0/3.0))-t14480*t14527*t14485*t14487*2.0+a1p*c*c0p*t14506*t14486*(4.0/3.0)+t14513*t14524*t14507*t14525*t14493*8.0-t14513*t14541*t14525*t14518*t14493*2.0-t14523*t14507*t14525*t14518*t14493*2.0-t14513*t14507*t14508*t14518*t14493*t14488*8.0-Aa*t14523*t14564*t14493*t14574*t14576*t14588*4.0+Aa*t14513*t14564*t14493*t14574*t14576*(t14584-t14596)*4.0-t14513*t14564*t14493*t14576*t14594*t14586*t14588*t14589*2.0+Aa*a1a*c*t14513*t14564*t14493*t14574*t14486*t14588*(4.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t14598 = rho_a+rho_b;
                double t14599 = 1.0/pow(t14598,1.0/3.0);
                double t14600 = c*t14599;
                double t14601 = sqrt(t14600);
                double t14602 = b1p*t14601;
                double t14603 = pow(t14600,3.0/2.0);
                double t14604 = b3p*t14603;
                double t14605 = c*c;
                double t14606 = 1.0/pow(t14598,2.0/3.0);
                double t14607 = b4p*t14605*t14606;
                double t14608 = b2p*c*t14599;
                double t14609 = t14602+t14604+t14607+t14608;
                double t14610 = 1.0/pow(t14598,7.0/3.0);
                double t14611 = 1.0/pow(t14598,8.0/3.0);
                double t14612 = 1.0/sqrt(t14600);
                double t14613 = a1p*c*t14599;
                double t14614 = t14613+1.0;
                double t14615 = 1.0/c0p;
                double t14616 = 1.0/t14609;
                double t14617 = t14615*t14616*(1.0/2.0);
                double t14618 = t14617+1.0;
                double t14619 = 1.0/t14618;
                double t14620 = 1.0/pow(t14598,4.0/3.0);
                double t14663 = 1.0/pow(t14598,5.0/3.0);
                double t14665 = b4p*t14605*t14663*(2.0/3.0);
                double t14666 = b2p*c*t14620*(1.0/3.0);
                double t14667 = b1p*c*t14620*t14612*(1.0/6.0);
                double t14668 = b3p*c*t14601*t14620*(1.0/2.0);
                double t14621 = t14665+t14666+t14667+t14668;
                double t14622 = rho_a-rho_b;
                double t14623 = 1.0/t14598;
                double t14624 = t14622*t14623;
                double t14625 = two_13*2.0;
                double t14626 = t14625-2.0;
                double t14627 = 1.0/t14626;
                double t14628 = 1.0/c0f;
                double t14629 = b1f*t14601;
                double t14630 = b3f*t14603;
                double t14631 = b4f*t14605*t14606;
                double t14632 = b2f*c*t14599;
                double t14633 = t14630+t14631+t14632+t14629;
                double t14634 = 1.0/t14633;
                double t14635 = t14634*t14628*(1.0/2.0);
                double t14636 = t14635+1.0;
                double t14637 = log(t14636);
                double t14638 = a1f*c*t14599;
                double t14639 = t14638+1.0;
                double t14640 = log(t14618);
                double t14652 = c0f*t14637*t14639*2.0;
                double t14653 = c0p*t14640*t14614*2.0;
                double t14641 = t14652-t14653;
                double t14642 = t14622*t14622;
                double t14643 = t14624+1.0;
                double t14644 = pow(t14643,4.0/3.0);
                double t14645 = -t14624+1.0;
                double t14646 = pow(t14645,4.0/3.0);
                double t14647 = t14644+t14646-2.0;
                double t14648 = 1.0/(t14598*t14598);
                double t14649 = t14622*t14648;
                double t14650 = t14623+t14649;
                double t14651 = 1.0/(t14598*t14598*t14598*t14598);
                double t14654 = pow(t14643,1.0/3.0);
                double t14655 = t14650*t14654*(4.0/3.0);
                double t14656 = pow(t14645,1.0/3.0);
                double t14690 = t14650*t14656*(4.0/3.0);
                double t14657 = -t14690+t14655;
                double t14658 = t14642*t14642;
                double t14659 = t14623-t14649;
                double t14660 = t14654*t14659*(4.0/3.0);
                double t14694 = t14656*t14659*(4.0/3.0);
                double t14661 = t14660-t14694;
                double t14662 = 1.0/(t14598*t14598*t14598*t14598*t14598);
                double t14664 = 1.0/(t14609*t14609);
                double t14669 = t14621*t14621;
                double t14670 = 1.0/t14636;
                double t14671 = 1.0/(t14633*t14633);
                double t14672 = 1.0/pow(t14600,3.0/2.0);
                double t14673 = b4p*t14611*t14605*(1.0E1/9.0);
                double t14674 = b2p*c*t14610*(4.0/9.0);
                double t14675 = b1p*c*t14610*t14612*(2.0/9.0);
                double t14676 = b3p*c*t14601*t14610*(2.0/3.0);
                double t14677 = b3p*t14611*t14612*t14605*(1.0/1.2E1);
                double t14678 = t14673+t14674+t14675+t14676+t14677-b1p*t14611*t14605*t14672*(1.0/3.6E1);
                double t14679 = t14614*t14619*t14664*t14678;
                double t14680 = b4f*t14605*t14663*(2.0/3.0);
                double t14681 = b2f*c*t14620*(1.0/3.0);
                double t14682 = b1f*c*t14620*t14612*(1.0/6.0);
                double t14683 = b3f*c*t14601*t14620*(1.0/2.0);
                double t14684 = t14680+t14681+t14682+t14683;
                double t14685 = 1.0/(t14609*t14609*t14609);
                double t14686 = t14684*t14684;
                double t14687 = 1.0/(t14618*t14618);
                double t14688 = 1.0/(t14609*t14609*t14609*t14609);
                double t14689 = t14614*t14615*t14669*t14687*t14688*(1.0/2.0);
                double t14691 = t14621*t14614*t14619*t14664;
                double t14692 = a1f*c*c0f*t14620*t14637*(2.0/3.0);
                double t14695 = t14670*t14671*t14639*t14684;
                double t14696 = a1p*c*c0p*t14620*t14640*(2.0/3.0);
                double t14693 = t14691+t14692-t14695-t14696;
                double t14697 = 1.0/(t14598*t14598*t14598);
                double t14698 = a1p*c*t14620*t14621*t14619*t14664*(2.0/3.0);
                double t14699 = 1.0/(t14598*t14598*t14598*t14598*t14598*t14598);
                double t14700 = 1.0/d2fz0;
                double t14701 = 1.0/Aa;
                double t14702 = b1a*t14601;
                double t14703 = b3a*t14603;
                double t14704 = b4a*t14605*t14606;
                double t14705 = b2a*c*t14599;
                double t14706 = t14702+t14703+t14704+t14705;
                double t14707 = 1.0/t14706;
                double t14708 = t14701*t14707*(1.0/2.0);
                double t14709 = t14708+1.0;
                double t14710 = log(t14709);
                double t14711 = a1a*c*t14599;
                double t14712 = t14711+1.0;
                double t14713 = t14662*t14658*4.0;
                double t14714 = t14622*t14656*t14697*(8.0/3.0);
                double t14715 = 1.0/pow(t14643,2.0/3.0);
                double t14716 = t14650*t14715*t14659*(4.0/9.0);
                double t14717 = 1.0/pow(t14645,2.0/3.0);
                double t14718 = t14650*t14717*t14659*(4.0/9.0);
                double t14719 = t14714+t14716+t14718-t14622*t14654*t14697*(8.0/3.0);
                double t14720 = t14622*t14642*t14651*4.0;
                double t14721 = 1.0/t14709;
                double t14722 = t14720+t14713;
                double t14723 = 1.0/(t14706*t14706);
                double t14724 = b4a*t14605*t14663*(2.0/3.0);
                double t14725 = b2a*c*t14620*(1.0/3.0);
                double t14726 = b1a*c*t14620*t14612*(1.0/6.0);
                double t14727 = b3a*c*t14601*t14620*(1.0/2.0);
                double t14728 = t14724+t14725+t14726+t14727;
                double t14729 = t14651*t14658;
                double t14730 = t14729-1.0;
                double t14731 = t14728*t14728;
                double t14732 = t14720-t14713;
                v_rho_a_rho_b[Q] += scale * (t14598*(t14679+t14689+t14698-t14614*t14619*t14685*t14669*2.0+t14651*t14661*t14627*t14658*(t14691+t14692-t14695-t14696)-a1p*c*c0p*t14610*t14640*(8.0/9.0)-t14651*t14627*t14647*t14658*(t14679+t14689+t14698-t14614*t14619*t14685*t14669*2.0+1.0/(t14633*t14633*t14633)*t14670*t14639*t14686*2.0-t14670*t14671*t14639*(b2f*c*t14610*(4.0/9.0)+b4f*t14611*t14605*(1.0E1/9.0)+b1f*c*t14610*t14612*(2.0/9.0)+b3f*c*t14601*t14610*(2.0/3.0)-b1f*t14611*t14605*t14672*(1.0/3.6E1)+b3f*t14611*t14612*t14605*(1.0/1.2E1))-1.0/(t14633*t14633*t14633*t14633)*1.0/(t14636*t14636)*t14628*t14639*t14686*(1.0/2.0)+a1f*c*c0f*t14610*t14637*(8.0/9.0)-a1p*c*c0p*t14610*t14640*(8.0/9.0)-a1f*c*t14620*t14670*t14671*t14684*(2.0/3.0))+t14641*t14642*t14651*t14627*t14647*1.2E1+t14641*t14661*t14662*t14627*t14658*4.0+t14641*t14651*t14627*t14719*t14658-t14641*t14662*t14627*t14657*t14658*4.0-t14651*t14627*t14657*t14693*t14658-t14662*t14627*t14647*t14693*t14658*8.0-t14641*t14627*t14647*t14658*t14699*2.0E1+t14622*t14641*t14642*t14651*t14661*t14627*4.0+t14622*t14641*t14642*t14651*t14627*t14657*4.0+Aa*t14700*t14710*t14712*t14722*t14661*t14627*2.0+Aa*t14700*t14710*t14712*t14730*t14627*t14719*2.0-Aa*t14700*t14710*t14712*t14627*t14657*(t14713-t14622*t14642*t14651*4.0)*2.0+Aa*t14700*t14710*t14712*t14627*t14647*(t14642*t14651*1.2E1-t14658*t14699*2.0E1)*2.0-t14700*t14712*t14721*t14730*t14723*t14661*t14627*t14728+t14700*t14712*t14721*t14730*t14723*t14627*t14728*t14657+t14700*t14712*t14721*t14722*t14723*t14627*t14647*t14728-t14700*t14712*t14721*t14723*t14732*t14627*t14647*t14728-t14700*t14712*t14721*t14730*t14731*1.0/(t14706*t14706*t14706)*t14627*t14647*2.0+t14700*t14712*t14721*t14730*t14723*t14627*t14647*(b2a*c*t14610*(4.0/9.0)+b4a*t14611*t14605*(1.0E1/9.0)+b1a*c*t14610*t14612*(2.0/9.0)+b3a*c*t14601*t14610*(2.0/3.0)-b1a*t14611*t14605*t14672*(1.0/3.6E1)+b3a*t14611*t14612*t14605*(1.0/1.2E1))+Aa*a1a*c*t14700*t14620*t14710*t14730*t14661*t14627*(2.0/3.0)-Aa*a1a*c*t14610*t14700*t14710*t14730*t14627*t14647*(8.0/9.0)-Aa*a1a*c*t14700*t14620*t14710*t14730*t14627*t14657*(2.0/3.0)-Aa*a1a*c*t14700*t14620*t14710*t14722*t14627*t14647*(2.0/3.0)+Aa*a1a*c*t14700*t14620*t14710*t14732*t14627*t14647*(2.0/3.0)+t14700*t14701*t14712*t14730*t14731*1.0/(t14706*t14706*t14706*t14706)*t14627*1.0/(t14709*t14709)*t14647*(1.0/2.0)+a1a*c*t14700*t14620*t14721*t14730*t14723*t14627*t14647*t14728*(2.0/3.0))-t14621*t14614*t14619*t14664*2.0+t14651*t14627*t14647*t14658*(t14691+t14692-t14695-t14696)*2.0+a1p*c*c0p*t14620*t14640*(4.0/3.0)-t14641*t14651*t14661*t14627*t14658+t14641*t14651*t14627*t14657*t14658+t14641*t14662*t14627*t14647*t14658*8.0-Aa*t14700*t14710*t14712*t14730*t14661*t14627*2.0+Aa*t14700*t14710*t14712*t14730*t14627*t14657*2.0+Aa*t14700*t14710*t14712*t14722*t14627*t14647*2.0-Aa*t14700*t14710*t14712*t14732*t14627*t14647*2.0-t14700*t14712*t14721*t14730*t14723*t14627*t14647*t14728*2.0+Aa*a1a*c*t14700*t14620*t14710*t14730*t14627*t14647*(4.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t14734 = rho_a+rho_b;
                double t14735 = 1.0/pow(t14734,1.0/3.0);
                double t14736 = c*t14735;
                double t14737 = sqrt(t14736);
                double t14738 = b1p*t14737;
                double t14739 = pow(t14736,3.0/2.0);
                double t14740 = b3p*t14739;
                double t14741 = c*c;
                double t14742 = 1.0/pow(t14734,2.0/3.0);
                double t14743 = b4p*t14741*t14742;
                double t14744 = b2p*c*t14735;
                double t14745 = t14740+t14743+t14744+t14738;
                double t14746 = 1.0/pow(t14734,7.0/3.0);
                double t14747 = 1.0/pow(t14734,8.0/3.0);
                double t14748 = 1.0/sqrt(t14736);
                double t14749 = a1p*c*t14735;
                double t14750 = t14749+1.0;
                double t14751 = 1.0/c0p;
                double t14752 = 1.0/t14745;
                double t14753 = t14751*t14752*(1.0/2.0);
                double t14754 = t14753+1.0;
                double t14755 = 1.0/t14754;
                double t14756 = 1.0/pow(t14734,4.0/3.0);
                double t14796 = 1.0/pow(t14734,5.0/3.0);
                double t14798 = b4p*t14741*t14796*(2.0/3.0);
                double t14799 = b2p*c*t14756*(1.0/3.0);
                double t14800 = b1p*c*t14756*t14748*(1.0/6.0);
                double t14801 = b3p*c*t14737*t14756*(1.0/2.0);
                double t14757 = t14800+t14801+t14798+t14799;
                double t14758 = rho_a-rho_b;
                double t14759 = 1.0/t14734;
                double t14760 = t14758*t14759;
                double t14761 = two_13*2.0;
                double t14762 = t14761-2.0;
                double t14763 = 1.0/t14762;
                double t14764 = 1.0/c0f;
                double t14765 = b1f*t14737;
                double t14766 = b3f*t14739;
                double t14767 = b4f*t14741*t14742;
                double t14768 = b2f*c*t14735;
                double t14769 = t14765+t14766+t14767+t14768;
                double t14770 = 1.0/t14769;
                double t14771 = t14770*t14764*(1.0/2.0);
                double t14772 = t14771+1.0;
                double t14773 = log(t14772);
                double t14774 = a1f*c*t14735;
                double t14775 = t14774+1.0;
                double t14776 = log(t14754);
                double t14784 = c0f*t14773*t14775*2.0;
                double t14785 = c0p*t14750*t14776*2.0;
                double t14777 = t14784-t14785;
                double t14778 = t14758*t14758;
                double t14779 = t14760+1.0;
                double t14780 = pow(t14779,4.0/3.0);
                double t14781 = -t14760+1.0;
                double t14782 = pow(t14781,4.0/3.0);
                double t14783 = t14780+t14782-2.0;
                double t14786 = 1.0/(t14734*t14734);
                double t14787 = t14758*t14786;
                double t14788 = t14759+t14787;
                double t14789 = 1.0/(t14734*t14734*t14734*t14734);
                double t14790 = pow(t14779,1.0/3.0);
                double t14791 = t14790*t14788*(4.0/3.0);
                double t14792 = pow(t14781,1.0/3.0);
                double t14825 = t14792*t14788*(4.0/3.0);
                double t14793 = -t14825+t14791;
                double t14794 = 1.0/(t14734*t14734*t14734*t14734*t14734);
                double t14795 = t14778*t14778;
                double t14797 = 1.0/(t14745*t14745);
                double t14802 = 1.0/t14772;
                double t14803 = 1.0/(t14769*t14769);
                double t14804 = b4f*t14741*t14796*(2.0/3.0);
                double t14805 = b2f*c*t14756*(1.0/3.0);
                double t14806 = b1f*c*t14756*t14748*(1.0/6.0);
                double t14807 = b3f*c*t14737*t14756*(1.0/2.0);
                double t14808 = t14804+t14805+t14806+t14807;
                double t14809 = t14802*t14803*t14808*t14775;
                double t14810 = a1p*c*c0p*t14756*t14776*(2.0/3.0);
                double t14826 = t14750*t14755*t14757*t14797;
                double t14827 = a1f*c*c0f*t14773*t14756*(2.0/3.0);
                double t14811 = t14810-t14826+t14809-t14827;
                double t14812 = t14757*t14757;
                double t14813 = 1.0/pow(t14736,3.0/2.0);
                double t14814 = b4p*t14741*t14747*(1.0E1/9.0);
                double t14815 = b2p*c*t14746*(4.0/9.0);
                double t14816 = b1p*c*t14746*t14748*(2.0/9.0);
                double t14817 = b3p*c*t14737*t14746*(2.0/3.0);
                double t14818 = b3p*t14741*t14747*t14748*(1.0/1.2E1);
                double t14819 = t14814+t14815+t14816+t14817+t14818-b1p*t14741*t14813*t14747*(1.0/3.6E1);
                double t14820 = 1.0/(t14745*t14745*t14745);
                double t14821 = t14820*t14812*t14750*t14755*2.0;
                double t14822 = t14808*t14808;
                double t14823 = 1.0/(t14754*t14754);
                double t14824 = 1.0/(t14745*t14745*t14745*t14745);
                double t14828 = t14786*2.0;
                double t14829 = 1.0/(t14734*t14734*t14734);
                double t14830 = t14829*t14758*2.0;
                double t14831 = t14830+t14828;
                double t14832 = t14788*t14788;
                double t14833 = a1p*c*c0p*t14746*t14776*(8.0/9.0);
                double t14834 = 1.0/(t14734*t14734*t14734*t14734*t14734*t14734);
                double t14835 = 1.0/d2fz0;
                double t14836 = 1.0/Aa;
                double t14837 = b1a*t14737;
                double t14838 = b3a*t14739;
                double t14839 = b4a*t14741*t14742;
                double t14840 = b2a*c*t14735;
                double t14841 = t14840+t14837+t14838+t14839;
                double t14842 = 1.0/t14841;
                double t14843 = t14842*t14836*(1.0/2.0);
                double t14844 = t14843+1.0;
                double t14845 = log(t14844);
                double t14846 = a1a*c*t14735;
                double t14847 = t14846+1.0;
                double t14848 = t14831*t14790*(4.0/3.0);
                double t14849 = 1.0/pow(t14779,2.0/3.0);
                double t14850 = t14832*t14849*(4.0/9.0);
                double t14851 = 1.0/pow(t14781,2.0/3.0);
                double t14852 = t14832*t14851*(4.0/9.0);
                double t14853 = t14850+t14852+t14848-t14831*t14792*(4.0/3.0);
                double t14854 = t14758*t14778*t14789*4.0;
                double t14855 = t14794*t14795*4.0;
                double t14856 = t14854+t14855;
                double t14857 = 1.0/t14844;
                double t14858 = t14795*t14789;
                double t14859 = t14858-1.0;
                double t14860 = 1.0/(t14841*t14841);
                double t14861 = b4a*t14741*t14796*(2.0/3.0);
                double t14862 = b2a*c*t14756*(1.0/3.0);
                double t14863 = b1a*c*t14756*t14748*(1.0/6.0);
                double t14864 = b3a*c*t14737*t14756*(1.0/2.0);
                double t14865 = t14861+t14862+t14863+t14864;
                double t14866 = t14865*t14865;
                v_rho_b_rho_b[Q] += scale * (-t14734*(t14821+t14833-t14750*t14755*t14819*t14797-t14763*t14783*t14795*t14789*(t14821+t14833-t14750*t14755*t14819*t14797-t14802*t14822*t14775*1.0/(t14769*t14769*t14769)*2.0+t14802*t14803*t14775*(b2f*c*t14746*(4.0/9.0)+b4f*t14741*t14747*(1.0E1/9.0)+b1f*c*t14746*t14748*(2.0/9.0)+b3f*c*t14737*t14746*(2.0/3.0)-b1f*t14741*t14813*t14747*(1.0/3.6E1)+b3f*t14741*t14747*t14748*(1.0/1.2E1))+t14822*1.0/(t14772*t14772)*t14764*t14775*1.0/(t14769*t14769*t14769*t14769)*(1.0/2.0)-a1f*c*c0f*t14746*t14773*(8.0/9.0)-t14812*t14750*t14751*t14823*t14824*(1.0/2.0)+a1f*c*t14802*t14803*t14808*t14756*(2.0/3.0)-a1p*c*t14755*t14756*t14757*t14797*(2.0/3.0))-t14812*t14750*t14751*t14823*t14824*(1.0/2.0)-t14811*t14763*t14783*t14794*t14795*8.0-t14811*t14763*t14793*t14795*t14789*2.0+t14834*t14763*t14783*t14777*t14795*2.0E1+t14763*t14793*t14794*t14777*t14795*8.0+t14763*t14853*t14777*t14795*t14789+t14763*t14783*t14777*t14778*t14789*1.2E1-t14811*t14763*t14783*t14758*t14778*t14789*8.0+t14763*t14783*t14758*t14794*t14777*t14778*3.2E1+t14763*t14793*t14758*t14777*t14778*t14789*8.0-a1p*c*t14755*t14756*t14757*t14797*(2.0/3.0)+Aa*t14763*t14835*t14853*t14845*t14847*t14859*2.0+Aa*t14763*t14835*t14845*t14793*t14847*t14856*4.0+Aa*t14763*t14835*t14845*t14783*t14847*(t14834*t14795*2.0E1+t14778*t14789*1.2E1+t14758*t14794*t14778*3.2E1)*2.0-t14860*t14763*t14835*t14783*t14847*t14856*t14865*t14857*2.0-t14860*t14763*t14835*t14793*t14847*t14865*t14857*t14859*2.0+1.0/(t14841*t14841*t14841)*t14763*t14835*t14783*t14847*t14857*t14866*t14859*2.0-t14860*t14763*t14835*t14783*t14847*t14857*t14859*(b2a*c*t14746*(4.0/9.0)+b4a*t14741*t14747*(1.0E1/9.0)+b1a*c*t14746*t14748*(2.0/9.0)+b3a*c*t14737*t14746*(2.0/3.0)-b1a*t14741*t14813*t14747*(1.0/3.6E1)+b3a*t14741*t14747*t14748*(1.0/1.2E1))+Aa*a1a*c*t14763*t14835*t14845*t14756*t14783*t14856*(4.0/3.0)+Aa*a1a*c*t14763*t14835*t14746*t14845*t14783*t14859*(8.0/9.0)+Aa*a1a*c*t14763*t14835*t14845*t14756*t14793*t14859*(4.0/3.0)-1.0/(t14841*t14841*t14841*t14841)*t14763*t14835*1.0/(t14844*t14844)*t14836*t14783*t14847*t14866*t14859*(1.0/2.0)-a1a*c*t14860*t14763*t14835*t14756*t14783*t14865*t14857*t14859*(2.0/3.0))-t14750*t14755*t14757*t14797*2.0+a1p*c*c0p*t14756*t14776*(4.0/3.0)-t14811*t14763*t14783*t14795*t14789*2.0+t14763*t14783*t14794*t14777*t14795*8.0+t14763*t14793*t14777*t14795*t14789*2.0+t14763*t14783*t14758*t14777*t14778*t14789*8.0+Aa*t14763*t14835*t14845*t14783*t14847*t14856*4.0+Aa*t14763*t14835*t14845*t14793*t14847*t14859*4.0-t14860*t14763*t14835*t14783*t14847*t14865*t14857*t14859*2.0+Aa*a1a*c*t14763*t14835*t14845*t14756*t14783*t14859*(4.0/3.0));
            }
            
        }
    }
}

}