#include <libmints/vector.h>
#include "LYP_Cfunctional.h"
#include "utility.h"
#include <cmath>

using namespace psi;

namespace psi {

LYP_CFunctional::LYP_CFunctional()
{
    name_ = "LYP_C";
    description_ = "    LYP Correlation\n";
    citation_ = "    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 (1989)\n";
    alpha_ = 1.0;
    omega_ = 0.0;
    lrc_ = false;
    gga_ = true;
    meta_ = false;
    parameters_["A"] =   4.9180000000000001E-02;
    parameters_["B"] =   1.3200000000000001E-01;
    parameters_["C"] =   2.5330000000000003E-01;
    parameters_["Dd"] =   3.4899999999999998E-01;
    parameters_["CFext"] =   3.6462398978764774E+01;
}
LYP_CFunctional::~LYP_CFunctional()
{
}
void LYP_CFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    double A = parameters_["A"];
    double B = parameters_["B"];
    double C = parameters_["C"];
    double Dd = parameters_["Dd"];
    double CFext = parameters_["CFext"];

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
        } else if (rho_b < lsda_cutoff_) {
        } else {
            // v
            if (deriv >= 0) {
                double t45993 = rho_a+rho_b;
                double t45994 = 1.0/pow(t45993,1.0/3.0);
                double t45995 = Dd*t45994;
                double t45996 = t45995+1.0;
                double t45997 = 1.0/t45996;
                double t45998 = t45993*t45993;
                double t45999 = t45998*(2.0/3.0);
                double t46000 = gamma_ab*2.0;
                double t46001 = gamma_aa+gamma_bb+t46000;
                double t46002 = 1.0/t45993;
                v[Q] += scale * A*rho_a*rho_b*t45997*t46002*-4.0-A*B*1.0/pow(t45993,1.1E1/3.0)*t45997*exp(-C*t45994)*(t45998*t46001*(-2.0/3.0)+gamma_aa*(t45999-rho_b*rho_b)+gamma_bb*(t45999-rho_a*rho_a)+rho_a*rho_b*((gamma_aa+gamma_bb)*(C*t45994*(1.0/1.8E1)+Dd*t45994*t45997*(1.0/1.8E1)-5.0/2.0)+CFext*(pow(rho_a,8.0/3.0)+pow(rho_b,8.0/3.0))-t46001*(C*t45994*(7.0/1.8E1)+Dd*t45994*t45997*(7.0/1.8E1)-4.7E1/1.8E1)-t46002*(gamma_aa*rho_a+gamma_bb*rho_b)*(C*t45994*(1.0/9.0)+Dd*t45994*t45997*(1.0/9.0)-1.1E1/9.0)));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t46004 = rho_a+rho_b;
                double t46005 = 1.0/pow(t46004,1.0/3.0);
                double t46006 = Dd*t46005;
                double t46007 = t46006+1.0;
                double t46008 = 1.0/t46007;
                double t46009 = t46004*t46004;
                double t46010 = t46009*(2.0/3.0);
                double t46011 = gamma_ab*2.0;
                double t46012 = gamma_aa+gamma_bb+t46011;
                double t46013 = 1.0/t46004;
                double t46040 = C*t46005;
                double t46014 = exp(-t46040);
                double t46015 = C*t46005*(7.0/1.8E1);
                double t46016 = Dd*t46005*t46008*(7.0/1.8E1);
                double t46017 = t46015+t46016-4.7E1/1.8E1;
                double t46018 = t46012*t46017;
                double t46019 = gamma_aa+gamma_bb;
                double t46020 = C*t46005*(1.0/1.8E1);
                double t46021 = Dd*t46005*t46008*(1.0/1.8E1);
                double t46022 = t46020+t46021-5.0/2.0;
                double t46023 = pow(rho_a,8.0/3.0);
                double t46024 = pow(rho_b,8.0/3.0);
                double t46025 = t46023+t46024;
                double t46026 = gamma_aa*rho_a;
                double t46027 = gamma_bb*rho_b;
                double t46028 = t46026+t46027;
                double t46029 = C*t46005*(1.0/9.0);
                double t46030 = Dd*t46005*t46008*(1.0/9.0);
                double t46031 = t46030+t46029-1.1E1/9.0;
                double t46032 = t46013*t46031*t46028;
                double t46047 = t46022*t46019;
                double t46048 = CFext*t46025;
                double t46033 = t46032+t46018-t46047-t46048;
                double t46034 = rho_b*(4.0/3.0);
                double t46035 = 1.0/pow(t46004,4.0/3.0);
                double t46036 = 1.0/(t46007*t46007);
                double t46037 = Dd*Dd;
                double t46038 = 1.0/pow(t46004,5.0/3.0);
                double t46039 = 1.0/(t46004*t46004);
                double t46041 = rho_b*rho_b;
                double t46042 = t46010-t46041;
                double t46043 = gamma_aa*t46042;
                double t46044 = rho_a*rho_a;
                double t46045 = t46010-t46044;
                double t46046 = gamma_bb*t46045;
                double t46051 = t46012*t46009*(2.0/3.0);
                double t46052 = rho_a*rho_b*t46033;
                double t46049 = -t46051+t46043-t46052+t46046;
                double t46050 = 1.0/(t46004*t46004*t46004*t46004*t46004);
                v_rho_a[Q] += scale * A*rho_b*t46013*t46008*-4.0+A*rho_a*rho_b*t46008*t46039*4.0+A*B*1.0/pow(t46004,1.1E1/3.0)*t46014*t46008*(rho_b*t46033-gamma_aa*(rho_a*(4.0/3.0)+t46034)+gamma_bb*(rho_a*(2.0/3.0)-t46034)+t46012*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_a,5.0/3.0)*(8.0/3.0)-t46019*(C*t46035*(1.0/5.4E1)+Dd*t46008*t46035*(1.0/5.4E1)-t46036*t46037*t46038*(1.0/5.4E1))+t46012*(C*t46035*(7.0/5.4E1)+Dd*t46008*t46035*(7.0/5.4E1)-t46036*t46037*t46038*(7.0/5.4E1))+t46013*t46028*(C*t46035*(1.0/2.7E1)+Dd*t46008*t46035*(1.0/2.7E1)-t46036*t46037*t46038*(1.0/2.7E1))-gamma_aa*t46013*t46031+t46031*t46028*t46039))-A*Dd*rho_a*rho_b*1.0/pow(t46004,7.0/3.0)*t46036*(4.0/3.0)+A*B*1.0/pow(t46004,1.4E1/3.0)*t46014*t46008*t46049*(1.1E1/3.0)-A*B*C*t46014*t46050*t46008*t46049*(1.0/3.0)-A*B*Dd*t46014*t46050*t46036*t46049*(1.0/3.0);
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t46054 = rho_a+rho_b;
                double t46055 = 1.0/pow(t46054,1.0/3.0);
                double t46056 = Dd*t46055;
                double t46057 = t46056+1.0;
                double t46058 = 1.0/t46057;
                double t46059 = t46054*t46054;
                double t46060 = t46059*(2.0/3.0);
                double t46061 = gamma_ab*2.0;
                double t46062 = gamma_aa+gamma_bb+t46061;
                double t46063 = 1.0/t46054;
                double t46090 = C*t46055;
                double t46064 = exp(-t46090);
                double t46065 = C*t46055*(7.0/1.8E1);
                double t46066 = Dd*t46055*t46058*(7.0/1.8E1);
                double t46067 = t46065+t46066-4.7E1/1.8E1;
                double t46068 = t46062*t46067;
                double t46069 = gamma_aa+gamma_bb;
                double t46070 = C*t46055*(1.0/1.8E1);
                double t46071 = Dd*t46055*t46058*(1.0/1.8E1);
                double t46072 = t46070+t46071-5.0/2.0;
                double t46073 = pow(rho_a,8.0/3.0);
                double t46074 = pow(rho_b,8.0/3.0);
                double t46075 = t46073+t46074;
                double t46076 = gamma_aa*rho_a;
                double t46077 = gamma_bb*rho_b;
                double t46078 = t46076+t46077;
                double t46079 = C*t46055*(1.0/9.0);
                double t46080 = Dd*t46055*t46058*(1.0/9.0);
                double t46081 = t46080+t46079-1.1E1/9.0;
                double t46082 = t46063*t46081*t46078;
                double t46097 = t46072*t46069;
                double t46098 = CFext*t46075;
                double t46083 = t46082+t46068-t46097-t46098;
                double t46084 = rho_a*(4.0/3.0);
                double t46085 = 1.0/pow(t46054,4.0/3.0);
                double t46086 = 1.0/(t46057*t46057);
                double t46087 = Dd*Dd;
                double t46088 = 1.0/pow(t46054,5.0/3.0);
                double t46089 = 1.0/(t46054*t46054);
                double t46091 = rho_b*rho_b;
                double t46092 = t46060-t46091;
                double t46093 = gamma_aa*t46092;
                double t46094 = rho_a*rho_a;
                double t46095 = t46060-t46094;
                double t46096 = gamma_bb*t46095;
                double t46101 = t46062*t46059*(2.0/3.0);
                double t46102 = rho_a*rho_b*t46083;
                double t46099 = -t46101-t46102+t46093+t46096;
                double t46100 = 1.0/(t46054*t46054*t46054*t46054*t46054);
                v_rho_b[Q] += scale * A*rho_a*t46063*t46058*-4.0+A*rho_a*rho_b*t46058*t46089*4.0+A*B*1.0/pow(t46054,1.1E1/3.0)*t46064*t46058*(rho_a*t46083-gamma_bb*(rho_b*(4.0/3.0)+t46084)+gamma_aa*(rho_b*(2.0/3.0)-t46084)+t46062*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_b,5.0/3.0)*(8.0/3.0)-t46069*(C*t46085*(1.0/5.4E1)+Dd*t46058*t46085*(1.0/5.4E1)-t46086*t46087*t46088*(1.0/5.4E1))+t46062*(C*t46085*(7.0/5.4E1)+Dd*t46058*t46085*(7.0/5.4E1)-t46086*t46087*t46088*(7.0/5.4E1))+t46063*t46078*(C*t46085*(1.0/2.7E1)+Dd*t46058*t46085*(1.0/2.7E1)-t46086*t46087*t46088*(1.0/2.7E1))-gamma_bb*t46063*t46081+t46081*t46078*t46089))-A*Dd*rho_a*rho_b*1.0/pow(t46054,7.0/3.0)*t46086*(4.0/3.0)+A*B*1.0/pow(t46054,1.4E1/3.0)*t46064*t46058*t46099*(1.1E1/3.0)-A*B*C*t46100*t46064*t46058*t46099*(1.0/3.0)-A*B*Dd*t46100*t46064*t46086*t46099*(1.0/3.0);
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t46104 = rho_a+rho_b;
                double t46105 = 1.0/pow(t46104,1.0/3.0);
                double t46106 = Dd*t46105;
                double t46107 = t46106+1.0;
                double t46108 = 1.0/t46107;
                v_gamma_aa[Q] += scale * A*B*1.0/pow(t46104,1.1E1/3.0)*t46108*exp(-C*t46105)*(rho_b*rho_b+rho_a*rho_b*(C*t46105*(1.0/3.0)+Dd*t46105*t46108*(1.0/3.0)+(rho_a*(C*t46105*(1.0/9.0)+Dd*t46105*t46108*(1.0/9.0)-1.1E1/9.0))/t46104-1.0/9.0));
            }
            
            // v_gamma_ab
            if (deriv >= 1) {
                double t46110 = rho_a+rho_b;
                double t46111 = 1.0/pow(t46110,1.0/3.0);
                double t46112 = Dd*t46111;
                double t46113 = t46112+1.0;
                double t46114 = 1.0/t46113;
                v_gamma_ab[Q] += scale * A*B*1.0/pow(t46110,1.1E1/3.0)*t46114*exp(-C*t46111)*((t46110*t46110)*(4.0/3.0)+rho_a*rho_b*(C*t46111*(7.0/9.0)+Dd*t46111*t46114*(7.0/9.0)-4.7E1/9.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t46116 = rho_a+rho_b;
                double t46117 = 1.0/pow(t46116,1.0/3.0);
                double t46118 = Dd*t46117;
                double t46119 = t46118+1.0;
                double t46120 = 1.0/t46119;
                v_gamma_bb[Q] += scale * A*B*t46120*1.0/pow(t46116,1.1E1/3.0)*exp(-C*t46117)*(rho_a*rho_a+rho_a*rho_b*(C*t46117*(1.0/3.0)+Dd*t46120*t46117*(1.0/3.0)+(rho_b*(C*t46117*(1.0/9.0)+Dd*t46120*t46117*(1.0/9.0)-1.1E1/9.0))/t46116-1.0/9.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t46124 = rho_a+rho_b;
                double t46125 = 1.0/pow(t46124,1.0/3.0);
                double t46126 = Dd*t46125;
                double t46127 = t46126+1.0;
                double t46128 = 1.0/t46127;
                double t46129 = t46124*t46124;
                double t46130 = t46129*(2.0/3.0);
                double t46131 = gamma_ab*2.0;
                double t46132 = gamma_aa+gamma_bb+t46131;
                double t46133 = 1.0/(t46127*t46127);
                double t46160 = C*t46125;
                double t46134 = exp(-t46160);
                double t46135 = C*t46125*(7.0/1.8E1);
                double t46136 = Dd*t46125*t46128*(7.0/1.8E1);
                double t46137 = t46135+t46136-4.7E1/1.8E1;
                double t46138 = t46132*t46137;
                double t46139 = gamma_aa+gamma_bb;
                double t46140 = C*t46125*(1.0/1.8E1);
                double t46141 = Dd*t46125*t46128*(1.0/1.8E1);
                double t46142 = t46140+t46141-5.0/2.0;
                double t46143 = pow(rho_a,8.0/3.0);
                double t46144 = pow(rho_b,8.0/3.0);
                double t46145 = t46143+t46144;
                double t46146 = 1.0/t46124;
                double t46147 = gamma_aa*rho_a;
                double t46148 = gamma_bb*rho_b;
                double t46149 = t46147+t46148;
                double t46150 = C*t46125*(1.0/9.0);
                double t46151 = Dd*t46125*t46128*(1.0/9.0);
                double t46152 = t46150+t46151-1.1E1/9.0;
                double t46153 = t46152*t46146*t46149;
                double t46186 = t46142*t46139;
                double t46187 = CFext*t46145;
                double t46154 = t46153+t46138-t46186-t46187;
                double t46155 = rho_b*(4.0/3.0);
                double t46156 = 1.0/pow(t46124,4.0/3.0);
                double t46157 = Dd*Dd;
                double t46158 = 1.0/pow(t46124,5.0/3.0);
                double t46159 = 1.0/(t46124*t46124);
                double t46161 = 1.0/pow(t46124,1.1E1/3.0);
                double t46162 = C*t46156*(1.0/5.4E1);
                double t46163 = Dd*t46128*t46156*(1.0/5.4E1);
                double t46194 = t46133*t46157*t46158*(1.0/5.4E1);
                double t46164 = t46162+t46163-t46194;
                double t46165 = pow(rho_a,5.0/3.0);
                double t46166 = CFext*t46165*(8.0/3.0);
                double t46167 = C*t46156*(7.0/5.4E1);
                double t46168 = Dd*t46128*t46156*(7.0/5.4E1);
                double t46196 = t46133*t46157*t46158*(7.0/5.4E1);
                double t46169 = t46167+t46168-t46196;
                double t46170 = t46132*t46169;
                double t46171 = C*t46156*(1.0/2.7E1);
                double t46172 = Dd*t46128*t46156*(1.0/2.7E1);
                double t46181 = t46133*t46157*t46158*(1.0/2.7E1);
                double t46173 = t46171+t46172-t46181;
                double t46174 = t46146*t46173*t46149;
                double t46175 = t46152*t46149*t46159;
                double t46195 = t46164*t46139;
                double t46197 = gamma_aa*t46152*t46146;
                double t46176 = t46170+t46174+t46166+t46175-t46195-t46197;
                double t46177 = 1.0/pow(t46124,7.0/3.0);
                double t46178 = 1.0/(t46124*t46124*t46124);
                double t46179 = 1.0/(t46127*t46127*t46127);
                double t46180 = 1.0/pow(t46124,8.0/3.0);
                double t46182 = rho_a*2.0;
                double t46183 = rho_b*2.0;
                double t46184 = t46182+t46183;
                double t46185 = t46132*t46184*(2.0/3.0);
                double t46188 = rho_b*t46154;
                double t46189 = rho_a*(4.0/3.0);
                double t46190 = t46155+t46189;
                double t46191 = rho_a*(2.0/3.0);
                double t46192 = t46155-t46191;
                double t46193 = gamma_bb*t46192;
                double t46198 = rho_a*rho_b*t46176;
                double t46199 = gamma_aa*t46190;
                double t46200 = t46193-t46185-t46188+t46198+t46199;
                double t46201 = 1.0/(t46124*t46124*t46124*t46124*t46124);
                double t46202 = rho_b*rho_b;
                double t46203 = t46130-t46202;
                double t46204 = gamma_aa*t46203;
                double t46205 = rho_a*rho_a;
                double t46206 = t46130-t46205;
                double t46207 = gamma_bb*t46206;
                double t46210 = t46132*t46129*(2.0/3.0);
                double t46211 = rho_a*rho_b*t46154;
                double t46208 = -t46210-t46211+t46204+t46207;
                double t46209 = 1.0/(t46124*t46124*t46124*t46124*t46124*t46124);
                double t46212 = 1.0/pow(t46124,1.9E1/3.0);
                v_rho_a_rho_a[Q] += scale * A*rho_b*t46128*t46159*8.0-A*Dd*rho_b*t46133*t46177*(8.0/3.0)-A*rho_a*rho_b*t46128*t46178*8.0+A*Dd*rho_a*rho_b*1.0/pow(t46124,1.0E1/3.0)*t46133*(4.0E1/9.0)+A*B*t46200*1.0/pow(t46124,1.4E1/3.0)*t46134*t46128*(2.2E1/3.0)-A*B*1.0/pow(t46124,1.7E1/3.0)*t46134*t46208*t46128*(1.54E2/9.0)+A*B*t46134*t46161*t46128*(gamma_ab*(8.0/3.0)+gamma_bb*2.0-rho_b*t46176*2.0-rho_a*rho_b*(CFext*pow(rho_a,2.0/3.0)*(4.0E1/9.0)+t46139*(C*t46177*(2.0/8.1E1)+Dd*t46128*t46177*(2.0/8.1E1)-t46133*t46180*t46157*(1.0/2.7E1)+Dd*t46157*t46178*t46179*(1.0/8.1E1))-t46132*(C*t46177*(1.4E1/8.1E1)+Dd*t46128*t46177*(1.4E1/8.1E1)-t46133*t46180*t46157*(7.0/2.7E1)+Dd*t46157*t46178*t46179*(7.0/8.1E1))-t46146*t46149*(C*t46177*(4.0/8.1E1)+Dd*t46128*t46177*(4.0/8.1E1)-t46133*t46180*t46157*(2.0/2.7E1)+Dd*t46157*t46178*t46179*(2.0/8.1E1))+gamma_aa*t46146*t46173*2.0+gamma_aa*t46152*t46159*2.0-t46152*t46149*t46178*2.0-t46173*t46149*t46159*2.0))-A*rho_a*rho_b*t46161*t46157*t46179*(8.0/9.0)-A*B*t46212*t46134*t46208*t46157*t46179*(2.0/9.0)-A*B*(C*C)*t46212*t46134*t46208*t46128*(1.0/9.0)-A*B*C*t46200*t46201*t46134*t46128*(2.0/3.0)+A*B*C*t46134*t46208*t46128*t46209*(2.6E1/9.0)-A*B*Dd*t46200*t46201*t46133*t46134*(2.0/3.0)+A*B*Dd*t46133*t46134*t46208*t46209*(2.6E1/9.0)-A*B*C*Dd*t46212*t46133*t46134*t46208*(2.0/9.0);
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t46214 = rho_a+rho_b;
                double t46215 = 1.0/pow(t46214,1.0/3.0);
                double t46216 = Dd*t46215;
                double t46217 = t46216+1.0;
                double t46218 = 1.0/t46217;
                double t46219 = 1.0/(t46214*t46214);
                double t46220 = 1.0/pow(t46214,7.0/3.0);
                double t46221 = 1.0/(t46217*t46217);
                double t46222 = t46214*t46214;
                double t46223 = t46222*(2.0/3.0);
                double t46224 = gamma_ab*2.0;
                double t46225 = gamma_aa+gamma_bb+t46224;
                double t46226 = 1.0/t46214;
                double t46265 = C*t46215;
                double t46227 = exp(-t46265);
                double t46228 = C*t46215*(7.0/1.8E1);
                double t46229 = Dd*t46215*t46218*(7.0/1.8E1);
                double t46230 = t46228+t46229-4.7E1/1.8E1;
                double t46231 = t46230*t46225;
                double t46232 = gamma_aa+gamma_bb;
                double t46233 = 1.0/pow(t46214,4.0/3.0);
                double t46234 = Dd*Dd;
                double t46235 = 1.0/pow(t46214,5.0/3.0);
                double t46236 = C*t46215*(1.0/9.0);
                double t46237 = Dd*t46215*t46218*(1.0/9.0);
                double t46238 = t46236+t46237-1.1E1/9.0;
                double t46239 = gamma_aa*rho_a;
                double t46240 = gamma_bb*rho_b;
                double t46241 = t46240+t46239;
                double t46242 = C*t46233*(1.0/5.4E1);
                double t46243 = Dd*t46233*t46218*(1.0/5.4E1);
                double t46269 = t46221*t46234*t46235*(1.0/5.4E1);
                double t46244 = t46242+t46243-t46269;
                double t46245 = C*t46233*(7.0/5.4E1);
                double t46246 = Dd*t46233*t46218*(7.0/5.4E1);
                double t46272 = t46221*t46234*t46235*(7.0/5.4E1);
                double t46247 = t46245-t46272+t46246;
                double t46248 = t46225*t46247;
                double t46249 = C*t46233*(1.0/2.7E1);
                double t46250 = Dd*t46233*t46218*(1.0/2.7E1);
                double t46263 = t46221*t46234*t46235*(1.0/2.7E1);
                double t46251 = t46250-t46263+t46249;
                double t46252 = t46241*t46251*t46226;
                double t46253 = t46241*t46219*t46238;
                double t46254 = C*t46215*(1.0/1.8E1);
                double t46255 = Dd*t46215*t46218*(1.0/1.8E1);
                double t46256 = t46254+t46255-5.0/2.0;
                double t46257 = pow(rho_a,8.0/3.0);
                double t46258 = pow(rho_b,8.0/3.0);
                double t46259 = t46257+t46258;
                double t46260 = 1.0/(t46214*t46214*t46214);
                double t46261 = 1.0/(t46217*t46217*t46217);
                double t46262 = 1.0/pow(t46214,8.0/3.0);
                double t46264 = t46241*t46226*t46238;
                double t46266 = t46232*t46256;
                double t46267 = CFext*t46259;
                double t46268 = rho_b*(4.0/3.0);
                double t46270 = pow(rho_a,5.0/3.0);
                double t46271 = CFext*t46270*(8.0/3.0);
                double t46282 = t46232*t46244;
                double t46290 = gamma_aa*t46226*t46238;
                double t46273 = t46252+t46253+t46271-t46290-t46282+t46248;
                double t46274 = 1.0/pow(t46214,1.4E1/3.0);
                double t46275 = rho_a*2.0;
                double t46276 = rho_b*2.0;
                double t46277 = t46275+t46276;
                double t46278 = t46225*t46277*(2.0/3.0);
                double t46279 = t46231+t46264-t46266-t46267;
                double t46280 = rho_a*(4.0/3.0);
                double t46281 = t46280+t46268;
                double t46283 = pow(rho_b,5.0/3.0);
                double t46284 = CFext*t46283*(8.0/3.0);
                double t46285 = 1.0/pow(t46214,1.1E1/3.0);
                double t46286 = rho_b*t46279;
                double t46287 = rho_a*(2.0/3.0);
                double t46288 = t46268-t46287;
                double t46289 = gamma_bb*t46288;
                double t46291 = rho_a*rho_b*t46273;
                double t46292 = gamma_aa*t46281;
                double t46293 = t46291+t46292-t46286-t46278+t46289;
                double t46294 = 1.0/(t46214*t46214*t46214*t46214*t46214);
                double t46295 = rho_a*t46279;
                double t46296 = rho_b*(2.0/3.0);
                double t46297 = t46280-t46296;
                double t46298 = gamma_aa*t46297;
                double t46299 = gamma_bb*t46281;
                double t46303 = gamma_bb*t46226*t46238;
                double t46300 = -t46303+t46252+t46253-t46282+t46248+t46284;
                double t46301 = rho_a*rho_b*t46300;
                double t46302 = t46301-t46295-t46278+t46298+t46299;
                double t46304 = rho_b*rho_b;
                double t46305 = t46223-t46304;
                double t46306 = gamma_aa*t46305;
                double t46307 = rho_a*rho_a;
                double t46308 = t46223-t46307;
                double t46309 = gamma_bb*t46308;
                double t46310 = 1.0/(t46214*t46214*t46214*t46214*t46214*t46214);
                double t46312 = t46222*t46225*(2.0/3.0);
                double t46313 = rho_a*rho_b*t46279;
                double t46311 = -t46312-t46313+t46306+t46309;
                double t46314 = 1.0/pow(t46214,1.9E1/3.0);
                v_rho_a_rho_b[Q] += scale * A*t46226*t46218*-4.0+A*rho_a*t46218*t46219*4.0+A*rho_b*t46218*t46219*4.0-A*Dd*rho_a*t46220*t46221*(4.0/3.0)-A*Dd*rho_b*t46220*t46221*(4.0/3.0)-A*rho_a*rho_b*t46260*t46218*8.0+A*Dd*rho_a*rho_b*t46221*1.0/pow(t46214,1.0E1/3.0)*(4.0E1/9.0)-A*B*1.0/pow(t46214,1.7E1/3.0)*t46218*t46227*(t46306+t46309-t46222*t46225*(2.0/3.0)-rho_a*rho_b*(t46231+t46264-CFext*t46259-t46232*t46256))*(1.54E2/9.0)-A*B*t46218*t46227*t46285*(gamma_ab*(-8.0/3.0)-t46231-t46264+t46266+t46267+rho_a*t46273+rho_b*(t46252+t46253+t46248+t46284-t46232*t46244-gamma_bb*t46226*t46238)+rho_a*rho_b*(t46232*(C*t46220*(2.0/8.1E1)+Dd*t46220*t46218*(2.0/8.1E1)-t46221*t46234*t46262*(1.0/2.7E1)+Dd*t46260*t46234*t46261*(1.0/8.1E1))-t46225*(C*t46220*(1.4E1/8.1E1)+Dd*t46220*t46218*(1.4E1/8.1E1)-t46221*t46234*t46262*(7.0/2.7E1)+Dd*t46260*t46234*t46261*(7.0/8.1E1))-t46241*t46226*(C*t46220*(4.0/8.1E1)+Dd*t46220*t46218*(4.0/8.1E1)-t46221*t46234*t46262*(2.0/2.7E1)+Dd*t46260*t46234*t46261*(2.0/8.1E1))+gamma_aa*t46251*t46226+gamma_bb*t46251*t46226+gamma_aa*t46219*t46238+gamma_bb*t46219*t46238-t46241*t46251*t46219*2.0-t46241*t46260*t46238*2.0))+A*B*t46302*t46218*t46227*t46274*(1.1E1/3.0)+A*B*t46218*t46227*t46274*t46293*(1.1E1/3.0)-A*rho_a*rho_b*t46234*t46261*t46285*(8.0/9.0)-A*B*t46311*t46314*t46234*t46261*t46227*(2.0/9.0)-A*B*(C*C)*t46311*t46314*t46218*t46227*(1.0/9.0)+A*B*C*t46310*t46311*t46218*t46227*(2.6E1/9.0)-A*B*C*t46302*t46218*t46227*t46294*(1.0/3.0)-A*B*C*t46218*t46227*t46293*t46294*(1.0/3.0)+A*B*Dd*t46310*t46221*t46311*t46227*(2.6E1/9.0)-A*B*Dd*t46221*t46302*t46227*t46294*(1.0/3.0)-A*B*Dd*t46221*t46227*t46293*t46294*(1.0/3.0)-A*B*C*Dd*t46221*t46311*t46314*t46227*(2.0/9.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t46316 = rho_a+rho_b;
                double t46317 = 1.0/pow(t46316,1.0/3.0);
                double t46318 = Dd*t46317;
                double t46319 = t46318+1.0;
                double t46320 = 1.0/t46319;
                double t46321 = t46316*t46316;
                double t46322 = t46321*(2.0/3.0);
                double t46323 = gamma_ab*2.0;
                double t46324 = gamma_aa+gamma_bb+t46323;
                double t46325 = 1.0/(t46319*t46319);
                double t46352 = C*t46317;
                double t46326 = exp(-t46352);
                double t46327 = C*t46317*(7.0/1.8E1);
                double t46328 = Dd*t46320*t46317*(7.0/1.8E1);
                double t46329 = t46327+t46328-4.7E1/1.8E1;
                double t46330 = t46324*t46329;
                double t46331 = gamma_aa+gamma_bb;
                double t46332 = C*t46317*(1.0/1.8E1);
                double t46333 = Dd*t46320*t46317*(1.0/1.8E1);
                double t46334 = t46332+t46333-5.0/2.0;
                double t46335 = pow(rho_a,8.0/3.0);
                double t46336 = pow(rho_b,8.0/3.0);
                double t46337 = t46335+t46336;
                double t46338 = 1.0/t46316;
                double t46339 = gamma_aa*rho_a;
                double t46340 = gamma_bb*rho_b;
                double t46341 = t46340+t46339;
                double t46342 = C*t46317*(1.0/9.0);
                double t46343 = Dd*t46320*t46317*(1.0/9.0);
                double t46344 = t46342+t46343-1.1E1/9.0;
                double t46345 = t46341*t46344*t46338;
                double t46378 = t46331*t46334;
                double t46379 = CFext*t46337;
                double t46346 = t46330+t46345-t46378-t46379;
                double t46347 = rho_a*(4.0/3.0);
                double t46348 = 1.0/pow(t46316,4.0/3.0);
                double t46349 = Dd*Dd;
                double t46350 = 1.0/pow(t46316,5.0/3.0);
                double t46351 = 1.0/(t46316*t46316);
                double t46353 = 1.0/pow(t46316,1.1E1/3.0);
                double t46354 = C*t46348*(1.0/5.4E1);
                double t46355 = Dd*t46320*t46348*(1.0/5.4E1);
                double t46386 = t46350*t46325*t46349*(1.0/5.4E1);
                double t46356 = t46354+t46355-t46386;
                double t46357 = pow(rho_b,5.0/3.0);
                double t46358 = CFext*t46357*(8.0/3.0);
                double t46359 = C*t46348*(7.0/5.4E1);
                double t46360 = Dd*t46320*t46348*(7.0/5.4E1);
                double t46388 = t46350*t46325*t46349*(7.0/5.4E1);
                double t46361 = t46360+t46359-t46388;
                double t46362 = t46324*t46361;
                double t46363 = C*t46348*(1.0/2.7E1);
                double t46364 = Dd*t46320*t46348*(1.0/2.7E1);
                double t46373 = t46350*t46325*t46349*(1.0/2.7E1);
                double t46365 = t46363+t46364-t46373;
                double t46366 = t46341*t46338*t46365;
                double t46367 = t46341*t46351*t46344;
                double t46387 = t46331*t46356;
                double t46389 = gamma_bb*t46344*t46338;
                double t46368 = t46362+t46366+t46358+t46367-t46387-t46389;
                double t46369 = 1.0/pow(t46316,7.0/3.0);
                double t46370 = 1.0/(t46316*t46316*t46316);
                double t46371 = 1.0/(t46319*t46319*t46319);
                double t46372 = 1.0/pow(t46316,8.0/3.0);
                double t46374 = rho_a*2.0;
                double t46375 = rho_b*2.0;
                double t46376 = t46374+t46375;
                double t46377 = t46324*t46376*(2.0/3.0);
                double t46380 = rho_a*t46346;
                double t46381 = rho_b*(2.0/3.0);
                double t46382 = t46381-t46347;
                double t46383 = gamma_aa*t46382;
                double t46384 = rho_b*(4.0/3.0);
                double t46385 = t46347+t46384;
                double t46390 = 1.0/(t46316*t46316*t46316*t46316*t46316);
                double t46391 = rho_b*rho_b;
                double t46392 = t46322-t46391;
                double t46393 = gamma_aa*t46392;
                double t46394 = rho_a*rho_a;
                double t46395 = t46322-t46394;
                double t46396 = gamma_bb*t46395;
                double t46399 = t46321*t46324*(2.0/3.0);
                double t46400 = rho_a*rho_b*t46346;
                double t46397 = -t46400+t46393+t46396-t46399;
                double t46398 = 1.0/(t46316*t46316*t46316*t46316*t46316*t46316);
                double t46401 = 1.0/pow(t46316,1.9E1/3.0);
                v_rho_b_rho_b[Q] += scale * A*rho_a*t46320*t46351*8.0-A*Dd*rho_a*t46325*t46369*(8.0/3.0)-A*rho_a*rho_b*t46320*t46370*8.0+A*Dd*rho_a*rho_b*1.0/pow(t46316,1.0E1/3.0)*t46325*(4.0E1/9.0)-A*B*t46320*1.0/pow(t46316,1.7E1/3.0)*t46326*t46397*(1.54E2/9.0)-A*B*t46320*1.0/pow(t46316,1.4E1/3.0)*t46326*(t46380+t46383+t46377-gamma_bb*t46385-rho_a*rho_b*t46368)*(2.2E1/3.0)+A*B*t46320*t46326*t46353*(gamma_aa*2.0+gamma_ab*(8.0/3.0)-rho_a*t46368*2.0-rho_a*rho_b*(CFext*pow(rho_b,2.0/3.0)*(4.0E1/9.0)+t46331*(C*t46369*(2.0/8.1E1)+Dd*t46320*t46369*(2.0/8.1E1)-t46325*t46372*t46349*(1.0/2.7E1)+Dd*t46370*t46371*t46349*(1.0/8.1E1))-t46324*(C*t46369*(1.4E1/8.1E1)+Dd*t46320*t46369*(1.4E1/8.1E1)-t46325*t46372*t46349*(7.0/2.7E1)+Dd*t46370*t46371*t46349*(7.0/8.1E1))-t46341*t46338*(C*t46369*(4.0/8.1E1)+Dd*t46320*t46369*(4.0/8.1E1)-t46325*t46372*t46349*(2.0/2.7E1)+Dd*t46370*t46371*t46349*(2.0/8.1E1))+gamma_bb*t46351*t46344*2.0+gamma_bb*t46338*t46365*2.0-t46341*t46370*t46344*2.0-t46341*t46351*t46365*2.0))-A*rho_a*rho_b*t46353*t46371*t46349*(8.0/9.0)-A*B*t46401*t46326*t46371*t46349*t46397*(2.0/9.0)+A*B*C*t46320*t46326*t46390*(t46380+t46383+t46377-gamma_bb*t46385-rho_a*rho_b*t46368)*(2.0/3.0)+A*B*Dd*t46325*t46326*t46390*(t46380+t46383+t46377-gamma_bb*t46385-rho_a*rho_b*t46368)*(2.0/3.0)-A*B*(C*C)*t46320*t46401*t46326*t46397*(1.0/9.0)+A*B*C*t46320*t46326*t46397*t46398*(2.6E1/9.0)+A*B*Dd*t46325*t46326*t46397*t46398*(2.6E1/9.0)-A*B*C*Dd*t46401*t46325*t46326*t46397*(2.0/9.0);
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t46403 = rho_a+rho_b;
                double t46404 = 1.0/pow(t46403,1.0/3.0);
                double t46405 = Dd*t46404;
                double t46406 = t46405+1.0;
                double t46407 = 1.0/t46406;
                double t46421 = C*t46404;
                double t46408 = exp(-t46421);
                double t46409 = C*t46404*(1.0/3.0);
                double t46410 = Dd*t46404*t46407*(1.0/3.0);
                double t46411 = 1.0/t46403;
                double t46412 = C*t46404*(1.0/9.0);
                double t46413 = Dd*t46404*t46407*(1.0/9.0);
                double t46414 = t46412+t46413-1.1E1/9.0;
                double t46415 = rho_a*t46411*t46414;
                double t46416 = t46410+t46415+t46409-1.0/9.0;
                double t46417 = 1.0/pow(t46403,4.0/3.0);
                double t46418 = Dd*Dd;
                double t46419 = 1.0/pow(t46403,5.0/3.0);
                double t46420 = 1.0/(t46406*t46406);
                double t46422 = rho_b*rho_b;
                double t46423 = rho_a*rho_b*t46416;
                double t46424 = t46422+t46423;
                double t46425 = 1.0/(t46403*t46403*t46403*t46403*t46403);
                v_rho_a_gamma_aa[Q] += scale * A*B*1.0/pow(t46403,1.4E1/3.0)*t46424*t46407*t46408*(-1.1E1/3.0)+A*B*1.0/pow(t46403,1.1E1/3.0)*t46407*t46408*(rho_b*t46416-rho_a*rho_b*(C*t46417*(1.0/9.0)-t46411*t46414+rho_a*t46411*(C*t46417*(1.0/2.7E1)+Dd*t46407*t46417*(1.0/2.7E1)-t46420*t46418*t46419*(1.0/2.7E1))+rho_a*1.0/(t46403*t46403)*t46414+Dd*t46407*t46417*(1.0/9.0)-t46420*t46418*t46419*(1.0/9.0)))+A*B*C*t46424*t46407*t46425*t46408*(1.0/3.0)+A*B*Dd*t46420*t46424*t46425*t46408*(1.0/3.0);
            }
            
            // v_rho_a_gamma_ab
            if (deriv >= 2) {
                double t46427 = rho_a+rho_b;
                double t46428 = 1.0/pow(t46427,1.0/3.0);
                double t46429 = Dd*t46428;
                double t46430 = t46429+1.0;
                double t46431 = 1.0/t46430;
                double t46437 = C*t46428;
                double t46432 = exp(-t46437);
                double t46433 = C*t46428*(7.0/9.0);
                double t46434 = Dd*t46431*t46428*(7.0/9.0);
                double t46435 = t46433+t46434-4.7E1/9.0;
                double t46436 = 1.0/pow(t46427,4.0/3.0);
                double t46438 = t46427*t46427;
                double t46439 = t46438*(4.0/3.0);
                double t46440 = rho_a*rho_b*t46435;
                double t46441 = t46440+t46439;
                double t46442 = 1.0/(t46427*t46427*t46427*t46427*t46427);
                double t46443 = 1.0/(t46430*t46430);
                v_rho_a_gamma_ab[Q] += scale * A*B*t46431*t46432*t46441*1.0/pow(t46427,1.4E1/3.0)*(-1.1E1/3.0)+A*B*t46431*t46432*1.0/pow(t46427,1.1E1/3.0)*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_b*t46435-rho_a*rho_b*(C*t46436*(7.0/2.7E1)-(Dd*Dd)*t46443*1.0/pow(t46427,5.0/3.0)*(7.0/2.7E1)+Dd*t46431*t46436*(7.0/2.7E1)))+A*B*C*t46431*t46432*t46441*t46442*(1.0/3.0)+A*B*Dd*t46432*t46441*t46442*t46443*(1.0/3.0);
            }
            
            // v_rho_a_gamma_bb
            if (deriv >= 2) {
                double t46445 = rho_a+rho_b;
                double t46446 = 1.0/pow(t46445,1.0/3.0);
                double t46447 = Dd*t46446;
                double t46448 = t46447+1.0;
                double t46449 = 1.0/t46448;
                double t46463 = C*t46446;
                double t46450 = exp(-t46463);
                double t46451 = C*t46446*(1.0/3.0);
                double t46452 = Dd*t46446*t46449*(1.0/3.0);
                double t46453 = 1.0/t46445;
                double t46454 = C*t46446*(1.0/9.0);
                double t46455 = Dd*t46446*t46449*(1.0/9.0);
                double t46456 = t46454+t46455-1.1E1/9.0;
                double t46457 = rho_b*t46453*t46456;
                double t46458 = t46451+t46452+t46457-1.0/9.0;
                double t46459 = 1.0/pow(t46445,4.0/3.0);
                double t46460 = Dd*Dd;
                double t46461 = 1.0/pow(t46445,5.0/3.0);
                double t46462 = 1.0/(t46448*t46448);
                double t46464 = rho_a*rho_a;
                double t46465 = rho_a*rho_b*t46458;
                double t46466 = t46464+t46465;
                double t46467 = 1.0/(t46445*t46445*t46445*t46445*t46445);
                v_rho_a_gamma_bb[Q] += scale * A*B*t46450*1.0/pow(t46445,1.4E1/3.0)*t46466*t46449*(-1.1E1/3.0)+A*B*t46450*1.0/pow(t46445,1.1E1/3.0)*t46449*(rho_a*2.0+rho_b*t46458-rho_a*rho_b*(C*t46459*(1.0/9.0)+rho_b*t46453*(C*t46459*(1.0/2.7E1)+Dd*t46449*t46459*(1.0/2.7E1)-t46460*t46461*t46462*(1.0/2.7E1))+rho_b*1.0/(t46445*t46445)*t46456+Dd*t46449*t46459*(1.0/9.0)-t46460*t46461*t46462*(1.0/9.0)))+A*B*C*t46450*t46466*t46449*t46467*(1.0/3.0)+A*B*Dd*t46450*t46462*t46466*t46467*(1.0/3.0);
            }
            
            // v_rho_b_gamma_aa
            if (deriv >= 2) {
                double t46469 = rho_a+rho_b;
                double t46470 = 1.0/pow(t46469,1.0/3.0);
                double t46471 = Dd*t46470;
                double t46472 = t46471+1.0;
                double t46473 = 1.0/t46472;
                double t46487 = C*t46470;
                double t46474 = exp(-t46487);
                double t46475 = C*t46470*(1.0/3.0);
                double t46476 = Dd*t46470*t46473*(1.0/3.0);
                double t46477 = 1.0/t46469;
                double t46478 = C*t46470*(1.0/9.0);
                double t46479 = Dd*t46470*t46473*(1.0/9.0);
                double t46480 = t46478+t46479-1.1E1/9.0;
                double t46481 = rho_a*t46480*t46477;
                double t46482 = t46481+t46475+t46476-1.0/9.0;
                double t46483 = 1.0/pow(t46469,4.0/3.0);
                double t46484 = Dd*Dd;
                double t46485 = 1.0/pow(t46469,5.0/3.0);
                double t46486 = 1.0/(t46472*t46472);
                double t46488 = rho_b*rho_b;
                double t46489 = rho_a*rho_b*t46482;
                double t46490 = t46488+t46489;
                double t46491 = 1.0/(t46469*t46469*t46469*t46469*t46469);
                v_rho_b_gamma_aa[Q] += scale * A*B*t46490*t46473*t46474*1.0/pow(t46469,1.4E1/3.0)*(-1.1E1/3.0)+A*B*t46473*t46474*1.0/pow(t46469,1.1E1/3.0)*(rho_b*2.0+rho_a*t46482-rho_a*rho_b*(C*t46483*(1.0/9.0)+rho_a*t46477*(C*t46483*(1.0/2.7E1)+Dd*t46473*t46483*(1.0/2.7E1)-t46484*t46485*t46486*(1.0/2.7E1))+rho_a*t46480*1.0/(t46469*t46469)+Dd*t46473*t46483*(1.0/9.0)-t46484*t46485*t46486*(1.0/9.0)))+A*B*C*t46490*t46473*t46491*t46474*(1.0/3.0)+A*B*Dd*t46490*t46491*t46474*t46486*(1.0/3.0);
            }
            
            // v_rho_b_gamma_ab
            if (deriv >= 2) {
                double t46493 = rho_a+rho_b;
                double t46494 = 1.0/pow(t46493,1.0/3.0);
                double t46495 = Dd*t46494;
                double t46496 = t46495+1.0;
                double t46497 = 1.0/t46496;
                double t46503 = C*t46494;
                double t46498 = exp(-t46503);
                double t46499 = C*t46494*(7.0/9.0);
                double t46500 = Dd*t46494*t46497*(7.0/9.0);
                double t46501 = t46500+t46499-4.7E1/9.0;
                double t46502 = 1.0/pow(t46493,4.0/3.0);
                double t46504 = t46493*t46493;
                double t46505 = t46504*(4.0/3.0);
                double t46506 = rho_a*rho_b*t46501;
                double t46507 = t46505+t46506;
                double t46508 = 1.0/(t46493*t46493*t46493*t46493*t46493);
                double t46509 = 1.0/(t46496*t46496);
                v_rho_b_gamma_ab[Q] += scale * A*B*t46507*1.0/pow(t46493,1.4E1/3.0)*t46497*t46498*(-1.1E1/3.0)+A*B*1.0/pow(t46493,1.1E1/3.0)*t46497*t46498*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_a*t46501-rho_a*rho_b*(C*t46502*(7.0/2.7E1)-(Dd*Dd)*t46509*1.0/pow(t46493,5.0/3.0)*(7.0/2.7E1)+Dd*t46502*t46497*(7.0/2.7E1)))+A*B*C*t46507*t46508*t46497*t46498*(1.0/3.0)+A*B*Dd*t46507*t46508*t46509*t46498*(1.0/3.0);
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t46511 = rho_a+rho_b;
                double t46512 = 1.0/pow(t46511,1.0/3.0);
                double t46513 = Dd*t46512;
                double t46514 = t46513+1.0;
                double t46515 = 1.0/t46514;
                double t46529 = C*t46512;
                double t46516 = exp(-t46529);
                double t46517 = C*t46512*(1.0/3.0);
                double t46518 = Dd*t46512*t46515*(1.0/3.0);
                double t46519 = 1.0/t46511;
                double t46520 = C*t46512*(1.0/9.0);
                double t46521 = Dd*t46512*t46515*(1.0/9.0);
                double t46522 = t46520+t46521-1.1E1/9.0;
                double t46523 = rho_b*t46522*t46519;
                double t46524 = t46523+t46517+t46518-1.0/9.0;
                double t46525 = 1.0/pow(t46511,4.0/3.0);
                double t46526 = Dd*Dd;
                double t46527 = 1.0/pow(t46511,5.0/3.0);
                double t46528 = 1.0/(t46514*t46514);
                double t46530 = rho_a*rho_a;
                double t46531 = rho_a*rho_b*t46524;
                double t46532 = t46530+t46531;
                double t46533 = 1.0/(t46511*t46511*t46511*t46511*t46511);
                v_rho_b_gamma_bb[Q] += scale * A*B*1.0/pow(t46511,1.4E1/3.0)*t46532*t46515*t46516*(-1.1E1/3.0)+A*B*1.0/pow(t46511,1.1E1/3.0)*t46515*t46516*(rho_a*t46524-rho_a*rho_b*(C*t46525*(1.0/9.0)-t46522*t46519+rho_b*t46519*(C*t46525*(1.0/2.7E1)+Dd*t46515*t46525*(1.0/2.7E1)-t46526*t46527*t46528*(1.0/2.7E1))+rho_b*1.0/(t46511*t46511)*t46522+Dd*t46515*t46525*(1.0/9.0)-t46526*t46527*t46528*(1.0/9.0)))+A*B*C*t46532*t46515*t46533*t46516*(1.0/3.0)+A*B*Dd*t46532*t46533*t46516*t46528*(1.0/3.0);
            }
            
        }
    }
}

}
