#include <libmints/vector.h>
#include "PZ81_Cfunctional.h"
#include "utility.h"
#include <cmath>

using namespace psi;

namespace psi {

PZ81_CFunctional::PZ81_CFunctional()
{
    name_ = "PZ81_C";
    description_ = "    PZ81 Correlation\n";
    citation_ = "    J.P. Perdew, A. Zunger, Phys. Rev. B., 23, 5048-5079, 1981\n";
    alpha_ = 1.0;
    omega_ = 0.0;
    lrc_ = false;
    gga_ = false;
    meta_ = false;
    parameters_["c"] =   6.2035049089939986E-01;
    parameters_["two_13"] =   1.2599210498948732E+00;
    parameters_["EcPld_1"] =  -1.4230000000000001E-01;
    parameters_["EcPld_2"] =   1.0528999999999999E+00;
    parameters_["EcPld_3"] =   3.3339999999999997E-01;
    parameters_["EcFld_1"] =  -8.4300000000000000E-02;
    parameters_["EcFld_2"] =   1.3980999999999999E+00;
    parameters_["EcFld_3"] =   2.6110000000000000E-01;
    parameters_["EcPhd_1"] =   3.1099999999999999E-02;
    parameters_["EcPhd_2"] =  -4.8000000000000001E-02;
    parameters_["EcPhd_3"] =   2.0000000000000000E-03;
    parameters_["EcPhd_4"] =  -1.1599999999999999E-02;
    parameters_["EcFhd_1"] =   1.5550000000000000E-02;
    parameters_["EcFhd_2"] =  -2.6900000000000000E-02;
    parameters_["EcFhd_3"] =   6.9999999999999999E-04;
    parameters_["EcFhd_4"] =  -4.7999999999999996E-03;
}
PZ81_CFunctional::~PZ81_CFunctional()
{
}
void PZ81_CFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    double c = parameters_["c"];
    double two_13 = parameters_["two_13"];
    double EcPld_1 = parameters_["EcPld_1"];
    double EcPld_2 = parameters_["EcPld_2"];
    double EcPld_3 = parameters_["EcPld_3"];
    double EcFld_1 = parameters_["EcFld_1"];
    double EcFld_2 = parameters_["EcFld_2"];
    double EcFld_3 = parameters_["EcFld_3"];
    double EcPhd_1 = parameters_["EcPhd_1"];
    double EcPhd_2 = parameters_["EcPhd_2"];
    double EcPhd_3 = parameters_["EcPhd_3"];
    double EcPhd_4 = parameters_["EcPhd_4"];
    double EcFhd_1 = parameters_["EcFhd_1"];
    double EcFhd_2 = parameters_["EcFhd_2"];
    double EcFhd_3 = parameters_["EcFhd_3"];
    double EcFhd_4 = parameters_["EcFhd_4"];

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
                double t39204 = rho_a+rho_b;
                double t39205 = 1.0/pow(t39204,1.0/3.0);
                double t39206 = c*t39205;
                double t39207 = log(t39206);
                double t39208 = EcPhd_1*t39207;
                double t39209 = pow(2.0,1.0/3.0);
                double t39210 = t39209*2.0;
                double t39211 = t39210-2.0;
                double t39212 = two_13*2.0;
                double t39213 = t39212-2.0;
                double t39214 = 1.0/t39213;
                double t39215 = sqrt(t39206);
                double t39216 = EcPld_2*t39215;
                double t39217 = EcPld_3*c*t39205;
                double t39218 = t39216+t39217+1.0;
                double t39219 = 1.0/t39218;
                double t39220 = EcPld_1*t39219;
                v[Q] += scale * t39204*(heaviside(-c*t39205+1.0)*(EcPhd_2+t39208+t39211*t39214*(EcFhd_2-EcPhd_2-t39208+EcFhd_1*t39207+EcFhd_4*c*t39205-EcPhd_4*c*t39205+EcFhd_3*c*t39205*t39207-EcPhd_3*c*t39205*t39207)+EcPhd_4*c*t39205+EcPhd_3*c*t39205*t39207)+heaviside(t39206-1.0)*(t39220-t39211*t39214*(t39220-EcFld_1/(EcFld_2*t39215+EcFld_3*c*t39205+1.0))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t39222 = rho_a+rho_b;
                double t39223 = 1.0/pow(t39222,4.0/3.0);
                double t39224 = 1.0/pow(t39222,1.0/3.0);
                double t39225 = c*t39224;
                double t39226 = 1.0/sqrt(t39225);
                double t39227 = sqrt(t39225);
                double t39228 = EcPld_3*c*t39223*(1.0/3.0);
                double t39229 = EcPld_2*c*t39223*t39226*(1.0/6.0);
                double t39230 = t39228+t39229;
                double t39231 = EcPld_2*t39227;
                double t39232 = EcPld_3*c*t39224;
                double t39233 = t39231+t39232+1.0;
                double t39234 = 1.0/(t39233*t39233);
                double t39235 = EcPld_1*t39230*t39234;
                double t39236 = pow(2.0,1.0/3.0);
                double t39237 = t39236*2.0;
                double t39238 = t39237-2.0;
                double t39239 = two_13*2.0;
                double t39240 = t39239-2.0;
                double t39241 = 1.0/t39240;
                double t39242 = 1.0/t39222;
                double t39243 = EcPhd_1*t39242*(1.0/3.0);
                double t39244 = log(t39225);
                double t39245 = EcPhd_3*c*t39223*(1.0/3.0);
                double t39246 = EcPhd_4*c*t39223*(1.0/3.0);
                double t39247 = EcPhd_3*c*t39223*t39244*(1.0/3.0);
                double t39248 = t39225-1.0;
                double t39249 = EcPhd_1*t39244;
                double t39250 = dirac(t39248);
                double t39251 = EcFld_2*t39227;
                double t39252 = EcFld_3*c*t39224;
                double t39253 = t39251+t39252+1.0;
                double t39254 = 1.0/t39233;
                double t39255 = EcPld_1*t39254;
                double t39256 = -t39225+1.0;
                double t39257 = heaviside(t39256);
                double t39258 = EcFhd_1*t39244;
                double t39259 = EcFhd_4*c*t39224;
                double t39260 = EcPhd_4*c*t39224;
                double t39261 = EcFhd_3*c*t39224*t39244;
                double t39262 = EcPhd_3*c*t39224*t39244;
                double t39263 = heaviside(t39248);
                double t39264 = 1.0/t39253;
                double t39265 = t39255-EcFld_1*t39264;
                double t39266 = t39255-t39241*t39238*t39265;
                v_rho_a[Q] += scale * t39257*(EcPhd_2+t39260+t39262+t39249+t39241*t39238*(EcFhd_2-EcPhd_2-t39260+t39261-t39262-t39249+t39258+t39259))+t39263*t39266+t39222*(t39263*(t39235-t39241*t39238*(t39235-EcFld_1*1.0/(t39253*t39253)*(EcFld_3*c*t39223*(1.0/3.0)+EcFld_2*c*t39223*t39226*(1.0/6.0))))-t39257*(t39243+t39245+t39246+t39247-t39241*t39238*(t39243+t39245+t39246+t39247-EcFhd_1*t39242*(1.0/3.0)-EcFhd_3*c*t39223*(1.0/3.0)-EcFhd_4*c*t39223*(1.0/3.0)-EcFhd_3*c*t39223*t39244*(1.0/3.0)))-c*t39223*t39250*t39266*(1.0/3.0)+c*t39223*t39250*(EcPhd_2+t39260+t39262+t39249+t39241*t39238*(EcFhd_2-EcPhd_2+t39261-t39249+t39258+t39259-EcPhd_4*c*t39224-EcPhd_3*c*t39224*t39244))*(1.0/3.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t39268 = rho_a+rho_b;
                double t39269 = 1.0/pow(t39268,4.0/3.0);
                double t39270 = 1.0/pow(t39268,1.0/3.0);
                double t39271 = c*t39270;
                double t39272 = 1.0/sqrt(t39271);
                double t39273 = sqrt(t39271);
                double t39274 = EcPld_3*c*t39269*(1.0/3.0);
                double t39275 = EcPld_2*c*t39272*t39269*(1.0/6.0);
                double t39276 = t39274+t39275;
                double t39277 = EcPld_2*t39273;
                double t39278 = EcPld_3*c*t39270;
                double t39279 = t39277+t39278+1.0;
                double t39280 = 1.0/(t39279*t39279);
                double t39281 = EcPld_1*t39280*t39276;
                double t39282 = pow(2.0,1.0/3.0);
                double t39283 = t39282*2.0;
                double t39284 = t39283-2.0;
                double t39285 = two_13*2.0;
                double t39286 = t39285-2.0;
                double t39287 = 1.0/t39286;
                double t39288 = 1.0/t39268;
                double t39289 = EcPhd_1*t39288*(1.0/3.0);
                double t39290 = log(t39271);
                double t39291 = EcPhd_3*c*t39269*(1.0/3.0);
                double t39292 = EcPhd_4*c*t39269*(1.0/3.0);
                double t39293 = EcPhd_3*c*t39290*t39269*(1.0/3.0);
                double t39294 = t39271-1.0;
                double t39295 = EcPhd_1*t39290;
                double t39296 = dirac(t39294);
                double t39297 = EcFld_2*t39273;
                double t39298 = EcFld_3*c*t39270;
                double t39299 = t39297+t39298+1.0;
                double t39300 = 1.0/t39279;
                double t39301 = EcPld_1*t39300;
                double t39302 = -t39271+1.0;
                double t39303 = heaviside(t39302);
                double t39304 = EcFhd_1*t39290;
                double t39305 = EcFhd_4*c*t39270;
                double t39306 = EcPhd_4*c*t39270;
                double t39307 = EcFhd_3*c*t39270*t39290;
                double t39308 = EcPhd_3*c*t39270*t39290;
                double t39309 = heaviside(t39294);
                double t39310 = 1.0/t39299;
                double t39311 = t39301-EcFld_1*t39310;
                double t39312 = t39301-t39311*t39284*t39287;
                v_rho_b[Q] += scale * t39303*(EcPhd_2+t39306+t39308+t39295+t39284*t39287*(EcFhd_2-EcPhd_2+t39304+t39305-t39306+t39307-t39308-t39295))+t39312*t39309+t39268*(t39309*(t39281-t39284*t39287*(t39281-EcFld_1*1.0/(t39299*t39299)*(EcFld_3*c*t39269*(1.0/3.0)+EcFld_2*c*t39272*t39269*(1.0/6.0))))-t39303*(t39291+t39292+t39293+t39289-t39284*t39287*(t39291+t39292+t39293+t39289-EcFhd_1*t39288*(1.0/3.0)-EcFhd_3*c*t39269*(1.0/3.0)-EcFhd_4*c*t39269*(1.0/3.0)-EcFhd_3*c*t39290*t39269*(1.0/3.0)))-c*t39312*t39269*t39296*(1.0/3.0)+c*t39269*t39296*(EcPhd_2+t39306+t39308+t39295+t39284*t39287*(EcFhd_2-EcPhd_2+t39304+t39305+t39307-t39295-EcPhd_4*c*t39270-EcPhd_3*c*t39270*t39290))*(1.0/3.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t39319 = rho_a+rho_b;
                double t39320 = 1.0/pow(t39319,4.0/3.0);
                double t39321 = 1.0/pow(t39319,1.0/3.0);
                double t39322 = c*t39321;
                double t39323 = 1.0/sqrt(t39322);
                double t39324 = sqrt(t39322);
                double t39325 = EcPld_3*c*t39320*(1.0/3.0);
                double t39326 = EcPld_2*c*t39320*t39323*(1.0/6.0);
                double t39327 = t39325+t39326;
                double t39328 = EcPld_2*t39324;
                double t39329 = EcPld_3*c*t39321;
                double t39330 = t39328+t39329+1.0;
                double t39331 = 1.0/(t39330*t39330);
                double t39332 = EcPld_1*t39331*t39327;
                double t39333 = t39322-1.0;
                double t39334 = heaviside(t39333);
                double t39335 = pow(2.0,1.0/3.0);
                double t39336 = t39335*2.0;
                double t39337 = t39336-2.0;
                double t39338 = two_13*2.0;
                double t39339 = t39338-2.0;
                double t39340 = 1.0/t39339;
                double t39341 = EcFld_3*c*t39320*(1.0/3.0);
                double t39342 = EcFld_2*c*t39320*t39323*(1.0/6.0);
                double t39343 = t39341+t39342;
                double t39344 = EcFld_2*t39324;
                double t39345 = EcFld_3*c*t39321;
                double t39346 = t39344+t39345+1.0;
                double t39347 = t39327*t39327;
                double t39348 = 1.0/(t39330*t39330*t39330);
                double t39349 = EcPld_1*t39347*t39348*2.0;
                double t39350 = 1.0/pow(t39319,7.0/3.0);
                double t39351 = 1.0/(t39346*t39346);
                double t39352 = c*c;
                double t39353 = 1.0/pow(t39319,8.0/3.0);
                double t39354 = 1.0/pow(t39322,3.0/2.0);
                double t39355 = EcPld_3*c*t39350*(4.0/9.0);
                double t39356 = EcPld_2*c*t39323*t39350*(2.0/9.0);
                double t39357 = t39355+t39356-EcPld_2*t39352*t39353*t39354*(1.0/3.6E1);
                double t39358 = EcPld_1*t39331*t39357;
                double t39359 = 1.0/(t39319*t39319);
                double t39360 = EcPhd_1*t39359*(1.0/3.0);
                double t39361 = log(t39322);
                double t39362 = EcPhd_3*c*t39350*(5.0/9.0);
                double t39363 = EcPhd_4*c*t39350*(4.0/9.0);
                double t39364 = EcPhd_3*c*t39350*t39361*(4.0/9.0);
                double t39365 = 1.0/t39330;
                double t39366 = EcPld_1*t39365;
                double t39367 = t39332-EcFld_1*t39351*t39343;
                double t39368 = t39332-t39340*t39337*t39367;
                double t39369 = dirac(t39333);
                double t39370 = EcPhd_1*t39361;
                double t39371 = 1.0/t39319;
                double t39372 = EcPhd_1*t39371*(1.0/3.0);
                double t39373 = EcPhd_3*c*t39320*(1.0/3.0);
                double t39374 = EcPhd_4*c*t39320*(1.0/3.0);
                double t39375 = EcPhd_3*c*t39320*t39361*(1.0/3.0);
                double t39376 = 1.0/t39346;
                double t39392 = EcFld_1*t39376;
                double t39377 = -t39392+t39366;
                double t39378 = t39366-t39340*t39337*t39377;
                double t39379 = dirac(t39333,1.0);
                double t39380 = EcFhd_1*t39361;
                double t39381 = EcFhd_4*c*t39321;
                double t39382 = EcPhd_4*c*t39321;
                double t39383 = EcFhd_3*c*t39321*t39361;
                double t39384 = EcPhd_3*c*t39321*t39361;
                double t39385 = EcFhd_2-EcPhd_2-t39370+t39380+t39381-t39382+t39383-t39384;
                double t39386 = t39340*t39337*t39385;
                double t39387 = EcPhd_2+t39370+t39382+t39384+t39386;
                double t39388 = -t39322+1.0;
                double t39389 = heaviside(t39388);
                double t39390 = t39372+t39373+t39374+t39375-EcFhd_1*t39371*(1.0/3.0)-EcFhd_3*c*t39320*(1.0/3.0)-EcFhd_4*c*t39320*(1.0/3.0)-EcFhd_3*c*t39320*t39361*(1.0/3.0);
                double t39391 = t39372+t39373+t39374+t39375-t39340*t39390*t39337;
                v_rho_a_rho_a[Q] += scale * -t39319*(-t39389*(t39360+t39362+t39363+t39364-t39340*t39337*(t39360+t39362+t39363+t39364-EcFhd_1*t39359*(1.0/3.0)-EcFhd_3*c*t39350*(5.0/9.0)-EcFhd_4*c*t39350*(4.0/9.0)-EcFhd_3*c*t39350*t39361*(4.0/9.0)))+t39334*(-t39349+t39358+t39340*t39337*(t39349-t39358-EcFld_1*(t39343*t39343)*1.0/(t39346*t39346*t39346)*2.0+EcFld_1*t39351*(EcFld_3*c*t39350*(4.0/9.0)-EcFld_2*t39352*t39353*t39354*(1.0/3.6E1)+EcFld_2*c*t39323*t39350*(2.0/9.0))))+c*t39320*t39391*t39369*(2.0/3.0)+c*t39320*t39368*t39369*(2.0/3.0)-c*t39350*t39369*t39378*(4.0/9.0)+c*t39350*t39369*t39387*(4.0/9.0)-t39352*t39353*t39378*t39379*(1.0/9.0)+t39352*t39353*t39387*t39379*(1.0/9.0))+t39334*t39368*2.0-t39391*t39389*2.0+c*t39320*t39369*t39387*(2.0/3.0)-c*t39320*t39369*(t39366+t39340*t39337*(t39392-t39366))*(2.0/3.0);
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t39394 = rho_a+rho_b;
                double t39395 = 1.0/pow(t39394,4.0/3.0);
                double t39396 = 1.0/pow(t39394,1.0/3.0);
                double t39397 = c*t39396;
                double t39398 = 1.0/sqrt(t39397);
                double t39399 = sqrt(t39397);
                double t39400 = EcPld_3*c*t39395*(1.0/3.0);
                double t39401 = EcPld_2*c*t39395*t39398*(1.0/6.0);
                double t39402 = t39400+t39401;
                double t39403 = EcPld_2*t39399;
                double t39404 = EcPld_3*c*t39396;
                double t39405 = t39403+t39404+1.0;
                double t39406 = 1.0/(t39405*t39405);
                double t39407 = EcPld_1*t39402*t39406;
                double t39408 = t39397-1.0;
                double t39409 = heaviside(t39408);
                double t39410 = pow(2.0,1.0/3.0);
                double t39411 = t39410*2.0;
                double t39412 = t39411-2.0;
                double t39413 = two_13*2.0;
                double t39414 = t39413-2.0;
                double t39415 = 1.0/t39414;
                double t39416 = EcFld_3*c*t39395*(1.0/3.0);
                double t39417 = EcFld_2*c*t39395*t39398*(1.0/6.0);
                double t39418 = t39416+t39417;
                double t39419 = EcFld_2*t39399;
                double t39420 = EcFld_3*c*t39396;
                double t39421 = t39420+t39419+1.0;
                double t39422 = t39402*t39402;
                double t39423 = 1.0/(t39405*t39405*t39405);
                double t39424 = EcPld_1*t39422*t39423*2.0;
                double t39425 = 1.0/pow(t39394,7.0/3.0);
                double t39426 = 1.0/(t39421*t39421);
                double t39427 = c*c;
                double t39428 = 1.0/pow(t39394,8.0/3.0);
                double t39429 = 1.0/pow(t39397,3.0/2.0);
                double t39430 = EcPld_3*c*t39425*(4.0/9.0);
                double t39431 = EcPld_2*c*t39425*t39398*(2.0/9.0);
                double t39432 = t39430+t39431-EcPld_2*t39427*t39428*t39429*(1.0/3.6E1);
                double t39433 = EcPld_1*t39432*t39406;
                double t39434 = 1.0/(t39394*t39394);
                double t39435 = EcPhd_1*t39434*(1.0/3.0);
                double t39436 = log(t39397);
                double t39437 = EcPhd_3*c*t39425*(5.0/9.0);
                double t39438 = EcPhd_4*c*t39425*(4.0/9.0);
                double t39439 = EcPhd_3*c*t39425*t39436*(4.0/9.0);
                double t39440 = 1.0/t39405;
                double t39441 = EcPld_1*t39440;
                double t39442 = t39407-EcFld_1*t39426*t39418;
                double t39443 = t39407-t39412*t39415*t39442;
                double t39444 = dirac(t39408);
                double t39445 = EcPhd_1*t39436;
                double t39446 = 1.0/t39394;
                double t39447 = EcPhd_1*t39446*(1.0/3.0);
                double t39448 = EcPhd_3*c*t39395*(1.0/3.0);
                double t39449 = EcPhd_4*c*t39395*(1.0/3.0);
                double t39450 = EcPhd_3*c*t39436*t39395*(1.0/3.0);
                double t39451 = 1.0/t39421;
                double t39467 = EcFld_1*t39451;
                double t39452 = t39441-t39467;
                double t39468 = t39412*t39415*t39452;
                double t39453 = t39441-t39468;
                double t39454 = dirac(t39408,1.0);
                double t39455 = EcFhd_1*t39436;
                double t39456 = EcFhd_4*c*t39396;
                double t39457 = EcPhd_4*c*t39396;
                double t39458 = EcFhd_3*c*t39436*t39396;
                double t39459 = EcPhd_3*c*t39436*t39396;
                double t39460 = EcFhd_2-EcPhd_2-t39445+t39455+t39456-t39457+t39458-t39459;
                double t39461 = t39412*t39415*t39460;
                double t39462 = EcPhd_2+t39461+t39445+t39457+t39459;
                double t39463 = -t39397+1.0;
                double t39464 = heaviside(t39463);
                double t39465 = t39450+t39447+t39448+t39449-EcFhd_1*t39446*(1.0/3.0)-EcFhd_3*c*t39395*(1.0/3.0)-EcFhd_4*c*t39395*(1.0/3.0)-EcFhd_3*c*t39436*t39395*(1.0/3.0);
                double t39466 = t39450+t39447+t39448+t39449-t39412*t39415*t39465;
                v_rho_a_rho_b[Q] += scale * -t39394*(-t39464*(t39435+t39437+t39438+t39439-t39412*t39415*(t39435+t39437+t39438+t39439-EcFhd_1*t39434*(1.0/3.0)-EcFhd_3*c*t39425*(5.0/9.0)-EcFhd_4*c*t39425*(4.0/9.0)-EcFhd_3*c*t39425*t39436*(4.0/9.0)))+t39409*(-t39424+t39433+t39412*t39415*(t39424-t39433-EcFld_1*1.0/(t39421*t39421*t39421)*(t39418*t39418)*2.0+EcFld_1*t39426*(EcFld_3*c*t39425*(4.0/9.0)-EcFld_2*t39427*t39428*t39429*(1.0/3.6E1)+EcFld_2*c*t39425*t39398*(2.0/9.0))))-c*t39425*t39444*t39453*(4.0/9.0)+c*t39425*t39444*t39462*(4.0/9.0)+c*t39443*t39444*t39395*(2.0/3.0)+c*t39444*t39466*t39395*(2.0/3.0)-t39453*t39427*t39454*t39428*(1.0/9.0)+t39462*t39427*t39454*t39428*(1.0/9.0))+t39443*t39409*2.0-t39464*t39466*2.0-c*t39444*t39453*t39395*(2.0/3.0)+c*t39444*t39462*t39395*(2.0/3.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t39470 = rho_a+rho_b;
                double t39471 = 1.0/pow(t39470,4.0/3.0);
                double t39472 = 1.0/pow(t39470,1.0/3.0);
                double t39473 = c*t39472;
                double t39474 = 1.0/sqrt(t39473);
                double t39475 = sqrt(t39473);
                double t39476 = EcPld_3*c*t39471*(1.0/3.0);
                double t39477 = EcPld_2*c*t39471*t39474*(1.0/6.0);
                double t39478 = t39476+t39477;
                double t39479 = EcPld_2*t39475;
                double t39480 = EcPld_3*c*t39472;
                double t39481 = t39480+t39479+1.0;
                double t39482 = 1.0/(t39481*t39481);
                double t39483 = EcPld_1*t39482*t39478;
                double t39484 = t39473-1.0;
                double t39485 = heaviside(t39484);
                double t39486 = pow(2.0,1.0/3.0);
                double t39487 = t39486*2.0;
                double t39488 = t39487-2.0;
                double t39489 = two_13*2.0;
                double t39490 = t39489-2.0;
                double t39491 = 1.0/t39490;
                double t39492 = EcFld_3*c*t39471*(1.0/3.0);
                double t39493 = EcFld_2*c*t39471*t39474*(1.0/6.0);
                double t39494 = t39492+t39493;
                double t39495 = EcFld_2*t39475;
                double t39496 = EcFld_3*c*t39472;
                double t39497 = t39495+t39496+1.0;
                double t39498 = t39478*t39478;
                double t39499 = 1.0/(t39481*t39481*t39481);
                double t39500 = EcPld_1*t39498*t39499*2.0;
                double t39501 = 1.0/pow(t39470,7.0/3.0);
                double t39502 = 1.0/(t39497*t39497);
                double t39503 = c*c;
                double t39504 = 1.0/pow(t39470,8.0/3.0);
                double t39505 = 1.0/pow(t39473,3.0/2.0);
                double t39506 = EcPld_3*c*t39501*(4.0/9.0);
                double t39507 = EcPld_2*c*t39501*t39474*(2.0/9.0);
                double t39508 = t39506+t39507-EcPld_2*t39503*t39504*t39505*(1.0/3.6E1);
                double t39509 = EcPld_1*t39508*t39482;
                double t39510 = 1.0/(t39470*t39470);
                double t39511 = EcPhd_1*t39510*(1.0/3.0);
                double t39512 = log(t39473);
                double t39513 = EcPhd_3*c*t39501*(5.0/9.0);
                double t39514 = EcPhd_4*c*t39501*(4.0/9.0);
                double t39515 = EcPhd_3*c*t39501*t39512*(4.0/9.0);
                double t39516 = 1.0/t39481;
                double t39517 = EcPld_1*t39516;
                double t39518 = t39483-EcFld_1*t39502*t39494;
                double t39519 = t39483-t39491*t39518*t39488;
                double t39520 = dirac(t39484);
                double t39521 = EcPhd_1*t39512;
                double t39522 = 1.0/t39470;
                double t39523 = EcPhd_1*t39522*(1.0/3.0);
                double t39524 = EcPhd_3*c*t39471*(1.0/3.0);
                double t39525 = EcPhd_4*c*t39471*(1.0/3.0);
                double t39526 = EcPhd_3*c*t39512*t39471*(1.0/3.0);
                double t39527 = 1.0/t39497;
                double t39543 = EcFld_1*t39527;
                double t39528 = -t39543+t39517;
                double t39529 = t39517-t39491*t39528*t39488;
                double t39530 = dirac(t39484,1.0);
                double t39531 = EcFhd_1*t39512;
                double t39532 = EcFhd_4*c*t39472;
                double t39533 = EcPhd_4*c*t39472;
                double t39534 = EcFhd_3*c*t39512*t39472;
                double t39535 = EcPhd_3*c*t39512*t39472;
                double t39536 = EcFhd_2-EcPhd_2-t39521+t39531+t39532-t39533+t39534-t39535;
                double t39537 = t39491*t39536*t39488;
                double t39538 = EcPhd_2+t39521+t39533+t39535+t39537;
                double t39539 = -t39473+1.0;
                double t39540 = heaviside(t39539);
                double t39541 = t39523+t39524+t39525+t39526-EcFhd_1*t39522*(1.0/3.0)-EcFhd_3*c*t39471*(1.0/3.0)-EcFhd_4*c*t39471*(1.0/3.0)-EcFhd_3*c*t39512*t39471*(1.0/3.0);
                double t39542 = t39523+t39524+t39525+t39526-t39541*t39491*t39488;
                v_rho_b_rho_b[Q] += scale * -t39470*(-t39540*(t39511+t39513+t39514+t39515-t39491*t39488*(t39511+t39513+t39514+t39515-EcFhd_1*t39510*(1.0/3.0)-EcFhd_3*c*t39501*(5.0/9.0)-EcFhd_4*c*t39501*(4.0/9.0)-EcFhd_3*c*t39501*t39512*(4.0/9.0)))+t39485*(-t39500+t39509+t39491*t39488*(t39500-t39509-EcFld_1*(t39494*t39494)*1.0/(t39497*t39497*t39497)*2.0+EcFld_1*t39502*(EcFld_3*c*t39501*(4.0/9.0)-EcFld_2*t39503*t39504*t39505*(1.0/3.6E1)+EcFld_2*c*t39501*t39474*(2.0/9.0))))-c*t39501*t39520*t39529*(4.0/9.0)+c*t39501*t39520*t39538*(4.0/9.0)+c*t39520*t39542*t39471*(2.0/3.0)+c*t39520*t39471*t39519*(2.0/3.0)-t39503*t39530*t39504*t39529*(1.0/9.0)+t39503*t39530*t39504*t39538*(1.0/9.0))-t39540*t39542*2.0+t39519*t39485*2.0+c*t39520*t39471*t39538*(2.0/3.0)-c*t39520*t39471*(t39517+t39491*t39488*(t39543-t39517))*(2.0/3.0);
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t39570 = rho_a+rho_b;
                double t39571 = 1.0/pow(t39570,1.0/3.0);
                double t39572 = c*t39571;
                double t39573 = log(t39572);
                double t39574 = EcPhd_1*t39573;
                double t39575 = pow(2.0,1.0/3.0);
                double t39576 = t39575*2.0;
                double t39577 = t39576-2.0;
                double t39578 = two_13*2.0;
                double t39579 = t39578-2.0;
                double t39580 = 1.0/t39579;
                double t39581 = sqrt(t39572);
                double t39582 = EcPld_2*t39581;
                double t39583 = EcPld_3*c*t39571;
                double t39584 = t39582+t39583+1.0;
                double t39585 = 1.0/t39584;
                double t39586 = EcPld_1*t39585;
                v[Q] += scale * t39570*(heaviside(-c*t39571+1.0)*(EcPhd_2+t39574+t39580*t39577*(EcFhd_2-EcPhd_2-t39574+EcFhd_1*t39573+EcFhd_4*c*t39571-EcPhd_4*c*t39571+EcFhd_3*c*t39571*t39573-EcPhd_3*c*t39571*t39573)+EcPhd_4*c*t39571+EcPhd_3*c*t39571*t39573)+heaviside(t39572-1.0)*(t39586-t39580*t39577*(t39586-EcFld_1/(EcFld_2*t39581+EcFld_3*c*t39571+1.0))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t39588 = rho_a+rho_b;
                double t39589 = 1.0/pow(t39588,4.0/3.0);
                double t39590 = 1.0/pow(t39588,1.0/3.0);
                double t39591 = c*t39590;
                double t39592 = 1.0/sqrt(t39591);
                double t39593 = sqrt(t39591);
                double t39594 = EcPld_3*c*t39589*(1.0/3.0);
                double t39595 = EcPld_2*c*t39592*t39589*(1.0/6.0);
                double t39596 = t39594+t39595;
                double t39597 = EcPld_2*t39593;
                double t39598 = EcPld_3*c*t39590;
                double t39599 = t39597+t39598+1.0;
                double t39600 = 1.0/(t39599*t39599);
                double t39601 = EcPld_1*t39600*t39596;
                double t39602 = pow(2.0,1.0/3.0);
                double t39603 = t39602*2.0;
                double t39604 = t39603-2.0;
                double t39605 = two_13*2.0;
                double t39606 = t39605-2.0;
                double t39607 = 1.0/t39606;
                double t39608 = 1.0/t39588;
                double t39609 = EcPhd_1*t39608*(1.0/3.0);
                double t39610 = log(t39591);
                double t39611 = EcPhd_3*c*t39589*(1.0/3.0);
                double t39612 = EcPhd_4*c*t39589*(1.0/3.0);
                double t39613 = EcPhd_3*c*t39610*t39589*(1.0/3.0);
                double t39614 = t39591-1.0;
                double t39615 = EcPhd_1*t39610;
                double t39616 = dirac(t39614);
                double t39617 = EcFld_2*t39593;
                double t39618 = EcFld_3*c*t39590;
                double t39619 = t39617+t39618+1.0;
                double t39620 = 1.0/t39599;
                double t39621 = EcPld_1*t39620;
                double t39622 = -t39591+1.0;
                double t39623 = heaviside(t39622);
                double t39624 = EcFhd_1*t39610;
                double t39625 = EcFhd_4*c*t39590;
                double t39626 = EcPhd_4*c*t39590;
                double t39627 = EcFhd_3*c*t39610*t39590;
                double t39628 = EcPhd_3*c*t39610*t39590;
                double t39629 = heaviside(t39614);
                double t39630 = 1.0/t39619;
                double t39631 = t39621-EcFld_1*t39630;
                double t39632 = t39621-t39604*t39631*t39607;
                v_rho_a[Q] += scale * t39623*(EcPhd_2+t39615+t39626+t39628+t39604*t39607*(EcFhd_2-EcPhd_2-t39615+t39624+t39625-t39626+t39627-t39628))+t39632*t39629+t39588*(t39629*(t39601-t39604*t39607*(t39601-EcFld_1*1.0/(t39619*t39619)*(EcFld_3*c*t39589*(1.0/3.0)+EcFld_2*c*t39592*t39589*(1.0/6.0))))-t39623*(t39611+t39612+t39613+t39609-t39604*t39607*(t39611+t39612+t39613+t39609-EcFhd_1*t39608*(1.0/3.0)-EcFhd_3*c*t39589*(1.0/3.0)-EcFhd_4*c*t39589*(1.0/3.0)-EcFhd_3*c*t39610*t39589*(1.0/3.0)))-c*t39632*t39616*t39589*(1.0/3.0)+c*t39616*t39589*(EcPhd_2+t39615+t39626+t39628+t39604*t39607*(EcFhd_2-EcPhd_2-t39615+t39624+t39625+t39627-EcPhd_4*c*t39590-EcPhd_3*c*t39610*t39590))*(1.0/3.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t39634 = rho_a+rho_b;
                double t39635 = 1.0/pow(t39634,4.0/3.0);
                double t39636 = 1.0/pow(t39634,1.0/3.0);
                double t39637 = c*t39636;
                double t39638 = 1.0/sqrt(t39637);
                double t39639 = sqrt(t39637);
                double t39640 = EcPld_3*c*t39635*(1.0/3.0);
                double t39641 = EcPld_2*c*t39635*t39638*(1.0/6.0);
                double t39642 = t39640+t39641;
                double t39643 = EcPld_2*t39639;
                double t39644 = EcPld_3*c*t39636;
                double t39645 = t39643+t39644+1.0;
                double t39646 = 1.0/(t39645*t39645);
                double t39647 = EcPld_1*t39642*t39646;
                double t39648 = pow(2.0,1.0/3.0);
                double t39649 = t39648*2.0;
                double t39650 = t39649-2.0;
                double t39651 = two_13*2.0;
                double t39652 = t39651-2.0;
                double t39653 = 1.0/t39652;
                double t39654 = 1.0/t39634;
                double t39655 = EcPhd_1*t39654*(1.0/3.0);
                double t39656 = log(t39637);
                double t39657 = EcPhd_3*c*t39635*(1.0/3.0);
                double t39658 = EcPhd_4*c*t39635*(1.0/3.0);
                double t39659 = EcPhd_3*c*t39635*t39656*(1.0/3.0);
                double t39660 = t39637-1.0;
                double t39661 = EcPhd_1*t39656;
                double t39662 = dirac(t39660);
                double t39663 = EcFld_2*t39639;
                double t39664 = EcFld_3*c*t39636;
                double t39665 = t39663+t39664+1.0;
                double t39666 = 1.0/t39645;
                double t39667 = EcPld_1*t39666;
                double t39668 = -t39637+1.0;
                double t39669 = heaviside(t39668);
                double t39670 = EcFhd_1*t39656;
                double t39671 = EcFhd_4*c*t39636;
                double t39672 = EcPhd_4*c*t39636;
                double t39673 = EcFhd_3*c*t39636*t39656;
                double t39674 = EcPhd_3*c*t39636*t39656;
                double t39675 = heaviside(t39660);
                double t39676 = 1.0/t39665;
                double t39677 = t39667-EcFld_1*t39676;
                double t39678 = t39667-t39650*t39653*t39677;
                v_rho_b[Q] += scale * t39669*(EcPhd_2+t39661+t39672+t39674+t39650*t39653*(EcFhd_2-EcPhd_2-t39661+t39670+t39671-t39672+t39673-t39674))+t39675*t39678+t39634*(t39675*(t39647-t39650*t39653*(t39647-EcFld_1*1.0/(t39665*t39665)*(EcFld_3*c*t39635*(1.0/3.0)+EcFld_2*c*t39635*t39638*(1.0/6.0))))-t39669*(t39655+t39657+t39658+t39659-t39650*t39653*(t39655+t39657+t39658+t39659-EcFhd_1*t39654*(1.0/3.0)-EcFhd_3*c*t39635*(1.0/3.0)-EcFhd_4*c*t39635*(1.0/3.0)-EcFhd_3*c*t39635*t39656*(1.0/3.0)))-c*t39635*t39662*t39678*(1.0/3.0)+c*t39635*t39662*(EcPhd_2+t39661+t39672+t39674+t39650*t39653*(EcFhd_2-EcPhd_2-t39661+t39670+t39671+t39673-EcPhd_4*c*t39636-EcPhd_3*c*t39636*t39656))*(1.0/3.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t39685 = rho_a+rho_b;
                double t39686 = 1.0/pow(t39685,4.0/3.0);
                double t39687 = 1.0/pow(t39685,1.0/3.0);
                double t39688 = c*t39687;
                double t39689 = 1.0/sqrt(t39688);
                double t39690 = sqrt(t39688);
                double t39691 = EcPld_3*c*t39686*(1.0/3.0);
                double t39692 = EcPld_2*c*t39686*t39689*(1.0/6.0);
                double t39693 = t39691+t39692;
                double t39694 = EcPld_2*t39690;
                double t39695 = EcPld_3*c*t39687;
                double t39696 = t39694+t39695+1.0;
                double t39697 = 1.0/(t39696*t39696);
                double t39698 = EcPld_1*t39693*t39697;
                double t39699 = t39688-1.0;
                double t39700 = heaviside(t39699);
                double t39701 = pow(2.0,1.0/3.0);
                double t39702 = t39701*2.0;
                double t39703 = t39702-2.0;
                double t39704 = two_13*2.0;
                double t39705 = t39704-2.0;
                double t39706 = 1.0/t39705;
                double t39707 = EcFld_3*c*t39686*(1.0/3.0);
                double t39708 = EcFld_2*c*t39686*t39689*(1.0/6.0);
                double t39709 = t39707+t39708;
                double t39710 = EcFld_2*t39690;
                double t39711 = EcFld_3*c*t39687;
                double t39712 = t39710+t39711+1.0;
                double t39713 = t39693*t39693;
                double t39714 = 1.0/(t39696*t39696*t39696);
                double t39715 = EcPld_1*t39713*t39714*2.0;
                double t39716 = 1.0/pow(t39685,7.0/3.0);
                double t39717 = 1.0/(t39712*t39712);
                double t39718 = c*c;
                double t39719 = 1.0/pow(t39685,8.0/3.0);
                double t39720 = 1.0/pow(t39688,3.0/2.0);
                double t39721 = EcPld_3*c*t39716*(4.0/9.0);
                double t39722 = EcPld_2*c*t39716*t39689*(2.0/9.0);
                double t39723 = t39721+t39722-EcPld_2*t39720*t39718*t39719*(1.0/3.6E1);
                double t39724 = EcPld_1*t39723*t39697;
                double t39725 = 1.0/(t39685*t39685);
                double t39726 = EcPhd_1*t39725*(1.0/3.0);
                double t39727 = log(t39688);
                double t39728 = EcPhd_3*c*t39716*(5.0/9.0);
                double t39729 = EcPhd_4*c*t39716*(4.0/9.0);
                double t39730 = EcPhd_3*c*t39716*t39727*(4.0/9.0);
                double t39731 = 1.0/t39696;
                double t39732 = EcPld_1*t39731;
                double t39733 = t39698-EcFld_1*t39717*t39709;
                double t39734 = t39698-t39703*t39706*t39733;
                double t39735 = dirac(t39699);
                double t39736 = EcPhd_1*t39727;
                double t39737 = 1.0/t39685;
                double t39738 = EcPhd_1*t39737*(1.0/3.0);
                double t39739 = EcPhd_3*c*t39686*(1.0/3.0);
                double t39740 = EcPhd_4*c*t39686*(1.0/3.0);
                double t39741 = EcPhd_3*c*t39727*t39686*(1.0/3.0);
                double t39742 = 1.0/t39712;
                double t39758 = EcFld_1*t39742;
                double t39743 = t39732-t39758;
                double t39759 = t39703*t39706*t39743;
                double t39744 = t39732-t39759;
                double t39745 = dirac(t39699,1.0);
                double t39746 = EcFhd_1*t39727;
                double t39747 = EcFhd_4*c*t39687;
                double t39748 = EcPhd_4*c*t39687;
                double t39749 = EcFhd_3*c*t39727*t39687;
                double t39750 = EcPhd_3*c*t39727*t39687;
                double t39751 = EcFhd_2-EcPhd_2-t39750-t39736+t39746+t39747-t39748+t39749;
                double t39752 = t39703*t39706*t39751;
                double t39753 = EcPhd_2+t39750+t39752+t39736+t39748;
                double t39754 = -t39688+1.0;
                double t39755 = heaviside(t39754);
                double t39756 = t39740+t39741+t39738+t39739-EcFhd_1*t39737*(1.0/3.0)-EcFhd_3*c*t39686*(1.0/3.0)-EcFhd_4*c*t39686*(1.0/3.0)-EcFhd_3*c*t39727*t39686*(1.0/3.0);
                double t39757 = t39740+t39741+t39738+t39739-t39703*t39706*t39756;
                v_rho_a_rho_a[Q] += scale * -t39685*(-t39755*(t39730+t39726+t39728+t39729-t39703*t39706*(t39730+t39726+t39728+t39729-EcFhd_1*t39725*(1.0/3.0)-EcFhd_3*c*t39716*(5.0/9.0)-EcFhd_4*c*t39716*(4.0/9.0)-EcFhd_3*c*t39716*t39727*(4.0/9.0)))+t39700*(-t39715+t39724+t39703*t39706*(t39715-t39724-EcFld_1*1.0/(t39712*t39712*t39712)*(t39709*t39709)*2.0+EcFld_1*t39717*(EcFld_3*c*t39716*(4.0/9.0)-EcFld_2*t39720*t39718*t39719*(1.0/3.6E1)+EcFld_2*c*t39716*t39689*(2.0/9.0))))-c*t39716*t39735*t39744*(4.0/9.0)+c*t39716*t39735*t39753*(4.0/9.0)+c*t39734*t39735*t39686*(2.0/3.0)+c*t39735*t39757*t39686*(2.0/3.0)-t39744*t39718*t39745*t39719*(1.0/9.0)+t39753*t39718*t39745*t39719*(1.0/9.0))+t39700*t39734*2.0-t39755*t39757*2.0-c*t39735*t39744*t39686*(2.0/3.0)+c*t39735*t39753*t39686*(2.0/3.0);
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t39761 = rho_a+rho_b;
                double t39762 = 1.0/pow(t39761,4.0/3.0);
                double t39763 = 1.0/pow(t39761,1.0/3.0);
                double t39764 = c*t39763;
                double t39765 = 1.0/sqrt(t39764);
                double t39766 = sqrt(t39764);
                double t39767 = EcPld_3*c*t39762*(1.0/3.0);
                double t39768 = EcPld_2*c*t39762*t39765*(1.0/6.0);
                double t39769 = t39767+t39768;
                double t39770 = EcPld_2*t39766;
                double t39771 = EcPld_3*c*t39763;
                double t39772 = t39770+t39771+1.0;
                double t39773 = 1.0/(t39772*t39772);
                double t39774 = EcPld_1*t39773*t39769;
                double t39775 = t39764-1.0;
                double t39776 = heaviside(t39775);
                double t39777 = pow(2.0,1.0/3.0);
                double t39778 = t39777*2.0;
                double t39779 = t39778-2.0;
                double t39780 = two_13*2.0;
                double t39781 = t39780-2.0;
                double t39782 = 1.0/t39781;
                double t39783 = EcFld_3*c*t39762*(1.0/3.0);
                double t39784 = EcFld_2*c*t39762*t39765*(1.0/6.0);
                double t39785 = t39783+t39784;
                double t39786 = EcFld_2*t39766;
                double t39787 = EcFld_3*c*t39763;
                double t39788 = t39786+t39787+1.0;
                double t39789 = t39769*t39769;
                double t39790 = 1.0/(t39772*t39772*t39772);
                double t39791 = EcPld_1*t39790*t39789*2.0;
                double t39792 = 1.0/pow(t39761,7.0/3.0);
                double t39793 = 1.0/(t39788*t39788);
                double t39794 = c*c;
                double t39795 = 1.0/pow(t39761,8.0/3.0);
                double t39796 = 1.0/pow(t39764,3.0/2.0);
                double t39797 = EcPld_3*c*t39792*(4.0/9.0);
                double t39798 = EcPld_2*c*t39765*t39792*(2.0/9.0);
                double t39799 = t39797+t39798-EcPld_2*t39794*t39795*t39796*(1.0/3.6E1);
                double t39800 = EcPld_1*t39773*t39799;
                double t39801 = 1.0/(t39761*t39761);
                double t39802 = EcPhd_1*t39801*(1.0/3.0);
                double t39803 = log(t39764);
                double t39804 = EcPhd_3*c*t39792*(5.0/9.0);
                double t39805 = EcPhd_4*c*t39792*(4.0/9.0);
                double t39806 = EcPhd_3*c*t39803*t39792*(4.0/9.0);
                double t39807 = 1.0/t39772;
                double t39808 = EcPld_1*t39807;
                double t39809 = t39774-EcFld_1*t39793*t39785;
                double t39810 = t39774-t39782*t39809*t39779;
                double t39811 = dirac(t39775);
                double t39812 = EcPhd_1*t39803;
                double t39813 = 1.0/t39761;
                double t39814 = EcPhd_1*t39813*(1.0/3.0);
                double t39815 = EcPhd_3*c*t39762*(1.0/3.0);
                double t39816 = EcPhd_4*c*t39762*(1.0/3.0);
                double t39817 = EcPhd_3*c*t39803*t39762*(1.0/3.0);
                double t39818 = 1.0/t39788;
                double t39834 = EcFld_1*t39818;
                double t39819 = -t39834+t39808;
                double t39820 = t39808-t39782*t39819*t39779;
                double t39821 = dirac(t39775,1.0);
                double t39822 = EcFhd_1*t39803;
                double t39823 = EcFhd_4*c*t39763;
                double t39824 = EcPhd_4*c*t39763;
                double t39825 = EcFhd_3*c*t39803*t39763;
                double t39826 = EcPhd_3*c*t39803*t39763;
                double t39827 = -t39764+1.0;
                double t39828 = heaviside(t39827);
                double t39829 = t39814+t39815+t39816+t39817-EcFhd_1*t39813*(1.0/3.0)-EcFhd_3*c*t39762*(1.0/3.0)-EcFhd_4*c*t39762*(1.0/3.0)-EcFhd_3*c*t39803*t39762*(1.0/3.0);
                double t39830 = t39814+t39815+t39816+t39817-t39782*t39829*t39779;
                double t39831 = EcFhd_2-EcPhd_2-t39812+t39822+t39823-t39824+t39825-t39826;
                double t39832 = t39831*t39782*t39779;
                double t39833 = EcPhd_2+t39812+t39832+t39824+t39826;
                v_rho_a_rho_b[Q] += scale * t39810*t39776*2.0-t39830*t39828*2.0+t39761*(t39828*(t39802+t39804+t39805+t39806-t39782*t39779*(t39802+t39804+t39805+t39806-EcFhd_1*t39801*(1.0/3.0)-EcFhd_3*c*t39792*(5.0/9.0)-EcFhd_4*c*t39792*(4.0/9.0)-EcFhd_3*c*t39803*t39792*(4.0/9.0)))+t39776*(-t39800+t39791+t39782*t39779*(t39800-t39791+EcFld_1*(t39785*t39785)*1.0/(t39788*t39788*t39788)*2.0-EcFld_1*t39793*(EcFld_3*c*t39792*(4.0/9.0)-EcFld_2*t39794*t39795*t39796*(1.0/3.6E1)+EcFld_2*c*t39765*t39792*(2.0/9.0))))-c*t39810*t39811*t39762*(2.0/3.0)-c*t39811*t39830*t39762*(2.0/3.0)+c*t39811*t39820*t39792*(4.0/9.0)+t39820*t39821*t39794*t39795*(1.0/9.0)-t39821*t39833*t39794*t39795*(1.0/9.0)-c*t39811*t39792*(EcPhd_2+t39812+t39824+t39826+t39782*t39779*(EcFhd_2-EcPhd_2-t39812+t39822+t39823+t39825-EcPhd_4*c*t39763-EcPhd_3*c*t39803*t39763))*(4.0/9.0))+c*t39811*t39833*t39762*(2.0/3.0)-c*t39811*t39762*(t39808+t39782*t39779*(t39834-t39808))*(2.0/3.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t39836 = rho_a+rho_b;
                double t39837 = 1.0/pow(t39836,4.0/3.0);
                double t39838 = 1.0/pow(t39836,1.0/3.0);
                double t39839 = c*t39838;
                double t39840 = 1.0/sqrt(t39839);
                double t39841 = sqrt(t39839);
                double t39842 = EcPld_3*c*t39837*(1.0/3.0);
                double t39843 = EcPld_2*c*t39840*t39837*(1.0/6.0);
                double t39844 = t39842+t39843;
                double t39845 = EcPld_2*t39841;
                double t39846 = EcPld_3*c*t39838;
                double t39847 = t39845+t39846+1.0;
                double t39848 = 1.0/(t39847*t39847);
                double t39849 = EcPld_1*t39844*t39848;
                double t39850 = t39839-1.0;
                double t39851 = heaviside(t39850);
                double t39852 = pow(2.0,1.0/3.0);
                double t39853 = t39852*2.0;
                double t39854 = t39853-2.0;
                double t39855 = two_13*2.0;
                double t39856 = t39855-2.0;
                double t39857 = 1.0/t39856;
                double t39858 = EcFld_3*c*t39837*(1.0/3.0);
                double t39859 = EcFld_2*c*t39840*t39837*(1.0/6.0);
                double t39860 = t39858+t39859;
                double t39861 = EcFld_2*t39841;
                double t39862 = EcFld_3*c*t39838;
                double t39863 = t39861+t39862+1.0;
                double t39864 = t39844*t39844;
                double t39865 = 1.0/(t39847*t39847*t39847);
                double t39866 = EcPld_1*t39864*t39865*2.0;
                double t39867 = 1.0/pow(t39836,7.0/3.0);
                double t39868 = 1.0/(t39863*t39863);
                double t39869 = c*c;
                double t39870 = 1.0/pow(t39836,8.0/3.0);
                double t39871 = 1.0/pow(t39839,3.0/2.0);
                double t39872 = EcPld_3*c*t39867*(4.0/9.0);
                double t39873 = EcPld_2*c*t39840*t39867*(2.0/9.0);
                double t39874 = t39872+t39873-EcPld_2*t39870*t39871*t39869*(1.0/3.6E1);
                double t39875 = EcPld_1*t39874*t39848;
                double t39876 = 1.0/(t39836*t39836);
                double t39877 = EcPhd_1*t39876*(1.0/3.0);
                double t39878 = log(t39839);
                double t39879 = EcPhd_3*c*t39867*(5.0/9.0);
                double t39880 = EcPhd_4*c*t39867*(4.0/9.0);
                double t39881 = EcPhd_3*c*t39867*t39878*(4.0/9.0);
                double t39882 = 1.0/t39847;
                double t39883 = EcPld_1*t39882;
                double t39884 = t39849-EcFld_1*t39860*t39868;
                double t39885 = t39849-t39854*t39857*t39884;
                double t39886 = dirac(t39850);
                double t39887 = EcPhd_1*t39878;
                double t39888 = 1.0/t39836;
                double t39889 = EcPhd_1*t39888*(1.0/3.0);
                double t39890 = EcPhd_3*c*t39837*(1.0/3.0);
                double t39891 = EcPhd_4*c*t39837*(1.0/3.0);
                double t39892 = EcPhd_3*c*t39837*t39878*(1.0/3.0);
                double t39893 = 1.0/t39863;
                double t39909 = EcFld_1*t39893;
                double t39894 = -t39909+t39883;
                double t39895 = t39883-t39854*t39857*t39894;
                double t39896 = dirac(t39850,1.0);
                double t39897 = EcFhd_1*t39878;
                double t39898 = EcFhd_4*c*t39838;
                double t39899 = EcPhd_4*c*t39838;
                double t39900 = EcFhd_3*c*t39838*t39878;
                double t39901 = EcPhd_3*c*t39838*t39878;
                double t39902 = EcFhd_2-EcPhd_2+t39900-t39901-t39887+t39897+t39898-t39899;
                double t39903 = t39902*t39854*t39857;
                double t39904 = EcPhd_2+t39901+t39903+t39887+t39899;
                double t39905 = -t39839+1.0;
                double t39906 = heaviside(t39905);
                double t39907 = t39890+t39891+t39892+t39889-EcFhd_1*t39888*(1.0/3.0)-EcFhd_3*c*t39837*(1.0/3.0)-EcFhd_4*c*t39837*(1.0/3.0)-EcFhd_3*c*t39837*t39878*(1.0/3.0);
                double t39908 = t39890+t39891+t39892+t39889-t39907*t39854*t39857;
                v_rho_b_rho_b[Q] += scale * -t39836*(-t39906*(t39880+t39881+t39877+t39879-t39854*t39857*(t39880+t39881+t39877+t39879-EcFhd_1*t39876*(1.0/3.0)-EcFhd_3*c*t39867*(5.0/9.0)-EcFhd_4*c*t39867*(4.0/9.0)-EcFhd_3*c*t39867*t39878*(4.0/9.0)))+t39851*(-t39866+t39875+t39854*t39857*(t39866-t39875-EcFld_1*(t39860*t39860)*1.0/(t39863*t39863*t39863)*2.0+EcFld_1*t39868*(EcFld_3*c*t39867*(4.0/9.0)-EcFld_2*t39870*t39871*t39869*(1.0/3.6E1)+EcFld_2*c*t39840*t39867*(2.0/9.0))))+c*t39904*t39867*t39886*(4.0/9.0)+c*t39908*t39837*t39886*(2.0/3.0)+c*t39837*t39885*t39886*(2.0/3.0)-c*t39867*t39886*t39895*(4.0/9.0)+t39904*t39870*t39869*t39896*(1.0/9.0)-t39870*t39895*t39869*t39896*(1.0/9.0))-t39906*t39908*2.0+t39851*t39885*2.0+c*t39904*t39837*t39886*(2.0/3.0)-c*t39837*t39886*(t39883+t39854*t39857*(t39909-t39883))*(2.0/3.0);
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t38725 = rho_a+rho_b;
                double t38726 = 1.0/pow(t38725,1.0/3.0);
                double t38727 = 1.0/t38725;
                double t38728 = rho_a-rho_b;
                double t38729 = t38727*t38728;
                double t38730 = c*t38726;
                double t38731 = log(t38730);
                double t38732 = EcPhd_1*t38731;
                double t38733 = two_13*2.0;
                double t38734 = t38733-2.0;
                double t38735 = 1.0/t38734;
                double t38736 = sqrt(t38730);
                double t38737 = EcPld_2*t38736;
                double t38738 = EcPld_3*c*t38726;
                double t38739 = t38737+t38738+1.0;
                double t38740 = 1.0/t38739;
                double t38741 = EcPld_1*t38740;
                double t38742 = t38729+1.0;
                double t38743 = pow(t38742,4.0/3.0);
                double t38744 = -t38729+1.0;
                double t38745 = pow(t38744,4.0/3.0);
                double t38746 = t38743+t38745-2.0;
                v[Q] += scale * t38725*(heaviside(-c*t38726+1.0)*(EcPhd_2+t38732+t38735*t38746*(EcFhd_2-EcPhd_2-t38732+EcFhd_1*t38731+EcFhd_4*c*t38726-EcPhd_4*c*t38726+EcFhd_3*c*t38731*t38726-EcPhd_3*c*t38731*t38726)+EcPhd_4*c*t38726+EcPhd_3*c*t38731*t38726)+heaviside(t38730-1.0)*(t38741-t38735*t38746*(t38741-EcFld_1/(EcFld_2*t38736+EcFld_3*c*t38726+1.0))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t38748 = rho_a+rho_b;
                double t38749 = 1.0/pow(t38748,1.0/3.0);
                double t38750 = 1.0/t38748;
                double t38751 = rho_a-rho_b;
                double t38752 = t38750*t38751;
                double t38753 = c*t38749;
                double t38754 = log(t38753);
                double t38755 = EcPhd_1*t38754;
                double t38756 = 1.0/pow(t38748,4.0/3.0);
                double t38757 = two_13*2.0;
                double t38758 = t38757-2.0;
                double t38759 = 1.0/t38758;
                double t38760 = t38752+1.0;
                double t38761 = pow(t38760,4.0/3.0);
                double t38762 = -t38752+1.0;
                double t38763 = pow(t38762,4.0/3.0);
                double t38764 = t38761+t38763-2.0;
                double t38765 = EcPhd_1*t38750*(1.0/3.0);
                double t38766 = EcPhd_3*c*t38756*(1.0/3.0);
                double t38767 = EcPhd_4*c*t38756*(1.0/3.0);
                double t38768 = 1.0/(t38748*t38748);
                double t38779 = t38751*t38768;
                double t38769 = t38750-t38779;
                double t38770 = EcFhd_1*t38754;
                double t38771 = EcFhd_4*c*t38749;
                double t38772 = EcPhd_4*c*t38749;
                double t38773 = EcFhd_3*c*t38754*t38749;
                double t38774 = EcPhd_3*c*t38754*t38749;
                double t38775 = EcPhd_3*c*t38754*t38756*(1.0/3.0);
                double t38776 = 1.0/sqrt(t38753);
                double t38777 = sqrt(t38753);
                double t38778 = pow(t38760,1.0/3.0);
                double t38780 = t38769*t38778*(4.0/3.0);
                double t38781 = pow(t38762,1.0/3.0);
                double t38782 = t38780-t38781*t38769*(4.0/3.0);
                double t38783 = EcFld_2*t38777;
                double t38784 = EcFld_3*c*t38749;
                double t38785 = t38783+t38784+1.0;
                double t38786 = EcPld_2*t38777;
                double t38787 = EcPld_3*c*t38749;
                double t38788 = t38786+t38787+1.0;
                double t38789 = EcPld_3*c*t38756*(1.0/3.0);
                double t38790 = EcPld_2*c*t38756*t38776*(1.0/6.0);
                double t38791 = t38790+t38789;
                double t38792 = 1.0/(t38788*t38788);
                double t38793 = t38753-1.0;
                double t38794 = EcFhd_2-EcPhd_2+t38770+t38771-t38772-t38755+t38773-t38774;
                double t38795 = dirac(t38793);
                double t38796 = 1.0/t38788;
                double t38797 = EcPld_1*t38796;
                double t38798 = 1.0/t38785;
                double t38801 = EcFld_1*t38798;
                double t38799 = -t38801+t38797;
                double t38800 = heaviside(t38793);
                v_rho_a[Q] += scale * -t38748*(heaviside(-t38753+1.0)*(t38765+t38766+t38775+t38767-t38764*t38759*(t38765+t38766+t38775+t38767-EcFhd_1*t38750*(1.0/3.0)-EcFhd_3*c*t38756*(1.0/3.0)-EcFhd_4*c*t38756*(1.0/3.0)-EcFhd_3*c*t38754*t38756*(1.0/3.0))-t38782*t38794*t38759)-t38800*(t38764*t38759*(EcFld_1*1.0/(t38785*t38785)*(EcFld_3*c*t38756*(1.0/3.0)+EcFld_2*c*t38756*t38776*(1.0/6.0))-EcPld_1*t38791*t38792)+EcPld_1*t38791*t38792-t38782*t38759*t38799)+c*t38756*t38795*(t38797-t38764*t38759*t38799)*(1.0/3.0)-c*t38756*t38795*(EcPhd_2+t38772+t38755+t38774+t38764*t38794*t38759)*(1.0/3.0))+heaviside(-c*t38749+1.0)*(EcPhd_2+t38772+t38755+t38774+t38764*t38759*(EcFhd_2-EcPhd_2+t38770+t38771-t38755+t38773-EcPhd_4*c*t38749-EcPhd_3*c*t38754*t38749))+t38800*(t38797+t38764*t38759*(t38801-t38797));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t38803 = rho_a+rho_b;
                double t38804 = 1.0/pow(t38803,1.0/3.0);
                double t38805 = 1.0/t38803;
                double t38806 = rho_a-rho_b;
                double t38807 = t38805*t38806;
                double t38808 = c*t38804;
                double t38809 = log(t38808);
                double t38810 = EcPhd_1*t38809;
                double t38811 = 1.0/pow(t38803,4.0/3.0);
                double t38812 = two_13*2.0;
                double t38813 = t38812-2.0;
                double t38814 = 1.0/t38813;
                double t38815 = t38807+1.0;
                double t38816 = pow(t38815,4.0/3.0);
                double t38817 = -t38807+1.0;
                double t38818 = pow(t38817,4.0/3.0);
                double t38819 = t38816+t38818-2.0;
                double t38820 = EcPhd_1*t38805*(1.0/3.0);
                double t38821 = EcPhd_3*c*t38811*(1.0/3.0);
                double t38822 = EcPhd_4*c*t38811*(1.0/3.0);
                double t38823 = 1.0/(t38803*t38803);
                double t38824 = t38823*t38806;
                double t38825 = t38805+t38824;
                double t38826 = EcFhd_1*t38809;
                double t38827 = EcFhd_4*c*t38804;
                double t38828 = EcPhd_4*c*t38804;
                double t38829 = EcFhd_3*c*t38804*t38809;
                double t38830 = EcPhd_3*c*t38804*t38809;
                double t38831 = EcPhd_3*c*t38811*t38809*(1.0/3.0);
                double t38832 = 1.0/sqrt(t38808);
                double t38833 = sqrt(t38808);
                double t38834 = pow(t38815,1.0/3.0);
                double t38835 = t38825*t38834*(4.0/3.0);
                double t38836 = pow(t38817,1.0/3.0);
                double t38837 = t38835-t38825*t38836*(4.0/3.0);
                double t38838 = EcFld_2*t38833;
                double t38839 = EcFld_3*c*t38804;
                double t38840 = t38838+t38839+1.0;
                double t38841 = EcPld_2*t38833;
                double t38842 = EcPld_3*c*t38804;
                double t38843 = t38841+t38842+1.0;
                double t38844 = EcPld_3*c*t38811*(1.0/3.0);
                double t38845 = EcPld_2*c*t38811*t38832*(1.0/6.0);
                double t38846 = t38844+t38845;
                double t38847 = 1.0/(t38843*t38843);
                double t38848 = t38808-1.0;
                double t38849 = EcFhd_2-EcPhd_2-t38810-t38830+t38826+t38827-t38828+t38829;
                double t38850 = dirac(t38848);
                double t38851 = 1.0/t38843;
                double t38852 = EcPld_1*t38851;
                double t38853 = 1.0/t38840;
                double t38856 = EcFld_1*t38853;
                double t38854 = t38852-t38856;
                double t38855 = heaviside(t38848);
                double t38857 = t38852-t38814*t38854*t38819;
                v_rho_b[Q] += scale * -t38803*(heaviside(-t38808+1.0)*(t38820+t38821+t38822+t38831-t38814*t38819*(t38820+t38821+t38822+t38831-EcFhd_1*t38805*(1.0/3.0)-EcFhd_3*c*t38811*(1.0/3.0)-EcFhd_4*c*t38811*(1.0/3.0)-EcFhd_3*c*t38811*t38809*(1.0/3.0))+t38814*t38837*t38849)-t38855*(t38814*t38819*(EcFld_1*1.0/(t38840*t38840)*(EcFld_3*c*t38811*(1.0/3.0)+EcFld_2*c*t38811*t38832*(1.0/6.0))-EcPld_1*t38846*t38847)+EcPld_1*t38846*t38847+t38814*t38854*t38837)+c*t38811*t38850*t38857*(1.0/3.0)-c*t38811*t38850*(EcPhd_2+t38810+t38830+t38828+t38814*t38819*t38849)*(1.0/3.0))+t38855*t38857+heaviside(-c*t38804+1.0)*(EcPhd_2+t38810+t38830+t38828+t38814*t38819*(EcFhd_2-EcPhd_2-t38810+t38826+t38827+t38829-EcPhd_4*c*t38804-EcPhd_3*c*t38804*t38809));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t38864 = rho_a+rho_b;
                double t38865 = 1.0/pow(t38864,4.0/3.0);
                double t38866 = 1.0/pow(t38864,1.0/3.0);
                double t38867 = c*t38866;
                double t38874 = 1.0/sqrt(t38867);
                double t38876 = EcPld_3*c*t38865*(1.0/3.0);
                double t38877 = EcPld_2*c*t38865*t38874*(1.0/6.0);
                double t38868 = t38876+t38877;
                double t38869 = 1.0/t38864;
                double t38870 = rho_a-rho_b;
                double t38871 = t38870*t38869;
                double t38872 = 1.0/(t38864*t38864);
                double t38909 = t38870*t38872;
                double t38873 = -t38909+t38869;
                double t38875 = sqrt(t38867);
                double t38878 = EcPld_2*t38875;
                double t38879 = EcPld_3*c*t38866;
                double t38880 = t38878+t38879+1.0;
                double t38881 = two_13*2.0;
                double t38882 = t38881-2.0;
                double t38883 = 1.0/t38882;
                double t38884 = t38871+1.0;
                double t38885 = -t38871+1.0;
                double t38886 = EcFld_3*c*t38865*(1.0/3.0);
                double t38887 = EcFld_2*c*t38865*t38874*(1.0/6.0);
                double t38888 = t38886+t38887;
                double t38889 = EcFld_2*t38875;
                double t38890 = EcFld_3*c*t38866;
                double t38891 = t38890+t38889+1.0;
                double t38892 = t38868*t38868;
                double t38893 = 1.0/(t38880*t38880*t38880);
                double t38894 = EcPld_1*t38892*t38893*2.0;
                double t38895 = 1.0/pow(t38864,7.0/3.0);
                double t38896 = 1.0/(t38891*t38891);
                double t38897 = c*c;
                double t38898 = 1.0/pow(t38864,8.0/3.0);
                double t38899 = 1.0/pow(t38867,3.0/2.0);
                double t38900 = 1.0/(t38880*t38880);
                double t38901 = EcPld_3*c*t38895*(4.0/9.0);
                double t38902 = EcPld_2*c*t38874*t38895*(2.0/9.0);
                double t38903 = t38901+t38902-EcPld_2*t38897*t38898*t38899*(1.0/3.6E1);
                double t38904 = pow(t38884,1.0/3.0);
                double t38905 = pow(t38885,1.0/3.0);
                double t38906 = t38872*2.0;
                double t38907 = 1.0/(t38864*t38864*t38864);
                double t38911 = t38870*t38907*2.0;
                double t38908 = -t38911+t38906;
                double t38910 = t38873*t38873;
                double t38912 = t38905*t38908*(4.0/3.0);
                double t38913 = 1.0/pow(t38884,2.0/3.0);
                double t38914 = t38910*t38913*(4.0/9.0);
                double t38915 = 1.0/pow(t38885,2.0/3.0);
                double t38916 = t38910*t38915*(4.0/9.0);
                double t38917 = t38912+t38914+t38916-t38904*t38908*(4.0/3.0);
                double t38918 = log(t38867);
                double t38919 = pow(t38884,4.0/3.0);
                double t38920 = pow(t38885,4.0/3.0);
                double t38921 = t38920+t38919-2.0;
                double t38922 = EcPhd_1*t38872*(1.0/3.0);
                double t38923 = EcPhd_3*c*t38895*(5.0/9.0);
                double t38924 = EcPhd_4*c*t38895*(4.0/9.0);
                double t38925 = t38904*t38873*(4.0/3.0);
                double t38945 = t38905*t38873*(4.0/3.0);
                double t38926 = t38925-t38945;
                double t38927 = EcPhd_3*c*t38918*t38895*(4.0/9.0);
                double t38928 = t38867-1.0;
                double t38929 = 1.0/t38880;
                double t38930 = EcPld_1*t38929;
                double t38931 = 1.0/t38891;
                double t38950 = EcFld_1*t38931;
                double t38932 = t38930-t38950;
                double t38933 = EcFhd_1*t38918;
                double t38934 = EcPhd_1*t38918;
                double t38935 = EcFhd_4*c*t38866;
                double t38936 = EcFhd_3*c*t38918*t38866;
                double t38937 = dirac(t38928);
                double t38938 = EcFhd_1*t38869*(1.0/3.0);
                double t38939 = EcPhd_1*t38869*(1.0/3.0);
                double t38940 = EcFhd_3*c*t38865*(1.0/3.0);
                double t38941 = EcFhd_4*c*t38865*(1.0/3.0);
                double t38942 = EcPhd_3*c*t38865*(1.0/3.0);
                double t38943 = EcPhd_4*c*t38865*(1.0/3.0);
                double t38944 = EcFhd_3*c*t38918*t38865*(1.0/3.0);
                double t38946 = EcPhd_4*c*t38866;
                double t38947 = EcPhd_3*c*t38918*t38866;
                double t38948 = EcFld_1*t38896*t38888;
                double t38951 = EcPld_1*t38900*t38868;
                double t38949 = -t38951+t38948;
                double t38967 = t38921*t38932*t38883;
                double t38952 = t38930-t38967;
                double t38953 = dirac(t38928,1.0);
                double t38954 = EcFhd_2-EcPhd_2+t38933-t38934+t38935+t38936-t38946-t38947;
                double t38955 = EcPld_1*t38900*t38903;
                double t38956 = t38921*t38954*t38883;
                double t38957 = EcPhd_2+t38934+t38946+t38947+t38956;
                double t38959 = EcPhd_3*c*t38918*t38865*(1.0/3.0);
                double t38958 = t38940+t38941-t38942-t38943+t38944+t38938-t38939-t38959;
                double t38960 = -t38867+1.0;
                double t38961 = heaviside(t38960);
                double t38962 = t38921*t38883*t38958;
                double t38963 = t38942+t38943+t38962+t38939+t38959-t38926*t38954*t38883;
                double t38964 = heaviside(t38928);
                double t38965 = t38921*t38883*t38949;
                double t38966 = t38951+t38965-t38932*t38926*t38883;
                v_rho_a_rho_a[Q] += scale * t38864*(t38964*(-t38955+t38894+t38921*t38883*(t38955-t38894+EcFld_1*1.0/(t38891*t38891*t38891)*(t38888*t38888)*2.0-EcFld_1*t38896*(EcFld_3*c*t38895*(4.0/9.0)-EcFld_2*t38897*t38898*t38899*(1.0/3.6E1)+EcFld_2*c*t38874*t38895*(2.0/9.0)))-t38932*t38917*t38883+t38926*t38883*t38949*2.0)+t38961*(t38922+t38923+t38924+t38927-t38921*t38883*(t38922+t38923+t38924+t38927-EcFhd_1*t38872*(1.0/3.0)-EcFhd_3*c*t38895*(5.0/9.0)-EcFhd_4*c*t38895*(4.0/9.0)-EcFhd_3*c*t38918*t38895*(4.0/9.0))+t38917*t38954*t38883-t38926*t38883*t38958*2.0)-c*t38963*t38865*t38937*(2.0/3.0)+c*t38952*t38937*t38895*(4.0/9.0)-c*t38865*t38937*t38966*(2.0/3.0)-c*t38937*t38957*t38895*(4.0/9.0)+t38952*t38953*t38897*t38898*(1.0/9.0)-t38953*t38957*t38897*t38898*(1.0/9.0))-t38961*t38963*2.0+t38964*t38966*2.0-c*t38952*t38865*t38937*(2.0/3.0)+c*t38865*t38937*t38957*(2.0/3.0);
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t38969 = rho_a+rho_b;
                double t38970 = rho_a-rho_b;
                double t38971 = 1.0/t38969;
                double t38972 = t38970*t38971;
                double t38973 = 1.0/(t38969*t38969*t38969);
                double t38974 = t38972+1.0;
                double t38975 = 1.0/(t38969*t38969);
                double t38976 = t38970*t38975;
                double t38977 = -t38972+1.0;
                double t38978 = t38971+t38976;
                double t38979 = t38971-t38976;
                double t38980 = 1.0/pow(t38969,1.0/3.0);
                double t38981 = c*t38980;
                double t38982 = log(t38981);
                double t38983 = 1.0/pow(t38969,7.0/3.0);
                double t38984 = two_13*2.0;
                double t38985 = t38984-2.0;
                double t38986 = 1.0/t38985;
                double t38987 = EcPhd_1*t38975*(1.0/3.0);
                double t38988 = EcPhd_3*c*t38983*(5.0/9.0);
                double t38989 = EcPhd_4*c*t38983*(4.0/9.0);
                double t38990 = pow(t38974,1.0/3.0);
                double t38991 = pow(t38977,1.0/3.0);
                double t38992 = 1.0/pow(t38969,4.0/3.0);
                double t38993 = EcFhd_1*t38971*(1.0/3.0);
                double t38994 = EcFhd_3*c*t38992*(1.0/3.0);
                double t38995 = EcFhd_4*c*t38992*(1.0/3.0);
                double t38996 = EcFhd_3*c*t38982*t38992*(1.0/3.0);
                double t39049 = EcPhd_1*t38971*(1.0/3.0);
                double t39050 = EcPhd_3*c*t38992*(1.0/3.0);
                double t39051 = EcPhd_4*c*t38992*(1.0/3.0);
                double t39052 = EcPhd_3*c*t38982*t38992*(1.0/3.0);
                double t38997 = t38993+t38994+t38995+t38996-t39050-t39051-t39052-t39049;
                double t38998 = EcPhd_3*c*t38982*t38983*(4.0/9.0);
                double t39002 = 1.0/sqrt(t38981);
                double t39004 = EcPld_3*c*t38992*(1.0/3.0);
                double t39005 = EcPld_2*c*t38992*t39002*(1.0/6.0);
                double t38999 = t39004+t39005;
                double t39000 = t38990*t38978*(4.0/3.0);
                double t39053 = t38991*t38978*(4.0/3.0);
                double t39001 = t39000-t39053;
                double t39003 = sqrt(t38981);
                double t39006 = EcPld_2*t39003;
                double t39007 = EcPld_3*c*t38980;
                double t39008 = t39006+t39007+1.0;
                double t39009 = t38990*t38979*(4.0/3.0);
                double t39057 = t38991*t38979*(4.0/3.0);
                double t39010 = t39009-t39057;
                double t39011 = EcFld_3*c*t38992*(1.0/3.0);
                double t39012 = EcFld_2*c*t38992*t39002*(1.0/6.0);
                double t39013 = t39011+t39012;
                double t39014 = EcFld_2*t39003;
                double t39015 = EcFld_3*c*t38980;
                double t39016 = t39014+t39015+1.0;
                double t39017 = 1.0/(t39016*t39016);
                double t39018 = EcFld_1*t39013*t39017;
                double t39019 = 1.0/(t39008*t39008);
                double t39059 = EcPld_1*t38999*t39019;
                double t39020 = t39018-t39059;
                double t39021 = pow(t38974,4.0/3.0);
                double t39022 = pow(t38977,4.0/3.0);
                double t39023 = t39021+t39022-2.0;
                double t39024 = t38999*t38999;
                double t39025 = 1.0/(t39008*t39008*t39008);
                double t39026 = EcPld_1*t39024*t39025*2.0;
                double t39027 = c*c;
                double t39028 = 1.0/pow(t38969,8.0/3.0);
                double t39029 = 1.0/pow(t38981,3.0/2.0);
                double t39030 = EcPld_3*c*t38983*(4.0/9.0);
                double t39031 = EcPld_2*c*t38983*t39002*(2.0/9.0);
                double t39032 = t39030+t39031-EcPld_2*t39027*t39028*t39029*(1.0/3.6E1);
                double t39033 = t38970*t38973*t38991*(8.0/3.0);
                double t39034 = 1.0/pow(t38974,2.0/3.0);
                double t39035 = t38978*t38979*t39034*(4.0/9.0);
                double t39036 = 1.0/pow(t38977,2.0/3.0);
                double t39037 = t38978*t38979*t39036*(4.0/9.0);
                double t39038 = t39033+t39035+t39037-t38970*t38990*t38973*(8.0/3.0);
                double t39039 = t38981-1.0;
                double t39040 = 1.0/t39008;
                double t39041 = EcPld_1*t39040;
                double t39042 = 1.0/t39016;
                double t39060 = EcFld_1*t39042;
                double t39043 = t39041-t39060;
                double t39044 = EcFhd_1*t38982;
                double t39045 = EcPhd_1*t38982;
                double t39046 = EcFhd_4*c*t38980;
                double t39047 = EcFhd_3*c*t38980*t38982;
                double t39048 = dirac(t39039);
                double t39054 = EcPhd_4*c*t38980;
                double t39055 = EcPhd_3*c*t38980*t38982;
                double t39056 = t38986*t38997*t39023;
                double t39058 = EcFhd_2-EcPhd_2+t39044-t39045-t39054+t39046-t39055+t39047;
                double t39061 = t38986*t39020*t39023;
                double t39073 = t38986*t39023*t39043;
                double t39062 = t39041-t39073;
                double t39063 = dirac(t39039,1.0);
                double t39064 = t38986*t39023*t39058;
                double t39065 = EcPhd_2+t39045+t39054+t39055+t39064;
                double t39066 = -t38981+1.0;
                double t39067 = heaviside(t39066);
                double t39068 = t38986*t39001*t39058;
                double t39069 = t39050+t39051+t39052+t39056+t39049+t39068;
                double t39070 = t39050+t39051+t39052+t39056+t39049-t38986*t39010*t39058;
                double t39071 = heaviside(t39039);
                double t39072 = t39061+t39059-t38986*t39010*t39043;
                v_rho_a_rho_b[Q] += scale * t39071*(t39061+t39059+t38986*t39001*(t39041-t39060))+t39071*t39072-t39070*t39067-t39067*t39069-t38969*(-t39071*(t39026-EcPld_1*t39032*t39019-t38986*t39001*t39020+t38986*t39010*t39020+t38986*t39043*t39038-t38986*t39023*(t39026-EcFld_1*(t39013*t39013)*1.0/(t39016*t39016*t39016)*2.0-EcPld_1*t39032*t39019+EcFld_1*t39017*(EcFld_3*c*t38983*(4.0/9.0)-EcFld_2*t39027*t39028*t39029*(1.0/3.6E1)+EcFld_2*c*t38983*t39002*(2.0/9.0))))-t39067*(t38987+t38988+t38989+t38998-t38986*t39023*(t38987+t38988+t38989+t38998-EcFhd_1*t38975*(1.0/3.0)-EcFhd_3*c*t38983*(5.0/9.0)-EcFhd_4*c*t38983*(4.0/9.0)-EcFhd_3*c*t38982*t38983*(4.0/9.0))+t38986*t38997*t39001-t38986*t38997*t39010-t38986*t39038*t39058)+c*t38992*t39070*t39048*(1.0/3.0)-c*t38983*t39062*t39048*(4.0/9.0)+c*t38992*t39072*t39048*(1.0/3.0)+c*t38983*t39065*t39048*(4.0/9.0)+c*t38992*t39048*t39069*(1.0/3.0)-t39062*t39027*t39063*t39028*(1.0/9.0)+t39027*t39063*t39028*t39065*(1.0/9.0)+c*t38992*t39048*(t39061+t39059+t38986*t39001*t39043)*(1.0/3.0))-c*t38992*t39062*t39048*(2.0/3.0)+c*t38992*t39065*t39048*(2.0/3.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t39075 = rho_a+rho_b;
                double t39076 = 1.0/pow(t39075,4.0/3.0);
                double t39077 = 1.0/pow(t39075,1.0/3.0);
                double t39078 = c*t39077;
                double t39086 = 1.0/sqrt(t39078);
                double t39088 = EcPld_3*c*t39076*(1.0/3.0);
                double t39089 = EcPld_2*c*t39076*t39086*(1.0/6.0);
                double t39079 = t39088+t39089;
                double t39080 = 1.0/t39075;
                double t39081 = rho_a-rho_b;
                double t39082 = t39080*t39081;
                double t39083 = 1.0/(t39075*t39075);
                double t39084 = t39081*t39083;
                double t39085 = t39080+t39084;
                double t39087 = sqrt(t39078);
                double t39090 = EcPld_2*t39087;
                double t39091 = EcPld_3*c*t39077;
                double t39092 = t39090+t39091+1.0;
                double t39093 = two_13*2.0;
                double t39094 = t39093-2.0;
                double t39095 = 1.0/t39094;
                double t39096 = t39082+1.0;
                double t39097 = -t39082+1.0;
                double t39098 = EcFld_3*c*t39076*(1.0/3.0);
                double t39099 = EcFld_2*c*t39076*t39086*(1.0/6.0);
                double t39100 = t39098+t39099;
                double t39101 = EcFld_2*t39087;
                double t39102 = EcFld_3*c*t39077;
                double t39103 = t39101+t39102+1.0;
                double t39104 = t39079*t39079;
                double t39105 = 1.0/(t39092*t39092*t39092);
                double t39106 = EcPld_1*t39104*t39105*2.0;
                double t39107 = 1.0/pow(t39075,7.0/3.0);
                double t39108 = 1.0/(t39103*t39103);
                double t39109 = c*c;
                double t39110 = 1.0/pow(t39075,8.0/3.0);
                double t39111 = 1.0/pow(t39078,3.0/2.0);
                double t39112 = 1.0/(t39092*t39092);
                double t39113 = EcPld_3*c*t39107*(4.0/9.0);
                double t39114 = EcPld_2*c*t39107*t39086*(2.0/9.0);
                double t39115 = t39113+t39114-EcPld_2*t39110*t39111*t39109*(1.0/3.6E1);
                double t39116 = pow(t39096,1.0/3.0);
                double t39117 = pow(t39097,1.0/3.0);
                double t39118 = t39083*2.0;
                double t39119 = 1.0/(t39075*t39075*t39075);
                double t39120 = t39081*t39119*2.0;
                double t39121 = t39120+t39118;
                double t39122 = t39085*t39085;
                double t39123 = t39121*t39116*(4.0/3.0);
                double t39124 = 1.0/pow(t39096,2.0/3.0);
                double t39125 = t39122*t39124*(4.0/9.0);
                double t39126 = 1.0/pow(t39097,2.0/3.0);
                double t39127 = t39122*t39126*(4.0/9.0);
                double t39128 = t39123+t39125+t39127-t39121*t39117*(4.0/3.0);
                double t39129 = log(t39078);
                double t39130 = pow(t39096,4.0/3.0);
                double t39131 = pow(t39097,4.0/3.0);
                double t39132 = t39130+t39131-2.0;
                double t39133 = EcPhd_1*t39083*(1.0/3.0);
                double t39134 = EcPhd_3*c*t39107*(5.0/9.0);
                double t39135 = EcPhd_4*c*t39107*(4.0/9.0);
                double t39136 = t39116*t39085*(4.0/3.0);
                double t39156 = t39117*t39085*(4.0/3.0);
                double t39137 = t39136-t39156;
                double t39138 = EcPhd_3*c*t39107*t39129*(4.0/9.0);
                double t39139 = t39078-1.0;
                double t39140 = 1.0/t39092;
                double t39141 = EcPld_1*t39140;
                double t39142 = 1.0/t39103;
                double t39162 = EcFld_1*t39142;
                double t39143 = t39141-t39162;
                double t39144 = EcFhd_1*t39129;
                double t39145 = EcPhd_1*t39129;
                double t39146 = EcFhd_4*c*t39077;
                double t39147 = EcFhd_3*c*t39129*t39077;
                double t39148 = dirac(t39139);
                double t39149 = EcFhd_1*t39080*(1.0/3.0);
                double t39150 = EcPhd_1*t39080*(1.0/3.0);
                double t39151 = EcFhd_3*c*t39076*(1.0/3.0);
                double t39152 = EcFhd_4*c*t39076*(1.0/3.0);
                double t39153 = EcPhd_3*c*t39076*(1.0/3.0);
                double t39154 = EcPhd_4*c*t39076*(1.0/3.0);
                double t39155 = EcFhd_3*c*t39129*t39076*(1.0/3.0);
                double t39157 = EcPhd_4*c*t39077;
                double t39158 = EcPhd_3*c*t39129*t39077;
                double t39159 = EcPhd_3*c*t39129*t39076*(1.0/3.0);
                double t39160 = EcFld_1*t39100*t39108;
                double t39163 = EcPld_1*t39112*t39079;
                double t39161 = t39160-t39163;
                double t39177 = t39132*t39143*t39095;
                double t39164 = t39141-t39177;
                double t39165 = dirac(t39139,1.0);
                double t39166 = EcFhd_2-EcPhd_2+t39144-t39145+t39146+t39147-t39157-t39158;
                double t39167 = EcPld_1*t39112*t39115;
                double t39168 = t39150-t39151-t39152+t39153+t39154-t39155-t39149+t39159;
                double t39169 = t39132*t39166*t39095;
                double t39170 = EcPhd_2+t39145+t39157+t39158+t39169;
                double t39171 = -t39078+1.0;
                double t39172 = heaviside(t39171);
                double t39173 = t39137*t39166*t39095;
                double t39174 = t39150+t39153+t39154+t39173+t39159-t39132*t39095*t39168;
                double t39175 = heaviside(t39139);
                double t39176 = t39132*t39161*t39095;
                v_rho_b_rho_b[Q] += scale * t39175*(t39163+t39176+t39137*t39095*(t39141-t39162))*2.0-t39172*t39174*2.0-t39075*(t39175*(-t39106+t39167+t39132*t39095*(t39106-t39167-EcFld_1*(t39100*t39100)*1.0/(t39103*t39103*t39103)*2.0+EcFld_1*t39108*(EcFld_3*c*t39107*(4.0/9.0)-EcFld_2*t39110*t39111*t39109*(1.0/3.6E1)+EcFld_2*c*t39107*t39086*(2.0/9.0)))+t39143*t39128*t39095+t39161*t39137*t39095*2.0)-t39172*(t39133+t39134+t39135+t39138-t39132*t39095*(t39133+t39134+t39135+t39138-EcFhd_1*t39083*(1.0/3.0)-EcFhd_3*c*t39107*(5.0/9.0)-EcFhd_4*c*t39107*(4.0/9.0)-EcFhd_3*c*t39107*t39129*(4.0/9.0))+t39128*t39166*t39095-t39137*t39095*t39168*2.0)+c*t39107*t39170*t39148*(4.0/9.0)-c*t39107*t39164*t39148*(4.0/9.0)+c*t39174*t39076*t39148*(2.0/3.0)+t39110*t39170*t39109*t39165*(1.0/9.0)-t39110*t39109*t39164*t39165*(1.0/9.0)+c*t39076*t39148*(t39163+t39176+t39143*t39137*t39095)*(2.0/3.0))+c*t39170*t39076*t39148*(2.0/3.0)-c*t39164*t39076*t39148*(2.0/3.0);
            }
            
        }
    }
}

}
