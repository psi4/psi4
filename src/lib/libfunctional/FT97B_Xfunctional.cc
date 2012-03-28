#include <libmints/vector.h>
#include "FT97B_Xfunctional.h"
#include "utility.h"
#include <cmath>

using namespace psi;

namespace psi {

FT97B_XFunctional::FT97B_XFunctional()
{
    name_ = "FT97B_X";
    description_ = "    Filitov and Theil 1997 Exchange\n";
    citation_ = "    M. Filatov and W. Theil, Mol. Phys., 91(5), 847-859, 1997.\n";
    alpha_ = 1.0;
    omega_ = 0.0;
    lrc_ = false;
    gga_ = true;
    meta_ = false;
    parameters_["c"] =  -9.3052573634909974E-01;
    parameters_["d0"] =   2.9136440000000000E-03;
    parameters_["d1"] =   9.4741689999999995E-04;
    parameters_["d2"] =   2.5011489999999999E+03;
}
FT97B_XFunctional::~FT97B_XFunctional()
{
}
void FT97B_XFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    double c = parameters_["c"];
    double d0 = parameters_["d0"];
    double d1 = parameters_["d1"];
    double d2 = parameters_["d2"];

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
                double t12250 = 1.0/pow(rho_b,8.0/3.0);
                double t12251 = log(gamma_bb*t12250+sqrt((gamma_bb*gamma_bb)*(t12250*t12250)+1.0));
                double t12253 = d2*d2;
                double t12254 = gamma_bb+t12253;
                double t12255 = 1.0/t12254;
                double t12256 = d1*gamma_bb*t12255;
                double t12252 = d0+t12256;
                v[Q] += scale * -c*pow(rho_b,4.0/3.0)*((gamma_bb*t12250*t12252*1.0/sqrt(gamma_bb*t12250*(t12251*t12251)*(t12252*t12252)*9.0+1.0))/c-1.0);
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t12259 = 1.0/pow(rho_b,8.0/3.0);
                double t12267 = gamma_bb*t12259;
                double t12260 = log(t12267+sqrt(t12267*t12267+1.0));
                double t12262 = d2*d2;
                double t12263 = gamma_bb+t12262;
                double t12264 = 1.0/t12263;
                double t12265 = d1*gamma_bb*t12264;
                double t12261 = d0+t12265;
                double t12266 = 1.0/c;
                double t12268 = t12260*t12260;
                double t12269 = t12261*t12261;
                double t12270 = gamma_bb*t12259*t12268*t12269*9.0;
                double t12271 = t12270+1.0;
                double t12272 = 1.0/sqrt(t12271);
                double t12273 = 1.0/pow(rho_b,1.1E1/3.0);
                double t12274 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * c*pow(rho_b,1.0/3.0)*(gamma_bb*t12261*t12272*t12266*t12259-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t12261*t12272*t12273*t12266*(8.0/3.0)-gamma_bb*t12261*1.0/pow(t12271,3.0/2.0)*t12266*t12259*(gamma_bb*t12273*t12268*t12269*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t12260*t12274*t12269*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t12274+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t12278 = 1.0/pow(rho_b,8.0/3.0);
                double t12286 = gamma_bb*t12278;
                double t12279 = log(t12286+sqrt(t12286*t12286+1.0));
                double t12281 = d2*d2;
                double t12282 = gamma_bb+t12281;
                double t12283 = 1.0/t12282;
                double t12284 = d1*gamma_bb*t12283;
                double t12280 = d0+t12284;
                double t12285 = 1.0/c;
                double t12287 = t12279*t12279;
                double t12288 = t12280*t12280;
                double t12289 = gamma_bb*t12278*t12287*t12288*9.0;
                double t12290 = t12289+1.0;
                double t12291 = 1.0/sqrt(t12290);
                double t12292 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12293 = d1*t12283;
                double t12294 = 1.0/(t12282*t12282);
                double t12295 = t12293-d1*gamma_bb*t12294;
                v_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t12280*t12291*t12285*t12278+gamma_bb*t12291*t12285*t12295*t12278-gamma_bb*t12280*1.0/pow(t12290,3.0/2.0)*t12285*t12278*(t12278*t12287*t12288*9.0+gamma_bb*t12292*t12279*t12288*1.0/sqrt((gamma_bb*gamma_bb)*t12292+1.0)*1.8E1+gamma_bb*t12280*t12295*t12278*t12287*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t12301 = 1.0/pow(rho_b,8.0/3.0);
                double t12309 = gamma_bb*t12301;
                double t12302 = log(t12309+sqrt(t12309*t12309+1.0));
                double t12304 = d2*d2;
                double t12305 = gamma_bb+t12304;
                double t12306 = 1.0/t12305;
                double t12307 = d1*gamma_bb*t12306;
                double t12303 = d0+t12307;
                double t12308 = 1.0/c;
                double t12310 = t12302*t12302;
                double t12311 = t12303*t12303;
                double t12312 = gamma_bb*t12301*t12310*t12311*9.0;
                double t12313 = t12312+1.0;
                double t12314 = 1.0/sqrt(t12313);
                double t12315 = 1.0/pow(rho_b,1.1E1/3.0);
                double t12316 = gamma_bb*gamma_bb;
                double t12317 = 1.0/pow(t12313,3.0/2.0);
                double t12318 = 1.0/pow(rho_b,1.4E1/3.0);
                double t12319 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12320 = t12316*t12319;
                double t12321 = t12320+1.0;
                double t12322 = 1.0/sqrt(t12321);
                double t12323 = gamma_bb*t12310*t12311*t12315*2.4E1;
                double t12324 = 1.0/pow(rho_b,1.9E1/3.0);
                double t12325 = t12302*t12311*t12322*t12324*t12316*4.8E1;
                double t12326 = t12323+t12325;
                v_rho_b_rho_b[Q] += scale * c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t12301*t12303*t12314*t12308-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t12303*t12314*t12308*t12318*(8.8E1/9.0)-gamma_bb*t12303*t12315*t12308*t12317*t12326*(8.0/3.0)-gamma_bb*t12301*t12303*t12308*t12317*(gamma_bb*t12310*t12311*t12318*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t12302*t12311*t12322*t12316*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t12311*t12316*1.28E2)/t12321-1.0/pow(rho_b,3.8E1/3.0)*t12302*t12311*1.0/pow(t12321,3.0/2.0)*(t12316*t12316)*1.28E2)*(1.0/2.0)+gamma_bb*t12301*t12303*1.0/pow(t12313,5.0/2.0)*t12308*(t12326*t12326)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t12303*t12314*t12315*t12308*(8.0/3.0)-gamma_bb*t12301*t12303*t12308*t12317*t12326*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t12372 = d2*d2;
                double t12373 = gamma_bb+t12372;
                double t12374 = 1.0/pow(rho_b,8.0/3.0);
                double t12379 = gamma_bb*t12374;
                double t12375 = log(t12379+sqrt(t12379*t12379+1.0));
                double t12376 = 1.0/t12373;
                double t12381 = d1*gamma_bb*t12376;
                double t12377 = d0+t12381;
                double t12378 = 1.0/c;
                double t12380 = t12375*t12375;
                double t12382 = t12377*t12377;
                double t12383 = gamma_bb*t12380*t12382*t12374*9.0;
                double t12384 = t12383+1.0;
                double t12385 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12386 = d1*t12376;
                double t12387 = 1.0/(t12373*t12373);
                double t12396 = d1*gamma_bb*t12387;
                double t12388 = t12386-t12396;
                double t12389 = 1.0/sqrt(t12384);
                double t12390 = t12380*t12382*t12374*9.0;
                double t12391 = gamma_bb*gamma_bb;
                double t12392 = t12391*t12385;
                double t12393 = t12392+1.0;
                double t12394 = 1.0/sqrt(t12393);
                double t12395 = gamma_bb*t12382*t12375*t12385*t12394*1.8E1;
                double t12397 = gamma_bb*t12380*t12374*t12377*t12388*1.8E1;
                double t12398 = t12390+t12395+t12397;
                double t12399 = 1.0/pow(t12384,3.0/2.0);
                double t12400 = d1*t12387*2.0;
                double t12401 = 1.0/(t12373*t12373*t12373);
                double t12402 = t12400-d1*gamma_bb*t12401*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * c*pow(rho_b,4.0/3.0)*(t12374*t12378*t12388*t12389*-2.0+gamma_bb*t12402*t12374*t12378*t12389+t12374*t12377*t12378*t12398*t12399+gamma_bb*t12374*t12378*t12388*t12398*t12399+gamma_bb*t12374*t12377*t12378*t12399*(t12382*t12375*t12385*t12394*3.6E1+t12380*t12374*t12377*t12388*3.6E1+gamma_bb*t12380*t12374*(t12388*t12388)*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t12382*1.8E1)/t12393-1.0/pow(rho_b,3.2E1/3.0)*t12382*t12391*t12375*1.0/pow(t12393,3.0/2.0)*1.8E1-gamma_bb*t12402*t12380*t12374*t12377*1.8E1+gamma_bb*t12375*t12385*t12394*t12377*t12388*7.2E1)*(1.0/2.0)-gamma_bb*t12374*1.0/pow(t12384,5.0/2.0)*t12377*t12378*(t12398*t12398)*(3.0/4.0));
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t12333 = 1.0/pow(rho_b,8.0/3.0);
                double t12341 = gamma_bb*t12333;
                double t12334 = log(t12341+sqrt(t12341*t12341+1.0));
                double t12336 = d2*d2;
                double t12337 = gamma_bb+t12336;
                double t12338 = 1.0/t12337;
                double t12339 = d1*gamma_bb*t12338;
                double t12335 = d0+t12339;
                double t12340 = 1.0/c;
                double t12342 = t12334*t12334;
                double t12343 = t12335*t12335;
                double t12344 = gamma_bb*t12333*t12342*t12343*9.0;
                double t12345 = t12344+1.0;
                double t12346 = 1.0/sqrt(t12345);
                double t12347 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12348 = d1*t12338;
                double t12349 = 1.0/(t12337*t12337);
                double t12357 = d1*gamma_bb*t12349;
                double t12350 = t12348-t12357;
                double t12351 = 1.0/pow(rho_b,1.1E1/3.0);
                double t12352 = gamma_bb*gamma_bb;
                double t12353 = t12352*t12347;
                double t12354 = t12353+1.0;
                double t12355 = 1.0/sqrt(t12354);
                double t12356 = 1.0/pow(t12345,3.0/2.0);
                double t12358 = 1.0/pow(rho_b,1.9E1/3.0);
                double t12359 = t12333*t12342*t12343*9.0;
                double t12360 = gamma_bb*t12334*t12343*t12355*t12347*1.8E1;
                double t12361 = gamma_bb*t12350*t12333*t12342*t12335*1.8E1;
                double t12362 = t12360+t12361+t12359;
                double t12363 = gamma_bb*t12342*t12351*t12343*2.4E1;
                double t12364 = t12334*t12343*t12352*t12355*t12358*4.8E1;
                double t12365 = t12363+t12364;
                v_rho_b_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t12340*t12351*t12335*t12346*(-8.0/3.0)-gamma_bb*t12340*t12350*t12351*t12346*(8.0/3.0)+t12340*t12333*t12335*t12356*t12365*(1.0/2.0)+gamma_bb*t12340*t12350*t12333*t12356*t12365*(1.0/2.0)+gamma_bb*t12340*t12351*t12335*t12362*t12356*(4.0/3.0)+gamma_bb*t12340*t12333*t12335*t12356*(t12342*t12351*t12343*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t12343*t12352*4.8E1)/t12354+gamma_bb*t12350*t12342*t12351*t12335*4.8E1+gamma_bb*t12334*t12343*t12355*t12358*1.44E2+t12350*t12334*t12352*t12335*t12355*t12358*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t12334*t12343*t12352*1.0/pow(t12354,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t12340*t12333*t12335*t12362*1.0/pow(t12345,5.0/2.0)*t12365*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t12340*t12333*t12335*t12346+gamma_bb*t12340*t12350*t12333*t12346-gamma_bb*t12340*t12333*t12335*t12362*t12356*(1.0/2.0))*(4.0/3.0);
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t12417 = 1.0/pow(rho_a,8.0/3.0);
                double t12418 = log(gamma_aa*t12417+sqrt((gamma_aa*gamma_aa)*(t12417*t12417)+1.0));
                double t12420 = d2*d2;
                double t12421 = gamma_aa+t12420;
                double t12422 = 1.0/t12421;
                double t12423 = d1*gamma_aa*t12422;
                double t12419 = d0+t12423;
                v[Q] += scale * -c*pow(rho_a,4.0/3.0)*((gamma_aa*t12417*t12419*1.0/sqrt(gamma_aa*t12417*(t12418*t12418)*(t12419*t12419)*9.0+1.0))/c-1.0);
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t12425 = 1.0/pow(rho_a,8.0/3.0);
                double t12433 = gamma_aa*t12425;
                double t12426 = log(t12433+sqrt(t12433*t12433+1.0));
                double t12428 = d2*d2;
                double t12429 = gamma_aa+t12428;
                double t12430 = 1.0/t12429;
                double t12431 = d1*gamma_aa*t12430;
                double t12427 = d0+t12431;
                double t12432 = 1.0/c;
                double t12434 = t12426*t12426;
                double t12435 = t12427*t12427;
                double t12436 = gamma_aa*t12425*t12434*t12435*9.0;
                double t12437 = t12436+1.0;
                double t12438 = 1.0/sqrt(t12437);
                double t12439 = 1.0/pow(rho_a,1.1E1/3.0);
                double t12440 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * c*pow(rho_a,1.0/3.0)*(gamma_aa*t12432*t12425*t12427*t12438-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t12432*t12427*t12438*t12439*(8.0/3.0)-gamma_aa*t12432*t12425*t12427*1.0/pow(t12437,3.0/2.0)*(gamma_aa*t12434*t12435*t12439*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t12440*t12426*t12435*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t12440+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t12443 = 1.0/pow(rho_a,8.0/3.0);
                double t12451 = gamma_aa*t12443;
                double t12444 = log(t12451+sqrt(t12451*t12451+1.0));
                double t12446 = d2*d2;
                double t12447 = gamma_aa+t12446;
                double t12448 = 1.0/t12447;
                double t12449 = d1*gamma_aa*t12448;
                double t12445 = d0+t12449;
                double t12450 = 1.0/c;
                double t12452 = t12444*t12444;
                double t12453 = t12445*t12445;
                double t12454 = gamma_aa*t12443*t12452*t12453*9.0;
                double t12455 = t12454+1.0;
                double t12456 = 1.0/sqrt(t12455);
                double t12457 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12458 = d1*t12448;
                double t12459 = 1.0/(t12447*t12447);
                double t12460 = t12458-d1*gamma_aa*t12459;
                v_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t12450*t12443*t12445*t12456+gamma_aa*t12450*t12460*t12443*t12456-gamma_aa*t12450*t12443*t12445*1.0/pow(t12455,3.0/2.0)*(t12443*t12452*t12453*9.0+gamma_aa*t12444*t12453*t12457*1.0/sqrt((gamma_aa*gamma_aa)*t12457+1.0)*1.8E1+gamma_aa*t12460*t12443*t12452*t12445*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t12466 = 1.0/pow(rho_a,8.0/3.0);
                double t12474 = gamma_aa*t12466;
                double t12467 = log(t12474+sqrt(t12474*t12474+1.0));
                double t12469 = d2*d2;
                double t12470 = gamma_aa+t12469;
                double t12471 = 1.0/t12470;
                double t12472 = d1*gamma_aa*t12471;
                double t12468 = d0+t12472;
                double t12473 = 1.0/c;
                double t12475 = t12467*t12467;
                double t12476 = t12468*t12468;
                double t12477 = gamma_aa*t12466*t12475*t12476*9.0;
                double t12478 = t12477+1.0;
                double t12479 = 1.0/sqrt(t12478);
                double t12480 = 1.0/pow(rho_a,1.1E1/3.0);
                double t12481 = gamma_aa*gamma_aa;
                double t12482 = 1.0/pow(t12478,3.0/2.0);
                double t12483 = 1.0/pow(rho_a,1.4E1/3.0);
                double t12484 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12485 = t12481*t12484;
                double t12486 = t12485+1.0;
                double t12487 = 1.0/sqrt(t12486);
                double t12488 = gamma_aa*t12480*t12475*t12476*2.4E1;
                double t12489 = 1.0/pow(rho_a,1.9E1/3.0);
                double t12490 = t12481*t12467*t12476*t12487*t12489*4.8E1;
                double t12491 = t12490+t12488;
                v_rho_a_rho_a[Q] += scale * c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t12473*t12466*t12468*t12479-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t12473*t12483*t12468*t12479*(8.8E1/9.0)-gamma_aa*t12480*t12473*t12482*t12491*t12468*(8.0/3.0)-gamma_aa*t12473*t12482*t12466*t12468*(gamma_aa*t12483*t12475*t12476*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t12481*t12467*t12476*t12487*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t12481*t12476*1.28E2)/t12486-1.0/pow(rho_a,3.8E1/3.0)*(t12481*t12481)*t12467*t12476*1.0/pow(t12486,3.0/2.0)*1.28E2)*(1.0/2.0)+gamma_aa*t12473*(t12491*t12491)*t12466*t12468*1.0/pow(t12478,5.0/2.0)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t12480*t12473*t12468*t12479*(8.0/3.0)-gamma_aa*t12473*t12482*t12491*t12466*t12468*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t12534 = d2*d2;
                double t12535 = gamma_aa+t12534;
                double t12536 = 1.0/pow(rho_a,8.0/3.0);
                double t12541 = gamma_aa*t12536;
                double t12537 = log(t12541+sqrt(t12541*t12541+1.0));
                double t12538 = 1.0/t12535;
                double t12543 = d1*gamma_aa*t12538;
                double t12539 = d0+t12543;
                double t12540 = 1.0/c;
                double t12542 = t12537*t12537;
                double t12544 = t12539*t12539;
                double t12545 = gamma_aa*t12542*t12544*t12536*9.0;
                double t12546 = t12545+1.0;
                double t12547 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12548 = d1*t12538;
                double t12549 = 1.0/(t12535*t12535);
                double t12558 = d1*gamma_aa*t12549;
                double t12550 = t12548-t12558;
                double t12551 = 1.0/sqrt(t12546);
                double t12552 = t12542*t12544*t12536*9.0;
                double t12553 = gamma_aa*gamma_aa;
                double t12554 = t12553*t12547;
                double t12555 = t12554+1.0;
                double t12556 = 1.0/sqrt(t12555);
                double t12557 = gamma_aa*t12544*t12537*t12547*t12556*1.8E1;
                double t12559 = gamma_aa*t12550*t12542*t12536*t12539*1.8E1;
                double t12560 = t12552+t12557+t12559;
                double t12561 = 1.0/pow(t12546,3.0/2.0);
                double t12562 = d1*t12549*2.0;
                double t12563 = 1.0/(t12535*t12535*t12535);
                double t12564 = t12562-d1*gamma_aa*t12563*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * c*pow(rho_a,4.0/3.0)*(t12540*t12550*t12551*t12536*-2.0+gamma_aa*t12540*t12551*t12536*t12564+t12540*t12560*t12561*t12536*t12539+gamma_aa*t12540*t12550*t12560*t12561*t12536+gamma_aa*t12540*t12561*t12536*t12539*(t12550*t12542*t12536*t12539*3.6E1+t12544*t12537*t12547*t12556*3.6E1+gamma_aa*(t12550*t12550)*t12542*t12536*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t12544*1.8E1)/t12555-1.0/pow(rho_a,3.2E1/3.0)*t12544*t12553*t12537*1.0/pow(t12555,3.0/2.0)*1.8E1-gamma_aa*t12542*t12536*t12564*t12539*1.8E1+gamma_aa*t12550*t12537*t12547*t12556*t12539*7.2E1)*(1.0/2.0)-gamma_aa*t12540*(t12560*t12560)*t12536*1.0/pow(t12546,5.0/2.0)*t12539*(3.0/4.0));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t12495 = 1.0/pow(rho_a,8.0/3.0);
                double t12503 = gamma_aa*t12495;
                double t12496 = log(t12503+sqrt(t12503*t12503+1.0));
                double t12498 = d2*d2;
                double t12499 = gamma_aa+t12498;
                double t12500 = 1.0/t12499;
                double t12501 = d1*gamma_aa*t12500;
                double t12497 = d0+t12501;
                double t12502 = 1.0/c;
                double t12504 = t12496*t12496;
                double t12505 = t12497*t12497;
                double t12506 = gamma_aa*t12504*t12505*t12495*9.0;
                double t12507 = t12506+1.0;
                double t12508 = 1.0/sqrt(t12507);
                double t12509 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12510 = d1*t12500;
                double t12511 = 1.0/(t12499*t12499);
                double t12519 = d1*gamma_aa*t12511;
                double t12512 = t12510-t12519;
                double t12513 = 1.0/pow(rho_a,1.1E1/3.0);
                double t12514 = gamma_aa*gamma_aa;
                double t12515 = t12514*t12509;
                double t12516 = t12515+1.0;
                double t12517 = 1.0/sqrt(t12516);
                double t12518 = 1.0/pow(t12507,3.0/2.0);
                double t12520 = 1.0/pow(rho_a,1.9E1/3.0);
                double t12521 = t12504*t12505*t12495*9.0;
                double t12522 = gamma_aa*t12505*t12517*t12509*t12496*1.8E1;
                double t12523 = gamma_aa*t12512*t12504*t12495*t12497*1.8E1;
                double t12524 = t12521+t12522+t12523;
                double t12525 = gamma_aa*t12504*t12513*t12505*2.4E1;
                double t12526 = t12520*t12505*t12514*t12517*t12496*4.8E1;
                double t12527 = t12525+t12526;
                v_rho_a_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t12502*t12513*t12508*t12497*(-8.0/3.0)-gamma_aa*t12502*t12512*t12513*t12508*(8.0/3.0)+t12502*t12518*t12527*t12495*t12497*(1.0/2.0)+gamma_aa*t12502*t12512*t12518*t12527*t12495*(1.0/2.0)+gamma_aa*t12502*t12513*t12524*t12518*t12497*(4.0/3.0)+gamma_aa*t12502*t12518*t12495*t12497*(t12504*t12513*t12505*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t12505*t12514*4.8E1)/t12516+gamma_aa*t12512*t12504*t12513*t12497*4.8E1+gamma_aa*t12520*t12505*t12517*t12496*1.44E2+t12520*t12512*t12514*t12517*t12496*t12497*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t12505*t12514*1.0/pow(t12516,3.0/2.0)*t12496*4.8E1)*(1.0/2.0)-gamma_aa*t12502*t12524*1.0/pow(t12507,5.0/2.0)*t12527*t12495*t12497*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t12502*t12508*t12495*t12497+gamma_aa*t12502*t12512*t12508*t12495-gamma_aa*t12502*t12524*t12518*t12495*t12497*(1.0/2.0))*(4.0/3.0);
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t11952 = 1.0/pow(rho_a,8.0/3.0);
                double t11953 = log(gamma_aa*t11952+sqrt((gamma_aa*gamma_aa)*(t11952*t11952)+1.0));
                double t11955 = d2*d2;
                double t11956 = gamma_aa+t11955;
                double t11957 = 1.0/t11956;
                double t11958 = d1*gamma_aa*t11957;
                double t11954 = d0+t11958;
                double t11959 = 1.0/c;
                double t11960 = 1.0/pow(rho_b,8.0/3.0);
                double t11961 = log(gamma_bb*t11960+sqrt((gamma_bb*gamma_bb)*(t11960*t11960)+1.0));
                double t11963 = gamma_bb+t11955;
                double t11964 = 1.0/t11963;
                double t11965 = d1*gamma_bb*t11964;
                double t11962 = d0+t11965;
                v[Q] += scale * -c*pow(rho_a,4.0/3.0)*(gamma_aa*t11952*t11954*t11959*1.0/sqrt(gamma_aa*t11952*(t11953*t11953)*(t11954*t11954)*9.0+1.0)-1.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t11960*t11962*t11959*1.0/sqrt(gamma_bb*t11960*(t11961*t11961)*(t11962*t11962)*9.0+1.0)-1.0);
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t11967 = 1.0/pow(rho_a,8.0/3.0);
                double t11975 = gamma_aa*t11967;
                double t11968 = log(t11975+sqrt(t11975*t11975+1.0));
                double t11970 = d2*d2;
                double t11971 = gamma_aa+t11970;
                double t11972 = 1.0/t11971;
                double t11973 = d1*gamma_aa*t11972;
                double t11969 = d0+t11973;
                double t11974 = 1.0/c;
                double t11976 = t11968*t11968;
                double t11977 = t11969*t11969;
                double t11978 = gamma_aa*t11967*t11976*t11977*9.0;
                double t11979 = t11978+1.0;
                double t11980 = 1.0/sqrt(t11979);
                double t11981 = 1.0/pow(rho_a,1.1E1/3.0);
                double t11982 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * c*pow(rho_a,1.0/3.0)*(gamma_aa*t11980*t11974*t11967*t11969-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t11980*t11981*t11974*t11969*(8.0/3.0)-gamma_aa*t11974*t11967*t11969*1.0/pow(t11979,3.0/2.0)*(gamma_aa*t11981*t11976*t11977*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t11982*t11968*t11977*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t11982+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t11984 = 1.0/pow(rho_b,8.0/3.0);
                double t11992 = gamma_bb*t11984;
                double t11985 = log(t11992+sqrt(t11992*t11992+1.0));
                double t11987 = d2*d2;
                double t11988 = gamma_bb+t11987;
                double t11989 = 1.0/t11988;
                double t11990 = d1*gamma_bb*t11989;
                double t11986 = d0+t11990;
                double t11991 = 1.0/c;
                double t11993 = t11985*t11985;
                double t11994 = t11986*t11986;
                double t11995 = gamma_bb*t11984*t11993*t11994*9.0;
                double t11996 = t11995+1.0;
                double t11997 = 1.0/sqrt(t11996);
                double t11998 = 1.0/pow(rho_b,1.1E1/3.0);
                double t11999 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * c*pow(rho_b,1.0/3.0)*(gamma_bb*t11991*t11984*t11986*t11997-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t11991*t11986*t11997*t11998*(8.0/3.0)-gamma_bb*t11991*t11984*t11986*1.0/pow(t11996,3.0/2.0)*(gamma_bb*t11993*t11994*t11998*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t11985*t11994*t11999*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t11999+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t12001 = 1.0/pow(rho_a,8.0/3.0);
                double t12009 = gamma_aa*t12001;
                double t12002 = log(t12009+sqrt(t12009*t12009+1.0));
                double t12004 = d2*d2;
                double t12005 = gamma_aa+t12004;
                double t12006 = 1.0/t12005;
                double t12007 = d1*gamma_aa*t12006;
                double t12003 = d0+t12007;
                double t12008 = 1.0/c;
                double t12010 = t12002*t12002;
                double t12011 = t12003*t12003;
                double t12012 = gamma_aa*t12001*t12010*t12011*9.0;
                double t12013 = t12012+1.0;
                double t12014 = 1.0/sqrt(t12013);
                double t12015 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12016 = d1*t12006;
                double t12017 = 1.0/(t12005*t12005);
                double t12018 = t12016-d1*gamma_aa*t12017;
                v_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t12001*t12003*t12014*t12008+gamma_aa*t12001*t12014*t12008*t12018-gamma_aa*t12001*t12003*1.0/pow(t12013,3.0/2.0)*t12008*(t12001*t12010*t12011*9.0+gamma_aa*t12002*t12011*t12015*1.0/sqrt((gamma_aa*gamma_aa)*t12015+1.0)*1.8E1+gamma_aa*t12001*t12010*t12003*t12018*1.8E1)*(1.0/2.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t12021 = 1.0/pow(rho_b,8.0/3.0);
                double t12029 = gamma_bb*t12021;
                double t12022 = log(t12029+sqrt(t12029*t12029+1.0));
                double t12024 = d2*d2;
                double t12025 = gamma_bb+t12024;
                double t12026 = 1.0/t12025;
                double t12027 = d1*gamma_bb*t12026;
                double t12023 = d0+t12027;
                double t12028 = 1.0/c;
                double t12030 = t12022*t12022;
                double t12031 = t12023*t12023;
                double t12032 = gamma_bb*t12021*t12030*t12031*9.0;
                double t12033 = t12032+1.0;
                double t12034 = 1.0/sqrt(t12033);
                double t12035 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12036 = d1*t12026;
                double t12037 = 1.0/(t12025*t12025);
                double t12038 = t12036-d1*gamma_bb*t12037;
                v_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t12021*t12023*t12034*t12028+gamma_bb*t12021*t12034*t12028*t12038-gamma_bb*t12021*t12023*1.0/pow(t12033,3.0/2.0)*t12028*(t12021*t12030*t12031*9.0+gamma_bb*t12022*t12031*t12035*1.0/sqrt((gamma_bb*gamma_bb)*t12035+1.0)*1.8E1+gamma_bb*t12021*t12030*t12023*t12038*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t12042 = 1.0/pow(rho_a,8.0/3.0);
                double t12050 = gamma_aa*t12042;
                double t12043 = log(t12050+sqrt(t12050*t12050+1.0));
                double t12045 = d2*d2;
                double t12046 = gamma_aa+t12045;
                double t12047 = 1.0/t12046;
                double t12048 = d1*gamma_aa*t12047;
                double t12044 = d0+t12048;
                double t12049 = 1.0/c;
                double t12051 = t12043*t12043;
                double t12052 = t12044*t12044;
                double t12053 = gamma_aa*t12042*t12051*t12052*9.0;
                double t12054 = t12053+1.0;
                double t12055 = 1.0/sqrt(t12054);
                double t12056 = 1.0/pow(rho_a,1.1E1/3.0);
                double t12057 = gamma_aa*gamma_aa;
                double t12058 = 1.0/pow(t12054,3.0/2.0);
                double t12059 = 1.0/pow(rho_a,1.4E1/3.0);
                double t12060 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12061 = t12060*t12057;
                double t12062 = t12061+1.0;
                double t12063 = 1.0/sqrt(t12062);
                double t12064 = gamma_aa*t12051*t12052*t12056*2.4E1;
                double t12065 = 1.0/pow(rho_a,1.9E1/3.0);
                double t12066 = t12043*t12052*t12063*t12065*t12057*4.8E1;
                double t12067 = t12064+t12066;
                v_rho_a_rho_a[Q] += scale * c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t12042*t12044*t12055*t12049-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t12044*t12055*t12049*t12059*(8.8E1/9.0)-gamma_aa*t12044*t12056*t12049*t12058*t12067*(8.0/3.0)-gamma_aa*t12042*t12044*t12049*t12058*(gamma_aa*t12051*t12052*t12059*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t12043*t12052*t12063*t12057*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t12052*t12057*1.28E2)/t12062-1.0/pow(rho_a,3.8E1/3.0)*t12043*t12052*1.0/pow(t12062,3.0/2.0)*(t12057*t12057)*1.28E2)*(1.0/2.0)+gamma_aa*t12042*t12044*1.0/pow(t12054,5.0/2.0)*t12049*(t12067*t12067)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t12044*t12055*t12056*t12049*(8.0/3.0)-gamma_aa*t12042*t12044*t12049*t12058*t12067*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t12070 = 1.0/pow(rho_b,8.0/3.0);
                double t12078 = gamma_bb*t12070;
                double t12071 = log(t12078+sqrt(t12078*t12078+1.0));
                double t12073 = d2*d2;
                double t12074 = gamma_bb+t12073;
                double t12075 = 1.0/t12074;
                double t12076 = d1*gamma_bb*t12075;
                double t12072 = d0+t12076;
                double t12077 = 1.0/c;
                double t12079 = t12071*t12071;
                double t12080 = t12072*t12072;
                double t12081 = gamma_bb*t12070*t12080*t12079*9.0;
                double t12082 = t12081+1.0;
                double t12083 = 1.0/sqrt(t12082);
                double t12084 = 1.0/pow(rho_b,1.1E1/3.0);
                double t12085 = gamma_bb*gamma_bb;
                double t12086 = 1.0/pow(t12082,3.0/2.0);
                double t12087 = 1.0/pow(rho_b,1.4E1/3.0);
                double t12088 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12089 = t12085*t12088;
                double t12090 = t12089+1.0;
                double t12091 = 1.0/sqrt(t12090);
                double t12092 = gamma_bb*t12080*t12084*t12079*2.4E1;
                double t12093 = 1.0/pow(rho_b,1.9E1/3.0);
                double t12094 = t12071*t12080*t12091*t12093*t12085*4.8E1;
                double t12095 = t12092+t12094;
                v_rho_b_rho_b[Q] += scale * c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t12070*t12072*t12083*t12077-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t12072*t12083*t12077*t12087*(8.8E1/9.0)-gamma_bb*t12072*t12084*t12077*t12086*t12095*(8.0/3.0)-gamma_bb*t12070*t12072*t12077*t12086*(gamma_bb*t12080*t12087*t12079*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t12071*t12080*t12091*t12085*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t12080*t12085*1.28E2)/t12090-1.0/pow(rho_b,3.8E1/3.0)*t12071*t12080*1.0/pow(t12090,3.0/2.0)*(t12085*t12085)*1.28E2)*(1.0/2.0)+gamma_bb*t12070*t12072*1.0/pow(t12082,5.0/2.0)*t12077*(t12095*t12095)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t12072*t12083*t12084*t12077*(8.0/3.0)-gamma_bb*t12070*t12072*t12077*t12086*t12095*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t12169 = d2*d2;
                double t12170 = gamma_aa+t12169;
                double t12171 = 1.0/pow(rho_a,8.0/3.0);
                double t12176 = gamma_aa*t12171;
                double t12172 = log(t12176+sqrt(t12176*t12176+1.0));
                double t12173 = 1.0/t12170;
                double t12178 = d1*gamma_aa*t12173;
                double t12174 = d0+t12178;
                double t12175 = 1.0/c;
                double t12177 = t12172*t12172;
                double t12179 = t12174*t12174;
                double t12180 = gamma_aa*t12171*t12177*t12179*9.0;
                double t12181 = t12180+1.0;
                double t12182 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12183 = d1*t12173;
                double t12184 = 1.0/(t12170*t12170);
                double t12193 = d1*gamma_aa*t12184;
                double t12185 = t12183-t12193;
                double t12186 = 1.0/sqrt(t12181);
                double t12187 = t12171*t12177*t12179*9.0;
                double t12188 = gamma_aa*gamma_aa;
                double t12189 = t12182*t12188;
                double t12190 = t12189+1.0;
                double t12191 = 1.0/sqrt(t12190);
                double t12192 = gamma_aa*t12172*t12182*t12191*t12179*1.8E1;
                double t12194 = gamma_aa*t12171*t12174*t12185*t12177*1.8E1;
                double t12195 = t12192+t12194+t12187;
                double t12196 = 1.0/pow(t12181,3.0/2.0);
                double t12197 = d1*t12184*2.0;
                double t12198 = 1.0/(t12170*t12170*t12170);
                double t12199 = t12197-d1*gamma_aa*t12198*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * c*pow(rho_a,4.0/3.0)*(t12171*t12175*t12185*t12186*-2.0+gamma_aa*t12171*t12175*t12186*t12199+t12171*t12174*t12175*t12195*t12196+gamma_aa*t12171*t12175*t12185*t12195*t12196+gamma_aa*t12171*t12174*t12175*t12196*(t12172*t12182*t12191*t12179*3.6E1+t12171*t12174*t12185*t12177*3.6E1+gamma_aa*t12171*(t12185*t12185)*t12177*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t12179*1.8E1)/t12190-1.0/pow(rho_a,3.2E1/3.0)*t12172*1.0/pow(t12190,3.0/2.0)*t12179*t12188*1.8E1-gamma_aa*t12171*t12174*t12177*t12199*1.8E1+gamma_aa*t12172*t12182*t12191*t12174*t12185*7.2E1)*(1.0/2.0)-gamma_aa*t12171*1.0/pow(t12181,5.0/2.0)*t12174*t12175*(t12195*t12195)*(3.0/4.0));
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t12205 = d2*d2;
                double t12206 = gamma_bb+t12205;
                double t12207 = 1.0/pow(rho_b,8.0/3.0);
                double t12212 = gamma_bb*t12207;
                double t12208 = log(t12212+sqrt(t12212*t12212+1.0));
                double t12209 = 1.0/t12206;
                double t12214 = d1*gamma_bb*t12209;
                double t12210 = d0+t12214;
                double t12211 = 1.0/c;
                double t12213 = t12208*t12208;
                double t12215 = t12210*t12210;
                double t12216 = gamma_bb*t12213*t12215*t12207*9.0;
                double t12217 = t12216+1.0;
                double t12218 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12219 = d1*t12209;
                double t12220 = 1.0/(t12206*t12206);
                double t12229 = d1*gamma_bb*t12220;
                double t12221 = t12219-t12229;
                double t12222 = 1.0/sqrt(t12217);
                double t12223 = t12213*t12215*t12207*9.0;
                double t12224 = gamma_bb*gamma_bb;
                double t12225 = t12224*t12218;
                double t12226 = t12225+1.0;
                double t12227 = 1.0/sqrt(t12226);
                double t12228 = gamma_bb*t12215*t12208*t12218*t12227*1.8E1;
                double t12230 = gamma_bb*t12210*t12221*t12213*t12207*1.8E1;
                double t12231 = t12230+t12223+t12228;
                double t12232 = 1.0/pow(t12217,3.0/2.0);
                double t12233 = d1*t12220*2.0;
                double t12234 = 1.0/(t12206*t12206*t12206);
                double t12235 = t12233-d1*gamma_bb*t12234*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * c*pow(rho_b,4.0/3.0)*(t12211*t12221*t12222*t12207*-2.0+gamma_bb*t12211*t12222*t12207*t12235+t12210*t12211*t12231*t12232*t12207+gamma_bb*t12211*t12221*t12231*t12232*t12207+gamma_bb*t12210*t12211*t12232*t12207*(t12210*t12221*t12213*t12207*3.6E1+t12215*t12208*t12218*t12227*3.6E1+gamma_bb*(t12221*t12221)*t12213*t12207*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t12215*1.8E1)/t12226-1.0/pow(rho_b,3.2E1/3.0)*t12215*t12224*t12208*1.0/pow(t12226,3.0/2.0)*1.8E1-gamma_bb*t12210*t12213*t12207*t12235*1.8E1+gamma_bb*t12210*t12221*t12208*t12218*t12227*7.2E1)*(1.0/2.0)-gamma_bb*t12210*t12211*(t12231*t12231)*t12207*1.0/pow(t12217,5.0/2.0)*(3.0/4.0));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t12097 = 1.0/pow(rho_a,8.0/3.0);
                double t12105 = gamma_aa*t12097;
                double t12098 = log(t12105+sqrt(t12105*t12105+1.0));
                double t12100 = d2*d2;
                double t12101 = gamma_aa+t12100;
                double t12102 = 1.0/t12101;
                double t12103 = d1*gamma_aa*t12102;
                double t12099 = d0+t12103;
                double t12104 = 1.0/c;
                double t12106 = t12098*t12098;
                double t12107 = t12099*t12099;
                double t12108 = gamma_aa*t12106*t12107*t12097*9.0;
                double t12109 = t12108+1.0;
                double t12110 = 1.0/sqrt(t12109);
                double t12111 = 1.0/pow(rho_a,1.6E1/3.0);
                double t12112 = d1*t12102;
                double t12113 = 1.0/(t12101*t12101);
                double t12121 = d1*gamma_aa*t12113;
                double t12114 = t12112-t12121;
                double t12115 = 1.0/pow(rho_a,1.1E1/3.0);
                double t12116 = gamma_aa*gamma_aa;
                double t12117 = t12111*t12116;
                double t12118 = t12117+1.0;
                double t12119 = 1.0/sqrt(t12118);
                double t12120 = 1.0/pow(t12109,3.0/2.0);
                double t12122 = 1.0/pow(rho_a,1.9E1/3.0);
                double t12123 = t12106*t12107*t12097*9.0;
                double t12124 = gamma_aa*t12111*t12107*t12119*t12098*1.8E1;
                double t12125 = gamma_aa*t12114*t12106*t12097*t12099*1.8E1;
                double t12126 = t12123+t12124+t12125;
                double t12127 = gamma_aa*t12106*t12115*t12107*2.4E1;
                double t12128 = t12122*t12107*t12116*t12119*t12098*4.8E1;
                double t12129 = t12127+t12128;
                v_rho_a_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t12110*t12104*t12115*t12099*(-8.0/3.0)-gamma_aa*t12110*t12104*t12114*t12115*(8.0/3.0)+t12120*t12104*t12129*t12097*t12099*(1.0/2.0)+gamma_aa*t12120*t12104*t12114*t12129*t12097*(1.0/2.0)+gamma_aa*t12120*t12104*t12115*t12126*t12099*(4.0/3.0)+gamma_aa*t12120*t12104*t12097*t12099*(t12106*t12115*t12107*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t12107*t12116*4.8E1)/t12118+gamma_aa*t12114*t12106*t12115*t12099*4.8E1+gamma_aa*t12122*t12107*t12119*t12098*1.44E2+t12122*t12114*t12116*t12119*t12098*t12099*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t12107*t12116*1.0/pow(t12118,3.0/2.0)*t12098*4.8E1)*(1.0/2.0)-gamma_aa*t12104*t12126*1.0/pow(t12109,5.0/2.0)*t12129*t12097*t12099*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t12110*t12104*t12097*t12099+gamma_aa*t12110*t12104*t12114*t12097-gamma_aa*t12120*t12104*t12126*t12097*t12099*(1.0/2.0))*(4.0/3.0);
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t12135 = 1.0/pow(rho_b,8.0/3.0);
                double t12143 = gamma_bb*t12135;
                double t12136 = log(t12143+sqrt(t12143*t12143+1.0));
                double t12138 = d2*d2;
                double t12139 = gamma_bb+t12138;
                double t12140 = 1.0/t12139;
                double t12141 = d1*gamma_bb*t12140;
                double t12137 = d0+t12141;
                double t12142 = 1.0/c;
                double t12144 = t12136*t12136;
                double t12145 = t12137*t12137;
                double t12146 = gamma_bb*t12135*t12144*t12145*9.0;
                double t12147 = t12146+1.0;
                double t12148 = 1.0/sqrt(t12147);
                double t12149 = 1.0/pow(rho_b,1.6E1/3.0);
                double t12150 = d1*t12140;
                double t12151 = 1.0/(t12139*t12139);
                double t12159 = d1*gamma_bb*t12151;
                double t12152 = t12150-t12159;
                double t12153 = 1.0/pow(rho_b,1.1E1/3.0);
                double t12154 = gamma_bb*gamma_bb;
                double t12155 = t12154*t12149;
                double t12156 = t12155+1.0;
                double t12157 = 1.0/sqrt(t12156);
                double t12158 = 1.0/pow(t12147,3.0/2.0);
                double t12160 = 1.0/pow(rho_b,1.9E1/3.0);
                double t12161 = t12135*t12144*t12145*9.0;
                double t12162 = gamma_bb*t12136*t12145*t12157*t12149*1.8E1;
                double t12163 = gamma_bb*t12152*t12135*t12144*t12137*1.8E1;
                double t12164 = t12161+t12162+t12163;
                double t12165 = gamma_bb*t12144*t12153*t12145*2.4E1;
                double t12166 = t12160*t12136*t12145*t12154*t12157*4.8E1;
                double t12167 = t12165+t12166;
                v_rho_b_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t12142*t12153*t12137*t12148*(-8.0/3.0)-gamma_bb*t12142*t12152*t12153*t12148*(8.0/3.0)+t12142*t12135*t12137*t12158*t12167*(1.0/2.0)+gamma_bb*t12142*t12152*t12135*t12158*t12167*(1.0/2.0)+gamma_bb*t12142*t12153*t12137*t12164*t12158*(4.0/3.0)+gamma_bb*t12142*t12135*t12137*t12158*(t12144*t12153*t12145*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t12145*t12154*4.8E1)/t12156+gamma_bb*t12152*t12144*t12153*t12137*4.8E1+gamma_bb*t12160*t12136*t12145*t12157*1.44E2+t12160*t12152*t12136*t12154*t12137*t12157*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t12136*t12145*t12154*1.0/pow(t12156,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t12142*t12135*t12137*t12164*1.0/pow(t12147,5.0/2.0)*t12167*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t12142*t12135*t12137*t12148+gamma_bb*t12142*t12152*t12135*t12148-gamma_bb*t12142*t12135*t12137*t12164*t12158*(1.0/2.0))*(4.0/3.0);
            }
            
        }
    }
}

}
