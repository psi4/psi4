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
                double t38391 = 1.0/pow(rho_b,8.0/3.0);
                double t38392 = log(gamma_bb*t38391+sqrt((gamma_bb*gamma_bb)*(t38391*t38391)+1.0));
                double t38394 = d2*d2;
                double t38395 = gamma_bb+t38394;
                double t38396 = 1.0/t38395;
                double t38397 = d1*gamma_bb*t38396;
                double t38393 = d0+t38397;
                v[Q] += scale * -c*pow(rho_b,4.0/3.0)*((gamma_bb*t38391*t38393*1.0/sqrt(gamma_bb*t38391*(t38392*t38392)*(t38393*t38393)*9.0+1.0))/c-1.0);
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t38400 = 1.0/pow(rho_b,8.0/3.0);
                double t38408 = gamma_bb*t38400;
                double t38401 = log(t38408+sqrt(t38408*t38408+1.0));
                double t38403 = d2*d2;
                double t38404 = gamma_bb+t38403;
                double t38405 = 1.0/t38404;
                double t38406 = d1*gamma_bb*t38405;
                double t38402 = d0+t38406;
                double t38407 = 1.0/c;
                double t38409 = t38401*t38401;
                double t38410 = t38402*t38402;
                double t38411 = gamma_bb*t38400*t38410*t38409*9.0;
                double t38412 = t38411+1.0;
                double t38413 = 1.0/sqrt(t38412);
                double t38414 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38415 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * c*pow(rho_b,1.0/3.0)*(gamma_bb*t38400*t38402*t38413*t38407-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t38402*t38413*t38414*t38407*(8.0/3.0)-gamma_bb*t38400*t38402*1.0/pow(t38412,3.0/2.0)*t38407*(gamma_bb*t38410*t38414*t38409*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t38401*t38410*t38415*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t38415+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t38419 = 1.0/pow(rho_b,8.0/3.0);
                double t38427 = gamma_bb*t38419;
                double t38420 = log(t38427+sqrt(t38427*t38427+1.0));
                double t38422 = d2*d2;
                double t38423 = gamma_bb+t38422;
                double t38424 = 1.0/t38423;
                double t38425 = d1*gamma_bb*t38424;
                double t38421 = d0+t38425;
                double t38426 = 1.0/c;
                double t38428 = t38420*t38420;
                double t38429 = t38421*t38421;
                double t38430 = gamma_bb*t38419*t38428*t38429*9.0;
                double t38431 = t38430+1.0;
                double t38432 = 1.0/sqrt(t38431);
                double t38433 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38434 = d1*t38424;
                double t38435 = 1.0/(t38423*t38423);
                double t38436 = t38434-d1*gamma_bb*t38435;
                v_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t38421*t38432*t38426*t38419+gamma_bb*t38432*t38426*t38436*t38419-gamma_bb*t38421*1.0/pow(t38431,3.0/2.0)*t38426*t38419*(t38419*t38428*t38429*9.0+gamma_bb*t38420*t38433*t38429*1.0/sqrt((gamma_bb*gamma_bb)*t38433+1.0)*1.8E1+gamma_bb*t38421*t38436*t38419*t38428*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t38442 = 1.0/pow(rho_b,8.0/3.0);
                double t38450 = gamma_bb*t38442;
                double t38443 = log(t38450+sqrt(t38450*t38450+1.0));
                double t38445 = d2*d2;
                double t38446 = gamma_bb+t38445;
                double t38447 = 1.0/t38446;
                double t38448 = d1*gamma_bb*t38447;
                double t38444 = d0+t38448;
                double t38449 = 1.0/c;
                double t38451 = t38443*t38443;
                double t38452 = t38444*t38444;
                double t38453 = gamma_bb*t38442*t38451*t38452*9.0;
                double t38454 = t38453+1.0;
                double t38455 = 1.0/sqrt(t38454);
                double t38456 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38457 = gamma_bb*gamma_bb;
                double t38458 = 1.0/pow(t38454,3.0/2.0);
                double t38459 = 1.0/pow(rho_b,1.4E1/3.0);
                double t38460 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38461 = t38460*t38457;
                double t38462 = t38461+1.0;
                double t38463 = 1.0/sqrt(t38462);
                double t38464 = gamma_bb*t38451*t38452*t38456*2.4E1;
                double t38465 = 1.0/pow(rho_b,1.9E1/3.0);
                double t38466 = t38443*t38452*t38463*t38465*t38457*4.8E1;
                double t38467 = t38464+t38466;
                v_rho_b_rho_b[Q] += scale * c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t38442*t38444*t38455*t38449-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t38444*t38455*t38449*t38459*(8.8E1/9.0)-gamma_bb*t38444*t38456*t38449*t38458*t38467*(8.0/3.0)-gamma_bb*t38442*t38444*t38449*t38458*(gamma_bb*t38451*t38452*t38459*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t38443*t38452*t38463*t38457*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t38452*t38457*1.28E2)/t38462-1.0/pow(rho_b,3.8E1/3.0)*t38443*t38452*1.0/pow(t38462,3.0/2.0)*(t38457*t38457)*1.28E2)*(1.0/2.0)+gamma_bb*t38442*t38444*1.0/pow(t38454,5.0/2.0)*t38449*(t38467*t38467)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t38444*t38455*t38456*t38449*(8.0/3.0)-gamma_bb*t38442*t38444*t38449*t38458*t38467*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t38513 = d2*d2;
                double t38514 = gamma_bb+t38513;
                double t38515 = 1.0/pow(rho_b,8.0/3.0);
                double t38520 = gamma_bb*t38515;
                double t38516 = log(t38520+sqrt(t38520*t38520+1.0));
                double t38517 = 1.0/t38514;
                double t38522 = d1*gamma_bb*t38517;
                double t38518 = d0+t38522;
                double t38519 = 1.0/c;
                double t38521 = t38516*t38516;
                double t38523 = t38518*t38518;
                double t38524 = gamma_bb*t38521*t38523*t38515*9.0;
                double t38525 = t38524+1.0;
                double t38526 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38527 = d1*t38517;
                double t38528 = 1.0/(t38514*t38514);
                double t38537 = d1*gamma_bb*t38528;
                double t38529 = t38527-t38537;
                double t38530 = 1.0/sqrt(t38525);
                double t38531 = t38521*t38523*t38515*9.0;
                double t38532 = gamma_bb*gamma_bb;
                double t38533 = t38532*t38526;
                double t38534 = t38533+1.0;
                double t38535 = 1.0/sqrt(t38534);
                double t38536 = gamma_bb*t38523*t38516*t38526*t38535*1.8E1;
                double t38538 = gamma_bb*t38521*t38515*t38518*t38529*1.8E1;
                double t38539 = t38531+t38536+t38538;
                double t38540 = 1.0/pow(t38525,3.0/2.0);
                double t38541 = d1*t38528*2.0;
                double t38542 = 1.0/(t38514*t38514*t38514);
                double t38543 = t38541-d1*gamma_bb*t38542*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * c*pow(rho_b,4.0/3.0)*(t38530*t38515*t38519*t38529*-2.0+gamma_bb*t38530*t38515*t38543*t38519+t38540*t38515*t38518*t38519*t38539+gamma_bb*t38540*t38515*t38519*t38529*t38539+gamma_bb*t38540*t38515*t38518*t38519*(t38523*t38516*t38526*t38535*3.6E1+t38521*t38515*t38518*t38529*3.6E1+gamma_bb*t38521*t38515*(t38529*t38529)*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t38523*1.8E1)/t38534-1.0/pow(rho_b,3.2E1/3.0)*t38523*t38532*t38516*1.0/pow(t38534,3.0/2.0)*1.8E1-gamma_bb*t38521*t38515*t38543*t38518*1.8E1+gamma_bb*t38516*t38526*t38535*t38518*t38529*7.2E1)*(1.0/2.0)-gamma_bb*t38515*1.0/pow(t38525,5.0/2.0)*t38518*t38519*(t38539*t38539)*(3.0/4.0));
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t38474 = 1.0/pow(rho_b,8.0/3.0);
                double t38482 = gamma_bb*t38474;
                double t38475 = log(t38482+sqrt(t38482*t38482+1.0));
                double t38477 = d2*d2;
                double t38478 = gamma_bb+t38477;
                double t38479 = 1.0/t38478;
                double t38480 = d1*gamma_bb*t38479;
                double t38476 = d0+t38480;
                double t38481 = 1.0/c;
                double t38483 = t38475*t38475;
                double t38484 = t38476*t38476;
                double t38485 = gamma_bb*t38474*t38483*t38484*9.0;
                double t38486 = t38485+1.0;
                double t38487 = 1.0/sqrt(t38486);
                double t38488 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38489 = d1*t38479;
                double t38490 = 1.0/(t38478*t38478);
                double t38498 = d1*gamma_bb*t38490;
                double t38491 = t38489-t38498;
                double t38492 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38493 = gamma_bb*gamma_bb;
                double t38494 = t38493*t38488;
                double t38495 = t38494+1.0;
                double t38496 = 1.0/sqrt(t38495);
                double t38497 = 1.0/pow(t38486,3.0/2.0);
                double t38499 = 1.0/pow(rho_b,1.9E1/3.0);
                double t38500 = t38474*t38483*t38484*9.0;
                double t38501 = gamma_bb*t38475*t38484*t38496*t38488*1.8E1;
                double t38502 = gamma_bb*t38491*t38474*t38483*t38476*1.8E1;
                double t38503 = t38500+t38501+t38502;
                double t38504 = gamma_bb*t38483*t38492*t38484*2.4E1;
                double t38505 = t38475*t38484*t38493*t38496*t38499*4.8E1;
                double t38506 = t38504+t38505;
                v_rho_b_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t38481*t38492*t38476*t38487*(-8.0/3.0)-gamma_bb*t38481*t38491*t38492*t38487*(8.0/3.0)+t38506*t38481*t38474*t38476*t38497*(1.0/2.0)+gamma_bb*t38503*t38481*t38492*t38476*t38497*(4.0/3.0)+gamma_bb*t38506*t38481*t38491*t38474*t38497*(1.0/2.0)+gamma_bb*t38481*t38474*t38476*t38497*(t38483*t38492*t38484*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t38484*t38493*4.8E1)/t38495+gamma_bb*t38491*t38483*t38492*t38476*4.8E1+gamma_bb*t38475*t38484*t38496*t38499*1.44E2+t38491*t38475*t38493*t38476*t38496*t38499*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t38475*t38484*t38493*1.0/pow(t38495,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t38503*t38506*t38481*t38474*t38476*1.0/pow(t38486,5.0/2.0)*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t38481*t38474*t38476*t38487+gamma_bb*t38481*t38491*t38474*t38487-gamma_bb*t38503*t38481*t38474*t38476*t38497*(1.0/2.0))*(4.0/3.0);
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t38558 = 1.0/pow(rho_a,8.0/3.0);
                double t38559 = log(gamma_aa*t38558+sqrt((gamma_aa*gamma_aa)*(t38558*t38558)+1.0));
                double t38561 = d2*d2;
                double t38562 = gamma_aa+t38561;
                double t38563 = 1.0/t38562;
                double t38564 = d1*gamma_aa*t38563;
                double t38560 = d0+t38564;
                v[Q] += scale * -c*pow(rho_a,4.0/3.0)*((gamma_aa*t38560*t38558*1.0/sqrt(gamma_aa*(t38560*t38560)*t38558*(t38559*t38559)*9.0+1.0))/c-1.0);
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t38566 = 1.0/pow(rho_a,8.0/3.0);
                double t38574 = gamma_aa*t38566;
                double t38567 = log(t38574+sqrt(t38574*t38574+1.0));
                double t38569 = d2*d2;
                double t38570 = gamma_aa+t38569;
                double t38571 = 1.0/t38570;
                double t38572 = d1*gamma_aa*t38571;
                double t38568 = d0+t38572;
                double t38573 = 1.0/c;
                double t38575 = t38567*t38567;
                double t38576 = t38568*t38568;
                double t38577 = gamma_aa*t38566*t38575*t38576*9.0;
                double t38578 = t38577+1.0;
                double t38579 = 1.0/sqrt(t38578);
                double t38580 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38581 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * c*pow(rho_a,1.0/3.0)*(gamma_aa*t38573*t38566*t38568*t38579-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t38580*t38573*t38568*t38579*(8.0/3.0)-gamma_aa*t38573*t38566*t38568*1.0/pow(t38578,3.0/2.0)*(gamma_aa*t38580*t38575*t38576*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t38581*t38567*t38576*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t38581+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t38584 = 1.0/pow(rho_a,8.0/3.0);
                double t38592 = gamma_aa*t38584;
                double t38585 = log(t38592+sqrt(t38592*t38592+1.0));
                double t38587 = d2*d2;
                double t38588 = gamma_aa+t38587;
                double t38589 = 1.0/t38588;
                double t38590 = d1*gamma_aa*t38589;
                double t38586 = d0+t38590;
                double t38591 = 1.0/c;
                double t38593 = t38585*t38585;
                double t38594 = t38586*t38586;
                double t38595 = gamma_aa*t38584*t38593*t38594*9.0;
                double t38596 = t38595+1.0;
                double t38597 = 1.0/sqrt(t38596);
                double t38598 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38599 = d1*t38589;
                double t38600 = 1.0/(t38588*t38588);
                double t38601 = t38599-d1*gamma_aa*t38600;
                v_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t38591*t38584*t38586*t38597+gamma_aa*t38601*t38591*t38584*t38597-gamma_aa*t38591*t38584*t38586*1.0/pow(t38596,3.0/2.0)*(t38584*t38593*t38594*9.0+gamma_aa*t38585*t38594*t38598*1.0/sqrt((gamma_aa*gamma_aa)*t38598+1.0)*1.8E1+gamma_aa*t38601*t38584*t38593*t38586*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t38607 = 1.0/pow(rho_a,8.0/3.0);
                double t38615 = gamma_aa*t38607;
                double t38608 = log(t38615+sqrt(t38615*t38615+1.0));
                double t38610 = d2*d2;
                double t38611 = gamma_aa+t38610;
                double t38612 = 1.0/t38611;
                double t38613 = d1*gamma_aa*t38612;
                double t38609 = d0+t38613;
                double t38614 = 1.0/c;
                double t38616 = t38608*t38608;
                double t38617 = t38609*t38609;
                double t38618 = gamma_aa*t38607*t38616*t38617*9.0;
                double t38619 = t38618+1.0;
                double t38620 = 1.0/sqrt(t38619);
                double t38621 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38622 = gamma_aa*gamma_aa;
                double t38623 = 1.0/pow(t38619,3.0/2.0);
                double t38624 = 1.0/pow(rho_a,1.4E1/3.0);
                double t38625 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38626 = t38622*t38625;
                double t38627 = t38626+1.0;
                double t38628 = 1.0/sqrt(t38627);
                double t38629 = gamma_aa*t38621*t38616*t38617*2.4E1;
                double t38630 = 1.0/pow(rho_a,1.9E1/3.0);
                double t38631 = t38630*t38622*t38608*t38617*t38628*4.8E1;
                double t38632 = t38631+t38629;
                v_rho_a_rho_a[Q] += scale * c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t38620*t38614*t38607*t38609-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t38620*t38614*t38624*t38609*(8.8E1/9.0)-gamma_aa*t38621*t38614*t38623*t38632*t38609*(8.0/3.0)-gamma_aa*t38614*t38623*t38607*t38609*(gamma_aa*t38624*t38616*t38617*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t38622*t38608*t38617*t38628*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t38622*t38617*1.28E2)/t38627-1.0/pow(rho_a,3.8E1/3.0)*(t38622*t38622)*t38608*t38617*1.0/pow(t38627,3.0/2.0)*1.28E2)*(1.0/2.0)+gamma_aa*t38614*(t38632*t38632)*t38607*t38609*1.0/pow(t38619,5.0/2.0)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t38620*t38621*t38614*t38609*(8.0/3.0)-gamma_aa*t38614*t38623*t38632*t38607*t38609*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t38675 = d2*d2;
                double t38676 = gamma_aa+t38675;
                double t38677 = 1.0/pow(rho_a,8.0/3.0);
                double t38682 = gamma_aa*t38677;
                double t38678 = log(t38682+sqrt(t38682*t38682+1.0));
                double t38679 = 1.0/t38676;
                double t38684 = d1*gamma_aa*t38679;
                double t38680 = d0+t38684;
                double t38681 = 1.0/c;
                double t38683 = t38678*t38678;
                double t38685 = t38680*t38680;
                double t38686 = gamma_aa*t38683*t38685*t38677*9.0;
                double t38687 = t38686+1.0;
                double t38688 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38689 = d1*t38679;
                double t38690 = 1.0/(t38676*t38676);
                double t38699 = d1*gamma_aa*t38690;
                double t38691 = t38689-t38699;
                double t38692 = 1.0/sqrt(t38687);
                double t38693 = t38683*t38685*t38677*9.0;
                double t38694 = gamma_aa*gamma_aa;
                double t38695 = t38694*t38688;
                double t38696 = t38695+1.0;
                double t38697 = 1.0/sqrt(t38696);
                double t38698 = gamma_aa*t38685*t38678*t38688*t38697*1.8E1;
                double t38700 = gamma_aa*t38680*t38691*t38683*t38677*1.8E1;
                double t38701 = t38700+t38693+t38698;
                double t38702 = 1.0/pow(t38687,3.0/2.0);
                double t38703 = d1*t38690*2.0;
                double t38704 = 1.0/(t38676*t38676*t38676);
                double t38705 = t38703-d1*gamma_aa*t38704*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * c*pow(rho_a,4.0/3.0)*(t38681*t38691*t38692*t38677*-2.0+gamma_aa*t38705*t38681*t38692*t38677+t38701*t38702*t38680*t38681*t38677+gamma_aa*t38701*t38702*t38681*t38691*t38677+gamma_aa*t38702*t38680*t38681*t38677*(t38680*t38691*t38683*t38677*3.6E1+t38685*t38678*t38688*t38697*3.6E1+gamma_aa*(t38691*t38691)*t38683*t38677*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t38685*1.8E1)/t38696-1.0/pow(rho_a,3.2E1/3.0)*t38685*t38694*t38678*1.0/pow(t38696,3.0/2.0)*1.8E1-gamma_aa*t38705*t38680*t38683*t38677*1.8E1+gamma_aa*t38680*t38691*t38678*t38688*t38697*7.2E1)*(1.0/2.0)-gamma_aa*(t38701*t38701)*t38680*t38681*t38677*1.0/pow(t38687,5.0/2.0)*(3.0/4.0));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t38636 = 1.0/pow(rho_a,8.0/3.0);
                double t38644 = gamma_aa*t38636;
                double t38637 = log(t38644+sqrt(t38644*t38644+1.0));
                double t38639 = d2*d2;
                double t38640 = gamma_aa+t38639;
                double t38641 = 1.0/t38640;
                double t38642 = d1*gamma_aa*t38641;
                double t38638 = d0+t38642;
                double t38643 = 1.0/c;
                double t38645 = t38637*t38637;
                double t38646 = t38638*t38638;
                double t38647 = gamma_aa*t38636*t38645*t38646*9.0;
                double t38648 = t38647+1.0;
                double t38649 = 1.0/sqrt(t38648);
                double t38650 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38651 = d1*t38641;
                double t38652 = 1.0/(t38640*t38640);
                double t38660 = d1*gamma_aa*t38652;
                double t38653 = t38651-t38660;
                double t38654 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38655 = gamma_aa*gamma_aa;
                double t38656 = t38650*t38655;
                double t38657 = t38656+1.0;
                double t38658 = 1.0/sqrt(t38657);
                double t38659 = 1.0/pow(t38648,3.0/2.0);
                double t38661 = 1.0/pow(rho_a,1.9E1/3.0);
                double t38662 = t38636*t38645*t38646*9.0;
                double t38663 = gamma_aa*t38650*t38637*t38646*t38658*1.8E1;
                double t38664 = gamma_aa*t38653*t38636*t38645*t38638*1.8E1;
                double t38665 = t38662+t38663+t38664;
                double t38666 = gamma_aa*t38645*t38654*t38646*2.4E1;
                double t38667 = t38661*t38637*t38646*t38655*t38658*4.8E1;
                double t38668 = t38666+t38667;
                v_rho_a_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t38643*t38654*t38638*t38649*(-8.0/3.0)-gamma_aa*t38643*t38653*t38654*t38649*(8.0/3.0)+t38643*t38636*t38638*t38659*t38668*(1.0/2.0)+gamma_aa*t38643*t38653*t38636*t38659*t38668*(1.0/2.0)+gamma_aa*t38643*t38654*t38638*t38665*t38659*(4.0/3.0)+gamma_aa*t38643*t38636*t38638*t38659*(t38645*t38654*t38646*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t38646*t38655*4.8E1)/t38657+gamma_aa*t38653*t38645*t38654*t38638*4.8E1+gamma_aa*t38661*t38637*t38646*t38658*1.44E2+t38661*t38653*t38637*t38655*t38638*t38658*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t38637*t38646*t38655*1.0/pow(t38657,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_aa*t38643*t38636*t38638*t38665*1.0/pow(t38648,5.0/2.0)*t38668*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t38643*t38636*t38638*t38649+gamma_aa*t38643*t38653*t38636*t38649-gamma_aa*t38643*t38636*t38638*t38665*t38659*(1.0/2.0))*(4.0/3.0);
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t38093 = 1.0/pow(rho_a,8.0/3.0);
                double t38094 = log(gamma_aa*t38093+sqrt((gamma_aa*gamma_aa)*(t38093*t38093)+1.0));
                double t38096 = d2*d2;
                double t38097 = gamma_aa+t38096;
                double t38098 = 1.0/t38097;
                double t38099 = d1*gamma_aa*t38098;
                double t38095 = d0+t38099;
                double t38100 = 1.0/c;
                double t38101 = 1.0/pow(rho_b,8.0/3.0);
                double t38102 = log(gamma_bb*t38101+sqrt((gamma_bb*gamma_bb)*(t38101*t38101)+1.0));
                double t38104 = gamma_bb+t38096;
                double t38105 = 1.0/t38104;
                double t38106 = d1*gamma_bb*t38105;
                double t38103 = d0+t38106;
                v[Q] += scale * -c*pow(rho_b,4.0/3.0)*(gamma_bb*t38100*t38101*t38103*1.0/sqrt(gamma_bb*t38101*(t38102*t38102)*(t38103*t38103)*9.0+1.0)-1.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t38100*t38093*t38095*1.0/sqrt(gamma_aa*t38093*(t38094*t38094)*(t38095*t38095)*9.0+1.0)-1.0);
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t38108 = 1.0/pow(rho_a,8.0/3.0);
                double t38116 = gamma_aa*t38108;
                double t38109 = log(t38116+sqrt(t38116*t38116+1.0));
                double t38111 = d2*d2;
                double t38112 = gamma_aa+t38111;
                double t38113 = 1.0/t38112;
                double t38114 = d1*gamma_aa*t38113;
                double t38110 = d0+t38114;
                double t38115 = 1.0/c;
                double t38117 = t38109*t38109;
                double t38118 = t38110*t38110;
                double t38119 = gamma_aa*t38108*t38117*t38118*9.0;
                double t38120 = t38119+1.0;
                double t38121 = 1.0/sqrt(t38120);
                double t38122 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38123 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * c*pow(rho_a,1.0/3.0)*(gamma_aa*t38110*t38121*t38115*t38108-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t38110*t38121*t38122*t38115*(8.0/3.0)-gamma_aa*t38110*1.0/pow(t38120,3.0/2.0)*t38115*t38108*(gamma_aa*t38122*t38117*t38118*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t38123*t38109*t38118*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t38123+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t38125 = 1.0/pow(rho_b,8.0/3.0);
                double t38133 = gamma_bb*t38125;
                double t38126 = log(t38133+sqrt(t38133*t38133+1.0));
                double t38128 = d2*d2;
                double t38129 = gamma_bb+t38128;
                double t38130 = 1.0/t38129;
                double t38131 = d1*gamma_bb*t38130;
                double t38127 = d0+t38131;
                double t38132 = 1.0/c;
                double t38134 = t38126*t38126;
                double t38135 = t38127*t38127;
                double t38136 = gamma_bb*t38125*t38134*t38135*9.0;
                double t38137 = t38136+1.0;
                double t38138 = 1.0/sqrt(t38137);
                double t38139 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38140 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * c*pow(rho_b,1.0/3.0)*(gamma_bb*t38132*t38125*t38127*t38138-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t38132*t38127*t38138*t38139*(8.0/3.0)-gamma_bb*t38132*t38125*t38127*1.0/pow(t38137,3.0/2.0)*(gamma_bb*t38134*t38135*t38139*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t38140*t38126*t38135*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t38140+1.0)*4.8E1)*(1.0/2.0));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t38142 = 1.0/pow(rho_a,8.0/3.0);
                double t38150 = gamma_aa*t38142;
                double t38143 = log(t38150+sqrt(t38150*t38150+1.0));
                double t38145 = d2*d2;
                double t38146 = gamma_aa+t38145;
                double t38147 = 1.0/t38146;
                double t38148 = d1*gamma_aa*t38147;
                double t38144 = d0+t38148;
                double t38149 = 1.0/c;
                double t38151 = t38143*t38143;
                double t38152 = t38144*t38144;
                double t38153 = gamma_aa*t38142*t38151*t38152*9.0;
                double t38154 = t38153+1.0;
                double t38155 = 1.0/sqrt(t38154);
                double t38156 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38157 = d1*t38147;
                double t38158 = 1.0/(t38146*t38146);
                double t38159 = t38157-d1*gamma_aa*t38158;
                v_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t38142*t38144*t38155*t38149+gamma_aa*t38142*t38155*t38149*t38159-gamma_aa*t38142*t38144*1.0/pow(t38154,3.0/2.0)*t38149*(t38142*t38151*t38152*9.0+gamma_aa*t38143*t38152*t38156*1.0/sqrt((gamma_aa*gamma_aa)*t38156+1.0)*1.8E1+gamma_aa*t38142*t38151*t38144*t38159*1.8E1)*(1.0/2.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t38162 = 1.0/pow(rho_b,8.0/3.0);
                double t38170 = gamma_bb*t38162;
                double t38163 = log(t38170+sqrt(t38170*t38170+1.0));
                double t38165 = d2*d2;
                double t38166 = gamma_bb+t38165;
                double t38167 = 1.0/t38166;
                double t38168 = d1*gamma_bb*t38167;
                double t38164 = d0+t38168;
                double t38169 = 1.0/c;
                double t38171 = t38163*t38163;
                double t38172 = t38164*t38164;
                double t38173 = gamma_bb*t38162*t38171*t38172*9.0;
                double t38174 = t38173+1.0;
                double t38175 = 1.0/sqrt(t38174);
                double t38176 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38177 = d1*t38167;
                double t38178 = 1.0/(t38166*t38166);
                double t38179 = t38177-d1*gamma_bb*t38178;
                v_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t38162*t38164*t38175*t38169+gamma_bb*t38162*t38175*t38169*t38179-gamma_bb*t38162*t38164*1.0/pow(t38174,3.0/2.0)*t38169*(t38162*t38171*t38172*9.0+gamma_bb*t38163*t38172*t38176*1.0/sqrt((gamma_bb*gamma_bb)*t38176+1.0)*1.8E1+gamma_bb*t38162*t38171*t38164*t38179*1.8E1)*(1.0/2.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t38183 = 1.0/pow(rho_a,8.0/3.0);
                double t38191 = gamma_aa*t38183;
                double t38184 = log(t38191+sqrt(t38191*t38191+1.0));
                double t38186 = d2*d2;
                double t38187 = gamma_aa+t38186;
                double t38188 = 1.0/t38187;
                double t38189 = d1*gamma_aa*t38188;
                double t38185 = d0+t38189;
                double t38190 = 1.0/c;
                double t38192 = t38184*t38184;
                double t38193 = t38185*t38185;
                double t38194 = gamma_aa*t38183*t38192*t38193*9.0;
                double t38195 = t38194+1.0;
                double t38196 = 1.0/sqrt(t38195);
                double t38197 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38198 = gamma_aa*gamma_aa;
                double t38199 = 1.0/pow(t38195,3.0/2.0);
                double t38200 = 1.0/pow(rho_a,1.4E1/3.0);
                double t38201 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38202 = t38201*t38198;
                double t38203 = t38202+1.0;
                double t38204 = 1.0/sqrt(t38203);
                double t38205 = gamma_aa*t38192*t38193*t38197*2.4E1;
                double t38206 = 1.0/pow(rho_a,1.9E1/3.0);
                double t38207 = t38204*t38206*t38184*t38193*t38198*4.8E1;
                double t38208 = t38205+t38207;
                v_rho_a_rho_a[Q] += scale * c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t38190*t38183*t38185*t38196-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t38200*t38190*t38185*t38196*(8.8E1/9.0)-gamma_aa*t38190*t38208*t38185*t38197*t38199*(8.0/3.0)-gamma_aa*t38190*t38183*t38185*t38199*(gamma_aa*t38200*t38192*t38193*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t38204*t38184*t38193*t38198*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t38193*t38198*1.28E2)/t38203-1.0/pow(rho_a,3.8E1/3.0)*1.0/pow(t38203,3.0/2.0)*t38184*t38193*(t38198*t38198)*1.28E2)*(1.0/2.0)+gamma_aa*t38190*(t38208*t38208)*t38183*t38185*1.0/pow(t38195,5.0/2.0)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t38190*t38185*t38196*t38197*(8.0/3.0)-gamma_aa*t38190*t38208*t38183*t38185*t38199*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t38211 = 1.0/pow(rho_b,8.0/3.0);
                double t38219 = gamma_bb*t38211;
                double t38212 = log(t38219+sqrt(t38219*t38219+1.0));
                double t38214 = d2*d2;
                double t38215 = gamma_bb+t38214;
                double t38216 = 1.0/t38215;
                double t38217 = d1*gamma_bb*t38216;
                double t38213 = d0+t38217;
                double t38218 = 1.0/c;
                double t38220 = t38212*t38212;
                double t38221 = t38213*t38213;
                double t38222 = gamma_bb*t38211*t38220*t38221*9.0;
                double t38223 = t38222+1.0;
                double t38224 = 1.0/sqrt(t38223);
                double t38225 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38226 = gamma_bb*gamma_bb;
                double t38227 = 1.0/pow(t38223,3.0/2.0);
                double t38228 = 1.0/pow(rho_b,1.4E1/3.0);
                double t38229 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38230 = t38226*t38229;
                double t38231 = t38230+1.0;
                double t38232 = 1.0/sqrt(t38231);
                double t38233 = gamma_bb*t38220*t38221*t38225*2.4E1;
                double t38234 = 1.0/pow(rho_b,1.9E1/3.0);
                double t38235 = t38212*t38221*t38232*t38234*t38226*4.8E1;
                double t38236 = t38233+t38235;
                v_rho_b_rho_b[Q] += scale * c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t38211*t38213*t38224*t38218-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t38213*t38224*t38218*t38228*(8.8E1/9.0)-gamma_bb*t38213*t38225*t38218*t38227*t38236*(8.0/3.0)-gamma_bb*t38211*t38213*t38218*t38227*(gamma_bb*t38220*t38221*t38228*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t38212*t38221*t38232*t38226*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t38221*t38226*1.28E2)/t38231-1.0/pow(rho_b,3.8E1/3.0)*t38212*t38221*1.0/pow(t38231,3.0/2.0)*(t38226*t38226)*1.28E2)*(1.0/2.0)+gamma_bb*t38211*t38213*1.0/pow(t38223,5.0/2.0)*t38218*(t38236*t38236)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t38213*t38224*t38225*t38218*(8.0/3.0)-gamma_bb*t38211*t38213*t38218*t38227*t38236*(1.0/2.0))*(8.0/3.0);
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t38310 = d2*d2;
                double t38311 = gamma_aa+t38310;
                double t38312 = 1.0/pow(rho_a,8.0/3.0);
                double t38317 = gamma_aa*t38312;
                double t38313 = log(t38317+sqrt(t38317*t38317+1.0));
                double t38314 = 1.0/t38311;
                double t38319 = d1*gamma_aa*t38314;
                double t38315 = d0+t38319;
                double t38316 = 1.0/c;
                double t38318 = t38313*t38313;
                double t38320 = t38315*t38315;
                double t38321 = gamma_aa*t38320*t38312*t38318*9.0;
                double t38322 = t38321+1.0;
                double t38323 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38324 = d1*t38314;
                double t38325 = 1.0/(t38311*t38311);
                double t38334 = d1*gamma_aa*t38325;
                double t38326 = t38324-t38334;
                double t38327 = 1.0/sqrt(t38322);
                double t38328 = t38320*t38312*t38318*9.0;
                double t38329 = gamma_aa*gamma_aa;
                double t38330 = t38323*t38329;
                double t38331 = t38330+1.0;
                double t38332 = 1.0/sqrt(t38331);
                double t38333 = gamma_aa*t38320*t38313*t38323*t38332*1.8E1;
                double t38335 = gamma_aa*t38312*t38315*t38326*t38318*1.8E1;
                double t38336 = t38333+t38335+t38328;
                double t38337 = 1.0/pow(t38322,3.0/2.0);
                double t38338 = d1*t38325*2.0;
                double t38339 = 1.0/(t38311*t38311*t38311);
                double t38340 = t38338-d1*gamma_aa*t38339*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * c*pow(rho_a,4.0/3.0)*(t38312*t38316*t38326*t38327*-2.0+gamma_aa*t38312*t38340*t38316*t38327+t38312*t38315*t38316*t38336*t38337+gamma_aa*t38312*t38316*t38326*t38336*t38337+gamma_aa*t38312*t38315*t38316*t38337*(t38320*t38313*t38323*t38332*3.6E1+t38312*t38315*t38326*t38318*3.6E1+gamma_aa*t38312*(t38326*t38326)*t38318*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t38320*1.8E1)/t38331-1.0/pow(rho_a,3.2E1/3.0)*t38320*t38313*1.0/pow(t38331,3.0/2.0)*t38329*1.8E1-gamma_aa*t38312*t38340*t38315*t38318*1.8E1+gamma_aa*t38313*t38323*t38332*t38315*t38326*7.2E1)*(1.0/2.0)-gamma_aa*t38312*1.0/pow(t38322,5.0/2.0)*t38315*t38316*(t38336*t38336)*(3.0/4.0));
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t38346 = d2*d2;
                double t38347 = gamma_bb+t38346;
                double t38348 = 1.0/pow(rho_b,8.0/3.0);
                double t38353 = gamma_bb*t38348;
                double t38349 = log(t38353+sqrt(t38353*t38353+1.0));
                double t38350 = 1.0/t38347;
                double t38355 = d1*gamma_bb*t38350;
                double t38351 = d0+t38355;
                double t38352 = 1.0/c;
                double t38354 = t38349*t38349;
                double t38356 = t38351*t38351;
                double t38357 = gamma_bb*t38354*t38356*t38348*9.0;
                double t38358 = t38357+1.0;
                double t38359 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38360 = d1*t38350;
                double t38361 = 1.0/(t38347*t38347);
                double t38370 = d1*gamma_bb*t38361;
                double t38362 = t38360-t38370;
                double t38363 = 1.0/sqrt(t38358);
                double t38364 = t38354*t38356*t38348*9.0;
                double t38365 = gamma_bb*gamma_bb;
                double t38366 = t38365*t38359;
                double t38367 = t38366+1.0;
                double t38368 = 1.0/sqrt(t38367);
                double t38369 = gamma_bb*t38356*t38349*t38359*t38368*1.8E1;
                double t38371 = gamma_bb*t38351*t38362*t38354*t38348*1.8E1;
                double t38372 = t38371+t38364+t38369;
                double t38373 = 1.0/pow(t38358,3.0/2.0);
                double t38374 = d1*t38361*2.0;
                double t38375 = 1.0/(t38347*t38347*t38347);
                double t38376 = t38374-d1*gamma_bb*t38375*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * c*pow(rho_b,4.0/3.0)*(t38352*t38362*t38363*t38348*-2.0+gamma_bb*t38352*t38363*t38348*t38376+t38351*t38352*t38372*t38373*t38348+gamma_bb*t38352*t38362*t38372*t38373*t38348+gamma_bb*t38351*t38352*t38373*t38348*(t38351*t38362*t38354*t38348*3.6E1+t38356*t38349*t38359*t38368*3.6E1+gamma_bb*(t38362*t38362)*t38354*t38348*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t38356*1.8E1)/t38367-1.0/pow(rho_b,3.2E1/3.0)*t38356*t38365*t38349*1.0/pow(t38367,3.0/2.0)*1.8E1-gamma_bb*t38351*t38354*t38348*t38376*1.8E1+gamma_bb*t38351*t38362*t38349*t38359*t38368*7.2E1)*(1.0/2.0)-gamma_bb*t38351*t38352*(t38372*t38372)*t38348*1.0/pow(t38358,5.0/2.0)*(3.0/4.0));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t38238 = 1.0/pow(rho_a,8.0/3.0);
                double t38246 = gamma_aa*t38238;
                double t38239 = log(t38246+sqrt(t38246*t38246+1.0));
                double t38241 = d2*d2;
                double t38242 = gamma_aa+t38241;
                double t38243 = 1.0/t38242;
                double t38244 = d1*gamma_aa*t38243;
                double t38240 = d0+t38244;
                double t38245 = 1.0/c;
                double t38247 = t38239*t38239;
                double t38248 = t38240*t38240;
                double t38249 = gamma_aa*t38238*t38247*t38248*9.0;
                double t38250 = t38249+1.0;
                double t38251 = 1.0/sqrt(t38250);
                double t38252 = 1.0/pow(rho_a,1.6E1/3.0);
                double t38253 = d1*t38243;
                double t38254 = 1.0/(t38242*t38242);
                double t38262 = d1*gamma_aa*t38254;
                double t38255 = t38253-t38262;
                double t38256 = 1.0/pow(rho_a,1.1E1/3.0);
                double t38257 = gamma_aa*gamma_aa;
                double t38258 = t38252*t38257;
                double t38259 = t38258+1.0;
                double t38260 = 1.0/sqrt(t38259);
                double t38261 = 1.0/pow(t38250,3.0/2.0);
                double t38263 = 1.0/pow(rho_a,1.9E1/3.0);
                double t38264 = t38238*t38247*t38248*9.0;
                double t38265 = gamma_aa*t38260*t38252*t38239*t38248*1.8E1;
                double t38266 = gamma_aa*t38240*t38255*t38238*t38247*1.8E1;
                double t38267 = t38264+t38265+t38266;
                double t38268 = gamma_aa*t38247*t38256*t38248*2.4E1;
                double t38269 = t38260*t38263*t38239*t38248*t38257*4.8E1;
                double t38270 = t38268+t38269;
                v_rho_a_gamma_aa[Q] += scale * -c*pow(rho_a,4.0/3.0)*(t38240*t38251*t38245*t38256*(-8.0/3.0)-gamma_aa*t38251*t38245*t38255*t38256*(8.0/3.0)+t38240*t38261*t38270*t38245*t38238*(1.0/2.0)+gamma_aa*t38240*t38261*t38245*t38256*t38267*(4.0/3.0)+gamma_aa*t38261*t38270*t38245*t38255*t38238*(1.0/2.0)+gamma_aa*t38240*t38261*t38245*t38238*(t38247*t38256*t38248*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t38248*t38257*4.8E1)/t38259+gamma_aa*t38240*t38255*t38247*t38256*4.8E1+gamma_aa*t38260*t38263*t38239*t38248*1.44E2+t38240*t38260*t38263*t38255*t38239*t38257*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t38239*t38248*t38257*1.0/pow(t38259,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_aa*t38240*1.0/pow(t38250,5.0/2.0)*t38270*t38245*t38238*t38267*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t38240*t38251*t38245*t38238+gamma_aa*t38251*t38245*t38255*t38238-gamma_aa*t38240*t38261*t38245*t38238*t38267*(1.0/2.0))*(4.0/3.0);
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t38276 = 1.0/pow(rho_b,8.0/3.0);
                double t38284 = gamma_bb*t38276;
                double t38277 = log(t38284+sqrt(t38284*t38284+1.0));
                double t38279 = d2*d2;
                double t38280 = gamma_bb+t38279;
                double t38281 = 1.0/t38280;
                double t38282 = d1*gamma_bb*t38281;
                double t38278 = d0+t38282;
                double t38283 = 1.0/c;
                double t38285 = t38277*t38277;
                double t38286 = t38278*t38278;
                double t38287 = gamma_bb*t38276*t38285*t38286*9.0;
                double t38288 = t38287+1.0;
                double t38289 = 1.0/sqrt(t38288);
                double t38290 = 1.0/pow(rho_b,1.6E1/3.0);
                double t38291 = d1*t38281;
                double t38292 = 1.0/(t38280*t38280);
                double t38300 = d1*gamma_bb*t38292;
                double t38293 = -t38300+t38291;
                double t38294 = 1.0/pow(rho_b,1.1E1/3.0);
                double t38295 = gamma_bb*gamma_bb;
                double t38296 = t38290*t38295;
                double t38297 = t38296+1.0;
                double t38298 = 1.0/sqrt(t38297);
                double t38299 = 1.0/pow(t38288,3.0/2.0);
                double t38301 = 1.0/pow(rho_b,1.9E1/3.0);
                double t38302 = t38276*t38285*t38286*9.0;
                double t38303 = gamma_bb*t38290*t38277*t38286*t38298*1.8E1;
                double t38304 = gamma_bb*t38293*t38276*t38285*t38278*1.8E1;
                double t38305 = t38302+t38303+t38304;
                double t38306 = gamma_bb*t38285*t38294*t38286*2.4E1;
                double t38307 = t38301*t38277*t38286*t38295*t38298*4.8E1;
                double t38308 = t38306+t38307;
                v_rho_b_gamma_bb[Q] += scale * -c*pow(rho_b,4.0/3.0)*(t38283*t38294*t38278*t38289*(-8.0/3.0)-gamma_bb*t38283*t38293*t38294*t38289*(8.0/3.0)+t38308*t38283*t38276*t38278*t38299*(1.0/2.0)+gamma_bb*t38305*t38283*t38294*t38278*t38299*(4.0/3.0)+gamma_bb*t38308*t38283*t38293*t38276*t38299*(1.0/2.0)+gamma_bb*t38283*t38276*t38278*t38299*(t38285*t38294*t38286*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t38286*t38295*4.8E1)/t38297+gamma_bb*t38301*t38277*t38286*t38298*1.44E2+gamma_bb*t38293*t38285*t38294*t38278*4.8E1+t38301*t38293*t38277*t38295*t38278*t38298*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t38277*t38286*t38295*1.0/pow(t38297,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t38305*t38308*t38283*t38276*t38278*1.0/pow(t38288,5.0/2.0)*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t38283*t38276*t38278*t38289+gamma_bb*t38283*t38293*t38276*t38289-gamma_bb*t38305*t38283*t38276*t38278*t38299*(1.0/2.0))*(4.0/3.0);
            }
            
        }
    }
}

}
