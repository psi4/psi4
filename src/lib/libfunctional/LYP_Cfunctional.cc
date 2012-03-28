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
            {
            double t9423 = rho_a+rho_b;
            double t9424 = 1.0/pow(t9423,1.0/3.0);
            double t9425 = Dd*t9424;
            double t9426 = t9425+1.0;
            double t9427 = 1.0/t9426;
            double t9428 = t9423*t9423;
            double t9429 = t9428*(2.0/3.0);
            double t9430 = gamma_ab*2.0;
            double t9431 = gamma_aa+gamma_bb+t9430;
            double t9432 = 1.0/t9423;
            v[Q] += scale * A*rho_a*rho_b*t9432*t9427*-4.0-A*B*1.0/pow(t9423,1.1E1/3.0)*t9427*exp(-C*t9424)*(t9431*t9428*(-2.0/3.0)+gamma_aa*(t9429-rho_b*rho_b)+gamma_bb*(t9429-rho_a*rho_a)+rho_a*rho_b*((gamma_aa+gamma_bb)*(C*t9424*(1.0/1.8E1)+Dd*t9424*t9427*(1.0/1.8E1)-5.0/2.0)+CFext*(pow(rho_a,8.0/3.0)+pow(rho_b,8.0/3.0))-t9431*(C*t9424*(7.0/1.8E1)+Dd*t9424*t9427*(7.0/1.8E1)-4.7E1/1.8E1)-t9432*(gamma_aa*rho_a+gamma_bb*rho_b)*(C*t9424*(1.0/9.0)+Dd*t9424*t9427*(1.0/9.0)-1.1E1/9.0)));
            }
            
            // v_rho_a
            {
            double t9434 = rho_a+rho_b;
            double t9435 = 1.0/pow(t9434,1.0/3.0);
            double t9436 = Dd*t9435;
            double t9437 = t9436+1.0;
            double t9438 = 1.0/t9437;
            double t9439 = t9434*t9434;
            double t9440 = t9439*(2.0/3.0);
            double t9441 = gamma_ab*2.0;
            double t9442 = gamma_aa+gamma_bb+t9441;
            double t9443 = 1.0/t9434;
            double t9470 = C*t9435;
            double t9444 = exp(-t9470);
            double t9445 = C*t9435*(7.0/1.8E1);
            double t9446 = Dd*t9435*t9438*(7.0/1.8E1);
            double t9447 = t9445+t9446-4.7E1/1.8E1;
            double t9448 = t9442*t9447;
            double t9449 = gamma_aa+gamma_bb;
            double t9450 = C*t9435*(1.0/1.8E1);
            double t9451 = Dd*t9435*t9438*(1.0/1.8E1);
            double t9452 = t9450+t9451-5.0/2.0;
            double t9453 = pow(rho_a,8.0/3.0);
            double t9454 = pow(rho_b,8.0/3.0);
            double t9455 = t9453+t9454;
            double t9456 = gamma_aa*rho_a;
            double t9457 = gamma_bb*rho_b;
            double t9458 = t9456+t9457;
            double t9459 = C*t9435*(1.0/9.0);
            double t9460 = Dd*t9435*t9438*(1.0/9.0);
            double t9461 = t9460+t9459-1.1E1/9.0;
            double t9462 = t9443*t9461*t9458;
            double t9477 = t9452*t9449;
            double t9478 = CFext*t9455;
            double t9463 = t9462+t9448-t9477-t9478;
            double t9464 = rho_b*(4.0/3.0);
            double t9465 = 1.0/pow(t9434,4.0/3.0);
            double t9466 = 1.0/(t9437*t9437);
            double t9467 = Dd*Dd;
            double t9468 = 1.0/pow(t9434,5.0/3.0);
            double t9469 = 1.0/(t9434*t9434);
            double t9471 = rho_b*rho_b;
            double t9472 = t9440-t9471;
            double t9473 = gamma_aa*t9472;
            double t9474 = rho_a*rho_a;
            double t9475 = t9440-t9474;
            double t9476 = gamma_bb*t9475;
            double t9481 = t9442*t9439*(2.0/3.0);
            double t9482 = rho_a*rho_b*t9463;
            double t9479 = -t9481+t9473-t9482+t9476;
            double t9480 = 1.0/(t9434*t9434*t9434*t9434*t9434);
            v_rho_a[Q] += scale * A*rho_b*t9443*t9438*-4.0+A*rho_a*rho_b*t9438*t9469*4.0+A*B*1.0/pow(t9434,1.1E1/3.0)*t9444*t9438*(rho_b*t9463-gamma_aa*(rho_a*(4.0/3.0)+t9464)+gamma_bb*(rho_a*(2.0/3.0)-t9464)+t9442*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_a,5.0/3.0)*(8.0/3.0)-t9449*(C*t9465*(1.0/5.4E1)+Dd*t9438*t9465*(1.0/5.4E1)-t9466*t9467*t9468*(1.0/5.4E1))+t9442*(C*t9465*(7.0/5.4E1)+Dd*t9438*t9465*(7.0/5.4E1)-t9466*t9467*t9468*(7.0/5.4E1))+t9443*t9458*(C*t9465*(1.0/2.7E1)+Dd*t9438*t9465*(1.0/2.7E1)-t9466*t9467*t9468*(1.0/2.7E1))-gamma_aa*t9443*t9461+t9461*t9458*t9469))-A*Dd*rho_a*rho_b*1.0/pow(t9434,7.0/3.0)*t9466*(4.0/3.0)+A*B*1.0/pow(t9434,1.4E1/3.0)*t9444*t9438*t9479*(1.1E1/3.0)-A*B*C*t9444*t9480*t9438*t9479*(1.0/3.0)-A*B*Dd*t9444*t9480*t9466*t9479*(1.0/3.0);
            }
            
            // v_rho_b
            {
            double t9484 = rho_a+rho_b;
            double t9485 = 1.0/pow(t9484,1.0/3.0);
            double t9486 = Dd*t9485;
            double t9487 = t9486+1.0;
            double t9488 = 1.0/t9487;
            double t9489 = t9484*t9484;
            double t9490 = t9489*(2.0/3.0);
            double t9491 = gamma_ab*2.0;
            double t9492 = gamma_aa+gamma_bb+t9491;
            double t9493 = 1.0/t9484;
            double t9520 = C*t9485;
            double t9494 = exp(-t9520);
            double t9495 = C*t9485*(7.0/1.8E1);
            double t9496 = Dd*t9485*t9488*(7.0/1.8E1);
            double t9497 = t9495+t9496-4.7E1/1.8E1;
            double t9498 = t9492*t9497;
            double t9499 = gamma_aa+gamma_bb;
            double t9500 = C*t9485*(1.0/1.8E1);
            double t9501 = Dd*t9485*t9488*(1.0/1.8E1);
            double t9502 = t9500+t9501-5.0/2.0;
            double t9503 = pow(rho_a,8.0/3.0);
            double t9504 = pow(rho_b,8.0/3.0);
            double t9505 = t9503+t9504;
            double t9506 = gamma_aa*rho_a;
            double t9507 = gamma_bb*rho_b;
            double t9508 = t9506+t9507;
            double t9509 = C*t9485*(1.0/9.0);
            double t9510 = Dd*t9485*t9488*(1.0/9.0);
            double t9511 = t9510+t9509-1.1E1/9.0;
            double t9512 = t9493*t9511*t9508;
            double t9527 = t9499*t9502;
            double t9528 = CFext*t9505;
            double t9513 = t9498+t9512-t9527-t9528;
            double t9514 = rho_a*(4.0/3.0);
            double t9515 = 1.0/pow(t9484,4.0/3.0);
            double t9516 = 1.0/(t9487*t9487);
            double t9517 = Dd*Dd;
            double t9518 = 1.0/pow(t9484,5.0/3.0);
            double t9519 = 1.0/(t9484*t9484);
            double t9521 = rho_b*rho_b;
            double t9522 = t9490-t9521;
            double t9523 = gamma_aa*t9522;
            double t9524 = rho_a*rho_a;
            double t9525 = t9490-t9524;
            double t9526 = gamma_bb*t9525;
            double t9531 = t9492*t9489*(2.0/3.0);
            double t9532 = rho_a*rho_b*t9513;
            double t9529 = -t9531+t9523-t9532+t9526;
            double t9530 = 1.0/(t9484*t9484*t9484*t9484*t9484);
            v_rho_b[Q] += scale * A*rho_a*t9493*t9488*-4.0+A*rho_a*rho_b*t9488*t9519*4.0+A*B*1.0/pow(t9484,1.1E1/3.0)*t9494*t9488*(rho_a*t9513-gamma_bb*(rho_b*(4.0/3.0)+t9514)+gamma_aa*(rho_b*(2.0/3.0)-t9514)+t9492*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_b,5.0/3.0)*(8.0/3.0)-t9499*(C*t9515*(1.0/5.4E1)+Dd*t9488*t9515*(1.0/5.4E1)-t9516*t9517*t9518*(1.0/5.4E1))+t9492*(C*t9515*(7.0/5.4E1)+Dd*t9488*t9515*(7.0/5.4E1)-t9516*t9517*t9518*(7.0/5.4E1))+t9493*t9508*(C*t9515*(1.0/2.7E1)+Dd*t9488*t9515*(1.0/2.7E1)-t9516*t9517*t9518*(1.0/2.7E1))-gamma_bb*t9493*t9511+t9511*t9508*t9519))-A*Dd*rho_a*rho_b*1.0/pow(t9484,7.0/3.0)*t9516*(4.0/3.0)+A*B*1.0/pow(t9484,1.4E1/3.0)*t9494*t9488*t9529*(1.1E1/3.0)-A*B*C*t9494*t9488*t9530*t9529*(1.0/3.0)-A*B*Dd*t9494*t9530*t9516*t9529*(1.0/3.0);
            }
            
            // v_gamma_aa
            {
            double t9534 = rho_a+rho_b;
            double t9535 = 1.0/pow(t9534,1.0/3.0);
            double t9536 = Dd*t9535;
            double t9537 = t9536+1.0;
            double t9538 = 1.0/t9537;
            v_gamma_aa[Q] += scale * A*B*1.0/pow(t9534,1.1E1/3.0)*t9538*exp(-C*t9535)*(rho_b*rho_b+rho_a*rho_b*(C*t9535*(1.0/3.0)+Dd*t9535*t9538*(1.0/3.0)+(rho_a*(C*t9535*(1.0/9.0)+Dd*t9535*t9538*(1.0/9.0)-1.1E1/9.0))/t9534-1.0/9.0));
            }
            
            // v_gamma_ab
            {
            double t9540 = rho_a+rho_b;
            double t9541 = 1.0/pow(t9540,1.0/3.0);
            double t9542 = Dd*t9541;
            double t9543 = t9542+1.0;
            double t9544 = 1.0/t9543;
            v_gamma_ab[Q] += scale * A*B*1.0/pow(t9540,1.1E1/3.0)*t9544*exp(-C*t9541)*((t9540*t9540)*(4.0/3.0)+rho_a*rho_b*(C*t9541*(7.0/9.0)+Dd*t9541*t9544*(7.0/9.0)-4.7E1/9.0));
            }
            
            // v_gamma_bb
            {
            double t9546 = rho_a+rho_b;
            double t9547 = 1.0/pow(t9546,1.0/3.0);
            double t9548 = Dd*t9547;
            double t9549 = t9548+1.0;
            double t9550 = 1.0/t9549;
            v_gamma_bb[Q] += scale * A*B*t9550*1.0/pow(t9546,1.1E1/3.0)*exp(-C*t9547)*(rho_a*rho_a+rho_a*rho_b*(C*t9547*(1.0/3.0)+Dd*t9550*t9547*(1.0/3.0)+(rho_b*(C*t9547*(1.0/9.0)+Dd*t9550*t9547*(1.0/9.0)-1.1E1/9.0))/t9546-1.0/9.0));
            }
            
            // v_rho_a_rho_a
            {
            double t9554 = rho_a+rho_b;
            double t9555 = 1.0/pow(t9554,1.0/3.0);
            double t9556 = Dd*t9555;
            double t9557 = t9556+1.0;
            double t9558 = 1.0/t9557;
            double t9559 = t9554*t9554;
            double t9560 = t9559*(2.0/3.0);
            double t9561 = gamma_ab*2.0;
            double t9562 = gamma_aa+gamma_bb+t9561;
            double t9563 = 1.0/(t9557*t9557);
            double t9590 = C*t9555;
            double t9564 = exp(-t9590);
            double t9565 = C*t9555*(7.0/1.8E1);
            double t9566 = Dd*t9555*t9558*(7.0/1.8E1);
            double t9567 = t9565+t9566-4.7E1/1.8E1;
            double t9568 = t9562*t9567;
            double t9569 = gamma_aa+gamma_bb;
            double t9570 = C*t9555*(1.0/1.8E1);
            double t9571 = Dd*t9555*t9558*(1.0/1.8E1);
            double t9572 = t9570+t9571-5.0/2.0;
            double t9573 = pow(rho_a,8.0/3.0);
            double t9574 = pow(rho_b,8.0/3.0);
            double t9575 = t9573+t9574;
            double t9576 = 1.0/t9554;
            double t9577 = gamma_aa*rho_a;
            double t9578 = gamma_bb*rho_b;
            double t9579 = t9577+t9578;
            double t9580 = C*t9555*(1.0/9.0);
            double t9581 = Dd*t9555*t9558*(1.0/9.0);
            double t9582 = t9580+t9581-1.1E1/9.0;
            double t9583 = t9582*t9576*t9579;
            double t9616 = t9572*t9569;
            double t9617 = CFext*t9575;
            double t9584 = t9583+t9568-t9616-t9617;
            double t9585 = rho_b*(4.0/3.0);
            double t9586 = 1.0/pow(t9554,4.0/3.0);
            double t9587 = Dd*Dd;
            double t9588 = 1.0/pow(t9554,5.0/3.0);
            double t9589 = 1.0/(t9554*t9554);
            double t9591 = 1.0/pow(t9554,1.1E1/3.0);
            double t9592 = C*t9586*(1.0/5.4E1);
            double t9593 = Dd*t9558*t9586*(1.0/5.4E1);
            double t9624 = t9563*t9587*t9588*(1.0/5.4E1);
            double t9594 = t9592+t9593-t9624;
            double t9595 = pow(rho_a,5.0/3.0);
            double t9596 = CFext*t9595*(8.0/3.0);
            double t9597 = C*t9586*(7.0/5.4E1);
            double t9598 = Dd*t9558*t9586*(7.0/5.4E1);
            double t9626 = t9563*t9587*t9588*(7.0/5.4E1);
            double t9599 = t9597+t9598-t9626;
            double t9600 = t9562*t9599;
            double t9601 = C*t9586*(1.0/2.7E1);
            double t9602 = Dd*t9558*t9586*(1.0/2.7E1);
            double t9611 = t9563*t9587*t9588*(1.0/2.7E1);
            double t9603 = t9601+t9602-t9611;
            double t9604 = t9576*t9579*t9603;
            double t9605 = t9582*t9579*t9589;
            double t9625 = t9594*t9569;
            double t9627 = gamma_aa*t9582*t9576;
            double t9606 = t9596+t9600+t9604+t9605-t9625-t9627;
            double t9607 = 1.0/pow(t9554,7.0/3.0);
            double t9608 = 1.0/(t9554*t9554*t9554);
            double t9609 = 1.0/(t9557*t9557*t9557);
            double t9610 = 1.0/pow(t9554,8.0/3.0);
            double t9612 = rho_a*2.0;
            double t9613 = rho_b*2.0;
            double t9614 = t9612+t9613;
            double t9615 = t9562*t9614*(2.0/3.0);
            double t9618 = rho_b*t9584;
            double t9619 = rho_a*(4.0/3.0);
            double t9620 = t9585+t9619;
            double t9621 = rho_a*(2.0/3.0);
            double t9622 = t9585-t9621;
            double t9623 = gamma_bb*t9622;
            double t9628 = rho_a*rho_b*t9606;
            double t9629 = gamma_aa*t9620;
            double t9630 = t9623-t9615-t9618+t9628+t9629;
            double t9631 = 1.0/(t9554*t9554*t9554*t9554*t9554);
            double t9632 = rho_b*rho_b;
            double t9633 = t9560-t9632;
            double t9634 = gamma_aa*t9633;
            double t9635 = rho_a*rho_a;
            double t9636 = t9560-t9635;
            double t9637 = gamma_bb*t9636;
            double t9640 = t9562*t9559*(2.0/3.0);
            double t9641 = rho_a*rho_b*t9584;
            double t9638 = -t9640-t9641+t9634+t9637;
            double t9639 = 1.0/(t9554*t9554*t9554*t9554*t9554*t9554);
            double t9642 = 1.0/pow(t9554,1.9E1/3.0);
            v_rho_a_rho_a[Q] += scale * A*rho_b*t9558*t9589*8.0-A*Dd*rho_b*t9563*t9607*(8.0/3.0)-A*rho_a*rho_b*t9558*t9608*8.0+A*Dd*rho_a*rho_b*1.0/pow(t9554,1.0E1/3.0)*t9563*(4.0E1/9.0)+A*B*1.0/pow(t9554,1.4E1/3.0)*t9564*t9558*t9630*(2.2E1/3.0)-A*B*1.0/pow(t9554,1.7E1/3.0)*t9564*t9558*t9638*(1.54E2/9.0)+A*B*t9564*t9591*t9558*(gamma_ab*(8.0/3.0)+gamma_bb*2.0-rho_b*t9606*2.0-rho_a*rho_b*(CFext*pow(rho_a,2.0/3.0)*(4.0E1/9.0)+t9569*(C*t9607*(2.0/8.1E1)+Dd*t9558*t9607*(2.0/8.1E1)-t9563*t9587*t9610*(1.0/2.7E1)+Dd*t9587*t9608*t9609*(1.0/8.1E1))-t9562*(C*t9607*(1.4E1/8.1E1)+Dd*t9558*t9607*(1.4E1/8.1E1)-t9563*t9587*t9610*(7.0/2.7E1)+Dd*t9587*t9608*t9609*(7.0/8.1E1))-t9576*t9579*(C*t9607*(4.0/8.1E1)+Dd*t9558*t9607*(4.0/8.1E1)-t9563*t9587*t9610*(2.0/2.7E1)+Dd*t9587*t9608*t9609*(2.0/8.1E1))+gamma_aa*t9582*t9589*2.0+gamma_aa*t9576*t9603*2.0-t9582*t9579*t9608*2.0-t9579*t9589*t9603*2.0))-A*rho_a*rho_b*t9591*t9587*t9609*(8.0/9.0)-A*B*t9564*t9587*t9642*t9609*t9638*(2.0/9.0)-A*B*(C*C)*t9564*t9558*t9642*t9638*(1.0/9.0)-A*B*C*t9564*t9558*t9630*t9631*(2.0/3.0)+A*B*C*t9564*t9558*t9638*t9639*(2.6E1/9.0)-A*B*Dd*t9563*t9564*t9630*t9631*(2.0/3.0)+A*B*Dd*t9563*t9564*t9638*t9639*(2.6E1/9.0)-A*B*C*Dd*t9563*t9564*t9642*t9638*(2.0/9.0);
            }
            
            // v_rho_a_rho_b
            {
            double t9644 = rho_a+rho_b;
            double t9645 = 1.0/pow(t9644,1.0/3.0);
            double t9646 = Dd*t9645;
            double t9647 = t9646+1.0;
            double t9648 = 1.0/t9647;
            double t9649 = 1.0/(t9644*t9644);
            double t9650 = 1.0/pow(t9644,7.0/3.0);
            double t9651 = 1.0/(t9647*t9647);
            double t9652 = t9644*t9644;
            double t9653 = t9652*(2.0/3.0);
            double t9654 = gamma_ab*2.0;
            double t9655 = gamma_aa+gamma_bb+t9654;
            double t9656 = 1.0/t9644;
            double t9695 = C*t9645;
            double t9657 = exp(-t9695);
            double t9658 = C*t9645*(7.0/1.8E1);
            double t9659 = Dd*t9645*t9648*(7.0/1.8E1);
            double t9660 = t9658+t9659-4.7E1/1.8E1;
            double t9661 = t9660*t9655;
            double t9662 = gamma_aa+gamma_bb;
            double t9663 = 1.0/pow(t9644,4.0/3.0);
            double t9664 = Dd*Dd;
            double t9665 = 1.0/pow(t9644,5.0/3.0);
            double t9666 = C*t9645*(1.0/9.0);
            double t9667 = Dd*t9645*t9648*(1.0/9.0);
            double t9668 = t9666+t9667-1.1E1/9.0;
            double t9669 = gamma_aa*rho_a;
            double t9670 = gamma_bb*rho_b;
            double t9671 = t9670+t9669;
            double t9672 = C*t9663*(1.0/5.4E1);
            double t9673 = Dd*t9663*t9648*(1.0/5.4E1);
            double t9699 = t9651*t9664*t9665*(1.0/5.4E1);
            double t9674 = t9672+t9673-t9699;
            double t9675 = C*t9663*(7.0/5.4E1);
            double t9676 = Dd*t9663*t9648*(7.0/5.4E1);
            double t9702 = t9651*t9664*t9665*(7.0/5.4E1);
            double t9677 = t9675+t9676-t9702;
            double t9678 = t9655*t9677;
            double t9679 = C*t9663*(1.0/2.7E1);
            double t9680 = Dd*t9663*t9648*(1.0/2.7E1);
            double t9693 = t9651*t9664*t9665*(1.0/2.7E1);
            double t9681 = t9680-t9693+t9679;
            double t9682 = t9671*t9681*t9656;
            double t9683 = t9671*t9649*t9668;
            double t9684 = C*t9645*(1.0/1.8E1);
            double t9685 = Dd*t9645*t9648*(1.0/1.8E1);
            double t9686 = t9684+t9685-5.0/2.0;
            double t9687 = pow(rho_a,8.0/3.0);
            double t9688 = pow(rho_b,8.0/3.0);
            double t9689 = t9687+t9688;
            double t9690 = 1.0/(t9644*t9644*t9644);
            double t9691 = 1.0/(t9647*t9647*t9647);
            double t9692 = 1.0/pow(t9644,8.0/3.0);
            double t9694 = t9671*t9656*t9668;
            double t9696 = t9662*t9686;
            double t9697 = CFext*t9689;
            double t9698 = rho_b*(4.0/3.0);
            double t9700 = pow(rho_a,5.0/3.0);
            double t9701 = CFext*t9700*(8.0/3.0);
            double t9712 = t9662*t9674;
            double t9720 = gamma_aa*t9656*t9668;
            double t9703 = t9682+t9683+t9678+t9701-t9720-t9712;
            double t9704 = 1.0/pow(t9644,1.4E1/3.0);
            double t9705 = rho_a*2.0;
            double t9706 = rho_b*2.0;
            double t9707 = t9705+t9706;
            double t9708 = t9655*t9707*(2.0/3.0);
            double t9709 = t9661+t9694-t9696-t9697;
            double t9710 = rho_a*(4.0/3.0);
            double t9711 = t9698+t9710;
            double t9713 = pow(rho_b,5.0/3.0);
            double t9714 = CFext*t9713*(8.0/3.0);
            double t9715 = 1.0/pow(t9644,1.1E1/3.0);
            double t9716 = rho_b*t9709;
            double t9717 = rho_a*(2.0/3.0);
            double t9718 = t9698-t9717;
            double t9719 = gamma_bb*t9718;
            double t9721 = rho_a*rho_b*t9703;
            double t9722 = gamma_aa*t9711;
            double t9723 = t9721+t9722-t9716-t9708+t9719;
            double t9724 = 1.0/(t9644*t9644*t9644*t9644*t9644);
            double t9725 = rho_a*t9709;
            double t9726 = rho_b*(2.0/3.0);
            double t9727 = t9710-t9726;
            double t9728 = gamma_aa*t9727;
            double t9729 = gamma_bb*t9711;
            double t9733 = gamma_bb*t9656*t9668;
            double t9730 = t9682+t9683+t9678-t9712+t9714-t9733;
            double t9731 = rho_a*rho_b*t9730;
            double t9732 = t9731-t9725-t9708+t9728+t9729;
            double t9734 = rho_b*rho_b;
            double t9735 = t9653-t9734;
            double t9736 = gamma_aa*t9735;
            double t9737 = rho_a*rho_a;
            double t9738 = t9653-t9737;
            double t9739 = gamma_bb*t9738;
            double t9740 = 1.0/(t9644*t9644*t9644*t9644*t9644*t9644);
            double t9742 = t9652*t9655*(2.0/3.0);
            double t9743 = rho_a*rho_b*t9709;
            double t9741 = -t9742-t9743+t9736+t9739;
            double t9744 = 1.0/pow(t9644,1.9E1/3.0);
            v_rho_a_rho_b[Q] += scale * A*t9656*t9648*-4.0+A*rho_a*t9648*t9649*4.0+A*rho_b*t9648*t9649*4.0-A*Dd*rho_a*t9650*t9651*(4.0/3.0)-A*Dd*rho_b*t9650*t9651*(4.0/3.0)-A*rho_a*rho_b*t9690*t9648*8.0+A*Dd*rho_a*rho_b*t9651*1.0/pow(t9644,1.0E1/3.0)*(4.0E1/9.0)-A*B*1.0/pow(t9644,1.7E1/3.0)*t9648*t9657*(t9736+t9739-t9652*t9655*(2.0/3.0)-rho_a*rho_b*(t9661+t9694-CFext*t9689-t9662*t9686))*(1.54E2/9.0)-A*B*t9648*t9657*t9715*(gamma_ab*(-8.0/3.0)-t9661-t9694+t9696+t9697+rho_a*t9703+rho_b*(t9682+t9683+t9678+t9714-t9662*t9674-gamma_bb*t9656*t9668)+rho_a*rho_b*(t9662*(C*t9650*(2.0/8.1E1)+Dd*t9650*t9648*(2.0/8.1E1)-t9651*t9664*t9692*(1.0/2.7E1)+Dd*t9690*t9664*t9691*(1.0/8.1E1))-t9655*(C*t9650*(1.4E1/8.1E1)+Dd*t9650*t9648*(1.4E1/8.1E1)-t9651*t9664*t9692*(7.0/2.7E1)+Dd*t9690*t9664*t9691*(7.0/8.1E1))-t9671*t9656*(C*t9650*(4.0/8.1E1)+Dd*t9650*t9648*(4.0/8.1E1)-t9651*t9664*t9692*(2.0/2.7E1)+Dd*t9690*t9664*t9691*(2.0/8.1E1))+gamma_aa*t9681*t9656+gamma_bb*t9681*t9656+gamma_aa*t9649*t9668+gamma_bb*t9649*t9668-t9671*t9681*t9649*2.0-t9671*t9690*t9668*2.0))+A*B*t9648*t9657*t9704*t9723*(1.1E1/3.0)+A*B*t9648*t9657*t9704*t9732*(1.1E1/3.0)-A*rho_a*rho_b*t9664*t9691*t9715*(8.0/9.0)-A*B*t9664*t9691*t9657*t9741*t9744*(2.0/9.0)-A*B*(C*C)*t9648*t9657*t9741*t9744*(1.0/9.0)+A*B*C*t9648*t9657*t9740*t9741*(2.6E1/9.0)-A*B*C*t9648*t9657*t9723*t9724*(1.0/3.0)-A*B*C*t9648*t9657*t9732*t9724*(1.0/3.0)+A*B*Dd*t9651*t9657*t9740*t9741*(2.6E1/9.0)-A*B*Dd*t9651*t9657*t9723*t9724*(1.0/3.0)-A*B*Dd*t9651*t9657*t9732*t9724*(1.0/3.0)-A*B*C*Dd*t9651*t9657*t9741*t9744*(2.0/9.0);
            }
            
            // v_rho_b_rho_b
            {
            double t9746 = rho_a+rho_b;
            double t9747 = 1.0/pow(t9746,1.0/3.0);
            double t9748 = Dd*t9747;
            double t9749 = t9748+1.0;
            double t9750 = 1.0/t9749;
            double t9751 = t9746*t9746;
            double t9752 = t9751*(2.0/3.0);
            double t9753 = gamma_ab*2.0;
            double t9754 = gamma_aa+gamma_bb+t9753;
            double t9755 = 1.0/(t9749*t9749);
            double t9782 = C*t9747;
            double t9756 = exp(-t9782);
            double t9757 = C*t9747*(7.0/1.8E1);
            double t9758 = Dd*t9750*t9747*(7.0/1.8E1);
            double t9759 = t9757+t9758-4.7E1/1.8E1;
            double t9760 = t9754*t9759;
            double t9761 = gamma_aa+gamma_bb;
            double t9762 = C*t9747*(1.0/1.8E1);
            double t9763 = Dd*t9750*t9747*(1.0/1.8E1);
            double t9764 = t9762+t9763-5.0/2.0;
            double t9765 = pow(rho_a,8.0/3.0);
            double t9766 = pow(rho_b,8.0/3.0);
            double t9767 = t9765+t9766;
            double t9768 = 1.0/t9746;
            double t9769 = gamma_aa*rho_a;
            double t9770 = gamma_bb*rho_b;
            double t9771 = t9770+t9769;
            double t9772 = C*t9747*(1.0/9.0);
            double t9773 = Dd*t9750*t9747*(1.0/9.0);
            double t9774 = t9772+t9773-1.1E1/9.0;
            double t9775 = t9771*t9774*t9768;
            double t9808 = t9761*t9764;
            double t9809 = CFext*t9767;
            double t9776 = t9760+t9775-t9808-t9809;
            double t9777 = rho_a*(4.0/3.0);
            double t9778 = 1.0/pow(t9746,4.0/3.0);
            double t9779 = Dd*Dd;
            double t9780 = 1.0/pow(t9746,5.0/3.0);
            double t9781 = 1.0/(t9746*t9746);
            double t9783 = 1.0/pow(t9746,1.1E1/3.0);
            double t9784 = C*t9778*(1.0/5.4E1);
            double t9785 = Dd*t9750*t9778*(1.0/5.4E1);
            double t9817 = t9780*t9755*t9779*(1.0/5.4E1);
            double t9786 = t9784+t9785-t9817;
            double t9787 = pow(rho_b,5.0/3.0);
            double t9788 = CFext*t9787*(8.0/3.0);
            double t9789 = C*t9778*(7.0/5.4E1);
            double t9790 = Dd*t9750*t9778*(7.0/5.4E1);
            double t9819 = t9780*t9755*t9779*(7.0/5.4E1);
            double t9791 = t9790+t9789-t9819;
            double t9792 = t9754*t9791;
            double t9793 = C*t9778*(1.0/2.7E1);
            double t9794 = Dd*t9750*t9778*(1.0/2.7E1);
            double t9803 = t9780*t9755*t9779*(1.0/2.7E1);
            double t9795 = t9793+t9794-t9803;
            double t9796 = t9771*t9768*t9795;
            double t9797 = t9771*t9781*t9774;
            double t9818 = t9761*t9786;
            double t9820 = gamma_bb*t9774*t9768;
            double t9798 = t9792+t9796+t9788+t9797-t9820-t9818;
            double t9799 = 1.0/pow(t9746,7.0/3.0);
            double t9800 = 1.0/(t9746*t9746*t9746);
            double t9801 = 1.0/(t9749*t9749*t9749);
            double t9802 = 1.0/pow(t9746,8.0/3.0);
            double t9804 = rho_a*2.0;
            double t9805 = rho_b*2.0;
            double t9806 = t9804+t9805;
            double t9807 = t9754*t9806*(2.0/3.0);
            double t9810 = rho_a*t9776;
            double t9811 = rho_b*(2.0/3.0);
            double t9812 = t9777-t9811;
            double t9813 = gamma_aa*t9812;
            double t9814 = rho_b*(4.0/3.0);
            double t9815 = t9777+t9814;
            double t9816 = gamma_bb*t9815;
            double t9821 = rho_a*rho_b*t9798;
            double t9822 = -t9810+t9821+t9813-t9807+t9816;
            double t9823 = 1.0/(t9746*t9746*t9746*t9746*t9746);
            double t9824 = rho_b*rho_b;
            double t9825 = t9752-t9824;
            double t9826 = gamma_aa*t9825;
            double t9827 = rho_a*rho_a;
            double t9828 = t9752-t9827;
            double t9829 = gamma_bb*t9828;
            double t9832 = t9751*t9754*(2.0/3.0);
            double t9833 = rho_a*rho_b*t9776;
            double t9830 = -t9832-t9833+t9826+t9829;
            double t9831 = 1.0/(t9746*t9746*t9746*t9746*t9746*t9746);
            double t9834 = 1.0/pow(t9746,1.9E1/3.0);
            v_rho_b_rho_b[Q] += scale * A*rho_a*t9750*t9781*8.0-A*Dd*rho_a*t9755*t9799*(8.0/3.0)-A*rho_a*rho_b*t9750*t9800*8.0+A*Dd*rho_a*rho_b*1.0/pow(t9746,1.0E1/3.0)*t9755*(4.0E1/9.0)+A*B*t9750*1.0/pow(t9746,1.4E1/3.0)*t9756*t9822*(2.2E1/3.0)-A*B*t9750*1.0/pow(t9746,1.7E1/3.0)*t9756*t9830*(1.54E2/9.0)+A*B*t9750*t9756*t9783*(gamma_aa*2.0+gamma_ab*(8.0/3.0)-rho_a*t9798*2.0-rho_a*rho_b*(CFext*pow(rho_b,2.0/3.0)*(4.0E1/9.0)+t9761*(C*t9799*(2.0/8.1E1)+Dd*t9750*t9799*(2.0/8.1E1)-t9755*t9779*t9802*(1.0/2.7E1)+Dd*t9779*t9800*t9801*(1.0/8.1E1))-t9754*(C*t9799*(1.4E1/8.1E1)+Dd*t9750*t9799*(1.4E1/8.1E1)-t9755*t9779*t9802*(7.0/2.7E1)+Dd*t9779*t9800*t9801*(7.0/8.1E1))-t9771*t9768*(C*t9799*(4.0/8.1E1)+Dd*t9750*t9799*(4.0/8.1E1)-t9755*t9779*t9802*(2.0/2.7E1)+Dd*t9779*t9800*t9801*(2.0/8.1E1))+gamma_bb*t9781*t9774*2.0+gamma_bb*t9768*t9795*2.0-t9771*t9781*t9795*2.0-t9771*t9774*t9800*2.0))-A*rho_a*rho_b*t9783*t9779*t9801*(8.0/9.0)-A*B*t9756*t9779*t9801*t9830*t9834*(2.0/9.0)-A*B*(C*C)*t9750*t9756*t9830*t9834*(1.0/9.0)+A*B*C*t9750*t9756*t9830*t9831*(2.6E1/9.0)-A*B*C*t9750*t9756*t9822*t9823*(2.0/3.0)+A*B*Dd*t9755*t9756*t9830*t9831*(2.6E1/9.0)-A*B*Dd*t9755*t9756*t9822*t9823*(2.0/3.0)-A*B*C*Dd*t9755*t9756*t9830*t9834*(2.0/9.0);
            }
            
            // v_rho_a_gamma_aa
            {
            double t9836 = rho_a+rho_b;
            double t9837 = 1.0/pow(t9836,1.0/3.0);
            double t9838 = Dd*t9837;
            double t9839 = t9838+1.0;
            double t9840 = 1.0/t9839;
            double t9854 = C*t9837;
            double t9841 = exp(-t9854);
            double t9842 = C*t9837*(1.0/3.0);
            double t9843 = Dd*t9840*t9837*(1.0/3.0);
            double t9844 = 1.0/t9836;
            double t9845 = C*t9837*(1.0/9.0);
            double t9846 = Dd*t9840*t9837*(1.0/9.0);
            double t9847 = t9845+t9846-1.1E1/9.0;
            double t9848 = rho_a*t9844*t9847;
            double t9849 = t9842+t9843+t9848-1.0/9.0;
            double t9850 = 1.0/pow(t9836,4.0/3.0);
            double t9851 = Dd*Dd;
            double t9852 = 1.0/pow(t9836,5.0/3.0);
            double t9853 = 1.0/(t9839*t9839);
            double t9855 = rho_b*rho_b;
            double t9856 = rho_a*rho_b*t9849;
            double t9857 = t9855+t9856;
            double t9858 = 1.0/(t9836*t9836*t9836*t9836*t9836);
            v_rho_a_gamma_aa[Q] += scale * A*B*t9840*t9841*1.0/pow(t9836,1.4E1/3.0)*t9857*(-1.1E1/3.0)+A*B*t9840*t9841*1.0/pow(t9836,1.1E1/3.0)*(rho_b*t9849-rho_a*rho_b*(C*t9850*(1.0/9.0)-t9844*t9847+rho_a*t9844*(C*t9850*(1.0/2.7E1)+Dd*t9840*t9850*(1.0/2.7E1)-t9851*t9852*t9853*(1.0/2.7E1))+rho_a*1.0/(t9836*t9836)*t9847+Dd*t9840*t9850*(1.0/9.0)-t9851*t9852*t9853*(1.0/9.0)))+A*B*C*t9840*t9841*t9857*t9858*(1.0/3.0)+A*B*Dd*t9841*t9853*t9857*t9858*(1.0/3.0);
            }
            
            // v_rho_a_gamma_ab
            {
            double t9860 = rho_a+rho_b;
            double t9861 = 1.0/pow(t9860,1.0/3.0);
            double t9862 = Dd*t9861;
            double t9863 = t9862+1.0;
            double t9864 = 1.0/t9863;
            double t9870 = C*t9861;
            double t9865 = exp(-t9870);
            double t9866 = C*t9861*(7.0/9.0);
            double t9867 = Dd*t9861*t9864*(7.0/9.0);
            double t9868 = t9866+t9867-4.7E1/9.0;
            double t9869 = 1.0/pow(t9860,4.0/3.0);
            double t9871 = t9860*t9860;
            double t9872 = t9871*(4.0/3.0);
            double t9873 = rho_a*rho_b*t9868;
            double t9874 = t9872+t9873;
            double t9875 = 1.0/(t9860*t9860*t9860*t9860*t9860);
            double t9876 = 1.0/(t9863*t9863);
            v_rho_a_gamma_ab[Q] += scale * A*B*1.0/pow(t9860,1.4E1/3.0)*t9864*t9865*t9874*(-1.1E1/3.0)+A*B*1.0/pow(t9860,1.1E1/3.0)*t9864*t9865*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_b*t9868-rho_a*rho_b*(C*t9869*(7.0/2.7E1)-(Dd*Dd)*1.0/pow(t9860,5.0/3.0)*t9876*(7.0/2.7E1)+Dd*t9864*t9869*(7.0/2.7E1)))+A*B*C*t9864*t9865*t9874*t9875*(1.0/3.0)+A*B*Dd*t9865*t9874*t9875*t9876*(1.0/3.0);
            }
            
            // v_rho_a_gamma_bb
            {
            double t9878 = rho_a+rho_b;
            double t9879 = 1.0/pow(t9878,1.0/3.0);
            double t9880 = Dd*t9879;
            double t9881 = t9880+1.0;
            double t9882 = 1.0/t9881;
            double t9896 = C*t9879;
            double t9883 = exp(-t9896);
            double t9884 = C*t9879*(1.0/3.0);
            double t9885 = Dd*t9882*t9879*(1.0/3.0);
            double t9886 = 1.0/t9878;
            double t9887 = C*t9879*(1.0/9.0);
            double t9888 = Dd*t9882*t9879*(1.0/9.0);
            double t9889 = t9887+t9888-1.1E1/9.0;
            double t9890 = rho_b*t9886*t9889;
            double t9891 = t9890+t9884+t9885-1.0/9.0;
            double t9892 = 1.0/pow(t9878,4.0/3.0);
            double t9893 = Dd*Dd;
            double t9894 = 1.0/pow(t9878,5.0/3.0);
            double t9895 = 1.0/(t9881*t9881);
            double t9897 = rho_a*rho_a;
            double t9898 = rho_a*rho_b*t9891;
            double t9899 = t9897+t9898;
            double t9900 = 1.0/(t9878*t9878*t9878*t9878*t9878);
            v_rho_a_gamma_bb[Q] += scale * A*B*t9882*t9883*1.0/pow(t9878,1.4E1/3.0)*t9899*(-1.1E1/3.0)+A*B*t9882*t9883*1.0/pow(t9878,1.1E1/3.0)*(rho_a*2.0+rho_b*t9891-rho_a*rho_b*(C*t9892*(1.0/9.0)+rho_b*t9886*(C*t9892*(1.0/2.7E1)+Dd*t9882*t9892*(1.0/2.7E1)-t9893*t9894*t9895*(1.0/2.7E1))+rho_b*1.0/(t9878*t9878)*t9889+Dd*t9882*t9892*(1.0/9.0)-t9893*t9894*t9895*(1.0/9.0)))+A*B*C*t9882*t9883*t9899*t9900*(1.0/3.0)+A*B*Dd*t9883*t9895*t9899*t9900*(1.0/3.0);
            }
            
            // v_rho_b_gamma_aa
            {
            double t9902 = rho_a+rho_b;
            double t9903 = 1.0/pow(t9902,1.0/3.0);
            double t9904 = Dd*t9903;
            double t9905 = t9904+1.0;
            double t9906 = 1.0/t9905;
            double t9920 = C*t9903;
            double t9907 = exp(-t9920);
            double t9908 = C*t9903*(1.0/3.0);
            double t9909 = Dd*t9903*t9906*(1.0/3.0);
            double t9910 = 1.0/t9902;
            double t9911 = C*t9903*(1.0/9.0);
            double t9912 = Dd*t9903*t9906*(1.0/9.0);
            double t9913 = t9911+t9912-1.1E1/9.0;
            double t9914 = rho_a*t9910*t9913;
            double t9915 = t9914+t9908+t9909-1.0/9.0;
            double t9916 = 1.0/pow(t9902,4.0/3.0);
            double t9917 = Dd*Dd;
            double t9918 = 1.0/pow(t9902,5.0/3.0);
            double t9919 = 1.0/(t9905*t9905);
            double t9921 = rho_b*rho_b;
            double t9922 = rho_a*rho_b*t9915;
            double t9923 = t9921+t9922;
            double t9924 = 1.0/(t9902*t9902*t9902*t9902*t9902);
            v_rho_b_gamma_aa[Q] += scale * A*B*1.0/pow(t9902,1.4E1/3.0)*t9923*t9906*t9907*(-1.1E1/3.0)+A*B*1.0/pow(t9902,1.1E1/3.0)*t9906*t9907*(rho_b*2.0+rho_a*t9915-rho_a*rho_b*(C*t9916*(1.0/9.0)+rho_a*t9910*(C*t9916*(1.0/2.7E1)+Dd*t9906*t9916*(1.0/2.7E1)-t9917*t9918*t9919*(1.0/2.7E1))+rho_a*1.0/(t9902*t9902)*t9913+Dd*t9906*t9916*(1.0/9.0)-t9917*t9918*t9919*(1.0/9.0)))+A*B*C*t9923*t9906*t9924*t9907*(1.0/3.0)+A*B*Dd*t9923*t9924*t9907*t9919*(1.0/3.0);
            }
            
            // v_rho_b_gamma_ab
            {
            double t9926 = rho_a+rho_b;
            double t9927 = 1.0/pow(t9926,1.0/3.0);
            double t9928 = Dd*t9927;
            double t9929 = t9928+1.0;
            double t9930 = 1.0/t9929;
            double t9936 = C*t9927;
            double t9931 = exp(-t9936);
            double t9932 = C*t9927*(7.0/9.0);
            double t9933 = Dd*t9930*t9927*(7.0/9.0);
            double t9934 = t9932+t9933-4.7E1/9.0;
            double t9935 = 1.0/pow(t9926,4.0/3.0);
            double t9937 = t9926*t9926;
            double t9938 = t9937*(4.0/3.0);
            double t9939 = rho_a*rho_b*t9934;
            double t9940 = t9938+t9939;
            double t9941 = 1.0/(t9926*t9926*t9926*t9926*t9926);
            double t9942 = 1.0/(t9929*t9929);
            v_rho_b_gamma_ab[Q] += scale * A*B*t9930*t9931*t9940*1.0/pow(t9926,1.4E1/3.0)*(-1.1E1/3.0)+A*B*t9930*t9931*1.0/pow(t9926,1.1E1/3.0)*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_a*t9934-rho_a*rho_b*(C*t9935*(7.0/2.7E1)-(Dd*Dd)*t9942*1.0/pow(t9926,5.0/3.0)*(7.0/2.7E1)+Dd*t9930*t9935*(7.0/2.7E1)))+A*B*C*t9930*t9931*t9940*t9941*(1.0/3.0)+A*B*Dd*t9931*t9940*t9941*t9942*(1.0/3.0);
            }
            
            // v_rho_b_gamma_bb
            {
            double t9944 = rho_a+rho_b;
            double t9945 = 1.0/pow(t9944,1.0/3.0);
            double t9946 = Dd*t9945;
            double t9947 = t9946+1.0;
            double t9948 = 1.0/t9947;
            double t9962 = C*t9945;
            double t9949 = exp(-t9962);
            double t9950 = C*t9945*(1.0/3.0);
            double t9951 = Dd*t9945*t9948*(1.0/3.0);
            double t9952 = 1.0/t9944;
            double t9953 = C*t9945*(1.0/9.0);
            double t9954 = Dd*t9945*t9948*(1.0/9.0);
            double t9955 = t9953+t9954-1.1E1/9.0;
            double t9956 = rho_b*t9952*t9955;
            double t9957 = t9950+t9951+t9956-1.0/9.0;
            double t9958 = 1.0/pow(t9944,4.0/3.0);
            double t9959 = Dd*Dd;
            double t9960 = 1.0/pow(t9944,5.0/3.0);
            double t9961 = 1.0/(t9947*t9947);
            double t9963 = rho_a*rho_a;
            double t9964 = rho_a*rho_b*t9957;
            double t9965 = t9963+t9964;
            double t9966 = 1.0/(t9944*t9944*t9944*t9944*t9944);
            v_rho_b_gamma_bb[Q] += scale * A*B*1.0/pow(t9944,1.4E1/3.0)*t9965*t9948*t9949*(-1.1E1/3.0)+A*B*1.0/pow(t9944,1.1E1/3.0)*t9948*t9949*(rho_a*t9957-rho_a*rho_b*(C*t9958*(1.0/9.0)-t9952*t9955+rho_b*t9952*(C*t9958*(1.0/2.7E1)+Dd*t9948*t9958*(1.0/2.7E1)-t9960*t9961*t9959*(1.0/2.7E1))+rho_b*1.0/(t9944*t9944)*t9955+Dd*t9948*t9958*(1.0/9.0)-t9960*t9961*t9959*(1.0/9.0)))+A*B*C*t9965*t9948*t9966*t9949*(1.0/3.0)+A*B*Dd*t9961*t9965*t9966*t9949*(1.0/3.0);
            }
            
        }
    }
}

}
