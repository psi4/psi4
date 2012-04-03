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
                double t11709 = rho_a+rho_b;
                double t11710 = 1.0/pow(t11709,1.0/3.0);
                double t11711 = Dd*t11710;
                double t11712 = t11711+1.0;
                double t11713 = 1.0/t11712;
                double t11714 = t11709*t11709;
                double t11715 = t11714*(2.0/3.0);
                double t11716 = gamma_ab*2.0;
                double t11717 = gamma_aa+gamma_bb+t11716;
                double t11718 = 1.0/t11709;
                v[Q] += scale * (A*rho_a*rho_b*t11713*t11718*-4.0-A*B*t11713*1.0/pow(t11709,1.1E1/3.0)*exp(-C*t11710)*(t11714*t11717*(-2.0/3.0)+gamma_aa*(t11715-rho_b*rho_b)+gamma_bb*(t11715-rho_a*rho_a)+rho_a*rho_b*((gamma_aa+gamma_bb)*(C*t11710*(1.0/1.8E1)+Dd*t11710*t11713*(1.0/1.8E1)-5.0/2.0)+CFext*(pow(rho_a,8.0/3.0)+pow(rho_b,8.0/3.0))-t11717*(C*t11710*(7.0/1.8E1)+Dd*t11710*t11713*(7.0/1.8E1)-4.7E1/1.8E1)-t11718*(gamma_aa*rho_a+gamma_bb*rho_b)*(C*t11710*(1.0/9.0)+Dd*t11710*t11713*(1.0/9.0)-1.1E1/9.0))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t11720 = rho_a+rho_b;
                double t11721 = 1.0/pow(t11720,1.0/3.0);
                double t11722 = Dd*t11721;
                double t11723 = t11722+1.0;
                double t11724 = 1.0/t11723;
                double t11725 = t11720*t11720;
                double t11726 = t11725*(2.0/3.0);
                double t11727 = gamma_ab*2.0;
                double t11728 = gamma_aa+gamma_bb+t11727;
                double t11729 = 1.0/t11720;
                double t11756 = C*t11721;
                double t11730 = exp(-t11756);
                double t11731 = C*t11721*(7.0/1.8E1);
                double t11732 = Dd*t11721*t11724*(7.0/1.8E1);
                double t11733 = t11731+t11732-4.7E1/1.8E1;
                double t11734 = t11733*t11728;
                double t11735 = gamma_aa+gamma_bb;
                double t11736 = C*t11721*(1.0/1.8E1);
                double t11737 = Dd*t11721*t11724*(1.0/1.8E1);
                double t11738 = t11736+t11737-5.0/2.0;
                double t11739 = pow(rho_a,8.0/3.0);
                double t11740 = pow(rho_b,8.0/3.0);
                double t11741 = t11740+t11739;
                double t11742 = gamma_aa*rho_a;
                double t11743 = gamma_bb*rho_b;
                double t11744 = t11742+t11743;
                double t11745 = C*t11721*(1.0/9.0);
                double t11746 = Dd*t11721*t11724*(1.0/9.0);
                double t11747 = t11745+t11746-1.1E1/9.0;
                double t11748 = t11744*t11729*t11747;
                double t11764 = t11735*t11738;
                double t11765 = CFext*t11741;
                double t11749 = t11734-t11764-t11765+t11748;
                double t11750 = rho_b*(4.0/3.0);
                double t11751 = 1.0/pow(t11720,4.0/3.0);
                double t11752 = 1.0/(t11723*t11723);
                double t11753 = Dd*Dd;
                double t11754 = 1.0/pow(t11720,5.0/3.0);
                double t11755 = 1.0/(t11720*t11720);
                double t11757 = rho_b*rho_b;
                double t11758 = t11726-t11757;
                double t11759 = gamma_aa*t11758;
                double t11760 = rho_a*rho_a;
                double t11761 = t11760-t11726;
                double t11762 = gamma_bb*t11761;
                double t11763 = t11725*t11728*(2.0/3.0);
                double t11766 = rho_a*rho_b*t11749;
                double t11767 = 1.0/(t11720*t11720*t11720*t11720*t11720);
                v_rho_a[Q] += scale * (A*rho_b*t11724*t11729*-4.0+A*rho_a*rho_b*t11724*t11755*4.0-A*Dd*rho_a*rho_b*1.0/pow(t11720,7.0/3.0)*t11752*(4.0/3.0)-A*B*1.0/pow(t11720,1.4E1/3.0)*t11730*t11724*(t11762+t11763+t11766-t11759)*(1.1E1/3.0)+A*B*1.0/pow(t11720,1.1E1/3.0)*t11730*t11724*(rho_b*t11749-gamma_aa*(rho_a*(4.0/3.0)+t11750)+gamma_bb*(rho_a*(2.0/3.0)-t11750)+t11728*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_a,5.0/3.0)*(8.0/3.0)-t11735*(C*t11751*(1.0/5.4E1)+Dd*t11724*t11751*(1.0/5.4E1)-t11752*t11753*t11754*(1.0/5.4E1))+t11728*(C*t11751*(7.0/5.4E1)+Dd*t11724*t11751*(7.0/5.4E1)-t11752*t11753*t11754*(7.0/5.4E1))+t11744*t11729*(C*t11751*(1.0/2.7E1)+Dd*t11724*t11751*(1.0/2.7E1)-t11752*t11753*t11754*(1.0/2.7E1))-gamma_aa*t11729*t11747+t11744*t11755*t11747))+A*B*C*t11730*t11724*t11767*(t11762+t11763+t11766-t11759)*(1.0/3.0)+A*B*Dd*t11730*t11752*t11767*(t11762+t11763+t11766-t11759)*(1.0/3.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t11769 = rho_a+rho_b;
                double t11770 = 1.0/pow(t11769,1.0/3.0);
                double t11771 = Dd*t11770;
                double t11772 = t11771+1.0;
                double t11773 = 1.0/t11772;
                double t11774 = t11769*t11769;
                double t11775 = t11774*(2.0/3.0);
                double t11776 = gamma_ab*2.0;
                double t11777 = gamma_aa+gamma_bb+t11776;
                double t11778 = 1.0/t11769;
                double t11805 = C*t11770;
                double t11779 = exp(-t11805);
                double t11780 = C*t11770*(7.0/1.8E1);
                double t11781 = Dd*t11770*t11773*(7.0/1.8E1);
                double t11782 = t11780+t11781-4.7E1/1.8E1;
                double t11783 = t11782*t11777;
                double t11784 = gamma_aa+gamma_bb;
                double t11785 = C*t11770*(1.0/1.8E1);
                double t11786 = Dd*t11770*t11773*(1.0/1.8E1);
                double t11787 = t11785+t11786-5.0/2.0;
                double t11788 = pow(rho_a,8.0/3.0);
                double t11789 = pow(rho_b,8.0/3.0);
                double t11790 = t11788+t11789;
                double t11791 = gamma_aa*rho_a;
                double t11792 = gamma_bb*rho_b;
                double t11793 = t11791+t11792;
                double t11794 = C*t11770*(1.0/9.0);
                double t11795 = Dd*t11770*t11773*(1.0/9.0);
                double t11796 = t11794+t11795-1.1E1/9.0;
                double t11797 = t11793*t11778*t11796;
                double t11813 = t11784*t11787;
                double t11814 = CFext*t11790;
                double t11798 = -t11813-t11814+t11783+t11797;
                double t11799 = rho_a*(4.0/3.0);
                double t11800 = 1.0/pow(t11769,4.0/3.0);
                double t11801 = 1.0/(t11772*t11772);
                double t11802 = Dd*Dd;
                double t11803 = 1.0/pow(t11769,5.0/3.0);
                double t11804 = 1.0/(t11769*t11769);
                double t11806 = rho_b*rho_b;
                double t11807 = t11806-t11775;
                double t11808 = gamma_aa*t11807;
                double t11809 = rho_a*rho_a;
                double t11810 = t11809-t11775;
                double t11811 = gamma_bb*t11810;
                double t11812 = t11774*t11777*(2.0/3.0);
                double t11815 = rho_a*rho_b*t11798;
                double t11816 = 1.0/(t11769*t11769*t11769*t11769*t11769);
                v_rho_b[Q] += scale * (A*rho_a*t11773*t11778*-4.0+A*rho_a*rho_b*t11804*t11773*4.0-A*Dd*rho_a*rho_b*t11801*1.0/pow(t11769,7.0/3.0)*(4.0/3.0)-A*B*t11773*1.0/pow(t11769,1.4E1/3.0)*t11779*(t11811+t11812+t11815+t11808)*(1.1E1/3.0)+A*B*t11773*1.0/pow(t11769,1.1E1/3.0)*t11779*(rho_a*t11798-gamma_bb*(rho_b*(4.0/3.0)+t11799)+gamma_aa*(rho_b*(2.0/3.0)-t11799)+t11777*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_b,5.0/3.0)*(8.0/3.0)-t11784*(C*t11800*(1.0/5.4E1)+Dd*t11800*t11773*(1.0/5.4E1)-t11801*t11802*t11803*(1.0/5.4E1))+t11777*(C*t11800*(7.0/5.4E1)+Dd*t11800*t11773*(7.0/5.4E1)-t11801*t11802*t11803*(7.0/5.4E1))+t11793*t11778*(C*t11800*(1.0/2.7E1)+Dd*t11800*t11773*(1.0/2.7E1)-t11801*t11802*t11803*(1.0/2.7E1))-gamma_bb*t11778*t11796+t11804*t11793*t11796))+A*B*C*t11816*t11773*t11779*(t11811+t11812+t11815+t11808)*(1.0/3.0)+A*B*Dd*t11801*t11816*t11779*(t11811+t11812+t11815+t11808)*(1.0/3.0));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t11818 = rho_a+rho_b;
                double t11819 = 1.0/pow(t11818,1.0/3.0);
                double t11820 = Dd*t11819;
                double t11821 = t11820+1.0;
                double t11822 = 1.0/t11821;
                v_gamma_aa[Q] += scale * (A*B*t11822*1.0/pow(t11818,1.1E1/3.0)*exp(-C*t11819)*(rho_b*rho_b+rho_a*rho_b*(C*t11819*(1.0/3.0)+Dd*t11822*t11819*(1.0/3.0)+(rho_a*(C*t11819*(1.0/9.0)+Dd*t11822*t11819*(1.0/9.0)-1.1E1/9.0))/t11818-1.0/9.0)));
            }
            
            // v_gamma_ab
            if (deriv >= 1) {
                double t11824 = rho_a+rho_b;
                double t11825 = 1.0/pow(t11824,1.0/3.0);
                double t11826 = Dd*t11825;
                double t11827 = t11826+1.0;
                double t11828 = 1.0/t11827;
                v_gamma_ab[Q] += scale * (A*B*1.0/pow(t11824,1.1E1/3.0)*t11828*exp(-C*t11825)*((t11824*t11824)*(4.0/3.0)+rho_a*rho_b*(C*t11825*(7.0/9.0)+Dd*t11825*t11828*(7.0/9.0)-4.7E1/9.0)));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t11830 = rho_a+rho_b;
                double t11831 = 1.0/pow(t11830,1.0/3.0);
                double t11832 = Dd*t11831;
                double t11833 = t11832+1.0;
                double t11834 = 1.0/t11833;
                v_gamma_bb[Q] += scale * (A*B*1.0/pow(t11830,1.1E1/3.0)*t11834*exp(-C*t11831)*(rho_a*rho_a+rho_a*rho_b*(C*t11831*(1.0/3.0)+Dd*t11831*t11834*(1.0/3.0)+(rho_b*(C*t11831*(1.0/9.0)+Dd*t11831*t11834*(1.0/9.0)-1.1E1/9.0))/t11830-1.0/9.0)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t11838 = rho_a+rho_b;
                double t11839 = 1.0/pow(t11838,1.0/3.0);
                double t11840 = Dd*t11839;
                double t11841 = t11840+1.0;
                double t11842 = 1.0/t11841;
                double t11843 = t11838*t11838;
                double t11844 = t11843*(2.0/3.0);
                double t11845 = gamma_ab*2.0;
                double t11846 = gamma_aa+gamma_bb+t11845;
                double t11847 = 1.0/(t11841*t11841);
                double t11874 = C*t11839;
                double t11848 = exp(-t11874);
                double t11849 = C*t11839*(7.0/1.8E1);
                double t11850 = Dd*t11842*t11839*(7.0/1.8E1);
                double t11851 = t11850+t11849-4.7E1/1.8E1;
                double t11852 = t11851*t11846;
                double t11853 = gamma_aa+gamma_bb;
                double t11854 = C*t11839*(1.0/1.8E1);
                double t11855 = Dd*t11842*t11839*(1.0/1.8E1);
                double t11856 = t11854+t11855-5.0/2.0;
                double t11857 = pow(rho_a,8.0/3.0);
                double t11858 = pow(rho_b,8.0/3.0);
                double t11859 = t11857+t11858;
                double t11860 = 1.0/t11838;
                double t11861 = gamma_aa*rho_a;
                double t11862 = gamma_bb*rho_b;
                double t11863 = t11861+t11862;
                double t11864 = C*t11839*(1.0/9.0);
                double t11865 = Dd*t11842*t11839*(1.0/9.0);
                double t11866 = t11864+t11865-1.1E1/9.0;
                double t11867 = t11860*t11863*t11866;
                double t11900 = t11853*t11856;
                double t11901 = CFext*t11859;
                double t11868 = -t11900-t11901+t11852+t11867;
                double t11869 = rho_b*(4.0/3.0);
                double t11870 = 1.0/pow(t11838,4.0/3.0);
                double t11871 = Dd*Dd;
                double t11872 = 1.0/pow(t11838,5.0/3.0);
                double t11873 = 1.0/(t11838*t11838);
                double t11875 = 1.0/pow(t11838,1.1E1/3.0);
                double t11876 = C*t11870*(1.0/5.4E1);
                double t11877 = Dd*t11842*t11870*(1.0/5.4E1);
                double t11908 = t11871*t11872*t11847*(1.0/5.4E1);
                double t11878 = -t11908+t11876+t11877;
                double t11879 = pow(rho_a,5.0/3.0);
                double t11880 = CFext*t11879*(8.0/3.0);
                double t11881 = C*t11870*(7.0/5.4E1);
                double t11882 = Dd*t11842*t11870*(7.0/5.4E1);
                double t11910 = t11871*t11872*t11847*(7.0/5.4E1);
                double t11883 = -t11910+t11881+t11882;
                double t11884 = t11846*t11883;
                double t11885 = C*t11870*(1.0/2.7E1);
                double t11886 = Dd*t11842*t11870*(1.0/2.7E1);
                double t11895 = t11871*t11872*t11847*(1.0/2.7E1);
                double t11887 = t11885+t11886-t11895;
                double t11888 = t11860*t11863*t11887;
                double t11889 = t11863*t11873*t11866;
                double t11909 = t11853*t11878;
                double t11911 = gamma_aa*t11860*t11866;
                double t11890 = -t11911+t11880-t11909+t11884+t11888+t11889;
                double t11891 = 1.0/pow(t11838,7.0/3.0);
                double t11892 = 1.0/(t11838*t11838*t11838);
                double t11893 = 1.0/(t11841*t11841*t11841);
                double t11894 = 1.0/pow(t11838,8.0/3.0);
                double t11896 = rho_a*2.0;
                double t11897 = rho_b*2.0;
                double t11898 = t11896+t11897;
                double t11899 = t11846*t11898*(2.0/3.0);
                double t11902 = rho_b*t11868;
                double t11903 = rho_a*(4.0/3.0);
                double t11904 = t11903+t11869;
                double t11905 = rho_a*(2.0/3.0);
                double t11906 = t11905-t11869;
                double t11907 = gamma_bb*t11906;
                double t11914 = gamma_aa*t11904;
                double t11915 = rho_a*rho_b*t11890;
                double t11912 = t11902-t11914-t11915+t11907+t11899;
                double t11913 = 1.0/(t11838*t11838*t11838*t11838*t11838);
                double t11916 = rho_b*rho_b;
                double t11917 = t11844-t11916;
                double t11918 = gamma_aa*t11917;
                double t11919 = rho_a*rho_a;
                double t11920 = t11844-t11919;
                double t11921 = gamma_bb*t11920;
                double t11924 = t11843*t11846*(2.0/3.0);
                double t11925 = rho_a*rho_b*t11868;
                double t11922 = t11921-t11924-t11925+t11918;
                double t11923 = 1.0/(t11838*t11838*t11838*t11838*t11838*t11838);
                double t11926 = 1.0/pow(t11838,1.9E1/3.0);
                v_rho_a_rho_a[Q] += scale * (A*rho_b*t11842*t11873*8.0-A*Dd*rho_b*t11891*t11847*(8.0/3.0)-A*rho_a*rho_b*t11842*t11892*8.0+A*Dd*rho_a*rho_b*1.0/pow(t11838,1.0E1/3.0)*t11847*(4.0E1/9.0)-A*B*t11912*t11842*1.0/pow(t11838,1.4E1/3.0)*t11848*(2.2E1/3.0)-A*B*t11922*t11842*1.0/pow(t11838,1.7E1/3.0)*t11848*(1.54E2/9.0)+A*B*t11842*t11848*t11875*(gamma_ab*(8.0/3.0)+gamma_bb*2.0-rho_b*t11890*2.0-rho_a*rho_b*(CFext*pow(rho_a,2.0/3.0)*(4.0E1/9.0)+t11853*(C*t11891*(2.0/8.1E1)+Dd*t11842*t11891*(2.0/8.1E1)-t11871*t11847*t11894*(1.0/2.7E1)+Dd*t11871*t11892*t11893*(1.0/8.1E1))-t11846*(C*t11891*(1.4E1/8.1E1)+Dd*t11842*t11891*(1.4E1/8.1E1)-t11871*t11847*t11894*(7.0/2.7E1)+Dd*t11871*t11892*t11893*(7.0/8.1E1))-t11860*t11863*(C*t11891*(4.0/8.1E1)+Dd*t11842*t11891*(4.0/8.1E1)-t11871*t11847*t11894*(2.0/2.7E1)+Dd*t11871*t11892*t11893*(2.0/8.1E1))+gamma_aa*t11860*t11887*2.0+gamma_aa*t11873*t11866*2.0-t11863*t11892*t11866*2.0-t11863*t11873*t11887*2.0))-A*rho_a*rho_b*t11871*t11875*t11893*(8.0/9.0)-A*B*t11922*t11871*t11926*t11848*t11893*(2.0/9.0)-A*B*(C*C)*t11922*t11842*t11926*t11848*(1.0/9.0)+A*B*C*t11912*t11913*t11842*t11848*(2.0/3.0)+A*B*C*t11922*t11842*t11923*t11848*(2.6E1/9.0)+A*B*Dd*t11912*t11913*t11847*t11848*(2.0/3.0)+A*B*Dd*t11922*t11923*t11847*t11848*(2.6E1/9.0)-A*B*C*Dd*t11922*t11926*t11847*t11848*(2.0/9.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t11928 = rho_a+rho_b;
                double t11929 = 1.0/pow(t11928,1.0/3.0);
                double t11930 = Dd*t11929;
                double t11931 = t11930+1.0;
                double t11932 = 1.0/t11931;
                double t11933 = 1.0/(t11928*t11928);
                double t11934 = 1.0/pow(t11928,7.0/3.0);
                double t11935 = 1.0/(t11931*t11931);
                double t11936 = t11928*t11928;
                double t11937 = t11936*(2.0/3.0);
                double t11938 = gamma_ab*2.0;
                double t11939 = gamma_aa+gamma_bb+t11938;
                double t11940 = 1.0/t11928;
                double t11979 = C*t11929;
                double t11941 = exp(-t11979);
                double t11942 = C*t11929*(7.0/1.8E1);
                double t11943 = Dd*t11932*t11929*(7.0/1.8E1);
                double t11944 = t11942+t11943-4.7E1/1.8E1;
                double t11945 = t11944*t11939;
                double t11946 = gamma_aa+gamma_bb;
                double t11947 = 1.0/pow(t11928,4.0/3.0);
                double t11948 = Dd*Dd;
                double t11949 = 1.0/pow(t11928,5.0/3.0);
                double t11950 = C*t11929*(1.0/9.0);
                double t11951 = Dd*t11932*t11929*(1.0/9.0);
                double t11952 = t11950+t11951-1.1E1/9.0;
                double t11953 = gamma_aa*rho_a;
                double t11954 = gamma_bb*rho_b;
                double t11955 = t11953+t11954;
                double t11956 = C*t11947*(1.0/5.4E1);
                double t11957 = Dd*t11932*t11947*(1.0/5.4E1);
                double t11983 = t11935*t11948*t11949*(1.0/5.4E1);
                double t11958 = t11956-t11983+t11957;
                double t11959 = C*t11947*(7.0/5.4E1);
                double t11960 = Dd*t11932*t11947*(7.0/5.4E1);
                double t11986 = t11935*t11948*t11949*(7.0/5.4E1);
                double t11961 = t11960+t11959-t11986;
                double t11962 = t11961*t11939;
                double t11963 = C*t11947*(1.0/2.7E1);
                double t11964 = Dd*t11932*t11947*(1.0/2.7E1);
                double t11977 = t11935*t11948*t11949*(1.0/2.7E1);
                double t11965 = t11963+t11964-t11977;
                double t11966 = t11940*t11955*t11965;
                double t11967 = t11933*t11952*t11955;
                double t11968 = C*t11929*(1.0/1.8E1);
                double t11969 = Dd*t11932*t11929*(1.0/1.8E1);
                double t11970 = t11968+t11969-5.0/2.0;
                double t11971 = pow(rho_a,8.0/3.0);
                double t11972 = pow(rho_b,8.0/3.0);
                double t11973 = t11971+t11972;
                double t11974 = 1.0/(t11928*t11928*t11928);
                double t11975 = 1.0/(t11931*t11931*t11931);
                double t11976 = 1.0/pow(t11928,8.0/3.0);
                double t11978 = t11940*t11952*t11955;
                double t11980 = t11970*t11946;
                double t11981 = CFext*t11973;
                double t11982 = rho_b*(4.0/3.0);
                double t11984 = pow(rho_a,5.0/3.0);
                double t11985 = CFext*t11984*(8.0/3.0);
                double t11994 = t11946*t11958;
                double t12004 = gamma_aa*t11940*t11952;
                double t11987 = t11962+t11966+t11967+t11985-t11994-t12004;
                double t11988 = 1.0/pow(t11928,1.4E1/3.0);
                double t11989 = rho_a*2.0;
                double t11990 = rho_b*2.0;
                double t11991 = t11990+t11989;
                double t11992 = rho_a*(4.0/3.0);
                double t11993 = t11982+t11992;
                double t11995 = pow(rho_b,5.0/3.0);
                double t11996 = CFext*t11995*(8.0/3.0);
                double t11997 = 1.0/pow(t11928,1.1E1/3.0);
                double t11998 = t11980-t11945+t11981-t11978;
                double t11999 = rho_b*t11998;
                double t12000 = gamma_aa*t11993;
                double t12001 = rho_a*(2.0/3.0);
                double t12002 = t11982-t12001;
                double t12003 = gamma_bb*t12002;
                double t12005 = rho_a*rho_b*t11987;
                double t12008 = t11991*t11939*(2.0/3.0);
                double t12006 = t11999+t12000+t12003+t12005-t12008;
                double t12007 = 1.0/(t11928*t11928*t11928*t11928*t11928);
                double t12009 = rho_b*(2.0/3.0);
                double t12010 = t11992-t12009;
                double t12011 = gamma_aa*t12010;
                double t12012 = gamma_bb*t11993;
                double t12016 = gamma_bb*t11940*t11952;
                double t12013 = t11962+t11966+t11967-t11994+t11996-t12016;
                double t12014 = rho_a*rho_b*t12013;
                double t12015 = rho_a*t11998;
                double t12017 = t12011+t12012+t12014+t12015-t12008;
                double t12018 = rho_b*rho_b;
                double t12019 = t11937-t12018;
                double t12020 = gamma_aa*t12019;
                double t12021 = rho_a*rho_a;
                double t12022 = t11937-t12021;
                double t12023 = gamma_bb*t12022;
                double t12024 = 1.0/(t11928*t11928*t11928*t11928*t11928*t11928);
                double t12025 = rho_a*rho_b*(t11980-t11945+t11981-t11978);
                double t12027 = t11936*t11939*(2.0/3.0);
                double t12026 = t12020+t12023+t12025-t12027;
                double t12028 = 1.0/pow(t11928,1.9E1/3.0);
                v_rho_a_rho_b[Q] += scale * (A*t11940*t11932*-4.0+A*rho_a*t11932*t11933*4.0+A*rho_b*t11932*t11933*4.0-A*Dd*rho_a*t11934*t11935*(4.0/3.0)-A*Dd*rho_b*t11934*t11935*(4.0/3.0)-A*rho_a*rho_b*t11932*t11974*8.0+A*Dd*rho_a*rho_b*t11935*1.0/pow(t11928,1.0E1/3.0)*(4.0E1/9.0)-A*B*t11932*t11941*1.0/pow(t11928,1.7E1/3.0)*(t12020+t12023-t11936*t11939*(2.0/3.0)-rho_a*rho_b*(t11945+t11978-CFext*t11973-t11970*t11946))*(1.54E2/9.0)-A*B*t11932*t11941*t11997*(gamma_ab*(-8.0/3.0)+t11980-t11945+t11981-t11978+rho_a*t11987+rho_b*(t11962+t11966+t11967+t11996-t11946*t11958-gamma_bb*t11940*t11952)+rho_a*rho_b*(t11946*(C*t11934*(2.0/8.1E1)+Dd*t11932*t11934*(2.0/8.1E1)-t11935*t11948*t11976*(1.0/2.7E1)+Dd*t11974*t11948*t11975*(1.0/8.1E1))-t11939*(C*t11934*(1.4E1/8.1E1)+Dd*t11932*t11934*(1.4E1/8.1E1)-t11935*t11948*t11976*(7.0/2.7E1)+Dd*t11974*t11948*t11975*(7.0/8.1E1))-t11940*t11955*(C*t11934*(4.0/8.1E1)+Dd*t11932*t11934*(4.0/8.1E1)-t11935*t11948*t11976*(2.0/2.7E1)+Dd*t11974*t11948*t11975*(2.0/8.1E1))+gamma_aa*t11933*t11952+gamma_aa*t11940*t11965+gamma_bb*t11933*t11952+gamma_bb*t11940*t11965-t11933*t11955*t11965*2.0-t11952*t11955*t11974*2.0))+A*B*t11932*t11941*t11988*t12006*(1.1E1/3.0)+A*B*t11932*t11941*t11988*(t12011+t12012+t12014+t12015-t11991*t11939*(2.0/3.0))*(1.1E1/3.0)-A*rho_a*rho_b*t11948*t11975*t11997*(8.0/9.0)-A*B*t11941*t11948*t11975*t12026*t12028*(2.0/9.0)-A*B*(C*C)*t11932*t11941*t12026*t12028*(1.0/9.0)-A*B*C*t11932*t11941*t12006*t12007*(1.0/3.0)+A*B*C*t11932*t11941*t12024*t12026*(2.6E1/9.0)-A*B*C*t11932*t11941*t12007*t12017*(1.0/3.0)-A*B*Dd*t11941*t11935*t12006*t12007*(1.0/3.0)+A*B*Dd*t11941*t11935*t12024*t12026*(2.6E1/9.0)-A*B*Dd*t11941*t11935*t12007*t12017*(1.0/3.0)-A*B*C*Dd*t11941*t11935*t12026*t12028*(2.0/9.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t12030 = rho_a+rho_b;
                double t12031 = 1.0/pow(t12030,1.0/3.0);
                double t12032 = Dd*t12031;
                double t12033 = t12032+1.0;
                double t12034 = 1.0/t12033;
                double t12035 = t12030*t12030;
                double t12036 = t12035*(2.0/3.0);
                double t12037 = gamma_ab*2.0;
                double t12038 = gamma_aa+gamma_bb+t12037;
                double t12039 = 1.0/(t12033*t12033);
                double t12066 = C*t12031;
                double t12040 = exp(-t12066);
                double t12041 = C*t12031*(7.0/1.8E1);
                double t12042 = Dd*t12031*t12034*(7.0/1.8E1);
                double t12043 = t12041+t12042-4.7E1/1.8E1;
                double t12044 = t12043*t12038;
                double t12045 = gamma_aa+gamma_bb;
                double t12046 = C*t12031*(1.0/1.8E1);
                double t12047 = Dd*t12031*t12034*(1.0/1.8E1);
                double t12048 = t12046+t12047-5.0/2.0;
                double t12049 = pow(rho_a,8.0/3.0);
                double t12050 = pow(rho_b,8.0/3.0);
                double t12051 = t12050+t12049;
                double t12052 = 1.0/t12030;
                double t12053 = gamma_aa*rho_a;
                double t12054 = gamma_bb*rho_b;
                double t12055 = t12053+t12054;
                double t12056 = C*t12031*(1.0/9.0);
                double t12057 = Dd*t12031*t12034*(1.0/9.0);
                double t12058 = t12056+t12057-1.1E1/9.0;
                double t12059 = t12052*t12055*t12058;
                double t12092 = t12045*t12048;
                double t12093 = CFext*t12051;
                double t12060 = t12044-t12092-t12093+t12059;
                double t12061 = rho_a*(4.0/3.0);
                double t12062 = 1.0/pow(t12030,4.0/3.0);
                double t12063 = Dd*Dd;
                double t12064 = 1.0/pow(t12030,5.0/3.0);
                double t12065 = 1.0/(t12030*t12030);
                double t12067 = 1.0/pow(t12030,1.1E1/3.0);
                double t12068 = C*t12062*(1.0/5.4E1);
                double t12069 = Dd*t12034*t12062*(1.0/5.4E1);
                double t12101 = t12063*t12064*t12039*(1.0/5.4E1);
                double t12070 = -t12101+t12068+t12069;
                double t12071 = pow(rho_b,5.0/3.0);
                double t12072 = CFext*t12071*(8.0/3.0);
                double t12073 = C*t12062*(7.0/5.4E1);
                double t12074 = Dd*t12034*t12062*(7.0/5.4E1);
                double t12103 = t12063*t12064*t12039*(7.0/5.4E1);
                double t12075 = -t12103+t12073+t12074;
                double t12076 = t12038*t12075;
                double t12077 = C*t12062*(1.0/2.7E1);
                double t12078 = Dd*t12034*t12062*(1.0/2.7E1);
                double t12087 = t12063*t12064*t12039*(1.0/2.7E1);
                double t12079 = t12077+t12078-t12087;
                double t12080 = t12052*t12055*t12079;
                double t12081 = t12055*t12065*t12058;
                double t12102 = t12070*t12045;
                double t12104 = gamma_bb*t12052*t12058;
                double t12082 = -t12102-t12104+t12080+t12072+t12081+t12076;
                double t12083 = 1.0/pow(t12030,7.0/3.0);
                double t12084 = 1.0/(t12030*t12030*t12030);
                double t12085 = 1.0/(t12033*t12033*t12033);
                double t12086 = 1.0/pow(t12030,8.0/3.0);
                double t12088 = rho_a*2.0;
                double t12089 = rho_b*2.0;
                double t12090 = t12088+t12089;
                double t12091 = t12090*t12038*(2.0/3.0);
                double t12094 = rho_a*t12060;
                double t12095 = rho_b*(2.0/3.0);
                double t12096 = t12061-t12095;
                double t12097 = gamma_aa*t12096;
                double t12098 = rho_b*(4.0/3.0);
                double t12099 = t12061+t12098;
                double t12100 = gamma_bb*t12099;
                double t12105 = rho_a*rho_b*t12082;
                double t12106 = t12100+t12105-t12091-t12094+t12097;
                double t12107 = 1.0/(t12030*t12030*t12030*t12030*t12030);
                double t12108 = rho_b*rho_b;
                double t12109 = t12036-t12108;
                double t12110 = gamma_aa*t12109;
                double t12111 = rho_a*rho_a;
                double t12112 = t12111-t12036;
                double t12113 = gamma_bb*t12112;
                double t12114 = t12035*t12038*(2.0/3.0);
                double t12115 = rho_a*rho_b*t12060;
                double t12116 = -t12110+t12113+t12114+t12115;
                double t12117 = 1.0/(t12030*t12030*t12030*t12030*t12030*t12030);
                double t12118 = 1.0/pow(t12030,1.9E1/3.0);
                v_rho_b_rho_b[Q] += scale * (A*rho_a*t12034*t12065*8.0-A*Dd*rho_a*t12083*t12039*(8.0/3.0)-A*rho_a*rho_b*t12034*t12084*8.0+A*Dd*rho_a*rho_b*1.0/pow(t12030,1.0E1/3.0)*t12039*(4.0E1/9.0)+A*B*1.0/pow(t12030,1.4E1/3.0)*t12040*t12034*t12106*(2.2E1/3.0)+A*B*1.0/pow(t12030,1.7E1/3.0)*t12040*t12034*t12116*(1.54E2/9.0)+A*B*t12040*t12034*t12067*(gamma_aa*2.0+gamma_ab*(8.0/3.0)-rho_a*t12082*2.0-rho_a*rho_b*(CFext*pow(rho_b,2.0/3.0)*(4.0E1/9.0)+t12045*(C*t12083*(2.0/8.1E1)+Dd*t12034*t12083*(2.0/8.1E1)-t12063*t12039*t12086*(1.0/2.7E1)+Dd*t12063*t12084*t12085*(1.0/8.1E1))-t12038*(C*t12083*(1.4E1/8.1E1)+Dd*t12034*t12083*(1.4E1/8.1E1)-t12063*t12039*t12086*(7.0/2.7E1)+Dd*t12063*t12084*t12085*(7.0/8.1E1))-t12052*t12055*(C*t12083*(4.0/8.1E1)+Dd*t12034*t12083*(4.0/8.1E1)-t12063*t12039*t12086*(2.0/2.7E1)+Dd*t12063*t12084*t12085*(2.0/8.1E1))+gamma_bb*t12052*t12079*2.0+gamma_bb*t12065*t12058*2.0-t12055*t12084*t12058*2.0-t12055*t12065*t12079*2.0))-A*rho_a*rho_b*t12063*t12067*t12085*(8.0/9.0)+A*B*t12040*t12063*t12118*t12085*(-t12110+t12113+t12114+t12115)*(2.0/9.0)+A*B*(C*C)*t12040*t12034*t12118*(-t12110+t12113+t12114+t12115)*(1.0/9.0)-A*B*C*t12040*t12034*t12106*t12107*(2.0/3.0)-A*B*C*t12040*t12034*t12116*t12117*(2.6E1/9.0)-A*B*Dd*t12040*t12106*t12107*t12039*(2.0/3.0)-A*B*Dd*t12040*t12116*t12117*t12039*(2.6E1/9.0)+A*B*C*Dd*t12040*t12118*t12039*(-t12110+t12113+t12114+t12115)*(2.0/9.0));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t12120 = rho_a+rho_b;
                double t12121 = 1.0/pow(t12120,1.0/3.0);
                double t12122 = Dd*t12121;
                double t12123 = t12122+1.0;
                double t12124 = 1.0/t12123;
                double t12138 = C*t12121;
                double t12125 = exp(-t12138);
                double t12126 = C*t12121*(1.0/3.0);
                double t12127 = Dd*t12121*t12124*(1.0/3.0);
                double t12128 = 1.0/t12120;
                double t12129 = C*t12121*(1.0/9.0);
                double t12130 = Dd*t12121*t12124*(1.0/9.0);
                double t12131 = t12130+t12129-1.1E1/9.0;
                double t12132 = rho_a*t12131*t12128;
                double t12133 = t12132+t12126+t12127-1.0/9.0;
                double t12134 = 1.0/pow(t12120,4.0/3.0);
                double t12135 = Dd*Dd;
                double t12136 = 1.0/pow(t12120,5.0/3.0);
                double t12137 = 1.0/(t12123*t12123);
                double t12139 = rho_b*rho_b;
                double t12140 = rho_a*rho_b*t12133;
                double t12141 = t12140+t12139;
                double t12142 = 1.0/(t12120*t12120*t12120*t12120*t12120);
                v_rho_a_gamma_aa[Q] += scale * (A*B*1.0/pow(t12120,1.4E1/3.0)*t12141*t12124*t12125*(-1.1E1/3.0)+A*B*1.0/pow(t12120,1.1E1/3.0)*t12124*t12125*(rho_b*t12133-rho_a*rho_b*(C*t12134*(1.0/9.0)-t12131*t12128+rho_a*t12128*(C*t12134*(1.0/2.7E1)+Dd*t12124*t12134*(1.0/2.7E1)-t12135*t12136*t12137*(1.0/2.7E1))+rho_a*1.0/(t12120*t12120)*t12131+Dd*t12124*t12134*(1.0/9.0)-t12135*t12136*t12137*(1.0/9.0)))+A*B*C*t12141*t12124*t12142*t12125*(1.0/3.0)+A*B*Dd*t12141*t12142*t12125*t12137*(1.0/3.0));
            }
            
            // v_rho_a_gamma_ab
            if (deriv >= 2) {
                double t12144 = rho_a+rho_b;
                double t12145 = 1.0/pow(t12144,1.0/3.0);
                double t12146 = Dd*t12145;
                double t12147 = t12146+1.0;
                double t12148 = 1.0/t12147;
                double t12154 = C*t12145;
                double t12149 = exp(-t12154);
                double t12150 = C*t12145*(7.0/9.0);
                double t12151 = Dd*t12145*t12148*(7.0/9.0);
                double t12152 = t12150+t12151-4.7E1/9.0;
                double t12153 = 1.0/pow(t12144,4.0/3.0);
                double t12155 = t12144*t12144;
                double t12156 = t12155*(4.0/3.0);
                double t12157 = rho_a*rho_b*t12152;
                double t12158 = t12156+t12157;
                double t12159 = 1.0/(t12144*t12144*t12144*t12144*t12144);
                double t12160 = 1.0/(t12147*t12147);
                v_rho_a_gamma_ab[Q] += scale * (A*B*1.0/pow(t12144,1.4E1/3.0)*t12148*t12149*t12158*(-1.1E1/3.0)+A*B*1.0/pow(t12144,1.1E1/3.0)*t12148*t12149*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_b*t12152-rho_a*rho_b*(C*t12153*(7.0/2.7E1)-(Dd*Dd)*t12160*1.0/pow(t12144,5.0/3.0)*(7.0/2.7E1)+Dd*t12153*t12148*(7.0/2.7E1)))+A*B*C*t12148*t12149*t12158*t12159*(1.0/3.0)+A*B*Dd*t12160*t12149*t12158*t12159*(1.0/3.0));
            }
            
            // v_rho_a_gamma_bb
            if (deriv >= 2) {
                double t12162 = rho_a+rho_b;
                double t12163 = 1.0/pow(t12162,1.0/3.0);
                double t12164 = Dd*t12163;
                double t12165 = t12164+1.0;
                double t12166 = 1.0/t12165;
                double t12180 = C*t12163;
                double t12167 = exp(-t12180);
                double t12168 = C*t12163*(1.0/3.0);
                double t12169 = Dd*t12163*t12166*(1.0/3.0);
                double t12170 = 1.0/t12162;
                double t12171 = C*t12163*(1.0/9.0);
                double t12172 = Dd*t12163*t12166*(1.0/9.0);
                double t12173 = t12171+t12172-1.1E1/9.0;
                double t12174 = rho_b*t12170*t12173;
                double t12175 = t12174+t12168+t12169-1.0/9.0;
                double t12176 = 1.0/pow(t12162,4.0/3.0);
                double t12177 = Dd*Dd;
                double t12178 = 1.0/pow(t12162,5.0/3.0);
                double t12179 = 1.0/(t12165*t12165);
                double t12181 = rho_a*rho_a;
                double t12182 = rho_a*rho_b*t12175;
                double t12183 = t12181+t12182;
                double t12184 = 1.0/(t12162*t12162*t12162*t12162*t12162);
                v_rho_a_gamma_bb[Q] += scale * (A*B*1.0/pow(t12162,1.4E1/3.0)*t12183*t12166*t12167*(-1.1E1/3.0)+A*B*1.0/pow(t12162,1.1E1/3.0)*t12166*t12167*(rho_a*2.0+rho_b*t12175-rho_a*rho_b*(C*t12176*(1.0/9.0)+rho_b*t12170*(C*t12176*(1.0/2.7E1)+Dd*t12166*t12176*(1.0/2.7E1)-t12177*t12178*t12179*(1.0/2.7E1))+rho_b*1.0/(t12162*t12162)*t12173+Dd*t12166*t12176*(1.0/9.0)-t12177*t12178*t12179*(1.0/9.0)))+A*B*C*t12183*t12166*t12184*t12167*(1.0/3.0)+A*B*Dd*t12183*t12184*t12167*t12179*(1.0/3.0));
            }
            
            // v_rho_b_gamma_aa
            if (deriv >= 2) {
                double t12186 = rho_a+rho_b;
                double t12187 = 1.0/pow(t12186,1.0/3.0);
                double t12188 = Dd*t12187;
                double t12189 = t12188+1.0;
                double t12190 = 1.0/t12189;
                double t12204 = C*t12187;
                double t12191 = exp(-t12204);
                double t12192 = C*t12187*(1.0/3.0);
                double t12193 = Dd*t12190*t12187*(1.0/3.0);
                double t12194 = 1.0/t12186;
                double t12195 = C*t12187*(1.0/9.0);
                double t12196 = Dd*t12190*t12187*(1.0/9.0);
                double t12197 = t12195+t12196-1.1E1/9.0;
                double t12198 = rho_a*t12194*t12197;
                double t12199 = t12192+t12193+t12198-1.0/9.0;
                double t12200 = 1.0/pow(t12186,4.0/3.0);
                double t12201 = Dd*Dd;
                double t12202 = 1.0/pow(t12186,5.0/3.0);
                double t12203 = 1.0/(t12189*t12189);
                double t12205 = rho_b*rho_b;
                double t12206 = rho_a*rho_b*t12199;
                double t12207 = t12205+t12206;
                double t12208 = 1.0/(t12186*t12186*t12186*t12186*t12186);
                v_rho_b_gamma_aa[Q] += scale * (A*B*t12207*t12190*t12191*1.0/pow(t12186,1.4E1/3.0)*(-1.1E1/3.0)+A*B*t12190*t12191*1.0/pow(t12186,1.1E1/3.0)*(rho_b*2.0+rho_a*t12199-rho_a*rho_b*(C*t12200*(1.0/9.0)+rho_a*t12194*(C*t12200*(1.0/2.7E1)+Dd*t12200*t12190*(1.0/2.7E1)-t12201*t12202*t12203*(1.0/2.7E1))+rho_a*1.0/(t12186*t12186)*t12197+Dd*t12200*t12190*(1.0/9.0)-t12201*t12202*t12203*(1.0/9.0)))+A*B*C*t12207*t12190*t12208*t12191*(1.0/3.0)+A*B*Dd*t12203*t12207*t12208*t12191*(1.0/3.0));
            }
            
            // v_rho_b_gamma_ab
            if (deriv >= 2) {
                double t12210 = rho_a+rho_b;
                double t12211 = 1.0/pow(t12210,1.0/3.0);
                double t12212 = Dd*t12211;
                double t12213 = t12212+1.0;
                double t12214 = 1.0/t12213;
                double t12220 = C*t12211;
                double t12215 = exp(-t12220);
                double t12216 = C*t12211*(7.0/9.0);
                double t12217 = Dd*t12211*t12214*(7.0/9.0);
                double t12218 = t12216+t12217-4.7E1/9.0;
                double t12219 = 1.0/pow(t12210,4.0/3.0);
                double t12221 = t12210*t12210;
                double t12222 = t12221*(4.0/3.0);
                double t12223 = rho_a*rho_b*t12218;
                double t12224 = t12222+t12223;
                double t12225 = 1.0/(t12210*t12210*t12210*t12210*t12210);
                double t12226 = 1.0/(t12213*t12213);
                v_rho_b_gamma_ab[Q] += scale * (A*B*1.0/pow(t12210,1.4E1/3.0)*t12214*t12215*t12224*(-1.1E1/3.0)+A*B*1.0/pow(t12210,1.1E1/3.0)*t12214*t12215*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_a*t12218-rho_a*rho_b*(C*t12219*(7.0/2.7E1)-(Dd*Dd)*1.0/pow(t12210,5.0/3.0)*t12226*(7.0/2.7E1)+Dd*t12214*t12219*(7.0/2.7E1)))+A*B*C*t12214*t12215*t12224*t12225*(1.0/3.0)+A*B*Dd*t12215*t12224*t12225*t12226*(1.0/3.0));
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t12228 = rho_a+rho_b;
                double t12229 = 1.0/pow(t12228,1.0/3.0);
                double t12230 = Dd*t12229;
                double t12231 = t12230+1.0;
                double t12232 = 1.0/t12231;
                double t12246 = C*t12229;
                double t12233 = exp(-t12246);
                double t12234 = C*t12229*(1.0/3.0);
                double t12235 = Dd*t12232*t12229*(1.0/3.0);
                double t12236 = 1.0/t12228;
                double t12237 = C*t12229*(1.0/9.0);
                double t12238 = Dd*t12232*t12229*(1.0/9.0);
                double t12239 = t12237+t12238-1.1E1/9.0;
                double t12240 = rho_b*t12236*t12239;
                double t12241 = t12240+t12234+t12235-1.0/9.0;
                double t12242 = 1.0/pow(t12228,4.0/3.0);
                double t12243 = Dd*Dd;
                double t12244 = 1.0/pow(t12228,5.0/3.0);
                double t12245 = 1.0/(t12231*t12231);
                double t12247 = rho_a*rho_a;
                double t12248 = rho_a*rho_b*t12241;
                double t12249 = t12247+t12248;
                double t12250 = 1.0/(t12228*t12228*t12228*t12228*t12228);
                v_rho_b_gamma_bb[Q] += scale * (A*B*t12232*t12233*1.0/pow(t12228,1.4E1/3.0)*t12249*(-1.1E1/3.0)+A*B*t12232*t12233*1.0/pow(t12228,1.1E1/3.0)*(rho_a*t12241-rho_a*rho_b*(C*t12242*(1.0/9.0)-t12236*t12239+rho_b*t12236*(C*t12242*(1.0/2.7E1)+Dd*t12232*t12242*(1.0/2.7E1)-t12243*t12244*t12245*(1.0/2.7E1))+rho_b*1.0/(t12228*t12228)*t12239+Dd*t12232*t12242*(1.0/9.0)-t12243*t12244*t12245*(1.0/9.0)))+A*B*C*t12232*t12250*t12233*t12249*(1.0/3.0)+A*B*Dd*t12250*t12233*t12245*t12249*(1.0/3.0));
            }
            
        }
    }
}

}
