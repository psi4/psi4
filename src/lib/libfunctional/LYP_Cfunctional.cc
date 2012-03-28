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
                double t19853 = rho_a+rho_b;
                double t19854 = 1.0/pow(t19853,1.0/3.0);
                double t19855 = Dd*t19854;
                double t19856 = t19855+1.0;
                double t19857 = 1.0/t19856;
                double t19858 = t19853*t19853;
                double t19859 = t19858*(2.0/3.0);
                double t19860 = gamma_ab*2.0;
                double t19861 = gamma_aa+gamma_bb+t19860;
                double t19862 = 1.0/t19853;
                v[Q] += scale * A*rho_a*rho_b*t19862*t19857*-4.0-A*B*1.0/pow(t19853,1.1E1/3.0)*t19857*exp(-C*t19854)*(t19861*t19858*(-2.0/3.0)+gamma_aa*(t19859-rho_b*rho_b)+gamma_bb*(t19859-rho_a*rho_a)+rho_a*rho_b*((gamma_aa+gamma_bb)*(C*t19854*(1.0/1.8E1)+Dd*t19854*t19857*(1.0/1.8E1)-5.0/2.0)+CFext*(pow(rho_a,8.0/3.0)+pow(rho_b,8.0/3.0))-t19861*(C*t19854*(7.0/1.8E1)+Dd*t19854*t19857*(7.0/1.8E1)-4.7E1/1.8E1)-t19862*(gamma_aa*rho_a+gamma_bb*rho_b)*(C*t19854*(1.0/9.0)+Dd*t19854*t19857*(1.0/9.0)-1.1E1/9.0)));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t19864 = rho_a+rho_b;
                double t19865 = 1.0/pow(t19864,1.0/3.0);
                double t19866 = Dd*t19865;
                double t19867 = t19866+1.0;
                double t19868 = 1.0/t19867;
                double t19869 = t19864*t19864;
                double t19870 = t19869*(2.0/3.0);
                double t19871 = gamma_ab*2.0;
                double t19872 = gamma_aa+gamma_bb+t19871;
                double t19873 = 1.0/t19864;
                double t19900 = C*t19865;
                double t19874 = exp(-t19900);
                double t19875 = C*t19865*(7.0/1.8E1);
                double t19876 = Dd*t19865*t19868*(7.0/1.8E1);
                double t19877 = t19875+t19876-4.7E1/1.8E1;
                double t19878 = t19872*t19877;
                double t19879 = gamma_aa+gamma_bb;
                double t19880 = C*t19865*(1.0/1.8E1);
                double t19881 = Dd*t19865*t19868*(1.0/1.8E1);
                double t19882 = t19880+t19881-5.0/2.0;
                double t19883 = pow(rho_a,8.0/3.0);
                double t19884 = pow(rho_b,8.0/3.0);
                double t19885 = t19883+t19884;
                double t19886 = gamma_aa*rho_a;
                double t19887 = gamma_bb*rho_b;
                double t19888 = t19886+t19887;
                double t19889 = C*t19865*(1.0/9.0);
                double t19890 = Dd*t19865*t19868*(1.0/9.0);
                double t19891 = t19890+t19889-1.1E1/9.0;
                double t19892 = t19873*t19891*t19888;
                double t19908 = t19882*t19879;
                double t19909 = CFext*t19885;
                double t19893 = -t19908-t19909+t19892+t19878;
                double t19894 = rho_b*(4.0/3.0);
                double t19895 = 1.0/pow(t19864,4.0/3.0);
                double t19896 = 1.0/(t19867*t19867);
                double t19897 = Dd*Dd;
                double t19898 = 1.0/pow(t19864,5.0/3.0);
                double t19899 = 1.0/(t19864*t19864);
                double t19901 = rho_b*rho_b;
                double t19902 = t19901-t19870;
                double t19903 = gamma_aa*t19902;
                double t19904 = rho_a*rho_a;
                double t19905 = t19904-t19870;
                double t19906 = gamma_bb*t19905;
                double t19907 = t19872*t19869*(2.0/3.0);
                double t19910 = rho_a*rho_b*t19893;
                double t19911 = 1.0/(t19864*t19864*t19864*t19864*t19864);
                v_rho_a[Q] += scale * A*rho_b*t19873*t19868*-4.0+A*rho_a*rho_b*t19868*t19899*4.0-A*Dd*rho_a*rho_b*1.0/pow(t19864,7.0/3.0)*t19896*(4.0/3.0)-A*B*1.0/pow(t19864,1.4E1/3.0)*t19874*t19868*(t19910+t19903+t19906+t19907)*(1.1E1/3.0)+A*B*1.0/pow(t19864,1.1E1/3.0)*t19874*t19868*(rho_b*t19893-gamma_aa*(rho_a*(4.0/3.0)+t19894)+gamma_bb*(rho_a*(2.0/3.0)-t19894)+t19872*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_a,5.0/3.0)*(8.0/3.0)-t19879*(C*t19895*(1.0/5.4E1)+Dd*t19868*t19895*(1.0/5.4E1)-t19896*t19897*t19898*(1.0/5.4E1))+t19872*(C*t19895*(7.0/5.4E1)+Dd*t19868*t19895*(7.0/5.4E1)-t19896*t19897*t19898*(7.0/5.4E1))+t19873*t19888*(C*t19895*(1.0/2.7E1)+Dd*t19868*t19895*(1.0/2.7E1)-t19896*t19897*t19898*(1.0/2.7E1))-gamma_aa*t19873*t19891+t19891*t19888*t19899))+A*B*C*t19911*t19874*t19868*(t19910+t19903+t19906+t19907)*(1.0/3.0)+A*B*Dd*t19911*t19874*t19896*(t19910+t19903+t19906+t19907)*(1.0/3.0);
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t19913 = rho_a+rho_b;
                double t19914 = 1.0/pow(t19913,1.0/3.0);
                double t19915 = Dd*t19914;
                double t19916 = t19915+1.0;
                double t19917 = 1.0/t19916;
                double t19918 = t19913*t19913;
                double t19919 = t19918*(2.0/3.0);
                double t19920 = gamma_ab*2.0;
                double t19921 = gamma_aa+gamma_bb+t19920;
                double t19922 = 1.0/t19913;
                double t19949 = C*t19914;
                double t19923 = exp(-t19949);
                double t19924 = C*t19914*(7.0/1.8E1);
                double t19925 = Dd*t19914*t19917*(7.0/1.8E1);
                double t19926 = t19924+t19925-4.7E1/1.8E1;
                double t19927 = t19921*t19926;
                double t19928 = gamma_aa+gamma_bb;
                double t19929 = C*t19914*(1.0/1.8E1);
                double t19930 = Dd*t19914*t19917*(1.0/1.8E1);
                double t19931 = t19930+t19929-5.0/2.0;
                double t19932 = pow(rho_a,8.0/3.0);
                double t19933 = pow(rho_b,8.0/3.0);
                double t19934 = t19932+t19933;
                double t19935 = gamma_aa*rho_a;
                double t19936 = gamma_bb*rho_b;
                double t19937 = t19935+t19936;
                double t19938 = C*t19914*(1.0/9.0);
                double t19939 = Dd*t19914*t19917*(1.0/9.0);
                double t19940 = t19938+t19939-1.1E1/9.0;
                double t19941 = t19922*t19940*t19937;
                double t19957 = t19931*t19928;
                double t19958 = CFext*t19934;
                double t19942 = t19941+t19927-t19957-t19958;
                double t19943 = rho_a*(4.0/3.0);
                double t19944 = 1.0/pow(t19913,4.0/3.0);
                double t19945 = 1.0/(t19916*t19916);
                double t19946 = Dd*Dd;
                double t19947 = 1.0/pow(t19913,5.0/3.0);
                double t19948 = 1.0/(t19913*t19913);
                double t19950 = rho_b*rho_b;
                double t19951 = t19950-t19919;
                double t19952 = gamma_aa*t19951;
                double t19953 = rho_a*rho_a;
                double t19954 = t19953-t19919;
                double t19955 = gamma_bb*t19954;
                double t19956 = t19921*t19918*(2.0/3.0);
                double t19959 = rho_a*rho_b*t19942;
                double t19960 = 1.0/(t19913*t19913*t19913*t19913*t19913);
                v_rho_b[Q] += scale * A*rho_a*t19922*t19917*-4.0+A*rho_a*rho_b*t19917*t19948*4.0-A*Dd*rho_a*rho_b*1.0/pow(t19913,7.0/3.0)*t19945*(4.0/3.0)-A*B*1.0/pow(t19913,1.4E1/3.0)*t19923*t19917*(t19952+t19955+t19956+t19959)*(1.1E1/3.0)+A*B*1.0/pow(t19913,1.1E1/3.0)*t19923*t19917*(rho_a*t19942-gamma_bb*(rho_b*(4.0/3.0)+t19943)+gamma_aa*(rho_b*(2.0/3.0)-t19943)+t19921*(rho_a*2.0+rho_b*2.0)*(2.0/3.0)-rho_a*rho_b*(CFext*pow(rho_b,5.0/3.0)*(8.0/3.0)-t19928*(C*t19944*(1.0/5.4E1)+Dd*t19917*t19944*(1.0/5.4E1)-t19945*t19946*t19947*(1.0/5.4E1))+t19921*(C*t19944*(7.0/5.4E1)+Dd*t19917*t19944*(7.0/5.4E1)-t19945*t19946*t19947*(7.0/5.4E1))+t19922*t19937*(C*t19944*(1.0/2.7E1)+Dd*t19917*t19944*(1.0/2.7E1)-t19945*t19946*t19947*(1.0/2.7E1))-gamma_bb*t19922*t19940+t19940*t19937*t19948))+A*B*C*t19923*t19960*t19917*(t19952+t19955+t19956+t19959)*(1.0/3.0)+A*B*Dd*t19923*t19960*t19945*(t19952+t19955+t19956+t19959)*(1.0/3.0);
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t19962 = rho_a+rho_b;
                double t19963 = 1.0/pow(t19962,1.0/3.0);
                double t19964 = Dd*t19963;
                double t19965 = t19964+1.0;
                double t19966 = 1.0/t19965;
                v_gamma_aa[Q] += scale * A*B*1.0/pow(t19962,1.1E1/3.0)*t19966*exp(-C*t19963)*(rho_b*rho_b+rho_a*rho_b*(C*t19963*(1.0/3.0)+Dd*t19963*t19966*(1.0/3.0)+(rho_a*(C*t19963*(1.0/9.0)+Dd*t19963*t19966*(1.0/9.0)-1.1E1/9.0))/t19962-1.0/9.0));
            }
            
            // v_gamma_ab
            if (deriv >= 1) {
                double t19968 = rho_a+rho_b;
                double t19969 = 1.0/pow(t19968,1.0/3.0);
                double t19970 = Dd*t19969;
                double t19971 = t19970+1.0;
                double t19972 = 1.0/t19971;
                v_gamma_ab[Q] += scale * A*B*t19972*1.0/pow(t19968,1.1E1/3.0)*exp(-C*t19969)*((t19968*t19968)*(4.0/3.0)+rho_a*rho_b*(C*t19969*(7.0/9.0)+Dd*t19972*t19969*(7.0/9.0)-4.7E1/9.0));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t19974 = rho_a+rho_b;
                double t19975 = 1.0/pow(t19974,1.0/3.0);
                double t19976 = Dd*t19975;
                double t19977 = t19976+1.0;
                double t19978 = 1.0/t19977;
                v_gamma_bb[Q] += scale * A*B*1.0/pow(t19974,1.1E1/3.0)*t19978*exp(-C*t19975)*(rho_a*rho_a+rho_a*rho_b*(C*t19975*(1.0/3.0)+Dd*t19975*t19978*(1.0/3.0)+(rho_b*(C*t19975*(1.0/9.0)+Dd*t19975*t19978*(1.0/9.0)-1.1E1/9.0))/t19974-1.0/9.0));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t19982 = rho_a+rho_b;
                double t19983 = 1.0/pow(t19982,1.0/3.0);
                double t19984 = Dd*t19983;
                double t19985 = t19984+1.0;
                double t19986 = 1.0/t19985;
                double t19987 = t19982*t19982;
                double t19988 = t19987*(2.0/3.0);
                double t19989 = gamma_ab*2.0;
                double t19990 = gamma_aa+gamma_bb+t19989;
                double t19991 = 1.0/(t19985*t19985);
                double t20018 = C*t19983;
                double t19992 = exp(-t20018);
                double t19993 = C*t19983*(7.0/1.8E1);
                double t19994 = Dd*t19983*t19986*(7.0/1.8E1);
                double t19995 = t19993+t19994-4.7E1/1.8E1;
                double t19996 = t19990*t19995;
                double t19997 = gamma_aa+gamma_bb;
                double t19998 = C*t19983*(1.0/1.8E1);
                double t19999 = Dd*t19983*t19986*(1.0/1.8E1);
                double t20000 = t19998+t19999-5.0/2.0;
                double t20001 = pow(rho_a,8.0/3.0);
                double t20002 = pow(rho_b,8.0/3.0);
                double t20003 = t20001+t20002;
                double t20004 = 1.0/t19982;
                double t20005 = gamma_aa*rho_a;
                double t20006 = gamma_bb*rho_b;
                double t20007 = t20005+t20006;
                double t20008 = C*t19983*(1.0/9.0);
                double t20009 = Dd*t19983*t19986*(1.0/9.0);
                double t20010 = t20008+t20009-1.1E1/9.0;
                double t20011 = t20010*t20004*t20007;
                double t20044 = t19997*t20000;
                double t20045 = CFext*t20003;
                double t20012 = t19996+t20011-t20044-t20045;
                double t20013 = rho_b*(4.0/3.0);
                double t20014 = 1.0/pow(t19982,4.0/3.0);
                double t20015 = Dd*Dd;
                double t20016 = 1.0/pow(t19982,5.0/3.0);
                double t20017 = 1.0/(t19982*t19982);
                double t20019 = 1.0/pow(t19982,1.1E1/3.0);
                double t20020 = C*t20014*(1.0/5.4E1);
                double t20021 = Dd*t19986*t20014*(1.0/5.4E1);
                double t20052 = t19991*t20015*t20016*(1.0/5.4E1);
                double t20022 = t20020+t20021-t20052;
                double t20023 = pow(rho_a,5.0/3.0);
                double t20024 = CFext*t20023*(8.0/3.0);
                double t20025 = C*t20014*(7.0/5.4E1);
                double t20026 = Dd*t19986*t20014*(7.0/5.4E1);
                double t20054 = t19991*t20015*t20016*(7.0/5.4E1);
                double t20027 = t20025+t20026-t20054;
                double t20028 = t19990*t20027;
                double t20029 = C*t20014*(1.0/2.7E1);
                double t20030 = Dd*t19986*t20014*(1.0/2.7E1);
                double t20039 = t19991*t20015*t20016*(1.0/2.7E1);
                double t20031 = t20030+t20029-t20039;
                double t20032 = t20004*t20031*t20007;
                double t20033 = t20010*t20007*t20017;
                double t20053 = t19997*t20022;
                double t20055 = gamma_aa*t20010*t20004;
                double t20034 = t20032+t20024+t20033-t20053+t20028-t20055;
                double t20035 = 1.0/pow(t19982,7.0/3.0);
                double t20036 = 1.0/(t19982*t19982*t19982);
                double t20037 = 1.0/(t19985*t19985*t19985);
                double t20038 = 1.0/pow(t19982,8.0/3.0);
                double t20040 = rho_a*2.0;
                double t20041 = rho_b*2.0;
                double t20042 = t20040+t20041;
                double t20043 = t19990*t20042*(2.0/3.0);
                double t20046 = rho_b*t20012;
                double t20047 = rho_a*(4.0/3.0);
                double t20048 = t20013+t20047;
                double t20049 = rho_a*(2.0/3.0);
                double t20050 = t20013-t20049;
                double t20051 = gamma_bb*t20050;
                double t20056 = rho_a*rho_b*t20034;
                double t20057 = gamma_aa*t20048;
                double t20058 = t20051-t20043-t20046+t20056+t20057;
                double t20059 = 1.0/(t19982*t19982*t19982*t19982*t19982);
                double t20060 = rho_b*rho_b;
                double t20061 = t19988-t20060;
                double t20062 = gamma_aa*t20061;
                double t20063 = rho_a*rho_a;
                double t20064 = t19988-t20063;
                double t20065 = gamma_bb*t20064;
                double t20068 = t19990*t19987*(2.0/3.0);
                double t20069 = rho_a*rho_b*t20012;
                double t20066 = t20062+t20065-t20068-t20069;
                double t20067 = 1.0/(t19982*t19982*t19982*t19982*t19982*t19982);
                double t20070 = 1.0/pow(t19982,1.9E1/3.0);
                v_rho_a_rho_a[Q] += scale * A*rho_b*t19986*t20017*8.0-A*Dd*rho_b*t19991*t20035*(8.0/3.0)-A*rho_a*rho_b*t19986*t20036*8.0+A*Dd*rho_a*rho_b*1.0/pow(t19982,1.0E1/3.0)*t19991*(4.0E1/9.0)+A*B*1.0/pow(t19982,1.4E1/3.0)*t19992*t19986*t20058*(2.2E1/3.0)-A*B*1.0/pow(t19982,1.7E1/3.0)*t19992*t19986*t20066*(1.54E2/9.0)+A*B*t19992*t19986*t20019*(gamma_ab*(8.0/3.0)+gamma_bb*2.0-rho_b*t20034*2.0-rho_a*rho_b*(CFext*pow(rho_a,2.0/3.0)*(4.0E1/9.0)+t19997*(C*t20035*(2.0/8.1E1)+Dd*t19986*t20035*(2.0/8.1E1)-t19991*t20015*t20038*(1.0/2.7E1)+Dd*t20015*t20036*t20037*(1.0/8.1E1))-t19990*(C*t20035*(1.4E1/8.1E1)+Dd*t19986*t20035*(1.4E1/8.1E1)-t19991*t20015*t20038*(7.0/2.7E1)+Dd*t20015*t20036*t20037*(7.0/8.1E1))-t20004*t20007*(C*t20035*(4.0/8.1E1)+Dd*t19986*t20035*(4.0/8.1E1)-t19991*t20015*t20038*(2.0/2.7E1)+Dd*t20015*t20036*t20037*(2.0/8.1E1))+gamma_aa*t20004*t20031*2.0+gamma_aa*t20010*t20017*2.0-t20010*t20007*t20036*2.0-t20031*t20007*t20017*2.0))-A*rho_a*rho_b*t20015*t20019*t20037*(8.0/9.0)-A*B*t19992*t20015*t20070*t20037*t20066*(2.0/9.0)-A*B*(C*C)*t19992*t19986*t20070*t20066*(1.0/9.0)+A*B*C*t19992*t19986*t20066*t20067*(2.6E1/9.0)-A*B*C*t19992*t19986*t20058*t20059*(2.0/3.0)+A*B*Dd*t19991*t19992*t20066*t20067*(2.6E1/9.0)-A*B*Dd*t19991*t19992*t20058*t20059*(2.0/3.0)-A*B*C*Dd*t19991*t19992*t20070*t20066*(2.0/9.0);
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t20072 = rho_a+rho_b;
                double t20073 = 1.0/pow(t20072,1.0/3.0);
                double t20074 = Dd*t20073;
                double t20075 = t20074+1.0;
                double t20076 = 1.0/t20075;
                double t20077 = 1.0/(t20072*t20072);
                double t20078 = 1.0/pow(t20072,7.0/3.0);
                double t20079 = 1.0/(t20075*t20075);
                double t20080 = t20072*t20072;
                double t20081 = t20080*(2.0/3.0);
                double t20082 = gamma_ab*2.0;
                double t20083 = gamma_aa+gamma_bb+t20082;
                double t20084 = 1.0/t20072;
                double t20123 = C*t20073;
                double t20085 = exp(-t20123);
                double t20086 = C*t20073*(7.0/1.8E1);
                double t20087 = Dd*t20073*t20076*(7.0/1.8E1);
                double t20088 = t20086+t20087-4.7E1/1.8E1;
                double t20089 = t20083*t20088;
                double t20090 = gamma_aa+gamma_bb;
                double t20091 = 1.0/pow(t20072,4.0/3.0);
                double t20092 = Dd*Dd;
                double t20093 = 1.0/pow(t20072,5.0/3.0);
                double t20094 = C*t20073*(1.0/9.0);
                double t20095 = Dd*t20073*t20076*(1.0/9.0);
                double t20096 = t20094+t20095-1.1E1/9.0;
                double t20097 = gamma_aa*rho_a;
                double t20098 = gamma_bb*rho_b;
                double t20099 = t20097+t20098;
                double t20100 = C*t20091*(1.0/5.4E1);
                double t20101 = Dd*t20091*t20076*(1.0/5.4E1);
                double t20127 = t20092*t20093*t20079*(1.0/5.4E1);
                double t20102 = t20100+t20101-t20127;
                double t20103 = C*t20091*(7.0/5.4E1);
                double t20104 = Dd*t20091*t20076*(7.0/5.4E1);
                double t20130 = t20092*t20093*t20079*(7.0/5.4E1);
                double t20105 = t20103-t20130+t20104;
                double t20106 = t20105*t20083;
                double t20107 = C*t20091*(1.0/2.7E1);
                double t20108 = Dd*t20091*t20076*(1.0/2.7E1);
                double t20121 = t20092*t20093*t20079*(1.0/2.7E1);
                double t20109 = -t20121+t20107+t20108;
                double t20110 = t20109*t20084*t20099;
                double t20111 = t20077*t20096*t20099;
                double t20112 = C*t20073*(1.0/1.8E1);
                double t20113 = Dd*t20073*t20076*(1.0/1.8E1);
                double t20114 = t20112+t20113-5.0/2.0;
                double t20115 = pow(rho_a,8.0/3.0);
                double t20116 = pow(rho_b,8.0/3.0);
                double t20117 = t20115+t20116;
                double t20118 = 1.0/(t20072*t20072*t20072);
                double t20119 = 1.0/(t20075*t20075*t20075);
                double t20120 = 1.0/pow(t20072,8.0/3.0);
                double t20122 = t20084*t20096*t20099;
                double t20124 = t20114*t20090;
                double t20125 = CFext*t20117;
                double t20126 = rho_b*(4.0/3.0);
                double t20128 = pow(rho_a,5.0/3.0);
                double t20129 = CFext*t20128*(8.0/3.0);
                double t20140 = t20102*t20090;
                double t20148 = gamma_aa*t20084*t20096;
                double t20131 = t20110+t20111-t20140+t20106+t20129-t20148;
                double t20132 = 1.0/pow(t20072,1.4E1/3.0);
                double t20133 = rho_a*2.0;
                double t20134 = rho_b*2.0;
                double t20135 = t20133+t20134;
                double t20136 = t20135*t20083*(2.0/3.0);
                double t20137 = t20122-t20124-t20125+t20089;
                double t20138 = rho_a*(4.0/3.0);
                double t20139 = t20126+t20138;
                double t20141 = pow(rho_b,5.0/3.0);
                double t20142 = CFext*t20141*(8.0/3.0);
                double t20143 = 1.0/pow(t20072,1.1E1/3.0);
                double t20144 = rho_b*t20137;
                double t20145 = rho_a*(2.0/3.0);
                double t20146 = t20126-t20145;
                double t20147 = gamma_bb*t20146;
                double t20149 = rho_a*rho_b*t20131;
                double t20150 = gamma_aa*t20139;
                double t20151 = t20150-t20144-t20136+t20147+t20149;
                double t20152 = 1.0/(t20072*t20072*t20072*t20072*t20072);
                double t20153 = rho_a*t20137;
                double t20154 = rho_b*(2.0/3.0);
                double t20155 = t20154-t20138;
                double t20156 = gamma_aa*t20155;
                double t20158 = gamma_bb*t20084*t20096;
                double t20157 = t20110+t20111-t20140+t20106+t20142-t20158;
                double t20159 = rho_b*rho_b;
                double t20160 = t20081-t20159;
                double t20161 = gamma_aa*t20160;
                double t20162 = rho_a*rho_a;
                double t20163 = t20081-t20162;
                double t20164 = gamma_bb*t20163;
                double t20165 = 1.0/(t20072*t20072*t20072*t20072*t20072*t20072);
                double t20167 = t20080*t20083*(2.0/3.0);
                double t20168 = rho_a*rho_b*t20137;
                double t20166 = t20161+t20164-t20167-t20168;
                double t20169 = 1.0/pow(t20072,1.9E1/3.0);
                v_rho_a_rho_b[Q] += scale * A*t20084*t20076*-4.0+A*rho_a*t20076*t20077*4.0+A*rho_b*t20076*t20077*4.0-A*Dd*rho_a*t20078*t20079*(4.0/3.0)-A*Dd*rho_b*t20078*t20079*(4.0/3.0)-A*rho_a*rho_b*t20118*t20076*8.0+A*Dd*rho_a*rho_b*1.0/pow(t20072,1.0E1/3.0)*t20079*(4.0E1/9.0)-A*B*1.0/pow(t20072,1.7E1/3.0)*t20076*t20085*(t20161+t20164-t20080*t20083*(2.0/3.0)-rho_a*rho_b*(t20122+t20089-CFext*t20117-t20114*t20090))*(1.54E2/9.0)-A*B*t20143*t20076*t20085*(gamma_ab*(-8.0/3.0)-t20122+t20124+t20125-t20089+rho_a*t20131+rho_b*(t20110+t20111+t20106+t20142-t20102*t20090-gamma_bb*t20084*t20096)+rho_a*rho_b*(t20090*(C*t20078*(2.0/8.1E1)+Dd*t20076*t20078*(2.0/8.1E1)-t20120*t20092*t20079*(1.0/2.7E1)+Dd*t20118*t20092*t20119*(1.0/8.1E1))-t20083*(C*t20078*(1.4E1/8.1E1)+Dd*t20076*t20078*(1.4E1/8.1E1)-t20120*t20092*t20079*(7.0/2.7E1)+Dd*t20118*t20092*t20119*(7.0/8.1E1))-t20084*t20099*(C*t20078*(4.0/8.1E1)+Dd*t20076*t20078*(4.0/8.1E1)-t20120*t20092*t20079*(2.0/2.7E1)+Dd*t20118*t20092*t20119*(2.0/8.1E1))+gamma_aa*t20109*t20084+gamma_bb*t20109*t20084+gamma_aa*t20077*t20096+gamma_bb*t20077*t20096-t20109*t20077*t20099*2.0-t20118*t20096*t20099*2.0))+A*B*t20132*t20151*t20076*t20085*(1.1E1/3.0)-A*rho_a*rho_b*t20143*t20092*t20119*(8.0/9.0)-A*B*t20132*t20076*t20085*(t20153+t20136+t20156-gamma_bb*t20139-rho_a*rho_b*t20157)*(1.1E1/3.0)-A*B*t20092*t20119*t20085*t20166*t20169*(2.0/9.0)+A*B*C*t20152*t20076*t20085*(t20153+t20136+t20156-gamma_bb*t20139-rho_a*rho_b*t20157)*(1.0/3.0)+A*B*Dd*t20152*t20085*t20079*(t20153+t20136+t20156-gamma_bb*t20139-rho_a*rho_b*t20157)*(1.0/3.0)-A*B*(C*C)*t20076*t20085*t20166*t20169*(1.0/9.0)-A*B*C*t20151*t20152*t20076*t20085*(1.0/3.0)+A*B*C*t20165*t20076*t20085*t20166*(2.6E1/9.0)-A*B*Dd*t20151*t20152*t20085*t20079*(1.0/3.0)+A*B*Dd*t20165*t20085*t20166*t20079*(2.6E1/9.0)-A*B*C*Dd*t20085*t20166*t20079*t20169*(2.0/9.0);
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t20171 = rho_a+rho_b;
                double t20172 = 1.0/pow(t20171,1.0/3.0);
                double t20173 = Dd*t20172;
                double t20174 = t20173+1.0;
                double t20175 = 1.0/t20174;
                double t20176 = t20171*t20171;
                double t20177 = t20176*(2.0/3.0);
                double t20178 = gamma_ab*2.0;
                double t20179 = gamma_aa+gamma_bb+t20178;
                double t20180 = 1.0/(t20174*t20174);
                double t20207 = C*t20172;
                double t20181 = exp(-t20207);
                double t20182 = C*t20172*(7.0/1.8E1);
                double t20183 = Dd*t20172*t20175*(7.0/1.8E1);
                double t20184 = t20182+t20183-4.7E1/1.8E1;
                double t20185 = t20184*t20179;
                double t20186 = gamma_aa+gamma_bb;
                double t20187 = C*t20172*(1.0/1.8E1);
                double t20188 = Dd*t20172*t20175*(1.0/1.8E1);
                double t20189 = t20187+t20188-5.0/2.0;
                double t20190 = pow(rho_a,8.0/3.0);
                double t20191 = pow(rho_b,8.0/3.0);
                double t20192 = t20190+t20191;
                double t20193 = 1.0/t20171;
                double t20194 = gamma_aa*rho_a;
                double t20195 = gamma_bb*rho_b;
                double t20196 = t20194+t20195;
                double t20197 = C*t20172*(1.0/9.0);
                double t20198 = Dd*t20172*t20175*(1.0/9.0);
                double t20199 = t20197+t20198-1.1E1/9.0;
                double t20200 = t20193*t20196*t20199;
                double t20233 = t20186*t20189;
                double t20234 = CFext*t20192;
                double t20201 = t20200-t20233-t20234+t20185;
                double t20202 = rho_a*(4.0/3.0);
                double t20203 = 1.0/pow(t20171,4.0/3.0);
                double t20204 = Dd*Dd;
                double t20205 = 1.0/pow(t20171,5.0/3.0);
                double t20206 = 1.0/(t20171*t20171);
                double t20208 = 1.0/pow(t20171,1.1E1/3.0);
                double t20209 = C*t20203*(1.0/5.4E1);
                double t20210 = Dd*t20203*t20175*(1.0/5.4E1);
                double t20242 = t20204*t20205*t20180*(1.0/5.4E1);
                double t20211 = t20210-t20242+t20209;
                double t20212 = pow(rho_b,5.0/3.0);
                double t20213 = CFext*t20212*(8.0/3.0);
                double t20214 = C*t20203*(7.0/5.4E1);
                double t20215 = Dd*t20203*t20175*(7.0/5.4E1);
                double t20244 = t20204*t20205*t20180*(7.0/5.4E1);
                double t20216 = t20214+t20215-t20244;
                double t20217 = t20216*t20179;
                double t20218 = C*t20203*(1.0/2.7E1);
                double t20219 = Dd*t20203*t20175*(1.0/2.7E1);
                double t20228 = t20204*t20205*t20180*(1.0/2.7E1);
                double t20220 = t20218+t20219-t20228;
                double t20221 = t20220*t20193*t20196;
                double t20222 = t20206*t20196*t20199;
                double t20243 = t20211*t20186;
                double t20245 = gamma_bb*t20193*t20199;
                double t20223 = t20221+t20213+t20222-t20243+t20217-t20245;
                double t20224 = 1.0/pow(t20171,7.0/3.0);
                double t20225 = 1.0/(t20171*t20171*t20171);
                double t20226 = 1.0/(t20174*t20174*t20174);
                double t20227 = 1.0/pow(t20171,8.0/3.0);
                double t20229 = rho_a*2.0;
                double t20230 = rho_b*2.0;
                double t20231 = t20230+t20229;
                double t20232 = t20231*t20179*(2.0/3.0);
                double t20235 = rho_a*t20201;
                double t20236 = rho_b*(2.0/3.0);
                double t20237 = t20202-t20236;
                double t20238 = gamma_aa*t20237;
                double t20239 = rho_b*(4.0/3.0);
                double t20240 = t20202+t20239;
                double t20241 = gamma_bb*t20240;
                double t20246 = rho_a*rho_b*t20223;
                double t20247 = -t20232+t20241-t20235+t20246+t20238;
                double t20248 = 1.0/(t20171*t20171*t20171*t20171*t20171);
                double t20249 = rho_b*rho_b;
                double t20250 = t20177-t20249;
                double t20251 = gamma_aa*t20250;
                double t20252 = rho_a*rho_a;
                double t20253 = t20252-t20177;
                double t20254 = gamma_bb*t20253;
                double t20255 = t20176*t20179*(2.0/3.0);
                double t20256 = rho_a*rho_b*t20201;
                double t20257 = -t20251+t20254+t20255+t20256;
                double t20258 = 1.0/(t20171*t20171*t20171*t20171*t20171*t20171);
                double t20259 = 1.0/pow(t20171,1.9E1/3.0);
                v_rho_b_rho_b[Q] += scale * A*rho_a*t20206*t20175*8.0-A*Dd*rho_a*t20224*t20180*(8.0/3.0)-A*rho_a*rho_b*t20225*t20175*8.0+A*Dd*rho_a*rho_b*1.0/pow(t20171,1.0E1/3.0)*t20180*(4.0E1/9.0)+A*B*1.0/pow(t20171,1.4E1/3.0)*t20181*t20175*t20247*(2.2E1/3.0)+A*B*1.0/pow(t20171,1.7E1/3.0)*t20181*t20175*t20257*(1.54E2/9.0)+A*B*t20181*t20208*t20175*(gamma_aa*2.0+gamma_ab*(8.0/3.0)-rho_a*t20223*2.0-rho_a*rho_b*(CFext*pow(rho_b,2.0/3.0)*(4.0E1/9.0)+t20186*(C*t20224*(2.0/8.1E1)+Dd*t20224*t20175*(2.0/8.1E1)-t20204*t20180*t20227*(1.0/2.7E1)+Dd*t20204*t20225*t20226*(1.0/8.1E1))-t20179*(C*t20224*(1.4E1/8.1E1)+Dd*t20224*t20175*(1.4E1/8.1E1)-t20204*t20180*t20227*(7.0/2.7E1)+Dd*t20204*t20225*t20226*(7.0/8.1E1))-t20193*t20196*(C*t20224*(4.0/8.1E1)+Dd*t20224*t20175*(4.0/8.1E1)-t20204*t20180*t20227*(2.0/2.7E1)+Dd*t20204*t20225*t20226*(2.0/8.1E1))+gamma_bb*t20220*t20193*2.0+gamma_bb*t20206*t20199*2.0-t20220*t20206*t20196*2.0-t20225*t20196*t20199*2.0))-A*rho_a*rho_b*t20204*t20208*t20226*(8.0/9.0)+A*B*t20204*t20181*t20226*t20259*(-t20251+t20254+t20255+t20256)*(2.0/9.0)+A*B*(C*C)*t20181*t20175*t20259*(-t20251+t20254+t20255+t20256)*(1.0/9.0)-A*B*C*t20181*t20175*t20247*t20248*(2.0/3.0)-A*B*C*t20181*t20175*t20257*t20258*(2.6E1/9.0)-A*B*Dd*t20180*t20181*t20247*t20248*(2.0/3.0)-A*B*Dd*t20180*t20181*t20257*t20258*(2.6E1/9.0)+A*B*C*Dd*t20180*t20181*t20259*(-t20251+t20254+t20255+t20256)*(2.0/9.0);
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t20261 = rho_a+rho_b;
                double t20262 = 1.0/pow(t20261,1.0/3.0);
                double t20263 = Dd*t20262;
                double t20264 = t20263+1.0;
                double t20265 = 1.0/t20264;
                double t20279 = C*t20262;
                double t20266 = exp(-t20279);
                double t20267 = C*t20262*(1.0/3.0);
                double t20268 = Dd*t20262*t20265*(1.0/3.0);
                double t20269 = 1.0/t20261;
                double t20270 = C*t20262*(1.0/9.0);
                double t20271 = Dd*t20262*t20265*(1.0/9.0);
                double t20272 = t20270+t20271-1.1E1/9.0;
                double t20273 = rho_a*t20272*t20269;
                double t20274 = t20273+t20267+t20268-1.0/9.0;
                double t20275 = 1.0/pow(t20261,4.0/3.0);
                double t20276 = Dd*Dd;
                double t20277 = 1.0/pow(t20261,5.0/3.0);
                double t20278 = 1.0/(t20264*t20264);
                double t20280 = rho_b*rho_b;
                double t20281 = rho_a*rho_b*t20274;
                double t20282 = t20280+t20281;
                double t20283 = 1.0/(t20261*t20261*t20261*t20261*t20261);
                v_rho_a_gamma_aa[Q] += scale * A*B*1.0/pow(t20261,1.4E1/3.0)*t20282*t20265*t20266*(-1.1E1/3.0)+A*B*1.0/pow(t20261,1.1E1/3.0)*t20265*t20266*(rho_b*t20274-rho_a*rho_b*(C*t20275*(1.0/9.0)-t20272*t20269+rho_a*t20269*(C*t20275*(1.0/2.7E1)+Dd*t20265*t20275*(1.0/2.7E1)-t20276*t20277*t20278*(1.0/2.7E1))+rho_a*1.0/(t20261*t20261)*t20272+Dd*t20265*t20275*(1.0/9.0)-t20276*t20277*t20278*(1.0/9.0)))+A*B*C*t20282*t20265*t20283*t20266*(1.0/3.0)+A*B*Dd*t20282*t20283*t20266*t20278*(1.0/3.0);
            }
            
            // v_rho_a_gamma_ab
            if (deriv >= 2) {
                double t20285 = rho_a+rho_b;
                double t20286 = 1.0/pow(t20285,1.0/3.0);
                double t20287 = Dd*t20286;
                double t20288 = t20287+1.0;
                double t20289 = 1.0/t20288;
                double t20295 = C*t20286;
                double t20290 = exp(-t20295);
                double t20291 = C*t20286*(7.0/9.0);
                double t20292 = Dd*t20286*t20289*(7.0/9.0);
                double t20293 = t20291+t20292-4.7E1/9.0;
                double t20294 = 1.0/pow(t20285,4.0/3.0);
                double t20296 = t20285*t20285;
                double t20297 = t20296*(4.0/3.0);
                double t20298 = rho_a*rho_b*t20293;
                double t20299 = t20297+t20298;
                double t20300 = 1.0/(t20285*t20285*t20285*t20285*t20285);
                double t20301 = 1.0/(t20288*t20288);
                v_rho_a_gamma_ab[Q] += scale * A*B*t20290*1.0/pow(t20285,1.4E1/3.0)*t20289*t20299*(-1.1E1/3.0)+A*B*t20290*1.0/pow(t20285,1.1E1/3.0)*t20289*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_b*t20293-rho_a*rho_b*(C*t20294*(7.0/2.7E1)-(Dd*Dd)*t20301*1.0/pow(t20285,5.0/3.0)*(7.0/2.7E1)+Dd*t20294*t20289*(7.0/2.7E1)))+A*B*C*t20300*t20290*t20289*t20299*(1.0/3.0)+A*B*Dd*t20300*t20301*t20290*t20299*(1.0/3.0);
            }
            
            // v_rho_a_gamma_bb
            if (deriv >= 2) {
                double t20303 = rho_a+rho_b;
                double t20304 = 1.0/pow(t20303,1.0/3.0);
                double t20305 = Dd*t20304;
                double t20306 = t20305+1.0;
                double t20307 = 1.0/t20306;
                double t20321 = C*t20304;
                double t20308 = exp(-t20321);
                double t20309 = C*t20304*(1.0/3.0);
                double t20310 = Dd*t20304*t20307*(1.0/3.0);
                double t20311 = 1.0/t20303;
                double t20312 = C*t20304*(1.0/9.0);
                double t20313 = Dd*t20304*t20307*(1.0/9.0);
                double t20314 = t20312+t20313-1.1E1/9.0;
                double t20315 = rho_b*t20311*t20314;
                double t20316 = t20310+t20315+t20309-1.0/9.0;
                double t20317 = 1.0/pow(t20303,4.0/3.0);
                double t20318 = Dd*Dd;
                double t20319 = 1.0/pow(t20303,5.0/3.0);
                double t20320 = 1.0/(t20306*t20306);
                double t20322 = rho_a*rho_a;
                double t20323 = rho_a*rho_b*t20316;
                double t20324 = t20322+t20323;
                double t20325 = 1.0/(t20303*t20303*t20303*t20303*t20303);
                v_rho_a_gamma_bb[Q] += scale * A*B*1.0/pow(t20303,1.4E1/3.0)*t20324*t20307*t20308*(-1.1E1/3.0)+A*B*1.0/pow(t20303,1.1E1/3.0)*t20307*t20308*(rho_a*2.0+rho_b*t20316-rho_a*rho_b*(C*t20317*(1.0/9.0)+rho_b*t20311*(C*t20317*(1.0/2.7E1)+Dd*t20307*t20317*(1.0/2.7E1)-t20320*t20318*t20319*(1.0/2.7E1))+rho_b*1.0/(t20303*t20303)*t20314+Dd*t20307*t20317*(1.0/9.0)-t20320*t20318*t20319*(1.0/9.0)))+A*B*C*t20324*t20307*t20325*t20308*(1.0/3.0)+A*B*Dd*t20320*t20324*t20325*t20308*(1.0/3.0);
            }
            
            // v_rho_b_gamma_aa
            if (deriv >= 2) {
                double t20327 = rho_a+rho_b;
                double t20328 = 1.0/pow(t20327,1.0/3.0);
                double t20329 = Dd*t20328;
                double t20330 = t20329+1.0;
                double t20331 = 1.0/t20330;
                double t20345 = C*t20328;
                double t20332 = exp(-t20345);
                double t20333 = C*t20328*(1.0/3.0);
                double t20334 = Dd*t20331*t20328*(1.0/3.0);
                double t20335 = 1.0/t20327;
                double t20336 = C*t20328*(1.0/9.0);
                double t20337 = Dd*t20331*t20328*(1.0/9.0);
                double t20338 = t20336+t20337-1.1E1/9.0;
                double t20339 = rho_a*t20335*t20338;
                double t20340 = t20333+t20334+t20339-1.0/9.0;
                double t20341 = 1.0/pow(t20327,4.0/3.0);
                double t20342 = Dd*Dd;
                double t20343 = 1.0/pow(t20327,5.0/3.0);
                double t20344 = 1.0/(t20330*t20330);
                double t20346 = rho_b*rho_b;
                double t20347 = rho_a*rho_b*t20340;
                double t20348 = t20346+t20347;
                double t20349 = 1.0/(t20327*t20327*t20327*t20327*t20327);
                v_rho_b_gamma_aa[Q] += scale * A*B*t20331*t20332*1.0/pow(t20327,1.4E1/3.0)*t20348*(-1.1E1/3.0)+A*B*t20331*t20332*1.0/pow(t20327,1.1E1/3.0)*(rho_b*2.0+rho_a*t20340-rho_a*rho_b*(C*t20341*(1.0/9.0)+rho_a*t20335*(C*t20341*(1.0/2.7E1)+Dd*t20331*t20341*(1.0/2.7E1)-t20342*t20343*t20344*(1.0/2.7E1))+rho_a*1.0/(t20327*t20327)*t20338+Dd*t20331*t20341*(1.0/9.0)-t20342*t20343*t20344*(1.0/9.0)))+A*B*C*t20331*t20332*t20348*t20349*(1.0/3.0)+A*B*Dd*t20332*t20344*t20348*t20349*(1.0/3.0);
            }
            
            // v_rho_b_gamma_ab
            if (deriv >= 2) {
                double t20351 = rho_a+rho_b;
                double t20352 = 1.0/pow(t20351,1.0/3.0);
                double t20353 = Dd*t20352;
                double t20354 = t20353+1.0;
                double t20355 = 1.0/t20354;
                double t20361 = C*t20352;
                double t20356 = exp(-t20361);
                double t20357 = C*t20352*(7.0/9.0);
                double t20358 = Dd*t20352*t20355*(7.0/9.0);
                double t20359 = t20357+t20358-4.7E1/9.0;
                double t20360 = 1.0/pow(t20351,4.0/3.0);
                double t20362 = t20351*t20351;
                double t20363 = t20362*(4.0/3.0);
                double t20364 = rho_a*rho_b*t20359;
                double t20365 = t20363+t20364;
                double t20366 = 1.0/(t20351*t20351*t20351*t20351*t20351);
                double t20367 = 1.0/(t20354*t20354);
                v_rho_b_gamma_ab[Q] += scale * A*B*1.0/pow(t20351,1.4E1/3.0)*t20355*t20356*t20365*(-1.1E1/3.0)+A*B*1.0/pow(t20351,1.1E1/3.0)*t20355*t20356*(rho_a*(8.0/3.0)+rho_b*(8.0/3.0)+rho_a*t20359-rho_a*rho_b*(C*t20360*(7.0/2.7E1)-(Dd*Dd)*1.0/pow(t20351,5.0/3.0)*t20367*(7.0/2.7E1)+Dd*t20360*t20355*(7.0/2.7E1)))+A*B*C*t20355*t20356*t20365*t20366*(1.0/3.0)+A*B*Dd*t20356*t20365*t20366*t20367*(1.0/3.0);
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t20369 = rho_a+rho_b;
                double t20370 = 1.0/pow(t20369,1.0/3.0);
                double t20371 = Dd*t20370;
                double t20372 = t20371+1.0;
                double t20373 = 1.0/t20372;
                double t20387 = C*t20370;
                double t20374 = exp(-t20387);
                double t20375 = C*t20370*(1.0/3.0);
                double t20376 = Dd*t20370*t20373*(1.0/3.0);
                double t20377 = 1.0/t20369;
                double t20378 = C*t20370*(1.0/9.0);
                double t20379 = Dd*t20370*t20373*(1.0/9.0);
                double t20380 = t20378+t20379-1.1E1/9.0;
                double t20381 = rho_b*t20380*t20377;
                double t20382 = t20381+t20375+t20376-1.0/9.0;
                double t20383 = 1.0/pow(t20369,4.0/3.0);
                double t20384 = Dd*Dd;
                double t20385 = 1.0/pow(t20369,5.0/3.0);
                double t20386 = 1.0/(t20372*t20372);
                double t20388 = rho_a*rho_a;
                double t20389 = rho_a*rho_b*t20382;
                double t20390 = t20388+t20389;
                double t20391 = 1.0/(t20369*t20369*t20369*t20369*t20369);
                v_rho_b_gamma_bb[Q] += scale * A*B*t20390*t20373*t20374*1.0/pow(t20369,1.4E1/3.0)*(-1.1E1/3.0)+A*B*t20373*t20374*1.0/pow(t20369,1.1E1/3.0)*(rho_a*t20382-rho_a*rho_b*(C*t20383*(1.0/9.0)-t20380*t20377+rho_b*t20377*(C*t20383*(1.0/2.7E1)+Dd*t20373*t20383*(1.0/2.7E1)-t20384*t20385*t20386*(1.0/2.7E1))+rho_b*t20380*1.0/(t20369*t20369)+Dd*t20373*t20383*(1.0/9.0)-t20384*t20385*t20386*(1.0/9.0)))+A*B*C*t20390*t20373*t20391*t20374*(1.0/3.0)+A*B*Dd*t20390*t20391*t20374*t20386*(1.0/3.0);
            }
            
        }
    }
}

}
