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
                double t1970 = 1.0/pow(rho_b,8.0/3.0);
                double t1971 = log(gamma_bb*t1970+sqrt((gamma_bb*gamma_bb)*(t1970*t1970)+1.0));
                double t1973 = d2*d2;
                double t1974 = gamma_bb+t1973;
                double t1975 = 1.0/t1974;
                double t1976 = d1*gamma_bb*t1975;
                double t1972 = d0+t1976;
                v[Q] += scale * (-c*pow(rho_b,4.0/3.0)*((gamma_bb*t1970*t1972*1.0/sqrt(gamma_bb*t1970*(t1971*t1971)*(t1972*t1972)*9.0+1.0))/c-1.0));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t1979 = 1.0/pow(rho_b,8.0/3.0);
                double t1987 = gamma_bb*t1979;
                double t1980 = log(t1987+sqrt(t1987*t1987+1.0));
                double t1982 = d2*d2;
                double t1983 = gamma_bb+t1982;
                double t1984 = 1.0/t1983;
                double t1985 = d1*gamma_bb*t1984;
                double t1981 = d0+t1985;
                double t1986 = 1.0/c;
                double t1988 = t1980*t1980;
                double t1989 = t1981*t1981;
                double t1990 = gamma_bb*t1979*t1988*t1989*9.0;
                double t1991 = t1990+1.0;
                double t1992 = 1.0/sqrt(t1991);
                double t1993 = 1.0/pow(rho_b,1.1E1/3.0);
                double t1994 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * (c*pow(rho_b,1.0/3.0)*(gamma_bb*t1981*t1992*t1986*t1979-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t1981*t1992*t1993*t1986*(8.0/3.0)-gamma_bb*t1981*1.0/pow(t1991,3.0/2.0)*t1986*t1979*(gamma_bb*t1993*t1988*t1989*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t1980*t1994*t1989*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t1994+1.0)*4.8E1)*(1.0/2.0)));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t1998 = 1.0/pow(rho_b,8.0/3.0);
                double t2006 = gamma_bb*t1998;
                double t1999 = log(t2006+sqrt(t2006*t2006+1.0));
                double t2001 = d2*d2;
                double t2002 = gamma_bb+t2001;
                double t2003 = 1.0/t2002;
                double t2004 = d1*gamma_bb*t2003;
                double t2000 = d0+t2004;
                double t2005 = 1.0/c;
                double t2007 = t1999*t1999;
                double t2008 = t2000*t2000;
                double t2009 = gamma_bb*t1998*t2007*t2008*9.0;
                double t2010 = t2009+1.0;
                double t2011 = 1.0/sqrt(t2010);
                double t2012 = 1.0/pow(rho_b,1.6E1/3.0);
                double t2013 = d1*t2003;
                double t2014 = 1.0/(t2002*t2002);
                double t2015 = t2013-d1*gamma_bb*t2014;
                v_gamma_bb[Q] += scale * (-c*pow(rho_b,4.0/3.0)*(t1998*t2000*t2011*t2005+gamma_bb*t1998*t2011*t2005*t2015-gamma_bb*t1998*t2000*1.0/pow(t2010,3.0/2.0)*t2005*(t1998*t2007*t2008*9.0+gamma_bb*t1999*t2012*t2008*1.0/sqrt((gamma_bb*gamma_bb)*t2012+1.0)*1.8E1+gamma_bb*t1998*t2000*t2015*t2007*1.8E1)*(1.0/2.0)));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t2021 = 1.0/pow(rho_b,8.0/3.0);
                double t2029 = gamma_bb*t2021;
                double t2022 = log(t2029+sqrt(t2029*t2029+1.0));
                double t2024 = d2*d2;
                double t2025 = gamma_bb+t2024;
                double t2026 = 1.0/t2025;
                double t2027 = d1*gamma_bb*t2026;
                double t2023 = d0+t2027;
                double t2028 = 1.0/c;
                double t2030 = t2022*t2022;
                double t2031 = t2023*t2023;
                double t2032 = gamma_bb*t2021*t2030*t2031*9.0;
                double t2033 = t2032+1.0;
                double t2034 = 1.0/sqrt(t2033);
                double t2035 = 1.0/pow(rho_b,1.1E1/3.0);
                double t2036 = gamma_bb*gamma_bb;
                double t2037 = 1.0/pow(t2033,3.0/2.0);
                double t2038 = 1.0/pow(rho_b,1.4E1/3.0);
                double t2039 = 1.0/pow(rho_b,1.6E1/3.0);
                double t2040 = t2036*t2039;
                double t2041 = t2040+1.0;
                double t2042 = 1.0/sqrt(t2041);
                double t2043 = gamma_bb*t2030*t2031*t2035*2.4E1;
                double t2044 = 1.0/pow(rho_b,1.9E1/3.0);
                double t2045 = t2022*t2031*t2042*t2044*t2036*4.8E1;
                double t2046 = t2043+t2045;
                v_rho_b_rho_b[Q] += scale * (c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t2021*t2023*t2034*t2028-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t2023*t2034*t2028*t2038*(8.8E1/9.0)-gamma_bb*t2023*t2035*t2028*t2037*t2046*(8.0/3.0)-gamma_bb*t2021*t2023*t2028*t2037*(gamma_bb*t2030*t2031*t2038*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t2022*t2031*t2042*t2036*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t2031*t2036*1.28E2)/t2041-1.0/pow(rho_b,3.8E1/3.0)*t2022*t2031*1.0/pow(t2041,3.0/2.0)*(t2036*t2036)*1.28E2)*(1.0/2.0)+gamma_bb*t2021*t2023*1.0/pow(t2033,5.0/2.0)*t2028*(t2046*t2046)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t2023*t2034*t2035*t2028*(8.0/3.0)-gamma_bb*t2021*t2023*t2028*t2037*t2046*(1.0/2.0))*(8.0/3.0));
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t2092 = d2*d2;
                double t2093 = gamma_bb+t2092;
                double t2094 = 1.0/pow(rho_b,8.0/3.0);
                double t2099 = gamma_bb*t2094;
                double t2095 = log(t2099+sqrt(t2099*t2099+1.0));
                double t2096 = 1.0/t2093;
                double t2101 = d1*gamma_bb*t2096;
                double t2097 = d0+t2101;
                double t2098 = 1.0/c;
                double t2100 = t2095*t2095;
                double t2102 = t2097*t2097;
                double t2103 = gamma_bb*t2094*t2100*t2102*9.0;
                double t2104 = t2103+1.0;
                double t2105 = 1.0/pow(rho_b,1.6E1/3.0);
                double t2106 = d1*t2096;
                double t2107 = 1.0/(t2093*t2093);
                double t2116 = d1*gamma_bb*t2107;
                double t2108 = t2106-t2116;
                double t2109 = 1.0/sqrt(t2104);
                double t2110 = t2094*t2100*t2102*9.0;
                double t2111 = gamma_bb*gamma_bb;
                double t2112 = t2111*t2105;
                double t2113 = t2112+1.0;
                double t2114 = 1.0/sqrt(t2113);
                double t2115 = gamma_bb*t2095*t2102*t2105*t2114*1.8E1;
                double t2117 = gamma_bb*t2094*t2097*t2100*t2108*1.8E1;
                double t2118 = t2110+t2115+t2117;
                double t2119 = 1.0/pow(t2104,3.0/2.0);
                double t2120 = d1*t2107*2.0;
                double t2121 = 1.0/(t2093*t2093*t2093);
                double t2122 = t2120-d1*gamma_bb*t2121*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * (c*pow(rho_b,4.0/3.0)*(t2094*t2098*t2108*t2109*-2.0+gamma_bb*t2094*t2098*t2122*t2109+t2094*t2097*t2098*t2118*t2119+gamma_bb*t2094*t2098*t2108*t2118*t2119+gamma_bb*t2094*t2097*t2098*t2119*(t2094*t2097*t2100*t2108*3.6E1+t2095*t2102*t2105*t2114*3.6E1+gamma_bb*t2094*t2100*(t2108*t2108)*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t2102*1.8E1)/t2113-1.0/pow(rho_b,3.2E1/3.0)*t2095*t2102*t2111*1.0/pow(t2113,3.0/2.0)*1.8E1-gamma_bb*t2094*t2097*t2100*t2122*1.8E1+gamma_bb*t2095*t2097*t2105*t2114*t2108*7.2E1)*(1.0/2.0)-gamma_bb*t2094*t2097*t2098*1.0/pow(t2104,5.0/2.0)*(t2118*t2118)*(3.0/4.0)));
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t2053 = 1.0/pow(rho_b,8.0/3.0);
                double t2061 = gamma_bb*t2053;
                double t2054 = log(t2061+sqrt(t2061*t2061+1.0));
                double t2056 = d2*d2;
                double t2057 = gamma_bb+t2056;
                double t2058 = 1.0/t2057;
                double t2059 = d1*gamma_bb*t2058;
                double t2055 = d0+t2059;
                double t2060 = 1.0/c;
                double t2062 = t2054*t2054;
                double t2063 = t2055*t2055;
                double t2064 = gamma_bb*t2053*t2062*t2063*9.0;
                double t2065 = t2064+1.0;
                double t2066 = 1.0/sqrt(t2065);
                double t2067 = 1.0/pow(rho_b,1.6E1/3.0);
                double t2068 = d1*t2058;
                double t2069 = 1.0/(t2057*t2057);
                double t2077 = d1*gamma_bb*t2069;
                double t2070 = t2068-t2077;
                double t2071 = 1.0/pow(rho_b,1.1E1/3.0);
                double t2072 = gamma_bb*gamma_bb;
                double t2073 = t2072*t2067;
                double t2074 = t2073+1.0;
                double t2075 = 1.0/sqrt(t2074);
                double t2076 = 1.0/pow(t2065,3.0/2.0);
                double t2078 = 1.0/pow(rho_b,1.9E1/3.0);
                double t2079 = t2053*t2062*t2063*9.0;
                double t2080 = gamma_bb*t2054*t2063*t2075*t2067*1.8E1;
                double t2081 = gamma_bb*t2070*t2053*t2062*t2055*1.8E1;
                double t2082 = t2080+t2081+t2079;
                double t2083 = gamma_bb*t2062*t2071*t2063*2.4E1;
                double t2084 = t2054*t2063*t2072*t2075*t2078*4.8E1;
                double t2085 = t2083+t2084;
                v_rho_b_gamma_bb[Q] += scale * (-c*pow(rho_b,4.0/3.0)*(t2060*t2071*t2055*t2066*(-8.0/3.0)-gamma_bb*t2060*t2070*t2071*t2066*(8.0/3.0)+t2060*t2053*t2055*t2076*t2085*(1.0/2.0)+gamma_bb*t2060*t2070*t2053*t2076*t2085*(1.0/2.0)+gamma_bb*t2060*t2071*t2055*t2082*t2076*(4.0/3.0)+gamma_bb*t2060*t2053*t2055*t2076*(t2062*t2071*t2063*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t2063*t2072*4.8E1)/t2074+gamma_bb*t2070*t2062*t2071*t2055*4.8E1+gamma_bb*t2054*t2063*t2075*t2078*1.44E2+t2070*t2054*t2072*t2055*t2075*t2078*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t2054*t2063*t2072*1.0/pow(t2074,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t2060*t2053*t2055*t2082*1.0/pow(t2065,5.0/2.0)*t2085*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t2060*t2053*t2055*t2066+gamma_bb*t2060*t2070*t2053*t2066-gamma_bb*t2060*t2053*t2055*t2082*t2076*(1.0/2.0))*(4.0/3.0));
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t2137 = 1.0/pow(rho_a,8.0/3.0);
                double t2138 = log(gamma_aa*t2137+sqrt((gamma_aa*gamma_aa)*(t2137*t2137)+1.0));
                double t2140 = d2*d2;
                double t2141 = gamma_aa+t2140;
                double t2142 = 1.0/t2141;
                double t2143 = d1*gamma_aa*t2142;
                double t2139 = d0+t2143;
                v[Q] += scale * (-c*pow(rho_a,4.0/3.0)*((gamma_aa*t2137*t2139*1.0/sqrt(gamma_aa*t2137*(t2138*t2138)*(t2139*t2139)*9.0+1.0))/c-1.0));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t2145 = 1.0/pow(rho_a,8.0/3.0);
                double t2153 = gamma_aa*t2145;
                double t2146 = log(t2153+sqrt(t2153*t2153+1.0));
                double t2148 = d2*d2;
                double t2149 = gamma_aa+t2148;
                double t2150 = 1.0/t2149;
                double t2151 = d1*gamma_aa*t2150;
                double t2147 = d0+t2151;
                double t2152 = 1.0/c;
                double t2154 = t2146*t2146;
                double t2155 = t2147*t2147;
                double t2156 = gamma_aa*t2145*t2154*t2155*9.0;
                double t2157 = t2156+1.0;
                double t2158 = 1.0/sqrt(t2157);
                double t2159 = 1.0/pow(rho_a,1.1E1/3.0);
                double t2160 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * (c*pow(rho_a,1.0/3.0)*(gamma_aa*t2152*t2145*t2147*t2158-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t2152*t2147*t2158*t2159*(8.0/3.0)-gamma_aa*t2152*t2145*t2147*1.0/pow(t2157,3.0/2.0)*(gamma_aa*t2154*t2155*t2159*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t2160*t2146*t2155*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t2160+1.0)*4.8E1)*(1.0/2.0)));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t2163 = 1.0/pow(rho_a,8.0/3.0);
                double t2171 = gamma_aa*t2163;
                double t2164 = log(t2171+sqrt(t2171*t2171+1.0));
                double t2166 = d2*d2;
                double t2167 = gamma_aa+t2166;
                double t2168 = 1.0/t2167;
                double t2169 = d1*gamma_aa*t2168;
                double t2165 = d0+t2169;
                double t2170 = 1.0/c;
                double t2172 = t2164*t2164;
                double t2173 = t2165*t2165;
                double t2174 = gamma_aa*t2163*t2172*t2173*9.0;
                double t2175 = t2174+1.0;
                double t2176 = 1.0/sqrt(t2175);
                double t2177 = 1.0/pow(rho_a,1.6E1/3.0);
                double t2178 = d1*t2168;
                double t2179 = 1.0/(t2167*t2167);
                double t2180 = t2178-d1*gamma_aa*t2179;
                v_gamma_aa[Q] += scale * (-c*pow(rho_a,4.0/3.0)*(t2170*t2163*t2165*t2176+gamma_aa*t2170*t2180*t2163*t2176-gamma_aa*t2170*t2163*t2165*1.0/pow(t2175,3.0/2.0)*(t2163*t2172*t2173*9.0+gamma_aa*t2164*t2173*t2177*1.0/sqrt((gamma_aa*gamma_aa)*t2177+1.0)*1.8E1+gamma_aa*t2180*t2163*t2172*t2165*1.8E1)*(1.0/2.0)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t2186 = 1.0/pow(rho_a,8.0/3.0);
                double t2194 = gamma_aa*t2186;
                double t2187 = log(t2194+sqrt(t2194*t2194+1.0));
                double t2189 = d2*d2;
                double t2190 = gamma_aa+t2189;
                double t2191 = 1.0/t2190;
                double t2192 = d1*gamma_aa*t2191;
                double t2188 = d0+t2192;
                double t2193 = 1.0/c;
                double t2195 = t2187*t2187;
                double t2196 = t2188*t2188;
                double t2197 = gamma_aa*t2186*t2195*t2196*9.0;
                double t2198 = t2197+1.0;
                double t2199 = 1.0/sqrt(t2198);
                double t2200 = 1.0/pow(rho_a,1.1E1/3.0);
                double t2201 = gamma_aa*gamma_aa;
                double t2202 = 1.0/pow(t2198,3.0/2.0);
                double t2203 = 1.0/pow(rho_a,1.4E1/3.0);
                double t2204 = 1.0/pow(rho_a,1.6E1/3.0);
                double t2205 = t2201*t2204;
                double t2206 = t2205+1.0;
                double t2207 = 1.0/sqrt(t2206);
                double t2208 = gamma_aa*t2195*t2196*t2200*2.4E1;
                double t2209 = 1.0/pow(rho_a,1.9E1/3.0);
                double t2210 = t2187*t2196*t2201*t2207*t2209*4.8E1;
                double t2211 = t2210+t2208;
                v_rho_a_rho_a[Q] += scale * (c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t2193*t2186*t2188*t2199-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t2193*t2188*t2199*t2203*(8.8E1/9.0)-gamma_aa*t2193*t2188*t2200*t2202*t2211*(8.0/3.0)-gamma_aa*t2193*t2186*t2188*t2202*(gamma_aa*t2195*t2196*t2203*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t2187*t2196*t2201*t2207*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t2196*t2201*1.28E2)/t2206-1.0/pow(rho_a,3.8E1/3.0)*t2187*t2196*(t2201*t2201)*1.0/pow(t2206,3.0/2.0)*1.28E2)*(1.0/2.0)+gamma_aa*t2193*t2186*t2188*1.0/pow(t2198,5.0/2.0)*(t2211*t2211)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t2193*t2188*t2199*t2200*(8.0/3.0)-gamma_aa*t2193*t2186*t2188*t2202*t2211*(1.0/2.0))*(8.0/3.0));
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t2254 = d2*d2;
                double t2255 = gamma_aa+t2254;
                double t2256 = 1.0/pow(rho_a,8.0/3.0);
                double t2261 = gamma_aa*t2256;
                double t2257 = log(t2261+sqrt(t2261*t2261+1.0));
                double t2258 = 1.0/t2255;
                double t2263 = d1*gamma_aa*t2258;
                double t2259 = d0+t2263;
                double t2260 = 1.0/c;
                double t2262 = t2257*t2257;
                double t2264 = t2259*t2259;
                double t2265 = gamma_aa*t2262*t2264*t2256*9.0;
                double t2266 = t2265+1.0;
                double t2267 = 1.0/pow(rho_a,1.6E1/3.0);
                double t2268 = d1*t2258;
                double t2269 = 1.0/(t2255*t2255);
                double t2278 = d1*gamma_aa*t2269;
                double t2270 = t2268-t2278;
                double t2271 = 1.0/sqrt(t2266);
                double t2272 = t2262*t2264*t2256*9.0;
                double t2273 = gamma_aa*gamma_aa;
                double t2274 = t2273*t2267;
                double t2275 = t2274+1.0;
                double t2276 = 1.0/sqrt(t2275);
                double t2277 = gamma_aa*t2264*t2257*t2267*t2276*1.8E1;
                double t2279 = gamma_aa*t2270*t2262*t2256*t2259*1.8E1;
                double t2280 = t2272+t2277+t2279;
                double t2281 = 1.0/pow(t2266,3.0/2.0);
                double t2282 = d1*t2269*2.0;
                double t2283 = 1.0/(t2255*t2255*t2255);
                double t2284 = t2282-d1*gamma_aa*t2283*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * (c*pow(rho_a,4.0/3.0)*(t2260*t2270*t2271*t2256*-2.0+gamma_aa*t2260*t2271*t2256*t2284+t2260*t2280*t2281*t2256*t2259+gamma_aa*t2260*t2270*t2280*t2281*t2256+gamma_aa*t2260*t2281*t2256*t2259*(t2270*t2262*t2256*t2259*3.6E1+t2264*t2257*t2267*t2276*3.6E1+gamma_aa*(t2270*t2270)*t2262*t2256*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t2264*1.8E1)/t2275-1.0/pow(rho_a,3.2E1/3.0)*t2264*t2273*t2257*1.0/pow(t2275,3.0/2.0)*1.8E1-gamma_aa*t2262*t2256*t2284*t2259*1.8E1+gamma_aa*t2270*t2257*t2267*t2276*t2259*7.2E1)*(1.0/2.0)-gamma_aa*t2260*(t2280*t2280)*t2256*1.0/pow(t2266,5.0/2.0)*t2259*(3.0/4.0)));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t2215 = 1.0/pow(rho_a,8.0/3.0);
                double t2223 = gamma_aa*t2215;
                double t2216 = log(t2223+sqrt(t2223*t2223+1.0));
                double t2218 = d2*d2;
                double t2219 = gamma_aa+t2218;
                double t2220 = 1.0/t2219;
                double t2221 = d1*gamma_aa*t2220;
                double t2217 = d0+t2221;
                double t2222 = 1.0/c;
                double t2224 = t2216*t2216;
                double t2225 = t2217*t2217;
                double t2226 = gamma_aa*t2215*t2224*t2225*9.0;
                double t2227 = t2226+1.0;
                double t2228 = 1.0/sqrt(t2227);
                double t2229 = 1.0/pow(rho_a,1.6E1/3.0);
                double t2230 = d1*t2220;
                double t2231 = 1.0/(t2219*t2219);
                double t2239 = d1*gamma_aa*t2231;
                double t2232 = t2230-t2239;
                double t2233 = 1.0/pow(rho_a,1.1E1/3.0);
                double t2234 = gamma_aa*gamma_aa;
                double t2235 = t2234*t2229;
                double t2236 = t2235+1.0;
                double t2237 = 1.0/sqrt(t2236);
                double t2238 = 1.0/pow(t2227,3.0/2.0);
                double t2240 = 1.0/pow(rho_a,1.9E1/3.0);
                double t2241 = t2215*t2224*t2225*9.0;
                double t2242 = gamma_aa*t2216*t2225*t2237*t2229*1.8E1;
                double t2243 = gamma_aa*t2232*t2215*t2224*t2217*1.8E1;
                double t2244 = t2241+t2242+t2243;
                double t2245 = gamma_aa*t2224*t2233*t2225*2.4E1;
                double t2246 = t2240*t2216*t2225*t2234*t2237*4.8E1;
                double t2247 = t2245+t2246;
                v_rho_a_gamma_aa[Q] += scale * (-c*pow(rho_a,4.0/3.0)*(t2222*t2233*t2217*t2228*(-8.0/3.0)-gamma_aa*t2222*t2232*t2233*t2228*(8.0/3.0)+t2222*t2215*t2217*t2238*t2247*(1.0/2.0)+gamma_aa*t2222*t2232*t2215*t2238*t2247*(1.0/2.0)+gamma_aa*t2222*t2233*t2217*t2244*t2238*(4.0/3.0)+gamma_aa*t2222*t2215*t2217*t2238*(t2224*t2233*t2225*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t2225*t2234*4.8E1)/t2236+gamma_aa*t2232*t2224*t2233*t2217*4.8E1+gamma_aa*t2240*t2216*t2225*t2237*1.44E2+t2240*t2232*t2216*t2234*t2217*t2237*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t2216*t2225*t2234*1.0/pow(t2236,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_aa*t2222*t2215*t2217*t2244*1.0/pow(t2227,5.0/2.0)*t2247*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t2222*t2215*t2217*t2228+gamma_aa*t2222*t2232*t2215*t2228-gamma_aa*t2222*t2215*t2217*t2244*t2238*(1.0/2.0))*(4.0/3.0));
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t1672 = 1.0/pow(rho_a,8.0/3.0);
                double t1673 = log(gamma_aa*t1672+sqrt((gamma_aa*gamma_aa)*(t1672*t1672)+1.0));
                double t1675 = d2*d2;
                double t1676 = gamma_aa+t1675;
                double t1677 = 1.0/t1676;
                double t1678 = d1*gamma_aa*t1677;
                double t1674 = d0+t1678;
                double t1679 = 1.0/c;
                double t1680 = 1.0/pow(rho_b,8.0/3.0);
                double t1681 = log(gamma_bb*t1680+sqrt((gamma_bb*gamma_bb)*(t1680*t1680)+1.0));
                double t1683 = gamma_bb+t1675;
                double t1684 = 1.0/t1683;
                double t1685 = d1*gamma_bb*t1684;
                double t1682 = d0+t1685;
                v[Q] += scale * (-c*pow(rho_a,4.0/3.0)*(gamma_aa*t1672*t1674*t1679*1.0/sqrt(gamma_aa*t1672*(t1673*t1673)*(t1674*t1674)*9.0+1.0)-1.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t1680*t1682*t1679*1.0/sqrt(gamma_bb*t1680*(t1681*t1681)*(t1682*t1682)*9.0+1.0)-1.0));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t1687 = 1.0/pow(rho_a,8.0/3.0);
                double t1695 = gamma_aa*t1687;
                double t1688 = log(t1695+sqrt(t1695*t1695+1.0));
                double t1690 = d2*d2;
                double t1691 = gamma_aa+t1690;
                double t1692 = 1.0/t1691;
                double t1693 = d1*gamma_aa*t1692;
                double t1689 = d0+t1693;
                double t1694 = 1.0/c;
                double t1696 = t1688*t1688;
                double t1697 = t1689*t1689;
                double t1698 = gamma_aa*t1687*t1696*t1697*9.0;
                double t1699 = t1698+1.0;
                double t1700 = 1.0/sqrt(t1699);
                double t1701 = 1.0/pow(rho_a,1.1E1/3.0);
                double t1702 = gamma_aa*gamma_aa;
                v_rho_a[Q] += scale * (c*pow(rho_a,1.0/3.0)*(gamma_aa*t1694*t1687*t1689*t1700-1.0)*(-4.0/3.0)+c*pow(rho_a,4.0/3.0)*(gamma_aa*t1694*t1689*t1700*t1701*(8.0/3.0)-gamma_aa*t1694*t1687*t1689*1.0/pow(t1699,3.0/2.0)*(gamma_aa*t1696*t1697*t1701*2.4E1+1.0/pow(rho_a,1.9E1/3.0)*t1688*t1697*t1702*1.0/sqrt(1.0/pow(rho_a,1.6E1/3.0)*t1702+1.0)*4.8E1)*(1.0/2.0)));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t1704 = 1.0/pow(rho_b,8.0/3.0);
                double t1712 = gamma_bb*t1704;
                double t1705 = log(t1712+sqrt(t1712*t1712+1.0));
                double t1707 = d2*d2;
                double t1708 = gamma_bb+t1707;
                double t1709 = 1.0/t1708;
                double t1710 = d1*gamma_bb*t1709;
                double t1706 = d0+t1710;
                double t1711 = 1.0/c;
                double t1713 = t1705*t1705;
                double t1714 = t1706*t1706;
                double t1715 = gamma_bb*t1704*t1713*t1714*9.0;
                double t1716 = t1715+1.0;
                double t1717 = 1.0/sqrt(t1716);
                double t1718 = 1.0/pow(rho_b,1.1E1/3.0);
                double t1719 = gamma_bb*gamma_bb;
                v_rho_b[Q] += scale * (c*pow(rho_b,1.0/3.0)*(gamma_bb*t1711*t1704*t1706*t1717-1.0)*(-4.0/3.0)+c*pow(rho_b,4.0/3.0)*(gamma_bb*t1711*t1706*t1717*t1718*(8.0/3.0)-gamma_bb*t1711*t1704*t1706*1.0/pow(t1716,3.0/2.0)*(gamma_bb*t1713*t1714*t1718*2.4E1+1.0/pow(rho_b,1.9E1/3.0)*t1705*t1714*t1719*1.0/sqrt(1.0/pow(rho_b,1.6E1/3.0)*t1719+1.0)*4.8E1)*(1.0/2.0)));
            }
            
            // v_gamma_aa
            if (deriv >= 1) {
                double t1721 = 1.0/pow(rho_a,8.0/3.0);
                double t1729 = gamma_aa*t1721;
                double t1722 = log(t1729+sqrt(t1729*t1729+1.0));
                double t1724 = d2*d2;
                double t1725 = gamma_aa+t1724;
                double t1726 = 1.0/t1725;
                double t1727 = d1*gamma_aa*t1726;
                double t1723 = d0+t1727;
                double t1728 = 1.0/c;
                double t1730 = t1722*t1722;
                double t1731 = t1723*t1723;
                double t1732 = gamma_aa*t1721*t1730*t1731*9.0;
                double t1733 = t1732+1.0;
                double t1734 = 1.0/sqrt(t1733);
                double t1735 = 1.0/pow(rho_a,1.6E1/3.0);
                double t1736 = d1*t1726;
                double t1737 = 1.0/(t1725*t1725);
                double t1738 = t1736-d1*gamma_aa*t1737;
                v_gamma_aa[Q] += scale * (-c*pow(rho_a,4.0/3.0)*(t1721*t1723*t1734*t1728+gamma_aa*t1721*t1734*t1728*t1738-gamma_aa*t1721*t1723*1.0/pow(t1733,3.0/2.0)*t1728*(t1721*t1730*t1731*9.0+gamma_aa*t1722*t1731*t1735*1.0/sqrt((gamma_aa*gamma_aa)*t1735+1.0)*1.8E1+gamma_aa*t1721*t1730*t1723*t1738*1.8E1)*(1.0/2.0)));
            }
            
            // v_gamma_bb
            if (deriv >= 1) {
                double t1741 = 1.0/pow(rho_b,8.0/3.0);
                double t1749 = gamma_bb*t1741;
                double t1742 = log(t1749+sqrt(t1749*t1749+1.0));
                double t1744 = d2*d2;
                double t1745 = gamma_bb+t1744;
                double t1746 = 1.0/t1745;
                double t1747 = d1*gamma_bb*t1746;
                double t1743 = d0+t1747;
                double t1748 = 1.0/c;
                double t1750 = t1742*t1742;
                double t1751 = t1743*t1743;
                double t1752 = gamma_bb*t1741*t1750*t1751*9.0;
                double t1753 = t1752+1.0;
                double t1754 = 1.0/sqrt(t1753);
                double t1755 = 1.0/pow(rho_b,1.6E1/3.0);
                double t1756 = d1*t1746;
                double t1757 = 1.0/(t1745*t1745);
                double t1758 = t1756-d1*gamma_bb*t1757;
                v_gamma_bb[Q] += scale * (-c*pow(rho_b,4.0/3.0)*(t1741*t1743*t1754*t1748+gamma_bb*t1741*t1754*t1748*t1758-gamma_bb*t1741*t1743*1.0/pow(t1753,3.0/2.0)*t1748*(t1741*t1750*t1751*9.0+gamma_bb*t1742*t1751*t1755*1.0/sqrt((gamma_bb*gamma_bb)*t1755+1.0)*1.8E1+gamma_bb*t1741*t1750*t1743*t1758*1.8E1)*(1.0/2.0)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t1762 = 1.0/pow(rho_a,8.0/3.0);
                double t1770 = gamma_aa*t1762;
                double t1763 = log(t1770+sqrt(t1770*t1770+1.0));
                double t1765 = d2*d2;
                double t1766 = gamma_aa+t1765;
                double t1767 = 1.0/t1766;
                double t1768 = d1*gamma_aa*t1767;
                double t1764 = d0+t1768;
                double t1769 = 1.0/c;
                double t1771 = t1763*t1763;
                double t1772 = t1764*t1764;
                double t1773 = gamma_aa*t1762*t1771*t1772*9.0;
                double t1774 = t1773+1.0;
                double t1775 = 1.0/sqrt(t1774);
                double t1776 = 1.0/pow(rho_a,1.1E1/3.0);
                double t1777 = gamma_aa*gamma_aa;
                double t1778 = 1.0/pow(t1774,3.0/2.0);
                double t1779 = 1.0/pow(rho_a,1.4E1/3.0);
                double t1780 = 1.0/pow(rho_a,1.6E1/3.0);
                double t1781 = t1780*t1777;
                double t1782 = t1781+1.0;
                double t1783 = 1.0/sqrt(t1782);
                double t1784 = gamma_aa*t1771*t1772*t1776*2.4E1;
                double t1785 = 1.0/pow(rho_a,1.9E1/3.0);
                double t1786 = t1763*t1772*t1783*t1785*t1777*4.8E1;
                double t1787 = t1784+t1786;
                v_rho_a_rho_a[Q] += scale * (c*1.0/pow(rho_a,2.0/3.0)*(gamma_aa*t1762*t1764*t1775*t1769-1.0)*(-4.0/9.0)-c*pow(rho_a,4.0/3.0)*(gamma_aa*t1764*t1775*t1769*t1779*(8.8E1/9.0)-gamma_aa*t1764*t1776*t1769*t1778*t1787*(8.0/3.0)-gamma_aa*t1762*t1764*t1769*t1778*(gamma_aa*t1771*t1772*t1779*8.8E1+1.0/pow(rho_a,2.2E1/3.0)*t1763*t1772*t1783*t1777*4.32E2+(gamma_aa*1.0/pow(rho_a,1.0E1)*t1772*t1777*1.28E2)/t1782-1.0/pow(rho_a,3.8E1/3.0)*t1763*t1772*1.0/pow(t1782,3.0/2.0)*(t1777*t1777)*1.28E2)*(1.0/2.0)+gamma_aa*t1762*t1764*1.0/pow(t1774,5.0/2.0)*t1769*(t1787*t1787)*(3.0/4.0))+c*pow(rho_a,1.0/3.0)*(gamma_aa*t1764*t1775*t1776*t1769*(8.0/3.0)-gamma_aa*t1762*t1764*t1769*t1778*t1787*(1.0/2.0))*(8.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t1790 = 1.0/pow(rho_b,8.0/3.0);
                double t1798 = gamma_bb*t1790;
                double t1791 = log(t1798+sqrt(t1798*t1798+1.0));
                double t1793 = d2*d2;
                double t1794 = gamma_bb+t1793;
                double t1795 = 1.0/t1794;
                double t1796 = d1*gamma_bb*t1795;
                double t1792 = d0+t1796;
                double t1797 = 1.0/c;
                double t1799 = t1791*t1791;
                double t1800 = t1792*t1792;
                double t1801 = gamma_bb*t1790*t1799*t1800*9.0;
                double t1802 = t1801+1.0;
                double t1803 = 1.0/sqrt(t1802);
                double t1804 = 1.0/pow(rho_b,1.1E1/3.0);
                double t1805 = gamma_bb*gamma_bb;
                double t1806 = 1.0/pow(t1802,3.0/2.0);
                double t1807 = 1.0/pow(rho_b,1.4E1/3.0);
                double t1808 = 1.0/pow(rho_b,1.6E1/3.0);
                double t1809 = t1805*t1808;
                double t1810 = t1809+1.0;
                double t1811 = 1.0/sqrt(t1810);
                double t1812 = gamma_bb*t1799*t1800*t1804*2.4E1;
                double t1813 = 1.0/pow(rho_b,1.9E1/3.0);
                double t1814 = t1791*t1800*t1811*t1813*t1805*4.8E1;
                double t1815 = t1812+t1814;
                v_rho_b_rho_b[Q] += scale * (c*1.0/pow(rho_b,2.0/3.0)*(gamma_bb*t1790*t1792*t1797*t1803-1.0)*(-4.0/9.0)-c*pow(rho_b,4.0/3.0)*(gamma_bb*t1792*t1797*t1803*t1807*(8.8E1/9.0)-gamma_bb*t1792*t1797*t1804*t1806*t1815*(8.0/3.0)-gamma_bb*t1790*t1792*t1797*t1806*(gamma_bb*t1799*t1800*t1807*8.8E1+1.0/pow(rho_b,2.2E1/3.0)*t1791*t1800*t1811*t1805*4.32E2+(gamma_bb*1.0/pow(rho_b,1.0E1)*t1800*t1805*1.28E2)/t1810-1.0/pow(rho_b,3.8E1/3.0)*t1791*t1800*1.0/pow(t1810,3.0/2.0)*(t1805*t1805)*1.28E2)*(1.0/2.0)+gamma_bb*t1790*t1792*t1797*1.0/pow(t1802,5.0/2.0)*(t1815*t1815)*(3.0/4.0))+c*pow(rho_b,1.0/3.0)*(gamma_bb*t1792*t1797*t1803*t1804*(8.0/3.0)-gamma_bb*t1790*t1792*t1797*t1806*t1815*(1.0/2.0))*(8.0/3.0));
            }
            
            // v_gamma_aa_gamma_aa
            if (deriv >= 2) {
                double t1889 = d2*d2;
                double t1890 = gamma_aa+t1889;
                double t1891 = 1.0/pow(rho_a,8.0/3.0);
                double t1896 = gamma_aa*t1891;
                double t1892 = log(t1896+sqrt(t1896*t1896+1.0));
                double t1893 = 1.0/t1890;
                double t1898 = d1*gamma_aa*t1893;
                double t1894 = d0+t1898;
                double t1895 = 1.0/c;
                double t1897 = t1892*t1892;
                double t1899 = t1894*t1894;
                double t1900 = gamma_aa*t1891*t1897*t1899*9.0;
                double t1901 = t1900+1.0;
                double t1902 = 1.0/pow(rho_a,1.6E1/3.0);
                double t1903 = d1*t1893;
                double t1904 = 1.0/(t1890*t1890);
                double t1913 = d1*gamma_aa*t1904;
                double t1905 = t1903-t1913;
                double t1906 = 1.0/sqrt(t1901);
                double t1907 = t1891*t1897*t1899*9.0;
                double t1908 = gamma_aa*gamma_aa;
                double t1909 = t1902*t1908;
                double t1910 = t1909+1.0;
                double t1911 = 1.0/sqrt(t1910);
                double t1912 = gamma_aa*t1892*t1899*t1902*t1911*1.8E1;
                double t1914 = gamma_aa*t1891*t1894*t1897*t1905*1.8E1;
                double t1915 = t1912+t1914+t1907;
                double t1916 = 1.0/pow(t1901,3.0/2.0);
                double t1917 = d1*t1904*2.0;
                double t1918 = 1.0/(t1890*t1890*t1890);
                double t1919 = t1917-d1*gamma_aa*t1918*2.0;
                v_gamma_aa_gamma_aa[Q] += scale * (c*pow(rho_a,4.0/3.0)*(t1891*t1895*t1905*t1906*-2.0+gamma_aa*t1891*t1895*t1906*t1919+t1891*t1894*t1895*t1915*t1916+gamma_aa*t1891*t1895*t1905*t1915*t1916+gamma_aa*t1891*t1894*t1895*t1916*(t1891*t1894*t1897*t1905*3.6E1+t1892*t1899*t1902*t1911*3.6E1+gamma_aa*t1891*t1897*(t1905*t1905)*1.8E1+(gamma_aa*1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t1899*1.8E1)/t1910-1.0/pow(rho_a,3.2E1/3.0)*t1892*t1899*1.0/pow(t1910,3.0/2.0)*t1908*1.8E1-gamma_aa*t1891*t1894*t1897*t1919*1.8E1+gamma_aa*t1892*t1894*t1902*t1911*t1905*7.2E1)*(1.0/2.0)-gamma_aa*t1891*t1894*t1895*1.0/pow(t1901,5.0/2.0)*(t1915*t1915)*(3.0/4.0)));
            }
            
            // v_gamma_bb_gamma_bb
            if (deriv >= 2) {
                double t1925 = d2*d2;
                double t1926 = gamma_bb+t1925;
                double t1927 = 1.0/pow(rho_b,8.0/3.0);
                double t1932 = gamma_bb*t1927;
                double t1928 = log(t1932+sqrt(t1932*t1932+1.0));
                double t1929 = 1.0/t1926;
                double t1934 = d1*gamma_bb*t1929;
                double t1930 = d0+t1934;
                double t1931 = 1.0/c;
                double t1933 = t1928*t1928;
                double t1935 = t1930*t1930;
                double t1936 = gamma_bb*t1933*t1935*t1927*9.0;
                double t1937 = t1936+1.0;
                double t1938 = 1.0/pow(rho_b,1.6E1/3.0);
                double t1939 = d1*t1929;
                double t1940 = 1.0/(t1926*t1926);
                double t1949 = d1*gamma_bb*t1940;
                double t1941 = t1939-t1949;
                double t1942 = 1.0/sqrt(t1937);
                double t1943 = t1933*t1935*t1927*9.0;
                double t1944 = gamma_bb*gamma_bb;
                double t1945 = t1944*t1938;
                double t1946 = t1945+1.0;
                double t1947 = 1.0/sqrt(t1946);
                double t1948 = gamma_bb*t1935*t1928*t1938*t1947*1.8E1;
                double t1950 = gamma_bb*t1930*t1941*t1933*t1927*1.8E1;
                double t1951 = t1950+t1943+t1948;
                double t1952 = 1.0/pow(t1937,3.0/2.0);
                double t1953 = d1*t1940*2.0;
                double t1954 = 1.0/(t1926*t1926*t1926);
                double t1955 = t1953-d1*gamma_bb*t1954*2.0;
                v_gamma_bb_gamma_bb[Q] += scale * (c*pow(rho_b,4.0/3.0)*(t1931*t1941*t1942*t1927*-2.0+gamma_bb*t1931*t1942*t1927*t1955+t1930*t1931*t1951*t1952*t1927+gamma_bb*t1931*t1941*t1951*t1952*t1927+gamma_bb*t1930*t1931*t1952*t1927*(t1930*t1941*t1933*t1927*3.6E1+t1935*t1928*t1938*t1947*3.6E1+gamma_bb*(t1941*t1941)*t1933*t1927*1.8E1+(gamma_bb*1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t1935*1.8E1)/t1946-1.0/pow(rho_b,3.2E1/3.0)*t1935*t1944*t1928*1.0/pow(t1946,3.0/2.0)*1.8E1-gamma_bb*t1930*t1933*t1927*t1955*1.8E1+gamma_bb*t1930*t1941*t1928*t1938*t1947*7.2E1)*(1.0/2.0)-gamma_bb*t1930*t1931*(t1951*t1951)*t1927*1.0/pow(t1937,5.0/2.0)*(3.0/4.0)));
            }
            
            // v_rho_a_gamma_aa
            if (deriv >= 2) {
                double t1817 = 1.0/pow(rho_a,8.0/3.0);
                double t1825 = gamma_aa*t1817;
                double t1818 = log(t1825+sqrt(t1825*t1825+1.0));
                double t1820 = d2*d2;
                double t1821 = gamma_aa+t1820;
                double t1822 = 1.0/t1821;
                double t1823 = d1*gamma_aa*t1822;
                double t1819 = d0+t1823;
                double t1824 = 1.0/c;
                double t1826 = t1818*t1818;
                double t1827 = t1819*t1819;
                double t1828 = gamma_aa*t1817*t1826*t1827*9.0;
                double t1829 = t1828+1.0;
                double t1830 = 1.0/sqrt(t1829);
                double t1831 = 1.0/pow(rho_a,1.6E1/3.0);
                double t1832 = d1*t1822;
                double t1833 = 1.0/(t1821*t1821);
                double t1841 = d1*gamma_aa*t1833;
                double t1834 = t1832-t1841;
                double t1835 = 1.0/pow(rho_a,1.1E1/3.0);
                double t1836 = gamma_aa*gamma_aa;
                double t1837 = t1831*t1836;
                double t1838 = t1837+1.0;
                double t1839 = 1.0/sqrt(t1838);
                double t1840 = 1.0/pow(t1829,3.0/2.0);
                double t1842 = 1.0/pow(rho_a,1.9E1/3.0);
                double t1843 = t1817*t1826*t1827*9.0;
                double t1844 = gamma_aa*t1831*t1818*t1827*t1839*1.8E1;
                double t1845 = gamma_aa*t1834*t1817*t1826*t1819*1.8E1;
                double t1846 = t1843+t1844+t1845;
                double t1847 = gamma_aa*t1826*t1835*t1827*2.4E1;
                double t1848 = t1842*t1818*t1827*t1836*t1839*4.8E1;
                double t1849 = t1847+t1848;
                v_rho_a_gamma_aa[Q] += scale * (-c*pow(rho_a,4.0/3.0)*(t1830*t1824*t1835*t1819*(-8.0/3.0)-gamma_aa*t1830*t1824*t1834*t1835*(8.0/3.0)+t1840*t1824*t1817*t1819*t1849*(1.0/2.0)+gamma_aa*t1840*t1824*t1834*t1817*t1849*(1.0/2.0)+gamma_aa*t1840*t1824*t1835*t1819*t1846*(4.0/3.0)+gamma_aa*t1840*t1824*t1817*t1819*(t1826*t1835*t1827*2.4E1+(1.0/(rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a*rho_a)*t1827*t1836*4.8E1)/t1838+gamma_aa*t1834*t1826*t1835*t1819*4.8E1+gamma_aa*t1842*t1818*t1827*t1839*1.44E2+t1842*t1834*t1818*t1836*t1819*t1839*9.6E1-gamma_aa*1.0/pow(rho_a,3.5E1/3.0)*t1818*t1827*t1836*1.0/pow(t1838,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_aa*t1824*t1817*t1819*t1846*1.0/pow(t1829,5.0/2.0)*t1849*(3.0/4.0))-c*pow(rho_a,1.0/3.0)*(t1830*t1824*t1817*t1819+gamma_aa*t1830*t1824*t1834*t1817-gamma_aa*t1840*t1824*t1817*t1819*t1846*(1.0/2.0))*(4.0/3.0));
            }
            
            // v_rho_b_gamma_bb
            if (deriv >= 2) {
                double t1855 = 1.0/pow(rho_b,8.0/3.0);
                double t1863 = gamma_bb*t1855;
                double t1856 = log(t1863+sqrt(t1863*t1863+1.0));
                double t1858 = d2*d2;
                double t1859 = gamma_bb+t1858;
                double t1860 = 1.0/t1859;
                double t1861 = d1*gamma_bb*t1860;
                double t1857 = d0+t1861;
                double t1862 = 1.0/c;
                double t1864 = t1856*t1856;
                double t1865 = t1857*t1857;
                double t1866 = gamma_bb*t1855*t1864*t1865*9.0;
                double t1867 = t1866+1.0;
                double t1868 = 1.0/sqrt(t1867);
                double t1869 = 1.0/pow(rho_b,1.6E1/3.0);
                double t1870 = d1*t1860;
                double t1871 = 1.0/(t1859*t1859);
                double t1879 = d1*gamma_bb*t1871;
                double t1872 = t1870-t1879;
                double t1873 = 1.0/pow(rho_b,1.1E1/3.0);
                double t1874 = gamma_bb*gamma_bb;
                double t1875 = t1874*t1869;
                double t1876 = t1875+1.0;
                double t1877 = 1.0/sqrt(t1876);
                double t1878 = 1.0/pow(t1867,3.0/2.0);
                double t1880 = 1.0/pow(rho_b,1.9E1/3.0);
                double t1881 = t1855*t1864*t1865*9.0;
                double t1882 = gamma_bb*t1856*t1865*t1877*t1869*1.8E1;
                double t1883 = gamma_bb*t1872*t1855*t1864*t1857*1.8E1;
                double t1884 = t1881+t1882+t1883;
                double t1885 = gamma_bb*t1864*t1873*t1865*2.4E1;
                double t1886 = t1880*t1856*t1865*t1874*t1877*4.8E1;
                double t1887 = t1885+t1886;
                v_rho_b_gamma_bb[Q] += scale * (-c*pow(rho_b,4.0/3.0)*(t1862*t1873*t1857*t1868*(-8.0/3.0)-gamma_bb*t1862*t1872*t1873*t1868*(8.0/3.0)+t1862*t1855*t1857*t1878*t1887*(1.0/2.0)+gamma_bb*t1862*t1872*t1855*t1878*t1887*(1.0/2.0)+gamma_bb*t1862*t1873*t1857*t1884*t1878*(4.0/3.0)+gamma_bb*t1862*t1855*t1857*t1878*(t1864*t1873*t1865*2.4E1+(1.0/(rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b*rho_b)*t1865*t1874*4.8E1)/t1876+gamma_bb*t1872*t1864*t1873*t1857*4.8E1+gamma_bb*t1880*t1856*t1865*t1877*1.44E2+t1880*t1872*t1856*t1874*t1857*t1877*9.6E1-gamma_bb*1.0/pow(rho_b,3.5E1/3.0)*t1856*t1865*t1874*1.0/pow(t1876,3.0/2.0)*4.8E1)*(1.0/2.0)-gamma_bb*t1862*t1855*t1857*t1884*1.0/pow(t1867,5.0/2.0)*t1887*(3.0/4.0))-c*pow(rho_b,1.0/3.0)*(t1862*t1855*t1857*t1868+gamma_bb*t1862*t1872*t1855*t1868-gamma_bb*t1862*t1855*t1857*t1884*t1878*(1.0/2.0))*(4.0/3.0));
            }
            
        }
    }
}

}
