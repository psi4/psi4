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
                double t2784 = rho_a+rho_b;
                double t2785 = 1.0/pow(t2784,1.0/3.0);
                double t2786 = c*t2785;
                double t2787 = log(t2786);
                double t2788 = EcPhd_1*t2787;
                double t2789 = pow(2.0,1.0/3.0);
                double t2790 = t2789*2.0;
                double t2791 = t2790-2.0;
                double t2792 = two_13*2.0;
                double t2793 = t2792-2.0;
                double t2794 = 1.0/t2793;
                double t2795 = sqrt(t2786);
                double t2796 = EcPld_2*t2795;
                double t2797 = EcPld_3*c*t2785;
                double t2798 = t2796+t2797+1.0;
                double t2799 = 1.0/t2798;
                double t2800 = EcPld_1*t2799;
                v[Q] += scale * (t2784*(heaviside(-c*t2785+1.0)*(EcPhd_2+t2788+t2791*t2794*(EcFhd_2-EcPhd_2-t2788+EcFhd_1*t2787+EcFhd_4*c*t2785-EcPhd_4*c*t2785+EcFhd_3*c*t2785*t2787-EcPhd_3*c*t2785*t2787)+EcPhd_4*c*t2785+EcPhd_3*c*t2785*t2787)+heaviside(t2786-1.0)*(t2800-t2791*t2794*(t2800-EcFld_1/(EcFld_2*t2795+EcFld_3*c*t2785+1.0)))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t2802 = rho_a+rho_b;
                double t2803 = 1.0/pow(t2802,4.0/3.0);
                double t2804 = 1.0/pow(t2802,1.0/3.0);
                double t2805 = c*t2804;
                double t2806 = 1.0/sqrt(t2805);
                double t2807 = sqrt(t2805);
                double t2808 = EcPld_3*c*t2803*(1.0/3.0);
                double t2809 = EcPld_2*c*t2803*t2806*(1.0/6.0);
                double t2810 = t2808+t2809;
                double t2811 = EcPld_2*t2807;
                double t2812 = EcPld_3*c*t2804;
                double t2813 = t2811+t2812+1.0;
                double t2814 = 1.0/(t2813*t2813);
                double t2815 = EcPld_1*t2810*t2814;
                double t2816 = pow(2.0,1.0/3.0);
                double t2817 = t2816*2.0;
                double t2818 = t2817-2.0;
                double t2819 = two_13*2.0;
                double t2820 = t2819-2.0;
                double t2821 = 1.0/t2820;
                double t2822 = 1.0/t2802;
                double t2823 = EcPhd_1*t2822*(1.0/3.0);
                double t2824 = log(t2805);
                double t2825 = EcPhd_3*c*t2803*(1.0/3.0);
                double t2826 = EcPhd_4*c*t2803*(1.0/3.0);
                double t2827 = EcPhd_3*c*t2803*t2824*(1.0/3.0);
                double t2828 = t2805-1.0;
                double t2829 = EcPhd_1*t2824;
                double t2830 = dirac(t2828);
                double t2831 = EcFld_2*t2807;
                double t2832 = EcFld_3*c*t2804;
                double t2833 = t2831+t2832+1.0;
                double t2834 = 1.0/t2813;
                double t2835 = EcPld_1*t2834;
                double t2836 = -t2805+1.0;
                double t2837 = heaviside(t2836);
                double t2838 = EcFhd_1*t2824;
                double t2839 = EcFhd_4*c*t2804;
                double t2840 = EcPhd_4*c*t2804;
                double t2841 = EcFhd_3*c*t2804*t2824;
                double t2842 = EcPhd_3*c*t2804*t2824;
                double t2843 = heaviside(t2828);
                double t2844 = 1.0/t2833;
                double t2845 = t2835-EcFld_1*t2844;
                double t2846 = t2835-t2821*t2818*t2845;
                v_rho_a[Q] += scale * (t2837*(EcPhd_2+t2840+t2842+t2829+t2821*t2818*(EcFhd_2-EcPhd_2-t2840+t2841-t2842-t2829+t2838+t2839))+t2843*t2846+t2802*(t2843*(t2815-t2821*t2818*(t2815-EcFld_1*1.0/(t2833*t2833)*(EcFld_3*c*t2803*(1.0/3.0)+EcFld_2*c*t2803*t2806*(1.0/6.0))))-t2837*(t2823+t2825+t2826+t2827-t2821*t2818*(t2823+t2825+t2826+t2827-EcFhd_1*t2822*(1.0/3.0)-EcFhd_3*c*t2803*(1.0/3.0)-EcFhd_4*c*t2803*(1.0/3.0)-EcFhd_3*c*t2803*t2824*(1.0/3.0)))-c*t2803*t2830*t2846*(1.0/3.0)+c*t2803*t2830*(EcPhd_2+t2840+t2842+t2829+t2821*t2818*(EcFhd_2-EcPhd_2+t2841-t2829+t2838+t2839-EcPhd_4*c*t2804-EcPhd_3*c*t2804*t2824))*(1.0/3.0)));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t2848 = rho_a+rho_b;
                double t2849 = 1.0/pow(t2848,4.0/3.0);
                double t2850 = 1.0/pow(t2848,1.0/3.0);
                double t2851 = c*t2850;
                double t2852 = 1.0/sqrt(t2851);
                double t2853 = sqrt(t2851);
                double t2854 = EcPld_3*c*t2849*(1.0/3.0);
                double t2855 = EcPld_2*c*t2852*t2849*(1.0/6.0);
                double t2856 = t2854+t2855;
                double t2857 = EcPld_2*t2853;
                double t2858 = EcPld_3*c*t2850;
                double t2859 = t2857+t2858+1.0;
                double t2860 = 1.0/(t2859*t2859);
                double t2861 = EcPld_1*t2860*t2856;
                double t2862 = pow(2.0,1.0/3.0);
                double t2863 = t2862*2.0;
                double t2864 = t2863-2.0;
                double t2865 = two_13*2.0;
                double t2866 = t2865-2.0;
                double t2867 = 1.0/t2866;
                double t2868 = 1.0/t2848;
                double t2869 = EcPhd_1*t2868*(1.0/3.0);
                double t2870 = log(t2851);
                double t2871 = EcPhd_3*c*t2849*(1.0/3.0);
                double t2872 = EcPhd_4*c*t2849*(1.0/3.0);
                double t2873 = EcPhd_3*c*t2870*t2849*(1.0/3.0);
                double t2874 = t2851-1.0;
                double t2875 = EcPhd_1*t2870;
                double t2876 = dirac(t2874);
                double t2877 = EcFld_2*t2853;
                double t2878 = EcFld_3*c*t2850;
                double t2879 = t2877+t2878+1.0;
                double t2880 = 1.0/t2859;
                double t2881 = EcPld_1*t2880;
                double t2882 = -t2851+1.0;
                double t2883 = heaviside(t2882);
                double t2884 = EcFhd_1*t2870;
                double t2885 = EcFhd_4*c*t2850;
                double t2886 = EcPhd_4*c*t2850;
                double t2887 = EcFhd_3*c*t2850*t2870;
                double t2888 = EcPhd_3*c*t2850*t2870;
                double t2889 = heaviside(t2874);
                double t2890 = 1.0/t2879;
                double t2891 = t2881-EcFld_1*t2890;
                double t2892 = t2881-t2864*t2891*t2867;
                v_rho_b[Q] += scale * (t2883*(EcPhd_2+t2875+t2886+t2888+t2864*t2867*(EcFhd_2-EcPhd_2-t2875+t2884+t2885-t2886+t2887-t2888))+t2892*t2889+t2848*(t2889*(t2861-t2864*t2867*(t2861-EcFld_1*1.0/(t2879*t2879)*(EcFld_3*c*t2849*(1.0/3.0)+EcFld_2*c*t2852*t2849*(1.0/6.0))))-t2883*(t2871+t2872+t2873+t2869-t2864*t2867*(t2871+t2872+t2873+t2869-EcFhd_1*t2868*(1.0/3.0)-EcFhd_3*c*t2849*(1.0/3.0)-EcFhd_4*c*t2849*(1.0/3.0)-EcFhd_3*c*t2870*t2849*(1.0/3.0)))-c*t2892*t2849*t2876*(1.0/3.0)+c*t2849*t2876*(EcPhd_2+t2875+t2886+t2888+t2864*t2867*(EcFhd_2-EcPhd_2-t2875+t2884+t2885+t2887-EcPhd_4*c*t2850-EcPhd_3*c*t2850*t2870))*(1.0/3.0)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t2899 = rho_a+rho_b;
                double t2900 = 1.0/pow(t2899,4.0/3.0);
                double t2901 = 1.0/pow(t2899,1.0/3.0);
                double t2902 = c*t2901;
                double t2903 = 1.0/sqrt(t2902);
                double t2904 = sqrt(t2902);
                double t2905 = EcPld_3*c*t2900*(1.0/3.0);
                double t2906 = EcPld_2*c*t2900*t2903*(1.0/6.0);
                double t2907 = t2905+t2906;
                double t2908 = EcPld_2*t2904;
                double t2909 = EcPld_3*c*t2901;
                double t2910 = t2908+t2909+1.0;
                double t2911 = 1.0/(t2910*t2910);
                double t2912 = EcPld_1*t2911*t2907;
                double t2913 = t2902-1.0;
                double t2914 = heaviside(t2913);
                double t2915 = pow(2.0,1.0/3.0);
                double t2916 = t2915*2.0;
                double t2917 = t2916-2.0;
                double t2918 = two_13*2.0;
                double t2919 = t2918-2.0;
                double t2920 = 1.0/t2919;
                double t2921 = EcFld_3*c*t2900*(1.0/3.0);
                double t2922 = EcFld_2*c*t2900*t2903*(1.0/6.0);
                double t2923 = t2921+t2922;
                double t2924 = EcFld_2*t2904;
                double t2925 = EcFld_3*c*t2901;
                double t2926 = t2924+t2925+1.0;
                double t2927 = t2907*t2907;
                double t2928 = 1.0/(t2910*t2910*t2910);
                double t2929 = EcPld_1*t2927*t2928*2.0;
                double t2930 = 1.0/pow(t2899,7.0/3.0);
                double t2931 = 1.0/(t2926*t2926);
                double t2932 = c*c;
                double t2933 = 1.0/pow(t2899,8.0/3.0);
                double t2934 = 1.0/pow(t2902,3.0/2.0);
                double t2935 = EcPld_3*c*t2930*(4.0/9.0);
                double t2936 = EcPld_2*c*t2903*t2930*(2.0/9.0);
                double t2937 = t2935+t2936-EcPld_2*t2932*t2933*t2934*(1.0/3.6E1);
                double t2938 = EcPld_1*t2911*t2937;
                double t2939 = 1.0/(t2899*t2899);
                double t2940 = EcPhd_1*t2939*(1.0/3.0);
                double t2941 = log(t2902);
                double t2942 = EcPhd_3*c*t2930*(5.0/9.0);
                double t2943 = EcPhd_4*c*t2930*(4.0/9.0);
                double t2944 = EcPhd_3*c*t2930*t2941*(4.0/9.0);
                double t2945 = 1.0/t2910;
                double t2946 = EcPld_1*t2945;
                double t2947 = t2912-EcFld_1*t2931*t2923;
                double t2948 = t2912-t2920*t2917*t2947;
                double t2949 = dirac(t2913);
                double t2950 = EcPhd_1*t2941;
                double t2951 = 1.0/t2899;
                double t2952 = EcPhd_1*t2951*(1.0/3.0);
                double t2953 = EcPhd_3*c*t2900*(1.0/3.0);
                double t2954 = EcPhd_4*c*t2900*(1.0/3.0);
                double t2955 = EcPhd_3*c*t2900*t2941*(1.0/3.0);
                double t2956 = 1.0/t2926;
                double t2972 = EcFld_1*t2956;
                double t2957 = -t2972+t2946;
                double t2958 = t2946-t2920*t2917*t2957;
                double t2959 = dirac(t2913,1.0);
                double t2960 = EcFhd_1*t2941;
                double t2961 = EcFhd_4*c*t2901;
                double t2962 = EcPhd_4*c*t2901;
                double t2963 = EcFhd_3*c*t2901*t2941;
                double t2964 = EcPhd_3*c*t2901*t2941;
                double t2965 = EcFhd_2-EcPhd_2-t2950+t2960+t2961-t2962+t2963-t2964;
                double t2966 = t2920*t2917*t2965;
                double t2967 = EcPhd_2+t2950+t2962+t2964+t2966;
                double t2968 = -t2902+1.0;
                double t2969 = heaviside(t2968);
                double t2970 = t2952+t2953+t2954+t2955-EcFhd_1*t2951*(1.0/3.0)-EcFhd_3*c*t2900*(1.0/3.0)-EcFhd_4*c*t2900*(1.0/3.0)-EcFhd_3*c*t2900*t2941*(1.0/3.0);
                double t2971 = t2952+t2953+t2954+t2955-t2920*t2970*t2917;
                v_rho_a_rho_a[Q] += scale * (-t2899*(-t2969*(t2940+t2942+t2943+t2944-t2920*t2917*(t2940+t2942+t2943+t2944-EcFhd_1*t2939*(1.0/3.0)-EcFhd_3*c*t2930*(5.0/9.0)-EcFhd_4*c*t2930*(4.0/9.0)-EcFhd_3*c*t2930*t2941*(4.0/9.0)))+t2914*(-t2929+t2938+t2920*t2917*(t2929-t2938-EcFld_1*(t2923*t2923)*1.0/(t2926*t2926*t2926)*2.0+EcFld_1*t2931*(EcFld_3*c*t2930*(4.0/9.0)-EcFld_2*t2932*t2933*t2934*(1.0/3.6E1)+EcFld_2*c*t2903*t2930*(2.0/9.0))))+c*t2900*t2971*t2949*(2.0/3.0)+c*t2900*t2948*t2949*(2.0/3.0)-c*t2930*t2949*t2958*(4.0/9.0)+c*t2930*t2949*t2967*(4.0/9.0)-t2932*t2933*t2958*t2959*(1.0/9.0)+t2932*t2933*t2967*t2959*(1.0/9.0))+t2914*t2948*2.0-t2971*t2969*2.0+c*t2900*t2949*t2967*(2.0/3.0)-c*t2900*t2949*(t2946+t2920*t2917*(t2972-t2946))*(2.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t2974 = rho_a+rho_b;
                double t2975 = 1.0/pow(t2974,4.0/3.0);
                double t2976 = 1.0/pow(t2974,1.0/3.0);
                double t2977 = c*t2976;
                double t2978 = 1.0/sqrt(t2977);
                double t2979 = sqrt(t2977);
                double t2980 = EcPld_3*c*t2975*(1.0/3.0);
                double t2981 = EcPld_2*c*t2975*t2978*(1.0/6.0);
                double t2982 = t2980+t2981;
                double t2983 = EcPld_2*t2979;
                double t2984 = EcPld_3*c*t2976;
                double t2985 = t2983+t2984+1.0;
                double t2986 = 1.0/(t2985*t2985);
                double t2987 = EcPld_1*t2982*t2986;
                double t2988 = t2977-1.0;
                double t2989 = heaviside(t2988);
                double t2990 = pow(2.0,1.0/3.0);
                double t2991 = t2990*2.0;
                double t2992 = t2991-2.0;
                double t2993 = two_13*2.0;
                double t2994 = t2993-2.0;
                double t2995 = 1.0/t2994;
                double t2996 = EcFld_3*c*t2975*(1.0/3.0);
                double t2997 = EcFld_2*c*t2975*t2978*(1.0/6.0);
                double t2998 = t2996+t2997;
                double t2999 = EcFld_2*t2979;
                double t3000 = EcFld_3*c*t2976;
                double t3001 = t2999+t3000+1.0;
                double t3002 = t2982*t2982;
                double t3003 = 1.0/(t2985*t2985*t2985);
                double t3004 = EcPld_1*t3002*t3003*2.0;
                double t3005 = 1.0/pow(t2974,7.0/3.0);
                double t3006 = 1.0/(t3001*t3001);
                double t3007 = c*c;
                double t3008 = 1.0/pow(t2974,8.0/3.0);
                double t3009 = 1.0/pow(t2977,3.0/2.0);
                double t3010 = EcPld_3*c*t3005*(4.0/9.0);
                double t3011 = EcPld_2*c*t2978*t3005*(2.0/9.0);
                double t3012 = t3010+t3011-EcPld_2*t3007*t3008*t3009*(1.0/3.6E1);
                double t3013 = EcPld_1*t2986*t3012;
                double t3014 = 1.0/(t2974*t2974);
                double t3015 = EcPhd_1*t3014*(1.0/3.0);
                double t3016 = log(t2977);
                double t3017 = EcPhd_3*c*t3005*(5.0/9.0);
                double t3018 = EcPhd_4*c*t3005*(4.0/9.0);
                double t3019 = EcPhd_3*c*t3005*t3016*(4.0/9.0);
                double t3020 = 1.0/t2985;
                double t3021 = EcPld_1*t3020;
                double t3022 = t2987-EcFld_1*t2998*t3006;
                double t3023 = t2987-t2992*t2995*t3022;
                double t3024 = dirac(t2988);
                double t3025 = EcPhd_1*t3016;
                double t3026 = 1.0/t2974;
                double t3027 = EcPhd_1*t3026*(1.0/3.0);
                double t3028 = EcPhd_3*c*t2975*(1.0/3.0);
                double t3029 = EcPhd_4*c*t2975*(1.0/3.0);
                double t3030 = EcPhd_3*c*t2975*t3016*(1.0/3.0);
                double t3031 = 1.0/t3001;
                double t3047 = EcFld_1*t3031;
                double t3032 = t3021-t3047;
                double t3048 = t2992*t2995*t3032;
                double t3033 = t3021-t3048;
                double t3034 = dirac(t2988,1.0);
                double t3035 = EcFhd_1*t3016;
                double t3036 = EcFhd_4*c*t2976;
                double t3037 = EcPhd_4*c*t2976;
                double t3038 = EcFhd_3*c*t2976*t3016;
                double t3039 = EcPhd_3*c*t2976*t3016;
                double t3040 = EcFhd_2-EcPhd_2-t3025+t3035+t3036-t3037+t3038-t3039;
                double t3041 = t2992*t2995*t3040;
                double t3042 = EcPhd_2+t3041+t3025+t3037+t3039;
                double t3043 = -t2977+1.0;
                double t3044 = heaviside(t3043);
                double t3045 = t3030+t3027+t3028+t3029-EcFhd_1*t3026*(1.0/3.0)-EcFhd_3*c*t2975*(1.0/3.0)-EcFhd_4*c*t2975*(1.0/3.0)-EcFhd_3*c*t2975*t3016*(1.0/3.0);
                double t3046 = t3030+t3027+t3028+t3029-t2992*t2995*t3045;
                v_rho_a_rho_b[Q] += scale * (-t2974*(-t3044*(t3015+t3017+t3018+t3019-t2992*t2995*(t3015+t3017+t3018+t3019-EcFhd_1*t3014*(1.0/3.0)-EcFhd_3*c*t3005*(5.0/9.0)-EcFhd_4*c*t3005*(4.0/9.0)-EcFhd_3*c*t3005*t3016*(4.0/9.0)))+t2989*(-t3004+t3013+t2992*t2995*(t3004-t3013-EcFld_1*(t2998*t2998)*1.0/(t3001*t3001*t3001)*2.0+EcFld_1*t3006*(EcFld_3*c*t3005*(4.0/9.0)-EcFld_2*t3007*t3008*t3009*(1.0/3.6E1)+EcFld_2*c*t2978*t3005*(2.0/9.0))))+c*t2975*t3023*t3024*(2.0/3.0)+c*t2975*t3024*t3046*(2.0/3.0)-c*t3005*t3024*t3033*(4.0/9.0)+c*t3005*t3024*t3042*(4.0/9.0)-t3033*t3007*t3034*t3008*(1.0/9.0)+t3042*t3007*t3034*t3008*(1.0/9.0))+t2989*t3023*2.0-t3044*t3046*2.0-c*t2975*t3024*t3033*(2.0/3.0)+c*t2975*t3024*t3042*(2.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t3050 = rho_a+rho_b;
                double t3051 = 1.0/pow(t3050,4.0/3.0);
                double t3052 = 1.0/pow(t3050,1.0/3.0);
                double t3053 = c*t3052;
                double t3054 = 1.0/sqrt(t3053);
                double t3055 = sqrt(t3053);
                double t3056 = EcPld_3*c*t3051*(1.0/3.0);
                double t3057 = EcPld_2*c*t3051*t3054*(1.0/6.0);
                double t3058 = t3056+t3057;
                double t3059 = EcPld_2*t3055;
                double t3060 = EcPld_3*c*t3052;
                double t3061 = t3060+t3059+1.0;
                double t3062 = 1.0/(t3061*t3061);
                double t3063 = EcPld_1*t3062*t3058;
                double t3064 = t3053-1.0;
                double t3065 = heaviside(t3064);
                double t3066 = pow(2.0,1.0/3.0);
                double t3067 = t3066*2.0;
                double t3068 = t3067-2.0;
                double t3069 = two_13*2.0;
                double t3070 = t3069-2.0;
                double t3071 = 1.0/t3070;
                double t3072 = EcFld_3*c*t3051*(1.0/3.0);
                double t3073 = EcFld_2*c*t3051*t3054*(1.0/6.0);
                double t3074 = t3072+t3073;
                double t3075 = EcFld_2*t3055;
                double t3076 = EcFld_3*c*t3052;
                double t3077 = t3075+t3076+1.0;
                double t3078 = t3058*t3058;
                double t3079 = 1.0/(t3061*t3061*t3061);
                double t3080 = EcPld_1*t3078*t3079*2.0;
                double t3081 = 1.0/pow(t3050,7.0/3.0);
                double t3082 = 1.0/(t3077*t3077);
                double t3083 = c*c;
                double t3084 = 1.0/pow(t3050,8.0/3.0);
                double t3085 = 1.0/pow(t3053,3.0/2.0);
                double t3086 = EcPld_3*c*t3081*(4.0/9.0);
                double t3087 = EcPld_2*c*t3054*t3081*(2.0/9.0);
                double t3088 = t3086+t3087-EcPld_2*t3083*t3084*t3085*(1.0/3.6E1);
                double t3089 = EcPld_1*t3062*t3088;
                double t3090 = 1.0/(t3050*t3050);
                double t3091 = EcPhd_1*t3090*(1.0/3.0);
                double t3092 = log(t3053);
                double t3093 = EcPhd_3*c*t3081*(5.0/9.0);
                double t3094 = EcPhd_4*c*t3081*(4.0/9.0);
                double t3095 = EcPhd_3*c*t3081*t3092*(4.0/9.0);
                double t3096 = 1.0/t3061;
                double t3097 = EcPld_1*t3096;
                double t3098 = t3063-EcFld_1*t3082*t3074;
                double t3099 = t3063-t3071*t3068*t3098;
                double t3100 = dirac(t3064);
                double t3101 = EcPhd_1*t3092;
                double t3102 = 1.0/t3050;
                double t3103 = EcPhd_1*t3102*(1.0/3.0);
                double t3104 = EcPhd_3*c*t3051*(1.0/3.0);
                double t3105 = EcPhd_4*c*t3051*(1.0/3.0);
                double t3106 = EcPhd_3*c*t3051*t3092*(1.0/3.0);
                double t3107 = 1.0/t3077;
                double t3123 = EcFld_1*t3107;
                double t3108 = t3097-t3123;
                double t3124 = t3071*t3068*t3108;
                double t3109 = t3097-t3124;
                double t3110 = dirac(t3064,1.0);
                double t3111 = EcFhd_1*t3092;
                double t3112 = EcFhd_4*c*t3052;
                double t3113 = EcPhd_4*c*t3052;
                double t3114 = EcFhd_3*c*t3052*t3092;
                double t3115 = EcPhd_3*c*t3052*t3092;
                double t3116 = EcFhd_2-EcPhd_2-t3101+t3111+t3112-t3113+t3114-t3115;
                double t3117 = t3071*t3068*t3116;
                double t3118 = EcPhd_2+t3101+t3113+t3115+t3117;
                double t3119 = -t3053+1.0;
                double t3120 = heaviside(t3119);
                double t3121 = t3103+t3104+t3105+t3106-EcFhd_1*t3102*(1.0/3.0)-EcFhd_3*c*t3051*(1.0/3.0)-EcFhd_4*c*t3051*(1.0/3.0)-EcFhd_3*c*t3051*t3092*(1.0/3.0);
                double t3122 = t3103+t3104+t3105+t3106-t3071*t3068*t3121;
                v_rho_b_rho_b[Q] += scale * (-t3050*(-t3120*(t3091+t3093+t3094+t3095-t3071*t3068*(t3091+t3093+t3094+t3095-EcFhd_1*t3090*(1.0/3.0)-EcFhd_3*c*t3081*(5.0/9.0)-EcFhd_4*c*t3081*(4.0/9.0)-EcFhd_3*c*t3081*t3092*(4.0/9.0)))+t3065*(-t3080+t3089+t3071*t3068*(t3080-t3089-EcFld_1*(t3074*t3074)*1.0/(t3077*t3077*t3077)*2.0+EcFld_1*t3082*(EcFld_3*c*t3081*(4.0/9.0)-EcFld_2*t3083*t3084*t3085*(1.0/3.6E1)+EcFld_2*c*t3054*t3081*(2.0/9.0))))+c*t3051*t3099*t3100*(2.0/3.0)+c*t3051*t3100*t3122*(2.0/3.0)-c*t3081*t3100*t3109*(4.0/9.0)+c*t3081*t3100*t3118*(4.0/9.0)-t3083*t3084*t3110*t3109*(1.0/9.0)+t3083*t3084*t3110*t3118*(1.0/9.0))+t3065*t3099*2.0-t3120*t3122*2.0-c*t3051*t3100*t3109*(2.0/3.0)+c*t3051*t3100*t3118*(2.0/3.0));
            }
            
        } else if (rho_b < lsda_cutoff_) {
            // v
            if (deriv >= 0) {
                double t3151 = rho_a+rho_b;
                double t3152 = 1.0/pow(t3151,1.0/3.0);
                double t3153 = c*t3152;
                double t3154 = log(t3153);
                double t3155 = EcPhd_1*t3154;
                double t3156 = pow(2.0,1.0/3.0);
                double t3157 = t3156*2.0;
                double t3158 = t3157-2.0;
                double t3159 = two_13*2.0;
                double t3160 = t3159-2.0;
                double t3161 = 1.0/t3160;
                double t3162 = sqrt(t3153);
                double t3163 = EcPld_2*t3162;
                double t3164 = EcPld_3*c*t3152;
                double t3165 = t3163+t3164+1.0;
                double t3166 = 1.0/t3165;
                double t3167 = EcPld_1*t3166;
                v[Q] += scale * (t3151*(heaviside(-c*t3152+1.0)*(EcPhd_2+t3155+t3161*t3158*(EcFhd_2-EcPhd_2-t3155+EcFhd_1*t3154+EcFhd_4*c*t3152-EcPhd_4*c*t3152+EcFhd_3*c*t3152*t3154-EcPhd_3*c*t3152*t3154)+EcPhd_4*c*t3152+EcPhd_3*c*t3152*t3154)+heaviside(t3153-1.0)*(t3167-t3161*t3158*(t3167-EcFld_1/(EcFld_2*t3162+EcFld_3*c*t3152+1.0)))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t3169 = rho_a+rho_b;
                double t3170 = 1.0/pow(t3169,4.0/3.0);
                double t3171 = 1.0/pow(t3169,1.0/3.0);
                double t3172 = c*t3171;
                double t3173 = 1.0/sqrt(t3172);
                double t3174 = sqrt(t3172);
                double t3175 = EcPld_3*c*t3170*(1.0/3.0);
                double t3176 = EcPld_2*c*t3170*t3173*(1.0/6.0);
                double t3177 = t3175+t3176;
                double t3178 = EcPld_2*t3174;
                double t3179 = EcPld_3*c*t3171;
                double t3180 = t3178+t3179+1.0;
                double t3181 = 1.0/(t3180*t3180);
                double t3182 = EcPld_1*t3181*t3177;
                double t3183 = pow(2.0,1.0/3.0);
                double t3184 = t3183*2.0;
                double t3185 = t3184-2.0;
                double t3186 = two_13*2.0;
                double t3187 = t3186-2.0;
                double t3188 = 1.0/t3187;
                double t3189 = 1.0/t3169;
                double t3190 = EcPhd_1*t3189*(1.0/3.0);
                double t3191 = log(t3172);
                double t3192 = EcPhd_3*c*t3170*(1.0/3.0);
                double t3193 = EcPhd_4*c*t3170*(1.0/3.0);
                double t3194 = EcPhd_3*c*t3170*t3191*(1.0/3.0);
                double t3195 = t3172-1.0;
                double t3196 = EcPhd_1*t3191;
                double t3197 = dirac(t3195);
                double t3198 = EcFld_2*t3174;
                double t3199 = EcFld_3*c*t3171;
                double t3200 = t3198+t3199+1.0;
                double t3201 = 1.0/t3180;
                double t3202 = EcPld_1*t3201;
                double t3203 = -t3172+1.0;
                double t3204 = heaviside(t3203);
                double t3205 = EcFhd_1*t3191;
                double t3206 = EcFhd_4*c*t3171;
                double t3207 = EcPhd_4*c*t3171;
                double t3208 = EcFhd_3*c*t3171*t3191;
                double t3209 = EcPhd_3*c*t3171*t3191;
                double t3210 = heaviside(t3195);
                double t3211 = 1.0/t3200;
                double t3212 = t3202-EcFld_1*t3211;
                double t3213 = t3202-t3185*t3188*t3212;
                v_rho_a[Q] += scale * (t3204*(EcPhd_2+t3196+t3207+t3209+t3185*t3188*(EcFhd_2-EcPhd_2-t3196+t3205+t3206-t3207+t3208-t3209))+t3210*t3213+t3169*(t3210*(t3182-t3185*t3188*(t3182-EcFld_1*1.0/(t3200*t3200)*(EcFld_3*c*t3170*(1.0/3.0)+EcFld_2*c*t3170*t3173*(1.0/6.0))))-t3204*(t3190+t3192+t3193+t3194-t3185*t3188*(t3190+t3192+t3193+t3194-EcFhd_1*t3189*(1.0/3.0)-EcFhd_3*c*t3170*(1.0/3.0)-EcFhd_4*c*t3170*(1.0/3.0)-EcFhd_3*c*t3170*t3191*(1.0/3.0)))-c*t3170*t3197*t3213*(1.0/3.0)+c*t3170*t3197*(EcPhd_2+t3196+t3207+t3209+t3185*t3188*(EcFhd_2-EcPhd_2-t3196+t3205+t3206+t3208-EcPhd_4*c*t3171-EcPhd_3*c*t3171*t3191))*(1.0/3.0)));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t3215 = rho_a+rho_b;
                double t3216 = 1.0/pow(t3215,4.0/3.0);
                double t3217 = 1.0/pow(t3215,1.0/3.0);
                double t3218 = c*t3217;
                double t3219 = 1.0/sqrt(t3218);
                double t3220 = sqrt(t3218);
                double t3221 = EcPld_3*c*t3216*(1.0/3.0);
                double t3222 = EcPld_2*c*t3216*t3219*(1.0/6.0);
                double t3223 = t3221+t3222;
                double t3224 = EcPld_2*t3220;
                double t3225 = EcPld_3*c*t3217;
                double t3226 = t3224+t3225+1.0;
                double t3227 = 1.0/(t3226*t3226);
                double t3228 = EcPld_1*t3223*t3227;
                double t3229 = pow(2.0,1.0/3.0);
                double t3230 = t3229*2.0;
                double t3231 = t3230-2.0;
                double t3232 = two_13*2.0;
                double t3233 = t3232-2.0;
                double t3234 = 1.0/t3233;
                double t3235 = 1.0/t3215;
                double t3236 = EcPhd_1*t3235*(1.0/3.0);
                double t3237 = log(t3218);
                double t3238 = EcPhd_3*c*t3216*(1.0/3.0);
                double t3239 = EcPhd_4*c*t3216*(1.0/3.0);
                double t3240 = EcPhd_3*c*t3216*t3237*(1.0/3.0);
                double t3241 = t3218-1.0;
                double t3242 = EcPhd_1*t3237;
                double t3243 = dirac(t3241);
                double t3244 = EcFld_2*t3220;
                double t3245 = EcFld_3*c*t3217;
                double t3246 = t3244+t3245+1.0;
                double t3247 = 1.0/t3226;
                double t3248 = EcPld_1*t3247;
                double t3249 = -t3218+1.0;
                double t3250 = heaviside(t3249);
                double t3251 = EcFhd_1*t3237;
                double t3252 = EcFhd_4*c*t3217;
                double t3253 = EcPhd_4*c*t3217;
                double t3254 = EcFhd_3*c*t3217*t3237;
                double t3255 = EcPhd_3*c*t3217*t3237;
                double t3256 = heaviside(t3241);
                double t3257 = 1.0/t3246;
                double t3258 = t3248-EcFld_1*t3257;
                double t3259 = t3248-t3231*t3234*t3258;
                v_rho_b[Q] += scale * (t3250*(EcPhd_2+t3242+t3253+t3255+t3231*t3234*(EcFhd_2-EcPhd_2-t3242+t3251+t3252-t3253+t3254-t3255))+t3256*t3259+t3215*(t3256*(t3228-t3231*t3234*(t3228-EcFld_1*1.0/(t3246*t3246)*(EcFld_3*c*t3216*(1.0/3.0)+EcFld_2*c*t3216*t3219*(1.0/6.0))))-t3250*(t3240+t3236+t3238+t3239-t3231*t3234*(t3240+t3236+t3238+t3239-EcFhd_1*t3235*(1.0/3.0)-EcFhd_3*c*t3216*(1.0/3.0)-EcFhd_4*c*t3216*(1.0/3.0)-EcFhd_3*c*t3216*t3237*(1.0/3.0)))-c*t3216*t3243*t3259*(1.0/3.0)+c*t3216*t3243*(EcPhd_2+t3242+t3253+t3255+t3231*t3234*(EcFhd_2-EcPhd_2-t3242+t3251+t3252+t3254-EcPhd_4*c*t3217-EcPhd_3*c*t3217*t3237))*(1.0/3.0)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t3266 = rho_a+rho_b;
                double t3267 = 1.0/pow(t3266,4.0/3.0);
                double t3268 = 1.0/pow(t3266,1.0/3.0);
                double t3269 = c*t3268;
                double t3270 = 1.0/sqrt(t3269);
                double t3271 = sqrt(t3269);
                double t3272 = EcPld_3*c*t3267*(1.0/3.0);
                double t3273 = EcPld_2*c*t3270*t3267*(1.0/6.0);
                double t3274 = t3272+t3273;
                double t3275 = EcPld_2*t3271;
                double t3276 = EcPld_3*c*t3268;
                double t3277 = t3275+t3276+1.0;
                double t3278 = 1.0/(t3277*t3277);
                double t3279 = EcPld_1*t3274*t3278;
                double t3280 = t3269-1.0;
                double t3281 = heaviside(t3280);
                double t3282 = pow(2.0,1.0/3.0);
                double t3283 = t3282*2.0;
                double t3284 = t3283-2.0;
                double t3285 = two_13*2.0;
                double t3286 = t3285-2.0;
                double t3287 = 1.0/t3286;
                double t3288 = EcFld_3*c*t3267*(1.0/3.0);
                double t3289 = EcFld_2*c*t3270*t3267*(1.0/6.0);
                double t3290 = t3288+t3289;
                double t3291 = EcFld_2*t3271;
                double t3292 = EcFld_3*c*t3268;
                double t3293 = t3291+t3292+1.0;
                double t3294 = t3274*t3274;
                double t3295 = 1.0/(t3277*t3277*t3277);
                double t3296 = EcPld_1*t3294*t3295*2.0;
                double t3297 = 1.0/pow(t3266,7.0/3.0);
                double t3298 = 1.0/(t3293*t3293);
                double t3299 = c*c;
                double t3300 = 1.0/pow(t3266,8.0/3.0);
                double t3301 = 1.0/pow(t3269,3.0/2.0);
                double t3302 = EcPld_3*c*t3297*(4.0/9.0);
                double t3303 = EcPld_2*c*t3270*t3297*(2.0/9.0);
                double t3304 = t3302+t3303-EcPld_2*t3299*t3300*t3301*(1.0/3.6E1);
                double t3305 = EcPld_1*t3278*t3304;
                double t3306 = 1.0/(t3266*t3266);
                double t3307 = EcPhd_1*t3306*(1.0/3.0);
                double t3308 = log(t3269);
                double t3309 = EcPhd_3*c*t3297*(5.0/9.0);
                double t3310 = EcPhd_4*c*t3297*(4.0/9.0);
                double t3311 = EcPhd_3*c*t3297*t3308*(4.0/9.0);
                double t3312 = 1.0/t3277;
                double t3313 = EcPld_1*t3312;
                double t3314 = t3279-EcFld_1*t3290*t3298;
                double t3315 = t3279-t3284*t3287*t3314;
                double t3316 = dirac(t3280);
                double t3317 = EcPhd_1*t3308;
                double t3318 = 1.0/t3266;
                double t3319 = EcPhd_1*t3318*(1.0/3.0);
                double t3320 = EcPhd_3*c*t3267*(1.0/3.0);
                double t3321 = EcPhd_4*c*t3267*(1.0/3.0);
                double t3322 = EcPhd_3*c*t3267*t3308*(1.0/3.0);
                double t3323 = 1.0/t3293;
                double t3339 = EcFld_1*t3323;
                double t3324 = t3313-t3339;
                double t3340 = t3284*t3287*t3324;
                double t3325 = t3313-t3340;
                double t3326 = dirac(t3280,1.0);
                double t3327 = EcFhd_1*t3308;
                double t3328 = EcFhd_4*c*t3268;
                double t3329 = EcPhd_4*c*t3268;
                double t3330 = EcFhd_3*c*t3268*t3308;
                double t3331 = EcPhd_3*c*t3268*t3308;
                double t3332 = EcFhd_2-EcPhd_2+t3330-t3331-t3317+t3327+t3328-t3329;
                double t3333 = t3284*t3287*t3332;
                double t3334 = EcPhd_2+t3331+t3333+t3317+t3329;
                double t3335 = -t3269+1.0;
                double t3336 = heaviside(t3335);
                double t3337 = t3320+t3321+t3322+t3319-EcFhd_1*t3318*(1.0/3.0)-EcFhd_3*c*t3267*(1.0/3.0)-EcFhd_4*c*t3267*(1.0/3.0)-EcFhd_3*c*t3267*t3308*(1.0/3.0);
                double t3338 = t3320+t3321+t3322+t3319-t3284*t3287*t3337;
                v_rho_a_rho_a[Q] += scale * (-t3266*(-t3336*(t3310+t3311+t3307+t3309-t3284*t3287*(t3310+t3311+t3307+t3309-EcFhd_1*t3306*(1.0/3.0)-EcFhd_3*c*t3297*(5.0/9.0)-EcFhd_4*c*t3297*(4.0/9.0)-EcFhd_3*c*t3297*t3308*(4.0/9.0)))+t3281*(-t3296+t3305+t3284*t3287*(t3296-t3305-EcFld_1*(t3290*t3290)*1.0/(t3293*t3293*t3293)*2.0+EcFld_1*t3298*(EcFld_3*c*t3297*(4.0/9.0)-EcFld_2*t3299*t3300*t3301*(1.0/3.6E1)+EcFld_2*c*t3270*t3297*(2.0/9.0))))+c*t3267*t3315*t3316*(2.0/3.0)-c*t3297*t3316*t3325*(4.0/9.0)+c*t3297*t3316*t3334*(4.0/9.0)+c*t3267*t3316*t3338*(2.0/3.0)-t3299*t3300*t3325*t3326*(1.0/9.0)+t3299*t3300*t3334*t3326*(1.0/9.0))+t3281*t3315*2.0-t3336*t3338*2.0-c*t3267*t3316*t3325*(2.0/3.0)+c*t3267*t3316*t3334*(2.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t3342 = rho_a+rho_b;
                double t3343 = 1.0/pow(t3342,4.0/3.0);
                double t3344 = 1.0/pow(t3342,1.0/3.0);
                double t3345 = c*t3344;
                double t3346 = 1.0/sqrt(t3345);
                double t3347 = sqrt(t3345);
                double t3348 = EcPld_3*c*t3343*(1.0/3.0);
                double t3349 = EcPld_2*c*t3343*t3346*(1.0/6.0);
                double t3350 = t3348+t3349;
                double t3351 = EcPld_2*t3347;
                double t3352 = EcPld_3*c*t3344;
                double t3353 = t3351+t3352+1.0;
                double t3354 = 1.0/(t3353*t3353);
                double t3355 = EcPld_1*t3350*t3354;
                double t3356 = t3345-1.0;
                double t3357 = heaviside(t3356);
                double t3358 = pow(2.0,1.0/3.0);
                double t3359 = t3358*2.0;
                double t3360 = t3359-2.0;
                double t3361 = two_13*2.0;
                double t3362 = t3361-2.0;
                double t3363 = 1.0/t3362;
                double t3364 = EcFld_3*c*t3343*(1.0/3.0);
                double t3365 = EcFld_2*c*t3343*t3346*(1.0/6.0);
                double t3366 = t3364+t3365;
                double t3367 = EcFld_2*t3347;
                double t3368 = EcFld_3*c*t3344;
                double t3369 = t3367+t3368+1.0;
                double t3370 = t3350*t3350;
                double t3371 = 1.0/(t3353*t3353*t3353);
                double t3372 = EcPld_1*t3370*t3371*2.0;
                double t3373 = 1.0/pow(t3342,7.0/3.0);
                double t3374 = 1.0/(t3369*t3369);
                double t3375 = c*c;
                double t3376 = 1.0/pow(t3342,8.0/3.0);
                double t3377 = 1.0/pow(t3345,3.0/2.0);
                double t3378 = EcPld_3*c*t3373*(4.0/9.0);
                double t3379 = EcPld_2*c*t3346*t3373*(2.0/9.0);
                double t3380 = t3378+t3379-EcPld_2*t3375*t3376*t3377*(1.0/3.6E1);
                double t3381 = EcPld_1*t3380*t3354;
                double t3382 = 1.0/(t3342*t3342);
                double t3383 = EcPhd_1*t3382*(1.0/3.0);
                double t3384 = log(t3345);
                double t3385 = EcPhd_3*c*t3373*(5.0/9.0);
                double t3386 = EcPhd_4*c*t3373*(4.0/9.0);
                double t3387 = EcPhd_3*c*t3373*t3384*(4.0/9.0);
                double t3388 = 1.0/t3353;
                double t3389 = EcPld_1*t3388;
                double t3390 = t3355-EcFld_1*t3374*t3366;
                double t3391 = t3355-t3360*t3363*t3390;
                double t3392 = dirac(t3356);
                double t3393 = EcPhd_1*t3384;
                double t3394 = 1.0/t3342;
                double t3395 = EcPhd_1*t3394*(1.0/3.0);
                double t3396 = EcPhd_3*c*t3343*(1.0/3.0);
                double t3397 = EcPhd_4*c*t3343*(1.0/3.0);
                double t3398 = EcPhd_3*c*t3343*t3384*(1.0/3.0);
                double t3399 = 1.0/t3369;
                double t3415 = EcFld_1*t3399;
                double t3400 = t3389-t3415;
                double t3416 = t3360*t3363*t3400;
                double t3401 = t3389-t3416;
                double t3402 = dirac(t3356,1.0);
                double t3403 = EcFhd_1*t3384;
                double t3404 = EcFhd_4*c*t3344;
                double t3405 = EcPhd_4*c*t3344;
                double t3406 = EcFhd_3*c*t3344*t3384;
                double t3407 = EcPhd_3*c*t3344*t3384;
                double t3408 = EcFhd_2-EcPhd_2-t3393+t3403+t3404-t3405+t3406-t3407;
                double t3409 = t3360*t3363*t3408;
                double t3410 = EcPhd_2+t3393+t3405+t3407+t3409;
                double t3411 = -t3345+1.0;
                double t3412 = heaviside(t3411);
                double t3413 = t3395+t3396+t3397+t3398-EcFhd_1*t3394*(1.0/3.0)-EcFhd_3*c*t3343*(1.0/3.0)-EcFhd_4*c*t3343*(1.0/3.0)-EcFhd_3*c*t3343*t3384*(1.0/3.0);
                double t3414 = t3395+t3396+t3397+t3398-t3360*t3363*t3413;
                v_rho_a_rho_b[Q] += scale * (-t3342*(-t3412*(t3383+t3385+t3386+t3387-t3360*t3363*(t3383+t3385+t3386+t3387-EcFhd_1*t3382*(1.0/3.0)-EcFhd_3*c*t3373*(5.0/9.0)-EcFhd_4*c*t3373*(4.0/9.0)-EcFhd_3*c*t3373*t3384*(4.0/9.0)))+t3357*(-t3372+t3381+t3360*t3363*(t3372-t3381-EcFld_1*(t3366*t3366)*1.0/(t3369*t3369*t3369)*2.0+EcFld_1*t3374*(EcFld_3*c*t3373*(4.0/9.0)-EcFld_2*t3375*t3376*t3377*(1.0/3.6E1)+EcFld_2*c*t3346*t3373*(2.0/9.0))))+c*t3343*t3391*t3392*(2.0/3.0)-c*t3373*t3392*t3401*(4.0/9.0)+c*t3373*t3392*t3410*(4.0/9.0)+c*t3343*t3392*t3414*(2.0/3.0)-t3375*t3376*t3401*t3402*(1.0/9.0)+t3375*t3376*t3410*t3402*(1.0/9.0))+t3391*t3357*2.0-t3412*t3414*2.0-c*t3343*t3392*t3401*(2.0/3.0)+c*t3343*t3392*t3410*(2.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t3418 = rho_a+rho_b;
                double t3419 = 1.0/pow(t3418,4.0/3.0);
                double t3420 = 1.0/pow(t3418,1.0/3.0);
                double t3421 = c*t3420;
                double t3422 = 1.0/sqrt(t3421);
                double t3423 = sqrt(t3421);
                double t3424 = EcPld_3*c*t3419*(1.0/3.0);
                double t3425 = EcPld_2*c*t3422*t3419*(1.0/6.0);
                double t3426 = t3424+t3425;
                double t3427 = EcPld_2*t3423;
                double t3428 = EcPld_3*c*t3420;
                double t3429 = t3427+t3428+1.0;
                double t3430 = 1.0/(t3429*t3429);
                double t3431 = EcPld_1*t3430*t3426;
                double t3432 = t3421-1.0;
                double t3433 = heaviside(t3432);
                double t3434 = pow(2.0,1.0/3.0);
                double t3435 = t3434*2.0;
                double t3436 = t3435-2.0;
                double t3437 = two_13*2.0;
                double t3438 = t3437-2.0;
                double t3439 = 1.0/t3438;
                double t3440 = EcFld_3*c*t3419*(1.0/3.0);
                double t3441 = EcFld_2*c*t3422*t3419*(1.0/6.0);
                double t3442 = t3440+t3441;
                double t3443 = EcFld_2*t3423;
                double t3444 = EcFld_3*c*t3420;
                double t3445 = t3443+t3444+1.0;
                double t3446 = t3426*t3426;
                double t3447 = 1.0/(t3429*t3429*t3429);
                double t3448 = EcPld_1*t3446*t3447*2.0;
                double t3449 = 1.0/pow(t3418,7.0/3.0);
                double t3450 = 1.0/(t3445*t3445);
                double t3451 = c*c;
                double t3452 = 1.0/pow(t3418,8.0/3.0);
                double t3453 = 1.0/pow(t3421,3.0/2.0);
                double t3454 = EcPld_3*c*t3449*(4.0/9.0);
                double t3455 = EcPld_2*c*t3422*t3449*(2.0/9.0);
                double t3456 = t3454+t3455-EcPld_2*t3451*t3452*t3453*(1.0/3.6E1);
                double t3457 = EcPld_1*t3430*t3456;
                double t3458 = 1.0/(t3418*t3418);
                double t3459 = EcPhd_1*t3458*(1.0/3.0);
                double t3460 = log(t3421);
                double t3461 = EcPhd_3*c*t3449*(5.0/9.0);
                double t3462 = EcPhd_4*c*t3449*(4.0/9.0);
                double t3463 = EcPhd_3*c*t3460*t3449*(4.0/9.0);
                double t3464 = 1.0/t3429;
                double t3465 = EcPld_1*t3464;
                double t3466 = t3431-EcFld_1*t3450*t3442;
                double t3467 = t3431-t3436*t3439*t3466;
                double t3468 = dirac(t3432);
                double t3469 = EcPhd_1*t3460;
                double t3470 = 1.0/t3418;
                double t3471 = EcPhd_1*t3470*(1.0/3.0);
                double t3472 = EcPhd_3*c*t3419*(1.0/3.0);
                double t3473 = EcPhd_4*c*t3419*(1.0/3.0);
                double t3474 = EcPhd_3*c*t3460*t3419*(1.0/3.0);
                double t3475 = 1.0/t3445;
                double t3491 = EcFld_1*t3475;
                double t3476 = -t3491+t3465;
                double t3477 = t3465-t3436*t3439*t3476;
                double t3478 = dirac(t3432,1.0);
                double t3479 = EcFhd_1*t3460;
                double t3480 = EcFhd_4*c*t3420;
                double t3481 = EcPhd_4*c*t3420;
                double t3482 = EcFhd_3*c*t3420*t3460;
                double t3483 = EcPhd_3*c*t3420*t3460;
                double t3484 = EcFhd_2-EcPhd_2+t3480-t3481+t3482-t3483-t3469+t3479;
                double t3485 = t3436*t3439*t3484;
                double t3486 = EcPhd_2+t3481+t3483+t3485+t3469;
                double t3487 = -t3421+1.0;
                double t3488 = heaviside(t3487);
                double t3489 = t3471+t3472+t3473+t3474-EcFhd_1*t3470*(1.0/3.0)-EcFhd_3*c*t3419*(1.0/3.0)-EcFhd_4*c*t3419*(1.0/3.0)-EcFhd_3*c*t3460*t3419*(1.0/3.0);
                double t3490 = t3471+t3472+t3473+t3474-t3436*t3439*t3489;
                v_rho_b_rho_b[Q] += scale * (-t3418*(-t3488*(t3461+t3462+t3463+t3459-t3436*t3439*(t3461+t3462+t3463+t3459-EcFhd_1*t3458*(1.0/3.0)-EcFhd_3*c*t3449*(5.0/9.0)-EcFhd_4*c*t3449*(4.0/9.0)-EcFhd_3*c*t3460*t3449*(4.0/9.0)))+t3433*(-t3448+t3457+t3436*t3439*(t3448-t3457-EcFld_1*(t3442*t3442)*1.0/(t3445*t3445*t3445)*2.0+EcFld_1*t3450*(EcFld_3*c*t3449*(4.0/9.0)-EcFld_2*t3451*t3452*t3453*(1.0/3.6E1)+EcFld_2*c*t3422*t3449*(2.0/9.0))))+c*t3490*t3419*t3468*(2.0/3.0)+c*t3419*t3467*t3468*(2.0/3.0)-c*t3449*t3468*t3477*(4.0/9.0)+c*t3449*t3468*t3486*(4.0/9.0)-t3451*t3452*t3477*t3478*(1.0/9.0)+t3451*t3452*t3486*t3478*(1.0/9.0))+t3433*t3467*2.0-t3490*t3488*2.0+c*t3419*t3468*t3486*(2.0/3.0)-c*t3419*t3468*(t3465+t3436*t3439*(t3491-t3465))*(2.0/3.0));
            }
            
        } else {
            // v
            if (deriv >= 0) {
                double t2304 = rho_a+rho_b;
                double t2305 = 1.0/pow(t2304,1.0/3.0);
                double t2306 = 1.0/t2304;
                double t2307 = rho_a-rho_b;
                double t2308 = t2306*t2307;
                double t2309 = c*t2305;
                double t2310 = log(t2309);
                double t2311 = EcPhd_1*t2310;
                double t2312 = two_13*2.0;
                double t2313 = t2312-2.0;
                double t2314 = 1.0/t2313;
                double t2315 = sqrt(t2309);
                double t2316 = EcPld_2*t2315;
                double t2317 = EcPld_3*c*t2305;
                double t2318 = t2316+t2317+1.0;
                double t2319 = 1.0/t2318;
                double t2320 = EcPld_1*t2319;
                double t2321 = t2308+1.0;
                double t2322 = pow(t2321,4.0/3.0);
                double t2323 = -t2308+1.0;
                double t2324 = pow(t2323,4.0/3.0);
                double t2325 = t2322+t2324-2.0;
                v[Q] += scale * (t2304*(heaviside(-c*t2305+1.0)*(EcPhd_2+t2311+t2314*t2325*(EcFhd_2-EcPhd_2-t2311+EcFhd_1*t2310+EcFhd_4*c*t2305-EcPhd_4*c*t2305+EcFhd_3*c*t2310*t2305-EcPhd_3*c*t2310*t2305)+EcPhd_4*c*t2305+EcPhd_3*c*t2310*t2305)+heaviside(t2309-1.0)*(t2320-t2314*t2325*(t2320-EcFld_1/(EcFld_2*t2315+EcFld_3*c*t2305+1.0)))));
            }
            
            // v_rho_a
            if (deriv >= 1) {
                double t2327 = rho_a+rho_b;
                double t2328 = 1.0/pow(t2327,1.0/3.0);
                double t2329 = 1.0/t2327;
                double t2330 = rho_a-rho_b;
                double t2331 = t2330*t2329;
                double t2332 = c*t2328;
                double t2333 = log(t2332);
                double t2334 = EcPhd_1*t2333;
                double t2335 = 1.0/pow(t2327,4.0/3.0);
                double t2336 = two_13*2.0;
                double t2337 = t2336-2.0;
                double t2338 = 1.0/t2337;
                double t2339 = t2331+1.0;
                double t2340 = pow(t2339,4.0/3.0);
                double t2341 = -t2331+1.0;
                double t2342 = pow(t2341,4.0/3.0);
                double t2343 = t2340+t2342-2.0;
                double t2344 = EcPhd_1*t2329*(1.0/3.0);
                double t2345 = EcPhd_3*c*t2335*(1.0/3.0);
                double t2346 = EcPhd_4*c*t2335*(1.0/3.0);
                double t2347 = 1.0/(t2327*t2327);
                double t2358 = t2330*t2347;
                double t2348 = t2329-t2358;
                double t2349 = EcFhd_1*t2333;
                double t2350 = EcFhd_4*c*t2328;
                double t2351 = EcPhd_4*c*t2328;
                double t2352 = EcFhd_3*c*t2333*t2328;
                double t2353 = EcPhd_3*c*t2333*t2328;
                double t2354 = EcPhd_3*c*t2333*t2335*(1.0/3.0);
                double t2355 = 1.0/sqrt(t2332);
                double t2356 = sqrt(t2332);
                double t2357 = pow(t2339,1.0/3.0);
                double t2359 = t2348*t2357*(4.0/3.0);
                double t2360 = pow(t2341,1.0/3.0);
                double t2361 = t2359-t2360*t2348*(4.0/3.0);
                double t2362 = EcFld_2*t2356;
                double t2363 = EcFld_3*c*t2328;
                double t2364 = t2362+t2363+1.0;
                double t2365 = EcPld_2*t2356;
                double t2366 = EcPld_3*c*t2328;
                double t2367 = t2365+t2366+1.0;
                double t2368 = EcPld_3*c*t2335*(1.0/3.0);
                double t2369 = EcPld_2*c*t2335*t2355*(1.0/6.0);
                double t2370 = t2368+t2369;
                double t2371 = 1.0/(t2367*t2367);
                double t2372 = t2332-1.0;
                double t2373 = EcFhd_2-EcPhd_2+t2350-t2351-t2334+t2352-t2353+t2349;
                double t2374 = dirac(t2372);
                double t2375 = 1.0/t2367;
                double t2376 = EcPld_1*t2375;
                double t2377 = 1.0/t2364;
                double t2380 = EcFld_1*t2377;
                double t2378 = -t2380+t2376;
                double t2379 = heaviside(t2372);
                v_rho_a[Q] += scale * (-t2327*(heaviside(-t2332+1.0)*(t2344+t2345+t2354+t2346-t2343*t2338*(t2344+t2345+t2354+t2346-EcFhd_1*t2329*(1.0/3.0)-EcFhd_3*c*t2335*(1.0/3.0)-EcFhd_4*c*t2335*(1.0/3.0)-EcFhd_3*c*t2333*t2335*(1.0/3.0))-t2361*t2373*t2338)-t2379*(t2343*t2338*(EcFld_1*1.0/(t2364*t2364)*(EcFld_3*c*t2335*(1.0/3.0)+EcFld_2*c*t2335*t2355*(1.0/6.0))-EcPld_1*t2370*t2371)+EcPld_1*t2370*t2371-t2361*t2338*t2378)+c*t2335*t2374*(t2376-t2343*t2338*t2378)*(1.0/3.0)-c*t2335*t2374*(EcPhd_2+t2351+t2334+t2353+t2343*t2373*t2338)*(1.0/3.0))+heaviside(-c*t2328+1.0)*(EcPhd_2+t2351+t2334+t2353+t2343*t2338*(EcFhd_2-EcPhd_2+t2350-t2334+t2352+t2349-EcPhd_4*c*t2328-EcPhd_3*c*t2333*t2328))+t2379*(t2376+t2343*t2338*(t2380-t2376)));
            }
            
            // v_rho_b
            if (deriv >= 1) {
                double t2382 = rho_a+rho_b;
                double t2383 = 1.0/pow(t2382,1.0/3.0);
                double t2384 = 1.0/t2382;
                double t2385 = rho_a-rho_b;
                double t2386 = t2384*t2385;
                double t2387 = c*t2383;
                double t2388 = log(t2387);
                double t2389 = EcPhd_1*t2388;
                double t2390 = 1.0/pow(t2382,4.0/3.0);
                double t2391 = two_13*2.0;
                double t2392 = t2391-2.0;
                double t2393 = 1.0/t2392;
                double t2394 = t2386+1.0;
                double t2395 = pow(t2394,4.0/3.0);
                double t2396 = -t2386+1.0;
                double t2397 = pow(t2396,4.0/3.0);
                double t2398 = t2395+t2397-2.0;
                double t2399 = EcPhd_1*t2384*(1.0/3.0);
                double t2400 = EcPhd_3*c*t2390*(1.0/3.0);
                double t2401 = EcPhd_4*c*t2390*(1.0/3.0);
                double t2402 = 1.0/(t2382*t2382);
                double t2403 = t2385*t2402;
                double t2404 = t2384+t2403;
                double t2405 = EcFhd_1*t2388;
                double t2406 = EcFhd_4*c*t2383;
                double t2407 = EcPhd_4*c*t2383;
                double t2408 = EcFhd_3*c*t2383*t2388;
                double t2409 = EcPhd_3*c*t2383*t2388;
                double t2410 = EcPhd_3*c*t2390*t2388*(1.0/3.0);
                double t2411 = 1.0/sqrt(t2387);
                double t2412 = sqrt(t2387);
                double t2413 = pow(t2394,1.0/3.0);
                double t2414 = t2404*t2413*(4.0/3.0);
                double t2415 = pow(t2396,1.0/3.0);
                double t2416 = t2414-t2404*t2415*(4.0/3.0);
                double t2417 = EcFld_2*t2412;
                double t2418 = EcFld_3*c*t2383;
                double t2419 = t2417+t2418+1.0;
                double t2420 = EcPld_2*t2412;
                double t2421 = EcPld_3*c*t2383;
                double t2422 = t2420+t2421+1.0;
                double t2423 = EcPld_3*c*t2390*(1.0/3.0);
                double t2424 = EcPld_2*c*t2390*t2411*(1.0/6.0);
                double t2425 = t2423+t2424;
                double t2426 = 1.0/(t2422*t2422);
                double t2427 = t2387-1.0;
                double t2428 = EcFhd_2-EcPhd_2-t2389+t2405+t2406-t2407+t2408-t2409;
                double t2429 = dirac(t2427);
                double t2430 = 1.0/t2422;
                double t2431 = EcPld_1*t2430;
                double t2432 = 1.0/t2419;
                double t2435 = EcFld_1*t2432;
                double t2433 = t2431-t2435;
                double t2434 = heaviside(t2427);
                double t2436 = t2431-t2393*t2398*t2433;
                v_rho_b[Q] += scale * (-t2382*(heaviside(-t2387+1.0)*(t2399+t2400+t2401+t2410-t2393*t2398*(t2399+t2400+t2401+t2410-EcFhd_1*t2384*(1.0/3.0)-EcFhd_3*c*t2390*(1.0/3.0)-EcFhd_4*c*t2390*(1.0/3.0)-EcFhd_3*c*t2390*t2388*(1.0/3.0))+t2393*t2416*t2428)-t2434*(t2393*t2398*(EcFld_1*1.0/(t2419*t2419)*(EcFld_3*c*t2390*(1.0/3.0)+EcFld_2*c*t2390*t2411*(1.0/6.0))-EcPld_1*t2425*t2426)+EcPld_1*t2425*t2426+t2393*t2433*t2416)+c*t2390*t2436*t2429*(1.0/3.0)-c*t2390*t2429*(EcPhd_2+t2389+t2407+t2409+t2393*t2398*t2428)*(1.0/3.0))+t2434*t2436+heaviside(-c*t2383+1.0)*(EcPhd_2+t2389+t2407+t2409+t2393*t2398*(EcFhd_2-EcPhd_2-t2389+t2405+t2406+t2408-EcPhd_4*c*t2383-EcPhd_3*c*t2383*t2388)));
            }
            
            // v_rho_a_rho_a
            if (deriv >= 2) {
                double t2443 = rho_a+rho_b;
                double t2444 = 1.0/pow(t2443,4.0/3.0);
                double t2445 = 1.0/pow(t2443,1.0/3.0);
                double t2446 = c*t2445;
                double t2453 = 1.0/sqrt(t2446);
                double t2455 = EcPld_3*c*t2444*(1.0/3.0);
                double t2456 = EcPld_2*c*t2444*t2453*(1.0/6.0);
                double t2447 = t2455+t2456;
                double t2448 = 1.0/t2443;
                double t2449 = rho_a-rho_b;
                double t2450 = t2448*t2449;
                double t2451 = 1.0/(t2443*t2443);
                double t2488 = t2451*t2449;
                double t2452 = t2448-t2488;
                double t2454 = sqrt(t2446);
                double t2457 = EcPld_2*t2454;
                double t2458 = EcPld_3*c*t2445;
                double t2459 = t2457+t2458+1.0;
                double t2460 = two_13*2.0;
                double t2461 = t2460-2.0;
                double t2462 = 1.0/t2461;
                double t2463 = t2450+1.0;
                double t2464 = -t2450+1.0;
                double t2465 = EcFld_3*c*t2444*(1.0/3.0);
                double t2466 = EcFld_2*c*t2444*t2453*(1.0/6.0);
                double t2467 = t2465+t2466;
                double t2468 = EcFld_2*t2454;
                double t2469 = EcFld_3*c*t2445;
                double t2470 = t2468+t2469+1.0;
                double t2471 = t2447*t2447;
                double t2472 = 1.0/(t2459*t2459*t2459);
                double t2473 = EcPld_1*t2471*t2472*2.0;
                double t2474 = 1.0/pow(t2443,7.0/3.0);
                double t2475 = 1.0/(t2470*t2470);
                double t2476 = c*c;
                double t2477 = 1.0/pow(t2443,8.0/3.0);
                double t2478 = 1.0/pow(t2446,3.0/2.0);
                double t2479 = 1.0/(t2459*t2459);
                double t2480 = EcPld_3*c*t2474*(4.0/9.0);
                double t2481 = EcPld_2*c*t2453*t2474*(2.0/9.0);
                double t2482 = t2480+t2481-EcPld_2*t2476*t2477*t2478*(1.0/3.6E1);
                double t2483 = pow(t2463,1.0/3.0);
                double t2484 = pow(t2464,1.0/3.0);
                double t2485 = t2451*2.0;
                double t2486 = 1.0/(t2443*t2443*t2443);
                double t2490 = t2449*t2486*2.0;
                double t2487 = -t2490+t2485;
                double t2489 = t2452*t2452;
                double t2491 = t2484*t2487*(4.0/3.0);
                double t2492 = 1.0/pow(t2463,2.0/3.0);
                double t2493 = t2492*t2489*(4.0/9.0);
                double t2494 = 1.0/pow(t2464,2.0/3.0);
                double t2495 = t2494*t2489*(4.0/9.0);
                double t2496 = t2491+t2493+t2495-t2483*t2487*(4.0/3.0);
                double t2497 = log(t2446);
                double t2498 = pow(t2463,4.0/3.0);
                double t2499 = pow(t2464,4.0/3.0);
                double t2500 = t2498+t2499-2.0;
                double t2501 = EcPhd_1*t2451*(1.0/3.0);
                double t2502 = EcPhd_3*c*t2474*(5.0/9.0);
                double t2503 = EcPhd_4*c*t2474*(4.0/9.0);
                double t2504 = t2452*t2483*(4.0/3.0);
                double t2524 = t2452*t2484*(4.0/3.0);
                double t2505 = t2504-t2524;
                double t2506 = EcPhd_3*c*t2474*t2497*(4.0/9.0);
                double t2507 = t2446-1.0;
                double t2508 = 1.0/t2459;
                double t2509 = EcPld_1*t2508;
                double t2510 = 1.0/t2470;
                double t2529 = EcFld_1*t2510;
                double t2511 = t2509-t2529;
                double t2512 = EcFhd_1*t2497;
                double t2513 = EcPhd_1*t2497;
                double t2514 = EcFhd_4*c*t2445;
                double t2515 = EcFhd_3*c*t2445*t2497;
                double t2516 = dirac(t2507);
                double t2517 = EcFhd_1*t2448*(1.0/3.0);
                double t2518 = EcPhd_1*t2448*(1.0/3.0);
                double t2519 = EcFhd_3*c*t2444*(1.0/3.0);
                double t2520 = EcFhd_4*c*t2444*(1.0/3.0);
                double t2521 = EcPhd_3*c*t2444*(1.0/3.0);
                double t2522 = EcPhd_4*c*t2444*(1.0/3.0);
                double t2523 = EcFhd_3*c*t2444*t2497*(1.0/3.0);
                double t2525 = EcPhd_4*c*t2445;
                double t2526 = EcPhd_3*c*t2445*t2497;
                double t2527 = EcFld_1*t2475*t2467;
                double t2530 = EcPld_1*t2447*t2479;
                double t2528 = -t2530+t2527;
                double t2546 = t2462*t2500*t2511;
                double t2531 = t2509-t2546;
                double t2532 = dirac(t2507,1.0);
                double t2533 = EcFhd_2-EcPhd_2+t2512-t2513+t2514+t2515-t2525-t2526;
                double t2534 = EcPld_1*t2482*t2479;
                double t2535 = t2462*t2500*t2533;
                double t2536 = EcPhd_2+t2513+t2525+t2526+t2535;
                double t2538 = EcPhd_3*c*t2444*t2497*(1.0/3.0);
                double t2537 = t2520-t2521-t2522+t2523+t2517-t2518+t2519-t2538;
                double t2539 = -t2446+1.0;
                double t2540 = heaviside(t2539);
                double t2541 = t2462*t2500*t2537;
                double t2542 = t2521+t2522+t2541+t2518+t2538-t2462*t2505*t2533;
                double t2543 = heaviside(t2507);
                double t2544 = t2462*t2500*t2528;
                double t2545 = t2530+t2544-t2462*t2511*t2505;
                v_rho_a_rho_a[Q] += scale * (-t2443*(t2543*(-t2473+t2534+t2462*t2500*(t2473-t2534-EcFld_1*1.0/(t2470*t2470*t2470)*(t2467*t2467)*2.0+EcFld_1*t2475*(EcFld_3*c*t2474*(4.0/9.0)-EcFld_2*t2476*t2477*t2478*(1.0/3.6E1)+EcFld_2*c*t2453*t2474*(2.0/9.0)))+t2462*t2496*t2511-t2462*t2505*t2528*2.0)-t2540*(t2501+t2502+t2503+t2506-t2462*t2500*(t2501+t2502+t2503+t2506-EcFhd_1*t2451*(1.0/3.0)-EcFhd_3*c*t2474*(5.0/9.0)-EcFhd_4*c*t2474*(4.0/9.0)-EcFhd_3*c*t2474*t2497*(4.0/9.0))+t2462*t2496*t2533-t2462*t2505*t2537*2.0)+c*t2444*t2542*t2516*(2.0/3.0)-c*t2474*t2531*t2516*(4.0/9.0)+c*t2444*t2516*t2545*(2.0/3.0)+c*t2474*t2516*t2536*(4.0/9.0)-t2476*t2477*t2531*t2532*(1.0/9.0)+t2476*t2477*t2532*t2536*(1.0/9.0))-t2540*t2542*2.0+t2543*t2545*2.0-c*t2444*t2531*t2516*(2.0/3.0)+c*t2444*t2516*t2536*(2.0/3.0));
            }
            
            // v_rho_a_rho_b
            if (deriv >= 2) {
                double t2548 = rho_a+rho_b;
                double t2549 = rho_a-rho_b;
                double t2550 = 1.0/t2548;
                double t2551 = t2550*t2549;
                double t2552 = 1.0/(t2548*t2548*t2548);
                double t2553 = t2551+1.0;
                double t2554 = 1.0/(t2548*t2548);
                double t2555 = t2554*t2549;
                double t2556 = -t2551+1.0;
                double t2557 = t2550+t2555;
                double t2558 = t2550-t2555;
                double t2559 = 1.0/pow(t2548,1.0/3.0);
                double t2560 = c*t2559;
                double t2561 = log(t2560);
                double t2562 = 1.0/pow(t2548,7.0/3.0);
                double t2563 = two_13*2.0;
                double t2564 = t2563-2.0;
                double t2565 = 1.0/t2564;
                double t2566 = EcPhd_1*t2554*(1.0/3.0);
                double t2567 = EcPhd_3*c*t2562*(5.0/9.0);
                double t2568 = EcPhd_4*c*t2562*(4.0/9.0);
                double t2569 = pow(t2553,1.0/3.0);
                double t2570 = pow(t2556,1.0/3.0);
                double t2571 = 1.0/pow(t2548,4.0/3.0);
                double t2572 = EcFhd_1*t2550*(1.0/3.0);
                double t2573 = EcFhd_3*c*t2571*(1.0/3.0);
                double t2574 = EcFhd_4*c*t2571*(1.0/3.0);
                double t2575 = EcFhd_3*c*t2561*t2571*(1.0/3.0);
                double t2627 = EcPhd_1*t2550*(1.0/3.0);
                double t2628 = EcPhd_3*c*t2571*(1.0/3.0);
                double t2629 = EcPhd_4*c*t2571*(1.0/3.0);
                double t2630 = EcPhd_3*c*t2561*t2571*(1.0/3.0);
                double t2576 = t2572+t2573+t2574+t2575-t2630-t2627-t2628-t2629;
                double t2577 = EcPhd_3*c*t2561*t2562*(4.0/9.0);
                double t2580 = 1.0/sqrt(t2560);
                double t2582 = EcPld_3*c*t2571*(1.0/3.0);
                double t2583 = EcPld_2*c*t2571*t2580*(1.0/6.0);
                double t2578 = t2582+t2583;
                double t2579 = t2570*t2557*(4.0/3.0);
                double t2581 = sqrt(t2560);
                double t2584 = EcPld_2*t2581;
                double t2585 = EcPld_3*c*t2559;
                double t2586 = t2584+t2585+1.0;
                double t2587 = t2570*t2558*(4.0/3.0);
                double t2588 = EcFld_3*c*t2571*(1.0/3.0);
                double t2589 = EcFld_2*c*t2571*t2580*(1.0/6.0);
                double t2590 = t2588+t2589;
                double t2591 = EcFld_2*t2581;
                double t2592 = EcFld_3*c*t2559;
                double t2593 = t2591+t2592+1.0;
                double t2594 = 1.0/(t2593*t2593);
                double t2595 = EcFld_1*t2590*t2594;
                double t2596 = 1.0/(t2586*t2586);
                double t2637 = EcPld_1*t2578*t2596;
                double t2597 = t2595-t2637;
                double t2635 = t2558*t2569*(4.0/3.0);
                double t2598 = t2587-t2635;
                double t2599 = pow(t2553,4.0/3.0);
                double t2600 = pow(t2556,4.0/3.0);
                double t2601 = t2599+t2600-2.0;
                double t2602 = t2578*t2578;
                double t2603 = 1.0/(t2586*t2586*t2586);
                double t2604 = EcPld_1*t2602*t2603*2.0;
                double t2605 = c*c;
                double t2606 = 1.0/pow(t2548,8.0/3.0);
                double t2607 = 1.0/pow(t2560,3.0/2.0);
                double t2608 = EcPld_3*c*t2562*(4.0/9.0);
                double t2609 = EcPld_2*c*t2562*t2580*(2.0/9.0);
                double t2610 = t2608+t2609-EcPld_2*t2605*t2606*t2607*(1.0/3.6E1);
                double t2611 = t2552*t2570*t2549*(8.0/3.0);
                double t2612 = 1.0/pow(t2553,2.0/3.0);
                double t2613 = t2557*t2558*t2612*(4.0/9.0);
                double t2614 = 1.0/pow(t2556,2.0/3.0);
                double t2615 = t2557*t2558*t2614*(4.0/9.0);
                double t2616 = t2611+t2613+t2615-t2552*t2549*t2569*(8.0/3.0);
                double t2617 = t2560-1.0;
                double t2618 = 1.0/t2586;
                double t2619 = EcPld_1*t2618;
                double t2620 = 1.0/t2593;
                double t2639 = EcFld_1*t2620;
                double t2621 = t2619-t2639;
                double t2622 = EcFhd_1*t2561;
                double t2623 = EcPhd_1*t2561;
                double t2624 = EcFhd_4*c*t2559;
                double t2625 = EcFhd_3*c*t2561*t2559;
                double t2626 = dirac(t2617);
                double t2631 = EcPhd_4*c*t2559;
                double t2632 = EcPhd_3*c*t2561*t2559;
                double t2638 = t2557*t2569*(4.0/3.0);
                double t2633 = t2579-t2638;
                double t2634 = t2565*t2576*t2601;
                double t2636 = EcFhd_2-EcPhd_2+t2622-t2631-t2623-t2632+t2624+t2625;
                double t2640 = t2565*t2597*t2601;
                double t2652 = t2565*t2601*t2621;
                double t2641 = -t2652+t2619;
                double t2642 = dirac(t2617,1.0);
                double t2643 = t2565*t2601*t2636;
                double t2644 = EcPhd_2+t2631+t2623+t2632+t2643;
                double t2645 = -t2560+1.0;
                double t2646 = heaviside(t2645);
                double t2647 = t2630+t2634+t2627+t2628+t2629-t2565*t2633*t2636;
                double t2648 = t2565*t2598*t2636;
                double t2649 = t2630+t2634+t2627+t2628+t2629+t2648;
                double t2650 = heaviside(t2617);
                double t2651 = t2640+t2637-t2565*t2621*t2633;
                v_rho_a_rho_b[Q] += scale * (t2650*(t2640+t2637+t2565*(t2587-t2635)*(t2619-t2639))+t2650*t2651-t2646*t2647-t2646*t2649-t2548*(-t2650*(t2604-EcPld_1*t2596*t2610-t2565*t2597*t2598+t2565*t2597*t2633+t2565*t2621*t2616-t2565*t2601*(t2604-EcFld_1*(t2590*t2590)*1.0/(t2593*t2593*t2593)*2.0-EcPld_1*t2596*t2610+EcFld_1*t2594*(EcFld_3*c*t2562*(4.0/9.0)-EcFld_2*t2605*t2606*t2607*(1.0/3.6E1)+EcFld_2*c*t2562*t2580*(2.0/9.0))))-t2646*(t2566+t2567+t2568+t2577-t2565*t2601*(t2566+t2567+t2568+t2577-EcFhd_1*t2554*(1.0/3.0)-EcFhd_3*c*t2562*(5.0/9.0)-EcFhd_4*c*t2562*(4.0/9.0)-EcFhd_3*c*t2561*t2562*(4.0/9.0))+t2565*t2576*t2598-t2565*t2576*t2633-t2565*t2616*t2636)-c*t2562*t2641*t2626*(4.0/9.0)+c*t2571*t2651*t2626*(1.0/3.0)+c*t2562*t2626*t2644*(4.0/9.0)+c*t2571*t2626*t2647*(1.0/3.0)+c*t2571*t2626*t2649*(1.0/3.0)-t2605*t2641*t2606*t2642*(1.0/9.0)+t2605*t2606*t2642*t2644*(1.0/9.0)+c*t2571*t2626*(t2640+t2637+t2565*t2598*t2621)*(1.0/3.0))-c*t2571*t2641*t2626*(2.0/3.0)+c*t2571*t2626*t2644*(2.0/3.0));
            }
            
            // v_rho_b_rho_b
            if (deriv >= 2) {
                double t2654 = rho_a+rho_b;
                double t2655 = 1.0/pow(t2654,4.0/3.0);
                double t2656 = 1.0/pow(t2654,1.0/3.0);
                double t2657 = c*t2656;
                double t2665 = 1.0/sqrt(t2657);
                double t2667 = EcPld_3*c*t2655*(1.0/3.0);
                double t2668 = EcPld_2*c*t2655*t2665*(1.0/6.0);
                double t2658 = t2667+t2668;
                double t2659 = 1.0/t2654;
                double t2660 = rho_a-rho_b;
                double t2661 = t2660*t2659;
                double t2662 = 1.0/(t2654*t2654);
                double t2663 = t2660*t2662;
                double t2664 = t2663+t2659;
                double t2666 = sqrt(t2657);
                double t2669 = EcPld_2*t2666;
                double t2670 = EcPld_3*c*t2656;
                double t2671 = t2670+t2669+1.0;
                double t2672 = two_13*2.0;
                double t2673 = t2672-2.0;
                double t2674 = 1.0/t2673;
                double t2675 = t2661+1.0;
                double t2676 = -t2661+1.0;
                double t2677 = EcFld_3*c*t2655*(1.0/3.0);
                double t2678 = EcFld_2*c*t2655*t2665*(1.0/6.0);
                double t2679 = t2677+t2678;
                double t2680 = EcFld_2*t2666;
                double t2681 = EcFld_3*c*t2656;
                double t2682 = t2680+t2681+1.0;
                double t2683 = t2658*t2658;
                double t2684 = 1.0/(t2671*t2671*t2671);
                double t2685 = EcPld_1*t2683*t2684*2.0;
                double t2686 = 1.0/pow(t2654,7.0/3.0);
                double t2687 = 1.0/(t2682*t2682);
                double t2688 = c*c;
                double t2689 = 1.0/pow(t2654,8.0/3.0);
                double t2690 = 1.0/pow(t2657,3.0/2.0);
                double t2691 = 1.0/(t2671*t2671);
                double t2692 = EcPld_3*c*t2686*(4.0/9.0);
                double t2693 = EcPld_2*c*t2665*t2686*(2.0/9.0);
                double t2694 = t2692+t2693-EcPld_2*t2690*t2688*t2689*(1.0/3.6E1);
                double t2695 = pow(t2675,1.0/3.0);
                double t2696 = pow(t2676,1.0/3.0);
                double t2697 = t2662*2.0;
                double t2698 = 1.0/(t2654*t2654*t2654);
                double t2699 = t2660*t2698*2.0;
                double t2700 = t2697+t2699;
                double t2701 = t2664*t2664;
                double t2702 = t2695*t2700*(4.0/3.0);
                double t2703 = 1.0/pow(t2675,2.0/3.0);
                double t2704 = t2701*t2703*(4.0/9.0);
                double t2705 = 1.0/pow(t2676,2.0/3.0);
                double t2706 = t2701*t2705*(4.0/9.0);
                double t2707 = t2702+t2704+t2706-t2696*t2700*(4.0/3.0);
                double t2708 = log(t2657);
                double t2709 = pow(t2675,4.0/3.0);
                double t2710 = pow(t2676,4.0/3.0);
                double t2711 = t2710+t2709-2.0;
                double t2712 = EcPhd_1*t2662*(1.0/3.0);
                double t2713 = EcPhd_3*c*t2686*(5.0/9.0);
                double t2714 = EcPhd_4*c*t2686*(4.0/9.0);
                double t2715 = t2664*t2695*(4.0/3.0);
                double t2735 = t2664*t2696*(4.0/3.0);
                double t2716 = t2715-t2735;
                double t2717 = EcPhd_3*c*t2686*t2708*(4.0/9.0);
                double t2718 = t2657-1.0;
                double t2719 = 1.0/t2671;
                double t2720 = EcPld_1*t2719;
                double t2721 = 1.0/t2682;
                double t2740 = EcFld_1*t2721;
                double t2722 = t2720-t2740;
                double t2723 = EcFhd_1*t2708;
                double t2724 = EcPhd_1*t2708;
                double t2725 = EcFhd_4*c*t2656;
                double t2726 = EcFhd_3*c*t2656*t2708;
                double t2727 = dirac(t2718);
                double t2728 = EcFhd_1*t2659*(1.0/3.0);
                double t2729 = EcPhd_1*t2659*(1.0/3.0);
                double t2730 = EcFhd_3*c*t2655*(1.0/3.0);
                double t2731 = EcFhd_4*c*t2655*(1.0/3.0);
                double t2732 = EcPhd_3*c*t2655*(1.0/3.0);
                double t2733 = EcPhd_4*c*t2655*(1.0/3.0);
                double t2734 = EcFhd_3*c*t2655*t2708*(1.0/3.0);
                double t2736 = EcPhd_4*c*t2656;
                double t2737 = EcPhd_3*c*t2656*t2708;
                double t2738 = EcFld_1*t2687*t2679;
                double t2741 = EcPld_1*t2691*t2658;
                double t2739 = -t2741+t2738;
                double t2757 = t2674*t2711*t2722;
                double t2742 = t2720-t2757;
                double t2743 = dirac(t2718,1.0);
                double t2744 = EcFhd_2-EcPhd_2+t2723-t2724+t2725+t2726-t2736-t2737;
                double t2745 = EcPld_1*t2691*t2694;
                double t2746 = t2674*t2711*t2744;
                double t2747 = EcPhd_2+t2724+t2736+t2737+t2746;
                double t2749 = EcPhd_3*c*t2655*t2708*(1.0/3.0);
                double t2748 = t2730+t2731-t2732-t2733+t2734+t2728-t2729-t2749;
                double t2750 = -t2657+1.0;
                double t2751 = heaviside(t2750);
                double t2752 = t2674*t2711*t2748;
                double t2753 = t2674*t2716*t2744;
                double t2754 = t2732+t2733+t2752+t2753+t2729+t2749;
                double t2755 = heaviside(t2718);
                double t2756 = t2674*t2711*t2739;
                v_rho_b_rho_b[Q] += scale * (t2755*(t2741+t2756+t2674*t2716*(t2720-t2740))*2.0-t2751*t2754*2.0-t2654*(t2755*(-t2685+t2745+t2674*t2711*(t2685-t2745-EcFld_1*1.0/(t2682*t2682*t2682)*(t2679*t2679)*2.0+EcFld_1*t2687*(EcFld_3*c*t2686*(4.0/9.0)-EcFld_2*t2690*t2688*t2689*(1.0/3.6E1)+EcFld_2*c*t2665*t2686*(2.0/9.0)))+t2674*t2722*t2707+t2674*t2716*t2739*2.0)-t2751*(t2712+t2713+t2714+t2717-t2674*t2711*(t2712+t2713+t2714+t2717-EcFhd_1*t2662*(1.0/3.0)-EcFhd_3*c*t2686*(5.0/9.0)-EcFhd_4*c*t2686*(4.0/9.0)-EcFhd_3*c*t2686*t2708*(4.0/9.0))+t2674*t2707*t2744+t2674*t2716*t2748*2.0)+c*t2655*t2727*t2754*(2.0/3.0)-c*t2686*t2742*t2727*(4.0/9.0)+c*t2686*t2727*t2747*(4.0/9.0)-t2688*t2689*t2742*t2743*(1.0/9.0)+t2688*t2689*t2743*t2747*(1.0/9.0)+c*t2655*t2727*(t2741+t2756+t2674*t2722*t2716)*(2.0/3.0))-c*t2655*t2742*t2727*(2.0/3.0)+c*t2655*t2727*t2747*(2.0/3.0));
            }
            
        }
    }
}

}