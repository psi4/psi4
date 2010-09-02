/**********************************************************
* S_functional.cc: definiitions for S_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include <libmints/properties.h>
#include <libciomr/libciomr.h>
#include "S_functional.h"
#include <stdlib.h>
#include <cmath>
#include <string>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

S_Functional::S_Functional(int npoints, int deriv) : Functional(npoints, deriv)
{
}
S_Functional::~S_Functional()
{
}
void S_Functional::init()
{
    
    name_ = "S";
    description_ = "Simple Slater LSDA exchange functional";
    citation_ = "J.C. Slater, A simplification of the hartree-fock method, Phys. Rev., 81(3):385-390, 1951";
    
    double c = 9.305257363490997E-01;
    params_.push_back(make_pair("c",c));

    is_gga_ = false;
    is_meta_ = false;
    is_hybrid_ = false;

    exact_exchange_ = 0.0;

}
void S_Functional::computeRKSFunctional(shared_ptr<Properties> prop)
{
    int ntrue = prop->getTrueSize();

    const double* rho_a;
    const double* gamma_aa;
    const double* tau_a;

    rho_a = prop->getDensity();
    if (is_gga_) {
        gamma_aa = prop->getDensityGradientSquared();
    }
    if (is_meta_) {
        tau_a = prop->getKEDensity();
    }

    double c = params_[0].second;
    
    //Functional
    for (int index = 0; index < ntrue; index++) {
        functional_[index] = -c*pow(rho_a[index],4.0/3.0); 
    }
    //First Partials
    for (int index = 0; index < ntrue && deriv_ >= 1; index++) {
        v_rho_a_[index] = -4.0/3.0*c*pow(rho_a[index],1.0/3.0); 
        v_rho_b_[index] = v_rho_a_[index];
        if (is_gga_) {
            //Nothing for S
        }
    }
    //Second Partials
    for (int index = 0; index < ntrue && deriv_ >= 2; index++) {
        v_rho_a_rho_a_[index] = -4.0/9.0*c*pow(rho_a[index],-2.0/3.0); 
        v_rho_b_rho_b_[index] = v_rho_a_rho_a_[index];
        v_rho_a_rho_b_[index] = 0.0;
        if (is_gga_) {
            //Nothing for S
        }
    }
}
void S_Functional::computeUKSFunctional(shared_ptr<Properties> prop)
{
    int ntrue = prop->getTrueSize();

    const double* rho_a;
    const double* rho_b;
    const double* gamma_aa;
    const double* gamma_ab;
    const double* gamma_bb;
    const double* tau_a;
    const double* tau_b;

    rho_a = prop->getDensityA();
    rho_a = prop->getDensityB();
    if (is_gga_) {
        gamma_aa = prop->getDensityGradientSquaredAA();
        gamma_ab = prop->getDensityGradientSquaredAB();
        gamma_bb = prop->getDensityGradientSquaredBB();
    }
    if (is_meta_) {
        tau_a = prop->getKEDensityA();
        tau_b = prop->getKEDensityB();
    }

    double c = params_[0].second;
    
    //Functional
    for (int index = 0; index < ntrue; index++) {
        functional_[index] = -c*(pow(rho_a[index],4.0/3.0)+pow(rho_a[index],4.0/3.0)); 
    }
    //First Partials
    for (int index = 0; index < ntrue && deriv_ >= 1; index++) {
        v_rho_a_[index] = -4.0/3.0*c*pow(rho_a[index],1.0/3.0); 
        v_rho_a_[index] = -4.0/3.0*c*pow(rho_b[index],1.0/3.0); 
        if (is_gga_) {
            //Nothing for S
        }
    }
    //Second Partials
    for (int index = 0; index < ntrue && deriv_ >= 2; index++) {
        v_rho_a_rho_a_[index] = -4.0/9.0*c*pow(rho_a[index],-2.0/3.0); 
        v_rho_b_rho_b_[index] = -4.0/9.0*c*pow(rho_b[index],-2.0/3.0); 
        v_rho_a_rho_b_[index] = 0.0;
        if (is_gga_) {
            //Nothing for S
        }
    }
}

}}

