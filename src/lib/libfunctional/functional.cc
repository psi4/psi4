/**********************************************************
* functional.cc: definiitions for functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include <libciomr/libciomr.h>
#include "functional.h"
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {


Functional::Functional(int npoints, int deriv) : npoints_(npoints), deriv_(deriv)
{
    init();
    allocate();
}
void Functional::allocate()
{
    functional_ = init_array(npoints_);

    if (deriv_ >= 1) {
        v_rho_a_ = init_array(npoints_);
        v_rho_b_ = init_array(npoints_);
        if (is_gga_) {
            v_gamma_aa_ = init_array(npoints_);
            v_gamma_ab_ = init_array(npoints_);
            v_gamma_bb_ = init_array(npoints_);
        }
    }      
    if (deriv_ >= 2) {
        v_rho_a_rho_a_ = init_array(npoints_);
        v_rho_a_rho_b_ = init_array(npoints_);
        v_rho_b_rho_b_ = init_array(npoints_);
        if (is_gga_) {
            v_rho_a_gamma_aa_ = init_array(npoints_);
            v_rho_a_gamma_ab_ = init_array(npoints_);
            v_rho_a_gamma_bb_ = init_array(npoints_);
            v_rho_b_gamma_aa_ = init_array(npoints_);
            v_rho_b_gamma_ab_ = init_array(npoints_);
            v_rho_b_gamma_bb_ = init_array(npoints_);
            v_gamma_aa_gamma_aa_ = init_array(npoints_);
            v_gamma_aa_gamma_ab_ = init_array(npoints_);
            v_gamma_aa_gamma_bb_ = init_array(npoints_);
            v_gamma_ab_gamma_ab_ = init_array(npoints_);
            v_gamma_ab_gamma_bb_ = init_array(npoints_);
            v_gamma_bb_gamma_bb_ = init_array(npoints_);
        }
    }      
} 
Functional::~Functional()
{
    free(functional_);

    if (deriv_ >= 1) {
        free(v_rho_a_);
        free(v_rho_b_);
        if (is_gga_) {
            free(v_gamma_aa_);
            free(v_gamma_ab_);
            free(v_gamma_bb_);
        }
    }      
    if (deriv_ >= 2) {
        free(v_rho_a_rho_a_);
        free(v_rho_a_rho_b_);
        free(v_rho_b_rho_b_);
        if (is_gga_) {
            free(v_rho_a_gamma_aa_);
            free(v_rho_a_gamma_ab_);
            free(v_rho_a_gamma_bb_);
            free(v_rho_b_gamma_aa_);
            free(v_rho_b_gamma_ab_);
            free(v_rho_b_gamma_bb_);
            free(v_gamma_aa_gamma_aa_);
            free(v_gamma_aa_gamma_ab_);
            free(v_gamma_aa_gamma_bb_);
            free(v_gamma_ab_gamma_ab_);
            free(v_gamma_ab_gamma_bb_);
            free(v_gamma_bb_gamma_bb_);
        }
    }      
}
std::string Functional::getParametersString() 
{
    std::string params = "";

    params += "     Parameter      Value      \n";
    params += "-------------------------------\n";

    for (int A = 0; A<params_.size(); A++)
        params += sprintf("  %20s   %20.16E\n",params_[A].first.c_str(),params_[A].second); 

    return params;
}

}}
