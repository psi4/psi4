/**********************************************************
* functional.cc: definiitions for functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include <boost/algorithm/string.hpp> 
#include <libciomr/libciomr.h>
#include "functional.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {


Functional::Functional(int npoints, int deriv) : npoints_(npoints), deriv_(deriv)
{
    //Default opposite spin density cutoff
    cutoff_ = 1E-12;
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
        if (is_meta_) {
            v_tau_a_ = init_array(npoints_);
            v_tau_b_ = init_array(npoints_);
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
        if (is_meta_) {
            v_rho_a_tau_a_ = init_array(npoints_);
            v_rho_a_tau_b_ = init_array(npoints_);
            v_rho_b_tau_a_ = init_array(npoints_);
            v_rho_b_tau_b_ = init_array(npoints_);
            v_tau_a_tau_b_ = init_array(npoints_);
            v_tau_a_tau_b_ = init_array(npoints_);
            v_tau_b_tau_b_ = init_array(npoints_);
            if (is_gga_) {
                v_gamma_aa_tau_a_ = init_array(npoints_);
                v_gamma_aa_tau_b_ = init_array(npoints_);
                v_gamma_ab_tau_a_ = init_array(npoints_);
                v_gamma_ab_tau_b_ = init_array(npoints_);
                v_gamma_bb_tau_a_ = init_array(npoints_);
                v_gamma_bb_tau_b_ = init_array(npoints_);
            }
        }
    }      
}
void Functional::reallocate(int npoints, int deriv)
{
    if (npoints == npoints_ && deriv == deriv_)
        return;

    release();
    npoints_ = npoints;
    deriv_ = deriv;
    allocate();
}
void Functional::release()
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
        if (is_meta_) {
            free(v_tau_a_);
            free(v_tau_b_);
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
        if (is_meta_) {
            free(v_rho_a_tau_a_);
            free(v_rho_a_tau_b_);
            free(v_rho_b_tau_a_);
            free(v_rho_b_tau_b_);
            free(v_tau_a_tau_b_);
            free(v_tau_a_tau_b_);
            free(v_tau_b_tau_b_);
            if (is_gga_) {
                free(v_gamma_aa_tau_a_);
                free(v_gamma_aa_tau_b_);
                free(v_gamma_ab_tau_a_);
                free(v_gamma_ab_tau_b_);
                free(v_gamma_bb_tau_a_);
                free(v_gamma_bb_tau_b_);
            }
        }
    }      
}
Functional::~Functional()
{
    release();
}
std::string Functional::getParametersString() 
{
    std::stringstream params;

    params << "     Parameter           Value         \n";
    params << "  -------------------------------------\n";

    params.precision(16);

    for (int A = 0; A<params_.size(); A++) {
        params.width(2);
        params << "  ";
        params.width(15);
        params << left << params_[A].first; 
        params.width(24);
        params << scientific << left << params_[A].second << endl; 
    }
    
    return params.str();
}
void Functional::setParameter(const std::string & key, double val)
{
    for (int A = 0; A<params_.size(); A++) {
        if (boost::to_upper(params_[A].first) == boost::to_upper(key.upper()))
            params_[A] = make_pair(key,val);
    }

}

}}
