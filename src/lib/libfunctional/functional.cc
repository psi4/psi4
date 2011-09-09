/**********************************************************
* functional.cc: definitions for functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include <boost/algorithm/string.hpp> 
#include <libmints/points.h>
#include <libciomr/libciomr.h>
#include "functional.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

#if 0
std::string Functional::testFunctionals()
{
    std::vector<std::string> names = Functional::availableNames();
    std::stringstream s;

    s << "  Testing Functionals:" << endl;

    boost::shared_ptr<Properties> props = Properties::get_testbed();
    int npoints = props->npoints();

    for (int A = 0; A < names.size(); A++ ) { 
        boost::shared_ptr<Functional> func = Functional::createFunctional(names[A],npoints,2);
        s << func->testFunctional(props);
    }

    return s.str();
}
std::string Functional::testFunctional(boost::shared_ptr<Properties> props)
{    
    std::stringstream s;
    
    int npoints = props->npoints();
    double* rho_a = props->property_value("RHO_A")->pointer();
    double* rho_b = props->property_value("RHO_B")->pointer();
    double* gamma_aa = props->property_value("GAMMA_AA")->pointer();
    double* gamma_ab = props->property_value("GAMMA_AB")->pointer();
    double* gamma_bb = props->property_value("GAMMA_BB")->pointer();
    double* tau_a = props->property_value("TAU_A")->pointer();
    double* tau_b = props->property_value("TAU_B")->pointer();

    double* functional = getFunctional();

    double* v_rho_a = getV_RhoA();
    double* v_rho_b = getV_RhoB();
    double* v_gamma_aa = getV_GammaAA();
    double* v_gamma_ab = getV_GammaAB();
    double* v_gamma_bb = getV_GammaBB();
    double* v_tau_a = getV_TauA();
    double* v_tau_b = getV_TauB();

    double* v_rho_a_rho_a = getV_RhoA_RhoA(); 
    double* v_rho_a_rho_b = getV_RhoA_RhoB(); 
    double* v_rho_b_rho_b = getV_RhoB_RhoB(); 
    
    double* v_rho_a_gamma_aa = getV_RhoA_GammaAA(); 
    double* v_rho_a_gamma_ab = getV_RhoA_GammaAB(); 
    double* v_rho_a_gamma_bb = getV_RhoA_GammaBB(); 
    double* v_rho_b_gamma_aa = getV_RhoB_GammaAA(); 
    double* v_rho_b_gamma_ab = getV_RhoB_GammaAB(); 
    double* v_rho_b_gamma_bb = getV_RhoB_GammaBB(); 
    
    double* v_gamma_aa_gamma_aa = getV_GammaAA_GammaAA(); 
    double* v_gamma_aa_gamma_ab = getV_GammaAA_GammaAB(); 
    double* v_gamma_aa_gamma_bb = getV_GammaAA_GammaBB(); 
    double* v_gamma_ab_gamma_ab = getV_GammaAB_GammaAB(); 
    double* v_gamma_ab_gamma_bb = getV_GammaAB_GammaBB(); 
    double* v_gamma_bb_gamma_bb = getV_GammaBB_GammaBB(); 
    
    double* v_rho_a_tau_a = getV_RhoA_TauA(); 
    double* v_rho_a_tau_b = getV_RhoA_TauB(); 
    double* v_rho_b_tau_a = getV_RhoB_TauA(); 
    double* v_rho_b_tau_b = getV_RhoB_TauB(); 
    
    double* v_gamma_aa_tau_a = getV_GammaAA_TauA(); 
    double* v_gamma_ab_tau_a = getV_GammaAB_TauA(); 
    double* v_gamma_bb_tau_a = getV_GammaBB_TauA(); 
    double* v_gamma_aa_tau_b = getV_GammaAA_TauB(); 
    double* v_gamma_ab_tau_b = getV_GammaAB_TauB(); 
    double* v_gamma_bb_tau_b = getV_GammaBB_TauB(); 
 
    double* v_tau_a_tau_a = getV_TauA_TauA(); 
    double* v_tau_a_tau_b = getV_TauA_TauB(); 
    double* v_tau_b_tau_b = getV_TauB_TauB(); 
    
    s << endl;
    s << "  Testing Functional " << getName() << ":" << endl;

    s << endl;       
 
    s << "  RKS Results:" << endl;
    s << " ---------------" << endl;

    computeRKSFunctional(props);


    for (int Q = 0; Q < npoints; Q++) {

        //What conditions are we testing?
        s.setf(ios::scientific);
        s.precision(2);
        s << endl;
        s << " **rho_a= " << rho_a[Q] << " rho_b= " << rho_b[Q] << " gamma_aa= " << gamma_aa[Q] << " gamma_ab= " << gamma_ab[Q] << \
            " gamma_bb= " << gamma_bb[Q] << " tau_a= " << tau_a[Q] << " tau_b= " << tau_b[Q] << endl;             

        s.setf(ios::scientific);
        s.precision(11);

        s << endl;

        //Functional
        s << "   " << left << setw(20) << "functional" << "= " << right << setw(30) << functional[Q] << endl;            
        
        s << endl;
        
        //First Partials
        s << "   " << left << setw(20) << "v_rho_a" << "= " << right << setw(30) << v_rho_a[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_b" << "= " << right << setw(30) << v_rho_b[Q] << endl;            
 
        if (isGGA()) {
            s << "   " << left << setw(20) << "v_gamma_aa" << "= " << right << setw(30) << v_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab" << "= " << right << setw(30) << v_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb" << "= " << right << setw(30) << v_gamma_bb[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
        }
        if (isMeta()) {
            s << "   " << left << setw(20) << "v_tau_a" << "= " << right << setw(30) << v_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_b" << "= " << right << setw(30) << v_tau_b[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
        }

        //Second Partials
        s << endl;

        s << "   " << left << setw(20) << "v_rho_a_rho_a" << "= " << right << setw(30) << v_rho_a_rho_a[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_a_rho_b" << "= " << right << setw(30) << v_rho_a_rho_b[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_b_rho_b" << "= " << right << setw(30) << v_rho_b_rho_b[Q] << endl;            
    
        if (isGGA()) {
            s << "   " << left << setw(20) << "v_rho_a_gamma_aa" << "= " << right << setw(30) << v_rho_a_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_ab" << "= " << right << setw(30) << v_rho_a_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_bb" << "= " << right << setw(30) << v_rho_a_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_aa" << "= " << right << setw(30) << v_rho_b_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_ab" << "= " << right << setw(30) << v_rho_b_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_bb" << "= " << right << setw(30) << v_rho_b_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_aa" << "= " << right << setw(30) << v_gamma_aa_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_ab" << "= " << right << setw(30) << v_gamma_aa_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_bb" << "= " << right << setw(30) << v_gamma_aa_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_ab" << "= " << right << setw(30) << v_gamma_ab_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_bb" << "= " << right << setw(30) << v_gamma_ab_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_gamma_bb" << "= " << right << setw(30) << v_gamma_bb_gamma_bb[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_rho_a_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
        }
        if (isMeta()) {
            s << "   " << left << setw(20) << "v_rho_a_tau_a" << "= " << right << setw(30) << v_rho_a_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_tau_b" << "= " << right << setw(30) << v_rho_a_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_a" << "= " << right << setw(30) << v_rho_b_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_b" << "= " << right << setw(30) << v_rho_b_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_a" << "= " << right << setw(30) << v_tau_a_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_b" << "= " << right << setw(30) << v_tau_a_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_b_tau_b" << "= " << right << setw(30) << v_tau_b_tau_b[Q] << endl;            
            if (isGGA()) {
                s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << v_gamma_aa_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << v_gamma_aa_tau_b[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << v_gamma_ab_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << v_gamma_ab_tau_b[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << v_gamma_bb_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << v_gamma_bb_tau_b[Q] << endl;            
            } else {
                s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            }
        } else {
            s << "   " << left << setw(20) << "v_rho_a_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_b_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
        }
    }
    s << endl;       
 
    s << "  UKS Results:" << endl;
    s << " ---------------" << endl;

    computeUKSFunctional(props);

    for (int Q = 0; Q < npoints; Q++) {

        //What conditions are we testing?
        s.setf(ios::scientific);
        s.precision(2);
        s << endl;
        s << " **rho_a= " << rho_a[Q] << " rho_b= " << rho_b[Q] << " gamma_aa= " << gamma_aa[Q] << " gamma_ab= " << gamma_ab[Q] << \
            " gamma_bb= " << gamma_bb[Q] << " tau_a= " << tau_a[Q] << " tau_b= " << tau_b[Q] << endl;             

        s.setf(ios::scientific);
        s.precision(11);

        s << endl;

        //Functional
        s << "   " << left << setw(20) << "functional" << "= " << right << setw(30) << functional[Q] << endl;            
        
        s << endl;
        
        //First Partials
        s << "   " << left << setw(20) << "v_rho_a" << "= " << right << setw(30) << v_rho_a[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_b" << "= " << right << setw(30) << v_rho_b[Q] << endl;            
 
        if (isGGA()) {
            s << "   " << left << setw(20) << "v_gamma_aa" << "= " << right << setw(30) << v_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab" << "= " << right << setw(30) << v_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb" << "= " << right << setw(30) << v_gamma_bb[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
        }
        if (isMeta()) {
            s << "   " << left << setw(20) << "v_tau_a" << "= " << right << setw(30) << v_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_b" << "= " << right << setw(30) << v_tau_b[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
        }

        //Second Partials
        s << endl;

        s << "   " << left << setw(20) << "v_rho_a_rho_a" << "= " << right << setw(30) << v_rho_a_rho_a[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_a_rho_b" << "= " << right << setw(30) << v_rho_a_rho_b[Q] << endl;            
        s << "   " << left << setw(20) << "v_rho_b_rho_b" << "= " << right << setw(30) << v_rho_b_rho_b[Q] << endl;            
    
        if (isGGA()) {
            s << "   " << left << setw(20) << "v_rho_a_gamma_aa" << "= " << right << setw(30) << v_rho_a_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_ab" << "= " << right << setw(30) << v_rho_a_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_bb" << "= " << right << setw(30) << v_rho_a_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_aa" << "= " << right << setw(30) << v_rho_b_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_ab" << "= " << right << setw(30) << v_rho_b_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_bb" << "= " << right << setw(30) << v_rho_b_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_aa" << "= " << right << setw(30) << v_gamma_aa_gamma_aa[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_ab" << "= " << right << setw(30) << v_gamma_aa_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_bb" << "= " << right << setw(30) << v_gamma_aa_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_ab" << "= " << right << setw(30) << v_gamma_ab_gamma_ab[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_bb" << "= " << right << setw(30) << v_gamma_ab_gamma_bb[Q] << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_gamma_bb" << "= " << right << setw(30) << v_gamma_bb_gamma_bb[Q] << endl;            
        } else {
            s << "   " << left << setw(20) << "v_rho_a_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_aa" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_ab" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_gamma_bb" << "= " << right << setw(30) << 0.0 << endl;            
        }
        if (isMeta()) {
            s << "   " << left << setw(20) << "v_rho_a_tau_a" << "= " << right << setw(30) << v_rho_a_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_a_tau_b" << "= " << right << setw(30) << v_rho_a_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_a" << "= " << right << setw(30) << v_rho_b_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_b" << "= " << right << setw(30) << v_rho_b_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_a" << "= " << right << setw(30) << v_tau_a_tau_a[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_b" << "= " << right << setw(30) << v_tau_a_tau_b[Q] << endl;            
            s << "   " << left << setw(20) << "v_tau_b_tau_b" << "= " << right << setw(30) << v_tau_b_tau_b[Q] << endl;            
            if (isGGA()) {
                s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << v_gamma_aa_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << v_gamma_aa_tau_b[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << v_gamma_ab_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << v_gamma_ab_tau_b[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << v_gamma_bb_tau_a[Q] << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << v_gamma_bb_tau_b[Q] << endl;            
            } else {
                s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
                s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            }
        } else {
            s << "   " << left << setw(20) << "v_rho_a_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_a_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_rho_b_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_a_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_tau_b_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_aa_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_ab_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_tau_a" << "= " << right << setw(30) << 0.0 << endl;            
            s << "   " << left << setw(20) << "v_gamma_bb_tau_b" << "= " << right << setw(30) << 0.0 << endl;            
        }
    }
    return s.str();
}
#endif
Functional::Functional(int npoints, int deriv) : npoints_(npoints), deriv_(deriv)
{
    //Default opposite spin density cutoff
    cutoff_ = 1E-20;
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
            v_tau_a_tau_a_ = init_array(npoints_);
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
            free(v_tau_a_tau_a_);
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
        if (boost::to_upper_copy(params_[A].first) == boost::to_upper_copy(key))
            params_[A] = make_pair(key,val);
    }

}

}}
