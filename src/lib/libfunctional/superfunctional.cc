/**********************************************************
* superfunctional.cc: definitions for superfunctionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include "superfunctional.h"
#include <boost/algorithm/string.hpp>
#include <libmints/mints.h>
#include <libfock/points.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

#if 0
std::string SuperFunctional::testSuperFunctionals()
{
    std::vector<std::string> names = SuperFunctional::availableNames();
    std::stringstream s;

    s << "  Testing SuperFunctionals:" << endl;

    boost::shared_ptr<Properties> props = Properties::get_testbed();
    int npoints = props->npoints();

    for (int A = 0; A < names.size(); A++ ) {
        boost::shared_ptr<SuperFunctional> func = SuperFunctional::createSuperFunctional(names[A],npoints,2);
        s << func->testSuperFunctional(props);
    }

    return s.str();
}
std::string SuperFunctional::testSuperFunctional(boost::shared_ptr<Properties> props)
{
    std::stringstream s;
    int npoints = props->npoints();
    const double* rho_a = props->property_value("RHO_A")->pointer();
    const double* rho_b = props->property_value("RHO_B")->pointer();
    const double* gamma_aa = props->property_value("GAMMA_AA")->pointer();
    const double* gamma_ab = props->property_value("GAMMA_AB")->pointer();
    const double* gamma_bb = props->property_value("GAMMA_BB")->pointer();
    const double* tau_a = props->property_value("TAU_A")->pointer();
    const double* tau_b = props->property_value("TAU_B")->pointer();

    double* functional = getFunctionalValue();

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
    s << "  Testing SuperFunctional " << getName() << ":" << endl;

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
        s.precision(12);

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
        s.precision(12);

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
SuperFunctional::SuperFunctional(int npoints, int deriv) : npoints_(npoints), deriv_(deriv),
    exact_exchange_(0.0), pt2_(0.0), dashD_weight_(0.0), omega_(0.0), oldGGA_(false), oldMeta_(false)
{
    allocate();
}
SuperFunctional::~SuperFunctional()
{
    release();
}
boost::shared_ptr<SuperFunctional> SuperFunctional::buildSuperFunctional(const std::string & build, int npoints, int deriv)
{
    boost::shared_ptr<SuperFunctional> superfun = (boost::shared_ptr<SuperFunctional>) new SuperFunctional(npoints,deriv);
    //TODO
    return superfun;
}
void SuperFunctional::allocate()
{
    functional_ = init_array(npoints_);

    if (deriv_ >= 1) {
        v_rho_a_ = init_array(npoints_);
        v_rho_b_ = init_array(npoints_);
        if (isGGA()) {
            v_gamma_aa_ = init_array(npoints_);
            v_gamma_ab_ = init_array(npoints_);
            v_gamma_bb_ = init_array(npoints_);
        }
        if (isMeta()) {
            v_tau_a_ = init_array(npoints_);
            v_tau_b_ = init_array(npoints_);
        }
    }
    if (deriv_ >= 2) {
        v_rho_a_rho_a_ = init_array(npoints_);
        v_rho_a_rho_b_ = init_array(npoints_);
        v_rho_b_rho_b_ = init_array(npoints_);
        if (isGGA()) {
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
        if (isMeta()) {
            v_rho_a_tau_a_ = init_array(npoints_);
            v_rho_a_tau_b_ = init_array(npoints_);
            v_rho_b_tau_a_ = init_array(npoints_);
            v_rho_b_tau_b_ = init_array(npoints_);
            v_tau_a_tau_b_ = init_array(npoints_);
            v_tau_a_tau_b_ = init_array(npoints_);
            v_tau_b_tau_b_ = init_array(npoints_);
            if (isGGA()) {
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
void SuperFunctional::reallocate(int npoints, int deriv)
{
    if (npoints == npoints_ && deriv == deriv_ && isGGA() == oldGGA_ && isMeta() == oldMeta_)
        return;

    //Clear the total registers
    release();
    npoints_ = npoints;
    deriv_ = deriv;
    oldGGA_ = isGGA();
    oldMeta_ = isMeta();
    allocate();

    //Clear the individual registers
    for (int A = 0; A < functionals_.size(); A++) {
        functionals_[A].first->setNPoints(npoints_);
        functionals_[A].first->setDeriv(deriv_);
    }
}
void SuperFunctional::release()
{
    free(functional_);

    if (deriv_ >= 1) {
        free(v_rho_a_);
        free(v_rho_b_);
        if (oldGGA_) {
            free(v_gamma_aa_);
            free(v_gamma_ab_);
            free(v_gamma_bb_);
        }
        if (oldMeta_) {
            free(v_tau_a_);
            free(v_tau_b_);
        }
    }
    if (deriv_ >= 2) {
        free(v_rho_a_rho_a_);
        free(v_rho_a_rho_b_);
        free(v_rho_b_rho_b_);
        if (oldGGA_) {
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
        if (oldMeta_) {
            free(v_rho_a_tau_a_);
            free(v_rho_a_tau_b_);
            free(v_rho_b_tau_a_);
            free(v_rho_b_tau_b_);
            free(v_tau_a_tau_b_);
            free(v_tau_a_tau_b_);
            free(v_tau_b_tau_b_);
            if (oldGGA_) {
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
void SuperFunctional::setParameter(const std::string & name, const std::string & param, double val)
{
    for (int A = 0; A < functionals_.size(); A++)
        if (to_upper_copy(functionals_[A].first->getName()) == to_upper_copy(name)) {
            functionals_[A].first->setParameter(param,val);
            break;
        }
}
void SuperFunctional::setOmega(double omega)
{
    omega_ = omega;
    for (int A = 0; A < functionals_.size(); A++) {
        functionals_[A].first->setParameter("omega", omega_);
    }
}
int SuperFunctional::getLocalAnsatz()
{
    if (isMeta()) return 2;
    else if (isGGA()) return 1;
    else return 0;
}
bool SuperFunctional::isGGA()
{
    for (int A = 0; A < functionals_.size(); A++)
        if (functionals_[A].first->isGGA())
            return true;

    return false;
}
bool SuperFunctional::isMeta()
{
    for (int A = 0; A < functionals_.size(); A++)
        if (functionals_[A].first->isMeta())
            return true;

    return false;
}
void SuperFunctional::computeRKSFunctional(boost::shared_ptr<RKSFunctions> props)
{
    for (int A = 0; A < functionals_.size(); A++) {
        functionals_[A].first->computeRKSFunctional(props);
    }
    //Add up the results
    collectResults();
}
void SuperFunctional::computeUKSFunctional(boost::shared_ptr<UKSFunctions> props)
{
    for (int A = 0; A < functionals_.size(); A++) {
        functionals_[A].first->computeUKSFunctional(props);
    }
    //Add up the results
    collectResults();
}
void SuperFunctional::collectResults()
{
    //Zero stuff out
    memset(functional_,'\0',npoints_*sizeof(double));

    if (deriv_ >= 1) {
        memset(v_rho_a_,'\0',npoints_*sizeof(double));
        memset(v_rho_b_,'\0',npoints_*sizeof(double));
        if (isGGA()) {
            memset(v_gamma_aa_,'\0',npoints_*sizeof(double));
            memset(v_gamma_ab_,'\0',npoints_*sizeof(double));
            memset(v_gamma_bb_,'\0',npoints_*sizeof(double));
        }
        if (isMeta()) {
            memset(v_tau_a_,'\0',npoints_*sizeof(double));
            memset(v_tau_b_,'\0',npoints_*sizeof(double));
        }
    }
    if (deriv_ >= 2) {
        memset(v_rho_a_rho_a_,'\0',npoints_*sizeof(double));
        memset(v_rho_a_rho_b_,'\0',npoints_*sizeof(double));
        memset(v_rho_b_rho_b_,'\0',npoints_*sizeof(double));
        if (isGGA()) {
            memset(v_rho_a_gamma_aa_,'\0',npoints_*sizeof(double));
            memset(v_rho_a_gamma_ab_,'\0',npoints_*sizeof(double));
            memset(v_rho_a_gamma_bb_,'\0',npoints_*sizeof(double));
            memset(v_rho_b_gamma_aa_,'\0',npoints_*sizeof(double));
            memset(v_rho_b_gamma_ab_,'\0',npoints_*sizeof(double));
            memset(v_rho_b_gamma_bb_,'\0',npoints_*sizeof(double));
            memset(v_gamma_aa_gamma_aa_,'\0',npoints_*sizeof(double));
            memset(v_gamma_aa_gamma_ab_,'\0',npoints_*sizeof(double));
            memset(v_gamma_aa_gamma_bb_,'\0',npoints_*sizeof(double));
            memset(v_gamma_ab_gamma_ab_,'\0',npoints_*sizeof(double));
            memset(v_gamma_ab_gamma_bb_,'\0',npoints_*sizeof(double));
            memset(v_gamma_bb_gamma_bb_,'\0',npoints_*sizeof(double));
        }
        if (isMeta()) {
            memset(v_rho_a_tau_a_,'\0',npoints_*sizeof(double));
            memset(v_rho_a_tau_b_,'\0',npoints_*sizeof(double));
            memset(v_rho_b_tau_a_,'\0',npoints_*sizeof(double));
            memset(v_rho_b_tau_b_,'\0',npoints_*sizeof(double));
            memset(v_tau_a_tau_b_,'\0',npoints_*sizeof(double));
            memset(v_tau_a_tau_b_,'\0',npoints_*sizeof(double));
            memset(v_tau_b_tau_b_,'\0',npoints_*sizeof(double));
            if (isGGA()) {
                memset(v_gamma_aa_tau_a_,'\0',npoints_*sizeof(double));
                memset(v_gamma_aa_tau_b_,'\0',npoints_*sizeof(double));
                memset(v_gamma_ab_tau_a_,'\0',npoints_*sizeof(double));
                memset(v_gamma_ab_tau_b_,'\0',npoints_*sizeof(double));
                memset(v_gamma_bb_tau_a_,'\0',npoints_*sizeof(double));
                memset(v_gamma_bb_tau_b_,'\0',npoints_*sizeof(double));
            }
        }
    }
    for (int A = 0; A < functionals_.size(); A++) {

        double weight = functionals_[A].second;

        //Get register handles
        double* functional = functionals_[A].first->getFunctional();

        double* v_rho_a = functionals_[A].first->getV_RhoA();
        double* v_rho_b = functionals_[A].first->getV_RhoB();
        double* v_gamma_aa = functionals_[A].first->getV_GammaAA();
        double* v_gamma_ab = functionals_[A].first->getV_GammaAB();
        double* v_gamma_bb = functionals_[A].first->getV_GammaBB();
        double* v_tau_a = functionals_[A].first->getV_TauA();
        double* v_tau_b = functionals_[A].first->getV_TauB();

        double* v_rho_a_rho_a = functionals_[A].first->getV_RhoA_RhoA();
        double* v_rho_a_rho_b = functionals_[A].first->getV_RhoA_RhoB();
        double* v_rho_b_rho_b = functionals_[A].first->getV_RhoB_RhoB();

        double* v_rho_a_gamma_aa = functionals_[A].first->getV_RhoA_GammaAA();
        double* v_rho_a_gamma_ab = functionals_[A].first->getV_RhoA_GammaAB();
        double* v_rho_a_gamma_bb = functionals_[A].first->getV_RhoA_GammaBB();
        double* v_rho_b_gamma_aa = functionals_[A].first->getV_RhoB_GammaAA();
        double* v_rho_b_gamma_ab = functionals_[A].first->getV_RhoB_GammaAB();
        double* v_rho_b_gamma_bb = functionals_[A].first->getV_RhoB_GammaBB();

        double* v_gamma_aa_gamma_aa = functionals_[A].first->getV_GammaAA_GammaAA();
        double* v_gamma_aa_gamma_ab = functionals_[A].first->getV_GammaAA_GammaAB();
        double* v_gamma_aa_gamma_bb = functionals_[A].first->getV_GammaAA_GammaBB();
        double* v_gamma_ab_gamma_ab = functionals_[A].first->getV_GammaAB_GammaAB();
        double* v_gamma_ab_gamma_bb = functionals_[A].first->getV_GammaAB_GammaBB();
        double* v_gamma_bb_gamma_bb = functionals_[A].first->getV_GammaBB_GammaBB();

        double* v_rho_a_tau_a = functionals_[A].first->getV_RhoA_TauA();
        double* v_rho_a_tau_b = functionals_[A].first->getV_RhoA_TauB();
        double* v_rho_b_tau_a = functionals_[A].first->getV_RhoB_TauA();
        double* v_rho_b_tau_b = functionals_[A].first->getV_RhoB_TauB();

        double* v_gamma_aa_tau_a = functionals_[A].first->getV_GammaAA_TauA();
        double* v_gamma_ab_tau_a = functionals_[A].first->getV_GammaAB_TauA();
        double* v_gamma_bb_tau_a = functionals_[A].first->getV_GammaBB_TauA();
        double* v_gamma_aa_tau_b = functionals_[A].first->getV_GammaAA_TauB();
        double* v_gamma_ab_tau_b = functionals_[A].first->getV_GammaAB_TauB();
        double* v_gamma_bb_tau_b = functionals_[A].first->getV_GammaBB_TauB();

        double* v_tau_a_tau_a = functionals_[A].first->getV_TauA_TauA();
        double* v_tau_a_tau_b = functionals_[A].first->getV_TauA_TauB();
        double* v_tau_b_tau_b = functionals_[A].first->getV_TauB_TauB();

        //Add up the results
        for (int Q = 0; Q < npoints_; Q++) {
            functional_[Q] += weight*functional[Q];

            if (deriv_ >= 1) {
                v_rho_a_[Q] += weight*v_rho_a[Q];
                v_rho_b_[Q] += weight*v_rho_b[Q];
                if (functionals_[A].first->isGGA()) {
                    v_gamma_aa_[Q] += weight*v_gamma_aa[Q];
                    v_gamma_ab_[Q] += weight*v_gamma_ab[Q];
                    v_gamma_bb_[Q] += weight*v_gamma_bb[Q];
                }
                if (functionals_[A].first->isMeta()) {
                    v_tau_a_[Q] += weight*v_tau_a[Q];
                    v_tau_b_[Q] += weight*v_tau_b[Q];
                }
            }
            if (deriv_ >= 2) {
                v_rho_a_rho_a_[Q] += weight*v_rho_a_rho_a[Q];
                v_rho_a_rho_b_[Q] += weight*v_rho_a_rho_b[Q];
                v_rho_b_rho_b_[Q] += weight*v_rho_b_rho_b[Q];
                if (functionals_[A].first->isGGA()) {
                    v_rho_a_gamma_aa_[Q] += weight*v_rho_a_gamma_aa[Q];
                    v_rho_a_gamma_ab_[Q] += weight*v_rho_a_gamma_ab[Q];
                    v_rho_a_gamma_bb_[Q] += weight*v_rho_a_gamma_bb[Q];
                    v_rho_b_gamma_aa_[Q] += weight*v_rho_b_gamma_aa[Q];
                    v_rho_b_gamma_ab_[Q] += weight*v_rho_b_gamma_ab[Q];
                    v_rho_b_gamma_bb_[Q] += weight*v_rho_b_gamma_bb[Q];
                    v_gamma_aa_gamma_aa_[Q] += weight*v_gamma_aa_gamma_aa[Q];
                    v_gamma_aa_gamma_ab_[Q] += weight*v_gamma_aa_gamma_ab[Q];
                    v_gamma_aa_gamma_bb_[Q] += weight*v_gamma_aa_gamma_bb[Q];
                    v_gamma_ab_gamma_ab_[Q] += weight*v_gamma_ab_gamma_ab[Q];
                    v_gamma_ab_gamma_bb_[Q] += weight*v_gamma_ab_gamma_bb[Q];
                    v_gamma_bb_gamma_bb_[Q] += weight*v_gamma_bb_gamma_bb[Q];
                }
                if (functionals_[A].first->isMeta()) {
                    v_rho_a_tau_a_[Q] += weight*v_rho_a_tau_a[Q];
                    v_rho_a_tau_b_[Q] += weight*v_rho_a_tau_b[Q];
                    v_rho_b_tau_a_[Q] += weight*v_rho_b_tau_a[Q];
                    v_rho_b_tau_b_[Q] += weight*v_rho_b_tau_b[Q];
                    v_tau_a_tau_b_[Q] += weight*v_tau_a_tau_b[Q];
                    v_tau_a_tau_b_[Q] += weight*v_tau_a_tau_b[Q];
                    v_tau_b_tau_b_[Q] += weight*v_tau_b_tau_b[Q];
                    if (functionals_[A].first->isGGA()) {
                        v_gamma_aa_tau_a_[Q] += weight*v_gamma_aa_tau_a[Q];
                        v_gamma_aa_tau_b_[Q] += weight*v_gamma_aa_tau_b[Q];
                        v_gamma_ab_tau_a_[Q] += weight*v_gamma_ab_tau_a[Q];
                        v_gamma_ab_tau_b_[Q] += weight*v_gamma_ab_tau_b[Q];
                        v_gamma_bb_tau_a_[Q] += weight*v_gamma_bb_tau_a[Q];
                        v_gamma_bb_tau_b_[Q] += weight*v_gamma_bb_tau_b[Q];
                    }
                }
            }
        }
    }
}
void SuperFunctional::addFunctional(const boost::shared_ptr<Functional> & f, double weight)
{
    functionals_.push_back(make_pair(f, weight));
    reallocate(npoints_,deriv_);
}
void SuperFunctional::addFunctional(const std::string & name, double weight)
{
    functionals_.push_back(make_pair(Functional::createFunctional(name, npoints_, deriv_), weight));
}
void SuperFunctional::setDashD(boost::shared_ptr<Dispersion> disp, double weight)
{
    dashD_ = disp;
    dashD_weight_ = weight;
}
std::string SuperFunctional::getComposition()
{
    std::stringstream s;
    s << "  SuperFunctional " << name_ << endl;
    s << "  Composed of:" << endl;
    for (int A = 0; A < functionals_.size(); A++) {
        s.precision(4);
        s << "  " << functionals_[A].second*100.0 << "% " << functionals_[A].first->getName() << endl;
    }
    if (isHybrid() || isDoubleHybrid() || isDashD() || isRangeCorrected())
        s << "  Additionally composed of: " << endl;
    if (isHybrid())
        s << "  " << exact_exchange_*100.0 << "% Exact Exchange" << endl;
    if (isDoubleHybrid())
        s << "  " << pt2_*100.0 << "% PT2 Long Range Correlation" << endl;
    if (isDashD())
        s << "  " << dashD_weight_*100.0 << "% " << dashD_->getName() << endl;
    if (isRangeCorrected())
        s << "  Range-Correction with omega of " << omega_ << endl;

    return s.str();
}

}}
