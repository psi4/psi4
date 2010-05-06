/*
 *  x_lda_functional.cc
 *  Definition of class X_B88_Functional for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "functional.h"
#include "integrator.h"
#include "x_b88_functional.h"
#include <libmints/properties.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
X_B88_Functional::X_B88_Functional(int _block_size): Functional(_block_size)
{
    value_ = init_array(_block_size);
    gradA_ = init_array(_block_size);
    gradB_ = init_array(_block_size);
    gradAA_ = init_array(_block_size);
    gradAB_ = init_array(_block_size);
    gradBB_ = init_array(_block_size);
}
X_B88_Functional::~X_B88_Functional()
{
    free(value_);
    free(gradA_);
    free(gradB_);
    free(gradAA_);
    free(gradAB_);
    free(gradBB_);
}
void X_B88_Functional::computeFunctionalRKS(shared_ptr<Properties> prop)
{
	double tol = 1.0E-20;
	const double *rho = prop->getDensity();
        const double *del_rho_2 = prop->getDensityGradientSquared();
        int ntrue = prop->getTrueSize();
        double c = 3.0/8.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(M_PI,-1.0/3.0);
        double d = 0.0042;
        double xs, rho13, rho1, rho83, sqrt_deal, drho2, arcsinhxs, den;

        for (int grid_index = 0; grid_index<ntrue; grid_index++) {

            //LDA Contribution:
            rho1 = rho[grid_index];
            rho13 = pow(rho[grid_index],1.0/3.0);
            rho83 = rho1*rho1*rho13*rho13;	    

            value_[grid_index] = -c*rho13*rho1;
	    gradA_[grid_index] = -4.0/3.0*c*rho13;

            //Gradient Correction:
            drho2 = del_rho_2[grid_index];
            xs = sqrt(drho2)/(rho13*rho1); //x_s = |\nabla \rho_s |/\rho^{4/3}
            //Note that $\arcsinh(x) = \ln(x+\sqrt{x^2+1})$
            arcsinhxs = log(xs+sqrt(xs*xs+1.0));
            den = 1.0+6.0*d*xs*arcsinhxs;
            sqrt_deal = sqrt(1.0+drho2/rho83);
            
            value_[grid_index] += -d*rho13*rho1*xs*xs/den;
            gradA_[grid_index] += 4.0/3.0*d*drho2/(rho1*rho1*rho13*den);
            gradA_[grid_index] += d*drho2/(rho1*rho13*den*den)*(-8.0*d*sqrt(drho2)*arcsinhxs/(rho1*rho1*rho13)-8.0*d*drho2/(rho1*rho83*sqrt_deal));

	    gradB_[grid_index] = gradA_[grid_index];

            gradAA_[grid_index] = -d/(rho1*rho13*den);
            gradAA_[grid_index] += d*drho2/(rho1*rho13*den*den)*(3.0*d*arcsinhxs/(sqrt(drho2)*rho13*rho1)+3.0*d/(rho83*sqrt_deal));

            gradBB_[grid_index] = gradAA_[grid_index]; //It's called symmetry.
            //gradAB_[grid_index] = 0.0; //Already there

            //fprintf(outfile,"\n  FUNCTIONAL TEST:\n");
            //fprintf(outfile,"  Value: %14.10f\n",2.0*value_[grid_index]);
            //fprintf(outfile,"  Grad_A: %14.10f\n",gradA_[grid_index]);
            //fprintf(outfile,"  Grad_B: %14.10f\n",gradB_[grid_index]);
            //fprintf(outfile,"  Grad_AA: %14.10f\n",gradAA_[grid_index]);
            //fprintf(outfile,"  Grad_AB: %14.10f\n",gradAB_[grid_index]);
            //fprintf(outfile,"  Grad_BB: %14.10f\n",gradBB_[grid_index]);
            
        }
}
void X_B88_Functional::computeFunctionalUKS(shared_ptr<Properties> prop)
{
    //TODO
}
string X_B88_Functional::getName()
{
	string s("");
	s+="B88";
	return s;
}
string X_B88_Functional::getDescription()
{
	string s("  ");//Leave the spaces
	s+="Becke Exchange Functional with Correct Asymptotic Behavior";
	return s;
}
string X_B88_Functional::getCitation()
{
	string s("  "); //Leave the spaces
	s+="A. D. Becke, Phys. Rev. A, 38(6), pp. 3098-3100, 1988";
	return s;
}


}}
