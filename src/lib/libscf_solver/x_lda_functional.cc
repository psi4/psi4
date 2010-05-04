/*
 *  x_lda_functional.cc
 *  Definition of class X_LDA_Functional for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "functional.h"
#include "integrator.h"
#include "x_lda_functional.h"
#include <libmints/properties.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
X_LDA_Functional::X_LDA_Functional(int _block_size): Functional(_block_size)
{
    value_ = init_array(_block_size);
    gradA_ = init_array(_block_size);
    gradB_ = init_array(_block_size);
}
X_LDA_Functional::~X_LDA_Functional()
{
    free(value_);
    free(gradA_);
    free(gradB_);
}
void X_LDA_Functional::computeFunctionalRKS(shared_ptr<Properties> prop)
{
	double tol = 1.0E-20;
	const double *rho = prop->getDensity();
        int ntrue = prop->getTrueSize();
        double c = 3.0/8.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(M_PI,-1.0/3.0);

        for (int grid_index = 0; grid_index<ntrue; grid_index++) {
	    value_[grid_index] = -c*pow(rho[grid_index],4.0/3.0);
	    gradA_[grid_index] = -4.0/3.0*c*pow(rho[grid_index],1.0/3.0);
	    gradB_[grid_index] = gradA_[grid_index];
        }
}
void X_LDA_Functional::computeFunctionalUKS(shared_ptr<Properties> prop)
{
    //TODO
}
string X_LDA_Functional::getName()
{
	string s("");
	s+="X_LDA";
	return s;
}
string X_LDA_Functional::getDescription()
{
	string s("  ");//Leave the spaces
	s+="Exchange Local-Density Approximation";
	return s;
}
string X_LDA_Functional::getCitation()
{
	string s("  "); //Leave the spaces
	s+="Slater, J.C., Phys. Rev., 81(3), 1951, pp. 385-390";
	return s;
}


}}
