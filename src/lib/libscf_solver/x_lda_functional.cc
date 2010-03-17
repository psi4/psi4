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

double X_LDA_Functional::getValue(shared_ptr<Properties> prop)
{
	double tol = 1.0E-20;
	double rho = 1.0*prop->getDensity();
	if (rho<tol)
		return 0.0;
	double c = 3.0/8.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(PI,-1.0/3.0);  
	double val = -c*pow(rho,4.0/3.0);
	return val;
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
