/*
 *  null_functional.cc
 *  Definition of class NULL_Functional for use in KS-DFT (in the case where no correlation functional is present)
 *
 *  Created by Robert Parrish on 05/02/10.
 *
 */

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "functional.h"
#include "integrator.h"
#include "null_functional.h"
#include <libmints/mints.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
NULL_Functional::NULL_Functional(int _block_size): Functional(_block_size)
{
    value_ = init_array(_block_size);
    gradA_ = init_array(_block_size);
    gradB_ = init_array(_block_size);
}
NULL_Functional::~NULL_Functional()
{
    free(value_);
    free(gradA_);
    free(gradB_);
}
void NULL_Functional::computeFunctionalRKS(shared_ptr<Properties> prop)
{
}
void NULL_Functional::computeFunctionalUKS(shared_ptr<Properties> prop)
{
}
string NULL_Functional::getName()
{
	string s("");
	s+="NULL";
	return s;
}
string NULL_Functional::getDescription()
{
	string s("  ");//Leave the spaces
	s+="Blank functional for use when no correlation is requested";
	return s;
}
string NULL_Functional::getCitation()
{
	string s("  "); //Leave the spaces
	s+="Parrish, R. P., Nature, 2010.";
	return s;
}


}}
