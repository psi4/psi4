/*
 *  functionalfactory.cc
 *  Definition of class FunctionalFactory for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "functional.h"
#include "functionalfactory.h"
#include "x_lda_functional.h"
#include <string>
#include <cstring>
using namespace psi;

namespace psi { namespace scf {

Functional * FunctionalFactory::getFunctional(string id, int block_size)
{
	if (id == "X_LDA") {
		//X_LDA_Functional f(1);
		//printf("%s",f.getName().c_str());
		return new X_LDA_Functional(block_size);
	}
	else {
            throw std::domain_error("Requested Functional does not exist!");		
	}
}

}}
