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
#include "x_b88_functional.h"
#include "x_b3_functional.h"
#include "null_functional.h"
#include <libmints/mints.h>
#include <string>
#include <cstring>
using namespace psi;
using namespace std;

namespace psi { namespace scf {

const double Functional::tol = 1.0e-20;

Functional * FunctionalFactory::getExchangeFunctional(string x_id,string c_id, int block_size)
{
    if (x_id == "X_LDA") {
        return new X_LDA_Functional(block_size);
    }
    if (x_id == "X_B88") {
        return new X_B88_Functional(block_size);
    }
    if (x_id == "X_B3") {
        return new X_B3_Functional(block_size);
    }
    else {
        throw std::domain_error("Requested Exchange Functional does not exist!");		
    }
}
Functional * FunctionalFactory::getCorrelationFunctional(string x_id,string c_id, int block_size)
{
    if (c_id == "") {
        return new NULL_Functional(block_size);
    }
    else {
        throw std::domain_error("Requested Correlation Functional does not exist!");		
    }
}

}}
