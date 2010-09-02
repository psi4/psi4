/**********************************************************
* S_functional.cc: definiitions for S_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include <libmints/properties.h>
#include <libciomr/libciomr.h>
#include "S_functional.h"
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

S_Functional::S_Functional(int npoints, int deriv) : Functional(npoints, deriv)
{
}
S_Functional::~S_Functional()
{
}
void S_Functional::init()
{
}
void S_Functional::computeRKSFunctional(shared_ptr<Properties> prop)
{
}
void S_Functional::computeUKSFunctional(shared_ptr<Properties> prop)
{
}

}}

