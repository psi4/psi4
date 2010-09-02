/**********************************************************
* factory.cc: defines for functional factory for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include "functional.h"
#include "S_functional.h"

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

boost::shared_ptr<Functional> Functional::createFunctional(const std::string & name, int npoints, int deriv)
{
 //   if (name == "S")
 //       return boost::shared_ptr<Functional>( new S_Functional(npoints,deriv));
}

}}

