/**********************************************************
* factory.cc: defines for functional factory for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include "functional.h"
#include "superfunctional.h"
#include "S_functional.h"
#include "B88_functional.h"
#include <boost/algorithm/string.hpp>
#include <string>
#include <sstream>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

boost::shared_ptr<Functional> Functional::createFunctional(const std::string & name, int npoints, int deriv)
{
    if (boost::to_upper_copy(name) == "S")
        return boost::shared_ptr<Functional> (new S_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == "B88")
        return boost::shared_ptr<Functional> (new B88_Functional(npoints,deriv));
}
std::string Functional::availableFunctionals()
{
    std::stringstream f;
    f << "   Available Exchange Functionals:   " << endl;
    f << "    Name:        LSDA:     GGA:     Meta: " << endl;
    f << "  ------------ --------- -------- ---------" << endl;
    
    f << "   S               X                   " << endl;
    f << "   B88             X        X          " << endl;
   
    f << endl; 
    f << "   Available Correlation Functionals:   " << endl;
    f << "    Name:        LSDA:     GGA:     Meta: " << endl;
    f << "  ------------ --------- -------- ---------" << endl;
    
        
    return f.str();
    
}
boost::shared_ptr<SuperFunctional> SuperFunctional::createSuperFunctional(const std::string & name, int npoints, int deriv)
{
    shared_ptr<SuperFunctional> superfun = (shared_ptr<SuperFunctional>) new SuperFunctional(npoints, deriv);

    //Alias Table
    if (boost::to_upper_copy(name) == "S") {
        superfun->addFunctional(Functional::createFunctional("S",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("S");
        superfun->setDescription("Simple Slater LSDA exchange functional");
        superfun->setCitation("J.C. Slater, Phys. Rev., 81(3):385–390, 1951");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion(""),  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == "B88") {
        superfun->addFunctional(Functional::createFunctional("B88",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B88");
        superfun->setDescription("Becke 88 Exchange");
        superfun->setCitation("A.D. Becke, Phys. Rev. A, 38(6):3098–3100, 1988");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion(""),  0.0000000000000000E+00);
    }

    return superfun;
}

}}

