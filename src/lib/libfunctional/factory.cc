/**********************************************************
* factory.cc: defines for functional factory for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include "functional.h"
#include "superfunctional.h"
#include "S_X_functional.h"
#include "B_X_functional.h"
#include "B88_X_functional.h"
#include "PBE_X_functional.h"
#include "PW91_X_functional.h"
#include "FT97B_X_functional.h"
#include "LYP_C_functional.h"
#include "VWN5RPA_C_functional.h"
#include "VWN5_C_functional.h"
#include "PZ81_C_functional.h"
#include "P86_C_functional.h"
#include "PW91_C_functional.h"
#include "PW92_C_functional.h"
#include "PBE_C_functional.h"
#include "FT97_C_functional.h"
#include "EDF1_functional.h"
#include "EDF1_X_functional.h"
#include "EDF1_C_functional.h"
#include "B97_0_functional.h"
#include "B97_1_functional.h"
#include "B97_2_functional.h"
#include "B97_D2_functional.h"
#include "HCTH_functional.h"
#include "HCTH407_functional.h"
#include "HCTH147_functional.h"
#include "HCTH120_functional.h"
#include "M05_functional.h"
#include "M05_2X_functional.h"
#include "TauHCTH_functional.h"
#include "TauHCTH0_functional.h"
#include "wS_X_functional.h"
#include "wBf_X_functional.h"
#include "wPBEf_X_functional.h"
#include "wB97_functional.h"
#include "wB97X_functional.h"
#include <boost/algorithm/string.hpp>
#include <string>
#include <sstream>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {

boost::shared_ptr<Functional> Functional::createFunctional(const std::string & name, int npoints, int deriv)
{
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("S_X")))
        return boost::shared_ptr<Functional> (new S_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B_X")))
        return boost::shared_ptr<Functional> (new B_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B88_X")))
        return boost::shared_ptr<Functional> (new B88_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_X")))
        return boost::shared_ptr<Functional> (new PBE_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW91_X")))
        return boost::shared_ptr<Functional> (new PW91_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("FT97B_X")))
        return boost::shared_ptr<Functional> (new FT97B_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("LYP_C")))
        return boost::shared_ptr<Functional> (new LYP_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("VWN5RPA_C")))
        return boost::shared_ptr<Functional> (new VWN5RPA_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("VWN5_C")))
        return boost::shared_ptr<Functional> (new VWN5_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PZ81_C")))
        return boost::shared_ptr<Functional> (new PZ81_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("P86_C")))
        return boost::shared_ptr<Functional> (new P86_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW91_C")))
        return boost::shared_ptr<Functional> (new PW91_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW92_C")))
        return boost::shared_ptr<Functional> (new PW92_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_C")))
        return boost::shared_ptr<Functional> (new PBE_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("FT97_C")))
        return boost::shared_ptr<Functional> (new FT97_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("EDF1")))
        return boost::shared_ptr<Functional> (new EDF1_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("EDF1_X")))
        return boost::shared_ptr<Functional> (new EDF1_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("EDF1_C")))
        return boost::shared_ptr<Functional> (new EDF1_C_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_0")))
        return boost::shared_ptr<Functional> (new B97_0_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_1")))
        return boost::shared_ptr<Functional> (new B97_1_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_2")))
        return boost::shared_ptr<Functional> (new B97_2_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_D2")))
        return boost::shared_ptr<Functional> (new B97_D2_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH")))
        return boost::shared_ptr<Functional> (new HCTH_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH407")))
        return boost::shared_ptr<Functional> (new HCTH407_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH147")))
        return boost::shared_ptr<Functional> (new HCTH147_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH120")))
        return boost::shared_ptr<Functional> (new HCTH120_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("M05")))
        return boost::shared_ptr<Functional> (new M05_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("M05_2X")))
        return boost::shared_ptr<Functional> (new M05_2X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("TauHCTH")))
        return boost::shared_ptr<Functional> (new TauHCTH_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("TauHCTH0")))
        return boost::shared_ptr<Functional> (new TauHCTH0_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wS_X")))
        return boost::shared_ptr<Functional> (new wS_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wBf_X")))
        return boost::shared_ptr<Functional> (new wBf_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wPBEf_X")))
        return boost::shared_ptr<Functional> (new wPBEf_X_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wB97")))
        return boost::shared_ptr<Functional> (new wB97_Functional(npoints,deriv));
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wB97X")))
        return boost::shared_ptr<Functional> (new wB97X_Functional(npoints,deriv));
}
std::string Functional::availableFunctionals()
{
    std::stringstream f;
    f << "   Available Exchange Functionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta: " << endl;
    f << "  ------------ --------- --------- ---------" << endl;
    f << "   S_X             X                     " << endl;
    f << "   B_X             X         X           " << endl;
    f << "   B88_X           X         X           " << endl;
    f << "   PBE_X           X         X           " << endl;
    f << "   PW91_X          X         X           " << endl;
    f << "   FT97B_X         X         X           " << endl;
    f << "   EDF1_X          X         X           " << endl;
    f << "   wS_X            X                     " << endl;
    f << "   wBf_X           X         X           " << endl;
    f << "   wPBEf_X         X         X           " << endl;
 
    f << endl; 
    f << "   Available Correlation Functionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta: " << endl;
    f << "  ------------ --------- --------- ---------" << endl;
    f << "   LYP_C           X         X           " << endl;
    f << "   VWN5RPA_C       X                     " << endl;
    f << "   VWN5_C          X                     " << endl;
    f << "   PZ81_C          X                     " << endl;
    f << "   P86_C           X         X           " << endl;
    f << "   PW91_C          X         X           " << endl;
    f << "   PW92_C          X                     " << endl;
    f << "   PBE_C           X         X           " << endl;
    f << "   FT97_C          X         X           " << endl;
    f << "   EDF1_C          X         X           " << endl;
    
    f << endl; 
    f << "   Available Combined Functionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta: " << endl;
    f << "  ------------ --------- --------- ---------" << endl;
    f << "   EDF1            X         X           " << endl;
    f << "   B97_0           X         X           " << endl;
    f << "   B97_1           X         X           " << endl;
    f << "   B97_2           X         X           " << endl;
    f << "   B97_D2          X         X           " << endl;
    f << "   HCTH            X         X           " << endl;
    f << "   HCTH407         X         X           " << endl;
    f << "   HCTH147         X         X           " << endl;
    f << "   HCTH120         X         X           " << endl;
    f << "   M05             X         X         X " << endl;
    f << "   M05_2X          X         X         X " << endl;
    f << "   TauHCTH         X         X         X " << endl;
    f << "   TauHCTH0        X         X         X " << endl;
    f << "   wB97            X         X           " << endl;
    f << "   wB97X           X         X           " << endl;
    
    return f.str();
    
}
std::vector<std::string> Functional::availableNames()
{
    std::vector<std::string> names;
    names.push_back("S_X");
    names.push_back("B_X");
    names.push_back("B88_X");
    names.push_back("PBE_X");
    names.push_back("PW91_X");
    names.push_back("FT97B_X");
    names.push_back("LYP_C");
    names.push_back("VWN5RPA_C");
    names.push_back("VWN5_C");
    names.push_back("PZ81_C");
    names.push_back("P86_C");
    names.push_back("PW91_C");
    names.push_back("PW92_C");
    names.push_back("PBE_C");
    names.push_back("FT97_C");
    names.push_back("EDF1");
    names.push_back("EDF1_X");
    names.push_back("EDF1_C");
    names.push_back("B97_0");
    names.push_back("B97_1");
    names.push_back("B97_2");
    names.push_back("B97_D2");
    names.push_back("HCTH");
    names.push_back("HCTH407");
    names.push_back("HCTH147");
    names.push_back("HCTH120");
    names.push_back("M05");
    names.push_back("M05_2X");
    names.push_back("TauHCTH");
    names.push_back("TauHCTH0");
    names.push_back("wS_X");
    names.push_back("wBf_X");
    names.push_back("wPBEf_X");
    names.push_back("wB97");
    names.push_back("wB97X");
    return names;
}
boost::shared_ptr<SuperFunctional> SuperFunctional::createSuperFunctional(const std::string & name, int npoints, int deriv)
{
    boost::shared_ptr<SuperFunctional> superfun = (boost::shared_ptr<SuperFunctional>) new SuperFunctional(npoints, deriv);

    //Alias Table
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("S_X"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("S_X");
        superfun->setDescription("Slater Exchange");
        superfun->setCitation("J.C. Slater, Phys. Rev., 81(3):385-390, 1951");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B88_X"))) {
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B88_X");
        superfun->setDescription("Becke 88 Exchange (GGA Only)");
        superfun->setCitation("A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B_X"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B_X");
        superfun->setDescription("Becke Exchange (S+B88)");
        superfun->setCitation("A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_X"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE_X");
        superfun->setDescription("PBE Exchange");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW91_X"))) {
        superfun->addFunctional(Functional::createFunctional("PW91_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PW91_X");
        superfun->setDescription("Perdew 91 Exchange Functional");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("FT97B_X"))) {
        superfun->addFunctional(Functional::createFunctional("FT97B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("FT97B_X");
        superfun->setDescription("Filitov and Theil 1997 Exchange B");
        superfun->setCitation("M. Filatov and W. Theil, Mol. Phys., 91(5), 847-859, 1997");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B3_X"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  8.0000000000000004E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  7.1999999999999997E-01);
        superfun->setName("B3_X");
        superfun->setDescription("B3 (80% S + 72% B88 + 20% HF)");
        superfun->setCitation("A.D. Becke, J. Chem. Phys., 98(7):5648-5652, 1993");
        superfun->setExactExchange(  2.0000000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PZ81_C"))) {
        superfun->addFunctional(Functional::createFunctional("PZ81_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PZ81_C");
        superfun->setDescription("PZ81 Correlation");
        superfun->setCitation("J.P. Perdew, A. Zunger, Phys. Rev. B., 23, 5048-5079, 1981");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("P86_C"))) {
        superfun->addFunctional(Functional::createFunctional("P86_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("P86_C");
        superfun->setDescription("P86 Correlation (PZ81 LSDA + P86 GGA)");
        superfun->setCitation("J.P. Perdew, Phys. Rev. B., 33, 8822-8824, 1986");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("VWN5_C"))) {
        superfun->addFunctional(Functional::createFunctional("VWN5_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("VWN5_C");
        superfun->setDescription("VWN5 Correlation Functional");
        superfun->setCitation("S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("VWN5RPA_C"))) {
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("VWN5RPA_C");
        superfun->setDescription("VWN5 (RPA) Correlation Functional");
        superfun->setCitation("S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("LYP_C"))) {
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("LYP_C");
        superfun->setDescription("LYP Correlation");
        superfun->setCitation("B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW91_C"))) {
        superfun->addFunctional(Functional::createFunctional("PW91_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PW91_C");
        superfun->setDescription("PW91 Correlation");
        superfun->setCitation("J.P. Perdew, et. al., Phys. Rev. B., 46(11), 6671-6687, 1992");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW92_C"))) {
        superfun->addFunctional(Functional::createFunctional("PW92_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PW92_C");
        superfun->setDescription("PW92 LSDA Correlation");
        superfun->setCitation("J.P. Perdew and Y. Wang, Phys. Rev. B., 45(23), 13244, 1992");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_C"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE_C");
        superfun->setDescription("PBE Correlation");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("FT97_C"))) {
        superfun->addFunctional(Functional::createFunctional("FT97_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("FT97_C");
        superfun->setDescription("FT97 Correlation (Involves Ei functions)");
        superfun->setCitation("M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PW91"))) {
        superfun->addFunctional(Functional::createFunctional("PW91_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("PW91_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PW91");
        superfun->setDescription("PW91 XC ");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("FT97"))) {
        superfun->addFunctional(Functional::createFunctional("FT97B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("FT97_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("FT97");
        superfun->setDescription("FT97 XC ");
        superfun->setCitation("M. Filatov and W. Theil, Int. J. Quant. Chem., 62, 603-616, 1997");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("EDF1"))) {
        superfun->addFunctional(Functional::createFunctional("EDF1",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("EDF1");
        superfun->setDescription("Empirical Density Function #1");
        superfun->setCitation("R.D. Adamson, P.M.W. Gill, and J.A. Pople, Chem. Phys. Lett., 284, 6-11, 1998");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("EDF2"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  2.8110000000000002E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  6.2270000000000003E-01);
        superfun->addFunctional(Functional::createFunctional("EDF1_X",npoints,deriv), -5.5100000000000003E-02);
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  3.0290000000000000E-01);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  5.9980000000000000E-01);
        superfun->addFunctional(Functional::createFunctional("EDF1_C",npoints,deriv), -5.3000000000000000E-03);
        superfun->setName("EDF2");
        superfun->setDescription("Empirical Density Function #2");
        superfun->setCitation("C.Y. Lin, M.W. George, and P.M.W. Gill, Aust. J. Chem., 57, 365-370, 2004");
        superfun->setExactExchange(  1.6950000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BLYP"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BLYP");
        superfun->setDescription("BLYP GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BLYP_D1"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BLYP_D1");
        superfun->setDescription("BLYP-D1 GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D1",  1.3999999999999999E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BLYP_D2"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BLYP_D2");
        superfun->setDescription("BLYP-D2 GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  1.2000000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B3LYP"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  8.0000000000000004E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  7.1999999999999997E-01);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  8.1000000000000005E-01);
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  1.9000000000000000E-01);
        superfun->setName("B3LYP");
        superfun->setDescription("B3LYP (with VWN3 Correlation)");
        superfun->setCitation("P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994");
        superfun->setExactExchange(  2.0000000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B3LYP5"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  8.0000000000000004E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  7.1999999999999997E-01);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  8.1000000000000005E-01);
        superfun->addFunctional(Functional::createFunctional("VWN5_C",npoints,deriv),  1.9000000000000000E-01);
        superfun->setName("B3LYP5");
        superfun->setDescription("B3LYP (with VWN5 Correlation)");
        superfun->setCitation("P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994");
        superfun->setExactExchange(  2.0000000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B3LYP_D2"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  8.0000000000000004E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  7.1999999999999997E-01);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  8.1000000000000005E-01);
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  1.9000000000000000E-01);
        superfun->setName("B3LYP_D2");
        superfun->setDescription("B3LYP (with VWN3 Correlation)");
        superfun->setCitation("P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994");
        superfun->setExactExchange(  2.0000000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  1.0500000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B3LYP5_D2"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  8.0000000000000004E-01);
        superfun->addFunctional(Functional::createFunctional("B88_X",npoints,deriv),  7.1999999999999997E-01);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  8.1000000000000005E-01);
        superfun->addFunctional(Functional::createFunctional("VWN5_C",npoints,deriv),  1.9000000000000000E-01);
        superfun->setName("B3LYP5_D2");
        superfun->setDescription("B3LYP (with VWN5 Correlation)");
        superfun->setCitation("P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994");
        superfun->setExactExchange(  2.0000000000000001E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  1.0500000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BP86"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("P86_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BP86");
        superfun->setDescription("BP86 GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BP86_D1"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("P86_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BP86_D1");
        superfun->setDescription("BP86-D1 GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D1",  1.3000000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("BP86_D2"))) {
        superfun->addFunctional(Functional::createFunctional("B_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("P86_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("BP86_D2");
        superfun->setDescription("BP86-D2 GGA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  1.0500000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE");
        superfun->setDescription("PBE GGA");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_D1"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE_D1");
        superfun->setDescription("PBE-D1 GGA");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D1",  6.9999999999999996E-01),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE_D2"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE_D2");
        superfun->setDescription("PBE-D2 GGA");
        superfun->setCitation("J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  7.5000000000000000E-01),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_0"))) {
        superfun->addFunctional(Functional::createFunctional("B97_0",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B97_0");
        superfun->setDescription("B97-0 Hybrid GGA");
        superfun->setCitation("A.D. Becke, J. Chem. Phys., 107(20), 8554-8560, 1997");
        superfun->setExactExchange(  1.9430000000000000E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_1"))) {
        superfun->addFunctional(Functional::createFunctional("B97_1",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B97_1");
        superfun->setDescription("B97-1 Hybrid GGA");
        superfun->setCitation("F.A. Hamprecht et. al., J. Chem. Phys., 109(15), 6264-6271, 1998");
        superfun->setExactExchange(  2.0999999999999999E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_2"))) {
        superfun->addFunctional(Functional::createFunctional("B97_2",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B97_2");
        superfun->setDescription("B97-2 Hybrid GGA");
        superfun->setCitation("P.J. Wilson et. al., J. Chem. Phys., 115(20), 9233-9242, 2001");
        superfun->setExactExchange(  2.0999999999999999E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("B97_D2"))) {
        superfun->addFunctional(Functional::createFunctional("B97_D2",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("B97_D2");
        superfun->setDescription("B97-D2 Pure GGA");
        superfun->setCitation("S. Grimme, J. Comput. Chem., 27, 1787-1799, 2006");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
        superfun->setDashD(Dispersion::createDispersion("-D2",  1.2500000000000000E+00),1.0);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH"))) {
        superfun->addFunctional(Functional::createFunctional("HCTH",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("HCTH");
        superfun->setDescription("HCTH Pure GGA");
        superfun->setCitation("F.A. Hamprecht et. al., J. Chem. Phys., 109(15), 6264-6271");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH120"))) {
        superfun->addFunctional(Functional::createFunctional("HCTH120",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("HCTH120");
        superfun->setDescription("HCTH120 Pure GGA");
        superfun->setCitation("A.D. Boese, et. al., J. Chem. Phys., 112(4), 1670-1678, 2000");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH147"))) {
        superfun->addFunctional(Functional::createFunctional("HCTH147",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("HCTH147");
        superfun->setDescription("HCTH147 Pure GGA");
        superfun->setCitation("A.D. Boese, et. al., J. Chem. Phys., 112(4), 1670-1678, 2000");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("HCTH407"))) {
        superfun->addFunctional(Functional::createFunctional("HCTH407",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("HCTH407");
        superfun->setDescription("HCTH407 Pure GGA");
        superfun->setCitation("A.D. Boese and N.C. Handy, J. Chem. Phys., 114(13), 5497-5503, 2001");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("M05"))) {
        superfun->addFunctional(Functional::createFunctional("M05",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("M05");
        superfun->setDescription("M05 Hybrid Meta-GGA");
        superfun->setCitation("Yan Zhao, Nathan E. Schultz, and D. G. Truhlar, J. Chem. Phys., 123, 161103, 2005");
        superfun->setExactExchange(  2.8000000000000003E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("M05_2X"))) {
        superfun->addFunctional(Functional::createFunctional("M05_2X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("M05_2X");
        superfun->setDescription("M05-2X Hybrid Meta-GGA");
        superfun->setCitation("Zhao, Y., Schultz, N. E., Truhlar, D. G., J. Chem. Theory Comput. 2, 364, 2006");
        superfun->setExactExchange(  5.6000000000000005E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("TauHCTH"))) {
        superfun->addFunctional(Functional::createFunctional("TauHCTH",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("TauHCTH");
        superfun->setDescription("TauHCTH Power Series Pure Meta-GGA");
        superfun->setCitation("A.D. Boese and N.C. Handy, J. Chem. Phys., 116(22), 9559, 2002");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("TauHCTH0"))) {
        superfun->addFunctional(Functional::createFunctional("TauHCTH0",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("TauHCTH0");
        superfun->setDescription("TauHCTH0 Power Series Hybrid Meta-GGA");
        superfun->setCitation("A.D. Boese and N.C. Handy, J. Chem. Phys., 116(22), 9559, 2002");
        superfun->setExactExchange(  1.4999999999999999E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wB97"))) {
        superfun->addFunctional(Functional::createFunctional("wB97",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wB97");
        superfun->setDescription("Range-Corrected B97 Pure GGA");
        superfun->setCitation("J. Chai and M. Head-Gordon, J. Chem. Phys., 128, 084106, 2008");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  4.0000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wB97X"))) {
        superfun->addFunctional(Functional::createFunctional("wB97X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wB97X");
        superfun->setDescription("Range-Corrected B97 Hybrid GGA");
        superfun->setCitation("J. Chai and M. Head-Gordon, J. Chem. Phys., 128, 084106, 2008");
        superfun->setExactExchange(  1.5770599999999999E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.0000000000000004E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wS_X"))) {
        superfun->addFunctional(Functional::createFunctional("wS_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wS_X");
        superfun->setDescription("Range-Corrected Slater LSDA");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("SVWN"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("SVWN");
        superfun->setDescription("S+VWN");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("SVWN5"))) {
        superfun->addFunctional(Functional::createFunctional("S_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("VWN5_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("SVWN5");
        superfun->setDescription("S+VWN5");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wSVWN"))) {
        superfun->addFunctional(Functional::createFunctional("wS_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("VWN5RPA_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wSVWN");
        superfun->setDescription("S+VWN");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wSVWN5"))) {
        superfun->addFunctional(Functional::createFunctional("wS_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("VWN5_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wSVWN5");
        superfun->setDescription("S+VWN5");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wBf_X"))) {
        superfun->addFunctional(Functional::createFunctional("wBf_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wBf_X");
        superfun->setDescription("Becke Exchange (S+B88), SR (Fake)");
        superfun->setCitation("A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wBLYPf"))) {
        superfun->addFunctional(Functional::createFunctional("wBf_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("LYP_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wBLYPf");
        superfun->setDescription("BLYP SR (Fake)");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wPBEf_X"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wPBEf_X");
        superfun->setDescription("wPBE_X (Fake)");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wPBEf"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  1.0000000000000000E+00);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wPBEf");
        superfun->setDescription("wPBE (Fake)");
        superfun->setCitation("Null");
        superfun->setExactExchange(  0.0000000000000000E+00);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("PBE0"))) {
        superfun->addFunctional(Functional::createFunctional("PBE_X",npoints,deriv),  7.5000000000000000E-01);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("PBE0");
        superfun->setDescription("PBE0");
        superfun->setCitation("Null");
        superfun->setExactExchange(  2.5000000000000000E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  0.0000000000000000E+00);
    }
    if (boost::to_upper_copy(name) == boost::to_upper_copy(std::string("wPBE0f"))) {
        superfun->addFunctional(Functional::createFunctional("wPBEf_X",npoints,deriv),  7.5000000000000000E-01);
        superfun->addFunctional(Functional::createFunctional("PBE_C",npoints,deriv),  1.0000000000000000E+00);
        superfun->setName("wPBE0f");
        superfun->setDescription("wPBE0 (Fake)");
        superfun->setCitation("Null");
        superfun->setExactExchange(  2.5000000000000000E-01);
        superfun->setPT2(  0.0000000000000000E+00);
        superfun->setOmega(  3.3000000000000002E-01);
    }

    return superfun;
}
std::string SuperFunctional::availableSuperFunctionals()
{
    std::stringstream f;
    
    f << "   Available Exchange SuperFunctionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta:    Hybrid:     PT2:      RC:       -D:" << endl;
    f << "  ------------ --------- --------- --------- --------- --------- --------- ---------" << endl;
    f << "   S_X             X                                                             " << endl;
    f << "   B88_X           X         X                                                   " << endl;
    f << "   B_X             X         X                                                   " << endl;
    f << "   PBE_X           X         X                                                   " << endl;
    f << "   PW91_X          X         X                                                   " << endl;
    f << "   FT97B_X         X         X                                                   " << endl;
    f << "   B3_X            X         X                   X                               " << endl;
    f << "   wS_X            X                                                 X           " << endl;
    f << "   wBf_X           X         X                                       X           " << endl;
 
    f << endl; 
    f << "   Available Correlation SuperFunctionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta:    Hybrid:     PT2:      RC:       -D:" << endl;
    f << "  ------------ --------- --------- --------- --------- --------- --------- ---------" << endl;
    f << "   PZ81_C          X                                                             " << endl;
    f << "   P86_C           X         X                                                   " << endl;
    f << "   VWN5_C          X                                                             " << endl;
    f << "   VWN5RPA_C       X                                                             " << endl;
    f << "   LYP_C           X         X                                                   " << endl;
    f << "   PW91_C          X         X                                                   " << endl;
    f << "   PW92_C          X                                                             " << endl;
    f << "   PBE_C           X         X                                                   " << endl;
    f << "   FT97_C          X         X                                                   " << endl;

    f << endl; 
    f << "   Available Combined SuperFunctionals:   " << endl;
    f << "    Name:        LSDA:      GGA:     Meta:    Hybrid:     PT2:      RC:       -D:" << endl;
    f << "  ------------ --------- --------- --------- --------- --------- --------- ---------" << endl;
    f << "   PW91            X         X                                                   " << endl;
    f << "   FT97            X         X                                                   " << endl;
    f << "   EDF1            X         X                                                   " << endl;
    f << "   EDF2            X         X                   X                               " << endl;
    f << "   BLYP            X         X                                                   " << endl;
    f << "   BLYP_D1         X         X                                                 X " << endl;
    f << "   BLYP_D2         X         X                                                 X " << endl;
    f << "   B3LYP           X         X                   X                               " << endl;
    f << "   B3LYP5          X         X                   X                               " << endl;
    f << "   B3LYP_D2        X         X                   X                             X " << endl;
    f << "   B3LYP5_D2       X         X                   X                             X " << endl;
    f << "   BP86            X         X                                                   " << endl;
    f << "   BP86_D1         X         X                                                 X " << endl;
    f << "   BP86_D2         X         X                                                 X " << endl;
    f << "   PBE             X         X                                                   " << endl;
    f << "   PBE_D1          X         X                                                 X " << endl;
    f << "   PBE_D2          X         X                                                 X " << endl;
    f << "   B97_0           X         X                   X                               " << endl;
    f << "   B97_1           X         X                   X                               " << endl;
    f << "   B97_2           X         X                   X                               " << endl;
    f << "   B97_D2          X         X                                                 X " << endl;
    f << "   HCTH            X         X                                                   " << endl;
    f << "   HCTH120         X         X                                                   " << endl;
    f << "   HCTH147         X         X                                                   " << endl;
    f << "   HCTH407         X         X                                                   " << endl;
    f << "   M05             X         X         X         X                               " << endl;
    f << "   M05_2X          X         X         X         X                               " << endl;
    f << "   TauHCTH         X         X         X                                         " << endl;
    f << "   TauHCTH0        X         X         X         X                               " << endl;
    f << "   wB97            X         X                                       X           " << endl;
    f << "   wB97X           X         X                   X                   X           " << endl;
    f << "   SVWN            X                                                             " << endl;
    f << "   SVWN5           X                                                             " << endl;
    f << "   wSVWN           X                                                 X           " << endl;
    f << "   wSVWN5          X                                                 X           " << endl;
    f << "   wBLYPf          X         X                                       X           " << endl;
    f << "   wPBEf_X         X         X                                       X           " << endl;
    f << "   wPBEf           X         X                                       X           " << endl;
    f << "   PBE0            X         X                   X                               " << endl;
    f << "   wPBE0f          X         X                   X                   X           " << endl;

    return f.str();
}
std::vector<std::string> SuperFunctional::availableNames()
{
    std::vector<std::string> names;
    names.push_back("S_X");
    names.push_back("B88_X");
    names.push_back("B_X");
    names.push_back("PBE_X");
    names.push_back("PW91_X");
    names.push_back("FT97B_X");
    names.push_back("B3_X");
    names.push_back("PZ81_C");
    names.push_back("P86_C");
    names.push_back("VWN5_C");
    names.push_back("VWN5RPA_C");
    names.push_back("LYP_C");
    names.push_back("PW91_C");
    names.push_back("PW92_C");
    names.push_back("PBE_C");
    names.push_back("FT97_C");
    names.push_back("PW91");
    names.push_back("FT97");
    names.push_back("EDF1");
    names.push_back("EDF2");
    names.push_back("BLYP");
    names.push_back("BLYP_D1");
    names.push_back("BLYP_D2");
    names.push_back("B3LYP");
    names.push_back("B3LYP5");
    names.push_back("B3LYP_D2");
    names.push_back("B3LYP5_D2");
    names.push_back("BP86");
    names.push_back("BP86_D1");
    names.push_back("BP86_D2");
    names.push_back("PBE");
    names.push_back("PBE_D1");
    names.push_back("PBE_D2");
    names.push_back("B97_0");
    names.push_back("B97_1");
    names.push_back("B97_2");
    names.push_back("B97_D2");
    names.push_back("HCTH");
    names.push_back("HCTH120");
    names.push_back("HCTH147");
    names.push_back("HCTH407");
    names.push_back("M05");
    names.push_back("M05_2X");
    names.push_back("TauHCTH");
    names.push_back("TauHCTH0");
    names.push_back("wB97");
    names.push_back("wB97X");
    names.push_back("wS_X");
    names.push_back("SVWN");
    names.push_back("SVWN5");
    names.push_back("wSVWN");
    names.push_back("wSVWN5");
    names.push_back("wBf_X");
    names.push_back("wBLYPf");
    names.push_back("wPBEf_X");
    names.push_back("wPBEf");
    names.push_back("PBE0");
    names.push_back("wPBE0f");
    return names;
}

}}

