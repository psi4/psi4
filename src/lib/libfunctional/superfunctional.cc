/**********************************************************
* superfunctional.cc: definitions for superfunctionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include "superfunctional.h"
#include <boost/algorithm/string.hpp>
#include <libmints/properties.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {


SuperFunctional::SuperFunctional(int npoints, int deriv) : npoints_(npoints), deriv_(deriv),
    exact_exchange_(0.0), pt2_(0.0), dashD_weight_(0.0), omega_(0.0)
{
}
SuperFunctional::~SuperFunctional()
{
}
boost::shared_ptr<SuperFunctional> SuperFunctional::buildSuperFunctional(const std::string & build, int npoints, int deriv)
{
    boost::shared_ptr<SuperFunctional> superfun = (boost::shared_ptr<SuperFunctional>) new SuperFunctional(npoints,deriv);
    return superfun;
}
void SuperFunctional::reallocate()
{
    for (int A = 0; A < functionals_.size(); A++) {
        functionals_[A].first->setNPoints(npoints_);
        functionals_[A].first->setDeriv(deriv_);
    }
}
void SuperFunctional::setParameter(const std::string & name, const std::string & param, double val)
{
    for (int A = 0; A < functionals_.size(); A++)
        if (to_upper_copy(functionals_[A].first->getName()) == to_upper_copy(name)) {
            functionals_[A].first->setParameter(param,val);
            break;
        }
}
bool SuperFunctional::isGGA()
{
    for (int A = 0; A < functionals_.size(); A++)
        if (functionals_[A].first->isGGA())
            return true;

    return false;
}
bool SuperFunctional::isMeta()
{
    for (int A = 0; A < functionals_.size(); A++)
        if (functionals_[A].first->isMeta())
            return true;

    return false;
}
void SuperFunctional::computeRKSFunctional(boost::shared_ptr<Properties> props) 
{
    for (int A = 0; A < functionals_.size(); A++)
        functionals_[A].first->computeRKSFunctional(props);
}
void SuperFunctional::computeUKSFunctional(boost::shared_ptr<Properties> props) 
{
    for (int A = 0; A < functionals_.size(); A++)
        functionals_[A].first->computeRKSFunctional(props);
}
void SuperFunctional::addFunctional(const boost::shared_ptr<Functional> & f, double weight)
{
    functionals_.push_back(make_pair(f, weight));
}
void SuperFunctional::addFunctional(const std::string & name, double weight)
{
    functionals_.push_back(make_pair(Functional::createFunctional(name, npoints_, deriv_), weight));
}
void SuperFunctional::setDashD(shared_ptr<Dispersion> disp, double weight)
{
    dashD_ = disp;
    dashD_weight_ = weight;
}
std::string SuperFunctional::getComposition()
{
    std::stringstream s;
    s << "  SuperFunctional " << name_ << endl;
    s << "  Composed of:" << endl;
    for (int A = 0; A < functionals_.size(); A++) {
        s.precision(4);
        s << "  " << functionals_[A].second*100.0 << "% " << functionals_[A].first->getName() << endl;
    }
    if (isHybrid() || isDoubleHybrid() || isDashD() || isRangeCorrected())
        s << "  Additionally composed of: " << endl;
    if (isHybrid())
        s << "  " << exact_exchange_*100.0 << "% Exact Exchange" << endl;
    if (isDoubleHybrid())
        s << "  " << pt2_*100.0 << "% PT2 Long Range Correlation" << endl;
    if (isDashD())
        s << "  " << dashD_weight_*100.0 << "% " << dashD_->getName() << endl;
    if (isRangeCorrected())
        s << "  Range-Correction with omega of " << omega_ << endl; 
        
    return s.str();
}

}}
