/**********************************************************
* dispersion.cc: definitions for -D for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include "dispersion.h"
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {


Dispersion::Dispersion() 
{
}
Dispersion::~Dispersion()
{
}
boost::shared_ptr<Dispersion> Dispersion::createDispersion(const std::string & name)
{
    boost::shared_ptr<Dispersion> disp;
    return disp;
}
D1::D1() 
{
    name_ = "-D1";
    description_ = "Grimme's -D1 Dispersion Correction";
    citation_ = "Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473"; 
}
D1::~D1() 
{
}
double D1::computeEnergy(shared_ptr<Molecule> mol)
{
    double energy = 0.0;
    return energy;
}
shared_ptr<Matrix> D1::computeGradient(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> grad;
    return grad;
}
shared_ptr<Matrix> D1::computeHessian(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> hess;
    return hess;
}
D2::D2() 
{
    name_ = "-D2";
    description_ = "Grimme's -D2 Dispersion Correction";
    citation_ = "Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799";
}
D2::~D2() 
{
}
double D2::computeEnergy(shared_ptr<Molecule> mol)
{
    double energy = 0.0;
    return energy;
}
shared_ptr<Matrix> D2::computeGradient(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> grad;
    return grad;
}
shared_ptr<Matrix> D2::computeHessian(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> hess;
    return hess;
}
D3::D3() 
{
    name_ = "-D3";
    description_ = "Grimme's -D3 Dispersion Correction";
    citation_ = "Grimme, S. (2010),  J.C.P., 132: 154104";
}
D3::~D3() 
{
}
double D3::computeEnergy(shared_ptr<Molecule> mol)
{
    double energy = 0.0;
    return energy;
}
shared_ptr<Matrix> D3::computeGradient(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> grad;
    return grad;
}
shared_ptr<Matrix> D3::computeHessian(shared_ptr<Molecule> mol)
{
    shared_ptr<Matrix> hess;
    return hess;
}

}}
