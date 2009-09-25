/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include <cmath>
#include "cartesian.h"

using namespace psi::intder;

double Cartesian::Distance(Cartesian& B)
{
  double result = 0.0;

  result = sqrt(pow(B.x - x, 2.0) + pow(B.y - y, 2.0) + pow(B.z - z, 2.0));
  return result;
}

double Cartesian::Distance(double bx, double by, double bz)
{
  double result = 0.0;

  result = sqrt(pow(bx - x, 2.0) + pow(by - y, 2.0) + pow(bz - z, 2.0));
  return result;
}

