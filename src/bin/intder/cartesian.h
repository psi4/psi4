/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_cartesian_h_
#define _psi_bin_intder_cartesian_h_

namespace psi { namespace intder {

class Cartesian
{
  double x, y, z;

public:
  Cartesian(double x1, double y1, double z1)
    { x = x1; y = y1; z = z1; }
  Cartesian()
    { x = 0.0; y = 0.0; z = 0.0; }

  double Distance(Cartesian& B);
  double Distance(double bx, double by, double bz);

  double& getX()
    { return x; }
  double& getY()
    { return y; }
  double& getZ()
    { return z; }
};

}} // namespace psi::intder

#endif // header guard

