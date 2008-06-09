/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libbasis_osrecur_h_
#define _psi_src_lib_libbasis_osrecur_h_

#include <psitypes.h>

namespace psi {

class OI_OSRecursor {

  int maxam1_;
  int maxam2_;

  PSI_FLOAT** OIX_;
  PSI_FLOAT** OIY_;
  PSI_FLOAT** OIZ_;

  // No default constructor
  OI_OSRecursor();
  // No assignment operator
  OI_OSRecursor& operator=(const OI_OSRecursor&);

 public:
  OI_OSRecursor(int, int);
  ~OI_OSRecursor();

  PSI_FLOAT** OIX() const { return OIX_; };
  PSI_FLOAT** OIY() const { return OIY_; };
  PSI_FLOAT** OIZ() const { return OIZ_; };
  void compute(PSI_FLOAT PA[3], PSI_FLOAT PB[3], PSI_FLOAT gamma, int am1, int am2);
};

} // end of namespace psi

#endif
