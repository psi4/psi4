/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/

#include <string.h>

namespace psi { namespace response {

struct Params {
  int print;             /* Output level control */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;          /* cacheing level for libdpd */
  int ref;               /* reference determinant (0=RHF, 1=ROHF, 2=UHF) */
  double *omega;         /* energy of applied field (a.u) for dynamic properties */
  int nomega;            /* number of field energies desired */
  std::string prop;      /* desired property */
};

}} // namespace psi::response
