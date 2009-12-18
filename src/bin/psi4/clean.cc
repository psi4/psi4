/*! \file clean
    \defgroup PSI4
*/

#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

namespace psi {

/*!
** psiclean():
**
** Remove all psio files.  This should serve as a more elegant replacement
** for the PSI3 utility psiclean, which executed system calls to delete 
** the files.
*/

void psiclean(void) {

  ip_cwk_clear();
  ip_cwk_add(const_cast<char*>(":DEFAULT"));
  ip_cwk_add(const_cast<char*>(":PSI"));
  ip_set_uppercase(1);

  for(int i=0; i < PSIO_MAXUNIT; i++) {
    psio_open(i, PSIO_OPEN_OLD);
    psio_close(i, 0);
  }
}

} //end ::psi

