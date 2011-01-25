/*! \file clean
    \defgroup PSI4
*/

#include <libpsio/psio.h>

namespace psi {

    extern FILE *outfile;

/*!
** psiclean():
**
** Remove all psio files.  This should serve as a more elegant replacement
** for the PSI3 utility psiclean, which executed system calls to delete
** the files.
*/

void psiclean(void) {

  fprintf(outfile, "psiclean: Cleaning up files.\n");
//  ip_cwk_clear();
//  ip_cwk_add(const_cast<char*>(":DEFAULT"));
//  ip_cwk_add(const_cast<char*>(":PSI"));
//  ip_set_uppercase(1);

  //for(int i=0; i < 32; i++) { // skip the checkpoint file by default
  // skip checkpoint and optimization data files by default
  psio_open(0, PSIO_OPEN_OLD);
  psio_close(0, 0);
  for(int i=2; i < 32; i++) {
    psio_open(i, PSIO_OPEN_OLD);
    psio_close(i, 0);
  }
  for(int i=33; i < PSIO_MAXUNIT; i++) {
    psio_open(i, PSIO_OPEN_OLD);
    psio_close(i, 0);
  }
}

} //end ::psi

