/*!
 \file
 \ingroup PSIO
 */

#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols) {

      double* buf = new double[cols];
      memset(static_cast<void*>(buf),'\0',cols*sizeof(double));

      psio_address next_psio = PSIO_ZERO;
      for (int i=0; i<rows; i++) {
          PSIO::write(unit,key,(char *) (buf),
          sizeof(double)*cols,next_psio,&next_psio);
      }

      delete[] buf;
}

  /*!
   ** PSIO_ZERO_DISK()
   **
   ** \ingroup PSIO
   */

  int psio_zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols) {
    _default_psio_lib_->zero_disk(unit, key, rows, cols);
    return 1;
  }

}

