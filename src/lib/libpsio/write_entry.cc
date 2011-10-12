/*!
 \file
 \ingroup PSIO
 */

#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::write_entry(unsigned int unit, const char *key, char *buffer, ULI size) {
  psio_address end = PSIO_ZERO;
  write(unit, key, buffer, size, PSIO_ZERO, &end);
}

  /*!
   ** PSIO_WRITE_ENTRY()
   **
   ** \ingroup PSIO
   */

  int psio_write_entry(unsigned int unit, const char *key, char *buffer, ULI size) {
    _default_psio_lib_->write_entry(unit, key, buffer, size);
    return 1;
  }

}

