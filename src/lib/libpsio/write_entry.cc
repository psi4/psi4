/*!
 \file write_entry.cc
 \ingroup (PSIO)
 */

#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

void PSIO::write_entry(unsigned int unit, char *key, char *buffer, ULI size) {
  psio_address end;
  write(unit, key, buffer, size, PSIO_ZERO, &end);
}

extern "C" {
  
  /*!
   ** PSIO_WRITE_ENTRY()
   **
   ** \ingroup (PSIO)
   */
  int psio_write_entry(unsigned int unit, char *key, char *buffer, ULI size) {
    psio_address end;
    return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
  }
}

