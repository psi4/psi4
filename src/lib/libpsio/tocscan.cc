/*!
 \file
 \ingroup PSIO
 */

#include <cstring>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

psio_tocentry*PSIO::tocscan(unsigned int unit, char *key) {
  psio_tocentry *this_entry;
  
  if (key == NULL)
    return (NULL);
  
  if ((strlen(key)+1) > PSIO_KEYLEN)
    psio_error(unit, PSIO_ERROR_KEYLEN);
  
  this_entry = psio_unit[unit].toc;
  
  while (this_entry != NULL) {
    if (!strcmp(this_entry->key, key))
      return (this_entry);
    this_entry = this_entry->next;
  }
  
  return (NULL);
}

  /*!
   ** PSIO_TOCSCAN(): Scans the TOC for a particular keyword and returns either
   ** a pointer to the entry or NULL to the caller.
   **
   ** \ingroup PSIO
   */

  psio_tocentry *psio_tocscan(unsigned int unit, char *key) {
    return _default_psio_lib_->tocscan(unit, key);
  }

}

