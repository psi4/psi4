/*!
 \file
 \ingroup PSIO
 */

#include <libpsio/psio.h>

namespace psi {
  /*!
   ** PSIO_GET_LENGTH(): Given a start page and offset for two data sets,
   ** compute the number of bytes between them.  Note that eadd denotes the
   ** beginning of the next entry and not the end of the current entry.
   **
   ** \ingroup PSIO
   */

  ULI psio_get_length(psio_address sadd, psio_address eadd) {
    
    ULI full_page_bytes;
    
    /* Number of bytes on fullpages */
    full_page_bytes = (eadd.page - sadd.page- 1)*PSIO_PAGELEN;
    
    if (full_page_bytes < 0) { /* We're on a single page */
      return (eadd.offset - sadd.offset);
    } else if (full_page_bytes == 0) { /* We're on the next page */
      return ((PSIO_PAGELEN - sadd.offset) + eadd.offset);
    } else {
      return ((PSIO_PAGELEN - sadd.offset) + full_page_bytes + eadd.offset);
    }
  }

}

