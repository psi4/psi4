/*!
 \file
 \ingroup PSIO
 */

#include <cstring>
#include <cstdlib>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::tocclean(unsigned int unit, const char *key) {
  psio_tocentry *this_entry, *last_entry, *prev_entry;
  psio_ud *this_unit;

  this_unit = &(psio_unit[unit]);

  this_entry = tocscan(unit, key);
  if (this_entry == NULL) {
    if (!strcmp(key, ""))
      this_entry = this_unit->toc;
    else {
      fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s in unit %d\n", key, unit);
      psio_error(unit, PSIO_ERROR_NOTOCENT);
    }
  } else
    this_entry = this_entry->next;

  /* Get the end of the TOC and work backwards */
  last_entry = toclast(unit);

  while ((last_entry != this_entry) && (last_entry != NULL)) {
    /* Now free all the remaining members */
    prev_entry = last_entry->last;
    free(last_entry);
    last_entry = prev_entry;
    this_unit->toclen--;
  }

  /* Update on disk */
  wt_toclen(unit, this_unit->toclen);
  tocwrite(unit);
}

  /*!
   ** PSIO_TOCCLEAN(): Delete all TOC entries after the given key.
   ** If a blank key is given, the entire TOC will be wiped.
   **
   ** \ingroup PSIO
   */

  int psio_tocclean(unsigned int unit, const char *key) {
    _default_psio_lib_->tocclean(unit, key);
    return 0;
  }

}

