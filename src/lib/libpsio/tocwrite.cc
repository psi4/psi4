/*!
 \file
 \ingroup PSIO
 */

#include <cstdlib>
#include <unistd.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::tocwrite(unsigned int unit) {
  unsigned int i;
  psio_ud *this_unit;
  psio_tocentry *this_entry;
  ULI entry_size;
  psio_address address;

  this_unit = &(psio_unit[unit]);
  entry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);

  if (!open_check(unit))
    return;

  wt_toclen(unit, this_unit->toclen);

  this_entry = this_unit->toc;
  address = psio_get_address(PSIO_ZERO, sizeof(ULI));
  for (i=0; i < this_unit->toclen; i++) {
    rw(unit, (char *) this_entry, address, entry_size, 1);
    this_entry = this_entry->next;
    if (this_entry != NULL)
      address = this_entry->sadd;
  }
}

  /*!
   ** PSIO_TOCWRITE(): Write the table of contents for file number 'unit'.
   **
   ** \param unit  = The PSI unit to which we will write the TOC.
   **
   ** NB: This function should NOT call psio_error because the latter calls it!
   **
   ** \ingroup PSIO
   */

  int psio_tocwrite(unsigned int unit) {
    _default_psio_lib_->tocwrite(unit);
    return 1;
  }

}

