/*!
 \file toclast.cc
 \ingroup (PSIO)
 */

#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

psio_tocentry*PSIO::toclast(unsigned int unit) {
  psio_tocentry *this_entry = psio_unit[unit].toc;
  
  while (this_entry->next != NULL)
    this_entry = this_entry->next;
  
  return (this_entry);
}
