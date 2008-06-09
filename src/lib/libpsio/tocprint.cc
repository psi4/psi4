/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::tocprint(unsigned int unit, FILE *output) {
  psio_tocentry *this_entry;
  
  this_entry = psio_unit[unit].toc;
  
  fprintf(output, "\nTable of Contents for Unit %5u\n", unit);
  fprintf(
          output,
          "----------------------------------------------------------------------------\n");
  fprintf(
          output,
          "Key                                   Spage    Soffset      Epage    Eoffset\n");
  fprintf(
          output,
          "----------------------------------------------------------------------------\n");
  
  while (this_entry != NULL) {
    fprintf(output, "%-32s %10lu %10lu %10lu %10lu\n", this_entry->key,
            this_entry->sadd.page, this_entry->sadd.offset,
            this_entry->eadd.page, this_entry->eadd.offset);
    this_entry = this_entry->next;
  }
  fprintf(output, "\n");
  fflush(output);
}

  /*!
   ** PSIO_TOCPRINT(): Print the table of contents for the given unit
   **
   ** \ingroup PSIO
   */

  void psio_tocprint(unsigned int unit, FILE *output) {
    return _default_psio_lib_->tocprint(unit, output);
  }


}

