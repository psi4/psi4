/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

    extern FILE *outfile;

void PSIO::tocprint(unsigned int unit) {
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  fprintf(outfile, "\nTable of Contents for Unit %5u\n", unit);
  fprintf(
          outfile,
          "----------------------------------------------------------------------------\n");
  fprintf(
          outfile,
          "Key                                   Spage    Soffset      Epage    Eoffset\n");
  fprintf(
          outfile,
          "----------------------------------------------------------------------------\n");

  while (this_entry != NULL) {
    fprintf(outfile, "%-32s %10lu %10lu %10lu %10lu\n", this_entry->key,
            this_entry->sadd.page, this_entry->sadd.offset,
            this_entry->eadd.page, this_entry->eadd.offset);
    this_entry = this_entry->next;
  }
  fprintf(outfile, "\n");
  fflush(outfile);
}

  /*!
   ** PSIO_TOCPRINT(): Print the table of contents for the given unit
   **
   ** \ingroup PSIO
   */

  void psio_tocprint(unsigned int unit) {
    return _default_psio_lib_->tocprint(unit);
  }


}

