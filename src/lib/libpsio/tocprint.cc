/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

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

