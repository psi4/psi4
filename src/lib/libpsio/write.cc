/*!
 \file
 \ingroup PSIO
 */

#include <cstdlib>
#include <cstring>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::write(unsigned int unit, const char *key, char *buffer, ULI size,
                 psio_address start, psio_address *end) {
  psio_ud *this_unit;
  psio_tocentry *this_entry, *last_entry;
  psio_address start_toc, start_data, end_data; /* global addresses */
  ULI tocentry_size;
  int dirty = 0;

  this_unit = &(psio_unit[unit]);

  /* Find the entry in the TOC */
  this_entry = tocscan(unit, key);

  tocentry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);

  if (this_entry == NULL) { /* New TOC entry */
    if (start.page||start.offset)
      psio_error(unit, PSIO_ERROR_BLKSTART);

    dirty = 1; /* set flag for writing the TOC header */

    this_entry = (psio_tocentry *) malloc(sizeof(psio_tocentry));
    ::strncpy(this_entry->key, key, PSIO_KEYLEN);
    this_entry->key[PSIO_KEYLEN-1] = '\0';
    this_entry->next = NULL;
    this_entry->last = NULL;

    /* Compute the global address of the new entry */
    if (!(this_unit->toclen)) { /* First TOC entry */
      this_entry->sadd.page = 0;
      this_entry->sadd.offset = sizeof(ULI); /* offset for the toclen value stored first */
      this_unit->toc = this_entry;
    } else { /* Use ending address from last TOC entry */
      last_entry = toclast(unit);
      this_entry->sadd = last_entry->eadd;
      last_entry->next = this_entry;
      this_entry->last = last_entry;
    }

    /* compute important global addresses for the entry */
    start_toc = this_entry->sadd;
    start_data = psio_get_address(start_toc, tocentry_size);
    start_data = psio_get_global_address(start_data, start);
    end_data = psio_get_address(start_data, size);

    /* Set the end address for this_entry */
    this_entry->eadd = end_data;

    /* Update the unit's TOC stats */
    this_unit->toclen++;
    wt_toclen(unit, this_unit->toclen);

    /* Update end (an entry-relative address) for the caller */
    *end = psio_get_address(start, size);
  } else { /* Old TOC entry */

    /* Compute the global starting page and offset for the block */
    start_toc = this_entry->sadd;
    start_data = psio_get_address(start_toc, tocentry_size);
    start_data = psio_get_global_address(start_data, start);

    /* Make sure this block doesn't start past the end of the entry */
    if (start_data.page > this_entry->eadd.page)
      psio_error(unit, PSIO_ERROR_BLKSTART);
    else if ((start_data.page == this_entry->eadd.page) &&(start_data.offset
        > this_entry->eadd.offset))
      psio_error(unit, PSIO_ERROR_BLKSTART);

    /* Compute the new global ending address for the entry, if necessary */
    end_data = psio_get_address(start_data, size);
    if (end_data.page > this_entry->eadd.page) {
      if (this_entry->next != NULL) {
        fprintf(stderr, "PSIO_ERROR: Attempt to write into next entry: %d, %s\n", unit, key);
        psio_error(unit, PSIO_ERROR_BLKEND);
      }
      this_entry->eadd = end_data;
      dirty = 1; /* set flag for writing the TOC header */
    } else if ((end_data.page == this_entry->eadd.page) &&(end_data.offset
        > this_entry->eadd.offset)) {
      if (this_entry->next != NULL) {
        fprintf(stderr, "PSIO_ERROR: Attempt to write into next entry: %d, %s\n", unit, key);
        psio_error(unit, PSIO_ERROR_BLKEND);
      }
      this_entry->eadd = end_data;
      dirty = 1; /* set flag for writing the TOC header */
    }

    /* Update end (an entry-relative address) for the caller */
    *end = psio_get_address(start, size);
  }

  if (dirty) /* Need to first write/update the TOC header for this record */
    rw(unit, (char *) this_entry, start_toc, tocentry_size, 1);

  /* Now write the actual data to the unit */
  rw(unit, buffer, start_data, size, 1);

#ifdef PSIO_STATS
  psio_writlen[unit] += size;
#endif
}

  /*!
   ** PSIO_WRITE(): Writes data to a TOC entry in a PSI file.
   **
   **  \param unit    = The PSI unit number used to identify the file to all read
   **                   and write functions.
   **  \param key     = The TOC keyword identifying the desired entry.
   **  \param buffer  = The buffer from which the data is written.
   **  \param size    = The number of bytes to write.
   **  \param start   = The entry-relative starting page/offset to write the data.
   **  \param end     = A pointer to the entry-relative page/offset for the next
   **                   byte after the end of the write request.
   **
   ** \ingroup PSIO
   */

  int psio_write(unsigned int unit, const char *key, char *buffer, ULI size,
                 psio_address start, psio_address *end) {
    _default_psio_lib_->write(unit, key, buffer, size, start, end);
    return 1;
  }

}

