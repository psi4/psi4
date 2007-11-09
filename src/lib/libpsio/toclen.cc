/*!
 \file toclen.cc
 \ingroup (PSIO)
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

unsigned int PSIO::toclen(unsigned int unit) {
  unsigned int toclen=0;
  psio_tocentry *this_entry;
  
  this_entry = psio_unit[unit].toc;
  
  while (this_entry != NULL) {
    ++toclen;
    this_entry = this_entry->next;
  }
  
  return (toclen);
}

ULI PSIO::rd_toclen(unsigned int unit) {
  int errcod, stream;
  psio_ud *this_unit;
  ULI toclen;
  
  this_unit = &(psio_unit[unit]);
  
  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if (errcod == -1)
    psio_error(unit, PSIO_ERROR_LSEEK);
  
  /* Read the value */
  errcod =:: read(stream, (char *) &toclen, sizeof(ULI));
  if(errcod != sizeof(ULI)) return(0); /* assume that all is well (see comments above) */

  return(toclen);
}

void PSIO::wt_toclen(unsigned int unit, ULI toclen) {
  int errcod, stream;
  psio_ud *this_unit;
  
  this_unit = &(psio_unit[unit]);
  
  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if (errcod == -1) {
    fprintf(stderr, "Error in PSIO_WT_TOCLEN()!\n");
    exit(_error_exit_code_);
  }
  
  /* Write the value */
  errcod =:: write(stream, (char *) &toclen, sizeof(ULI));
  if(errcod != sizeof(ULI)) {
    fprintf(stderr, "PSIO_ERROR: Failed to write toclen to unit %d.\n", unit);
    psio_error(unit,PSIO_ERROR_WRITE);
  }
}

#if 0
      /*!
       ** PSIO_TOCLEN(): Compute the length of the TOC for a given unit using the in-core TOC list.
       **
       ** \ingroup (PSIO)
       */

      unsigned int psio_toclen(unsigned int unit)
      {
        return __psio_toclen(_default_psio_lib_,unit);
      }

      /*!
       ** PSIO_RD_TOCLEN(): Read the length of the TOC for a given unit directly from the file.
       **
       ** \param unit = PSI unit number from which to read the toclen.
       **
       ** NB: Note that we do not exit if the read request of the toclen from
       ** the file fails. This is because the request may be to an new file
       ** for which the toclen has not yet been written.  (We allow the user
       ** to open files with status PSIO_OPEN_OLD even if they don't exist,
       ** because sometimes you can't know this in advance.)
       **
       ** \ingroup (PSIO)
       */
      ULI psio_rd_toclen(unsigned int unit)
      {
        return __psio_rd_toclen(_default_psio_lib_,unit);
      }

      /*!
       ** PSIO_WT_TOCLEN(): Write the length of the TOC for a given unit directly to the file.
       **
       ** \param unit = PSI unit number to which to write the toclen.
       **
       ** \ingroup (PSIO)
       */
      void psio_wt_toclen(unsigned int unit, ULI toclen)
      {
        return __psio_wt_toclen(_default_psio_lib_,unit,toclen);
      }
#endif
      
