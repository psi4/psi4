/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <unistd.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>

namespace psi {

void PSIO::rw(unsigned int unit, char *buffer, psio_address address, ULI size,
              int wrt) {
  int errcod;
  unsigned int i;
  ULI errcod_uli;
  ULI page, offset;
  ULI buf_offset;
  ULI this_page, this_page_max, this_page_total;
  unsigned int first_vol, this_vol, numvols;
  ULI bytes_left, num_full_pages;
  psio_ud *this_unit;
  
  this_unit = &(psio_unit[unit]);
  numvols = this_unit->numvols;
  page = address.page;
  offset = address.offset;
  
  /* Seek all volumes to correct starting positions */
  first_vol = page % numvols;
  errcod = psio_volseek(&(this_unit->vol[first_vol]), page, offset, numvols);
  if (errcod == -1)
    psio_error(unit, PSIO_ERROR_LSEEK);
  for (i=1, this_page=page+1; i < numvols; i++, this_page++) {
    this_vol = this_page % numvols;
    errcod = psio_volseek(&(this_unit->vol[this_vol]), this_page, (ULI) 0,
                          numvols);
    if (errcod == -1)
      psio_error(unit, PSIO_ERROR_LSEEK);
  }
  
  /* Number of bytes left on the first page */
  this_page_max = PSIO_PAGELEN - offset;
  
  /* If we have enough room on this page, use it */
  if (size <= this_page_max)
    this_page_total = size;
  else
    this_page_total = this_page_max;
  
  buf_offset = 0;
  if (wrt) {
	if (Communicator::world->me() == 0) {
      errcod_uli =:: write(this_unit->vol[first_vol].stream, &(buffer[buf_offset]),
          this_page_total);
	}
        Communicator::world->bcast(&errcod_uli, 1, 0);
        //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
    if(errcod_uli != this_page_total) psio_error(unit,PSIO_ERROR_WRITE);
  }
  else {
	if (Communicator::world->me() == 0) {
      errcod_uli = ::read(this_unit->vol[first_vol].stream, &(buffer[buf_offset]),
          this_page_total);
        }
        Communicator::world->bcast(&errcod_uli, 1, 0);
        //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
    if(errcod_uli != this_page_total)
      psio_error(unit,PSIO_ERROR_READ);
    else
      Communicator::world->bcast(&(buffer[buf_offset]), this_page_total, 0);
      //Communicator::world->raw_bcast(&(buffer[buf_offset]), this_page_total, 0);

  }

  /* Total number of bytes remaining to be read/written */
  bytes_left = size - this_page_total;

  /* Read/Write all the full pages */
  num_full_pages = bytes_left/PSIO_PAGELEN;
  buf_offset += this_page_total;
  for(i=0,this_page=page+1; i < num_full_pages; i++,this_page++) {
    this_vol = this_page % numvols;
    this_page_total = PSIO_PAGELEN;
    if(wrt) {
      if (Communicator::world->me() == 0) {
        errcod_uli = ::write(this_unit->vol[this_vol].stream, &(buffer[buf_offset]),
            this_page_total);
      }
      Communicator::world->bcast(&errcod_uli, 1, 0);
      //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
      if(errcod_uli != this_page_total) psio_error(unit,PSIO_ERROR_WRITE);
    }
    else {
      if (Communicator::world->me() == 0) {
        errcod_uli = ::read(this_unit->vol[this_vol].stream, &(buffer[buf_offset]),
            this_page_total);
      }
      Communicator::world->bcast(&errcod_uli, 1, 0);
      //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
      if(errcod_uli != this_page_total)
        psio_error(unit,PSIO_ERROR_READ);
      else
        Communicator::world->bcast(&(buffer[buf_offset]), this_page_total, 0);
        //Communicator::world->raw_bcast(&(buffer[buf_offset]), this_page_total, 0);

    }
    buf_offset += this_page_total;
  }

  /* Read/Write the final partial page */
  bytes_left -= num_full_pages * PSIO_PAGELEN;
  this_vol = this_page % numvols;
  if(bytes_left) {
    if(wrt) {
      if (Communicator::world->me() == 0) {
        errcod_uli = ::write(this_unit->vol[this_vol].stream, &(buffer[buf_offset]),
            bytes_left);
      }
      Communicator::world->bcast(&errcod_uli, 1, 0);
      //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
      if(errcod_uli != bytes_left) psio_error(unit,PSIO_ERROR_WRITE);
    }
    else {
      if (Communicator::world->me() == 0) {
        errcod_uli = ::read(this_unit->vol[this_vol].stream, &(buffer[buf_offset]),
            bytes_left);
      }
      Communicator::world->bcast(&errcod_uli, 1, 0);
      //Communicator::world->raw_bcast(&errcod_uli, sizeof(ULI), 0);
      if(errcod_uli != bytes_left)
        psio_error(unit,PSIO_ERROR_READ);
      else
        Communicator::world->bcast(&(buffer[buf_offset]), bytes_left, 0);
        //Communicator::world->raw_bcast(&(buffer[buf_offset]), bytes_left, 0);
    }
  }
}

  /*!
   ** PSIO_RW(): Central function for all reads and writes on a PSIO unit.
   **
   ** \params unit    = The PSI unit number.
   ** \params buffer  = The buffer containing the bytes for the read/write event.
   ** \params address = the PSIO global address for the start of the read/write.
   ** \params size    = The number of bytes to read/write.
   ** \params         = Indicates if the call is to read (0) or write (0) the input data.
   **
   ** \ingroup PSIO
   */

  int psio_rw(unsigned int unit, char *buffer, psio_address address, ULI size,
              int wrt) {
    _default_psio_lib_->rw(unit, buffer, address, size, wrt);
    return 1;
  }

}

