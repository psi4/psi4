/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <boost/shared_ptr.hpp>

namespace psi {

/* Definition of global data */
shared_ptr<PSIO> _default_psio_lib_;

int PSIO::_error_exit_code_ = 1;
psio_address PSIO_ZERO = { 0, 0 };

PSIO::PSIO() {
  int i, j;
  
  psio_unit = (psio_ud *) malloc(sizeof(psio_ud)*PSIO_MAXUNIT);
#ifdef PSIO_STATS
  psio_readlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
  psio_writlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
#endif
  state_ = 1;
  
  if (psio_unit == NULL) {
    fprintf(stderr, "Error in PSIO_INIT()!\n");
    exit(_error_exit_code_);
  }
  
  for (i=0; i < PSIO_MAXUNIT; i++) {
#ifdef PSIO_STATS
    psio_readlen[i] = psio_writlen[i] = 0;
#endif      
    psio_unit[i].numvols = 0;
    for (j=0; j < PSIO_MAXVOL; j++) {
      psio_unit[i].vol[j].path = NULL;
      psio_unit[i].vol[j].stream = -1;
    }
    psio_unit[i].toclen = 0;
    psio_unit[i].toc = NULL;
  }
}

  int psio_init(void) {
    if (_default_psio_lib_.get() == 0) {
      shared_ptr<PSIO> temp(new PSIO);
      _default_psio_lib_ = temp;
      if (_default_psio_lib_ == 0) {
        fprintf(stderr,"LIBPSIO::init() -- failed to allocate the memory");
        exit(PSIO::_error_exit_code_);
      }
    }

    return 1;
  }

  int psio_state() {
    return _default_psio_lib_->state();
  }

}

