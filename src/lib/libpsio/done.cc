/*!
 ** \file
 ** \ingroup PSIO
 */

#include <cstdlib>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

#ifdef PSIO_STATS
#include <ctime>
#endif

namespace psi {

PSIO::~PSIO() {
#ifdef PSIO_STATS
  int i;
  ULI total_read=0, total_write=0;
  FILE *io_out;
  time_t my_time;
  my_time = time(NULL);
  io_out = fopen("io.dat", "a+");
  fprintf(io_out, "\nLIBPSIO Read/Write Statistics\n\n");
  fprintf(io_out, "Run at: %s\n", ctime(&my_time));
  fprintf(io_out, "Unit      Read(kB)    Write(kB)\n");
  fprintf(io_out, "-------------------------------\n");
  for(i=0; i < PSIO_MAXUNIT; i++) {
    total_read += psio_readlen[i];
    total_write += psio_writlen[i];

    if(psio_readlen[i] || psio_writlen[i])
    fprintf(io_out, "%3d   %10.1f   %10.1f\n",i,
        ((double) psio_readlen[i])/((double) 1024),
        ((double) psio_writlen[i])/((double) 1024));
  }
  fprintf(io_out, "-------------------------------\n");
  fprintf(io_out, "Total %10.1f   %10.1f\n",
      ((double) total_read)/((double) 1024),
      ((double) total_write)/((double) 1024));
  fclose(io_out);
  free(psio_readlen);
  free(psio_writlen);
#endif

  free(psio_unit);
  state_ = 0;
  files_keywords_.clear();
}

int psio_done(void) {
  if(_default_psio_lib_){
      // The old pointer implementation of this used to set the pointer to zero for
      // the test used in psio_init.  This is not necessary with smart pointers
      _default_psio_lib_.reset();
  }
  return true;
}

}

