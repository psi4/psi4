/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

void start_io(int argc, char *argv[])
{

  tstart();
  chkpt_init(PSIO_OPEN_OLD);

  return;
}

void stop_io()
{
  tstop();

  return;
}


//void punt(const char *errmsg)
//{
//  fprintf(stderr, "Error: %s\n", errmsg);
//  stop_io();
//  abort();
//}

}} // namespace psi::oeprop
