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
  int errcod;

  errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    abort();
  ip_cwk_add(":OEPROP");
  tstart(outfile);
  psio_init(); psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);

  return;
}

void stop_io()
{
  tstop(outfile);
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);

  return;
}


void punt(const char *errmsg)
{
  fprintf(stderr, "Error: %s\n", errmsg);
  stop_io();
  abort();
}

}} // namespace psi::oeprop
