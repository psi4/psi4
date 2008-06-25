/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void start_io(int argc, char *argv[])
{
  //  abort();
  //ip_cwk_add(":INPUT");
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_cwk_add(":INPUT");
  
  tstart(outfile);
  psio_init();
  psio_ipv1_config();

  //free(extra_args);

  return;
}

void stop_io()
{
  tstop(outfile);
  psio_done();
//  psi_stop(infile,outfile,psi_file_prefix);
}

}} // namespace psi::input

