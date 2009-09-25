/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** SETUP_IO.C
**
** Contains some setup/shutdown I/O stuff
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "globals.h"

namespace psi { namespace detcas {

/*
** init_io(): Function opens input and output files
*/
void init_io(int argc, char *argv[])
{
  int i;
  int num_extra_args;
  char **extra_args;

  extra_args = (char **) malloc(argc*sizeof(char *)); 
  for (i=1,num_extra_args=0; i<argc; i++) {
    if (strcmp(argv[i], "--quiet") == 0) {
      Params.print_lvl = 0;
    }
    else {
      extra_args[num_extra_args++] = argv[i];
    }
  }
 
  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args, extra_args, 0); 

  /*
  init_in_out(argc-parsed,argv+parsed);
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  */

  if (Params.print_lvl) tstart(outfile);
  ip_cwk_add(":DETCAS");
  psio_init(); psio_ipv1_config();

  free(extra_args);
}



/*
** close_io(): Function closes down I/O and exits
*/
void close_io(void)
{
   psio_done();
   if (Params.print_lvl) tstop(outfile);
   psi_stop(infile,outfile,psi_file_prefix);
}



/*
** check() acts as an abort function if the condition 'a' is not true;
**   it shuts down i/o and returns an error
*/
void check(int a, const char *errmsg)
{
  if (!a) {
    fprintf(outfile, "%s\n", errmsg);
    close_io();
    exit(1);
  }
}

}} // end namespace psi::detcas


/*
** The gprgid() function required by the PSI I/O libraries
*/
extern "C" {
  const char *gprgid()
  {
    const char *prgid = "DETCAS";
    return(prgid);
  }
}


