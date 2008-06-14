/*! \defgroup PSICLEAN psiclean: Delete scratch files */

/*! 
** \file
** \ingroup PSICLEAN
** \brief Delete scratch files
**
** Utility program to delete scratch files.  Generalization of earlier
** PSI2.0 shell script which was limited to scratch files being put
** in /tmp[0-9]/$user/$name.* .  Here we will search the default path
** instead.
**
** C. David Sherrill
** 
*/ 

#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psi4-dec.h>

#define MAX_STRING 300

namespace psi { namespace psiclean {

int psiclean(int argc, char *argv[])
{
  ULI i, nvol;
  int errcod;
  char *vpath;
  char *basename;
  char fileslist[MAX_STRING];
  char cmdstring[MAX_STRING];

  // psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  
  /* Initialize the I/O system */
  psio_init(); psio_ipv1_config();

  /* Get the number of volumes */
  nvol = psio_get_numvols_default();

  errcod = psio_get_filename_default(&basename);

  for (i=0; i<nvol; i++) {
    errcod = psio_get_volpath_default(i, &vpath);

    sprintf(fileslist,"%s%s.*",vpath,basename);
    sprintf(cmdstring,"echo Removing files %s%s",vpath,basename);
    system(cmdstring);
    sprintf(cmdstring,"ls -l %s",fileslist);
    system(cmdstring);
    sprintf(cmdstring,"/bin/rm %s",fileslist);
    system(cmdstring);
    free(vpath);
  }
  free(basename);

  /* we're done, clean up */
  //psio_done();
  //psi_stop(infile,outfile,psi_file_prefix);

  return PSI_RETURN_SUCCESS;
}

}}

