/*!
** \file
** \brief Close input and output, stop input parser
** \ingroup
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libplugin/plugin.h>
namespace psi {

/*!
** psi_stop(): This function closes input and output files and
** deinitializes Input Parsing library.
**
** Arguments: none
**
** Returns: one of standard PSI error codes
** \ingroup CIOMR
*/

int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix)
{
  free(psi_file_prefix);
  fflush(outfile);
  fclose(outfile);
  fclose(infile);

  infile = NULL;
  outfile = NULL;
  psi_file_prefix = NULL;

  psi::yetiEnv.free();

  return(PSI_RETURN_SUCCESS);
}

}

