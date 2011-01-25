/*!
** \file
** \brief Close input and output, stop input parser
** \ingroup
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libciomr/libciomr.h>

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
  fclose(outfile);
  fclose(infile);

  infile = NULL;
  outfile = NULL;
  psi_file_prefix = NULL;

  return(PSI_RETURN_SUCCESS);
}

}

