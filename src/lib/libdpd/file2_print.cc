/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include "dpd.h"

extern "C" {

/* dpd_file2_print(): Prints out data for all irreps of a two-index dpdfile.
**
** Arguments:
**   struct dpdfile2 *File: A pointer to the dpdfile to be printed.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_file2_print(dpdfile2 *File, FILE *outfile)
{
  int i, my_irrep;
  dpdparams2 *Params;

  my_irrep = File->my_irrep;
  Params = File->params;

  fprintf(outfile, "\n\tDPD File2: %s\n", File->label);
  fprintf(outfile,   "\tDPD Parameters:\n");
  fprintf(outfile,   "\t------------------\n");
  fprintf(outfile,   "\tpnum = %d   qnum = %d   irrep = %d \n",
      Params->pnum, Params->qnum, File->my_irrep);
  fprintf(outfile,   "\tIrreps = %1d\n\n", Params->nirreps);
  fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(outfile, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
      fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	      Params->rowtot[i], Params->coltot[i^my_irrep]);
  fflush(outfile);

  dpd_file2_mat_init(File);
  dpd_file2_mat_rd(File);
  dpd_file2_mat_print(File, outfile);
  dpd_file2_mat_close(File);

  return 0;

}

} /* extern "C" */
