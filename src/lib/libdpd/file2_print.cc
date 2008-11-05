/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/* dpd_file2_print(): Prints out data for all irreps of a two-index dpdfile.
**
** Arguments:
**   struct dpdfile2 *File: A pointer to the dpdfile to be printed.
**   FILE *out: The formatted output file stream.
*/

int dpd_file2_print(dpdfile2 *File, FILE *out)
{
  int i, my_irrep;
  dpdparams2 *Params;

  my_irrep = File->my_irrep;
  Params = File->params;

  fprintf(out, "\n\tDPD File2: %s\n", File->label);
  fprintf(out,   "\tDPD Parameters:\n");
  fprintf(out,   "\t------------------\n");
  fprintf(out,   "\tpnum = %d   qnum = %d   irrep = %d \n",
      Params->pnum, Params->qnum, File->my_irrep);
  fprintf(out,   "\tIrreps = %1d\n\n", Params->nirreps);
  fprintf(out, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(out, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
      fprintf(out,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	      Params->rowtot[i], Params->coltot[i^my_irrep]);
  fflush(out);

  dpd_file2_mat_init(File);
  dpd_file2_mat_rd(File);
  dpd_file2_mat_print(File, out);
  dpd_file2_mat_close(File);

  return 0;

}

} // namespace psi
