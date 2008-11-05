/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

int dpd_file4_print(dpdfile4 *File, FILE *out)
{
  int i, h, my_irrep;
  dpdparams4 *Params;

  my_irrep = File->my_irrep;
  Params = File->params;

  fprintf(out, "\n\tDPD File4: %s\n", File->label);
  fprintf(out, "\n\tDPD Parameters:\n");
  fprintf(out,   "\t---------------\n");
  fprintf(out,   "\tpqnum = %d   rsnum = %d\n",
	  Params->pqnum, Params->rsnum);
  fprintf(out, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(out, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
      fprintf(out,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	      Params->rowtot[i], Params->coltot[i^my_irrep]);
  fflush(outfile);

  for(h=0; h < File->params->nirreps; h++) {
      fprintf(out, "\n\tFile %3d DPD File4: %s\n", File->filenum,
	      File->label);
      fprintf(out,   "\tMatrix for Irrep %1d\n", h);
      fprintf(out,   "\t----------------------------------------\n");
      dpd_file4_mat_irrep_init(File, h);
      dpd_file4_mat_irrep_rd(File, h);
      dpd_4mat_irrep_print(File->matrix[h], File->params, h, my_irrep, outfile);
      dpd_file4_mat_irrep_close(File, h);
    }

  return 0;

}

} // namespace psi
