/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/* dpd_buf4_print(): Prints out data for all irreps of a dpd
** four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf to be printed.
**   FILE *outfile: The formatted output file stream.
**   int print_data: 0 = print buf4 parameters only; 1 = print matrices
*/

int dpd_buf4_print(dpdbuf4 *Buf, FILE *outfile, int print_data)
{
  int h, i, all_buf_irrep;
  dpdparams4 *Params;

  all_buf_irrep = Buf->file.my_irrep;
  Params = Buf->params;

  fprintf(outfile, "\n\tDPD Buf4 for file4: %s\n", Buf->file.label);
  fprintf(outfile, "\n\tDPD Parameters:\n");
  fprintf(outfile,   "\t---------------\n");
  fprintf(outfile,   "\tpqnum = %d   rsnum = %d\n",
	  Params->pqnum, Params->rsnum);
  fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(outfile, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
    fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	    Params->rowtot[i], Params->coltot[i^all_buf_irrep]);
  fflush(outfile);

  if(print_data) {
    for(h=0; h < Buf->params->nirreps; h++) {
      fprintf(outfile, "\n\tFile %3d DPD Buf4: %s\n", Buf->file.filenum,
	      Buf->file.label);
      fprintf(outfile,   "\tMatrix for Irrep %1d\n", h);
      fprintf(outfile,   "\t----------------------------------------\n");
      dpd_buf4_mat_irrep_init(Buf, h);
      dpd_buf4_mat_irrep_rd(Buf, h);
      dpd_4mat_irrep_print(Buf->matrix[h], Buf->params, h, all_buf_irrep, outfile);
      dpd_buf4_mat_irrep_close(Buf, h);
    }
  }

  return 0;

}

}
