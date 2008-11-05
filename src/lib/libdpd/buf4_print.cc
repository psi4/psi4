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

int dpd_buf4_print(dpdbuf4 *Buf, FILE *out, int print_data)
{
  int h, i, all_buf_irrep;
  dpdparams4 *Params;

  all_buf_irrep = Buf->file.my_irrep;
  Params = Buf->params;

  fprintf(out, "\n\tDPD Buf4 for file4: %s\n", Buf->file.label);
  fprintf(out, "\n\tDPD Parameters:\n");
  fprintf(out,   "\t---------------\n");
  fprintf(out,   "\tpqnum = %d   rsnum = %d\n",
	  Params->pqnum, Params->rsnum);
  fprintf(out, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(out, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
    fprintf(out,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	    Params->rowtot[i], Params->coltot[i^all_buf_irrep]);
  fflush(out);

  if(print_data) {
    for(h=0; h < Buf->params->nirreps; h++) {
      fprintf(out, "\n\tFile %3d DPD Buf4: %s\n", Buf->file.filenum,
	      Buf->file.label);
      fprintf(out,   "\tMatrix for Irrep %1d\n", h);
      fprintf(out,   "\t----------------------------------------\n");
      dpd_buf4_mat_irrep_init(Buf, h);
      dpd_buf4_mat_irrep_rd(Buf, h);
      dpd_4mat_irrep_print(Buf->matrix[h], Buf->params, h, all_buf_irrep, outfile);
      dpd_buf4_mat_irrep_close(Buf, h);
    }
  }

  return 0;

}

} // namespace psi
