/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {
	
/* dpd_buf4_init(): Initializes a dpd four-index buffer for reading or writing
**   data.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf to be initialized.
**   int inputfile: The PSI unit number for the data on disk.
**   int pqnum: The index combination for the bra indices for the
**              data as it will be used in memory.
**   int rsnum: The index combination for the ket indices for the
**              data as it will be used in memory.
**   int file_pqnum: The index combination for the bra indices for the
**                   data as it will be stored on disk.
**   int file_rsnum: The index combination for the ket indices for the
**                   data as it will be stored on disk.
**   int anti: Boolean flag which indicates whether the data needs to
**             be antisymmetrized as it is read from disk.
**   char *label: The string labelling the PSIO TOC entry on disk.
*/

int dpd_buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum,
		 int file_pqnum, int file_rsnum, int anti, char *label)
{
  int h, nirreps, nump, nrows, p, Gp, Gr, offset;

  Buf->dpdnum = dpd_default;
  Buf->params = &(dpd_list[dpd_default].params4[pqnum][rsnum]);

  Buf->anti = anti;

  dpd_file4_init(&(Buf->file), inputfile, irrep, file_pqnum, file_rsnum, label);

  Buf->matrix = (double ***) malloc(Buf->params->nirreps*sizeof(double **));

  /* Set up shifted matrix info */
  nirreps = Buf->params->nirreps;
  Buf->shift.shift_type = 0;
  Buf->shift.rowtot = init_int_matrix(nirreps, nirreps);
  Buf->shift.coltot = init_int_matrix(nirreps, nirreps);
  Buf->shift.matrix = (double ****) malloc(nirreps * sizeof(double ***));

  /* row_offset lookup array */
  /* For a (pq,rs) buffer (assuming p and q are NOT packed), on which
     row of the irrep submatrix h does a given value of the orbital
     index, p, first appear?  This should work for non-totally
     symmetric quantities. */
  for(h=0,nump=0; h < nirreps; h++) nump += Buf->params->ppi[h];
  Buf->row_offset = init_int_matrix(nirreps, nump);
  for(h=0; h < nirreps; h++) {
    for(p=0; p < nump; p++) Buf->row_offset[h][p] = -1;
    for(Gp=0,nrows=0; Gp < nirreps; Gp++) {
      for(p=0; p < Buf->params->ppi[Gp]; p++) {
	if(Buf->params->qpi[Gp^h])
	  Buf->row_offset[h][Buf->params->poff[Gp]+p] = nrows;
	nrows += Buf->params->qpi[Gp^h];
      }
    }
  }

  /* col_offset lookup array */
  /* For a (pq,rs) buffer (assuming r and s are NOT packed), in which
  column of the irrep submatrix h does a given irrep of the orbital
  index, r, first appear?  This should work for non-totally-symmetric
  quantities.  */
  Buf->col_offset = init_int_matrix(nirreps, nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gr=0,offset=0; Gr < nirreps; Gr++) {
      Buf->col_offset[h][Gr] = offset;
      offset += Buf->params->rpi[Gr] * Buf->params->spi[Gr^h^Buf->file.my_irrep];
    }
  }

  return 0;
}

} // namespace psi

