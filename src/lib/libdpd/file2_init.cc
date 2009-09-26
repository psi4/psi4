/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {

/* dpd_file2_init(): Initializes a dpd two-index file for reading
** or writing data.
**
** Arguments:
**   dpdfile2 *File: A pointer to the two-index dpdfile.
**   int filenum: The PSI unit number for this file.
**   int irrep: The symmetry of the data (=0 for totally-symmetric)
**   int pnum: The orbital subspace number for the left index [see
**             dpd_init()].
**   int qnum: The orbital subspace number for the right index [see
**             dpd_init()].
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_file2_init(dpdfile2 *File, int filenum, int irrep, int pnum,
		   int qnum, const char *label)
{
  int i, q, rs, nirreps;
  struct dpd_file2_cache_entry *this_entry;
  int *coltot, *colidx, **colorb, *qpi, *qoff, *qsym;

  File->dpdnum = dpd_default;
  File->params = &(dpd_list[dpd_default].params2[pnum][qnum]);
  strcpy(File->label,label);
  File->filenum = filenum;
  File->my_irrep = irrep;

  nirreps = File->params->nirreps;

  this_entry = dpd_file2_cache_scan(filenum, irrep, pnum, qnum, label, dpd_default);
  if(this_entry != NULL) {
      File->incore = 1;
      File->matrix = this_entry->matrix;
    }
  else {
      File->incore = 0;
      File->matrix = (double ***) malloc(File->params->nirreps*sizeof(double **));
    }

  /* Construct logical subfile pointers */
  File->lfiles = (psio_address *) malloc(File->params->nirreps *
					 sizeof(psio_address));
  File->lfiles[0] = PSIO_ZERO;
  for(i=1; i < File->params->nirreps; i++)
    File->lfiles[i] = psio_get_address(File->lfiles[i-1],
				       (File->params->rowtot[i-1] *
					File->params->coltot[(i-1)^irrep] *
					sizeof(double)));

  /* Force all two-index files into cache */
  /*  dpd_file2_cache_add(File); */

  return 0;
}

}
