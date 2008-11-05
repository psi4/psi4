/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "dpd.h"

namespace psi {

int dpd_file2_mat_print(dpdfile2 *File, FILE *out)
{
  div_t fraction;
  int i,j;
  int rows, cols, cols_per_page, num_pages, last_page, page, first_col;
  int h, my_irrep;
  dpdparams2 *Params;

  Params = File->params;
  my_irrep = File->my_irrep;

  cols_per_page = 9;

  for(h=0; h < File->params->nirreps; h++) {

      fprintf(out, "\n\tFile %3d DPD File2: %s\n", File->filenum,
	      File->label);
      fprintf(out,   "\tMatrix for Irrep %1d\n", h);
      fprintf(out,   "\t----------------------------------------\n");
      
      rows = Params->rowtot[h];
      cols = Params->coltot[h^my_irrep];

      /* Determine the number of cols_per_page groups */
      fraction = div(cols,cols_per_page);
      num_pages = fraction.quot;  /* Number of complete column groups */
      last_page = fraction.rem;  /* Number of columns in last group */

      /* Loop over the complete column groups */
      for(page=0; page < num_pages; page++) {
	  first_col = page*cols_per_page;

	  fprintf(out,"\n            ");
	  for(i=first_col; i < first_col+cols_per_page; i++) 
	      fprintf(out,"         %5d     ",i);

	  fprintf(out,"\n            ");
	  for(i=first_col; i < first_col+cols_per_page; i++) 
	      fprintf(out,"          (%3d)    ",
		      Params->colorb[h^my_irrep][i]);

	  fprintf (out,"\n");
	  for(i=0; i < rows; i++) {
	      fprintf(out,"\n%5d  (%3d)",i, Params->roworb[h][i]);

	      for(j=first_col; j < first_col+cols_per_page; j++)        
		  fprintf (out,"%19.15f",File->matrix[h][i][j]);
	    }

	  fprintf (out,"\n");
	}

      /* Now print the remaining columns */
      if(last_page) {
	  first_col = page*cols_per_page;
	  
	  fprintf(out,"\n            ");
	  for(i=first_col; i < first_col+last_page; i++) 
	      fprintf(out,"         %5d     ",i);
      
	  fprintf(out,"\n            ");
	  for(i=first_col; i < first_col+last_page; i++) 
	      fprintf(out,"          (%3d)    ",
		      Params->colorb[h^my_irrep][i]);
	  
	  fprintf (out,"\n");
	  for(i=0; i < rows; i++) {
	      fprintf(out,"\n%5d  (%3d)",i, Params->roworb[h][i]);
	  
	      for(j=first_col; j < first_col+last_page; j++)
		  fprintf (out,"%19.15f", File->matrix[h][i][j]);
	    }
	  fprintf (out,"\n");
	}
    }

  return 0;
}

} // namespace psi
