/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "dpd.h"

namespace psi {

int dpd_file2_mat_print(dpdfile2 *File, FILE *outfile)
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

      fprintf(outfile, "\n\tFile %3d DPD File2: %s\n", File->filenum,
	      File->label);
      fprintf(outfile,   "\tMatrix for Irrep %1d\n", h);
      fprintf(outfile,   "\t----------------------------------------\n");
      
      rows = Params->rowtot[h];
      cols = Params->coltot[h^my_irrep];

      /* Determine the number of cols_per_page groups */
      fraction = div(cols,cols_per_page);
      num_pages = fraction.quot;  /* Number of complete column groups */
      last_page = fraction.rem;  /* Number of columns in last group */

      /* Loop over the complete column groups */
      for(page=0; page < num_pages; page++) {
	  first_col = page*cols_per_page;

	  fprintf(outfile,"\n            ");
	  for(i=first_col; i < first_col+cols_per_page; i++) 
	      fprintf(outfile,"         %5d     ",i);

	  fprintf(outfile,"\n            ");
	  for(i=first_col; i < first_col+cols_per_page; i++) 
	      fprintf(outfile,"          (%3d)    ",
		      Params->colorb[h^my_irrep][i]);

	  fprintf (outfile,"\n");
	  for(i=0; i < rows; i++) {
	      fprintf(outfile,"\n%5d  (%3d)",i, Params->roworb[h][i]);

	      for(j=first_col; j < first_col+cols_per_page; j++)        
		  fprintf (outfile,"%19.15f",File->matrix[h][i][j]);
	    }

	  fprintf (outfile,"\n");
	}

      /* Now print the remaining columns */
      if(last_page) {
	  first_col = page*cols_per_page;
	  
	  fprintf(outfile,"\n            ");
	  for(i=first_col; i < first_col+last_page; i++) 
	      fprintf(outfile,"         %5d     ",i);
      
	  fprintf(outfile,"\n            ");
	  for(i=first_col; i < first_col+last_page; i++) 
	      fprintf(outfile,"          (%3d)    ",
		      Params->colorb[h^my_irrep][i]);
	  
	  fprintf (outfile,"\n");
	  for(i=0; i < rows; i++) {
	      fprintf(outfile,"\n%5d  (%3d)",i, Params->roworb[h][i]);
	  
	      for(j=first_col; j < first_col+last_page; j++)
		  fprintf (outfile,"%19.15f", File->matrix[h][i][j]);
	    }
	  fprintf (outfile,"\n");
	}
    }

  return 0;
}

}
