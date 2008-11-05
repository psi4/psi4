/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "dpd.h"

namespace psi {
	
int dpd_4mat_irrep_print(double **matrix, dpdparams4 *Params,
			 int block, int my_irrep, FILE *out)
{
  div_t fraction;
  int i,j,r_irrep;
  int rows, cols, cols_per_page, num_pages, last_page, page, first_col;

  cols_per_page = 5;

  r_irrep = block^my_irrep;

  rows = Params->rowtot[block];
  cols = Params->coltot[r_irrep];

  /* Determine the number of cols_per_page groups */
  fraction = div(cols,cols_per_page);
  num_pages = fraction.quot;  /* Number of complete column groups */
  last_page = fraction.rem;  /* Number of columns in last group */

  /* Loop over the complete column groups */
  for(page=0; page < num_pages; page++) {
      first_col = page*cols_per_page;

      fprintf(out,"\n           ");
      for(i=first_col; i < first_col+cols_per_page; i++) 
          fprintf(out,"              %5d",i);

      fprintf(out,"\n               ");
      for(i=first_col; i < first_col+cols_per_page; i++) 
          fprintf(out,"          (%3d,%3d)",
                  Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

      fprintf (out,"\n");
      for(i=0; i < rows; i++) {
          fprintf(out,"\n%5d  (%3d,%3d)",i,
                  Params->roworb[block][i][0], Params->roworb[block][i][1]);

          for(j=first_col; j < first_col+cols_per_page; j++)        
              fprintf (out,"%19.15f",matrix[i][j]);
        }

      fprintf (out,"\n");
    }

  /* Now print the remaining columns */
  if(last_page) {
      first_col = page*cols_per_page;

      fprintf(out,"\n           ");
      for(i=first_col; i < first_col+last_page; i++) 
	  fprintf(out,"              %5d",i);
      
      fprintf(out,"\n               ");
      for(i=first_col; i < first_col+last_page; i++) 
	  fprintf(out,"          (%3d,%3d)",
		  Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

      fprintf (out,"\n");
      for(i=0; i < rows; i++) {
	  fprintf(out,"\n%5d  (%3d,%3d)",i,
		  Params->roworb[block][i][0], Params->roworb[block][i][1]);

	  for(j=first_col; j < first_col+last_page; j++)
	      fprintf (out,"%19.15f",matrix[i][j]);
	}

      fprintf (out,"\n");
    }

  return 0;

}

} // namespace psi

