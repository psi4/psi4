#include <cstdio>
#include <cstdlib>

/*!
** \file
** \brief Print a matrix to a file in a formatted style
** \ingroup QT
*/

namespace psi {
	
/*!
** mat_print(): Prints a matrix to a file in a formatted style
**
** \param matrix  = matrix to print
** \param rows    = number of rows
** \param cols    = number of columns
** \param outfile = output file pointer for printing
**
** Returns: Always returns zero...
**
** \ingroup QT
*/
int mat_print(double **matrix, int rows, int cols, FILE *outfile)
{
  div_t fraction;
  int i,j;
  int cols_per_page, num_pages, last_page, page, first_col;

  cols_per_page = 5;

  /* Determine the number of cols_per_page groups */
  fraction = div(cols,cols_per_page);
  num_pages = fraction.quot;  /* Number of complete column groups */
  last_page = fraction.rem;  /* Number of columns in last group */

  /* Loop over the complete column groups */
  for(page=0; page < num_pages; page++) {
      first_col = page*cols_per_page;

      fprintf(outfile,"\n      ");
      for(i=first_col; i < first_col+cols_per_page; i++) 
          fprintf(outfile,"      %5d        ",i);

      fprintf (outfile,"\n");
      for(i=0; i < rows; i++) {
          fprintf(outfile,"\n%5d ",i);

          for(j=first_col; j < first_col+cols_per_page; j++)        
              fprintf (outfile,"%19.15f",matrix[i][j]);
        }

      fprintf (outfile,"\n");
    }

  /* Now print the remaining columns */
  if(last_page) {
      first_col = page*cols_per_page;

      fprintf(outfile,"\n      ");
      for(i=first_col; i < first_col+last_page; i++) 
          fprintf(outfile,"      %5d        ",i);
      
      fprintf (outfile,"\n");
      for(i=0; i < rows; i++) {
	  fprintf(outfile,"\n%5d ",i);

	  for(j=first_col; j < first_col+last_page; j++)
	      fprintf (outfile,"%19.15f",matrix[i][j]);
	}

      fprintf (outfile,"\n");
    }

  return 0;

}

}

