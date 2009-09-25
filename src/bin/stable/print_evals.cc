/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void print_evals(double **evals, int *rank)
{
  int h, i;

  fprintf(outfile, "\t  #  ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "    %3s  ",moinfo.labels[h]);
  fprintf(outfile, "\n");
  fprintf(outfile, "\t---- ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "---------");
  fprintf(outfile, "\n");

  for(i=0; i < 5; i++) {
    fprintf(outfile, "\t %2d  ", i);
      for(h=0; h < moinfo.nirreps; h++) {
	  if(rank[h] <= i) fprintf(outfile, "         ");
	  else fprintf(outfile, " %7.4f ", evals[h][i]);
	}
      fprintf(outfile, "\n");
    }

  fprintf(outfile, "\n");
}

}} // namespace psi::stable
