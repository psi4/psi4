/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void read_basset_info()
{
  int i;

  exps  = chkpt_rd_exps();
  contr = chkpt_rd_contr();
  sprim = chkpt_rd_sprim();
  snuc  = chkpt_rd_snuc();
  stype = chkpt_rd_stype();
  snumg = chkpt_rd_snumg();
  sloc  = chkpt_rd_sloc();
  sloc_new  = chkpt_rd_sloc_new();
  
  lmax = 0;
  for(i=0;i<nshell;i++)
    lmax = (lmax >= stype[i]) ? lmax : stype[i];
  lmax--;

  if (print_lvl >= PRINTBASSETLEVEL) {
    fprintf(outfile,"    Prim#  \t   Exponent \tContr. coeff.\n");
    for(i=0;i<nprim;i++)
      fprintf(outfile,"  %5d  \t%12.5lf\t%13.11lf\n",i+1,exps[i],contr[i]);
    fprintf(outfile,"\n\n  Shell# Nuc#\t  L\tSPRIM\t SLOC\tSNUMG\n");
    for(i=0;i<nshell;i++)
      fprintf(outfile,"  %3d  \t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\n",i+1,
              snuc[i],stype[i]-1,sprim[i],sloc[i],snumg[i]);
    fprintf(outfile,"\n\n");
  }

}

}} // namespace psi::oeprop
