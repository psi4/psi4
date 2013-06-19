/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_opdm(void);
void uhf_opdm(void);
void rhf_sf_opdm(void);
void uhf_sf_opdm(void);

void opdm(void)
{
  if(params.gradient) {
    if(params.ref == 0) rhf_opdm();
    else if(params.ref == 2) uhf_sf_opdm();
  }
  else {
    if(params.ref == 0) rhf_opdm();
    else if(params.ref == 2) uhf_opdm();
  }
}

void rhf_opdm(void)
{
  int h, i, j, a, b;
  int I, J, A, B;
  double trace_IJ=0.0;
  double trace_AB=0.0;
  double alpha;
  dpdfile2 DIJ, DAB, D;
  dpdbuf4 T2, T2A, T2B;

  if(params.gradient) alpha = 1.0;
  else alpha = 2.0;

  dpd_->file2_init(&DIJ, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_->contract442(&T2A, &T2B, &DIJ, 0, 0, -alpha, 0);
  dpd_->buf4_close(&T2B);
  dpd_->buf4_close(&T2A);
  trace_IJ = dpd_->file2_trace(&DIJ);
  if(params.gradient) dpd_->file2_copy(&DIJ, PSIF_CC_OEI, "Dij");
  dpd_->file2_close(&DIJ);

  dpd_->file2_init(&DAB, PSIF_CC_OEI, 0, 1, 1, "DAB");
  dpd_->file2_scm(&DAB,0.0);
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_->contract442(&T2A, &T2B, &DAB, 3, 3, alpha, 0);
  dpd_->buf4_close(&T2B);
  dpd_->buf4_close(&T2A);
  trace_AB = dpd_->file2_trace(&DAB);
  if(params.gradient) dpd_->file2_copy(&DAB, PSIF_CC_OEI, "Dab");
  dpd_->file2_close(&DAB);

  if(params.gradient) {
    /* zero out Dia DIA Dai DAI */
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
    dpd_->file2_scm(&D,0.0);
    dpd_->file2_close(&D);
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
    dpd_->file2_scm(&D,0.0);
    dpd_->file2_close(&D);
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
    dpd_->file2_scm(&D,0.0);
    dpd_->file2_close(&D);
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
    dpd_->file2_scm(&D,0.0);
  }
  
/*
  fprintf(outfile, "\n\tTrace of OPDM(2)_IJ     = %20.15f\n", trace_IJ);
  fprintf(outfile, "\n\tTrace of OPDM(2)_AB     = %20.15f\n", trace_AB);
*/
  fprintf(outfile, "\tTrace of OPDM(2)        = %20.15f\n", trace_IJ+trace_AB);
  fflush(outfile);

  return;
}

void uhf_opdm(void)
{
  int h, i, j, a, b;
  int I, J, A, B;
  double traceA=0.0;
  double traceB=0.0;
  dpdfile2 D;
  dpdfile2 T1A, T1B;
  dpdbuf4 T2A, T2B;

  if(params.semicanonical) {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
    dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&T1A, &T1B, &D, 0, 0, -1, 0);
    dpd_->file2_close(&T1A);
    dpd_->file2_close(&T1B);
    
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_->contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  else {	
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_->contract442(&T2A, &T2B, &D, 0, 0, -1, 0);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
  traceA += dpd_->file2_trace(&D);
  dpd_->buf4_close(&T2A);
  dpd_->buf4_close(&T2B);
  dpd_->file2_close(&D);
    
  if(params.semicanonical) {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, "Dij");
    dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&T1A, &T1B, &D, 0, 0, -1, 0);
    dpd_->file2_close(&T1A);
    dpd_->file2_close(&T1B);
    
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_->contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  else { 
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, "Dij");
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_->contract442(&T2A, &T2B, &D, 0, 0, -1, 0);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->contract442(&T2A, &T2B, &D, 1, 1, -1, 1);
  traceB += dpd_->file2_trace(&D);
  dpd_->buf4_close(&T2A);
  dpd_->buf4_close(&T2B);
  dpd_->file2_close(&D);
  
  if(params.semicanonical) {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
    dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&T1A, &T1B, &D, 1, 1, 1, 0);
    dpd_->file2_close(&T1A);
    dpd_->file2_close(&T1B);
    
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_->contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  else {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_->contract442(&T2A, &T2B, &D, 3, 3, 1, 0);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->contract442(&T2A, &T2B, &D, 2, 2, 1, 1);
  traceA += dpd_->file2_trace(&D);
  dpd_->buf4_close(&T2A);
  dpd_->buf4_close(&T2B);
  dpd_->file2_close(&D);
    
  if(params.semicanonical) {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, "Dab");
    dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&T1A, &T1B, &D, 1, 1, 1, 0);
    dpd_->file2_close(&T1A);
    dpd_->file2_close(&T1B);
    
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_->contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  else {
    dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, "Dab");
    dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_->contract442(&T2A, &T2B, &D, 3, 3, 1, 0);
    dpd_->buf4_close(&T2A);
    dpd_->buf4_close(&T2B);
  }
  
  dpd_->buf4_init(&T2A, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->buf4_init(&T2B, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
  traceB += dpd_->file2_trace(&D);
  dpd_->buf4_close(&T2A);
  dpd_->buf4_close(&T2B);
  dpd_->file2_close(&D);
   
  /*
  fprintf(outfile,"\n");
  fprintf(outfile,"\tTrace of Alpha OPDM(2)  = %20.15f\n", fabs(traceA));
  fprintf(outfile,"\tTrace of Beta OPDM(2)   = %20.15f\n", fabs(traceB));
  */
    
  return;
}


/* This code isn't used anymore, and if you try to resurrect it, be sure to
note that the contract244 calls that point to the same file4 are not safe
in the dpd_cache structure.  -TDC, 11/21/08
*/
void rhf_sf_opdm(void) 
{
  dpdfile2 D;
  dpdbuf4 TA, TB;
  double trace_IJ = 0.0; 
  double trace_ij = 0.0;
  double trace_AB = 0.0; 
  double trace_ab = 0.0;

  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_->buf4_copy(&TA, PSIF_CC_TMP0, "tIJAB");
  dpd_->buf4_close(&TA);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  dpd_->buf4_copy(&TA, PSIF_CC_TMP0, "tijab");
  dpd_->buf4_close(&TA);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->buf4_copy(&TA, PSIF_CC_TMP0, "tIjAb");
  dpd_->buf4_close(&TA);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_->buf4_copy(&TA, PSIF_CC_TMP0, "tiJaB");
  dpd_->buf4_close(&TA);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_->file2_scm(&D, 0.0);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_->buf4_init(&TB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_->contract442(&TA, &TB, &D, 0, 0, -1.0, 0.0);
  dpd_->buf4_close(&TA);
  dpd_->buf4_close(&TB);

  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->buf4_init(&TB, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->contract442(&TA, &TB, &D, 0, 0, -1.0, 1.0);
  dpd_->buf4_close(&TA);
  dpd_->buf4_close(&TB);
  trace_IJ += dpd_->file2_trace(&D);
  dpd_->file2_close(&D);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  dpd_->file2_scm(&D, 0.0);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  dpd_->buf4_init(&TB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  dpd_->contract442(&TA, &TB, &D, 0, 0, -1.0, 0.0);
  dpd_->buf4_close(&TB);
  dpd_->buf4_close(&TA);

  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_->buf4_init(&TB, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_->contract442(&TA, &TB, &D, 0, 0, -1.0, 1.0);
  dpd_->buf4_close(&TB);
  dpd_->buf4_close(&TA);
  trace_ij += dpd_->file2_trace(&D);
  dpd_->file2_close(&D);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  dpd_->file2_scm(&D, 0.0);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_->buf4_init(&TB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_->contract442(&TA, &TB, &D, 3, 3, 1.0, 0.0);
  dpd_->buf4_close(&TA);
  dpd_->buf4_close(&TB);

  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_->buf4_init(&TB, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_->contract442(&TA, &TB, &D, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&TA);
  dpd_->buf4_close(&TB);
  trace_AB += dpd_->file2_trace(&D);
  dpd_->file2_close(&D);

  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  dpd_->file2_scm(&D, 0.0);
  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_->buf4_init(&TB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_->contract442(&TA, &TB, &D, 3, 3, 1.0, 0.0);
  dpd_->buf4_close(&TB);
  dpd_->buf4_close(&TA);

  dpd_->buf4_init(&TA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->buf4_init(&TB, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_->contract442(&TA, &TB, &D, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&TB);
  dpd_->buf4_close(&TA);
  trace_ab += dpd_->file2_trace(&D);
  dpd_->file2_close(&D);

  /* zero out Dia DIA Dai DAI */
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  dpd_->file2_scm(&D,0.0);
  dpd_->file2_close(&D);
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  dpd_->file2_scm(&D,0.0);
  dpd_->file2_close(&D);
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  dpd_->file2_scm(&D,0.0);
  dpd_->file2_close(&D);
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  dpd_->file2_scm(&D,0.0);
  dpd_->file2_close(&D);
  

  fprintf(outfile, "\n\tTrace of IJ onepdm = %20.15f\n", trace_IJ);
  fprintf(outfile, "\tTrace of ij onepdm = %20.15f\n", trace_ij);
  fprintf(outfile, "\tTrace of oo onepdm = %20.15f\n", trace_IJ+trace_ij);
  fprintf(outfile, "\tTrace of AB onepdm = %20.15f\n", trace_AB);
  fprintf(outfile, "\tTrace of ab onepdm = %20.15f\n", trace_ab);
  fprintf(outfile, "\tTrace of vv onepdm = %20.15f\n", (trace_AB+trace_ab));
  fprintf(outfile, "\tTrace of total onepdm = %20.15f\n", trace_IJ+trace_ij+trace_AB+trace_ab);
}

void uhf_sf_opdm(void)
{
  fprintf(outfile,"\n\tNot yet yo! -MLA\n");
  exit(PSI_RETURN_FAILURE);
}

}} /* End namespaces */
