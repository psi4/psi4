/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** INTS.C
**
** Return values of one and two-electron integrals
**
** C. David Sherrill
** University of California, Berkeley
**
** Based on code from the DETCI program
** April 1998
*/

#include <cstdlib>
#include <cstdio>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

void read_integrals()
{
  int i, j, ij, k, l, kl, ijkl;
  int nbstri;
  double value;

  /* allocate memory for one and two electron integrals */
  nbstri = CalcInfo.nbstri;
  CalcInfo.onel_ints = init_array(nbstri);
  CalcInfo.twoel_ints = init_array(nbstri * (nbstri + 1) / 2);

  /* now read them in */

  if (Params.use_fzc_h) {
    if (Params.print_lvl > 3) 
      fprintf(outfile, "\n\tOne-electron integrals (frozen core operator):\n");
    iwl_rdone(Params.oei_file, PSIF_MO_FZC, CalcInfo.onel_ints, nbstri, 
              Params.oei_erase, (Params.print_lvl>3), outfile);
  }
  else {
    if (Params.print_lvl > 3) 
      fprintf(outfile, "\n\tOne-electron integrals (bare):\n");
    iwl_rdone(Params.oei_file, PSIF_MO_OEI, CalcInfo.onel_ints, nbstri, 
              Params.oei_erase, (Params.print_lvl>3), outfile);
  }

  if (Params.print_lvl > 4) 
    fprintf(outfile, "\n\tTwo-electron integrals:\n");

  iwl_rdtwo(Params.tei_file, CalcInfo.twoel_ints, ioff, 
     CalcInfo.nmo, Params.filter_ints ? CalcInfo.num_fzc_orbs : 0, 
     Params.filter_ints ? CalcInfo.num_fzv_orbs : 0, 
     (Params.print_lvl>4), outfile);

} 



double get_onel(int i, int j)
{
  int ij;

  ij = INDEX(i,j);
  return(CalcInfo.onel_ints[ij]);
}


double get_twoel(int i, int j, int k, int l)
{
  int ij, kl, ijkl;

  ij = INDEX(i,j);
  kl = INDEX(k,l);
  ijkl = INDEX(ij,kl);

  return(CalcInfo.twoel_ints[ijkl]);
}


/*
** get_mat_block()
**
** This function gets an irrep block of a full matrix
**
** C. David Sherrill
** May 1998
*/
void get_mat_block(double **src, double **dst, int dst_dim, int dst_offset,
                   int *dst2src)
{

  int P, Q, p, q;

  for (P=0; P<dst_dim; P++) {
    p = dst2src[P+dst_offset];
    for (Q=0; Q<dst_dim; Q++) {
      q = dst2src[Q+dst_offset];
      dst[P][Q] = src[p][q];
    } 
  }

}

}} // end namespace psi::detcas

