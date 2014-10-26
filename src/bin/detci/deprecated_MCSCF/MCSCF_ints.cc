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
#include "psi4-dec.h"

namespace psi { namespace detci {

void mcscf_read_integrals()
{
  int i, j, ij, k, l, kl, ijkl;
  int nbstri;
  double value;

  /* allocate memory for one and two electron integrals */
  nbstri = CalcInfo.nmotri;
  CalcInfo.onel_ints = init_array(nbstri);
  CalcInfo.onel_ints_bare = init_array(nbstri);
  CalcInfo.twoel_ints = init_array(nbstri * (nbstri + 1) / 2);

  /* now read them in */

  if (MCSCF_Parameters.use_fzc_h) {
    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf("\n\tOne-electron integrals (frozen core operator):\n");
    }
    // can't erase file, re-reading it below
    iwl_rdone(MCSCF_Parameters.oei_file, PSIF_MO_FZC, CalcInfo.onel_ints, nbstri, 
              0, (MCSCF_Parameters.print_lvl>3), "outfile");
  }
  else {
    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf("\n\tOne-electron integrals (bare):\n");
    }
    // can't erase file, re-reading it below
    iwl_rdone(MCSCF_Parameters.oei_file, PSIF_MO_OEI, CalcInfo.onel_ints, nbstri, 
              0, (MCSCF_Parameters.print_lvl>3), "outfile");
  }

  /* even if we utilize frozen core operator for some terms, the
     current Lagrangian code is forced to use bare h, so let's grab
     that, too (both should be available from the transformation code)
     -CDS 10/3/14
  */
  if (MCSCF_Parameters.print_lvl > 3) { 
    outfile->Printf("\n\tOne-electron integrals (bare):\n");
    }
  iwl_rdone(MCSCF_Parameters.oei_file, PSIF_MO_OEI, CalcInfo.onel_ints_bare, nbstri,
            MCSCF_Parameters.oei_erase ? 0 : 1, (MCSCF_Parameters.print_lvl>3), "outfile");

  
  if (MCSCF_Parameters.print_lvl > 6) 
    outfile->Printf("\n\tTwo-electron integrals:\n");

  iwl_rdtwo(MCSCF_Parameters.tei_file, CalcInfo.twoel_ints, ioff, 
     CalcInfo.nmo, MCSCF_Parameters.filter_ints ? CalcInfo.num_fzc_orbs : 0, 
     MCSCF_Parameters.filter_ints ? CalcInfo.num_fzv_orbs : 0, 
     (MCSCF_Parameters.print_lvl>6), "outfile");

} 



double mcscf_get_onel(int i, int j)
{
  int ij;

  ij = INDEX(i,j);
  return(CalcInfo.onel_ints[ij]);
}


double mcscf_get_twoel(int i, int j, int k, int l)
{
  int ij, kl, ijkl;

  ij = INDEX(i,j);
  kl = INDEX(k,l);
  ijkl = INDEX(ij,kl);

  return(CalcInfo.twoel_ints[ijkl]);
}


/*
** mcscf_get_mat_block()
**
** This function gets an irrep block of a full matrix
**
** C. David Sherrill
** May 1998
*/
void mcscf_get_mat_block(double **src, double **dst, int dst_dim, int dst_offset,
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

}} // end namespace psi::detci

