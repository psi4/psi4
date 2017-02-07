/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
**  \ingroup DETCI
**  \brief Some code associated with the H0block structure, which does
**    the preconditioning for the Olsen/Davidson procedure.
**
** David Sherrill
** July 1995
**
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/ciwave.h"
#include "psi4/detci/slaterd.h"

namespace psi { namespace detci {

#define SMALL_DET 1e-10


/*
** H0block_init()
**
** initialize everything but buf_num and buf_member, which depend on the
** CIvector structure
*/
void CIWavefunction::H0block_init(unsigned int size) {

   unsigned int size2;

   if (size > Parameters_->h0blocksize) H0block_->size = Parameters_->h0blocksize;
   else H0block_->size = size;
   H0block_->coupling_size = Parameters_->h0block_coupling_size;

   if (H0block_->coupling_size)
     size2 = H0block_->size + H0block_->coupling_size;
   else size2 = H0block_->size;

   if (print_ > 1)
     outfile->Printf("    Total H0block size (including coupling): %d\n",size2);

   H0block_->osize = H0block_->size;
   H0block_->guess_size = Parameters_->h0guess_size;
   H0block_->oguess_size = H0block_->guess_size;
   H0block_->ocoupling_size = H0block_->coupling_size;

   if (H0block_->size) {
      H0block_->H0b = init_matrix(H0block_->size, H0block_->size);
      if (Parameters_->precon == PRECON_GEN_DAVIDSON)
        H0block_->H0b_diag_transpose = init_array(H0block_->size);
      H0block_->H0b_diag = init_matrix(H0block_->size, H0block_->size);
      H0block_->H0b_eigvals = init_array(H0block_->size);
      H0block_->tmp1 = init_matrix(H0block_->size, H0block_->size);
      H0block_->H00 = init_array(size2);
      H0block_->c0b = init_array(size2);
      H0block_->c0bp = init_array(size2);
      H0block_->s0b = init_array(size2);
      H0block_->s0bp = init_array(size2);
      H0block_->alplist = init_int_array(size2);
      H0block_->betlist = init_int_array(size2);
      H0block_->alpidx = init_int_array(size2);
      H0block_->betidx = init_int_array(size2);
      H0block_->blknum = init_int_array(size2);
      H0block_->pair = init_int_array(size2);
      if (Parameters_->precon == PRECON_H0BLOCK_INVERT)
        H0block_->H0b_inv = init_matrix(H0block_->size, H0block_->size);

      if (Parameters_->h0block_coupling) {
        H0block_->tmp_array1 = init_array(size2);
        H0block_->tmp_array2 = init_array(size2);
        }
     }
}


void CIWavefunction::H0block_free(void)
{
   int i;

   if (H0block_->osize) {
      free_matrix(H0block_->H0b, H0block_->osize);
      if (Parameters_->precon == PRECON_GEN_DAVIDSON)
         free(H0block_->H0b_diag_transpose);
      free_matrix(H0block_->H0b_diag, H0block_->osize);
      free_matrix(H0block_->tmp1, H0block_->osize);
      free(H0block_->H00);
      free(H0block_->c0b);
      free(H0block_->c0bp);
      free(H0block_->s0b);
      free(H0block_->s0bp);
      free(H0block_->alplist);
      free(H0block_->betlist);
      free(H0block_->alpidx);
      free(H0block_->betidx);
      free(H0block_->blknum);
      free(H0block_->pair);
      if (Parameters_->precon == PRECON_H0BLOCK_INVERT)
         free_matrix(H0block_->H0b_inv, H0block_->osize);
      if (Parameters_->h0block_coupling) {
         free(H0block_->tmp_array1);
         free(H0block_->tmp_array2);
         }
      if (H0block_->nbuf) {
         free(H0block_->buf_num);
         for (i=0; i<H0block_->nbuf; i++) free(H0block_->buf_member[i]);
         free(H0block_->buf_member);
         }
      }
}

void CIWavefunction::H0block_print(void)
{
   int i;

   outfile->Printf( "\nMembers of H0 block:\n\n");
   for (i=0; i<H0block_->size; i++) {
     std::string configstring(print_config(CalcInfo_->num_ci_orbs, CalcInfo_->num_alp_expl,
         CalcInfo_->num_bet_expl, alplist_[H0block_->alplist[i]] +
         H0block_->alpidx[i], betlist_[H0block_->betlist[i]] +
         H0block_->betidx[i], CalcInfo_->num_drc_orbs));
      outfile->Printf( "  %3d [%3d] %10.6lf  Block %2d (%4d,%4d)  %s\n",
         i+1, H0block_->pair[i] + 1, H0block_->H00[i], H0block_->blknum[i],
         H0block_->alpidx[i], H0block_->betidx[i], configstring.c_str());
      }
}


int CIWavefunction::H0block_calc(double E)
{
   int i, j, size;
   double detH0 = -1.0;
   double c_tmp, s_tmp, tval1, tval2, tval3;
   double *H0xc0, *H0xs0;

   size = H0block_->size;

   if (print_ > 4) {
      outfile->Printf( "\nc0b = \n");
      print_mat(&(H0block_->c0b), 1, H0block_->size, "outfile");
      outfile->Printf( "\ns0b = \n");
      print_mat(&(H0block_->s0b), 1, H0block_->size, "outfile");
      }

   if (Parameters_->precon == PRECON_GEN_DAVIDSON) {
     H0xc0 = init_array(size);
     H0xs0 = init_array(size);
     for (i=0; i<size; i++) {
        for (j=0; j<size; j++)
           H0block_->H0b_diag_transpose[j] = H0block_->H0b_diag[j][i];
        dot_arr(H0block_->H0b_diag_transpose, H0block_->c0b, size, &H0xc0[i]);
        dot_arr(H0block_->H0b_diag_transpose, H0block_->s0b, size, &H0xs0[i]);
        }
     for (i=0; i<size; i++) {
        c_tmp = s_tmp = 0.0;
        for (j=0; j<size; j++) {
           tval1 = H0xc0[j] * H0block_->H0b_diag[i][j];
           tval2 = H0xs0[j] * H0block_->H0b_diag[i][j];
           tval3 = H0block_->H0b_eigvals[j] - E;
           if (fabs(tval3) < HD_MIN) tval3 = 0.0;
           else tval3 = 1.0/tval3;
           tval1 *= tval3;
           tval2 *= tval3;
           c_tmp += tval1;
           s_tmp += tval2;
           }
        H0block_->c0bp[i] = c_tmp;
        H0block_->s0bp[i] = s_tmp;
        }

     if (print_ > 4) {
        outfile->Printf( "\nc0b = \n");
        print_mat(&(H0block_->c0b),1,H0block_->size+H0block_->coupling_size,"outfile");
        outfile->Printf( "\nc0bp = \n");
        print_mat(&(H0block_->c0bp),1,H0block_->size+H0block_->coupling_size,"outfile");
        outfile->Printf( "\ns0b = \n");
        print_mat(&(H0block_->s0b), 1, H0block_->size, "outfile");
        outfile->Printf( "\ns0bp = \n");
        print_mat(&(H0block_->s0bp), 1, H0block_->size, "outfile");
        }
     free(H0xc0);
     free(H0xs0);
     return(1);
     }
   else if (Parameters_->precon == PRECON_H0BLOCK_INVERT ||
            Parameters_->precon == PRECON_H0BLOCK_ITER_INVERT) {

     /* form H0b-E and take its inverse */
     /* subtract E from the diagonal */
     for (i=0; i<size; i++) {
        H0block_->c0bp[i] = H0block_->c0b[i]; /* necessary for pople */
        H0block_->s0bp[i] = H0block_->s0b[i]; /* also necessary for pople */
        for (j=0; j<size; j++) {
           H0block_->tmp1[i][j] = H0block_->H0b[i][j];
           if (i==j) H0block_->tmp1[i][i] -= E;
           }
        }

     if (print_ > 4) {
        outfile->Printf( "\n E = %lf\n", E);
        outfile->Printf( " H0 - E\n");
        print_mat(H0block_->tmp1, H0block_->size, H0block_->size, "outfile");
        }

     if (Parameters_->precon == PRECON_H0BLOCK_ITER_INVERT) {
       pople(H0block_->tmp1, H0block_->c0bp, size, 1, 1e-9, "outfile",
             print_);
       if (Parameters_->update == UPDATE_OLSEN) {
         for (i=0; i<size; i++)
            for (j=0; j<size; j++) {
               H0block_->tmp1[i][j] = H0block_->H0b[i][j];
               if (i==j) H0block_->tmp1[i][i] -= E;
               }
         pople(H0block_->tmp1,H0block_->s0bp,size,1,1e-9,"outfile",
               print_);
         }
       detH0 = 1.0;
       }
     else {
       detH0 = invert_matrix(H0block_->tmp1, H0block_->H0b_inv, size, "outfile");
       if (print_ > 4) {
         outfile->Printf( "\nINV(H0 - E)\n");
         print_mat(H0block_->H0b_inv, H0block_->size, H0block_->size, "outfile");
         }

       /* get c0bp = (H0b - E)^{-1} * c0b */
       mmult(H0block_->H0b_inv, 0, &(H0block_->c0b), 1, &(H0block_->c0bp), 1,
             size, size, 1, 0);

       /* get s0bp = (H0b - E)^{-1} * s0b */
       mmult(H0block_->H0b_inv, 0, &(H0block_->s0b), 1, &(H0block_->s0bp), 1,
             size, size, 1, 0);
       }

     if (print_ > 4) {
        outfile->Printf( "\nc0b = \n");
        print_mat(&(H0block_->c0b), 1, H0block_->size, "outfile");
        outfile->Printf( "\nc0bp = \n");
        print_mat(&(H0block_->c0bp), 1, H0block_->size, "outfile");
        outfile->Printf( "\ns0b = \n");
        print_mat(&(H0block_->s0b), 1, H0block_->size, "outfile");
        outfile->Printf( "\ns0bp = \n");
        print_mat(&(H0block_->s0bp), 1, H0block_->size, "outfile");
        outfile->Printf("DET H0 = %5.4E\n", detH0);
        }

      if (detH0 < SMALL_DET) return(0);
      else  return(1);
      }

   return 0;
}


/*
** eventually replace this with a somewhat more efficient CIvect member
** function which employs the new H0block_->buf_member matrix.
**
** cscode == 0 refers to c0b, while cscode == 1 refers to s0b
**
*/
void CIWavefunction::H0block_gather(double **mat, int al, int bl, int cscode, int mscode,
   int phase)
{

   double *target;
   int i, aidx, bidx;

   if (cscode == 0)
      target = H0block_->c0b;
   else if (cscode == 1)
      target = H0block_->s0b;
   else {
     outfile->Printf("(H0block_gather): invalid cscode\n");
      return;
      }

   for (i=0; i<(H0block_->size+H0block_->coupling_size); i++) {
      if (H0block_->alplist[i] == al && H0block_->betlist[i] == bl) {
         aidx = H0block_->alpidx[i];
         bidx = H0block_->betidx[i];
         target[i] = mat[aidx][bidx];
         }
      if (mscode && H0block_->alplist[i] == bl && H0block_->betlist[i] == al) {
         aidx = H0block_->alpidx[i];
         bidx = H0block_->betidx[i];
         if (phase==1) target[i] = mat[bidx][aidx] ;
         else target[i] = -mat[bidx][aidx];
         }
      }
}


/*
** Calculate the contributions to x and y due to the H0 block
*/
void CIWavefunction::H0block_xy(double *x, double *y, double E)
{
   int i;
   double tx=0.0, ty=0.0, tval, c;


   for (i=0; i<H0block_->size; i++) {
      tval = H0block_->H00[i] - E;
      c = H0block_->c0b[i];
      if (fabs(tval) < HD_MIN) tval = HD_MIN; /* prevent /0 */
      tval = 1.0 / tval;
      tx += c * c * tval;
      ty += c * H0block_->s0b[i] * tval;
      }

   *x -= tx;
   *y -= ty;

  /*
   outfile->Printf("-tx = %lf -ty = %lf\n",tx,ty);
   for (i=0; i<H0block_->size; i++)
      outfile->Printf("H0block_->c0b[%d] = %lf\n",i,H0block_->c0b[i]);
   for (i=0; i<H0block_->size; i++)
      outfile->Printf("H0block_->c0bp[%d] = %lf\n",i,H0block_->c0bp[i]);
   for (i=0; i<H0block_->size; i++)
      outfile->Printf("H0block_->s0b[%d] = %lf\n",i,H0block_->s0b[i]);
   for (i=0; i<H0block_->size; i++)
      outfile->Printf("H0block_->s0bp[%d] = %lf\n",i,H0block_->s0bp[i]);
  */

   dot_arr(H0block_->c0b, H0block_->c0bp, H0block_->size, &tx);
   *x += tx;
   dot_arr(H0block_->s0b, H0block_->c0bp, H0block_->size, &ty);
 /*
   dot_arr(H0block_->c0b, H0block_->s0bp, H0block_->size, &ty);
 */
   *y += ty;
   /* outfile->Printf("+tx = %lf +ty = %lf\n",tx,ty); */

}

void CIWavefunction::H0block_setup(int num_blocks, int *Ia_code, int *Ib_code)
{
   int member;
   int ac, bc, ai, bi;
   int q, found, size;

   size = H0block_->size + H0block_->coupling_size;
   for (member=0; member < size; member++) {
      ac = H0block_->alplist[member];
      ai = H0block_->alpidx[member];
      bc = H0block_->betlist[member];
      bi = H0block_->betidx[member];

      /* set the "pairs" entry */
      for (q=0,found=0; q<size && !found; q++) {
         if (H0block_->alplist[q] == bc && H0block_->betlist[q] == ac &&
               H0block_->alpidx[q] == bi && H0block_->betidx[q] == ai) {
            H0block_->pair[member] = q;
            found=1;
            }
         }
      if (!found) H0block_->pair[member] = -1;

      /* set the blknum array */
      for (q=0,found=0; q<num_blocks && !found; q++) {
         if (Ia_code[q] == ac && Ib_code[q] == bc) {
            H0block_->blknum[member] = q;
            found = 1;
            }
         }
      if (!found) {
         H0block_->blknum[member] = -1;
        outfile->Printf("(H0block_setup): Can't find CI block!\n");
         }

      } /* end loop over members of H0block */


}

/*
** H0block_pairup
**
** This function makes sure that the H0 block members are all paired up.
** If a member does not have a pair in the block, it is removed from the
** block in order to ensure that the proper spin symmetry is present
** in the eigenvectors of the H0 block.  This is only appropriate for
** Ms=0 cases; otherwise, the function should not be called.  This
** turns out to be nontrivial, as the pairing structure can be in
** principle rather mixed-up.  Thus I'll opt for a recursive solution.
**
** All parameters are taken from the H0block structure.
**
** Author: David Sherrill, November 1996
** Modified by Matt Leininger July 1998
**
*/
void CIWavefunction::H0block_pairup(int guess)
{
  int i,first,newsize,size;

  size = 0;
  if (guess==2) {
    size = H0block_->size + H0block_->coupling_size;
    if (H0block_->coupling_size == 0) return;
    }
  else if (guess==1) size = H0block_->guess_size;
  else if (guess==0) size = H0block_->size;

  /* nothing to be done if no H0block */
  if (size < 1) return;

  for (first=-1,i=0; i<size; i++) {
    if (H0block_->pair[i] == -1) {
      first = i;
      break;
    }
  }

  if (first == -1) return; /* all paired up */
  else {
    newsize = first;
    for (i=0; i<newsize; i++) {
      if (H0block_->pair[i] >= newsize) H0block_->pair[i] = -1;
    }
  }

  if (first == 0) {
    outfile->Printf( "    Warning!  H0block size reduced to zero by ");
    outfile->Printf( "    H0block_pairup!\n");
  }

  if (guess==2) H0block_->coupling_size = newsize - H0block_->size;
  else if (guess==1) H0block_->guess_size = newsize;
  else if (guess==0) H0block_->size = newsize;

  H0block_pairup(guess);

}

/*
** H0block_spin_cpl_chk
**
** This function makes sure that all element of a spin-coupling set
** are included in H0 block. If a member does not have all determinants
** of the spin-coupling set in the block, it is removed from the
** block in order to ensure that the proper spin symmetry is present
** in the eigenvectors of the H0 block. This only works if the
** averaged diagonal elements over a spin-coupling set are used.
**
** All parameters are taken from the H0block structure.
**
** Author: Matt Leininger, February 1998
** Substantial modifications by David Sherrill, October 1998
**
*/
void CIWavefunction::H0block_spin_cpl_chk(void)
{
  int i, newsize;
  double zero = 1E-13;
  double diff = 0.0, spin_cpl_vals2;

  /* nothing to be done if no H0block */
  if (H0block_->size > 0) {

    if (H0block_->coupling_size > 0)
      spin_cpl_vals2 = H0block_->H00[H0block_->size];
    else
      spin_cpl_vals2 = H0block_->spin_cp_vals;

    i = H0block_->size-1;
    diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
   /* outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && diff < zero) {
      i--;
      diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
    /*  outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    }

    newsize = i+1;

    if (newsize == 0) {
      outfile->Printf( "    Warning!  H0block size reduced to zero by ");
      outfile->Printf( "    H0block_spin_cpl_chk!\n");
    }
    H0block_->size = newsize;
  }

  /****************************************************************
  ** Also need to check the H0 block of the initial guess which may
  ** be smaller than the H0block_->size
  ** CDS: I am assuming that guess_size has to be <= H0block_->size
  *******************************************************************/
  if (H0block_->guess_size > 0) {

    if (H0block_->guess_size >= H0block_->osize) {
      if (H0block_->coupling_size > 0)
        spin_cpl_vals2 = H0block_->H00[H0block_->size];
      else
        spin_cpl_vals2 = H0block_->spin_cp_vals;
      newsize = H0block_->osize;
    }
    else {
      spin_cpl_vals2 = H0block_->H00[H0block_->guess_size];
      newsize = H0block_->guess_size;
    }

    i = newsize - 1;
    diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
   /* outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && fabs(diff) < zero) {
      i--;
      diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
     /* outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    }

    newsize = i+1;
    if (newsize == 0) {
      outfile->Printf( "    Warning!  H0block guess size reduced to zero by ");
      outfile->Printf( "    H0block_spin_cpl_chk!\n");
    }

    H0block_->guess_size = newsize;
  }

  /****************************************************************
  ** Also need to check the H0 block of the h0block coupling which
  ** will be larger than the H0block_->size
  *****************************************************************/
  if (H0block_->coupling_size > 0) {

    spin_cpl_vals2 = H0block_->spin_cp_vals;
    newsize = H0block_->size + H0block_->coupling_size;

    i = newsize - 1;
    diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
   /* outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && fabs(diff) < zero) {
      i--;
      diff = fabs(H0block_->H00[i] - spin_cpl_vals2);
     /* outfile->Printf("diff[%d] = %20.15f\n", i, diff); */
    }

    newsize = i+1;

    if (newsize < H0block_->size) {
      outfile->Printf( "    H0block coupling size reduced below 0 ???\n");
      newsize = H0block_->size;
    }

    if (newsize == H0block_->size) {
      outfile->Printf(
        "    Warning! H0block coupling size reduced to H0block size by ");
      outfile->Printf( "    H0block_spin_cpl_chk!\n");
    }

    H0block_->coupling_size = newsize - H0block_->size;
   }

}

/*
** H0block_filter_setup()
**
** This function takes a pair of user-specified
** determinants and adds them to the most H0block space (at the expense
** of two previously determined determinants, if they aren't already
** present) and stores which H0block determinant numbers they are.
**
** C. David Sherrill, July 2003
*/
void CIWavefunction::H0block_filter_setup(void)
{
  int Iac, Ibc, Iaridx, Ibridx;
  int Jac, Jbc, Jaridx, Jbridx;
  int i, found1, found2, replace;

  Iac = Parameters_->filter_guess_Iac;
  Ibc = Parameters_->filter_guess_Ibc;
  Iaridx = Parameters_->filter_guess_Iaridx;
  Ibridx = Parameters_->filter_guess_Ibridx;
  Jac = Parameters_->filter_guess_Jac;
  Jbc = Parameters_->filter_guess_Jbc;
  Jaridx = Parameters_->filter_guess_Jaridx;
  Jbridx = Parameters_->filter_guess_Jbridx;

  /* figure out if these determinants are in the list already */
  found1 = 0;
  for (i=0; i<H0block_->size && !found1; i++) {
    if (H0block_->alplist[i] == Iac && H0block_->alpidx[i] == Iaridx &&
        H0block_->betlist[i] == Ibc && H0block_->betidx[i] == Ibridx) {
      found1 = 1;
      Parameters_->filter_guess_H0_det1 = i;
    }
  }

  found2 = 0;
  for (i=0; i<H0block_->size && !found2; i++) {
    if (H0block_->alplist[i] == Jac && H0block_->alpidx[i] == Jaridx &&
        H0block_->betlist[i] == Jbc && H0block_->betidx[i] == Jbridx) {
      found2 = 1;
      Parameters_->filter_guess_H0_det2 = i;
    }
  }

  /* if not in the H0block list, we need to force it there.  remove the
     last one in the list and replace it with the determinant I
   */
  if (!found1) {
    replace = H0block_->size-1;
    H0block_->alplist[replace] = Iac;
    H0block_->alpidx[replace] = Iaridx;
    H0block_->betlist[replace] = Ibc;
    H0block_->betidx[replace] = Ibridx;
    Parameters_->filter_guess_H0_det1 = replace;
  }
  if (!found2) {
    if (found1) replace = H0block_->size-1;
    else replace = H0block_->size-2;
    H0block_->alplist[replace] = Jac;
    H0block_->alpidx[replace] = Jaridx;
    H0block_->betlist[replace] = Jbc;
    H0block_->betidx[replace] = Jbridx;
    Parameters_->filter_guess_H0_det2 = replace;
  }

}

void CIWavefunction::H0block_fill()
{
   int i, j, size;
   int Ia, Ib, Ja, Jb;
   int Ialist, Iblist;
   SlaterDeterminant I, J;
   double *evals, **evecs;

   /* fill lower triangle */
   for (i=0; i<H0block_->size; i++) {

      Ialist = H0block_->alplist[i];
      Iblist = H0block_->betlist[i];
      Ia = H0block_->alpidx[i];
      Ib = H0block_->betidx[i];
      I.set(CalcInfo_->num_alp_expl,
          alplist_[Ialist][Ia].occs, CalcInfo_->num_bet_expl,
          betlist_[Iblist][Ib].occs);
      for (j=0; j<=i; j++) {
         Ialist = H0block_->alplist[j];
         Iblist = H0block_->betlist[j];
         Ia = H0block_->alpidx[j];
         Ib = H0block_->betidx[j];
         J.set(CalcInfo_->num_alp_expl,
            alplist_[Ialist][Ia].occs, CalcInfo_->num_bet_expl,
            betlist_[Iblist][Ib].occs);

         /* pointers in next line avoids copying structures I and J */
         H0block_->H0b[i][j] = matrix_element(&I, &J);
         if (i==j) H0block_->H0b[i][i] += CalcInfo_->edrc;
         /* outfile->Printf(" i = %d   j = %d\n",i,j); */
         }

      H0block_->H00[i] = H0block_->H0b[i][i];
      }

   /* fill upper triangle */
   fill_sym_matrix(H0block_->H0b, H0block_->size);

   if (Parameters_->precon == PRECON_GEN_DAVIDSON)
     size = H0block_->size;
   else
     size = H0block_->guess_size;

   if (print_ > 2) {
     outfile->Printf("H0block size = %d in H0block_fill\n",H0block_->size);
     outfile->Printf(
             "H0block guess size = %d in H0block_fill\n",H0block_->guess_size);
     outfile->Printf(
             "H0block coupling size = %d in H0block_fill\n",
             H0block_->coupling_size);
     outfile->Printf("Diagonalizing H0block_->H0b size %d in h0block_fill in"
                     " detci.cc ... ", size);

   }

   sq_rsp(size, size, H0block_->H0b, H0block_->H0b_eigvals, 1,
          H0block_->H0b_diag, 1.0E-14);

   if (print_) {
      outfile->Printf( "    H0 Block Eigenvalue = %12.8lf\n",
             H0block_->H0b_eigvals[0] + CalcInfo_->enuc);

      }

   if (print_ > 5 && size < 1000) {
      for (i=0; i<size; i++) H0block_->H0b_eigvals[i] += CalcInfo_->enuc;
      outfile->Printf( "\nH0 Block Eigenvectors\n");
      eivout(H0block_->H0b_diag, H0block_->H0b_eigvals,
             size, size, "outfile");
      outfile->Printf( "\nH0b matrix\n");
      print_mat(H0block_->H0b, size, size, "outfile");
      }
}

void CIWavefunction::H0block_coupling_calc(double E)
{
   int i, j, size, size2;
   double tval1, tval2;
   double *delta_2, *gamma_1, *gamma_2, *H_12, *delta_1;
   SlaterDeterminant I, J;
   int Ia, Ib;
   int Ialist, Iblist;

   size = H0block_->size;
   size2 = H0block_->size + H0block_->coupling_size;

   H_12 = init_array(H0block_->coupling_size);
   delta_1 = init_array(H0block_->size);
   delta_2 = init_array(H0block_->coupling_size);
   gamma_1 = init_array(H0block_->size);
   gamma_2 = init_array(H0block_->coupling_size);

   if (print_ > 5) {
      outfile->Printf( "\nc0b in H0block_coupling_calc = \n");
      print_mat(&(H0block_->c0b), 1, size2, "outfile");
      outfile->Printf( "\nc0bp in H0block_coupling_calc = \n");
      print_mat(&(H0block_->c0bp), 1, size2, "outfile");
      }

     /* copy to delta_1 */
     for (i=0; i<size; i++)
        delta_1[i] = H0block_->c0bp[i];

     /* form delta_2 array  (D-E)^-1 r_2 */
     for (i=size; i<size2; i++) {
        tval1 = H0block_->H00[i] - E;
        if (fabs(tval1) > HD_MIN)
          H0block_->c0bp[i] = H0block_->c0b[i]/tval1;
        else H0block_->c0bp[i] = 0.0;
        delta_2[i-size] = H0block_->c0bp[i];
        }
/*
     for (i=0; i<size2; i++)
        outfile->Printf("In Hcc H0block_->c0bp[%d] = %lf\n", i, H0block_->c0bp[i]);
*/

     zero_arr(gamma_2, size);
     /* Construct H_12 coupling block on-the-fly */
     for (i=0; i<size; i++) {
        Ialist = H0block_->alplist[i];
        Iblist = H0block_->betlist[i];
        Ia = H0block_->alpidx[i];
        Ib = H0block_->betidx[i];
        I.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Ia].occs,
              CalcInfo_->num_bet_expl, betlist_[Iblist][Ib].occs);
        for (j=size; j<size2; j++) {
           Ialist = H0block_->alplist[j];
           Iblist = H0block_->betlist[j];
           Ia = H0block_->alpidx[j];
           Ib = H0block_->betidx[j];
           J.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Ia].occs,
                 CalcInfo_->num_bet_expl, betlist_[Iblist][Ib].occs);
           H_12[j-size] = matrix_element(&I, &J);
           } /* end loop over j */

        dot_arr(H_12, delta_2, H0block_->coupling_size, &tval2);
        gamma_1[i] = tval2;
        for (j=0; j<H0block_->coupling_size; j++)
           gamma_2[j] += H_12[j] * delta_1[i];

        } /* end loop over i */


     /* Construct delta_1 = (H_11)^-1 gamma_1, delta_2 = (D_2-E)^-1 * gamma_2 */
     /* First delta_2 */
     for (i=size; i<size2; i++) {
        tval1 = H0block_->H00[i] - E;
        if (fabs(tval1) > HD_MIN)
          delta_2[i-size] = gamma_2[i-size]/tval1;
        else delta_2[i-size] = 0.0;
        }

     /* Now delta_1 */

     /* form H0b-E and take its inverse */
     for (i=0; i<size; i++) {
        delta_1[i] = gamma_1[i];
        for (j=0; j<size; j++) {
           H0block_->tmp1[i][j] = H0block_->H0b[i][j];
           if (i==j) H0block_->tmp1[i][i] -= E;
           }
        }

     if (print_ > 4) {
        outfile->Printf( "\n E = %lf\n", E);
        outfile->Printf( " H0 - E\n");
        print_mat(H0block_->tmp1, H0block_->size, H0block_->size, "outfile");
        }

/*
       for (i=0; i<size; i++)
          outfile->Printf("gamma_1[%d] = %lf\n", i, gamma_1[i]);

       pople(H0block_->tmp1, delta_1, size, 1, 1e-9, outfile,
             print_);
*/
       flin(H0block_->tmp1, delta_1, size, 1, &tval1);

     /*
       detH0 = invert_matrix(H0block_->tmp1, H0block_->H0b_inv, size, outfile);
       mmult(H0block_->H0b_inv,0,&(gamma_1),1,&(delta_1),1,size,size,1,0);
     */

    /*
       if (Parameters_->update == UPDATE_OLSEN) {
         for (i=0; i<size; i++)
            for (j=0; j<size; j++) {
               H0block_->tmp1[i][j] = H0block_->H0b[i][j];
               if (i==j) H0block_->tmp1[i][i] -= E;
               }
         pople(H0block_->tmp1,H0block_->s0bp,size,1,1e-9,outfile,
               print_);
         }
    */

      /* Construction of delta_1 and delta_2 completed */
      /* Now modify correction vectors in H0block structure */

      for (i=0; i<size; i++) H0block_->c0bp[i] -= delta_1[i];
      for (i=size; i<size2; i++) H0block_->c0bp[i] -= delta_2[i-size];

     /*
      for (i=0; i<size2; i++) {
         if (i>=H0block_->coupling_size)
           H0block_->c0bp[i] -= delta_2[i-size];
         else H0block_->c0bp[i] -= delta_1[i];
         }
    */

   free(H_12); free(delta_1);
   free(delta_2); free(gamma_1);
   free(gamma_2);
}


}} // namespace psi::detci
