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
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "ci_tol.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

#define SMALL_DET 1e-10

extern struct stringwr **alplist;
extern struct stringwr **betlist;

#define CONFIG_STRING_MAX 200

extern void print_config(int nbf, int num_alp_el, int num_bet_el,
   struct stringwr *stralp, struct stringwr *strbet,
   int num_fzc_orbs, char *outstring);


/*
** H0block_init()
**
** initialize everything but buf_num and buf_member, which depend on the
** CIvector structure
*/
void H0block_init(unsigned int size) {

   unsigned int size2;

   if (size > Parameters.h0blocksize) H0block.size = Parameters.h0blocksize;
   else H0block.size = size;
   H0block.coupling_size = Parameters.h0block_coupling_size;

   if (H0block.coupling_size) 
     size2 = H0block.size + H0block.coupling_size;
   else size2 = H0block.size;

   if (Parameters.print_lvl > 1)
     fprintf(outfile,"Total H0block size (including coupling): %d\n",size2);
   
   H0block.osize = H0block.size;
   H0block.guess_size = Parameters.h0guess_size;
   H0block.oguess_size = H0block.guess_size;
   H0block.ocoupling_size = H0block.coupling_size;

   if (H0block.size) {
      H0block.H0b = init_matrix(H0block.size, H0block.size);
      if (Parameters.precon == PRECON_GEN_DAVIDSON) 
        H0block.H0b_diag_transpose = init_array(H0block.size);
    /*  H0block.H0b_diag_transpose = init_matrix(H0block.size, H0block.size); */
      H0block.H0b_diag = init_matrix(H0block.size, H0block.size);
      H0block.H0b_eigvals = init_array(H0block.size);
      /* if (Parameters.precon == PRECON_H0BLOCK_INVERT || 
          Parameters.precon == PRECON_H0BLOCK_ITER_INVERT) */
      H0block.tmp1 = init_matrix(H0block.size, H0block.size); 
      if (Parameters.precon == PRECON_H0BLOCK_INVERT)
        H0block.H0b_inv = init_matrix(H0block.size, H0block.size);
      H0block.H00 = init_array(size2);
      H0block.c0b = init_array(size2);
      H0block.c0bp = init_array(size2);
      H0block.s0b = init_array(size2);
      H0block.s0bp = init_array(size2);
      H0block.alplist = init_int_array(size2);
      H0block.betlist = init_int_array(size2);
      H0block.alpidx = init_int_array(size2);
      H0block.betidx = init_int_array(size2);
      H0block.blknum = init_int_array(size2);
      H0block.pair = init_int_array(size2);

      if (Parameters.h0block_coupling) {
        H0block.tmp_array1 = init_array(size2);
        H0block.tmp_array2 = init_array(size2);
        }
     }
}


void H0block_free(void) 
{
   int i;

   if (H0block.osize) {
      free(H0block.H00);
      free(H0block.c0b);
      free(H0block.c0bp);
      free(H0block.s0b);
      free(H0block.s0bp);
      free(H0block.alplist);
      free(H0block.betlist);
      free(H0block.alpidx);
      free(H0block.betidx);
      free(H0block.blknum);
      if (Parameters.precon == PRECON_GEN_DAVIDSON)
        free(H0block.H0b_diag_transpose);
      free_matrix(H0block.H0b, H0block.osize);
      if (Parameters.precon == PRECON_H0BLOCK_INVERT)
        free_matrix(H0block.H0b_inv, H0block.osize);
    /*  if (Parameters.precon == PRECON_H0BLOCK_INVERT || 
          Parameters.precon == PRECON_H0BLOCK_ITER_INVERT) */
        free_matrix(H0block.tmp1, H0block.osize); 
      free(H0block.pair);
      if (H0block.nbuf) {
         free(H0block.buf_num);
         for (i=0; i<H0block.nbuf; i++) free(H0block.buf_member[i]);
         free(H0block.buf_member);
         }
      }
}

void H0block_print(void)
{
   int i;
   char configstring[CONFIG_STRING_MAX];

   fprintf(outfile, "\nMembers of H0 block:\n\n");
   for (i=0; i<H0block.size; i++) {
      print_config(CalcInfo.num_ci_orbs, CalcInfo.num_alp_expl, 
         CalcInfo.num_bet_expl, alplist[H0block.alplist[i]] + 
         H0block.alpidx[i], betlist[H0block.betlist[i]] + 
         H0block.betidx[i], CalcInfo.num_fzc_orbs,configstring);
      fprintf(outfile, "  %3d [%3d] %10.6lf  Block %2d (%4d,%4d)  %s\n", 
         i+1, H0block.pair[i] + 1, H0block.H00[i], H0block.blknum[i],
         H0block.alpidx[i], H0block.betidx[i], configstring);
      }
}


int H0block_calc(double E)
{
   static int first_call = 1;
   int i, j, size;
   double detH0 = -1.0;
   double c_tmp, s_tmp, tval1, tval2, tval3;
   double *H0xc0, *H0xs0;

   size = H0block.size;

   if (Parameters.print_lvl > 4) {
      fprintf(outfile, "\nc0b = \n");
      print_mat(&(H0block.c0b), 1, H0block.size, outfile);
      fprintf(outfile, "\ns0b = \n");
      print_mat(&(H0block.s0b), 1, H0block.size, outfile); 
      }

   if (Parameters.precon == PRECON_GEN_DAVIDSON) {
     if (first_call) {
       first_call = 0;
     /*  for (i=0; i<size; i++)
          for (j=0; j<size; j++)  
             H0block.H0b_diag_transpose[i][j] = H0block.H0b_diag[j][i];
     */
       }
     H0xc0 = init_array(size);
     H0xs0 = init_array(size);
     for (i=0; i<size; i++) {
        for (j=0; j<size; j++) 
           H0block.H0b_diag_transpose[j] = H0block.H0b_diag[j][i];
        dot_arr(H0block.H0b_diag_transpose, H0block.c0b, size, &H0xc0[i]);
        dot_arr(H0block.H0b_diag_transpose, H0block.s0b, size, &H0xs0[i]);
        }
     for (i=0; i<size; i++) {
        c_tmp = s_tmp = 0.0;
        for (j=0; j<size; j++) {
           tval1 = H0xc0[j] * H0block.H0b_diag[i][j];
           tval2 = H0xs0[j] * H0block.H0b_diag[i][j];
           tval3 = H0block.H0b_eigvals[j] - E; 
           if (fabs(tval3) < HD_MIN) tval3 = 0.0;
           else tval3 = 1.0/tval3;
           tval1 *= tval3;
           tval2 *= tval3;
           c_tmp += tval1;
           s_tmp += tval2;
           }
        H0block.c0bp[i] = c_tmp;
        H0block.s0bp[i] = s_tmp;
        }

     if (Parameters.print_lvl > 4) {
        fprintf(outfile, "\nc0b = \n");
        print_mat(&(H0block.c0b),1,H0block.size+H0block.coupling_size,outfile);
        fprintf(outfile, "\nc0bp = \n");
        print_mat(&(H0block.c0bp),1,H0block.size+H0block.coupling_size,outfile);
        fprintf(outfile, "\ns0b = \n");
        print_mat(&(H0block.s0b), 1, H0block.size, outfile); 
        fprintf(outfile, "\ns0bp = \n");
        print_mat(&(H0block.s0bp), 1, H0block.size, outfile); 
        }
     free(H0xc0);
     free(H0xs0);
     return(1);
     }
   else if (Parameters.precon == PRECON_H0BLOCK_INVERT || 
            Parameters.precon == PRECON_H0BLOCK_ITER_INVERT) {

     /* form H0b-E and take its inverse */
     /* subtract E from the diagonal */
     for (i=0; i<size; i++) {
        H0block.c0bp[i] = H0block.c0b[i]; /* necessary for pople */
        H0block.s0bp[i] = H0block.s0b[i]; /* also necessary for pople */
        for (j=0; j<size; j++) {
           H0block.tmp1[i][j] = H0block.H0b[i][j];
           if (i==j) H0block.tmp1[i][i] -= E; 
           }
        }

     if (Parameters.print_lvl > 4) {
        fprintf(outfile, "\n E = %lf\n", E);
        fprintf(outfile, " H0 - E\n");
        print_mat(H0block.tmp1, H0block.size, H0block.size, outfile);
        }

     if (Parameters.precon == PRECON_H0BLOCK_ITER_INVERT) {
       pople(H0block.tmp1, H0block.c0bp, size, 1, 1e-9, outfile, 
             Parameters.print_lvl);
       if (Parameters.update == UPDATE_OLSEN) {
         for (i=0; i<size; i++) 
            for (j=0; j<size; j++) {
               H0block.tmp1[i][j] = H0block.H0b[i][j];
               if (i==j) H0block.tmp1[i][i] -= E;
               }
         pople(H0block.tmp1,H0block.s0bp,size,1,1e-9,outfile,
               Parameters.print_lvl);
         }
       detH0 = 1.0;
       }
     else {
       detH0 = invert_matrix(H0block.tmp1, H0block.H0b_inv, size, outfile);
       if (Parameters.print_lvl > 4) {
         fprintf(outfile, "\nINV(H0 - E)\n");
         print_mat(H0block.H0b_inv, H0block.size, H0block.size, outfile);
         }
 
       /* get c0bp = (H0b - E)^{-1} * c0b */
       mmult(H0block.H0b_inv, 0, &(H0block.c0b), 1, &(H0block.c0bp), 1,
             size, size, 1, 0);

       /* get s0bp = (H0b - E)^{-1} * s0b */
       mmult(H0block.H0b_inv, 0, &(H0block.s0b), 1, &(H0block.s0bp), 1,
             size, size, 1, 0);
       } 

     if (Parameters.print_lvl > 4) {
        fprintf(outfile, "\nc0b = \n");
        print_mat(&(H0block.c0b), 1, H0block.size, outfile);
        fprintf(outfile, "\nc0bp = \n");
        print_mat(&(H0block.c0bp), 1, H0block.size, outfile);
        fprintf(outfile, "\ns0b = \n");
        print_mat(&(H0block.s0b), 1, H0block.size, outfile); 
        fprintf(outfile, "\ns0bp = \n");
        print_mat(&(H0block.s0bp), 1, H0block.size, outfile); 
        fprintf(outfile,"DET H0 = %5.4E\n", detH0);
        }
  
      if (detH0 < SMALL_DET) return(0);
      else  return(1);
      }

}


/*
** eventually replace this with a somewhat more efficient CIvect member
** function which employs the new H0block.buf_member matrix.
** 
** cscode == 0 refers to c0b, while cscode == 1 refers to s0b 
**
*/
void H0block_gather(double **mat, int al, int bl, int cscode, int mscode,
   int phase)
{

   double *target;
   int i, aidx, bidx;

   if (cscode == 0) 
      target = H0block.c0b;
   else if (cscode == 1) 
      target = H0block.s0b;
   else {
      printf("(H0block_gather): invalid cscode\n"); 
      return;
      }

   for (i=0; i<(H0block.size+H0block.coupling_size); i++) {
      if (H0block.alplist[i] == al && H0block.betlist[i] == bl) {
         aidx = H0block.alpidx[i];
         bidx = H0block.betidx[i];
         target[i] = mat[aidx][bidx];
         }
      if (mscode && H0block.alplist[i] == bl && H0block.betlist[i] == al) {
         aidx = H0block.alpidx[i];
         bidx = H0block.betidx[i];
         if (phase==1) target[i] = mat[bidx][aidx] ;
         else target[i] = -mat[bidx][aidx];
         }
      }
}


/*
** Calculate the contributions to x and y due to the H0 block
*/
void H0block_xy(double *x, double *y, double E)
{
   int i;
   double tx=0.0, ty=0.0, tval, c;


   for (i=0; i<H0block.size; i++) {
      tval = H0block.H00[i] - E;
      c = H0block.c0b[i];
      if (fabs(tval) < HD_MIN) tval = HD_MIN; /* prevent /0 */
      tval = 1.0 / tval;
      tx += c * c * tval;
      ty += c * H0block.s0b[i] * tval;
      }

   *x -= tx;
   *y -= ty;

  /*
   fprintf(outfile,"-tx = %lf -ty = %lf\n",tx,ty);
   for (i=0; i<H0block.size; i++) 
      fprintf(outfile,"H0block.c0b[%d] = %lf\n",i,H0block.c0b[i]);
   for (i=0; i<H0block.size; i++) 
      fprintf(outfile,"H0block.c0bp[%d] = %lf\n",i,H0block.c0bp[i]);
   for (i=0; i<H0block.size; i++) 
      fprintf(outfile,"H0block.s0b[%d] = %lf\n",i,H0block.s0b[i]);
   for (i=0; i<H0block.size; i++) 
      fprintf(outfile,"H0block.s0bp[%d] = %lf\n",i,H0block.s0bp[i]);
  */

   dot_arr(H0block.c0b, H0block.c0bp, H0block.size, &tx);
   *x += tx;
   dot_arr(H0block.s0b, H0block.c0bp, H0block.size, &ty);
 /*
   dot_arr(H0block.c0b, H0block.s0bp, H0block.size, &ty); 
 */
   *y += ty;
   /* fprintf(outfile,"+tx = %lf +ty = %lf\n",tx,ty); */

}

void H0block_setup(int num_blocks, int *Ia_code, int *Ib_code)
{
   int member;
   int ac, bc, ai, bi;
   int q, found, size;

   size = H0block.size + H0block.coupling_size;
   for (member=0; member < size; member++) {
      ac = H0block.alplist[member];
      ai = H0block.alpidx[member];
      bc = H0block.betlist[member];
      bi = H0block.betidx[member];

      /* set the "pairs" entry */
      for (q=0,found=0; q<size && !found; q++) {
         if (H0block.alplist[q] == bc && H0block.betlist[q] == ac &&
               H0block.alpidx[q] == bi && H0block.betidx[q] == ai) {
            H0block.pair[member] = q;
            found=1;
            }
         }
      if (!found) H0block.pair[member] = -1;

      /* set the blknum array */
      for (q=0,found=0; q<num_blocks && !found; q++) {
         if (Ia_code[q] == ac && Ib_code[q] == bc) {
            H0block.blknum[member] = q;
            found = 1;
            }
         }
      if (!found) {
         H0block.blknum[member] = -1;
         printf("(H0block_setup): Can't find CI block!\n");
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
void H0block_pairup(int guess)
{
  int i,first,newsize,size;

  size = 0;
  if (guess==2) {
    size = H0block.size + H0block.coupling_size;
    if (H0block.coupling_size == 0) return;
    }
  else if (guess==1) size = H0block.guess_size;
  else if (guess==0) size = H0block.size;

  /* nothing to be done if no H0block */
  if (size < 1) return;

  for (first=-1,i=0; i<size; i++) {
    if (H0block.pair[i] == -1) {
      first = i;
      break;
    }
  }

  if (first == -1) return; /* all paired up */
  else {
    newsize = first;
    for (i=0; i<newsize; i++) {
      if (H0block.pair[i] >= newsize) H0block.pair[i] = -1;
    }
  }

  if (first == 0) {
    fprintf(outfile, "Warning!  H0block size reduced to zero by ");
    fprintf(outfile, "H0block_pairup!\n");
  }

  if (guess==2) H0block.coupling_size = newsize - H0block.size;
  else if (guess==1) H0block.guess_size = newsize;
  else if (guess==0) H0block.size = newsize;

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
void H0block_spin_cpl_chk(void) 
{
  int i,newsize,tmpsize;
  double zero = 1E-13;
  double diff = 0.0, spin_cpl_vals2;

  /* nothing to be done if no H0block */
  if (H0block.size > 0) { 

    if (H0block.coupling_size > 0)
      spin_cpl_vals2 = H0block.H00[H0block.size];
    else
      spin_cpl_vals2 = H0block.spin_cp_vals;

    i = H0block.size-1;
    diff = fabs(H0block.H00[i] - spin_cpl_vals2);
   /* fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && diff < zero) {
      i--;
      diff = fabs(H0block.H00[i] - spin_cpl_vals2);
    /*  fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    } 

    newsize = i+1;

    if (newsize == 0) {
      fprintf(outfile, "Warning!  H0block size reduced to zero by ");
      fprintf(outfile, "H0block_spin_cpl_chk!\n");
    }
    H0block.size = newsize;
  }  

  /****************************************************************
  ** Also need to check the H0 block of the initial guess which may
  ** be smaller than the H0block.size
  ** CDS: I am assuming that guess_size has to be <= H0block.size
  *******************************************************************/
  if (H0block.guess_size > 0) { 

    if (H0block.guess_size >= H0block.osize) { 
      if (H0block.coupling_size > 0)
        spin_cpl_vals2 = H0block.H00[H0block.size];
      else
        spin_cpl_vals2 = H0block.spin_cp_vals;
      newsize = H0block.osize;
    }
    else {
      spin_cpl_vals2 = H0block.H00[H0block.guess_size];
      newsize = H0block.guess_size;
    }

    i = newsize - 1;
    diff = fabs(H0block.H00[i] - spin_cpl_vals2);
   /* fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && fabs(diff) < zero) {
      i--;
      diff = fabs(H0block.H00[i] - spin_cpl_vals2);
     /* fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    }

    newsize = i+1;
    if (newsize == 0) {
      fprintf(outfile, "Warning!  H0block guess size reduced to zero by ");
      fprintf(outfile, "H0block_spin_cpl_chk!\n");
    }

    H0block.guess_size = newsize;
  }
 
  /****************************************************************
  ** Also need to check the H0 block of the h0block coupling which
  ** will be larger than the H0block.size
  *****************************************************************/
  if (H0block.coupling_size > 0) { 

    spin_cpl_vals2 = H0block.spin_cp_vals;
    newsize = H0block.size + H0block.coupling_size;

    i = newsize - 1;
    diff = fabs(H0block.H00[i] - spin_cpl_vals2);
   /* fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    while (i > 0 && fabs(diff) < zero) {
      i--;
      diff = fabs(H0block.H00[i] - spin_cpl_vals2);
     /* fprintf(outfile,"diff[%d] = %20.15f\n", i, diff); */
    }

    newsize = i+1;

    if (newsize < H0block.size) {
      fprintf(outfile, "H0block coupling size reduced below 0 ???\n");
      newsize = H0block.size;
    }

    if (newsize == H0block.size) {
      fprintf(outfile, 
        "Warning! H0block coupling size reduced to H0block size by ");
      fprintf(outfile, "H0block_spin_cpl_chk!\n");
    }

    H0block.coupling_size = newsize - H0block.size;
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
void H0block_filter_setup(void)
{
  int Iac, Ibc, Iaridx, Ibridx;
  int Jac, Jbc, Jaridx, Jbridx;
  int i, found1, found2, replace;

  Iac = Parameters.filter_guess_Iac;
  Ibc = Parameters.filter_guess_Ibc;
  Iaridx = Parameters.filter_guess_Iaridx;
  Ibridx = Parameters.filter_guess_Ibridx;
  Jac = Parameters.filter_guess_Jac;
  Jbc = Parameters.filter_guess_Jbc;
  Jaridx = Parameters.filter_guess_Jaridx;
  Jbridx = Parameters.filter_guess_Jbridx;

  /* figure out if these determinants are in the list already */
  found1 = 0;
  for (i=0; i<H0block.size && !found1; i++) {
    if (H0block.alplist[i] == Iac && H0block.alpidx[i] == Iaridx &&
        H0block.betlist[i] == Ibc && H0block.betidx[i] == Ibridx) {
      found1 = 1;
      Parameters.filter_guess_H0_det1 = i;
    }
  }

  found2 = 0;
  for (i=0; i<H0block.size && !found2; i++) {
    if (H0block.alplist[i] == Jac && H0block.alpidx[i] == Jaridx &&
        H0block.betlist[i] == Jbc && H0block.betidx[i] == Jbridx) {
      found2 = 1;
      Parameters.filter_guess_H0_det2 = i;
    }
  }

  /* if not in the H0block list, we need to force it there.  remove the
     last one in the list and replace it with the determinant I
   */
  if (!found1) {
    replace = H0block.size-1;
    H0block.alplist[replace] = Iac;
    H0block.alpidx[replace] = Iaridx;
    H0block.betlist[replace] = Ibc;
    H0block.betidx[replace] = Ibridx;
    Parameters.filter_guess_H0_det1 = replace;
  }
  if (!found2) {
    if (found1) replace = H0block.size-1;
    else replace = H0block.size-2;
    H0block.alplist[replace] = Jac;
    H0block.alpidx[replace] = Jaridx;
    H0block.betlist[replace] = Jbc;
    H0block.betidx[replace] = Jbridx;
    Parameters.filter_guess_H0_det2 = replace;
  }

}

}} // namespace psi::detci

