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
**  \brief Mitrushenkov iterative scheme for RAS CI's
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** 19 June 1995
**
** Matt L. Leininger
** Center for Computational Quantum Chemistry
** 18 February 1998
**
** Last updated 2/18/98 Added civect class calls and worked out numerous bugs
**
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

extern void xeaxmy(double *x, double *y, double a, int size);
extern void xeaxpby(double *x, double *y, double a, double b, int size);
extern void xexy(double *x, double *y, int size);
extern void buf_ols_denom(double *a, double *hd, double E, int len);
extern void buf_ols_updt(double *a, double *c, double *norm, double *ovrlap,
   double *c1norm, int len);
extern double buf_xy1(double *c, double *hd, double E, int len);


#define MITRUSH_E_DIFF_MIN 5.0E-6


/*
** mitrush_iter()
**
** Adapted from ci_iter() in the FCI program (5/95 and before).
** Performs Olsen and Mitrushenkov iterations to converge on the CI vector.
** New version uses the CIvect class.
**
** David Sherrill
** 19 June 1995
**
*/
void CIWavefunction::mitrush_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
      **betlist, int nroots, double *evals, double conv_rms, double conv_e,
      double enuc, double edrc, int maxiter, int maxnvect)
{

   int i, j, ij, k, l, curr, last, iter=0, L, tmpi;
   int num_alp_str, num_bet_str, detH0 = -1;
   double *oei, *tei;
   double E, E_curr, E_last, E_est, E12, norm=1.0, S;
   double **H2x2, *evals2x2, **evecs2x2, alast, acur;
   double x, y, c1norm = 0.0;
   int sm_tridim, buf;
   double *sm_mat, *sm_evals, **sm_evecs;
   int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
   double testS = 0.0;
   double tval, *mi_coeff, *buffer1, *buffer2;
   double **alpha, chknorm;
   int diag_method;
   CIvect Cvec;
   CIvect Sigma;

   Cvec.set(Parameters_->icore, maxnvect, 1,
            Parameters_->s_filenum, CIblks_);
   Sigma.set(Parameters_->icore, maxnvect, 1,
             Parameters_->s_filenum, CIblks_);

   // Open I/O files but not with OPEN_OLD
   Cvec.init_io_files(false);
   Sigma.init_io_files(false);

   /* set up the vector pointers/info */
   if (Cvec.read_new_first_buf() == -1) Cvec.write_new_first_buf();
   if (Sigma.read_new_first_buf() == -1) Sigma.write_new_first_buf();
   if (Cvec.read_num_vecs() == -1) Cvec.write_num_vecs(0);
   if (Sigma.read_num_vecs() == -1) Sigma.write_num_vecs(0);

   Cvec.h0block_buf_init();

   curr = 0;
   last = (Parameters_->diag_method == METHOD_MITRUSHENKOV) ? 1 : 0;
   diag_method = Parameters_->diag_method;

   /* set buffer pointers */

   buffer1 = *(Hd.blockptr(0));
   buffer2 = Hd.buf_malloc();
   Hd.buf_unlock();

   /* get some of the stuff from CalcInfo for easier access */

   num_alp_str = CalcInfo_->num_alp_str;
   num_bet_str = CalcInfo_->num_bet_str;
   if (Parameters_->fci) oei = CalcInfo_->tf_onel_ints->pointer();
   else oei = CalcInfo_->gmat->pointer();
   tei = CalcInfo_->twoel_ints->pointer();

   /* small arrays to hold most important config information */

   mi_iac = init_int_array(Parameters_->nprint);
   mi_ibc = init_int_array(Parameters_->nprint);
   mi_iaidx = init_int_array(Parameters_->nprint);
   mi_ibidx = init_int_array(Parameters_->nprint);
   mi_coeff = init_array(Parameters_->nprint);


   /* stuff for the 2x2 Davidson procedure */

   H2x2 = init_matrix(2,2);
   evals2x2 = init_array(2);
   evecs2x2 = init_matrix(2,2);
   alpha = init_matrix(1,1);
    /* not used but is necessary for call to olsen_iter_xy */

   /* setup initial guess vector */

   if (Parameters_->restart) {
     outfile->Printf("\nAttempting Restart with 1 vector\n");
     if ((i=Cvec.read_num_vecs())< 1) {
       throw PsiException("CI vector file should contain only 1 vector.",__FILE__,__LINE__);
     }

     Cvec.buf_lock(buffer1);
     Cvec.read(0, 0);
     tval = Cvec * Cvec;
     if ((tval - 1.0) > ZERO) outfile->Printf("CI vector may be corrupted."
     " Attempting to correct by renormalizing.\n");
     Cvec.symnorm(tval,CI_VEC,TRUE);
     }
   else if (Parameters_->guess_vector == PARM_GUESS_VEC_UNIT && nroots == 1 &&
      Parameters_->num_init_vecs == 1) { /* use unit vector */
      tval = 1.0;
      Cvec.buf_lock(buffer1);
      Cvec.init_vals(0, 1, &(CalcInfo_->ref_alp_list), &(CalcInfo_->ref_alp_rel),
                     &(CalcInfo_->ref_bet_list), &(CalcInfo_->ref_bet_rel),
                     H0block_->blknum, &tval);
      Cvec.write_num_vecs(1);
      k = 1;
      Cvec.read(0, 0);
      Cvec.symnorm(1.0,CI_VEC,TRUE);
      Sigma.set_zero_blocks_all();
      }

   else { /* use H0BLOCK eigenvector guess */
      if (Parameters_->precon == PRECON_GEN_DAVIDSON) L = H0block_->size;
      else L = H0block_->guess_size;
      sm_evals = init_array(L);
      /* need to fill out sm_evecs into b (pad w/ 0's) */
      outfile->Printf( "Using %d initial trial vectors\n",
         Parameters_->num_init_vecs);

      Cvec.buf_lock(buffer1);
      for (i=0,k=0; i<L && k<1; i++) {

         /* check sm_evecs[i] to see if it has the correct spin symmetry */
         /* if Ms=0 ... */
         for (j=0,tmpi=0; j<L && !tmpi && Parameters_->Ms0; j++) {
            l = H0block_->pair[j];
            if (l == -1) { tmpi = 1; break; }
            tval = H0block_->H0b_diag[l][i];
            if ((int) Parameters_->S % 2) tval = -tval;
            if (H0block_->H0b_diag[j][i] - tval > 1.0E-8) tmpi = 1;
            }
         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = H0block_->H0b_diag[j][i];

         Cvec.init_vals(k, L, H0block_->alplist, H0block_->alpidx,
            H0block_->betlist, H0block_->betidx, H0block_->blknum, sm_evals);
         Cvec.write_num_vecs(1);
         k++;
         /* MLL added 6-15-98 */
         Cvec.read(0, 0);
         Cvec.symnorm(1.0,CI_VEC,TRUE);
         }

      free(sm_evals);
      Sigma.set_zero_blocks_all();
      }

   if (print_ > 4) {
      outfile->Printf( "\nC(0) vector = \n");
      Cvec.print();
      }

   Sigma.buf_lock(buffer2);
   Cvec.read(0, 0);
   sigma(Cvec, Sigma, oei, tei, 0);
   Sigma.write_num_vecs(1);

   Cvec.copy_zero_blocks(Sigma);

   if (print_ > 4) {
      outfile->Printf( "\nSigma vector\n");
      Sigma.print();

      }

   /* get H00 */
   E = Cvec * Sigma;
   E += edrc;
   E_last = CalcInfo_->escf - CalcInfo_->enuc;

   /* Cvec.print(outfile); */
   Cvec.buf_unlock();
   Sigma.buf_unlock();


   /*
    * get y = C(0) * (Hd - E)^-1 * sigma(0)
    * and x = C(0) * (Hd - E)^-1 * C(0)
    * Mental note: x and y are doubles not vectors
    */
   olsen_iter_xy(Cvec,Sigma,Hd,&x,&y,buffer1,buffer2,E,curr,1,alpha,alplist,
                 betlist);

   if (print_ > 3) {
     outfile->Printf( "Straight x = %12.6lf\n", x);
     outfile->Printf( "Straight y = %12.6lf\n", y);
     }

   if (Parameters_->precon >= PRECON_GEN_DAVIDSON && H0block_->size) {
       detH0 = 1;
       /*
       detH0 = H0block_calc(E);
       */

       if (detH0 == 0) {
         outfile->Printf( "H0block inverse is nearly nonsingular:");
         outfile->Printf(" initiating DAVIDSON preconditioner\n");

         Parameters_->precon = PRECON_DAVIDSON;
         }
       if (Parameters_->precon >= PRECON_GEN_DAVIDSON) {
         /*
         H0block_xy(&x,&y,E);
         */
         if (print_ > 3) {
            outfile->Printf( "x = %12.6lf\n", x);
            outfile->Printf( "y = %12.6lf\n", y);
            }
         }
     }

   E_est = y / x;  /* should I add fzc here? */
   if (print_ > 2) {
      outfile->Printf( "E_est = %12.6lf E-edrc = %12.6lf E = %12.6lf\n",
           E_est,E-edrc,E);
       outfile->Printf( "x = %lf  y = %lf\n",x,y);
       }
    /* calculate delta_C and C(1) */
   olsen_update(Cvec, Sigma, Hd, E, E_est, &norm, &c1norm, &S, buffer1, buffer2,
                 curr, last, "outfile", iter, alplist, betlist);

   norm = sqrt(1.0 / norm);
   Cvec.buf_lock(buffer1);
   Cvec.read(last, 0);
   Cvec.symnorm(norm,CI_VEC,TRUE);
   S *= norm;
   Cvec.buf_unlock();
   if (Parameters_->calc_ssq && Parameters_->icore==1)
     Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, 0);

  /*
   Cvec.buf_lock(buffer1);
   Cvec.read(0,0);
   outfile->Printf(" Cvec[0] = \n");
   Cvec.print(outfile);
   Cvec.read(1,0);
   outfile->Printf(" Cvec[1] = \n");
   Cvec.print(outfile);
   Cvec.buf_unlock();
  */

   /* S is the overlap of the C_(i-1) and C_i */

   outfile->Printf("Iter  0  ROOT 1 ECI = %14.9lf", enuc + E);
   outfile->Printf("    Delta_E %10.3E   Delta_C %10.3E\n", E - E_last, c1norm);
   Process::environment.globals["DETCI AVG DVEC NORM"] = c1norm;


   iter = 1;
   if (Parameters_->diag_method == METHOD_MITRUSHENKOV) {
      curr = 1; last = 0;
      }

   while (1) {

      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);
      Cvec.read(curr,0);
      /* chknorm = Cvec.checknorm();
      outfile->Printf("Norm of curr CI vect = %lf\n",chknorm);
     */
      if (print_ > 4) {
         outfile->Printf( "\nC(%2d) vector (symm'd norm'd)\n", iter) ;
         Cvec.print();

         }

      Cvec.read(curr,0);
      /* Sigma.read(curr,0);
      */
      sigma(Cvec, Sigma, oei, tei, curr);
      Cvec.copy_zero_blocks(Sigma);
      Cvec.read(curr,0);
      if (print_ > 4) {
         outfile->Printf("\nC(%2d) vector (symm'd norm'd) second time\n",iter);
         Cvec.print();

         }

      if (print_ > 4) {
         Sigma.read(curr,0);
         outfile->Printf( "\n curr = %d\n", curr);
         outfile->Printf( "\nSigma(%2d) vector\n", iter);
         Sigma.print();

         }

      Cvec.read(curr,0);
      Sigma.read(curr,0);
      /* calculate H(ii) */
      E_curr = Cvec * Sigma;
      E_curr += edrc;
      E = E_curr;

      /* check for convergence and exit if reached */
      if (iter == maxiter)
         outfile->Printf( "Maximum number of iterations reached!\n");

      /* if the 2x2 matrix diagonalization can be done, take Mitrush Step */

      if (Parameters_->diag_method==METHOD_MITRUSHENKOV &&
          diag_method==METHOD_MITRUSHENKOV && S < S_MAX &&
          fabs(E_last-E_curr) > MITRUSH_E_DIFF_MIN) {
        outfile->Printf( "Taking Mitrushenkov step (S =%10.6lf <%10.6lf)\n",
                 S, S_MAX);
         /* calculate H(i,i-1) = H(i-1,i) */
        Cvec.read(last, 0);
        E12 = Cvec * Sigma;
        E12 += edrc * S;

        /* fill up little H matrix and solve 2x2 eigenvalue problem */
        H2x2[0][0] = E_last;
        H2x2[0][1] = H2x2[1][0] = E12;
        H2x2[1][1] = E_curr;

        solve_2x2_pep(H2x2, S, evals2x2, evecs2x2);

        /* recalculate C(i) = alpha(i-1) * C(i-1) + alpha(i) * C(i)     */
        /* and sigma(i) = alpha(i-1) * sigma(i-1) + alpha(i) * sigma(i) */
        alast = evecs2x2[0][0];
        acur = evecs2x2[0][1];
        norm = 1.0/sqrt(acur * acur + alast * alast + 2.0 * acur * alast * S);

        Cvec.buf_unlock();
        Sigma.buf_unlock();

        /*
        outfile->Printf("alpha_0 = %lf  alpha_1 = %lf\n",
                evecs2x2[0][0],evecs2x2[0][1]);
        Cvec.buf_lock(buffer1);
        Cvec.read(0,0);
        outfile->Printf(" Cvec[%d] = \n",0);
        Cvec.print(outfile);
        Cvec.read(1,0);
        outfile->Printf(" Cvec[%d] = \n",1);
        Cvec.print(outfile);
        Cvec.buf_unlock();
        */

        /* Construct C(i) = alpha(i)*C(i) + alpha(i-1)*C(i-1) */
        mitrush_update(Cvec,Sigma,norm,acur,alast,buffer1,buffer2,curr,last);

        /* put H(ii) = E(i) */
        E = evals2x2[0];

        } /* end Mitrushenkov Step */

      else {
         Cvec.buf_unlock();
         Sigma.buf_unlock();
         if (Parameters_->diag_method==METHOD_MITRUSHENKOV &&
             diag_method==METHOD_MITRUSHENKOV)
           {
            diag_method = METHOD_OLSEN;
            last = curr;
           }
         }

      /*
      ** get y = C(0) * (Hd - E)^-1 * sigma(0)
      ** and x = C(0) * (Hd - E)^-1 * C(0)
      */
      olsen_iter_xy(Cvec,Sigma,Hd,&x,&y,buffer1,buffer2,E,curr,1,alpha,
                    alplist, betlist);

      if (print_ > 3) {
        outfile->Printf( "Straight x = %12.6lf\n", x);
        outfile->Printf( "Straight y = %12.6lf\n", y);
        }

      if (Parameters_->precon >= PRECON_GEN_DAVIDSON && H0block_->size) {
        detH0 = H0block_calc(E);
        /* outfile->Printf("detH0 = %d\n", detH0);
         */
        if (detH0 == 0) {
          outfile->Printf( "H0block inverse is nearly singular:");
          outfile->Printf(" initiating DAVIDSON preconditioner\n");
          Parameters_->precon = PRECON_DAVIDSON;
          }
        if (Parameters_->precon >= PRECON_GEN_DAVIDSON) {
          H0block_xy(&x,&y,E);
          if (print_ > 3) {
            outfile->Printf( "x = %12.6lf\n", x);
            outfile->Printf( "y = %12.6lf\n", y);
            }
          }
        }

      E_est = y / x;
      if (print_ > 2) {
        outfile->Printf( "E_est = %12.6lf E = %12.6lf\n",
                E_est+edrc+enuc, E+enuc);
        /* outfile->Printf( "x = %lf  y = %lf\n",x,y); */
        }

      /* calculate delta_C and C(next) */
      olsen_update(Cvec,Sigma,Hd,E,E_est,&norm,&c1norm,&S,buffer1,buffer2,
                   curr,last,"outfile",iter, alplist, betlist);

      norm = sqrt(1.0 / norm);
      Cvec.buf_lock(buffer1);
      Cvec.read(last, 0);
      S *= norm;
      if (Parameters_->precon >= PRECON_GEN_DAVIDSON &&
          diag_method==METHOD_MITRUSHENKOV && S<S_MAX)
        Cvec.symnorm(norm,CI_VEC,FALSE);
      else Cvec.symnorm(norm,CI_VEC,TRUE);

     /*
      Cvec.buf_lock(buffer1);
      Cvec.buf_unlock();
      Cvec.read(0,0);
      outfile->Printf(" Cvec[0] = \n");
      Cvec.print(outfile);
      Cvec.read(1,0);
      outfile->Printf(" Cvec[1] = \n");
      Cvec.print(outfile);
     */

      if (diag_method==METHOD_MITRUSHENKOV) {
         curr = !curr;
         last = !last;
        }

      if ((fabs(E - E_last) < conv_e && c1norm < conv_rms) || iter >=maxiter) {
        outfile->Printf( "Iter %2d  ROOT 1 ECI = %14.9lf", iter, E + enuc);
        outfile->Printf( "    Delta_E %10.3E   Delta_C %10.3E %c\n"
          ,E-E_last,c1norm,(fabs(E - E_last) < conv_e && c1norm < conv_rms)
          ? 'c' : ' ');
        evals[0] = E;
        free_matrix(H2x2,2);
        free(evals2x2);
        free_matrix(evecs2x2,2);
        outfile->Printf( "\n\n* ROOT 1 CI total energy = %19.15lf\n", E + enuc);

        Cvec.max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx,
           mi_ibidx, mi_coeff, Parameters_->neg_only);
        print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
           mi_coeff);
        free(mi_iac);  free(mi_ibc); free(mi_iaidx);  free(mi_ibidx);
        free(mi_coeff);
        Parameters_->diag_h_converged = true;
        return;
        }

      Cvec.buf_unlock();
      Sigma.buf_unlock();
      outfile->Printf("Iter %2d  ROOT 1 ECI = %14.9lf", iter, enuc+E);
      outfile->Printf("    Delta_E %10.3E   Delta_C %10.3E\n",E-E_last,c1norm);

      iter++;
      Parameters_->diag_iters_taken = iter;
      E_last = E;
      if (Parameters_->calc_ssq && Parameters_->icore==1)
        Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, 0);
      } /* end while (1) */

}

/*
** olsen_update()
**
** Update the CI vector according to the Olsen method.
** Try to assume storage of only one C, one S, and Hd on disk.  Let
** the mitrush_update() function assume two C and two S.  No, can't
** assume that or we'll never get Mitrush started.  Just fix it so that
** Cn can be identical to Cc.
** The update formula is c0' = c0 + (Ho-Eo)^-1(E_est*Co - sigma_o)
** this is the same formula as the davidson procedure except lambda_i
** has been substituted by E_est
*/
void CIWavefunction::olsen_update(CIvect &C, CIvect &S, CIvect &Hd, double E, double E_est,
      double *norm, double *c1norm, double *ovrlap, double *buffer1,
      double *buffer2,
      int curr, int next, std::string out, int iter, struct stringwr **alplist,
      struct stringwr **betlist)
{

   int buf;
   double nx=0.0, ox=0.0, tmp1, tmp2, normc1=0.0, tmpnorm=0.0;
   double rnorm=0.0, rnormtmp=0.0;

   for (buf=0; buf<C.buf_per_vect_; buf++) {
      tmp1 = 0.0;
      tmp2 = 0.0;
      C.buf_lock(buffer1);
      S.buf_lock(buffer2);
      C.read(curr, buf);
      S.read(curr, buf);
      /* C = E_est * C - S, C is buffer1*/
      xeaxmy(buffer1, buffer2, E_est, C.buf_size_[buf]);
      C.buf_unlock();
      S.buf_unlock();
      Hd.buf_lock(buffer2);
      if (Parameters_->hd_otf == FALSE) Hd.read(0,buf);
      else Hd.diag_mat_els_otf(alplist, betlist, CalcInfo_->onel_ints->pointer(),
           CalcInfo_->twoel_ints->pointer(), CalcInfo_->edrc, CalcInfo_->num_alp_expl,
           CalcInfo_->num_bet_expl, CalcInfo_->nmo, buf, Parameters_->hd_ave);
      /* Check norm of residual vector i.e. before preconditioning */
      dot_arr(buffer1, buffer1, C.buf_size_[buf], &rnormtmp);
      /* C = C/(Hd - E) */
      buf_ols_denom(buffer1, buffer2, E, S.buf_size_[buf]);
      /* buffer1 is now equal to C^1, i.e. the correction to C_i
      ** without the H0block correction
      */
      Hd.buf_unlock();
      /* C_new = C_i + C^1 */
      C.buf_lock(buffer2);
      C.read(curr, buf);
      buf_ols_updt(buffer1,buffer2,&tmp1,&tmp2,&tmpnorm,C.buf_size_[buf]);
      if (Parameters_->precon >= PRECON_GEN_DAVIDSON)
        C.h0block_buf_ols(&tmp1,&tmp2,&tmpnorm,E_est);
      if (C.buf_offdiag_[buf]) {
         tmp1 *= 2.0;
         tmp2 *= 2.0;
         tmpnorm *= 2.0;
         rnormtmp *= 2.0;
         }
      normc1 += tmpnorm;
      nx += tmp1;
      ox += tmp2;
      rnorm += rnormtmp;
      C.write(next, buf);
      C.buf_unlock();
 }

   *norm = nx;
   /* outfile->Printf("\n ovrlap(ox) = %20.16f\n", ox); */
   *ovrlap = ox;
   if (normc1 <= 1.0E-13) {
     outfile->Printf("Norm of correction vector = %5.4e\n", normc1);
     outfile->Printf("This may cause numerical errors which would" \
         " deteriorate the diagonalization procedure.\n");
     }
   *c1norm = sqrt(rnorm);
   normc1 = sqrt(normc1);
}

/*
** olsen_iter_xy(): This function evaluates the quantities x and y defined by
**    x = C^(i) * (Hd - E)^-1 * C^(i)
**    y = C^(i) * (Hd - E)^-1 * sigma^(i)
**    The diagonal elements of H (Hd) are used in this routine.
**    These quantities are subtracted out and substituted in later
**    routines if the gen_davidson or h0block_inv preconditioners
**    are used. A very intelligent approach originally implimented
**    by a fellow Jedi.
**
** Parameters:
**   C       =  reference of current iteration's ci vector
**   S       =  reference of current iteration's sigma vector
**   Hd      =  reference of vector of diagonal elements of H
**   x       =  pointer to double to hold x
**   y       =  pointer to double to hold y
**   buffer1 = pointer to first I/O buffer
**   buffer2 = pointer to second I/O buffer
**   E       =  current iteration's energy
**   curvect = current vector
**   L       = number of vectors in b or sigma files
**   alplist = alpha string list for use with diag energies on the fly
**   betlist = beta string list for use with diag energies on the fly
**
** Returns: none
*/
void CIWavefunction::olsen_iter_xy(CIvect &C, CIvect &S, CIvect &Hd, double *x, double *y,
      double *buffer1, double *buffer2, double E, int curvect, int L,
      double **alpha, struct stringwr **alplist, struct stringwr **betlist)
{

   int buf, i,j;
   double tx = 0.0, ty = 0.0, tmpy = 0.0;
   double *sigma0b1, *sigma0b2;
   *x = 0.0;
   *y = 0.0;

   Hd.buf_lock(buffer2);
   if (Parameters_->diag_method==METHOD_DAVIDSON_LIU_SEM) {
     sigma0b1 = init_array(H0block_->size);
     sigma0b2 = init_array(H0block_->size);
     }
   for (buf=0; buf<C.buf_per_vect_; buf++) {
      tx = ty = 0.0;
      C.buf_lock(buffer1);
      C.read(curvect,buf);
      if (Parameters_->diag_method==METHOD_DAVIDSON_LIU_SEM)
        C.h0block_gather_vec(CI_VEC);
      if (Parameters_->hd_otf == FALSE) Hd.read(0,buf);
      else Hd.diag_mat_els_otf(alplist, betlist, CalcInfo_->onel_ints->pointer(),
           CalcInfo_->twoel_ints->pointer(), CalcInfo_->edrc, CalcInfo_->num_alp_expl,
           CalcInfo_->num_bet_expl, CalcInfo_->nmo, buf, Parameters_->hd_ave);
      tx = buf_xy1(buffer1, buffer2, E, Hd.buf_size_[buf]);
      /* buffer2 = Hd * Ci */
      C.buf_unlock();
      S.buf_lock(buffer1);
      if (Parameters_->diag_method <= METHOD_MITRUSHENKOV) {
        /* Olsen and Mitrushenkov iterators */
        S.read(curvect,buf);
        dot_arr(buffer1, buffer2, C.buf_size_[buf], &ty);
        }
      else { /* Dot buffer2 with all Sigma vectors on disk */
        for (i=0; i<L; i++) {
           S.read(i,buf);
           dot_arr(buffer1, buffer2, C.buf_size_[buf], &tmpy);
           ty += tmpy * alpha[i][curvect];
           zero_arr(sigma0b1,H0block_->size);
           S.h0block_gather_multivec(sigma0b1);
           for (j=0; j<H0block_->size; j++)
              sigma0b2[j] += alpha[i][curvect] * sigma0b1[j];
           }

       }
      if (C.buf_offdiag_[buf]) {
         *x += 2.0 * tx;
         *y += 2.0 * ty;
         }
      else {
         *x += tx;
         *y += ty;
         }
      S.buf_unlock();
      }

   Hd.buf_unlock();
   if (Parameters_->diag_method==METHOD_DAVIDSON_LIU_SEM) {
     for (j=0; j<H0block_->size; j++)
        H0block_->s0b[j] = sigma0b2[j];
     free(sigma0b1);
     free(sigma0b2);
     }

}

/*
** mitrush_update()
** Perform the Mitrushenkov update.  New version 3/96
**
*/
void CIWavefunction::mitrush_update(CIvect &C, CIvect &S, double norm, double acur,
   double alast, double *buffer1, double *buffer2, int curr, int next)
{
   int i, j, k, buf, blk, al, bl;
   double phase, tval;

   if (!Parameters_->Ms0) phase = 1.0;
   else phase = ((int) Parameters_->S % 2) ? -1.0 : 1.0;

   for (buf=0; buf<C.buf_per_vect_; buf++) {
      C.buf_lock(buffer1);
      C.read(curr, buf);
      C.buf_unlock();
      C.buf_lock(buffer2);
      C.read(next, buf);
      xeaxpby(buffer2, buffer1, alast, acur, C.buf_size_[buf]);
      /* buffer2 is the new C(i) vector constructed from
      ** alpha(i-1)*C(i-1) + alpha(i)*C(i).
      ** However, this vector is not normalized or symmetrized yet.
      ** Still need to construct the olsen_update of this new
      ** C(i) vector.
      */
      C.write(curr, buf);
      C.buf_unlock();
      }
   C.buf_lock(buffer1);
   C.read(curr,0);
   C.symnorm(norm,0,1);
   C.buf_unlock();

   for (buf=0; buf<S.buf_per_vect_; buf++) {
      S.buf_lock(buffer1);
      S.read(curr, buf);
      S.buf_unlock();
      S.buf_lock(buffer2);
      S.read(next, buf);
      xeaxpby(buffer2, buffer1, alast, acur, S.buf_size_[buf]);
      S.write(curr, buf);
      S.buf_unlock();
      }
   S.buf_lock(buffer1);
   S.read(curr,0);
   S.symnorm(norm,1,1);
   S.buf_unlock();
}







}} // namespace psi::detci
