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

#define EXTERN 
/* #define DEBUG */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "globals.h"
#include "civect.h"
#include "ci_tol.h"

namespace psi { namespace detci {

extern int H0block_calc(double E);
extern int H0block_coupling_calc(double E, struct stringwr *alplist,
   struct stringwr *betlist);
extern void H0block_xy(double *x, double *y, double E);
extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode, 
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG, 
   struct stringwr **alplist, struct stringwr **betlist,
   FILE *outfile);


/* #define DEBUG */
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
void mitrush_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
      **betlist, int nroots, double *evals, double conv_rms, double conv_e, 
      double enuc, double efzc, int maxiter, int maxnvect, FILE *outfile, 
      int print_lvl)
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

   Cvec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit,
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
   Sigma.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore,Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_s_tmp_units, Parameters.first_s_tmp_unit,
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

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
   last = (Parameters.diag_method == METHOD_MITRUSHENKOV) ? 1 : 0;
   diag_method = Parameters.diag_method;

   /* set buffer pointers */

   buffer1 = *(Hd.blockptr(0));
   buffer2 = Hd.buf_malloc();
   Hd.buf_unlock();

   /* get some of the stuff from CalcInfo for easier access */

   num_alp_str = CalcInfo.num_alp_str;
   num_bet_str = CalcInfo.num_bet_str;
   if (Parameters.fci) oei = CalcInfo.tf_onel_ints;
   else oei = CalcInfo.gmat[0];
   tei = CalcInfo.twoel_ints;

   /* small arrays to hold most important config information */

   mi_iac = init_int_array(Parameters.nprint);
   mi_ibc = init_int_array(Parameters.nprint);
   mi_iaidx = init_int_array(Parameters.nprint);
   mi_ibidx = init_int_array(Parameters.nprint);
   mi_coeff = init_array(Parameters.nprint);


   /* stuff for the 2x2 Davidson procedure */

   H2x2 = init_matrix(2,2);
   evals2x2 = init_array(2);
   evecs2x2 = init_matrix(2,2);
   alpha = init_matrix(1,1); 
    /* not used but is necessary for call to olsen_iter_xy */ 

   /* setup initial guess vector */

   if (Parameters.restart) {
     fprintf(outfile,"\nAttempting Restart with 1 vector\n");
     if ((i=Cvec.read_num_vecs())< 1) {
       fprintf(outfile, "CI vector file contains %d vectors, need 1.\n", i);
       exit(0);
     }
     fflush(outfile);
     Cvec.buf_lock(buffer1);
     Cvec.read(0, 0);
     tval = Cvec * Cvec;
     if ((tval - 1.0) > ZERO) fprintf(outfile,"CI vector may be corrupted."
     " Attempting to correct by renormalizing.\n"); 
     Cvec.symnorm(tval,CI_VEC,TRUE);
     }
   else if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT && nroots == 1 &&
      Parameters.num_init_vecs == 1) { /* use unit vector */
      tval = 1.0;
      Cvec.buf_lock(buffer1);
      Cvec.init_vals(0, 1, &(CalcInfo.ref_alp_list), &(CalcInfo.ref_alp_rel),
                     &(CalcInfo.ref_bet_list), &(CalcInfo.ref_bet_rel), 
                     H0block.blknum, &tval);
      Cvec.write_num_vecs(1);
      k = 1;
      Cvec.read(0, 0);
      Cvec.symnorm(1.0,CI_VEC,TRUE);
      Sigma.set_zero_blocks_all();
      }

   else { /* use H0BLOCK eigenvector guess */
      if (Parameters.precon == PRECON_GEN_DAVIDSON) L = H0block.size;
      else L = H0block.guess_size;
      sm_evals = init_array(L);
      /* need to fill out sm_evecs into b (pad w/ 0's) */
      fprintf(outfile, "Using %d initial trial vectors\n",
         Parameters.num_init_vecs);

      Cvec.buf_lock(buffer1);
      for (i=0,k=0; i<L && k<1; i++) {

         /* check sm_evecs[i] to see if it has the correct spin symmetry */
         /* if Ms=0 ... */
         for (j=0,tmpi=0; j<L && !tmpi && Parameters.Ms0; j++) {
            l = H0block.pair[j];
            if (l == -1) { tmpi = 1; break; }
            tval = H0block.H0b_diag[l][i];
            if ((int) Parameters.S % 2) tval = -tval;
            if (H0block.H0b_diag[j][i] - tval > 1.0E-8) tmpi = 1;
            }
         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = H0block.H0b_diag[j][i];

         Cvec.init_vals(k, L, H0block.alplist, H0block.alpidx,
            H0block.betlist, H0block.betidx, H0block.blknum, sm_evals);
         Cvec.write_num_vecs(1);
         k++;
         /* MLL added 6-15-98 */
         Cvec.read(0, 0);
         Cvec.symnorm(1.0,CI_VEC,TRUE);
         }

      free(sm_evals);
      Sigma.set_zero_blocks_all();
      }

   if (print_lvl > 4) {
      fprintf(outfile, "\nC(0) vector = \n");
      Cvec.print(outfile);
      }

   Sigma.buf_lock(buffer2);
   Cvec.read(0, 0);
   sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, 0);
   Sigma.write_num_vecs(1);

   Cvec.copy_zero_blocks(Sigma);

   if (print_lvl > 4) {
      fprintf(outfile, "\nSigma vector\n");
      Sigma.print(outfile);
      fflush(outfile);
      }

   /* get H00 */
   E = Cvec * Sigma; 
   E += efzc;
   E_last = CalcInfo.escf - CalcInfo.enuc;

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

   if (Parameters.print_lvl > 3) {
     fprintf(outfile, "Straight x = %12.6lf\n", x);
     fprintf(outfile, "Straight y = %12.6lf\n", y);
     }

   if (Parameters.precon >= PRECON_GEN_DAVIDSON && H0block.size) {
       detH0 = 1;
       /* 
       detH0 = H0block_calc(E); 
       */ 
       fflush(outfile); 
       if (detH0 == 0) {
         fprintf(outfile, "H0block inverse is nearly nonsingular:");
         fprintf(outfile," initiating DAVIDSON preconditioner\n");
         fflush(outfile);
         Parameters.precon = PRECON_DAVIDSON;
         }
       if (Parameters.precon >= PRECON_GEN_DAVIDSON) {
         /* 
         H0block_xy(&x,&y,E); 
         */
         if (Parameters.print_lvl > 3) {
            fprintf(outfile, "x = %12.6lf\n", x);
            fprintf(outfile, "y = %12.6lf\n", y);
            }
         }
     }

   E_est = y / x;  /* should I add fzc here? */
   if (Parameters.print_lvl > 2) {
      fprintf(outfile, "E_est = %12.6lf E-efzc = %12.6lf E = %12.6lf\n", 
           E_est,E-efzc,E);
       fprintf(outfile, "x = %lf  y = %lf\n",x,y);
       }
    /* calculate delta_C and C(1) */
   olsen_update(Cvec, Sigma, Hd, E, E_est, &norm, &c1norm, &S, buffer1, buffer2,
                 curr, last, outfile, iter, alplist, betlist);
 
   norm = sqrt(1.0 / norm);
   Cvec.buf_lock(buffer1);
   Cvec.read(last, 0);
   Cvec.symnorm(norm,CI_VEC,TRUE); 
   S *= norm;
   Cvec.buf_unlock();
   if (Parameters.calc_ssq && Parameters.icore==1)
     Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, 0);

  /*
   Cvec.buf_lock(buffer1);
   Cvec.read(0,0);
   fprintf(outfile," Cvec[0] = \n");
   Cvec.print(outfile);
   Cvec.read(1,0);
   fprintf(outfile," Cvec[1] = \n");
   Cvec.print(outfile);
   Cvec.buf_unlock();
  */

   /* S is the overlap of the C_(i-1) and C_i */

   fprintf(outfile,"Iter  0  ROOT 1 ECI = %14.9lf", enuc + E);
   fprintf(outfile,"    Delta_E %10.3E   Delta_C %10.3E\n", E - E_last, c1norm);
   fflush(outfile);

   iter = 1;
   if (Parameters.diag_method == METHOD_MITRUSHENKOV) {
      curr = 1; last = 0;
      }

   while (1) {

      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);
      Cvec.read(curr,0);
      /* chknorm = Cvec.checknorm();
      fprintf(outfile,"Norm of curr CI vect = %lf\n",chknorm);
     */ 
      if (print_lvl > 4) {
         fprintf(outfile, "\nC(%2d) vector (symm'd norm'd)\n", iter) ;
         Cvec.print(outfile);
         fflush(outfile);
         }

      Cvec.read(curr,0);
      /* Sigma.read(curr,0);
      */
      sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, curr);
      Cvec.copy_zero_blocks(Sigma);
      Cvec.read(curr,0);
      if (print_lvl > 4) {
         fprintf(outfile,"\nC(%2d) vector (symm'd norm'd) second time\n",iter);
         Cvec.print(outfile);
         fflush(outfile);
         }

      if (print_lvl > 4) {
         Sigma.read(curr,0);
         fprintf(outfile, "\n curr = %d\n", curr);
         fprintf(outfile, "\nSigma(%2d) vector\n", iter);
         Sigma.print(outfile);
         fflush(outfile);
         }

      Cvec.read(curr,0);
      Sigma.read(curr,0);
      /* calculate H(ii) */
      E_curr = Cvec * Sigma;      
      E_curr += efzc;
      E = E_curr;

      /* check for convergence and exit if reached */
      if (iter == maxiter)
         fprintf(outfile, "Maximum number of iterations reached!\n");

      /* if the 2x2 matrix diagonalization can be done, take Mitrush Step */

      if (Parameters.diag_method==METHOD_MITRUSHENKOV && 
          diag_method==METHOD_MITRUSHENKOV && S < S_MAX && 
          fabs(E_last-E_curr) > MITRUSH_E_DIFF_MIN) {
        fprintf(outfile, "Taking Mitrushenkov step (S =%10.6lf <%10.6lf)\n",
                 S, S_MAX);
         /* calculate H(i,i-1) = H(i-1,i) */
        Cvec.read(last, 0);         
        E12 = Cvec * Sigma;
        E12 += efzc * S;

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
        fprintf(outfile,"alpha_0 = %lf  alpha_1 = %lf\n",
                evecs2x2[0][0],evecs2x2[0][1]);
        Cvec.buf_lock(buffer1);
        Cvec.read(0,0);
        fprintf(outfile," Cvec[%d] = \n",0);
        Cvec.print(outfile);
        Cvec.read(1,0);
        fprintf(outfile," Cvec[%d] = \n",1);
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
         if (Parameters.diag_method==METHOD_MITRUSHENKOV && 
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

      if (Parameters.print_lvl > 3) {
        fprintf(outfile, "Straight x = %12.6lf\n", x);
        fprintf(outfile, "Straight y = %12.6lf\n", y);
        }

      if (Parameters.precon >= PRECON_GEN_DAVIDSON && H0block.size) {
        detH0 = H0block_calc(E);
        /* fprintf(outfile,"detH0 = %d\n", detH0);
        fflush(outfile); */
        if (detH0 == 0) {
          fprintf(outfile, "H0block inverse is nearly singular:");
          fprintf(outfile," initiating DAVIDSON preconditioner\n");
          Parameters.precon = PRECON_DAVIDSON;
          }
        if (Parameters.precon >= PRECON_GEN_DAVIDSON) {
          H0block_xy(&x,&y,E);
          if (Parameters.print_lvl > 3) {
            fprintf(outfile, "x = %12.6lf\n", x);
            fprintf(outfile, "y = %12.6lf\n", y);
            }
          }
        }

      E_est = y / x;
      if (Parameters.print_lvl > 2) {
        fprintf(outfile, "E_est = %12.6lf E = %12.6lf\n", 
                E_est+efzc+enuc, E+enuc);
        /* fprintf(outfile, "x = %lf  y = %lf\n",x,y); */
        }

      /* calculate delta_C and C(next) */
      olsen_update(Cvec,Sigma,Hd,E,E_est,&norm,&c1norm,&S,buffer1,buffer2,
                   curr,last,outfile,iter, alplist, betlist);

      norm = sqrt(1.0 / norm);
      Cvec.buf_lock(buffer1);
      Cvec.read(last, 0);
      S *= norm; 
      if (Parameters.precon >= PRECON_GEN_DAVIDSON && 
          diag_method==METHOD_MITRUSHENKOV && S<S_MAX) 
        Cvec.symnorm(norm,CI_VEC,FALSE); 
      else Cvec.symnorm(norm,CI_VEC,TRUE);

     /*
      Cvec.buf_lock(buffer1);
      Cvec.buf_unlock();
      Cvec.read(0,0);
      fprintf(outfile," Cvec[0] = \n");
      Cvec.print(outfile);
      Cvec.read(1,0);
      fprintf(outfile," Cvec[1] = \n");
      Cvec.print(outfile);
     */

      if (diag_method==METHOD_MITRUSHENKOV) {
         curr = !curr;
         last = !last;
        }

      if ((fabs(E - E_last) < conv_e && c1norm < conv_rms) || iter >=maxiter) {
        fprintf(outfile, "Iter %2d  ROOT 1 ECI = %14.9lf", iter, E + enuc);
        fprintf(outfile, "    Delta_E %10.3E   Delta_C %10.3E %c\n"
          ,E-E_last,c1norm,(fabs(E - E_last) < conv_e && c1norm < conv_rms)
          ? 'c' : ' ');
        evals[0] = E;
        free_matrix(H2x2,2);
        free(evals2x2);
        free_matrix(evecs2x2,2);
        fprintf(outfile, "\n\n* ROOT 1 CI total energy = %19.15lf\n", E + enuc);

        Cvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx,
           mi_ibidx, mi_coeff, Parameters.neg_only);
        print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
           mi_coeff, AlphaG, BetaG, alplist, betlist, outfile);
        free(mi_iac);  free(mi_ibc); free(mi_iaidx);  free(mi_ibidx);
        free(mi_coeff);
        return;
        }

      Cvec.buf_unlock();
      Sigma.buf_unlock();
      fprintf(outfile,"Iter %2d  ROOT 1 ECI = %14.9lf", iter, enuc+E);
      fprintf(outfile,"    Delta_E %10.3E   Delta_C %10.3E\n",E-E_last,c1norm);
      fflush(outfile);
      iter++;
      E_last = E;
      if (Parameters.calc_ssq && Parameters.icore==1)
        Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, 0);
      } /* end while (1) */

}

}} // namespace psi::detci

