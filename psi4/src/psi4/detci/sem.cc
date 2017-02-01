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
    \ingroup DETCI
    \brief Enter brief description of file here
*/
/*
** Simultaneous Expansion Method Iterator
**
** C. David Sherrill
** August 29, 1995
**
** Last modified by MLL on 25 November 1997
*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libqt/slaterdset.h"
#include "psi4/libmints/vector.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"
#include "psi4/physconst.h"

namespace psi { namespace detci {


#define MALPHA_TOLERANCE 1E-15

void CIWavefunction::sem_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
      **betlist, double *evals, double conv_e,
      double conv_rms, double enuc, double edrc,
      int nroots, int maxiter, int maxnvect)
{
   int i, j, k, l, ij, I, L, L2=0, L3=0, tmpi, detH0;
   unsigned long det1, N;
   int num_alp_str, num_bet_str, Llast;
   int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx, *root_converged;
   int *Lvec, *did_root, num_root_converged;
   double *mi_coeff, *clpse_norm, **clpse_dot, **tmpmat;
   double *oei, *tei, **G, ***alpha, **lambda, ****m_alpha, ***m_lambda;
   int sm_tridim;
   double *sm_mat, *sm_evals, **sm_evecs;
   int iter = 0, converged = 0;
   int iter2 = 0; /* iterations since last collapse */
   double tval, tval2, *lastroot, *dvecnorm, *buffer1, *buffer2;
   double ***M, **sigma_overlap, Mtmp, **cmp_cncoe, **tr_cmp_cncoe;
   int *lse_do_arr, lse_do = 0, collapse_num = 0, iter_tmp = 0;
   int form_M = 0, tmpval;
   int last_lse_collapse_num = -Parameters_->lse_collapse;
   double *x, *y, tmpx, tmpy;
   double lse_tolerance, *renorm_c, *E_est, ovlpmax=0.0;
   double cknorm, tvalmatt=0.0, tmp; /* Add by CDS for debugging purposes */
   int errcod;
   std::string str;
   bool dvec_read_fail = false;

   CIvect Cvec(Parameters_->icore, maxnvect, 1,
          Parameters_->c_filenum, CIblks_, CalcInfo_, Parameters_,
          H0block_, false);
   CIvect Cvec2(Parameters_->icore, maxnvect, 1,
          Parameters_->c_filenum, CIblks_, CalcInfo_, Parameters_,
          H0block_, false);
   CIvect Sigma(Parameters_->icore, maxnvect, 1,
                Parameters_->s_filenum, CIblks_, CalcInfo_, Parameters_,
                H0block_, false);
   CIvect Sigma2(Parameters_->icore, maxnvect, 1,
                 Parameters_->s_filenum, CIblks_, CalcInfo_, Parameters_,
                 H0block_, false);

   int dvec_file, dvec2_file;
   if (!Parameters_->nodfile){
     dvec_file = Parameters_->d_filenum;
     dvec2_file = Parameters_->d_filenum;
   }
   else{
     dvec_file = Parameters_->c_filenum;
     dvec2_file = Parameters_->s_filenum;
   }

   CIvect Dvec(Parameters_->icore, maxnvect, 1,
               dvec_file, CIblks_, CalcInfo_, Parameters_,
               H0block_, false);
   CIvect Dvec2(Parameters_->icore, maxnvect, 1,
                dvec2_file, CIblks_, CalcInfo_, Parameters_,
                H0block_, false);

   /* open the files: some of these CIvectors are logical vectors that
      point to the same files ... don't need to repeat the file opens
      for those
   */
   bool open_old;
   if (Parameters_->restart) open_old = true;
   else open_old = false;
   Cvec.init_io_files(open_old);
   Sigma.init_io_files(open_old);

   if (Parameters_->guess_vector == PARM_GUESS_VEC_DFILE) open_old = true;
   else open_old = false;
   Dvec.init_io_files(open_old);


   /* set up the vector pointers/info */
   if (Cvec.read_new_first_buf() == -1) Cvec.write_new_first_buf();
   if (Sigma.read_new_first_buf() == -1) Sigma.write_new_first_buf();
   /* should not need to have first_buf info in D file, it never changes
    * unless it is only a logical D file, which isn't controlled here
    * anyway. */
   if (Cvec.read_num_vecs() == -1) Cvec.write_num_vecs(0);
   if (Sigma.read_num_vecs() == -1) Sigma.write_num_vecs(0);
   if (Dvec.read_num_vecs() == -1) Dvec.write_num_vecs(0);

   /* allocate memory */
   Dvec.h0block_buf_init();
   buffer1 = *(Hd.blockptr(0));
   buffer2 = Hd.buf_malloc();
   Hd.buf_unlock();


   /* get some of the stuff from CalcInfo for easier access */
   num_alp_str = CalcInfo_->num_alp_str;
   num_bet_str = CalcInfo_->num_bet_str;
   if (Parameters_->fci) oei = CalcInfo_->tf_onel_ints->pointer();
   else oei = CalcInfo_->gmat->pointer();
   tei = CalcInfo_->twoel_ints->pointer();

   lastroot = init_array(nroots);
   dvecnorm = init_array(nroots);
   root_converged = init_int_array(nroots);
   did_root = init_int_array(nroots);
   clpse_norm = init_array(maxnvect);
   clpse_dot = init_matrix(maxnvect, maxnvect);
   tmpmat = init_matrix(maxnvect, nroots);
   lse_do_arr = init_int_array(nroots);
   renorm_c = init_array(nroots);
   x = init_array(nroots);
   y = init_array(nroots);
   E_est = init_array(nroots);

   /* small arrays to hold most important config information */
   mi_iac = init_int_array(Parameters_->nprint);
   mi_ibc = init_int_array(Parameters_->nprint);
   mi_iaidx = init_int_array(Parameters_->nprint);
   mi_ibidx = init_int_array(Parameters_->nprint);
   mi_coeff = init_array(Parameters_->nprint);

   G = init_matrix(maxnvect, maxnvect);
   cmp_cncoe = init_matrix(maxnvect, maxnvect);
   tr_cmp_cncoe = init_matrix(maxnvect, maxnvect);
   sigma_overlap = init_matrix(maxnvect, maxnvect);
   Lvec = init_int_array(maxnvect);
   lambda = init_matrix(maxnvect, maxnvect);

   m_lambda = (double ***) malloc (sizeof(double **) * maxnvect);
   for (i=0; i<maxnvect; i++)
      m_lambda[i] = init_matrix(nroots, maxnvect);

   M = (double ***) malloc (sizeof(double **) * nroots);
   for (i=0; i<nroots; i++)
         M[i] = init_matrix(maxnvect, maxnvect);

   alpha = (double ***) malloc (sizeof(double **) * maxnvect);
   for (i=0; i<maxnvect; i++) {
      alpha[i] = init_matrix(maxnvect, maxnvect);
      }

   m_alpha = (double ****) malloc (sizeof(double ***) * maxiter);
   for (i=0; i<maxiter; i++) {
      m_alpha[i] = (double ***) malloc (sizeof(double **) * nroots);
      for (j=0; j<nroots; j++)
         m_alpha[i][j] = init_matrix(maxnvect, maxnvect);
      }

   if (Parameters_->lse) lse_tolerance = Parameters_->lse_tolerance;

   if (Parameters_->nodfile == FALSE &&
     Parameters_->guess_vector == PARM_GUESS_VEC_DFILE) {
     if ((i = Dvec.read_num_vecs()) != nroots) {
       if (print_) outfile->Printf( "    D file contains %d not %d vectors.  Trying another guess.\n", i, nroots);
       dvec_read_fail = true;
       /*
       if (Parameters_->h0blocksize == 0) {
         Parameters_->guess_vector = PARM_GUESS_VEC_UNIT;
         outfile->Printf( "unit vector guess.\n");
       }
       else {
         Parameters_->guess_vector = PARM_GUESS_VEC_H0_BLOCK;
         outfile->Printf( "H0block guess.\n");
       }
       */
     }
   }

   if (Parameters_->restart) {  /* restart option! */
      // L = Parameters_->restart_vecs;
      L = Cvec.read_num_vecs();
      i = Sigma.read_num_vecs();
      if (i != L) {
        outfile->Printf( "    %d C vectors and %d Sigma vectors.\n", i, L);
        if (i < L) {
          L = i;
          Cvec.write_num_vecs(L);
        }
        if (print_) outfile->Printf( "    Using %d vectors \n", L);
      }
      if (L < nroots) {
        str = "Restart failed...  ";
        str += std::to_string( L) ;
        str += " vectors for ";
        str += std::to_string( nroots) ;
        str += " roots";
        throw PsiException(str,__FILE__,__LINE__);
      }

      if (print_) outfile->Printf( "\n    Attempting Restart with %d vectors\n\n", L);

   /* open detci.dat and write file_offset and file_number array out to
      detci.dat */

      //Cvec.reset_detfile(CI_VEC);
      //Cvec2.reset_detfile(CI_VEC);
      i = Cvec.read_new_first_buf();
      Cvec.set_new_first_buf(i);
      Cvec2.set_new_first_buf(i);
      //Sigma.reset_detfile(SIGMA_VEC);
      //Sigma2.reset_detfile(SIGMA_VEC);
      j = Sigma.read_new_first_buf();
      Sigma.set_new_first_buf(j);
      Sigma2.set_new_first_buf(j);
      /* the first buffer of D file should not change unless there
       * is only a logical D file not a physical D file */
      if (Parameters_->nodfile) {
        //Dvec.reset_detfile(CI_VEC);
    Dvec.set_new_first_buf(i);
        //Dvec2.reset_detfile(SIGMA_VEC);
        Dvec2.set_new_first_buf(j);
    /* I hope I don't double-reset anything w/ next 2 lines */
        Dvec.restart_reord_fp(maxnvect-1);
        Dvec2.restart_reord_fp(maxnvect-1);
        }

      // The next few lines might help debugging restarts
      // Cvec.civect_psio_debug();
      // Sigma.civect_psio_debug();
      // Dvec.civect_psio_debug();
      //

      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);

      for (i=0; i<L; i++) {
         Sigma.read(i, 0);
         if (print_ > 4) {
            outfile->Printf( "    Sigma[%d] =\n", i);
            Sigma.print();
            }
         for (j=0; j<=i; j++) {
            Cvec.read(j, 0);
            if (print_ > 4) {
               outfile->Printf( "    C[%d] =\n", j);
               Cvec.print();
               }
            G[j][i] = G[i][j] = Cvec * Sigma;
            }
         }
      Cvec.buf_unlock();
      Sigma.buf_unlock();

      if (print_ > 3) {
         outfile->Printf( "\n    G matrix (%2d) = \n", iter);
         print_mat(G, L, L, "outfile");
         }

      /* solve the L x L eigenvalue problem G a = lambda a for M roots */
      sq_rsp(L, L, G, lambda[iter2], 1, alpha[iter2], 1.0E-14);
      if (print_ > 4) {
         outfile->Printf( "\n     G eigenvectors and eigenvalues:\n");
         eivout(alpha[iter2], lambda[iter2], L, L, "outfile");
         }

      /* loop over roots and write out the required number of vects */
      if (nroots + L > maxnvect) {
         throw PsiException("    Error: Can't do restart if maxnvect < nroots + L",__FILE__,__LINE__);
         }


      /* gather C and Sigma */
      Cvec.buf_lock(buffer1);
      Dvec.buf_lock(buffer2);
      for (i=0; i<nroots; i++) {
         Dvec.gather(i, L, i, alpha[iter2], Cvec);
         }
      for (i=0; i<nroots; i++) {
         Cvec.copy(Dvec, i, i);
         }
      Cvec.buf_unlock();
      Dvec.buf_unlock();
      Sigma.buf_lock(buffer1);
      Dvec2.buf_lock(buffer2);
      for (i=0; i<nroots; i++) {
         Dvec2.gather(i, L, i, alpha[iter2], Sigma);
         }
      for (i=0; i<nroots; i++) {
         Sigma.copy(Dvec2, i, i);
         }
      Dvec2.buf_unlock();
      Sigma.buf_unlock();


      /*
      for (i=L; i<nroots+L; i++) {
         Cvec.restart_gather(i, L, i-L, alpha[iter2], buffer1, buffer2);
         Sigma.restart_gather(i, L, i-L, alpha[iter2], buffer1, buffer2);
         }

      Cvec.restart_reord_fp(L);
      Cvec2.restart_reord_fp(L);
      Sigma.restart_reord_fp(L);
      Sigma2.restart_reord_fp(L);
      */

      zero_mat(G, L, L);
      k = nroots;
   }

   /* previous-run d vector */
   else if (Parameters_->guess_vector==PARM_GUESS_VEC_DFILE && !dvec_read_fail) {
     if (print_) outfile->Printf( "    Attempting to use %d previous converged vectors\n\n",
        nroots);

     if (Parameters_->nodfile) {
       i = Cvec.read_new_first_buf();
       Cvec.set_new_first_buf(i);
       Cvec2.set_new_first_buf(i);
       j = Sigma.read_new_first_buf();
       Sigma.set_new_first_buf(j);
       Sigma2.set_new_first_buf(j);
       /* the first buffer of D file should not change unless there
        * is only a logical D file not a physical D file */
       Dvec.set_new_first_buf(i);
       Dvec2.set_new_first_buf(j);
       Dvec.restart_reord_fp(maxnvect-1);
       Dvec2.restart_reord_fp(maxnvect-1);
     }
     Cvec.buf_lock(buffer1);
     Dvec.buf_lock(buffer2);
     if ((i = Dvec.read_num_vecs()) < nroots) {
       str = "Only ";
       str += std::to_string( i) ;
       str += " vectors available in D file for ";
       str += std::to_string( nroots) ;
       str += " roots!";
       throw PsiException(str,__FILE__,__LINE__);
     }


     for (i=0; i<nroots; i++) {
       Cvec.copy(Dvec, i, i);
       }
     Cvec.buf_unlock();
     Dvec.buf_unlock();

     k = nroots;
   }

   /* unit vector */
   else if (Parameters_->guess_vector == PARM_GUESS_VEC_UNIT ||
            (dvec_read_fail && Parameters_->h0blocksize==0)) {
     tval = 1.0;
     Cvec.buf_lock(buffer1);
     Cvec.init_vals(0, 1, &(CalcInfo_->ref_alp_list), &(CalcInfo_->ref_alp_rel),
        &(CalcInfo_->ref_bet_list), &(CalcInfo_->ref_bet_rel), H0block_->blknum,
        &tval);
     Cvec.buf_unlock();
     Cvec.write_num_vecs(1);
     Sigma.set_zero_blocks_all();
     k = 1;
   }

   else { /* use H0BLOCK eigenvector guess */
      if (Parameters_->precon == PRECON_GEN_DAVIDSON) L = H0block_->size;
      else L = H0block_->guess_size;

      /* outfile->Printf( " L = %d in sem.cc line 345\n", L); */
      /* N = CIblks.vectlen; The variable N is never used */
      sm_evals = init_array(L);

      /* need to fill out sm_evecs into b (pad w/ 0's) */
      if (print_) outfile->Printf( "    Using %d initial trial vectors\n\n", Parameters_->num_init_vecs);

      Cvec.buf_lock(buffer1);
      for (i=0,k=0; i<L && k < Parameters_->num_init_vecs; i++) {

         /* if Ms=0 check sm_evecs[i] to see if it has the correct
          * spin symmetry
          */
         tmpi=0;
         for (j=0; Parameters_->Ms0 && j<L && !tmpi; j++) {
            l = H0block_->pair[j];
            if (l == -1) {
              outfile->Printf("(sem_iter): Warning: unpaired h0block member!\n");
               tmpi = 1;
               }
            tval = H0block_->H0b_diag[l][i];
            if ((int) Parameters_->S % 2) tval = -tval;
            if (fabs(H0block_->H0b_diag[j][i] - tval) > 1.0E-8) {
              tmpi = 1;
              outfile->Printf("(sem_iter): H0block_->H0b_diag[%d][%d]"
                      " - H0block_->H0b_diag[%d][%d] = %lf - %lf = %lf"
                      " > 1.0E-8\n", j, i, l, i, H0block_->H0b_diag[j][i],
                     tval, (H0block_->H0b_diag[j][i] - tval));
              }
            }

         /* also check that it satisfies any user-specified properties */
         if (!tmpi && Parameters_->filter_guess) {
       j = Parameters_->filter_guess_H0_det1;
       l = Parameters_->filter_guess_H0_det2;
           tval = H0block_->H0b_diag[l][i];
       if (Parameters_->filter_guess_sign == -1) tval = -tval;
       if (fabs(H0block_->H0b_diag[j][i] - tval) > 1.0E-8) {
         tmpi = 1;
         outfile->Printf( "(sem_iter): Guess vector failed user-specified"
                          " criterion.\n");
         outfile->Printf( "(sem_iter): H0block_->H0b_diag[%d][%d]"
                 " - H0block_->H0b_diag[%d][%d] = %lf - %lf = %lf"
             " > 1.0E-8\n", j, i, l, i, H0block_->H0b_diag[j][i],
             tval, (H0block_->H0b_diag[j][i] - tval));
       }
     }

         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = H0block_->H0b_diag[j][i];

         Cvec.init_vals(k, L, H0block_->alplist, H0block_->alpidx,
            H0block_->betlist, H0block_->betidx, H0block_->blknum, sm_evals);

         if (Parameters_->calc_ssq && Parameters_->icore==1) {
            Cvec.buf_unlock();
            tval = Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, k);
            Cvec.buf_lock(buffer1);
            if (fabs(tval - (Parameters_->S*(Parameters_->S+1.0))) > 1.0E-3) {
              outfile->Printf(
                 "Computed <S^2> not as desired, discarding guess\n");
              }
            else k++;
            }
         else k++;
         }
      Cvec.buf_unlock();
      Cvec.write_num_vecs(k);
      Sigma.set_zero_blocks_all();

      free(sm_evals);
      /* free(sm_mat);
      free_matrix(sm_evecs, L); */
      }

   /* Gather h0block elemts of initial CI_VEC - could move this to init_vals */
   Cvec.buf_lock(buffer1);
   Cvec.read(0,0);
   Cvec.symnorm(1.0,CI_VEC,TRUE);
   Cvec.buf_unlock();

   /* write file_offset and file_number array out to detci.dat */
   //Cvec.write_detfile(CI_VEC);
   //Sigma.write_detfile(SIGMA_VEC);
   //if (print_ > 1)
   //  outfile->Printf("Restart info written.\n");


   if (k < nroots) {
      str = "(sem_iter): Failure to get required number of guess vects.\n";
      str += "  Got ";
      str += std::to_string( k) ;
      str += ", need nroots=";
      str += std::to_string( nroots) ;
      str += " to proceed.  Aborting";
      throw PsiException(str,__FILE__,__LINE__);
   }

   L = k;
   Llast = 0;

   if (Parameters_->nodfile) {
     Dvec.restart_reord_fp(maxnvect-1);
     Dvec2.restart_reord_fp(maxnvect-1);
     }

   if (print_){
     outfile->Printf("     Iter   Root       Total Energy       Delta E      C RMS\n\n");
   }

   /* begin iteration */
   while (!converged && iter <= maxiter) {

      Lvec[iter2] = L;

      if (print_ > 2) {
        outfile->Printf( "L[cur] = %3d, L[last] = %3d\n", L,
          iter2 > 0 ? Lvec[iter2-1] : 999);
      }

      /* form contributions to the G matrix */
      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);

      for (i=Llast; i<L; i++) {
         Cvec.read(i, 0);
         if (print_ > 3) {
            outfile->Printf( "b[%d] =\n", i);
            Cvec.print();
            }

         sigma(Cvec, Sigma, oei, tei, i);

         if (Parameters_->z_scale_H) {
           Cvec.buf_unlock();
           Sigma.buf_unlock();
           Sigma.scale_sigma(Hd, Cvec, alplist, betlist, i, buffer1, buffer2);
           Cvec.buf_lock(buffer1);
           Sigma.buf_lock(buffer2);
           Sigma.read(i,0);
           }

         if (print_ > 3) { /* and this as well */
            outfile->Printf( "H * b[%d] = \n", i);
            Sigma.print();
            }

         for (j=0; j<L; j++) {
            Cvec.read(j, 0);
            G[j][i] = G[i][j] = Cvec * Sigma;
            }

         }
       Sigma.write_num_vecs(L);
       Llast = L;


      /* solve the L x L eigenvalue problem G a = lambda a for M roots */
      sq_rsp(L, L, G, lambda[iter2], 1, alpha[iter2], 1.0E-14);

      if (print_ > 4) {
         outfile->Printf( "\n G eigenvectors and eigenvalues:\n");
         eivout(alpha[iter2], lambda[iter2], L, L, "outfile");
         }

      if (print_ > 3) {
        outfile->Printf( "\nG matrix (%2d) = \n", iter);
        print_mat(G, L, L, "outfile");
        }

      Cvec.buf_unlock();
      Sigma.buf_unlock();

     if (Parameters_->lse && (maxnvect-L <= nroots*Parameters_->collapse_size)
         && L>2) form_M = 1;
     else form_M = 0;

     if (form_M) {
       /* Form sigma_overlap matrix */
       Sigma.buf_lock(buffer1);
       Sigma2.buf_lock(buffer2);
       zero_mat(sigma_overlap,maxnvect, maxnvect);
       for (i=0; i<L; i++) {
          Sigma.read(i, 0);
          if (print_ > 2) {
            outfile->Printf("Sigma[%d] = ", i);
            Sigma.print();

          }
          for (j=i; j<L; j++) {
             Sigma2.read(j, 0);
             if (print_ > 2) {
               outfile->Printf("Sigma2[%d] = ", j);
               Sigma2.print();

             }
             sigma_overlap[i][j] = sigma_overlap[j][i] = Sigma * Sigma2;
          }
       }
       Sigma.buf_unlock();
       Sigma2.buf_unlock();

       /* Form Mij matrix */
       /* This formula assumes the b vectors are orthogonal */
       for (k=0; k<nroots; k++) {
          zero_mat(M[k],maxnvect,maxnvect);
          for (i=0; i<L; i++) {
             for (j=i; j<L; j++) {
                M[k][i][j] = M[k][j][i] =  sigma_overlap[i][j]
                             -2.0 * lambda[iter2][k]* G[i][j];
                if (i==j) M[k][i][i] += lambda[iter2][k] * lambda[iter2][k];
                }
             }
          } /* end loop over k (nroots) */

       if (print_ > 2) {
         outfile->Printf( "\nsigma_overlap matrix (%2d) = \n", iter);
         print_mat(sigma_overlap, L, L, "outfile");

         for (k=0; k<nroots; k++) {
            outfile->Printf( "\nM matrix (%2d) for root %d = \n", iter, k);
            print_mat(M[k], L, L, "outfile");
            outfile->Printf( "\n");
            }
         }

       /* solve the L x L eigenvalue problem M a = lambda a for M roots */
       for (k=0; k<nroots; k++) {
          sq_rsp(L, L, M[k], m_lambda[iter2][k], 1, m_alpha[iter2][k], 1.0E-14);
          if (print_ > 2) {
            outfile->Printf( "\n M eigenvectors and eigenvalues root %d:\n",k);
            eivout(m_alpha[iter2][k], m_lambda[iter2][k], L, L, "outfile");
            }
          }

       }

     if (Parameters_->print_sigma_overlap) {
       /* Form sigma_overlap matrix */
       Sigma.buf_lock(buffer1);
       Sigma2.buf_lock(buffer2);
       zero_mat(sigma_overlap,maxnvect, maxnvect);
       for (i=0; i<L; i++) {
          Sigma.read(i, 0);
          if (print_ > 2) {
            outfile->Printf("Sigma[%d] = ", i);
            Sigma.print();

            }
          for (j=i; j<L; j++) {
             Sigma2.read(j, 0);
             if (print_ > 2) {
               outfile->Printf("Sigma2[%d] = ", j);
               Sigma2.print();

               }
             sigma_overlap[i][j] = sigma_overlap[j][i] = Sigma * Sigma2;
             }
          }
       Sigma.buf_unlock();
       Sigma2.buf_unlock();

       outfile->Printf( "\nsigma_overlap matrix (%2d) = \n", iter);
       print_mat(sigma_overlap, L, L, "outfile");

       for (i=0; i<L; i++) {
         outfile->Printf( "\nGuess energy #%d = %15.9lf\n", i,
           -1.0 * sqrt(sigma_overlap[i][i]) + CalcInfo_->enuc + CalcInfo_->edrc);
         }

       /* diagonalize sigma_overlap to see what that does
          should be the same as diagonalizing H^2 in the basis of the
          Davidson subspace vectors
        */

       /* Form Mij matrix */
       /* This formula assumes the b vectors are orthogonal */
       zero_mat(M[0],maxnvect,maxnvect);
       for (i=0; i<L; i++) {
          for (j=i; j<L; j++) {
             M[0][i][j] = M[0][j][i] =  sigma_overlap[i][j];
          }
       }

       /* solve the L x L eigenvalue problem M a = lambda a for M roots */
       sq_rsp(L, L, M[0], m_lambda[0][0], 1, m_alpha[0][0], 1.0E-14);
       for (i=0; i<L; i++) {
         m_lambda[0][0][i] = -1.0 * sqrt(m_lambda[0][0][i]) +
           CalcInfo_->enuc + CalcInfo_->edrc;
       }
       outfile->Printf( "\n Guess energy from H^2 = %15.9lf\n",
         m_lambda[0][0][L]);
       outfile->Printf( "\n M eigenvectors and eigenvalues root %d:\n",0);
       eivout(m_alpha[0][0], m_lambda[0][0], L, L, "outfile");
       }

      /* before we form correction vectors see if enough room to
       * append new b vectors to Davidson subspace.
       */

      if ((iter2 - Parameters_->collapse_size + 1 >= 0) && (Lvec[iter2 -
           Parameters_->collapse_size + 1] + nroots * Parameters_->collapse_size
           > maxnvect) && iter != maxiter) {

         Cvec.set_nvect(maxnvect);
         Cvec2.set_nvect(maxnvect);
         Sigma.set_nvect(maxnvect);
         if (Parameters_->lse) Sigma2.set_nvect(maxnvect);
         collapse_num++;
         lse_do = 0;

         Cvec.buf_lock(buffer1);
         Dvec.buf_lock(buffer2);
         zero_int_array(lse_do_arr, nroots);
         for (i=0; i<nroots; i++) {
            if (form_M && ((collapse_num-last_lse_collapse_num)
               >= Parameters_->lse_collapse) &&
               (fabs(lambda[iter2][i]-lastroot[i]) < lse_tolerance)
                && (m_lambda[iter2][i][0] > MALPHA_TOLERANCE)) {
               lse_do_arr[i] = 1; lse_do++;
              }
            else {
              lse_do_arr[i] = 0;
              if (i==0) break;
              }
            }

         for (i=0; i<nroots; i++) {
            if (lse_do_arr[i])
              for (j=0; j<L; j++) cmp_cncoe[j][i] = m_alpha[iter2][i][j][0];
            else
              for (j=0; j<L; j++) cmp_cncoe[j][i] = alpha[iter2][j][i];
            }

        /* transpose the cmp_cncoe matrix to prepare for schmidt orthog */
        zero_mat(tr_cmp_cncoe, maxnvect, maxnvect);
        for (i=0; i<L; i++)
           for (j=0; j<L; j++)
              tr_cmp_cncoe[i][j] = cmp_cncoe[j][i];

        tmp = 0.0;
        for (i=0; i<maxnvect; i++)
           tmp += cmp_cncoe[i][1] * tr_cmp_cncoe[0][i];

        if (lse_do && nroots>1) schmidt(tr_cmp_cncoe, nroots, L, "outfile");

        /* transpose the cmp_cncoe matrix to prepare for schmidt orthog */
        for (i=0; i<maxnvect; i++)
           for (j=0; j<maxnvect; j++)
              cmp_cncoe[i][j] = tr_cmp_cncoe[j][i];

      if (lse_do) last_lse_collapse_num = collapse_num;
      for (i=0; i<nroots; i++)
         Dvec.gather(i, L, i, cmp_cncoe, Cvec);

       for (i=0; i<nroots; i++)
          Cvec.copy(Dvec, maxnvect-nroots+i, i);

       Cvec.buf_unlock();
       Dvec.buf_unlock();
       Sigma.buf_lock(buffer1);
       Dvec2.buf_lock(buffer2);

       for (i=0; i<nroots; i++)
          Dvec2.gather(i, L, i, cmp_cncoe, Sigma);

       for (i=0; i<nroots; i++)
          Sigma.copy(Dvec2, maxnvect-nroots+i, i);

       Sigma.buf_unlock();
       Dvec2.buf_unlock();
       L2 = L3 = nroots;

       if (print_ > 2) {
         outfile->Printf( "Gathered vectors 0 to %d and wrote to positions \
            %d to %d\n", L-1, maxnvect-nroots, maxnvect-1);
       }

       for (i=1; i<Parameters_->collapse_size; i++) {

          if (Parameters_->nodfile) {
            Dvec.restart_reord_fp(maxnvect-2);
            Dvec2.restart_reord_fp(maxnvect-2);
            }

            /* do all the C's */
            zero_int_array(did_root, nroots);
            zero_mat(clpse_dot, nroots, maxnvect);
            Cvec.buf_lock(buffer1);
            Dvec.buf_lock(buffer2);
            for (j=0; j<nroots; j++)
               Dvec.gather(j, Lvec[iter2-i], j, alpha[iter2-i], Cvec);

            for (j=0; j<nroots; j++) {
               if (root_converged[j]) continue;
               if (Dvec.schmidt_add2(Cvec, maxnvect-L2, maxnvect-1,
                  j, maxnvect-L2-1, clpse_dot[j], &(clpse_norm[j]),&tval)) {
                  did_root[j] = 1;
                  L2++;
                  if (tval>ovlpmax) ovlpmax = tval;
                  }
               }
            Cvec.buf_unlock();
            Dvec.buf_unlock();

            /* do all the Sigmas */
            Sigma.buf_lock(buffer1);
            Dvec2.buf_lock(buffer2);
            for (j=0; j<nroots; j++) {
               if (!did_root[j]) continue;
                 for (k=0; k<Lvec[iter2-i]; k++)
                    tmpmat[k][j] = alpha[iter2-i][k][j] * clpse_norm[j];
               Dvec2.gather(j, Lvec[iter2-i], j, tmpmat, Sigma);
               }
            for (j=0; j<nroots; j++) {
               if (!did_root[j]) continue;
               for (k=maxnvect-L3; k<maxnvect; k++) {
                  Dvec2.civ_xpeay(-clpse_norm[j] * clpse_dot[j][k], Sigma, j,k);
                  }
               Sigma.copy(Dvec2, maxnvect-L3-1, j);
               L3++;
               }
            Sigma.buf_unlock();
            Dvec2.buf_unlock();
            } /* end collapse_size > 1 */

        if (L2 != L3) {
         outfile->Printf("(sem_iter): L2 != L3.  Bad. \n");
          outfile->Printf("(sem_iter): L2 != L3.  Bad. \n");
          }

        Cvec.restart_reord_fp(maxnvect-L2);
        Cvec2.restart_reord_fp(maxnvect-L2);

    Cvec.write_new_first_buf();
        Sigma.restart_reord_fp(maxnvect-L3);
        Sigma2.restart_reord_fp(maxnvect-L3);
    Sigma.write_new_first_buf();
        Cvec.set_nvect(L2);
        Cvec2.set_nvect(L2);
        Sigma.set_nvect(L3);
        Sigma2.set_nvect(L3);
    Cvec.write_num_vecs(L2);
    Sigma.write_num_vecs(L3);
        L = L2;
        Llast = L;
        iter2 = 0;  Lvec[0] = L;
        if (print_ > 2) {
          outfile->Printf( "L = %d, L2 = %d, L3 = %d\n",L, L2, L3);
        }

        /* write file_offset and file_number array out to detci.dat */
    /*
        Cvec.write_detfile(CI_VEC);
        Sigma.write_detfile(SIGMA_VEC);
        if (print_ > 1)
          outfile->Printf("Restart info written.\n");

    */

        if (Parameters_->nodfile) {
          Dvec.set_nvect(L2);
          Dvec2.set_nvect(L3);
          //Dvec.reset_detfile(CI_VEC);
          //Dvec2.reset_detfile(SIGMA_VEC);
          Dvec.restart_reord_fp(maxnvect-1);
          Dvec2.restart_reord_fp(maxnvect-1);
          }

        /* Schmidt-Orthogonalize again to ensure numerical stability */
        if (ovlpmax > S_MAX) {
           Cvec.buf_lock(buffer1);
           Cvec2.buf_lock(buffer2);
           L2 = L3 = 1;
           zero_mat(clpse_dot,maxnvect,maxnvect);
           zero_arr(clpse_norm,nroots);
           for (j=1; j<L; j++) {
              if (Cvec.schmidt_add2(Cvec2,0, j-1, j, j,
                clpse_dot[j],&(clpse_norm[j]),&ovlpmax)) L2++;
              }
           if (ovlpmax > S_MAX) outfile->Printf("Near degeneracy in b space\n");
           Cvec.buf_unlock();
           Cvec2.buf_unlock();
           Sigma.buf_lock(buffer1);
           Sigma2.buf_lock(buffer2);
           for (j=1; j<L; j++) {
              for (k=0; k<L3; k++)
                 Sigma.civ_xpeay(-clpse_norm[j]*clpse_dot[j][k],Sigma2,j,k);
              L3++;
              }
           L = L2;
           Sigma.buf_unlock();
           Sigma2.buf_unlock();
           outfile->Printf("  Second Schmidt-Orthogonalization performed.\n");
           }

        /* need to re-form G here and re-diagonalize it */
        zero_mat(G, maxnvect, maxnvect);
        Cvec.buf_lock(buffer1);
        Sigma.buf_lock(buffer2);


         /* MLL debug */
        /*
         if (1) {
           for (i=0; i<L; i++) {
              Cvec.read(i,0);
              sigma(Cvec, Sigma, oei, tei, i);
              if (print_ > 1) {
                outfile->Printf(
                  "Exact Sigma: (redid multiplication) H * b[%d] = \n", i);
                Sigma.print(outfile);
                }
              }
           }
         */

        /* Reforming G matrix after collapse */
        for (i=0; i<L; i++) {
           Cvec.read(i, 0);
           for (j=0; j<=i; j++) {
               Sigma.read(j, 0);
               G[j][i] = G[i][j] = Cvec * Sigma;
               }
           }

        /* solve the L x L eigenvalue problem G a = lambda a for M roots */
        sq_rsp(L, L, G, lambda[iter2], 1, alpha[iter2], 1.0E-14);

        if (print_ > 4) {
           outfile->Printf( "\n G eigenvectors and eigenvalues:\n");
           eivout(alpha[iter2], lambda[iter2], L, L, "outfile");
           }

        Cvec.buf_unlock();
        Sigma.buf_unlock();

        if (print_ > 1) {
           outfile->Printf("  Collapsed Davidson subspace to %d vectors\n",L);
           if (lse_do) {
             outfile->Printf("  Least Squares Extrapolation for Root%c",
                     (lse_do>1) ? 's' : ' ');
             for (i=0; i<nroots; i++)
                if (lse_do_arr[i]) outfile->Printf(" %d", i);
             outfile->Printf("\n");
             }
          }

        } /* end collapse routine */

        if (Parameters_->update == UPDATE_DAVIDSON) {
          /* form the d part of the correction vector */
          Dvec.dcalc(nroots, L, alpha[iter2], lambda[iter2], dvecnorm, Cvec,
                     Sigma, buffer1, buffer2, root_converged, (print_ > 4),
                     E_est);
          }
        else if (Parameters_->update == UPDATE_OLSEN) {
          /* Compute x and y values for E_est */
          Cvec.buf_lock(buffer1);
          Dvec.buf_lock(buffer2);
          for (i=0; i<nroots; i++) Dvec.gather(i,L,i,alpha[iter2],Cvec);
          Dvec.buf_unlock();
          Cvec.buf_unlock();
          for (i=0; i<nroots; i++) {
             olsen_iter_xy(Dvec,Sigma,Hd,&tmpx,&tmpy,buffer1,buffer2,
               lambda[iter2][i]+edrc,i,L,alpha[iter2], alplist, betlist);
             x[i] = tmpx;
             y[i] = tmpy;
             /* outfile->Printf("x[%d] = %lf    y[%d] = %lf\n",i,x[i],i,y[i]);
               E_est[i] += edrc; */
             errcod = H0block_calc(lambda[iter2][i]);
             if (!errcod)
               outfile->Printf("Determinant of H0block is too small.\n");
             if (Parameters_->precon>=PRECON_GEN_DAVIDSON)
               H0block_xy(&x[i],&y[i],lambda[iter2][i]);
        /*
             outfile->Printf(
                     "Modified x[%d] = %lf y[%d] = %lf\n",i,x[i],i,y[i]);
        */
             E_est[i] = y[i]/x[i];
             /*
               outfile->Printf("E_est[%d] = %20.12f lambda[%d] = %20.12f\n",i,
                    E_est[i]+edrc+enuc,i,lambda[iter2][i]+edrc+enuc);
             */
             }
         Dvec.dcalc(nroots,L,alpha[iter2],lambda[iter2],dvecnorm,Cvec,Sigma,
          buffer1,buffer2,root_converged,(print_ > 4),E_est);
         }
        else {
          throw PsiException("UPDATE option not recognized.  Choose DAVIDSON or OLSEN",__FILE__,__LINE__);
          }


      Cvec.copy_zero_blocks(Sigma);


      /* check for convergence */
      converged = 1;
      for (i=0; i<nroots; i++) {
         if (dvecnorm[i] <= conv_rms && fabs(lambda[iter2][i] - lastroot[i])
            <= conv_e) root_converged[i] = 1;
         else {
            root_converged[i] = 0;
            converged = 0;
            }
         if (print_) {
            outfile->Printf( "   @CI %2d:    %2d  %18.12lf   %10.4E   %10.4E %c\n",
               iter, i, (lambda[iter2][i] + enuc + edrc),
               (lambda[iter2][i] - lastroot[i]), dvecnorm[i],
               root_converged[i] ? 'c' : ' ');
         }
      }

      if ((nroots > 1) && print_) outfile->Printf( "\n");

      if (iter == maxiter && !Parameters_->mcscf) {
         outfile->Printf( "\nMaximum number of CI iterations reached\n");
         }

      if (converged){
          double avg_vec_norm = 0.0;
          for (i=0; i<nroots; i++){
            avg_vec_norm += dvecnorm[i] * Parameters_->average_weights[i];
          }
          Process::environment.globals["DETCI AVG DVEC NORM"] = avg_vec_norm;
          Parameters_->diag_h_converged = true;
      }

      if (converged || iter == maxiter) {

         Cvec.buf_lock(buffer1);
         Dvec.buf_lock(buffer2);
         //if (Parameters_->nodfile) Dvec.reset_detfile(CI_VEC);
         for (i=0; i<nroots; i++) {
            evals[i] = lambda[iter2][i];
            tval = alpha[iter2][0][i];
            Dvec.civ_xeay(tval, Cvec, i, 0);

            for (j=1; j<L; j++) {
               tval = alpha[iter2][j][i];
               Dvec.civ_xpeay(tval, Cvec, i, j);
               }

            // if (print_) {
            //   outfile->Printf( "\n* ROOT %d CI total energy = %17.13lf", i+1,
            //      evals[i] + enuc + edrc);

            //   if (nroots > 1) {
            //      outfile->Printf( "  (%6.4lf eV, %9.2lf 1/cm)\n",
            //        (evals[i] - evals[0]) * pc_hartree2ev,
            //        (evals[i] - evals[0]) * pc_hartree2wavenumbers);
            //      }
            //   else outfile->Printf( "\n");

            //    zero_arr(mi_coeff, Parameters_->nprint);
            //    Dvec.max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc,
            //       mi_iaidx, mi_ibidx, mi_coeff, Parameters_->neg_only);
            //    print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
            //       mi_coeff);
            //    outfile->Printf( "\n");
            // }
            Dvec.write_num_vecs(i+1);  // only if nodfile ?
         }
         Cvec.buf_unlock();
         Dvec.buf_unlock();
         break;
         }

      else {
         for (i=0; i<nroots; i++) lastroot[i] = lambda[iter2][i];
         }

      /*
        Cvec.buf_lock(buffer1);
        Cvec2.buf_lock(buffer2);
        for (i=0; i<L; i++) {
           Cvec.read(i,0);
           cknorm = Cvec.checknorm();
           outfile->Printf("\ncknorm for b vector %d = %20.15f\n\n", i,cknorm);
           outfile->Printf("before correctn vecs added to list of b vecs\n");
           outfile->Printf("\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile);
           for (j=i; j<L; j++) {
              Cvec2.read(j,0);
              outfile->Printf("Cvec2 (b vector) %d =\n", j);
              Cvec2.print(outfile);
              tvalmatt = Cvec * Cvec2;
              outfile->Printf("Cvec[%d] * Cvec[%d] = %20.15f\n",i,j,tvalmatt);
              }
           }
         Cvec.buf_unlock();
         Cvec2.buf_unlock();
       */

      /* form the correction vector and normalize */

      Dvec.buf_lock(buffer1);
      for (k=0; k<nroots; k++) {
         if (root_converged[k]) continue;
         Hd.buf_lock(buffer2);
         if (Parameters_->precon == PRECON_EVANGELISTI)
           tval = Dvec.dcalc_evangelisti(k, L, lambda[iter2][k]+edrc, Hd, Cvec,
                buffer1, buffer2, Parameters_->precon, L, alplist,
                betlist, alpha[iter2]);
         else{
             tval = Dvec.dcalc2(k, lambda[iter2][k]+edrc, Hd,
                Parameters_->precon, alplist, betlist);
             errcod = H0block_calc(lambda[iter2][k]+edrc); /* MLL */
          }
         if (Parameters_->precon >= PRECON_GEN_DAVIDSON && (iter >= 1)) {
           if (Parameters_->h0block_coupling && (iter >= 2))
             H0block_coupling_calc(lambda[iter2][k]+edrc);
           Dvec.h0block_buf_precon(&tval, k);
           }
         if (tval < 1.0E-13 && print_ > 0) {
           outfile->Printf("    Warning: Norm of "
                  "correction (root %d) is < 1.0E-13\n", k);
           }
         Dvec.read(k,0);
         if (Parameters_->filter_zero_det) {
           tval -= Dvec.zero_det(Parameters_->filter_zero_det_Iac,
                     Parameters_->filter_zero_det_Iaridx,
                     Parameters_->filter_zero_det_Ibc,
                     Parameters_->filter_zero_det_Ibridx);
         }
         tval = sqrt(1.0 / tval);
         Dvec.symnorm(tval,0,0);

         if (print_ > 4) {
            outfile->Printf( "\nsecond d matrix root %d\n", k);
            Dvec.print();
            }

         Hd.buf_unlock();
         Cvec.buf_lock(buffer2);
         /* Schmidt orthog and append d's to b */
         if (Dvec.schmidt_add(Cvec, L)) L++;
         Cvec.buf_unlock();
         if (L > maxnvect) {
            str = "(sem_iter): L(";
            str += std::to_string( L) ;
            str += ") > maxnvect(";
            str += std::to_string( maxnvect) ;
            str += "!  Aborting!";
            throw PsiException(str,__FILE__,__LINE__);
            }
         } /* end loop over roots for new expansion vectors */

        Dvec.buf_unlock();
    Cvec.write_num_vecs(L);

        /* MLL Debug 1-8-98 If CI vector is converged too tight the
        ** norm of the correction vector i.e. the residual vector
        ** may become exceedingly small resulting in numerical instabilities
        */


     /*
        Cvec.buf_lock(buffer1);
        for (i=0; i<L; i++) {
           Cvec.read(i,0);
           outfile->Printf("\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile);
           }
        Cvec.buf_unlock();

        Cvec.buf_lock(buffer1);
        Cvec2.buf_lock(buffer2);
        for (i=0; i<L; i++) {
           Cvec.read(i,0);
           cknorm = Cvec.checknorm();
           outfile->Printf("\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile);
           outfile->Printf("\ncknorm for b vector %d after"
                   " correction vectors = %20.15f\n\n",i,cknorm);
           for (j=i; j<L; j++) {
              Cvec2.read(j,0);
              tvalmatt = Cvec * Cvec2;
              outfile->Printf("Cvec[%d] * Cvec[%d] = %20.15f\n",i,j,tvalmatt);
              }
           }
         Cvec.buf_unlock();
         Cvec2.buf_unlock();
     */

        /* not sure that you really want to do this every iter...  CDS
        if (Parameters_->calc_ssq && Parameters_->icore==1) {
          for (k=0; k<L; k++)
            tval = Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, k);
          }
        */

        iter++;
        iter2++;
  } /* end iteration */
  Parameters_->diag_iters_taken = iter;

   /* PT correction */
   /*
   if (Parameters_->calc_pt_corr) {
     Dvec.buf_lock(buffer1);
     Dvec.read(0,0);
     Dvec.pt_correction();
     Dvec.buf_unlock();
   }
   */

   /* Compute S^2 */
   if (Parameters_->calc_ssq && Parameters_->icore==1) {
     for (k=0; k<nroots; k++)
       Dvec.calc_ssq(buffer1, buffer2, alplist, betlist, k);
   }

   Cvec.close_io_files(1);
   Sigma.close_io_files(1);
   if (Parameters_->nodfile == FALSE) Dvec.close_io_files(1);

   // Free N-D arrays
   for (i=0; i<maxnvect; i++)
      free_matrix(m_lambda[i], nroots);
   free(m_lambda);

   for (i=0; i<nroots; i++)
        free_matrix(M[i], maxnvect);
   free(M);

   for (i=0; i<maxnvect; i++) {
      free_matrix(alpha[i], maxnvect);
      }
   free(alpha);

   for (i=0; i<maxiter; i++) {
      for (j=0; j<nroots; j++) {
         free_matrix(m_alpha[i][j], maxnvect);
      }
     free(m_alpha[i]);
   }
   free(m_alpha);

   // Free buffers
   free(buffer1);
   free(buffer2);

   // Free arrays
   free(lastroot);  free(dvecnorm);     free(root_converged);
   free(did_root);  free(clpse_norm);   free(lse_do_arr);
   free(renorm_c);  free(x);            free(y);
   free(E_est);     free(mi_iac);       free(Lvec);
   free(mi_ibc);    free(mi_iaidx);     free(mi_ibidx);
   free(mi_coeff);

   // Free matrices
   free_matrix(clpse_dot, maxnvect);    free_matrix(G, maxnvect);
   free_matrix(tmpmat, maxnvect);       free_matrix(cmp_cncoe, maxnvect);
   free_matrix(tr_cmp_cncoe, maxnvect); free_matrix(sigma_overlap, maxnvect);
   free_matrix(lambda, maxnvect);
}

}} // namespace psi::detci
