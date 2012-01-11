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

/* #define DEBUG */ 

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterdset.h>
#include <physconst.h>
#include "structs.h"
#include "ci_tol.h"
#define EXTERN
#include "globals.h"
#include "civect.h"

namespace psi { namespace detci {

extern int H0block_calc(double E);
extern void H0block_xy(double *x, double *y, double E);
extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
   struct stringwr **alplist, struct stringwr **betlist,
   FILE *outfile);
extern void parse_import_vector(SlaterDetSet *sdset, int *i_alplist, 
   int *i_alpidx, int *i_betlist, int *i_betidx, int *i_blknums);

extern void H0block_coupling_calc(double E, struct stringwr **alplist,
   struct stringwr **betlist);

#define MALPHA_TOLERANCE 1E-15

void sem_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
      **betlist, double *evals, double conv_e, 
      double conv_rms, double enuc, double efzc, 
      int nroots, int maxiter, int maxnvect, FILE *outfile, int print_lvl)
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
   int last_lse_collapse_num = -Parameters.lse_collapse;
   double *x, *y, tmpx, tmpy;
   double lse_tolerance, *renorm_c, *E_est, ovlpmax=0.0;
   double cknorm, tvalmatt=0.0, tmp; /* Add by CDS for debugging purposes */
   int errcod; 
 
   CIvect Cvec;
   CIvect Cvec2;
   CIvect Sigma;
   CIvect Sigma2;
   CIvect Dvec;
   CIvect Dvec2;


   Cvec.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit, 
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
   Cvec2.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit, 
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
   Sigma.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_s_tmp_units, Parameters.first_s_tmp_unit,
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
   Sigma2.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
      CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
      CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
      CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
      Parameters.num_s_tmp_units, Parameters.first_s_tmp_unit,
      CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
   if (!Parameters.nodfile) {
     Dvec.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
        CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
        CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
        CalcInfo.nirreps, AlphaG->subgr_per_irrep, nroots,
        Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
        CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
     Dvec2.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
        CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
        CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
        CalcInfo.nirreps, AlphaG->subgr_per_irrep, nroots,
        Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
        CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
     }
   else {
     Dvec.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
        CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
        CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
        CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
        Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit,
        CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
     Dvec2.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
        CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
        CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
        CalcInfo.nirreps, AlphaG->subgr_per_irrep, maxnvect,
        Parameters.num_s_tmp_units, Parameters.first_s_tmp_unit,
        CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
     }

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
   num_alp_str = CalcInfo.num_alp_str;
   num_bet_str = CalcInfo.num_bet_str;
   if (Parameters.fci) oei = CalcInfo.tf_onel_ints;
   else oei = CalcInfo.gmat[0];
   tei = CalcInfo.twoel_ints;

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
   mi_iac = init_int_array(Parameters.nprint);
   mi_ibc = init_int_array(Parameters.nprint);
   mi_iaidx = init_int_array(Parameters.nprint);
   mi_ibidx = init_int_array(Parameters.nprint);
   mi_coeff = init_array(Parameters.nprint);

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

   if (Parameters.lse) lse_tolerance = Parameters.lse_tolerance;

   if (Parameters.nodfile == FALSE) {
     if (Parameters.guess_vector == PARM_GUESS_VEC_DFILE && 
         (i = Dvec.read_num_vecs()) != nroots) {
       fprintf(outfile, "D file contains %d not %d vectors.  Attempting ",
               i, nroots);
       if (Parameters.h0blocksize == 0) {
         Parameters.guess_vector == PARM_GUESS_VEC_UNIT;
         fprintf(outfile, "unit vector guess.\n");
       }
       else {
         Parameters.guess_vector = PARM_GUESS_VEC_H0_BLOCK;
         fprintf(outfile, "H0block guess.\n");
       }
     }
   }

   if (Parameters.restart) {  /* restart option! */
      // L = Parameters.restart_vecs;
      L = Cvec.read_num_vecs();
      i = Sigma.read_num_vecs();
      if (i != L) {
        fprintf(outfile, "%d C vectors and %d Sigma vectors.\n", i, L);
	if (i < L) {
          L = i;
	  Cvec.write_num_vecs(L);
        }
	fprintf(outfile, "Using %d vectors \n", L);
      }
      if (L < nroots) {
        fprintf(outfile, "\nRestart failed... %d vectors for %d roots\n", 
                L, nroots);
        exit(0);
      }
        
      fprintf(outfile, "\nAttempting Restart with %d vectors\n", L);

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
      if (Parameters.nodfile) {
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
      // fflush(outfile);

      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);    

      for (i=0; i<L; i++) {
         Sigma.read(i, 0);
         if (print_lvl > 4) { 
            fprintf(outfile, "Sigma[%d] =\n", i);
            Sigma.print(outfile);
            }
         for (j=0; j<=i; j++) {
            Cvec.read(j, 0);
            if (print_lvl > 4) { 
               fprintf(outfile, "C[%d] =\n", j);
               Cvec.print(outfile);
               }
            G[j][i] = G[i][j] = Cvec * Sigma;
            }
         }
      Cvec.buf_unlock();
      Sigma.buf_unlock();

      if (print_lvl > 3) {
         fprintf(outfile, "\nG matrix (%2d) = \n", iter);
         print_mat(G, L, L, outfile);
         }

      /* solve the L x L eigenvalue problem G a = lambda a for M roots */
      sq_rsp(L, L, G, lambda[iter2], 1, alpha[iter2], 1.0E-14);
      if (print_lvl > 4) {
         fprintf(outfile, "\n G eigenvectors and eigenvalues:\n");
         eivout(alpha[iter2], lambda[iter2], L, L, outfile);
         }

      /* loop over roots and write out the required number of vects */
      if (nroots + L > maxnvect) {
         fprintf(outfile,"Error: Can't do restart if maxnvect < nroots + L\n");
         exit(0);
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
   else if (Parameters.guess_vector == PARM_GUESS_VEC_DFILE) {
     fprintf(outfile, "Attempting to use %d previous converged vectors\n", 
        nroots);
     if (Parameters.nodfile) {
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
       fprintf(outfile, "Only %d vectors available in D file for %d roots!\n",
               i, nroots); 
       exit(0);
     }
     

     for (i=0; i<nroots; i++) {
       Cvec.copy(Dvec, i, i);
       }
     Cvec.buf_unlock();
     Dvec.buf_unlock();

     k = nroots;
   }

   /* unit vector */
   else if (Parameters.guess_vector == PARM_GUESS_VEC_UNIT) {
     tval = 1.0;
     Cvec.buf_lock(buffer1);
     Cvec.init_vals(0, 1, &(CalcInfo.ref_alp_list), &(CalcInfo.ref_alp_rel),
        &(CalcInfo.ref_bet_list), &(CalcInfo.ref_bet_rel), H0block.blknum, 
        &tval);
     Cvec.buf_unlock();
     Cvec.write_num_vecs(1);
     Sigma.set_zero_blocks_all();
     k = 1;
   } 

   /* import a previously exported CI vector */
   else if (Parameters.guess_vector == PARM_GUESS_VEC_IMPORT) {

     SlaterDetSet *dets;
     int *import_alplist, *import_alpidx, *import_betlist, *import_betidx;
     int *import_blknums;

     slaterdetset_read(PSIF_CIVECT, "CI vector", &dets);

     // store the alpha graph, relative alpha index, beta graph, relative
     // beta index, and CI block number for each imported determinant
     import_alplist = init_int_array(dets->size);
     import_alpidx  = init_int_array(dets->size);
     import_betlist = init_int_array(dets->size);
     import_betidx  = init_int_array(dets->size); 
     import_blknums = init_int_array(dets->size); 

     parse_import_vector(dets, import_alplist, import_alpidx, import_betlist,
       import_betidx, import_blknums);
     
     k=0;
     for (i=0; i<nroots; i++) {

       zero_arr(buffer2, dets->size);
       slaterdetset_read_vect(PSIF_CIVECT, "CI vector", buffer2, 
         dets->size, i);

       // initialize the values in Cvec
       Cvec.buf_lock(buffer1);
       Cvec.init_vals(i, dets->size, import_alplist, import_alpidx,
         import_betlist, import_betidx, import_blknums, buffer2);
       Cvec.buf_unlock();
       k++; // increment number of vectors
     }

     Cvec.write_num_vecs(k);
     Sigma.set_zero_blocks_all();

     // when we're done, free the memory
     slaterdetset_delete_full(dets); 
     free(import_alplist);  free(import_alpidx);
     free(import_betlist);  free(import_betidx);
     free(import_blknums); 
   }
 
   else { /* use H0BLOCK eigenvector guess */
      if (Parameters.precon == PRECON_GEN_DAVIDSON) L = H0block.size;
      else L = H0block.guess_size;

      /* fprintf(outfile, " L = %d in sem.cc line 345\n", L); */
      /* N = CIblks.vectlen; The variable N is never used */
      sm_evals = init_array(L);
     
      /* need to fill out sm_evecs into b (pad w/ 0's) */
      fprintf(outfile, "Using %d initial trial vectors\n", 
         Parameters.num_init_vecs);

      Cvec.buf_lock(buffer1);
      for (i=0,k=0; i<L && k < Parameters.num_init_vecs; i++) {

         /* if Ms=0 check sm_evecs[i] to see if it has the correct 
          * spin symmetry 
          */
         tmpi=0;
         for (j=0; Parameters.Ms0 && j<L && !tmpi; j++) {
            l = H0block.pair[j];
            if (l == -1) { 
               printf("(sem_iter): Warning: unpaired h0block member!\n");
               tmpi = 1; 
               }
            tval = H0block.H0b_diag[l][i];
            if ((int) Parameters.S % 2) tval = -tval;
            if (fabs(H0block.H0b_diag[j][i] - tval) > 1.0E-8) {
              tmpi = 1;
              fprintf(outfile,"(sem_iter): H0block.H0b_diag[%d][%d]" 
                      " - H0block.H0b_diag[%d][%d] = %lf - %lf = %lf"
                      " > 1.0E-8\n", j, i, l, i, H0block.H0b_diag[j][i], 
                     tval, (H0block.H0b_diag[j][i] - tval));
              }
            }

         /* also check that it satisfies any user-specified properties */
         if (!tmpi && Parameters.filter_guess) {
	   j = Parameters.filter_guess_H0_det1;
	   l = Parameters.filter_guess_H0_det2;
           tval = H0block.H0b_diag[l][i];  
	   if (Parameters.filter_guess_sign == -1) tval = -tval;
	   if (fabs(H0block.H0b_diag[j][i] - tval) > 1.0E-8) {
	     tmpi = 1;
	     fprintf(outfile, "(sem_iter): Guess vector failed user-specified"
	                      " criterion.\n");
	     fprintf(outfile, "(sem_iter): H0block.H0b_diag[%d][%d]"
	             " - H0block.H0b_diag[%d][%d] = %lf - %lf = %lf"
		     " > 1.0E-8\n", j, i, l, i, H0block.H0b_diag[j][i],
		     tval, (H0block.H0b_diag[j][i] - tval));
	   }
	 }

         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = H0block.H0b_diag[j][i];

         Cvec.init_vals(k, L, H0block.alplist, H0block.alpidx, 
            H0block.betlist, H0block.betidx, H0block.blknum, sm_evals);

         if (Parameters.calc_ssq && Parameters.icore==1) {
            Cvec.buf_unlock();
            tval = Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, k);
            Cvec.buf_lock(buffer1);  
            if (fabs(tval - (Parameters.S*(Parameters.S+1.0))) > 1.0E-3) {
              fprintf(outfile, 
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
   //if (Parameters.print_lvl > 1) 
   //  fprintf(outfile,"Restart info written.\n");
   fflush(outfile);

   if (k < nroots) {
      printf("(sem_iter): Failure to get required number of guess vects.\n");
      printf("  Got %d, need nroots=%d to proceed.  Aborting\n", k, nroots);
      exit(1);
   }

   L = k;
   Llast = 0;
 
   if (Parameters.nodfile) {
     Dvec.restart_reord_fp(maxnvect-1);
     Dvec2.restart_reord_fp(maxnvect-1);
     }

   /* begin iteration */
   while (!converged && iter <= maxiter) {

      Lvec[iter2] = L;

      #ifdef DEBUG
      fprintf(outfile, "L[cur] = %3d, L[last] = %3d\n", L, 
        iter2 > 0 ? Lvec[iter2-1] : 999);
      #endif

      /* form contributions to the G matrix */
      Cvec.buf_lock(buffer1);
      Sigma.buf_lock(buffer2);

      for (i=Llast; i<L; i++) {
         Cvec.read(i, 0);
         if (print_lvl > 3) { 
            fprintf(outfile, "b[%d] =\n", i);
            Cvec.print(outfile);
            }

         sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, i);

         if (Parameters.z_scale_H) {
           Cvec.buf_unlock();
           Sigma.buf_unlock();
           Sigma.scale_sigma(Hd, Cvec, alplist, betlist, i, buffer1, buffer2);
           Cvec.buf_lock(buffer1);
           Sigma.buf_lock(buffer2);
           Sigma.read(i,0);
           }

         if (print_lvl > 3) { /* and this as well */
            fprintf(outfile, "H * b[%d] = \n", i);
            Sigma.print(outfile);
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

      if (print_lvl > 4) {
         fprintf(outfile, "\n G eigenvectors and eigenvalues:\n");
         eivout(alpha[iter2], lambda[iter2], L, L, outfile);
         }

      if (print_lvl > 3) {
        fprintf(outfile, "\nG matrix (%2d) = \n", iter);
        print_mat(G, L, L, outfile);
        }

      Cvec.buf_unlock();
      Sigma.buf_unlock();

     if (Parameters.lse && (maxnvect-L <= nroots*Parameters.collapse_size) 
         && L>2) form_M = 1; 
     else form_M = 0;

     if (form_M) {
       /* Form sigma_overlap matrix */
       Sigma.buf_lock(buffer1);
       Sigma2.buf_lock(buffer2);
       zero_mat(sigma_overlap,maxnvect, maxnvect);
       for (i=0; i<L; i++) {
          Sigma.read(i, 0);
          if (print_lvl > 2) {
            fprintf(outfile,"Sigma[%d] = ", i);
            Sigma.print(outfile);
            fflush(outfile);
            }
          for (j=i; j<L; j++) {
             Sigma2.read(j, 0);
             if (print_lvl > 2) {
               fprintf(outfile,"Sigma2[%d] = ", j);
               Sigma2.print(outfile);
               fflush(outfile);
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

       if (print_lvl > 2) {
         fprintf(outfile, "\nsigma_overlap matrix (%2d) = \n", iter);
         print_mat(sigma_overlap, L, L, outfile);
 
         for (k=0; k<nroots; k++) {
            fprintf(outfile, "\nM matrix (%2d) for root %d = \n", iter, k);
            print_mat(M[k], L, L, outfile);
            fprintf(outfile, "\n");
            }
         }
 
       /* solve the L x L eigenvalue problem M a = lambda a for M roots */
       for (k=0; k<nroots; k++) {
          sq_rsp(L, L, M[k], m_lambda[iter2][k], 1, m_alpha[iter2][k], 1.0E-14);
          if (print_lvl > 2) {
            fprintf(outfile, "\n M eigenvectors and eigenvalues root %d:\n",k);
            eivout(m_alpha[iter2][k], m_lambda[iter2][k], L, L, outfile);
            }
          } 
 
       }

     if (Parameters.print_sigma_overlap) {
       /* Form sigma_overlap matrix */
       Sigma.buf_lock(buffer1);
       Sigma2.buf_lock(buffer2);
       zero_mat(sigma_overlap,maxnvect, maxnvect);
       for (i=0; i<L; i++) {
          Sigma.read(i, 0);
          if (print_lvl > 2) {
            fprintf(outfile,"Sigma[%d] = ", i);
            Sigma.print(outfile);
            fflush(outfile);
            }
          for (j=i; j<L; j++) {
             Sigma2.read(j, 0);
             if (print_lvl > 2) {
               fprintf(outfile,"Sigma2[%d] = ", j);
               Sigma2.print(outfile);
               fflush(outfile);
               }
             sigma_overlap[i][j] = sigma_overlap[j][i] = Sigma * Sigma2;
             }
          }
       Sigma.buf_unlock();
       Sigma2.buf_unlock();

       fprintf(outfile, "\nsigma_overlap matrix (%2d) = \n", iter);
       print_mat(sigma_overlap, L, L, outfile);

       for (i=0; i<L; i++) {
         fprintf(outfile, "\nGuess energy #%d = %15.9lf\n", i,
           -1.0 * sqrt(sigma_overlap[i][i]) + CalcInfo.enuc + CalcInfo.efzc);
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
           CalcInfo.enuc + CalcInfo.efzc;
       }
       fprintf(outfile, "\n Guess energy from H^2 = %15.9lf\n", 
         m_lambda[0][0][L]);
       fprintf(outfile, "\n M eigenvectors and eigenvalues root %d:\n",0);
       eivout(m_alpha[0][0], m_lambda[0][0], L, L, outfile);
       }

      /* before we form correction vectors see if enough room to 
       * append new b vectors to Davidson subspace.
       */
 
      if ((iter2 - Parameters.collapse_size + 1 >= 0) && (Lvec[iter2 - 
           Parameters.collapse_size + 1] + nroots * Parameters.collapse_size 
           > maxnvect) && iter != maxiter) {

         Cvec.set_nvect(maxnvect);
         Cvec2.set_nvect(maxnvect);
         Sigma.set_nvect(maxnvect);
         if (Parameters.lse) Sigma2.set_nvect(maxnvect);
         collapse_num++;
         lse_do = 0;

         Cvec.buf_lock(buffer1);
         Dvec.buf_lock(buffer2);
         zero_int_array(lse_do_arr, nroots);
         for (i=0; i<nroots; i++) {
            if (form_M && ((collapse_num-last_lse_collapse_num) 
               >= Parameters.lse_collapse) && 
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

        if (lse_do && nroots>1) schmidt(tr_cmp_cncoe, nroots, L, outfile);

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
 
       #ifdef DEBUG
       fprintf(outfile, "Gathered vectors 0 to %d and wrote to positions \
          %d to %d\n", L-1, maxnvect-nroots, maxnvect-1);
       #endif

       for (i=1; i<Parameters.collapse_size; i++) {

          if (Parameters.nodfile) {
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
          printf("(sem_iter): L2 != L3.  Bad. \n");
          fprintf(outfile,"(sem_iter): L2 != L3.  Bad. \n");
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
        #ifdef DEBUG
        fprintf(outfile, "L = %d, L2 = %d, L3 = %d\n",L, L2, L3);
        #endif

        /* write file_offset and file_number array out to detci.dat */
	/*
        Cvec.write_detfile(CI_VEC);
        Sigma.write_detfile(SIGMA_VEC);
        if (Parameters.print_lvl > 1) 
          fprintf(outfile,"Restart info written.\n");
        fflush(outfile);
	*/

        if (Parameters.nodfile) {
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
           if (ovlpmax > S_MAX) fprintf(outfile,"Near degeneracy in b space\n");
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
           fprintf(outfile,"  Second Schmidt-Orthogonalization performed.\n");
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
              sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, i);
              if (print_lvl > 1) { 
                fprintf(outfile, 
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

        if (print_lvl > 4) {
           fprintf(outfile, "\n G eigenvectors and eigenvalues:\n");
           eivout(alpha[iter2], lambda[iter2], L, L, outfile);
           }

        Cvec.buf_unlock();
        Sigma.buf_unlock();

        if (print_lvl > 1) {
           fprintf(outfile,"  Collapsed Davidson subspace to %d vectors\n",L);
           if (lse_do) { 
             fprintf(outfile,"  Least Squares Extrapolation for Root%c",
                     (lse_do>1) ? 's' : ' ');
             for (i=0; i<nroots; i++) 
                if (lse_do_arr[i]) fprintf(outfile," %d", i); 
             fprintf(outfile,"\n");   
             }
          }

        } /* end collapse routine */

        if (Parameters.update == UPDATE_DAVIDSON) {
          /* form the d part of the correction vector */
          Dvec.dcalc(nroots, L, alpha[iter2], lambda[iter2], dvecnorm, Cvec, 
                     Sigma, buffer1, buffer2, root_converged, (print_lvl > 4), 
                     outfile, E_est);
          }
        else if (Parameters.update == UPDATE_OLSEN) {
          /* Compute x and y values for E_est */ 
          Cvec.buf_lock(buffer1);
          Dvec.buf_lock(buffer2);
          for (i=0; i<nroots; i++) Dvec.gather(i,L,i,alpha[iter2],Cvec); 
          Dvec.buf_unlock();
          Cvec.buf_unlock();
          for (i=0; i<nroots; i++) {
             olsen_iter_xy(Dvec,Sigma,Hd,&tmpx,&tmpy,buffer1,buffer2,
               lambda[iter2][i]+efzc,i,L,alpha[iter2], alplist, betlist);
             x[i] = tmpx;
             y[i] = tmpy; 
             /* fprintf(outfile,"x[%d] = %lf    y[%d] = %lf\n",i,x[i],i,y[i]);  
               E_est[i] += efzc; */ 
             errcod = H0block_calc(lambda[iter2][i]);
             if (!errcod) 
               fprintf(outfile,"Determinant of H0block is too small.\n");
             if (Parameters.precon>=PRECON_GEN_DAVIDSON) 
               H0block_xy(&x[i],&y[i],lambda[iter2][i]);
        /*   
             fprintf(outfile,
                     "Modified x[%d] = %lf y[%d] = %lf\n",i,x[i],i,y[i]); 
        */
             E_est[i] = y[i]/x[i]; 
             /* 
               fprintf(outfile,"E_est[%d] = %20.12f lambda[%d] = %20.12f\n",i,
                    E_est[i]+efzc+enuc,i,lambda[iter2][i]+efzc+enuc);
             */
             }
         Dvec.dcalc(nroots,L,alpha[iter2],lambda[iter2],dvecnorm,Cvec,Sigma,
          buffer1,buffer2,root_converged,(print_lvl > 4),outfile,E_est);
         }
        else { 
          fprintf(outfile,
                  "UPDATE option not recognized.  Choose DAVIDSON or OLSEN\n");
          exit(0);
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
         fprintf(outfile, "Iter %2d  Root %2d = %13.9lf",
            iter, i+1, (lambda[iter2][i] + enuc + efzc));
         fprintf(outfile, "   Delta_E %10.3E   Delta_C %10.3E %c\n",
            lambda[iter2][i] - lastroot[i], dvecnorm[i],
            root_converged[i] ? 'c' : ' ');
         fflush(outfile);
         }

      if (nroots > 1) fprintf(outfile, "\n");

      if (iter == maxiter) {
         fprintf(outfile, "\nMaximum number of iterations reached\n");
         }

      if (converged || iter == maxiter) {
         fflush(outfile);
         Cvec.buf_lock(buffer1);
         Dvec.buf_lock(buffer2);
         //if (Parameters.nodfile) Dvec.reset_detfile(CI_VEC);
         for (i=0; i<nroots; i++) {
            evals[i] = lambda[iter2][i];
            tval = alpha[iter2][0][i];
            Dvec.civ_xeay(tval, Cvec, i, 0);

            for (j=1; j<L; j++) {
               tval = alpha[iter2][j][i];
               Dvec.civ_xpeay(tval, Cvec, i, j);
               }

            fprintf(outfile, "\n* ROOT %d CI total energy = %17.13lf", i+1, 
               evals[i] + enuc + efzc);

            if (nroots > 1) {
               fprintf(outfile, "  (%6.4lf eV, %9.2lf 1/cm)\n", 
                 (evals[i] - evals[0]) * _hartree2ev,
                 (evals[i] - evals[0]) * _hartree2wavenumbers);
               }
            else fprintf(outfile, "\n");

            if (Parameters.print_lvl) {
               zero_arr(mi_coeff, Parameters.nprint);
               Dvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc,
                  mi_iaidx, mi_ibidx, mi_coeff, Parameters.neg_only);
               print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx, 
                  mi_coeff, AlphaG, BetaG, alplist, betlist, outfile);
               fprintf(outfile, "\n");
               }
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
           fprintf(outfile,"\ncknorm for b vector %d = %20.15f\n\n", i,cknorm);
           fprintf(outfile,"before correctn vecs added to list of b vecs\n");
           fprintf(outfile,"\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile); 
           for (j=i; j<L; j++) {
              Cvec2.read(j,0);
              fprintf(outfile,"Cvec2 (b vector) %d =\n", j);
              Cvec2.print(outfile);
              tvalmatt = Cvec * Cvec2;
              fprintf(outfile,"Cvec[%d] * Cvec[%d] = %20.15f\n",i,j,tvalmatt);
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
         if (Parameters.precon == PRECON_EVANGELISTI)
           tval = Dvec.dcalc_evangelisti(k, L, lambda[iter2][k]+efzc, Hd, Cvec, 
                buffer1, buffer2, Parameters.precon, L, alplist, 
                betlist, alpha[iter2]); 
         else tval = Dvec.dcalc2(k, lambda[iter2][k]+efzc, Hd,
                Parameters.precon, alplist, betlist);
         if (Parameters.precon >= PRECON_GEN_DAVIDSON && (iter >= 1)) { 
           if (Parameters.h0block_coupling && (iter >= 2)) 
             H0block_coupling_calc(lambda[iter2][k]+efzc, alplist, betlist);
           Dvec.h0block_buf_precon(&tval, k);
           }
         if (tval < 1.0E-13 && print_lvl > 0) { 
           fprintf(outfile,"Warning: Norm of "  
                  "correction (root %d) is < 1.0E-13\n", k);  
           }
         Dvec.read(k,0); 
         if (Parameters.filter_zero_det) {
           tval -= Dvec.zero_det(Parameters.filter_zero_det_Iac,
                     Parameters.filter_zero_det_Iaridx, 
                     Parameters.filter_zero_det_Ibc,
                     Parameters.filter_zero_det_Ibridx);
         }
         tval = sqrt(1.0 / tval);
         Dvec.symnorm(tval,0,0);

         if (print_lvl > 4) {
            fprintf(outfile, "\nsecond d matrix root %d\n", k);
            Dvec.print(outfile);
            }

         Hd.buf_unlock();  
         Cvec.buf_lock(buffer2);
         /* Schmidt orthog and append d's to b */
         if (Dvec.schmidt_add(Cvec, L)) L++;
         Cvec.buf_unlock();
         if (L > maxnvect) {
            fprintf(outfile, "(sem_iter): L(%2d) > maxnvect(%2d)!",L,maxnvect);
            fprintf(outfile, " Aborting!\n");
            exit(0);
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
           fprintf(outfile,"\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile);
           }
        Cvec.buf_unlock();

        Cvec.buf_lock(buffer1);
        Cvec2.buf_lock(buffer2);
        for (i=0; i<L; i++) {
           Cvec.read(i,0); 
           cknorm = Cvec.checknorm();
           fprintf(outfile,"\nCvec (b vector) %d =\n", i);
           Cvec.print(outfile); 
           fprintf(outfile,"\ncknorm for b vector %d after" 
                   " correction vectors = %20.15f\n\n",i,cknorm);
           for (j=i; j<L; j++) {
              Cvec2.read(j,0);
              tvalmatt = Cvec * Cvec2;
              fprintf(outfile,"Cvec[%d] * Cvec[%d] = %20.15f\n",i,j,tvalmatt);
              } 
           }
         Cvec.buf_unlock();
         Cvec2.buf_unlock();
     */

        /* not sure that you really want to do this every iter...  CDS
        if (Parameters.calc_ssq && Parameters.icore==1) {
          for (k=0; k<L; k++) 
            tval = Cvec.calc_ssq(buffer1, buffer2, alplist, betlist, k); 
          }
        */

        iter++;
        iter2++;
      } /* end iteration */


   /* Dump the vector to a PSIO file
      Added by Edward valeev (August 2002) */
   if (Parameters.export_ci_vector && Parameters.icore==1) {
     StringSet alphastrings, betastrings;
     SlaterDetSet dets;
     //SlaterDetVector vec;
     short int *fzc_occ;
     unsigned char *newocc;
     int irrep, gr, l, n;

     if (CalcInfo.num_fzc_orbs > 0) {
       fzc_occ = (short int *) malloc(CalcInfo.num_fzc_orbs*sizeof(short int));
       for (int l=0; l<CalcInfo.num_fzc_orbs; l++) {
         fzc_occ[l] = CalcInfo.order[l]; /* put it in Pitzer order */
       }
     }

     newocc = (unsigned char *) malloc(((AlphaG->num_el > BetaG->num_el) ? 
       AlphaG->num_el : BetaG->num_el)*sizeof(unsigned char));

     stringset_init(&alphastrings,AlphaG->num_str,AlphaG->num_el,
                    CalcInfo.num_fzc_orbs, fzc_occ);
     int list_gr = 0;
     int offset = 0;
     for(irrep=0; irrep<AlphaG->nirreps; irrep++) {
       for(gr=0; gr<AlphaG->subgr_per_irrep; gr++,list_gr++) {
         int nlists_per_gr = AlphaG->sg[irrep][gr].num_strings;
         for(l=0; l<nlists_per_gr; l++) {
           /* convert occs to Pitzer order */
           for (n=0; n<AlphaG->num_el; n++) {
             newocc[n] = (unsigned char) 
               CalcInfo.order[alplist[list_gr][l].occs[n] + 
               CalcInfo.num_fzc_orbs];
           }
	   stringset_add(&alphastrings,l+offset,newocc);
         }
	 offset += nlists_per_gr;
       }
     }
   
     stringset_init(&betastrings,BetaG->num_str,BetaG->num_el,
                    CalcInfo.num_fzc_orbs, fzc_occ);
     list_gr = 0;
     offset = 0;
     for(irrep=0; irrep<BetaG->nirreps; irrep++) {
       for(gr=0; gr<BetaG->subgr_per_irrep; gr++,list_gr++) {
         int nlists_per_gr = BetaG->sg[irrep][gr].num_strings;
         for(l=0; l<nlists_per_gr; l++) {
           /* convert occs to Pitzer order */
           for (n=0; n<BetaG->num_el; n++) {
             newocc[n] = (unsigned char) 
               CalcInfo.order[betlist[list_gr][l].occs[n] +
               CalcInfo.num_fzc_orbs];
           }
	   stringset_add(&betastrings,l+offset,newocc);
         }
	 offset += nlists_per_gr;
       }
     }
     free(newocc);
     if (CalcInfo.num_fzc_orbs > 0)
       free(fzc_occ);

     int ii;
     int size = CIblks.vectlen;
     int Iarel, Ialist, Ibrel, Iblist;
     slaterdetset_init(&dets,size,&alphastrings,&betastrings);
     for (ii=0; ii<size; ii++) {
       Dvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
       int irrep = Ialist/AlphaG->subgr_per_irrep;
       int gr = Ialist%AlphaG->subgr_per_irrep;
       int Ia = Iarel + AlphaG->list_offset[Ialist];
       irrep = Iblist/BetaG->subgr_per_irrep;
       gr = Iblist%BetaG->subgr_per_irrep;
       int Ib = Ibrel + BetaG->list_offset[Iblist];
       slaterdetset_add(&dets, ii, Ia, Ib);
     }

     // Don't init, don't need the memory allocated
     // slaterdetvector_init(&vec, &dets);

     Dvec.buf_lock(buffer1);
     for (ii=0; ii<Parameters.num_export; ii++) {
       zero_arr(buffer1, size);
       Dvec.read(ii,0);
       // slaterdetvector_set(&vec, buffer1);
       // slaterdetvector_write(PSIF_CIVECT,"CI vector",&vec);
       slaterdetset_write(PSIF_CIVECT,"CI vector",&dets);
       slaterdetset_write_vect(PSIF_CIVECT,"CI vector",buffer1,size,ii);
     }

     Dvec.buf_unlock();
     slaterdetset_delete_full(&dets);
   }
   else if (Parameters.export_ci_vector && Parameters.icore != 1) {
     fprintf(outfile, "\nWarning: requested CI vector export, unavailable " \
       "for icore = %d\n", Parameters.icore);
   }

   /* PT correction */
   /*
   if (Parameters.calc_pt_corr) {
     Dvec.buf_lock(buffer1);
     Dvec.read(0,0);
     Dvec.pt_correction();
     Dvec.buf_unlock();
   }
   */

   /* Compute S^2 */
   if (Parameters.calc_ssq && Parameters.icore==1) {
     for (k=0; k<nroots; k++) 
       Dvec.calc_ssq(buffer1, buffer2, alplist, betlist, k); 
   }

   Cvec.close_io_files(1);
   Sigma.close_io_files(1);
   if (Parameters.nodfile == FALSE) Dvec.close_io_files(1);

   free(mi_iac); free(mi_ibc); free(mi_iaidx); free(mi_ibidx); free(mi_coeff);
   free(dvecnorm);  free(lastroot);  free(root_converged);
   free(clpse_norm);  free(did_root);  free_matrix(clpse_dot, nroots);
   free_matrix(tmpmat, maxnvect);  free(Lvec);
   free(buffer1);
   free(buffer2);
}

}} // namespace psi::detci

