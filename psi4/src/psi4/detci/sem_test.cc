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
** Simultaneous Expansion Method for the Iterative Solution of
** Several of the Lowest Eigenvalues and Corresponding Eivenvectors of
** Large Real-Symmetric Matrices
**
** Algorithm due to Bowen Liu
** IBM Research Laboratory
**
** Implemented for Schaefer Group by David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
** In-core version for now!
** February 1994
**
** Updated 12/7/95 for testing within rasci code
** Updated 11/21/97 for least squares extrapolation and debugging of sem
**                  code
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

#define MAX_B_ROWS 200
#define MIN_F_DENOM 1.0E-3


/*** do a little test routine
main()
{
double **A ;
double *evals, **evecs ;
int i, j, used;
void sem() ;
std::string OutFileRMR ;

   ffile(&outfile, "output.dat", 0) ;
   tstart(outfile) ;

   A = init_matrix(50,50) ;
   evals = init_array(4) ;
   evecs = init_matrix(4, 50) ;
   for (i=0; i<50; i++) {
      for (j=0; j<=i; j++) {
         if (i!=j) {A[i][j] = 1.0; A[j][i] = 1.0; }
         else if (i<5) A[i][j] = 1.0 + 0.1 * (double) i ;
         else A[i][j] = 2.0 * (double) i + 1.0;
         }
      }
   sem(A, 50, 4, evecs, evals, 1.0E-10, 6, &used);


   outfile->Printf( "Ok, the eigenvectors are sideways!\n");
   eivout(evecs, evals, 4, 50, outfile);
   outfile->Printf( "\nused %d expansion vectors\n", used);

   tstop(outfile);
   fclose(outfile);
}
***/


/*
** sem(): Use Liu's Simultaneous Expansion Method to find the lowest
**    few eigenvalues of a real matrix.
**
** Arguments:
**   A        =  matrix to find eigenvalues of
**   N        =  size of matrix
**   M        =  number of eigenvalues to solve for
**   L        =  number of initial vectors in subspace
**   evecs    =  matrix for eigenvectors
**   evals    =  array for eigenvalues
**   b        =  set of subspace vectors (dimensions maxiter x N)
**                 space for rows i < L should not yet be allocated!!
**   conv_e   =  convergence tolerance.  The lowest energy eigenvalue must
**                 be converged to within this range.  It is interesting
**                 that higher roots _may_ converge faster.
**   conv_rms =  the required tolerance for convergence of the CI correction
**                 vector
**   maxiter  =  max number of iterations allowed
**   offst    =  offset to add to eigenvalues in printing (e.g. enuc)
**   vu       =  pointer to int to hold how many expansion vectors used
**   outfile  =  output file
**
** Returns: none
*/
void CIWavefunction::sem_test(double **A, int N, int M, int L, double **evecs, double *evals,
      double **b, double conv_e, double conv_rms, int maxiter, double offst,
      int *vu, int maxnvect)
{
   double *tmp_vec, **tmp_mat ;
   double **jnk;
   int sm_tridim;
   double *sm_mat;
   double *sm_evals;
   double **sm_evecs;
   int i, j, ij, k, I;
   double **G, **d;
   double *lambda, **alpha, **f ;
   double *converged_root;
   double **m_lambda, ***m_alpha;
   double tval, *dvecnorm;
   int converged=0, iter=1;
   int iter2=0; /* iterations since last collapse */
   double *lastroot;
   int lse_do=0, last_lse_collapse_num=-Parameters_->lse_collapse, collapse_num=0;
   double lse_tolerance=Parameters_->lse_tolerance;
   double **sigma_overlap, ***Mmatrix;
   int *Lvec;

   /* check parameters */
   if (evecs == NULL || evals == NULL) {
     outfile->Printf("(sem): passed uncallocated pointers for evecs or evals\n") ;
      return ;
      }


   for (I=0; I<N; I++)
      A[I][I] -= CalcInfo_->edrc;


   /* make space for temp vector */
   tmp_vec = init_array(N);
   lastroot = init_array(N);
   converged_root = init_array(M);
   dvecnorm = init_array(M);

   /* allocate other arrays with ~fixed dimensions during iteration */
   d = init_matrix(M, N);    /* M and N are both fixed */
   f = init_matrix(M, N);
   G = init_matrix(maxnvect, maxnvect);
   tmp_mat = init_matrix(maxnvect, N);
   jnk = init_matrix(maxnvect, N);
   lambda = init_array(maxnvect);
   alpha = init_matrix(maxnvect, maxnvect);
   sigma_overlap = init_matrix(maxnvect, maxnvect);
   Mmatrix = (double ***) malloc (sizeof(double **) * M);
   for (i=0; i<M; i++)
      Mmatrix[i] = init_matrix(maxnvect, maxnvect);
   m_lambda = init_matrix(M, maxnvect);
   m_alpha = (double ***) malloc (sizeof(double **) * maxnvect);
   for (i=0; i<maxnvect; i++) {
      m_alpha[i] = init_matrix(maxnvect, maxnvect);
      }
   Lvec = init_int_array(maxnvect);


   /* ITERATE */
   while (!converged && iter <= maxiter) {

      Lvec[iter2] = L;
      /* form G matrix */
      mmult(b, 0, A, 0, tmp_mat, 0, L, N, N, 0); /* tmp = B * A    */
      mmult(tmp_mat, 0, b, 1, G, 0, L, N, L, 0); /* G = tmp * B(T) */

      /* solve the L x L eigenvalue problem G a = lambda a for M roots */
      sq_rsp(L, L, G, lambda, 1, alpha, 1E-14);

      if (N<100 && print_ >=3) {
        outfile->Printf("\n b matrix\n");
        print_mat(b,L,N,"outfile");
        outfile->Printf("\n sigma matrix\n");
        print_mat(tmp_mat,L,N,"outfile");
        outfile->Printf("\n G matrix (%d)\n", iter-1);
        print_mat(G,L,L,"outfile");
        outfile->Printf("\n Eigenvectors and eigenvalues of G matrix (%d)\n", iter-1);
        eivout(alpha, lambda, L, L, "outfile");
        }

      lse_do = 0;
      if (Parameters_->lse && (maxnvect-L <= M*Parameters_->collapse_size) && L>2 &&
         (lse_tolerance > fabs(lambda[0]-lastroot[0])) && iter>=3 &&
         ((collapse_num-last_lse_collapse_num)>= Parameters_->lse_collapse))
        lse_do = 1;
      if (lse_do) {
        /* Form sigma_overlap matrix */
        zero_mat(sigma_overlap,maxnvect,maxnvect);
        mmult(b, 0, A, 0, tmp_mat, 0, L, N, N, 0);
        mmult(tmp_mat, 0, tmp_mat, 1, sigma_overlap, 0, L, N, L, 0);

       /* Form Mij matrix */
       for (k=0; k<M; k++) {
          zero_mat(Mmatrix[k],maxnvect,maxnvect);
          for (i=0; i<L; i++) {
             for (j=i; j<L; j++) {
                Mmatrix[k][i][j] = Mmatrix[k][j][i] =  sigma_overlap[i][j]
                             -2.0 * lambda[k] * G[i][j];
                if (i==j) Mmatrix[k][i][i] += pow(lambda[k],2.0);
                }
             }
          } /* end loop over k (nroots) */

         if (print_ > 2) {
           outfile->Printf( "\nsigma_overlap matrix (%2d) = \n", iter-1);
           print_mat(sigma_overlap, L, L, "outfile");

           for (k=0; k<M; k++) {
              outfile->Printf( "\nM matrix (%2d) for root %d = \n", iter, k);
              print_mat(Mmatrix[k], L, L, "outfile");
              outfile->Printf( "\n");
              }
           }

        /* solve the L x L eigenvalue problem M a = lambda a for M roots */
       for (k=0; k<M; k++) {
          sq_rsp(L, L, Mmatrix[k], m_lambda[k], 1, m_alpha[k], 1.0E-14);
          if (print_ > 2) {
            outfile->Printf( "\n M eigenvectors and eigenvalues root %d:\n",k);
            eivout(m_alpha[k], m_lambda[k], L, L, "outfile");
            }
          }

        } /* end if lse_do */

      if ((Parameters_->collapse_size>0) && (iter2-Parameters_->collapse_size+1 > 0)
         && (Lvec[iter2-Parameters_->collapse_size+1]+M*Parameters_->collapse_size
         > maxnvect) && iter!=maxiter) {

        collapse_num++;
        if (lse_do) last_lse_collapse_num = collapse_num;
        /* copy ci vector into d matrix */
        zero_mat(d, M, N);
        if (lse_do)
          for (k=0; k<M; k++)
             for (i=0; i<L; i++)
                for (I=0; I<N; I++)
                   d[k][I] += m_alpha[k][i][0] * b[i][I];

        else
          for (k=0; k<M; k++)
             for (i=0; i<L; i++)
                for (I=0; I<N; I++)
                   d[k][I] += alpha[i][k] * b[i][I];

        /* copy d matrix to end of b matrix */
        for (i=0; i<M; i++)
           for (I=0; I<N; I++)
              b[maxnvect-1-i][I] = d[i][I];

        /* reorder b matrix pointers */
        for (i=0; i<L; i++) jnk[i] = b[i];
        for (i=0; i<L; i++) b[i] = jnk[maxnvect-1-i];

        L = M;
        iter2 = 0;

        /* zero out old parts of b matrix */
        for (i=L; i<maxnvect; i++) zero_arr(b[i], N);

        /* reform G matrix */
        mmult(b, 0, A, 0, tmp_mat, 0, L, N, N, 0); /* tmp = B * A    */
        mmult(tmp_mat, 0, b, 1, G, 0, L, N, L, 0); /* G = tmp * B(T) */

        /* solve the L x L eigenvalue problem G a = lambda a for M roots */
        sq_rsp(L, L, G, lambda, 1, alpha, 1E-14);

        if (N<100 && print_ >= 3) {
        outfile->Printf(" Reformed G matrix (%d)\n",iter-1);
        print_mat(G,L,L,"outfile");
        outfile->Printf("\n");
        }

        if (lse_do) outfile->Printf(" Least Squares Extrapolation\n");
        outfile->Printf(" Collapse Davidson subspace to %d vectors\n", L);
       } /* end collapse */

      /* form the d part of the correction vector */
      zero_mat(d, M, N);
      for (k=0; k<M; k++) {
         for (i=0; i<L; i++) {
            mmult(A,0,&(b[i]),1,&(tmp_vec),1,N,N,1,0); /* tmp=A*b[i] */
            for (I=0; I<N; I++) {
               d[k][I] += alpha[i][k] * (tmp_vec[I] - lambda[k] * b[i][I]);
               }
            }
         }

      if (N<100 && print_ >= 3) {
        outfile->Printf(" D vectors for iter (%d)\n",iter-1);
        print_mat(d,M,N,"outfile");
        }

      /* check for convergence */
      converged = 1;
      for (i=0; i<M; i++) {
         dot_arr(d[i], d[i], N, &tval);
         tval = sqrt(tval);
         dvecnorm[i] = tval;
         if (dvecnorm[i] <= conv_rms && fabs(lambda[i] - lastroot[i]) <= conv_e) converged_root[i] = 1;
         else {
          converged_root[i] = 0;
          converged = 0;
          }
         outfile->Printf( "Iter %3d  Root %d = %13.9lf",
            iter-1, i+1, (lambda[i] + offst));
         outfile->Printf( "    Delta_E %10.3E   Delta_C %10.3E %c\n",
            lambda[i] - lastroot[i], tval, converged_root[i] ? 'c' : ' ');

         }

      if (M>1) {
        outfile->Printf( "\n");

        }

      if (converged || iter == maxiter) {
         for (i=0; i<M; i++) {
            evals[i] = lambda[i];
            for (j=0; j<L; j++) {
               tval = alpha[j][i];
               for (I=0; I<N; I++)
                  evecs[i][I] += tval * b[j][I];
               }
            }
         break;
         }
      else {
         for (i=0; i<M; i++) lastroot[i] = lambda[i];
         }

      /* form the correction vector and normalize */
      for (k=0; k<M; k++) {
         for (I=0; I<N; I++) {
            tval = lambda[k] - A[I][I];
            /* make sure denom isn't 0.  If so, make some arbitrary val  */
            /* It might be interesting to figure the best way to do this */

            /* previous way to do it
            if (fabs(tval) < MIN_F_DENOM) tval = 0.1;
            f[k][I] = d[k][I] / tval;
            */

            /* the way GUGA does it */
            if (fabs(tval) < 1.0E-8) f[k][I] = 0.0;
            else f[k][I] = d[k][I] / tval;
            }
         }
      normalize(f, M, N);

      /* Schmidt orthog and append f's to b */
      for (i=0; i<M; i++)
         if (converged_root[i] == 0)
           if (schmidt_add(b, L, N, f[i])) L++;
      outfile->Printf(" Number of b vectors = %d\n", L);

      if (L > maxnvect) {
         std::string str = "(test_sem): L(";
         str += std::to_string( L) ;
         str += ") > maxnvect(";
         str += std::to_string( maxnvect) ;
         str += ")! Aborting!";
         throw PsiException(str,__FILE__,__LINE__);
         }

      /* Again Schmidt orthog b's (minimize numerical error) */
      /* Doesn't this mess up the sigma vectors slightly */
      schmidt(b, L, N, "outfile");
      iter++ ;
      iter2++;
      }

   *vu = L;

   free(lambda);
   free(tmp_vec);
   free_matrix(d, M);
   free_matrix(f, M);
   free_matrix(b, maxnvect);
   free_matrix(G, maxnvect);
   free_matrix(tmp_mat, maxnvect);
   free_matrix(alpha, maxnvect);
   free_matrix(sigma_overlap, maxnvect);
   for (i=0; i<M; i++) free_matrix(Mmatrix[i], maxnvect);
   free(Mmatrix);
   for (i=0; i<maxnvect; i++) free_matrix(m_alpha[i], maxnvect);
   free(m_alpha);
}

}} // namespace psi::detci
