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
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

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
FILE *outfile ;

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


   fprintf(outfile, "Ok, the eigenvectors are sideways!\n");
   eivout(evecs, evals, 4, 50, outfile);
   fprintf(outfile, "\nused %d expansion vectors\n", used);

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
void sem_test(double **A, int N, int M, int L, double **evecs, double *evals, 
      double **b, double conv_e, double conv_rms, int maxiter, double offst, 
      int *vu, int maxnvect, FILE *outfile)
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
   int lse_do=0, last_lse_collapse_num=-Parameters.lse_collapse, collapse_num=0;
   double lse_tolerance=Parameters.lse_tolerance;
   double **sigma_overlap, ***Mmatrix;
   int *Lvec;

   /* check parameters */
   if (evecs == NULL || evals == NULL) {
      printf("(sem): passed uncallocated pointers for evecs or evals\n") ;
      return ;
      }


   for (I=0; I<N; I++)
      A[I][I] -= CalcInfo.efzc;
   

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
 
      if (N<100 && Parameters.print_lvl >=3) {
        fprintf(outfile,"\n b matrix\n");
        print_mat(b,L,N,outfile);
        fprintf(outfile,"\n sigma matrix\n");
        print_mat(tmp_mat,L,N,outfile); 
        fprintf(outfile,"\n G matrix (%d)\n", iter-1);
        print_mat(G,L,L,outfile);
        fprintf(outfile,"\n Eigenvectors and eigenvalues of G matrix (%d)\n", iter-1);
        eivout(alpha, lambda, L, L, outfile);
        }

      lse_do = 0;
      if (Parameters.lse && (maxnvect-L <= M*Parameters.collapse_size) && L>2 &&
         (lse_tolerance > fabs(lambda[0]-lastroot[0])) && iter>=3 &&
         ((collapse_num-last_lse_collapse_num)>= Parameters.lse_collapse)) 
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

         if (Parameters.print_lvl > 2) {
           fprintf(outfile, "\nsigma_overlap matrix (%2d) = \n", iter-1);
           print_mat(sigma_overlap, L, L, outfile);

           for (k=0; k<M; k++) {
              fprintf(outfile, "\nM matrix (%2d) for root %d = \n", iter, k);
              print_mat(Mmatrix[k], L, L, outfile);
              fprintf(outfile, "\n");
              }
           } 

        /* solve the L x L eigenvalue problem M a = lambda a for M roots */
       for (k=0; k<M; k++) {
          sq_rsp(L, L, Mmatrix[k], m_lambda[k], 1, m_alpha[k], 1.0E-14);
          if (Parameters.print_lvl > 2) {
            fprintf(outfile, "\n M eigenvectors and eigenvalues root %d:\n",k);
            eivout(m_alpha[k], m_lambda[k], L, L, outfile);
            }
          }

        } /* end if lse_do */ 

      if ((Parameters.collapse_size>0) && (iter2-Parameters.collapse_size+1 > 0)
         && (Lvec[iter2-Parameters.collapse_size+1]+M*Parameters.collapse_size 
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
        
        if (N<100 && Parameters.print_lvl >= 3) {
        fprintf(outfile," Reformed G matrix (%d)\n",iter-1);
        print_mat(G,L,L,outfile);
        fprintf(outfile,"\n");
        }

        if (lse_do) fprintf(outfile," Least Squares Extrapolation\n");
        fprintf(outfile," Collapse Davidson subspace to %d vectors\n", L);
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

      if (N<100 && Parameters.print_lvl >= 3) {
        fprintf(outfile," D vectors for iter (%d)\n",iter-1);
        print_mat(d,M,N,outfile);
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
         fprintf(outfile, "Iter %3d  Root %d = %13.9lf",
            iter-1, i+1, (lambda[i] + offst));
         fprintf(outfile, "    Delta_E %10.3E   Delta_C %10.3E %c\n",
            lambda[i] - lastroot[i], tval, converged_root[i] ? 'c' : ' ');
         fflush(outfile);
         }

      if (M>1) {
        fprintf(outfile, "\n");
        fflush(outfile);
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
      fprintf(outfile," Number of b vectors = %d\n", L);

      if (L > maxnvect) {
         std::string str = "(test_sem): L(";
         str += static_cast<std::ostringstream*>( &(std::ostringstream() << L) )->str();
         str += ") > maxnvect(";
         str += static_cast<std::ostringstream*>( &(std::ostringstream() << maxnvect) )->str();
         str += ")! Aborting!";
         throw PsiException(str,__FILE__,__LINE__);
         }

      /* Again Schmidt orthog b's (minimize numerical error) */
      /* Doesn't this mess up the sigma vectors slightly */
      schmidt(b, L, N, outfile); 
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
}

}} // namespace psi::detci

