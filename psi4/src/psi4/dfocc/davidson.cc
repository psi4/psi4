/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

#define BIGNUM 1E100
//#define maxiter 1000
#define maxiter cc_maxiter

/*!
** davidson(): Computes the lowest few eigenvalues and eigenvectors of a
** symmetric matrix, A, using the Davidson-Liu algorithm.  
**
** The matrix must be small enough to fit entirely in core.  This algorithm 
** is useful if one is interested in only a few roots of the matrix
** rather than the whole spectrum.
**
** NB: This implementation will keep up to eight guess vectors for each
** root desired before collapsing to one vector per root.  In
** addition, if smart_guess=1 (the default), guess vectors are
** constructed by diagonalization of a sub-matrix of A; otherwise,
** unit vectors are used.
**
** TDC, July-August 2002
** UB, November 2017
**
** \param A      = matrix to diagonalize
** \param N      = dimension of A
** \param M      = number of roots desired
** \param eps    = eigenvalues
** \param v      = eigenvectors
** \param cutoff = tolerance for convergence of eigenvalues
** \param print  = Boolean for printing additional information
**
** Returns: number of converged roots
*/

int DFOCC::davidson(SharedTensor2d A, int M, SharedTensor1d eps, SharedTensor2d v, double cutoff, int print) {
  int i, j, k, L, I;
  double minimum;
  int min_pos, numf, iter, converged, maxdim, skip_check;
  int init_dim;
  int smart_guess = 1;
  double norm, denom, diff;

  SharedTensor2d B, Bnew, sigma, G, alpha, Res;
  SharedTensor1d diag, lambda, lambda_old;
  SharedTensor1i small2big, conv;

  // get dimension of matrix A
  int N = A->dim1();

  // define maximum dimension
  maxdim = 8 * M;

  // Current set of guess vectors, stored by row
  B = std::make_shared<Tensor2d>("B matrix", maxdim, N);
  // Guess vectors formed from old vectors, stored by row
  Bnew = std::make_shared<Tensor2d>("New B matrix", M, N);
  // Sigma vectors, stored by column
  sigma = std::make_shared<Tensor2d>("sigma matrix", N, maxdim);
  // Davidson mini-Hamitonian
  G = std::make_shared<Tensor2d>("G matrix", maxdim, maxdim);
  // Residual eigenvectors, stored by row
  Res = std::make_shared<Tensor2d>("f matrix", maxdim, N);
  // Eigenvectors of G
  alpha = std::make_shared<Tensor2d>("alpha matrix", maxdim, maxdim);
  // Eigenvalues of G
  lambda = std::make_shared<Tensor1d>("lambda vector", maxdim);
  // Approximate roots from previous iteration
  lambda_old = std::make_shared<Tensor1d>("lambda vector", maxdim);

  //===================================
  //======== Initial guess ============
  //===================================
  // Use eigenvectors of a sub-matrix as initial guesses
  if (smart_guess) { 

    // init dim
    if (N > 7*M) init_dim = 7*M;
    else init_dim = M;

    // Diagonals of A
    diag = std::make_shared<Tensor1d>("A diag vector", N);
    // Addressing
    small2big = std::make_shared<Tensor1i>("A diag vector", 7*M);

    // diag of A
    for (i=0; i < N; i++) diag->set(i, A->get(i,i)); 

    // Form addressing array 
    for (i=0; i < init_dim; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if(diag->get(j) < minimum) {	
                 minimum = diag->get(j);
	         min_pos = j; 
		 small2big->set(i,j);
	      } 
         }
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }

    // Form G
    for (i=0; i < init_dim; i++) {
	 int ii = small2big->get(i);
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      G->set(i, j, A->get(ii,jj));

         }
    }

    // Diagonalize G
    G->diagonalize(init_dim, alpha, lambda, 1e-12, true);

    // Form B 
    for (i=0; i < init_dim; i++) {
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      B->set(i, jj, alpha->get(j,i));
	 }
    }

    diag.reset();
    small2big.reset();
  }// end if(smart_guess)

  // Use unit vectors as initial guesses
  else { 
    // Diagonals of A
    diag = std::make_shared<Tensor1d>("A diag vector", N);
    for (i=0; i < N; i++) diag->set(i, A->get(i,i)); 

    // Form B
    for (i=0; i < M; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if (diag->get(j) < minimum) { 
		  minimum = diag->get(j); 
	          min_pos = j; 
	      }
         }

	 B->set(i, min_pos, 1.0);
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }
    diag.reset();
  } // end else


  //===================================
  //======== Loop =====================
  //===================================
  // start
  L = init_dim;
  iter = 0;
  converged = 0;

  // boolean array for convergence of each root
  conv = std::make_shared<Tensor1i>("conv vector", M);

  // Head of Loop
  while(converged < M && iter < maxiter) {

    skip_check = 0;
    if(print) printf("\niter = %d\n", iter); 

    // Form G
    // sigma = A*B'
    sigma->gemm(false, true, A, B, 1.0, 0.0);

    // G = B*sigma
    G->gemm(false, false, B, sigma, 1.0, 0.0);

    // Diagonalize G
    G->diagonalize(L, alpha, lambda, 1e-12, true);

    // Form preconditioned residue vectors
    for (k=0; k < M; k++) {
         for (I=0; I < N; I++) {
	      Res->set(k,I, 0.0);
	      double value = 0.0;
	      for(i=0; i < L; i++) {
		  value += alpha->get(i,k) * (sigma->get(I,i) - lambda->get(k) * B->get(i,I));
	      }
              Res->add(k,I,value);
	      denom = lambda->get(k) - A->get(I,I);
	      if (fabs(denom) > 1e-6) Res->set(k,I, Res->get(k,I)/denom);
	      else Res->set(k,I,0.0);
         }
    }

    // Normalize each residual
    for (k=0; k < M; k++) {
         norm = 0.0;
         for (I=0; I < N; I++) {
	     norm += Res->get(k,I) * Res->get(k,I);
         }
      
	 norm = std::sqrt(norm);
         for (I=0; I < N; I++) {
	      if (norm > 1e-6) Res->set(k,I, Res->get(k,I)/norm);
	      else Res->set(k,I,0.0);
         }
    }

    // Schmidt orthogonalize the Res[k] against the set of B[i] and add new vectors
    for (k=0,numf=0; k < M; k++) {
         SharedTensor1d Rk = std::make_shared<Tensor1d>("Res[k] vector", N);
	 Rk->row_vector(Res,k);
         if (B->gs_add(L,Rk)) { 
             L++; 
	     numf++; 
	 }
    }

    // If L is close to maxdim, collapse to one guess per root
    if (maxdim - L < M) {
        if (print) {
	    printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
	    printf("Collapsing eigenvectors.\n");
        }

	// form Bnew
        for (i=0; i < M; i++) {
	     for (k=0; k < N; k++) {
		  double sum = 0.0;
	          for (j=0; j < L; j++) {
	               sum += alpha->get(j,i) * B->get(j,k);
	          } 
		  Bnew->set(i,k,sum);
	     }
        }

        // Copy new vectors into place
        Bnew->zero();
        for (int i=0; i < M; i++) {
	     for (int k=0; k < N; k++) {
                  B->set(i,k, Bnew->get(i,k));
             }
        }
        skip_check = 1;
        L = M;
    } // end if

    // check convergence on all roots
    if (!skip_check) {
        converged = 0;
	conv->zero();
        if (print) {
	    printf("Root      Eigenvalue       Delta  Converged?\n");
	    printf("---- -------------------- ------- ----------\n");
        }  

        for (k=0; k < M; k++) {
	     diff = std::fabs(lambda->get(k) - lambda_old->get(k));
	     if (diff < cutoff) {
		 conv->set(k,1);
	         converged++;
	     }
	     lambda_old->set(k,lambda->get(k));
	     if (print) {
	         printf("%3d  %20.14f %4.3e    %1s\n", k, lambda->get(k), diff, conv->get(k) == 1 ? "Y" : "N");
	     }
        }
    }// end if (!skip_check)

    iter++;
  }// end of loop

  // Generate final eigenvalues and eigenvectors
  if (converged == M) {
      for (i=0; i < M; i++) {
           eps->set(i,lambda->get(i));
	   for (I=0; I < N; I++) {
	        double sum = 0.0;
                for (j=0; j < L; j++) {
		     sum += alpha->get(j,i) * B->get(j,I);
	        }  
	        v->add(I,i,sum);
           }
      }
      if (print) printf("Davidson algorithm converged in %d iterations.\n", iter);
  }

  conv.reset();
  B.reset();
  Bnew.reset();
  sigma.reset();
  G.reset();
  Res.reset();
  alpha.reset();
  lambda.reset();
  lambda_old.reset();

  return converged;
}

/*!
** cis_davidson(): Computes the lowest few eigenvalues and eigenvectors of a
** symmetric matrix, A, using the Davidson-Liu algorithm.  
**
** The matrix must be small enough to fit entirely in core.  This algorithm 
** is useful if one is interested in only a few roots of the matrix
** rather than the whole spectrum.
**
** NB: This implementation will keep up to eight guess vectors for each
** root desired before collapsing to one vector per root.  In
** addition, if smart_guess=1 (the default), guess vectors are
** constructed by diagonalization of a sub-matrix of A; otherwise,
** unit vectors are used.
**
** TDC, July-August 2002
** UB, November 2018
**
** \param A      = matrix to diagonalize
** \param N      = dimension of A
** \param M      = number of roots desired
** \param eps    = eigenvalues
** \param v      = eigenvectors
** \param cutoff = tolerance for convergence of eigenvalues
**
** Returns: number of converged roots
*/

int DFOCC::cis_davidson(int M, SharedTensor1d eps, SharedTensor2d v, double cutoff) {
  int L;
  double minimum;
  int min_pos, numf, iter, converged, maxdim, skip_check;
  int init_dim;
  int smart_guess = 1;
  double norm, denom, diff;

  SharedTensor2d B, Bnew, sigma, G, alpha, Res, G2, alpha2;
  SharedTensor1d diag, lambda, lambda_old, lambda2;
  SharedTensor1i small2big, conv;

  // get dimension of matrix A
  int N = naoccA * navirA;

  // define maximum dimension
  maxdim = 8 * M;

  // malloc
  auto Hdiag = std::make_shared<Tensor1d>("H <IA>", N);
  auto Cci = std::make_shared<Tensor1d>("C <IA>", N);
  auto Sci = std::make_shared<Tensor1d>("Sigma <IA>", N);

  // diagonal H
  if (cis_alg == "MO_BASIS") cis_diagonal(Hdiag);
  else if (cis_alg == "AO_BASIS") cis_diagonal_approx(Hdiag);

  // Current set of guess vectors, stored by row
  B = std::make_shared<Tensor2d>("B matrix", maxdim, N);
  // Guess vectors formed from old vectors, stored by row
  Bnew = std::make_shared<Tensor2d>("New B matrix", M, N);
  // Sigma vectors, stored by column
  sigma = std::make_shared<Tensor2d>("sigma matrix", N, maxdim);
  // Davidson mini-Hamitonian
  G = std::make_shared<Tensor2d>("G matrix", maxdim, maxdim);
  // Residual eigenvectors, stored by row
  Res = std::make_shared<Tensor2d>("f matrix", maxdim, N);
  // Eigenvectors of G
  alpha = std::make_shared<Tensor2d>("alpha matrix", maxdim, maxdim);
  // Eigenvalues of G
  lambda = std::make_shared<Tensor1d>("lambda vector", maxdim);
  // Approximate roots from previous iteration
  lambda_old = std::make_shared<Tensor1d>("old lambda vector", maxdim);

  //===================================
  //======== Initial guess ============
  //===================================
  // Use eigenvectors of a sub-matrix as initial guesses
  /*
  if (smart_guess) { 

    // init dim
    if (N > 7*M) init_dim = 7*M;
    else init_dim = M;

    // Diagonals of A
    diag = std::make_shared<Tensor1d>("A diag vector", N);
    // Addressing
    small2big = std::make_shared<Tensor1i>("A diag vector", 7*M);

    // diag of A
    for (i=0; i < N; i++) diag->set(i, A->get(i,i)); 

    // Form addressing array 
    for (i=0; i < init_dim; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if(diag->get(j) < minimum) {	
                 minimum = diag->get(j);
	         min_pos = j; 
		 small2big->set(i,j);
	      } 
         }
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }

    // Form G
    for (i=0; i < init_dim; i++) {
	 int ii = small2big->get(i);
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      G->set(i, j, A->get(ii,jj));

         }
    }

    // Diagonalize G
    G->diagonalize(init_dim, alpha, lambda, 1e-12, true);

    // Form B 
    for (i=0; i < init_dim; i++) {
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      B->set(i, jj, alpha->get(j,i));
	 }
    }

    diag.reset();
    small2big.reset();
  }// end if(smart_guess)
  */

  // Use unit vectors as initial guesses
  //else { 
    // Diagonals of A
    diag = std::make_shared<Tensor1d>("A diag vector", N);
    diag->axpy(Hdiag, 1.0);

    // Form B
    for (int i=0; i < M; i++) {
         //outfile->Printf("\troot: %1d\n",i);
         minimum = diag->get(0);
         min_pos = 0;
         for (int j=1; j < N; j++) {
	      if (diag->get(j) < minimum) { 
		  minimum = diag->get(j); 
	          min_pos = j; 
		  //outfile->Printf("\tmin_pos: %2d\n",min_pos);
	      }
         }

	 B->set(i, min_pos, 1.0);
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }
    diag.reset();
    //B->print();
  //} // end else


  //===================================
  //======== Loop =====================
  //===================================
  // start
  //L = init_dim;
  L = M;
  iter = 0;
  converged = 0;

  // boolean array for convergence of each root
  conv = std::make_shared<Tensor1i>("conv vector", M);

  // Head of Loop
  while(converged < M && iter < maxiter) {

    skip_check = 0;
    outfile->Printf("\n\tDavidson iter: %1d\n",iter);

    // Get Coefficients
    for (int i=0; i < L; i++) {
	 // Get C for the current root
         for (int J=0; J < N; J++) Cci->set(J, B->get(i,J));
	 
	 // Form sigma 
	 if (cis_alg == "MO_BASIS") cis_sigma(Cci, Sci);
	 else if (cis_alg == "AO_BASIS") cis_sigma_ao(Cci, Sci);

	 // Put into overall sigma
         for (int J=0; J < N; J++) sigma->set(J, i, Sci->get(J));
    }
    //sigma->print();

    // Form G
    // G = B*sigma
    G->gemm(false, false, B, sigma, 1.0, 0.0);

    // Diagonalize G
    G->diagonalize(L, alpha, lambda, 1e-12, true);

    // Form preconditioned residue vectors
    for (int k=0; k < M; k++) {
         for (int I=0; I < N; I++) {
	      Res->set(k,I, 0.0);
	      double value = 0.0;
	      for(int i=0; i < L; i++) {
		  value += alpha->get(i,k) * (sigma->get(I,i) - lambda->get(k) * B->get(i,I));
	      }
              Res->add(k,I,value);
	      denom = lambda->get(k) - Hdiag->get(I);
	      if (fabs(denom) > 1e-6) Res->set(k,I, Res->get(k,I)/denom);
	      else Res->set(k,I,0.0);
         }
    }

    // Normalize each residual
    for (int k=0; k < M; k++) {
         norm = 0.0;
         for (int I=0; I < N; I++) {
	     norm += Res->get(k,I) * Res->get(k,I);
         }
      
	 norm = std::sqrt(norm);
         for (int I=0; I < N; I++) {
	      if (norm > 1e-6) Res->set(k,I, Res->get(k,I)/norm);
	      else Res->set(k,I,0.0);
         }
    }

    // Schmidt orthogonalize the Res[k] against the set of B[i] and add new vectors
    for (int k=0,numf=0; k < M; k++) {
        SharedTensor1d Rk = std::make_shared<Tensor1d>("Res[k] vector", N);
        Rk->row_vector(Res,k);
        if (B->gs_add(L,Rk)) { 
            L++; 
            numf++; 
        }
    }

    // If L is close to maxdim, collapse to one guess per root
    if (maxdim - L < M) {
	outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
	outfile->Printf("Collapsing eigenvectors.\n");

	// form Bnew
        Bnew->zero();
        for (int i=0; i < M; i++) {
	     for (int k=0; k < N; k++) {
		  double sum = 0.0;
	          for (int j=0; j < L; j++) {
	               sum += alpha->get(j,i) * B->get(j,k);
	          } 
                  Bnew->set(i,k,sum);
	     }
        }

        // Copy new vectors into place
        for (int i=0; i < M; i++) {
	     for (int k=0; k < N; k++) {
                  B->set(i,k, Bnew->get(i,k));
             }
        }
        skip_check = 1;
        L = M;
    } // end if

    // check convergence on all roots
    if (!skip_check) {
        converged = 0;
	conv->zero();
	outfile->Printf("\tRoot      Eigenvalue       Delta  Converged?\n");
	outfile->Printf("\t---- -------------------- ------- ----------\n");

        for (int k=0; k < M; k++) {
	     diff = std::fabs(lambda->get(k) - lambda_old->get(k));
	     if (diff < cutoff) {
		 conv->set(k,1);
	         converged++;
	     }
	     lambda_old->set(k,lambda->get(k));
	     outfile->Printf("\t%3d  %20.14f %4.3e    %1s\n", k+1, lambda->get(k), diff, conv->get(k) == 1 ? "Y" : "N");
        }
    }// end if (!skip_check)

    iter++;
  }// end of loop

  // Generate final eigenvalues and eigenvectors
  if (converged == M) {
      for (int i=0; i < M; i++) {
           eps->set(i,lambda->get(i));
	   for (int I=0; I < N; I++) {
	        double sum = 0.0;
                for (int j=0; j < L; j++) {
		     sum += alpha->get(j,i) * B->get(j,I);
	        }  
	        v->add(I,i,sum);
           }
      }
      outfile->Printf("\n\tDavidson algorithm converged in %d iterations.\n", iter);
  }

  conv.reset();
  B.reset();
  Bnew.reset();
  sigma.reset();
  G.reset();
  Res.reset();
  alpha.reset();
  lambda.reset();
  lambda_old.reset();
  Hdiag.reset();
  Cci.reset();
  Sci.reset();

  return converged;
}

/*!
** qdpt_davidson(): Computes the lowest few eigenvalues and eigenvectors of a
** symmetric matrix, A, using the Davidson-Liu algorithm.  
**
** The matrix must be small enough to fit entirely in core.  This algorithm 
** is useful if one is interested in only a few roots of the matrix
** rather than the whole spectrum.
**
** NB: This implementation will keep up to eight guess vectors for each
** root desired before collapsing to one vector per root.  In
** addition, if smart_guess=1 (the default), guess vectors are
** constructed by diagonalization of a sub-matrix of A; otherwise,
** unit vectors are used.
**
** UB, November 2018
**
** \param A      = matrix to diagonalize
** \param N      = dimension of A
** \param M      = number of roots desired
** \param eps    = eigenvalues
** \param v      = eigenvectors
** \param cutoff = tolerance for convergence of eigenvalues
**
** Returns: number of converged roots
*/

int DFOCC::qdpt_davidson(int M, SharedTensor1d eps, SharedTensor2d v, double cutoff) {
  int L;
  double minimum;
  int min_pos, numf, iter, converged, maxdim, skip_check;
  int init_dim;
  int smart_guess = 1;
  double norm, denom, diff;

  SharedTensor2d B, Bnew, sigma, G, alpha, Res, G2, alpha2;
  SharedTensor1d diag, lambda, lambda_old, lambda2;
  SharedTensor1i small2big, conv;

  // get dimension of matrix A
  int N = nconfig;

  // define maximum dimension
  maxdim = 8 * M;

  // malloc
  auto Hdiag = std::make_shared<Tensor1d>("Diagonal H", N);
  auto Cci = std::make_shared<Tensor1d>("CI Coeff", N);
  auto Sci = std::make_shared<Tensor1d>("Sigma for a root", N);

  // diagonal H
  if (wfn_type_ == "DF-CIS") cis_diagonal(Hdiag);
  //if (wfn_type_ == "QDPT2") qdpt2_diagonal(Hdiag);
  //else if (wfn_type_ == "FCI" || wfn_type_ == "CAS") fci_diagonal(Hdiag);

  // Current set of guess vectors, stored by row
  B = std::make_shared<Tensor2d>("B matrix", maxdim, N);
  // Guess vectors formed from old vectors, stored by row
  Bnew = std::make_shared<Tensor2d>("New B matrix", M, N);
  // Sigma vectors, stored by column
  sigma = std::make_shared<Tensor2d>("sigma matrix", maxdim, N);
  // Davidson mini-Hamitonian
  G = std::make_shared<Tensor2d>("G matrix", maxdim, maxdim);
  // Residual eigenvectors, stored by row
  Res = std::make_shared<Tensor2d>("f matrix", maxdim, N);
  // Eigenvectors of G
  alpha = std::make_shared<Tensor2d>("alpha matrix", maxdim, maxdim);
  // Eigenvalues of G
  lambda = std::make_shared<Tensor1d>("lambda vector", maxdim);
  // Approximate roots from previous iteration
  lambda_old = std::make_shared<Tensor1d>("old lambda vector", maxdim);

  //===================================
  //======== Initial guess ============
  //===================================
  // Use unit vectors as initial guesses
  //else { 
    // Diagonals of A
    diag = std::make_shared<Tensor1d>("A diag vector", N);
    diag->axpy(Hdiag, 1.0);

    // Form B
    for (int i=0; i < M; i++) {
         //outfile->Printf("\troot: %1d\n",i);
         minimum = diag->get(0);
         min_pos = 0;
         for (int j=1; j < N; j++) {
	      if (diag->get(j) < minimum) { 
		  minimum = diag->get(j); 
	          min_pos = j; 
		  //outfile->Printf("\tmin_pos: %2d\n",min_pos);
	      }
         }

	 B->set(i, min_pos, 1.0);
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }
    diag.reset();
    //B->print();
  //} // end else

  //===================================
  //======== Loop =====================
  //===================================
  // start
  //L = init_dim;
  L = M;
  iter = 0;
  converged = 0;

  // boolean array for convergence of each root
  conv = std::make_shared<Tensor1i>("conv vector", M);

  // Head of Loop
  while(converged < M && iter < maxiter) {

    skip_check = 0;
    outfile->Printf("\n\tDavidson iter: %1d\n",iter);

    // Get Coefficients
    for (int i=0; i < L; i++) {
	 // Get C for the current root
         for (int J=0; J < N; J++) Cci->set(J, B->get(i,J));
	 
	 // Form sigma 
	 if (wfn_type_ == "DF-CIS") cis_sigma(Cci, Sci);
	 //if (wfn_type_ == "QDPT2") qdpt2_sigma(Hdiag, Cci, Sci);
	 //else if (wfn_type_ == "FCI" || wfn_type_ == "CAS")  fci_sigma(Cci, Sci);

	 // Put into overall sigma
         for (int J=0; J < N; J++) sigma->set(i, J, Sci->get(J));
    }
    //sigma->print();

    // Form G
    // G = sigma*B'
    G->gemm(false, true, sigma, B, 1.0, 0.0);

    // Diagonalize G
    G->diagonalize(L, alpha, lambda, 1e-12, true);

    // Form preconditioned residue vectors
    for (int k=0; k < M; k++) {
         for (int I=0; I < N; I++) {
	      Res->set(k,I, 0.0);
	      double value = 0.0;
	      for(int i=0; i < L; i++) {
		  value += alpha->get(i,k) * (sigma->get(i,I) - lambda->get(k) * B->get(i,I));
	      }
              Res->add(k,I,value);
	      denom = lambda->get(k) - Hdiag->get(I);
	      if (fabs(denom) > 1e-6) Res->set(k,I, Res->get(k,I)/denom);
	      else Res->set(k,I,0.0);
         }
    }

    // Normalize each residual
    for (int k=0; k < M; k++) {
         norm = 0.0;
         for (int I=0; I < N; I++) {
	     norm += Res->get(k,I) * Res->get(k,I);
         }
      
	 norm = std::sqrt(norm);
         for (int I=0; I < N; I++) {
	      if (norm > 1e-6) Res->set(k,I, Res->get(k,I)/norm);
	      else Res->set(k,I,0.0);
         }
    }

    // Schmidt orthogonalize the Res[k] against the set of B[i] and add new vectors
    for (int k=0,numf=0; k < M; k++) {
         SharedTensor1d Rk = std::make_shared<Tensor1d>("Res[k] vector", N);
	 Rk->row_vector(Res,k);
         if (B->gs_add(L,Rk)) { 
             L++; 
	     numf++; 
	 }
    }

    // If L is close to maxdim, collapse to one guess per root
    if (maxdim - L < M) {
	outfile->Printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
	outfile->Printf("Collapsing eigenvectors.\n");

	// form Bnew
        for (int i=0; i < M; i++) {
	     for (int k=0; k < N; k++) {
		  double sum = 0.0;
	          for (int j=0; j < L; j++) {
	               sum += alpha->get(j,i) * B->get(j,k);
	          } 
		  Bnew->set(i,k,sum);
	     }
        }

        // Copy new vectors into place
        Bnew->zero();
        for (int i=0; i < M; i++) {
	     for (int k=0; k < N; k++) {
                  B->set(i,k, Bnew->get(i,k));
             }
        }
        skip_check = 1;
        L = M;
    } // end if

    // check convergence on all roots
    if (!skip_check) {
        converged = 0;
	conv->zero();
	outfile->Printf("\n\tRoot      Eigenvalue       Delta  Converged?\n");
	outfile->Printf("\t---- -------------------- ------- ----------\n");

        for (int k=0; k < M; k++) {
	     diff = std::fabs(lambda->get(k) - lambda_old->get(k));
	     if (diff < cutoff) {
		 conv->set(k,1);
	         converged++;
	     }
	     lambda_old->set(k,lambda->get(k));
	     outfile->Printf("\t%3d  %20.14f %4.3e    %1s\n", k+1, lambda->get(k), diff, conv->get(k) == 1 ? "Y" : "N");
        }
    }// end if (!skip_check)

    iter++;
  }// end of loop

  // Generate final eigenvalues and eigenvectors
  if (converged == M) {
      for (int i=0; i < M; i++) {
           eps->set(i,lambda->get(i));
	   for (int I=0; I < N; I++) {
	        double sum = 0.0;
                for (int j=0; j < L; j++) {
		     sum += alpha->get(j,i) * B->get(j,I);
	        }  
	        v->add(I,i,sum);
           }
      }
      outfile->Printf("\n\tDavidson algorithm converged in %d iterations.\n", iter);
  }

  conv.reset();
  B.reset();
  Bnew.reset();
  sigma.reset();
  G.reset();
  Res.reset();
  alpha.reset();
  lambda.reset();
  lambda_old.reset();
  Hdiag.reset();
  Cci.reset();
  Sci.reset();

  return converged;
}

}  // namespace dfoccwave
}  // namespace psi


