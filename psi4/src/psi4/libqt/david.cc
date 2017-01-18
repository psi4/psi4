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

/*!
  \file
  \brief In-core Davidson-Liu diagonalization of symm matrices
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

namespace psi {

#define BIGNUM 1E100
#define MAXIT 1000

/*!
** david(): Computes the lowest few eigenvalues and eigenvectors of a
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
** \ingroup QT
*/

int david(double **A, int N, int M, double *eps, double **v,
	  double cutoff, int print)
{
  int i, j, k, L, I;
  double minimum;
  int min_pos, numf, iter, *conv, converged, maxdim, skip_check;
  int *small2big, init_dim;
  int smart_guess = 1;
  double *Adiag, **b, **bnew, **sigma, **G;
  double *lambda, **alpha, **f, *lambda_old;
  double norm, denom, diff;

  maxdim = 8 * M;

  b = block_matrix(maxdim, N);  /* current set of guess vectors,
				   stored by row */
  bnew = block_matrix(M, N); /* guess vectors formed from old vectors,
				stored by row*/
  sigma = block_matrix(N, maxdim); /* sigma vectors, stored by column */
  G = block_matrix(maxdim, maxdim); /* Davidson mini-Hamitonian */
  f = block_matrix(maxdim, N); /* residual eigenvectors, stored by row */
  alpha = block_matrix(maxdim, maxdim); /* eigenvectors of G */
  lambda = init_array(maxdim); /* eigenvalues of G */
  lambda_old = init_array(maxdim); /* approximate roots from previous
				      iteration */

  if(smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */

    if(N > 7*M) init_dim = 7*M;
    else init_dim = M;
    Adiag = init_array(N);
    small2big = init_int_array(7*M);
    for(i=0; i < N; i++) { Adiag[i] = A[i][i]; }
    for(i=0; i < init_dim; i++) {
      minimum = Adiag[0];
      min_pos = 0;
      for(j=1; j < N; j++)
	if(Adiag[j] < minimum) {
	  minimum = Adiag[j];
	  min_pos = j;
	  small2big[i] = j;
	}

      Adiag[min_pos] = BIGNUM;
      lambda_old[i] = minimum;
    }
    for(i=0; i < init_dim; i++) {
      for(j=0; j < init_dim; j++)
	G[i][j] = A[small2big[i]][small2big[j]];
    }

    sq_rsp(init_dim, init_dim, G, lambda, 1, alpha, 1e-12);

    for(i=0; i < init_dim; i++) {
      for(j=0; j < init_dim; j++)
	b[i][small2big[j]] = alpha[j][i];
    }

    free(Adiag);
    free(small2big);
  }
  else { /* Use unit vectors as initial guesses */
    Adiag = init_array(N);
    for(i=0; i < N; i++) { Adiag[i] = A[i][i]; }
    for(i=0; i < M; i++) {
      minimum = Adiag[0];
      min_pos = 0;
      for(j=1; j < N; j++)
	if(Adiag[j] < minimum) { minimum = Adiag[j]; min_pos = j; }

      b[i][min_pos] = 1.0;
      Adiag[min_pos] = BIGNUM;
      lambda_old[i] = minimum;
    }
    free(Adiag);
  }

  L = init_dim;
  iter =0;
  converged = 0;
  conv = init_int_array(M); /* boolean array for convergence of each
			       root */
  while(converged < M && iter < MAXIT) {

    skip_check = 0;
    if(print) printf("\niter = %d\n", iter);

    /* form mini-matrix */
    C_DGEMM('n','t', N, L, N, 1.0, &(A[0][0]), N, &(b[0][0]), N,
	    0.0, &(sigma[0][0]), maxdim);
    C_DGEMM('n','n', L, L, N, 1.0, &(b[0][0]), N,
	    &(sigma[0][0]), maxdim, 0.0, &(G[0][0]), maxdim);

    /* diagonalize mini-matrix */
    sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);

    /* form preconditioned residue vectors */
    for(k=0; k < M; k++)
      for(I=0; I < N; I++) {
	f[k][I] = 0.0;
	for(i=0; i < L; i++) {
	  f[k][I] += alpha[i][k] * (sigma[I][i] - lambda[k] * b[i][I]);
	}
	denom = lambda[k] - A[I][I];
	if(fabs(denom) > 1e-6) f[k][I] /= denom;
	else f[k][I] = 0.0;
      }

    /* normalize each residual */
    for(k=0; k < M; k++) {
      norm = 0.0;
      for(I=0; I < N; I++) {
	norm += f[k][I] * f[k][I];
      }
      norm = sqrt(norm);
      for(I=0; I < N; I++) {
	if(norm > 1e-6) f[k][I] /= norm;
	else f[k][I] = 0.0;
      }
    }

    /* schmidt orthogonalize the f[k] against the set of b[i] and add
       new vectors */
    for(k=0,numf=0; k < M; k++)
      if(schmidt_add(b, L, N, f[k])) { L++; numf++; }

    /* If L is close to maxdim, collapse to one guess per root */
    if(maxdim - L < M) {
      if(print) {
	printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
	printf("Collapsing eigenvectors.\n");
      }
      for(i=0; i < M; i++) {
	memset((void *) bnew[i], 0, N*sizeof(double));
	for(j=0; j < L; j++) {
	  for(k=0; k < N; k++) {
	    bnew[i][k] += alpha[j][i] * b[j][k];
	  }
	}
      }

      /* copy new vectors into place */
      for(i=0; i < M; i++)
	for(k=0; k < N; k++)
	  b[i][k] = bnew[i][k];

      skip_check = 1;

      L = M;
    }

    /* check convergence on all roots */
    if(!skip_check) {
      converged = 0;
      zero_int_array(conv, M);
      if(print) {
	printf("Root      Eigenvalue       Delta  Converged?\n");
	printf("---- -------------------- ------- ----------\n");
      }
      for(k=0; k < M; k++) {
	diff = fabs(lambda[k] - lambda_old[k]);
	if(diff < cutoff) {
	  conv[k] = 1;
	  converged++;
	}
	lambda_old[k] = lambda[k];
	if(print) {
	  printf("%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
		 conv[k] == 1 ? "Y" : "N");
	}
      }
    }

    iter++;
  }

  /* generate final eigenvalues and eigenvectors */
  if(converged == M) {
    for(i=0; i < M; i++) {
      eps[i] = lambda[i];
      for(j=0; j < L; j++) {
	for(I=0; I < N; I++) {
	  v[I][i] += alpha[j][i] * b[j][I];
	}
      }
    }
    if(print) printf("Davidson algorithm converged in %d iterations.\n", iter);
  }

  free(conv);
  free_block(b);
  free_block(bnew);
  free_block(sigma);
  free_block(G);
  free_block(f);
  free_block(alpha);
  free(lambda);
  free(lambda_old);

  return converged;
}

}
