/*! \file
    \ingroup LMP2
    \brief compute the projection matrix 
*/
#include <iostream>
#include <fstream>              // file I/O support
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;

void LMP2::projection() {

  int i, j, ij, cnt, v, K, L, I;
  int r, s, k, l;
  int *ij_owner, *ij_local, count;
  int num_zero;
  double **Xt;
  double **St;
  double **Ft;
  double **Fbar;
  double **temp;
  double **evecs, *evals;

  /* Set up ij_owner and ij_local */
  ij_owner = get_ij_owner();
  ij_local = get_ij_local();

/*  ij_owner = init_int_array(ij_pairs);
  ij_local = init_int_array(ij_pairs);


  v=0;
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++) {
      ij_owner[ij] = v%nprocs;
      v++;
    }
  }

  v=0;
  count=0;
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++) {
      ij_local[ij] = count;
      if(v%nprocs == nprocs-1) count++;
      v++;
    }
  }
*/

  /* Compute the complete virtual space projector */
  Rt_full = block_matrix(nso,nso);
  for(i=0; i < nso; i++) Rt_full[i][i] = 1.0;

  C_DGEMM('n','n',nso,nso,nso,-1.0,&(D[0][0]),nso,&(aoovlp[0][0]),nso,
          1.0,&(Rt_full[0][0]),nso);
  free_block(D);

  // Transform the overlap matrix from the AO to PAO basis
  paoovlp = block_matrix(nso,nso);
  temp = block_matrix(nso,nso);
  C_DGEMM('t', 'n', nso, nso, nso, 1, &(Rt_full[0][0]), nso, &(aoovlp[0][0]), nso, 0, &(temp[0][0]), nso);
  C_DGEMM('n', 'n', nso, nso, nso, 1, &(temp[0][0]), nso, &(Rt_full[0][0]), nso, 0, &(paoovlp[0][0]), nso);

  // Transform the Fock matrix from the AO to PAO basis
  paoF = block_matrix(nso,nso);
  C_DGEMM('t', 'n', nso, nso, nso, 1, &(Rt_full[0][0]), nso, &(aoF[0][0]), nso, 0, &(temp[0][0]), nso);
  C_DGEMM('n', 'n', nso, nso, nso, 1, &(temp[0][0]), nso, &(Rt_full[0][0]), nso, 0, &(paoF[0][0]), nso);
  free_block(temp);

  free_block(aoovlp);
  free_block(aoF);

  /* Build the virtual metric and W transforms for each pair domain */
  if(ij_pairs%nprocs == 0) {
    W = (double ***) malloc((ij_pairs/nprocs) * sizeof(double **));
  }
  else {
    if(myid < ij_pairs%nprocs) {
      W = (double ***) malloc(((ij_pairs/nprocs) + 1) * sizeof(double **));
    }
    else {
      W = (double ***) malloc(ij_pairs/nprocs * sizeof(double **));
    }
  }

//  W = (double ***) malloc(ij_pairs * sizeof(double **));
  loevals = (double **) malloc(ij_pairs * sizeof(double*));
  pairdom_nrlen = init_int_array(ij_pairs); /* dimension of non-redundant basis */
  num_zero = 0;

  v=0;
  for(i=0, ij=0; i < nocc; i++) {
    for(j=0; j <= i; j++, ij++, v++) {

      if(v%nprocs != myid)
        continue;

    /* Build the virtual space overlap matrix for this pair */
    St = block_matrix(pairdom_len[ij],pairdom_len[ij]);
    for(r=0, K=0; r < natom; r++) {
      if(pairdomain[ij][r]) {
        for(k=aostart[r]; k <= aostop[r]; k++,K++) {
          for(s=0, L=0; s < natom; s++) {
            if(pairdomain[ij][s]) {
              for(l=aostart[s]; l <= aostop[s]; l++,L++) {

                St[K][L] = paoovlp[k][l];
              }
            }
          }
        }
      }
    }

//    fprintf(outfile, "pairdom_len[ij] = %d\n", pairdom_len[ij]);
//    print_mat(St, pairdom_len[ij], pairdom_len[ij], outfile);

    /* Diagonalize metric */
    evals = init_array(pairdom_len[ij]);
    evecs = block_matrix(pairdom_len[ij],pairdom_len[ij]);
    sq_rsp(pairdom_len[ij], pairdom_len[ij], St, evals, 1, evecs, 1e-12);
    free_block(St);

//    fprintf(outfile, "\nevecs and evals for St diag for ij = %d:\n", ij);
//    fprintf(outfile,   "======\n");
//    eivout(evecs, evals, pairdom_len[ij], pairdom_len[ij], outfile);


    /* Count the number of zero eigenvalues */
    for(k=0, cnt=0; k < pairdom_len[ij]; k++) if(evals[k] <= 1e-6) cnt++;

    pairdom_nrlen[ij] = pairdom_len[ij]-cnt;

    /*
      fprintf(outfile, "\n\tS-tilde eigenvalues for ij = %d:\n", ij);
      for(i=0; i < pairdom_len[ij]; i++) fprintf(outfile, "\t%d %20.12f\n", i, evals[i]);

      fprintf(outfile, "\n\tS-tilde eigenvectors for ij = %d:\n", ij);
      print_mat(evecs,pairdom_len[ij],pairdom_len[ij],outfile);
    */

    /* Build the projected, non-redundant transform (X-tilde) */
    Xt = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
    for(k=0, I=0; k < pairdom_len[ij]; k++) {
      if(evals[k] > 1e-6) {
	for(l=0; l < pairdom_len[ij]; l++)
	  Xt[l][I] = evecs[l][k]/sqrt(evals[k]);
	I++;
      }
      else num_zero++;
    }
    
//    fprintf(outfile, "\nXt[%d]:\n", ij);
//    fprintf(outfile,   "======\n");
//    print_mat(Xt, pairdom_len[ij], pairdom_nrlen[ij], outfile);
	    
    free_block(evecs);
    free(evals);

//    fprintf(outfile, "\n\taoF matrix for ij = %d:\n", ij);
//    print_mat(aoF,pairdom_len[ij],pairdom_len[ij],outfile);

    /* Build the virtual space Fock matrix for this pair */
    Ft = block_matrix(pairdom_len[ij],pairdom_len[ij]);
    for(r=0, K=0; r < natom; r++) {
      if(pairdomain[ij][r]) {
        for(k=aostart[r]; k <= aostop[r]; k++,K++) {
          for(s=0, L=0; s < natom; s++) {
            if(pairdomain[ij][s]) {
              for(l=aostart[s]; l <= aostop[s]; l++,L++) {

                Ft[K][L] = paoF[k][l];
              }
            }
          }
        }
      }
    }


//    free_block(temp);

//    fprintf(outfile, "\nFt[%d]:\n", ij);
//    fprintf(outfile,   "======\n");
//    print_mat(Ft,pairdom_len[ij],pairdom_len[ij],outfile);

    // Project the Fock matrix into the non-redundant virtual space 
    Fbar = block_matrix(pairdom_nrlen[ij],pairdom_nrlen[ij]);
    temp = block_matrix(pairdom_nrlen[ij],pairdom_len[ij]);
    C_DGEMM('t','n',pairdom_nrlen[ij],pairdom_len[ij],pairdom_len[ij],1.0,
	    &(Xt[0][0]),pairdom_nrlen[ij],&(Ft[0][0]),pairdom_len[ij],0.0,&(temp[0][0]),
            pairdom_len[ij]);
    C_DGEMM('n','n',pairdom_nrlen[ij],pairdom_nrlen[ij],pairdom_len[ij],1.0,
	    &(temp[0][0]),pairdom_len[ij],&(Xt[0][0]),pairdom_nrlen[ij],0.0,
            &(Fbar[0][0]),pairdom_nrlen[ij]);
    free_block(temp);
    free_block(Ft);

//    fprintf(outfile, "\nFbar[%d] before diag:\n", ij);
//    fprintf(outfile,   "======\n");
//    print_mat(Fbar,pairdom_nrlen[ij],pairdom_nrlen[ij],outfile);
 
    // Diagonalize Fbar 
    loevals[ij] = init_array(pairdom_nrlen[ij]);
    evecs = block_matrix(pairdom_nrlen[ij],pairdom_nrlen[ij]);
    zero_arr(loevals[ij], pairdom_nrlen[ij]);
    sq_rsp(pairdom_nrlen[ij],pairdom_nrlen[ij],Fbar,loevals[ij],1,evecs,1e-12);

//    fprintf(outfile, "\nFbar[%d]:\n", ij);
//    fprintf(outfile,   "======\n");
//    print_mat(evecs,pairdom_nrlen[ij],pairdom_nrlen[ij],outfile);
    

    // Finally, build the W matrix 
    if(myid == ij_owner[ij]) {
      if(ij_pairs%nprocs == 0) {
        W[ij_local[ij]] = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
      }
      else {
        if(myid < ij_pairs%nprocs)
            W[ij_local[ij]] = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
        else 
            W[ij_local[ij]] = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
      }
    }

//    W[ij] = block_matrix(pairdom_len[ij],pairdom_nrlen[ij]);
   
    C_DGEMM('n','n',pairdom_len[ij],pairdom_nrlen[ij],pairdom_nrlen[ij],1.0,
	    &(Xt[0][0]),pairdom_nrlen[ij],&(evecs[0][0]),pairdom_nrlen[ij],
	    0.0,&(W[ij_local[ij]][0][0]),pairdom_nrlen[ij]);

    free_block(evecs);

    if(print > 1) {
      fprintf(outfile, "\nW-Trans[%d]:\n", ij);
      fprintf(outfile,   "======\n");
      print_mat(W[ij_local[ij]],pairdom_len[ij],pairdom_nrlen[ij],outfile);
    }

//    fprintf(outfile, "\nEigenvalues[%d]:\n", ij);
//    fprintf(outfile,   "======\n");
//    for(i=0; i < pairdom_nrlen[ij]; i++)
//    fprintf(outfile, "%d %20.12f\n", i, loevals[ij][i]);


//    fprintf(outfile, "pairdom_nrlen[%d] = %d\n", ij, pairdom_nrlen[ij]);
  
    free_block(Xt);
    free_block(Fbar);

    } // j loop
  } // i loop

  fflush(outfile);


}

}} // namespace psi::lmp2
