/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace cphf {

void out_of_core(double ***F, double ***S, double ***B, double **A)
{
  int coord;
  int i, j, ij;
  double value;
  double *inbuf;
  char *label;
  double cutoff = 0.0;
  int k, l, n;
  int ii, jj, kk;
  int ji, ik, ki, il, li;
  int jk, kj, kl, lk, jl, lj; 
  int iq, jq, kq, lq;
  int tmp = 0;
  int endflag = 0;
  int bufsize = 0;
  int labelindex = 0;
  double intvalue = 0;
  struct iwlbuf buffer;
  
  /* Grab the MO-basis overlap and Fock derivative integrals from disk */
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  inbuf = init_array(ntri);
  
  for(coord=0; coord < natom*3; coord++) {

    sprintf(label, "MO-basis Fock Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    
    for(i=0; i < PSIO_KEYLEN; i++) {
      label[i] = '\0';
    }

    for(i=0, ij=0; i < nmo; i++) {
      for(j=0; j <= i; j++, ij++) {
	     F[coord][i][j] = F[coord][j][i] = inbuf[ij];
      }
    }
    if(print_lvl > 5) {
      fprintf(outfile, "F[%d] Deriv (MO):\n", coord);
      print_mat(F[coord], nmo, nmo, outfile);
    }
  }

  for(coord=0; coord < natom*3; coord++) {
    
    sprintf(label, "MO-basis Overlap Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    
    for(i=0; i < PSIO_KEYLEN; i++) {
      label[i] = '\0';
    }

    for(i=0, ij=0; i < nmo; i++) {
      for(j=0; j <= i; j++, ij++) {
        S[coord][i][j] = S[coord][j][i] = inbuf[ij];
      }
    }

    if(print_lvl > 5) {
      fprintf(outfile, "S[%d] Deriv (MO):\n", coord);
      print_mat(S[coord], nmo, nmo, outfile);
    }
  }

  /* Fock and overlap derivative components */
  for(coord=0; coord < natom*3; coord++) {
    for(i=0; i < nmo; i++) {
      for(j=0; j < nmo; j++) {
        B[coord][i][j] = F[coord][i][j] - S[coord][i][j] * evals[j];
      }
    }

    if(print_lvl > 5) {
      fprintf(outfile, "B[%d]_ij = F_ij-S_ij*e_j:\n", coord);
      print_mat(B[coord], nmo, nmo, outfile);
    }
  }

  free(inbuf);
  free(label);

  /* Two-Electron Integrals Out-Of-Core 
   *
   * Out-of-Core Cases
   * Case 1  : (pp|pp)
   * Case 2  : (pp|pq)
   * Case 3  : (pq|qq)
   * Case 4  : (pp|qq)
   * Case 5  : (pq|pq)
   * Case 6  : (pp|qr) 
   * Case 7  : (pq|qr)
   * Case 8  : (pq|pr)
   * Case 9  : (pq|rq)
   * Case 10 : (pq|rr)
   * Case 11 : (pq|rs)
   *
   * MO basis
   * B_a_pq = -S_a_rs*(2*(pq|rs)-(pr|qs))
   * A_pq,rs = -4*(pq|rs)+(pr|qs)+(ps|qr)
   * 
  */

  /* Prepare arrays in QTS ordering so  
   * the integrals can be screened according
   * to orbital space; this is needed for 
   * the B-matrix.
  */

  qtsorder = init_int_array(nmo);
  
  reorder_qt(clsdpi, openpi, frdoccpi, fruoccpi, qtsorder, orbspi, nirreps);

  iwl_buf_init(&buffer, PSIF_MO_TEI, cutoff, 1, 1);
  
  while (endflag == 0)
  {
    bufsize = buffer.inbuf;
    endflag = buffer.lastbuf;
    labelindex = 0;

    for (n = 0; n < bufsize; n++) {
      /* Pitzer ordering */ 
      i = abs(buffer.labels[labelindex++]);   
      j = buffer.labels[labelindex++];
      k = buffer.labels[labelindex++];
      l = buffer.labels[labelindex++];
      intvalue = buffer.values[n];

      /* QTS ordering */
      iq = qtsorder[i];
      jq = qtsorder[j];
      kq = qtsorder[k];
      lq = qtsorder[l];

      /* Case 1 : (pp|pp) */
      if ((i == j) && (j == k) && (k == l)) {
        /* (oo|oo) */	
	if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][i] -= S[coord][i][i] * intvalue;
          }
	}

	ii = INDEX(i,i);
	
	A[ii][ii] -= 2 * intvalue;
      }  
     
      /* Case 2 : (pp|pq) */
      else if ((i == j) && (j == k)) {
        /* (oo|ov) */
	if (iq < ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][l] -= S[coord][i][i] * intvalue;
            B[coord][l][i] -= S[coord][i][i] * intvalue;	    
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][i] -= 2 * S[coord][i][l] * intvalue;
            B[coord][i][l] -= S[coord][i][i] * intvalue;
            B[coord][l][i] -= S[coord][i][i] * intvalue;
          }
        }
	
	ii = INDEX(i,i);
	il = INDEX(i,l);
	
	A[ii][il] -= 2 * intvalue;
	A[il][ii] -= 2 * intvalue;
      }  
        
      /* Case 3 : (pq|qq) */
      else if ((j == k) && (k == l)) {
	/* (vo|oo) */      
	if (iq >= ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][j] -= S[coord][j][j] * intvalue;
	    B[coord][j][i] -= S[coord][j][j] * intvalue;
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) { 
            B[coord][i][j] -= S[coord][j][j] * intvalue;
	    B[coord][j][i] -= S[coord][j][j] * intvalue;
	    B[coord][j][j] -= 2 * S[coord][i][j] * intvalue;
          }
	}

	ij = INDEX(i,j);
	jj = INDEX(j,j);
        
	A[ij][jj] -= 2 * intvalue;
	A[jj][ij] -= 2 * intvalue;
      } 
        
      /* Case 4 : (pp|qq) */
      else if (i == j && k == l) {
	/* (vv|oo) */      
        if (iq >= ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][i] -= 2 * S[coord][k][k] * intvalue;
	  }
	}
	/* (oo|vv) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][k] -= 2 * S[coord][i][i] * intvalue;
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) { 
	    B[coord][i][i] -= 2 * S[coord][k][k] * intvalue;
            B[coord][k][k] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][i][k] += S[coord][i][k] * intvalue;
	    B[coord][k][i] += S[coord][i][k] * intvalue;
	  }
	}

	ii = INDEX(i,i);
	kk = INDEX(k,k);
	ik = INDEX(i,k);
        
	A[ik][ik] += 1 * intvalue;
	A[ii][kk] -= 4 * intvalue;
	A[kk][ii] -= 4 * intvalue;
      } 
	    
      /* Case 5 : (pq|pq) */
      else if ((i == k) && (j == l)) {
	/* (vo|vo) */      
        if (iq >= ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (ov|ov) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][j] += S[coord][i][i] * intvalue;
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][i][j] * intvalue;
            B[coord][i][i] += S[coord][j][j] * intvalue;
            B[coord][j][j] += S[coord][i][i] * intvalue;
          }
	}

	ij = INDEX(i,j);
        ii = INDEX(i,i);
	jj = INDEX(j,j);
        
	A[ij][ij] -= 3 * intvalue;
	A[ii][jj] += 2 * intvalue;
	A[jj][ii] += 2 * intvalue;
      }
      
      /* Case 6 : (pp|qr) */
      else if (k != l && i == j) {
	/* (vv|oo) */      
        if (iq >= ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {	
            B[coord][i][i] -= 4 * S[coord][k][l] * intvalue;
	  }
	}
	/* (oo|vv) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
            B[coord][k][l] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][l][k] -= 2 * S[coord][i][i] * intvalue;
	  }
	}
        /* (oo|vo) */	
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][l] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][l][k] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][i][k] += S[coord][i][l] * intvalue;
	    B[coord][k][i] += S[coord][i][l] * intvalue;
	  }
	}
	/* (oo|ov) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][l] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][l][k] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][i][l] += S[coord][i][k] * intvalue;
	    B[coord][l][i] += S[coord][i][k] * intvalue;
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][i] -= 4 * S[coord][k][l] * intvalue;
	    B[coord][k][l] -= 2 * S[coord][i][i] * intvalue;
	    B[coord][l][k] -= 2 * S[coord][i][i] * intvalue;
            B[coord][i][k] += S[coord][i][l] * intvalue;
	    B[coord][k][i] += S[coord][i][l] * intvalue;
	    B[coord][i][l] += S[coord][i][k] * intvalue;
	    B[coord][l][i] += S[coord][i][k] * intvalue;
	  }
	}
	
	ii = INDEX(i,i);
	kl = INDEX(k,l);
	ik = INDEX(i,k);
	il = INDEX(i,l);
        
 	A[ii][kl] -= 4 * intvalue;
	A[kl][ii] -= 4 * intvalue;
	A[ik][il] += 1 * intvalue;
        A[il][ik] += 1 * intvalue;
      }
     
      /* Case 7 : (pq|qr) */
      else if (i != l && j == k) {
	/* (vo|ov) */      
        if (iq >= ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][l] += S[coord][j][j] * intvalue;
            B[coord][l][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (ov|vo) */
        else if (iq < ndocc && jq >= ndocc && kq >= ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][j] += 2 * S[coord][i][l] * intvalue;
	  }
	}
	/* (oo|ov) */
        else if (iq < ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][l] -= 3 * S[coord][i][j] * intvalue;
            B[coord][l][j] -= 3 * S[coord][i][j] * intvalue;	    
	    B[coord][i][l] += S[coord][j][j] * intvalue;
	    B[coord][l][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (vo|oo) */
	else if (iq >= ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 3 * S[coord][j][l] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][j][l] * intvalue;
	    B[coord][i][l] += S[coord][j][j] * intvalue;
	    B[coord][l][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (oo|oo) */
        else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 3 * S[coord][j][l] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][j][l] * intvalue;
	    B[coord][j][l] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][l][j] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][j][j] += 2 * S[coord][i][l] * intvalue;
	    B[coord][i][l] += S[coord][j][j] * intvalue;
	    B[coord][l][i] += S[coord][j][j] * intvalue;
	  }
	}
	
        ij = INDEX(i,j);
        jl = INDEX(j,l);
	il = INDEX(i,l);
	jj = INDEX(j,j);

	A[ij][jl] -= 3 * intvalue;
	A[jl][ij] -= 3 * intvalue;
	A[il][jj] += 2 * intvalue;
	A[jj][il] += 2 * intvalue;
      }
      
      /* Case 8 : (pq|pr) */
      else if (j != l && i == k) {
	/* (vo|vo) */
	if (iq >= ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][i] += 2 * S[coord][j][l] * intvalue;
	  }
	}
	/* (ov|ov) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][l] += S[coord][i][i] * intvalue;
	    B[coord][l][j] += S[coord][i][i] * intvalue;
	  }
	}
	/* (oo|ov) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][l] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][l][i] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][j][l] += S[coord][i][i] * intvalue;
	    B[coord][l][j] += S[coord][i][i] * intvalue;
	  }
	}
	/* (ov|oo) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 3 * S[coord][i][l] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][i][l] * intvalue;
	    B[coord][j][l] += S[coord][i][i] * intvalue;
	    B[coord][l][j] += S[coord][i][i] * intvalue;
	  }
	}
	/* (oo|oo) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 3 * S[coord][i][l] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][i][l] * intvalue;
	    B[coord][i][l] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][l][i] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][i][i] += 2 * S[coord][j][l] * intvalue;
	    B[coord][j][l] += S[coord][i][i] * intvalue;
	    B[coord][l][j] += S[coord][i][i] * intvalue;
	  }
	}
	      
	ij = INDEX(i,j);
	il = INDEX(i,l);
	ii = INDEX(i,i);
	jl = INDEX(j,l);
        
	A[ij][il] -= 3 * intvalue;
	A[il][ij] -= 3 * intvalue;
	A[ii][jl] += 2 * intvalue;
	A[jl][ii] += 2 * intvalue;
      }

      /* Case 9 : (pq|rq) */
      else if (i !=k && j == l) {
        /* (vo|vo) */
	if (iq >= ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][k] += S[coord][j][j] * intvalue;
	    B[coord][k][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (ov|ov) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq >= ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][j] += 2 * S[coord][i][k] * intvalue;
	  }
	}
	/* (oo|vo) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][k] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][k][j] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][i][k] += S[coord][j][j] * intvalue;
	    B[coord][k][i] += S[coord][j][j] * intvalue;
	  }
	}
	/* (vo|oo) */
	else if (iq >= ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][j] -= 3 * S[coord][j][k] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][j][k] * intvalue;
	    B[coord][i][k] += S[coord][j][j] * intvalue;
	    B[coord][k][i] += S[coord][j][j] * intvalue;
	  }
	}
        /* (oo|oo) */	
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][i][j] -= 3 * S[coord][j][k] * intvalue;
	    B[coord][j][i] -= 3 * S[coord][j][k] * intvalue;
	    B[coord][j][k] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][k][j] -= 3 * S[coord][i][j] * intvalue;
	    B[coord][j][j] += 2 * S[coord][i][k] * intvalue;
	    B[coord][i][k] += S[coord][j][j] * intvalue;
	    B[coord][k][i] += S[coord][j][j] * intvalue;
	  }
	}
	
	ij = INDEX(i,j);
	kj = INDEX(k,j);
	jj = INDEX(j,j);
	ik = INDEX(i,k);
        
	A[ij][kj] -= 3 * intvalue;
	A[kj][ij] -= 3 * intvalue;
	A[jj][ik] += 2 * intvalue;
	A[ik][jj] += 2 * intvalue;
      }

      /* Case 10 : (pq|rr) */
      else if (i != j && k == l) {
        /* (vv|oo) */
	if (iq >= ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 2 * S[coord][k][k] * intvalue;
            B[coord][j][i] -= 2 * S[coord][k][k] * intvalue;
	  }
	}
	/* (oo|vv) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][k] -= 4 * S[coord][i][j] * intvalue;
	  }
	}
	/* (ov|oo) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 2 * S[coord][k][k] * intvalue;
            B[coord][j][i] -= 2 * S[coord][k][k] * intvalue;
	    B[coord][j][k] += S[coord][i][k] * intvalue;
            B[coord][k][j] += S[coord][i][k] * intvalue;
	  }
	}
	/* (vo|oo) */
	else if (iq >= ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 2 * S[coord][k][k] * intvalue;
            B[coord][j][i] -= 2 * S[coord][k][k] * intvalue;
	    B[coord][i][k] += S[coord][j][k] * intvalue;
            B[coord][k][i] += S[coord][j][k] * intvalue;
	  }
	}
        /* (oo|oo) */	
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
            B[coord][k][k] -= 4 * S[coord][i][j] * intvalue;
	    B[coord][i][j] -= 2 * S[coord][k][k] * intvalue;
	    B[coord][j][i] -= 2 * S[coord][k][k] * intvalue;
	    B[coord][i][k] += S[coord][j][k] * intvalue;
	    B[coord][k][i] += S[coord][j][k] * intvalue;
            B[coord][j][k] += S[coord][i][k] * intvalue;
            B[coord][k][j] += S[coord][i][k] * intvalue;
	  }
	}

	ij = INDEX(i,j);
	kk = INDEX(k,k);
	ik = INDEX(i,k);
	jk = INDEX(j,k);
        
	A[ij][kk] -= 4 * intvalue;
	A[kk][ij] -= 4 * intvalue;
	A[ik][jk] += 1 * intvalue;
	A[jk][ik] += 1 * intvalue;
      }

      /* Case 11 : (pq|rs) */
      else if (i != j && j != k && k != l) {
        /* (vv|oo) */
	if (iq >= ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 4 * S[coord][k][l] * intvalue;
            B[coord][j][i] -= 4 * S[coord][k][l] * intvalue;
	  }
	}
	/* (vo|ov) */
	else if (iq >= ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][l] += S[coord][j][k] * intvalue;
            B[coord][l][i] += S[coord][k][j] * intvalue;
	  }
	}
	/* (vo|vo) */
	else if (iq >= ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][k] += S[coord][j][l] * intvalue;
            B[coord][k][i] += S[coord][l][j] * intvalue;
	  }
	}
	/* (ov|ov) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][l] += S[coord][i][k] * intvalue;
            B[coord][l][j] += S[coord][i][k] * intvalue;
	  }
	}
	/* (ov|vo) */
	else if (iq < ndocc && jq >= ndocc && kq >= ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][j][k] += S[coord][i][l] * intvalue;
            B[coord][k][j] += S[coord][i][l] * intvalue;
	  }
	}
	/* (oo|vv) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][l] -= 4 * S[coord][i][j] * intvalue;
            B[coord][l][k] -= 4 * S[coord][i][j] * intvalue;
	  }
	}
	/* (vo|oo) */
	else if (iq >= ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 4 * S[coord][k][l] * intvalue;
            B[coord][j][i] -= 4 * S[coord][k][l] * intvalue;
	    B[coord][i][l] += S[coord][j][k] * intvalue;
            B[coord][l][i] += S[coord][k][j] * intvalue;
	    B[coord][i][k] += S[coord][j][l] * intvalue;
	    B[coord][k][i] += S[coord][l][j] * intvalue;
	  }
	}
	/* (ov|oo) */
	else if (iq < ndocc && jq >= ndocc && kq < ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 4 * S[coord][k][l] * intvalue;
            B[coord][j][i] -= 4 * S[coord][k][l] * intvalue;
	    B[coord][j][l] += S[coord][i][k] * intvalue;
            B[coord][l][j] += S[coord][i][k] * intvalue;
	    B[coord][j][k] += S[coord][i][l] * intvalue;
            B[coord][k][j] += S[coord][i][l] * intvalue;
	  }
	}
	/* (oo|vo) */
	else if (iq < ndocc && jq < ndocc && kq >= ndocc && lq < ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][l] -= 4 * S[coord][i][j] * intvalue;
            B[coord][l][k] -= 4 * S[coord][i][j] * intvalue;
            B[coord][j][k] += S[coord][i][l] * intvalue;
            B[coord][k][j] += S[coord][i][l] * intvalue;
	    B[coord][i][k] += S[coord][j][l] * intvalue;
            B[coord][k][i] += S[coord][l][j] * intvalue;
	  }
	}
	/* (oo|ov) */
	else if (iq < ndocc && jq < ndocc && kq < ndocc && lq >= ndocc) {
          for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][k][l] -= 4 * S[coord][i][j] * intvalue;
            B[coord][l][k] -= 4 * S[coord][i][j] * intvalue;
	    B[coord][j][l] += S[coord][i][k] * intvalue;
            B[coord][l][j] += S[coord][i][k] * intvalue;
	    B[coord][i][l] += S[coord][j][k] * intvalue;
            B[coord][l][i] += S[coord][k][j] * intvalue;
	  }
	}
	/* (oo|oo) */      
	else if(iq < ndocc && jq < ndocc && kq < ndocc && lq < ndocc) {
	  for (coord = 0; coord < natom * 3; coord++) {
	    B[coord][i][j] -= 4 * S[coord][k][l] * intvalue;
	    B[coord][j][i] -= 4 * S[coord][k][l] * intvalue;
	    B[coord][k][l] -= 4 * S[coord][i][j] * intvalue;
	    B[coord][l][k] -= 4 * S[coord][i][j] * intvalue;
	    B[coord][i][k] += S[coord][j][l] * intvalue;
	    B[coord][k][i] += S[coord][j][l] * intvalue;
	    B[coord][i][l] += S[coord][j][k] * intvalue;
	    B[coord][l][i] += S[coord][j][k] * intvalue;
	    B[coord][j][k] += S[coord][i][l] * intvalue;
	    B[coord][k][j] += S[coord][i][l] * intvalue;
	    B[coord][j][l] += S[coord][i][k] * intvalue;
	    B[coord][l][j] += S[coord][i][k] * intvalue;
	  }
	}
	
	ij = INDEX(i,j);
	kl = INDEX(k,l);
	ik = INDEX(i,k);
	jl = INDEX(j,l);
	il = INDEX(i,l);
	jk = INDEX(j,k);
        
	A[ij][kl] -= 4 * intvalue;
	A[kl][ij] -= 4 * intvalue;
	A[ik][jl] += 1 * intvalue;
	A[jl][ik] += 1 * intvalue;
	A[il][jk] += 1 * intvalue;
	A[jk][il] += 1 * intvalue;
      }

    } /* end of for loop: buffer */
    
    if (endflag == 0) {
      iwl_buf_fetch(&buffer);
    }
  } /* End of while loop: out-of-core */
  
  iwl_buf_close(&buffer,1);

  if(print_lvl > 5) { 
    for(coord=0; coord < 3*natom; coord++) {
      fprintf(outfile, "B(p,q)[%d] (MO):\n", coord);
      print_mat(B[coord], nmo, nmo, outfile);
    }
  }
}

void sort_B(double ***B, double **Baijk)
{
  int coord, a, i, AI;
  int asym, isym;
  int afirst, alast;
  int ifirst, ilast;
   
  for(coord=0; coord < natom*3; coord++) {
     
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

        for(isym = 0; isym < nirreps; isym++) {
         
          ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    Baijk[coord][AI] = B[coord][a][i];
	  }
	}
      }
    }
  }
  
  return;
}

void sort_A(double **A, double **Aaibj)
{
  int asym, isym, bsym, jsym, psym, qsym;
  int a, i, b, j, p, q;
  int ai, bj, pi, qj;
  int AI, BJ, PI, QJ;
  int afirst, alast, ifirst, ilast;
  int bfirst, blast, jfirst, jlast;
  int pfirst, plast, qfirst, qlast;

  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	   
        ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {
          
	  ai = INDEX(a,i);
  	  
	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
               
	      for(jsym=0; jsym < nirreps; jsym++) {

	        jfirst = ofirst[jsym];
	        jlast = olast[jsym];

	        for(j=jfirst; j <= jlast; j++,BJ++) {
                  
		  bj = INDEX(b,j);
	  	  
	          Aaibj[AI][BJ] = A[ai][bj];
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}

}} // namespace psi::cphf
