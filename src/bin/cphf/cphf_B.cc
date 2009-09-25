/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

/* Convert atomic unit of magnetic flux density to tesla */
#define _aumfd2tesla 2.35051742E5

namespace psi { namespace cphf {

void zval_to_symbol(double,char*);

/* cphf_B(): Solve the first-order CPHF equations for a magnetic 
** field perturbation.
**
** The CPHF equations are:
**
** A U^a = B0^a
**
** where A is the MO hessian for a complex perburbation, U^a is the
** orbital response (CPHF) coefficient, and B0^a is the
** perturbation-dependent inhomogenous factor.  For magnetic field
** perturbations, B0^a is given by (MO basis):
**
** (B0^a)_ai = (mu^a)_ai
**
** where (mu^a)_pq is a magnetic dipole moment integral.
**
*/

void cphf_B(double ***UX, double **lx)
{
  int stat, coord;
  double *scratch;
  double ***B;
  double **B0;
  double **TMP;
  double **X;
  char *label;
  int a, asym, afirst, alast;
  int i, isym, ifirst, ilast;
  int b, bsym, bfirst, blast;
  int j, jsym, jfirst, jlast;
  int ab, ij, ib, aj, abij, ajib;
  int AI, BJ;
  int *ipiv, error;
  double **A, ***UB;
  double **I, **I_q;
  double ***S;
  double cutoff = 0.0;
  int endflag = 0;
  int bufsize = 0;
  int labelindex = 0;
  double intvalue = 0;
  struct iwlbuf buffer;
  int p, q, r, s;
  int pq, rs;
  int p_, q_, r_, s_;
  int n;
  int ai, bj, bi;
  double **A0;
  double tmp;
  double *eval;
  double **evec;
  double *rotstr;

  /*
  Compute MO Hessian for a complex perturbation (IC) 
    Hessian eigenvalues will be negative:
      A[AI][BJ] = (a==b) * (i==j) * (evals[i] - evals[a]);
      A[AI][BJ] += ints[abij] - ints[ajib];
    Hessian eigenvalues will be positive:
      A[AI][BJ] = (a==b) * (i==j) * (evals[a] - evals[i]);
      A[AI][BJ] += ints[ajib] - ints[abij];
  */
 
  ints = init_array(ntei);
  iwl_rdtwo(PSIF_MO_TEI, ints, ioff, nmo, 0, 0, 0, outfile);

  A = block_matrix(num_ai, num_ai);
  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {

	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
	      ab = INDEX(a,b);
	      ib = INDEX(i,b);

	      for(jsym=0; jsym < nirreps; jsym++) {

		jfirst = ofirst[jsym];
		jlast = olast[jsym];

		for(j=jfirst; j <= jlast; j++,BJ++) {
		  ij = INDEX(i,j);
		  aj = INDEX(a,j);

		  abij = INDEX(ab,ij);
		  ajib = INDEX(aj,ib);

		  A[AI][BJ] = (a==b) * (i==j) * (evals[a] - evals[i]);
		  A[AI][BJ] += ints[ajib] - ints[abij];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  free(ints);

  //fprintf(outfile, "\nMO Hessian (IC):\n");
  //print_mat(A, num_ai, num_ai, outfile);

  /* Compute MO Hessian for a complex perturbation (OC) */

  /*
  A = block_matrix(num_ai, num_ai);
  iwl_buf_init(&buffer, PSIF_MO_TEI, cutoff, 1, 1);
  
  while (endflag == 0)
  {
    bufsize = buffer.inbuf;
    endflag = buffer.lastbuf;
    labelindex = 0;

    for (n = 0; n < bufsize; n++) {
      p = abs(buffer.labels[labelindex++]);   
      q = buffer.labels[labelindex++];
      r = buffer.labels[labelindex++];
      s = buffer.labels[labelindex++];
      intvalue = buffer.values[n];

      pq = INDEX(p,q);
      rs = INDEX(r,s);

      p_ = qtsorder[p];
      q_ = qtsorder[q];
      r_ = qtsorder[r];
      s_ = qtsorder[s];

      if (p!=q && r!=s && pq!=rs) {
	if (p_ < ndocc && q_ < ndocc && r_ >= ndocc && s_ >= ndocc) {
          ai = (r_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + q_;
          aj = (r_-ndocc) * ndocc + q_;
          bi = (s_-ndocc) * ndocc + p_;
          A[ai][bj] -= intvalue;
          A[aj][bi] -= intvalue;
          A[bi][aj] -= intvalue;
          A[bj][ai] -= intvalue;
	}
	else if (p_ < ndocc && q_ >= ndocc && r_ < ndocc && s_ >= ndocc) {
          ai = (q_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + r_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ < ndocc && s_ >= ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (s_-ndocc) * ndocc + r_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ >= ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (r_-ndocc) * ndocc + s_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ >= ndocc && r_ < ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + r_;
          bj = (q_-ndocc) * ndocc + s_;
          aj = (p_-ndocc) * ndocc + s_;
          bi = (q_-ndocc) * ndocc + r_;
          A[ai][bj] -= intvalue;
          A[aj][bi] -= intvalue;
          A[bi][aj] -= intvalue;
          A[bj][ai] -= intvalue;
	}
      }  
      else if (p!=q && r!=s && pq==rs) {
	if (p_ < ndocc && q_ >= ndocc && r_ < ndocc && s_ >= ndocc) {
          ai = (q_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + r_;
          A[ai][bj] += intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ < ndocc && s_ >= ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (s_-ndocc) * ndocc + r_;
          A[ai][bj] += intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ >= ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (r_-ndocc) * ndocc + s_;
          A[ai][bj] += intvalue;
	}
	else if (p_ >= ndocc && q_ >= ndocc && r_ < ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + r_;
          bj = (q_-ndocc) * ndocc + s_;
          aj = (p_-ndocc) * ndocc + s_;
          bi = (q_-ndocc) * ndocc + r_;
          A[ai][bj] -= intvalue;
          A[aj][bi] -= intvalue;
          A[bi][aj] -= intvalue;
          A[bj][ai] -= intvalue;
	}
      }  
      else if (p!=q && r==s) {
	if (p_ < ndocc && q_ < ndocc && r_ >= ndocc && s_ >= ndocc) {
          ai = (r_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + q_;
          aj = (r_-ndocc) * ndocc + q_;
          bi = (s_-ndocc) * ndocc + p_;
          A[ai][bj] -= intvalue;
          A[aj][bi] -= intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ >= ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (r_-ndocc) * ndocc + s_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ >= ndocc && r_ < ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + r_;
          bj = (q_-ndocc) * ndocc + s_;
          aj = (p_-ndocc) * ndocc + s_;
          bi = (q_-ndocc) * ndocc + r_;
          A[ai][bj] -= intvalue;
          A[bi][aj] -= intvalue;
	}
      }  
      else if (p==q && r!=s) {
	if (p_ < ndocc && q_ < ndocc && r_ >= ndocc && s_ >= ndocc) {
          ai = (r_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + q_;
          aj = (r_-ndocc) * ndocc + q_;
          bi = (s_-ndocc) * ndocc + p_;
          A[ai][bj] -= intvalue;
          A[bi][aj] -= intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ < ndocc && s_ >= ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (s_-ndocc) * ndocc + r_;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ >= ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (r_-ndocc) * ndocc + s_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ >= ndocc && r_ < ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + r_;
          bj = (q_-ndocc) * ndocc + s_;
          aj = (p_-ndocc) * ndocc + s_;
          bi = (q_-ndocc) * ndocc + r_;
          A[ai][bj] -= intvalue;
          A[aj][bi] -= intvalue;
	}
      }  
      else if (p==q && r==s && pq!=rs) {
	if (p_ < ndocc && q_ < ndocc && r_ >= ndocc && s_ >= ndocc) {
          ai = (r_-ndocc) * ndocc + p_;
          bj = (s_-ndocc) * ndocc + q_;
          A[ai][bj] -= intvalue;
	}
	else if (p_ >= ndocc && q_ < ndocc && r_ >= ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + q_;
          bj = (r_-ndocc) * ndocc + s_;
          A[ai][bj] += intvalue;
          A[bj][ai] += intvalue;
	}
	else if (p_ >= ndocc && q_ >= ndocc && r_ < ndocc && s_ < ndocc) {
          ai = (p_-ndocc) * ndocc + r_;
          bj = (q_-ndocc) * ndocc + s_;
          A[ai][bj] -= intvalue;
	}
      }  
    } 
    if (endflag == 0) {
      iwl_buf_fetch(&buffer);
    }
  } 
  
  iwl_buf_close(&buffer,1);

  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	   
        ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {
          
	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
               
	      for(jsym=0; jsym < nirreps; jsym++) {

	        jfirst = ofirst[jsym];
	        jlast = olast[jsym];

	        for(j=jfirst; j <= jlast; j++,BJ++) {
                  
	          A[AI][BJ] += (a==b) * (i==j) * (evals[a] - evals[i]);
                }
              }
            }
          }
        }
      }
    }
  }*/

  /* Diagonalize Magnetic Hessian */
  /*
  eval = init_array(num_ai);
  evec = block_matrix(num_ai,num_ai);
  sq_rsp(num_ai,num_ai,A,eval,0,evec,1E-14);
  fprintf(outfile, "\n\tMagnetic Hessian Eigenvalues\n\n");
  for(i=0; i<num_ai; i++) fprintf(outfile,"\t%d\t%12.8lf\n",i,eval[i]);
  */
  
  /*
  fprintf(outfile, "\nMO Hessian:\n");
  print_mat(A, num_ai, num_ai, outfile);
  */

  /*
  for(i=0; i<nuocc*ndocc; i++) 
    for(j=0; j<nuocc*ndocc; j++)
      tmp += A[i][j] - A0[i][j];

  fprintf(outfile,"\nDiff %20.10lf\n",tmp);

  fflush(outfile); 
  */

  /* Allocate space for the B vectors */
  B = (double ***)malloc(3 * sizeof(double **));
  for(coord=0; coord < 3; coord++) B[coord] = block_matrix(nmo,nmo);

  /* Transform magnetic dipole integrals to the MO basis */
  noei_ao = nao*(nao+1)/2;
  scratch = init_array(noei_ao);
  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);

  iwl_rdone(PSIF_OEI, PSIF_AO_LX, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * scratch[ij];
      TMP[j][i] = +0.5 * scratch[ij];
    }

  /*
  fprintf(outfile, "\tAO-basis LX Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  /*
  fprintf(outfile, "\n\tSO-basis LX Integrals:\n");
  print_mat(TMP, nso, nso, outfile);
  */

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[0][0][0]),nmo);

  zero_arr(scratch,noei_ao);
  zero_mat(TMP,nao,nao);
  zero_mat(X,nao,nao);

  /*
  fprintf(outfile, "\n\tMO-basis LX Integrals:\n");
  print_mat(B[0], nmo, nmo, outfile);
  */

  iwl_rdone(PSIF_OEI, PSIF_AO_LY, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * scratch[ij];
      TMP[j][i] = +0.5 * scratch[ij];
    }

  /*
  fprintf(outfile, "\n\tAO-basis LY Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  /*
  fprintf(outfile, "\n\tSO-basis LY Integrals:\n");
  print_mat(TMP, nso, nso, outfile);
  */

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[1][0][0]),nmo);

  zero_arr(scratch,noei_ao);
  zero_mat(TMP,nao,nao);
  zero_mat(X,nao,nao);

  /*
  fprintf(outfile, "\n\tMO-basis LY Integrals:\n");
  print_mat(B[1], nmo, nmo, outfile);
  */

  iwl_rdone(PSIF_OEI, PSIF_AO_LZ, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * scratch[ij];
      TMP[j][i] = +0.5 * scratch[ij];
    }

  /*
  fprintf(outfile, "\n\tAO-basis LZ Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  /*
  fprintf(outfile, "\n\tSO-basis LZ Integrals:\n");
  print_mat(TMP, nso, nso, outfile);
  */

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[2][0][0]),nmo);

  /*
  fprintf(outfile, "\n\tMO-basis LZ Integrals:\n");
  print_mat(B[2], nmo, nmo, outfile);
  */

  free(scratch);
  free_block(TMP);
  free_block(X);

  /*
  fprintf(outfile, "\n\tSCF MO's:\n");
  print_mat(scf, nso, nmo, outfile);
  fprintf(outfile, "\n\tMO-basis LX Integrals:\n");
  print_mat(B[0], nmo, nmo, outfile);
  fprintf(outfile, "\n\tMO-basis LY Integrals:\n");
  print_mat(B[1], nmo, nmo, outfile);
  fprintf(outfile, "\n\tMO-basis LZ Integrals:\n");
  print_mat(B[2], nmo, nmo, outfile);
  */

  /* Sort the B's into vector storage */
  B0 = block_matrix(3, num_ai);
  for(coord=0; coord < 3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym = 0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    B0[coord][AI] = B[coord][a][i];
	  }
	}
      }
    }
  }
  for(coord=0; coord<3; coord++)
    free_block(B[coord]);
  free(B);

  /* Solve the CPHF equations */
  ipiv = init_int_array(num_ai);
  for(coord=0; coord < 3; coord++) 
    error = C_DGESV(num_ai,1,&(A[0][0]),num_ai,&(ipiv[0]),&(B0[coord][0]),num_ai);

  free_block(A);

  /* Sort the UB matrices to matrix form */
  UB = (double ***) malloc(3 * sizeof(double **));
  for(coord=0; coord < 3; coord++)  
    UB[coord] = block_matrix(nmo,nmo);

  for(coord=0; coord < 3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym=0; isym < nirreps; isym++) {

	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    UB[coord][a][i] = B0[coord][AI];
	  }
	}
      }
    }
  }
  free_block(B0);

  /*
  fprintf(outfile,"\n\tUB[X]\n\n");
    print_mat(UB[0],nmo,nmo,outfile);

  fprintf(outfile,"\n\tUB[Y]\n\n");
    print_mat(UB[1],nmo,nmo,outfile);

  fprintf(outfile,"\n\tUB[Z]\n\n");
    print_mat(UB[2],nmo,nmo,outfile);
  */

  S = (double ***)malloc(natom*3*sizeof(double**));
  for(coord=0; coord<natom*3; coord++) {
    S[coord] = block_matrix(nmo,nmo);
  }
  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  label = (char *)malloc(PSIO_KEYLEN*sizeof(char)); 

  for(coord=0; coord < natom*3; coord++) {

    sprintf(label, "AO-basis Half-Diff Overlap (%d)", coord);
    psio_open(PSIF_OEI, PSIO_OPEN_OLD);
    psio_read_entry(PSIF_OEI, label, (char *) &(TMP[0][0]), 
                    nao*nao*sizeof(double));
    psio_close(PSIF_OEI, 1);

    /* Transform from AO -> SO */
    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    zero_mat(TMP,nao,nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    /* Transform from SO -> MO */
    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    zero_mat(TMP,nao,nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(S[coord][0][0]),nmo);
    
    //fprintf(outfile,"\nHalf-Diff S[%d] (MO):\n",coord);
    //print_mat(S[coord],nmo,nmo,outfile);
  }

  free_block(TMP);
  free_block(X);

  /* Atomic Axial Tensor */
  I = block_matrix(3, natom*3);  

  for(coord=0; coord < natom*3; coord++) {
    for(asym=0; asym < nirreps; asym++) {
      afirst = vfirst[asym];
      alast = vlast[asym];
      for(a=afirst; a <= alast; a++) {
        for(isym=0; isym < nirreps; isym++) {
          ifirst = ofirst[isym];
          ilast = olast[isym];
          for(i=ifirst; i <= ilast; i++) {
            I[0][coord] += 2.0 * UX[coord][a][i] * UB[0][a][i];
            I[1][coord] += 2.0 * UX[coord][a][i] * UB[1][a][i];
            I[2][coord] += 2.0 * UX[coord][a][i] * UB[2][a][i];
          }
        }
      } 
    }
  }

  for(coord=0; coord < natom*3; coord++) {
    for(asym=0; asym < nirreps; asym++) {
      afirst = vfirst[asym];
      alast = vlast[asym];
      for(a=afirst; a <= alast; a++) {
        for(isym=0; isym < nirreps; isym++) {
          ifirst = ofirst[isym];
          ilast = olast[isym];
          for(i=ifirst; i <= ilast; i++) {
            I[0][coord] += 2.0 * UB[0][a][i] * S[coord][i][a];
            I[1][coord] += 2.0 * UB[1][a][i] * S[coord][i][a];
            I[2][coord] += 2.0 * UB[2][a][i] * S[coord][i][a];
          }
        }
      } 
    }
  }
  
  /*
  for(coord=0; coord < natom*3; coord++) {
    I[0][coord] *= _bohr2angstroms * _aumfd2tesla * 1E-6;
    I[1][coord] *= _bohr2angstroms * _aumfd2tesla * 1E-6;
    I[2][coord] *= _bohr2angstroms * _aumfd2tesla * 1E-6;
  }

  fprintf(outfile,"\nAtomic Axial Tensor I (10^-6 Angstrom^-1 T^-1):\n");
  print_mat(I,3,natom*3,outfile);
  */

  /* 
  Nuclear Contribution to the AAT:
    Levi-Civita tensor epsilon_{ijk} 
    i=0 j=1 k=2 epsilon_{012}=+1  
    i=0 j=2 k=1 epsilon_{021}=-1  
    i=1 j=0 k=2 epsilon_{102}=-1  
    i=2 j=1 k=0 epsilon_{210}=-1  
    i=1 j=2 k=0 epsilon_{120}=+1  
    i=2 j=0 k=1 epsilon_{201}=+1  
  */

  fprintf(outfile,"\n\tElectronic Contributions to AAT (au):\n");
  fprintf(outfile,"\n\t     Bx\t\t     By\t\t     Bz\n");
  for(i=0; i<natom*3; i++) { 
    if((i%3)==0)
      fprintf(outfile,"%3sx\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i%3)==1)
      fprintf(outfile,"%3sy\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i%3)==2)
      fprintf(outfile,"%3sz\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i+1)%3==0) fprintf(outfile,"\n");
  }

  for(coord=0; coord < natom; coord++) {
    I[1][coord*3+0] += 0.25 * zvals[coord] * geom[coord][2];
    I[2][coord*3+0] -= 0.25 * zvals[coord] * geom[coord][1];
    I[0][coord*3+1] -= 0.25 * zvals[coord] * geom[coord][2];
    I[1][coord*3+2] -= 0.25 * zvals[coord] * geom[coord][0];
    I[2][coord*3+1] += 0.25 * zvals[coord] * geom[coord][0];
    I[0][coord*3+2] += 0.25 * zvals[coord] * geom[coord][1];
  }

  fprintf(outfile,"\n\tAtomic Axial Tensor:\n");
  fprintf(outfile,"\n\tUnits: au\n");
  fprintf(outfile,"\n\tTerms: Electronic + Nuclear\n");
  fprintf(outfile,"\n\t     Bx\t\t     By\t\t     Bz\n");
  for(i=0; i<natom*3; i++) { 
    if((i%3)==0)
      fprintf(outfile,"%3sx\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i%3)==1)
      fprintf(outfile,"%3sy\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i%3)==2)
      fprintf(outfile,"%3sz\t%10.6lf\t%10.6lf\t%10.6lf\n",
              asymbol[i],I[0][i],I[1][i],I[2][i]);
    if((i+1)%3==0) fprintf(outfile,"\n");
  }

  /* Transform the atomic axial tensor to normal coordinates*/
  I_q = block_matrix(3, natom*3);  

  C_DGEMM('n','n',3,natom*3,natom*3,1,&(I[0][0]),natom*3,
          &(lx[0][0]),natom*3,0,&(I_q[0][0]),natom*3);

  fprintf(outfile,"\n\tAtomic Axial Tensor:\n");
  fprintf(outfile,"\n\tUnits: au\n");
  fprintf(outfile,"\n\tTerms: Electronic + Nuclear\n");
  fprintf(outfile,"\n\t     Bx\t\t     By\t\t     Bz\n");
  for(i=0; i<nnc; i++) { 
    fprintf(outfile,"  Q%d\t%10.6lf\t%10.6lf\t%10.6lf\n",i+1,
            I_q[0][i],I_q[1][i],I_q[2][i]);
  }

  /* Rotational Strength */

  rotstr = init_array(natom*3);

  for(i=(3*natom-1); i >= (3*natom-nnc); i--) 
    for(j=0; j<3; j++)
      rotstr[i] += dipder_q[j][i] * I_q[j][i];
  
  fprintf(outfile,"\n\tRotational Strength (au):\n\n");
  for(i=natom*3-1; i >= (natom*3-nnc); i--) {
    fprintf(outfile,"  Q%d\t%20.15lf\n",(3*natom-i),rotstr[i]);
  }

  fprintf(outfile,"\n\tRotational Strength (10^-44 esu^2 cm^2):\n\n");
  for(i=natom*3-1; i >= (natom*3-nnc); i--) {
    tmp = (-1e-36*_bohr2angstroms*_h*rotstr[i]*1e10*1e44)/(_amu2kg*2*_pi*_c); 
    fprintf(outfile,"  Q%d\t%20.15lf\n",(3*natom-i),tmp);
  }
  
  for(coord=0; coord < 3; coord++) {
    free_block(UB[coord]);
    free_block(S[coord]);
  }
  free_block(I);
  free_block(I_q);
  free_block(dipder_q);
  free(label);
}

}} // namespace psi::cphf
