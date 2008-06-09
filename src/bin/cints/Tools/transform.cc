/*! \file transform.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
    -------------------------------------------------------------------------------------------------
    Functions for performing four-index transformation from cartesian to pure angular momentum basis
    -------------------------------------------------------------------------------------------------*/

#include<cstdio>
#include<memory.h>
#include<libciomr/libciomr.h>
#include<libqt/qt.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

#if !USE_MM
void transform_i(double *data, double *puream_data, double **c2p, int am_i, int nj, int nk, int nl)
{
  register int i,j,k,l;
  int step = nj*nk*nl;
  int offset;
  double *ii;
  int ni, nsph, isph;
  double tmp;

  ni = ioff[am_i + 1];
  nsph = 2*am_i + 1;
  
  for(j=0;j<nj;j++)
    for(k=0;k<nk;k++) {
      offset = (j*nk+k)*nl;
      for(l=0;l<nl;l++,offset++)
	for(isph=0;isph<nsph;isph++) {
	  ii = data + offset; tmp = 0.0;
	  for(i=0;i<ni;i++,ii+=step)
	    tmp += c2p[isph][i]*(*ii);
	  puream_data[step*isph+offset] = tmp;
	}
    }

  return;
}
	  

void transform_j(double *data, double *puream_data, double **c2p, int am_j, int ni, int nk, int nl)
{
  register int i,j,k,l;
  int step = nk*nl;
  int offset,offset_sph;
  double *jj;
  int nj, nsph, jsph;
  double tmp;

  nj = ioff[am_j + 1];
  nsph = 2*am_j + 1;
  
  for(i=0;i<ni;i++)
    for(k=0;k<nk;k++,offset+=nl) {
      offset = (i*nj*nk + k)*nl;
      offset_sph = (i*nsph*nk + k)*nl;
      for(l=0;l<nl;l++,offset++,offset_sph++)
	for(jsph=0;jsph<nsph;jsph++) {
	  jj = data + offset; tmp = 0.0;
	  for(j=0;j<nj;j++,jj+=step)
	    tmp += c2p[jsph][j]*(*jj);
	  puream_data[step*jsph+offset_sph] = tmp;
	}
    }

  return;
}



void transform_k(double *data, double *puream_data, double **c2p, int am_k, int ni, int nj, int nl)
{
  register int i,j,k,l;
  int step = nl;
  int offset,offset_sph;
  double *kk;
  int nk, nsph, ksph;
  double tmp;

  nk = ioff[am_k + 1];
  nsph = 2*am_k + 1;
  
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++) {
      offset = (i*nj + j)*nk*nl;
      offset_sph = (i*nj + j)*nsph*nl;
      for(l=0;l<nl;l++,offset++,offset_sph++)
	for(ksph=0;ksph<nsph;ksph++) {
	  kk = data + offset; tmp = 0.0;
	  for(k=0;k<nk;k++,kk+=step)
	    tmp += c2p[ksph][k]*(*kk);
	  puream_data[step*ksph+offset_sph] = tmp;
	}
    }

  return;
}



void transform_l(double *data, double *puream_data, double **c2p, int am_l, int ni, int nj, int nk)
{
  register int i,j,k,l;
  int step = 1;
  int offset,offset_sph;
  double *ll;
  int nl, nsph, lsph;
  double tmp;

  nl = ioff[am_l + 1];
  nsph = 2*am_l + 1;
  
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++) {
      offset = (i*nj + j)*nk*nl;
      offset_sph = (i*nj + j)*nk*nsph;
      for(k=0;k<nk;k++,offset+=nl,offset_sph+=nsph)
	for(lsph=0;lsph<nsph;lsph++) {
	  ll = data + offset; tmp = 0;
	  for(l=0;l<nl;l++,ll+=step)
	    tmp += c2p[lsph][l]*(*ll);
	  puream_data[step*lsph+offset_sph] = tmp;
	}
    }

}

#elif USE_MM && !SPARSE_C2P

double *transform_ijkl(double *data, double *puream_data, double **cc2pp_ij, double **cc2pp_kl, int trans_ij,
		       int am_i, int am_j, int am_k, int am_l)
{
  int ij,kl,ijkl;
  int n_ij, n_kl, nsph_ij, nsph_kl;
  double tmp;
#if USE_BLAS
  int zero = 0;
  int one = 1;
  double zero_d = 0.0;
  double one_d = 1.0;
  char true_c = 't';
  char false_c = 'n';
#else
  double *tmp_ptr;
  static double *cij_ckl[CINTS_MAX_AM*CINTS_MAX_AM*(CINTS_MAX_AM+1)*(CINTS_MAX_AM+1)/4];
  static double *cij_pkl[CINTS_MAX_AM*(CINTS_MAX_AM+1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)/2];
  static double *pij_pkl[(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)];
#endif

  n_ij = ioff[am_i+1]*ioff[am_j+1];
  n_kl = ioff[am_k+1]*ioff[am_l+1];
#if !USE_BLAS
  for(ij=0, tmp_ptr = data; ij<n_ij; ij++, tmp_ptr += n_kl)
    cij_ckl[ij] = tmp_ptr;
#endif
  
  nsph_kl = (2*am_k+1)*(2*am_l+1);
#if !USE_BLAS
  for(ij=0, tmp_ptr = puream_data; ij<n_ij; ij++, tmp_ptr += nsph_kl)
    cij_pkl[ij] = tmp_ptr;
  /* Insert any C-based matrix multiplication routine here */
  newmm_rking(&cij_ckl[0],0,cc2pp_kl,1,&cij_pkl[0],n_ij,n_kl,nsph_kl,1.0,0.0);
#else
/*  F_DGEMM(&true_c,&false_c,&nsph_kl,&n_ij,&n_kl,&one_d,&(cc2pp_kl[0][0]),&n_kl,
          data,&n_kl,&zero_d,puream_data,&nsph_kl);*/
    C_DGEMM('n', 't', n_ij, nsph_kl, n_kl, 1.0, data, n_kl, cc2pp_kl[0], n_kl,
            0.0, puream_data, nsph_kl);
#endif

  if (trans_ij) {
    nsph_ij = (2*am_i+1)*(2*am_j+1);
#if !USE_BLAS
    for(ij=0, tmp_ptr = data; ij<nsph_ij; ij++, tmp_ptr += nsph_kl)
      pij_pkl[ij] = tmp_ptr;
    /* Insert any C-based matrix multiplication routine here */
    newmm_rking(cc2pp_ij,0,&cij_pkl[0],0,&pij_pkl[0],nsph_ij,n_ij,nsph_kl,1.0,0.0);
#else
/*  F_DGEMM(&false_c,&false_c,&nsph_kl,&nsph_ij,&n_ij,&one_d,cc2pp_ij[0],&n_ij,
          puream_data,&nsph_kl,&zero_d,data,&nsph_kl);*/
    C_DGEMM('n', 'n', nsph_ij, nsph_kl, n_ij, 1.0, cc2pp_ij[0], n_ij,
            puream_data, nsph_kl, 0.0, data, nsph_kl);

#endif

    return data;
  }
  else
    return puream_data;

}

#elif USE_MM && SPARSE_C2P

double *transform_ijkl_sparse(double *data, double *puream_data,
		       mat_elem **pp2cc_sparse, mat_elem **cc2pp_sparse,
		       int *pp2cc_rowlength, int *cc2pp_rowlength, int trans_ij,
		       int am_i, int am_j, int am_k, int am_l)
{
  int ij,kl,ijkl,pq;
  int n_ij, n_kl, nsph_ij, nsph_kl;
  double tmp;

  double *tmp_ptr;
/*  static double *cij_ckl[CINTS_MAX_AM*CINTS_MAX_AM*(CINTS_MAX_AM+1)*(CINTS_MAX_AM+1)/4];
  static double *cij_pkl[CINTS_MAX_AM*(CINTS_MAX_AM+1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)/2];
  static double *pij_pkl[(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)*(2*CINTS_MAX_AM-1)];*/
  double *cijckl_row, *cijpkl_row, *pijpkl_row;

  n_ij = ioff[am_i+1]*ioff[am_j+1];
  n_kl = ioff[am_k+1]*ioff[am_l+1];
  nsph_kl = (2*am_k+1)*(2*am_l+1);
  for(ij=0, cijckl_row = data, cijpkl_row = puream_data; ij<n_ij;
      ij++, cijckl_row += n_kl, cijpkl_row += nsph_kl) {
    memset(cijpkl_row,0,sizeof(double)*nsph_kl);
    for(pq=0; pq<n_kl; pq++) {
      tmp = cijckl_row[pq];
      for(kl=0; kl<cc2pp_rowlength[pq]; kl++) {
	cijpkl_row[cc2pp_sparse[pq][kl].column] += tmp*cc2pp_sparse[pq][kl].value;
      }
    }
  }
			 
  if (trans_ij) {
    nsph_ij = (2*am_i+1)*(2*am_j+1);
    for(ij=0, pijpkl_row = data; ij<nsph_ij; ij++, pijpkl_row += nsph_kl) {
      memset(pijpkl_row,0,sizeof(double)*nsph_kl);
      for(pq=0; pq<pp2cc_rowlength[ij]; pq++) {
	tmp = pp2cc_sparse[ij][pq].value;
	cijpkl_row = puream_data + nsph_kl*pp2cc_sparse[ij][pq].column;
	for(kl=0; kl<nsph_kl; kl++) {
	  pijpkl_row[kl] += tmp*cijpkl_row[kl];
	}
      }
    }
    return data;
  }

  return puream_data;

}
#endif
};};
