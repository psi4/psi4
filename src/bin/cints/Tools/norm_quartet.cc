/*! \file norm_quartet.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<cstdlib>
#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"transform.h"

namespace psi { namespace CINTS {

/*!------------------------------------------------------------------------------------
  This functions processes the computed integrals quartet:
  1) Normalizes cartesian components to unity;
  2) Transforms integrals from cartesian to pure angular momentum basis (if needed);
 ------------------------------------------------------------------------------------*/
double *norm_quartet(double *data, double *puream_data, int am[4], int puream)
{
  int ii, jj, kk, ll;
  int ni, nj, nk, nl;
  double *tmp_ptr;
  const double *ptr_i, *ptr_j, *ptr_k, *ptr_l;
  double norm_i, norm_ij, norm_ijk;
#pragma disjoint (*tmp_ptr,norm_i)
#pragma disjoint (*tmp_ptr,norm_ij)
#pragma disjoint (*tmp_ptr,norm_ijk)


  /*--- numbers of integrals in each shell ---*/
  ni = ioff[am[0] + 1];
  nj = ioff[am[1] + 1];
  nk = ioff[am[2] + 1];
  nl = ioff[am[3] + 1];

  /*------------------------------------------------------------------------------------
    Normalize contracted integrals - right here each cartesian component in the shell
    has the same normalization coefficient so that only components with radial parts of
    x^l, y^l, and z^l are normalized to unity. After this block of code all basis
    functions are normalized to unity. Needed this so that integrals in terms of
    puream i-functions were computed properly.

    N.B.: possibly, LIBINT should be rewritten so that all intermediates are computed in
    terms of normalized functions?!
   ------------------------------------------------------------------------------------*/
  ptr_i = GTOs.bf_norm[am[0]];
  ptr_j = GTOs.bf_norm[am[1]];
  ptr_k = GTOs.bf_norm[am[2]];
  ptr_l = GTOs.bf_norm[am[3]];
  
  tmp_ptr = data;
  for(ii=0; ii<ni; ii++) {
    norm_i = ptr_i[ii];
    for(jj=0; jj<nj; jj++) {
      norm_ij = norm_i*ptr_j[jj];
      for(kk=0; kk<nk; kk++) {
	norm_ijk = norm_ij*ptr_k[kk];
	for(ll=0; ll<nl; ll++) {
	  *(tmp_ptr++) *= norm_ijk*ptr_l[ll];
	}
      }
    }
  }


  /*--------------------------------------------------
    Cartesian to pure angular momentum transformation
   --------------------------------------------------*/
  if (puream) {
#if !USE_MM
    if (am[0]) {
      ni = 2*am[0]+1;
      transform_i(data,puream_data,GTOs.cart2pureang[am[0]],am[0],nj,nk,nl);
      tmp_ptr = data;
      data = puream_data;
      puream_data = tmp_ptr;
    }
    if (am[1]) {
      nj = 2*am[1]+1;
      transform_j(data,puream_data,GTOs.cart2pureang[am[1]],am[1],ni,nk,nl);
      tmp_ptr = data;
      data = puream_data;
      puream_data = tmp_ptr;
    }
    if (am[2]) {
      nk = 2*am[2]+1;
      transform_k(data,puream_data,GTOs.cart2pureang[am[2]],am[2],ni,nj,nl);
      tmp_ptr = data;
      data = puream_data;
      puream_data = tmp_ptr;
    }
    if (am[3]) {
      nl = 2*am[3]+1;
      transform_l(data,puream_data,GTOs.cart2pureang[am[3]],am[3],ni,nj,nk);
      tmp_ptr = data;
      data = puream_data;
      puream_data = tmp_ptr;
    }
#elif USE_MM && !SPARSE_C2P
    /*--- N.B. am[0] >= am[1] ---*/
    if (am[0])
      data = transform_ijkl(data,puream_data,GTOs.cc2pp[am[0]][am[1]],
			    GTOs.cc2pp[am[2]][am[3]],1,am[0],am[1],am[2],am[3]);
    else
      data = transform_ijkl(data,puream_data,NULL,GTOs.cc2pp[am[2]][am[3]],0,
			    am[0],am[1],am[2],am[3]);
#elif USE_MM && SPARSE_C2P
    /*--- N.B. am[0] >= am[1] ---*/
    if (am[0])
      data = transform_ijkl_sparse(data,puream_data,
				   GTOs.pp2cc_sparse[am[0]][am[1]],
				   GTOs.cc2pp_sparse[am[2]][am[3]],
				   GTOs.pp2cc_rowlength[am[0]][am[1]],
				   GTOs.cc2pp_rowlength[am[2]][am[3]],
				   1,am[0],am[1],am[2],am[3]);
    else
      data = transform_ijkl_sparse(data,puream_data,NULL,
				   GTOs.cc2pp_sparse[am[2]][am[3]],
				   NULL,GTOs.cc2pp_rowlength[am[2]][am[3]],
				   0,am[0],am[1],am[2],am[3]);
#endif
  }

  return data;
}

};};
