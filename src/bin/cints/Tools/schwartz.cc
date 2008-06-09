/*! \file schwartz.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<cstdlib>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>

#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"quartet_data.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include"schwartz.h"

namespace psi { namespace CINTS {

void schwartz_eri()
{
  /*--- Various data structures ---*/
  struct shell_pair *sp_ij, *sp_kl;

  Libint_t Libint;
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;
#endif
  
  int ij, kl, ik, jl, ijkl;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int g, i, j, k, l, m, ii, jj, kk, ll;
  int a, b, c, d;
  int amax, bmax, cmax, dmax;
  int si;                                  /* GCC compiler screwes up if static is left out */
  int sj, sk, sl, si_g, sj_g;
  int sii, sjj, skk, sll , slll;
  int sij, skl, sijkl;
  int pi, pj, pk, pl ;
  int max_pj, max_pl;
  int pii, pjj, pkk, pll;
  int upk, num_unique_pk;
  int max_num_unique_quartets;
  int max_num_prim_comb;
  int max_bf_per_shell;

  int size, class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, li, lj;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;
  int effective = 0;         /* Flag to check the effectiveness of prescreening */

  double AB2, CD2;
  double *raw_data;             /* pointer to the unnormalized taregt quartet of integrals */
  double *data;                 /* pointer to the transformed normalized target quartet of integrals */
#ifdef NONDOUBLE_INTS
  REALTYPE *target_ints;        /* Pointer to the location of the target quartet on the stack of
			 	   integrals quartets if libint.a is using other than regular doubles */
#endif
  double temp;
  double max_elem;
  

  BasisSet.schwartz_eri = block_matrix(BasisSet.num_shells,BasisSet.num_shells);
  
  /*---------------
    Initialization
   ---------------*/
/*  init_fjt(BasisSet.max_am*4);*/
#ifndef USE_TAYLOR_FM
  init_fjt_table(&fjt_table);
#endif
  max_bf_per_shell = ioff[BasisSet.max_am];
  max_cart_class_size = (max_bf_per_shell)*
                        (max_bf_per_shell)*
                        (max_bf_per_shell)*
                        (max_bf_per_shell);
#ifdef NONDOUBLE_INTS
  raw_data = init_array(max_cart_class_size);
#endif
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
		      (BasisSet.max_num_prims*BasisSet.max_num_prims);
  UserOptions.memory -= init_libint(&Libint,BasisSet.max_am-1,max_num_prim_comb);

/*-------------------------------------------------
  generate all shell quartets 
  suitable for building the PK-matrix
 -------------------------------------------------*/
  for (sii=0; sii<BasisSet.num_shells; sii++)
    for (sjj=0; sjj<=sii; sjj++) {
	
      si = sii;
      sj = sjj;
      sk = sii;
      sl = sjj;

      if(BasisSet.shells[sii].am < BasisSet.shells[sjj].am){
	  dum = si;
	  si = sj;
	  sj = dum;
	  dum = sk;
	  sk = sl;
	  sl = dum;
      }
      ni = ioff[BasisSet.shells[si].am];
      nj = ioff[BasisSet.shells[sj].am];
      nk = ioff[BasisSet.shells[sk].am];
      nl = ioff[BasisSet.shells[sl].am];
      class_size = ni*nj*nk*nl;
      np_i = BasisSet.shells[si].n_prims;
      np_j = BasisSet.shells[sj].n_prims;
      np_k = BasisSet.shells[sk].n_prims;
      np_l = BasisSet.shells[sl].n_prims;
      orig_am[0] = BasisSet.shells[si].am-1;
      orig_am[1] = BasisSet.shells[sj].am-1;
      orig_am[2] = BasisSet.shells[sk].am-1;
      orig_am[3] = BasisSet.shells[sl].am-1;
      am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];

      sp_ij = &(BasisSet.shell_pairs[si][sj]);
      sp_kl = &(BasisSet.shell_pairs[sk][sl]);
      
      Libint.AB[0] = sp_ij->AB[0];
      Libint.AB[1] = sp_ij->AB[1];
      Libint.AB[2] = sp_ij->AB[2];
      Libint.CD[0] = sp_kl->AB[0];
      Libint.CD[1] = sp_kl->AB[1];
      Libint.CD[2] = sp_kl->AB[2];
		  
      AB2 = Libint.AB[0]*Libint.AB[0]+
	    Libint.AB[1]*Libint.AB[1]+
	    Libint.AB[2]*Libint.AB[2];
      CD2 = Libint.CD[0]*Libint.CD[0]+
	    Libint.CD[1]*Libint.CD[1]+
	    Libint.CD[2]*Libint.CD[2];

      /*--- Compute data for primitive quartets here ---*/
      num_prim_comb = 0;
      for (pi = 0; pi < np_i; pi++) {
	max_pj = (si == sj) ? pi+1 : np_j;
	for (pj = 0; pj < max_pj; pj++) {
	  m = (1 + (si == sj && pi != pj));
	  for (pk = 0; pk < np_k; pk++) {
	    max_pl = (sk == sl) ? pk+1 : np_l;
	    for (pl = 0; pl < max_pl; pl++){
	      n = m * (1 + (sk == sl && pk != pl));
#ifdef USE_TAYLOR_FM
	      quartet_data(&(Libint.PrimQuartet[num_prim_comb++]), NULL, AB2, CD2,
			   sp_ij, sp_kl, am, pi, pj, pk, pl, (double)n);
#else
	      quartet_data(&(Libint.PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
			   sp_ij, sp_kl, am, pi, pj, pk, pl, (double)n);
#endif
	    }
	  }
	}
      }

	      /*-------------------------------------------------
		Compute the quartet and find the largest element
	       -------------------------------------------------*/
      if (am) {
#ifdef NONDOUBLE_INTS
	  size = ioff[BasisSet.shells[si].am]*ioff[BasisSet.shells[sj].am]*
		 ioff[BasisSet.shells[sk].am]*ioff[BasisSet.shells[sl].am];
	  target_ints = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
	  for(i=0;i<size;i++)
	      raw_data[i] = (double) target_ints[i];
#else
	  raw_data = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
#endif
	  /*                    zero here means no transformation to puream basis */
	  /*                                    |                                 */
	  data = norm_quartet(raw_data, NULL, orig_am, 0);

	  max_elem = -1.0;
	  iimax = ni - 1;
	  for(ii=0; ii <= iimax; ii++){
	    jjmax = nj - 1;
	    for(jj=0; jj <= jjmax; jj++){
	      ijkl = jj+nl*(ii+nk*(jj+nj*ii));
	      if(fabs(data[ijkl]) > max_elem){
		max_elem = fabs(data[ijkl]);
	      }
	    }
	  }
      }
      else {
	  temp = 0.0;
	  for(p=0;p<num_prim_comb;p++)
	      temp += (double) Libint.PrimQuartet[p].F[0];
	  
	  max_elem = fabs(temp);
      }

      BasisSet.schwartz_eri[sii][sjj] = BasisSet.schwartz_eri[sjj][sii] = sqrt(max_elem);
      if (max_elem < UserOptions.cutoff)
	/*--- At least one shell quartet will be eliminated ---*/
	effective = 1;

    } /* end getting shell combination */

  
  /*---------
    Clean-up
   ---------*/
  free_libint(&Libint);
#ifndef USE_TAYLOR_FM
  free_fjt_table(&fjt_table);
#endif
#ifdef NONDOUBLE_INTS
  free(raw_data);
#endif

  if (effective == 0)
    fprintf(outfile,"  Schwartz prescreening of ERIs will be ineffective\n");

  return;
}

};};
