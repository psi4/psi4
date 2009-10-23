/*! \file te_deriv2_scf.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <libqt/qt.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#include "moinfo.h"
#include "compute_scf_opdm.h"
#include "norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include "deriv1_quartet_data.h"
#include "symmetrize.h"
#include "small_fns.h"
#include <Tools/prints.h>

#define PK_ORDER 0
#define INCLUDE_1ST 1
#define INCLUDE_2ND 1

namespace psi { namespace cints {

void te_deriv2_scf(void)
{
  /*--- Various data structures ---*/
  struct shell_pair *sp_ij, *sp_kl;
  Libderiv_t Libderiv;                /* Integrals library object */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;               /* table of auxiliary function F_m(u) for each primitive combination */
#endif
  double **hess_te;

  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  register int i, j, k, l, m, ii, jj, kk, ll;
  register int si, sj, sk, sl ;
  register int sii, sjj, skk, sll, slll;
  register int pi, pj, pk, pl ;
  int max_pj, max_pl;
  register int pii, pjj, pkk, pll ;
  int switch_ij, switch_kl, switch_ijkl;
  int center_i, center_j, center_k, center_l, center[4];

  int class_size;
  int max_class_size;
  int max_cart_class_size;

  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, quartet_size;

  int num_prim_comb, p, max_num_prim_comb;
  double AB2, CD2;
  double *FourInd;
  int I, J, K, L, IJ, KL, this_quartet, num_unique_pk, unique;
  int si_fao, sj_fao, sk_fao, sl_fao;
  int coord, coord1, coord2;
  int di, dj;
  int si_arr[3], sj_arr[3], sk_arr[3], sl_arr[3];
  double *data, value[12], svalue;
  double fac1, fac2, fac3;
  double ffac1, ffac2, ffac3;
  double d2acc[12][12], contr, pfac, perm_pf;

  /*---------------
    Initialization
    ---------------*/
  hess_te = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+DERIV_LVL);
  init_fjt_table(&fjt_table);
#endif
  init_libderiv_base();
  /*  init_libint_base(); */

  max_cart_class_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  max_class_size = max_cart_class_size;
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
    (BasisSet.max_num_prims*BasisSet.max_num_prims);
  init_libderiv12(&Libderiv,BasisSet.max_am-1,max_num_prim_comb,max_cart_class_size);
  FourInd = init_array(max_cart_class_size);
  /*  init_libint(&Libint,BasisSet.max_am-1,max_num_prim_comb); */

#if !PK_ORDER
  for (sii=0; sii<BasisSet.num_shells; sii++)
    for (sjj=0; sjj<=sii; sjj++)
      for (skk=0; skk<=sii; skk++)
	for (sll=0; sll<= ((sii == skk) ? sjj : skk); sll++){
	  si = sii; sj = sjj; sk = skk; sl = sll;
#endif
#if PK_ORDER
  for (sii=0; sii<BasisSet.num_shells; sii++)
    for (sjj=0; sjj<=sii; sjj++)
      for (skk=0; skk<=sjj; skk++)
	for (sll=0; sll<= skk; sll++){

	  si_arr[0] = sii; sj_arr[0] = sjj; sk_arr[0] = skk; sl_arr[0] = sll;
	  if (sii == sjj && sii == skk || sjj == skk && sjj == sll)
	    num_unique_pk = 1;
	  else if (sii == skk || sjj == sll) {
	    num_unique_pk = 2;
	    si_arr[1] = sii; sj_arr[1] = skk; sk_arr[1] = sjj; sl_arr[1] = sll;
	  }
	  else if (sjj == skk) {
	    num_unique_pk = 2;
	    si_arr[1] = sii; sj_arr[1] = sll; sk_arr[1] = sjj; sl_arr[1] = skk;
	  }
	  else if (sii == sjj || skk == sll) {
	    num_unique_pk = 2;
	    si_arr[1] = sii; sj_arr[1] = skk; sk_arr[1] = sjj; sl_arr[1] = sll;
	  }
	  else {
	    num_unique_pk = 3;
	    si_arr[1] = sii; sj_arr[1] = skk; sk_arr[1] = sjj; sl_arr[1] = sll;
	    si_arr[2] = sii; sj_arr[2] = sll; sk_arr[2] = sjj; sl_arr[2] = skk;
	  }

	  for(unique=0; unique < num_unique_pk; unique++) {
	    si = si_arr[unique];
	    sj = sj_arr[unique];
	    sk = sk_arr[unique];
	    sl = sl_arr[unique];
#endif

	    fac1 = fac2 = fac3 = 1.0;
	    /* place in "ascending" angular mom-
	       my simple way of optimizing PHG recursion (VRR) */
	    /* these first two are good for the HRR */
	    if(BasisSet.shells[si].am < BasisSet.shells[sj].am){
	      dum = si;
	      si = sj;
	      sj = dum;
	    }
	    if(BasisSet.shells[sk].am < BasisSet.shells[sl].am){
	      dum = sk;
	      sk = sl;
	      sl = dum;
	    }
	    /* this should be /good/ for the VRR */
	    if(BasisSet.shells[si].am + BasisSet.shells[sj].am > BasisSet.shells[sk].am + BasisSet.shells[sl].am){
	      dum = si;
	      si = sk;
	      sk = dum;
	      dum = sj;
	      sj = sl;
	      sl = dum;
	    }

	    center_i = BasisSet.shells[si].center-1;
	    center_j = BasisSet.shells[sj].center-1;
	    center_k = BasisSet.shells[sk].center-1;
	    center_l = BasisSet.shells[sl].center-1;
#if 0
	    if (center_i == center_j ||
		center_i == center_k ||
		center_i == center_l ||
		center_j == center_k ||
		center_j == center_l ||
		center_k == center_l)
	      continue;
#endif

	    ni = ioff[BasisSet.shells[si].am];
	    nj = ioff[BasisSet.shells[sj].am];
	    nk = ioff[BasisSet.shells[sk].am];
	    nl = ioff[BasisSet.shells[sl].am];
	    quartet_size = ni*nj*nk*nl;

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

	    Libderiv.AB[0] = sp_ij->AB[0];
	    Libderiv.AB[1] = sp_ij->AB[1];
	    Libderiv.AB[2] = sp_ij->AB[2];
	    Libderiv.CD[0] = sp_kl->AB[0];
	    Libderiv.CD[1] = sp_kl->AB[1];
	    Libderiv.CD[2] = sp_kl->AB[2];
	    /*
	      Libint.AB[0] = sp_ij->AB[0];
	      Libint.AB[1] = sp_ij->AB[1];
	      Libint.AB[2] = sp_ij->AB[2];
	      Libint.CD[0] = sp_kl->AB[0];
	      Libint.CD[1] = sp_kl->AB[1];
	      Libint.CD[2] = sp_kl->AB[2];
	    */
	    
	    AB2 = Libderiv.AB[0]*Libderiv.AB[0]+
	      Libderiv.AB[1]*Libderiv.AB[1]+
	      Libderiv.AB[2]*Libderiv.AB[2];
	    CD2 = Libderiv.CD[0]*Libderiv.CD[0]+
	      Libderiv.CD[1]*Libderiv.CD[1]+
	      Libderiv.CD[2]*Libderiv.CD[2];

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
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb]),
					NULL, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n);
#endif    
		    num_prim_comb++;
		    
		  }
		}
	      }
	    }

	    /*	    for(i=0;i<40000;i++) */
	    build_deriv12_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libderiv,num_prim_comb);

#if INCLUDE_1ST
	    norm_quartet(Libderiv.ABCD[0], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[1], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[2], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[6], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[7], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[8], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[9], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[10], NULL, orig_am, 0);
	    norm_quartet(Libderiv.ABCD[11], NULL, orig_am, 0);

	    /* For all the integrals in this buffer, find the permutationally unique ones */
	    this_quartet = 0;
	    for(i=0; i < ni; i++) {
	      I = BasisSet.shells[si].fao + i - 1;
	      for(j=0; j < nj; j++) {
		J = BasisSet.shells[sj].fao + j - 1;
		IJ = INDEX(I,J);
		for(k=0; k < nk; k++) {
		  K = BasisSet.shells[sk].fao + k - 1;
		  for(l=0; l < nl; l++,this_quartet++) {
		    L = BasisSet.shells[sl].fao + l - 1;
		    KL = INDEX(K,L);

		    if(si == sj && I < J)
		      continue;
		    if(sk == sl && K < L)
		      continue;
		    if(INDEX(si,sj) == INDEX(sk,sl) && IJ < KL)
		      continue;

		    value[0] = Libderiv.ABCD[0][this_quartet];
		    value[1] = Libderiv.ABCD[1][this_quartet];
		    value[2] = Libderiv.ABCD[2][this_quartet];
		    value[6] = Libderiv.ABCD[6][this_quartet];
		    value[7] = Libderiv.ABCD[7][this_quartet];
		    value[8] = Libderiv.ABCD[8][this_quartet];
		    value[9] = Libderiv.ABCD[9][this_quartet];
		    value[10] = Libderiv.ABCD[10][this_quartet];
		    value[11] = Libderiv.ABCD[11][this_quartet];

		    /* translational invariance */
		    value[3] = - value[0] - value[6] - value[9];
		    value[4] = - value[1] - value[7] - value[10];
		    value[5] = - value[2] - value[8] - value[11];

		    ffac1 = ffac2 = ffac3 = 1.0;

		    if(I!=J) ffac1 *= 2.0;
		    if(K!=L) ffac1 *= 2.0;
		    if(I!=K) ffac2 *= 2.0;
		    if(J!=L) ffac2 *= 2.0;
		    if(I!=L) ffac3 *= 2.0;
		    if(J!=K) ffac3 *= 2.0;
		    if(INDEX(I,J) != INDEX(K,L)) ffac1 *= 2.0;
		    if(INDEX(I,K) != INDEX(J,L)) ffac2 *= 2.0;
		    if(INDEX(I,L) != INDEX(J,K)) ffac3 *= 2.0;

		    if((I==J && I==K) || (J==K && J==L) || (I==J && I==L) || (I==K && I==L)) {
		      F[center_i*3][I][J] += Dens[K][L] * ffac1 * 0.5 * value[0];
		      F[center_i*3][K][L] += Dens[I][J] * ffac1 * 0.5 * value[0];

		      F[center_i*3+1][I][J] += Dens[K][L] * ffac1 * 0.5 * value[1];
		      F[center_i*3+1][K][L] += Dens[I][J] * ffac1 * 0.5 * value[1];

		      F[center_i*3+2][I][J] += Dens[K][L] * ffac1 * 0.5 * value[2];
		      F[center_i*3+2][K][L] += Dens[I][J] * ffac1 * 0.5 * value[2];

		      F[center_j*3][I][J] += Dens[K][L] * ffac1 * 0.5 * value[3];
		      F[center_j*3][K][L] += Dens[I][J] * ffac1 * 0.5 * value[3];

		      F[center_j*3+1][I][J] += Dens[K][L] * ffac1 * 0.5 * value[4];
		      F[center_j*3+1][K][L] += Dens[I][J] * ffac1 * 0.5 * value[4];

		      F[center_j*3+2][I][J] += Dens[K][L] * ffac1 * 0.5 * value[5];
		      F[center_j*3+2][K][L] += Dens[I][J] * ffac1 * 0.5 * value[5];

		      F[center_k*3][I][J] += Dens[K][L] * ffac1 * 0.5 * value[6];
		      F[center_k*3][K][L] += Dens[I][J] * ffac1 * 0.5 * value[6];

		      F[center_k*3+1][I][J] += Dens[K][L] * ffac1 * 0.5 * value[7];
		      F[center_k*3+1][K][L] += Dens[I][J] * ffac1 * 0.5 * value[7];

		      F[center_k*3+2][I][J] += Dens[K][L] * ffac1 * 0.5 * value[8];
		      F[center_k*3+2][K][L] += Dens[I][J] * ffac1 * 0.5 * value[8];

		      F[center_l*3][I][J] += Dens[K][L] * ffac1 * 0.5 * value[9];
		      F[center_l*3][K][L] += Dens[I][J] * ffac1 * 0.5 * value[9];

		      F[center_l*3+1][I][J] += Dens[K][L] * ffac1 * 0.5 * value[10];
		      F[center_l*3+1][K][L] += Dens[I][J] * ffac1 * 0.5 * value[10];

		      F[center_l*3+2][I][J] += Dens[K][L] * ffac1 * 0.5 * value[11];
		      F[center_l*3+2][K][L] += Dens[I][J] * ffac1 * 0.5 * value[11];
		    }
		    else if(I==K || J==L) {
		      F[center_i*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[0];
		      F[center_i*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[0];
		      F[center_i*3][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[0];
		      F[center_i*3][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[0];

		      F[center_i*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[1];
		      F[center_i*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[1];
		      F[center_i*3+1][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[1];
		      F[center_i*3+1][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[1];

		      F[center_i*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[2];
		      F[center_i*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[2];
		      F[center_i*3+2][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[2];
		      F[center_i*3+2][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[2];

		      F[center_j*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[3];
		      F[center_j*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[3];
		      F[center_j*3][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[3];
		      F[center_j*3][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[3];

		      F[center_j*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[4];
		      F[center_j*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[4];
		      F[center_j*3+1][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[4];
		      F[center_j*3+1][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[4];

		      F[center_j*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[5];
		      F[center_j*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[5];
		      F[center_j*3+2][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[5];
		      F[center_j*3+2][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[5];

		      F[center_k*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[6];
		      F[center_k*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[6];
		      F[center_k*3][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[6];
		      F[center_k*3][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[6];

		      F[center_k*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[7];
		      F[center_k*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[7];
		      F[center_k*3+1][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[7];
		      F[center_k*3+1][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[7];

		      F[center_k*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[8];
		      F[center_k*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[8];
		      F[center_k*3+2][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[8];
		      F[center_k*3+2][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[8];

		      F[center_l*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[9];
		      F[center_l*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[9];
		      F[center_l*3][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[9];
		      F[center_l*3][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[9];

		      F[center_l*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[10];
		      F[center_l*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[10];
		      F[center_l*3+1][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[10];
		      F[center_l*3+1][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[10];

		      F[center_l*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[11];
		      F[center_l*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[11];
		      F[center_l*3+2][I][K] -= Dens[J][L] * ffac2 * 0.5 * value[11];
		      F[center_l*3+2][J][L] -= Dens[I][K] * ffac2 * 0.5 * value[11];
		    }
		    else if(J==K || I==L) {
		      F[center_i*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[0];
		      F[center_i*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[0];
		      F[center_i*3][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[0];
		      F[center_i*3][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[0];

		      F[center_i*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[1];
		      F[center_i*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[1];
		      F[center_i*3+1][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[1];
		      F[center_i*3+1][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[1];

		      F[center_i*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[2];
		      F[center_i*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[2];
		      F[center_i*3+2][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[2];
		      F[center_i*3+2][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[2];

		      F[center_j*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[3];
		      F[center_j*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[3];
		      F[center_j*3][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[3];
		      F[center_j*3][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[3];

		      F[center_j*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[4];
		      F[center_j*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[4];
		      F[center_j*3+1][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[4];
		      F[center_j*3+1][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[4];

		      F[center_j*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[5];
		      F[center_j*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[5];
		      F[center_j*3+2][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[5];
		      F[center_j*3+2][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[5];

		      F[center_k*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[6];
		      F[center_k*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[6];
		      F[center_k*3][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[6];
		      F[center_k*3][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[6];

		      F[center_k*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[7];
		      F[center_k*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[7];
		      F[center_k*3+1][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[7];
		      F[center_k*3+1][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[7];

		      F[center_k*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[8];
		      F[center_k*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[8];
		      F[center_k*3+2][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[8];
		      F[center_k*3+2][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[8];

		      F[center_l*3][I][J] += Dens[K][L] * ffac1 * 0.75 * value[9];
		      F[center_l*3][K][L] += Dens[I][J] * ffac1 * 0.75 * value[9];
		      F[center_l*3][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[9];
		      F[center_l*3][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[9];

		      F[center_l*3+1][I][J] += Dens[K][L] * ffac1 * 0.75 * value[10];
		      F[center_l*3+1][K][L] += Dens[I][J] * ffac1 * 0.75 * value[10];
		      F[center_l*3+1][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[10];
		      F[center_l*3+1][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[10];

		      F[center_l*3+2][I][J] += Dens[K][L] * ffac1 * 0.75 * value[11];
		      F[center_l*3+2][K][L] += Dens[I][J] * ffac1 * 0.75 * value[11];
		      F[center_l*3+2][I][L] -= Dens[J][K] * ffac3 * 0.5 * value[11];
		      F[center_l*3+2][J][K] -= Dens[I][L] * ffac3 * 0.5 * value[11];
		    }
		    else if(I==J || K==L) {
		      F[center_i*3][I][J] += Dens[K][L] * ffac1 * value[0];
		      F[center_i*3][K][L] += Dens[I][J] * ffac1 * value[0];
		      F[center_i*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[0];
		      F[center_i*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[0];

		      F[center_i*3+1][I][J] += Dens[K][L] * ffac1 * value[1];
		      F[center_i*3+1][K][L] += Dens[I][J] * ffac1 * value[1];
		      F[center_i*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[1];
		      F[center_i*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[1];

		      F[center_i*3+2][I][J] += Dens[K][L] * ffac1 * value[2];
		      F[center_i*3+2][K][L] += Dens[I][J] * ffac1 * value[2];
		      F[center_i*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[2];
		      F[center_i*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[2];

		      F[center_j*3][I][J] += Dens[K][L] * ffac1 * value[3];
		      F[center_j*3][K][L] += Dens[I][J] * ffac1 * value[3];
		      F[center_j*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[3];
		      F[center_j*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[3];

		      F[center_j*3+1][I][J] += Dens[K][L] * ffac1 * value[4];
		      F[center_j*3+1][K][L] += Dens[I][J] * ffac1 * value[4];
		      F[center_j*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[4];
		      F[center_j*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[4];

		      F[center_j*3+2][I][J] += Dens[K][L] * ffac1 * value[5];
		      F[center_j*3+2][K][L] += Dens[I][J] * ffac1 * value[5];
		      F[center_j*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[5];
		      F[center_j*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[5];

		      F[center_k*3][I][J] += Dens[K][L] * ffac1 * value[6];
		      F[center_k*3][K][L] += Dens[I][J] * ffac1 * value[6];
		      F[center_k*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[6];
		      F[center_k*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[6];

		      F[center_k*3+1][I][J] += Dens[K][L] * ffac1 * value[7];
		      F[center_k*3+1][K][L] += Dens[I][J] * ffac1 * value[7];
		      F[center_k*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[7];
		      F[center_k*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[7];

		      F[center_k*3+2][I][J] += Dens[K][L] * ffac1 * value[8];
		      F[center_k*3+2][K][L] += Dens[I][J] * ffac1 * value[8];
		      F[center_k*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[8];
		      F[center_k*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[8];

		      F[center_l*3][I][J] += Dens[K][L] * ffac1 * value[9];
		      F[center_l*3][K][L] += Dens[I][J] * ffac1 * value[9];
		      F[center_l*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[9];
		      F[center_l*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[9];

		      F[center_l*3+1][I][J] += Dens[K][L] * ffac1 * value[10];
		      F[center_l*3+1][K][L] += Dens[I][J] * ffac1 * value[10];
		      F[center_l*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[10];
		      F[center_l*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[10];

		      F[center_l*3+2][I][J] += Dens[K][L] * ffac1 * value[11];
		      F[center_l*3+2][K][L] += Dens[I][J] * ffac1 * value[11];
		      F[center_l*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[11];
		      F[center_l*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[11];
		    }
		    else {
		      F[center_i*3][I][J] += Dens[K][L] * ffac1 * value[0];
		      F[center_i*3][K][L] += Dens[I][J] * ffac1 * value[0];
		      F[center_i*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[0];
		      F[center_i*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[0];	
		      F[center_i*3][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[0];
		      F[center_i*3][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[0];

		      F[center_i*3+1][I][J] += Dens[K][L] * ffac1 * value[1];
		      F[center_i*3+1][K][L] += Dens[I][J] * ffac1 * value[1];
		      F[center_i*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[1];
		      F[center_i*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[1];	
		      F[center_i*3+1][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[1];
		      F[center_i*3+1][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[1];

		      F[center_i*3+2][I][J] += Dens[K][L] * ffac1 * value[2];
		      F[center_i*3+2][K][L] += Dens[I][J] * ffac1 * value[2];
		      F[center_i*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[2];
		      F[center_i*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[2];	
		      F[center_i*3+2][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[2];
		      F[center_i*3+2][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[2];

		      F[center_j*3][I][J] += Dens[K][L] * ffac1 * value[3];
		      F[center_j*3][K][L] += Dens[I][J] * ffac1 * value[3];
		      F[center_j*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[3];
		      F[center_j*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[3];	
		      F[center_j*3][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[3];
		      F[center_j*3][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[3];

		      F[center_j*3+1][I][J] += Dens[K][L] * ffac1 * value[4];
		      F[center_j*3+1][K][L] += Dens[I][J] * ffac1 * value[4];
		      F[center_j*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[4];
		      F[center_j*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[4];	
		      F[center_j*3+1][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[4];
		      F[center_j*3+1][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[4];

		      F[center_j*3+2][I][J] += Dens[K][L] * ffac1 * value[5];
		      F[center_j*3+2][K][L] += Dens[I][J] * ffac1 * value[5];
		      F[center_j*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[5];
		      F[center_j*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[5];	
		      F[center_j*3+2][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[5];
		      F[center_j*3+2][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[5];

		      F[center_k*3][I][J] += Dens[K][L] * ffac1 * value[6];
		      F[center_k*3][K][L] += Dens[I][J] * ffac1 * value[6];
		      F[center_k*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[6];
		      F[center_k*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[6];	
		      F[center_k*3][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[6];
		      F[center_k*3][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[6];

		      F[center_k*3+1][I][J] += Dens[K][L] * ffac1 * value[7];
		      F[center_k*3+1][K][L] += Dens[I][J] * ffac1 * value[7];
		      F[center_k*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[7];
		      F[center_k*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[7];	
		      F[center_k*3+1][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[7];
		      F[center_k*3+1][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[7];

		      F[center_k*3+2][I][J] += Dens[K][L] * ffac1 * value[8];
		      F[center_k*3+2][K][L] += Dens[I][J] * ffac1 * value[8];
		      F[center_k*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[8];
		      F[center_k*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[8];	
		      F[center_k*3+2][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[8];
		      F[center_k*3+2][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[8];

		      F[center_l*3][I][J] += Dens[K][L] * ffac1 * value[9];
		      F[center_l*3][K][L] += Dens[I][J] * ffac1 * value[9];
		      F[center_l*3][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[9];
		      F[center_l*3][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[9];	
		      F[center_l*3][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[9];
		      F[center_l*3][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[9];

		      F[center_l*3+1][I][J] += Dens[K][L] * ffac1 * value[10];
		      F[center_l*3+1][K][L] += Dens[I][J] * ffac1 * value[10];
		      F[center_l*3+1][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[10];
		      F[center_l*3+1][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[10];	
		      F[center_l*3+1][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[10];
		      F[center_l*3+1][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[10];

		      F[center_l*3+2][I][J] += Dens[K][L] * ffac1 * value[11];
		      F[center_l*3+2][K][L] += Dens[I][J] * ffac1 * value[11];
		      F[center_l*3+2][I][K] -= Dens[J][L] * ffac2 * 0.25 * value[11];
		      F[center_l*3+2][J][L] -= Dens[I][K] * ffac2 * 0.25 * value[11];	
		      F[center_l*3+2][I][L] -= Dens[J][K] * ffac3 * 0.25 * value[11];
		      F[center_l*3+2][J][K] -= Dens[I][L] * ffac3 * 0.25 * value[11];
		    }

		  } /* l */
		} /* k */
	      } /* j */
	    } /* i */
#endif

	    /*--------------------------------------
	      
	      Form the contributions to the Hessian
	      
	      --------------------------------------*/
	    /*--- Figure out the prefactor ---*/
	    pfac = 1.0;
	    if (si == sj)
	      pfac *= 0.5;
	    if (sk == sl)
	      pfac *= 0.5;
	    if (si == sk && sj == sl || si == sl && sj == sk)
	      pfac *= 0.5;

#if INCLUDE_2ND
	    /*-------------
	      Form FourInd
	     -------------*/
	    si_fao = BasisSet.shells[si].fao-1;
	    sj_fao = BasisSet.shells[sj].fao-1;
	    sk_fao = BasisSet.shells[sk].fao-1;
	    sl_fao = BasisSet.shells[sl].fao-1;
	    /*--- RHF or UHF case ---*/
	    count = 0;
	    for (I = si_fao; I < si_fao+ni; I++)
	      for (J = sj_fao; J < sj_fao+nj; J++)
		for (K = sk_fao; K < sk_fao+nk; K++)
		  for (L = sl_fao; L < sl_fao+nl; L++) {
		    FourInd[count] = pfac*
		      (4.0*Dens[I][J]*Dens[K][L] -
		       Dens[I][K]*Dens[J][L] -
		       Dens[I][L]*Dens[K][J])*
		      GTOs.bf_norm[orig_am[0]][I-si_fao]*
		      GTOs.bf_norm[orig_am[1]][J-sj_fao]*
		      GTOs.bf_norm[orig_am[2]][K-sk_fao]*
		      GTOs.bf_norm[orig_am[3]][L-sl_fao];
		    count++;
		  }
	    
	    center[0] = center_i;
	    center[1] = center_j;
	    center[2] = center_k;
	    center[3] = center_l;

	    /*--- d2/didj contractions, neither di nor dj equal Bxyz ---*/
	    for(di=0;di<12;di++) {
	      if (di<3 || di>5) {
		for(dj=di;dj<12;dj++) {
		  if (dj<3 || dj>5) {
		    data = Libderiv.ABCD[12+di*12+dj];
		    
		    contr = 0.0;
		    for(k=0;k<quartet_size;k++)
		      contr += data[k]*FourInd[k];
		    d2acc[di][dj] = contr;
		  }
		}
	      }
	    }
	    
	    /*--- d2/dAidBj contractions by translational invariance ---*/
	    for(di=0;di<3;di++) {
	      for(dj=3;dj<6;dj++) {
		/* d2/dAidBj = - d2/dAidCj - d2/dAidDj ... */
		contr = - d2acc[di][dj+3] - d2acc[di][dj+6];
		/*             - d2/dAidAj */
		if (di <= dj-3)
		  contr -= d2acc[di][dj-3];
		else
		  contr -= d2acc[dj-3][di];
		d2acc[di][dj] = contr;
	      }
	    }
	    
	    /*--- d2/dBidBj contractions by translational invariance ---*/
	    for(di=3;di<6;di++) {
	      for(dj=di;dj<6;dj++) {
		/* d2/dBidBj = d2/dAidAj + d2/dAidCj  + d2/dAidDj  +... */
		contr = d2acc[di-3][dj-3] + d2acc[di-3][dj+3] + d2acc[di-3][dj+6] +
		  d2acc[dj-3][di+3] + d2acc[di+3][dj+3] + d2acc[di+3][dj+6] + 
		  d2acc[dj-3][di+6] + d2acc[dj+3][di+6] + d2acc[di+6][dj+6];
		d2acc[di][dj] = contr;
	      }
	    }
	    
	    /*--- d2/dBidCDj contractions by translational invariance ---*/
	    for(di=3;di<6;di++) {
	      for(dj=6;dj<12;dj++) {
		/* d2/dBidCDj = - d2/dAidCDj ... */
		contr = - d2acc[di-3][dj];
		/*             - d2/dCidCDj */
		if (di+3 <= dj)
		  contr -= d2acc[di+3][dj];
		else
		  contr -= d2acc[dj][di+3];
		/*             - d2/dDidCDj */
		if (di+6 <= dj)
		  contr -= d2acc[di+6][dj];
		else
		  contr -= d2acc[dj][di+6];
		d2acc[di][dj] = contr;
	      }
	    }
	    
	    /*--- add contributions to the Hessian ---*/
	    for(di=0;di<12;di++) {
	      coord1 = center[di/3]*3 + di%3;
	      for(dj=di;dj<12;dj++) {
		coord2 = center[dj/3]*3 + dj%3;
		if (coord1 == coord2 && di != dj)
		  perm_pf = 2.0;
		else
		  perm_pf = 1.0;
		hess_te[coord1][coord2] += perm_pf*d2acc[di][dj];
	      }
	    }
#endif

#if PK_ORDER
	  } /* unique */
#endif
	} /* sll */


  /*--- symmetrize and print out Hessian */
  symmetrize_hessian(hess_te);
  if (UserOptions.print_lvl >= PRINT_OEDERIV)
    print_atommat("Two-electron component of the molecular Hessian (a.u.)",hess_te);

  add_mat(Hess,hess_te,Hess,Molecule.num_atoms*3,Molecule.num_atoms*3);

  /*---------
    Clean-up
    ---------*/
  free_libderiv(&Libderiv);
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt_table(&fjt_table);
  free_fjt();
#endif
  free_block(hess_te);

  return;
}
}}
