/*! \file te_deriv2_scf_symm.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory.h>
#include <psitypes.h>
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
#include "hash.h"
#include "small_fns.h"
#include <Tools/prints.h>
#include <stdexcept>

#define INCLUDE_1ST 1
#define INCLUDE_2ND 1

namespace psi { namespace CINTS {

void te_deriv2_scf_symm(void)
{
  /*--- Various data structures ---*/
  struct shell_pair *sp_ij, *sp_kl;
  Libderiv_t Libderiv;                /* Integrals library object */
  htable_t htable;                    /* hashing table */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;               /* table of auxiliary function F_m(u) for each primitive combination */
#endif
  double **hess_te;

  PSI_INT_LEAST64 quartet_index;
  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int i, j, k, l, m, ii, jj, kk, ll;
  int si, sj, sk, sl ;
  int sii, sjj, skk, sll, slll;
  int usi, usj, usk, usl;
  int usii, usjj, uskk, usll, uslll;
  int pi, pj, pk, pl ;
  int max_pj, max_pl;
  int pii, pjj, pkk, pll ;
  int switch_ij, switch_kl, switch_ijkl;
  int center_i, center_j, center_k, center_l, center[4];

  int *si_arr, *sj_arr, *sk_arr, *sl_arr, *key_arr;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int usi_eq_usj, usi_eq_usk, usi_eq_usl, usj_eq_usl, usk_eq_usj, usk_eq_usl;
  int usij_eq_uskl, usik_eq_usjl, usil_eq_uskj;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  double lambda_T = 0.5/Symmetry.nirreps;
  int num_unique_quartets;
  int max_num_unique_quartets;
  int plquartet;
  PSI_INT_LEAST64 key, key1, key2, key3;
  int new_quartet, htable_ptr, nstri;
  double q4ijkl;

  int class_size;
  int max_class_size;
  int max_cart_class_size;

  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, quartet_size;

  int num_prim_comb, p, max_num_prim_comb;
  double AB2, CD2;
  double *FourInd;
  int I, J, K, L, IJ, KL, this_quartet, upk, num_unique_pk, unique;
  int si_fao, sj_fao, sk_fao, sl_fao;
  int coord, coord1, coord2;
  int di, dj;
  double *data, value[12], svalue;
  double fac1, fac2, fac3;
  double ffac1, ffac2, ffac3;
  double d2acc[12][12], contr, pfac, perm_pf, temp;

  /*---------------
    Initialization
    ---------------*/
  if (Symmetry.nirreps > 1)
    init_htable( &htable, Symmetry.max_stab_index );
  hess_te = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+DERIV_LVL);
  init_fjt_table(&fjt_table);
#endif
  init_libderiv_base();

  max_cart_class_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  max_class_size = max_cart_class_size;
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
    (BasisSet.max_num_prims*BasisSet.max_num_prims);
  init_libderiv12(&Libderiv,BasisSet.max_am-1,max_num_prim_comb,max_cart_class_size);
  FourInd = init_array(max_cart_class_size);
  max_num_unique_quartets = Symmetry.max_stab_index*
                            Symmetry.max_stab_index*
                            Symmetry.max_stab_index;
  if (Symmetry.nirreps == 1) {
    si_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sj_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sk_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sl_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
  }
  key_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);

  /*--------------------------------------------
    generate all symmetry unique shell quartets
    with ordering for building the PK-matrix
   --------------------------------------------*/
  quartet_index = 0;
  for (usii=0; usii<Symmetry.num_unique_shells; usii++)
    for (usjj=0; usjj<=usii; usjj++)
      for (uskk=0; uskk<=usjj; uskk++)
	for (usll=0; usll<=uskk; usll++, quartet_index++){

/*--- Decide what unique shell quartets out of (ij|kl), (ik|jl), and (il|jk) are unique
        with respect to permutations ---*/
	  usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
	  if (usii == usjj && usii == uskk || usjj == uskk && usjj == usll)
	    num_unique_pk = 1;
	  else if (usii == uskk || usjj == usll) {
	    num_unique_pk = 2;
	    usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
	  }
	  else if (usjj == uskk) {
	    num_unique_pk = 2;
	    usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
	  }
	  else if (usii == usjj || uskk == usll) {
	    num_unique_pk = 2;
	    usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
	  }
	  else {
	    num_unique_pk = 3;
	    usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
	    usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
	  }

	  count = 0;
	  for(upk=0;upk<num_unique_pk;upk++) {
	    /*--- For each combination of unique shells generate "petit list" of shells ---*/
	    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

	    si = Symmetry.us2s[usi];
	    sjj = Symmetry.us2s[usj];
	    skk = Symmetry.us2s[usk];
	    sll = Symmetry.us2s[usl];
	    total_am = BasisSet.shells[si].am+BasisSet.shells[sjj].am+BasisSet.shells[skk].am+BasisSet.shells[sll].am;

	    if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
	      /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
	      stab_i = Symmetry.atom_positions[BasisSet.shells[si].center-1];
	      stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
	      stab_k = Symmetry.atom_positions[BasisSet.shells[skk].center-1];
	      stab_l = Symmetry.atom_positions[BasisSet.shells[sll].center-1];
	      stab_ij = Symmetry.GnG[stab_i][stab_j];
	      stab_kl = Symmetry.GnG[stab_k][stab_l];
	      R_list = Symmetry.dcr[stab_i][stab_j];
	      S_list = Symmetry.dcr[stab_k][stab_l];
	      T_list = Symmetry.dcr[stab_ij][stab_kl];
	      fac1 = Symmetry.nirreps/
		     Symmetry.dcr_deg[Symmetry.GnG[stab_i][stab_j]][Symmetry.GnG[stab_k][stab_l]];
	      usi_eq_usj = (usi == usj);
	      usi_eq_usk = (usi == usk);
	      usi_eq_usl = (usi == usl);
	      usj_eq_usl = (usj == usl);
	      usk_eq_usj = (usk == usj);
	      usk_eq_usl = (usk == usl);
	      usij_eq_uskl = (INDEX(usi,usj) == INDEX(usk,usl));
	      usik_eq_usjl = (INDEX(usi,usk) == INDEX(usj,usl));
	      usil_eq_uskj = (INDEX(usi,usl) == INDEX(usk,usj));
	      if (!usi_eq_usj)
		fac1 *= 2.0;
	      if (!usk_eq_usl)
		fac1 *= 2.0;
	      if (usij_eq_uskl)
		fac1 *= 0.5;
		

	      for(dcr_ij=0;dcr_ij<Symmetry.dcr_dim[stab_i][stab_j];dcr_ij++){
		R = R_list[dcr_ij];
		sj = BasisSet.shells[sjj].trans_vec[R]-1;
		for(dcr_ijkl=0;dcr_ijkl<Symmetry.dcr_dim[stab_ij][stab_kl];dcr_ijkl++){
		  T = T_list[dcr_ijkl];
		  sk = BasisSet.shells[skk].trans_vec[T]-1;
		  slll = BasisSet.shells[sll].trans_vec[T]-1;
		  for(dcr_kl=0;dcr_kl<Symmetry.dcr_dim[stab_k][stab_l];dcr_kl++) {
		    S = S_list[dcr_kl];
		    sl = BasisSet.shells[slll].trans_vec[S]-1;

		    q4ijkl = fac1;
		    
		    key1 = compute_key(si,sj,sk,sl);
		    key2 = compute_key(si,sk,sj,sl);
		    key3 = compute_key(si,sl,sk,sj);
		    new_quartet = put_entry(&htable,key1,si,sj,sk,sl,q4ijkl,0,0);
		    if (new_quartet > -1) {
		      key_arr[count] = new_quartet;
		      count++;
		    }
		    
		    if ( (key1 == key3 && key3 != key2) ||
			 (key2 == key3 && key3 != key1) ||
			 (key2 != key3 && key1 != key3 && key2 != key1)) {
		      new_quartet = put_entry(&htable,key2,si,sk,sj,sl,0,q4ijkl,0);
		      if (new_quartet > -1) {
			key_arr[count] = new_quartet;
			count++;
		      }
		    }
		    
		    if ( (key1 == key2 && key3 != key1) ||
			 (key2 != key3 && key1 != key3 && key1 != key2)) {
		      new_quartet = put_entry(&htable,key3,si,sl,sk,sj,0,0,q4ijkl);
		      if (new_quartet > -1) {
			key_arr[count] = new_quartet;
			count++;
		      }
		    }
		  }
		}
	      } /* petite list is ready to be used */
	    }
	    else {
	      si_arr[count] = si;
	      sj_arr[count] = sjj;
	      sk_arr[count] = skk;
	      sl_arr[count] = sll;
	      count++;
	    }
	  }
	  num_unique_quartets = count;
	  if (count > 3*max_num_unique_quartets)
	    throw std::domain_error("Problem with hashing?");

	    /*----------------------------------
	      Compute the nonredundant quartets
	     ----------------------------------*/
	    for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	      if (Symmetry.nirreps == 1) {
		si = si_arr[plquartet];
		sj = sj_arr[plquartet];
		sk = sk_arr[plquartet];
		sl = sl_arr[plquartet];
		
		fac1 = fac2 = fac3 = 1.0;
		if (si != sj)
		  fac1 *= 2.0;
		if (sk != sl)
		  fac1 *= 2.0;
		if (si != sk)
		  fac2 *= 2.0;
		if (sj != sl)
		  fac2 *= 2.0;
		if (si != sl)
		  fac3 *= 2.0;
		if (sj != sk)
		  fac3 *= 2.0;
		if (INDEX(si,sj) == INDEX(sk,sl))
		  fac1 *= 0.5;
		if (INDEX(si,sk) == INDEX(sj,sl))
		  fac2 *= 0.5;
		if (INDEX(si,sl) == INDEX(sj,sk))
		  fac3 *= 0.5;
	      }
	      else {
		htable_ptr = key_arr[plquartet];
		if (htable_ptr >= htable.size || htable_ptr < 0)
		    throw std::domain_error("Problem with hashing?");
		htable.table[htable_ptr].key = EMPTY_KEY;
		si = htable.table[htable_ptr].si;
		sj = htable.table[htable_ptr].sj;
		sk = htable.table[htable_ptr].sk;
		sl = htable.table[htable_ptr].sl;

		fac1 = htable.table[htable_ptr].q4ijkl;
		fac2 = htable.table[htable_ptr].q4ikjl;
		fac3 = htable.table[htable_ptr].q4ilkj;
		if (si == sj && si == sk ||
		    sj == sk && sj == sl ||
		    si == sj && si == sl ||
		    si == sk && si == sl) {
		  fac2 = fac3 = fac1;
		}
		else if (si == sk || sj == sl) {
		  fac1 = fac3 = (fac1 + fac3);
		}
		else if (sj == sk || si == sl) {
		  fac1 = fac2 = (fac1 + fac2);
		}
		else if (si == sj || sk == sl) {
		  fac3 = fac2 = (fac3 + fac2);
		}
	      }

#if DEBUG
	      if (si < 0 || si >= BasisSet.num_shells)
		  throw std::domain_error("Problem with shell indices");
	      if (sj < 0 || sj >= BasisSet.num_shells)
		  throw std::domain_error("Problem with shell indices");
	      if (sk < 0 || sk >= BasisSet.num_shells)
		  throw std::domain_error("Problem with shell indices");
	      if (sl < 0 || sl >= BasisSet.num_shells)
		  throw std::domain_error("Problem with shell indices");
#endif
	      
	      /* place in "ascending" angular mom-
	       my simple way of optimizing PHG recursion (VRR) */
	      /* these first two are good for the HRR */
	      if(BasisSet.shells[si].am < BasisSet.shells[sj].am){
		dum = si;
		si = sj;
		sj = dum;
		temp = fac2;
		fac2 = fac3;
		fac3 = temp;
	      }
	      if(BasisSet.shells[sk].am < BasisSet.shells[sl].am){
		dum = sk;
		sk = sl;
		sl = dum;
		temp = fac2;
		fac2 = fac3;
		fac3 = temp;
	      }
	      /* this should be /good/ for the VRR */
	      if(BasisSet.shells[si].am + BasisSet.shells[sj].am >
		 BasisSet.shells[sk].am + BasisSet.shells[sl].am){
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
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
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

		    ffac1 = fac1;
		    ffac2 = fac2;
		    ffac3 = fac3;

		    /*		    if(I!=J) ffac1 *= 2.0;
		    if(K!=L) ffac1 *= 2.0;
		    if(I!=K) ffac2 *= 2.0;
		    if(J!=L) ffac2 *= 2.0;
		    if(I!=L) ffac3 *= 2.0;
		    if(J!=K) ffac3 *= 2.0;
		    if(INDEX(I,J) != INDEX(K,L)) ffac1 *= 2.0;
		    if(INDEX(I,K) != INDEX(J,L)) ffac2 *= 2.0;
		    if(INDEX(I,L) != INDEX(J,K)) ffac3 *= 2.0;*/
		    if (I != J && si == sj)
		      ffac1 *= 2.0;
		    if (K != L && sk == sl)
		      ffac1 *= 2.0;
		    if (I != K && si == sk)
		      ffac2 *= 2.0;
		    if (J != L && sj == sl)
		      ffac2 *= 2.0;
		    if (I != L && si == sl)
		      ffac3 *= 2.0;
		    if (J != K && sj == sk)
		      ffac3 *= 2.0;
		    if (INDEX(I,J) != INDEX(K,L) &&
			INDEX(si,sj) == INDEX(sk,sl))
		      ffac1 *= 2.0;
		    if (INDEX(I,K) != INDEX(J,L) &&
			INDEX(si,sk) == INDEX(sj,sl))
		      ffac2 *= 2.0;
		    if (INDEX(I,L) != INDEX(J,K) &&
			INDEX(si,sl) == INDEX(sj,sk))
		      ffac3 *= 2.0;


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
	    if (usi == usj)
	      pfac *= 0.5;
	    if (usk == usl)
	      pfac *= 0.5;
	    if (usi == usk && usj == usl || usi == usl && usj == usk)
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
	    } /* end of loop over symmetry unique quartets */
	} /* end of loop over unique shells */

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
};};
