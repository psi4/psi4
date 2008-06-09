/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<memory.h>
#include<pthread.h>
#include<psitypes.h>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"quartet_data.h"      /* From Default_Ints */
#include"norm_quartet.h"
#include"hash.h"
#include"transmat.h"
#include"read_scf_opdm.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
  #include"fjt.h"
#endif
#include"schwartz.h"
#include"shell_block_matrix.h"

using namespace std;
      
extern pthread_mutex_t fock_mutex;    /* Used to lock Fock matrices during update */


namespace psi {
  namespace CINTS {

    extern double ****Gskel, ****Gskel_o;   /* Global skeleton Fock matrices, updated in critical sections */
void *hf_fock_thread(void *tnum_ptr)
{
  const long int thread_num = (long int) tnum_ptr;
  /*--- Various data structures ---*/
  struct tebuf *tot_data;             /* buffer for non-zero integrals */
  struct shell_pair *sp_ij, *sp_kl, *sp_ik, *sp_il, *sp_jk, *sp_jl;
  Libint_t Libint;                    /* Integrals library object */
  htable_t htable;                    /* hashing table */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;           /* table of auxiliary function F_m(u) for each primitive combination */
#endif
  
  int total_te_count = 0;
  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int g, i, j, k, l, m, ii, jj, kk, ll;
  int a, b, c, d;
  int si;                                  /* GCC compiler screwes up if static is left out */
  int sj, sk, sl, si_g, sj_g;
  int sii, sjj, skk, sll , slll;
  int sij, skl, sijkl;
  int pi, pj, pk, pl ;
  int max_pj, max_pl;
  int pii, pjj, pkk, pll;
  int upk, num_unique_pk;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int *si_arr, *sj_arr, *sk_arr, *sl_arr, *key_arr;
  int usii,usjj,uskk,usll,usi,usj,usk,usl;
  int usi_eq_usj, usi_eq_usk, usi_eq_usl, usj_eq_usl, usk_eq_usj, usk_eq_usl;
  int usij_eq_uskl, usik_eq_usjl, usil_eq_uskj;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, li, lj;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;
  PSI_INT_LEAST64 key, key1, key2, key3;
  int new_quartet, htable_ptr, nstri;
  PSI_INT_LEAST64 quartet_index;

  double so_int;
  double lambda_T = 0.5/Symmetry.nirreps;
  double AB2, CD2;
  double *data;
  double pkblock_end_value = 0.0;
  double temp;
  double **tmpmat1;
  double *qijkl_arr, *qikjl_arr, *qiljk_arr;
  double q4ijkl;
  double fac1, fac2, fac3;
  double ffac1, ffac2, ffac3;
  double c1, c2, c3, c4;
  double dmax;
  double ****G, ****G_o;       /* Shell-blocked skeleton G matrices from this thread */

  /*---
    init hashing table to store and retrieve quartet data
    init table for Fj(T)
    ---*/
  if (Symmetry.nirreps > 1)
    init_htable( &htable, Symmetry.max_stab_index );

#ifndef USE_TAYLOR_FM
  init_fjt_table(&fjt_table);
#endif

  max_cart_class_size = (ioff[BasisSet.max_am])*
			(ioff[BasisSet.max_am])*
			(ioff[BasisSet.max_am])*
			(ioff[BasisSet.max_am]);
  max_num_unique_quartets = Symmetry.max_stab_index*
			    Symmetry.max_stab_index*
			    Symmetry.max_stab_index;
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
		      (BasisSet.max_num_prims*BasisSet.max_num_prims);
  /*--- init a LIBINT object ---*/
  pthread_mutex_lock(&fock_mutex);
  UserOptions.memory -= init_libint(&Libint,BasisSet.max_am-1,max_num_prim_comb);
  pthread_mutex_unlock(&fock_mutex);

  /*--- Allocate this thread's shell-blocked skeleton G's ---*/
  G = init_shell_block_matrix();
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    G_o = init_shell_block_matrix();
  
  tot_data = (struct tebuf*) malloc(max_cart_class_size*sizeof(struct tebuf));
  if (Symmetry.nirreps == 1) {
    si_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sj_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sk_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
    sl_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);
  }
  key_arr = (int *)malloc(sizeof(int)*3*max_num_unique_quartets);

  c1 = 1.0 - 0.5*UserOptions.hf_exch;
  c2 = 1.0 - 0.25*UserOptions.hf_exch;
  c3 = -0.5*UserOptions.hf_exch;
  c4 = -0.25*UserOptions.hf_exch;
  
/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  quartet_index = 0;
  for (usii=0; usii<Symmetry.num_unique_shells; usii++)
    for (usjj=0; usjj<=usii; usjj++)
      for (uskk=0; uskk<=usjj; uskk++)
	for (usll=0; usll<=uskk; usll++, quartet_index++){

	  /*--- Decide if this thread will do this ---*/
	  if ( quartet_index%UserOptions.num_threads != thread_num )
	    continue;
	    
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

		    /*----------------------------------------------------------
		      Obviously redundant or zero cases are be eliminated here!
		      Redundancies arise in DCD approach when usk == usl etc.
		     -------------------------------------------------------------*/
		      
		    if(!(total_am%2)||
		       (BasisSet.shells[si].center!=BasisSet.shells[sj].center)||
		       (BasisSet.shells[sj].center!=BasisSet.shells[sk].center)||
		       (BasisSet.shells[sk].center!=BasisSet.shells[sl].center)) {

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
		}
	      } /* petite list is ready to be used */
	    }
	    else {
	      if(!(total_am%2)||
		 (BasisSet.shells[si].center!=BasisSet.shells[sjj].center)||
		 (BasisSet.shells[sjj].center!=BasisSet.shells[skk].center)||
		 (BasisSet.shells[skk].center!=BasisSet.shells[sll].center)) {
		si_arr[count] = si;
		sj_arr[count] = sjj;
		sk_arr[count] = skk;
		sl_arr[count] = sll;
		count++;
	      }
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
		dmax = MAX(BasisSet.shell_pairs[si][sj].Dmax, BasisSet.shell_pairs[sk][sl].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[si][sk].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[si][sl].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[sj][sk].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[sj][sl].Dmax);
		if (BasisSet.schwartz_eri[si][sj]*BasisSet.schwartz_eri[sk][sl]*dmax < UserOptions.cutoff)
		  continue;
		
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
		dmax = MAX(BasisSet.shell_pairs[si][sj].Dmax, BasisSet.shell_pairs[sk][sl].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[si][sk].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[si][sl].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[sj][sk].Dmax);
		dmax = MAX(dmax, 0.25*BasisSet.shell_pairs[sj][sl].Dmax);
		if (BasisSet.schwartz_eri[si][sj]*BasisSet.schwartz_eri[sk][sl]*dmax < UserOptions.cutoff)
		  continue;

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
	      ioffset = BasisSet.shells[si].fao - 1;
              joffset = BasisSet.shells[sj].fao - 1;
              koffset = BasisSet.shells[sk].fao - 1;
              loffset = BasisSet.shells[sl].fao - 1;
	      ni = ioff[BasisSet.shells[si].am];
	      nj = ioff[BasisSet.shells[sj].am];
	      nk = ioff[BasisSet.shells[sk].am];
	      nl = ioff[BasisSet.shells[sl].am];
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
	      sp_ik = &(BasisSet.shell_pairs[si][sk]);
	      sp_il = &(BasisSet.shell_pairs[si][sl]);
	      sp_jk = &(BasisSet.shell_pairs[sj][sk]);
	      sp_jl = &(BasisSet.shell_pairs[sj][sl]);
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
				   sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
#else
		      quartet_data(&(Libint.PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
				   sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
#endif
		    }
		  }
		}
	      }

	      /*-----------------------------------------------------
		Compute the quartet and compute its contributions to
		appropriate blocks of Gsh and Gsh_o
	       -----------------------------------------------------*/
	      if (am) {
		data = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint,num_prim_comb);
		/*                    zero here means no transformation to puream basis */
		/*                                    |                                 */
		data = norm_quartet(data, NULL, orig_am, 0);

		/*--- Here just put non-redundant integrals to tot_data ---*/
		num = 0;
		if(si==sj && sk==sl && si==sk) {
		  iimax = ni - 1;
		  for(ii=0; ii <= iimax; ii++){
		    jjmax = ii;
		    for(jj=0; jj <= jjmax; jj++){
		      kkmax = ii;
		      for(kk=0; kk <= kkmax; kk++){
			llmax = (kk==ii)? jj : kk ;
			for(ll=0; ll <= llmax; ll++){
			  ijkl = ll+nl*(kk+nk*(jj+nj*ii));
/*			  if (fabs(data[ijkl])>UserOptions.cutoff) {*/
			    tot_data[num].i = (short int) ii;
			    tot_data[num].j = (short int) jj;
			    tot_data[num].k = (short int) kk;
			    tot_data[num].l = (short int) ll;
			    tot_data[num].val = data[ijkl];
			    num++;
/*			  }*/
			}
		      }
		    }
		  }
		}
		else if(si==sk && sj==sl){
		  iimax = ni - 1;
		  for(ii=0; ii <= iimax; ii++){
		    jjmax = nj - 1;
		    for(jj=0; jj <= jjmax; jj++){
		      kkmax = ii;
		      for(kk=0; kk <= kkmax; kk++){
			llmax = (kk==ii)? jj : nl - 1;
			for(ll=0; ll <= llmax; ll++){
			  ijkl = ll+nl*(kk+nk*(jj+nj*ii));
/*			  if(fabs(data[ijkl])>UserOptions.cutoff){*/
			    tot_data[num].i = (short int) ii;
			    tot_data[num].j = (short int) jj;
			    tot_data[num].k = (short int) kk;
			    tot_data[num].l = (short int) ll;
			    tot_data[num].val = data[ijkl];
			    num++;
/*			  }*/
			}
		      }
		    }
		  }
		}
		else {
		  iimax = ni - 1;
		  kkmax = nk - 1;
		  for(ii=0; ii <= iimax; ii++){
		    jjmax = (si == sj) ? ii : nj - 1;
		    for(jj=0; jj <= jjmax; jj++){
		      for(kk=0; kk <= kkmax; kk++){
			llmax = (sk == sl) ? kk : nl - 1;
			for(ll=0; ll <= llmax; ll++){
			  ijkl = ll+nl*(kk+nk*(jj+nj*ii));
/*			  if(fabs(data[ijkl])>UserOptions.cutoff){*/
			    tot_data[num].i = (short int) ii;
			    tot_data[num].j = (short int) jj;
			    tot_data[num].k = (short int) kk;
			    tot_data[num].l = (short int) ll;
			    tot_data[num].val = data[ijkl];
			    num++;
/*			  }*/
			}
		      }
		    }
		  }
		}

	      }
	      else {
		temp = 0.0;
		for(p=0;p<num_prim_comb;p++)
		  temp += Libint.PrimQuartet[p].F[0];

		tot_data[0].i = tot_data[0].j = tot_data[0].k = tot_data[0].l = 0;
		tot_data[0].val = temp;
		num = 0;
/*		if (fabs(temp) > UserOptions.cutoff)*/
		  num = 1;

	      }

	      total_te_count += num;

	      if (UserOptions.reftype == rhf) {
		for(n=0;n<num;n++) {
		  i = tot_data[n].i;
		  j = tot_data[n].j;
		  k = tot_data[n].k;
		  l = tot_data[n].l;
		  ii = i + ioffset;
		  jj = j + joffset;
		  kk = k + koffset;
		  ll = l + loffset;
		  temp = tot_data[n].val;
		  
		  ffac1 = fac1*temp;
		  ffac2 = fac2*temp;
		  ffac3 = fac3*temp;
		  if (ii != jj && si == sj)
		      ffac1 *= 2.0;
		  if (kk != ll && sk == sl)
		      ffac1 *= 2.0;
		  if (ii != kk && si == sk)
		      ffac2 *= 2.0;
		  if (jj != ll && sj == sl)
		      ffac2 *= 2.0;
		  if (ii != ll && si == sl)
		      ffac3 *= 2.0;
		  if (jj != kk && sj == sk)
		      ffac3 *= 2.0;
		  if (INDEX(ii,jj) != INDEX(kk,ll) &&
		      INDEX(si,sj) == INDEX(sk,sl))
		      ffac1 *= 2.0;
		  if (INDEX(ii,kk) != INDEX(jj,ll) &&
		      INDEX(si,sk) == INDEX(sj,sl))
		      ffac2 *= 2.0;
		  if (INDEX(ii,ll) != INDEX(jj,kk) &&
		      INDEX(si,sl) == INDEX(sj,sk))
		      ffac3 *= 2.0;
		  
		  if (ii == jj && ii == kk ||
		      jj == kk && jj == ll ||
		      ii == jj && ii == ll ||
		      ii == kk && ii == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c1*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c1*ffac1);
		  }
		  else if (ii == kk || jj == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c2*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c2*ffac1);
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c3*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c3*ffac2);
		  }
		  else if (jj == kk || ii == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c2*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c2*ffac1);
		      G[si][sl][i][l] += sp_jk->dmat[j][k]*(c3*ffac3);
		      G[sj][sk][j][k] += sp_il->dmat[i][l]*(c3*ffac3);
		  }
		  else if (ii == jj || kk == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(ffac1);
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c4*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c4*ffac2);
		  }
		  else {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(ffac1);
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c4*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c4*ffac2);
		      G[si][sl][i][l] += sp_jk->dmat[j][k]*(c4*ffac3);
		      G[sj][sk][j][k] += sp_il->dmat[i][l]*(c4*ffac3);
		  }
		}
	      }
	      else if (UserOptions.reftype == uhf || UserOptions.reftype == rohf) {
		for(n=0;n<num;n++) {
		  i = tot_data[n].i;
		  j = tot_data[n].j;
		  k = tot_data[n].k;
		  l = tot_data[n].l;
		  ii = i + ioffset;
		  jj = j + joffset;
		  kk = k + koffset;
		  ll = l + loffset;
		  temp = tot_data[n].val;
		  
		  ffac1 = fac1*temp;
		  ffac2 = fac2*temp;
		  ffac3 = fac3*temp;
		  if (ii != jj && si == sj)
		      ffac1 *= 2.0;
		  if (kk != ll && sk == sl)
		      ffac1 *= 2.0;
		  if (ii != kk && si == sk)
		      ffac2 *= 2.0;
		  if (jj != ll && sj == sl)
		      ffac2 *= 2.0;
		  if (ii != ll && si == sl)
		      ffac3 *= 2.0;
		  if (jj != kk && sj == sk)
		      ffac3 *= 2.0;
		  if (INDEX(ii,jj) != INDEX(kk,ll) &&
		      INDEX(si,sj) == INDEX(sk,sl))
		      ffac1 *= 2.0;
		  if (INDEX(ii,kk) != INDEX(jj,ll) &&
		      INDEX(si,sk) == INDEX(sj,sl))
		      ffac2 *= 2.0;
		  if (INDEX(ii,ll) != INDEX(jj,kk) &&
		      INDEX(si,sl) == INDEX(sj,sk))
		      ffac3 *= 2.0;
		  
		  if (ii == jj && ii == kk ||
		      jj == kk && jj == ll ||
		      ii == jj && ii == ll ||
		      ii == kk && ii == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c1*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c1*ffac1);
		      G_o[si][sj][i][j] -= sp_kl->dmato[k][l]*(c3*ffac1);
		      G_o[sk][sl][k][l] -= sp_ij->dmato[i][j]*(c3*ffac1);
		  }
		  else if (ii == kk || jj == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c2*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c2*ffac1);
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c3*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c3*ffac2);
		      G_o[si][sj][i][j] -= sp_kl->dmato[k][l]*(c4*ffac1);
		      G_o[sk][sl][k][l] -= sp_ij->dmato[i][j]*(c4*ffac1);
		      G_o[si][sk][i][k] -= sp_jl->dmato[j][l]*(c3*ffac2);
		      G_o[sj][sl][j][l] -= sp_ik->dmato[i][k]*(c3*ffac2);
		  }
		  else if (jj == kk || ii == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*(c2*ffac1);
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*(c2*ffac1);
		      G[si][sl][i][l] += sp_jk->dmat[j][k]*(c3*ffac3);
		      G[sj][sk][j][k] += sp_il->dmat[i][l]*(c3*ffac3);
		      G_o[si][sj][i][j] -= sp_kl->dmato[k][l]*(c4*ffac1);
		      G_o[sk][sl][k][l] -= sp_ij->dmato[i][j]*(c4*ffac1);
		      G_o[si][sl][i][l] -= sp_jk->dmato[j][k]*(c3*ffac3);
		      G_o[sj][sk][j][k] -= sp_il->dmato[i][l]*(c3*ffac3);
		  }
		  else if (ii == jj || kk == ll) {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*ffac1;
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*ffac1;
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c4*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c4*ffac2);
		      G_o[si][sk][i][k] -= sp_jl->dmato[j][l]*(c4*ffac2);
		      G_o[sj][sl][j][l] -= sp_ik->dmato[i][k]*(c4*ffac2);
		  }
		  else {
		      G[si][sj][i][j] += sp_kl->dmat[k][l]*ffac1;
		      G[sk][sl][k][l] += sp_ij->dmat[i][j]*ffac1;
		      G[si][sk][i][k] += sp_jl->dmat[j][l]*(c4*ffac2);
		      G[sj][sl][j][l] += sp_ik->dmat[i][k]*(c4*ffac2);
		      G[si][sl][i][l] += sp_jk->dmat[j][k]*(c4*ffac3);
		      G[sj][sk][j][k] += sp_il->dmat[i][l]*(c4*ffac3);
		      G_o[si][sk][i][k] -= sp_jl->dmato[j][l]*(c4*ffac2);
		      G_o[sj][sl][j][l] -= sp_ik->dmato[i][k]*(c4*ffac2);
		      G_o[si][sl][i][l] -= sp_jk->dmato[j][k]*(c4*ffac3);
		      G_o[sj][sk][j][k] -= sp_il->dmato[i][l]*(c4*ffac3);
		  }
		}
	      }
	      
	    } /* end of computing "petit" list */
	} /* end getting unique shell combination */

  pthread_mutex_lock(&fock_mutex);
  {
      for(si=0;si<BasisSet.num_shells;si++) {
	ni = ioff[BasisSet.shells[si].am];
	for(sj=0;sj<BasisSet.num_shells;sj++) {
	  nj = ioff[BasisSet.shells[sj].am];
	  for(i=0;i<ni;i++)
	    for(j=0;j<nj;j++)
	      Gskel[si][sj][i][j] += G[si][sj][i][j];
	}
      }
      if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
	for(si=0;si<BasisSet.num_shells;si++) {
	  ni = ioff[BasisSet.shells[si].am];
	  for(sj=0;sj<BasisSet.num_shells;sj++) {
	    nj = ioff[BasisSet.shells[sj].am];
	    for(i=0;i<ni;i++)
	      for(j=0;j<nj;j++)
		Gskel_o[si][sj][i][j] += G_o[si][sj][i][j];
	  }
	}
  }
   pthread_mutex_unlock(&fock_mutex);
  

  /*---------
    Clean-up
   ---------*/
  free_shell_block_matrix(G);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    free_shell_block_matrix(G_o);
  }
  if (Symmetry.nirreps > 1) {
    free_htable(&htable);
  }
  free_libint(&Libint);
  free(tot_data);
  if (Symmetry.nirreps == 1) {
    free(si_arr);
    free(sj_arr);
    free(sk_arr);
    free(sl_arr);
  }
  free(key_arr);

#ifndef USE_TAYLOR_FM
  free_fjt_table(&fjt_table);
#endif

  return NULL;
}

};};
