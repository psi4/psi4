/*! \file cc_bt2_thread.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<pthread.h>
#include<stdio.h>
#include<memory.h>
#include<stdlib.h>
#include<cstring>
#include<cmath>
#include<libqt/qt.h>
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
#include"Tools/int_fjt.h"
#endif

namespace psi { namespace CINTS {
  void *cc_bt2_thread(void *tnum_ptr)
  {
    const long int thread_num = (long int) tnum_ptr;
    const double toler = UserOptions.cutoff;
    const double m_sqrt1_2 = 1/sqrt(2.0);
    
    /*--- Various data structures ---*/
    struct shell_pair *sp_ij, *sp_kl;
    struct unique_shell_pair *usp_ij,*usp_kl;
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
    int i, j, k, l, m, p, q, r, s;
    int ni, nj, nk, nl;
    int p_abs, q_abs, r_abs, s_abs;
    int si, sj, sk, sl ;
    int sii, sjj, skk, sll , slll;
    int num_ij, swap_ij_kl;
    int pi, pj, pk, pl;
    int max_pj, max_pl;
    int *sj_arr, *sk_arr, *sl_arr;
    int si_fao, sj_fao, sk_fao, sl_fao;
    int usii,usjj,uskk,usll,usi,usj,usk,usl,usij;
    int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
    int *R_list, *S_list, *T_list;
    int R,S,T;
    int dcr_ij, dcr_kl, dcr_ijkl;
    int lambda_T = 1;
    int num_unique_quartets;
    int quartet_index;
    int plquartet;
    int max_num_unique_quartets;
    int max_num_prim_comb;
    int switch_ij, switch_kl, switch_ijkl;
    
    int size;
    int max_class_size;
    int max_cart_class_size;
    
    int np_i, np_j, np_k, np_l;
    int nr, ns, np, nq;
    int num_prim_comb;
    
    int num_ibatch, num_i_per_ibatch, ibatch, ibatch_length;
    int imin, imax, jmin;
    int max_bf_per_shell;
    int mo_i, mo_j, mo_a, mo_b, mo_ij;
    int ia;
    int rs_offset, rsi_offset, rsp_offset;
    
    double AB2, CD2;
    double value;
    double *raw_data;             /* pointer to the unnormalized taregt quartet of integrals */
    double *data;                 /* pointer to the transformed normalized target quartet of integrals */
#ifdef NONDOUBLE_INTS
    REALTYPE *target_ints;            /* Pointer to the location of the target quartet on the stack of
					 integrals quartets if libint.a is using other than regular doubles */
#endif
    
    double *rspq_ptr;
    double temp;
    double *mo_vec;
    double *rsiq_buf;             /* buffer for (rs|iq) integrals, where r,s run over shell sets,
				     i runs over I-batch, q runs over all AOs */
    double *rsi_row, *i_row;
    double *ia_block_ptr;
    double *rsia_buf;             /* buffer for (rs|ia) integrals, where r,s run over shell sets,
				     i runs over I-batch, q runs over all AOs */
    double *jsi_row;
    double *jbi_row;
    double iajb, ibja, pfac, k0ijab, k1ijab, eijab, e0, e1;
    
    double temp1,temp2,*iq_row,*ip_row;
    int rs,qrs;
    
    double *scratch_buf;          /* scratch used in permuting bra and ket */
    
    int nao, nocc;
    int ji, ki, il, li, jk, kj, lj, lk;
    int ii, jj, kk, ll, IJ, KL;
    
    nao = BasisSet.num_ao;
    nocc = CCInfo.nocc;
    
    /*---------------
      Initialization
      ---------------*/
#ifndef USE_TAYLOR_FM
    init_fjt_table(&fjt_table);
#endif
    
    /*-------------------------
      Allocate data structures
    -------------------------*/
    max_bf_per_shell = ioff[BasisSet.max_am];
    max_cart_class_size = (max_bf_per_shell)*
    (max_bf_per_shell)*
      (max_bf_per_shell)*
      (max_bf_per_shell);
    max_num_unique_quartets = Symmetry.max_stab_index*
      Symmetry.max_stab_index*
      Symmetry.max_stab_index;
    sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    if (Symmetry.nirreps > 1)
      max_class_size = max_cart_class_size;
    max_num_prim_comb = (BasisSet.max_num_prims*
			 BasisSet.max_num_prims)*
      (BasisSet.max_num_prims*
       BasisSet.max_num_prims);
    init_libint(&Libint, BasisSet.max_am-1, max_num_prim_comb);
    
#ifdef NONDOUBLE_INTS
    raw_data = init_array(max_cart_class_size);
#endif
    
    
    /*--------------------------------------------
      generate all symmetry unique shell quartets
      --------------------------------------------*/
    for (usii=0; usii<Symmetry.num_unique_shells; usii++)
      for (usjj=0; usjj<=usii; usjj++)
	for (uskk=0; uskk<=usii; uskk++)
	  for (usll=0; usll<= ((usii == uskk) ? usjj : uskk); usll++, quartet_index++){
	    
	    usi = usii;
	    usj = usjj;
	    usk = uskk;
	    usl = usll;
	    
	    /*--- Decide if this thread will do this ---*/
	    if ( quartet_index%UserOptions.num_threads != thread_num )
	      continue;
	    
	    si = Symmetry.us2s[usii];
	    sjj = Symmetry.us2s[usjj];
	    skk = Symmetry.us2s[uskk];
	    sll = Symmetry.us2s[usll];
	    
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
	      lambda_T = Symmetry.nirreps/Symmetry.dcr_deg[stab_ij][stab_kl];
	      ni = (BasisSet.puream ? 2*BasisSet.shells[si].am - 1 : ioff[BasisSet.shells[si].am]);
	      nj = (BasisSet.puream ? 2*BasisSet.shells[sjj].am - 1 : ioff[BasisSet.shells[sjj].am]);
	      nk = (BasisSet.puream ? 2*BasisSet.shells[skk].am - 1 : ioff[BasisSet.shells[skk].am]);
	      nl = (BasisSet.puream ? 2*BasisSet.shells[sll].am - 1 : ioff[BasisSet.shells[sll].am]);
	      
	      memset(sj_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sk_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sl_arr,0,sizeof(int)*max_num_unique_quartets);
	      count = 0;
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
		    
		    sj_arr[count] = sj;
		    sk_arr[count] = sk;
		    sl_arr[count] = sl;
		    count++;
		  }
		}
	      } /* petite list is ready to be used */
	      num_unique_quartets = count;
	    }
	    else { /*--- C1 symmetry case ---*/
	      num_unique_quartets = 1;
	      sj_arr[0] = usj;
	      sk_arr[0] = usk;
	      sl_arr[0] = usl;
	    }
	    
	    
	    /*----------------------------------
	      Compute the nonredundant quartets
	      ----------------------------------*/
	    for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	      si = Symmetry.us2s[usii];
	      sj = sj_arr[plquartet];
	      sk = sk_arr[plquartet];
	      sl = sl_arr[plquartet];
	      
	      total_am = BasisSet.shells[si].am + BasisSet.shells[sj].am + 
		BasisSet.shells[sk].am + BasisSet.shells[sl].am;
	      /* parity selection */
	      if (total_am%2 &&
		  BasisSet.shells[si].center==BasisSet.shells[sj].center &&
		  BasisSet.shells[sj].center==BasisSet.shells[sk].center &&
		  BasisSet.shells[sk].center==BasisSet.shells[sl].center)
		continue;
	      
	      switch_ij = 0;
	      switch_kl = 0;
	      switch_ijkl = 0;
	      /* place in "ascending" angular mom-
		 my simple way of optimizing PHG recursion (VRR) */
	      /* these first two are good for the HRR */
	      if(BasisSet.shells[si].am < BasisSet.shells[sj].am){
		dum = si;
		si = sj;
		sj = dum;
		switch_ij = 1;
	      }
	      if(BasisSet.shells[sk].am < BasisSet.shells[sl].am){
		dum = sk;
		sk = sl;
		sl = dum;
		switch_kl = 1;
	      }
	      /* this should be /good/ for the VRR */
	      if(BasisSet.shells[si].am + BasisSet.shells[sj].am > BasisSet.shells[sk].am + BasisSet.shells[sl].am){
		dum = si;
		si = sk;
		sk = dum;
		dum = sj;
		sj = sl;
		sl = dum;
		switch_ijkl = 1;
	      }
	      
	      np_i = BasisSet.shells[si].n_prims;
	      np_j = BasisSet.shells[sj].n_prims;
	      np_k = BasisSet.shells[sk].n_prims;
	      np_l = BasisSet.shells[sl].n_prims;
	      ni = ioff[BasisSet.shells[si].am];
	      nj = ioff[BasisSet.shells[sj].am];
	      nk = ioff[BasisSet.shells[sk].am];
	      nl = ioff[BasisSet.shells[sl].am];
	      si_fao = BasisSet.shells[si].fao-1;
	      sj_fao = BasisSet.shells[sj].fao-1;
	      sk_fao = BasisSet.shells[sk].fao-1;
	      sl_fao = BasisSet.shells[sl].fao-1;
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
	      
	      AB2 = Libint.AB[0]*Libint.AB[0]+Libint.AB[1]*Libint.AB[1]+Libint.AB[2]*Libint.AB[2];
	      CD2 = Libint.CD[0]*Libint.CD[0]+Libint.CD[1]*Libint.CD[1]+Libint.CD[2]*Libint.CD[2];
	      
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
	      
	      /*--- Compute the integrals ---*/
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
		/* No need to transforms integrals to sph. harm. basis */
		data = norm_quartet(raw_data, NULL, orig_am, 0);
	      }
	      else {
		temp = 0.0;
		for(p=0;p<num_prim_comb;p++)
		  temp += (double) Libint.PrimQuartet[p].F[0];
#ifdef NONDOUBLE_INTS
		raw_data[0] = temp;
		data = raw_data;
#else
		Libint.int_stack[0] = temp;
		data = Libint.int_stack;
#endif
	      }
	      
	      ijkl = 0;
	      for(i=0;i<ni;i++) {
		ii = i + si_fao; 
		for(j=0;j<nj;j++) {
		  jj = j + sj_fao; 
		  IJ = INDEX(ii,jj);
		  for(k=0;k<nk;k++) {
		    kk = k + sk_fao; 
		    for(l=0;l<nl;l++,ijkl++) {
		      ll = l + sl_fao; 
		      KL = INDEX(kk,ll);
		      
		      if(si == sj && ii < jj) continue;
		      if(sk == sl && kk < ll) continue;
		      if(INDEX(si,sj) == INDEX(sk,sl) && IJ < KL) continue;     
		      
		      value = data[ijkl];
		      
		      /*
			if(fabs(value) > 1e-8)
			fprintf(outfile, "%d %d %d %d %20.14f\n", ii+1, jj+1, kk+1, ll+1, value);
		      */
		      
		      ij = ii*nao + jj; ji = jj*nao + ii;
		      ik = ii*nao + kk; ki = kk*nao + ii;
		      il = ii*nao + ll; li = ll*nao + ii;
		      jk = jj*nao + kk; kj = kk*nao + jj;
		      jl = jj*nao + ll; lj = ll*nao + jj;
		      kl = kk*nao + ll; lk = ll*nao + kk;
		      
		      /* (ij|kl) */
		      C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jl], 1, CCInfo.T2_t[ik],1);
		      
		      if(ii!=jj && kk!=ll && IJ!=KL) {
			/* (ij|lk) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			/* (ji|kl) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			/* (ji|lk) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ik], 1, CCInfo.T2_t[jl],1);
			/* (kl|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			/* (kl|ji) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[li], 1, CCInfo.T2_t[kj],1);
			/* (lk|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[kj], 1, CCInfo.T2_t[li],1);
			/* (lk|ji) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ki], 1, CCInfo.T2_t[lj],1);
		      }
		      else if(ii!=jj && kk!=ll && IJ==KL) {
			/* (ij|lk) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			/* (ji|kl) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			/* (ji|lk) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ik], 1, CCInfo.T2_t[jl],1);
		      }
		      else if(ii!=jj && kk==ll) {
			/* (ji|kl) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			/* (kl|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			/* (kl|ji) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[li], 1, CCInfo.T2_t[kj],1);
		      }
		      else if(ii==jj && kk!=ll) {
			/* (ij|lk) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			/* (kl|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			/* (lk|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[kj], 1, CCInfo.T2_t[li],1);
		      }
		      else if(ii==jj && kk==ll && IJ!=KL) {
			/* (kl|ij) */
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
		      }
		    }
		  }
		}
	      }
	      
	    } /* end of RSPQ loop */
	  } /* end of "unique" RSQP loop */
    
    /*---------
      Clean-up
      ---------*/
#ifdef NONDOUBLE_INTS
    free(raw_data);
#endif
    free_libint(&Libint);
    free(sj_arr);
    free(sk_arr);
    free(sl_arr);
#ifndef USE_TAYLOR_FM
    free_fjt_table(&fjt_table);
#endif
    
    
    return NULL;
  }
};
};
