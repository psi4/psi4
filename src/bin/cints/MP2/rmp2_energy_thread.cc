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
#include<libciomr/libciomr.h>
#include<libqt/qt.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"schwartz.h"
#include"quartet_data.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
  #include"fjt.h"
#endif
#include"quartet_permutations.h"
#include"rmp2_energy.h"

#define SWAP1 1

namespace psi { 
  namespace CINTS {

void *rmp2_energy_thread(void *tnum_ptr)
{
  const long int thread_num = (long int) tnum_ptr;
  const double toler = UserOptions.cutoff;
  const double m_sqrt1_2 = 1/sqrt(2.0);

  extern RMP2_Status_t RMP2_Status;
  extern pthread_mutex_t rmp2_energy_mutex;
  extern pthread_mutex_t *rmp2_sindex_mutex;
  extern pthread_cond_t rmp2_energy_cond;
  extern double *jsia_buf;             /* buffer for (js|ia) integrals, where j runs over all d.-o. MOs,
				   s runs over all AOs, i - over I-batch, a - over all virtuals */
  extern double *jbia_buf;             /* buffer contains all MP2-type integrals */

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
  int i, m, p, q, r, s;
  int p_abs, q_abs, r_abs, s_abs;
  int si, sj, sk, sl ;
  int sii, sjj, skk, sll , slll;
  int num_ij, swap_ij_kl;
  int pi, pj, pk, pl;
  int max_pj, max_pl;
  int *sj_arr, *sk_arr, *sl_arr;
  int sr_fao, ss_fao, sp_fao, sq_fao;
  int usii,usjj,uskk,usll,usi,usj,usk,usl,usij;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

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
  pthread_mutex_lock(&rmp2_energy_mutex);
  init_libint(&Libint, BasisSet.max_am-1, max_num_prim_comb);
  pthread_mutex_unlock(&rmp2_energy_mutex);

#ifdef NONDOUBLE_INTS
  raw_data = init_array(max_cart_class_size);
#endif

  num_ibatch = RMP2_Status.num_ibatch;
  num_i_per_ibatch = RMP2_Status.num_i_per_ibatch;
  rsiq_buf = init_array(num_i_per_ibatch*BasisSet.num_ao*
			max_bf_per_shell*max_bf_per_shell);
  rsia_buf = init_array(num_i_per_ibatch*MOInfo.nactuocc*
			  max_bf_per_shell*max_bf_per_shell);
  scratch_buf = init_array(MAX(max_cart_class_size,
			       num_i_per_ibatch*BasisSet.num_ao*
			       max_bf_per_shell*max_bf_per_shell));
  
/*-----------------------------------
  generate all unique shell quartets
 -----------------------------------*/
  /*--- I-batch loop ---*/
  for (ibatch=0;ibatch<num_ibatch;ibatch++) {
    imin = ibatch * num_i_per_ibatch + MOInfo.nfrdocc;
    imax = MIN( imin+num_i_per_ibatch , MOInfo.ndocc );
    ibatch_length = imax - imin;
    jmin = MOInfo.nfrdocc;
    if (thread_num == 0)
      fprintf(outfile,"  Pass #%d, MO %d through MO %d\n",ibatch,imin+1,imax);
    fflush(outfile);

    /*--- "unique" R,S loop ---*/
    usij = 0;
    for (usii=0; usii<Symmetry.num_unique_shells; usii++)
      for (usjj=0; usjj<=usii; usjj++, usij++) {
	/*--- Decide if this thread will do this ---*/
	if ( usij%UserOptions.num_threads != thread_num )
	  continue;
	usi = usii; usj = usjj;
	/*--- As usual, swap order usi and usj according to their angular momenta ---*/
	if(BasisSet.shells[Symmetry.us2s[usi]].am < BasisSet.shells[Symmetry.us2s[usj]].am){
	  dum = usi;
	  usi = usj;
	  usj = dum;
	}
	
	sii = Symmetry.us2s[usi];
	sjj = Symmetry.us2s[usj];
	if (Symmetry.nirreps > 1) {
	  stab_i = Symmetry.atom_positions[BasisSet.shells[sii].center-1];
	  stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
	  stab_ij = Symmetry.GnG[stab_i][stab_j];
	  R_list = Symmetry.dcr[stab_i][stab_j];
	  num_ij = Symmetry.dcr_dim[stab_i][stab_j];
	}
	else
	  num_ij = 1;

	/*--- R,S loop ---*/
	for(dcr_ij=0;dcr_ij<num_ij;dcr_ij++) {
	  if (Symmetry.nirreps > 1)
	    R = R_list[dcr_ij];
	  else
	    R = 0;
	  si = sii;
	  sj = BasisSet.shells[sjj].trans_vec[R]-1;

	  /*--- "Unique" P,Q loop ---*/
	  for (uskk=0; uskk<Symmetry.num_unique_shells; uskk++)
	    for (usll=0; usll<=uskk; usll++){

	      /*--- For each combination of unique shells generate "petit list" of shells ---*/
	      usk = uskk; usl = usll;
	      /*--- As usual, swap order usk and usl according to their angular momenta ---*/
	      if(BasisSet.shells[Symmetry.us2s[usk]].am < BasisSet.shells[Symmetry.us2s[usl]].am){
		dum = usk;
		usk = usl;
		usl = dum;
	      }
	      /*--- DO NOT SWAP bra and ket at this time. Do it later, in the main loop ---*/
	      if(BasisSet.shells[Symmetry.us2s[usi]].am + BasisSet.shells[Symmetry.us2s[usj]].am >
		 BasisSet.shells[Symmetry.us2s[usk]].am + BasisSet.shells[Symmetry.us2s[usl]].am)
		swap_ij_kl = 1;
	      else
		swap_ij_kl = 0;

	      skk = Symmetry.us2s[usk];
	      sll = Symmetry.us2s[usl];
	      if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
		/*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
		stab_k = Symmetry.atom_positions[BasisSet.shells[skk].center-1];
		stab_l = Symmetry.atom_positions[BasisSet.shells[sll].center-1];
		stab_kl = Symmetry.GnG[stab_k][stab_l];
		S_list = Symmetry.dcr[stab_k][stab_l];
		T_list = Symmetry.dcr[stab_ij][stab_kl];
		lambda_T = Symmetry.nirreps/Symmetry.dcr_deg[stab_ij][stab_kl];

		memset(sj_arr,0,sizeof(int)*max_num_unique_quartets);
		memset(sk_arr,0,sizeof(int)*max_num_unique_quartets);
		memset(sl_arr,0,sizeof(int)*max_num_unique_quartets);
		count = 0;
		
		for(dcr_ijkl=0;dcr_ijkl<Symmetry.dcr_dim[stab_ij][stab_kl];dcr_ijkl++){
		  T = T_list[dcr_ijkl];
		  sk = BasisSet.shells[skk].trans_vec[T]-1;
		  slll = BasisSet.shells[sll].trans_vec[T]-1;
		  for(dcr_kl=0;dcr_kl<Symmetry.dcr_dim[stab_k][stab_l];dcr_kl++) {
		    S = S_list[dcr_kl];
		    sl = BasisSet.shells[slll].trans_vec[S]-1;

		    total_am = BasisSet.shells[si].am +
			       BasisSet.shells[sj].am +
			       BasisSet.shells[sk].am +
			       BasisSet.shells[sl].am;
		    /*-------------------------------------------------------------
		      Obviously redundant or zero cases should be eliminated here!
		      Right now only zero case is eliminated. Redundancies arising
		      in DCD approach when usi == usj etc. may be eliminated too
		      but lambda_T will have to be replaced by an array (it won't
		      the same for every shell quartet in petite list anymore).
		     -------------------------------------------------------------*/
		    if(!(total_am%2)||
		       (BasisSet.shells[si].center!=BasisSet.shells[sj].center)||
		       (BasisSet.shells[sj].center!=BasisSet.shells[sk].center)||
		       (BasisSet.shells[sk].center!=BasisSet.shells[sl].center)) {
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
		total_am = BasisSet.shells[si].am +
			   BasisSet.shells[usj].am +
			   BasisSet.shells[usk].am +
			   BasisSet.shells[usl].am;
		if(!(total_am%2)||
		   (BasisSet.shells[si].center!=BasisSet.shells[usj].center)||
		   (BasisSet.shells[usj].center!=BasisSet.shells[usk].center)||
		   (BasisSet.shells[usk].center!=BasisSet.shells[usl].center)) {
		  num_unique_quartets = 1;
		  sj_arr[0] = usj;
		  sk_arr[0] = usk;
		  sl_arr[0] = usl;
		}
		else
		  num_unique_quartets = 0;
	      }

	      /*----------------------------------
	        Compute the nonredundant quartets
	       ----------------------------------*/
	      for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
		si = sii;
		sj = sj_arr[plquartet];
		sk = sk_arr[plquartet];
		sl = sl_arr[plquartet];
		if (BasisSet.schwartz_eri[si][sj]*BasisSet.schwartz_eri[sk][sl] < toler)
		  continue;
		/*--- As usual, we have to order bra-ket so that ket has the largest angular momentum */
		if (swap_ij_kl) {
		  dum = si;
		  si = sk;
		  sk = dum;
		  dum = sj;
		  sj = sl;
		  sl = dum;
		}
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
/*				if (fabs(sp_ij->Sovlp[pi][pj]*sp_kl->Sovlp[pk][pl]) > 1.0E-10)*/
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

		/*--- swap bra and ket back to the original order if needed ---*/
#if !SWAP1
		if (swap_ij_kl) {
		  dum = si;
		  si = sk;
		  sk = dum;
		  dum = sj;
		  sj = sl;
		  sl = dum;
		  sr_fao = BasisSet.shells[si].fao - 1;
		  ss_fao = BasisSet.shells[sj].fao - 1;
		  sp_fao = BasisSet.shells[sk].fao - 1;
		  sq_fao = BasisSet.shells[sl].fao - 1;
		  nr = ioff[BasisSet.shells[si].am];
		  ns = ioff[BasisSet.shells[sj].am];
		  np = ioff[BasisSet.shells[sk].am];
		  nq = ioff[BasisSet.shells[sl].am];
		  if (am /*&& (orig_am[0] + orig_am[1] != 0)*/) {
/*		    timer_on("Pre1Swap");*/
		      /*--- (pq|rs) -> (rs|pq) ---*/
		    ijkl_to_klij(data,scratch_buf,np*nq,nr*ns);
		    data = scratch_buf;
/*		    timer_off("Pre1Swap");*/
		  }
		}
		else {
		  sr_fao = BasisSet.shells[si].fao - 1;
		  ss_fao = BasisSet.shells[sj].fao - 1;
		  sp_fao = BasisSet.shells[sk].fao - 1;
		  sq_fao = BasisSet.shells[sl].fao - 1;
		  nr = ioff[BasisSet.shells[si].am];
		  ns = ioff[BasisSet.shells[sj].am];
		  np = ioff[BasisSet.shells[sk].am];
		  nq = ioff[BasisSet.shells[sl].am];
		}
#else
		if (swap_ij_kl) {
		  dum = si;
		  si = sk;
		  sk = dum;
		  dum = sj;
		  sj = sl;
		  sl = dum;
		  sr_fao = BasisSet.shells[si].fao - 1;
		  ss_fao = BasisSet.shells[sj].fao - 1;
		  sp_fao = BasisSet.shells[sk].fao - 1;
		  sq_fao = BasisSet.shells[sl].fao - 1;
		  nr = ioff[BasisSet.shells[si].am];
		  ns = ioff[BasisSet.shells[sj].am];
		  np = ioff[BasisSet.shells[sk].am];
		  nq = ioff[BasisSet.shells[sl].am];
		}
		else {
		  sr_fao = BasisSet.shells[si].fao - 1;
		  ss_fao = BasisSet.shells[sj].fao - 1;
		  sp_fao = BasisSet.shells[sk].fao - 1;
		  sq_fao = BasisSet.shells[sl].fao - 1;
		  nr = ioff[BasisSet.shells[si].am];
		  ns = ioff[BasisSet.shells[sj].am];
		  np = ioff[BasisSet.shells[sk].am];
		  nq = ioff[BasisSet.shells[sl].am];
		  if (am /*&& (orig_am[0] + orig_am[1] != 0)*/) {
/*		    timer_on("Pre1Swap");*/
		      /*--- (pq|rs) -> (rs|pq) ---*/
		    ijkl_to_klij(data,scratch_buf,nr*ns,np*nq);
		    data = scratch_buf;
/*		    timer_off("Pre1Swap");*/
		  }
		}
#endif

		/*--- step 1 of the transformation ---*/
/*		timer_on("Step 1");*/
#if !SWAP1
		if (usk != usl)
		for(r=0;r<nr;r++) {
		    for(s=0;s<ns;s++) {
			rs_offset = ((r*ns + s) * np * nq);
			rsi_offset = ((r*ns + s) * ibatch_length * BasisSet.num_ao);
			rsi_row = rsiq_buf + rsi_offset;
			for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=BasisSet.num_ao) {
			    mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
			    rspq_ptr = data + rs_offset;
			    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
				for(q=0,q_abs=sq_fao;
				    q<nq;
				    q++,q_abs++,rspq_ptr++) {
				    temp = (*rspq_ptr);
				    rsi_row[q_abs] += mo_vec[p_abs] * temp;
				    rsi_row[p_abs] += mo_vec[q_abs] * temp;
				}
			    }
			}
		    }
		}
		else
		for(r=0;r<nr;r++) {
		    for(s=0;s<ns;s++) {
			rs_offset = ((r*ns + s) * np * nq);
			rsi_offset = ((r*ns + s) * ibatch_length * BasisSet.num_ao);
			rsi_row = rsiq_buf + rsi_offset;
			for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=BasisSet.num_ao) {
			    mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
			    rspq_ptr = data + rs_offset;
			    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
				for(q=0,q_abs=sq_fao;
				    q<nq;
				    q++,q_abs++,rspq_ptr++) {
				    temp = (*rspq_ptr);
				    rsi_row[q_abs] += mo_vec[p_abs] * temp;
				}
			    }
			}
		    }
		}
#else
		if (usk != usl)
/*		for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		    mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
		    rspq_ptr = data;
		    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			for(q=0,q_abs=sq_fao;
			    q<nq;
			    q++,q_abs++) {
			    iq_row = rsiq_buf + ((mo_i*BasisSet.num_ao + q_abs) * nr * ns);
			    ip_row = rsiq_buf + ((mo_i*BasisSet.num_ao + p_abs) * nr * ns);
			    temp1 = mo_vec[q_abs];
			    temp2 = mo_vec[p_abs];
#if !USE_BLAS
			    for(rs=0;rs<nr*ns;rs++,rspq_ptr++,iq_row++,ip_row++) {
				(*iq_row) += temp2 * (*rspq_ptr);
				(*ip_row) += temp1 * (*rspq_ptr);
			    }
#else
			    C_DAXPY(nr*ns,temp2,rspq_ptr,1,iq_row,1);
			    C_DAXPY(nr*ns,temp1,rspq_ptr,1,ip_row,1);
			    rspq_ptr += nr*ns;
#endif
			}
		    }
		}*/
		for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		    mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
		    rspq_ptr = data;
		    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			for(q=0,q_abs=sq_fao;
			    q<nq;
			    q++,q_abs++) {
			    ip_row = rsiq_buf + ((mo_i*BasisSet.num_ao + p_abs) * nr * ns);
			    temp1 = mo_vec[q_abs];
#if !USE_BLAS
			    for(rs=0;rs<nr*ns;rs++,rspq_ptr++,iq_row++,ip_row++) {
				(*ip_row) += temp1 * (*rspq_ptr);
			    }
#else
			    C_DAXPY(nr*ns,temp1,rspq_ptr,1,ip_row,1);
			    rspq_ptr += nr*ns;
#endif
			}
		    }
		    rspq_ptr = data;
		    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
#if !USE_BLAS
			for(qrs=0;qrs<nq*nr*ns;qrs++,i_row++,rspq_ptr++) {
			    (*i_row) += temp * (*rspq_ptr);
			}
#else
			C_DAXPY(nq*nr*ns,temp,rspq_ptr,1,i_row,1);
			rspq_ptr += nq*nr*ns;
#endif
		    }
		}
		else
		for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		    mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
		    rspq_ptr = data;
		    for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
#if !USE_BLAS
			for(qrs=0;qrs<nq*nr*ns;qrs++,i_row++,rspq_ptr++) {
			    (*i_row) += temp * (*rspq_ptr);
			}
#else
			C_DAXPY(nq*nr*ns,temp,rspq_ptr,1,i_row,1);
			rspq_ptr += nq*nr*ns;
#endif
		    }
		}
#endif
/*		timer_off("Step 1");*/
	      
	      } /* end of computing "petit" list - end of P,Q loop */
	    } /* end of "unique" P,Q loop */
	  
	  /*--- step 2 of the transfromation ---*/
/*	  rsi_offset = ((r*ns + s) * ibatch_length * BasisSet.num_ao);*/
#if SWAP1
/*	  timer_on("Post1Swap");*/
	  ijkl_to_klij(rsiq_buf,scratch_buf,ibatch_length*BasisSet.num_ao,nr*ns);
	  rsi_row = scratch_buf;
/*	  timer_off("Post1Swap");*/
#else
	  rsi_row = rsiq_buf;
#endif
	  ia_block_ptr = rsia_buf;
	  for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
		  /*--- Can be done as a matrix multiply ---*/
/*		  timer_on("Step 2");*/
#if !USE_BLAS
		  ia = 0;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=BasisSet.num_ao) {
		      for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++,ia++) {
			  mo_vec = MOInfo.scf_evec_uocc[0][mo_a];
			  temp = 0.0;
			  for(q_abs=0;q_abs<BasisSet.num_ao;q_abs++) {
			      temp += mo_vec[q_abs] * rsi_row[q_abs];
			  }
			  ia_block_ptr[ia] = temp;
		      }
		  }
#else
		  C_DGEMM('n','t',ibatch_length,MOInfo.nactuocc,BasisSet.num_ao,1.0,
			  rsi_row,BasisSet.num_ao,MOInfo.scf_evec_uocc[0][0],BasisSet.num_ao,
			  0.0,ia_block_ptr,MOInfo.nactuocc);
		  rsi_row += BasisSet.num_ao*ibatch_length;
#endif
		  ia_block_ptr += MOInfo.nactuocc*ibatch_length;
/*		  timer_off("Step 2");*/
	      }
	  }

	  /*--- step 3 of the transformation ---*/
	  rsi_row = rsia_buf;
	  /*--- To update (JS|IA) need to lock mutex corresponding to the S and R indices ---*/
#if LOCK_RS_SHELL	  
	  pthread_mutex_lock(&rmp2_sindex_mutex[INDEX(si,sj)]);
#endif
	  for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
/*		  timer_on("Step 3");*/
		  r_abs = r + sr_fao;
		  s_abs = s + ss_fao;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=MOInfo.nactuocc) {
		      for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
#if !LOCK_RS_SHELL
			  pthread_mutex_lock(&rmp2_sindex_mutex[s_abs]);
#endif
			  jsi_row = jsia_buf + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.nactuocc;
			  temp = MOInfo.scf_evec_occ[0][mo_j+jmin][r_abs];
			  for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
			      jsi_row[mo_a] += temp * rsi_row[mo_a];
			  }
#if !LOCK_RS_SHELL
			  pthread_mutex_unlock(&rmp2_sindex_mutex[s_abs]);
#endif
			  if (usi != usj) {
#if !LOCK_RS_SHELL
			    pthread_mutex_lock(&rmp2_sindex_mutex[r_abs]);
#endif
			    jsi_row = jsia_buf + ((mo_j * BasisSet.num_ao + r_abs) * ibatch_length + mo_i) * MOInfo.nactuocc;
			    temp = MOInfo.scf_evec_occ[0][mo_j+jmin][s_abs];
			    for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
			      jsi_row[mo_a] += temp * rsi_row[mo_a];
			    }
#if !LOCK_RS_SHELL
			    pthread_mutex_unlock(&rmp2_sindex_mutex[r_abs]);
#endif
			  }
		      }
		  }
/*		  timer_off("Step 3");*/
	      }
	  }
#if LOCK_RS_SHELL
	  pthread_mutex_unlock(&rmp2_sindex_mutex[INDEX(si,sj)]);
#endif
	  memset(rsiq_buf,0,nr*ns*ibatch_length*BasisSet.num_ao*sizeof(double));
	  memset(rsia_buf,0,nr*ns*ibatch_length*MOInfo.nactuocc*sizeof(double));

	} /* end of R,S loop */
      } /* end of "unique" R,S loop */

    pthread_mutex_lock(&rmp2_energy_mutex);
    RMP2_Status.num_arrived++;
    if (RMP2_Status.num_arrived != UserOptions.num_threads) { /*--- there are some threads still working - wait here --*/
      pthread_cond_wait(&rmp2_energy_cond,&rmp2_energy_mutex);
      pthread_mutex_unlock(&rmp2_energy_mutex);
    }
    else { /*--- this is the last thread to get here - do the 4th step and energy calculation alone and wake everybody up ---*/
      /*--- step 4 of the transformation ---*/
/*    timer_on("Step 4");*/
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
	  for(mo_b=0;mo_b<MOInfo.nactuocc;mo_b++) {
	    jbi_row = jbia_buf + ((mo_j * MOInfo.nactuocc + mo_b) * ibatch_length + mo_i) * MOInfo.nactuocc;
	    for(s_abs=0;s_abs<BasisSet.num_ao;s_abs++) {
	      jsi_row = jsia_buf + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.nactuocc;
	      temp = MOInfo.scf_evec_uocc[0][mo_b][s_abs];
	      for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
		jbi_row[mo_a] += temp * jsi_row[mo_a];
	      }
	    }
	  }
	}
      }
/*    timer_off("Step 4");*/

      /*--- Finish the symmetrization step - zero out non-totally symmetric integrals in Abelian case */
      if (Symmetry.nirreps > 1)
	for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	  for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
	    for(mo_b=0;mo_b<MOInfo.nactuocc;mo_b++) {
	      for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
		if ((MOInfo.mo2symblk_occ[0][mo_i+imin] ^ MOInfo.mo2symblk_occ[0][mo_j+jmin]) ^
		    (MOInfo.mo2symblk_uocc[0][mo_a] ^ MOInfo.mo2symblk_uocc[0][mo_b]))
		    jbia_buf[((mo_j * MOInfo.nactuocc + mo_b) * ibatch_length + mo_i) * MOInfo.nactuocc + mo_a] = 0.0;
	      }
	    }
	  }
	}

#if PRINT
    /*--- Print them out if needed ---*/
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
	  for(mo_b=0;mo_b<MOInfo.nactuocc;mo_b++) {
	    for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
	      temp = jbia_buf[((mo_j * MOInfo.nactuocc + mo_b) * ibatch_length + mo_i) * MOInfo.nactuocc + mo_a];
	      if (fabs(temp) > ZERO)
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_j+jmin,mo_b,mo_i+imin,mo_a,
			mo_j*MOInfo.nactuocc+mo_b,
			mo_i*MOInfo.nactuocc+mo_a,
			temp);
	    }
	  }
	}
      }
#endif

    /*--- Compute a contribution to the MP2 energy from the current batch ---*/
/*    timer_on("Energy");*/
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
	  /*--- Evaluate pair energies ---*/
	  mo_ij = INDEX(mo_i + imin, mo_j + jmin);
	  e0 = 0.0;
	  e1 = 0.0;
	  for(mo_a=0;mo_a<MOInfo.nactuocc;mo_a++) {
	    for(mo_b=0;mo_b<=mo_a;mo_b++) {
	      iajb = jbia_buf[((mo_j * MOInfo.nactuocc + mo_b) * ibatch_length + mo_i) * MOInfo.nactuocc + mo_a];
	      ibja = jbia_buf[((mo_j * MOInfo.nactuocc + mo_a) * ibatch_length + mo_i) * MOInfo.nactuocc + mo_b];
	      pfac = (mo_i + imin != mo_j + jmin) ? 1.0 : m_sqrt1_2;
	      pfac *= (mo_a != mo_b) ? 1.0 : m_sqrt1_2;
	      k0ijab = pfac * (iajb + ibja);
	      eijab = (MOInfo.scf_evals_occ[0][mo_i+imin] +
		       MOInfo.scf_evals_occ[0][mo_j+jmin] -
		       MOInfo.scf_evals_uocc[0][mo_a] -
		       MOInfo.scf_evals_uocc[0][mo_b]);
	      e0 += k0ijab*k0ijab/eijab;
	      if (mo_i + imin != mo_j + jmin) {
		k1ijab = pfac * (iajb - ibja);
		e1 += 3*k1ijab*k1ijab/eijab;
	      }
	    }
	  }
	  RMP2_Status.emp2_0[mo_ij] = e0;
	  RMP2_Status.Emp2_0 += e0;
	  if (mo_i + imin != mo_j + jmin) {
	    RMP2_Status.emp2_1[mo_ij] = e1;
	    RMP2_Status.Emp2_1 += e1;
	  }
	}
      }
/*    timer_off("Energy");*/
      if (ibatch < num_ibatch-1) {
	fprintf(outfile,"  Energy after Pass #%d = %20.10lf\n",ibatch,RMP2_Status.Emp2_0+RMP2_Status.Emp2_1);
	fprintf(outfile,"  Singlet energy after Pass #%d = %20.10lf\n",ibatch,RMP2_Status.Emp2_0);
	fprintf(outfile,"  Triplet energy after Pass #%d = %20.10lf\n",ibatch,RMP2_Status.Emp2_1);
	fflush(outfile);
	memset(jsia_buf,0,MOInfo.nactdocc*BasisSet.num_ao*ibatch_length*MOInfo.nactuocc*sizeof(double));
	memset(jbia_buf,0,MOInfo.nactdocc*MOInfo.nactuocc*ibatch_length*MOInfo.nactuocc*sizeof(double));

	fprintf(outfile,"\n");
	fprintf(outfile,"  Singlet pair energies after Pass #%d:\n",ibatch);
	fprintf(outfile,"    i       j         e(ij)\n");
	fprintf(outfile,"  -----   -----   ------------\n");
	for(mo_i=0;mo_i<ibatch_length;mo_i++)
	  for(mo_j=0;mo_j<=mo_i+imin-jmin;mo_j++) {
	    mo_ij = INDEX(mo_i + imin, mo_j + jmin);
	    fprintf(outfile,"  %3d     %3d     %12.9lf\n",imin+mo_i+1,jmin+mo_j+1,RMP2_Status.emp2_0[mo_ij]);
	  }
	fprintf(outfile,"\n");
	fprintf(outfile,"  Triplet pair energies after Pass #%d:\n",ibatch);
	fprintf(outfile,"    i       j         e(ij)\n");
	fprintf(outfile,"  -----   -----   ------------\n");
	for(mo_i=0;mo_i<ibatch_length;mo_i++)
	  for(mo_j=0;mo_j<mo_i+imin-jmin;mo_j++) {
	    mo_ij = INDEX(mo_i + imin, mo_j + jmin);
	    fprintf(outfile,"  %3d     %3d     %12.9lf\n",imin+mo_i+1,jmin+mo_j+1,RMP2_Status.emp2_1[mo_ij]);
	  }
	fprintf(outfile,"\n");
	fflush(outfile);
      }
      /*--- Done with the non-threaded part - wake everybody up and prepare for the next ibatch ---*/
      RMP2_Status.num_arrived = 0;
      pthread_cond_broadcast(&rmp2_energy_cond);
      pthread_mutex_unlock(&rmp2_energy_mutex);
    }
  } /* End of "I"-loop */

  /*---------
    Clean-up
   ---------*/
  free(rsia_buf);
  free(rsiq_buf);
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

}
}
