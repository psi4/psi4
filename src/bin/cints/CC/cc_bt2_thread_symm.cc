/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#if HAVE_CMATH
# include <cmath>
#else
# include <cmath>
#endif

#include <cstring>
#include<cstdio>
#include<memory.h>
#include<cstdlib>
#include<pthread.h>

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
#include"int_fjt.h"
#endif

namespace psi { namespace CINTS {
  void *cc_bt2_thread_symm(void *tnum_ptr)
  {
    const long int thread_num = (long int) tnum_ptr;
    const double toler = UserOptions.cutoff;
    const double m_sqrt1_2 = 1/sqrt(2.0);
    int lambda_T = 1;
    
    int nao = BasisSet.num_ao;
    int nocc = CCInfo.nocc;
    
    /*---------------
      Initialization
      ---------------*/
#ifndef USE_TAYLOR_FM
    double_array_t fjt_table;
    init_fjt_table(&fjt_table);
#endif
    
    /*-------------------------
      Allocate data structures
      -------------------------*/
    int max_bf_per_shell = ioff[BasisSet.max_am];
    int max_cart_class_size = (max_bf_per_shell)*
      (max_bf_per_shell)*
      (max_bf_per_shell)*
      (max_bf_per_shell);
    int max_class_size;
    if (Symmetry.nirreps > 1)
      max_class_size = max_cart_class_size;
    int max_num_unique_quartets = Symmetry.max_stab_index*
      Symmetry.max_stab_index*
      Symmetry.max_stab_index;
    int* sj_arr = new int[max_num_unique_quartets];
    int* sk_arr = new int[max_num_unique_quartets];
    int* sl_arr = new int[max_num_unique_quartets];
    int max_num_prim_comb = BasisSet.max_num_prims*
      BasisSet.max_num_prims*
      BasisSet.max_num_prims*
      BasisSet.max_num_prims;
    Libint_t Libint;
    init_libint(&Libint, BasisSet.max_am-1, max_num_prim_comb);
    
    double* raw_data;
#ifdef NONDOUBLE_INTS
    raw_data = init_array(max_cart_class_size);
#endif
    
    
    /*--------------------------------------------
      generate all symmetry unique shell quartets
      --------------------------------------------*/
    for (int usii=0; usii<Symmetry.num_unique_shells; usii++)
      for (int uskk=0; uskk<Symmetry.num_unique_shells; uskk++) {
	int ik_index = usii*Symmetry.num_unique_shells + uskk;
	
	/*--- Decide if this thread will do this ---*/
	if ( ik_index%UserOptions.num_threads != thread_num )
	  continue;
	
	for (int usjj=0; usjj<=usii; usjj++)
	  for (int usll=0; usll<=uskk; usll++){
	    
	    int usi = usii;
	    int usj = usjj;
	    int usk = uskk;
	    int usl = usll;
	    
	    int si = Symmetry.us2s[usi];
	    int sjj = Symmetry.us2s[usj];
	    int skk = Symmetry.us2s[usk];
	    int sll = Symmetry.us2s[usl];
	    
	    int num_unique_quartets = 0;
	    
	    if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
	      /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
	      int stab_i = Symmetry.atom_positions[BasisSet.shells[si].center-1];
	      int stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
	      int stab_k = Symmetry.atom_positions[BasisSet.shells[skk].center-1];
	      int stab_l = Symmetry.atom_positions[BasisSet.shells[sll].center-1];
	      int stab_ij = Symmetry.GnG[stab_i][stab_j];
	      int stab_kl = Symmetry.GnG[stab_k][stab_l];
	      int* R_list = Symmetry.dcr[stab_i][stab_j];
	      int* S_list = Symmetry.dcr[stab_k][stab_l];
	      int* T_list = Symmetry.dcr[stab_ij][stab_kl];
	      lambda_T = Symmetry.nirreps/Symmetry.dcr_deg[stab_ij][stab_kl];
	      
	      memset(sj_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sk_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sl_arr,0,sizeof(int)*max_num_unique_quartets);
	      int count = 0;
	      for(int dcr_ij=0;dcr_ij<Symmetry.dcr_dim[stab_i][stab_j];dcr_ij++){
		int R = R_list[dcr_ij];
		int sj = BasisSet.shells[sjj].trans_vec[R]-1;
		for(int dcr_ijkl=0;dcr_ijkl<Symmetry.dcr_dim[stab_ij][stab_kl];dcr_ijkl++){
		  int T = T_list[dcr_ijkl];
		  int sk = BasisSet.shells[skk].trans_vec[T]-1;
		  int slll = BasisSet.shells[sll].trans_vec[T]-1;
		  for(int dcr_kl=0;dcr_kl<Symmetry.dcr_dim[stab_k][stab_l];dcr_kl++) {
		    int S = S_list[dcr_kl];
		    int sl = BasisSet.shells[slll].trans_vec[S]-1;
		    
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
	      lambda_T = 1;
	      sj_arr[0] = usj;
	      sk_arr[0] = usk;
	      sl_arr[0] = usl;
	    }
	    
	    
	    /*----------------------------------
	      Compute the nonredundant quartets
	      ----------------------------------*/
	    for(int plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	      int si = Symmetry.us2s[usii];
	      int sj = sj_arr[plquartet];
	      int sk = sk_arr[plquartet];
	      int sl = sl_arr[plquartet];
	      
	      int total_am = BasisSet.shells[si].am + BasisSet.shells[sj].am + 
		BasisSet.shells[sk].am + BasisSet.shells[sl].am;
	      /* parity selection */
	      if (total_am%2 &&
		  BasisSet.shells[si].center==BasisSet.shells[sj].center &&
		  BasisSet.shells[sj].center==BasisSet.shells[sk].center &&
		  BasisSet.shells[sk].center==BasisSet.shells[sl].center)
		continue;
	      
	      int switch_ij = 0;
	      int switch_kl = 0;
	      int switch_ijkl = 0;
	      /* place in "ascending" angular mom-
		 my simple way of optimizing PHG recursion (VRR) */
	      /* these first two are good for the HRR */
	      if(BasisSet.shells[si].am < BasisSet.shells[sj].am){
		int dum = si;
		si = sj;
		sj = dum;
		switch_ij = 1;
	      }
	      if(BasisSet.shells[sk].am < BasisSet.shells[sl].am){
		int dum = sk;
		sk = sl;
		sl = dum;
		switch_kl = 1;
	      }
	      /* this should be /good/ for the VRR */
	      if(BasisSet.shells[si].am + BasisSet.shells[sj].am > BasisSet.shells[sk].am + BasisSet.shells[sl].am){
		int dum = si;
		si = sk;
		sk = dum;
		dum = sj;
		sj = sl;
		sl = dum;
		switch_ijkl = 1;
	      }
	      
	      int np_i = BasisSet.shells[si].n_prims;
	      int np_j = BasisSet.shells[sj].n_prims;
	      int np_k = BasisSet.shells[sk].n_prims;
	      int np_l = BasisSet.shells[sl].n_prims;
	      int ni = ioff[BasisSet.shells[si].am];
	      int nj = ioff[BasisSet.shells[sj].am];
	      int nk = ioff[BasisSet.shells[sk].am];
	      int nl = ioff[BasisSet.shells[sl].am];
	      int si_fao = BasisSet.shells[si].fao-1;
	      int sj_fao = BasisSet.shells[sj].fao-1;
	      int sk_fao = BasisSet.shells[sk].fao-1;
	      int sl_fao = BasisSet.shells[sl].fao-1;
	      int orig_am[4];
	      orig_am[0] = BasisSet.shells[si].am-1;
	      orig_am[1] = BasisSet.shells[sj].am-1;
	      orig_am[2] = BasisSet.shells[sk].am-1;
	      orig_am[3] = BasisSet.shells[sl].am-1;
	      int am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];
	      
	      struct shell_pair* sp_ij = &(BasisSet.shell_pairs[si][sj]);
	      struct shell_pair* sp_kl = &(BasisSet.shell_pairs[sk][sl]);
	      
	      Libint.AB[0] = sp_ij->AB[0];
	      Libint.AB[1] = sp_ij->AB[1];
	      Libint.AB[2] = sp_ij->AB[2];
	      Libint.CD[0] = sp_kl->AB[0];
	      Libint.CD[1] = sp_kl->AB[1];
	      Libint.CD[2] = sp_kl->AB[2];
	      
	      double AB2 = Libint.AB[0]*Libint.AB[0]+Libint.AB[1]*Libint.AB[1]+Libint.AB[2]*Libint.AB[2];
	      double CD2 = Libint.CD[0]*Libint.CD[0]+Libint.CD[1]*Libint.CD[1]+Libint.CD[2]*Libint.CD[2];
	      
	      /*--- Compute data for primitive quartets here ---*/
	      int num_prim_comb = 0;
	      for (int pi = 0; pi < np_i; pi++) {
		int max_pj = (si == sj) ? pi+1 : np_j;
		for (int pj = 0; pj < max_pj; pj++) {
		  int m = (1 + (si == sj && pi != pj));
		  for (int pk = 0; pk < np_k; pk++) {
		    int max_pl = (sk == sl) ? pk+1 : np_l;
		    for (int pl = 0; pl < max_pl; pl++){
		      int n = m * (1 + (sk == sl && pk != pl));
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
	      double* data;
	      if (am) {
#ifdef NONDOUBLE_INTS
		int size = ni * nj * nk * nl;
		REALTYPE target_ints = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
	      for(i=0;i<size;i++)
		raw_data[i] = (double) target_ints[i];
#else
	      raw_data = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
#endif
	      /* No need to transforms integrals to sph. harm. basis */
	      data = norm_quartet(raw_data, NULL, orig_am, 0);
	      }
	      else {
		double temp = 0.0;
		for(int p=0;p<num_prim_comb;p++)
		  temp += (double) Libint.PrimQuartet[p].F[0];
#ifdef NONDOUBLE_INTS
		raw_data[0] = temp;
		data = raw_data;
#else
		Libint.int_stack[0] = temp;
		data = Libint.int_stack;
#endif
	      }
	      
	      int ijkl = 0;
	      for(int i=0;i<ni;i++) {
		int ii = i + si_fao; 
		for(int j=0;j<nj;j++) {
		  int jj = j + sj_fao; 
		  int IJ = INDEX(ii,jj);
		  for(int k=0;k<nk;k++) {
		    int kk = k + sk_fao; 
		    for(int l=0;l<nl;l++,ijkl++) {
		      int ll = l + sl_fao; 
		      int KL = INDEX(kk,ll);
		      
		      double value = data[ijkl];
		      
		      /*
			if(fabs(value) > 1e-8)
			fprintf(outfile, "%d %d %d %d %20.14f\n", ii+1, jj+1, kk+1, ll+1, value);
		      */
		      
		      int ij = ii*nao + jj; int ji = jj*nao + ii;
		      int ik = ii*nao + kk; int ki = kk*nao + ii;
		      int il = ii*nao + ll; int li = ll*nao + ii;
		      int jk = jj*nao + kk; int kj = kk*nao + jj;
		      int jl = jj*nao + ll; int lj = ll*nao + jj;
		      int kl = kk*nao + ll; int lk = ll*nao + kk;
		      
		      // (ij|kl)
		      C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jl], 1, CCInfo.T2_t[ik],1);
		      
		      /*
			if(ii!=jj && kk!=ll && IJ!=KL) {
			// (ij|lk)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			// (ji|kl)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			// (ji|lk)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ik], 1, CCInfo.T2_t[jl],1);
			// (kl|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			// (kl|ji)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[li], 1, CCInfo.T2_t[kj],1);
			// (lk|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[kj], 1, CCInfo.T2_t[li],1);
			// (lk|ji)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ki], 1, CCInfo.T2_t[lj],1);
			}
			else if(ii!=jj && kk!=ll && IJ==KL) {
			// (ij|lk)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			// (ji|kl)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			// (ji|lk)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[ik], 1, CCInfo.T2_t[jl],1);
			}
			else if(ii!=jj && kk==ll) {
			// (ji|kl)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[il], 1, CCInfo.T2_t[jk],1);
			// (kl|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			// (kl|ji)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[li], 1, CCInfo.T2_t[kj],1);
			}
			else if(ii==jj && kk!=ll) {
			// (ij|lk)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[jk], 1, CCInfo.T2_t[il],1);
			// (kl|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			// (lk|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[kj], 1, CCInfo.T2_t[li],1);
			}
			else if(ii==jj && kk==ll && IJ!=KL) {
			// (kl|ij)
			C_DAXPY(nocc*nocc, value, CCInfo.T2_s[lj], 1, CCInfo.T2_t[ki],1);
			}
		      */
		    }
		  }
		}
	      }
	      
	    } /* end of RSPQ loop */
	  } /* end of "unique" RSQP loop */
      }
    
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
  
};}; // namespace 
