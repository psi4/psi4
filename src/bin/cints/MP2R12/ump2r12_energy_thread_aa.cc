#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<cstdlib>
#include<libciomr/libciomr.h>
#include<libqt/qt.h>
#include<libpsio/psio.h>
#include<libint/libint.h>
#include<libr12/libr12.h>
#include<pthread.h>

#include"defines.h"
#include"psifiles.h"
#define EXTERN
#include"global.h"
#include"r12_quartet_data.h"
#include"norm_quartet.h"
#include"int_fjt.h"
#include"quartet_permutations.h"
#include"rmp2r12_energy.h"

namespace psi{ namespace CINTS{ namespace ump2r12_aa{

/*-------------------------------------------------------
  Algorithm

  ***Split into threads. Each thread do the following:
  
  Loop over I batches (batch size num_i_per_batch) of active DOCC

    Loop over all symmetry-unique shells UR, US<=UR
      ***if this UR, US is to be handled by this thread - proceed, else skip to next one
      Find all symmetry-distinct shell doublets resulting from (UR US|
      Loop over all resulting shell doublets R, S

        Loop over all symmetry-unique shells UP, UQ<=UP
	  Find all symmetry-distinct shell quartets resulting from (R S|UP UQ)
	  Loop over the resulting set of P, Q doublets
            Evaluate (RS|PQ), (RS|r12|PQ), and (RS|[r12,T1]|PQ)
	    Loop over p in P, q in Q, r in R, s in S, i in I
	      (rs|iq) += Cpi * (rs|pq)
	      (rs|ip) += Cqi * (rs|pq)
	      same for (rs|r12|pq) and (rs|[r12,T1]|pq)
	    End p, q, r, s, i loop
          End P,Q loop
        End UP, UQ loop

        Loop over r in R, s in S
          Loop over q < nao, x < num_mo, i in I
	    (rs|ix) += Cxs * (rs|is)
	    same for (rs|r12|is) and (rs|[r12,T1]|is)
	  End q, x, i
	End r, s loop
	***Lock (js| and (jr| (either individual or shell blocks depending on LOCK_RS_SHELL in rmp2r12_energy.h)
	Loop over r in R, s in S
	  Loop over i in I, x < num_mo, j <= i
	    (js|ix) += Cjr * (rs|ix)
	    (jr|ix) += Cjs * (rs|ix)
	    same for (rs|r12|ix), but
	    (js|[r12,T1]|ix) += Cjr * (rs|[r12,T1]|ix)
	    (jr|[r12,T1]|ix) -= Cjs * (rs|[r12,T1]|ix)                   <---- Note the minus sign here!!!
	  End i, x, j loop
        End r, s loop
	***Unlock (js| and (jr|

      End R, S loop
    End UR, US loop

    ***Barrier: threads wait until everyone is done
    ***Do the following in one thread only
    Loop over i in I, j <= i
      Loop over r < nao, x < num_mo, y < num_mo
        (jy|ix) += Cys * (js|ix)
	same for (js|r12|ix) and (js|[r12,T1]|ix)
      End r, x, y loop
    End i, j loop

  End I loop

  ***Merge all threads

 -------------------------------------------------------*/


void *ump2r12_energy_thread_aa(void *tnum_ptr)
{
  const long int thread_num = (long int) tnum_ptr;
  const double toler = UserOptions.cutoff;

  extern pthread_mutex_t rmp2r12_energy_mutex;
  extern pthread_mutex_t *rmp2r12_sindex_mutex;
  extern pthread_cond_t rmp2r12_energy_cond;
  extern RMP2R12_Status_t RMP2R12_Status;
  extern double *jsix_buf[NUM_TE_TYPES];             /* buffer for (js|ia) integrals, where j runs over all alpha occ act MOs,
						 s runs over all AOs, i - over I-batch, x - over all MOs */
  extern double *jyix_buf[NUM_TE_TYPES];             /* buffer contains all MP2-R12/A-type integrals (AA|AA) */
  extern double *jsix_buf2[NUM_TE_TYPES];             /* buffer for (js|ia) integrals, where j runs over all beta occ act MOs,
						 s runs over all AOs, i - over I-batch, x - over all MOs */
  extern double *jyix_buf2[NUM_TE_TYPES];             /* buffer contains all MP2-R12/A-type integrals (AA|BB)*/
  
  extern int *first, *last;                          /* first and last absolute (Pitzer) orbital indices in symblk */
  extern int *fstocc_alpha, *lstocc_alpha;           /* first and last occupied indices in Pitzer ordering for each symblk */
  extern int *fstocc_beta, *lstocc_beta;             /* first and last occupied indices in Pitzer ordering for each symblk */
  extern int *occ_alpha, *occ_beta;                  /* Pitzer to "full"(no frozen core) QTS index mapping */
  extern int *act2fullQTS_alpha;          /* Maps "active"(taking frozen core into account) QTS into "frozen" QTS index */
  extern int *act2fullQTS_beta;           /* Maps "active"(taking frozen core into account) QTS into "frozen" QTS index */
  extern int *ioff3;                                 /* returns pointers to the beginning of rows in a rectangular matrix */

  struct shell_pair *sp_ij, *sp_kl;
  Libr12_t Libr12;
  double_array_t fjt_table;

  int ij, kl, ik, jl, ijkl;
  int count, dum, status;
  int te_type;
  int n, num[NUM_TE_TYPES];
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int i, j, m, p, q, r, s, x, y;
  int isym, jsym, xsym, ysym;
  int p_abs, q_abs, r_abs, s_abs;
  int si, sj, sk, sl;
  int sii, sjj, skk, sll , slll;
  int num_ij, swap_ij_kl;
  int pi, pj, pk, pl ;
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

  int size, class_size;
  int max_class_size;
  int max_cart_class_size;

  int np_i, np_j, np_k, np_l;
  int nr, ns, np, nq;

  int num_ibatch, num_i_per_ibatch, ibatch, ibatch_first, ibatch_length;
  int imin, imax, jmin;
  int max_bf_per_shell;
  int mo_i, mo_j, mo_x, mo_y;
  int ix, xy;
  int rs, qrs;
  int rs_offset, rsi_offset, rsp_offset;
  int num_prim_comb;

  char ij_key_string[80];
  int iounits[NUM_TE_TYPES];                  /* integrals file numbers */
  int iounits2[NUM_TE_TYPES];                  /* integrals file numbers */
  double AB2, CD2;
  double *raw_data[NUM_TE_TYPES];             /* pointers to the unnormalized target quartet of integrals */
  double *data[NUM_TE_TYPES];                 /* pointers to the transformed normalized target quartet of integrals */

  double temp;
  double ssss, ss_r12_ss;
  double *rspq_ptr;
  double *mo_vec;
  double *rsiq_buf[NUM_TE_TYPES];             /* buffer for (rs|iq) integrals, where r,s run over shell sets,
						 i runs over I-batch, q runs over all AOs */
  double *rsi_row;
  double *ix_block_ptr;
  double *rsix_buf[NUM_TE_TYPES];                           /* buffer for (rs|iq) integrals, where r,s run over shell sets,
						 i runs over I-batch, x runs over all MOs */
  double *i_row;
  double **jy_buf;
  double *xy_buf;
  double *jsi_row;
  double *jyi_row;
  double temp1,temp2,*iq_row,*ip_row;
  double *scratch_buf;


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
#ifdef NONDOUBLE_INTS
  for(i=0;i<NUM_TE_TYPES;i++)
    raw_data[i] = init_array(max_cart_class_size);
#endif
  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  init_libr12(&Libr12,BasisSet.max_am-1,max_num_prim_comb);
  init_fjt_table(&fjt_table);

  num_i_per_ibatch = RMP2R12_Status.num_i_per_ibatch;
  num_ibatch = RMP2R12_Status.num_ibatch;
  ibatch_first = RMP2R12_Status.ibatch_first;
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    rsiq_buf[te_type] = init_array(num_i_per_ibatch*BasisSet.num_ao*
				   max_bf_per_shell*max_bf_per_shell);
    rsix_buf[te_type] = init_array(num_i_per_ibatch*MOInfo.num_mo*
				   max_bf_per_shell*max_bf_per_shell);
  }
  jy_buf = block_matrix(MOInfo.alpha_occ,MOInfo.num_mo);
  xy_buf = init_array(MOInfo.num_mo*MOInfo.num_mo);
  scratch_buf = init_array(MAX(MOInfo.alpha_act_occ*MOInfo.num_mo*
			       num_i_per_ibatch*MOInfo.num_mo,
			       MAX(max_cart_class_size,
				   num_i_per_ibatch*BasisSet.num_ao*
				   max_bf_per_shell*max_bf_per_shell)));

/*-----------------------------------
  generate all unique shell quartets
 -----------------------------------*/
  /*--- I-batch loop ---*/
  for (ibatch=ibatch_first;ibatch<num_ibatch;ibatch++) {
    imin = ibatch * num_i_per_ibatch + MOInfo.nfrdocc;
    imax = MIN( imin+num_i_per_ibatch , MOInfo.alpha_occ );
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

	/*----------------------------------------------------------------------
	  NOTE on swapping usi, usj, etc.:

	  For 2-electron integrals of Hermitian operators it does not matter
	  if we swap si and sj, or sk and sl, or even si,sj and sk,sl. It
	  matters here though since [r12,T1] is non-Hermitian! What we want
	  in the end are the integrals of the [r12,T1] operator of this type:
	  (si sj|[r12,T1]|sk sl). If we have to swap bra and ket in the process
	  we will end up with having to compute (sk sl|[r12,T2]|si sj) instead.
	  Therefore if we have to swap bra and ket (swap_ij_kl = 1), then we
	  have to take [r12,T2] integral instead and swap it's bra and ket back.
	 ----------------------------------------------------------------------*/
	
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

		Libr12.ShellQuartet.AB[0] = sp_ij->AB[0];
		Libr12.ShellQuartet.AB[1] = sp_ij->AB[1];
		Libr12.ShellQuartet.AB[2] = sp_ij->AB[2];
		Libr12.ShellQuartet.CD[0] = sp_kl->AB[0];
		Libr12.ShellQuartet.CD[1] = sp_kl->AB[1];
		Libr12.ShellQuartet.CD[2] = sp_kl->AB[2];
		Libr12.ShellQuartet.AC[0] = Molecule.centers[BasisSet.shells[si].center-1].x-
			Molecule.centers[BasisSet.shells[sk].center-1].x;
		Libr12.ShellQuartet.AC[1] = Molecule.centers[BasisSet.shells[si].center-1].y-
			Molecule.centers[BasisSet.shells[sk].center-1].y;
		Libr12.ShellQuartet.AC[2] = Molecule.centers[BasisSet.shells[si].center-1].z-
			Molecule.centers[BasisSet.shells[sk].center-1].z;
		Libr12.ShellQuartet.ABdotAC = Libr12.ShellQuartet.AB[0]*Libr12.ShellQuartet.AC[0]+
					      Libr12.ShellQuartet.AB[1]*Libr12.ShellQuartet.AC[1]+
					      Libr12.ShellQuartet.AB[2]*Libr12.ShellQuartet.AC[2];
		Libr12.ShellQuartet.CDdotCA = -1.0*(Libr12.ShellQuartet.CD[0]*Libr12.ShellQuartet.AC[0]+
						    Libr12.ShellQuartet.CD[1]*Libr12.ShellQuartet.AC[1]+
						    Libr12.ShellQuartet.CD[2]*Libr12.ShellQuartet.AC[2]);
		AB2 = Libr12.ShellQuartet.AB[0]*Libr12.ShellQuartet.AB[0]+
		      Libr12.ShellQuartet.AB[1]*Libr12.ShellQuartet.AB[1]+
		      Libr12.ShellQuartet.AB[2]*Libr12.ShellQuartet.AB[2];
		CD2 = Libr12.ShellQuartet.CD[0]*Libr12.ShellQuartet.CD[0]+
		      Libr12.ShellQuartet.CD[1]*Libr12.ShellQuartet.CD[1]+
		      Libr12.ShellQuartet.CD[2]*Libr12.ShellQuartet.CD[2];

		/*--------------------------------
		  contract by primitives out here
		 --------------------------------*/
		num_prim_comb = 0;
		for (pi = 0; pi < np_i; pi++)
		  for (pj = 0; pj < np_j; pj++)
		    for (pk = 0; pk < np_k; pk++)
		      for (pl = 0; pl < np_l; pl++){
			r12_quartet_data(&(Libr12.PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
					 sp_ij, sp_kl, am, pi, pj, pk, pl, lambda_T);
		      }

		if (am) {
		  build_r12_grt[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libr12,num_prim_comb);
		  if (swap_ij_kl)
		    /*--- (usi usj|[r12,T1]|usk usl) = (usk usl|[r12,T2]|usi usj) ---*/
		    Libr12.te_ptr[2] = Libr12.te_ptr[3];
#ifdef NONDOUBLE_INTS
		  size = ioff[BasisSet.shells[si].am]*ioff[BasisSet.shells[sj].am]*
			 ioff[BasisSet.shells[sk].am]*ioff[BasisSet.shells[sl].am];
		  for(i=0;i<NUM_TE_TYPES-1;i++)
		    for(j=0;j<size;j++)
		      raw_data[i][j] = (double) Libr12.te_ptr[i][j];
#else
		  for(i=0;i<NUM_TE_TYPES-1;i++)
		    raw_data[i] = Libr12.te_ptr[i];
#endif
		  /*--- Just normalize the integrals ---*/
		  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++)
		    data[te_type] = norm_quartet(raw_data[te_type], NULL, orig_am, 0);
		}
		else {
		  ssss = 0.0;
		  ss_r12_ss = 0.0;
		  for(p=0;p<num_prim_comb;p++) {
		    ssss += (double) Libr12.PrimQuartet[p].F[0];
		    ss_r12_ss += (double) Libr12.PrimQuartet[p].ss_r12_ss;
		  }
		  build_r12_grt[0][0][0][0](&Libr12,num_prim_comb);
#ifdef NONDOUBLE_INTS
		  raw_data[0][0] = ssss;
		  raw_data[1][0] = ss_r12_ss;
		  raw_data[2][0] = (double) Libr12.te_ptr[2][0];
		  raw_data[3][0] = (double) Libr12.te_ptr[3][0];
		  data[0] = raw_data[0];
		  data[1] = raw_data[1];
		  data[2] = raw_data[2];
		  data[3] = raw_data[3];
#else
		  Libr12.int_stack[2] = Libr12.te_ptr[2][0];
		  Libr12.int_stack[3] = Libr12.te_ptr[3][0];
		  Libr12.int_stack[0] = ssss;
		  Libr12.int_stack[1] = ss_r12_ss;
		  data[0] = Libr12.int_stack;
		  data[1] = Libr12.int_stack+1;
		  data[2] = Libr12.int_stack+2;
		  data[3] = Libr12.int_stack+3;
#endif
		}

		/*---
		  Swap bra and ket back to the reversed order (PQ|RS) if not done yet
		  Need this to make Step 1 a (fast) vector operation!
		  NOTE: This changes only the way integrals are stored! No need to
		  worry about non-hermiticity of [r12,T1] here!!!
		  ---*/
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
		  /*---
		    No need to swap bra and ket since the quartet
		    was already computed in reverse order
		   ---*/
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
		  /*---
		    Need to swap bra and ket, but do
		    the actual permutation in te_type-loop
		   ---*/
		}

		/*--- step 1 of the transformation ---*/
/*		timer_on("Step 1");*/
		for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
		  /*---
		    if bra and ket were not swapped before computing
		    the quartet, swap them now so that Step 1 is
		    a vector operation
		   ---*/
		  if ((!swap_ij_kl) && am) {
/*		    timer_on("Pre1Swap");*/
		    /*--- (rs|pq) -> (pq|rs) ---*/
		    ijkl_to_klij(data[te_type],scratch_buf,nr*ns,np*nq);
		    data[te_type] = scratch_buf;
/*		    timer_off("Pre1Swap");*/
		  }
		  if (usk != usl)
		    for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		      mo_vec = MOInfo.scf_evec_occ_alpha[mo_i+imin];
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			for(q=0,q_abs=sq_fao;
			    q<nq;
			    q++,q_abs++) {
			  ip_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + p_abs) * nr * ns);
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
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
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
		      mo_vec = MOInfo.scf_evec_occ_alpha[mo_i+imin];
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
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
		}
/*		timer_off("Step 1");*/
	      
	      } /* end of computing "petit" list - end of P,Q loop */
	    } /* end of "unique" P,Q loop */
	  
	  /*--- step 2 of the transfromation ---*/
	  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
/*	    timer_on("Post1Swap");*/
	    ijkl_to_klij(rsiq_buf[te_type],scratch_buf,ibatch_length*BasisSet.num_ao,nr*ns);
	    rsi_row = scratch_buf;
/*	    timer_off("Post1Swap");*/
	    ix_block_ptr = rsix_buf[te_type];
	    for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
		  /*--- Can be done as a matrix multiply now ---*/
/*		  timer_on("Step 2");*/
#if !USE_BLAS
		  ix = 0;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=BasisSet.num_ao) {
		    for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++,ix++) {
		      mo_vec = MOInfo.scf_evec_alpha[mo_x];
		      temp = ix_block_ptr[ix];
		      for(q_abs=0;q_abs<BasisSet.num_ao;q_abs++) {
			temp += mo_vec[q_abs] * rsi_row[q_abs];
		      }
		      ix_block_ptr[ix] = temp;
		    }
		  }
#else
		  C_DGEMM('n','t',ibatch_length,MOInfo.num_mo,BasisSet.num_ao,1.0,
			  rsi_row,BasisSet.num_ao,MOInfo.scf_evec_alpha[0],BasisSet.num_ao,
			  0.0,ix_block_ptr,MOInfo.num_mo);
		  rsi_row += BasisSet.num_ao*ibatch_length;
#endif
		  ix_block_ptr += MOInfo.num_mo*ibatch_length;
/*		  timer_off("Step 2");*/
	      }
	    }

	    /*--- step 3 of the transformation (AA|AA) spin case---*/
	    rsi_row = rsix_buf[te_type];
	    /*--- To update (JS|IA) need to lock mutex corresponding to the S and R indices ---*/
#if LOCK_RS_SHELL	  
	    pthread_mutex_lock(&rmp2r12_sindex_mutex[INDEX(si,sj)]);
#endif
	    for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
/*		  timer_on("Step 3");*/
		  r_abs = r + sr_fao;
		  s_abs = s + ss_fao;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=MOInfo.num_mo) {
		    for(mo_j=0;mo_j<(te_type==2 ? MOInfo.alpha_act_occ : mo_i);mo_j++) {
#if !LOCK_RS_SHELL
		      pthread_mutex_lock(&rmp2r12_sindex_mutex[s_abs]);
#endif
		      jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
		      temp = MOInfo.scf_evec_occ_alpha[mo_j+jmin][r_abs];
		      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			jsi_row[mo_x] += temp * rsi_row[mo_x];
		      }
#if !LOCK_RS_SHELL
		      pthread_mutex_unlock(&rmp2r12_sindex_mutex[s_abs]);
#endif
		      if (usi != usj) {
#if !LOCK_RS_SHELL
			pthread_mutex_lock(&rmp2r12_sindex_mutex[r_abs]);
#endif
			jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + r_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
			temp = MOInfo.scf_evec_occ_alpha[mo_j+jmin][s_abs];
			if (te_type == 2)   /*--- [r12,T1] integral - negative sign ---*/
			  temp *= (-1.0);
			for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			  jsi_row[mo_x] += temp * rsi_row[mo_x];
			}
#if !LOCK_RS_SHELL
			pthread_mutex_unlock(&rmp2r12_sindex_mutex[r_abs]);
#endif
		      }
		    }
		  }
/*		  timer_off("Step 3");*/
	      }
	    }
#if LOCK_RS_SHELL
	    pthread_mutex_unlock(&rmp2r12_sindex_mutex[INDEX(si,sj)]);
#endif


	    /*--- step 3 of the transformation (AA|BB) spin case---*/
	    rsi_row = rsix_buf[te_type];
	    /*--- To update (JS|IA) need to lock mutex corresponding to the S and R indices ---*/
#if LOCK_RS_SHELL	  
	    pthread_mutex_lock(&rmp2r12_sindex_mutex[INDEX(si,sj)]);
#endif
	    for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
/*		  timer_on("Step 3");*/
		  r_abs = r + sr_fao;
		  s_abs = s + ss_fao;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=MOInfo.num_mo) {
		    for(mo_j=0;mo_j<MOInfo.beta_act_occ;mo_j++) {
#if !LOCK_RS_SHELL
		      pthread_mutex_lock(&rmp2r12_sindex_mutex[s_abs]);
#endif
		      jsi_row = jsix_buf2[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
		      temp = MOInfo.scf_evec_occ_beta[mo_j+jmin][r_abs];
		      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			jsi_row[mo_x] += temp * rsi_row[mo_x];
		      }
#if !LOCK_RS_SHELL
		      pthread_mutex_unlock(&rmp2r12_sindex_mutex[s_abs]);
#endif
		      if (usi != usj) {
#if !LOCK_RS_SHELL
			pthread_mutex_lock(&rmp2r12_sindex_mutex[r_abs]);
#endif
			jsi_row = jsix_buf2[te_type] + ((mo_j * BasisSet.num_ao + r_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
			temp = MOInfo.scf_evec_occ_beta[mo_j+jmin][s_abs];
			if (te_type == 2)   /*--- [r12,T1] integral - negative sign ---*/
			  temp *= (-1.0);
			for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			  jsi_row[mo_x] += temp * rsi_row[mo_x];
			}
#if !LOCK_RS_SHELL
			pthread_mutex_unlock(&rmp2r12_sindex_mutex[r_abs]);
#endif
		      }
		    }
		  }
/*		  timer_off("Step 3");*/
	      }
	    }
#if LOCK_RS_SHELL
	    pthread_mutex_unlock(&rmp2r12_sindex_mutex[INDEX(si,sj)]);
#endif


	    memset(rsiq_buf[te_type],0,nr*ns*ibatch_length*BasisSet.num_ao*sizeof(double));
	    memset(rsix_buf[te_type],0,nr*ns*ibatch_length*MOInfo.num_mo*sizeof(double));
	  } /* end loop over te_type */
	  } /* end of R,S loop */
      } /* end of "unique" R,S loop */


    pthread_mutex_lock(&rmp2r12_energy_mutex);
    RMP2R12_Status.num_arrived++;
    if (RMP2R12_Status.num_arrived != UserOptions.num_threads) { /*--- there are some threads still working - wait here --*/
      pthread_cond_wait(&rmp2r12_energy_cond,&rmp2r12_energy_mutex);
      pthread_mutex_unlock(&rmp2r12_energy_mutex);
    }
    else { /*--- this is the last thread to get here - do the 4th step and dumping integrals alone and wake everybody up ---*/
    /*--- step 4 of the transformation ---*/
/*    timer_on("Step 4");*/
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
      /* make the alpha - alpha integrals */
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<(te_type==2 ? MOInfo.alpha_act_occ : mo_i+1);mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    jyi_row = jyix_buf[te_type] + ((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo;
	    for(s_abs=0;s_abs<BasisSet.num_ao;s_abs++) {
	      jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
	      temp = MOInfo.scf_evec_alpha[mo_y][s_abs];
	      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
		jyi_row[mo_x] += temp * jsi_row[mo_x];
	      }
	    }
	  }
	}
      }
      /* make the alpha - beta integrals */
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<MOInfo.beta_act_occ;mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    jyi_row = jyix_buf2[te_type] + ((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo;
	    for(s_abs=0;s_abs<BasisSet.num_ao;s_abs++) {
	      jsi_row = jsix_buf2[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
	      temp = MOInfo.scf_evec_beta[mo_y][s_abs];
	      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
		jyi_row[mo_x] += temp * jsi_row[mo_x];
	      }
	    }
	  }
	}
      }
    }
/*    timer_off("Step 4");*/

#if PRINT
    /*--- Print them out for now ---*/
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
      fprintf(outfile,"  Alpha Alpha Transformed integrals of type %d\n",te_type);
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<(te_type==2 ? MOInfo.alpha_act_occ : mo_i+1);mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
	      temp = jyix_buf[te_type][((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo + mo_x];
	      if (fabs(temp) > ZERO) {
		if (te_type < 2)
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_j+jmin,mo_y,mo_i+imin,mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			temp);
		else
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_i+imin,mo_x,mo_j+jmin,mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			temp);
	      }
	    }
	  }
	}
      }
      fprintf(outfile,"  Alpha Beta Transformed integrals of type %d\n",te_type);
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<MOInfo.beta_act_occ;mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
	      temp = jyix_buf2[te_type][((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo + mo_x];
	      if (fabs(temp) > ZERO) {
		if (te_type < 2)
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_j+jmin,mo_y,mo_i+imin,mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			temp);
		else
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_i+imin,mo_x,mo_j+jmin,mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			temp);
	      }
	    }
	  }
	}
      }
    }
#endif

    /*--- Files are opened and closed each pass to ensure integrity of TOCs
      if restart ever needed, these are the alpha-beta integrals ---*/
    iounits[0] = PSIF_MO_AA_TEI;
    iounits[1] = PSIF_MO_AA_R12;
    iounits[2] = PSIF_MO_AA_R12T1;
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++)
      psio_open(iounits[te_type], (ibatch != 0) ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);

    /*--------------------------------------------------------
      Dump all fully transformed integrals to disk. Zero out
      non-symmetrical integrals (what's left of the Pitzer's
      equal contribution theorem in Abelian case).
     --------------------------------------------------------*/
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	/*--------------------------------------------------------------------
	  Write integrals out in num_mo by num_mo batches corresponding to
	  each active ij pair.
	 --------------------------------------------------------------------*/
	for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	  i = act2fullQTS_alpha[mo_i+imin]; /*--- mo_i+imin is the index in QTS "active" indexing scheme
					          we need i to be in QTS "full" indexing scheme, that's
						  what the MP2R12 code expects ---*/
	  isym = MOInfo.mo2symblk_occ_alpha[mo_i+imin];
	  for(mo_j=0;mo_j<(te_type==2 ? MOInfo.alpha_act_occ : mo_i+1);mo_j++) {
	    jsym = MOInfo.mo2symblk_occ_alpha[mo_j+jmin];
	    j = act2fullQTS_alpha[mo_j+jmin];    /*--- Again, get the "full" QTS index ---*/
	    sprintf(ij_key_string,"Block_%d_x_%d_y",i,j);
	    memset(xy_buf,0,MOInfo.num_mo*MOInfo.num_mo*sizeof(double));
	    /*--- Put all integrals with common i and j into a buffer ---*/
	    for(mo_x=0,xy=0;mo_x<MOInfo.num_mo;mo_x++) {
	      x = mo_x;    /*--- The second index is a Pitzer index, that's fine ---*/
	      xsym = MOInfo.mo2symblk[x];
	      for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++,xy++) {
		y = mo_y;                    /*--- Again, the Pitzer index here is what we need ---*/
		ysym = MOInfo.mo2symblk[y];
		/*--- Skip this indegral if it's non-totally symmetric -
		  Pitzer's contribution theorem in Abelian case ---*/
		if ((isym ^ jsym) ^ (xsym ^ ysym))
		  continue;
		/*--- find the integral in ixjy_buf and put it in xy_buf ---*/
		m = (((mo_j*MOInfo.num_mo + mo_y)*ibatch_length + mo_i)*MOInfo.num_mo + mo_x);
		xy_buf[xy] = jyix_buf[te_type][m];
	      }
	    }
	    psio_write_entry(iounits[te_type], ij_key_string, (char *)xy_buf,
			     MOInfo.num_mo*MOInfo.num_mo*sizeof(double));
	  }
	}
	psio_close(iounits[te_type], 1);
    }

    /*--- Exactly like the above, but now we're writing out the alpha - beta integrals ---*/
    iounits2[0] = PSIF_MO_AB_TEI;
    iounits2[1] = PSIF_MO_AB_R12;
    iounits2[2] = PSIF_MO_AB_R12T1;
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++)
      psio_open(iounits2[te_type], (ibatch != 0) ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);

    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	/*--------------------------------------------------------------------
	  Write integrals out in num_mo by num_mo batches corresponding to
	  each active ij pair.
	 --------------------------------------------------------------------*/
	for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	  i = act2fullQTS_alpha[mo_i+imin]; /*--- mo_i+imin is the index in QTS "active" indexing scheme
					          we need i to be in QTS "full" indexing scheme, that's
						  what the MP2R12 code expects ---*/
	  isym = MOInfo.mo2symblk_occ_alpha[mo_i+imin];
	  for(mo_j=0;mo_j<MOInfo.beta_act_occ;mo_j++) {
	    jsym = MOInfo.mo2symblk_occ_beta[mo_j+jmin];
	    j = act2fullQTS_beta[mo_j+jmin];    /*--- Again, get the "full" QTS index ---*/
	    sprintf(ij_key_string,"Block_%d_x_%d_y",i,j);
	    memset(xy_buf,0,MOInfo.num_mo*MOInfo.num_mo*sizeof(double));
	    /*--- Put all integrals with common i and j into a buffer ---*/
	    for(mo_x=0,xy=0;mo_x<MOInfo.num_mo;mo_x++) {
	      x = mo_x;    /*--- The second index is a Pitzer index, that's fine ---*/
	      xsym = MOInfo.mo2symblk[x];
	      for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++,xy++) {
		y = mo_y;                    /*--- Again, the Pitzer index here is what we need ---*/
		ysym = MOInfo.mo2symblk[y];
		/*--- Skip this indegral if it's non-totally symmetric -
		  Pitzer's contribution theorem in Abelian case ---*/
		if ((isym ^ jsym) ^ (xsym ^ ysym))
		  continue;
		/*--- find the integral in ixjy_buf and put it in xy_buf ---*/
		m = (((mo_j*MOInfo.num_mo + mo_y)*ibatch_length + mo_i)*MOInfo.num_mo + mo_x);
		xy_buf[xy] = jyix_buf2[te_type][m];
	      }
	    }
	    psio_write_entry(iounits2[te_type], ij_key_string, (char *)xy_buf,
			     MOInfo.num_mo*MOInfo.num_mo*sizeof(double));
	  }
	}
	psio_close(iounits2[te_type], 1);
    }


    if (ibatch < num_ibatch-1)
      for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	memset(jsix_buf[te_type],0,MOInfo.alpha_act_occ*BasisSet.num_ao*ibatch_length*MOInfo.num_mo*sizeof(double));
	memset(jyix_buf[te_type],0,MOInfo.alpha_act_occ*MOInfo.num_mo*ibatch_length*MOInfo.num_mo*sizeof(double));
	memset(jsix_buf2[te_type],0,MOInfo.beta_act_occ*BasisSet.num_ao*ibatch_length*MOInfo.num_mo*sizeof(double));
	memset(jyix_buf2[te_type],0,MOInfo.beta_act_occ*MOInfo.num_mo*ibatch_length*MOInfo.num_mo*sizeof(double));
      }
    /*--- Done with the non-threaded part - wake everybody up and prepare for the next ibatch ---*/
    RMP2R12_Status.num_arrived = 0;
    pthread_cond_broadcast(&rmp2r12_energy_cond);
    pthread_mutex_unlock(&rmp2r12_energy_mutex);
    }
  } /* End of "I"-loop */

  /*---------
    Clean-up
   ---------*/
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    free(rsiq_buf[te_type]);
    free(rsix_buf[te_type]);
  }
  free(xy_buf);
  free_block(jy_buf);
  free(scratch_buf);
  free_libr12(&Libr12);
  free_fjt_table(&fjt_table);
#ifdef NONDOUBLE_INTS
  for(i=0;i<NUM_TE_TYPES;i++)
    free(raw_data[i]);
#endif
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);

  return NULL;
}

}}} /* End namespaces */
