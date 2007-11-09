/*! \file te_deriv1_scf_thread_symm.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory.h>
#include <pthread.h>
#include <psitypes.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include "deriv1_quartet_data.h"
#include "small_fns.h"

namespace psi { namespace CINTS {
void *te_deriv1_scf_thread_symm(void *tnum_ptr)
{
  const long int thread_num = (long int) tnum_ptr;

  extern double **grad_te;
  extern pthread_mutex_t deriv1_mutex;    /* Used to lock "global" gradient matrix grad_te during update */
  
  /*--- Various data structures ---*/
  struct shell_pair *sp_ij, *sp_kl;
  Libderiv_t Libderiv;                    /* Integrals library object */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;               /* table of auxiliary function F_m(u) for each primitive combination */
#endif

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
  register int pii, pjj, pkk, pll ;
  int switch_ij, switch_kl, switch_ijkl;
  int center_i, center_j, center_k, center_l;

  int *sj_arr, *sk_arr, *sl_arr;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int max_num_unique_quartets;
  int plquartet;

  int class_size;
  int max_class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, quartet_size;

  int si_fao, sj_fao, sk_fao, sl_fao;
  int sii_fao, sjj_fao, skk_fao, sll_fao;
  int ao_i, imax, ao_j, jmax, ao_k, kmax, ao_l, lmax;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p, max_num_prim_comb;

  int buf_offset, buf_4offset, buf_size, last_buf;
  int quartet_done, offset;

  int mosh_i, mosh_j;

  double AB2, CD2;
  double *FourInd;
  double **grad_te_local;
  double pfac;
  double temp;
  double alpha, beta;
  double **dens_i, **dens_j;
  double ddax, dday, ddaz, ddbx, ddby, ddbz,
    ddcx, ddcy, ddcz, dddx, dddy, dddz;

  /*---------------
    Initialization
    ---------------*/
#ifndef USE_TAYLOR_FM
  init_fjt_table(&fjt_table);
#endif

  max_cart_class_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  max_class_size = max_cart_class_size;
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
    (BasisSet.max_num_prims*BasisSet.max_num_prims);
  init_libderiv1(&Libderiv,BasisSet.max_am-1,max_num_prim_comb,max_class_size);
  FourInd = init_array(max_cart_class_size);
  max_num_unique_quartets = Symmetry.max_stab_index*
    Symmetry.max_stab_index*
    Symmetry.max_stab_index;
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);

  grad_te_local = block_matrix(Molecule.num_atoms,3);

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
	    class_size = ni*nj*nk*nl;

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
	      }
	    } /* petite list is ready to be used */
	    num_unique_quartets = count;
	  }
	  else { /*--- C1 symmetry case ---*/
	    num_unique_quartets = 1;
	    sj_arr[0] = usj;
	    sk_arr[0] = usk;
	    sl_arr[0] = usl;
	    ni = (BasisSet.puream ? 2*BasisSet.shells[usi].am - 1 : ioff[BasisSet.shells[usi].am]);
	    nj = (BasisSet.puream ? 2*BasisSet.shells[usj].am - 1 : ioff[BasisSet.shells[usj].am]);
	    nk = (BasisSet.puream ? 2*BasisSet.shells[usk].am - 1 : ioff[BasisSet.shells[usk].am]);
	    nl = (BasisSet.puream ? 2*BasisSet.shells[usl].am - 1 : ioff[BasisSet.shells[usl].am]);
	    ioffset = BasisSet.shells[usi].fbf - 1;
	    joffset = BasisSet.shells[usj].fbf - 1;
	    koffset = BasisSet.shells[usk].fbf - 1;
	    loffset = BasisSet.shells[usl].fbf - 1;
	  }


	  /*----------------------------------
	    Compute the nonredundant quartets
	    ----------------------------------*/
	  for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	    si = Symmetry.us2s[usii];
	    sj = sj_arr[plquartet];
	    sk = sk_arr[plquartet];
	    sl = sl_arr[plquartet];

	    /*--- Skip this quartet if all four centers are the same ---*/
	    if (BasisSet.shells[si].center == BasisSet.shells[sj].center &&
		BasisSet.shells[si].center == BasisSet.shells[sk].center &&
		BasisSet.shells[si].center == BasisSet.shells[sl].center &&
		DERIV_LVL == 1)
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
	    
	    AB2 = Libderiv.AB[0]*Libderiv.AB[0]+
	      Libderiv.AB[1]*Libderiv.AB[1]+
	      Libderiv.AB[2]*Libderiv.AB[2];
	    CD2 = Libderiv.CD[0]*Libderiv.CD[0]+
	      Libderiv.CD[1]*Libderiv.CD[1]+
	      Libderiv.CD[2]*Libderiv.CD[2];

	    /*-------------------------
	      Figure out the prefactor
	      -------------------------*/
	    pfac = 1.0;
	    if (usi == usj)
	      pfac *= 0.5;
	    if (usk == usl)
	      pfac *= 0.5;
	    if (usi == usk && usj == usl || usi == usl && usj == usk)
	      pfac *= 0.5;

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
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					NULL, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac*lambda_T);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac*lambda_T);
#endif		    
		  }
		}
	      }
	    }

	    /*-------------
	      Form FourInd
	      -------------*/
	    si_fao = BasisSet.shells[si].fao-1;
	    sj_fao = BasisSet.shells[sj].fao-1;
	    sk_fao = BasisSet.shells[sk].fao-1;
	    sl_fao = BasisSet.shells[sl].fao-1;
	    if (UserOptions.reftype == rhf || UserOptions.reftype == uhf) {           /*--- RHF or UHF ---*/
	      count = 0;
	      for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
		for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
		  for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
		    for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
		      FourInd[count] = (4.0*Dens[ao_i][ao_j]*Dens[ao_k][ao_l] -
					Dens[ao_i][ao_k]*Dens[ao_j][ao_l] -
					Dens[ao_i][ao_l]*Dens[ao_k][ao_j])*
			GTOs.bf_norm[orig_am[0]][ao_i-si_fao]*
			GTOs.bf_norm[orig_am[1]][ao_j-sj_fao]*
			GTOs.bf_norm[orig_am[2]][ao_k-sk_fao]*
			GTOs.bf_norm[orig_am[3]][ao_l-sl_fao];
		      count++;
		    }
	    }
	    else {                     /*--- ROHF or TCSCF ---*/
	      if (am)
		memset((char *) FourInd, 0, sizeof(double)*quartet_size);
	      else
		FourInd[0] = 0.0;
	      for(mosh_i=0;mosh_i<MOInfo.num_moshells;mosh_i++)
		for(mosh_j=0;mosh_j<MOInfo.num_moshells;mosh_j++) {
		  count = 0;
		  dens_i = ShDens[mosh_i];
		  dens_j = ShDens[mosh_j];
		  alpha = 8.0*MOInfo.Alpha[mosh_i][mosh_j];
		  beta = 4.0*MOInfo.Beta[mosh_i][mosh_j];
		  for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
		    for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
		      for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
			for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
			  FourInd[count] += (alpha*dens_i[ao_i][ao_j]*dens_j[ao_k][ao_l] +
					     beta*(dens_i[ao_i][ao_k]*dens_j[ao_j][ao_l] +
						   dens_i[ao_i][ao_l]*dens_j[ao_k][ao_j]));
			  count++;
			}
		}
	      /*--- Normalize it ---*/
	      count = 0;
	      for (ao_i = 0; ao_i < ni; ao_i++)
		for (ao_j = 0; ao_j < nj; ao_j++)
		  for (ao_k = 0; ao_k < nk; ao_k++)
		    for (ao_l = 0; ao_l < nl; ao_l++) {
		      FourInd[count] *= GTOs.bf_norm[orig_am[0]][ao_i]*
			GTOs.bf_norm[orig_am[1]][ao_j]*
			GTOs.bf_norm[orig_am[2]][ao_k]*
			GTOs.bf_norm[orig_am[3]][ao_l];
		      count++;
		    }
	    }
	      
	    build_deriv1_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libderiv,num_prim_comb);

	    center_i = BasisSet.shells[si].center-1;
	    center_j = BasisSet.shells[sj].center-1;
	    center_k = BasisSet.shells[sk].center-1;
	    center_l = BasisSet.shells[sl].center-1;
	    
	    ddax = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddax += Libderiv.ABCD[0][k]*FourInd[k];
	    grad_te_local[center_i][0] += ddax;

	    dday = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dday += Libderiv.ABCD[1][k]*FourInd[k];
	    grad_te_local[center_i][1] += dday;

	    ddaz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddaz += Libderiv.ABCD[2][k]*FourInd[k];
	    grad_te_local[center_i][2] += ddaz;

	    /*	    ddbx = 0.0;
		    for(k=0;k<quartet_size;k++)
		    ddbx += Libderiv.ABCD[3][k]*FourInd[k];
		    grad_te_local[center_j][0] += ddbx;

		    ddby = 0.0;
		    for(k=0;k<quartet_size;k++)
		    ddby += Libderiv.ABCD[4][k]*FourInd[k];
		    grad_te_local[center_j][1] += ddby;

		    ddbz = 0.0;
		    for(k=0;k<quartet_size;k++)
		    ddbz += Libderiv.ABCD[5][k]*FourInd[k];
		    grad_te_local[center_j][2] += ddbz;*/

	    ddcx = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddcx += Libderiv.ABCD[6][k]*FourInd[k];
	    grad_te_local[center_k][0] += ddcx;

	    ddcy = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddcy += Libderiv.ABCD[7][k]*FourInd[k];
	    grad_te_local[center_k][1] += ddcy;

	    ddcz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddcz += Libderiv.ABCD[8][k]*FourInd[k];
	    grad_te_local[center_k][2] += ddcz;

	    dddx = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddx += Libderiv.ABCD[9][k]*FourInd[k];
	    grad_te_local[center_l][0] += dddx;

	    dddy = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddy += Libderiv.ABCD[10][k]*FourInd[k];
	    grad_te_local[center_l][1] += dddy;

	    dddz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddz += Libderiv.ABCD[11][k]*FourInd[k];
	    grad_te_local[center_l][2] += dddz;
	      
	    grad_te_local[center_j][0] -= ddax + ddcx + dddx;
	    grad_te_local[center_j][1] -= dday + ddcy + dddy;
	    grad_te_local[center_j][2] -= ddaz + ddcz + dddz;

	  } /* end of loop over symmetry unique quartets */
	} /* end of unique shell loops */

  pthread_mutex_lock(&deriv1_mutex);
  for(i=0;i<Molecule.num_atoms;i++) {
    grad_te[i][0] += grad_te_local[i][0];
    grad_te[i][1] += grad_te_local[i][1];
    grad_te[i][2] += grad_te_local[i][2];
  }
  pthread_mutex_unlock(&deriv1_mutex);
  
  /*---------
    Clean-up
    ---------*/
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  free(FourInd);

  free_block(grad_te_local);
  free_libderiv(&Libderiv);
#ifndef USE_TAYLOR_FM
  free_fjt_table(&fjt_table);
#endif

  return NULL;
}
};
}
