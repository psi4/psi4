/*! \file te_deriv1_scf_thread.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
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
void *te_deriv1_scf_thread(void *tnum_ptr)
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
  register int i, j, k, l, m, ii, jj, kk, ll;
  register int si, sj, sk, sl ;
  register int sii, sjj, skk, sll, slll;
  register int pi, pj, pk, pl ;
  int max_pj, max_pl;
  register int pii, pjj, pkk, pll ;
  int switch_ij, switch_kl, switch_ijkl;
  int center_i, center_j, center_k, center_l;

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

  grad_te_local = block_matrix(Molecule.num_atoms,3);

/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  for (sii=0; sii<BasisSet.num_shells; sii++)
    for (sjj=0; sjj<=sii; sjj++)
      for (skk=0; skk<=sii; skk++)
	for (sll=0; sll<= ((sii == skk) ? sjj : skk); sll++, quartet_index++){

	    si = sii; sj = sjj; sk = skk; sl = sll;

	    /*--- Decide if this thread will do this ---*/
	    if ( quartet_index%UserOptions.num_threads != thread_num )
	      continue;

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
	    if (si == sj)
	      pfac *= 0.5;
	    if (sk == sl)
	      pfac *= 0.5;
	    if (si == sk && sj == sl || si == sl && sj == sk)
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
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
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
	}

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
