/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory.h>
#include<stdexcept>
#include <psitypes.h>
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
#include <Tools/prints.h>

namespace psi {
  namespace CINTS {
void te_deriv1_print(void)
{
  Libderiv_t Libderiv;                    /* Integrals library object */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;               /* table of auxiliary function F_m(u) for each primitive combination */
#endif

  int max_class_size;
  int max_cart_class_size;
  int max_num_prim_comb, num_prim_comb;
  int max_num_unique_quartets;
  int *sj_arr, *sk_arr, *sl_arr;
  int usii, usjj, uskk, usll;
  PSI_INT_LEAST64 quartet_index;
  int sii, sjj, skk, sll;
  int num_unique_quartets;
  int ni, nj, nk, nl;
  int si, sj, sk, sl;
  int plquartet, dum;
  int switch_ij, switch_kl, switch_ijkl;
  int np_i, np_j, np_k, np_l;
  int quartet_size;
  int orig_am[4], am;
  struct shell_pair *sp_ij, *sp_kl;
  double AB2, CD2;
  int pi, pj, pk, pl;
  int si_fao, sj_fao, sk_fao, sl_fao;
  int ao_i, ao_j, ao_k, ao_l;
  int count;
  int center_i, center_j, center_k, center_l;
  int iout, jout, kout, lout;
  int i, k, ij, kl;

  double *PrintFourInd;
  double **derivs;
  FILE **out;
  char lbl[20];

  /* some error checking up front */
  if(Symmetry.nirreps > 1 || BasisSet.puream) {
    fprintf(outfile, "\nThe derivative integral printing code will almost certainly fail\n");
    fprintf(outfile, "for non-C1 and/or puream cases.  Exiting.\n");
    throw std::domain_error("no symmetry or puream cases allowed when printing deriv. ints.");
  }

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
  PrintFourInd = init_array(max_cart_class_size);
  derivs = block_matrix(Molecule.num_atoms*3,max_cart_class_size);
  out = (FILE **) malloc(Molecule.num_atoms * 3 * sizeof(FILE *));

  for(i=0; i < Molecule.num_atoms*3; i++) {
    sprintf(lbl,"eri%d.dat", i);
    ffile(&out[i], lbl, 0);
  }

  max_num_unique_quartets = Symmetry.max_stab_index*
    Symmetry.max_stab_index*
    Symmetry.max_stab_index;
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);

  for (usii=0; usii<Symmetry.num_unique_shells; usii++)
    for (usjj=0; usjj<Symmetry.num_unique_shells; usjj++)
      for (uskk=0; uskk<Symmetry.num_unique_shells; uskk++)
	for (usll=0; usll<Symmetry.num_unique_shells; usll++, quartet_index++){

	  sii = Symmetry.us2s[usii];
	  sjj = Symmetry.us2s[usjj];
	  skk = Symmetry.us2s[uskk];
	  sll = Symmetry.us2s[usll];

	  num_unique_quartets = 1;
	  sj_arr[0] = usjj;
	  sk_arr[0] = uskk;
	  sl_arr[0] = usll;

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

	    num_prim_comb = 0;
	    for (pi = 0; pi < np_i; pi++) {
	      for (pj = 0; pj < np_j; pj++) {
		for (pk = 0; pk < np_k; pk++) {
		  for (pl = 0; pl < np_l; pl++){
#ifdef USE_TAYLOR_FM
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					NULL, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, 1);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, 1);
#endif		    
		  }
		}
	      }
	    }

	    si_fao = BasisSet.shells[si].fao-1;
	    sj_fao = BasisSet.shells[sj].fao-1;
	    sk_fao = BasisSet.shells[sk].fao-1;
	    sl_fao = BasisSet.shells[sl].fao-1;
	    count = 0;
	    for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
	      for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
		for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
		  for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
		    PrintFourInd[count] = GTOs.bf_norm[orig_am[0]][ao_i-si_fao] *
		      GTOs.bf_norm[orig_am[1]][ao_j-sj_fao] *
		      GTOs.bf_norm[orig_am[2]][ao_k-sk_fao] *
		      GTOs.bf_norm[orig_am[3]][ao_l-sl_fao];
		    count++;
		  }

	    build_deriv1_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libderiv,num_prim_comb);

	    center_i = BasisSet.shells[si].center-1;
	    center_j = BasisSet.shells[sj].center-1;
	    center_k = BasisSet.shells[sk].center-1;
	    center_l = BasisSet.shells[sl].center-1;

	    zero_mat(derivs, Molecule.num_atoms*3, max_cart_class_size);
	    for(k=0; k < quartet_size; k++) {
	      derivs[center_i*3+0][k] += Libderiv.ABCD[0][k] * PrintFourInd[k];
	      derivs[center_i*3+1][k] += Libderiv.ABCD[1][k] * PrintFourInd[k];
	      derivs[center_i*3+2][k] += Libderiv.ABCD[2][k] * PrintFourInd[k];

	      derivs[center_j*3+0][k] -= (Libderiv.ABCD[0][k] + Libderiv.ABCD[6][k] + Libderiv.ABCD[9][k]) * PrintFourInd[k];
	      derivs[center_j*3+1][k] -= (Libderiv.ABCD[1][k] + Libderiv.ABCD[7][k] + Libderiv.ABCD[10][k]) * PrintFourInd[k];
	      derivs[center_j*3+2][k] -= (Libderiv.ABCD[2][k] + Libderiv.ABCD[8][k] + Libderiv.ABCD[11][k]) * PrintFourInd[k];

	      derivs[center_k*3+0][k] += Libderiv.ABCD[6][k] * PrintFourInd[k];
	      derivs[center_k*3+1][k] += Libderiv.ABCD[7][k] * PrintFourInd[k];
	      derivs[center_k*3+2][k] += Libderiv.ABCD[8][k] * PrintFourInd[k];

	      derivs[center_l*3+0][k] += Libderiv.ABCD[9][k] * PrintFourInd[k];
	      derivs[center_l*3+1][k] += Libderiv.ABCD[10][k] * PrintFourInd[k];
	      derivs[center_l*3+2][k] += Libderiv.ABCD[11][k] * PrintFourInd[k];
	    }

	    for(k=0; k < Molecule.num_atoms*3; k++) {
	      count = 0;
	      for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
		for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
		  for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
		    for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
		      if(fabs(derivs[k][count]) > 1e-6) {
			iout = ao_i;  jout = ao_j; kout = ao_k; lout = ao_l;
			if(switch_ijkl) { dum = iout; iout = kout; kout = dum; dum = jout; jout = lout; lout = dum; }
			if(switch_kl) { dum = kout; kout = lout; lout = dum; }
			if(switch_ij) { dum = iout; iout = jout; jout = dum; }
			ij = INDEX(iout,jout); kl = INDEX(kout,lout);
			if(iout >= jout && kout >= lout && ij >= kl)
			  fprintf(out[k], "%3d %3d %3d %3d %20.14f\n", iout, jout, kout, lout, derivs[k][count]);
		      }

		      count++;
		    }
	    }

	  } /* end of loop over symmetry unique quartets */

	} /* end of unique shell loops */

  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  free(PrintFourInd);
  for(i=0; i < Molecule.num_atoms*3; i++) fclose(out[i]);
  free(out);
  free_block(derivs);
  free_libderiv(&Libderiv);
#ifndef USE_TAYLOR_FM
  free_fjt_table(&fjt_table);
#endif

  return;
}
  }
}
