/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<cstdlib>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<libr12/libr12.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"r12_quartet_data.h"
#include"iwl_tebuf.h"
#include"norm_quartet.h"
#include"int_fjt.h"


namespace psi { namespace CINTS {
void r12_te_ints()
{
  const double toler = UserOptions.cutoff;
  
  /*--- ASCII file to print integrals ---*/
  FILE *teout[NUM_TE_TYPES];
  char teout_filename[] = "teout0.dat";

  /*--- Various data structures ---*/
  struct iwlbuf TEOUT[NUM_TE_TYPES-1];             /* IWL buffer for target integrals */
  int itapTE[NUM_TE_TYPES-1];
  struct coordinates ericent[4];       /* coordinates of centers A, B, C, and D */
  struct tebuf *tot_data[NUM_TE_TYPES];              /* accum. for contracted integrals */
  struct shell_pair *sp_ij, *sp_kl;
  struct unique_shell_pair *usp_ij,*usp_kl;
  Libr12_t Libr12;
  double_array_t fjt_table;

  static char *te_operator[] = { "1/r12", "r12", "[r12,T1]" };
  int total_te_count[NUM_TE_TYPES] = {0, 0, 0, 0};
  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int te_type;
  int n, num[NUM_TE_TYPES];
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  register int i, j, k, l, ii, jj, kk, ll;
  register int si, sj, sk, sl ;
  register int sii, sjj, skk, sll , slll;
  register int pi, pj, pk, pl ;
  register int pii, pjj, pkk, pll ;
  int upk, num_unique_pk;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int *sj_arr, *sk_arr, *sl_arr;
  int *sj_fbf_arr, *sk_fbf_arr, *sl_fbf_arr;
  int usii,usjj,uskk,usll,usi,usj,usk,usl;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int class_size;
  int max_class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int so_ii, so_jj, so_kk, so_ll;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;

  double so_int[NUM_TE_TYPES];
  double AB2, CD2;
  double *data[NUM_TE_TYPES];
  double *puream_data[NUM_TE_TYPES];
  double **plist_data[NUM_TE_TYPES];
  double pkblock_end_value = 0.0;
  double value;
  double temp;
  double ssss, ss_r12_ss;

  /*---------------
    Initialization
   ---------------*/
#if PRINT
  for(te_type=0;te_type<NUM_TE_TYPES;te_type++) {
    teout_filename[5] = te_type + '0';
    teout[te_type] = fopen(teout_filename,"w");
  }
#endif
  itapTE[0] = IOUnits.itap33;
  itapTE[1] = IOUnits.itapR12;
  itapTE[2] = IOUnits.itapT1;
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    iwl_buf_init(&TEOUT[te_type], itapTE[te_type], toler, 0, 0);
  }
  init_fjt(BasisSet.max_am*4);
  init_libr12_base();

  
  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_cart_class_size = ioff[BasisSet.max_am]*
			ioff[BasisSet.max_am]*
			ioff[BasisSet.max_am]*
			ioff[BasisSet.max_am];
  max_num_unique_quartets = Symmetry.max_stab_index*Symmetry.max_stab_index*Symmetry.max_stab_index;
  for(te_type=0;te_type<NUM_TE_TYPES;te_type++) {
    tot_data[te_type] = (struct tebuf*) malloc(max_num_unique_quartets*max_cart_class_size*sizeof(struct tebuf));
    memset(tot_data[te_type], 0, (max_num_unique_quartets*max_cart_class_size)*sizeof(struct tebuf));
  }
  if (BasisSet.puream) {
    if (Symmetry.nirreps == 1) { /*--- allocate NUM_TE_TYPES scratch arrays for cart->puream transformation ---*/
      for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
	puream_data[te_type] = (double *) malloc(sizeof(double)*
						 (BasisSet.max_am*2-1)*
						 ioff[BasisSet.max_am]*
						 ioff[BasisSet.max_am]*
						 ioff[BasisSet.max_am]);
    }
    else {
      puream_data[0] = (double *) malloc(sizeof(double)*
					 (BasisSet.max_am*2-1)*
					 ioff[BasisSet.max_am]*
					 ioff[BasisSet.max_am]*
					 ioff[BasisSet.max_am]);
      for(te_type=1;te_type<NUM_TE_TYPES;te_type++)
	puream_data[te_type] = puream_data[0];
    }
  }
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  if (Symmetry.nirreps > 1) {
    if (BasisSet.puream)
      max_class_size = (2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1);
    else
      max_class_size = max_cart_class_size;
    for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
      plist_data[te_type] = block_matrix(max_num_unique_quartets,max_class_size);
    sj_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sk_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sl_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  }

  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  UserOptions.memory -= init_libr12(&Libr12,BasisSet.max_am-1,max_num_prim_comb);
  init_fjt_table(&fjt_table);

/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  for (usii=0; usii<Symmetry.num_unique_shells; usii++)
    for (usjj=0; usjj<=usii; usjj++)
      for (uskk=0; uskk<=usjj; uskk++)
	for (usll=0; usll<=uskk; usll++){

          /*--- Decide what shell quartets out of (ij|kl), (ik|jl), and (il|jk) are unique ---*/
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

	  for(upk=0;upk<num_unique_pk;upk++) {
	    /*--- For each combination of unique shells generate "petit list" of shells ---*/
	    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

	    /* place in "ascending" angular momentum order-
	       remember that [r12,T1] do not like it, I will have to worry about that later */
	    /* these first two are good for the HRR */
	    if(BasisSet.shells[Symmetry.us2s[usi]].am < BasisSet.shells[Symmetry.us2s[usj]].am){
	      dum = usi;
	      usi = usj;
	      usj = dum;
	    }
	    if(BasisSet.shells[Symmetry.us2s[usk]].am < BasisSet.shells[Symmetry.us2s[usl]].am){
	      dum = usk;
	      usk = usl;
	      usl = dum;
	    }
	    /* this should be /good/ for the VRR */
	    if(BasisSet.shells[Symmetry.us2s[usi]].am + BasisSet.shells[Symmetry.us2s[usj]].am >
	       BasisSet.shells[Symmetry.us2s[usk]].am + BasisSet.shells[Symmetry.us2s[usl]].am){
	      dum = usi;
	      usi = usk;
	      usk = dum;
	      dum = usj;
	      usj = usl;
	      usl = dum;
	    }

	    si = Symmetry.us2s[usi];
	    sjj = Symmetry.us2s[usj];
	    skk = Symmetry.us2s[usk];
	    sll = Symmetry.us2s[usl];
	    if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
	      /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
	      usp_ij = &(Symmetry.us_pairs[usi][usj]);
	      usp_kl = &(Symmetry.us_pairs[usk][usl]);
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
	      memset(sj_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sk_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sl_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
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

		    total_am = BasisSet.shells[si].am+BasisSet.shells[sj].am+BasisSet.shells[sk].am+BasisSet.shells[sl].am;
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
		      sj_fbf_arr[count] = BasisSet.shells[sj].fbf-1;
		      sk_fbf_arr[count] = BasisSet.shells[sk].fbf-1;
		      sl_fbf_arr[count] = BasisSet.shells[sl].fbf-1;
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

	    np_i = BasisSet.shells[si].n_prims;
	    np_j = BasisSet.shells[sjj].n_prims;
	    np_k = BasisSet.shells[skk].n_prims;
	    np_l = BasisSet.shells[sll].n_prims;
	    
	    orig_am[0] = BasisSet.shells[si].am-1;
	    orig_am[1] = BasisSet.shells[sjj].am-1;
	    orig_am[2] = BasisSet.shells[skk].am-1;
	    orig_am[3] = BasisSet.shells[sll].am-1;
	    am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];


	    /*----------------------------------
	      Compute the nonredundant quartets
	     ----------------------------------*/
	    for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	      sj = sj_arr[plquartet];
	      sk = sk_arr[plquartet];
	      sl = sl_arr[plquartet];

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
		build_r12_grt[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libr12, num_prim_comb);
		/*--- copy transformed data to plist_data to be used in the symmetrization step ---*/
		if (Symmetry.nirreps > 1)
		  for(te_type=0;te_type<NUM_TE_TYPES;te_type++) {
		    data[te_type] = norm_quartet(Libr12.te_ptr[te_type], puream_data[te_type], orig_am, BasisSet.puream);
		    memcpy(plist_data[te_type][plquartet],data[te_type],sizeof(double)*class_size);
		  }
		/*--- or just transform data ---*/
		else
		  for(te_type=0;te_type<NUM_TE_TYPES;te_type++) {
		    data[te_type] = norm_quartet(Libr12.te_ptr[te_type], puream_data[te_type], orig_am, BasisSet.puream);
		  }

	      }
	      else {
		ssss = 0.0;
		ss_r12_ss = 0.0;
		for(p=0;p<num_prim_comb;p++) {
		  ssss += Libr12.PrimQuartet[p].F[0];
		  ss_r12_ss += Libr12.PrimQuartet[p].ss_r12_ss;
		}
		build_r12_grt[0][0][0][0](&Libr12,num_prim_comb);
		if (Symmetry.nirreps > 1) {
		  plist_data[0][plquartet][0] = ssss;
		  plist_data[1][plquartet][0] = ss_r12_ss;
		  plist_data[2][plquartet][0] = Libr12.te_ptr[2][0];
		  plist_data[3][plquartet][0] = Libr12.te_ptr[3][0];
		}
		else {
		    /*--------------------------------------------------------
		      This was so freaking careless on my part. I forgot that
		      in hrr_grt_order_0000 te_ptr[2] and te_ptr[3] point to
		      int_stack[1] and int_stack[0] respectively. Therefore I
		      have to copy these into int_stack[2] and int_stack[3] first
		      before copying ssss and ss_r12_ss into those locations.
		     --------------------------------------------------------*/
		  Libr12.int_stack[2] = Libr12.te_ptr[2][0];
		  Libr12.int_stack[3] = Libr12.te_ptr[3][0];
		  Libr12.int_stack[0] = ssss;
		  Libr12.int_stack[1] = ss_r12_ss;
		  data[0] = Libr12.int_stack;
		  data[1] = Libr12.int_stack+1;
		  data[2] = Libr12.int_stack+2;
		  data[3] = Libr12.int_stack+3;
		}
	      }

	    } /* end of computing "petit" list */

	    memset(num,0,NUM_TE_TYPES*sizeof(int));
	    if (Symmetry.nirreps > 1) { /*--- Non-C1 case ---*/
	    /*------------------------------------------------------------------------
	      Now we have everything to build SO's. Need to distinguish several cases
	      that are slightly different from each other. To avoid extra if's inside
	      the loops I separated them:
	      1) usi == usj == usk == usl
	      2) usi == usj != usk == usl
	      3) usi == usk != usj == usl
	      4) usi == usj
	      5) usk == usl
	      6) general case
	      "Symmetrization" is based on Pitzer's equal contribution theorem.

	      The general procedure is as follows:
	      1) loop over irreps, if there are any SO products that arise from
	         this usp_ij pair and belong to this irrep, loop over them (ij)
	      2) from ij figure out product of which SOs this is, and what AOs
	         (or should I say, basis functions) contribute to these SO.
		 so_i and so_j are absolute indices, but i and j are relative
		 indices (i and j are indices of cart/puream components in usi
		 and usj)
	      3) if there are any products of SOs from usp_kl that belong to irrep,
	         loop over those (kl). Obviously, the whole point of introducing
		 restrictions on different cases was to improve efficiency, but now
		 it seems like a stupid idea. I should have left just one loop case.

		 For more comments see Default_Ints/te_ints.c
	     ------------------------------------------------------------------------*/
	    bf_i = BasisSet.shells[si].fbf-1;
	    if (usi == usj && usi == usk && usi == usl || usi == usk && usj == usl)
	      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
		if (npi_ij = usp_ij->SOpair_npi[irrep])
		  for(ij=0;ij<npi_ij;ij++) {
		    so_i = usp_ij->SOpair_so_i[irrep][ij];
		    so_j = usp_ij->SOpair_so_j[irrep][ij];
		    i = usp_ij->SOpair_bf_i[irrep][ij];
		    j = usp_ij->SOpair_bf_j[irrep][ij];
		    ind_offset = (i*nj + j)*nk*nl;
		    for(kl=0;kl<=ij;kl++) {
		      so_k = usp_kl->SOpair_so_i[irrep][kl];
		      so_l = usp_kl->SOpair_so_j[irrep][kl];
		      k = usp_kl->SOpair_bf_i[irrep][kl];
		      l = usp_kl->SOpair_bf_j[irrep][kl];
		      index = ind_offset + k*nl + l;
		      memset(so_int,0,NUM_TE_TYPES*sizeof(double));
		      for(s=0;s<num_unique_quartets;s++){
			ioffset = bf_i + i;
			bf_j = sj_fbf_arr[s]+j;
			bf_k = sk_fbf_arr[s]+k;
			bf_l = sl_fbf_arr[s]+l;
			temp = Symmetry.usotao[so_i][ioffset]*
			       Symmetry.usotao[so_j][bf_j]*
			       Symmetry.usotao[so_k][bf_k]*
			       Symmetry.usotao[so_l][bf_l];
			for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			  so_int[te_type] += temp*plist_data[te_type][s][index];
		      }
		      for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			if (fabs(so_int[te_type])>toler) {
			  so_ii = so_i;
			  so_jj = so_j;
			  so_kk = so_k;
			  so_ll = so_l;
			  if (so_ii < so_jj) {
			    SWAP(so_ii,so_jj);
			    if (te_type == 2) /* [r12,T1] is antisymmetric WRT i and j */
			      so_int[te_type] *= -1.0;
			  }
			  if (so_kk < so_ll) {
			    SWAP(so_kk,so_ll);
			    if (te_type == 3) /* [r12,T2] is antisymmetric WRT to k and l */
			      so_int[te_type] *= -1.0;
			  }
			  if ((so_ii < so_kk) || (so_ii == so_kk && so_jj < so_ll)) {
			    if (te_type < 2) { /* only eri and r12 ints have bra-ket symmetry */
			      SWAP(so_ii,so_kk);
			      SWAP(so_jj,so_ll);
			    }
			  }
			  tot_data[te_type][num[te_type]].i = (short int) so_ii;
			  tot_data[te_type][num[te_type]].j = (short int) so_jj;
			  tot_data[te_type][num[te_type]].k = (short int) so_kk;
			  tot_data[te_type][num[te_type]].l = (short int) so_ll;
			  tot_data[te_type][num[te_type]].val = so_int[te_type];
			  num[te_type]++;
			}
		    }
		  }
	      }
	    else
	      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
		if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
		  for(ij=0;ij<npi_ij;ij++) {
		    i = usp_ij->SOpair_bf_i[irrep][ij];
		    j = usp_ij->SOpair_bf_j[irrep][ij];
		    so_i = usp_ij->SOpair_so_i[irrep][ij];
		    so_j = usp_ij->SOpair_so_j[irrep][ij];
		    ind_offset = (i*nj + j)*nk*nl;
		    for(kl=0;kl<npi_kl;kl++) {
		      k = usp_kl->SOpair_bf_i[irrep][kl];
		      l = usp_kl->SOpair_bf_j[irrep][kl];
		      so_k = usp_kl->SOpair_so_i[irrep][kl];
		      so_l = usp_kl->SOpair_so_j[irrep][kl];
		      index = ind_offset + k*nl + l;
		      memset(so_int,0,NUM_TE_TYPES*sizeof(double));
		      for(s=0;s<num_unique_quartets;s++){
			ioffset = bf_i + i;
			bf_j = sj_fbf_arr[s]+j;
			bf_k = sk_fbf_arr[s]+k;
			bf_l = sl_fbf_arr[s]+l;
			temp = Symmetry.usotao[so_i][ioffset]*
			       Symmetry.usotao[so_j][bf_j]*
			       Symmetry.usotao[so_k][bf_k]*
			       Symmetry.usotao[so_l][bf_l];
			for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			  so_int[te_type] += temp*plist_data[te_type][s][index];
		      }
		      for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			if (fabs(so_int[te_type])>toler) {
			  so_ii = so_i;
			  so_jj = so_j;
			  so_kk = so_k;
			  so_ll = so_l;
			  if (so_ii < so_jj) {
			    SWAP(so_ii,so_jj);
			    if (te_type == 2) /* [r12,T1] is antisymmetric WRT i and j */
			      so_int[te_type] *= -1.0;
			  }
			  if (so_kk < so_ll) {
			    SWAP(so_kk,so_ll);
			    if (te_type == 3) /* [r12,T2] is antisymmetric WRT to k and l */
			      so_int[te_type] *= -1.0;
			  }
			  if ((so_ii < so_kk) || (so_ii == so_kk && so_jj < so_ll)) {
			    if (te_type < 2) { /* only eri and r12 ints have bra-ket symmetry */
			      SWAP(so_ii,so_kk);
			      SWAP(so_jj,so_ll);
			    }
			  }
			  tot_data[te_type][num[te_type]].i = (short int) so_ii;
			  tot_data[te_type][num[te_type]].j = (short int) so_jj;
			  tot_data[te_type][num[te_type]].k = (short int) so_kk;
			  tot_data[te_type][num[te_type]].l = (short int) so_ll;
			  tot_data[te_type][num[te_type]].val = so_int[te_type];
			  num[te_type]++;
			}
		    }
		  }
	      }
	    }
	    else { /*--- C1 symmetry ---*/
	      /*--- Here just put non-redundant integrals to tot_data ---*/
	      if(usi==usj&&usk==usl&&usi==usk) { /*--- All shells are the same - the (aa|aa) case ---*/
		iimax = ni - 1;
		for(ii=0; ii <= iimax; ii++){
		  jjmax = ii;
		  for(jj=0; jj <= jjmax; jj++){
		    kkmax = ii;
		    for(kk=0; kk <= kkmax; kk++){
		      llmax = (kk==ii)? jj : kk ;
		      for(ll=0; ll <= llmax; ll++){
			index = ll+nl*(kk+nk*(jj+nj*ii));
			for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			  if(fabs(data[te_type][index])>toler){
			    tot_data[te_type][num[te_type]].i = (short int) (ii+ioffset);
			    tot_data[te_type][num[te_type]].j = (short int) (jj+joffset);
			    tot_data[te_type][num[te_type]].k = (short int) (kk+koffset);
			    tot_data[te_type][num[te_type]].l = (short int) (ll+loffset);
			    tot_data[te_type][num[te_type]].val = data[te_type][index];
			    num[te_type]++;
			  }
		      }
		    }
		  }
		}
	      }
	      else if(usi==usk && usj==usl){    /*--- The (ab|ab) case ---*/
		iimax = ni - 1;
		for(ii=0; ii <= iimax; ii++){
		  jjmax = nj - 1;
		  for(jj=0; jj <= jjmax; jj++){
		    kkmax = ii;
		    for(kk=0; kk <= kkmax; kk++){
		      llmax = (kk==ii)? jj : nl - 1;
		      for(ll=0; ll <= llmax; ll++){
			index = ll+nl*(kk+nk*(jj+nj*ii));
			for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			  if(fabs(data[te_type][index])>toler){
			    i = ii + ioffset;
			    j = jj + joffset;
			    k = kk + koffset;
			    l = ll + loffset;
			    value = data[te_type][index];
			    if (i < j) {
			      SWAP(i,j);
			      SWAP(k,l);
			      if (te_type > 1) /* [r12,Ti] are antisymmetric WRT to i -> j, k -> l */
				value *= -1.0;
			    }
			    if (i < k) {
			      if (te_type < 2) { /* only eri and r12 ints have bra-ket symmetry */
				SWAP(i,k);
				SWAP(j,l);
			      }
			    }
			    tot_data[te_type][num[te_type]].i = (short int) i;
			    tot_data[te_type][num[te_type]].j = (short int) j;
			    tot_data[te_type][num[te_type]].k = (short int) k;
			    tot_data[te_type][num[te_type]].l = (short int) l;
			    tot_data[te_type][num[te_type]].val = value;
			    num[te_type]++;
			  }
		      }
		    }
		  }
		}
	      }
	      else {   /*--- The (ab|cd) case ---*/
		iimax = ni - 1;
		kkmax = nk - 1;
		for(ii=0; ii <= iimax; ii++){
		  jjmax = (usi==usj) ? ii : nj - 1;
		  for(jj=0; jj <= jjmax; jj++){
		    for(kk=0; kk <= kkmax; kk++){
		      llmax = (usk==usl) ? kk : nl - 1;
		      for(ll=0; ll <= llmax; ll++){
			index = ll+nl*(kk+nk*(jj+nj*ii));
			for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
			  if(fabs(data[te_type][index])>toler){
			    i = ii + ioffset;
			    j = jj + joffset;
			    k = kk + koffset;
			    l = ll + loffset;
			    value = data[te_type][index];
			    if (i < j) {
			      SWAP(i,j);
			      if (te_type == 2) /* [r12,T1] is antisymmetric WRT i and j */
				value *= -1.0;
			    }
			    if (k < l) {
			      SWAP(k,l);
			      if (te_type == 3) /* [r12,T2] is antisymmetric WRT k and l */
				value *= -1.0;
			    }
			    if ((i < k) || (i == k && j < l)) {
			      if (te_type < 2) { /* only eri and r12 ints have bra-ket symmetry */
				SWAP(i,k);
				SWAP(j,l);
			      }
			    }
			    tot_data[te_type][num[te_type]].i = (short int) i;
			    tot_data[te_type][num[te_type]].j = (short int) j;
			    tot_data[te_type][num[te_type]].k = (short int) k;
			    tot_data[te_type][num[te_type]].l = (short int) l;
			    tot_data[te_type][num[te_type]].val = value;
			    num[te_type]++;
			  }
		      }
		    }
		  }
		}
	      }
	    }


	    /*---------------
	      Write out ERIs
	     ---------------*/
	    if (num[0]) { /* Let's see if we need to write out something */
	      total_te_count[0] += num[0];
	      if (upk == num_unique_pk - 1) /* if this is the last quartet needed for a pk-block - let CSCF know
					       by setting index i of the last integral to negative of itself.
					       The only guy where this trick won't work will be (00|00).
					       But normally (00|00) is of (ss|ss) type and is enough for computing one
					       pk-matrix element. PK-buffer won't get overfull because of this one guy */
		tot_data[0][num[0]-1].i = -tot_data[0][num[0]-1].i;
	      iwl_buf_wrt_struct_nocut(&TEOUT[0], tot_data[0], num[0]);
	    }

	    /*----------------
	      Write out r12's
	     ----------------*/
	    if (num[1]) { /* Let's see if we need to write out something */
	      total_te_count[1] += num[1];
	      iwl_buf_wrt_struct_nocut(&TEOUT[1], tot_data[1], num[1]);
	    }

	    if (NUM_TE_TYPES > 2) {
	      /*---------------------
		Write out [r12,T1]'s
	       ---------------------*/
	      if (num[2]) { /* Let's see if we need to write out something */
		total_te_count[2] += num[2];
		iwl_buf_wrt_struct_nocut(&TEOUT[2], tot_data[2], num[2]);
	      }

	      /*-----------------------------------
		Convert [r12,T2]'s into [r12,T1]'s
	       -----------------------------------*/
	      for(i=0;i<num[3];i++)
		if (tot_data[3][i].i == tot_data[3][i].k && tot_data[3][i].j == tot_data[3][i].l) {
		  /* Do NOT write out if ij == kl because (ij|[r12,T1]|ij) = (ij|[r12,T2]|ij) */
		  tot_data[3][i].val = 0.0;
		  total_te_count[2]--;
		  continue;
		}
		else {
		  SWAP(tot_data[3][i].i,tot_data[3][i].k);
		  SWAP(tot_data[3][i].j,tot_data[3][i].l);
		}
	      if (num[3]) { /* Let's see if we need to write out something */
		total_te_count[2] += num[3];
		iwl_buf_wrt_struct(&TEOUT[2], tot_data[3], num[3], toler);
	      }
	    }
	    
#if PRINT
	    for(te_type=0;te_type<NUM_TE_TYPES;te_type++) {
	      if (te_type < 3) {
		for(n=0; n<num[te_type]; n++){
		  fprintf(teout[te_type], "%5d%5d%5d%5d%20.10lf\n",
			  abs(tot_data[te_type][n].i)+1, 
			  tot_data[te_type][n].j+1, 
			  tot_data[te_type][n].k+1, 
			  tot_data[te_type][n].l+1, 
			  tot_data[te_type][n].val);
		}
		fflush(teout[te_type]);
	      }
	      else { /* We've swapped ij and kl in [r12,T2] case, remember? */
		for(n=0; n<num[te_type]; n++){
		  fprintf(teout[te_type], "%5d%5d%5d%5d%20.10lf\n",
			  tot_data[te_type][n].k+1, 
			  tot_data[te_type][n].l+1, 
			  tot_data[te_type][n].i+1, 
			  tot_data[te_type][n].j+1, 
			  tot_data[te_type][n].val);
		}
		fflush(teout[te_type]);
	      }
	    }
#endif
	  } /* end of computing PK-quartets. */

	} /* end getting unique shell combination */

  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    iwl_buf_flush(&TEOUT[te_type], 1);
    iwl_buf_close(&TEOUT[te_type], 1);
  }

#if PRINT
  for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
    fclose(teout[te_type]);
#endif


  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++)
    fprintf(outfile,"    Wrote %d two-electron integrals of %s to IWL file %2d\n",
	    total_te_count[te_type],te_operator[te_type],itapTE[te_type]);
  fprintf(outfile,"\n");

  /*---------
    Clean-up
   ---------*/
  free_libr12(&Libr12);
  free_fjt(&fjt_table);
  for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
    free(tot_data[te_type]);
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  if (Symmetry.nirreps > 1) {
    for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
      free_block(plist_data[te_type]);
    free(sj_fbf_arr);
    free(sk_fbf_arr);
    free(sl_fbf_arr);
  }
  if (BasisSet.puream) {
    if (Symmetry.nirreps > 1)
      free(puream_data[0]);
    else
      for(te_type=0;te_type<NUM_TE_TYPES;te_type++)
	free(puream_data[te_type]);
  }

  return;
}

};};
