/*! \file
    \ingroup CINTS
    \brief Two-electron integrals.
*/
#include<cmath>
#include<cstring>
#include<cstdio>
#include<memory.h>
#include<cstdlib>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>

#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"schwartz.h"
#include"quartet_data.h"
#include"iwl_tebuf.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif

namespace psi { namespace cints {
  //! Two-electron integral computation.
void te_ints()
{
  const double toler = UserOptions.cutoff;

  /*--- ASCII file to print integrals ---*/
  FILE *eriout ;

  /*--- Various data structures ---*/
  struct iwlbuf ERIOUT;               /* IWL buffer for target integrals */
  struct tebuf *tot_data;              /* accum. for contracted integrals */
  struct shell_pair *sp_ij, *sp_kl;
  struct unique_shell_pair *usp_ij,*usp_kl;

  Libint_t Libint;
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;
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
  register int i, j, k, l, m, ii, jj, kk, ll;
  register int si, sj, sk, sl ;
  register int sii, sjj, skk, sll , slll;
  register int pi, pj, pk, pl;
  int max_pj, max_pl;
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

  int size, class_size;
  int max_class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;

  double so_int;
  double AB2, CD2;
  double *raw_data;             /* pointer to the unnormalized taregt quartet of integrals */
  double *data;                 /* pointer to the transformed normalized target quartet of integrals */
#ifdef NONDOUBLE_INTS
  REALTYPE *target_ints;        /* Pointer to the location of the target quartet on the stack of
                                   integrals quartets if libint.a is using other than regular doubles */
#endif
  double *puream_data;
  double **plist_data;
  double pkblock_end_value = 0.0;
  double temp;

  /*---------------
    Initialization
   ---------------*/
#if PRINT
  eriout = fopen("eriout.dat","w");
#endif
  iwl_buf_init(&ERIOUT, IOUnits.itap33, toler, 0, 0);
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4);
  init_fjt_table(&fjt_table);
#endif
  init_libint_base();

  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_cart_class_size = (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am]);
  max_num_unique_quartets = Symmetry.max_stab_index*
                            Symmetry.max_stab_index*
                            Symmetry.max_stab_index;
  tot_data = (struct tebuf*) malloc(max_num_unique_quartets*max_cart_class_size*sizeof(struct tebuf));
  memset(tot_data, 0, (max_num_unique_quartets*max_cart_class_size)*sizeof(struct tebuf));
#ifdef NONDOUBLE_INTS
  raw_data = init_array(max_cart_class_size);
#endif
  if (BasisSet.puream)
    puream_data = (double *) malloc(sizeof(double)*
                                    (BasisSet.max_am*2-1)*
                                    ioff[BasisSet.max_am]*
                                    ioff[BasisSet.max_am]*
                                    ioff[BasisSet.max_am]);
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  if (Symmetry.nirreps > 1) {
    if (BasisSet.puream)
      max_class_size = (2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1);
    else
      max_class_size = max_cart_class_size;
    plist_data = block_matrix(max_num_unique_quartets,max_class_size);
    sj_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sk_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sl_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  }

  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  init_libint(&Libint, BasisSet.max_am-1, max_num_prim_comb);

/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  if(UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"  -Electron repulsion integrals:\n\n");
  }
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

            /* place in "ascending" angular mom-
               my simple way of optimizing PHG recursion (VRR) */
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
                data = norm_quartet(raw_data, puream_data, orig_am, BasisSet.puream);
                /*--- copy data to plist_data to be used in the symmetrization step ---*/
                if (Symmetry.nirreps > 1)
                  memcpy(plist_data[plquartet],data,sizeof(double)*class_size);
              }
              else {
                temp = 0.0;
                for(p=0;p<num_prim_comb;p++)
                  temp += (double) Libint.PrimQuartet[p].F[0];
                if (Symmetry.nirreps > 1)
                  plist_data[plquartet][0] = temp;
                else {
#ifdef NONDOUBLE_INTS
                  raw_data[0] = temp;
                  data = raw_data;
#else
                  Libint.int_stack[0] = temp;
                  data = Libint.int_stack;
#endif
                }
              }

            } /* end of computing "petit" list */

            num = 0;
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

              NOTE: only the last case is commented.
             ------------------------------------------------------------------------*/
            bf_i = BasisSet.shells[si].fbf-1;
            if (usi == usj && usi == usk && usi == usl)
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
#if SCF_ONLY
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      index = ind_offset + k*nl + l;
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){
                        ioffset = bf_i + i;
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler)
                        if (so_i >= so_k) {
                          tot_data[num].i = (short int) so_i;
                          tot_data[num].j = (short int) so_j;
                          tot_data[num].k = (short int) so_k;
                          tot_data[num].l = (short int) so_l;
                          tot_data[num++].val = so_int;
                        }
                        else {
                          tot_data[num].i = (short int) so_k;
                          tot_data[num].j = (short int) so_l;
                          tot_data[num].k = (short int) so_i;
                          tot_data[num].l = (short int) so_j;
                          tot_data[num++].val = so_int;
                        }
                    }
                  }
              }
            else if (usi == usj && usk == usl)
              for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
                if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
                  for(ij=0;ij<npi_ij;ij++) {
                    so_i = usp_ij->SOpair_so_i[irrep][ij];
                    so_j = usp_ij->SOpair_so_j[irrep][ij];
                    i = usp_ij->SOpair_bf_i[irrep][ij];
                    j = usp_ij->SOpair_bf_j[irrep][ij];
                    ind_offset = (i*nj + j)*nk*nl;
                    for(kl=0;kl<npi_kl;kl++) {
                      so_k = usp_kl->SOpair_so_i[irrep][kl];
                      so_l = usp_kl->SOpair_so_j[irrep][kl];
#if SCF_ONLY
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      index = ind_offset + k*nl + l;
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){
                        ioffset = bf_i + i;
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler)
                        if (so_i >= so_k) {
                          tot_data[num].i = (short int) so_i;
                          tot_data[num].j = (short int) so_j;
                          tot_data[num].k = (short int) so_k;
                          tot_data[num].l = (short int) so_l;
                          tot_data[num++].val = so_int;
                        }
                        else {
                          tot_data[num].i = (short int) so_k;
                          tot_data[num].j = (short int) so_l;
                          tot_data[num].k = (short int) so_i;
                          tot_data[num].l = (short int) so_j;
                          tot_data[num++].val = so_int;
                        }
                    }
                  }
              }
            else if (usi == usk && usj == usl)
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
#if SCF_ONLY
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      index = ind_offset + k*nl + l;
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){
                        ioffset = bf_i + i;
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler)
                        if (so_i >= so_j)
                          if (so_k >= so_l)
                            if ((so_i > so_k) || (so_i == so_k && so_j >= so_l)) {
                              tot_data[num].i = (short int) so_i;
                              tot_data[num].j = (short int) so_j;
                              tot_data[num].k = (short int) so_k;
                              tot_data[num].l = (short int) so_l;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_k;
                              tot_data[num].j = (short int) so_l;
                              tot_data[num].k = (short int) so_i;
                              tot_data[num].l = (short int) so_j;
                              tot_data[num++].val = so_int;
                            }
                          else
                            if ((so_i > so_l) || (so_i == so_l && so_j >= so_k)) {
                              tot_data[num].i = (short int) so_i;
                              tot_data[num].j = (short int) so_j;
                              tot_data[num].k = (short int) so_l;
                              tot_data[num].l = (short int) so_k;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_l;
                              tot_data[num].j = (short int) so_k;
                              tot_data[num].k = (short int) so_i;
                              tot_data[num].l = (short int) so_j;
                              tot_data[num++].val = so_int;
                            }
                        else
                          if (so_k >= so_l)
                            if ((so_j > so_k) || (so_j == so_k && so_i >= so_l)) {
                              tot_data[num].i = (short int) so_j;
                              tot_data[num].j = (short int) so_i;
                              tot_data[num].k = (short int) so_k;
                              tot_data[num].l = (short int) so_l;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_k;
                              tot_data[num].j = (short int) so_l;
                              tot_data[num].k = (short int) so_j;
                              tot_data[num].l = (short int) so_i;
                              tot_data[num++].val = so_int;
                            }
                          else
                            if ((so_j > so_l) || (so_j == so_l && so_i >= so_k)) {
                              tot_data[num].i = (short int) so_j;
                              tot_data[num].j = (short int) so_i;
                              tot_data[num].k = (short int) so_l;
                              tot_data[num].l = (short int) so_k;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_l;
                              tot_data[num].j = (short int) so_k;
                              tot_data[num].k = (short int) so_j;
                              tot_data[num].l = (short int) so_i;
                              tot_data[num++].val = so_int;
                            }
                    }
                  }
              }
            else if (usi == usj)
              for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
                if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
                  for(ij=0;ij<npi_ij;ij++) {
                    so_i = usp_ij->SOpair_so_i[irrep][ij];
                    so_j = usp_ij->SOpair_so_j[irrep][ij];
                    i = usp_ij->SOpair_bf_i[irrep][ij];
                    j = usp_ij->SOpair_bf_j[irrep][ij];
                    ind_offset = (i*nj + j)*nk*nl;
                    for(kl=0;kl<npi_kl;kl++) {
                      so_k = usp_kl->SOpair_so_i[irrep][kl];
                      so_l = usp_kl->SOpair_so_j[irrep][kl];
#if SCF_ONLY
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      index = ind_offset + k*nl + l;
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){
                        ioffset = bf_i + i;
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler)
                        if (so_k >= so_l)
                          if ((so_i > so_k) || (so_i == so_k && so_j >= so_l)) {
                            tot_data[num].i = (short int) so_i;
                            tot_data[num].j = (short int) so_j;
                            tot_data[num].k = (short int) so_k;
                            tot_data[num].l = (short int) so_l;
                            tot_data[num++].val = so_int;
                          }
                          else {
                            tot_data[num].i = (short int) so_k;
                            tot_data[num].j = (short int) so_l;
                            tot_data[num].k = (short int) so_i;
                            tot_data[num].l = (short int) so_j;
                            tot_data[num++].val = so_int;
                          }
                        else
                          if ((so_i > so_l) || (so_i == so_l && so_j >= so_k)) {
                            tot_data[num].i = (short int) so_i;
                            tot_data[num].j = (short int) so_j;
                            tot_data[num].k = (short int) so_l;
                            tot_data[num].l = (short int) so_k;
                            tot_data[num++].val = so_int;
                          }
                          else {
                            tot_data[num].i = (short int) so_l;
                            tot_data[num].j = (short int) so_k;
                            tot_data[num].k = (short int) so_i;
                            tot_data[num].l = (short int) so_j;
                            tot_data[num++].val = so_int;
                          }
                    }
                  }
              }
            else if (usk == usl)
              for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
                if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
                  for(ij=0;ij<npi_ij;ij++) {
                    so_i = usp_ij->SOpair_so_i[irrep][ij];
                    so_j = usp_ij->SOpair_so_j[irrep][ij];
                    i = usp_ij->SOpair_bf_i[irrep][ij];
                    j = usp_ij->SOpair_bf_j[irrep][ij];
                    ind_offset = (i*nj + j)*nk*nl;
                    for(kl=0;kl<npi_kl;kl++) {
                      so_k = usp_kl->SOpair_so_i[irrep][kl];
                      so_l = usp_kl->SOpair_so_j[irrep][kl];
#if SCF_ONLY
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      index = ind_offset + k*nl + l;
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){
                        ioffset = bf_i + i;
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler)
                        if (so_i >= so_j)
                          if ((so_i > so_k) || (so_i == so_k && so_j >= so_l)){
                            tot_data[num].i = (short int) so_i;
                            tot_data[num].j = (short int) so_j;
                            tot_data[num].k = (short int) so_k;
                            tot_data[num].l = (short int) so_l;
                            tot_data[num++].val = so_int;
                          }
                          else {
                            tot_data[num].i = (short int) so_k;
                            tot_data[num].j = (short int) so_l;
                            tot_data[num].k = (short int) so_i;
                            tot_data[num].l = (short int) so_j;
                            tot_data[num++].val = so_int;
                          }
                        else
                          if ((so_j > so_k) || (so_j == so_k && so_i >= so_l)){
                            tot_data[num].i = (short int) so_j;
                            tot_data[num].j = (short int) so_i;
                            tot_data[num].k = (short int) so_k;
                            tot_data[num].l = (short int) so_l;
                            tot_data[num++].val = so_int;
                          }
                          else {
                            tot_data[num].i = (short int) so_k;
                            tot_data[num].j = (short int) so_l;
                            tot_data[num].k = (short int) so_j;
                            tot_data[num].l = (short int) so_i;
                            tot_data[num++].val = so_int;
                          }
                    }
                  }
              }
            else
              for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
                  /*---
                    npi_ij - number of pairs of SOs arising from the ij pair of unique shells
                             whose direct product transforms as irrep
                    check to see if for given irrep npi_ij and npi_kl are non-zero, i.e.
                    if there's a combination of SO's that thansforms as totally symm. irrep
                    ---*/
                if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
                  for(ij=0;ij<npi_ij;ij++) { /*--- Loop over SO pairs from usij ---*/
                    i = usp_ij->SOpair_bf_i[irrep][ij];       /*--- Basis function type to which this SO from usi corresponds,
                                                                    e.g. px, dxz, etc. ---*/
                    j = usp_ij->SOpair_bf_j[irrep][ij];
                    so_i = usp_ij->SOpair_so_i[irrep][ij];    /*--- Absolute index of this SO from usi ---*/
                    so_j = usp_ij->SOpair_so_j[irrep][ij];
                    ind_offset = (i*nj + j)*nk*nl;
                    for(kl=0;kl<npi_kl;kl++) { /*--- Loop over SO pairs from uskl ---*/
                      k = usp_kl->SOpair_bf_i[irrep][kl];
                      l = usp_kl->SOpair_bf_j[irrep][kl];
                      so_k = usp_kl->SOpair_so_i[irrep][kl];
                      so_l = usp_kl->SOpair_so_j[irrep][kl];
#if SCF_ONLY
                      /*--- Throw out the combinations which have SOs that belong to 4 different
                            irreps. Such integrals don't contribute to the Fock matrices ---*/
                      if (UserOptions.scf_only)
                        if (Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_j] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_k] &&
                            Symmetry.so2symblk[so_i] != Symmetry.so2symblk[so_l])
                          continue;
#endif
                      index = ind_offset + k*nl + l;          /*--- Based on basis function types figure
                                                                    out which particular integral from these
                                                                    quartets contributes ---*/
                      so_int = 0.0;
                      for(s=0;s<num_unique_quartets;s++){  /*--- Sum over petite list quartets ---*/
                        ioffset = bf_i + i;      /*--- Absolute basis function index ---*/
                        bf_j = sj_fbf_arr[s]+j;
                        bf_k = sk_fbf_arr[s]+k;
                        bf_l = sl_fbf_arr[s]+l;
                        so_int += Symmetry.usotao[so_i][ioffset]*
                                  Symmetry.usotao[so_j][bf_j]*
                                  Symmetry.usotao[so_k][bf_k]*
                                  Symmetry.usotao[so_l][bf_l]*
                                  plist_data[s][index];
                      }
                      if (fabs(so_int)>toler) /*--- For non-zero integrals pring
                                                    their indices into canonical order ---*/
                        if (so_i >= so_j)
                          if (so_k >= so_l)
                            if ((so_i > so_k) || (so_i == so_k && so_j >= so_l)) {
                              tot_data[num].i = (short int) so_i;
                              tot_data[num].j = (short int) so_j;
                              tot_data[num].k = (short int) so_k;
                              tot_data[num].l = (short int) so_l;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_k;
                              tot_data[num].j = (short int) so_l;
                              tot_data[num].k = (short int) so_i;
                              tot_data[num].l = (short int) so_j;
                              tot_data[num++].val = so_int;
                            }
                          else
                            if ((so_i > so_l) || (so_i == so_l && so_j >= so_k)) {
                              tot_data[num].i = (short int) so_i;
                              tot_data[num].j = (short int) so_j;
                              tot_data[num].k = (short int) so_l;
                              tot_data[num].l = (short int) so_k;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_l;
                              tot_data[num].j = (short int) so_k;
                              tot_data[num].k = (short int) so_i;
                              tot_data[num].l = (short int) so_j;
                              tot_data[num++].val = so_int;
                            }
                        else
                          if (so_k >= so_l)
                            if ((so_j > so_k) || (so_j == so_k && so_i >= so_l)) {
                              tot_data[num].i = (short int) so_j;
                              tot_data[num].j = (short int) so_i;
                              tot_data[num].k = (short int) so_k;
                              tot_data[num].l = (short int) so_l;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_k;
                              tot_data[num].j = (short int) so_l;
                              tot_data[num].k = (short int) so_j;
                              tot_data[num].l = (short int) so_i;
                              tot_data[num++].val = so_int;
                            }
                          else
                            if ((so_j > so_l) || (so_j == so_l && so_i >= so_k)) {
                              tot_data[num].i = (short int) so_j;
                              tot_data[num].j = (short int) so_i;
                              tot_data[num].k = (short int) so_l;
                              tot_data[num].l = (short int) so_k;
                              tot_data[num++].val = so_int;
                            }
                            else {
                              tot_data[num].i = (short int) so_l;
                              tot_data[num].j = (short int) so_k;
                              tot_data[num].k = (short int) so_j;
                              tot_data[num].l = (short int) so_i;
                              tot_data[num++].val = so_int;
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
                        if(fabs(data[index])>toler){
                            tot_data[num].i = (short int) (ii+ioffset);
                            tot_data[num].j = (short int) (jj+joffset);
                            tot_data[num].k = (short int) (kk+koffset);
                            tot_data[num].l = (short int) (ll+loffset);
                            tot_data[num].val = data[index];
                            num++;
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
                        if(fabs(data[index])>toler){
                          i = ii + ioffset;
                          j = jj + joffset;
                          k = kk + koffset;
                          l = ll + loffset;
                          if (i < j) {
                            SWAP(i,j);
                            SWAP(k,l);
                          }
                          if (i < k) {
                            SWAP(i,k);
                            SWAP(j,l);
                          }
                          tot_data[num].i = (short int) i;
                          tot_data[num].j = (short int) j;
                          tot_data[num].k = (short int) k;
                          tot_data[num].l = (short int) l;
                          tot_data[num].val = data[index];
                          num++;
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
                        if(fabs(data[index])>toler){
                          i = ii + ioffset;
                          j = jj + joffset;
                          k = kk + koffset;
                          l = ll + loffset;
                          if (i < j)
                            SWAP(i,j);
                          if (k < l)
                            SWAP(k,l);
                          if ((i < k) || (i == k && j < l)) {
                            SWAP(i,k);
                            SWAP(j,l);
                          }
                          tot_data[num].i = (short int) i;
                          tot_data[num].j = (short int) j;
                          tot_data[num].k = (short int) k;
                          tot_data[num].l = (short int) l;
                          tot_data[num].val = data[index];
                          num++;
                        }
                      }
                    }
                  }
                }
              }
            }


            if (num) { /* Let's see if we need to write out something */
              total_te_count += num;
              if (upk == num_unique_pk - 1) /* if this is the last quartet needed for a pk-block - let CSCF know
                                               by setting index i of the last integral to negative of itself.
                                               The only guy where this trick won't work will be (00|00).
                                               But normally (00|00) is of (ss|ss) type and is enough for computing one
                                               pk-matrix element. PK-buffer won't get overfull because of this one guy */
                tot_data[num-1].i = -tot_data[num-1].i;
              iwl_buf_wrt_struct_nocut(&ERIOUT, tot_data, num);
            }

            if(UserOptions.print_lvl >= PRINT_DEBUG) {
              /* fprintf(outfile,"  -Electron repulsion integrals:\n\n"); */
              for(n=0; n<num; n++){
                fprintf(outfile, "%5d%5d%5d%5d%20.10lf\n",
                        abs(tot_data[n].i),
                        tot_data[n].j,
                        tot_data[n].k,
                        tot_data[n].l,
                        tot_data[n].val);
              }
            }

#if PRINT
            for(n=0; n<num; n++){
              fprintf(eriout, "%5d%5d%5d%5d%20.10lf\n",
                      abs(tot_data[n].i)+1,
                      tot_data[n].j+1,
                      tot_data[n].k+1,
                      tot_data[n].l+1,
                      tot_data[n].val);
            }
            fflush(eriout);
#endif
          } /* end of computing PK-quartets. */

        } /* end getting unique shell combination */

  iwl_buf_flush(&ERIOUT, 1);
  iwl_buf_close(&ERIOUT, 1);

#if PRINT
  fclose(eriout);
#endif


  fprintf(outfile,"\n    Wrote %d two-electron integrals to IWL file %2d\n\n",total_te_count,IOUnits.itap33);

  /*---------
    Clean-up
   ---------*/
  free_libint(&Libint);
  free(tot_data);
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  if (Symmetry.nirreps > 1) {
    free_block(plist_data);
    free(sj_fbf_arr);
    free(sk_fbf_arr);
    free(sl_fbf_arr);
  }
  if (BasisSet.puream)
    free(puream_data);
#ifdef NONDOUBLE_INTS
  free(raw_data);
#endif
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt_table(&fjt_table);
  free_fjt();
#endif

  return;
}
}}
