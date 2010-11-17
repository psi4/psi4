/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory.h>
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
#include <stdexcept>
#include <Tools/prints.h>


namespace psi {
  namespace cints {

void te_deriv1_corr()
{
  //FILE *debugfile; // jtf
  /*--- Various data structures ---*/
  struct iwlbuf TPDM;                 /* IWL buffer for two-pdm matrix elements */
  struct shell_pair *sp_ij, *sp_kl;
  Libderiv_t Libderiv;                    /* Integrals library object */
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;               /* table of auxiliary function F_m(u) for each primitive combination */
#endif

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
  //debugfile = fopen("te_deriv1_debug", "w"); //jtf
  iwl_buf_init(&TPDM, IOUnits.itapG, 0.0, 1, 1);
  buf_offset = 0;
  buf_4offset = 0;
  buf_size = TPDM.inbuf;
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+DERIV_LVL);
  init_fjt_table(&fjt_table);
#endif
  init_libderiv_base();

  max_cart_class_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  max_class_size = max_cart_class_size;
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
		      (BasisSet.max_num_prims*BasisSet.max_num_prims);
  init_libderiv1(&Libderiv,BasisSet.max_am-1,max_num_prim_comb,max_cart_class_size);
  FourInd = init_array(max_cart_class_size);

  grad_te_local = block_matrix(Molecule.num_atoms,3);

  //fprintf(debugfile, "Symmetry Info:\n");
  //fprintf(debugfile, "nirrep = %d\n", Symmetry.nirreps);
  //fprintf(debugfile, "num_so = %d\n", Symmetry.num_so);
  //fprintf(debugfile, "num_unique_shells = %d\n", Symmetry.num_unique_shells);
  //fprintf(debugfile, "Symmetry.usotao matrix: \n");
  //for(i=0; i<Symmetry.num_so; i++){
    //for(j=0; j<Symmetry.num_so; j++){
      //fprintf(debugfile, "%lf  ", Symmetry.usotao[i][j]);
    //}
    //fprintf(debugfile, "\n");
  //}

/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  for (sii = 0; sii < BasisSet.num_shells; sii++)
    for (sjj = 0; sjj <= sii; sjj++)
      for (skk = 0; skk <= sii; skk++)
        for (sll = 0; sll <= ((sii == skk) ? sjj : skk); sll++) {

          si = sii;
          sj = sjj;
          sk = skk;
          sl = sll;

          /*--- Skip this quartet if all four centers are the same ---*/
          if (BasisSet.shells[si].center == BasisSet.shells[sj].center
              && BasisSet.shells[si].center == BasisSet.shells[sk].center
              && BasisSet.shells[si].center == BasisSet.shells[sl].center) {
            /*--- If reading in density - need to skip the appropriate block of that too ---*/
            /*	      if (read_dens) {
             last_buf = TPDM.lastbuf;
             quartet_done = 0;
             do {
             if (buf_offset < buf_size) {
             i = TPDM.labels[buf_4offset]   - si_fao;
             if (i < 0)
             quartet_done = 1;
             buf_offset++;
             buf_4offset += 4;
             }
             else if (!last_buf) {
             iwl_buf_fetch(&TPDM);
             buf_offset = 0;
             buf_4offset = 0;
             last_buf = TPDM.lastbuf;
             }
             else {
             punt(fpo,"  The last TPDM quartet not marked\n");
             }
             } while (!quartet_done);
             }*/

            continue;
          }

          switch_ij = 0;
          switch_kl = 0;
          switch_ijkl = 0;
          /* place in "ascending" angular mom-
           my simple way of optimizing PHG recursion (VRR) */
          /* these first two are good for the HRR */
          if (BasisSet.shells[si].am < BasisSet.shells[sj].am) {
            dum = si;
            si = sj;
            sj = dum;
            switch_ij = 1;
          }
          if (BasisSet.shells[sk].am < BasisSet.shells[sl].am) {
            dum = sk;
            sk = sl;
            sl = dum;
            switch_kl = 1;
          }
          /* this should be /good/ for the VRR */
          if (BasisSet.shells[si].am + BasisSet.shells[sj].am
              > BasisSet.shells[sk].am + BasisSet.shells[sl].am) {
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
          quartet_size = ni * nj * nk * nl;

          np_i = BasisSet.shells[si].n_prims;
          np_j = BasisSet.shells[sj].n_prims;
          np_k = BasisSet.shells[sk].n_prims;
          np_l = BasisSet.shells[sl].n_prims;

          orig_am[0] = BasisSet.shells[si].am - 1;
          orig_am[1] = BasisSet.shells[sj].am - 1;
          orig_am[2] = BasisSet.shells[sk].am - 1;
          orig_am[3] = BasisSet.shells[sl].am - 1;
          am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];

          sp_ij = &(BasisSet.shell_pairs[si][sj]);
          sp_kl = &(BasisSet.shell_pairs[sk][sl]);

          Libderiv.AB[0] = sp_ij->AB[0];
          Libderiv.AB[1] = sp_ij->AB[1];
          Libderiv.AB[2] = sp_ij->AB[2];
          Libderiv.CD[0] = sp_kl->AB[0];
          Libderiv.CD[1] = sp_kl->AB[1];
          Libderiv.CD[2] = sp_kl->AB[2];

          AB2 = Libderiv.AB[0] * Libderiv.AB[0] + Libderiv.AB[1]
              * Libderiv.AB[1] + Libderiv.AB[2] * Libderiv.AB[2];
          CD2 = Libderiv.CD[0] * Libderiv.CD[0] + Libderiv.CD[1]
              * Libderiv.CD[1] + Libderiv.CD[2] * Libderiv.CD[2];

          /*-------------------------
           Figure out the prefactor
           -------------------------*/
          pfac = 1.0;
          if (si == sj) pfac *= 0.5;
          if (sk == sl) pfac *= 0.5;
          if (si == sk && sj == sl || si == sl && sj == sk) pfac *= 0.5;
          pfac *= 8.0; /*--- The factor of 8 needed for correlated densities ---*/

          /*--- Compute data for primitive quartets here ---*/
          num_prim_comb = 0;
          for (pi = 0; pi < np_i; pi++) {
            max_pj = (si == sj) ? pi + 1 : np_j;
            for (pj = 0; pj < max_pj; pj++) {
              m = (1 + (si == sj && pi != pj));
              for (pk = 0; pk < np_k; pk++) {
                max_pl = (sk == sl) ? pk + 1 : np_l;
                for (pl = 0; pl < max_pl; pl++) {
                  n = m * (1 + (sk == sl && pk != pl));
#ifdef USE_TAYLOR_FM
                  deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
                      NULL, AB2, CD2,
                      sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
#else
                  deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
                      &fjt_table, AB2, CD2, sp_ij, sp_kl, am, pi, pj, pk, pl, n
                          * pfac);
#endif    

                }
              }
            }
          }

          /*--- Read in a shell quartet from disk ---*/
          memset(FourInd, 0, sizeof(double) * quartet_size);
          last_buf = TPDM.lastbuf;
          quartet_done = 0;
          sii_fao = BasisSet.shells[sii].fao - 1;
          sjj_fao = BasisSet.shells[sjj].fao - 1;
          skk_fao = BasisSet.shells[skk].fao - 1;
          sll_fao = BasisSet.shells[sll].fao - 1;
//          fprintf(debugfile, "\t Start quartet: %d %d %d %d\n", sii, sjj, skk, sll); // jtf
//          fprintf(debugfile, "Symm:\t %d, %d, %d, %d\n ", Symmetry. );
          do {
            if (buf_offset < buf_size) {
              i = TPDM.labels[buf_4offset] - sii_fao;
              if (i >= 0) {
                j = TPDM.labels[buf_4offset + 1] - sjj_fao;
                k = TPDM.labels[buf_4offset + 2] - skk_fao;
                l = TPDM.labels[buf_4offset + 3] - sll_fao;
                if (switch_ij) {
                  dum = i;
                  i = j;
                  j = dum;
                }
                if (switch_kl) {
                  dum = k;
                  k = l;
                  l = dum;
                }
                if (switch_ijkl) {
                  dum = i;
                  i = k;
                  k = dum;
                  dum = j;
                  j = l;
                  l = dum;
                }
                offset = ((i * nj + j) * nk + k) * nl + l;

                FourInd[offset] += TPDM.values[buf_offset]
                    * GTOs.bf_norm[orig_am[0]][i] * GTOs.bf_norm[orig_am[1]][j]
                    * GTOs.bf_norm[orig_am[2]][k] * GTOs.bf_norm[orig_am[3]][l];
                //fprintf(debugfile, "TPDM[%d]:\t %d, %d, %d, %d\t%lf\t%d %lf\n",
                    //buf_offset, i, j, k, l, TPDM.values[buf_offset], offset,
                    //FourInd[offset]); // jtf


              } else quartet_done = 1;
              buf_offset++;
              buf_4offset += 4;
            } else if (!last_buf) {
              iwl_buf_fetch(&TPDM);
              buf_offset = 0;
              buf_4offset = 0;
              last_buf = TPDM.lastbuf;
            } else {
              throw std::domain_error("The last TPDM quartet not marked");
            }
          } while (!quartet_done);

          build_deriv1_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](
              &Libderiv, num_prim_comb);

          center_i = BasisSet.shells[si].center - 1;
          center_j = BasisSet.shells[sj].center - 1;
          center_k = BasisSet.shells[sk].center - 1;
          center_l = BasisSet.shells[sl].center - 1;
          // jtf added debugging print statements
          //fprintf(debugfile, "Centers: %d, %d, %d, %d\n", center_i, center_j, center_k, center_l);
          //fprintf(debugfile, "Shells:  %d, %d, %d, %d\n", si, sj, sk, sl);
          ddax = 0.0;
          for (k = 0; k < quartet_size; k++){
            //fprintf(debugfile, "\t %lf", Libderiv.ABCD[0][k]);
            ddax += Libderiv.ABCD[0][k] * FourInd[k];
          }
          //fprintf(debugfile, "\n");
          grad_te_local[center_i][0] += ddax;

          dday = 0.0;
          for (k = 0; k < quartet_size; k++){
            //fprintf(debugfile, "\t %lf", Libderiv.ABCD[1][k]);
            dday += Libderiv.ABCD[1][k] * FourInd[k];
          }
          //fprintf(debugfile, "\n");
          grad_te_local[center_i][1] += dday;

          ddaz = 0.0;
          for (k = 0; k < quartet_size; k++){
            //fprintf(debugfile, "\t %lf", Libderiv.ABCD[2][k]);
            ddaz += Libderiv.ABCD[2][k] * FourInd[k];
          }
          //fprintf(debugfile, "\n");
          grad_te_local[center_i][2] += ddaz;

          /*ddbx = 0.0;
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
          for (k = 0; k < quartet_size; k++)
            ddcx += Libderiv.ABCD[6][k] * FourInd[k];
          grad_te_local[center_k][0] += ddcx;

          ddcy = 0.0;
          for (k = 0; k < quartet_size; k++)
            ddcy += Libderiv.ABCD[7][k] * FourInd[k];
          grad_te_local[center_k][1] += ddcy;

          ddcz = 0.0;
          for (k = 0; k < quartet_size; k++)
            ddcz += Libderiv.ABCD[8][k] * FourInd[k];
          grad_te_local[center_k][2] += ddcz;

          dddx = 0.0;
          for (k = 0; k < quartet_size; k++)
            dddx += Libderiv.ABCD[9][k] * FourInd[k];
          grad_te_local[center_l][0] += dddx;

          dddy = 0.0;
          for (k = 0; k < quartet_size; k++)
            dddy += Libderiv.ABCD[10][k] * FourInd[k];
          grad_te_local[center_l][1] += dddy;

          dddz = 0.0;
          for (k = 0; k < quartet_size; k++)
            dddz += Libderiv.ABCD[11][k] * FourInd[k];
          grad_te_local[center_l][2] += dddz;

          grad_te_local[center_j][0] -= ddax + ddcx + dddx;
          grad_te_local[center_j][1] -= dday + ddcy + dddy;
          grad_te_local[center_j][2] -= ddaz + ddcz + dddz;
        }

  if (UserOptions.print_lvl >= PRINT_TEDERIV)
    print_atomvec((char *)"Two-electron contribution to the forces (a.u.)", grad_te_local);



  for(i=0;i<Molecule.num_atoms;i++) {
    Grad[i][0] += grad_te_local[i][0];
    Grad[i][1] += grad_te_local[i][1];
    Grad[i][2] += grad_te_local[i][2];
  }
  
  /*---------
    Clean-up
   ---------*/
  free(FourInd);
  free_block(grad_te_local);
  free_libderiv(&Libderiv);
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt_table(&fjt_table);
  free_fjt();
#endif
  iwl_buf_close(&TPDM,1);
  //fclose(debugfile); // jtf
  return;
}
}}
