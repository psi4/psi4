/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory.h>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <libqt/qt.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#include "moinfo.h"
#include "compute_scf_opdm.h"
#include "norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif
#include "iwl_tebuf.h"
#include "symmetrize.h"
#include "small_fns.h"
#include "compute_eri.h"

namespace psi { namespace CINTS {
static inline int hash(int L, int M, int N);

void giao_te_deriv(void)
{
  /*--- Various data structures ---*/
  Libint_t Libint;                    /* Integrals library object */
  struct iwlbuf G, DGDBX, DGDBY, DGDBZ;  // IWL buffers

  /*---------------
    Initialization
    ---------------*/
  // NOTE: we are computing GIAO derivative integrals via ERIs in which one of functions
  // has angular momentum incremented by 1
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4+1,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+1);
//  init_fjt_table(&fjt_table);
#endif
  init_libint_base();

  int max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
    (BasisSet.max_num_prims*BasisSet.max_num_prims);
  init_libint(&Libint,BasisSet.max_am,max_num_prim_comb);

  iwl_buf_init(&G, 47, 0.0, 0, 0);
  iwl_buf_init(&DGDBX, IOUnits.itapdgdB[0], 0.0, 0, 0);
  iwl_buf_init(&DGDBY, IOUnits.itapdgdB[1], 0.0, 0, 0);
  iwl_buf_init(&DGDBZ, IOUnits.itapdgdB[2], 0.0, 0, 0);
  
  int max_cart_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  double* abcd_buf = new double[max_cart_size];
  double* dgdBx_buf = new double[max_cart_size];
  double* dgdBy_buf = new double[max_cart_size];
  double* dgdBz_buf = new double[max_cart_size];
  struct tebuf* g_tebuf = new struct tebuf[max_cart_size];
  struct tebuf* dgdBx_tebuf = new struct tebuf[max_cart_size];
  struct tebuf* dgdBy_tebuf = new struct tebuf[max_cart_size];
  struct tebuf* dgdBz_tebuf = new struct tebuf[max_cart_size];
  int max_cart_size_1000 = ioff[BasisSet.max_am+1]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  double* ap1bcd_buf = new double[max_cart_size_1000];
  double* abcp1d_buf = new double[max_cart_size_1000];

  const double dgdB_pfac = 0.5;
  
  for (int sii=0; sii<BasisSet.num_shells; sii++) {
    for (int sjj=0; sjj<BasisSet.num_shells; sjj++) {
      for (int skk=0; skk<BasisSet.num_shells; skk++) {
	for (int sll=0; sll<BasisSet.num_shells; sll++) {
	  int si = sii;
          int sj = sjj;
          int sk = skk;
          int sl = sll;

          if (si==0 && sj==4 && sk==8 && sl==9) {
            int stop = 0;
          }
          
          // These are angular momenta of GIAO Gaussians
          int am[4];
          am[0] = BasisSet.shells[si].am-1;
          am[1] = BasisSet.shells[sj].am-1;
          am[2] = BasisSet.shells[sk].am-1;
          am[3] = BasisSet.shells[sl].am-1;
          int total_am = am[0]+am[1]+am[2]+am[3];

          // If these shells are on the same origin -- GIAO derivatives will be 0 -- skip
          int center[4];
          center[0] = BasisSet.shells[si].center-1;
          center[1] = BasisSet.shells[sj].center-1;
          center[2] = BasisSet.shells[sk].center-1;
          center[3] = BasisSet.shells[sl].center-1;
//          if (center[0] == center[1] && center[0] == center[2] && center[0] == center[3])
//            continue;
            
          struct shell_pair* sp_ij = &(BasisSet.shell_pairs[si][sj]);
          struct shell_pair* sp_kl = &(BasisSet.shell_pairs[sk][sl]);
          struct coordinates A, C, AB, CD, ABxA_plus_CDxC;
          AB.x = sp_ij->AB[0];
          AB.y = sp_ij->AB[1];
          AB.z = sp_ij->AB[2];
          CD.x = sp_kl->AB[0];
          CD.y = sp_kl->AB[1];
          CD.z = sp_kl->AB[2];
          A = Molecule.centers[center[0]];
          C = Molecule.centers[center[2]];
          ABxA_plus_CDxC.x = (AB.y*A.z - AB.z*A.y) + (CD.y*C.z - CD.z*C.y);
          ABxA_plus_CDxC.y = (AB.z*A.x - AB.x*A.z) + (CD.z*C.x - CD.x*C.z);
          ABxA_plus_CDxC.z = (AB.x*A.y - AB.y*A.x) + (CD.x*C.y - CD.y*C.x);

          //
          // first compute (ab|cd) integrals
          //
          int inc[4] = {0,0,0,0};
          int num_ints = compute_eri(abcd_buf, &Libint, si, sj, sk, sl, inc[0], inc[1], inc[2], inc[3], false);

          // These are angular momenta of Gaussians over which we compute first quartet of ERIs
          int gam[4];
          for(int i=0; i<4; i++)
            gam[i] = am[i] + inc[i];

          int nao[4];
          nao[0] = ioff[gam[0]+1];
          nao[1] = ioff[gam[1]+1];
          nao[2] = ioff[gam[2]+1];
          nao[3] = ioff[gam[3]+1];
          int quartet_size = nao[0]*nao[1]*nao[2]*nao[3];
          if (num_ints)
          for(int i=0; i<quartet_size; i++) {
            double value = abcd_buf[i];
            dgdBx_buf[i] = dgdB_pfac * ABxA_plus_CDxC.x * value;
            dgdBy_buf[i] = dgdB_pfac * ABxA_plus_CDxC.y * value;
            dgdBz_buf[i] = dgdB_pfac * ABxA_plus_CDxC.z * value;
          }
          else
          for(int i=0; i<quartet_size; i++) {
            abcd_buf[i] = 0.0;
            dgdBx_buf[i] = 0.0;
            dgdBy_buf[i] = 0.0;
            dgdBz_buf[i] = 0.0;
          }

          //
          // then compute (a+1 b|cd) integrals
          //
          
          // if AB=0 then this contribution is 0
          if (center[0] != center[1]) {
            inc[0] = 1;
            num_ints = compute_eri(ap1bcd_buf, &Libint, si, sj, sk, sl, inc[0], inc[1], inc[2], inc[3], false);
            int gam_ap1[4];
            for(int i=0; i<4; i++)
              gam_ap1[i] = am[i] + inc[i];
            int nao_ap1[4];
            nao_ap1[0] = ioff[gam_ap1[0]+1];
            nao_ap1[1] = ioff[gam_ap1[1]+1];
            nao_ap1[2] = ioff[gam_ap1[2]+1];
            nao_ap1[3] = ioff[gam_ap1[3]+1];
            int quartet_size_ap1bcd = nao_ap1[0]*nao_ap1[1]*nao_ap1[2]*nao_ap1[3];
            inc[0] = 0;

            // (a+1x b|cd) contributes to dg/dBy and dg/dBz with these coefficients
            double ap1xpfac_y =  dgdB_pfac * AB.z;
            double ap1xpfac_z = -dgdB_pfac * AB.y;
            // (a+1y b|cd) contributes to dg/dBx and dg/dBz with these coefficients
            double ap1ypfac_x = -dgdB_pfac * AB.z;
            double ap1ypfac_z =  dgdB_pfac * AB.x;
            // (a+1z b|cd) contributes to dg/dBx and dg/dBy with these coefficients
            double ap1zpfac_x =  dgdB_pfac * AB.y;
            double ap1zpfac_y = -dgdB_pfac * AB.x;

            // Contract ERIs into the target GIAO derivative integral
            /*--- create all am components of si ---*/
            int a = 0;
            double* dgdBx_ptr = dgdBx_buf;
            double* dgdBy_ptr = dgdBy_buf;
            double* dgdBz_ptr = dgdBz_buf;
	    for(int ii = 0; ii <= am[0]; ii++){
	      int l1 = am[0] - ii;
	      for(int jj = 0; jj <= ii; jj++, a++){
	        int m1 = ii - jj;
	        int n1 = jj;
                
                int ap1x = hash(l1+1,m1,n1);
                int ap1y = hash(l1,m1+1,n1);
                int ap1z = hash(l1,m1,n1+1);

                int njkl = nao[1]*nao[2]*nao[3];
                double* ap1xbcd_ptr = ap1bcd_buf + njkl*ap1x;
                double* ap1ybcd_ptr = ap1bcd_buf + njkl*ap1y;
                double* ap1zbcd_ptr = ap1bcd_buf + njkl*ap1z;
                for(int jkl=0; jkl<njkl; jkl++, ap1xbcd_ptr++, ap1ybcd_ptr++, ap1zbcd_ptr++,
                                         dgdBx_ptr++, dgdBy_ptr++, dgdBz_ptr++) {
                  double ap1xbcd = *ap1xbcd_ptr;
                  double ap1ybcd = *ap1ybcd_ptr;
                  double ap1zbcd = *ap1zbcd_ptr;
                  *dgdBx_ptr += ap1ybcd * ap1ypfac_x + ap1zbcd * ap1zpfac_x;
                  *dgdBy_ptr += ap1zbcd * ap1zpfac_y + ap1xbcd * ap1xpfac_y;
                  *dgdBz_ptr += ap1xbcd * ap1xpfac_z + ap1ybcd * ap1ypfac_z;
                }
              }
            }
            
          }

          //
          // then compute (ab|c+1 d) integrals
          //

          // if CD=0 then this contribution is 0
          if (center[2] != center[3]) {
            inc[2] = 1;
            num_ints = compute_eri(abcp1d_buf, &Libint, si, sj, sk, sl, inc[0], inc[1], inc[2], inc[3], false);
            int gam_cp1[4];
            for(int i=0; i<4; i++)
              gam_cp1[i] = am[i] + inc[i];
            int nao_cp1[4];
            nao_cp1[0] = ioff[gam_cp1[0]+1];
            nao_cp1[1] = ioff[gam_cp1[1]+1];
            nao_cp1[2] = ioff[gam_cp1[2]+1];
            nao_cp1[3] = ioff[gam_cp1[3]+1];
            int quartet_size_abcp1d = nao_cp1[0]*nao_cp1[1]*nao_cp1[2]*nao_cp1[3];
            inc[2] = 0;
            
            // (ab|c+1x d) contributes to dg/dBy and dg/dBz with these coefficients
            double cp1xpfac_y =  dgdB_pfac * CD.z;
            double cp1xpfac_z = -dgdB_pfac * CD.y;
            // (ab|c+1y d) contributes to dg/dBx and dg/dBz with these coefficients
            double cp1ypfac_x = -dgdB_pfac * CD.z;
            double cp1ypfac_z =  dgdB_pfac * CD.x;
            // (ab|c+1z d) contributes to dg/dBx and dg/dBy with these coefficients
            double cp1zpfac_x =  dgdB_pfac * CD.y;
            double cp1zpfac_y = -dgdB_pfac * CD.x;

            // Contract ERIs into the target GIAO derivative integral
            int nij = nao[0]*nao[1];
            int nk1l = nao_cp1[2] * nao[3];
            int nl = nao[3];
            double* dgdBx_ptr = dgdBx_buf;
            double* dgdBy_ptr = dgdBy_buf;
            double* dgdBz_ptr = dgdBz_buf;
            double* ab_ptr = abcp1d_buf;
            for(int ij=0; ij<nij; ij++, ab_ptr+=nk1l) {
              /*--- create all am components of sk ---*/
              int c = 0;
              for(int ii = 0; ii <= am[2]; ii++){
                int l1 = am[2] - ii;
	        for(int jj = 0; jj <= ii; jj++, c++){
	          int m1 = ii - jj;
	          int n1 = jj;
                
                  int cp1x = hash(l1+1,m1,n1);
                  int cp1y = hash(l1,m1+1,n1);
                  int cp1z = hash(l1,m1,n1+1);

                  double* abcp1xd_ptr = ab_ptr + cp1x*nl;
                  double* abcp1yd_ptr = ab_ptr + cp1y*nl;
                  double* abcp1zd_ptr = ab_ptr + cp1z*nl;
                  for(int l=0; l<nl; l++, abcp1xd_ptr++, abcp1yd_ptr++, abcp1zd_ptr++,
                                          dgdBx_ptr++, dgdBy_ptr++, dgdBz_ptr++) {
                    double abcp1xd = *abcp1xd_ptr;
                    double abcp1yd = *abcp1yd_ptr;
                    double abcp1zd = *abcp1zd_ptr;
                    *dgdBx_ptr += abcp1yd * cp1ypfac_x + abcp1zd * cp1zpfac_x;
                    *dgdBy_ptr += abcp1zd * cp1zpfac_y + abcp1xd * cp1xpfac_y;
                    *dgdBz_ptr += abcp1xd * cp1xpfac_z + abcp1yd * cp1ypfac_z;
                  }
                }
              }
            }

          }

          //
          // Compute dg/dB over GIAO Gaussians
          //
          

          // Normalize integrals
          double* g_target_buf;
          double* dgdBx_target_buf;
          double* dgdBy_target_buf;
          double* dgdBz_target_buf;
          if (total_am) {
            // No transformation to spherical harmonics basis, only normalization
            g_target_buf = norm_quartet(abcd_buf, NULL, am, 0);
            dgdBx_target_buf = norm_quartet(dgdBx_buf, NULL, am, 0);
            dgdBy_target_buf = norm_quartet(dgdBy_buf, NULL, am, 0);
            dgdBz_target_buf = norm_quartet(dgdBz_buf, NULL, am, 0);
          }
          else {
            g_target_buf = abcd_buf;
            dgdBx_target_buf = dgdBx_buf;
            dgdBy_target_buf = dgdBy_buf;
            dgdBz_target_buf = dgdBz_buf;
          }

          // Print integrals out
          if(UserOptions.print_lvl >= PRINT_DEBUG) {
            int ijkl = 0;
            for(int i=0; i<nao[0]; i++) {
              int ii = i+BasisSet.shells[si].fao-1;
              for(int j=0; j<nao[1]; j++) {
                int jj = j+BasisSet.shells[sj].fao-1;
                for(int k=0; k<nao[2]; k++) {
                  int kk = k+BasisSet.shells[sk].fao-1;
                  for(int l=0; l<nao[3]; l++, ijkl++) {
                    int ll = l+BasisSet.shells[sl].fao-1;

                    if (fabs(abcd_buf[ijkl]) > UserOptions.cutoff) {
                      fprintf(stdout, "%5d%5d%5d%5d%20.10lf\n",
                              ii, jj, kk, ll, abcd_buf[ijkl]);
                    }
/*                    if (fabs(dgdBx_buf[ijkl]) > UserOptions.cutoff) {
                      fprintf(stdout, "%5d%5d%5d%5d%20.10lf\n",
                              ii, jj, kk, ll, dgdBx_buf[ijkl]);
                    }
                    if (fabs(dgdBy_buf[ijkl]) > UserOptions.cutoff) {
                      fprintf(stdout, "%5d%5d%5d%5d%20.10lf\n",
                              ii, jj, kk, ll, dgdBy_buf[ijkl]);
                    }
                    if (fabs(dgdBz_buf[ijkl]) > UserOptions.cutoff) {
                      fprintf(stdout, "%5d%5d%5d%5d%20.10lf\n",
                              ii, jj, kk, ll, dgdBz_buf[ijkl]);
                    }*/
                  }
                }
	      } 
	    }
          }

          // Write integrals out
          int fao[4];
          fao[0] = BasisSet.shells[si].fao - 1;
          fao[1] = BasisSet.shells[sj].fao - 1;
          fao[2] = BasisSet.shells[sk].fao - 1;
          fao[3] = BasisSet.shells[sl].fao - 1;
          int ijkl = 0;
          for(int i=0; i<nao[0]; i++) {
            int ii = i + fao[0];
            for(int j=0; j<nao[1]; j++) {
              int jj = j + fao[1];
              for(int k=0; k<nao[2]; k++) {
                int kk = k + fao[2];
                for(int l=0; l<nao[3]; l++, ijkl++) {
                  int ll = l + fao[3];
                    
                  g_tebuf[ijkl].i = ii;
                  dgdBx_tebuf[ijkl].i = ii;
                  dgdBy_tebuf[ijkl].i = ii;
                  dgdBz_tebuf[ijkl].i = ii;
                  g_tebuf[ijkl].j = jj;
                  dgdBx_tebuf[ijkl].j = jj;
                  dgdBy_tebuf[ijkl].j = jj;
                  dgdBz_tebuf[ijkl].j = jj;
                  g_tebuf[ijkl].k = kk;
                  dgdBx_tebuf[ijkl].k = kk;
                  dgdBy_tebuf[ijkl].k = kk;
                  dgdBz_tebuf[ijkl].k = kk;
                  g_tebuf[ijkl].l = ll;
                  dgdBx_tebuf[ijkl].l = ll;
                  dgdBy_tebuf[ijkl].l = ll;
                  dgdBz_tebuf[ijkl].l = ll;
                  g_tebuf[ijkl].val = g_target_buf[ijkl];
                  dgdBx_tebuf[ijkl].val = dgdBx_target_buf[ijkl];
                  dgdBy_tebuf[ijkl].val = dgdBy_target_buf[ijkl];
                  dgdBz_tebuf[ijkl].val = dgdBz_target_buf[ijkl];
                }
              }
            }
          }
          iwl_buf_wrt_struct(&G, g_tebuf, quartet_size, UserOptions.cutoff);
          iwl_buf_wrt_struct(&DGDBX, dgdBx_tebuf, quartet_size, UserOptions.cutoff);
          iwl_buf_wrt_struct(&DGDBY, dgdBy_tebuf, quartet_size, UserOptions.cutoff);
          iwl_buf_wrt_struct(&DGDBZ, dgdBz_tebuf, quartet_size, UserOptions.cutoff);

	} /* sll */
      } /* skk */
    } /* sjj */
  } /* sii */

  iwl_buf_flush(&G, 1);  iwl_buf_close(&G, 1);
  iwl_buf_flush(&DGDBX, 1);  iwl_buf_close(&DGDBX, 1);
  iwl_buf_flush(&DGDBY, 1);  iwl_buf_close(&DGDBY, 1);
  iwl_buf_flush(&DGDBZ, 1);  iwl_buf_close(&DGDBZ, 1);

  /*---------
    Clean-up
    ---------*/
  free_libint(&Libint);
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt();
#endif

  return;
}

static inline int hash(int L, int M, int N)
{
  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};

  int MpN = M+N;
  int result = N + io[MpN];
  return result;
}
};};
