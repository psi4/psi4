/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void setup_LR(struct RHO_Params rho_params)
{
  dpdfile2 L1,R1,Z1,F;
  dpdbuf4 L2,R2,Z2;
  int i,j,L_irr, R_irr, G_irr, L_root, R_root;
  char L1A_lbl[32], L1B_lbl[32], L2AA_lbl[32], L2BB_lbl[32], L2AB_lbl[32], L2RHF_lbl[32];
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32], R2RHF_lbl[32];
  char lbl[32];
  double tval, one_energy, this_energy, total_two_energy, two_energy, R0;

  L_irr = rho_params.L_irr;
  R_irr = rho_params.R_irr;
  G_irr = rho_params.G_irr;
  L_root = rho_params.L_root;
  R_root = rho_params.R_root;
  R0 = rho_params.R0;

  /*
  fprintf(outfile,"\n");
  if(L_root+1 == 0)
    fprintf(outfile,"\tDensity for GS %3s\n",moinfo.labels[L_irr]);
  else 
    fprintf(outfile,"\tDensity for ES %d%3s\n",L_root+1,moinfo.labels[L_irr]);
  fflush(outfile);
  */
  
  /*
  fprintf(outfile,"\n\tSetting up L and R to compute density\n");
  fprintf(outfile,"\tLeft-hand eigenvector: symmetry %s and excited root %d\n",
	  moinfo.labels[L_irr], L_root+1); 
  fprintf(outfile,"\tRight-hand eigenvector: symmetry %s and excited root %d\n",
	  moinfo.labels[R_irr], R_root+1); 
  fprintf(outfile,"\tR0 = %15.10lf\n",params.R0); 
  */

  /* form labels for the L to be copied */
  sprintf(L1A_lbl,"LIA %d %d", L_irr, L_root);
  sprintf(L1B_lbl,"Lia %d %d", L_irr, L_root);
  sprintf(L2AA_lbl,"LIJAB %d %d", L_irr, L_root );
  sprintf(L2BB_lbl,"Lijab %d %d", L_irr, L_root );
  sprintf(L2AB_lbl,"LIjAb %d %d", L_irr, L_root );
  sprintf(L2RHF_lbl,"2LIjAb - LIjbA %d %d", L_irr, L_root );

  /*** GLG is used by the ground state density code and contains R0*L + Zeta */
  /*** the first term disappears if L_irr != G_irr (excitation is asymmetric) */

  if ( (L_irr == G_irr) && (!params.calc_xi) ) {
    if (params.ref == 0) {
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2RHF_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "2LIjAb - LIjbA");
      dpd_buf4_close(&L2);
    }
    if ( (params.ref == 0) || (params.ref == 1) ) {
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1A_lbl);
      dpd_file2_copy(&L1, CC_GLG, "LIA");
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1B_lbl);
      dpd_file2_copy(&L1, CC_GLG, "Lia");
      dpd_file2_close(&L1);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2AA_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "LIJAB");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2BB_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "Lijab");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2AB_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "LIjAb");
      dpd_buf4_close(&L2);
    }
    else if (params.ref == 2) {
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1A_lbl);
      dpd_file2_copy(&L1, CC_GLG, "LIA");
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 2, 3, L1B_lbl);
      dpd_file2_copy(&L1, CC_GLG, "Lia");
      dpd_file2_close(&L1);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2AA_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "LIJAB");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 12, 17, 12, 17, 0, L2BB_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "Lijab");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, L2AB_lbl);
      dpd_buf4_copy(&L2, CC_GLG, "LIjAb");
      dpd_buf4_close(&L2);
    }
  }
 
  /* put copy of L in CC_GL for excited state parts of density */
  if (!params.ground) {
    if (params.ref == 0) {
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2RHF_lbl);
      dpd_buf4_copy(&L2, CC_GL, "2LIjAb - LIjbA");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
      dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "2LIjAb - LIjbA (IA,jb)");
      dpd_buf4_close(&L2);
    }
    if ( (params.ref==0) || (params.ref==1) ) {
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1A_lbl);
      dpd_file2_copy(&L1, CC_GL, "LIA");
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1B_lbl);
      dpd_file2_copy(&L1, CC_GL, "Lia");
      dpd_file2_close(&L1);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2AA_lbl);
      dpd_buf4_copy(&L2, CC_GL, "LIJAB");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2BB_lbl);
      dpd_buf4_copy(&L2, CC_GL, "Lijab");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2AB_lbl);
      dpd_buf4_copy(&L2, CC_GL, "LIjAb");
      dpd_buf4_close(&L2);
    }
    else if (params.ref == 2) {
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 0, 1, L1A_lbl);
      dpd_file2_copy(&L1, CC_GL, "LIA");
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_LAMPS, L_irr, 2, 3, L1B_lbl);
      dpd_file2_copy(&L1, CC_GL, "Lia");
      dpd_file2_close(&L1);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, L2AA_lbl);
      dpd_buf4_copy(&L2, CC_GL, "LIJAB");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 12, 17, 12, 17, 0, L2BB_lbl);
      dpd_buf4_copy(&L2, CC_GL, "Lijab");
      dpd_buf4_close(&L2);
      dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, L2AB_lbl);
      dpd_buf4_copy(&L2, CC_GL, "LIjAb");
      dpd_buf4_close(&L2);
    }
  }

  /* for ground-state density contributions L <- R0 L + Zeta */
  /* symmetry of L must be same as density */

  if ( (!rho_params.L_ground) && (!params.calc_xi)) {
    if (params.connect_xi) rho_params.R0 = 0.0;
    if ( (params.ref==0) || (params.ref==1) ) {
      if (L_irr == G_irr) {
        dpd_file2_init(&L1, CC_GLG, L_irr, 0, 1, "LIA");
        dpd_file2_scm(&L1, rho_params.R0);
        dpd_file2_close(&L1);
        dpd_file2_init(&L1, CC_GLG, L_irr, 0, 1, "Lia");
        dpd_file2_scm(&L1, rho_params.R0);
        dpd_file2_close(&L1);
        dpd_buf4_init(&L2, CC_GLG, L_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, L_irr, 2, 7, 2, 7, 0, "Lijab");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, L_irr, 0, 5, 0, 5, 0, "LIjAb");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
      }
      else {
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
        dpd_file2_scm(&L1, 0.0);
        dpd_file2_close(&L1);
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "Lia");
        dpd_file2_scm(&L1, 0.0);
        dpd_file2_close(&L1);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
      }
      if (params.use_zeta) {
        dpd_file2_init(&Z1, CC_LAMPS, G_irr, 0, 1, "ZIA");
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
        dpd_file2_axpy(&Z1, &L1, 1.0, 0);
        dpd_file2_close(&L1);
        dpd_file2_close(&Z1);
        dpd_file2_init(&Z1, CC_LAMPS, G_irr, 0, 1, "Zia");
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "Lia");
        dpd_file2_axpy(&Z1, &L1, 1.0, 0);
        dpd_file2_close(&L1);
        dpd_file2_close(&Z1);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 2, 7, 2, 7, 0, "ZIJAB");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 2, 7, 2, 7, 0, "Zijab");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 0, 5, 0, 5, 0, "ZIjAb");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
      }
    }
    else if (params.ref==2) {
      if (L_irr == G_irr) {
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
        dpd_file2_scm(&L1, rho_params.R0);
        dpd_file2_close(&L1);
        dpd_file2_init(&L1, CC_GLG, G_irr, 2, 3, "Lia");
        dpd_file2_scm(&L1, rho_params.R0);
        dpd_file2_close(&L1);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
        dpd_buf4_scm(&L2, rho_params.R0);
        dpd_buf4_close(&L2);
      }
      else {
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
        dpd_file2_scm(&L1, 0.0);
        dpd_file2_close(&L1);
        dpd_file2_init(&L1, CC_GLG, G_irr, 2, 3, "Lia");
        dpd_file2_scm(&L1, 0.0);
        dpd_file2_close(&L1);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
        dpd_buf4_init(&L2, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
        dpd_buf4_scm(&L2, 0.0);
        dpd_buf4_close(&L2);
      }
      /* check magnitude */
      dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
      tval = dpd_file2_dot_self(&L1);
      dpd_file2_close(&L1);
      fprintf(outfile,"Ro*L+Zeta in CC_GLG, LIA before zeta: %15.10lf\n",tval);

      if (params.use_zeta) {
        dpd_file2_init(&Z1, CC_LAMPS, G_irr, 0, 1, "ZIA");
        dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
        dpd_file2_axpy(&Z1, &L1, 1.0, 0);
        dpd_file2_close(&L1);
        dpd_file2_close(&Z1);
        dpd_file2_init(&Z1, CC_LAMPS, G_irr, 2, 3, "Zia");
        dpd_file2_init(&L1, CC_GLG, G_irr, 2, 3, "Lia");
        dpd_file2_axpy(&Z1, &L1, 1.0, 0);
        dpd_file2_close(&L1);
        dpd_file2_close(&Z1);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 2, 7, 2, 7, 0, "ZIJAB");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 12, 17, 12, 17, 0, "Zijab");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
        dpd_buf4_init(&Z2, CC_LAMPS, G_irr, 22, 28, 22, 28, 0, "ZIjAb");
        dpd_buf4_init(&L2, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
        dpd_buf4_axpy(&Z2, &L2, 1.0);
        dpd_buf4_close(&L2);
        dpd_buf4_close(&Z2);
	/* check magnitude */
	dpd_file2_init(&L1, CC_GLG, G_irr, 0, 1, "LIA");
	tval = dpd_file2_dot_self(&L1);
	dpd_file2_close(&L1);
	fprintf(outfile,"Ro*L+Zeta in CC_GLG, LIA: %15.10lf\n",tval);
      }
    }
  }

  /* now sort entries in CC_GLG */
  if (!params.calc_xi) {
    if( (params.ref==0) || (params.ref==1) ) { /** RHF/ROHF **/
      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GLG, qpsr, 0, 5, "LiJaB");
      dpd_buf4_close(&L2);
    
      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 2, 7, 0, "LIJAB");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LIAJB");
      dpd_buf4_close(&L2);
      
      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 2, 7, 0, "Lijab");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "Liajb");
      dpd_buf4_close(&L2);
     
      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LIAjb");
      dpd_buf4_close(&L2);
        
      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LiJaB");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LiaJB");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAjb");
      dpd_buf4_sort(&L2, CC_GLG, psrq, 10, 10, "LIbjA");
      dpd_buf4_sort(&L2, CC_GLG, rqps, 10, 10, "LjAIb");
      dpd_buf4_close(&L2);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_buf4_init(&L2, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GLG, qpsr, 23, 29, "LiJaB");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 0, 5, 2, 7, 0, "LIJAB");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 20, 20, "LIAJB");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 10, 15, 12, 17, 0, "Lijab");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 30, 30, "Liajb");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 20, 30, "LIAjb");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 23, 29, 23, 29, 0, "LiJaB");
      dpd_buf4_sort(&L2, CC_GLG, prqs, 30, 20, "LiaJB");
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, G_irr, 20, 30, 20, 30, 0, "LIAjb");
      dpd_buf4_sort(&L2, CC_GLG, psrq, 24, 27, "LIbjA");
      dpd_buf4_sort(&L2, CC_GLG, rqps, 27, 24, "LjAIb");
      dpd_buf4_close(&L2);
    }
  }

  if (!params.ground) { /* sort entries in CC_GL */
    if ( (params.ref==0) || (params.ref==1) ) { /** RHF/ROHF **/
      /* Build L2iJaB list */
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GL, qpsr, 0, 5, "LiJaB");
      dpd_buf4_sort(&L2, CC_GL, pqsr, 0, 5, "LIjaB");
      dpd_buf4_sort(&L2, CC_GL, qprs, 0, 5, "LiJAb");
      dpd_buf4_close(&L2);
      /* Build L2IAJB List */
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
      dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LIAJB");
      dpd_buf4_close(&L2);
      /* Build L2iajb List */
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
      dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "Liajb");
      dpd_buf4_close(&L2);
      /* Build L2IAjb List */
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LIAjb");
      dpd_buf4_close(&L2);
      /* Build L2iaJB List */
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
      dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LiaJB");
      dpd_buf4_close(&L2);
      /* Build L2IbjA and L2 jAIb Lists */
      dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
      dpd_buf4_sort(&L2, CC_GL, psrq, 10, 10, "LIbjA");
      dpd_buf4_sort(&L2, CC_GL, rqps, 10, 10, "LjAIb");
      dpd_buf4_close(&L2);
    }
    else if(params.ref == 2) { /** UHF **/
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GL, qpsr, 23, 29, "LiJaB");
      dpd_buf4_sort(&L2, CC_GL, pqsr, 22, 29, "LIjaB");
      dpd_buf4_sort(&L2, CC_GL, qprs, 23, 28, "LiJAb");
      dpd_buf4_close(&L2);
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
      dpd_buf4_sort(&L2, CC_GL, prqs, 20, 20, "LIAJB");
      dpd_buf4_close(&L2);
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
      dpd_buf4_sort(&L2, CC_GL, prqs, 30, 30, "Liajb");
      dpd_buf4_close(&L2);
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_sort(&L2, CC_GL, prqs, 20, 30, "LIAjb");
      dpd_buf4_close(&L2);
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
      dpd_buf4_sort(&L2, CC_GL, prqs, 30, 20, "LiaJB");
      dpd_buf4_close(&L2);
  
      dpd_buf4_init(&L2, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
      dpd_buf4_sort(&L2, CC_GL, psrq, 24, 27, "LIbjA");
      dpd_buf4_sort(&L2, CC_GL, rqps, 27, 24, "LjAIb");
      dpd_buf4_close(&L2);
    }
  }

  if (!rho_params.R_ground) { /* put copy of R into CC_GR */
    sprintf(R1A_lbl,"RIA %d %d", R_irr, R_root);
    sprintf(R1B_lbl,"Ria %d %d", R_irr, R_root);
    sprintf(R2AA_lbl,"RIJAB %d %d", R_irr, R_root );
    sprintf(R2BB_lbl,"Rijab %d %d", R_irr, R_root );
    sprintf(R2AB_lbl,"RIjAb %d %d", R_irr, R_root );
    sprintf(R2RHF_lbl,"2RIjAb - RIjbA %d %d", R_irr, R_root );
    if (params.ref == 0) {
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, R2RHF_lbl);
      dpd_buf4_copy(&R2, CC_GR, "2RIjAb - RIjbA");
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
      dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "2RIjAb - RIjbA (IA,jb)");
      dpd_buf4_close(&R2);
    }
    if ( (params.ref == 0) || (params.ref == 1) ) {
      dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, R1A_lbl);
      dpd_file2_copy(&R1, CC_GR, "RIA");
      dpd_file2_close(&R1);
      dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, R1B_lbl);
      dpd_file2_copy(&R1, CC_GR, "Ria");
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 2, 7, 2, 7, 0, R2AA_lbl);
      dpd_buf4_copy(&R2, CC_GR, "RIJAB");
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      dpd_buf4_copy(&R2, CC_GR, "Rijab");
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      dpd_buf4_copy(&R2, CC_GR, "RIjAb");
      dpd_buf4_close(&R2);
    }
    else if (params.ref == 2) {
      dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, R1A_lbl);
      dpd_file2_copy(&R1, CC_GR, "RIA");
      dpd_file2_close(&R1);
      dpd_file2_init(&R1, CC_RAMPS, R_irr, 2, 3, R1B_lbl);
      dpd_file2_copy(&R1, CC_GR, "Ria");
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 2, 7, 2, 7, 0, R2AA_lbl);
      dpd_buf4_copy(&R2, CC_GR, "RIJAB");
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      dpd_buf4_copy(&R2, CC_GR, "Rijab");
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, R_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      dpd_buf4_copy(&R2, CC_GR, "RIjAb");
      dpd_buf4_close(&R2);
    }
    /* Now sort copy of R in CC_GR */
    if( (params.ref==0) || (params.ref==1) ) { /** RHF/ROHF **/
      /* Build R2iJaB list */
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
      dpd_buf4_sort(&R2, CC_GR, qpsr, 0, 5, "RiJaB");
      dpd_buf4_sort(&R2, CC_GR, pqsr, 0, 5, "RIjaB");
      dpd_buf4_sort(&R2, CC_GR, qprs, 0, 5, "RiJAb");
      dpd_buf4_close(&R2);
      /* Build R2IAJB List */
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 2, 7, 0, "RIJAB");
      dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RIAJB");
      dpd_buf4_close(&R2);
      /* Build R2iajb List */
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 2, 7, 0, "Rijab");
      dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "Riajb");
      dpd_buf4_close(&R2);
      /* Build R2IAjb List */
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
      dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RIAjb");
      dpd_buf4_close(&R2);
      /* Build R2iaJB List */
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
      dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RiaJB");
      dpd_buf4_close(&R2);
      /* Build R2IbjA and R2 jAIb Lists */
      dpd_buf4_init(&R2, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
      dpd_buf4_sort(&R2, CC_GR, psrq, 10, 10, "RIbjA");
      dpd_buf4_sort(&R2, CC_GR, rqps, 10, 10, "RjAIb");
      dpd_buf4_close(&R2);
    }
    else if(params.ref == 2) { /** UHF **/
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
      dpd_buf4_sort(&R2, CC_GR, qpsr, 23, 29, "RiJaB");
      dpd_buf4_sort(&R2, CC_GR, pqsr, 22, 29, "RIjaB");
      dpd_buf4_sort(&R2, CC_GR, qprs, 23, 28, "RiJAb");
      dpd_buf4_close(&R2);
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 2, 7, 0, "RIJAB");
      dpd_buf4_sort(&R2, CC_GR, prqs, 20, 20, "RIAJB");
      dpd_buf4_close(&R2);
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 10, 15, 12, 17, 0, "Rijab");
      dpd_buf4_sort(&R2, CC_GR, prqs, 30, 30, "Riajb");
      dpd_buf4_close(&R2);
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
      dpd_buf4_sort(&R2, CC_GR, prqs, 20, 30, "RIAjb");
      dpd_buf4_close(&R2);
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
      dpd_buf4_sort(&R2, CC_GR, prqs, 30, 20, "RiaJB");
      dpd_buf4_close(&R2);
  
      dpd_buf4_init(&R2, CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
      dpd_buf4_sort(&R2, CC_GR, psrq, 24, 27, "RIbjA");
      dpd_buf4_sort(&R2, CC_GR, rqps, 27, 24, "RjAIb");
      dpd_buf4_close(&R2);
    }
  }

  return;
}

}} // namespace psi::ccdensity
