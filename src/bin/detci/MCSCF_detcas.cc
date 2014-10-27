/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup DETCAS
    \brief Orbital optimizer for detci
**
** DETCAS
**
** Orbital rotation program for determinant configuration interaction
** wavefunctions evaluated using the DETCI program
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
** Modified 05/28/99 by CDS to move BLAS interface to libqt
**
*/

#include "globaldefs.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
// #include <libipv1/ip_lib.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/wavefunction.h>
#include <libmints/molecule.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "psifiles.h"
#include "psi4-dec.h"
#include "structs.h"
#define EXTERN
#include "globals.h"
#include "MCSCF_setup_io.h"
#include "MCSCF_indpairs.h"

namespace psi { namespace detci {

extern void mcscf_get_mo_info(Options& options);
extern void get_parameters(Options &options);
extern void print_parameters(void);
extern void read_integrals(void);
extern void read_density_matrices(Options& options);
extern void read_lagrangian(void);
extern void form_independent_pairs(void);
extern void read_thetas(int npairs);
extern void write_thetas(int npairs);
extern int  read_ref_orbs(void);
extern int  write_ref_orbs(void);
extern void read_cur_orbs(void);
extern void form_F_act(void);
extern int  diis(int veclen, double *vec, double *errvec);
extern void mcscf_get_mat_block(double **src, double **dst, int dst_dim,
                          int dst_offset, int *dst2src);
extern void calc_dE_dT(int n, double **dEU, int npairs, int *ppair, 
                       int *qpair, double *theta, double *dET);
extern void form_appx_diag_mo_hess(int npairs, int *ppair, int *qpair, 
                              double *F_core, double *tei, double **opdm, 
                              double *tpdm, double *F_act, int firstact, 
                              int lastact, double *hess);
extern void form_diag_mo_hess(int npairs, int *ppair, int *qpair, 
                              double *F_core, double *tei, double **opdm, 
                              double *tpdm, double *F_act, int firstact, 
                              int lastact, double *hess);
extern void form_full_mo_hess(int npairs, int *ppair, int *qpair, 
                         double *oei, double *tei, double **opdm, 
                         double *tpdm, double **lag, double **hess);
extern void form_diag_mo_hess_yy(int npairs, int *ppair, int *qpair, 
                         double *oei, double *tei, double **opdm, 
                         double *tpdm, double **lag, double *hess);
extern void calc_orb_step(int npairs, double *grad, double *hess_diag, 
                          double *theta);
extern void calc_orb_step_full(int npairs, double *grad, double **hess, 
                          double *theta);
extern void calc_orb_step_bfgs(int npairs, double *grad, double **hess, 
                          double *theta);
extern void print_step(int npairs, int steptype);
extern void postmult_by_U(int irrep, int dim, double **mo_coeffs,
                          int npairs, int *p_arr, int *q_arr, 
                          double *theta_arr);
extern void premult_by_U(int irrep, int dim, double **mo_coeffs,
                         int npairs, int *p_arr, int *q_arr, 
                         double *theta_arr);
extern void postmult_by_exp_R(int irrep, int dim, double **mat,
                              int npairs, int *p_arr, int *q_arr, 
                              double *theta_arr);
extern void cleanup(void);


void mcscf_title(void);
void calc_gradient(void);
void bfgs_hessian(void);
void ds_hessian(void);
void calc_hessian(void);
void scale_gradient(void);
int check_conv(void);
int take_step(void);
void rotate_orbs(void);
double** lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec,
               int nmo, int npop, int print_lvl, int lag_file);

// struct calcinfo MCSCF_CalcInfo;
// struct params MCSCF_Parameters;
// int *ioff;
IndepPairs IndPairs;
// PsiReturnType detcas(Options &options);
PsiReturnType mcscf_update(Options &options);

#define MO_HESS_MIN 1.0E-1

}} // end namespace psi::detci


namespace psi { namespace detci {


PsiReturnType mcscf_update(Options &options)
{
  int converged = 0;
  int num_pairs = 0;
  int steptype = 0;

  MCSCF_Parameters.print_lvl = 1;
  MCSCF_CalcInfo.mo_hess = NULL;
  MCSCF_CalcInfo.mo_hess_diag = NULL;

  if (MCSCF_Parameters.print_lvl) tstart();
  get_parameters(options);     /* get running params (convergence, etc)    */
  mcscf_title();                     /* print program identification             */

  if (MCSCF_Parameters.print_lvl) print_parameters();

  mcscf_get_mo_info(options);               /* read DOCC, SOCC, frozen, nbfso, etc      */
  read_integrals();            /* get the 1 and 2 elec MO integrals        */
  read_density_matrices(options);

  MCSCF_CalcInfo.lag = lagcalc(MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.onel_ints_bare,
                     MCSCF_CalcInfo.twoel_ints, MCSCF_CalcInfo.nmo,
                     MCSCF_CalcInfo.npop, MCSCF_Parameters.print_lvl, PSIF_MO_LAG); 


  read_lagrangian();

  form_independent_pairs();
  num_pairs = IndPairs.get_num_pairs();

  read_thetas(num_pairs);
  if (MCSCF_Parameters.print_lvl > 2)
    IndPairs.print_vec(MCSCF_CalcInfo.theta_cur, "\n\tRotation Angles:");

  if (!read_ref_orbs()) {
    read_cur_orbs();
    write_ref_orbs();
    zero_arr(MCSCF_CalcInfo.theta_cur, num_pairs);
    write_thetas(num_pairs);
  }

  form_F_act();
  calc_gradient();
  converged = check_conv();

  if (MCSCF_Parameters.bfgs)
    bfgs_hessian();
  else if (MCSCF_Parameters.ds_hessian)
    ds_hessian();
  else
    calc_hessian();

  scale_gradient();

  if (!converged) {
    steptype = take_step();
    rotate_orbs();
    write_thetas(num_pairs);
  }
  else
    steptype = 0;

  print_step(num_pairs, steptype);

//  if (MCSCF_Parameters.print_lvl) quote();
  cleanup();
  //close_io();
  if (MCSCF_Parameters.print_lvl) tstop();

  if (converged){
    return EndLoop; 
  }
  else{
    return Success;
  }
}





/*
** mcscf_title(): Function prints a program identification
*/
void mcscf_title(void)
{
  if (MCSCF_Parameters.print_lvl) {
   outfile->Printf("\n");
   outfile->Printf("*******************************************************\n");
   outfile->Printf("                      D E T C A S \n");
   outfile->Printf("\n");
   outfile->Printf("                   C. David Sherrill\n") ;
   outfile->Printf("                     April 27 1998\n") ;
   outfile->Printf("*******************************************************\n");
   outfile->Printf("\n\n\n");
   }
  else {
   outfile->Printf( 
     "\nD E T C A S: C. David Sherrill, April 27 1998\n");
   }
  //fflush(outfile);
}




void form_independent_pairs(void)
{

  IndPairs.set(MCSCF_CalcInfo.nirreps, MAX_RAS_SPACES, MCSCF_CalcInfo.ras_opi,
               MCSCF_CalcInfo.ras_orbs, MCSCF_CalcInfo.frozen_docc, MCSCF_CalcInfo.fzc_orbs, 
               MCSCF_CalcInfo.rstr_docc, MCSCF_CalcInfo.cor_orbs,
               MCSCF_CalcInfo.rstr_uocc, MCSCF_CalcInfo.vir_orbs,
               MCSCF_CalcInfo.frozen_uocc, MCSCF_CalcInfo.fzv_orbs,
               MCSCF_CalcInfo.ci2relpitz, MCSCF_Parameters.ignore_ras_ras, MCSCF_Parameters.ignore_fz);

  if (MCSCF_Parameters.print_lvl > 3) IndPairs.print();

}


/*!
** calc_gradient()
**
** This function calculates the MO gradient from the MO Lagrangian
**
** \ingroup DETCAS
*/
void calc_gradient(void)
{
  int pair, npair, h, ir_npairs, ir_norbs, offset;
  double *ir_mo_grad, **ir_lag, *ir_theta_cur, value, rms;
  int *parr, *qarr, *ir_ppair, *ir_qpair;
  int p, q;

  npair = IndPairs.get_num_pairs();
  parr  = IndPairs.get_p_ptr();
  qarr  = IndPairs.get_q_ptr();

  MCSCF_CalcInfo.mo_grad = init_array(npair);

  /*
  calc_grad_1(npair, parr, qarr, MCSCF_CalcInfo.lag, MCSCF_CalcInfo.mo_grad);
  calc_grad_2(npair, parr, qarr, MCSCF_CalcInfo.onel_ints, MCSCF_CalcInfo.twoel_ints, 
              MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.F_act, 
              (MCSCF_CalcInfo.num_cor_orbs + MCSCF_CalcInfo.num_fzc_orbs), 
              MCSCF_CalcInfo.npop, MCSCF_CalcInfo.mo_grad); 
  */

  // scratch array for dEdTheta, big enough for any irrep
  ir_mo_grad = init_array(npair);
  
  // calculate dEdU, then dEdTheta
  for (h=0,offset=0; h<MCSCF_CalcInfo.nirreps; h++) {

    // Setup for this irrep
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    ir_norbs = MCSCF_CalcInfo.orbs_per_irr[h];
    if (h>0) offset += MCSCF_CalcInfo.orbs_per_irr[h-1];
    if (!ir_npairs) continue;
    ir_ppair = IndPairs.get_ir_prel_ptr(h);
    ir_qpair = IndPairs.get_ir_qrel_ptr(h);
    ir_lag = block_matrix(ir_norbs, ir_norbs);
    mcscf_get_mat_block(MCSCF_CalcInfo.lag, ir_lag, ir_norbs, offset, MCSCF_CalcInfo.pitz2ci);

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "Irrep %d of lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    ir_theta_cur = IndPairs.get_irrep_vec(h, MCSCF_CalcInfo.theta_cur); 

    // Need to mult the Lagrangian by 2 to get dEdU
    C_DSCAL(ir_norbs*ir_norbs, 2.0, ir_lag[0], 1);

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "Irrep %d of 2 * lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    if (MCSCF_Parameters.use_thetas) {
      // Calc dEdU
      premult_by_U(h, MCSCF_CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
                   ir_ppair, ir_qpair, ir_theta_cur);

      if (MCSCF_Parameters.print_lvl > 3) {
        outfile->Printf( "dE/dU:\n", h);
        print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
      }

      // Calculate dEdTheta
      calc_dE_dT(MCSCF_CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
                 ir_ppair, ir_qpair, ir_theta_cur, ir_mo_grad);
    }

    /* non-theta version */
    else {
      for (pair=0; pair<ir_npairs; pair++) {
        p = ir_ppair[pair];  q = ir_qpair[pair];
        ir_mo_grad[pair] = ir_lag[p][q] - ir_lag[q][p];
      }
    }

    // Put dEdTheta into the large gradient array
    IndPairs.put_irrep_vec(h, ir_mo_grad, MCSCF_CalcInfo.mo_grad);
    delete [] ir_theta_cur;
    free_block(ir_lag); 
  }


  rms = 0.0;
  for (pair=0; pair<npair; pair++) {
    value =  MCSCF_CalcInfo.mo_grad[pair];
    rms += value * value;
  }

  if (MCSCF_Parameters.print_lvl > 2) 
    IndPairs.print_vec(MCSCF_CalcInfo.mo_grad, "\n\tOrbital Gradient:");

  rms = sqrt(rms);
  MCSCF_CalcInfo.mo_grad_rms = rms;

  if (MCSCF_Parameters.print_lvl)
    outfile->Printf( "\n\tRMS Orbital Gradient: %6.4E\n", rms);
}


/*!
** bfgs_hessian()
**
** This function calculates a BFGS-updated MO/MO Hessian
**
** C. David Sherrill
** March 2004
**
** \ingroup DETCAS
*/
void bfgs_hessian(void)
{
  int i, j, npairs;
  double fac, fad, fae, tval;
  double *dx, *dg, *hdg, **hess_copy, **hess_copy2, hess_det;
  int *idx;

  npairs = IndPairs.get_num_pairs();

  psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);

  /* If no Hessian in the file */
  if (psio_tocscan(PSIF_DETCAS, "Hessian Inverse") == NULL) {
    calc_hessian();
    if (MCSCF_Parameters.hessian == "FULL") {
      MCSCF_CalcInfo.mo_hess = block_matrix(npairs,npairs);
      for (i=0; i<npairs; i++) {
        MCSCF_CalcInfo.mo_hess[i][i] = 1.0 / MCSCF_CalcInfo.mo_hess_diag[i];
      }
    }
    else {
      hess_copy = block_matrix(npairs, npairs);
      idx = init_int_array(npairs);
      for (i=0; i<npairs; i++) {
        for (j=0; j<npairs; j++) {
          hess_copy[i][j] = MCSCF_CalcInfo.mo_hess[i][j];     
        }
      }
      ludcmp(hess_copy,npairs,idx,&hess_det);

      for (i=0;i<npairs;i++) hess_det *= hess_copy[i][i];
      outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

      if (MCSCF_Parameters.level_shift) {
        while (hess_det < MCSCF_Parameters.determ_min) {
          outfile->Printf( "Level shifting hessian by %8.3E\n", MCSCF_Parameters.shift);
          for (i=0; i<npairs; i++) {
            MCSCF_CalcInfo.mo_hess[i][i] += MCSCF_Parameters.shift;
            for (j=0; j<npairs; j++) {
              hess_copy[i][j] = MCSCF_CalcInfo.mo_hess[i][j];
            }
          }
          ludcmp(hess_copy,npairs,idx,&hess_det);
          for (i=0;i<npairs;i++) hess_det *= hess_copy[i][i];
          outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

          for (i=0; i<npairs; i++) {
            for (j=0; j<npairs; j++) {
              hess_copy[i][j] = MCSCF_CalcInfo.mo_hess[i][j];
            }
          }
        }
      }

      /* n.b. this destroys hess_copy, and mo_hess now is INVERSE Hess */
      invert_matrix(hess_copy,MCSCF_CalcInfo.mo_hess,npairs,"outfile");

      if (MCSCF_Parameters.print_lvl > 3) {
        outfile->Printf( "\nInitial MO Hessian Inverse:\n");
        print_mat(MCSCF_CalcInfo.mo_hess,npairs,npairs,"outfile");
        outfile->Printf( "\n");
      }

      free(idx);
      free_block(hess_copy);
    } 

    /* write Hessian */
    psio_write_entry(PSIF_DETCAS, "Hessian Inverse", 
      (char *) MCSCF_CalcInfo.mo_hess[0], npairs*npairs*sizeof(double));
    /* write thetas */
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) MCSCF_CalcInfo.theta_cur,
      npairs*sizeof(double));
    /* write gradient */  
    psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) MCSCF_CalcInfo.mo_grad,
      npairs*sizeof(double));
    psio_close(PSIF_DETCAS, 1);


    return;
  } /* end initialization of BFGS data */

  dx = init_array(npairs);
  dg = init_array(npairs);
  hdg = init_array(npairs);

  /* read previous Hessian */
  MCSCF_CalcInfo.mo_hess = block_matrix(npairs,npairs);
  psio_read_entry(PSIF_DETCAS, "Hessian Inverse", (char *) MCSCF_CalcInfo.mo_hess[0], 
    npairs*npairs*sizeof(double));

  /* read previous thetas */
  psio_read_entry(PSIF_DETCAS, "Thetas", (char *) dx, npairs*sizeof(double));

  /* read previous gradient */
  psio_read_entry(PSIF_DETCAS, "MO Gradient", (char *) dg,
    npairs*sizeof(double));

  /* compute updated Hessian by BFGS procedure, see Numerical Recipies in C */
  
  /* get difference in thetas and gradient */
  for (i=0; i<npairs; i++) {
    dx[i] = MCSCF_CalcInfo.theta_cur[i] - dx[i];
    dg[i] = MCSCF_CalcInfo.mo_grad[i] - dg[i];
  }

  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf( "Delta Theta and Delta Grad arrays:\n");
    for (i=0; i<npairs; i++) {
      outfile->Printf( "%12.7lf   %12.7lf\n", dx[i], dg[i]);
    }
  }

  /* hdg = H^-1 * (grad(i+1) - grad(i)) */
  for (i=0; i<npairs; i++) {
    hdg[i] = 0.0;
    for (j=0; j<npairs; j++) {
      hdg[i] += MCSCF_CalcInfo.mo_hess[i][j] * dg[j];
    }
  }

  fac = fae = 0.0;

  for (i=0; i<npairs; i++) {
    fac += dg[i] * dx[i];
    fae += dg[i] * hdg[i];
  }
  fac = 1.0/fac;
  fad = 1.0/fae;

  /* the dg bit is the diff between BFGS and DFP */
  for (i=0; i<npairs; i++) dg[i] = fac*dx[i] - fad*hdg[i];

  for (i=0; i<npairs; i++) {
    for (j=0; j<npairs; j++) {
      MCSCF_CalcInfo.mo_hess[i][j] += fac*dx[i]*dx[j] - fad*hdg[i]*hdg[j]
         + fae*dg[i]*dg[j]; 
    }
  }

  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf( "\nBFGS MO Hessian Inverse:\n");
    print_mat(MCSCF_CalcInfo.mo_hess,npairs,npairs,"outfile");
    outfile->Printf( "\n");
  }

  /* 
     this doesn't work unless you fix that dg is not overwritten above 
     when that's ensured, it seems to match 
   */
  /*
  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf( "Check of dx = H dg\n");
    for (i=0; i<npairs; i++) {
      tval = 0.0;
      for (j=0; j<npairs; j++) {
        tval += MCSCF_CalcInfo.mo_hess[i][j] * dg[j];
      }
      outfile->Printf( "%12.7lf vs %12.7lf\n", dx[i], tval);
    }
  }
  */


  idx = init_int_array(npairs);
  hess_copy = block_matrix(npairs, npairs);
  for (i=0; i<npairs; i++) {
    for (j=0; j<npairs; j++) {
      hess_copy[i][j] = MCSCF_CalcInfo.mo_hess[i][j];     
    }
  }
  ludcmp(hess_copy,npairs,idx,&hess_det);
  for (i=0; i<npairs; i++) hess_det *= MCSCF_CalcInfo.mo_hess[i][i];
  outfile->Printf( "The determinant of the inverse Hessian is %8.3E\n",
    hess_det);

  /* just debug check */
  if (hess_det < 0.0) {
    hess_copy = block_matrix(npairs, npairs);
    hess_copy2 = block_matrix(npairs, npairs);
    for (i=0; i<npairs; i++) {
      for (j=0; j<npairs; j++) {
        hess_copy[i][j] = MCSCF_CalcInfo.mo_hess[i][j];     
      }
    }

    /* n.b. this destroys hess_copy */
    invert_matrix(hess_copy,hess_copy2,npairs,"outfile");

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "\nMO Hessian:\n");
      print_mat(hess_copy2,npairs,npairs,"outfile");
      outfile->Printf( "\n");
    }

    ludcmp(hess_copy2,npairs,idx,&hess_det);
    for (i=0;i<npairs;i++) hess_det *= hess_copy2[i][i];
    outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

      if (MCSCF_Parameters.level_shift) {
        while (hess_det < MCSCF_Parameters.determ_min) {
          outfile->Printf( "Level shifting hessian by %8.3E\n", MCSCF_Parameters.shift);
          for (i=0; i<npairs; i++) {
            hess_copy2[i][i] += MCSCF_Parameters.shift;
            for (j=0; j<npairs; j++) {
              hess_copy[i][j] = hess_copy2[i][j];
            }
          }
          ludcmp(hess_copy,npairs,idx,&hess_det);
          for (i=0;i<npairs;i++) hess_det *= hess_copy[i][i];
          outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);                                                                                
          for (i=0; i<npairs; i++) {
            for (j=0; j<npairs; j++) {
              hess_copy2[i][j] = hess_copy[i][j];
            }
          }
        }
      }
                                                                                
    /* n.b. this destroys hess_copy2, and mo_hess now is INVERSE Hess */
    invert_matrix(hess_copy2,MCSCF_CalcInfo.mo_hess,npairs,"outfile");

    free_block(hess_copy2);
  }

  free_block(hess_copy);

  /* write thetas */
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) MCSCF_CalcInfo.theta_cur,
    npairs*sizeof(double));
  /* write gradient */  
  psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) MCSCF_CalcInfo.mo_grad,
    npairs*sizeof(double));
  /* write updated Hessian */
  psio_write_entry(PSIF_DETCAS, "Hessian Inverse", (char *) MCSCF_CalcInfo.mo_hess[0],
    npairs*npairs*sizeof(double));
  psio_close(PSIF_DETCAS, 1);


  free(idx);
  free(dx);
  free(dg);
  free(hdg);
}


/*!
** ds_hessian()
**
** This function calculates a Hessian update by a difference of gradients
**
** C. David Sherrill
** March 2004
**
** \ingroup DETCAS
*/
void ds_hessian(void)
{
  int i, npairs;
  double tval;
  double *dx, *dg;

  npairs = IndPairs.get_num_pairs();

  psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);

  /* If no Hessian in the file */
  if (psio_tocscan(PSIF_DETCAS, "Hessian") == NULL) {
    calc_hessian();
    if (MCSCF_Parameters.hessian == "FULL") {
      MCSCF_CalcInfo.mo_hess_diag = init_array(npairs);
      for (i=0; i<npairs; i++) {
        MCSCF_CalcInfo.mo_hess_diag[i] = MCSCF_CalcInfo.mo_hess[i][i];
      }
    }
    for (i=0; i<npairs; i++) {
      if (MCSCF_CalcInfo.mo_hess_diag[i] < MO_HESS_MIN) {
        outfile->Printf( "Warning: MO Hessian denominator too small\n");
        MCSCF_CalcInfo.mo_hess_diag[i] = MO_HESS_MIN;
      }
    }

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "\nInitial MO Hessian:\n");
      for (i=0; i<npairs; i++) 
        outfile->Printf( "%12.6lf\n", MCSCF_CalcInfo.mo_hess_diag[i]);
      outfile->Printf( "\n");
    }

    /* write Hessian */
    psio_write_entry(PSIF_DETCAS, "Hessian", 
      (char *) MCSCF_CalcInfo.mo_hess_diag, npairs*sizeof(double));
    /* write thetas */
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) MCSCF_CalcInfo.theta_cur,
      npairs*sizeof(double));
    /* write gradient */  
    psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) MCSCF_CalcInfo.mo_grad,
      npairs*sizeof(double));
    psio_close(PSIF_DETCAS, 1);
    return;
  } /* end initialization of data */

  dx = init_array(npairs);
  dg = init_array(npairs);

  /* read previous Hessian */
  MCSCF_CalcInfo.mo_hess_diag = init_array(npairs);
  psio_read_entry(PSIF_DETCAS, "Hessian", (char *) MCSCF_CalcInfo.mo_hess_diag, 
    npairs*sizeof(double));

  /* read previous thetas */
  psio_read_entry(PSIF_DETCAS, "Thetas", (char *) dx, npairs*sizeof(double));

  /* read previous gradient */
  psio_read_entry(PSIF_DETCAS, "MO Gradient", (char *) dg,
    npairs*sizeof(double));

  /* compute updated Hessian by David's recipie */
  
  /* get difference in thetas and gradient */
  for (i=0; i<npairs; i++) {
    dx[i] = MCSCF_CalcInfo.theta_cur[i] - dx[i];
    dg[i] = MCSCF_CalcInfo.mo_grad[i] - dg[i];
  }

  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf( "Delta Theta and Delta Grad arrays:\n");
    for (i=0; i<npairs; i++) {
      outfile->Printf( "%12.7lf   %12.7lf\n", dx[i], dg[i]);
    }
  }

  /* H = 1/2 [ H + (grad(i+1) - grad(i))/(x(i+1)-x(i)) ] */
  for (i=0; i<npairs; i++) {
    tval = dg[i] / dx[i];
    if (tval < - MO_HESS_MIN) tval = - MO_HESS_MIN;
    if (tval > 500.0) tval = 500.0;
    MCSCF_CalcInfo.mo_hess_diag[i] = 0.5 * (MCSCF_CalcInfo.mo_hess_diag[i] + tval);
    if (MCSCF_CalcInfo.mo_hess_diag[i] < MO_HESS_MIN) 
      MCSCF_CalcInfo.mo_hess_diag[i] = MO_HESS_MIN;
  }

  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf( "\nDS MO Hessian:\n");
    for (i=0; i<npairs; i++) 
      outfile->Printf( "%12.6lf\n", MCSCF_CalcInfo.mo_hess_diag[i]);
    outfile->Printf( "\n");
  }

  /* write thetas */
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) MCSCF_CalcInfo.theta_cur,
    npairs*sizeof(double));
  /* write gradient */  
  psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) MCSCF_CalcInfo.mo_grad,
    npairs*sizeof(double));
  /* write updated Hessian */
  psio_write_entry(PSIF_DETCAS, "Hessian", (char *) MCSCF_CalcInfo.mo_hess_diag,
    npairs*sizeof(double));
  psio_close(PSIF_DETCAS, 1);


  free(dx);
  free(dg);
}



/*!
** calc_hessian()
**
** This function calculates an approximate MO Hessian from the 
** Fock matrix intermediates and two-electron integrals.
**
** C. David Sherrill
** April 1998
**
** \ingroup DETCAS
*/
void calc_hessian(void)
{

  int npairs, *ppair, *qpair, ncore;

  npairs = IndPairs.get_num_pairs();
  ppair = IndPairs.get_p_ptr();
  qpair = IndPairs.get_q_ptr();

  /* Now calculate the approximate diagonal MO Hessian */
  ncore = MCSCF_CalcInfo.num_fzc_orbs + MCSCF_CalcInfo.num_cor_orbs;

  if (MCSCF_Parameters.hessian == "DIAG") {
    MCSCF_CalcInfo.mo_hess_diag = init_array(npairs);
    
    if (MCSCF_Parameters.use_fzc_h == 1) 
      form_diag_mo_hess(npairs, ppair, qpair, MCSCF_CalcInfo.onel_ints, 
        MCSCF_CalcInfo.twoel_ints, MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.F_act, 
        ncore, MCSCF_CalcInfo.npop, MCSCF_CalcInfo.mo_hess_diag);
    else
      form_diag_mo_hess_yy(npairs, ppair, qpair, MCSCF_CalcInfo.onel_ints, 
        MCSCF_CalcInfo.twoel_ints, MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.lag,
        MCSCF_CalcInfo.mo_hess_diag);

    if (MCSCF_Parameters.print_lvl > 3)  
      IndPairs.print_vec(MCSCF_CalcInfo.mo_hess_diag,"\n\tDiagonal MO Hessian:");
  }
  else if (MCSCF_Parameters.hessian == "APPROX_DIAG") {
    MCSCF_CalcInfo.mo_hess_diag = init_array(npairs);
    form_appx_diag_mo_hess(npairs, ppair, qpair, MCSCF_CalcInfo.onel_ints, 
                      MCSCF_CalcInfo.twoel_ints, MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, 
                      MCSCF_CalcInfo.F_act, ncore, MCSCF_CalcInfo.npop, 
                      MCSCF_CalcInfo.mo_hess_diag);
    if (MCSCF_Parameters.print_lvl > 3)  
      IndPairs.print_vec(MCSCF_CalcInfo.mo_hess_diag,"\n\tAppx Diagonal MO Hessian:");
  }
  else if (MCSCF_Parameters.hessian == "FULL") {
    MCSCF_CalcInfo.mo_hess = block_matrix(npairs,npairs);
    form_full_mo_hess(npairs, ppair, qpair, MCSCF_CalcInfo.onel_ints, 
      MCSCF_CalcInfo.twoel_ints, MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.lag,
      MCSCF_CalcInfo.mo_hess);
    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "\nMO Hessian:\n");
      print_mat(MCSCF_CalcInfo.mo_hess,npairs,npairs,"outfile");
      outfile->Printf( "\n");
    }
  }
  else {
    outfile->Printf( "(detcas): Unrecognized Hessian option %s\n", 
      MCSCF_Parameters.hessian.c_str());
  }
 

}



/*!
** scale_gradient()
**
** Scales the orbital gradient by the approximate orbital Hessian
**
** \ingroup DETCAS
*/
void scale_gradient(void)
{
  int pair, npairs;
  double rms, value;

  npairs = IndPairs.get_num_pairs();

  MCSCF_CalcInfo.theta_step = init_array(npairs);

  // All this actually does is scale the gradient by the Hessian
  // If we have a diagonal (exact or approximate) Hessian

  // BFGS
  if (MCSCF_Parameters.scale_grad && MCSCF_Parameters.bfgs) {
    calc_orb_step_bfgs(npairs, MCSCF_CalcInfo.mo_grad, MCSCF_CalcInfo.mo_hess,
      MCSCF_CalcInfo.theta_step);
  }
  // non-BFGS diagonal Hessian
  else if (MCSCF_Parameters.scale_grad && ((MCSCF_Parameters.hessian == "DIAG") || 
       (MCSCF_Parameters.hessian == "APPROX_DIAG"))) {
    calc_orb_step(npairs, MCSCF_CalcInfo.mo_grad, MCSCF_CalcInfo.mo_hess_diag,
      MCSCF_CalcInfo.theta_step);
  }
  // non-BFGS full Hessian
  else if ((MCSCF_Parameters.scale_grad && (MCSCF_Parameters.hessian == "FULL")) ||
    MCSCF_Parameters.bfgs) {
    calc_orb_step_full(npairs, MCSCF_CalcInfo.mo_grad, MCSCF_CalcInfo.mo_hess,
      MCSCF_CalcInfo.theta_step);
  }
  // No Hessian available: take unit step down gradient direction
  else {
    for (pair=0; pair<npairs; pair++) 
      MCSCF_CalcInfo.theta_step[pair] = MCSCF_CalcInfo.mo_grad[pair];
  }
    

  if (MCSCF_Parameters.print_lvl > 3)  
    IndPairs.print_vec(MCSCF_CalcInfo.theta_step,"\n\tScaled Orbital Grad:");

  rms = 0.0;
  for (pair=0; pair<npairs; pair++) {
    value =  MCSCF_CalcInfo.theta_step[pair];
    rms += value * value;
  }

  rms = sqrt(rms);
  MCSCF_CalcInfo.scaled_mo_grad_rms = rms;
 
  if (MCSCF_Parameters.print_lvl)
    outfile->Printf( "\n\tScaled RMS Orbital Gradient: %6.4E\n", rms);

  if (MCSCF_Parameters.scale_step != 1.0) {
    for (pair=0; pair<npairs; pair++) 
      MCSCF_CalcInfo.theta_step[pair] *= MCSCF_Parameters.scale_step;
  }

}


/*!
** take_step()
**
** This function takes a step in orbital rotation (theta) space
**
** Returns: type of step taken; 1=regular (Newton-Raphson), 2=diis
**
** \ingroup DETCAS
*/
int take_step(void)
{
  int npairs, pair, took_diis;
  
  npairs = IndPairs.get_num_pairs();

  /* for debugging purposes */
  if (MCSCF_Parameters.force_step) {
    MCSCF_CalcInfo.theta_cur[MCSCF_Parameters.force_pair] = MCSCF_Parameters.force_value;
    outfile->Printf( "Forcing step for pair %d of size %8.3E\n",
      MCSCF_Parameters.force_pair, MCSCF_Parameters.force_value);
    return(1);
  }

  for (pair=0; pair<npairs; pair++) 
    MCSCF_CalcInfo.theta_cur[pair] += MCSCF_CalcInfo.theta_step[pair];
    //MCSCF_CalcInfo.theta_cur[pair] = MCSCF_CalcInfo.theta_step[pair];

  if (MCSCF_CalcInfo.iter >= MCSCF_Parameters.diis_start) 
    took_diis = diis(npairs, MCSCF_CalcInfo.theta_cur, MCSCF_CalcInfo.theta_step);
    //took_diis = diis(npairs, MCSCF_CalcInfo.theta_cur, MCSCF_CalcInfo.mo_grad);
  else 
    took_diis = 0;

  if (!took_diis) {
    if (MCSCF_Parameters.print_lvl) 
      outfile->Printf( "Taking regular step\n");
  }

  return(took_diis+1);

}


/*!
**
** rotate_orbs()
**
** Rotate the orbitals, irrep by irrep
**
** \ingroup DETCAS
*/
void rotate_orbs(void)
{
  double *ir_theta;
  int h, pair, ir_norbs, ir_npairs, *ir_ppair, *ir_qpair;

  // First, we need to come up with Theta vectors for each irrep
  ir_theta = init_array(IndPairs.get_num_pairs());  // always big enough
  for (h=0; h<MCSCF_CalcInfo.nirreps; h++) {
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    if (ir_npairs) {
      ir_norbs = MCSCF_CalcInfo.orbs_per_irr[h];
      ir_theta = IndPairs.get_irrep_vec(h, MCSCF_CalcInfo.theta_cur);
      ir_ppair  = IndPairs.get_ir_prel_ptr(h);
      ir_qpair  = IndPairs.get_ir_qrel_ptr(h);
      
      if (MCSCF_Parameters.print_lvl > 3) {
        outfile->Printf( "Thetas for irrep %d\n", h);
        for (pair=0; pair<ir_npairs; pair++) {
          outfile->Printf( "Pair (%2d,%2d) = %12.6lf\n",
                  ir_ppair[pair], ir_qpair[pair], ir_theta[pair]);
        }
        outfile->Printf( "\n");
        //fflush(outfile);
      }

      /* print old coefficients */
      if (MCSCF_Parameters.print_mos) {
        outfile->Printf( "\n\tOld molecular orbitals for irrep %s\n", 
          MCSCF_CalcInfo.labels[h]);
        print_mat(MCSCF_CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, "outfile");
      }

      if (MCSCF_Parameters.use_thetas)
        postmult_by_U(h, ir_norbs, MCSCF_CalcInfo.mo_coeffs[h], ir_npairs, 
          ir_ppair, ir_qpair, ir_theta);
      else
        postmult_by_exp_R(h, ir_norbs, MCSCF_CalcInfo.mo_coeffs[h], ir_npairs,
          ir_ppair, ir_qpair, ir_theta);

      /* print new coefficients */
      if (MCSCF_Parameters.print_mos) {
        outfile->Printf( "\n\tNew molecular orbitals for irrep %s\n", 
          MCSCF_CalcInfo.labels[h]);
        print_mat(MCSCF_CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, "outfile");
      }


      /* write the new block of MO coefficients to file30 */
      chkpt_init(PSIO_OPEN_OLD);
      chkpt_wt_scf_irrep(MCSCF_CalcInfo.mo_coeffs[h], h);
      chkpt_close();
      delete [] ir_theta;
    }
  }


}


/*!
** check_conv
**
** Check the summary file to see if we've converged
**
** Returns: 1 if converged, otherwise 0
**
** \ingroup DETCAS
*/
int check_conv(void)
{
  FILE *sumfile;
  char sumfile_name[] = "file14.dat";
  char comment[MAX_COMMENT];
  int i, entries, iter, nind;
  double rmsgrad, scaled_rmsgrad, energy, energy_last;
  int converged_energy=0, converged_grad=0, last_converged=0;

  ffile_noexit(&sumfile,sumfile_name,2);

  if (sumfile == NULL) {
    MCSCF_CalcInfo.iter = 0;
    return(0);
  }

  if (fscanf(sumfile, "%d", &entries) != 1) {
    outfile->Printf("(print_step): Trouble reading num entries in file %s\n",
            sumfile_name);
    fclose(sumfile);
    MCSCF_CalcInfo.iter = 0;
    return(0);
  } 

  MCSCF_CalcInfo.iter = entries;
  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &nind, &scaled_rmsgrad,
           &rmsgrad, &energy_last, comment);
  }
  fclose(sumfile);

  chkpt_init(PSIO_OPEN_OLD);
  energy = chkpt_rd_etot();
  chkpt_close();

  /* check for convergence */
  if (rmsgrad < MCSCF_Parameters.rms_grad_convergence)
    converged_grad = 1;
  if (fabs(energy_last - energy) < MCSCF_Parameters.energy_convergence)
    converged_energy = 1;
  if (strstr(comment, "CONV") != NULL)
    last_converged = 1;

  if (converged_grad && converged_energy && !last_converged) {
    outfile->Printf( "\n\t*** Calculation Converged ***\n");
    return(1);
  }
  else {
    outfile->Printf( "\n\t... calculation continuing ...\n");
    return(0);
  }
}

}} // end namespace psi::detci

