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
//#include <libipv1/ip_lib.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/wavefunction.h>
#include <libmints/molecule.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "psifiles.h"
#include "psi4-dec.h"
#include "globals.h"
#include "setup_io.h"
#include "indpairs.h"
#include "caswave.h"

namespace psi { namespace detcas {

extern void get_mo_info(Options& options);
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
extern void get_mat_block(double **src, double **dst, int dst_dim,
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

void title(void);
void quote(void);
void init_ioff(void);
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

struct calcinfo CalcInfo;
struct params Params;
int *ioff;
IndepPairs IndPairs;
PsiReturnType detcas(Options &options);

#define MO_HESS_MIN 1.0E-1

}} // end namespace psi::detcas

/* GLOBAL VARIABLES (other modules load these via globals.h) */
extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi { namespace detcas {

CASWavefunction::CASWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    set_reference_wavefunction(reference_wavefunction);
    init();
}

CASWavefunction::~CASWavefunction()
{

}

void CASWavefunction::init()
{
    // TODO-CDS:
    // The CC codes destroy the checkpoint object created by Wavefunction.
    // We'd like to be able to do the same here.  Need to convert everything
    // such that we don't explicitly need checkpoint
    // Destroy it. Otherwise we will see a "file already open" error.
    chkpt_.reset();

    // We're copying this stuff over like it's done in ccenergy, but
    // these Wavefunction member data are not actually used by the code
    // yet (Sept 2011 CDS).  Copying these only in case they're needed
    // by someone else who uses the CASWavefunction and expects them to
    // be available.

    // CDS-TODO: Note, some of these data should be updated to reflect what
    // the CI wavefunction is doing, not what the reference wavefunction
    // is doing, in case these are actually used elsewhere.
    nso_        = reference_wavefunction_->nso();
    nirrep_     = reference_wavefunction_->nirrep();
    nmo_        = reference_wavefunction_->nmo();
    for(int h = 0; h < nirrep_; ++h){
        soccpi_[h] = reference_wavefunction_->soccpi()[h];
        doccpi_[h] = reference_wavefunction_->doccpi()[h];
        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
        nmopi_[h]  = reference_wavefunction_->nmopi()[h];
        nsopi_[h]  = reference_wavefunction_->nsopi()[h];
    }
}

// double CASWavefunction::compute_energy()
// {
//     energy_ = 0.0;
//     PsiReturnType cas_return;
//     if ((cas_return = psi::detcas::detcas(options_)) == Success) {
//         // Get the total energy
//         energy_ = Process::environment.globals["CURRENT ENERGY"];
//     }
//     
//     return energy_;
// 
// }

PsiReturnType CASWavefunction::cas_update()
{
    PsiReturnType cas_return;
    cas_return = psi::detcas::detcas(options_);
    return cas_return;
}


//using namespace psi::detcas;

PsiReturnType detcas(Options &options)
{
  int converged = 0;
  int num_pairs = 0;
  int steptype = 0;

  Params.print_lvl = 1;
  CalcInfo.mo_hess = NULL;
  CalcInfo.mo_hess_diag = NULL;

  if (Params.print_lvl) tstart();
  get_parameters(options);     /* get running params (convergence, etc)    */
  init_ioff();                 /* set up the ioff array                    */
  title();                     /* print program identification             */

  if (Params.print_lvl) print_parameters();

  get_mo_info(options);               /* read DOCC, SOCC, frozen, nbfso, etc      */
  read_integrals();            /* get the 1 and 2 elec MO integrals        */
  read_density_matrices(options);

  CalcInfo.lag = lagcalc(CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.onel_ints_bare,
                     CalcInfo.twoel_ints, CalcInfo.nmo,
                     CalcInfo.npop, Params.print_lvl, PSIF_MO_LAG); 


  read_lagrangian();

  form_independent_pairs();
  num_pairs = IndPairs.get_num_pairs();

  read_thetas(num_pairs);
  if (Params.print_lvl > 2)
    IndPairs.print_vec(CalcInfo.theta_cur, "\n\tRotation Angles:");

  if (!read_ref_orbs()) {
    read_cur_orbs();
    write_ref_orbs();
    zero_arr(CalcInfo.theta_cur, num_pairs);
    write_thetas(num_pairs);
  }

  form_F_act();
  calc_gradient();
  converged = check_conv();

  if (Params.bfgs)
    bfgs_hessian();
  else if (Params.ds_hessian)
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

  if (Params.print_lvl) quote();
  cleanup();
  //close_io();
  if (Params.print_lvl) tstop();

  /* TODO  Fix this! Likely need to return non-Success flag */
  if (converged){
    return EndLoop; 
  }
  else{
    return Success;
  }
}


/*
** init_ioff(): Set up the ioff array for quick indexing
*/
void init_ioff(void)
{
  int i;

  /* set offsets for ij-type canonical ordering */
  ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
  ioff[0] = 0;
  for (i = 1; i < IOFF_MAX ; i++) ioff[i] = ioff[i-1] + i;
}



/*
** title(): Function prints a program identification
*/
void title(void)
{
  if (Params.print_lvl) {
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


void quote(void)
{
  outfile->Printf("\n\t\t \"Good bug ... dead bug\" \n\n");
  outfile->Printf("\t\t\t - Ed Valeev\n\n");
  //fflush(outfile);
}



void form_independent_pairs(void)
{

  IndPairs.set(CalcInfo.nirreps, MAX_RAS_SPACES, CalcInfo.ras_opi,
               CalcInfo.ras_orbs, CalcInfo.frozen_docc, CalcInfo.fzc_orbs, 
               CalcInfo.rstr_docc, CalcInfo.cor_orbs,
               CalcInfo.rstr_uocc, CalcInfo.vir_orbs,
               CalcInfo.frozen_uocc, CalcInfo.fzv_orbs,
               CalcInfo.ci2relpitz, Params.ignore_ras_ras, Params.ignore_fz);

  if (Params.print_lvl > 3) IndPairs.print();

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

  CalcInfo.mo_grad = init_array(npair);

  /*
  calc_grad_1(npair, parr, qarr, CalcInfo.lag, CalcInfo.mo_grad);
  calc_grad_2(npair, parr, qarr, CalcInfo.onel_ints, CalcInfo.twoel_ints, 
              CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.F_act, 
              (CalcInfo.num_cor_orbs + CalcInfo.num_fzc_orbs), 
              CalcInfo.npop, CalcInfo.mo_grad); 
  */

  // scratch array for dEdTheta, big enough for any irrep
  ir_mo_grad = init_array(npair);
  
  // calculate dEdU, then dEdTheta
  for (h=0,offset=0; h<CalcInfo.nirreps; h++) {

    // Setup for this irrep
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    ir_norbs = CalcInfo.orbs_per_irr[h];
    if (h>0) offset += CalcInfo.orbs_per_irr[h-1];
    if (!ir_npairs) continue;
    ir_ppair = IndPairs.get_ir_prel_ptr(h);
    ir_qpair = IndPairs.get_ir_qrel_ptr(h);
    ir_lag = block_matrix(ir_norbs, ir_norbs);
    get_mat_block(CalcInfo.lag, ir_lag, ir_norbs, offset, CalcInfo.pitz2ci);

    if (Params.print_lvl > 3) {
      outfile->Printf( "Irrep %d of lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    ir_theta_cur = IndPairs.get_irrep_vec(h, CalcInfo.theta_cur); 

    // Need to mult the Lagrangian by 2 to get dEdU
    C_DSCAL(ir_norbs*ir_norbs, 2.0, ir_lag[0], 1);

    if (Params.print_lvl > 3) {
      outfile->Printf( "Irrep %d of 2 * lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    if (Params.use_thetas) {
      // Calc dEdU
      premult_by_U(h, CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
                   ir_ppair, ir_qpair, ir_theta_cur);

      if (Params.print_lvl > 3) {
        outfile->Printf( "dE/dU:\n", h);
        print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
      }

      // Calculate dEdTheta
      calc_dE_dT(CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
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
    IndPairs.put_irrep_vec(h, ir_mo_grad, CalcInfo.mo_grad);
    delete [] ir_theta_cur;
    free_block(ir_lag); 
  }


  rms = 0.0;
  for (pair=0; pair<npair; pair++) {
    value =  CalcInfo.mo_grad[pair];
    rms += value * value;
  }

  if (Params.print_lvl > 2) 
    IndPairs.print_vec(CalcInfo.mo_grad, "\n\tOrbital Gradient:");

  rms = sqrt(rms);
  CalcInfo.mo_grad_rms = rms;

  if (Params.print_lvl)
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
    if (Params.hessian == "FULL") {
      CalcInfo.mo_hess = block_matrix(npairs,npairs);
      for (i=0; i<npairs; i++) {
        CalcInfo.mo_hess[i][i] = 1.0 / CalcInfo.mo_hess_diag[i];
      }
    }
    else {
      hess_copy = block_matrix(npairs, npairs);
      idx = init_int_array(npairs);
      for (i=0; i<npairs; i++) {
        for (j=0; j<npairs; j++) {
          hess_copy[i][j] = CalcInfo.mo_hess[i][j];     
        }
      }
      ludcmp(hess_copy,npairs,idx,&hess_det);

      for (i=0;i<npairs;i++) hess_det *= hess_copy[i][i];
      outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

      if (Params.level_shift) {
        while (hess_det < Params.determ_min) {
          outfile->Printf( "Level shifting hessian by %8.3E\n", Params.shift);
          for (i=0; i<npairs; i++) {
            CalcInfo.mo_hess[i][i] += Params.shift;
            for (j=0; j<npairs; j++) {
              hess_copy[i][j] = CalcInfo.mo_hess[i][j];
            }
          }
          ludcmp(hess_copy,npairs,idx,&hess_det);
          for (i=0;i<npairs;i++) hess_det *= hess_copy[i][i];
          outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

          for (i=0; i<npairs; i++) {
            for (j=0; j<npairs; j++) {
              hess_copy[i][j] = CalcInfo.mo_hess[i][j];
            }
          }
        }
      }

      /* n.b. this destroys hess_copy, and mo_hess now is INVERSE Hess */
      invert_matrix(hess_copy,CalcInfo.mo_hess,npairs,"outfile");

      if (Params.print_lvl > 3) {
        outfile->Printf( "\nInitial MO Hessian Inverse:\n");
        print_mat(CalcInfo.mo_hess,npairs,npairs,"outfile");
        outfile->Printf( "\n");
      }

      free(idx);
      free_block(hess_copy);
    } 

    /* write Hessian */
    psio_write_entry(PSIF_DETCAS, "Hessian Inverse", 
      (char *) CalcInfo.mo_hess[0], npairs*npairs*sizeof(double));
    /* write thetas */
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) CalcInfo.theta_cur,
      npairs*sizeof(double));
    /* write gradient */  
    psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) CalcInfo.mo_grad,
      npairs*sizeof(double));
    psio_close(PSIF_DETCAS, 1);


    return;
  } /* end initialization of BFGS data */

  dx = init_array(npairs);
  dg = init_array(npairs);
  hdg = init_array(npairs);

  /* read previous Hessian */
  CalcInfo.mo_hess = block_matrix(npairs,npairs);
  psio_read_entry(PSIF_DETCAS, "Hessian Inverse", (char *) CalcInfo.mo_hess[0], 
    npairs*npairs*sizeof(double));

  /* read previous thetas */
  psio_read_entry(PSIF_DETCAS, "Thetas", (char *) dx, npairs*sizeof(double));

  /* read previous gradient */
  psio_read_entry(PSIF_DETCAS, "MO Gradient", (char *) dg,
    npairs*sizeof(double));

  /* compute updated Hessian by BFGS procedure, see Numerical Recipies in C */
  
  /* get difference in thetas and gradient */
  for (i=0; i<npairs; i++) {
    dx[i] = CalcInfo.theta_cur[i] - dx[i];
    dg[i] = CalcInfo.mo_grad[i] - dg[i];
  }

  if (Params.print_lvl > 3) {
    outfile->Printf( "Delta Theta and Delta Grad arrays:\n");
    for (i=0; i<npairs; i++) {
      outfile->Printf( "%12.7lf   %12.7lf\n", dx[i], dg[i]);
    }
  }

  /* hdg = H^-1 * (grad(i+1) - grad(i)) */
  for (i=0; i<npairs; i++) {
    hdg[i] = 0.0;
    for (j=0; j<npairs; j++) {
      hdg[i] += CalcInfo.mo_hess[i][j] * dg[j];
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
      CalcInfo.mo_hess[i][j] += fac*dx[i]*dx[j] - fad*hdg[i]*hdg[j]
         + fae*dg[i]*dg[j]; 
    }
  }

  if (Params.print_lvl > 3) {
    outfile->Printf( "\nBFGS MO Hessian Inverse:\n");
    print_mat(CalcInfo.mo_hess,npairs,npairs,"outfile");
    outfile->Printf( "\n");
  }

  /* 
     this doesn't work unless you fix that dg is not overwritten above 
     when that's ensured, it seems to match 
   */
  /*
  if (Params.print_lvl > 3) {
    outfile->Printf( "Check of dx = H dg\n");
    for (i=0; i<npairs; i++) {
      tval = 0.0;
      for (j=0; j<npairs; j++) {
        tval += CalcInfo.mo_hess[i][j] * dg[j];
      }
      outfile->Printf( "%12.7lf vs %12.7lf\n", dx[i], tval);
    }
  }
  */


  idx = init_int_array(npairs);
  hess_copy = block_matrix(npairs, npairs);
  for (i=0; i<npairs; i++) {
    for (j=0; j<npairs; j++) {
      hess_copy[i][j] = CalcInfo.mo_hess[i][j];     
    }
  }
  ludcmp(hess_copy,npairs,idx,&hess_det);
  for (i=0; i<npairs; i++) hess_det *= CalcInfo.mo_hess[i][i];
  outfile->Printf( "The determinant of the inverse Hessian is %8.3E\n",
    hess_det);

  /* just debug check */
  if (hess_det < 0.0) {
    hess_copy = block_matrix(npairs, npairs);
    hess_copy2 = block_matrix(npairs, npairs);
    for (i=0; i<npairs; i++) {
      for (j=0; j<npairs; j++) {
        hess_copy[i][j] = CalcInfo.mo_hess[i][j];     
      }
    }

    /* n.b. this destroys hess_copy */
    invert_matrix(hess_copy,hess_copy2,npairs,"outfile");

    if (Params.print_lvl > 3) {
      outfile->Printf( "\nMO Hessian:\n");
      print_mat(hess_copy2,npairs,npairs,"outfile");
      outfile->Printf( "\n");
    }

    ludcmp(hess_copy2,npairs,idx,&hess_det);
    for (i=0;i<npairs;i++) hess_det *= hess_copy2[i][i];
    outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

      if (Params.level_shift) {
        while (hess_det < Params.determ_min) {
          outfile->Printf( "Level shifting hessian by %8.3E\n", Params.shift);
          for (i=0; i<npairs; i++) {
            hess_copy2[i][i] += Params.shift;
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
    invert_matrix(hess_copy2,CalcInfo.mo_hess,npairs,"outfile");

    free_block(hess_copy2);
  }

  free_block(hess_copy);

  /* write thetas */
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) CalcInfo.theta_cur,
    npairs*sizeof(double));
  /* write gradient */  
  psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) CalcInfo.mo_grad,
    npairs*sizeof(double));
  /* write updated Hessian */
  psio_write_entry(PSIF_DETCAS, "Hessian Inverse", (char *) CalcInfo.mo_hess[0],
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
    if (Params.hessian == "FULL") {
      CalcInfo.mo_hess_diag = init_array(npairs);
      for (i=0; i<npairs; i++) {
        CalcInfo.mo_hess_diag[i] = CalcInfo.mo_hess[i][i];
      }
    }
    for (i=0; i<npairs; i++) {
      if (CalcInfo.mo_hess_diag[i] < MO_HESS_MIN) {
        outfile->Printf( "Warning: MO Hessian denominator too small\n");
        CalcInfo.mo_hess_diag[i] = MO_HESS_MIN;
      }
    }

    if (Params.print_lvl > 3) {
      outfile->Printf( "\nInitial MO Hessian:\n");
      for (i=0; i<npairs; i++) 
        outfile->Printf( "%12.6lf\n", CalcInfo.mo_hess_diag[i]);
      outfile->Printf( "\n");
    }

    /* write Hessian */
    psio_write_entry(PSIF_DETCAS, "Hessian", 
      (char *) CalcInfo.mo_hess_diag, npairs*sizeof(double));
    /* write thetas */
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) CalcInfo.theta_cur,
      npairs*sizeof(double));
    /* write gradient */  
    psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) CalcInfo.mo_grad,
      npairs*sizeof(double));
    psio_close(PSIF_DETCAS, 1);
    return;
  } /* end initialization of data */

  dx = init_array(npairs);
  dg = init_array(npairs);

  /* read previous Hessian */
  CalcInfo.mo_hess_diag = init_array(npairs);
  psio_read_entry(PSIF_DETCAS, "Hessian", (char *) CalcInfo.mo_hess_diag, 
    npairs*sizeof(double));

  /* read previous thetas */
  psio_read_entry(PSIF_DETCAS, "Thetas", (char *) dx, npairs*sizeof(double));

  /* read previous gradient */
  psio_read_entry(PSIF_DETCAS, "MO Gradient", (char *) dg,
    npairs*sizeof(double));

  /* compute updated Hessian by David's recipie */
  
  /* get difference in thetas and gradient */
  for (i=0; i<npairs; i++) {
    dx[i] = CalcInfo.theta_cur[i] - dx[i];
    dg[i] = CalcInfo.mo_grad[i] - dg[i];
  }

  if (Params.print_lvl > 3) {
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
    CalcInfo.mo_hess_diag[i] = 0.5 * (CalcInfo.mo_hess_diag[i] + tval);
    if (CalcInfo.mo_hess_diag[i] < MO_HESS_MIN) 
      CalcInfo.mo_hess_diag[i] = MO_HESS_MIN;
  }

  if (Params.print_lvl > 3) {
    outfile->Printf( "\nDS MO Hessian:\n");
    for (i=0; i<npairs; i++) 
      outfile->Printf( "%12.6lf\n", CalcInfo.mo_hess_diag[i]);
    outfile->Printf( "\n");
  }

  /* write thetas */
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) CalcInfo.theta_cur,
    npairs*sizeof(double));
  /* write gradient */  
  psio_write_entry(PSIF_DETCAS, "MO Gradient", (char *) CalcInfo.mo_grad,
    npairs*sizeof(double));
  /* write updated Hessian */
  psio_write_entry(PSIF_DETCAS, "Hessian", (char *) CalcInfo.mo_hess_diag,
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
  ncore = CalcInfo.num_fzc_orbs + CalcInfo.num_cor_orbs;

  if (Params.hessian == "DIAG") {
    CalcInfo.mo_hess_diag = init_array(npairs);
    
    if (Params.use_fzc_h == 1) 
      form_diag_mo_hess(npairs, ppair, qpair, CalcInfo.onel_ints, 
        CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.F_act, 
        ncore, CalcInfo.npop, CalcInfo.mo_hess_diag);
    else
      form_diag_mo_hess_yy(npairs, ppair, qpair, CalcInfo.onel_ints, 
        CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.lag,
        CalcInfo.mo_hess_diag);

    if (Params.print_lvl > 3)  
      IndPairs.print_vec(CalcInfo.mo_hess_diag,"\n\tDiagonal MO Hessian:");
  }
  else if (Params.hessian == "APPROX_DIAG") {
    CalcInfo.mo_hess_diag = init_array(npairs);
    form_appx_diag_mo_hess(npairs, ppair, qpair, CalcInfo.onel_ints, 
                      CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, 
                      CalcInfo.F_act, ncore, CalcInfo.npop, 
                      CalcInfo.mo_hess_diag);
    if (Params.print_lvl > 3)  
      IndPairs.print_vec(CalcInfo.mo_hess_diag,"\n\tAppx Diagonal MO Hessian:");
  }
  else if (Params.hessian == "FULL") {
    CalcInfo.mo_hess = block_matrix(npairs,npairs);
    form_full_mo_hess(npairs, ppair, qpair, CalcInfo.onel_ints, 
      CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.lag,
      CalcInfo.mo_hess);
    if (Params.print_lvl > 3) {
      outfile->Printf( "\nMO Hessian:\n");
      print_mat(CalcInfo.mo_hess,npairs,npairs,"outfile");
      outfile->Printf( "\n");
    }
  }
  else {
    outfile->Printf( "(detcas): Unrecognized Hessian option %s\n", 
      Params.hessian.c_str());
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

  CalcInfo.theta_step = init_array(npairs);

  // All this actually does is scale the gradient by the Hessian
  // If we have a diagonal (exact or approximate) Hessian

  // BFGS
  if (Params.scale_grad && Params.bfgs) {
    calc_orb_step_bfgs(npairs, CalcInfo.mo_grad, CalcInfo.mo_hess,
      CalcInfo.theta_step);
  }
  // non-BFGS diagonal Hessian
  else if (Params.scale_grad && ((Params.hessian == "DIAG") || 
       (Params.hessian == "APPROX_DIAG"))) {
    calc_orb_step(npairs, CalcInfo.mo_grad, CalcInfo.mo_hess_diag,
      CalcInfo.theta_step);
  }
  // non-BFGS full Hessian
  else if ((Params.scale_grad && (Params.hessian == "FULL")) ||
    Params.bfgs) {
    calc_orb_step_full(npairs, CalcInfo.mo_grad, CalcInfo.mo_hess,
      CalcInfo.theta_step);
  }
  // No Hessian available: take unit step down gradient direction
  else {
    for (pair=0; pair<npairs; pair++) 
      CalcInfo.theta_step[pair] = CalcInfo.mo_grad[pair];
  }
    

  if (Params.print_lvl > 3)  
    IndPairs.print_vec(CalcInfo.theta_step,"\n\tScaled Orbital Grad:");

  rms = 0.0;
  for (pair=0; pair<npairs; pair++) {
    value =  CalcInfo.theta_step[pair];
    rms += value * value;
  }

  rms = sqrt(rms);
  CalcInfo.scaled_mo_grad_rms = rms;
 
  if (Params.print_lvl)
    outfile->Printf( "\n\tScaled RMS Orbital Gradient: %6.4E\n", rms);

  if (Params.scale_step != 1.0) {
    for (pair=0; pair<npairs; pair++) 
      CalcInfo.theta_step[pair] *= Params.scale_step;
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
  if (Params.force_step) {
    CalcInfo.theta_cur[Params.force_pair] = Params.force_value;
    outfile->Printf( "Forcing step for pair %d of size %8.3E\n",
      Params.force_pair, Params.force_value);
    return(1);
  }

  for (pair=0; pair<npairs; pair++) 
    CalcInfo.theta_cur[pair] += CalcInfo.theta_step[pair];
    //CalcInfo.theta_cur[pair] = CalcInfo.theta_step[pair];

  if (CalcInfo.iter >= Params.diis_start) 
    took_diis = diis(npairs, CalcInfo.theta_cur, CalcInfo.theta_step);
    //took_diis = diis(npairs, CalcInfo.theta_cur, CalcInfo.mo_grad);
  else 
    took_diis = 0;

  if (!took_diis) {
    if (Params.print_lvl) 
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
  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    if (ir_npairs) {
      ir_norbs = CalcInfo.orbs_per_irr[h];
      ir_theta = IndPairs.get_irrep_vec(h, CalcInfo.theta_cur);
      ir_ppair  = IndPairs.get_ir_prel_ptr(h);
      ir_qpair  = IndPairs.get_ir_qrel_ptr(h);
      
      if (Params.print_lvl > 3) {
        outfile->Printf( "Thetas for irrep %d\n", h);
        for (pair=0; pair<ir_npairs; pair++) {
          outfile->Printf( "Pair (%2d,%2d) = %12.6lf\n",
                  ir_ppair[pair], ir_qpair[pair], ir_theta[pair]);
        }
        outfile->Printf( "\n");
        //fflush(outfile);
      }

      /* print old coefficients */
      if (Params.print_mos) {
        outfile->Printf( "\n\tOld molecular orbitals for irrep %s\n", 
          CalcInfo.labels[h]);
        print_mat(CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, "outfile");
      }

      if (Params.use_thetas)
        postmult_by_U(h, ir_norbs, CalcInfo.mo_coeffs[h], ir_npairs, 
          ir_ppair, ir_qpair, ir_theta);
      else
        postmult_by_exp_R(h, ir_norbs, CalcInfo.mo_coeffs[h], ir_npairs,
          ir_ppair, ir_qpair, ir_theta);

      /* print new coefficients */
      if (Params.print_mos) {
        outfile->Printf( "\n\tNew molecular orbitals for irrep %s\n", 
          CalcInfo.labels[h]);
        print_mat(CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, "outfile");
      }


      /* write the new block of MO coefficients to file30 */
      chkpt_init(PSIO_OPEN_OLD);
      chkpt_wt_scf_irrep(CalcInfo.mo_coeffs[h], h);
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
    CalcInfo.iter = 0;
    return(0);
  }

  if (fscanf(sumfile, "%d", &entries) != 1) {
    outfile->Printf("(print_step): Trouble reading num entries in file %s\n",
            sumfile_name);
    fclose(sumfile);
    CalcInfo.iter = 0;
    return(0);
  } 

  CalcInfo.iter = entries;
  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &nind, &scaled_rmsgrad,
           &rmsgrad, &energy_last, comment);
  }
  fclose(sumfile);

  chkpt_init(PSIO_OPEN_OLD);
  energy = chkpt_rd_etot();
  chkpt_close();

  /* check for convergence */
  if (rmsgrad < Params.rms_grad_convergence)
    converged_grad = 1;
  if (fabs(energy_last - energy) < Params.energy_convergence)
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

}} // end namespace psi::detcas

