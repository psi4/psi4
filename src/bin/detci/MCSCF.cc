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

#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <psifiles.h>
#define EXTERN
#include "MCSCF.h"


namespace psi { namespace detci {

extern void mcscf_read_integrals(void);
extern void mcscf_get_mat_block(double **src, double **dst, int dst_dim,
                          int dst_offset, int *dst2src);

#define MO_HESS_MIN 1.0E-1

}} // end namespace psi::detci


namespace psi { namespace detci {

MCSCF::MCSCF(Options& options, OutFile& IterSummaryOut)
    : options_(options), IterSummaryOut_(IterSummaryOut)
{
    title();
    get_mo_info(options_);

    // Form independent pairs
    IndPairs.set(CalcInfo.nirreps, MAX_RAS_SPACES, CalcInfo.ras_opi,
                 MCSCF_CalcInfo.ras_orbs, MCSCF_CalcInfo.frozen_docc, MCSCF_CalcInfo.fzc_orbs,
                 MCSCF_CalcInfo.rstr_docc, MCSCF_CalcInfo.cor_orbs,
                 MCSCF_CalcInfo.rstr_uocc, MCSCF_CalcInfo.vir_orbs,
                 MCSCF_CalcInfo.frozen_uocc, MCSCF_CalcInfo.fzv_orbs,
                 MCSCF_CalcInfo.ci2relpitz, MCSCF_Parameters.ignore_ras_ras, MCSCF_Parameters.ignore_fz);

    if (MCSCF_Parameters.print_lvl > 3) IndPairs.print();
    num_indep_pairs_ = IndPairs.get_num_pairs();

    // Setup the DIIS manager
    diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(MCSCF_Parameters.diis_max_vecs,
                        "MCSCF DIIS", DIISManager::OldestAdded, DIISManager::InCore));
    diis_iter_ = 0; ndiis_vec_ = 0;
    diis_manager_->set_error_vector_size(1, DIISEntry::Pointer, num_indep_pairs_);
    diis_manager_->set_vector_size(1, DIISEntry::Pointer, num_indep_pairs_);

    // Setup general parameters
    theta_cur_ = init_array(num_indep_pairs_);
    zero_arr(theta_cur_, num_indep_pairs_);
}

MCSCF::~MCSCF()
{
}


int MCSCF::update(void)
{
  int converged = 0;
  int steptype = 0;

  mcscf_read_integrals();            /* get the 1 and 2 elec MO integrals        */
  read_density_matrices(options_);

  // Compute lagrangian
  MCSCF_CalcInfo.lag = lagcalc(MCSCF_CalcInfo.opdm, MCSCF_CalcInfo.tpdm, MCSCF_CalcInfo.onel_ints_bare,
                     MCSCF_CalcInfo.twoel_ints, CalcInfo.nmo,
                     MCSCF_CalcInfo.npop, MCSCF_Parameters.print_lvl, PSIF_MO_LAG); 

  if (MCSCF_Parameters.print_lvl > 2)
    IndPairs.print_vec(theta_cur_, "\n\tRotation Angles:");

  if (!read_ref_orbs()) {
    read_cur_orbs();
    write_ref_orbs();
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
  }
  else
    steptype = 0;

  print_step(MCSCF_CalcInfo.iter, num_indep_pairs_, steptype, IterSummaryOut_);

  // Cleanup the iteration
  iteration_clean();

  if (MCSCF_Parameters.print_lvl) tstop();

  return converged;
}


/*
** title(): Function prints a program identification
*/
void MCSCF::title(void)
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
}


/*!
** calc_gradient()
**
** This function calculates the MO gradient from the MO Lagrangian
**
** \ingroup DETCAS
*/
void MCSCF::calc_gradient(void)
{
  int pair, npair, h, ir_npairs, ir_norbs, offset;
  double *ir_mo_grad, **ir_lag, *ir_theta_cur, value, rms;
  int *parr, *qarr, *ir_ppair, *ir_qpair;
  int p, q;

  npair = IndPairs.get_num_pairs();
  parr  = IndPairs.get_p_ptr();
  qarr  = IndPairs.get_q_ptr();

  MCSCF_CalcInfo.mo_grad = init_array(npair);

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
    mcscf_get_mat_block(MCSCF_CalcInfo.lag, ir_lag, ir_norbs, offset, CalcInfo.reorder);

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "Irrep %d of lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    ir_theta_cur = IndPairs.get_irrep_vec(h, theta_cur_); 

    // Need to mult the Lagrangian by 2 to get dEdU
    C_DSCAL(ir_norbs*ir_norbs, 2.0, ir_lag[0], 1);

    if (MCSCF_Parameters.print_lvl > 3) {
      outfile->Printf( "Irrep %d of 2 * lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, "outfile");
    }

    if (MCSCF_Parameters.use_thetas) {
      // Calc dEdU
      premult_by_U(h, CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
                   ir_ppair, ir_qpair, ir_theta_cur);

      if (MCSCF_Parameters.print_lvl > 3) {
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
void MCSCF::bfgs_hessian(void)
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
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) theta_cur_,
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
    dx[i] = theta_cur_[i] - dx[i];
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
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) theta_cur_,
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
void MCSCF::ds_hessian(void)
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
    psio_write_entry(PSIF_DETCAS, "Thetas", (char *) theta_cur_,
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
    dx[i] = theta_cur_[i] - dx[i];
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
  psio_write_entry(PSIF_DETCAS, "Thetas", (char *) theta_cur_,
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
void MCSCF::calc_hessian(void)
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
void MCSCF::scale_gradient(void)
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
    outfile->Printf( "\n\tScaled RMS Orbital Gradient: %6.10E\n", rms);

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
int MCSCF::take_step(void)
{
  int took_diis;
  
  /* for debugging purposes */
  if (MCSCF_Parameters.force_step) {
    theta_cur_[MCSCF_Parameters.force_pair] = MCSCF_Parameters.force_value;
    outfile->Printf( "Forcing step for pair %d of size %8.3E\n",
      MCSCF_Parameters.force_pair, MCSCF_Parameters.force_value);
    return(1);
  }

  // Add new step to theta (Newton-Raphson)
  for (int pair=0; pair<num_indep_pairs_; pair++) 
    theta_cur_[pair] += MCSCF_CalcInfo.theta_step[pair];

  if (MCSCF_CalcInfo.iter >= MCSCF_Parameters.diis_start){ 
    // Add DIIS vector
    diis_manager_->add_entry(2, MCSCF_CalcInfo.theta_step, theta_cur_);
    ndiis_vec_++;

    // Check if we should skip DIIS
    if ((diis_iter_ % MCSCF_Parameters.diis_freq) || (ndiis_vec_ < MCSCF_Parameters.diis_min_vecs)){
        took_diis = 0;
    }
    else {
        // Extrapolate new thetas
        outfile->Printf("Taking a DIIS step\n");
        zero_arr(theta_cur_, num_indep_pairs_);
        diis_manager_->extrapolate(1, theta_cur_);
        diis_iter_++;
        took_diis = 1;
    }

  }
  else { 
    took_diis = 0;
  }
  if (!took_diis) {
    if (MCSCF_Parameters.print_lvl) 
      outfile->Printf("Taking regular step\n");
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
void MCSCF::rotate_orbs(void)
{
  double *ir_theta;
  int h, pair, ir_norbs, ir_npairs, *ir_ppair, *ir_qpair;

  // First, we need to come up with Theta vectors for each irrep
  ir_theta = init_array(IndPairs.get_num_pairs());  // always big enough
  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    if (ir_npairs) {
      ir_norbs = CalcInfo.orbs_per_irr[h];
      ir_theta = IndPairs.get_irrep_vec(h, theta_cur_);
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
          CalcInfo.labels[h]);
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
          CalcInfo.labels[h]);
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
int MCSCF::check_conv(void)
{
  int converged_energy=0, converged_grad=0, last_converged=0;


  /* check for convergence */
  if (MCSCF_CalcInfo.mo_grad_rms < MCSCF_Parameters.rms_grad_convergence)
    converged_grad = 1;
  if (fabs(MCSCF_CalcInfo.energy_old - MCSCF_CalcInfo.energy) < MCSCF_Parameters.energy_convergence)
    converged_energy = 1;

  if (converged_grad && converged_energy) {
    outfile->Printf( "\n\t*** Calculation Converged ***\n");
    return(1);
  }
  else {
    outfile->Printf( "\n\t... calculation continuing ...\n");
    return(0);
  }
}
/*
** finalize()
**
** This function frees any allocated global variables
**
*/
void MCSCF::finalize(void)
{
  int i;
  free(MCSCF_CalcInfo.frozen_docc);
  free(MCSCF_CalcInfo.frozen_uocc);
  free(MCSCF_CalcInfo.rstr_docc);
  free(MCSCF_CalcInfo.rstr_uocc);
  free(MCSCF_CalcInfo.ci2relpitz);
  free_int_matrix(MCSCF_CalcInfo.fzc_orbs);
  free_int_matrix(MCSCF_CalcInfo.fzv_orbs);
  for (i=0; i<MAX_RAS_SPACES; i++)
    free_int_matrix(MCSCF_CalcInfo.ras_orbs[i]);
  free(MCSCF_CalcInfo.ras_orbs);

  for (i=0; i<CalcInfo.nirreps; i++) {
    if (CalcInfo.orbs_per_irr[i])
      free_block(MCSCF_CalcInfo.mo_coeffs[i]);
  }
  free(MCSCF_CalcInfo.mo_coeffs);

  // DGAS updated
  diis_manager_->delete_diis_file();
  free(theta_cur_);
}

/*
** iteration_clean()
**
** Clean up intermediate quantities
**
*/
void MCSCF::iteration_clean(void)
{
  free(MCSCF_CalcInfo.onel_ints);
  free(MCSCF_CalcInfo.onel_ints_bare);
  free(MCSCF_CalcInfo.twoel_ints);
  free_block(MCSCF_CalcInfo.opdm);
  free(MCSCF_CalcInfo.tpdm);
  free_block(MCSCF_CalcInfo.lag);
  free(MCSCF_CalcInfo.F_act);
  free(MCSCF_CalcInfo.mo_grad);
  if (MCSCF_CalcInfo.mo_hess_diag != NULL) free(MCSCF_CalcInfo.mo_hess_diag);
  if (MCSCF_CalcInfo.mo_hess != NULL) free_block(MCSCF_CalcInfo.mo_hess);
  free(MCSCF_CalcInfo.theta_step);
}


}} // end namespace psi::detci

