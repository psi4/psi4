/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include "sapt2p3.h"



namespace psi { namespace sapt {

void SAPT2p3::ind_disp30()
{
  double **tAR = block_matrix(aoccA_,nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uAR Amplitudes", (char *)
    tAR[0], sizeof(double)*aoccA_*nvirA_);

  double inddisp_1 = 2.0*C_DDOT(aoccA_*nvirA_,tAR[0],1,wBAR_[foccA_],1);

  free_block(tAR);

  double **tBS = block_matrix(aoccB_,nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uBS Amplitudes", (char *)
    tBS[0], sizeof(double)*aoccB_*nvirB_);

  double inddisp_2 = 2.0*C_DDOT(aoccB_*nvirB_,tBS[0],1,wABS_[foccB_],1);

  free_block(tBS);

  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,0,nvirA_);
  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,0,nvirB_);

  double **vARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"IndDisp30 uARBS Amplitudes",(char *)
    tARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,0.0,vARBS[0],aoccB_*nvirB_);

  double inddisp_3 = 4.0*C_DDOT((long int) aoccA_*nvirA_*aoccB_*nvirB_,
    vARBS[0],1,tARBS[0],1);

  free_block(B_p_AR);
  free_block(B_p_BS);
  free_block(vARBS);
  free_block(tARBS);

  e_ind_disp30_ = inddisp_1 + inddisp_2 + inddisp_3;

  if (debug_) {
    outfile->Printf("\n    Ind-Disp30_1        = %18.12lf [Eh]\n",inddisp_1);
    outfile->Printf("    Ind-Disp30_2        = %18.12lf [Eh]\n",inddisp_2);
    outfile->Printf("    Ind-Disp30_3        = %18.12lf [Eh]\n",inddisp_3);
  }
  if (print_) {
    outfile->Printf("    Ind-Disp30          = %18.12lf [Eh]\n",e_ind_disp30_);

  }
}

}}
