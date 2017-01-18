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

#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp21()
{
  double e_disp210 = disp21_1(PSIF_SAPT_AMPS,"gARAR x tARBS",
    "tARBS Amplitudes",aoccA_,nvirA_,aoccB_,nvirB_);
  e_disp210 += disp21_2(PSIF_SAPT_AMPS,"T AR Intermediates",
    "Theta AR Intermediates",aoccA_,nvirA_);

  if (debug_) {
    outfile->Printf("    Disp210             = %18.12lf [Eh]\n",e_disp210);
    
  }

  double e_disp201 = disp21_1(PSIF_SAPT_AMPS,"gBSBS x tARBS",
    "tARBS Amplitudes",aoccA_,nvirA_,aoccB_,nvirB_);
  e_disp201 += disp21_2(PSIF_SAPT_AMPS,"T BS Intermediates",
    "Theta BS Intermediates",aoccB_,nvirB_);

  if (debug_) {
    outfile->Printf("    Disp201             = %18.12lf [Eh]\n\n",e_disp201);
    
  }

  e_disp21_ = e_disp210 + e_disp201;

  if (print_) {
    outfile->Printf("    Disp21              = %18.12lf [Eh]\n",e_disp21_);
    
  }
}

double SAPT2p::disp21_1(int ampfile, const char *glabel, const char *tlabel, 
  int aoccA, int nvirA, int aoccB, int nvirB)
{
  double **tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
  psio_->read_entry(ampfile,tlabel,(char *) tARBS[0],
    sizeof(double)*aoccA*nvirA*aoccB*nvirB);

  double **gARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
  psio_->read_entry(ampfile,glabel,(char *) gARBS[0],
    sizeof(double)*aoccA*nvirA*aoccB*nvirB);

  double energy = 4.0*C_DDOT((long int) aoccA*nvirA*aoccB*nvirB,tARBS[0],1,
    gARBS[0],1);

  free_block(tARBS);
  free_block(gARBS);

  if (debug_) {
    outfile->Printf("\n    Disp21_1            = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

double SAPT2p::disp21_2(int ampfile, const char *tlabel, 
  const char *thetalabel, int aoccA, int nvirA)
{
  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,tlabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **theta_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) theta_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double energy = 8.0*C_DDOT((long int) aoccA*nvirA*(ndf_+3),T_p_AR[0],1,
    theta_p_AR[0],1);

  free_block(T_p_AR);
  free_block(theta_p_AR);

  if (debug_) {
    outfile->Printf("    Disp21_2            = %18.12lf [Eh]\n",energy);
    
  }

  return(energy);
}

}}
