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

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

#include <iostream>

#include "psi4/psifiles.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

void SCF::read_so_oei()
{
  // Read all the SO one electron integrals in Pitzer order
  double* buffer = new double[nso*(nso+1)/2];

  // Read the kinetic energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;

  IWL::read_one(psio_.get(),PSIF_OEI, const_cast<char*>(PSIF_SO_T),buffer,nso*(nso+1)/2,0,0,"outfile");

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < H->get_rows(h); i++){
      for(int j=0; j < H->get_cols(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        H->set(h,i,j,buffer[ij]);
      }
    }
  }
  // Read the potential energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;

  IWL::read_one(psio_.get(),PSIF_OEI, const_cast<char*>(PSIF_SO_V),buffer,nso*(nso+1)/2,0,0,"outfile");
//  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_V),buffer,nso*(nso+1)/2,0,0,outfile);

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < H->get_rows(h); i++){
      for(int j=0; j < H->get_cols(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        H->add(h,i,j,buffer[ij]);
      }
    }
  }

  // Read the overlap integrals
  for(int k=0; k<nso*(nso+1)/2;++k) buffer[k] = 0.0;

  //  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_S),buffer,nso*(nso+1)/2,0,0,outfile);
  IWL::read_one(psio_.get(),PSIF_OEI, const_cast<char*>(PSIF_SO_S),buffer,nso*(nso+1)/2,0,0,"outfile");

  for(int h=0;h<nirreps;h++){
    for(int i=0; i < S->get_rows(h); i++){
      for(int j=0; j < S->get_rows(h); j++){
        int ij = INDEX(H->get_abs_row(h,i),H->get_abs_col(h,j));
        S->set(h,i,j, buffer[ij] );
      }
    }
  }

  delete[] buffer;

  if(options_.get_int("DEBUG") > 4){
    S->print();
    H->print();
  }
}

}} /* End Namespaces */
