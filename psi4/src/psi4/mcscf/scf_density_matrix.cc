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

#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::density_matrix()
{
  // Closed-shell density
  O->zero();
  for(int h = 0; h < nirreps; ++h){
    for(int i = 0; i < docc[h]; ++i){
      O->set(h, i, i, 1.0);
    }
  }
  transform(O,Dc,C_T);

  // ROHF density
  if(reference == rohf){
    O->zero();
    for(int h = 0; h < nirreps; ++h){
      for(int i = docc[h]; i < docc[h] + actv[h]; ++i){
        O->set(h, i, i, 1.0);
      }
    }
    transform(O,Do,C_T);
  }

  // TWOCON density
  if(reference == tcscf){
    for(int I = 0 ; I < nci; ++I){
      O->zero();
      O->set(tcscf_sym[I], tcscf_mos[I], tcscf_mos[I], 1.0);
      transform(O,Dtc[I],C_T);
    }
  }
}

}} /* End Namespaces */
