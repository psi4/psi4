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

#include <cstdio>
#include "psi4/liboptions/liboptions.h"
#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::canonicalize_MO()
{
  if(reference == tcscf){
    bool canonicalize_active_favg   = options_.get_bool("CANONICALIZE_ACTIVE_FAVG");
    bool canonicalize_inactive_favg = options_.get_bool("CANONICALIZE_INACTIVE_FAVG");
    if(canonicalize_active_favg || canonicalize_inactive_favg){
      outfile->Printf("\n\n  Forming Favg for final canonicalization");
      construct_Favg();
      transform(Favg,Favg_t,C);

      Feff_t->zero();

      for(int h =0; h < nirreps; ++h){
        // Set the (closed,closed) blocks
        for(int i = 0; i < sopi[h]; ++i){
          Feff_t->set(h,i,i,Favg_t->get(h,i,i));
        }
      }

      if(canonicalize_inactive_favg){
        // Set the diagonal blocks Fock
        for(int h =0; h < nirreps; ++h){
          // Set the (closed,closed) blocks
          for(int i = 0; i < docc[h]; ++i){
            for(int j = 0; j < docc[h]; ++j){
              Feff_t->set(h,i,j,Favg_t->get(h,i,j));
            }
          }
          // Set the (virtual,virtual) blocks
          for(int i = docc[h] + actv[h]; i < sopi[h]; ++i){
            for(int j = docc[h] + actv[h]; j < sopi[h]; ++j){
              Feff_t->set(h,i,j,Favg_t->get(h,i,j));
            }
          }
        }
      }
      if(canonicalize_active_favg){
        // Set the diagonal blocks Fock
        for(int h =0; h < nirreps; ++h){
          // Set the (active,active) blocks
          for(int i = docc[h]; i< docc[h] + actv[h]; ++i){
            for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
              Feff_t->set(h,i,j,Favg_t->get(h,i,j));
            }
          }
        }
      }
      Feff_t.diagonalize(C_t,epsilon);
      T.multiply(false,false,C,C_t);
      C = T;

      // For debugging purposes
      // construct_Favg();
      // transform(Favg,Favg_t,C);
      // Favg_t->print();
    }
  }

  outfile->Printf("\n\n  Orbitals are canonicalized as:");
  if(options_.get_bool("FAVG") || options_.get_bool("CANONICALIZE_INACTIVE_FAVG"))
    outfile->Printf("\n  inactive (docc + uocc) : Fock(avg)");
  else
    outfile->Printf("\n  inactive (docc + uocc) : Fock(core)");

  if(options_.get_bool("CANONICALIZE_ACTIVE_FAVG"))
    outfile->Printf("\n  active   (actv)        : Fock(avg)");
  else
    outfile->Printf("\n  active   (actv)        : Fock(core)");

}

}} /* End Namespaces */
