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

#include "psi4/liboptions/liboptions.h"
#include <cstdio>

#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::construct_Feff(int cycle)
{
  Feff_t_old = Feff_t;
  Feff_t = Fc_t;

  if(options_.get_bool("FAVG")){
    if(cycle >= options_.get_int("FAVG_START")){
      // Set the diagonal blocks Fock
      for(int h =0; h < nirreps; ++h){
        // Set the (closed,closed) blocks
        for(int i = 0; i < docc[h]; ++i){
          for(int j = 0; j < docc[h]; ++j){
            Feff_t->set(h,i,j,Favg_t->get(h,i,j));
          }
        }
        // Set the (active,active) blocks
        for(int i = docc[h]; i< docc[h] + actv[h]; ++i){
          for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
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
    if(cycle == options_.get_int("FAVG_START")){
      outfile->Printf("\n  *** Switching from Fc to F_avg ***");
    }
  }

  if(reference == rohf  && (cycle > turn_on_actv)){
    for(int h =0; h < nirreps; ++h){
      // Set the (closed,open) and (open,closed) blocks to 2(Fc-Fo)
      for(int i =0; i < docc[h]; ++i){
        for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
          double element = 2.0 * ( Fc_t->get(h,i,j) - Fo_t->get(h,i,j) );
          Feff_t->set(h,i,j,element);
          Feff_t->set(h,j,i,element);
        }
      }

      // Set the (virtual,open) and (open,virtual) blocks to 2 Fo
      for(int i = docc[h] + actv[h]; i < sopi[h]; ++i){
        for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
          double element = 2.0 * Fo_t->get(h,i,j);
          Feff_t->set(h,i,j,element);
          Feff_t->set(h,j,i,element);
        }
      }
    }
  }
  if(reference == tcscf && (cycle > turn_on_actv)){
    for(int I = 0 ; I < nci; ++I){
      int h = tcscf_sym[I];
      int i = tcscf_mos[I];
      // Set the (closed,tc) and (tc,closed) blocks to 2 (Fc - Ftc)
      for(int j = 0; j < docc[h]; ++j){
        double element = 2.0 * ( Fc_t->get(h,i,j) - Ftc_t[I]->get(h,i,j) );
        Feff_t->set(h,i,j,element);
        Feff_t->set(h,j,i,element);
      }
      // Set the (external,tc) and (tc,external) blocks to 2 Ftc
      for(int j = docc[h] + actv[h]; j < sopi[h]; ++j){
        double element = 2.0 * Ftc_t[I]->get(h,i,j);
        Feff_t->set(h,i,j,element);
        Feff_t->set(h,j,i,element);
      }
    }
    // Orbitals of the same symmetry
    if(options_.get_bool("INTERNAL_ROTATIONS")){
      for(int I = 0 ; I < nci; ++I){
        for(int J = I + 1 ; J < nci; ++J){
    	  if(tcscf_sym[I] == tcscf_sym[J]){
    		int h = tcscf_sym[I];
    		int i = tcscf_mos[I];
    		int j = tcscf_mos[J];
    		// Set the (tc,tc) and (tc,tc) blocks to 2 (Ftc_a - Ftc_b)
    		double element = 2.0 * (Ftc_t[I]->get(h,i,j) - Ftc_t[J]->get(h,i,j));
    		Feff_t->set(h,i,j,element);
    		Feff_t->set(h,j,i,element);
    	  }
    	}
      }
    }else{
      for(int I = 0 ; I < nci; ++I){
        for(int J = I + 1 ; J < nci; ++J){
          if(tcscf_sym[I] == tcscf_sym[J]){
            int h = tcscf_sym[I];
            int i = tcscf_mos[I];
            int j = tcscf_mos[J];
        	// Set the (tc,tc) and (tc,tc) blocks to 2 (Ftc_a - Ftc_b)
            double element = 0.0;
            Feff_t->set(h,i,j,element);
            Feff_t->set(h,j,i,element);
          }
        }
      }
    }
  }
  // Level shift
  double shift = options_.get_double("LEVEL_SHIFT");
  outfile->Printf("\n  Setting level shift to %.3f",shift);
  for(int h =0; h < nirreps; ++h){
    for(int i = docc[h] + actv[h]; i < sopi[h]; ++i){
      Feff_t->add(h,i,i,shift);
    }
  }
//
//  double dumping = static_cast<double>(options_.get_int("DUMPING")) / 100.0;
//  outfile->Printf("\n  Setting dumping to %.3f",dumping);
//  // Dumping
//  for(int h =0; h < nirreps; ++h){
//    // Set the (virtual,open) and (open,virtual) blocks to 2 Fo
//    for(int i = 0; i < sopi[h]; ++i){
//      for(int j = 0; j < sopi[h]; ++j){
//        double element = (1.0 - dumping) * Feff_t->get(h,i,j) + dumping * Feff_t_old->get(h,i,j);
//        Feff_t->set(h,i,j,element);
//      }
//    }
//  }
}

}} /* End Namespaces */
