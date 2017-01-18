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

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "mrcc.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCMRCC::print_mrccsd_energy(int cycle)
{
  delta_energy = current_energy-old_energy;
  if(cycle==0){
    print_method("\tMultireference Coupled Cluster\n\t\tUsing the DPD Library");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    outfile->Printf("\n  @CC Cycle      Energy          Delta E    ||DeltaT1|| ||DeltaT2|| Timing  DIIS");
    outfile->Printf("\n  @CC             [Eh]            [Eh]                               (Sec)");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
  }
  if(cycle>=0){
    outfile->Printf("\n  @CC %3d  %18.12f  %11.4e   %8.3e   %8.3e %7.0f",cycle,current_energy,delta_energy,delta_t1_amps,delta_t2_amps,total_time);

    bool is_converged = (delta_t1_amps < options_.get_double("R_CONVERGENCE") &&
                         delta_t2_amps < options_.get_double("R_CONVERGENCE") &&
                         fabs(delta_energy) < options_.get_double("E_CONVERGENCE"));

    if(is_converged && (cycle!=0)){
      char star = (options_.get_str("CORR_WFN") == "CCSD") ? '*' : ' ';
      outfile->Printf("\n  ------------------------------------------------------------------------------");
      outfile->Printf("\n\n%6c%1c Mk-MRCCSD total energy      = %20.12f\n",' ',star,current_energy);
    }
  }else if(cycle==-1){
    char star = ' ';
    if(options_.get_str("CORR_WFN") == "CCSD")
      star = '*';
    outfile->Printf("\n\n%6c%1c Mk-MRCCSD total energy      = %20.12f\n",' ',star,current_energy);
    print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
  }


}

}} /* End Namespaces */
