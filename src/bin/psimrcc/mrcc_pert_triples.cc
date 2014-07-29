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

/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <liboptions/liboptions.h>
#include <boost/shared_ptr.hpp>
#include <libchkpt/chkpt.hpp>
#include <psi4-dec.h>

#include "mrcc.h"
#include "mrccsd_t.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

void CCMRCC::compute_perturbative_triples()
{
  Timer timer;

  h_eff.set_eigenvalue(current_energy);
  h_eff.set_matrix(Heff,moinfo->get_nrefs());
  h_eff.set_right_eigenvector(right_eigenvector,moinfo->get_nrefs());
  h_eff.set_left_eigenvector(left_eigenvector,moinfo->get_nrefs());
  h_eff.set_zeroth_order_eigenvector(zeroth_order_eigenvector,moinfo->get_nrefs());

  MRCCSD_T mrccsd_t(options_,&h_eff);

  if(options_.get_bool("DIAGONALIZE_HEFF")){
    psi::fprintf(outfile,"\n\n  Diagonalizing Heff");
    current_energy = h_eff.diagonalize();
  }else{
    psi::fprintf(outfile,"\n\n  Computing the expectation value of Heff");
    current_energy = h_eff.expectation_value();
  }
  Process::environment.globals["CURRENT ENERGY"]    = current_energy;
  Process::environment.globals["MRCC TOTAL ENERGY"] = current_energy;

  psi::fprintf(outfile,"\n\n%6c* Mk-MRCCSD(T) total energy   =    %20.12f",' ',current_energy);
  psi::fprintf(outfile,"\n\n  Timing for triples:             %20.6f s",timer.get());
  fflush(outfile);
}


}}  /* End Namespaces */
