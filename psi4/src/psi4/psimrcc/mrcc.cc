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
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "algebra_interface.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"

namespace psi{ namespace psimrcc{

using namespace std;

CCMRCC::CCMRCC(SharedWavefunction ref_wfn, Options &options):
        CCManyBody(ref_wfn, options),
        options_(options)
{
  triples_type = ccsd;
  triples_coupling_type = cubic;
  ap_correction   = false; // Set tu true when computing the a posteriori correction
  current_energy  =  0.0;
  old_energy      = 10.0;

  // Parse the CORR_WFN parameter
  vector<string> theory_levels = split("PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT");
  for(size_t i=0;i<theory_levels.size();++i){
    if(options.get_str("CORR_WFN")==theory_levels[i])
      triples_type = TriplesType(i);
  }

  // Parse the COUPLING parameter
  vector<string> coupling_levels = split("NONE LINEAR QUADRATIC CUBIC");
  for(size_t i=0;i<coupling_levels.size();++i){
    if(options.get_str("COUPLING")==coupling_levels[i]){
      triples_coupling_type = TriplesCouplingType(i);
    }
  }

  // Parse the PERT_CBS parameter
  pert_cbs = options.get_bool("PERTURB_CBS");
  pert_cbs_coupling = options.get_bool("PERTURB_CBS_COUPLING");

  // Add the matrices that will store the intermediates
  add_matrices();

  // Generate the Fock matrices, Integrals and Denominators
  generate_integrals();
  generate_denominators();

  if(triples_type>ccsd)
    generate_triples_denominators();

  compute_reference_energy();

  DEBUGGING(1,
    blas->print_memory();
  )
}

CCMRCC::~CCMRCC()
{
}

}} /* End Namespaces */
