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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/**
 *  @file ccmrcc_compute.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains all the methods to compute the energy
*/

#include <cstdlib>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "mrcc.h"
#include "debugging.h"
#include "updater.h"

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;


/**
 * This is a generic coupled cluster cycle. By specifying the updater object you can get all the flavors of CC, single-reference, Mukherjee MRCC,...
 * @param updater the pointer to a Updater object
 */
void CCMRCC::compute_energy(Updater* updater)
{
  blas->diis_add("t1[o][v]{u}","t1_delta[o][v]{u}");
  blas->diis_add("t1[O][V]{u}","t1_delta[O][V]{u}");
  blas->diis_add("t2[oo][vv]{u}","t2_delta[oo][vv]{u}");
  blas->diis_add("t2[oO][vV]{u}","t2_delta[oO][vV]{u}");
  blas->diis_add("t2[OO][VV]{u}","t2_delta[OO][VV]{u}");

  Timer cc_timer;
  bool converged = false;
  // Start CC cycle
  int cycle = 0;
  while(!converged){

    updater->zero_internal_amps();

    synchronize_amps();
    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();
    build_t1_amplitudes();
    build_t2_amplitudes();
    blas->compute();
    if(triples_type>ccsd_t)
      build_t1_amplitudes_triples();
    if(triples_type>ccsd_t)
      build_t2_amplitudes_triples();

    converged = build_diagonalize_Heff(cycle,cc_timer.get());

    h_eff.set_eigenvalue(current_energy);
    h_eff.set_matrix(Heff,moinfo->get_nrefs());
    h_eff.set_right_eigenvector(right_eigenvector,moinfo->get_nrefs());
    h_eff.set_left_eigenvector(left_eigenvector,moinfo->get_nrefs());

    if(!converged){
      blas->diis_save_t_amps(cycle);
      updater->update(cycle,&h_eff);
      updater->zero_internal_delta_amps();
      compute_delta_amps();
      blas->diis(cycle,delta_energy,DiisCC);
    }

    if(cycle>options_.get_int("MAXITER")){
      outfile->Printf("\n\n\tThe calculation did not converge in %d cycles\n\tQuitting PSIMRCC\n",options_.get_int("MAXITER"));

      exit(1);
    }
    cycle++;
  }

  outfile->Printf("\n\n  Timing for singles and doubles: %20.6f s",cc_timer.get());

  if(options_.get_str("CORR_WFN")=="CCSD_T"){
    compute_perturbative_triples();
  }

  // Compute the apBWCCSD energy
  if(ap_correction){
    updater->zero_internal_amps();

    synchronize_amps();

    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_t1_amplitudes();
    build_t2_amplitudes();

    updater->update(cycle,&h_eff);

    updater->zero_internal_amps();

    synchronize_amps();

    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_t1_amplitudes();
    build_t2_amplitudes();

    updater->zero_internal_amps();

    converged = build_diagonalize_Heff(-1,cc_timer.get());
  }

  DEBUGGING(1,
    blas->print_memory();
  );

  CCOperation::print_timing();
}

}} /* End Namespaces */
