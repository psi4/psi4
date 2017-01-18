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

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "sort.h"
#include "mrcc.h"
#include "idmrpt2.h"
#include "mp2_ccsd.h"
#include "transform.h"
#include "debugging.h"
#include "psimrcc.h"
#include "updater.h"

namespace psi{ namespace psimrcc{

using namespace std;

/*!
 * Runs a MRPT2 and a MRCCSD computation
 * @todo move this code in the CCMRCC class
 */
void mrccsd(SharedWavefunction ref_wfn, Options & options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  CCMRCC        mrcc(ref_wfn, options);

  if(options.get_bool("PERTURB_CBS") && options.get_bool("PERTURB_CBS_COUPLING")){
    mrcc.compute_first_order_amps();
  }

  options.print();
  // Initialize the appropriate updater
  Updater* updater;
//  if(options_get_str("CORR_ANSATZ")=="SR")
//    updater = static_cast<Updater*>(new MkUpdater());
  if(options.get_str("CORR_ANSATZ")=="MK")
    updater = dynamic_cast<Updater*>(new MkUpdater(options));
  if(options.get_str("CORR_ANSATZ")=="BW")
    updater = dynamic_cast<Updater*>(new BWUpdater(options));

	// Compute the energy
  mrcc.compute_energy(updater);

  if(options.get_bool("PERTURB_CBS")){
    mrcc.perturbative_cbs();
  }

  delete updater;
}

/*!
 * Runs a CCSD_MP2 computation
 */
void mp2_ccsd(SharedWavefunction ref_wfn, Options &options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  MP2_CCSD        mp2_ccsd(ref_wfn, options);

  // Compute the initial amplitudes and CCSD_MP2 energy
  mp2_ccsd.compute_mp2_ccsd_energy();

  DEBUGGING(1,
    blas->print_memory();
  )
}

/*!
 * Runs a MRPT2 computation
 */
void mrpt2(SharedWavefunction ref_wfn, Options &options)
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  IDMRPT2        idmrpt2(ref_wfn, options);

  Updater* updater = dynamic_cast<Updater*>(new MkUpdater(options));

  // Compute the initial amplitudes and MP2 energy
  idmrpt2.compute_mrpt2_energy(updater);

  delete updater;

  DEBUGGING(1,
    blas->print_memory();
  )
}

/*!
 * Runs a integral transformation
 * @todo CCTransform is still unused in the code
 */
void transform_integrals()
{
//   CCTransform transf;
//   transf.read_so_integrals();
}

}} /* End Namespaces */
