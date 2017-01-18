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

//  MCSCF (2008) by Francesco Evangelista - frank@ccc.uga.edu
//  Center for Computational Chemistry, University of Georgia

/**
 *  @defgroup MCSCF MCSCF is a code for SCF/MCSCF computations
 *  @ingroup MCSCF
 *  @brief Contains main() and global variables
*/

// Standard libraries
#include <iostream>
#include <cstdlib>

// PSI C++ libraries
#include "psi4/libpsio/psio.hpp"

// PSI libraries
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmints/mintshelper.h"

#include "mcscf.h"
#include "scf.h"



namespace psi{
MOInfoSCF*     moinfo_scf = 0;

namespace mcscf{

MemoryManager* memory_manager = 0;

using namespace std;

/**
 * The main function
 * @param argc
 * @param argv[]
 * @return PSI_RETURN_SUCCESS if the program ran without any problem
 */
SharedWavefunction mcscf(SharedWavefunction ref_wfn, Options& options)
{
  using namespace psi;
  std::shared_ptr<PSIO> psio(new PSIO);

  tstart();

  memory_manager = new MemoryManager(Process::environment.get_memory());

  psio->open(PSIF_MCSCF, PSIO_OPEN_NEW);
  init_psi(options);

  SharedWavefunction wfn;
  if(options.get_str("REFERENCE") == "RHF"  ||
     options.get_str("REFERENCE") == "ROHF" ||
     options.get_str("REFERENCE") == "UHF"  ||
     options.get_str("REFERENCE") == "TWOCON")
  {
      // We need to start by generating some integrals
      MintsHelper* mints = new MintsHelper(ref_wfn->basisset(), options, 0);
      mints->integrals();
      delete mints;

      // Now, set the reference wavefunction for subsequent codes to use
      wfn = SharedWavefunction(new SCF(ref_wfn, options, psio));
      moinfo_scf      = new psi::MOInfoSCF(*(wfn.get()), options);
      wfn->compute_energy();

      Process::environment.globals["CURRENT ENERGY"] = wfn->reference_energy();
      Process::environment.globals["CURRENT REFERENCE ENERGY"] = wfn->reference_energy();
      Process::environment.globals["SCF TOTAL ENERGY"] = wfn->reference_energy();

  }else if(options.get_str("REFERENCE") == "MCSCF"){
      throw PSIEXCEPTION("REFERENCE = MCSCF not implemented yet");
  }

  if(moinfo_scf)     delete moinfo_scf;
  if(memory_manager) delete memory_manager;

  close_psi(options);
  psio->close(PSIF_MCSCF, 1);

  tstop();

  return wfn;
}

/**
 * Start Psi3, draw the program logo, version, and compilation details
 * @param argc
 * @param argv[]
 */
void init_psi(Options& options_)
{
  outfile->Printf("\n         ------------------------------------------");
  outfile->Printf("\n           MCSCF: a self-consistent field program");
  outfile->Printf("\n            written by Francesco A. Evangelista");
  outfile->Printf("\n         ------------------------------------------");
}

/**
 * Close psi by calling psio_done() and psi_stop()
 */
void close_psi(Options& options_)
{
  outfile->Printf("\n\n  MCSCF Execution Completed.\n\n");
}

}} /* End Namespaces */
