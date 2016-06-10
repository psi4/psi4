/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libmints/mints.h>
#include <libqt/qt.h>

#include <libsapt_solver/sapt0.h>
#include <libsapt_solver/sapt2.h>
#include <libsapt_solver/sapt2p.h>
#include <libsapt_solver/sapt2p3.h>
//#include <libsapt_solver/sapt_dft.h>
#include "wrapper.h"

namespace psi { namespace sapt {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType sapt(SharedWavefunction Dimer, SharedWavefunction MonomerA,
                   SharedWavefunction MonomerB, Options & options)
{
  tstart();

  boost::shared_ptr<PSIO> psio(new PSIO);

  if (options.get_str("SAPT_LEVEL") == "SAPT0") {
    SAPT0 sapt(Dimer, MonomerA, MonomerB, options, psio);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2") {
    SAPT2 sapt(Dimer, MonomerA, MonomerB, options, psio);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+") {
    SAPT2p sapt(Dimer, MonomerA, MonomerB, options, psio);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+3") {
    SAPT2p3 sapt(Dimer, MonomerA, MonomerB, options, psio);
    sapt.compute_energy();
  } else {
    throw PSIEXCEPTION("Unrecognized SAPT type");
  }
  // Shut down psi.

  tstop();

  return Success;
}

}}