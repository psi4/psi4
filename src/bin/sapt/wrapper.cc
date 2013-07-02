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

PsiReturnType sapt(Options & options)
{
  tstart();

  boost::shared_ptr<PSIO> psio(new PSIO);
  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  if (options.get_str("SAPT_LEVEL") == "SAPT0") {
    SAPT0 sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2") {
    SAPT2 sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+") {
    SAPT2p sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+3") {
    SAPT2p3 sapt(options, psio, chkpt);
    sapt.compute_energy();
//  } else if (options.get_str("SAPT_LEVEL") == "MP2C") {
//    MP2C sapt(options, psio, chkpt);
//    sapt.compute_energy();
  } else {
    throw PSIEXCEPTION("Unrecognized SAPT type");
  }

  // Shut down psi.

  tstop();

  return Success;
}

}}

