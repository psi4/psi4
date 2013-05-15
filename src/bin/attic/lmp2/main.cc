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

#include <boost/python.hpp>

#include <cstdlib>
#include <cstdio>
#include <string>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <psi4-dec.h>

#include <liblmp2_solver/lmp2.h>

using namespace boost;
using namespace std;

namespace psi { namespace lmp2 {

PsiReturnType lmp2(Options & options)
{
    tstart();

#ifdef HAVE_MADNESS

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    // Initialize the psi3 timer library.
    timer_init();

    string reference = options.get_str("REFERENCE");
    boost::shared_ptr<Wavefunction> lmp2;
    double energy;

    if (reference == "RHF") {
        lmp2 = boost::shared_ptr<Wavefunction>(new LMP2(options,
                                                        Process::environment.reference_wavefunction()));
    }
    else {
        throw InputException("Unknown reference " + reference, "REFERENCE", __FILE__, __LINE__);
        energy = 0.0;
    }

    // Set this early because the callback mechanism uses it.
    Process::environment.set_reference_wavefunction(lmp2);

    energy = lmp2->compute_energy();

    WorldComm->sync();

    // Set some environment variables
    Process::environment.globals["LMP2 ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;

    // Shut down psi.
    timer_done();
#else
    throw PSIEXCEPTION("LMP2 currently only works with the MADNESS parallel runtime.\n");
#endif

    tstop();

    return Success;
}

}}
