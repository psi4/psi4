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
#include <libmints/writer.h>
#include <libmints/writer_file_prefix.h>
#include <psi4-dec.h>

#include <libscf_solver/rhf.h>
#include <libscf_solver/rohf.h>
#include <libscf_solver/uhf.h>
#include <libscf_solver/cuhf.h>
#include <libscf_solver/ks.h>


using namespace boost;
using namespace std;

namespace psi { namespace scf {

PsiReturnType scf_dummy(Options & options)
{
    tstart();
    double energy;

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    string reference = options.get_str("REFERENCE");
    boost::shared_ptr<Wavefunction> scf;

    if (reference == "RHF") {
        scf = boost::shared_ptr<Wavefunction>(new RHF(options, psio));
    }
    else if (reference == "ROHF") {
        scf = boost::shared_ptr<Wavefunction>(new ROHF(options, psio));
    }
    else if (reference == "UHF") {
        scf = boost::shared_ptr<Wavefunction>(new UHF(options, psio));
    }
    else if (reference == "CUHF") {
        scf = boost::shared_ptr<Wavefunction>(new CUHF(options, psio));
    }
    else if (reference == "RKS") {
        scf = boost::shared_ptr<Wavefunction>(new RKS(options, psio));
    }
    else if (reference == "UKS") {
        scf = boost::shared_ptr<Wavefunction>(new UKS(options, psio));
    }
    else {
        throw InputException("Unknown reference " + reference, "REFERENCE", __FILE__, __LINE__);
        energy = 0.0;
    }

    outfile->Printf("\tWarning: Constructing empty scf wavefunction object.\n");
    outfile->Printf("\tLoading _some_ variables from chkpt (not all coded yet).\n");

    scf->load_values_from_chkpt();

    Process::environment.set_wavefunction(scf);

    if(psio->open_check(PSIF_CHKPT))
        psio->close(PSIF_CHKPT, 1);

    tstop();

    return Success;
}

}}
