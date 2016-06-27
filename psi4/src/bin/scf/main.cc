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

#include <boost/python.hpp>

#include <cstdlib>
#include <cstdio>
#include <string>

#include "psi4/include/psifiles.h"
#include "psi4/src/lib/libciomr/libciomr.h"
#include "psi4/src/lib/libpsio/psio.h"
#include "psi4/src/lib/libpsio/psio.hpp"
#include "psi4/src/lib/libparallel/parallel.h"
#include "psi4/src/lib/libiwl/iwl.h"
#include "psi4/src/lib/libqt/qt.h"


#include <libmints/writer.h>
#include <libmints/writer_file_prefix.h>
#include "psi4/include/psi4-dec.h"

#include "psi4/src/lib/libscf_solver/rhf.h"
#include <libscf_solver/rohf.h>
#include <libscf_solver/uhf.h>
#include <libscf_solver/cuhf.h>
#include <libscf_solver/ks.h>


using namespace boost;
using namespace std;

namespace psi { namespace scf {

SharedWavefunction scf(SharedWavefunction ref_wfn, Options & options, PyObject* pre, PyObject* post)
{
    tstart();

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    string reference = options.get_str("REFERENCE");
    boost::shared_ptr<Wavefunction> scf;
    double energy;
    //bool parallel = options.get_bool("PARALLEL");


    if (reference == "RHF") {
        scf = boost::shared_ptr<Wavefunction>(new RHF(ref_wfn, options, psio));
    }
    else if (reference == "ROHF") {
        scf = boost::shared_ptr<Wavefunction>(new ROHF(ref_wfn, options, psio));
    }
    else if (reference == "UHF") {
        scf = boost::shared_ptr<Wavefunction>(new UHF(ref_wfn, options, psio));
    }
    else if (reference == "CUHF") {
        scf = boost::shared_ptr<Wavefunction>(new CUHF(ref_wfn, options, psio));
    }
    else if (reference == "RKS") {
        scf = boost::shared_ptr<Wavefunction>(new RKS(ref_wfn, options, psio));
    }
    else if (reference == "UKS") {
        scf = boost::shared_ptr<Wavefunction>(new UKS(ref_wfn, options, psio));
    }
    else {
        throw InputException("Unknown reference " + reference, "REFERENCE", __FILE__, __LINE__);
        energy = 0.0;
    }

    // print the basis set
    if ( options.get_bool("PRINT_BASIS") ) {
        boost::shared_ptr<BasisSet> basisset = BasisSet::pyconstruct_orbital(scf->molecule(),
            "BASIS", options.get_str("BASIS"));
        basisset->print_detail();
    }

    // Set this early because the efp mechanism uses it.
    Process::environment.set_legacy_wavefunction(scf);

    if (pre)
        scf->add_preiteration_callback(pre);
    if (post)
        scf->add_postiteration_callback(post);

    energy = scf->compute_energy();

    // Print a molden file
    if ( options.get_bool("MOLDEN_WRITE") ) {
       boost::shared_ptr<MoldenWriter> molden(new MoldenWriter(scf));
       std::string filename = get_writer_file_prefix(scf->molecule()->name()) + ".molden";
       HF* hf = (HF*)scf.get();
       SharedVector occA = hf->occupation_a();
       SharedVector occB = hf->occupation_b();
       molden->write(filename, scf->Ca(), scf->Cb(), scf->epsilon_a(), scf->epsilon_b(),occA,occB);
    }

    // Print molecular orbitals
    if ( options.get_bool("PRINT_MOS") ) {
       boost::shared_ptr<MOWriter> mo(new MOWriter(scf,options));
       mo->write();
    }

    // Set some environment variables
    Process::environment.globals["SCF TOTAL ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = energy;

    // Shut down psi.
    tstop();

    return scf;
}

}}