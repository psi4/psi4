/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"

#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cctriples {

/*
 ** get_moinfo():  Routine to obtain basic orbital information from
 ** CC_INFO.
 **
 ** T. Daniel Crawford, October 1996.
 ** Modified by TDC, March 1999.
 */

void get_moinfo(std::shared_ptr<Wavefunction> wfn, Options &options) {
    int i, h, errcod, nactive, nirreps;
    std::string junk;
    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.labels = wfn->molecule()->irrep_labels();
    moinfo.enuc = wfn->molecule()->nuclear_repulsion_energy(wfn->get_dipole_field_strength());
    if (wfn->reference_wavefunction())
        moinfo.escf = wfn->reference_wavefunction()->energy();
    else
        moinfo.escf = wfn->energy();

    moinfo.orbspi = wfn->nmopi();
    moinfo.clsdpi = wfn->doccpi() - wfn->frzcpi();
    moinfo.openpi = wfn->soccpi();
    moinfo.frdocc = wfn->frzcpi();
    moinfo.fruocc = wfn->frzvpi();
    moinfo.uoccpi = moinfo.orbspi - moinfo.clsdpi - moinfo.openpi - moinfo.fruocc - moinfo.frdocc;

    nirreps = moinfo.nirreps;

    params.wfn = options.get_str("WFN");
    if (params.wfn != "CCSD" && params.wfn != "CCSD_T" && params.wfn != "CCSD_AT" && params.wfn != "BCCD" &&
        params.wfn != "BCCD_T") {
        throw PsiException("Invalid value of input keyword WFN", __FILE__, __LINE__);
    }

    params.nthreads = Process::environment.get_n_threads();
    if (options["CC_NUM_THREADS"].has_changed()) {
        params.nthreads = options.get_int("CC_NUM_THREADS");
    }

    params.semicanonical = 0;
    junk = options.get_str("REFERENCE");
    /* if no reference is given, assume rhf */
    if (junk == "RHF")
        params.ref = 0;
    else if (junk == "ROHF" && (params.wfn == "CCSD_T" || params.wfn == "BCCD_T")) {
        params.ref = 2;
        params.semicanonical = 1;
    } else if (junk == "ROHF")
        params.ref = 1;
    else if (junk == "UHF")
        params.ref = 2;
    else {
        throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
    }

    junk = options.get_str("DERTYPE");
    if (junk == "NONE")
        params.dertype = 0;
    else if (junk == "FIRST")
        params.dertype = 1;
    else {
        throw PsiException("Value of keyword DERTYPE is not applicable to CCSD(T)", __FILE__, __LINE__);
    }

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *)&(nactive), sizeof(int));

    if (params.ref == 2) { /** UHF **/

        moinfo.aoccpi = moinfo.clsdpi + wfn->soccpi();
        moinfo.boccpi = moinfo.clsdpi;
        moinfo.avirtpi = moinfo.uoccpi;
        moinfo.bvirtpi = moinfo.uoccpi + wfn->soccpi();

        moinfo.aocc_sym = init_int_array(nactive);
        moinfo.bocc_sym = init_int_array(nactive);
        moinfo.avir_sym = init_int_array(nactive);
        moinfo.bvir_sym = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry", (char *)moinfo.aocc_sym, sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry", (char *)moinfo.bocc_sym, sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry", (char *)moinfo.avir_sym, sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry", (char *)moinfo.bvir_sym, sizeof(int) * nactive);

        moinfo.aocc_off = init_int_array(moinfo.nirreps);
        moinfo.bocc_off = init_int_array(moinfo.nirreps);
        moinfo.avir_off = init_int_array(moinfo.nirreps);
        moinfo.bvir_off = init_int_array(moinfo.nirreps);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets", (char *)moinfo.aocc_off,
                        sizeof(int) * moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets", (char *)moinfo.bocc_off,
                        sizeof(int) * moinfo.nirreps);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets", (char *)moinfo.avir_off,
                        sizeof(int) * moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets", (char *)moinfo.bvir_off,
                        sizeof(int) * moinfo.nirreps);

    } else { /** RHF or ROHF **/

        moinfo.occpi = moinfo.clsdpi + wfn->soccpi();
        moinfo.virtpi = moinfo.uoccpi + wfn->soccpi();

        psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *)&(nactive), sizeof(int));

        moinfo.occ_sym = init_int_array(nactive);
        moinfo.vir_sym = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *)moinfo.occ_sym, sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *)moinfo.vir_sym, sizeof(int) * nactive);

        moinfo.occ_off = init_int_array(moinfo.nirreps);
        moinfo.vir_off = init_int_array(moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets", (char *)moinfo.occ_off, sizeof(int) * moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets", (char *)moinfo.vir_off, sizeof(int) * moinfo.nirreps);
    }

    outfile->Printf("\n\n");
    outfile->Printf("    Wave function   =    %6s\n", params.wfn.c_str());
    if (params.semicanonical) {
        outfile->Printf("    Reference wfn   =    ROHF changed to UHF for Semicanonical Orbitals\n");
    } else {
        outfile->Printf("    Reference wfn   =    %5s\n",
                        (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
    }
    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *)&(moinfo.eref), sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));

    outfile->Printf("\n    Nuclear Rep. energy (wfn)                = %20.15f\n", moinfo.enuc);
    outfile->Printf("    SCF energy          (wfn)                = %20.15f\n", moinfo.escf);
    outfile->Printf("    Reference energy    (file100)            = %20.15f\n", moinfo.eref);
    outfile->Printf("    CCSD energy         (file100)            = %20.15f\n", moinfo.ecc);
    outfile->Printf("    Total CCSD energy   (file100)            = %20.15f\n", moinfo.eref + moinfo.ecc);
}

/* Frees memory allocated in get_moinfo() and dumps some info. */
void cleanup() {
    if (params.ref == 2) {
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
    } else {
        free(moinfo.occ_sym);
        free(moinfo.vir_sym);
        free(moinfo.occ_off);
        free(moinfo.vir_off);
    }
}

}  // namespace cctriples
}  // namespace psi
