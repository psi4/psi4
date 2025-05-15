/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

/* get_moinfo(): Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified for ccresponse by TDC May, 2003
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn) {
    int i, j, h, p, q, errcod, nactive, nirreps, nfzc, nfzv;
    int *actpi, offset, act_offset;
    double **scf, ***C;
    psio_address next;

    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *)&(params.ref), sizeof(int));

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.nao = wfn->basisset()->nao();
    moinfo.labels = wfn->molecule()->irrep_labels();

    moinfo.sopi = wfn->nsopi();
    moinfo.orbspi = wfn->nmopi();
    moinfo.clsdpi = wfn->doccpi() - wfn->frzcpi();
    moinfo.openpi = wfn->soccpi();
    moinfo.frdocc = wfn->frzcpi();
    moinfo.fruocc = wfn->frzvpi();
    moinfo.uoccpi = moinfo.orbspi - moinfo.clsdpi - moinfo.openpi - moinfo.fruocc - moinfo.frdocc;

    moinfo.natom = wfn->molecule()->natom();

    nirreps = moinfo.nirreps;

    moinfo.ntri = moinfo.nmo * (moinfo.nmo + 1) / 2;
    moinfo.noei = moinfo.nso * (moinfo.nso + 1) / 2;
    moinfo.noei_ao = moinfo.nao * (moinfo.nao + 1) / 2;

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *)&(nactive), sizeof(int));
    moinfo.nactive = nactive;

    moinfo.nfzc = moinfo.frdocc.sum();

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
        moinfo.act_occpi = moinfo.occpi;

        moinfo.occ_sym = init_int_array(nactive);
        moinfo.vir_sym = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *)moinfo.occ_sym, sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *)moinfo.vir_sym, sizeof(int) * nactive);

        moinfo.occ_off = init_int_array(moinfo.nirreps);
        moinfo.vir_off = init_int_array(moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets", (char *)moinfo.occ_off, sizeof(int) * moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets", (char *)moinfo.vir_off, sizeof(int) * moinfo.nirreps);
    }

    moinfo.nvirt = moinfo.virtpi.sum();

    /*** arrange active SCF MO's ***/
    moinfo.actpi = moinfo.orbspi - moinfo.frdocc - moinfo.fruocc;
    moinfo.Ca = wfn->Ca_subset("SO", "ACTIVE");

    /* Get the active virtual orbitals */
    if (params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

        C = (double ***)malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for (h = 0; h < nirreps; h++) {
            if (moinfo.sopi[h] && moinfo.virtpi[h]) {
                C[h] = block_matrix(moinfo.sopi[h], moinfo.virtpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *)C[h][0],
                          sizeof(double) * moinfo.sopi[h] * moinfo.virtpi[h], next, &next);
            }
        }
        moinfo.C = C;
    }
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup() {
    int i;

    if (params.ref == 2) { /* UHF */
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
        free(moinfo.aocc_off);
        free(moinfo.bocc_off);
        free(moinfo.avir_off);
        free(moinfo.bvir_off);
    } else {
        for (i = 0; i < moinfo.nirreps; i++)
            if (moinfo.sopi[i] && moinfo.virtpi[i]) free_block(moinfo.C[i]);
        free(moinfo.C);
        free(moinfo.occ_sym);
        free(moinfo.vir_sym);
        free(moinfo.occ_off);
        free(moinfo.vir_off);
    }

    free(moinfo.mu_irreps);
    free(moinfo.l_irreps);
}

}  // namespace ccresponse
}  // namespace psi
