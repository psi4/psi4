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

namespace psi { namespace ccresponse {

/* get_moinfo(): Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified for ccresponse by TDC May, 2003
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn)
{
    int i, j, h, p, q, errcod, nactive, nirreps, nfzc, nfzv;
    int *actpi, offset, act_offset;
    double **scf, ***C;
    psio_address next;

    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(params.ref),
                    sizeof(int));

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.nao = wfn->basisset()->nao();
    moinfo.labels = wfn->molecule()->irrep_labels();

    moinfo.sopi = init_int_array(moinfo.nirreps);
    moinfo.orbspi = init_int_array(moinfo.nirreps);
    moinfo.clsdpi = init_int_array(moinfo.nirreps);
    moinfo.openpi = init_int_array(moinfo.nirreps);
    for(int h = 0; h < moinfo.nirreps; ++h){
        moinfo.sopi[h] = wfn->nsopi()[h];
        moinfo.orbspi[h] = wfn->nmopi()[h];
        moinfo.clsdpi[h] = wfn->doccpi()[h];
        moinfo.openpi[h] = wfn->soccpi()[h];
    }

    moinfo.natom = wfn->molecule()->natom();

    nirreps = moinfo.nirreps;

    moinfo.ntri = moinfo.nmo * (moinfo.nmo+1)/2;
    moinfo.noei = moinfo.nso * (moinfo.nso+1)/2;
    moinfo.noei_ao = moinfo.nao * (moinfo.nao+1)/2;

    /* Get frozen and active orbital lookups from CC_INFO */
    moinfo.frdocc = init_int_array(nirreps);
    moinfo.fruocc = init_int_array(nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep",
                    (char *) moinfo.frdocc, sizeof(int)*nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep",
                    (char *) moinfo.fruocc, sizeof(int)*nirreps);

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                    sizeof(int));
    moinfo.nactive = nactive;

    moinfo.nfzc = 0;
    for(h=0; h < nirreps; h++) moinfo.nfzc += moinfo.frdocc[h];

    if(params.ref == 2) { /** UHF **/

        moinfo.aoccpi = init_int_array(nirreps);
        moinfo.boccpi = init_int_array(nirreps);
        moinfo.avirtpi = init_int_array(nirreps);
        moinfo.bvirtpi = init_int_array(nirreps);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep",
                        (char *) moinfo.aoccpi, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep",
                        (char *) moinfo.boccpi, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep",
                        (char *) moinfo.avirtpi, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep",
                        (char *) moinfo.bvirtpi, sizeof(int)*moinfo.nirreps);

        moinfo.aocc_sym = init_int_array(nactive);
        moinfo.bocc_sym = init_int_array(nactive);
        moinfo.avir_sym = init_int_array(nactive);
        moinfo.bvir_sym = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry",
                        (char *) moinfo.aocc_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry",
                        (char *) moinfo.bocc_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry",
                        (char *) moinfo.avir_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry",
                        (char *) moinfo.bvir_sym, sizeof(int)*nactive);

        moinfo.aocc_off = init_int_array(moinfo.nirreps);
        moinfo.bocc_off = init_int_array(moinfo.nirreps);
        moinfo.avir_off = init_int_array(moinfo.nirreps);
        moinfo.bvir_off = init_int_array(moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets",
                        (char *) moinfo.aocc_off, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets",
                        (char *) moinfo.bocc_off, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets",
                        (char *) moinfo.avir_off, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets",
                        (char *) moinfo.bvir_off, sizeof(int)*moinfo.nirreps);

        moinfo.qt_aocc = init_int_array(nactive);
        moinfo.qt_bocc = init_int_array(nactive);
        moinfo.qt_avir = init_int_array(nactive);
        moinfo.qt_bvir = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order",
                        (char *) moinfo.qt_aocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order",
                        (char *) moinfo.qt_bocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order",
                        (char *) moinfo.qt_avir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order",
                        (char *) moinfo.qt_bvir, sizeof(int)*nactive);


    }
    else { /** RHF or ROHF **/

        moinfo.occpi = init_int_array(nirreps);
        moinfo.virtpi = init_int_array(nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep",
                        (char *) moinfo.occpi, sizeof(int)*nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep",
                        (char *) moinfo.virtpi, sizeof(int)*nirreps);

        moinfo.occ_sym = init_int_array(nactive);
        moinfo.vir_sym = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry",
                        (char *) moinfo.occ_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry",
                        (char *) moinfo.vir_sym, sizeof(int)*nactive);

        moinfo.occ_off = init_int_array(moinfo.nirreps);
        moinfo.vir_off = init_int_array(moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets",
                        (char *) moinfo.occ_off, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets",
                        (char *) moinfo.vir_off, sizeof(int)*moinfo.nirreps);

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo.qt_occ = init_int_array(nactive);
        moinfo.qt_vir = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Occ Order",
                        (char *) moinfo.qt_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Virt Order",
                        (char *) moinfo.qt_vir, sizeof(int)*nactive);

        moinfo.cc_occ = init_int_array(nactive);
        moinfo.cc_vir = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Occ Order",
                        (char *) moinfo.cc_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Virt Order",
                        (char *) moinfo.cc_vir, sizeof(int)*nactive);
    }

    /* Compute spatial-orbital reordering arrays */
    moinfo.pitzer2qt = init_int_array(moinfo.nmo);
    moinfo.qt2pitzer = init_int_array(moinfo.nmo);
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
               moinfo.pitzer2qt, moinfo.orbspi, moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++) {
        j = moinfo.pitzer2qt[i];
        moinfo.qt2pitzer[j] = i;
    }

    /* Adjust clsdpi array for frozen orbitals */
    for(i=0; i < nirreps; i++)
        moinfo.clsdpi[i] -= moinfo.frdocc[i];

    moinfo.uoccpi = init_int_array(moinfo.nirreps);
    for(i=0; i < nirreps; i++)
        moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
                moinfo.openpi[i] - moinfo.fruocc[i] -
                moinfo.frdocc[i];

    moinfo.nvirt = 0;
    for(i=0; i < nirreps; i++) moinfo.nvirt += moinfo.virtpi[i];

    /*** arrange active SCF MO's ***/
    actpi = init_int_array(nirreps);
    for(h=0; h < nirreps; h++)
        actpi[h] = moinfo.orbspi[h] - moinfo.frdocc[h] - moinfo.fruocc[h];
    moinfo.actpi = actpi;

    if(params.ref == 0 || params.ref == 1) /* RHF/ROHF */
        moinfo.scf = wfn->Ca()->to_block_matrix();
    else if(params.ref == 2) {  /* UHF */
        moinfo.scf_alpha = wfn->Ca()->to_block_matrix();
        moinfo.scf_beta = wfn->Cb()->to_block_matrix();
    }

    /* Get the active virtual orbitals */
    if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

        C = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo.sopi[h] && moinfo.virtpi[h]) {
                C[h] = block_matrix(moinfo.sopi[h],moinfo.virtpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
                        moinfo.sopi[h]*moinfo.virtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo.C = C;
    }

    /* Prepare memory for property integrals */
    moinfo.MU = (double ***) malloc(3 * sizeof(double **));
    moinfo.L = (double ***) malloc(3 * sizeof(double **));
    moinfo.Lcc = (double ***) malloc(3 * sizeof(double **));
    moinfo.P = (double ***) malloc(3 * sizeof(double **));
    moinfo.Pcc = (double ***) malloc(3 * sizeof(double **));
    moinfo.Q = (double ****) malloc(3 * sizeof(double ***));
    for(i=0; i < 3; i++)
        moinfo.Q[i] = (double ***) malloc(3 * sizeof(double **));
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup(void)
{
    int i;

    if(params.ref == 2) { /* UHF */
        free(moinfo.aoccpi);
        free(moinfo.boccpi);
        free(moinfo.avirtpi);
        free(moinfo.bvirtpi);
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
        free(moinfo.aocc_off);
        free(moinfo.bocc_off);
        free(moinfo.avir_off);
        free(moinfo.bvir_off);
        free(moinfo.qt_aocc);
        free(moinfo.qt_bocc);
        free(moinfo.qt_avir);
        free(moinfo.qt_bvir);
        free_block(moinfo.scf_alpha);
        free_block(moinfo.scf_beta);
    }
    else {
        for(i=0; i < moinfo.nirreps; i++)
            if(moinfo.sopi[i] && moinfo.virtpi[i]) free_block(moinfo.C[i]);
        free(moinfo.C);
        free(moinfo.occpi);
        free(moinfo.virtpi);
        free(moinfo.occ_sym);
        free(moinfo.vir_sym);
        free(moinfo.occ_off);
        free(moinfo.vir_off);
        free(moinfo.qt_occ);
        free(moinfo.qt_vir);
        free_block(moinfo.scf);
    }

    free(moinfo.sopi);
    free(moinfo.orbspi);
    free(moinfo.clsdpi);
    free(moinfo.openpi);
//    free(moinfo.uoccpi);
//    free(moinfo.fruocc);
//    free(moinfo.frdocc);
    free(moinfo.actpi);

    for(i=0; i < moinfo.nirreps; i++)
        free(moinfo.labels[i]);
    free(moinfo.labels);

    free(moinfo.MU);
    free(moinfo.L);
    free(moinfo.P);
    for(i=0; i < 3; i++)
        free(moinfo.Q[i]);
    free(moinfo.Q);
    free(moinfo.pitzer2qt);
    free(moinfo.qt2pitzer);

    free(moinfo.mu_irreps);
    free(moinfo.l_irreps);
}

}} // namespace psi::ccresponse
