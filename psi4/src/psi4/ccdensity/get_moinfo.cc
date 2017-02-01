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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <string.h>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn)
{
    int i, j, h, errcod;
    int nactive;
    double **scf_pitzer;

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.labels = wfn->molecule()->irrep_labels();
    moinfo.enuc = wfn->molecule()->nuclear_repulsion_energy();
    if(wfn->reference_wavefunction())
        moinfo.escf = wfn->reference_wavefunction()->reference_energy();
    else
        moinfo.escf = wfn->reference_energy();

    moinfo.orbspi = init_int_array(moinfo.nirreps);
    moinfo.clsdpi = init_int_array(moinfo.nirreps);
    moinfo.openpi = init_int_array(moinfo.nirreps);
    for(int h = 0; h < moinfo.nirreps; ++h){
        moinfo.orbspi[h] = wfn->nmopi()[h];
        moinfo.clsdpi[h] = wfn->doccpi()[h];
        moinfo.openpi[h] = wfn->soccpi()[h];
    }
    scf_pitzer = wfn->Ca()->to_block_matrix();

    moinfo.sym = 0;
    for (i=0;i<moinfo.nirreps;++i)
        for (j=0;j<moinfo.openpi[i];++j)
            moinfo.sym = moinfo.sym ^ i;

    /* Get frozen and active orbital lookups from CC_INFO */
    moinfo.frdocc = init_int_array(moinfo.nirreps);
    moinfo.fruocc = init_int_array(moinfo.nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep",
                    (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep",
                    (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                    sizeof(int));
    moinfo.nactive = nactive;

    if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

        moinfo.occpi = init_int_array(moinfo.nirreps);
        moinfo.virtpi = init_int_array(moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep",
                        (char *) moinfo.occpi, sizeof(int)*moinfo.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep",
                        (char *) moinfo.virtpi, sizeof(int)*moinfo.nirreps);

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

    }
    else if(params.ref == 2) { /** UHF **/
        moinfo.aoccpi = init_int_array(moinfo.nirreps);
        moinfo.boccpi = init_int_array(moinfo.nirreps);
        moinfo.avirtpi = init_int_array(moinfo.nirreps);
        moinfo.bvirtpi = init_int_array(moinfo.nirreps);

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
    for(i=0; i < moinfo.nirreps; i++)
        moinfo.clsdpi[i] -= moinfo.frdocc[i];

    moinfo.uoccpi = init_int_array(moinfo.nirreps);
    for(i=0; i < moinfo.nirreps; i++)
        moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
                moinfo.openpi[i] - moinfo.fruocc[i] -
                moinfo.frdocc[i];

    moinfo.nfzc = moinfo.nfzv = moinfo.nclsd = moinfo.nopen = moinfo.nuocc = 0;
    for(h=0; h < moinfo.nirreps; h++) {
        moinfo.nfzc += moinfo.frdocc[h];
        moinfo.nfzv += moinfo.fruocc[h];
        moinfo.nclsd += moinfo.clsdpi[h];
        moinfo.nopen += moinfo.openpi[h];
        moinfo.nuocc += moinfo.uoccpi[h];
    }

    if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo.qt_occ = init_int_array(nactive);
        moinfo.qt_vir = init_int_array(nactive);
        moinfo.cc_occ = init_int_array(nactive);
        moinfo.cc_vir = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Occ Order",
                        (char *) moinfo.qt_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Virt Order",
                        (char *) moinfo.qt_vir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Occ Order",
                        (char *) moinfo.cc_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Virt Order",
                        (char *) moinfo.cc_vir, sizeof(int)*nactive);

        /* Sort SCF MOs to QT order */
        moinfo.scf_qt = block_matrix(moinfo.nmo, moinfo.nmo);
        int I;
        for(i=0; i < moinfo.nmo; i++) {
            I = moinfo.pitzer2qt[i];  /* Pitzer --> QT */
            for(j=0; j < moinfo.nmo; j++) moinfo.scf_qt[j][I] = scf_pitzer[j][i];
        }
        free_block(scf_pitzer);
    }
    else if(params.ref == 2) { /** UHF **/

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo.qt_aocc = init_int_array(nactive);
        moinfo.qt_bocc = init_int_array(nactive);
        moinfo.qt_avir = init_int_array(nactive);
        moinfo.qt_bvir = init_int_array(nactive);
        moinfo.cc_aocc = init_int_array(nactive);
        moinfo.cc_bocc = init_int_array(nactive);
        moinfo.cc_avir = init_int_array(nactive);
        moinfo.cc_bvir = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order",
                        (char *) moinfo.qt_aocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order",
                        (char *) moinfo.qt_bocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order",
                        (char *) moinfo.qt_avir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order",
                        (char *) moinfo.qt_bvir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Occ Order",
                        (char *) moinfo.cc_aocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Occ Order",
                        (char *) moinfo.cc_bocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Virt Order",
                        (char *) moinfo.cc_avir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Virt Order",
                        (char *) moinfo.cc_bvir, sizeof(int)*nactive);

    }

    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
                    sizeof(double));

    outfile->Printf("\n\tNuclear Rep. energy (wfn)     = %20.15f\n",moinfo.enuc);
    outfile->Printf(  "\tSCF energy          (wfn)     = %20.15f\n",moinfo.escf);
    outfile->Printf(  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);

    if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
        psio_read_entry(PSIF_CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                        sizeof(double));
        outfile->Printf(  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
        outfile->Printf(  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
                moinfo.eref+moinfo.ecc);
    }
    else if( params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                        sizeof(double));
        outfile->Printf(  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
        outfile->Printf(  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
                moinfo.eref+moinfo.ecc);
    }
    else if(params.wfn == "CCSD_T") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc), sizeof(double));
        psio_read_entry(PSIF_CC_INFO, "(T) Energy", (char *) &(moinfo.et), sizeof(double));
        outfile->Printf(  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
        outfile->Printf(  "\t(T) energy          (CC_INFO) = %20.15f\n",moinfo.et);
        outfile->Printf(  "\tTotal CCSD(T) energy(CC_INFO) = %20.15f\n",
                moinfo.eref+moinfo.ecc+moinfo.et);
    }
    else if(params.wfn == "CC3" || params.wfn == "EOM_CC3") {
        psio_read_entry(PSIF_CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                        sizeof(double));
        outfile->Printf(  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
        outfile->Printf(  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
                moinfo.eref+moinfo.ecc);
    }


}

/* Frees memory allocated in get_moinfo(). */
void cleanup(void)
{
    int i;

    free(moinfo.orbspi);
    free(moinfo.clsdpi);
    free(moinfo.openpi);
    //  free(moinfo.uoccpi);
    //  free(moinfo.fruocc);
    //  free(moinfo.frdocc);
    for(i=0; i < moinfo.nirreps; i++)
        free(moinfo.labels[i]);
    free(moinfo.labels);

    if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
        free(moinfo.occ_sym);
        free(moinfo.vir_sym);
        free(moinfo.occ_off);
        free(moinfo.vir_off);
        free(moinfo.occpi);
        free(moinfo.virtpi);
        free(moinfo.qt_occ);
        free(moinfo.qt_vir);
        free(moinfo.cc_occ);
        free(moinfo.cc_vir);
        free(moinfo.pitzer2qt);
        free(moinfo.qt2pitzer);
        free_block(moinfo.scf_qt);
    }
    else if(params.ref == 2) { /** UHF **/
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
        free(moinfo.aocc_off);
        free(moinfo.bocc_off);
        free(moinfo.avir_off);
        free(moinfo.bvir_off);
        free(moinfo.aoccpi);
        free(moinfo.boccpi);
        free(moinfo.avirtpi);
        free(moinfo.bvirtpi);
        free(moinfo.qt_aocc);
        free(moinfo.qt_bocc);
        free(moinfo.qt_avir);
        free(moinfo.qt_bvir);
        free(moinfo.cc_aocc);
        free(moinfo.cc_bocc);
        free(moinfo.cc_avir);
        free(moinfo.cc_bvir);
    }
}


}} // namespace psi::ccdensity
