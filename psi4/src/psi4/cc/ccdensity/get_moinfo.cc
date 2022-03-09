/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

namespace psi {
namespace ccdensity {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CC_INFO.
**
** T. Daniel Crawford, October 1996.
** Modified by TDC, March 1999.
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn) {
    int i, j, h, errcod;
    int nactive;
    double **scf_pitzer;

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.labels = wfn->molecule()->irrep_labels();
    moinfo.enuc = wfn->molecule()->nuclear_repulsion_energy(wfn->get_dipole_field_strength());
    if (wfn->reference_wavefunction())
        moinfo.escf = wfn->reference_wavefunction()->energy();
    else
        moinfo.escf = wfn->energy();

    moinfo.orbspi = wfn->nmopi();
    moinfo.clsdpi = wfn->doccpi();
    moinfo.openpi = wfn->soccpi();

    moinfo.Ca = wfn->Ca();
    scf_pitzer = wfn->Ca()->to_block_matrix();

    moinfo.sym = 0;
    for (i = 0; i < moinfo.nirreps; ++i)
        for (j = 0; j < moinfo.openpi[i]; ++j) moinfo.sym = moinfo.sym ^ i;

    auto temp = std::vector<int>(moinfo.nirreps);

    /* Get frozen and active orbital lookups from CC_INFO */
    psio_read_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep", (char *)temp.data(), sizeof(int) * moinfo.nirreps);
    moinfo.frdocc = Dimension(temp);
    psio_read_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep", (char *)temp.data(), sizeof(int) * moinfo.nirreps);
    moinfo.fruocc = Dimension(temp);

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *)&(nactive), sizeof(int));
    moinfo.nactive = nactive;

    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

        psio_read_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep", (char *)temp.data(), sizeof(int) * moinfo.nirreps);
        moinfo.occpi = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.virtpi = Dimension(temp);

        moinfo.occ_sym = std::vector<int>(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *)moinfo.occ_sym.data(), sizeof(int) * nactive);
        moinfo.vir_sym = std::vector<int>(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *)moinfo.vir_sym.data(), sizeof(int) * nactive);

        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets", (char *)temp.data(), sizeof(int) * moinfo.nirreps);
        moinfo.occ_off = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets", (char *)temp.data(), sizeof(int) * moinfo.nirreps);
        moinfo.vir_off = Dimension(temp);

    } else if (params.ref == 2) { /** UHF **/

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.aoccpi = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.boccpi = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.avirtpi = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.bvirtpi = Dimension(temp);

        moinfo.aocc_sym = std::vector<int>(nactive);;
        moinfo.bocc_sym = std::vector<int>(nactive);;
        moinfo.avir_sym = std::vector<int>(nactive);;
        moinfo.bvir_sym = std::vector<int>(nactive);;

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry", (char *)moinfo.aocc_sym.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry", (char *)moinfo.bocc_sym.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry", (char *)moinfo.avir_sym.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry", (char *)moinfo.bvir_sym.data(), sizeof(int) * nactive);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.aocc_off = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.bocc_off = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.avir_off = Dimension(temp);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets", (char *)temp.data(),
                        sizeof(int) * moinfo.nirreps);
        moinfo.bvir_off = Dimension(temp);
    }

    /* Compute spatial-orbital reordering arrays */
    moinfo.pitzer2qt = std::vector<int>(moinfo.nmo);
    moinfo.qt2pitzer = std::vector<int>(moinfo.nmo);
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, moinfo.pitzer2qt.data(), moinfo.orbspi,
               moinfo.nirreps);
    for (i = 0; i < moinfo.nmo; i++) {
        j = moinfo.pitzer2qt[i];
        moinfo.qt2pitzer[j] = i;
    }

    /* Adjust clsdpi array for frozen orbitals */
    moinfo.clsdpi -= moinfo.frdocc;

    moinfo.uoccpi = moinfo.orbspi - moinfo.clsdpi - moinfo.openpi - moinfo.fruocc - moinfo.frdocc;

    moinfo.nfzc = moinfo.frdocc.sum();
    moinfo.nfzv = moinfo.fruocc.sum();
    moinfo.nclsd = moinfo.clsdpi.sum();
    moinfo.nopen = moinfo.openpi.sum();
    moinfo.nuocc = moinfo.uoccpi.sum();

    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo.qt_occ = std::vector<int>(nactive);
        moinfo.qt_vir = std::vector<int>(nactive);
        moinfo.cc_occ = std::vector<int>(nactive);
        moinfo.cc_vir = std::vector<int>(nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Occ Order", (char *)moinfo.qt_occ.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Virt Order", (char *)moinfo.qt_vir.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Occ Order", (char *)moinfo.cc_occ.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Virt Order", (char *)moinfo.cc_vir.data(), sizeof(int) * nactive);

        /* Sort SCF MOs to QT order */
        moinfo.scf_qt = block_matrix(moinfo.nmo, moinfo.nmo);
        int I;
        for (i = 0; i < moinfo.nmo; i++) {
            I = moinfo.pitzer2qt[i]; /* Pitzer --> QT */
            for (j = 0; j < moinfo.nmo; j++) moinfo.scf_qt[j][I] = scf_pitzer[j][i];
        }
        free_block(scf_pitzer);
    } else if (params.ref == 2) { /** UHF **/

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo.qt_aocc = std::vector<int>(nactive);
        moinfo.qt_bocc = std::vector<int>(nactive);
        moinfo.qt_avir = std::vector<int>(nactive);
        moinfo.qt_bvir = std::vector<int>(nactive);
        moinfo.cc_aocc = std::vector<int>(nactive);
        moinfo.cc_bocc = std::vector<int>(nactive);
        moinfo.cc_avir = std::vector<int>(nactive);
        moinfo.cc_bvir = std::vector<int>(nactive);

        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order", (char *)moinfo.qt_aocc.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order", (char *)moinfo.qt_bocc.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order", (char *)moinfo.qt_avir.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order", (char *)moinfo.qt_bvir.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Occ Order", (char *)moinfo.cc_aocc.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Occ Order", (char *)moinfo.cc_bocc.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Virt Order", (char *)moinfo.cc_avir.data(), sizeof(int) * nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Virt Order", (char *)moinfo.cc_bvir.data(), sizeof(int) * nactive);
    }

    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *)&(moinfo.eref), sizeof(double));

    outfile->Printf("\n\tNuclear Rep. energy (wfn)     = %20.15f\n", moinfo.enuc);
    outfile->Printf("\tSCF energy          (wfn)     = %20.15f\n", moinfo.escf);
    outfile->Printf("\tReference energy    (file100) = %20.15f\n", moinfo.eref);

    if (params.wfn == "CC2" || params.wfn == "EOM_CC2") {
        psio_read_entry(PSIF_CC_INFO, "CC2 Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCC2 energy          (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CC2 energy    (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    } else if (params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCCSD energy         (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CCSD energy   (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    } else if (params.wfn == "CCSD_T") {
        psio_read_entry(PSIF_CC_INFO, "CCSD Energy", (char *)&(moinfo.ecc), sizeof(double));
        psio_read_entry(PSIF_CC_INFO, "(T) Energy", (char *)&(moinfo.et), sizeof(double));
        outfile->Printf("\tCCSD energy         (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\t(T) energy          (CC_INFO) = %20.15f\n", moinfo.et);
        outfile->Printf("\tTotal CCSD(T) energy(CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc + moinfo.et);
    } else if (params.wfn == "CC3" || params.wfn == "EOM_CC3") {
        psio_read_entry(PSIF_CC_INFO, "CC3 Energy", (char *)&(moinfo.ecc), sizeof(double));
        outfile->Printf("\tCC3 energy          (CC_INFO) = %20.15f\n", moinfo.ecc);
        outfile->Printf("\tTotal CC3 energy    (CC_INFO) = %20.15f\n", moinfo.eref + moinfo.ecc);
    }
}

/* Frees memory allocated in get_moinfo(). */
void cleanup() {
    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
        free_block(moinfo.scf_qt);
    }
}

}  // namespace ccdensity
}  // namespace psi
