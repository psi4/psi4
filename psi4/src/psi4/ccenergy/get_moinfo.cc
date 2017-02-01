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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "MOInfo.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified by TDC, March 1999
*/

void CCEnergyWavefunction::get_moinfo(void)
{
    int i, j, h, p, q, errcod, nactive, nirreps;
    double ***Co, ***Cv, ***Ca, ***Cb;
    psio_address next;

    moinfo_.nirreps = nirrep_;
    moinfo_.nmo = nmo_;
    moinfo_.nso = nso_;
    moinfo_.nao = basisset_->nao();
    moinfo_.labels = molecule_->irrep_labels();
    moinfo_.enuc = molecule_->nuclear_repulsion_energy();
    moinfo_.conv = 0.0;
    if(reference_wavefunction_)
        moinfo_.escf = reference_wavefunction_->reference_energy();
    else
        moinfo_.escf = energy_;
    moinfo_.sopi = nsopi_;
    moinfo_.orbspi = nmopi_;
    moinfo_.openpi = soccpi_;
    moinfo_.clsdpi = init_int_array(moinfo_.nirreps);
    // TODO: figure out why I can't just directly assign doccpi_ here, the same as soccpi_.
    for(int h = 0; h < moinfo_.nirreps; ++h)
        moinfo_.clsdpi[h] = doccpi_[h];

    nirreps = moinfo_.nirreps;

    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(params_.ref),
                    sizeof(int));

    /* Get frozen and active orbital lookups from CC_INFO */
    moinfo_.frdocc = init_int_array(nirreps);
    moinfo_.fruocc = init_int_array(nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep",
                    (char *) moinfo_.frdocc, sizeof(int)*nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep",
                    (char *) moinfo_.fruocc, sizeof(int)*nirreps);

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                    sizeof(int));

    if(params_.ref == 2) { /** UHF **/

        moinfo_.aoccpi = init_int_array(nirreps);
        moinfo_.boccpi = init_int_array(nirreps);
        moinfo_.avirtpi = init_int_array(nirreps);
        moinfo_.bvirtpi = init_int_array(nirreps);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep",
                        (char *) moinfo_.aoccpi, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep",
                        (char *) moinfo_.boccpi, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep",
                        (char *) moinfo_.avirtpi, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep",
                        (char *) moinfo_.bvirtpi, sizeof(int)*moinfo_.nirreps);

        moinfo_.aocc_sym = init_int_array(nactive);
        moinfo_.bocc_sym = init_int_array(nactive);
        moinfo_.avir_sym = init_int_array(nactive);
        moinfo_.bvir_sym = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry",
                        (char *) moinfo_.aocc_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry",
                        (char *) moinfo_.bocc_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry",
                        (char *) moinfo_.avir_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry",
                        (char *) moinfo_.bvir_sym, sizeof(int)*nactive);

        moinfo_.aocc_off = init_int_array(moinfo_.nirreps);
        moinfo_.bocc_off = init_int_array(moinfo_.nirreps);
        moinfo_.avir_off = init_int_array(moinfo_.nirreps);
        moinfo_.bvir_off = init_int_array(moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets",
                        (char *) moinfo_.aocc_off, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets",
                        (char *) moinfo_.bocc_off, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets",
                        (char *) moinfo_.avir_off, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets",
                        (char *) moinfo_.bvir_off, sizeof(int)*moinfo_.nirreps);

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo_.qt_aocc = init_int_array(nactive);
        moinfo_.qt_bocc = init_int_array(nactive);
        moinfo_.qt_avir = init_int_array(nactive);
        moinfo_.qt_bvir = init_int_array(nactive);
        moinfo_.cc_aocc = init_int_array(nactive);
        moinfo_.cc_bocc = init_int_array(nactive);
        moinfo_.cc_avir = init_int_array(nactive);
        moinfo_.cc_bvir = init_int_array(nactive);

        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order",
                        (char *) moinfo_.qt_aocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order",
                        (char *) moinfo_.qt_bocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order",
                        (char *) moinfo_.qt_avir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order",
                        (char *) moinfo_.qt_bvir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Occ Order",
                        (char *) moinfo_.cc_aocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Occ Order",
                        (char *) moinfo_.cc_bocc, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Alpha Active Virt Order",
                        (char *) moinfo_.cc_avir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Beta Active Virt Order",
                        (char *) moinfo_.cc_bvir, sizeof(int)*nactive);


    }
    else { /** RHF or ROHF **/

        moinfo_.occpi = init_int_array(nirreps);
        moinfo_.virtpi = init_int_array(nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep",
                        (char *) moinfo_.occpi, sizeof(int)*nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep",
                        (char *) moinfo_.virtpi, sizeof(int)*nirreps);
        moinfo_.occ_sym = init_int_array(nactive);
        moinfo_.vir_sym = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry",
                        (char *) moinfo_.occ_sym, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry",
                        (char *) moinfo_.vir_sym, sizeof(int)*nactive);

        moinfo_.occ_off = init_int_array(moinfo_.nirreps);
        moinfo_.vir_off = init_int_array(moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Occ Orb Offsets",
                        (char *) moinfo_.occ_off, sizeof(int)*moinfo_.nirreps);
        psio_read_entry(PSIF_CC_INFO, "Active Virt Orb Offsets",
                        (char *) moinfo_.vir_off, sizeof(int)*moinfo_.nirreps);

        /* Get CC->QT and QT->CC active occupied and virtual reordering arrays */
        moinfo_.qt_occ = init_int_array(nactive);
        moinfo_.qt_vir = init_int_array(nactive);
        moinfo_.cc_occ = init_int_array(nactive);
        moinfo_.cc_vir = init_int_array(nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Occ Order",
                        (char *) moinfo_.qt_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "CC->QT Active Virt Order",
                        (char *) moinfo_.qt_vir, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Occ Order",
                        (char *) moinfo_.cc_occ, sizeof(int)*nactive);
        psio_read_entry(PSIF_CC_INFO, "QT->CC Active Virt Order",
                        (char *) moinfo_.cc_vir, sizeof(int)*nactive);

    }

    /* Build sosym array (for AO-basis BT2) */
    moinfo_.sosym = init_int_array(moinfo_.nso);
    for(h=0,q=0; h < nirreps; h++)
        for(p=0; p < moinfo_.sopi[h]; p++)
            moinfo_.sosym[q++] = h;

    /* Get the active virtual orbitals */
    if(params_.ref == 0 || params_.ref == 1) { /** RHF/ROHF **/

        Co = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo_.sopi[h] && moinfo_.occpi[h]) {
                Co[h] = block_matrix(moinfo_.sopi[h],moinfo_.occpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) Co[h][0],
                        moinfo_.sopi[h]*moinfo_.occpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo_.Co = Co;

        Cv = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo_.sopi[h] && moinfo_.virtpi[h]) {
                Cv[h] = block_matrix(moinfo_.sopi[h],moinfo_.virtpi[h]);
                psio_read(PSIF_CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) Cv[h][0],
                        moinfo_.sopi[h]*moinfo_.virtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo_.Cv = Cv;
    }
    else if(params_.ref == 2) { /** UHF **/

        Ca = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo_.sopi[h] && moinfo_.avirtpi[h]) {
                Ca[h] = block_matrix(moinfo_.sopi[h],moinfo_.avirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) Ca[h][0],
                        moinfo_.sopi[h]*moinfo_.avirtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo_.Cav = Ca;


        Cb = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo_.sopi[h] && moinfo_.bvirtpi[h]) {
                Cb[h] = block_matrix(moinfo_.sopi[h],moinfo_.bvirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Beta Virtual Orbs", (char *) Cb[h][0],
                        moinfo_.sopi[h]*moinfo_.bvirtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo_.Cbv = Cb;
    }

    /* Compute spatial-orbital reordering arrays */
    if(params_.ref == 0 || params_.ref == 1) {
        moinfo_.pitzer2qt = init_int_array(moinfo_.nmo);
        moinfo_.qt2pitzer = init_int_array(moinfo_.nmo);
        reorder_qt(moinfo_.clsdpi, moinfo_.openpi, moinfo_.frdocc, moinfo_.fruocc,
                   moinfo_.pitzer2qt, moinfo_.orbspi, moinfo_.nirreps);
        for(i=0; i < moinfo_.nmo; i++) {
            j = moinfo_.pitzer2qt[i];
            moinfo_.qt2pitzer[j] = i;
        }
    }
    else if(params_.ref == 2) {
        moinfo_.pitzer2qt_a = init_int_array(moinfo_.nmo);
        moinfo_.qt2pitzer_a = init_int_array(moinfo_.nmo);
        moinfo_.pitzer2qt_b = init_int_array(moinfo_.nmo);
        moinfo_.qt2pitzer_b = init_int_array(moinfo_.nmo);
        reorder_qt_uhf(moinfo_.clsdpi, moinfo_.openpi, moinfo_.frdocc, moinfo_.fruocc,
                       moinfo_.pitzer2qt_a, moinfo_.pitzer2qt_b, moinfo_.orbspi,
                       moinfo_.nirreps);
        for(i=0; i < moinfo_.nmo; i++) {
            j = moinfo_.pitzer2qt_a[i];
            moinfo_.qt2pitzer_a[j] = i;
            j = moinfo_.pitzer2qt_b[i];
            moinfo_.qt2pitzer_b[j] = i;
        }
    }

    /* Adjust clsdpi array for frozen orbitals */
    for(i=0; i < nirreps; i++)
        moinfo_.clsdpi[i] -= moinfo_.frdocc[i];

    moinfo_.uoccpi = init_int_array(moinfo_.nirreps);
    for(i=0; i < nirreps; i++)
        moinfo_.uoccpi[i] = moinfo_.orbspi[i] - moinfo_.clsdpi[i] -
                moinfo_.openpi[i] - moinfo_.fruocc[i] -
                moinfo_.frdocc[i];

    if(params_.ref == 0) {
        moinfo_.nvirt = 0;
        for(h=0; h < nirreps; h++) moinfo_.nvirt += moinfo_.virtpi[h];
    }

    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(moinfo_.eref),
                    sizeof(double));

    outfile->Printf("\n    Nuclear Rep. energy (wfn)     = %20.15f\n",moinfo_.enuc);
    outfile->Printf(  "    SCF energy          (wfn)     = %20.15f\n",moinfo_.escf);
    outfile->Printf(  "    Reference energy    (file100) = %20.15f\n",moinfo_.eref);
}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void CCEnergyWavefunction::cleanup(void)
{
    int i, h;
    char *keyw=NULL;

    if( params_.wfn == "CC2" || params_.wfn == "EOM_CC2" )
        psio_write_entry(PSIF_CC_INFO, "CC2 Energy", (char *) &(moinfo_.ecc),
                         sizeof(double));
    else if( params_.wfn == "CC3"  || params_.wfn == "EOM_CC3" )
        psio_write_entry(PSIF_CC_INFO, "CC3 Energy", (char *) &(moinfo_.ecc),
                         sizeof(double));
    else
        psio_write_entry(PSIF_CC_INFO, "CCSD Energy", (char *) &(moinfo_.ecc),
                         sizeof(double));

    if(params_.ref == 0 || params_.ref == 1) {
        for(h=0; h < moinfo_.nirreps; h++){
            if(moinfo_.sopi[h] && moinfo_.occpi[h]) free_block(moinfo_.Co[h]);
            if(moinfo_.sopi[h] && moinfo_.virtpi[h]) free_block(moinfo_.Cv[h]);
        }
        free(moinfo_.Cv);
        free(moinfo_.Co);
    }
    else if(params_.ref == 2) {
        for(h=0; h < moinfo_.nirreps; h++)
            if(moinfo_.sopi[h] && moinfo_.avirtpi[h]) free_block(moinfo_.Cav[h]);
        free(moinfo_.Cav);
        for(h=0; h < moinfo_.nirreps; h++)
            if(moinfo_.sopi[h] && moinfo_.bvirtpi[h]) free_block(moinfo_.Cbv[h]);
        free(moinfo_.Cbv);
    }

    // Wavefunction owns these arrays
//    free(moinfo.sopi);
//    free(moinfo.sosym);
//    free(moinfo.orbspi);
    free(moinfo_.clsdpi);
//    free(moinfo.openpi);
//    free(moinfo.uoccpi);
    //  free(moinfo.fruocc);
    //  free(moinfo.frdocc);
    for(i=0; i < moinfo_.nirreps; i++)
        free(moinfo_.labels[i]);
    free(moinfo_.labels);
    if(params_.ref == 2) {
        free(moinfo_.aoccpi);
        free(moinfo_.boccpi);
        free(moinfo_.avirtpi);
        free(moinfo_.bvirtpi);
        free(moinfo_.aocc_sym);
        free(moinfo_.bocc_sym);
        free(moinfo_.avir_sym);
        free(moinfo_.bvir_sym);
        free(moinfo_.aocc_off);
        free(moinfo_.bocc_off);
        free(moinfo_.avir_off);
        free(moinfo_.bvir_off);
        free(moinfo_.qt_aocc);
        free(moinfo_.qt_bocc);
        free(moinfo_.qt_avir);
        free(moinfo_.qt_bvir);
        free(moinfo_.cc_aocc);
        free(moinfo_.cc_bocc);
        free(moinfo_.cc_avir);
        free(moinfo_.cc_bvir);
    }
    else {
        free(moinfo_.occpi);
        free(moinfo_.virtpi);
        free(moinfo_.occ_sym);
        free(moinfo_.vir_sym);
        free(moinfo_.occ_off);
        free(moinfo_.vir_off);
        free(moinfo_.qt_occ);
        free(moinfo_.qt_vir);
        free(moinfo_.cc_occ);
        free(moinfo_.cc_vir);
    }

}

}} // namespace psi::ccenergy
