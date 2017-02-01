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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsio/psio.h"

#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/*
** get_moinfo():  Routine to obtain basic orbital information from
** CHKPT and CC_INFO.
**
** T. Daniel Crawford, October 1996
** Modified by TDC, March 1999
*/

void get_moinfo(std::shared_ptr<Wavefunction> wfn)
{
    int i, j, h, p, q, errcod, nactive, nirreps, sym;
    double ***C, ***Ca, ***Cb;
    psio_address next;

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.iopen = 0;
    for(int h=0; h < moinfo.nirreps; h++)
      moinfo.iopen += wfn->nsopi()[h];
    moinfo.irr_labs = wfn->molecule()->irrep_labels();
    moinfo.enuc = wfn->molecule()->nuclear_repulsion_energy();
    if(wfn->reference_wavefunction())
        moinfo.escf = wfn->reference_wavefunction()->reference_energy();
    else
        moinfo.escf = wfn->reference_energy();

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

    sym = 0;
    for (i=0;i<moinfo.nirreps;++i)
        for (j=0;j<moinfo.openpi[i];++j)
            sym = sym ^ i;
    moinfo.sym = sym;

    nirreps = moinfo.nirreps;

    moinfo.irr_labs_lowercase = (char **) malloc(sizeof(char *) * nirreps);
    for(i=0; i < nirreps; i++) {
      moinfo.irr_labs_lowercase[i] = (char *) malloc(4 * sizeof(char));
      moinfo.irr_labs_lowercase[i][0] = std::tolower(moinfo.irr_labs[i][0]);
      moinfo.irr_labs_lowercase[i][1] = std::tolower(moinfo.irr_labs[i][1]);
      moinfo.irr_labs_lowercase[i][2] = std::tolower(moinfo.irr_labs[i][2]);
      moinfo.irr_labs_lowercase[i][3] = '\0';
    }

    psio_read_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(params.ref),
                    sizeof(int));

    /* Get frozen and active orbital lookups from CC_INFO */
    moinfo.frdocc = init_int_array(nirreps);
    moinfo.fruocc = init_int_array(nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep",
                    (char *) moinfo.frdocc, sizeof(int)*nirreps);
    psio_read_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep",
                    (char *) moinfo.fruocc, sizeof(int)*nirreps);

    psio_read_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                    sizeof(int));

    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

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
    }

    else { /** UHF **/

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
    }

    /* Build sosym array (for AO-basis BT2) */
    moinfo.sosym = init_int_array(moinfo.nso);
    for(h=0,q=0; h < nirreps; h++)
        for(p=0; p < moinfo.sopi[h]; p++)
            moinfo.sosym[q++] = h;

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
    else if(params.ref == 2) { /** UHF **/

        Ca = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo.sopi[h] && moinfo.avirtpi[h]) {
                Ca[h] = block_matrix(moinfo.sopi[h],moinfo.avirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) Ca[h][0],
                        moinfo.sopi[h]*moinfo.avirtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo.Ca = Ca;


        Cb = (double ***) malloc(nirreps * sizeof(double **));
        next = PSIO_ZERO;
        for(h=0; h < nirreps; h++) {
            if(moinfo.sopi[h] && moinfo.bvirtpi[h]) {
                Cb[h] = block_matrix(moinfo.sopi[h],moinfo.bvirtpi[h]);
                psio_read(PSIF_CC_INFO, "UHF Active Beta Virtual Orbs", (char *) Cb[h][0],
                        moinfo.sopi[h]*moinfo.bvirtpi[h]*sizeof(double), next, &next);
            }
        }
        moinfo.Cb = Cb;
    }

    /* Adjust clsdpi array for frozen orbitals */
    for(i=0; i < nirreps; i++)
        moinfo.clsdpi[i] -= moinfo.frdocc[i];

    moinfo.uoccpi = init_int_array(moinfo.nirreps);
    for(i=0; i < nirreps; i++)
        moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
                moinfo.openpi[i] - moinfo.fruocc[i] -
                moinfo.frdocc[i];

    if(params.ref == 0) {
        moinfo.nvirt = 0;
        for(h=0; h < moinfo.nirreps; h++) moinfo.nvirt += moinfo.virtpi[h];
    }

    psio_read_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
                    sizeof(double));

    outfile->Printf("\n\tNuclear Rep. energy (wfn)     = %20.15f\n",moinfo.enuc);
    outfile->Printf(  "\tSCF energy          (wfn)     = %20.15f\n",moinfo.escf);
    outfile->Printf(  "\tReference energy    (file100) = %20.15f\n",moinfo.eref);


}

/* Frees memory allocated in get_moinfo() and dumps out the energy. */
void cleanup(void)
{
    int i, h;

    if(params.ref == 0 || params.ref == 1) {
        for(h=0; h < moinfo.nirreps; h++)
            if(moinfo.sopi[h] && moinfo.virtpi[h]) free_block(moinfo.C[h]);
        free(moinfo.C);
    }
    else if(params.ref == 2) {
        for(h=0; h < moinfo.nirreps; h++)
            if(moinfo.sopi[h] && moinfo.avirtpi[h]) free_block(moinfo.Ca[h]);
        free(moinfo.Ca);
        for(h=0; h < moinfo.nirreps; h++)
            if(moinfo.sopi[h] && moinfo.bvirtpi[h]) free_block(moinfo.Cb[h]);
        free(moinfo.Cb);
    }

    free(moinfo.sopi);
    free(moinfo.orbspi);
//    free(moinfo.sosym);
    free(moinfo.clsdpi);
    free(moinfo.openpi);
//    free(moinfo.uoccpi);
//    free(moinfo.fruocc);
//    free(moinfo.frdocc);
    for(i=0; i < moinfo.nirreps; i++)
        free(moinfo.irr_labs[i]);
    free(moinfo.irr_labs);
    for(i=0; i < moinfo.nirreps; i++)
        free(moinfo.irr_labs_lowercase[i]);
    free(moinfo.irr_labs_lowercase);
    if(params.ref == 2) {
        free(moinfo.aoccpi);
        free(moinfo.boccpi);
        free(moinfo.avirtpi);
        free(moinfo.bvirtpi);
        free(moinfo.aocc_sym);
        free(moinfo.bocc_sym);
        free(moinfo.avir_sym);
        free(moinfo.bvir_sym);
    }
    else {
        free(moinfo.occpi);
        free(moinfo.virtpi);
        free(moinfo.occ_sym);
        free(moinfo.vir_sym);
    }

}


}} // namespace psi::cceom
