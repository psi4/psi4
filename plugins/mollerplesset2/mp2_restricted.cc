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

#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include "psifiles.h"
#define EXTERN
#include "mp2.h"
#include "libdpd/dpd.gbl"

namespace psi{ namespace mollerplesset2{

double plugin_mp2_restricted(SharedWavefunction wfn, Options &options)
{
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted,
               IntegralTransform::DPDOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false);
    ints.set_dpd_id(0);
    ints.set_keep_iwl_so_ints(true);
    ints.set_keep_dpd_so_ints(true);
    ints.initialize();
    ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    double aaE2 = 0.0, abE2 = 0.0 , bbE2 = 0.0, e2 = 0.0;
    // The alpha-alpha and beta-beta spin cases
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), 
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    for(int h = 0; h < nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
                aaE2 += pow(K.matrix[h][ij][ab], 2) / 
                     (aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    bbE2 = aaE2;

    // The alpha-beta spin case
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), 
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
                abE2 += pow(K.matrix[h][ij][ab], 2) / 
                     ( aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    e2 = aaE2 + abE2 + bbE2;
    outfile->Printf("\n\n\t\tAA correlation energy = %20.16f\n", aaE2);
    outfile->Printf("\t\tAB correlation energy = %20.16f\n", abE2);
    outfile->Printf("\t\tBB correlation energy = %20.16f\n", bbE2);
  
    psio->close(PSIF_LIBTRANS_DPD, 1);

    return e2;
}

}} // End Namespaces
