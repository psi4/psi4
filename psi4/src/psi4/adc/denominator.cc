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

#include "psi4/libtrans/integraltransform.h"
#include <cmath>
#include "adc.h"

namespace psi{ namespace adc{

void
ADCWfn::shift_denom2(int root, int irrep, double omega)
{
    char lbl[32];
    dpdfile2 D, L;

    sprintf(lbl, "D_[%d]12", irrep);
    global_dpd_->file2_init(&D, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    global_dpd_->file2_mat_init(&D);
    global_dpd_->file2_mat_rd(&D);

    sprintf(lbl, "L^(%d)_[%d]12", root, irrep);
    global_dpd_->file2_init(&L, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    global_dpd_->file2_mat_init(&L);
    global_dpd_->file2_mat_rd(&L);

    for(int Isym = 0;Isym < nirrep_;Isym++){
        for(int i = 0;i < D.params->rowtot[Isym];i++){
            for(int a = 0;a < D.params->coltot[Isym^irrep];a++){
                double denom = omega - D.matrix[Isym][i][a];
                if(fabs(denom) > 1e-6)
                    L.matrix[Isym][i][a] = 1 / denom;
                else
                    L.matrix[Isym][i][a] = 0;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&L);
    global_dpd_->file2_mat_close(&L);
    global_dpd_->file2_close(&L);
    global_dpd_->file2_mat_close(&D);
    global_dpd_->file2_close(&D);
}

void
ADCWfn::shift_denom4(int irrep, double omega)
{
    char lbl[32];
    dpdbuf4 D;

    sprintf(lbl, "D_[%d]1234", irrep);
    global_dpd_->buf4_init(&D, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    for(int Gij = 0;Gij < nirrep_;Gij++){
        global_dpd_->buf4_mat_irrep_init(&D, Gij);

        for(int ij = 0;ij < D.params->rowtot[Gij];ij++){
            int i = D.params->roworb[Gij][ij][0];
            int j = D.params->roworb[Gij][ij][1];
            for(int ab = 0;ab < D.params->coltot[Gij^irrep];ab++){
                int a = D.params->colorb[Gij^irrep][ab][0];
                int b = D.params->colorb[Gij^irrep][ab][1];

                D.matrix[Gij][ij][ab] = 1.0 / (omega+aocce_[i]-avire_[a]+aocce_[j]-avire_[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, Gij);
        global_dpd_->buf4_mat_irrep_close(&D, Gij);
    }
    global_dpd_->buf4_close(&D);
}

}} // End Namespaces
