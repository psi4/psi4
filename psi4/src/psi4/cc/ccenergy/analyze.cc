/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::analyze() {
    double value, value2;

    dpdfile2 T1;
    dpdbuf4 T2;

    auto num_div = 500;
    auto max = 9.0;
    auto min = 0.0;
    auto width = (max - min) / (num_div);
    auto mode = std::ostream::app;
    auto printer = std::make_shared<PsiOutStream>("tamps.dat", mode);
    auto amp_array = init_array(num_div);

    auto nvir = moinfo_.virtpi[0];
    auto nocc = moinfo_.occpi[0];
    auto nso = moinfo_.nso;

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_mat_irrep_init(&T2, 0);
    global_dpd_->buf4_mat_irrep_rd(&T2, 0);
    auto T2trans = block_matrix(static_cast<size_t>(nocc) * nocc, static_cast<size_t>(nso) * nso);
    auto tmp = block_matrix(nvir, nso);
    auto tot1 = 0;
    auto tot2 = 0;
    auto position = 0;
    for (int ij = 0; ij < T2.params->rowtot[0]; ij++) {
        C_DGEMM('n', 't', nvir, nso, nvir, 1.0, &(T2.matrix[0][ij][0]), nvir, &(moinfo_.Cv[0][0][0]), nvir, 0.0,
                &(tmp[0][0]), nso);
        C_DGEMM('n', 'n', nso, nso, nvir, 1.0, &(moinfo_.Cv[0][0][0]), nvir, tmp[0], nso, 0.0, T2trans[ij], nso);

        for (int ab = 0; ab < nso * nso; ab++) {
            value = std::fabs(std::log10(std::fabs(T2trans[ij][ab])));
            tot2++;
            if ((value >= max) && (value <= (max + width))) {
                amp_array[num_div - 1]++;
                tot1++;
            } else if ((value <= min) && (value >= (min - width))) {
                amp_array[0]++;
                tot1++;
            } else if ((value < max) && (value > min)) {
                position = (int)floor((value - min) / width);
                amp_array[position]++;
                tot1++;
            }
        }
    }
    global_dpd_->buf4_mat_irrep_close(&T2, 0);
    global_dpd_->buf4_close(&T2);
    free_block(tmp);
    free_block(T2trans);

    value2 = 0;
    for (int i = num_div - 1; i >= 0; i--) {
        value = amp_array[i] / tot1;
        value2 += value;
        printer->Printf("%10.5lf %lf\n", -((i)*width) - min, value);
    }
    free(amp_array);
    printf("Total number of converged T2 amplitudes = %d\n", tot2);
    printf("Number of T2 amplitudes in analysis= %d\n", tot1);

    num_div = 40;
    max = 2;
    min = -5;
    width = (max - min) / (num_div);
    auto mode2 = std::ostream::app;
    auto printer2 = std::make_shared<PsiOutStream>("t1amps.dat", mode2);
    amp_array = init_array(num_div);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_print(&T1, "outfile");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    tot1 = tot2 = 0;
    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nso; a++) {
            value = std::log10(std::fabs(T1.matrix[0][i][a]));
            tot2++;
            if ((value >= max) && (value <= (max + width))) {
                amp_array[num_div - 1]++;
                tot1++;
            } else if ((value <= min) && (value >= (min - width))) {
                amp_array[0]++;
                tot1++;
            } else if ((value < max) && (value > min)) {
                position = (int)floor((value - min) / width);
                amp_array[position]++;
                tot1++;
            }
        }
    }

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);

    value2 = 0;
    for (int i = num_div - 1; i >= 0; i--) {
        value = amp_array[i] / tot1;
        value2 += value;
        printer2->Printf("%10.5lf %lf\n", ((i)*width) - min, value);
    }

    free(amp_array);
}

}  // namespace ccenergy
}  // namespace psi
