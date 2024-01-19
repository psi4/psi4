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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "ccdensity.h"
#include "Frozen.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {
#include "psi4/physconst.h"

void ex_oscillator_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S, struct TD_Params *U,
                            struct XTD_Params *xtd_data) {
    double lt_x, lt_y, lt_z;
    double rt_x, rt_y, rt_z;
    double ds_x, ds_y, ds_z;
    double f_x, f_y, f_z;
    double f;

    lt_x = lt_y = lt_z = 0.0;
    rt_x = rt_y = rt_z = 0.0;
    ds_x = ds_y = ds_z = 0.0;
    f_x = f_y = f_z = 0.0;
    f = 0;

    /*** Transform the SO dipole integrals to the MO basis ***/

    MintsHelper mints(wfn.basisset(), Process::environment.options, 0);
    auto dipole = mints.so_dipole();

    Matrix MUX_MO, MUY_MO, MUZ_MO, MUX_MO_A, MUY_MO_A, MUZ_MO_A, MUX_MO_B, MUY_MO_B, MUZ_MO_B;
    if ((params.ref == 0) || (params.ref == 1)) {
        MUX_MO = linalg::triplet(*wfn.Ca(), *dipole[0], *wfn.Ca(), true, false, false);
        MUY_MO = linalg::triplet(*wfn.Ca(), *dipole[1], *wfn.Ca(), true, false, false);
        MUZ_MO = linalg::triplet(*wfn.Ca(), *dipole[2], *wfn.Ca(), true, false, false);
    } else if (params.ref == 2) {
        MUX_MO_A = linalg::triplet(*wfn.Ca(), *dipole[0], *wfn.Ca(), true, false, false);
        MUY_MO_A = linalg::triplet(*wfn.Ca(), *dipole[1], *wfn.Ca(), true, false, false);
        MUZ_MO_A = linalg::triplet(*wfn.Ca(), *dipole[2], *wfn.Ca(), true, false, false);
        MUX_MO_B = linalg::triplet(*wfn.Cb(), *dipole[0], *wfn.Cb(), true, false, false);
        MUY_MO_B = linalg::triplet(*wfn.Cb(), *dipole[1], *wfn.Cb(), true, false, false);
        MUZ_MO_B = linalg::triplet(*wfn.Cb(), *dipole[2], *wfn.Cb(), true, false, false);
    }

    outfile->Printf("\n\tOscillator Strength for %d%3s to %d%3s\n", S->root + 1, moinfo.labels[S->irrep].c_str(),
                    U->root + 1, moinfo.labels[U->irrep].c_str());
    outfile->Printf("\t                              X    \t       Y    \t       Z\n");

    if ((params.ref == 0) || (params.ref == 1)) {
        lt_x = MUX_MO.vector_dot(moinfo.ltd_mat);
        lt_y = MUY_MO.vector_dot(moinfo.ltd_mat);
        lt_z = MUZ_MO.vector_dot(moinfo.ltd_mat);
        rt_x = MUX_MO.vector_dot(moinfo.rtd_mat);
        rt_y = MUY_MO.vector_dot(moinfo.rtd_mat);
        rt_z = MUZ_MO.vector_dot(moinfo.rtd_mat);

    } else if (params.ref == 2) {
        lt_x = MUX_MO_A.vector_dot(moinfo.ltd_a_mat);
        lt_y = MUY_MO_A.vector_dot(moinfo.ltd_a_mat);
        lt_z = MUZ_MO_A.vector_dot(moinfo.ltd_a_mat);
        rt_x = MUX_MO_A.vector_dot(moinfo.rtd_a_mat);
        rt_y = MUY_MO_A.vector_dot(moinfo.rtd_a_mat);
        rt_z = MUZ_MO_A.vector_dot(moinfo.rtd_a_mat);
        lt_x += MUX_MO_B.vector_dot(moinfo.ltd_b_mat);
        lt_y += MUY_MO_B.vector_dot(moinfo.ltd_b_mat);
        lt_z += MUZ_MO_B.vector_dot(moinfo.ltd_b_mat);
        rt_x += MUX_MO_B.vector_dot(moinfo.rtd_b_mat);
        rt_y += MUY_MO_B.vector_dot(moinfo.rtd_b_mat);
        rt_z += MUZ_MO_B.vector_dot(moinfo.rtd_b_mat);
    }

    ds_x = lt_x * rt_x;
    ds_y = lt_y * rt_y;
    ds_z = lt_z * rt_z;

    /* Use |w2 - w1| for oscillator strengths */
    // We view excitation energies as positive,
    // so we want to substract the lower state's energy from the
    // higher state's.
    // U should be the higher-energy excited state.
    double delta_e = U->cceom_energy - S->cceom_energy;

    f_x = (2 * delta_e * ds_x) / 3;
    f_y = (2 * delta_e * ds_y) / 3;
    f_z = (2 * delta_e * ds_z) / 3;

    f = f_x + f_y + f_z;

    /* Fill in XTD_Params for this Transition */
    xtd_data->root1 = S->root;
    xtd_data->root2 = U->root;
    xtd_data->irrep1 = S->irrep;
    xtd_data->irrep2 = U->irrep;
    xtd_data->cceom_energy = delta_e;
    xtd_data->OS = f;

    /* Compute Einstein A,B Coefficients */
    double hartree2Hz = pc_hartree2MHz * (1.0e6);
    double hbar = pc_h / (pc_twopi);
    /* SI Dipole Strength */
    double ds_si = (ds_x + ds_y + ds_z) * pc_dipmom_au2si * pc_dipmom_au2si;
    /* SI Transition Energy */
    double nu_si = delta_e * hartree2Hz;
    /* Einstein Coefficients */
    double einstein_b = (2.0 / 3.0) * (pc_pi / pow(hbar, 2.0)) * (1.0 / (4.0 * pc_pi * pc_e0)) * ds_si;
    double einstein_a = 8.0 * pc_pi * pc_h * pow((nu_si / pc_c), 3.0) * einstein_b;
    if (einstein_a < 1e-7) einstein_a = 0.0000000;
    xtd_data->einstein_a = einstein_a;
    xtd_data->einstein_b = einstein_b;

    outfile->Printf("\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n", lt_x, lt_y, lt_z);
    outfile->Printf("\t<n|mu_e|0>              %11.8lf \t %11.8lf \t %11.8lf\n", rt_x, rt_y, rt_z);
    outfile->Printf("\tDipole Strength         %11.8lf \n", ds_x + ds_y + ds_z);
    outfile->Printf("\tOscillator Strength     %11.8lf \n", f_x + f_y + f_z);
    outfile->Printf("\tEinstein A Coefficient   %11.8e \n", einstein_a);
    outfile->Printf("\tEinstein B Coefficient   %11.8e \n", einstein_b);

    // Save variables to wfn.
    // Process::environment.globals["CCname ROOT m -> ROOT n OSCILLATOR STRENGTH (LEN)"]
    // Process::environment.globals["CCname ROOT m -> ROOT n OSCILLATOR STRENGTH (LEN) - h TRANSITION"]
    // Process::environment.globals["CCname ROOT m (h) -> ROOT n (i) OSCILLATOR STRENGTH (LEN)"]
    // Process::environment.globals["CCname ROOT m (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (LEN)"]
    // Process::environment.globals["CCname ROOT m -> ROOT n EINSTEIN A (LEN)"]
    // Process::environment.globals["CCname ROOT m -> ROOT n EINSTEIN A (LEN) - h TRANSITION"]
    // Process::environment.globals["CCname ROOT m (h) -> ROOT n (i) EINSTEIN A (LEN)"]
    // Process::environment.globals["CCname ROOT m (IN h) -> ROOT n (IN i) EINSTEIN A (LEN)"]
    // Process::environment.globals["CCname ROOT m -> ROOT n EINSTEIN B (LEN)"]
    // Process::environment.globals["CCname ROOT m -> ROOT n EINSTEIN B (LEN) - h TRANSITION"]
    // Process::environment.globals["CCname ROOT m (h) -> ROOT n (i) EINSTEIN B (LEN)"]
    // Process::environment.globals["CCname ROOT m (IN h) -> ROOT n (IN i) EINSTEIN B (LEN)"]
    scalar_saver_excited(wfn, S, U, "OSCILLATOR STRENGTH (LEN)", f_x + f_y + f_z);
    scalar_saver_excited(wfn, S, U, "EINSTEIN A (LEN)", einstein_a);
    scalar_saver_excited(wfn, S, U, "EINSTEIN B (LEN)", einstein_b);

    return;
}
}
}  // namespace psi
