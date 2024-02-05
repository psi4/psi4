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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
/*
**  CCEOM: Program to calculate the EOM CCSD right-hand eigenvector and  energy
*/

#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "globals.h"

#include "psi4/cc/ccwave.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"

#include <cstdio>
#include <cstdlib>
#include <string>

namespace psi {
namespace cceom {

void init_io();
void get_moinfo(std::shared_ptr<Wavefunction>);
void cleanup();
void exit_io();
void diag(ccenergy::CCEnergyWavefunction&);
void get_params(Options &);
void get_eom_params(std::shared_ptr<Wavefunction>, Options &);
void form_dpd_dp();
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void sort_amps();
void hbar_norms();

/* local correlation functions */
void local_init();
void local_done();

}  // namespace cceom
}  // namespace psi

namespace psi {
namespace cceom {

PsiReturnType cceom(std::shared_ptr<ccenergy::CCEnergyWavefunction> ref_wfn, Options &options) {
    int i, h, done = 0, *cachefiles, **cachelist;
    init_io();
    outfile->Printf("\n\t**********************************************************\n");
    outfile->Printf("\t*  CCEOM: An Equation of Motion Coupled Cluster Program  *\n");
    outfile->Printf("\t**********************************************************\n");

    get_moinfo(ref_wfn);

    get_params(options);
    get_eom_params(ref_wfn, options);

    form_dpd_dp();

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if (params.ref == 2) { /* UHF */
        cachelist = cacheprep_uhf(params.cachelev, cachefiles);
        /* cachelist = init_int_matrix(32,32); */

        std::vector<int *> spaces;
        spaces.push_back(moinfo.aoccpi);
        spaces.push_back(moinfo.aocc_sym);
        spaces.push_back(moinfo.avirtpi);
        spaces.push_back(moinfo.avir_sym);
        spaces.push_back(moinfo.boccpi);
        spaces.push_back(moinfo.bocc_sym);
        spaces.push_back(moinfo.bvirtpi);
        spaces.push_back(moinfo.bvir_sym);
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 4, spaces);
    } else { /* RHF or ROHF */
        cachelist = cacheprep_rhf(params.cachelev, cachefiles);
        /* cachelist = init_int_matrix(12,12); */

        std::vector<int *> spaces;
        spaces.push_back(moinfo.occpi);
        spaces.push_back(moinfo.occ_sym);
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym);
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 2, spaces);
    }

    if (params.local) local_init();

    diag(*ref_wfn);

    dpd_close(0);
    if (params.local) local_done();
    cleanup();
    exit_io();
    return Success;
}

void init_io() {
    timer_on("cceom");
    for (int i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);
}

void exit_io() {
    int i;
    for (i = PSIF_CC_MIN; i <= PSIF_CC_DIIS_AMP; i++) psio_close(i, 1);
    for (i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i, 0);
    for (i = PSIF_CC_TMP11 + 1; i <= PSIF_CC_MAX; i++) psio_close(i, 1);
    timer_off("cceom");
}

void form_dpd_dp() {
    int h, h0, h1, cnt, nirreps;
    nirreps = moinfo.nirreps;

    dpd_dp = (int ***)malloc(nirreps * sizeof(int **));
    for (h = 0; h < nirreps; h++) {
        dpd_dp[h] = init_int_matrix(nirreps, 2);
        cnt = 0;
        for (h0 = 0; h0 < nirreps; h0++) {
            for (h1 = 0; h1 < nirreps; h1++) {
                if ((h0 ^ h1) == h) {
                    dpd_dp[h][cnt][0] = h0;
                    dpd_dp[h][cnt++][1] = h1;
                }
            }
        }
    }
}

}  // namespace cceom
}  // namespace psi
