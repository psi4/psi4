/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
/*
**  ccresponse: Program to compute CC linear response properties.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/psi4-dec.h"
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "globals.h"

namespace psi {
namespace ccresponse {

/* Function prototypes */
void init_io();
void title();
void get_moinfo(std::shared_ptr<Wavefunction>);
void get_params(std::shared_ptr<Wavefunction>, Options &);
void cleanup();
void exit_io();
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void hbar_extra();
void cc2_hbar_extra();
void sort_lamps();
void lambda_residuals();

void local_init();
void local_done();

void polar(std::shared_ptr<Wavefunction> ref_wfn);
void optrot(std::shared_ptr<Wavefunction> ref_wfn);
void roa(std::shared_ptr<Wavefunction> ref_wfn);

void preppert(std::shared_ptr<BasisSet> primary);

PsiReturnType ccresponse(std::shared_ptr<Wavefunction> ref_wfn, Options &options) {
    int **cachelist, *cachefiles;

    init_io();
    title();
    get_moinfo(ref_wfn);
    get_params(ref_wfn, options);

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if (params.ref == 2) { /*** UHF references ***/
        cachelist = cacheprep_uhf(params.cachelev, cachefiles);

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
    } else { /*** RHF/ROHF references ***/
        cachelist = cacheprep_rhf(params.cachelev, cachefiles);

        std::vector<int *> spaces;
        spaces.push_back(moinfo.occpi);
        spaces.push_back(moinfo.occ_sym);
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym);
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 2, spaces);
    }

    if (params.local) local_init();

    if (params.wfn == "CC2") {
        cc2_hbar_extra();
    } else {
        hbar_extra();
    }

    sort_lamps();                                /* should be removed sometime - provided by cclambda */
    if (params.wfn != "CC2") lambda_residuals(); /* don't do this for CC2 */

    preppert(ref_wfn->basisset());

    if (params.prop == "POLARIZABILITY") polar(ref_wfn);
    if (params.prop == "ROTATION") optrot(ref_wfn);
    if (params.prop == "ROA_TENSOR") roa(ref_wfn);

    if (params.local) local_done();

    dpd_close(0);

    if (params.ref == 2)
        cachedone_uhf(cachelist);
    else
        cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();

    exit_io();

    return PsiReturnType::Success;
}

void init_io() {
    int i;

    timer_on("ccresponse");

    for (i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);

    /* Clear out DIIS TOC Entries */
    psio_close(PSIF_CC_DIIS_AMP, 0);
    psio_close(PSIF_CC_DIIS_ERR, 0);

    psio_open(PSIF_CC_DIIS_AMP, 0);
    psio_open(PSIF_CC_DIIS_ERR, 0);
}

void title() {
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t*       ccresponse       *\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t**************************\n");
}

void exit_io() {
    int i;

    /* Close all dpd data files here */
    for (i = PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i, 1);
    for (i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i, 0); /* get rid of TMP files */
    for (i = PSIF_CC_TMP11 + 1; i <= PSIF_CC_MAX; i++) psio_close(i, 1);

    timer_off("ccresponse");
}

}
}  // namespace psi
