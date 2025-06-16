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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
/*
**  CCHBAR: Program to calculate the elements of the CCSD HBAR matrix.
*/

#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"

#include <cstdio>
#include <cstdlib>
#include <string>

namespace psi {
namespace cchbar {

void init_io();
void title();
void get_moinfo(std::shared_ptr<Wavefunction> ref_wfn, Options &);
void get_params(Options &);
void exit_io();
void F_build();
void Wmbej_build();
void Wmnie_build();
void Wmbij_build();
void Wabij_build();
void Wamef_build();
void Wabei_build();
void cc2_Zmbej_build();
void cc2_Wmbej_build();
void cc2_Wmbij_build();
void cc2_Wabei_build();
void purge();
void cleanup();
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void sort_amps();
void tau_build();
void taut_build();
void status(const char *, std::string);
void cc3_HET1();
void Fai_build();
void reference();
void norm_HET1();

using namespace psi;

PsiReturnType cchbar(std::shared_ptr<Wavefunction> ref_wfn, Options &options) {
    int **cachelist, *cachefiles;

    init_io();
    title();
    get_moinfo(ref_wfn, options);
    get_params(options);

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if (params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

        cachelist = cacheprep_rhf(params.cachelev, cachefiles);

        std::vector<int *> spaces;
        spaces.push_back(moinfo.occpi);
        spaces.push_back(moinfo.occ_sym);
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym);
        dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, nullptr, 2, spaces);
    } else if (params.ref == 2) { /** UHF **/

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
    }

    sort_amps();
    tau_build();
    taut_build();

    if (params.Tamplitude) {
        reference();
        Fai_build();
        Wabij_build();
    }

    F_build();
    if (params.print & 2) status("F elements", "outfile");

    Wamef_build();
    if (params.print & 2) status("Wamef elements", "outfile");
    Wmnie_build();
    if (params.print & 2) status("Wmnie elements", "outfile");

    if (params.wfn == "CC2" || params.wfn == "EOM_CC2") {
        cc2_Wmbej_build();
        if (params.print & 2) status("Wmbej elements", "outfile");
        cc2_Zmbej_build();
        if (params.print & 2) status("Zmbej elements", "outfile");
        cc2_Wmbij_build();
        if (params.print & 2) status("Wmbij elements", "outfile");
        cc2_Wabei_build();
        if (params.print & 2) status("Wabei elements", "outfile");
    } else {
        Wabei_build();
        if (params.print & 2) status("Wabei elements", "outfile");
        Wmbej_build();
        if (params.print & 2) status("Wmbej elements", "outfile");
        Wmbij_build();
        if (params.print & 2) status("Wmbij elements", "outfile");

        if (params.wfn == "CC3" || params.wfn == "EOM_CC3") {
            /* switch to ROHF to generate all spin cases of He^T1 elements */
            if ((params.dertype == 3 || params.dertype == 1) && params.ref == 0) {
                params.ref = 1;
                cc3_HET1(); /* compute remaining Wmbej [H,eT1] */
                norm_HET1();
                params.ref = 0;
            } else {
                cc3_HET1(); /* compute remaining Wmbej [H,eT1] */
                norm_HET1();
            }
        }
    }

    if (params.ref == 1) purge(); /** ROHF only **/
    dpd_close(0);

    if (params.ref == 2)
        cachedone_uhf(cachelist);
    else
        cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();
    exit_io();
    return Success;
}

void init_io() {
    timer_on("cchbar");

    /* Open all dpd data files */
    for (int i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);
}

void title() {
    outfile->Printf("\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t*         CCHBAR         *\n");
    outfile->Printf("\t\t\t*                        *\n");
    outfile->Printf("\t\t\t**************************\n");
    outfile->Printf("\n");
}

void exit_io() {
    /* Close all dpd data files here */
    for (int i = PSIF_CC_MIN; i < PSIF_CC_TMP; ++i) psio_close(i, 1);
    for (int i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; ++i) psio_close(i, 0); /* get rid of TMP files */
    for (int i = PSIF_CC_TMP11 + 1; i <= PSIF_CC_MAX; ++i) psio_close(i, 1);

    timer_off("cchbar");
}

}  // namespace cchbar
}  // namespace psi
