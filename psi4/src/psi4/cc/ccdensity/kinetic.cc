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
    \brief Computes the kinetic energy and the virial ratio for CC wave functions.
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

#define INDEX(i, j) ((i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i))

void kinetic(std::shared_ptr<Wavefunction> wfn) {
    int nmo, noei, stat, i, I, h, j, nclsd;
    int *order, *doccpi;
    double junk, tcorr, vcorr, tref, vref, ttot, vtot;
    double **T, **scf_pitzer, **scf_qt, **X;

    /* RHF/ROHF only for now */
    if (params.ref == 2) return;

    nmo = moinfo.nmo;
    noei = nmo * (nmo + 1) / 2;

    /*** Get the Pitzer -> QT reordering array ***/
    order = init_int_array(nmo);

    /* doccpi array must include frozen orbitals for reorder_qt() */
    doccpi = init_int_array(moinfo.nirreps);
    for (h = 0; h < moinfo.nirreps; h++) doccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];

    reorder_qt(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, order, moinfo.orbspi, moinfo.nirreps);

    /*** Reorder the SCF eigenvectors to QT ordering */
    scf_pitzer = wfn->Ca()->to_block_matrix();

    scf_qt = block_matrix(nmo, nmo);
    for (i = 0; i < nmo; i++) {
        I = order[i]; /* Pitzer --> QT */
        for (j = 0; j < nmo; j++) scf_qt[j][I] = scf_pitzer[j][i];
    }

    /*** Transform the kinetic energy integrals to the MO basis ***/
    T = wfn->mintshelper()->so_kinetic()->to_block_matrix();
    X = block_matrix(nmo, nmo);

    C_DGEMM('t', 'n', nmo, nmo, nmo, 1, &(scf_qt[0][0]), nmo, &(T[0][0]), nmo, 0, &(X[0][0]), nmo);
    C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(scf_qt[0][0]), nmo, 0, &(T[0][0]), nmo);

    /*** Contract the correlated kinetic energy ***/

    tcorr = 0.0;
    for (i = 0; i < nmo; i++)
        for (j = 0; j < nmo; j++) tcorr += T[i][j] * moinfo.opdm[i][j];

    /*** Compute the SCF kinetic energy ***/

    tref = 0.0;
    nclsd = moinfo.nfzc + moinfo.nclsd;
    for (i = 0; i < nclsd; i++) tref += T[i][i] * 2;
    for (i = nclsd; i < nclsd + moinfo.nopen; i++) tref += T[i][i];

    /*** Compute the virial ratios ***/
    ttot = tcorr + tref;
    vtot = moinfo.eref + moinfo.ecc - ttot;
    vref = moinfo.eref - tref;
    vcorr = moinfo.ecc - tcorr;

    outfile->Printf("\n\tVirial Theorem Data:\n");
    outfile->Printf("\t--------------------\n");
    outfile->Printf("\tKinetic energy (ref)   = %20.15f\n", tref);
    outfile->Printf("\tKinetic energy (corr)  = %20.15f\n", tcorr);
    outfile->Printf("\tKinetic energy (total) = %20.15f\n", ttot);

    outfile->Printf("\t-V/T (ref)             = %20.15f\n", -vref / tref);
    outfile->Printf("\t-V/T (corr)            = %20.15f\n", -vcorr / tcorr);
    outfile->Printf("\t-V/T (total)           = %20.15f\n", -vtot / ttot);

    wfn.set_scalar_variable("CC DELTA KINETIC ENERGY", tcorr);
    wfn.set_scalar_variable("CC DELTA POTENTIAL ENERGY", vcorr);
    wfn.set_scalar_variable("CC DELTA VIRIAL RATIO", -vtot / ttot + vref / tref);

    /*** Release memory ***/
    free_block(X);
    free_block(T);
    free_block(scf_qt);
    free_block(scf_pitzer);
    free(doccpi);
    free(order);
}
}  // namespace ccdensity
}  // namespace psi
