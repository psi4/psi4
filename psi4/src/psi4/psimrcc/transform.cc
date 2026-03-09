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

#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i, j) ((i > j) ? i : j)
#define MIN(i, j) ((i > j) ? j : i)
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))
#define four(i, j, k, l) INDEX(INDEX(i, j), INDEX(k, l))

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "matrix.h"
#include "transform.h"

namespace psi {
namespace psimrcc {

CCTransform::CCTransform(std::shared_ptr<PSIMRCCWfn> wfn) : fraction_of_memory_for_presorting(0.75), wfn_(wfn) {
    wfn_->blas()->add_index("[s>=s]");
    wfn_->blas()->add_index("[n>=n]");
    wfn_->blas()->add_index("[s]");
    tei_mo_indexing = wfn_->blas()->get_index("[n>=n]");
    tei_so_indexing = wfn_->blas()->get_index("[s>=s]");
    oei_so_indexing = wfn_->blas()->get_index("[s]");
    first_irrep_in_core = 0;
    last_irrep_in_core = 0;
}

CCTransform::~CCTransform() { free_memory(); }

/**
 * Read the one electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in oei_mo
 */
void CCTransform::read_oei_mo_integrals() {
    // Read all the (frozen + non-frozen) OEI in Pitzer order
    allocate_oei_mo();

    int nmo = wfn_->moinfo()->get_nmo();

    std::vector<double> H(INDEX(nmo - 1, nmo - 1) + 1, 0);

    iwl_rdone(PSIF_OEI, const_cast<char*>(PSIF_MO_FZC), H.data(), nmo * (nmo + 1) / 2, 0, 0, "outfile");
    //   else
    //     iwl_rdone(PSIF_OEI,PSIF_MO_FZC,H,norbs*(norbs+1)/2,0,1,outfile); //TODO fix it!

    for (int i = 0; i < nmo; i++)
        for (int j = 0; j < nmo; j++) oei_mo[i][j] = H[INDEX(i, j)];
}

/**
 * Free all the memory allocated by CCTransform
 */
void CCTransform::free_memory() { integral_map.clear(); }

/**
 * Allocate the oei_mo array
 */
void CCTransform::allocate_oei_mo() {
    if (oei_mo.size() == 0) {
        int nmo = wfn_->nmo();
        oei_mo = std::vector<std::vector<double>>(nmo, std::vector<double>(nmo, 0));
    }
}

/**
 * Allocate the oei_so array
 */
void CCTransform::allocate_oei_so() {
    if (oei_mo.size() == 0) {
        int nso = wfn_->nso();
        oei_so = std::vector<std::vector<double>>(nso, std::vector<double>(nso, 0));
    }
}

double CCTransform::oei(int p, int q) { return (oei_mo[p][q]); }

void CCTransform::transform_oei_so_integrals() {
    outfile->Printf("\n  CCTransform: transforming one-electron integrals");

    allocate_oei_mo();

    const int nso = wfn_->nso();
    const int nmo = wfn_->nmo();

    std::vector<std::vector<double>> A(nso, std::vector<double>(nmo, 0));
    const auto C = wfn_->moinfo()->get_scf_mos();

    // A(q,i) = H(q,p) * C(p,i)
    /*#ifdef CCTRANSFORM_USE_BLAS
      C_DGEMM_12(
    #else*/
    for (int q = 0; q < nso; q++)
        for (int j = 0; j < nmo; j++) {
            A[q][j] = 0.0;
            for (int p = 0; p < nso; p++) A[q][j] += oei_so[q][p] * C[p][j];
        }
    for (int i = 0; i < nmo; i++)
        for (int j = 0; j < nmo; j++) {
            oei_mo[i][j] = 0.0;
            for (int q = 0; q < nso; q++) oei_mo[i][j] += C[q][i] * A[q][j];
        }
    // #endif
}

}  // namespace psimrcc
}  // namespace psi
