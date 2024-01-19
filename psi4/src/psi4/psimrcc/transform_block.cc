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

#include <cmath>
#include <algorithm>

#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsi4util/libpsi4util.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i, j) ((i > j) ? i : j)
#define MIN(i, j) ((i > j) ? j : i)
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))
#define four(i, j, k, l) INDEX(INDEX(i, j), INDEX(k, l))

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include "algebra_interface.h"
#include "blas.h"
#include "matrix.h"
#include "index.h"
#include "transform.h"

namespace psi {
namespace psimrcc {

/**
 * Read at least one block of the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them
 * in the packed array tei_mo
 */
int CCTransform::read_tei_mo_integrals_block(int first_irrep) {
    std::vector<size_t> pairpi = tei_mo_indexing->get_pairpi();
    int last_irrep = allocate_tei_mo_block(first_irrep);

    // Write integrals to disk
    for (int h = first_irrep; h < last_irrep; ++h) {
        char data_label[80];
        sprintf(data_label, "PRESORTED_TEI_IRREP_%d", h);
        _default_psio_lib_->read_entry(
            PSIF_PSIMRCC_INTEGRALS, data_label, (char *)&(tei_mo[h][0]),
            static_cast<size_t>(INDEX(pairpi[h] - 1, pairpi[h] - 1) + 1) * static_cast<size_t>(sizeof(double)));
    }
    return (last_irrep);
}

/**
 * Allocate as many blocks of the tei_mo array and exit(EXIT_FAILURE) if there is not enough space
 */
int CCTransform::allocate_tei_mo_block(int first_irrep) {
    if (first_irrep > wfn_->nirrep()) {
        outfile->Printf("\n    Transform: allocate_tei_mo_block() was called with first_irrep > nirreps !");

        exit(EXIT_FAILURE);
    }

    size_t available_transform_memory =
        static_cast<size_t>(static_cast<double>(wfn_->free_memory_) * fraction_of_memory_for_presorting);

    int last_irrep = first_irrep;

    tei_mo = std::vector<std::vector<double>>(wfn_->nirrep());

    // Find how many irreps we can store in 95% of the free memory
    std::vector<size_t> pairpi = tei_mo_indexing->get_pairpi();
    for (int h = first_irrep; h < wfn_->nirrep(); ++h) {
        size_t required_memory = (INDEX(pairpi[h] - 1, pairpi[h] - 1) + 1) * static_cast<size_t>(sizeof(double));
        if (required_memory != 0) {
            if (required_memory < available_transform_memory) {
                available_transform_memory -= required_memory;
                wfn_->free_memory_ -= required_memory;
                tei_mo[h] = std::vector<double>(INDEX(pairpi[h] - 1, pairpi[h] - 1) + 1, 0);
                last_irrep++;
            }
        } else {
            last_irrep++;
        }
    }
    outfile->Printf("\n    Integrals from irreps %d -> %d will be read in core", first_irrep, last_irrep - 1);
    if (first_irrep == last_irrep) {
        outfile->Printf("\n    CCTransform: allocate_tei_mo_block() has not enough memory!");

        exit(EXIT_FAILURE);
    }
    first_irrep_in_core = first_irrep;
    last_irrep_in_core = last_irrep;
    return (last_irrep);
}

void CCTransform::free_tei_mo_block(int first_irrep, int last_irrep) {
    for (auto h = first_irrep; h < last_irrep; h++) {
        wfn_->free_memory_ += sizeof(double) * tei_mo[h].size();
    }
    tei_mo.clear();
}

double CCTransform::tei_block(int p, int q, int r, int s) {
    // Get the (pq|rs) integral
    int irrep(tei_mo_indexing->get_tuple_irrep(MAX(p, q), MIN(p, q)));
    if ((first_irrep_in_core <= irrep) && (irrep < last_irrep_in_core))
        return (tei_mo[tei_mo_indexing->get_tuple_irrep(MAX(p, q), MIN(p, q))]
                      [INDEX(tei_mo_indexing->get_tuple_rel_index(MAX(p, q), MIN(p, q)),
                             tei_mo_indexing->get_tuple_rel_index(MAX(r, s), MIN(r, s)))]);
    else
        return (0.0);
}

}  // namespace psimrcc
}  // namespace psi
