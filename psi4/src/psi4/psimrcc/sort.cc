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

/**
 *  @file ccsort.cpp
 *  @ingroup (PSIMRCC)
 */
#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include "blas.h"
#include "sort.h"
#include "transform.h"

namespace psi {
namespace psimrcc {

CCSort::CCSort(std::shared_ptr<PSIMRCCWfn> wfn, SortAlgorithm algorithm)
    : fraction_of_memory_for_sorting(0.5), nfzc(0), efzc(0.0), wfn_(wfn) {
    init();

    trans = std::make_shared<CCTransform>(wfn_);

    IntegralTransform* ints;
    // Use libtrans to generate MO basis integrals in Pitzer order
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    if (algorithm == mrpt2_sort) {
        // For MRPT2 calculations we need integrals of the form (OO|VV) and (OA|OA), where
        // O represents occ+act, V represents act+vir and A is all orbitals.
        std::shared_ptr<MOSpace> aocc;
        std::shared_ptr<MOSpace> avir;
        int nirrep = wfn_->nirrep();
        std::vector<int> aocc_orbs;
        std::vector<int> avir_orbs;
        std::vector<int> actv = wfn_->moinfo()->get_actv();
        std::vector<int> mopi = wfn_->moinfo()->get_mopi();
        std::vector<int> occ = wfn_->moinfo()->get_occ();
        int offset = 0;
        for (int h = 0; h < nirrep; ++h) {
            for (int i = 0; i < occ[h] + actv[h]; ++i) aocc_orbs.push_back(i + offset);
            for (int a = occ[h]; a < mopi[h]; ++a) avir_orbs.push_back(a + offset);
            offset += mopi[h];
        }
        aocc = std::make_shared<MOSpace>('M', aocc_orbs, aocc_orbs);
        avir = std::make_shared<MOSpace>('E', avir_orbs, avir_orbs);
        spaces.push_back(aocc);
        spaces.push_back(avir);
        ints = new IntegralTransform(wfn, spaces, IntegralTransform::TransformationType::Restricted,
                                     IntegralTransform::OutputType::DPDOnly, IntegralTransform::MOOrdering::PitzerOrder,
                                     IntegralTransform::FrozenOrbitals::None);
        ints->set_keep_dpd_so_ints(true);
        // Only transform the subclasses needed for MRPT2
        ints->transform_tei(aocc, MOSpace::all, aocc, MOSpace::all);
        ints->set_keep_dpd_so_ints(false);
        ints->transform_tei(aocc, aocc, avir, avir);
        build_integrals_mrpt2(ints);
    } else {
        ints = new IntegralTransform(wfn, spaces, IntegralTransform::TransformationType::Restricted,
                                     IntegralTransform::OutputType::IWLOnly, IntegralTransform::MOOrdering::PitzerOrder,
                                     IntegralTransform::FrozenOrbitals::None);
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        // Presort the integrals in the CCTransform class
        trans->presort_integrals();
        // Out-of-core algorithm: the transformed integrals and the CCMatrix(s) don't fit in
        // core or they are requested to be out-of-core
        build_integrals_out_of_core();
    }

    wfn_->moinfo()->set_fzcore_energy(efzc);
    outfile->Printf("\n\n    Frozen-core energy                     = %20.9f", efzc);
    delete ints;
}

CCSort::~CCSort() {}

/**
 * Initialize the CCSort class
 */
void CCSort::init() {
    // Find the frozen core orbitals in Pitzer ordering
    nfzc = wfn_->moinfo()->get_nfocc();
    intvec focc = wfn_->moinfo()->get_focc();
    intvec mopi = wfn_->moinfo()->get_mopi();
    int count1 = 0;
    int count2 = 0;
    frozen_core = std::vector<int>(nfzc, 0);
    for (int h = 0; h < wfn_->nirrep(); ++h) {
        for (int i = 0; i < focc[h]; ++i) frozen_core[count1++] = count2 + i;
        count2 += mopi[h];
    }
}

}  // namespace psimrcc
}  // namespace psi
