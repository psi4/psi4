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

/**
 *  @file ccsort.cpp
 *  @ingroup (PSIMRCC)
*/
#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/memory_manager.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include "blas.h"
#include "sort.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
extern MOInfo *moinfo;
extern MemoryManager* memory_manager;

using namespace std;

CCSort::CCSort(SharedWavefunction ref_wfn, SortAlgorithm algorithm):
    fraction_of_memory_for_sorting(0.5),nfzc(0),efzc(0.0),frozen_core(0)
{
    init();

    IntegralTransform *ints;
    // Use libtrans to generate MO basis integrals in Pitzer order
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    if(algorithm == mrpt2_sort){
        // For MRPT2 calculations we need integrals of the form (OO|VV) and (OA|OA), where
        // O represents occ+act, V represents act+vir and A is all orbitals.
        std::shared_ptr<MOSpace> aocc;
        std::shared_ptr<MOSpace> avir;
        int nirrep = ref_wfn->nirrep();
        std::vector<int> aocc_orbs;
        std::vector<int> avir_orbs;
        std::vector<int> actv = moinfo->get_actv();
        std::vector<int> mopi = moinfo->get_mopi();
        std::vector<int> occ  = moinfo->get_occ();
        int offset = 0;
        for(int h = 0; h < nirrep; ++h){
            for(int i = 0; i < occ[h] + actv[h]; ++i)
                aocc_orbs.push_back(i + offset);
            for(int a = occ[h]; a < mopi[h]; ++a)
                avir_orbs.push_back(a + offset);
            offset += mopi[h];
        }
        aocc = std::shared_ptr<MOSpace>(new MOSpace('M', aocc_orbs, aocc_orbs));
        avir = std::shared_ptr<MOSpace>(new MOSpace('E', avir_orbs, avir_orbs));
        spaces.push_back(aocc);
        spaces.push_back(avir);
        ints = new IntegralTransform(ref_wfn, spaces,
                                     IntegralTransform::Restricted,
                                     IntegralTransform::DPDOnly,
                                     IntegralTransform::PitzerOrder,
                                     IntegralTransform::None);
        ints->set_keep_dpd_so_ints(true);
        // Only transform the subclasses needed for MRPT2
        ints->transform_tei(aocc, MOSpace::all, aocc, MOSpace::all);
        ints->set_keep_dpd_so_ints(false);
        ints->transform_tei(aocc, aocc, avir, avir);
        build_integrals_mrpt2(ints);
    }else{
        ints = new IntegralTransform(ref_wfn, spaces,
                                     IntegralTransform::Restricted,
                                     IntegralTransform::IWLOnly,
                                     IntegralTransform::PitzerOrder,
                                     IntegralTransform::None);
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        // Presort the integrals in the CCTransform class
        trans->presort_integrals();
        // Out-of-core algorithm: the transformed integrals and the CCMatrix(s) don't fit in
        // core or they are requested to be out-of-core
        build_integrals_out_of_core();
    }

    moinfo->set_fzcore_energy(efzc);
    outfile->Printf("\n\n    Frozen-core energy                     = %20.9f",efzc);
    delete ints;

}

CCSort::~CCSort()
{
    cleanup();
}

/**
 * Initialize the CCSort class
 */
void CCSort::init()
{
    // Find the frozen core orbitals in Pitzer ordering
    nfzc        = moinfo->get_nfocc();
    intvec focc = moinfo->get_focc();
    intvec mopi = moinfo->get_mopi();
    allocate1(int,frozen_core,nfzc);
    int count1  = 0;
    int count2  = 0;
    for(int h = 0; h < moinfo->get_nirreps(); ++h){
        for(int i = 0; i < focc[h]; ++i)
            frozen_core[count1++] = count2 + i;
        count2 += mopi[h];
    }
}

/**
 * Clean the CCSort class
 */
void CCSort::cleanup()
{
    release1(frozen_core);
}

}} /* End Namespaces */
