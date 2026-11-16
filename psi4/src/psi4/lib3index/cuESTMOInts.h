/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef THREE_INDEX_CUESTMOINTS_H
#define THREE_INDEX_CUESTMOINTS_H

#ifdef USING_cuEST

#include <cstddef>
#include <cstdint>
#include <memory>

#include <cuest.h>

#include "psi4/libmints/typedefs.h"

namespace psi {

class BasisSet;
class Options;

/**
 * GPU density-fitted MO integrals via cuEST.
 *
 * Owns the cuEST AO pair list and DF integral plan built from the given
 * (primary, auxiliary) pair. The auxiliary basis decides what the integrals
 * are, so pass the basis for the method (e.g. DF_BASIS_MP2), not the SCF's.
 *
 * The cuestAOBasis_t handles themselves are owned by the BasisSet objects,
 * which this class keeps alive.
 */
class cuESTMOInts {
   protected:
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> auxiliary_;

    double pq_threshold_;
    int32_t integral_direct_;

    cuestAOBasis_t cuest_primary_basis_ = nullptr;
    cuestAOBasis_t cuest_auxiliary_basis_ = nullptr;
    cuestAOPairList_t cuest_pair_list_ = nullptr;
    cuestDFIntPlan_t cuest_df_plan_ = nullptr;

    cuestWorkspace_t* cuest_pair_list_ws_ptr_ = nullptr;
    cuestWorkspace_t* cuest_dfint_plan_ws_ptr_ = nullptr;

    void common_init();

   public:
    cuESTMOInts(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options,
                bool integral_direct = false);
    ~cuESTMOInts();

    cuESTMOInts(const cuESTMOInts&) = delete;
    cuESTMOInts& operator=(const cuESTMOInts&) = delete;

    /**
     * Fitted 3-index MO integrals B(Q|lr) = sum_mn B(Q|mn) C_ml C_nr.
     *
     * Cleft / Cright: Psi4 AO-basis coefficient matrices (nbf x nleft, nbf x nright, C1).
     * d_out:          caller-owned DEVICE buffer of naux * nleft * nright doubles,
     *                 written row-major as [Q][l][r].
     */
    void compute_3c_ints(SharedMatrix Cleft, SharedMatrix Cright, double* d_out);
};

}  // namespace psi

#endif  // USING_cuEST
#endif  // THREE_INDEX_CUESTMOINTS_H
