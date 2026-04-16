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
// The interface to cuEST was contributed by NVIDIA under the following terms:
// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CUESTJK_H
#define CUESTJK_H

#include "jk.h"

#ifdef USING_cuEST
#include <cuest.h>
#endif

namespace psi {

class BasisSet;
class Options;

#ifdef USING_cuEST

class PSI_API cuESTJK : public JK {
   protected:
    Options& options_;

    std::shared_ptr<BasisSet> auxiliary_;

    double condition_;
    double pq_threshold_;

    bool initialized_;

    cuestAOBasis_t cuest_primary_basis_;
    cuestAOBasis_t cuest_auxiliary_basis_;
    cuestAOPairList_t cuest_pair_list_;
    cuestDFIntPlan_t cuest_df_plan_;

    cuestDFCoulombComputeParameters_t cuest_coulomb_compute_params_;
    cuestDFSymmetricExchangeComputeParameters_t cuest_exchange_compute_params_;

    cuestWorkspace_t* cuest_pair_list_ws_ptr_;
    cuestWorkspace_t* cuest_dfint_plan_ws_ptr_;

    std::string name() override { return "cuESTJK"; }
    size_t memory_estimate() override;
    bool C1() const override { return true; }
    void preiterations() override;
    void compute_JK() override;
    void postiterations() override;

    cuestAOBasis_t build_cuest_basis(std::shared_ptr<BasisSet> basis,
                                     std::vector<cuestAOShell_t>& shells_out,
                                     cuestWorkspace_t& persistent_ws);
    void allocate_workspace(cuestWorkspaceDescriptor_t& desc, cuestWorkspace_t& ws);
    void free_workspace(cuestWorkspace_t& ws);
    void destroy_cuest_objects();

   public:
    cuESTJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    ~cuESTJK() override;

    void set_condition(double condition) { condition_ = condition; }
    void print_header() const override;

    cuestDFIntPlan_t cuest_df_plan() { return cuest_df_plan_; }
};

#endif // USING_cuEST
} // namespace psi
#endif // CUESTJK_H
