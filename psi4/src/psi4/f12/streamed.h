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

#pragma once

#include "mp2.h"

namespace psi {
namespace f12 {

struct PairDFScratch {
    einsums::Tensor<double, 2> aux_left;
    einsums::Tensor<double, 2> aux_right;
    einsums::Tensor<double, 2> aux_other;

    PairDFScratch(int workspace_id, int naux, int max_ri)
        : aux_left{"DF aux-left " + std::to_string(workspace_id), naux, max_ri},
          aux_right{"DF aux-right " + std::to_string(workspace_id), naux, max_ri},
          aux_other{"DF aux-other " + std::to_string(workspace_id), naux, max_ri} {}
};

class StreamedMP2F12 : public MP2F12 {
   public:
    StreamedMP2F12(SharedWavefunction ref_wfn, Options& options);
    double compute_energy() override;

   private:
    struct PairResult {
        double energy_s = 0.0;
        double energy_t = 0.0;
    };
    struct PairWorkspace {
        einsums::Tensor<double, 2> Gij;
        einsums::Tensor<double, 2> Gji;
        einsums::Tensor<double, 2> Fij;
        einsums::Tensor<double, 2> Fji;
        einsums::Tensor<double, 2> F2ij;
        einsums::Tensor<double, 2> F2ji;
        einsums::Tensor<double, 2> FGij;
        einsums::Tensor<double, 2> Ufij;
        einsums::Tensor<double, 2> Cij;
        einsums::Tensor<double, 2> Cji;
        einsums::Tensor<double, 2> Dij;
        einsums::Tensor<double, 2> product_full;
        einsums::Tensor<double, 2> product_occ_all;
        einsums::Tensor<double, 2> product_cabs_occ;
        einsums::Tensor<double, 2> product_vir_obs;
        PairDFScratch df_scratch;

        PairWorkspace(int workspace_id, int nobs, int nri, int nocc, int nact,
                      int nvir, int ncabs, int naux)
            : Gij{"ws Gij " + std::to_string(workspace_id), nobs, nri},
              Gji{"ws Gji " + std::to_string(workspace_id), nobs, nri},
              Fij{"ws Fij " + std::to_string(workspace_id), nri, nri},
              Fji{"ws Fji " + std::to_string(workspace_id), nri, nri},
              F2ij{"ws F2ij " + std::to_string(workspace_id), nact, nri},
              F2ji{"ws F2ji " + std::to_string(workspace_id), nact, nri},
              FGij{"ws FGij " + std::to_string(workspace_id), nact, nact},
              Ufij{"ws Ufij " + std::to_string(workspace_id), nact, nact},
              Cij{"ws Cij " + std::to_string(workspace_id), nvir, nvir},
              Cji{"ws Cji " + std::to_string(workspace_id), nvir, nvir},
              Dij{"ws Dij " + std::to_string(workspace_id), nvir, nvir},
              product_full{"ws product full " + std::to_string(workspace_id), nri, nri},
              product_occ_all{"ws product occ-all " + std::to_string(workspace_id), nocc, nri},
              product_cabs_occ{"ws product cabs-occ " + std::to_string(workspace_id), ncabs, nocc},
              product_vir_obs{"ws product vir-obs " + std::to_string(workspace_id), nvir, nobs},
              df_scratch{workspace_id, naux, nri} {}
    };

    int auxiliary_block_size_;

    int plan_pair_workers() const;

    void three_index_direct_j_ao_computer(
        const einsums::Tensor<double, 1>& j_weight,
        einsums::Tensor<double, 2>* j_ao,
        std::shared_ptr<BasisSet> bs1,
        std::shared_ptr<BasisSet> bs2);

    void three_index_mo_aux_blocked_pack(
        const std::vector<std::string>& int_types,
        const std::vector<einsums::Tensor<double, 3>*>& BPQ_outputs,
        std::shared_ptr<BasisSet> bs1,
        std::shared_ptr<BasisSet> bs2,
        const einsums::Tensor<double, 2>& C1,
        const einsums::Tensor<double, 2>& C2);

    void form_metric_inverse(einsums::Tensor<double, 2>* metric_inverse);

    void form_oper_ints_pack(
        const std::vector<std::string>& int_types,
        const std::vector<einsums::Tensor<double, 3>*>& DF_ERI_outputs);

    void form_oper_ints_pack(
        const std::vector<std::string>& int_types,
        const std::vector<einsums::Tensor<double, 2>*>& DF_ERI_outputs);

    void form_df_pair_ints_into(
        const std::string& int_type, einsums::Tensor<double, 3>* metric,
        einsums::Tensor<double, 3>* oper, einsums::Tensor<double, 2>* aux_oper,
        int pidx, int ridx, int rows, int cols, const std::vector<char>& order,
        einsums::Tensor<double, 2>* result, PairDFScratch* scratch);

    void form_df_fock_streamed(einsums::Tensor<double, 2>* f,
                              einsums::Tensor<double, 2>* k,
                              einsums::Tensor<double, 2>* fk,
                              const einsums::Tensor<double, 2>* metric_inverse,
                              einsums::Tensor<double, 3>* raw_g_occ_ri);
};

}  // namespace f12
}  // namespace psi
