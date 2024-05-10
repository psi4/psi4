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

#ifndef _psi_src_bin_psimrcc_ccsort_h
#define _psi_src_bin_psimrcc_ccsort_h

/**
 *  @file sort.h
 *  @ingroup (PSIMRCC)
 */

#include <iostream>
#include <map>
#include <vector>

#include "psimrcc_wfn.h"
#include "transform.h"

namespace psi {

class IntegralTransform;

namespace psimrcc {

class CCMatrix;

#include "matrix_types.h"

#ifndef INDEX
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))
#endif
#define four(i, j, k, l) INDEX(INDEX(i, j), INDEX(k, l))

enum SortAlgorithm { out_of_core_sort, mrpt2_sort };

/**
 *  @class CCSort
 *  @brief Grabs the MO integrals from CCTransform and fills the CCMatrix objects in CCBLAS
 */
class CCSort {
   public:
    CCSort(std::shared_ptr<PSIMRCCWfn> wfn, SortAlgorithm algorithm);
    ~CCSort();

   private:
    void init();

    //  // In-core algorithm
    //  void   build_integrals_in_core();
    //  void   frozen_core_energy_in_core();
    //  void   sort_integrals_in_core();
    //  void   form_two_electron_integrals_in_core(MatrixMap::iterator& iter);
    //  void   form_fock_in_core(MatrixMap::iterator& iter);
    //  double add_fock_two_in_core(int p, int q, int k, bool exchange);

    // Out-of-core algorithm
    void build_integrals_out_of_core();
    void frozen_core_energy_out_of_core();
    void sort_integrals_out_of_core(int first_irrep, int last_irrep, MatrixBlks& to_be_processed);
    void form_fock_one_out_of_core(MatrixBlks& to_be_processed);
    void form_fock_out_of_core(CCMatrix* Matrix, int h);
    void form_two_electron_integrals_out_of_core(CCMatrix* Matrix, int h);
    double add_fock_two_out_of_core(int p, int q, int k, bool exchange);
    void setup_out_of_core_list(MatMapIt& mat_it, int& mat_irrep, MatMapIt& mat_end, MatrixBlks& to_be_processed);
    void dump_integrals_to_disk(MatrixBlks& to_be_processed);

    // MRPT2 algorithm
    void build_integrals_mrpt2(IntegralTransform* ints);
    void frozen_core_energy_mrpt2();
    void allocate_and_sort_integrals_mrpt2();
    void allocate_amplitudes_mrpt2();
    void form_two_electron_integrals_mrpt2(MatrixMap::iterator& iter);
    void form_fock_mrpt2(MatrixMap::iterator& iter);
    double add_fock_two_mrpt2(int p, int q, int k, bool exchange);

    // Data
    double fraction_of_memory_for_sorting;
    int nfzc;
    double efzc;
    std::vector<int> frozen_core;
    // Transformer
    std::shared_ptr<CCTransform> trans;
    std::shared_ptr<PSIMRCCWfn> wfn_;
};

}  // namespace psimrcc
}  // namespace psi

#endif  // _psi_src_bin_psimrcc_ccsort_h
