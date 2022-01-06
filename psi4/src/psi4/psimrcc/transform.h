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

#ifndef _psi_src_bin_psimrcc_cctransform_h
#define _psi_src_bin_psimrcc_cctransform_h

#include <map>

#include "psimrcc_wfn.h"

namespace psi {

class IntegralTransform;

namespace psimrcc {

class CCIndex;

/**
        @author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCTransform {
   public:
    CCTransform(std::shared_ptr<PSIMRCCWfn> wfn);
    ~CCTransform();
    void print();
    // Presorting
    void presort_integrals();
    void read_oei_from_transqt() { read_oei_mo_integrals(); }
    void read_integrals_mrpt2(IntegralTransform* ints);
    int read_tei_mo_integrals_block(int first_irrep);
    void free_memory();
    void free_tei_mo_block(int first_irrep, int last_irrep);
    void transform_tei_integrals();
    double oei(int p, int q);
    double tei_block(int p, int q, int r, int s);
    double tei_mrpt2(int p, int q, int r, int s);

   private:
    std::vector<std::vector<double>> oei_mo;
    std::vector<std::vector<double>> oei_so;
    std::vector<std::vector<double>> tei_so;
    std::vector<std::vector<std::vector<double>>> tei_half_transformed;
    std::vector<std::vector<double>> tei_mo;
    CCIndex* oei_so_indexing;
    CCIndex* tei_so_indexing;
    CCIndex* tei_mo_indexing;
    std::shared_ptr<PSIMRCCWfn> wfn_;

    void read_oei_so_integrals();
    void read_oei_mo_integrals();
    void read_oei_mo_integrals_mrpt2();
    void read_tei_mo_integrals_mrpt2(IntegralTransform* ints);

    void transform_oei_so_integrals();

    void allocate_oei_so();
    void allocate_oei_mo();

    // Block
    int first_irrep_in_core;
    int last_irrep_in_core;
    int allocate_tei_mo_block(int first_irrep);
    std::map<size_t, double> integral_map;

    void presort_blocks(int first_irrep, int last_irrep);

    double fraction_of_memory_for_presorting;
};

}  // namespace psimrcc
}  // namespace psi

#endif  // _psi_src_bin_psimrcc_cctransform_h
