/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef _psi_src_bin_psimrcc_index_iterator_h_
#define _psi_src_bin_psimrcc_index_iterator_h_

/*! \file    index_iterator.h
    \ingroup (PSIMRCC)
    \brief   This class is used to iterate over n-tuples of MOs indices (p,q,r,..)
*/

#include <array>
#include <string>

#include "psimrcc_wfn.h"

namespace psi {
namespace psimrcc {

class CCIndex;

class CCIndexIterator {
   public:
    // Class Constructor and Destructor
    explicit CCIndexIterator(std::shared_ptr<PSIMRCCWfn> wfn, std::string str);
    explicit CCIndexIterator(std::shared_ptr<PSIMRCCWfn> wfn, std::string str, int select_irrep);
    explicit CCIndexIterator(std::shared_ptr<PSIMRCCWfn> wfn, CCIndex* index);
    explicit CCIndexIterator(std::shared_ptr<PSIMRCCWfn> wfn, CCIndex* index, int select_irrep);
    ~CCIndexIterator();

    // Class Public Methods
    bool first();
    void next();
    bool end() { return (absolute >= max_abs); }

    template <int N>
    short ind_abs() {
        return tuples[absolute][N];
    }

    int sym() { return symmetry; }
    size_t rel() { return relative; }
    size_t abs() { return absolute; }

   private:
    // Class private functions
    void startup(int min_sym, int max_sym);

    // Generica data
    int nirreps;

    // Index object
    CCIndex* ccindex;

    // Internal iterator
    size_t relative;  // Relative address of the current tuple
    size_t absolute;  // Absolute address of the current tuple
    size_t max_abs;   // Max absolute address of the current tuple
    size_t min_abs;   // Min absolute address of the current tuple
    int symmetry;     // Symmetry of the current tuple
    int current_block;

    // Properties of the tuples
    // These come from the Index object. Which should outlast the iterator, and which should
    // not be changing these values as long as this iterator is alive.
    int nelements;
    std::vector<std::vector<int>> element_irrep;
    const std::vector<std::array<short, 3>>& tuples;
    std::vector<size_t> block_last;
    std::vector<int> block_symmetry;
};

}  // namespace psimrcc
}  // namespace psi

#endif  // _psi_src_bin_psimrcc_index_iterator_h_
