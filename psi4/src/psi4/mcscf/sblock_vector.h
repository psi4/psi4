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

#ifndef _psi_src_lib_libmemtrix_sblock_vector_h_
#define _psi_src_lib_libmemtrix_sblock_vector_h_

#include <string>

#include "block_vector.h"

namespace psi {
namespace mcscf {

// Smart version of BlockVector
class SBlockVector {
   public:
    SBlockVector();
    SBlockVector(std::string label, int nirreps, int*& rows_size);
    SBlockVector(std::string label, int nirreps, vecint& rows_size);
    ~SBlockVector() {
        if (block_vector_)
            if (block_vector_->subtract_reference()) block_vector_ = nullptr;
    }

    // Manual allocation
    void allocate(std::string label, int nirreps, int*& rows_size);
    void allocate(std::string label, int nirreps, vecint& rows_size);

    void subtract_reference() {
        if (block_vector_) {
            if (block_vector_->subtract_reference()) block_vector_ = nullptr;
        }
    }

    // Copy constructor and assignment operator
    SBlockVector(const SBlockVector& src);
    SBlockVector& operator=(const SBlockVector& src);

    // Allow access to the implementation object
    const BlockVector* operator->() const { return block_vector_; }
    BlockVector* operator->() { return block_vector_; }

    // Access the implementation object
    BlockVector* getBlockVector() { return block_vector_; }

    // Checking functions
    bool is_allocated() { return (block_vector_); }
    void check(const char* cstr);

    void copy(SBlockVector& source);

   private:
    SBlockVector(BlockVector* block_vector);

    BlockVector* block_vector_;
};

}  // namespace mcscf
}  // namespace psi

#endif  // _psi_src_lib_libmemtrix_sblock_vector_h_
