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

#if !defined(psi_src_lib_libparallel_threaded_storage)
#define psi_src_lib_libparallel_threaded_storage

#include <vector>

namespace psi {

template<typename T>
class threaded_storage
{
    std::vector<T> storage_;

public:

    threaded_storage(const T& value = T())
        : storage_(Process::environment.get_n_threads(), value)
    { }

    void initialize(const T& value) {
        storage_.clear();
        for (int i=0; i<Process::environment.get_n_threads(); ++i) {
            storage_.push_back(value);
        }
    }

    T& operator*() {
        return storage_[0];
    }

    const T& operator[](size_t ind) const {
        return storage_[ind];
    }
    T& operator[](size_t ind) {
        return storage_[ind];
    }
};

}

#endif