/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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
        : storage_(WorldComm->nthread(), value)
    { }

    void initialize(const T& value) {
        storage_.clear();
        for (int i=0; i<WorldComm->nthread(); ++i) {
            storage_.push_back(value);
        }
    }

    T& operator*() {
        return storage_[WorldComm->thread_id(pthread_self())];
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

