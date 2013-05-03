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

#include <string>
#include <libdist_matrix/dist_mat.h>
#include <libmints/dimension.h>
#include "detail.h"

namespace psi { namespace libmatrix {

    struct libmints_matrix_wrapper;

    struct libdist_globals {

        static std::string interface_name;

        static void initialize(int argc, char** argv) {
            // Nothing to do. Psi4 already prepares the environment.
        }
    };

#if defined(HAVE_MADNESS)

    struct libdist_matrix_wrapper {
        libdist_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) : matrix_()
        {
            std::vector<int> dims(2);
            dims[0] = 2; dims[1] = 2;
            process_grid pgrid(dims);

            matrix_.initialize(pgrid, m.sum(), n.sum(), 64, name);
        }

        void print() const {
            matrix_.print();
        }

        void fill(double val) {
            matrix_.fill(val);
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libdist_matrix_wrapper& rhs) {
            matrix_.add(rhs.matrix_);
        }

        // Is it not unreasonable to require all matrix types to know how to work with the serial matrix object.
        void add(const libmints_matrix_wrapper& rhs) {
            // TODO: implement
        }
    private:
        Distributed_Matrix matrix_;
    };
#else
    // We are not configured with Madness, derive from is_not_available to
    // cause compiler error when libdist_matrix_wrapper is instantiated in
    // the code.
    UNAVAILABLE_MATRIX(libdist_matrix_wrapper);
#endif

}} // end namespaces

