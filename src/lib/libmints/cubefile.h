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

#ifndef _psi_src_lib_libmints_cubefile_h_
#define _psi_src_lib_libmints_cubefile_h_

#include <string>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

    class Wavefunction;
    class Matrix;

    class CubeFile{
        protected:
            /// The wavefunction object
            boost::shared_ptr<Wavefunction> wfn_;
            /// The AO basis (no spherical harmonics) density matrix
            boost::shared_ptr<Matrix> Dao_;
            /// The number of atomic orbitals
            int nao_;
            /// Padding (bohr) between the extrema of the molecule and the box walls
            double buffer_region_;
            /// The name of the cube file
            std::string filename_;
            /// Number of grid points in each dimension
            int nptsx_;
            int nptsy_;
            int nptsz_;
            /// The bounding box
            double minx_;
            double miny_;
            double minz_;
            double maxx_;
            double maxy_;
            double maxz_;
            /// The grid spacing in each dimension
            double incx_;
            double incy_;
            double incz_;
            /// Find the size of the grid needed to encapsulate the molecule
            void init_dimensions();
        public:
            CubeFile();
            void process_density();
            void set_npts(int x, int y, int z) { nptsx_ = x; nptsy_ = y; nptsz_ = z; }
            void set_buffer(double val) { buffer_region_ = val; }
            void set_filename(const std::string &str) { filename_ = str; }
    };

} // end namespace
#endif
