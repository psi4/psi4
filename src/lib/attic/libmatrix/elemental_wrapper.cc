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

#include "libmints_wrapper.h"
#include "elemental_wrapper.h"

namespace psi { namespace libmatrix {
#if defined(HAVE_ELEMENTAL)
    std::string libelemental_globals::interface_name = "libelemental";
    elem::mpi::Comm libelemental_globals::mpi_comm;
    int libelemental_globals::rank = -1;
#endif
}} // end namespaces

