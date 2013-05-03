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

/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/

namespace psi { namespace stable {

struct Params {
  int print_lvl;         /* Output level control */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;
  int ref;
  int follow_instab;     /* follow a UHF->UHF instability of same symm? */
  int num_evecs_print;   /* print n lowest eigenvectors of MO hessian */
  int rotation_method;   /* 0 = by angles, 1 = by antisymmetric matrix */
  double scale;          /* scale factor for orbital rotation step */
};

}} // namespace psi::stable
