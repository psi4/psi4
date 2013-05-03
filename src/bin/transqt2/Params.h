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
    \ingroup TRANSQT2
    \brief Enter brief description of file here 
*/

#include <string>

namespace psi {
  namespace transqt2 {

struct Params {
  std::string wfn;
  int ref;
  int cachelev;
  int dertype;
  int reset;          /* cmdline argument; if true, all CC-related
                         files are deleted at the beginning of the
                         run */
  int print_lvl;      /* Output level control */
  int print_tei;      /* Boolean for printing two-electron integrals */
  double tolerance;   /* Cutoff value for integrals in IWL Buffers */
  long int memory;    /* Memory available (in bytes) */
  int semicanonical;  /* Boolean for semicanonical orbitals */
  int delete_tei;     /* Boolean for the TEI integral file */
  int backtr;         /* Boolean for back-transforms (not yet implemented) */
};

  } // namespace transqt2
} // namespace psi
