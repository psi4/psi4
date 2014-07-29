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

#include <iostream>
#include <cstdio>

#include <libchkpt/chkpt.hpp>
#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

double SCF::compute_energy()
{
  psi::fprintf(outfile,"\n\n  Running an SCF calculation");

  startup();

  // Read the one electron integrals
  read_so_oei();

  // Read the two electron integrals
  // and construct the PK and K matrices
  read_so_tei();

  // Construct the S^-1/2 Matrix
  construct_S_inverse_sqrt();

  // Guess C
  initial_guess();

  // Iterate the SCF equations
  iterate_scf_equations();

  // Check the orthonormality of the MOs
  check_orthonormality();

  // Canonicalize MOs
  canonicalize_MO();

  // Print eigenvectors and MOs
  print_eigenvectors_and_MO();

  // Canonicalize MOs
  save_info();

  return(0.0);
}

}} /* End Namespaces */
