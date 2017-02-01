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
 * @END LICENSE
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "psi4/liboptions/liboptions.h"

#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::initial_guess()
{
  using namespace psi;

  SBlockMatrix H_t("H_t",nirreps,sopi,sopi);
  SBlockVector eigenvectors("H_t",nirreps,sopi);

  transform(H,H_t,S_sqrt_inv);

  H_t.diagonalize(C_t,eigenvectors);

  C.multiply(false,false,S_sqrt_inv,C_t);

  epsilon = eigenvectors;

  guess_occupation();
}

}} /* End Namespaces */
