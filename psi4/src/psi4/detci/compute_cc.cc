/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

/*! \file
**  \ingroup DETCI
**  \brief Arbitrary-order coupled-cluster code
**
** C. David Sherrill
** Center for Computational Molecular Science and Technology
** Georgia Institute of Technology
** March 2005
**
** Note: I think I need onel ints as g for formation of sigma
** in non-FCI cases, but make sure any CC parts don't try to get
** h and actually get g instead...
*/

#include <cstdio>
#include <cstdlib>
#include "psi4/psi4-dec.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi {
namespace detci {

// int cc_reqd_sblocks[CI_BLK_MAX];

/*
** compute_cc()
**
** This is the top-level function that controls the coupled-cluster
** computation
**
*/
void CIWavefunction::compute_cc() { outfile->Printf("compute_cc: Not yet available\n"); }

}  // namespace detci
}  // namespace psi
