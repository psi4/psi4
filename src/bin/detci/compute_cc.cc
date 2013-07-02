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

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterdset.h>
#include <libpsio/psio.h>
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

int cc_reqd_sblocks[CI_BLK_MAX];

/*
** compute_cc()
**
** This is the top-level function that controls the coupled-cluster
** computation
**
*/
void compute_cc(void)
{
  printf("compute_cc: Not yet available\n");
}

}} // namespace psi::detci

