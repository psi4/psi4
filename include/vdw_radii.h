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

/*------------------------------------------
  standard van der Waals radii for atoms
  (in angstrom)

  EV, March 30, 2001
 ------------------------------------------*/

#ifndef _psi_include_vdwradii_h_
#define _psi_include_vdwradii_h_

#define LAST_VDW_RADII_INDEX 9

#ifdef __cplusplus
extern "C" {
#endif

double atomic_vdw_radii[] = { 2.0,  /* default element or ghost */
                    1.2,  /* hydrogen */
		    1.4,  /* helium */
		    1.82, /* lithium */
		    1.8,  /* berrilium (there's no info in literature) */
		    1.6,  /* boron (there's no info in literature) */
		    1.70, /* carbon */
		    1.55, /* nitrogen */
		    1.52, /* oxygen */
		    1.47  /* fluorine */
};

#ifdef __cplusplus
}
#endif

#endif /* header guard */
