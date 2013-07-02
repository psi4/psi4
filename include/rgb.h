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
  Default RGB colors for atoms

  EV, March 30, 2001
 ------------------------------------------*/

#ifndef _psi_include_rgb_h_
#define _psi_include_rgb_h_

#define LAST_RGB_INDEX 9

#ifdef __cplusplus
extern "C" {
#endif

double atomic_rgb[][3] = { {0.40, 0.40, 0.40},  /* default element or ghost */
                    {1.00, 1.00, 1.00},  /* hydrogen */
		    {0.80, 0.80, 0.50},  /* helium */
		    {0.30, 0.80, 0.30},  /* lithium */
		    {0.65, 0.80, 0.25},  /* berrilium */
		    {0.50, 0.80, 0.15},  /* boron */
		    {0.25, 0.25, 0.25},  /* carbon */
		    {0.00, 0.00, 1.00},  /* nitrogen */
		    {1.00, 0.00, 0.00},  /* oxygen */
		    {0.75, 0.40, 0.15}  /* fluorine */
};

#ifdef __cplusplus
}
#endif

#endif /* header guard */
