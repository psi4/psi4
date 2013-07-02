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

#ifndef _psi_include_chkpt_params_h_
#define _psi_include_chkpt_params_h_

// Up to k-functions (L=7) currently supported
const int MAXANGMOM = 11;
const double LINDEP_CUTOFF = 1.0e-6;
//#define MAXANGMOM 11          /* Up to k-functions (L=7) currently supported */
//#define LINDEP_CUTOFF 1E-6    /* Threshold below which basis functions are considered
//                                 linearly dependent (see canonical orthogonalization in Szabo in Ostlund) */

#endif /* header guard */
