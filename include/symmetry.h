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

/*-------------------------------------------------------
  symmetry.h : Definitions of various symmetry constants

  Edward Valeev, Oct. 1999
 -------------------------------------------------------*/

#ifndef _psi_include_symmetry_h_
#define _psi_include_symmetry_h_

/*--------------------------------------------------------------
  Point groups at this moment are limited to D2h and sumbgroups
  Therefore the number of symmetry elements is at most 8. Each
  symmetry operation corresponds to a bit of a byte word. The
  correspondence is hardwired via defines "GFLAG" where G is
  the operation symbol.

  To describe nuclear stabilizers or subgroups in general I use
  a byte in which bits corresponding to the symmetry operations
  that constitute the group are set. The result is that each
  operation's contribution to the byte equals "GCODE".
 --------------------------------------------------------------*/
#define EFLAG 0
#define C2ZFLAG 1
#define C2YFLAG 2
#define C2XFLAG 3
#define IFLAG 4
#define SIGXYFLAG 5
#define SIGXZFLAG 6
#define SIGYZFLAG 7
#define C2XCODE 1<<C2XFLAG
#define C2YCODE 1<<C2YFLAG
#define C2ZCODE 1<<C2ZFLAG
#define ICODE 1<<IFLAG
#define SIGXYCODE 1<<SIGXYFLAG
#define SIGXZCODE 1<<SIGXZFLAG
#define SIGYZCODE 1<<SIGYZFLAG
#define ECODE 1<<EFLAG

/*-----------------
  Indices for axes
 -----------------*/
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

#endif /* header guard */
