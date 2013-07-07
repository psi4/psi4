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

#ifndef _DEFINES_H
#define	_DEFINES_H

// Set either MPQC_INTERFACE or PSI4_INTERFACE to 1
#define MPQC_INTERFACE 1

#if MPQC_INTERFACE
   #define OPEN_NAMESPACE namespace tiled_tensor{
   #define CLOSE_NAMESPACE }
#elif PSI4_INTERFACE
   #define OUT0 ExEnv::out0()
   #define OUT0 ExEnv::out0()
#else
   #error "You must select an interface in defines.h"
#endif

#define harikari(msg) std::cerr << msg << endl; abort()


#endif	/* _DEFINES_H */

