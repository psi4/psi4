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

/*!
  \file
  \brief Override strncpy to ensure strings are terminated
  \ingroup QT
  By Edward Valeev
*/

#include <cstring>
#include "psi4/libqt/qt.h"

namespace psi {

/*!
** strncpy(): Override default strncpy to ensure last byte is a string
**   terminating character.
**
** \param dest   = destination string
** \param source = source string
** \param n      = number of characters to copy
**
** Returns: pointer to destination string
*/
char* strncpy(char* dest, const char* source, size_t n) {
  if (n > 0) {
    ::strncpy(dest,source,n);
  }
  dest[n-1] = '\0';
  return dest;
}

}
