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

/*! \file    binary_io.h
    \ingroup optking
    \brief   header for binary reading and writing functions
*/

#ifndef _opt_io_h_
#define _opt_io_h_

#include "package.h"

namespace opt {

typedef unsigned long int ULI;

enum OPT_IO_FILE_STATUS {OPT_IO_OPEN_NEW, OPT_IO_OPEN_OLD} ;

bool opt_io_is_present(void);
void opt_io_remove(void);
void opt_io_open(OPT_IO_FILE_STATUS status);
void opt_io_close(int keep);
void opt_io_read_entry(const char *key, char *buffer, ULI size);
void opt_io_write_entry(const char *key, char *buffer, ULI size);
void opt_intco_dat_remove(void);
void opt_clean(void);

}

#endif