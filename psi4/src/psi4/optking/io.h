/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file    binary_io.h
    \ingroup optking
    \brief   header for binary reading and writing functions
*/

#ifndef _opt_io_h_
#define _opt_io_h_

#include "package.h"

#include <cstddef>

namespace opt {

enum OPT_IO_FILE_STATUS {OPT_IO_OPEN_NEW, OPT_IO_OPEN_OLD} ;

bool opt_io_is_present();
void opt_io_remove(bool force=false);
void opt_io_open(OPT_IO_FILE_STATUS status);
void opt_io_close(int keep);
void opt_io_read_entry(const char *key, char *buffer, size_t size);
void opt_io_write_entry(const char *key, char *buffer, size_t size);
void opt_intco_dat_remove();
void opt_clean();

}

#endif
