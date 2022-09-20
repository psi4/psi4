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

#ifndef PSIO_H
#define PSIO_H

#include <cstdio>
#include "psi4/libpsio/config.h"
#include <string>

namespace psi {

std::string decode_errno(const int errno_in);
int psio_init();
void psio_error(size_t unit, size_t errval, std::string prev_msg = "");
int psio_open(size_t unit, int status);
int psio_close(size_t unit, int keep);
std::string psio_getpid();

PSI_API psio_address psio_get_address(psio_address start, size_t shift);
psio_address psio_get_global_address(psio_address entry_start, psio_address rel_address);
int psio_volseek(psio_vol *vol, size_t page, size_t offset, size_t numvols);

int psio_tocwrite(size_t unit);
void psio_tocprint(size_t unit);
psio_tocentry *psio_tocscan(size_t unit, const char *key);
bool psio_tocentry_exists(size_t unit, const char *key);

int psio_write(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd, psio_address *eadd);
int psio_read(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd, psio_address *eadd);
int psio_write_entry(size_t unit, const char *key, char *buffer, size_t size);
int psio_read_entry(size_t unit, const char *key, char *buffer, size_t size);

int psio_open_check(size_t unit);
size_t psio_rd_toclen(size_t unit);
}  // namespace psi

#endif /* #ifndef PSIO_H */
