/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#include <optional>
#include "psi4/libpsio/config.h"
#include <string>

#ifdef _MSC_VER
#include <io.h>
#define SYSTEM_WRITE ::_write
#define SYSTEM_READ ::_read
#define SYSTEM_LSEEK ::_lseeki64
#define SYSTEM_OPEN ::_open
#define SYSTEM_CLOSE ::_close
#define SYSTEM_UNLINK ::_unlink
#define PSIO_OPEN_OLD_FLAGS _O_BINARY | _O_CREAT | _O_RDWR
#define PSIO_OPEN_NEW_FLAGS _O_BINARY | _O_CREAT | _O_RDWR | _O_TRUNC
#define PERMISSION_MODE _S_IWRITE
#else
#include <unistd.h>
#define SYSTEM_WRITE ::write
#define SYSTEM_READ ::read
#define SYSTEM_LSEEK ::lseek
#define SYSTEM_OPEN ::open
#define SYSTEM_CLOSE ::close
#define SYSTEM_UNLINK ::unlink
#define PSIO_OPEN_OLD_FLAGS O_CREAT | O_RDWR
#define PSIO_OPEN_NEW_FLAGS O_CREAT | O_RDWR | O_TRUNC
#define PERMISSION_MODE S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH
#endif

namespace psi {

std::string decode_errno(const int errno_in);
std::string psio_compose_err_msg(const std::string &beginning, const std::string &context, const size_t unit,
                                 const std::optional<int> errno_in = std::nullopt);

int psio_init();
void psio_error(size_t unit, size_t errval, std::string prev_msg = "");
int psio_open(size_t unit, int status);
int psio_close(size_t unit, int keep);
std::string psio_getpid();

PSI_API psio_address psio_get_address(psio_address start, size_t shift);
psio_address psio_get_global_address(psio_address entry_start, psio_address rel_address);
void psio_volseek(const psio_vol *vol, size_t page, const size_t offset, const size_t numvols, const size_t unit);

int psio_tocwrite(size_t unit);
void psio_tocprint(size_t unit);  // debug printing
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
