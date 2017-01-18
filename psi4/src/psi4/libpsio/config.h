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

#ifndef _psi_src_lib_libpsio_config_h_
#define _psi_src_lib_libpsio_config_h_

namespace psi {

#define PSIO_OPEN_NEW 0
#define PSIO_OPEN_OLD 1

#define PSIO_KEYLEN 80
#define PSIO_MAXVOL 8
#define PSIO_MAXUNIT 500
#define PSIO_PAGELEN 65536

#define PSIO_ERROR_INIT       1
#define PSIO_ERROR_DONE       2
#define PSIO_ERROR_MAXVOL     3
#define PSIO_ERROR_NOVOLPATH  4
#define PSIO_ERROR_OPEN       5
#define PSIO_ERROR_REOPEN     6
#define PSIO_ERROR_CLOSE      7
#define PSIO_ERROR_RECLOSE    8
#define PSIO_ERROR_OSTAT      9
#define PSIO_ERROR_LSEEK     10
#define PSIO_ERROR_READ      11
#define PSIO_ERROR_WRITE     12
#define PSIO_ERROR_NOTOCENT  13
#define PSIO_ERROR_TOCENTSZ  14
#define PSIO_ERROR_KEYLEN    15
#define PSIO_ERROR_BLKSIZ    16
#define PSIO_ERROR_BLKSTART  17
#define PSIO_ERROR_BLKEND    18
#define PSIO_ERROR_IDENTVOLPATH 19
#define PSIO_ERROR_MAXUNIT   20

typedef unsigned long int ULI; /* For convenience */

typedef struct {
    ULI page; /* First page of entry */
    ULI offset; /* Starting byte offset on fpage */
} psio_address;

typedef struct {
    char *path;
    int stream;
} psio_vol;

typedef struct psio_entry {
    char key[PSIO_KEYLEN];
    psio_address sadd;
    psio_address eadd;
    struct psio_entry *next;
    struct psio_entry *last;
} psio_tocentry;

typedef struct {
    ULI numvols;
    psio_vol vol[PSIO_MAXVOL];
    ULI toclen;
    psio_tocentry *toc;
} psio_ud;

/** A convenient address initialization struct */
extern psio_address PSIO_ZERO;

}

#endif /* header guard */